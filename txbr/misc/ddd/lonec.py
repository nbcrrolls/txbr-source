import sys
import os
import os.path
import getopt
import re
import numpy
import scipy
import scipy.misc
import scipy.ndimage.filters
import scipy.stats
import cv
import mrc


def usage():

    print
    print 'Usage: %s.py -d directory [Options]' %os.path.basename(sys.argv[0])
    print
    print "    --dark value [value in ('direct','evaluate')]"
    print "        Mode for evaluating the dark current (default is direct)"
    print "    --test"
    print "        Some testing"
    print "    --threshold value"
    print "        Set the threshold value for the noise in each image (default 5)"
    print "    --window size (default 3)"
    print "        The size of the window used to find the local maximum"
    print "    --DoPlot"
    print "        Some plotting"


def show( image ):
    '''Show an image with the openCV package'''

    cv.NamedWindow("Output")
    cv.ShowImage("output", image)
    cv.WaitKey()


def local_maxima( array2d, window_size=3 ):
    '''Find local maxima of a 2d array.
    Use the maximum filter with a window size of 3 for this purpose
    '''

    array2d_m = scipy.ndimage.filters.maximum_filter(array2d,window_size)

    return ((array2d==array2d_m) & (array2d!=0))


def connectivity( array2d ):

    array_int = numpy.roll(array2d,  1, 0) + array2d + numpy.roll(array2d,  -1, 0)

    array_int =  numpy.roll(array_int,  1, 1) + array_int + numpy.roll(array_int,  -1, 1)

    return array2d*array_int
    
    
def sharpen( u , noise_threshold = 5, window_size=3, test=False):

    doSharpen = True
    doGaussianBlur = True
    sigma = 2.0

    global sharpening_index

    if not doSharpen: return numpy.dstack((numpy.zeros_like(u),numpy.zeros_like(u),u))
    if not doSharpen: return u

    arr =  numpy.where(u<=noise_threshold,0.0,u)

    if doGaussianBlur: arr =  numpy.where(arr!=0,scipy.ndimage.gaussian_filter(u,sigma),0.0)

    final = numpy.where(local_maxima(arr,window_size=window_size),1.0,0.0)

    ne = numpy.sum(final)
    density = float(ne)/final.size
    print "%i: ne: %i    %f electron per pixel [ mean dist=%.1f - window size=%i]" %(sharpening_index,ne,density,1/density,window_size)

    sharpening_index = sharpening_index + 1

    if test:
        final = numpy.dstack((final*255,arr,u))
        #print numpy.min(numpy.min(final,axis=0),axis=0)
        return final
    else:
        return final


def rmOutliers():

    basename = "output"

    file_in = "%s.tif" %basename
    file_out = "%s.clean.tif" %basename

    u = scipy.misc.imread(file_in)

    score = 5.0*numpy.mean(u)

    print "Mean: %f" %numpy.mean(u)
    print "Mean: %f" %float(numpy.sum(u)/u.size)
    print scipy.stats.percentileofscore(numpy.ravel(u),score)

    u = scipy.stats.threshold(u,threshmax=score)

    print "Mean: %f" %numpy.mean(u)

    scipy.misc.imsave(file_out, u)



def getDarkImage( directory=None, files=None ):
    '''Get the dark current image.
    If a list of files containing low electron count images, the dark current image is
    evaluated by substracting the first two image contributions.
    Otherwise, check in the folder 'directory' for the reference image called "DarkReference.tif"
    '''

    if directory!=None:
	 dark_file = os.path.join(directory,"DarkReference.tif")

    if files!=None and len(files)>1:
        u1 = numpy.asarray(scipy.misc.imread(files[0]),dtype="float")
        u2 = numpy.asarray(scipy.misc.imread(files[1]),dtype="float")
        dark = u1 - numpy.where(u1-u2>=0,u1-u2,0.0)
    elif directory!=None and os.path.exists(dark_file):
        dark = numpy.asarray(scipy.misc.imread(dark_file),dtype="float")
    else:
        sys.exit("No dark image!")

    return dark
    

def main():

    '''Main routine for this flattening module'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hd:",[ 'dark=','threshold=','window=','test','doPlot','help','remove-outliers' ])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    global sharpening_index

    sharpening_index = 0

    directory = "."
    basename = "Image"
    dark = 'direct'
    test = False
    index0 = 0
    step = 50   # Dump intermediate output every 'step' images
    doPlot = False
    use_pil = True  # Use PIL module to save the image
    noise_threshold = 5
    window_size = 3
    removeOutliers = False

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-d"):
            directory = value_
        elif option_ in ("--dark"):
            dark = value_
        elif option_ in ("--threshold"):
            noise_threshold = int(value_)
        elif option_ in ("--window"):
            window_size = int(value_)
        elif option_ in ("--test"):
            test = True
        elif option_ in ("--doPlot"):
            doPlot = True
        elif option_ in ("--remove-outliers"):
            removeOutliers = True
        else:
            assert False, "unhandled option"

    print removeOutliers

    if removeOutliers:
        rmOutliers()
        sys.exit()

    doPlot = doPlot or test

    if not os.path.exists(directory):
        sys.exit("Please provide a directory!")

    files = os.listdir(directory)

    sequences = []
    for file in files:
        match = re.match('%s(\d*).tif$' %basename, file)
        if not match: continue
        sequences.append(int(match.group(1)))

    sequences.sort()

    if test: sequences = [ sequences[index0] ]

    files = [ os.path.join(directory,"%s%i.tif" %(basename,index)) for index in sequences if os.path.exists(os.path.join(directory,"%s%i.tif" %(basename,index))) ]

    n = len(files)
    print "There are %i image files." %(n)
    print "Window size: %i" %window_size
    print "Noise threshold: %i" %noise_threshold

    if n==0: sys.exit("There is no image file to process!")

    u0 = numpy.asarray(scipy.misc.imread(files[0]),dtype="float")
    nx,ny = u0.shape

    if dark=='direct':
        dark = getDarkImage( directory=directory )
    elif dark=='evaluate':
        dark = getDarkImage( files=files )
    else:
        dark = numpy.zeros((nx,ny))

    u = numpy.zeros_like(dark)

    if test:
        u = numpy.dstack((u,u,u))

    for index,file in enumerate(files[:]):
        slice = numpy.asarray(scipy.misc.imread(file),dtype="float") - dark
        u = u + sharpen(slice,noise_threshold=noise_threshold,window_size=window_size,test=test)
        if index%step==0:
            scipy.misc.imsave("output%i.tif" %index, u)

    ne = numpy.sum(u)
    density = float(ne)/u.size
    print "Final output: ne: %i    %f electron per pixel" %(ne,density)

    if use_pil:
        file_out = "output.tif"
        scipy.misc.imsave(file_out, u)
        if doPlot: os.system("/home/sphan/work/softwares/ImageJ/run %s &" %file_out)
    else:
        file_out = "output.mrc"
        u = u[::-1,:].T
        fout = mrc.MRCFile(file_out)
        fout.setHeader(u.shape[0],u.shape[1],1)
        fout.setZSliceAt(0,u)
        fout.updateHeader()
        if doPlot: os.system("imod %s" %file_out)



if __name__ == '__main__':

   main()
  
