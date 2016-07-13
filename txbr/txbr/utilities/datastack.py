import os.path
import glob
import re
import tempfile
import configobj
import scipy.ndimage.interpolation
import numpy
import mrc
import txbr.prefil
import txbr.utilities
import util

from txbr.txbrconf import log

global TIFF
global JPG

TIFF = "tif"
JPG = "jpg"

def loadStack( directory, item, workDirectory=".", scope=None, bin=None ):
    """
    *item*: either a MRCFile or a basename
    """

    if not os.path.lexists(directory):
        log.warning("%s does not exist!" %directory)
        return None

    log.debug("loadstack:  d: %s" %(directory))
    log.debug("loadstack: wd: %s" %(workDirectory))

    # MRC stack cases

    mrc_cases = [ item, \
                  item + txbr.utilities.ST_EXTENSION, \
                  item + txbr.utilities.MRC_EXTENSION, \
                  item + txbr.utilities.PREALI_EXTENSION ]

    for case in mrc_cases:
        if os.path.lexists(os.path.join(directory,case)):
             return MRCStack( case, directory, workDirectory=workDirectory, scope=scope, bin=bin )

    # Check tiff image stack cases. A metadata file with extension ".mdoc" must be detected in
    # the data directory
    
    img_cases = [ item, \
                  item + txbr.utilities.MDOC_EXTENSION ]

    for case in img_cases:
        if os.path.lexists(os.path.join(directory,case)):
            return ImageStack( item, directory, workDirectory=workDirectory, scope=scope, bin=bin )

    return None


class Stack:
    '''This class is used to describe stack of images that are used as inputs
    in a tomographic series.'''

    def __init__( self ):
        '''Initialize some variables of the stack'''

        self.basename = None
        self.angles = []
        self.subset = []
        self.workDirectory = "."


    def load( self, verbose=True ):

        raise NotImplementedError


    def getNumberOfImages( self ):
        '''Returns the number of images within the stack.'''

        return len(self.subset)


    def isEmpty(self):

        return self.getNumberOfImages()==0


    def getImageAt( self, index ):
        '''Returns a numpy array for an image at a given index within the stack.'''

        raise NotImplementedError


    def getImageAtAngle( self, angle ):
        '''Returns a numpy array for an image at a given angle within the stack.'''

        return self.getImageAt(self.indexOfAngle(angle))


    def getImageScale( self ):
        '''Returns the image scale.'''

        raise NotImplementedError


    def getImageSize( self ):
        '''Returns an image size.'''

        return NotImplementedError


    def indexOfAngle( self, angle ):

        index = self.subset.index(self.angles.index(angle))

        return index


    def getTiltAngles( self, unit='degree' ):
        '''Returns a list of the tilt angles available in that stack within the stack.'''

        angles = [ self.angles[index] for index in self.subset ]

        if unit=='degree':
            return angles
        else:
            return numpy.radians(angles)


    def getTiltAngleAt( self, index, unit='degree' ):
        '''Returns the tilt angle for a given index within the stack.'''

        index = self.subset[index]

        if unit=='degree':
            return self.angles[index]
        else:
            return numpy.radians(self.angles[index])


    def getReferenceIndex(self):
        '''Return the index corresponding to having the lowest tilt as possible'''

        #print "The reference index is: %i" %(numpy.argmin(numpy.abs(self.angles)))

        angles = [ self.angles[index] for index in self.subset ]

        return numpy.argmin(numpy.abs(angles))


    def getReferenceAngle(self):
        '''Return the angle corresponding to having the lowest tilt as possible'''

        return self.angles[self.getReferenceIndex()]


    def getReferenceImage(self):
        '''Return the image corresponding to having the lowest tilt as possible'''

        return self.getImageAt(self.getReferenceIndex())


    def filter( self, filt_angle, index2filter=None ):
        '''Filter the data'''

        return NotImplementedError


    def getFilteredImageAt( self, index ):

        return NotImplementedError


    def getFilteredImageAtAngle( self, angle ):
        '''Returns a numpy array for a filtered image at a given angle within the stack.'''

        return self.getFilteredImageAt(self.indexOfAngle(angle))


    def showImageAt( self, index ):
        '''Display an image of the stack at a given index as an openCV image'''

        u = self.getImageAt(index)
        u = (u-numpy.min(u))/numpy.max(u)

        util.showArray(u)


class MRCStack(Stack):
    '''A stack of images collected from an MRC file'''

    def __init__( self, filename, directory=".", workDirectory=".", scope=None, bin=1 ):

        self.filename = os.path.join(directory,filename )
        self.basename = os.path.splitext(filename)[0]
        self.directory = directory
        self.workDirectory = workDirectory
        self.bin = bin

        if self.bin==None: self.bin = 1
      
        self.mrc_file = mrc.MRCFile( self.filename )

        try:
            self.angles = mrc.tilts( self.filename )
        except:
            print "Could not load the tilt angles from the MRC stack header"
            pass

        if os.path.exists( "%s.rawtlt" %self.basename ):
            self.angles = numpy.loadtxt("%s.rawtlt" %self.basename).tolist()

        if len(self.angles)!=self.mrc_file.nz:
            print "The number of angles does not match the number of slices"

        self.subset = range(self.mrc_file.nz)   # subset of the images

        self.setBinFactor( bin )    # Set the binning and load the metadata information
        

    def setBinFactor( self, bin ):

        if bin==None: bin=1

        log.info("Set the bin factor of %s to %i" %(self,bin))

        self.bin = bin

        if self.bin==1:
            self.data_dir = self.directory    # If there is no binning the data should be pulled from the src directory
            self.filt_dir = os.path.join(self.workDirectory,"txbr-filter")
        else:
            self.data_dir = os.path.join(self.workDirectory,"data","bin%i" %self.bin)
            self.filt_dir = os.path.join(self.workDirectory,"txbr-filter","bin%i" %self.bin)

        self.load()


    def load( self ):

        self.slice4angle = {}

        for angle in self.angles:
            self.slice4angle[angle] = self.filename

        # Do binning if necessary

        if self.bin==1: return

        head,tail = os.path.split(self.filename)
        file_binned = os.path.join(self.data_dir,tail)

        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)

        if not os.path.exists(file_binned):
            os.system("newstack -bin %i %s %s" %(self.bin,self.filename,file_binned))

        self.mrc_file = mrc.MRCFile(file_binned)


    def getImageAt( self, index ):
        '''Returns a numpy array for an image at a given index within the stack.'''

        index = self.subset[index]
        u = self.mrc_file.getZSliceAt(index)

        return u.astype('float32')  # Otherwise openCV does not deal well with warping images


    def getImageSize( self ):
        '''Returns an image size.'''

        return (self.mrc_file.nx,self.mrc_file.ny)


    def getImageScale( self ):
        '''Returns the scale of the raw image stack.'''

        return (self.mrc_file.sx,self.mrc_file.sy)


    def filter( self, filt_angle, index2filter=None ):
        '''Filter the data'''

        if len(self.subset)==0: return

        # Create a temporary input file

        tempfilename = tempfile.mktemp()

        os.system("newstack -secs %s %s %s" %(",".join([str(s) for s in self.subset]), self.filename, tempfilename))

        # Output file

        file_out = "%s.SL" %( self.basename )

        # Other parameters

        filter_magnification = 1

        # Apply the filter

        txbr.prefil.rot_filter( tempfilename, file_out, filter_magnification, filt_angle)

        # Clean up

        os.remove(tempfilename)


    def getFilteredImageAt( self, index ):

        return mrc.MRCFile("%s.SL" %self.basename).getZSliceAt(index)


    def __repr__( self ):

        return "<%s: %s, bin=%i>" %( self.__class__, self.filename, self.bin )


    def __str__( self ):

        return "%s: %s, bin=%i" %( self.__class__, self.filename, self.bin )


class ImageStack(Stack):
    '''A stack of images as a collection of image stack using regular format (tiff,jpeg)'''

    def __init__( self , basename, directory, workDirectory=".", extension=None, bin=None, scope=None):

        self.directory = directory
        self.workDirectory = workDirectory
        self.basename = basename
        self.extension = extension
        
        try:
            self.bin = int(bin)
        except:
            self.bin = 1

        self.scope = scope

        mdoc_file = os.path.join(self.directory,"%s.mdoc" %self.basename)

        if os.path.exists(mdoc_file):
            self.metadata = configobj.ConfigObj(mdoc_file)
        else:
            self.metadata = None
        
        self.angles = []
        self.subset = []

        self.slice4angle = {}

        jpg_dir = os.path.join(self.directory,"src_jpegs")
        wib_dir = os.path.join(self.directory,"WIB")

        if self.extension==None and os.path.exists(wib_dir):
            self.extension=TIFF
        elif self.extension==None and os.path.exists(jpg_dir):
            self.extension=JPG

        self.src_dir = "."  # directory where the rawdata is

        if self.extension==JPG and os.path.exists(jpg_dir):
            self.src_dir = jpg_dir
        if self.extension==TIFF and os.path.exists(wib_dir):
            self.src_dir = wib_dir

        if self.bin==1:
            self.data_dir = self.src_dir    # If there is no binning the data should be pulled from the src directory
            self.filt_dir = os.path.join(self.workDirectory,"txbr-filter")
        else:
            self.data_dir = os.path.join(self.workDirectory,"data","bin%i" %self.bin)
            self.filt_dir = os.path.join(self.workDirectory,"txbr-filter","bin%i" %self.bin)

        self.load() # Load images from the hard drive


    def scanSourceDirectory( self, verbose=True ):

        # The names of the file to process is a concatenation of a basename and the tilt angle

        pattern = os.path.join( self.src_dir, '*%s*.%s' %(self.basename,self.extension))

        if verbose: print "Search for files like %s in %s" %(pattern,self.src_dir)

        slices = glob.glob( pattern )

        return slices
    

    def load( self, clear=False, verbose=True ):
        '''Scan the source directory for possible images and load some associated metadata (angles)
        either from the mdoc or the name of the file. This routine fills out the dictionary self.slice4angle
        which maps a tilt angle to an image (path) in either case. The two lists *self.angles* and *self.subset*
        are regenerated afterwards.'''

        slices = self.scanSourceDirectory( verbose )

        self.slice4angle = {}

        # Extract the tilt angles and map them to the images by creating the
        # dictionary self.slice4angle

        for tilt in slices:
            id_match = re.match('\S*_([+-.\d]*).%s' %self.extension, tilt)
            if self.metadata!=None and id_match:
                index = int(id_match.group(1))
                angle = float(self.metadata["ZValue = %i" %index]['TiltAngle'])
                self.slice4angle[angle] = tilt
                if verbose: print "Found file for angle %6s: %s" %(angle,tilt)
            elif id_match:  # The matching pattern corresponds to an angle
                angle = float(id_match.group(1))
                self.slice4angle[angle] = tilt
                if verbose: print "Found file for angle %6s: %s" %(angle,tilt)

        self.angles = sorted(self.slice4angle.keys())
        self.subset = range(len(self.angles))

        # Do binning if necessary

        if self.bin==1: return

        if not os.path.exists(self.data_dir): os.makedirs(self.data_dir)

        for angle,imgpath in self.slice4angle.iteritems():

            head,tail = os.path.split(imgpath)
            binned_output = os.path.join(self.data_dir,tail)

            if clear or not os.path.lexists(binned_output): 
                u = scipy.ndimage.interpolation.zoom(util.loadImage(imgpath),1.0/self.bin)
                util.saveImage( binned_output, u )

            self.slice4angle[angle] = binned_output


    def getImageAt( self, index, verbose=True ):
        '''Returns a numpy array for an image at a given index within the stack.'''

        index = self.subset[index]
        angle = self.angles[index]

        nameOfImage = self.slice4angle[angle]

        if verbose: print "Load Image: %s" %(nameOfImage)

        return util.loadImage(nameOfImage)


    def getImageScale( self ):
        '''Returns the image scale.'''

        if self.metadata!=None:
            s = float(self.metadata["PixelSpacing"])*self.bin
            return (s,s)
        else:
            return None
    

    def getImageSize( self ):
        '''Returns an image size.'''

        if self.metadata!=None:
            imgsize = self.metadata["ImageSize"].split()
            return [int(l)/self.bin for l in imgsize]

        if len(self.subset)==0 or len(self.angles)==0:
            return None
        
        return self.getImageAt(0).shape


    def showImageAt( self, index ):
        '''Display an image of the stack at a given index as an openCV image.
        This function is overriden from the parent class Stack for orientation reason (to
        match the MRC orientation)'''

        u = self.getImageAt(index)
        u = u[:,::-1].T.copy()
        u = (u-numpy.min(u))/numpy.max(u)

        util.showArray(u)


    def filter( self, filt_angle, index2filter=None ):
        '''Filter the data'''

        if index2filter==None: index2filter = range(len(self.subset))

        if len(self.subset)==0 or len(index2filter)==0: return

        if not os.path.lexists(self.filt_dir):
            os.makedirs(self.filt_dir)
        
        nx,ny = self.getImageSize()

        for index in index2filter:

            tempfilename = tempfile.mktemp()    # Create a temporary input file

            f = mrc.MRCFile(tempfilename)
            f.setHeader(nx,ny,1)
            f.setZSliceAt(0,self.getImageAt(index))
            f.updateHeader()

            # Output file

            head,tail = os.path.split(self.slice4angle[self.getTiltAngleAt(index)])
            file_out = os.path.join(self.filt_dir,"%s.SL" %tail)

            # Other parameters

            filter_magnification = 1

            # Apply the filter

            txbr.prefil.rot_filter( tempfilename, file_out, filter_magnification, filt_angle)

            # Clean up

            os.remove(tempfilename)


    def getFilteredImageAt( self, index ):

        head,tail = os.path.split(self.slice4angle[self.getTiltAngleAt(index)])
        filt_file = os.path.join(self.filt_dir,"%s.SL" %tail)

        f = mrc.MRCFile(filt_file)
       
        return f.getZSliceAt(0)


if __name__ == '__main__':

    options1 = {'directory':"/ncmirdata3/sphan/otf/ju", 'item':"Ju_cont2.mrc"}
    options2 = {'directory':"/ncmirdata3/sphan/otf/ju2", 'item':"CCDBid_74225_Ju_cont2.mrc"}

    stack = loadStack(**options1)
    stack.showImageAt(10)
