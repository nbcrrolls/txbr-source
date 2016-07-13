"""
.. `detectmarker-label`:

The module *detectmarker* allows to detect markers in an Electron Tomography tilt
series.

"""

import os.path
import math
import numpy
import numpy.fft
import scipy.ndimage.filters
import scipy.ndimage.measurements
import mrc
import modl
import util
import otfutil

import scipy.optimize

from txbr.txbrconf import log


def smooth(x,window_len=31,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y


def highpass(image, sigma=20):

    f = numpy.fft.rfft(image)
    g = scipy.ndimage.fourier_gaussian(f, sigma )
    h = numpy.fft.irfft(g)
    return h.real

def smooth_( x, y, order ):

    from util import polynomial_regression

    n = x.size
    x_ = numpy.resize(x,(n,1))

    y_shape= y.shape

    if y_shape[0]!=n: raise ValueError

    if len(y_shape)==1:
        y_ = numpy.resize(y,(n,1))
    else:
        y_ = numpy.resize(y,(n,y.size/n))

    order_ = numpy.array([order])

    smooth_y = numpy.empty_like(y_)

    for i in range(y.size/n):
        p = polynomial_regression(x_,y_[:,i],order_)
        smooth_y[:,i] = p.eval(x_)

    smooth_y.resize(y_shape)

    return smooth_y


def dumpMarkerPatchesIntoMRCFile( filename, u, doShow=True ):
    """This function dumps Marker patches into an MRC file (for checking purpose).

    :param filename: The name of the output MRC file.
    :param u: A list of *numpy* array of shape (nimg,n1,n2) containing the marker patches
        data. The length of the list corresponds to the total number of slices in the input
        MRC file. The variables *nimg*, *n1*, *n2* are respectively the number of patches,
        their width and height; *nimg* can be different for each microgaphs.

    """

    u = numpy.asarray(u)

    n,nimg,n1,n2 = u.shape

    nn = numpy.ceil(numpy.sqrt(n)).astype('int')

    img = numpy.zeros((nn*n1,nn*n2))

    f = mrc.MRCFile(filename)
    f.setHeader(nn*n1,nn*n2,nimg)

    for index1 in range(nimg):
        for index2 in range(n):
            x1 =(index2%nn)*n1
            x2 = x1 + n1
            y1 = (index2/nn)*n2
            y2 = y1 + n2
            slice = u[index2,index1]
            slice = (slice - numpy.mean(slice))/numpy.std(slice)
            img[x1:x2,y1:y2] = slice
        f.setZSliceAt(index1,img)

    f.updateHeader()

    if doShow: os.system("imod {0:s}".format(filename))


def crossCheck( uref, u, remove_background=False, verbose=True ):
    '''Cross-check two image patches.

    :param uref: The image reference (a numpy array)
    :param u: The image to compare with the reference (a numpy array)
    :param verbose: If true, the function outputs extra information.
    :returns: True if the two images are looking alike.
    '''

#    uref = numpy.asarray(uref)
#    u = numpy.asarray(u)
#
#    devref = numpy.max(uref) - numpy.min(uref)
#    dev = numpy.max(u) - numpy.min(u)
#
#    accept =  dev>0.25*devref
#
#    print " Accept: %s     %.2f / %.2f" %(accept, dev, 0.5*devref)

    if remove_background:
        uref = otfutil.removeBackground(uref)
        u = otfutil.removeBackground(u)

    uref = scipy.ndimage.filters.median_filter(uref,3)
    u = scipy.ndimage.filters.median_filter(u,3)

    d = numpy.sum(uref*u)/numpy.linalg.norm(uref)/numpy.linalg.norm(u)

    threshold = 0.95
    accept = d>threshold

    if verbose: log.info(" Accept: {0:s}     {1:.2f} / {2:.2f}".format(accept,d,threshold))

    return accept



    if u.shape!=uref.shape: return False

    CORR_RATIO_THRESHOLD = 0.45

    u1 = numpy.fft.fft2(uref)
    u2 = numpy.fft.fft2(u)

    phase =  u1*numpy.conj(u2)/numpy.abs( u1*u2 )
    peak = numpy.fft.ifft2(phase)

    peak_value = numpy.max(peak.real)

    accept = (peak_value>=CORR_RATIO_THRESHOLD)

    if verbose: log.info(" Accept: {0:s}     {1:.2f} / {1:.2f}".format(accept,peak_value,CORR_RATIO_THRESHOLD))

    return accept


def validatePeak( u, mean=None, verbose=True ):
    '''Validate a peak in an image..

    :param u: The peak image to valide (a numpy array)
    :param verbose: If true, the function outputs extra information.
    :returns: True if the two images are looking alike.
    '''


    log.info("min: {0:.1f}  max: {1:.1f}  mean: {2:.1f}".format(numpy.min(u),numpy.max(u),numpy.mean(u)))

    return True

  #  u = scipy.ndimage.filters.median_filter(u,3)

    SN_RATIO_THRESHOLD = 50.0
    SN_RATIO_THRESHOLD = 0.03

#    accept = math.fabs((numpy.max(u)-numpy.min(u))/numpy.std(u))>=SN_RATIO_THRESHOLD
    accept = math.fabs((numpy.max(u)-numpy.min(u))/numpy.mean(u))>=SN_RATIO_THRESHOLD
    accept = math.fabs(numpy.std(u)/mean)>=SN_RATIO_THRESHOLD
    accept = math.fabs(numpy.std(u))>=SN_RATIO_THRESHOLD
   # accept = numpy.max(u)-mean>

    if verbose:
        log.info("Accept: {0:s}  {1:.2f}  {2:.2e}".format(accept,math.fabs(numpy.std(u)),mean))

    return accept


def createHessianKernel( nx, ny, sigma=0.25 ):
    '''Creates the kernel images (to implement gaussian derivatives) to detect beads
    of a certain size.

    :param nx,ny: The width and heigh of the kernel image. Should be odd numbers;
        if not, the size is increased by one unit.
    :param sigma: Characteristic length of the kernel. Should be 1/5 of the bead
        diameter in practice.
    :returns: A tuple containing the 3 image kernel :math:`G_xx`, :math:`G_yy`
        and :math:`G_xy`
    '''

    log.info("Create Kernel for the detector image")

    # Make the box odd, so center corresponds to a discret point.
    
    if nx%2==0: nx+=1
    if ny%2==0: ny+=1

    origin = ( (nx-1)/2, (ny-1)/2 )
    y,x = numpy.meshgrid(range(ny),range(nx))

    x = x - origin[0]
    y = y - origin[1]

    r2 = x**2 + y**2

    Gxx = (x**2-sigma**2)/2.0/numpy.pi/sigma**4*numpy.exp(-r2/2.0/sigma**2)
    Gyy = (y**2-sigma**2)/2.0/numpy.pi/sigma**4*numpy.exp(-r2/2.0/sigma**2)
    Gxy = x*y/2.0/numpy.pi/sigma**4*numpy.exp(-r2/2.0/sigma**2)

    return Gxx, Gyy, Gxy


def findMarkersInStack( filename, slices=None, nmax=10000, size=5.0 ):
    '''Find markers in a stack of 2D images, currently a MRC file called *filename*.

    :param filename: The name of the MRC file to process (string).
    :param slices: A list containing the slices index to consider in the stack. If None, all the
        slices will be considered. Numbering convention; 0 corresponds to the first slice
        of the stack.
    :param nmax: The maximum number of beads to be detected (integer)
    :param size: The diameter (pixel size) of the markers (float).
    :returns: A list containing a numpy array (of shape (n,3)) for each of the slices in the
        stack where the detection occured
    '''

    if not os.path.exists(filename):
        log.warning("File {0:s} does not exists!".format(filename))
        return

    tokens = filename.split(".")
    if len(tokens)>1:
        basename = ".".join( token for token in tokens[:-1])
    else:
        basename = filename

    fin = mrc.MRCFile( filename )

    if slices==None:
        slices = range(fin.nz)

    Pts = []
    D = []

    for index,slice in enumerate(slices):
        points, detector, patches = findMarkersInImage( fin.getZSliceAt(slice), nmax=nmax, size=size)
        points = numpy.asarray(points)
        Z = numpy.ones((points.shape[0],1))*slice
        points = numpy.append(points,Z,axis=1)
        Pts.append(points)
        D.append(detector)
        showSignature = index==numpy.floor((len(slices)-1)/2.0)
        if showSignature:
            markerSignatureBasename = "{0:s}.bead.{1:}.mrc".format(basename,slice)
            dumpMarkerPatchesIntoMRCFile( markerSignatureBasename, patches, doShow=showSignature )

    mrc.createMRCFile( "{0:s}.detect.mrc".format(basename), numpy.dstack(D), view=True, scale=(fin.sx,fin.sy,fin.sz) )
    modl.saveAs3DModel( "{0:s}.bead.mod".format(basename), Pts, fin.nx, fin.ny, fin.nz, sx=fin.sx, sy=fin.sy, sz=fin.sz, doShow=True, mrcFileName=filename )

    return Pts


def findMarkersInImage( img, nmin = 0, nmax=10000, size=5, doValidateOnD=False, threshold=0.4, saveAllPointsInFile=None ):
    '''Find beads in a 2D image
    
    :param img: A numpy array containing images.
    
    '''
    
    # Detection mode in [ 'min', 'max', 'mix_mult']

    detection_mode = 'mix_mult'
  #  detection_mode = 'edge'

    # Other ossible options

    apply_gaussian_filter,rsmooth = True,1.0   # Gaussian filter on detection image
    apply_flattening,flat_order = True,1    # Flattening filter on detection image

    skipEdge = False
    secure_ratio = 0.8       # To avoid dealing with close peaks
    apply_cdm = size>10      # CDM correction for the big beads
    apply_cdm = True
    refine = 'micrograph'    # in [ 'detection', 'micrograph' ]

    # Prepare the template for detecting the peaks

    nx,ny = img.shape
    img = util.flattenImage(img, order=2)
    u = (img-numpy.mean(img))/numpy.std(img)    # work on normalized data

    lx,ly = int(size),int(size)
    Gxx,Gyy,Gxy = createHessianKernel( lx, ly, sigma=size/5.0 )

    log.info("Construct detection image...")

    Ixx = scipy.ndimage.filters.convolve( u, Gxx )
    Iyy = scipy.ndimage.filters.convolve( u, Gyy )
    Ixy = scipy.ndimage.filters.convolve( u, Gxy )

    beta = (Ixx+Iyy)/2.0
    delta = beta**2 - (Ixx*Iyy-Ixy**2)
    lambda1 = beta - numpy.sqrt(delta)
    lambda2 = beta + numpy.sqrt(delta)

    lambda_min = numpy.where(lambda1<lambda2,lambda1,lambda2)
    lambda_max = numpy.where(lambda1>lambda2,lambda1,lambda2)

    log.info("Detection image generated...")

    detections = {}

    detections['min'] = lambda_min
    detections['max'] = lambda_max
   # detections['mix_mult'] = lambda_min*lambda_max
   # detections['mix_mult'] = lambda_min*lambda_max*numpy.exp(-numpy.abs(lambda_max-lambda_min))
    detections['mix_mult'] = numpy.sign(lambda_min)*lambda_min*lambda_max*numpy.exp(-numpy.abs(lambda_max-lambda_min))
    detections['edge'] = numpy.abs(lambda_max/lambda_min)

    D = detections[detection_mode]
 
    if apply_flattening: D = util.flattenImage( D, order=flat_order )
    if apply_gaussian_filter:  D = scipy.ndimage.filters.gaussian_filter( D, rsmooth )

    u = D

    # Extract the peaks...

    m1 = numpy.min(u)
    m2 = numpy.max(u)
    m3 = numpy.mean(u)
    m4 = numpy.std(u)

    threshold = m3 + 3.5*m4
    log.info("min: {0:.1f}   threshold: {1:.1f}    max: {2:.1f}".format(numpy.min(u),threshold,numpy.max(u)))

    util.buildHistogram( u, doPlot=False  )

    sorted_ind = u.argsort( axis=None )

    mask = numpy.zeros((nx,ny)) # A mask eliminate closed points

    points = []
    patches = []

    l = int(numpy.ceil(secure_ratio*size))

    uref1 = None
    uref2 = None

    n = 0
    for index,t in enumerate(sorted_ind[-1::-1]):
        i,j = numpy.unravel_index(t,(nx,ny))
        if mask[i,j]==0:  # accepted
            x1 = max(0,i-l)
            x2 = min(nx,i+l)
            y1 = max(0,j-l)
            y2 = min(ny,j+l)
            borderEdge = x1==0 or y1==0 or x2==nx or y2==ny
            if skipEdge and borderEdge:
                mask[x1:x2,y1:y2] = 1
                continue
            if borderEdge:
                uc1 = numpy.zeros((2*l,2*l))
                uc2 = numpy.zeros((2*l,2*l))
                lx1 = max(0,l-i)
                lx2 = min(2*l,nx+l-i)
                ly1 = max(0,l-j)
                ly2 = min(2*l,ny+l-j)
                uc1[lx1:lx2,ly1:ly2] = u[x1:x2,y1:y2]
                uc2[lx1:lx2,ly1:ly2] = img[x1:x2,y1:y2]
            else:
                uc1 = u[x1:x2,y1:y2]
                uc2 = img[x1:x2,y1:y2]
            if uref1==None:
                uref1 = uc1
            if uref2==None:
                uref2 = uc2
         #   validatePeak(uc2)
         #   if not validatePeak(uc1,mean=m3): continue
         #   if not crossCheck( uref1, uc1 ): continue
         #   if not crossCheck( uref2, uc2, remove_background=True ): continue
         #   if not crossCheck( uref1, uc1 ) or not crossCheck( uref2, uc2, remove_background=True ): continue
            n = n + 1
            if borderEdge:
                uc1[0,:] = -0.5
                uc1[-1,:] = -0.5
                uc1[:,0] = -0.5
                uc1[:,-1] = -0.5
            patches.append((uc1,uc2))
            if apply_cdm:
                cdm_options = { 'detection':img[x1:x2,y1:y2], 'micrograph':-u[x1:x2,y1:y2] }
                patch = cdm_options[refine]
                patch = scipy.ndimage.filters.gaussian_filter( patch, rsmooth )
                patch = - otfutil.removeBackground(patch)
                tx,ty = scipy.ndimage.measurements.center_of_mass(patch)
                log.info("Residual Translation: (tx,ty)=({0:1f},{1:1f}) from patch of shape ({2:1f},{3:1f})".format(tx-(x2-x1)/2.0,ty-(y2-y1)/2.0,x2-x1,y2-y1))
             #   points.append([x1+tx,y1+ty])
                points.append([x1+tx+1,y1+ty+1])
            else:
            #    points.append([i,j])
                points.append([i+1,j+1])
            mask[x1:x2,y1:y2] = 1
            percent = (u[i,j]-m1)/(m2-m1)
            percent = 1.0*index/nx/ny
            alpha = 1.0*index/(m2-mask[i,j])
            log.info("#{0:}/{1:}: (i,j)=({2:4},{3:4})   value={4:.2f}   {5:.4f}  {6:.4f}   [min,mean,max]=[{7:.2f},{8:.2f},{9:.2f}]".format(n,index,i,j,u[i,j],percent,alpha,m1,m3,m2))
        if u[i,j]-threshold<=0:
            log.info("Reached threshold value {0:.1f} (>={1:.1f}) -- bead number: {2:}".format(u[i,j],threshold,n ))
            break
        if n>=nmax:
            log.info("Maximum number of beads reached... -- threshold: {0:.3f} > {1:.3f}".format( u[i,j], threshold ))
            break
            
    if saveAllPointsInFile!=None:
        numpy.savetxt(saveAllPointsInFile,points)

    # Finer validation on the first global pass

    doValidateOnD = False
  #  doValidateOnD = True

    p1_ref,p2_ref = patches[0]

    patches_ = []   # The new validated patches
    points_ = []   # The new validated points
    zz = []

    for index,(p1,p2) in enumerate(patches):
        if doValidateOnD:
            p1_ = util.flattenImage( p1.copy(), order=2 )
            corr = abs(numpy.sum(p1_ref*p1_)/numpy.linalg.norm(p1_ref)/numpy.linalg.norm(p1_))
            corr_th = 0.65
            corr_th = 0.6
#            corr_th = 0.58
#            corr_th = 0.55
            corr_th = 0.5
 #           corr_th = 0.0
        else:
            p2_ = util.flattenImage( p2.copy(), order=2 )
            corr = abs(numpy.sum(p2_ref*p2_)/numpy.linalg.norm(p2_ref)/numpy.linalg.norm(p2_))
            corr_th = 0.4
#            corr_th = 0.35
#            corr_th = 0.2
        zz.append(corr)
        if corr>corr_th or len(points_)<nmin:
            patches_.append((p1,p2))
            points_.append(points[index])


    patches = patches_
    if len(points_)!=0:
        points = numpy.row_stack(points_)
    else:
        points = numpy.zeros((0,2))

    log.info("Keep {0:} detected markers".format(len(points)))




#
#    hist = []
#    hist2 = []
#    hist3 = []
#    x = []
#    x2 = []
#
#    v0 = None
#    w0 = None
#    for i in range(1,len(vv),10):
#        pca = util.PCA(vv[:i],0.65)
#      #  pca = util.PCA(vv[:i],0.5)
#        print "n:%i -> %i,%.3f" %(i,pca.npc,pca.sumvariance[-1])
#        print pca.U.shape
#        v1 = pca.Vt[0,:]
#        w1 = pca.U[0,:]
#
#        if v0==None:
#            v0 = pca.Vt[0,:]
#            w0 = pca.U[0,:]
#        else:
#            a = numpy.sum(v0*v1)/numpy.linalg.norm(v0)/numpy.linalg.norm(v1)
#            b = numpy.sum(w0*w1)/numpy.linalg.norm(w0)/numpy.linalg.norm(w1)
#           # print "%.1f  %.1f" %(a,b)
#            print a
#
#        hist.append(pca.sumvariance[0])
#        x.append(i)
#        if len(pca.sumvariance)>1:
#            hist2.append((pca.sumvariance[1]-pca.sumvariance[0]))
#            hist3.append(pca.sumvariance[1])
#            x2.append(i)
#        if pca.npc>=2:
#            patches = patches[:i]
#            points = points[:i]
#            #break
#
#    der2 = numpy.diff(hist,n=1)


#    import pylab
#    pylab.figure()
#   # pylab.semilogy()
#    pylab.plot(zz)
# #   pylab.plot(zz_smooth)
##    pylab.plot(x,hist)
##    pylab.plot(x2,hist2)
##    pylab.plot(x2,hist3)
##    pylab.plot(x[1:],der2)
#    pylab.show()

    return points, D, patches


def rescaleModel( src, dest, stack, scale, size=20 ):
    '''
    :param scale: The scaling parameter. Equals to bin_dest/bin_src.
    '''

    print "Rescale Model {0:s} to {1:s}".format(src,dest)

    model_src = modl.Model(src)
    model_src.loadFromFile()

    model_dest = modl.Model(dest)

    sx,sy = stack.getImageScale()
    nx,ny = stack.getImageSize()
    nz = stack.getNumberOfImages()

    model_dest.nx = nx
    model_dest.ny = ny
    model_dest.nz = nz

    model_dest.imgref.xscale_new = sx
    model_dest.imgref.yscale_new = sy
    model_dest.imgref.zscale_new = sx/stack.bin # to fix...

    # First create a numpy array of the points

    XYZ = numpy.empty((0,5))

    for iobj,object_src in enumerate(model_src.objects):
        obj_dest = model_dest.addNewObject()
        for icont,contour_src in enumerate(object_src.contours):
            obj_dest.addNewContour()
            points = contour_src.points.copy()
            points[:,:2] = (points[:,:2]+(scale-1.0)/2.0)/scale
            ones = numpy.ones((len(points),1))
            XYZ = numpy.row_stack((XYZ,numpy.column_stack((ones*iobj,ones*icont,points))))
    
    # Second refine the points according to the slice

    n = stack.getNumberOfImages()

    npts = 0

    for z in range(n):
        index = numpy.where(XYZ[:,4]==z)
        npts += len(index[0])
        print "z={0:}  {1:} points".format(z,len(index[0]))
        xyz = XYZ[index]
        XYZ[index][:,2:] = refineMarkerPositions( xyz[:,2:], stack.getImageAt(z), width=size, height=size, clear=True )

    for iobj,icont,x,y,z in XYZ:
        iobj = int(iobj)
        icont = int(icont)
        if x!=numpy.nan and y!=numpy.nan and z!=numpy.nan:
            model_dest.objects[iobj].contours[icont].addPoint(x,y,z)

    print "Number of points: {0:} <-> {1:}".format(len(XYZ),npts)

    model_dest.save()


def refineMarkerPositions( pts, img, width=20, height=20, tol=5, clear=False, verbose=False ):
    '''Center the markers on an image stack'''

    half_width = width/2.0
    half_height = height/2.0

    for index,pt in enumerate(pts):
        
        ix,iy = pt[:2].astype('int')
        
        x1 = ix-half_width
        x2 = ix+half_width
        y1 = iy-half_height
        y2 = iy+half_height

        patch = img[x1:x2,y1:y2]
        patch = scipy.ndimage.filters.gaussian_filter( patch, 1.5 )
        patch = - otfutil.removeBackground(patch)

        tx,ty = scipy.ndimage.measurements.center_of_mass(patch)

        problem = tx+1-half_width>5 or ty+1-half_height>5
        if problem:
            problem = "Problem at index {0:}: pt=[{1:.1f},{2:.1f}]   t=[{3:.1f},{4:.1f}]".format(index,pt[0],pt[1],tx+1-half_width,ty+1-half_height)
        else:
            pt[0] = x1 + tx + 1
            pt[1] = y1 + ty + 1

        if problem and clear:
            pt[0] = numpy.nan
            pt[1] = numpy.nan

        if verbose or problem:
            print "tx={0:+6.2f}  ty={1:+6.2f}  {2:s}".format(tx+1-half_width,ty+1-half_height,problem)
        
    return pts


if __name__ == '__main__':

    option1 = { "directory":"/ncmirdata3/sphan/FHV-18/bin4", "fin":"FHV-18a.st", "size":6, "nmax":1500, "slices":None  }
    option2 = { "directory":"/ncmirdata3/sphan/FHV-18/bin4", "fin":"FHV-18a.st_low", "size":16, "nmax":500, "slices":[0] }
    option3 = { "directory":"/ncmirdata3/sphan/FHV-18/bin4", "fin":"FHV-18a.st", "size":6, "nmax":1500, "slices":None  }
    option4 = { "directory":"/ncmirdata3/sphan/FHV-18/bin4", "fin":"FHV-18a.st", "size":6, "nmax":1500, "slices":[60]  }
    option5 = { "directory":"/ncmirdata3/sphan/otf/mon", "fin":"SCN_night_OTO_cellchain001a.st", "size":4, "nmax":1500, "slices":[30,31]  }
    option6 = { "directory":"/ncmirdata3/sphan/otf/mon", "fin":"SCN_night_OTO_cellchain001a.st", "size":4, "nmax":1500, "slices":None  }
    option7 = { "directory":"/ncmirdata3/sphan/otf/mon2", "fin":"SCN_night_OTO_cellchain001a.st", "size":4, "nmax":1500, "slices":None  }
    option8 = { "directory":"/ncmirdata3/sphan/otf/25", "fin":"FHV-25-2a.st", "size":5, "nmax":1500, "slices":None  }
    option9 = { "directory":"/ncmirdata3/sphan/yongning", "fin":"test2a.st", "size":7, "nmax":1500, "slices":[60]  }
    option10 = { "directory":"/ncmirdata3/sphan/daniela/81780/data/bin4", "fin":"grid2_sect1_G4K_dish7area4_TS2-V2-mSOG_exp7_7_11_1a.st", "size":6, "nmax":1500, "slices":[60]  }
    option10 = { "directory":"/ncmirdata4/sphan/christine/Tumor_", "fin":"test-2.mrc", "size":6, "nmax":1500, "slices":[0]  }

    options = option10

#    fin = os.path.join( options["directory"], options["fin"])
#
#    findMarkersInStack( fin, nmax=options["nmax"], slices=options["slices"], size=options["size"] )


    option11 = { "directory":"/ccdbprod/ccdbprod2/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_20094/microscopy_83492/processed_data",
                 "basename":"grid2_sect3_G4K11K_2a",
                 "index":40 }
    options = option11

    import txbr.utilities

    model_path = os.path.join(options["directory"],"txbr-align","bin4","%s.mrk" %(options["basename"]))
    stack = txbr.utilities.loadStack( os.path.join(options["directory"],"CCDBID_83492"), "%s.st" %options["basename"])

    rescaleModel( model_path, "testa.modl", stack, 0.25 )


    options["basename"] = "grid2_sect3_G4K11K_2b"
    model_path = os.path.join(options["directory"],"txbr-align","bin4","%s.mrk" %(options["basename"]))
    stack = txbr.utilities.loadStack( os.path.join(options["directory"],"CCDBID_83492"), "%s.st" %options["basename"])

    rescaleModel( model_path, "testb.modl", stack, 0.25 )

    options["basename"] = "grid2_sect3_G4K11K_2c"

    model_path = os.path.join(options["directory"],"txbr-align","bin4","%s.mrk" %(options["basename"]))
    stack = txbr.utilities.loadStack( os.path.join(options["directory"],"CCDBID_83492"), "%s.st" %options["basename"])

    rescaleModel( model_path, "testc.modl", stack, 0.25 )







    import sys
    sys.exit(0)



    model_path = os.path.join(options["directory"],"txbr-detect","bin4","%s.bead.%03i.mod" %(options["basename"],options["index"]))
    model = modl.Model(model_path)
    model.loadFromFile()
    XY = model.points()



    #XY *= 4.0
    bin = 1.0/4.0
    XY[:,:2] = (XY[:,:2]+(bin-1.0)/2.0)/bin

    mrcfile = "/ccdbprod/ccdbprod2/home/CCDB_DATA_USER.portal/CCDB_DATA_USER/acquisition/project_20094/microscopy_83531/rawdata/CCDBID_83531/grid2_sect2_G4K11K_4a.st"
    m = mrc.MRCFile(mrcfile)
    img = m.getZSliceAt(options["index"])

    print model_path
    print XY

    refineMarkerPositions( XY, img )

    mnew = modl.Model("test.mod",model)
    obj = mnew.addNewObject()
    cont = obj.addNewContour()
    cont.points = XY
    mnew.save()

    print XY