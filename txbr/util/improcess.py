import numpy
import scipy.ndimage.measurements
import pylab

from PIL import Image

def __gaussianKernel__( nx, ny, sigma=0.25 ):
    '''Creates a gaussian kernel image.

    :param nx,ny: The width and heigh of the kernel image. Should be odd numbers;
        if not, the size is increased by one unit.
    :param sigma: Characteristic length of the kernel.
    :returns: The gaussian kernel :math:`G_sigma`
    '''

    print "Create Gaussian Kernel for the detector image"

    # Make the box odd, so center corresponds to a discret point.

    if nx%2==0: nx+=1
    if ny%2==0: ny+=1

    origin = ( (nx-1)/2, (ny-1)/2 )
    y,x = numpy.meshgrid(range(ny),range(nx))

    x = x - origin[0]
    y = y - origin[1]

    r2 = x**2 + y**2

    G_sigma = 1/2.0/numpy.pi/sigma*numpy.exp(-r2/2.0/sigma**2)

    return G_sigma



def loadImage( nameOfImage ):

    img = Image.open(nameOfImage).convert("I")    # Directly convert in the 32bit integer mode

    u = numpy.asarray(img)

    return u[::-1,:].T.astype('float')


def saveImage( nameOfImage, u, normalize=False, type=None):
    '''mode should be either int32 and float32'''

    print "Save image %s" %nameOfImage

    if type==None: type = u.dtype
    
    int_types = [numpy.int,numpy.int8,numpy.int16,numpy.int32,numpy.int64,numpy.uint8,numpy.uint16,numpy.uint32,numpy.uint64]
    float_types = [numpy.float,numpy.float16,numpy.float32,numpy.float64]

    if not type in int_types or not type in float_types:
        type = numpy.int32

    if type in int_types:
        if normalize:
            scale = float(numpy.iinfo('int16').max)/(numpy.max(u)-numpy.min(u))
            array = (u[:,::-1].T-numpy.min(u))*scale
            array = array.astype('int16')
        else:
            array = u[:,::-1].T.astype('int16')
    if type in float_types:
        if normalize:
            scale = 1.0/numpy.max(numpy.absolute(u))
            array = u[:,::-1].T*scale
            array = array.astype('float16')
        else:
            array = u[:,::-1].T.astype('float16')
    
    img = Image.fromarray(array)
    img.save(nameOfImage)
    

def buildHistogram( u, nbins=1000, doPlot=False ):
    '''Check the histogram of an image u'''

    umin = numpy.min(u)
    umax = numpy.max(u)

    x = numpy.arange(umin,umax,(umax-umin)/nbins)

    hist = scipy.ndimage.measurements.histogram( u, umin, umax, nbins )

    hist_der = numpy.diff(hist)
    hist_hess = numpy.diff(hist,n=2)

    print x[2+numpy.argmin(hist_hess)]

    if doPlot:
        pylab.figure()
        pylab.plot(x,hist)
        pylab.plot(x[1:],hist_der)
        pylab.plot(x[2:],hist_hess)
        pylab.show()


def flattenImage( img, order=4 ):
    '''Flatten image with a polynomial fit of order order'''

    nx,ny = img.shape

    x,y = numpy.meshgrid(range(ny),range(nx))

    v = [ ]
    for n in range(order+1):
        for i in range(n+1):
            v.append(x**i*y**(n-i))

    v = numpy.array(v)

   # v = numpy.array([numpy.ones((nx,ny)), x, y, x**2, x * y, y**2])

    img.resize((nx*ny))
    v.resize((v.shape[0],nx*ny))

    coefficients, residues, rank, singval = numpy.linalg.lstsq(v.T, img)

    flat_img = numpy.resize(img - numpy.dot(coefficients,v),(nx,ny))

    img.resize((nx,ny))

    return flat_img


def highpass(image, sigma=5):

    f = numpy.fft.rfft2(image)
    g = scipy.ndimage.fourier_gaussian(f, sigma, n=image.shape[1])
    h = numpy.fft.irfft2(g)
    return h.real


def histeq(im,nbr_bins=256):

   #get image histogram
   imhist,bins = numpy.histogram(im.flatten(),nbr_bins,normed=True)
   cdf = imhist.cumsum() #cumulative distribution function
   cdf = 255 * cdf / cdf[-1] #normalize

   #use linear interpolation of cdf to find new pixel values
   im2 = numpy.interp(im.flatten(),bins[:-1],cdf)

   return im2.reshape(im.shape), cdf
