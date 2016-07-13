import os.path
import sys
import numpy
import pylab
import mrc

from average import radialAverage


def getPowerSpectrumAt( filename, index, nx, ny, nbins = 300, filter=None ):
    '''Calculate the radially averaged power spectrum of an image. The image corresponds
    to the slice at a given 'index' from a MRC stack specified by its 'filenames'.'''

    file = mrc.MRCFile(filename)
    
    u = file.getZSliceAt(index)
    u = u[:nx,:ny]

    data = getPowerSpectrum(  u, nbins=nbins )

    numpy.savetxt( "%s.txt" %filename, data )

    return data



def getPowerSpectrum( u, filter=None, nbins = 300 ):
    '''Calculate the radially averaged power spectrum of an image.'''

    nx,ny = u.shape

    if filter=='wiener': u = scipy.signal.wiener(u)

    u = numpy.fft.fft2(u)

    if filter=="custom":
        x,y = numpy.meshgrid(range(nx),range(ny))
        kx = x - nx/2
        ky = y - ny/2
        k2 = kx**2 + ky**2
        f = 0.1*numpy.exp(-k2/500.0**2)
        f = numpy.ones(u.shape)*0.1
        f[0,0] = 1.0

        u = u*f

    pwspct = numpy.abs(u)**2
    pwspct = numpy.fft.fftshift(pwspct)

    k = numpy.linspace( 0, min(nx-1,ny-1)/2, nbins )
    f = radialAverage( pwspct, center=[nx/2,ny/2], nbins=nbins )

    data = numpy.column_stack((k,f))

    return data


def getScalingPowerLaw( scaling, valueAtOne=1.0 ):

    x = numpy.array([0.2,0.7])
    
    k = valueAtOne*x**scaling

    return k


def plotSpectrum( f, ref=None, title=None, legends=None, filename=None, loglog=True ):

	styles = [ '-v', ':.', '--o', '-*', '-+' ]

	pylab.figure()
	if loglog: pylab.loglog()
        else: pylab.semilogy()

        K = numpy.empty((0))

	p = []
	for index,t in enumerate(f):
		k = t[:,0]
		k = k/k[-1]
                K = numpy.concatenate((K,k))
		p += pylab.plot(k,t[:,1],styles[index%len(styles)],markersize=1.0)

	if ref!=None: p += pylab.plot( ref[0], ref[1] )
	if legends!=None: pylab.legend( p, legends )
        if loglog:
            K = K[numpy.argwhere(K>0)]
            xmin = numpy.min(K)
            xmax = numpy.max(K)
            pylab.xlim(xmin,xmax)

	pylab.xlabel(r"k/k$_{\rm Nyquist}$")
	pylab.ylabel(r"P(k)/P(0)")
	if title!=None: pylab.title(title)
	if filename!=None: pylab.savefig(filename)
	pylab.show()


if __name__=='__main__':
    
    nx = 4096
    ny = nx

    nx_de12 = 3000
    ny_de12 = nx_de12

    directory = "/ncmirdata3/sphan/noise"

    f_de12 = "de12_bright-2.mrc"
    f_tietz = "tietz_bright-2.mrc"
    f_gatan = "gatan_bright-2.mrc"

    f_de12 = os.path.join(directory,f_de12)
    f_tietz = os.path.join(directory,f_tietz)
    f_gatan = os.path.join(directory,f_gatan)

    calculate_data = 'calculate' in sys.argv

    nbins = 300

    if calculate_data:

            f1 = getPowerSpectrumAt( f_de12, 0, nx_de12, ny_de12, nbins=nbins )
            f2 = getPowerSpectrumAt( f_tietz, 0, nx, ny, nbins=nbins )
            f3 = getPowerSpectrumAt( f_gatan, 0, nx, ny, nbins=nbins )

    else:

            f1 = numpy.loadtxt( f_de12 + ".txt" )
            f2 = numpy.loadtxt( f_tietz + ".txt" )
            f3 = numpy.loadtxt( f_gatan + ".txt" )
            
    f1[:,1] = f1[:,1]/f1[0,1]
    f2[:,1] = f2[:,1]/f2[0,1]
    f3[:,1] = f3[:,1]/f3[0,1]

    f = (f1,f2,f3)
    ref = getScalingPowerLaw( -2, valueAtOne=numpy.nanmin(f[0][1]) )

    title = "Noise Power Spectrum"
    legends = ("DE12", "Tietz", "Gatan", r"1/k^2")
    filename="noise_spectrum.png"

    plotSpectrum( f , ref=ref, title=title, legends=legends, filename=filename )
