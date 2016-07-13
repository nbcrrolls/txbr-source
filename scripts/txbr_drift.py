#!/usr/bin/python

import os
import sys
import getopt
import multiprocessing
import numpy
import numpy.fft
import pylab
import scipy.ndimage.filters
import scipy.cluster.hierarchy
import mrc
import modl
import util

from txbr.utilities import loadStack
from txbr.onthefly import detect_task_star

from txbr.txbrconf import DETECT_DIR


def usage():
    """Usage for the command *txbr_sirt.py*."""

    print 'Usage: %s.py -b basename[,...]  [Options]' % os.path.basename(sys.argv[0])
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files (.fid and .txbr)"
    print "    --wd work_directory (default to directory)"
    print "        Name of the directory where the output files are stored"
    print "    --bin n"
    print "        Bin the data in the x and y directions"
    print "    --nproc value [default: 20 ]"
    print "        The number of CPUs to use."
    print "    -h or --help"
    print "        Help Information"


def hamming( u ):

    return u

    nx,ny = u.shape

    d = numpy.outer(numpy.hamming(nx),numpy.hamming(ny))
    u = u*d

    return u


def createDerKernel( nx, ny, sigma=0.25 ):
    '''Creates the kernel images (to implement gaussian derivatives) to detect beads
    of a certain size.

    :param nx,ny: The width and heigh of the kernel image. Should be odd numbers;
        if not, the size is increased by one unit.
    :param sigma: Characteristic length of the kernel. Should be 1/5 of the bead
        diameter in practice.
    :returns: A tuple containing the 3 image kernel :math:`G_xx`, :math:`G_yy`
        and :math:`G_xy`
    '''

    # Make the box odd, so center corresponds to a discret point.

    if nx%2==0: nx+=1
    if ny%2==0: ny+=1

    origin = ( (nx-1)/2, (ny-1)/2 )
    y,x = numpy.meshgrid(range(ny),range(nx))

    x = x - origin[0]
    y = y - origin[1]

    r2 = x**2 + y**2

    Gx = x/2.0/numpy.pi/sigma**2*numpy.exp(-r2/2.0/sigma**2)
    Gy = y/2.0/numpy.pi/sigma**2*numpy.exp(-r2/2.0/sigma**2)

    return Gx,Gy


def createDerivativeStacks( stack ):
    """Create the derivative stacks."""

    output_x = "%s.derx.mrc" %(stack.basename)
    output_y = "%s.dery.mrc" %(stack.basename)

    if os.path.lexists(output_x) and os.path.lexists(output_y): return

    nx,ny = stack.getImageSize()
    nz = stack.getNumberOfImages()

    size = 25
    lx,ly = int(size),int(size)
    Gx,Gy = createDerKernel( lx, ly, sigma=1.5 )

    fx = mrc.MRCFile(output_x)
    fx.setHeader( nx, ny, nz )

    fy = mrc.MRCFile(output_y)
    fy.setHeader( nx, ny, nz )

    for iz in range(nz):
        print "Derivative slice for index #%i" %iz
        u = stack.getImageAt(iz)
        Ix = scipy.ndimage.filters.convolve( u, Gx )
        Iy = scipy.ndimage.filters.convolve( u, Gy )
        fx.setZSliceAt( iz, Ix )
        fy.setZSliceAt( iz, Iy )

    fx.updateHeader()
    fy.updateHeader()



def __plot_transformations__( transformations ):

    tx = [ tf.T[0] for tf in transformations ]
    ty = [ tf.T[1] for tf in transformations ]
    axx = [ tf.M[0,0]-1 for tf in transformations ]
    axy = [ tf.M[0,1] for tf in transformations ]
    ayx = [ tf.M[1,0] for tf in transformations ]
    ayy = [ tf.M[1,1]-1 for tf in transformations ]

    pylab.figure()
    pylab.plot(tx,label="$t_x$")
    pylab.plot(ty,label="$t_y$")
    pylab.xlabel("N\%50")
    pylab.ylabel("Translation Coefficients")
    pylab.legend()
    pylab.title("Beam effect on an insect flight muscle (29k)")
    pylab.savefig("corr_tranlations.png")

    pylab.figure()
    pylab.plot(axx,label="$a_{xx}-1$")
    pylab.plot(axy,label="$a_{xy}$")
    pylab.plot(ayx,label="$a_{yx}$")
    pylab.plot(ayy,label="$a_{yy}-1$")
    pylab.xlabel("N\%50")
    pylab.ylabel("Linear Coefficients")
    pylab.legend()
    pylab.title("Beam effect on an insect flight muscle (29k)")
    pylab.savefig("corr_linear.png")


def __plot_corr_peaks__( corr_peaks ):

    pylab.figure()
    pylab.plot(corr_peaks)
    pylab.ylim(0.0,0.25)
    pylab.xlabel("N\%50")
    pylab.ylabel("Correlation Peaks")
    pylab.title("Correlation relative to First micrograph")
    pylab.savefig("corr_peaks.png")


def detect( stack, **keywords ):
    """Bead detection routine."""

    detect_dir = os.path.join( keywords['work_directory'], DETECT_DIR, "bin%i" %stack.bin )
    if not os.path.exists(detect_dir): os.makedirs(detect_dir)

    size = keywords['size']
    nmax = keywords['nmax']
    overwrite = keywords['overwrite']

    print "Bead detection (size=%i)" %size

    if keywords['test-detection']:
        overwrite = True
        views = [stack.getNumberOfImages()/2]
    else:
        views = range(stack.getNumberOfImages())

    pool = multiprocessing.Pool( processes=20 )

    tasks = [ [ stack, index, detect_dir, overwrite, nmax, size ] for index in views ]

    markers = pool.map( detect_task_star, tasks)

    peaks = numpy.array([ len(pts) for pts in markers ])
    npeak =  int(scipy.stats.cmedian(scipy.stats.trimboth(peaks,0.15)))

    print "npeak=%i [Median value of the markers versus tilt]" %(npeak)

    if keywords['test-detection']:
        sys.exit(0)

    final_detect_mod = os.path.join(detect_dir,"%s.bead.final.mod" %(stack.basename))
    nx,ny = stack.getImageSize()
    modl.saveAs3DModel( final_detect_mod , markers, nx, ny, len(markers), doShow=keywords["show-beads"], mrcFileName=stack.mrc_file.filename )

    return markers


def align( stack, markers, **keywords ):
    """Alignment routine."""

    nx,ny = stack.getImageSize()
    nz = stack.getNumberOfImages()

    views = range(nz)
    niter = keywords['niter']

    transformations = [ util.I2D for view in views ]
    
    for iter in range(niter):

        if iter!=0:
            print "Clusterize at iteration #%i" %iter
            XY = util.clusterize(numpy.row_stack(XY_[::20])[:2000],min_cluster_size=5)
            print "\niter %i  XY: %s" %(iter,XY.shape)
        else:
            print "Iteration #%i" %iter


        XY_ = []

        for index in range(nz):

            if iter==0:
                XY = markers[index>0 and index-1 or 0][:,:2]

            pts = markers[index]
            tf = transformations[index]

            index1,index2 = util.connectPairWise( XY, pts[:,:2], tf, dm = 2.5 )

            pts1_ = XY[index1]
            pts2_ = pts[index2][:,:2]
            
            if pts1_.size==0: continue

            pts1_ = pts1_[:,:2].copy()
            pts2_ = pts2_[:,:2].copy()

            if iter==0:
                p = util.polynomial_regression( pts2_, pts1_, [1,1] )
                warp2D = util.AffineTransformation2D( numpy.array([p[0].coeff,p[1].coeff]) )
            else:
                p = util.polynomial_regression( pts2_, pts1_, [3,3] )
                warp2D = util.PolynomialTransformation2D( 3, numpy.array([p[0].coeff,p[1].coeff]) )

            if iter==0 and index!=0:
                transformations[index] = transformations[index-1].compose(warp2D) # Only compose affine maps!!!
            else:
                transformations[index] = warp2D

            pts2_img = transformations[index].forward(pts2_)

            if iter!=0:
                diff = pts2_img-pts1_
                print "index:%i: n=%i diff=%s" %(index,len(pts1_),numpy.max(numpy.sqrt(diff[:,0]**2+diff[:,1]**2)))

            XY_.append(pts2_img)


    try:
        __plot_transformations__( [ tf.inv().compose(transformations[0]) for tf in transformations ] )
    except:
        pass

    xy_mod = os.path.join(keywords['work_directory'],"%s.xy.mod" %(stack.basename))
    XYZ = [numpy.column_stack((XY,[iz]*len(XY))) for iz in range(nz)]

    modl.saveAs3DModel( xy_mod , XYZ, nx, ny, nz, doShow=False, mrcFileName=stack.mrc_file.filename )

    return transformations


def createAlignedStack( stack, transformations, **keywords ):
    """Create the aligned stack"""

    if not keywords['create-stack']: return

    padx = 200
    pady = 200

    (nx,ny) = stack.getImageSize()
    nz = stack.getNumberOfImages() + 1

    nx -= 2*padx
    ny -= 2*pady

    u0 = 0
    umean = 0.0



    output = "output.mrc"
    print "Create Stack %s" %output
    
    corr_peaks = []
    
    iz0 = 0
    iz0 = nz-2

    u0 = stack.getImageAt(iz0)
    warp2D0 = transformations[iz0]
    u0 = util.warp( u0, warp2D0 )
    u0 = u0[padx:nx+padx,pady:ny+pady]
    u0 = hamming( u0 )

    f = mrc.MRCFile(output)
    f.setHeader( nx, ny, nz )

    for iz in range(nz-1):
        u = stack.getImageAt(iz)
        warp2D = transformations[iz]
#        warp2D.T[0] = 0
      #  warp2D = util.AffineTransformation2D( numpy.array([[0.0,1.0,0.0],[0.0,0.0,1.0]]) )
        u = util.warp( u, warp2D )
        u = u[padx:nx+padx,pady:ny+pady]
        u = hamming( u )
        f.setZSliceAt( iz, u )
        umean = u/nz + umean
        u1 = numpy.fft.fft2(u0)
        u2 = numpy.fft.fft2(u)
        u =  u1*numpy.conj(u2)/numpy.abs( u1*u2 )
        u = numpy.fft.ifft2(u)

        print "micrograph #%i:  correlation: %.2f" %(iz,numpy.max(u.real))
        print numpy.argmax(u.real)
        corr_peaks.append(numpy.max(u.real))

       # f.setZSliceAt( iz, numpy.fft.fftshift(u.real))


    f.setZSliceAt( nz-1, umean )

    f.updateHeader()



    output = "output_.mrc"

    print "Create Stack %s" %output

    der_x = "%s.derx.mrc" %(stack.basename)
    der_y = "%s.dery.mrc" %(stack.basename)

    fder_x = mrc.MRCFile(der_x)
    fder_y = mrc.MRCFile(der_y)

    output_x = "output.res.x.mrc"
    output_y = "output.res.y.mrc"

    f = mrc.MRCFile(output)
    f.setHeader( nx, ny, nz-1 )

    for iz in range(nz-1):
        u = stack.getImageAt(iz)
        derx = fder_x.getZSliceAt(iz)
        dery = fder_y.getZSliceAt(iz)
        warp2D = transformations[iz]
        u = util.warp( u, warp2D )
        derx = util.warp( derx, warp2D )
        dery = util.warp( dery, warp2D )
        u = u[padx:nx+padx,pady:ny+pady]
        derx = derx[padx:nx+padx,pady:ny+pady]
        dery = dery[padx:nx+padx,pady:ny+pady]
        u = hamming( u )
      #  f.setZSliceAt( iz, (umean-u)*derx/numpy.sqrt(derx**2+dery**2) )
        f.setZSliceAt( iz, (umean-u) )

    f.updateHeader()



#    umean = hamming( umean )
#
#    for iz in range(nz-1):
#        u = stack.getImageAt(iz)
#        warp2D = transformations[iz]
#        u = util.warp( u, warp2D )
#        u = u[padx:nx+padx,pady:ny+pady]
#        u = hamming( u )
#        u1 = numpy.fft.fft2(umean)
#        u2 = numpy.fft.fft2(u)
#        u =  u1*numpy.conj(u2)/numpy.abs( u1*u2 )
#        u = numpy.fft.ifft2(u)
#        print "micrograph #%i:  correlation: %.2f" %(iz,numpy.max(u.real))
#        print numpy.argmax(u.real)
#        corr_peaks.append(numpy.max(u.real))
#
#



    numpy.savetxt( "corr_peaks.txt", corr_peaks)

    __plot_corr_peaks__( corr_peaks )

    sx, sy, sz = mrc.scale(stack.mrc_file.filename)
    mrc.update_scale( output, sx, sy, sz )
    mrc.update_origin( output, (1-padx)*sx, (1-pady)*sy, 0.0)

    model = os.path.join(keywords['work_directory'],"%s.xy.mod" %(stack.basename))
    os.system("imod %s %s" %(output,model==None and "" or model))


def main():
    """The main routine for this iterative process"""

    try:

        flags1 = "hf:b:"
        flags2 = [ "help", "directory=", "wd=", "bin=", "size=", "niter=", "nproc=", "test-detection", "stack" ]

        opts, args = getopt.getopt( sys.argv[1:], flags1, flags2 )

    except getopt.GetoptError, err:

        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    keywords = {}

    keywords['directory'] = "."
    keywords['work_directory'] = "."
    keywords['basename'] = None
    keywords['bin'] = None
    keywords['nproc'] = 20
    keywords['size'] = 5
    keywords['nmax'] = 2000
    keywords['niter'] = 4
    keywords['overwrite'] = False
    keywords["show-beads"] = False
    keywords['test-detection'] = False
    keywords['create-stack'] = False
    
    for option,value in opts:
        if option in ('-h','--help'):
            usage()
            sys.exit()
        elif option in ('-d','--directory'):
            keywords['directory'] = value
        elif option in ('--wd'):
            keywords['work_directory'] = value
        elif option in ('-b'):
            keywords['basename'] = value
        elif option in ('--bin'):
            keywords['bin'] = int(value)
        elif option in ('--size'):
            keywords['size'] = float(value)
        elif option in ('--niter'):
            keywords['niter'] = int(value)
        elif option in ('--nproc'):
            keywords['nproc'] = int(value)
        elif option in ('--test-detection'):
            keywords['test-detection'] = True
        elif option in ('--stack'):
            keywords['create-stack'] = True
        else:
            assert False, "unhandled option"

    if keywords['basename']==None:
        usage()
        sys.exit()

    stack = loadStack( keywords['directory'], keywords['basename'], workDirectory=keywords['work_directory'], bin=keywords['bin'] )

    createDerivativeStacks( stack )


    markers = detect( stack, **keywords )

    transformations = align( stack, markers, **keywords )

    createAlignedStack( stack, transformations, **keywords )


if __name__ == '__main__': main()