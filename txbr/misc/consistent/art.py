import sys
import numpy
import numpy.linalg
import scipy.misc
import numpy.random

numpy.set_printoptions(precision=1)


def get2DTestArray( file=None, nx=21, ny=21):
    '''Generate a small 2D test array, or read it from a file'''

    if file!=None and os.path.exists(file):
        return scipy.misc.imread(file)

    ox = int(nx/2) + 1
    oz = int(ny/2) + 1
    half_width = 1

    im2D = numpy.zeros((nx,ny))
    im2D[ox-half_width:ox+half_width,oz-half_width:oz+half_width] = 1.

    return im2D


def generateTrajectoryTraces( npixel, tilts, dest_shape ):
    '''Generate the trajectory traces as they match in the reconstruction array.
    npixel: Number of pixels in the projections
    tilts: list of tilt angles in radians
    des_shape: Shape of the reconstruction array'''

    ntilt = len(tilts)
    nx,nz = dest_shape

    r = numpy.zeros((npixel,ntilt,nx,nz))
    ox = int(nx/2) + 1
    oz = int(nz/2) + 1

    for index,t in enumerate(tilts):
        for ix_im in range(npixel):
            for ix_vol in range(nx):
                z = oz + numpy.sin(t)*(ix_vol-ox) + (oz-ix_im)/numpy.cos(t)
                iz = int(z)
                #print "ix_im,index,ix_vol,iz= %i,%i,%i,%i" %(ix_im,index,ix_vol,iz)
                if iz<0 or iz>=nz: continue
                else: r[ix_im,index,ix_vol,iz] = 1

    return r


def project( r, vol):
    '''Generate the projection of a volume'''

    # r has a shape (npixel,ntilt,nx,nz)
    # vol has a shape (nx,nz)
    # the result should be of size should be of shape (npixel,ntilt)

    return numpy.tensordot(r,vol,axes=([2,3],[0,1]))



def directReconstruct( r, im ):

    # r has a shape (npixel,ntilt,nx,nz)
    # im has a shape (npixel,ntilt)
    
    nx = r.shape[-2]
    nz = r.shape[-1]

    print r.shape
    print im.shape

    r = numpy.resize(r,(im.size,nx*nz))

    im = im.ravel()

    rinv = numpy.linalg.pinv(r)

    vol = numpy.dot(rinv,im)

    vol.resize((nx,nz))

    return vol


def consistentReconstruction(  r, im, init = None ):

    # r has a shape (npixel,ntilt,nx,nz)
    # im has a shape (npixel,ntilt)

    if init==None:
        init = numpy.zeros((nx,nz))
        
    vol = init

    mix = 0.01
    n = 500

    norm = numpy.sum(numpy.sum(r**2,axis=0),axis=0)
    print norm.shape

    for i in range(n):

        vol_2 = vol + mix*numpy.tensordot((im - numpy.tensordot(r,vol)),r,axes=((0,1),(0,1)))/norm
        print "iter #%i: %s     %s" %(i,numpy.max(vol_2-vol),numpy.max(im - numpy.tensordot(r,vol)))
        vol = vol_2

    return vol





vol = get2DTestArray(nx=100,ny=100)
nx,nz = vol.shape

theta = numpy.radians(numpy.arange(-60,62,1))
#theta = numpy.radians(numpy.arange(-90,91,1))
ntilt = len(theta)

npixel = nx
tilts = theta

r = generateTrajectoryTraces( npixel, tilts, vol.shape )
im = numpy.tensordot(r,vol)



#im_grad_m =  numpy.roll(im,0,0) - numpy.roll(im,-1,0)
#im_grad_p =  numpy.roll(im,1,0) - numpy.roll(im,0,0)
#
#
#print im[:,0]
#
#i1 = numpy.where(im[:]==0.0)
#i2 = numpy.where(im[:]>=1)
#i2 = numpy.where((im_grad_p[:]!=0) | (im_grad_m[:]!=0))
#
#print i2
#
#im = im[i2]
#r = r[i2]

#print im.shape
#print r.shape
#
#sys.exit()


#r1 = numpy.sum(r[i1],axis=0)
#r2 = numpy.sum(r[i2],axis=0)
#
#r = numpy.array((r1,r2))
#im = numpy.array((0,1))
#
#print r.shape

#r.resize((im.size,nx*nz))











#r.resize((nx*ntilt,nx*nz))

#print "trajectory done"
#
#rinv = numpy.linalg.pinv(r)
#
#print rinv.shape
#
#im = im.ravel()
#
#X = numpy.dot(rinv,im)
#
#X.resize((nx,nz))

X = directReconstruct( r, im )
init = vol + 0.1*numpy.random.random(vol.shape)
#init=vol
#init=None
#X = consistentReconstruction(  r, im, init=init )


res = X - vol
print "Residual: %f" %(numpy.max(numpy.abs(res)))

#print numpy.max(X)
#index =  numpy.argmax(X)

#(imax,jmax) = numpy.unravel_index(index,X.shape)
#print "(imax,jmax)=(%i,%i)" %(imax,jmax)

scipy.misc.imshow(X)
