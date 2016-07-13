import numpy, pylab
import statsplot

from numpy.fft import fft2
from paircorrelation import PairCorrelationFunction_2D, \
                            PairCorrelationFunction_3D

def shiftArrayToCenter(a):
    
    u = numpy.empty_like(a)

    n1,n2 = a.shape
        
    n1 = n1/2
    n2 = n2/2
    
    u[:n1,:n2] = a[-n1:,-n2:]
    u[-n1:,:n2] = a[:n1,-n2:]
    u[:n1,-n2:] = a[-n1:,:n2]
    u[-n1:,-n2:] = a[:n1,:n2]
    
    return u


def getPairCorrelation( XYZ ):
    '''Return the pair correlation function for XYZ.
    XYZ is a numpy array of shape (n,dim); dim can be 2 or 3.
    '''
    
    n,dim = XYZ.shape

    rg = numpy.max(XYZ,axis=0)-numpy.min(XYZ,axis=0)
    
    S = numpy.max(rg)
    rMax = numpy.max(rg)
    dr = 0.5
    
    print 'S=%f    r=%f    dr=%f' %(S,rMax,dr)
    
    if dim==2:
        (gr, r, x, y) =  PairCorrelationFunction_2D(XYZ[:,0],XYZ[:,1],S,rMax,dr)
        return (r,gr)
    elif dim==3:
        (gr, r, x, y, z) =  PairCorrelationFunction_3D(XYZ[:,0],XYZ[:,1],XYZ[:,2],S,rMax,dr)
        return (r,gr)
    
    return None
    
    
#def subVol( XYZ, xlim=None, ylim=None, zlim=None):
#    '''Select points from XYZ that belongs to a box specified by x, y and z.
#    xlim = (xmin,xmax), ylim = (ymin,ymax) and zlim = (zmin,zmax)'''
#
#    points = XYZ
#
#    if xlim!=None:
#        indices, = numpy.where((points[:,0]>=xlim[0])*(points[:,0]<xlim[1]))
#        points = points[indices]
#
#    if ylim!=None:
#        indices, = numpy.where((points[:,1]>=ylim[0])*(points[:,1]<ylim[1]))
#        points = points[indices]
#
#    print points.shape
#
#    if zlim!=None:
#        indices, = numpy.where((points[:,2]>=zlim[0])*(points[:,2]<zlim[1]))
#        points = points[indices]
#
#    print points.shape
#
#    return points
#
#
#def subVol2( XYZ, u=[0.0,0.0,1.0], l=[100,200]):
#
#    points = XYZ
#
#    print "Min/Max points: %s/%s" %(numpy.min(XYZ,axis=0),numpy.max(XYZ,axis=0))
#
#    u = numpy.asarray(u)
#    l = numpy.asarray(l)
#
#    cote = numpy.dot(points,u)
#
#    print "Min/Max cote: %s/%s" %(numpy.min(cote,axis=0),numpy.max(cote,axis=0))
#
#    indices, = numpy.where((cote[:]>=l[0])*(cote[:]<l[1]))
#    points = points[indices]
#
#    return points
#
#
#def sliceVol( XYZ, u=[0.0,0.0,1.0], X0=[0,0,0], l=[100,100], xlim=None, ylim=None, zlim=None):
#    '''This routine returns a set of 3D points from XYZ. The set will
#    belong to a slice contained into two parallel planes whose normal
#    is given by u separated by a distance l=l[1]+l[0]. Points X0 will
#    be between the planes at a distance l[0] of one of them and l[1]
#    of the other.
#    '''
#
#    if xlim!=None or ylim!=None or zlim!=None:
#        XYZ = subVol(XYZ, xlim=xlim, ylim=ylim, zlim=zlim)
#
#    u = numpy.asarray(u)
#    X0 = numpy.asarray(X0)
#    l = numpy.asarray(l)
#
#    cote = numpy.dot(XYZ-X0,u)
#
#    indices, = numpy.where((cote[:]>=-l[0])*(cote[:]<l[1]))
#
#    return XYZ[indices]


def fourier( XYZ, nx, ny, nz ):
    
    print "nx=%i  nx=%i  nx=%i" %(nx, ny, nz)
    
    f_xy = numpy.zeros((nx,ny))
    f_xz = numpy.zeros((nx,nz))
    f_yz = numpy.zeros((ny,nz))
    
    for pt in XYZ:
        
        ix,iy,iz = pt.astype('int')
        
        f_xy[ix,iy] += 1
        f_xz[ix,iz] += 1
        f_yz[iy,iz] += 1
        
    f_xy = fft2(f_xy).imag**2
    f_xz = fft2(f_xz).imag**2
    f_yz = fft2(f_yz).imag**2
    
    f_xy = shiftArrayToCenter(f_xy)
    f_xz = shiftArrayToCenter(f_xz)
    f_yz = shiftArrayToCenter(f_yz)
    
    pylab.figure()
    pylab.imshow(f_xy)
    pylab.show()
    
    pylab.figure()
    pylab.imshow(f_xz)
    pylab.show()
    
    pylab.figure()
    pylab.imshow(f_yz)
    pylab.show()
        

def neighbours( XYZ, index, n=20 ):
    '''Get the n closest neighbours in XYZ from XYZ[index].'''
    
    point0 = XYZ[index]
    
    r = XYZ - point0
    distance = numpy.sqrt(r[:,0]**2 + r[:,1]**2 + r[:,2]**2)
    indices = numpy.argsort(distance)
    
    XYZ = XYZ[indices]
    
    return XYZ[:n]



        
    