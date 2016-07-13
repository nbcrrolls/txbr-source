import numpy
import scipy.ndimage.filters

def removeBackground( u ):

    nx,ny = u.shape

    u = scipy.ndimage.filters.gaussian_filter(u,1.5)
  #  u = scipy.ndimage.filters.maximum_filter(u,3)

    ratio = 0.20
    
    lx = max(2,int(ratio*nx/2.0))
    ly = max(2,int(ratio*ny/2.0))

    if 2*lx>=nx or 2*ly>=ny: return u

    background = (numpy.sum(u)-numpy.sum(u[lx:-lx,ly:-ly]))/(nx*ny-(nx-2.0*lx)*(ny-2.0*ly))

    #print "(lx,ly)=(%i,%i) min,max=%.2f,%.2f  background=%.2f" %(lx,ly,numpy.min(u),numpy.max(u),background)

    th = background + 0.7*(numpy.min(u)-background)

    v = numpy.where(u<=th,u-th,0.0)

    v = scipy.ndimage.filters.gaussian_filter(v,1.2)

    return v


def changeToPolarCoordinate( u, origin=None ):

    nx,ny = u.shape
    nr = min(nx-10,ny-10)/2.0
    nr = 250

    if origin==None:
        origin = numpy.array([nx/2,ny/2],dtype='int')

    nr_x = min(nx-10-origin[0],origin[0]-10)
    nr_y = min(ny-10-origin[1],origin[1]-10)
    nr = min(nr_x,nr_y)
    nr = min(250,nr)

    theta = numpy.linspace(0,2.0*numpy.pi,1000)
    r = numpy.linspace(nr/4,nr,1000)
    r = numpy.linspace(0,nr,500)
   # r = numpy.logspace(1,numpy.log10(nr),1000) - 1.0
##    r = numpy.logspace(1,numpy.log10(nr),1000) - 1.0
##    r = nr-r


    x = numpy.outer(numpy.cos(theta),r) + origin[0]
    y = numpy.outer(numpy.sin(theta),r) + origin[1]

    x0 = x.astype('int')
    y0 = y.astype('int')

#    print x0
#    print y0

    px = x - x0
    py = y - y0

#    v = numpy.zeros((len(theta),len(r)))
#
#    for i in range(len(theta)):
#        for j in range(len(r)):
#            v[i,j] = u[x0[i,j],y0[i,j]]

    v = (1.0-px)*(1.0-py)*u[x0,y0] + px*(1.0-py)*u[x0+1,y0] + (1.0-px)*py*u[x0,y0+1] + px*py*u[x0+1,y0+1]

#    v = numpy.empty_like(u)
#
#    for index,x_ in enumerate(x0[:-1]):
#        v[index,:] = px[index]*u[x_+1,:] + (1.0-px[index])*u[x_,:]
#    v[-1,:] = u[-1,:]
#
#    u = v.copy()
#
#    for index,y_ in enumerate(y0[:-1]):
#        v[:,index] = py[index]*u[:,y_+1] + (1.0-py[index])*u[:,y_]
#    v[:,-1] = u[:,-1]

    return v
