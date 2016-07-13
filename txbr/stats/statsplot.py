import numpy, pylab
from txbr.utilities import plot3DMarkers 

def plotCorrelation(r,g):
    '''Plot the correlation function g(r)'''
    
    pylab.figure()
    pylab.plot(r,g)
    pylab.xlabel("r")
    pylab.ylabel("g(r)")
    pylab.show()
    
    
def plot3D(XYZ):
            
    try:
        from enthought.mayavi import mlab
    except ImportError:
        log.info('Mlab is not available!')
        return

    scale = [1.0,1.0,1.0]
    
    x = XYZ[:,0].astype(float).ravel()
    y = XYZ[:,1].astype(float).ravel()
    z = XYZ[:,2].astype(float).ravel()
    
    #sf = numpy.max(XYZ)/50.0
    sf = numpy.max(XYZ)/100.0

    n = XYZ.shape[0]
    cs = numpy.arange(4*n,5*n)/5.0/n
    
    fig = mlab.figure()
    
    #mlab.points3d(x,y,z,cs,colormap='copper',scale_factor=sf)
    mlab.points3d(x,y,z,color=(0.9,0.9,0.9),scale_factor=sf)
    
    
def plotZFrame(x=[0.0,1.0],y=[0.0,1.0],z=[0.0,1.0]):
    '''Plot two boundary planes at z=zmin and z=zmax'''
    
    try:
        from enthought.mayavi import mlab
    except ImportError:
        log.info('Mlab is not available!')
        return

    xmin,xmax = x
    ymin,ymax = y
    zmin,zmax = z
    
    c = (0.1,0.1,0.1)

    box_zmin = numpy.array([[xmin,ymin,zmin],
                       [xmax,ymin,zmin],
                       [xmax,ymax,zmin],
                       [xmin,ymax,zmin],
                       [xmin,ymin,zmin]])
    
    mlab.plot3d(box_zmin[:,0], box_zmin[:,1], box_zmin[:,2], tube_radius=1, color=c)

    box_zmax = numpy.array([[xmin,ymin,zmax],
                       [xmax,ymin,zmax],
                       [xmax,ymax,zmax],
                       [xmin,ymax,zmax],
                       [xmin,ymin,zmax]])
    
    mlab.plot3d(box_zmax[:,0], box_zmax[:,1], box_zmax[:,2], tube_radius=1, color=c)
    

def plotReferenceAxes(X1,X2,X3,X4,tube_radius=5):
    
    try:
        from enthought.mayavi import mlab
    except ImportError:
        log.info('Mlab is not available!')
        return
    
    X1 = numpy.asarray(X1)
    X2 = numpy.asarray(X2)
    X3 = numpy.asarray(X3)
    X4 = numpy.asarray(X4)
    
    X = numpy.row_stack((X1,X2))
    mlab.plot3d(X[:,0], X[:,1], X[:,2], [1,5], tube_radius=tube_radius, colormap='Reds')
    
    X = numpy.row_stack((X1,X3))
    mlab.plot3d(X[:,0], X[:,1], X[:,2], [1,5], tube_radius=tube_radius, colormap='Greens')
    
    X = numpy.row_stack((X1,X4))
    mlab.plot3d(X[:,0], X[:,1], X[:,2], [1,5], tube_radius=tube_radius, colormap='Blues')
    
    u1 = X2 - X1
    u2 = X3 - X1
    u3 = X4 - X1
    
    u1 = u1/numpy.sqrt(numpy.dot(u1,u1))
    u2 = u2/numpy.sqrt(numpy.dot(u2,u2))
    u3 = u3/numpy.sqrt(numpy.dot(u3,u3))
    
#    print numpy.dot(u1,u2)
#    print numpy.dot(u1,u3)
#    print numpy.dot(u2,u3)
    
    print 'angle1: %.2f' %numpy.degrees(numpy.arccos(numpy.dot(u1,u2)))
    print 'angle2: %.2f' %numpy.degrees(numpy.arccos(numpy.dot(u1,u3)))
    print 'angle3: %.2f' %numpy.degrees(numpy.arccos(numpy.dot(u2,u3)))
