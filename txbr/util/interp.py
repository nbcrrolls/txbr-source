import sys
import logging
import getopt
import numpy
import pylab
import scipy.optimize

from plot import plot_contour


def f(X,constraints,verbose=False):

    shape0 = X.shape
    n1,n2 = constraints.shape

    X.resize((n1,n2))

    X_ = numpy.where(numpy.isfinite(constraints),constraints,X)

    f = numpy.sum(X_**2)
    f -= numpy.sum(X_[0:-1,:]*X_[1:,:])/2.0
    f -= numpy.sum(X_[:,0:-1]*X_[:,1:])/2.0

    X.resize(shape0)

    if verbose: print f

    return f


def dfdX(X,constraints,verbose=False):

    shape0 = X.shape
    n1,n2 = constraints.shape

    X.resize((n1,n2))

    X_ = numpy.where(numpy.isfinite(constraints),constraints,X)

    dfdX = 2.0*X_
    dfdX[0:-1,:] -= X_[1:,:]/2.0
    dfdX[1:,:] -= X_[0:-1,:]/2.0
    dfdX[:,0:-1] -= X_[:,1:]/2.0
    dfdX[:,1:] -= X_[:,0:-1]/2.0

    dfdX = numpy.where(numpy.isnan(constraints),dfdX,0)

    dfdX.resize((n1*n2))

    X.resize(shape0)

    return dfdX


def d2fdX2(X,constraints,verbose=False):

    shape0 = X.shape

    n1,n2 = constraints.shape

    X.resize((n1,n2))

    d2fdX2 = numpy.zeros((n1,n2,n1,n2))

    for position1,value1 in numpy.ndenumerate(X):
        if numpy.isnan(constraints[position1]):
            d2fdX2[position1][position1] = 2.0

    for i in range(n1-1):
        for j in range(n2):
            if numpy.isnan(constraints[i,j]) and numpy.isnan(constraints[i+1,j]):
                d2fdX2[i,j,i+1,j] -= 1.0/2.0
                d2fdX2[i+1,j,i,j] -= 1.0/2.0
    for i in range(n1):
        for j in range(n2-1):
            if numpy.isnan(constraints[i,j]) and numpy.isnan(constraints[i,j+1]):
                d2fdX2[i,j,i,j+1] -= 1.0/2.0
                d2fdX2[i,j+1,i,j] -= 1.0/2.0

    d2fdX2.resize((n1*n2,n1*n2))

    X.resize(shape0)

    return d2fdX2


def inter2D( constraints, mode, init=None ):
    '''This routine performs an interpolation scheme from a set of scattered points (ideally contours).
    Those points corresponds to all the positions where the 2D map constraints is defined (Nan where
    we need to provide an interpolation).ss
    mode: variable that specify how interpolation is calculated, either through the minimization
    of a certain energy (defined in this module as f) or through a recursive approach.
    '''

    if init==None or init.shape!=constraints.shape:
        map0 = numpy.zeros_like(constraints)
    else:
        map0 = init

    map0 = numpy.where(numpy.isfinite(constraints),constraints,map0)

    if mode=='minimization':

        EPS = 1.e-6
        MAX_ITER = 50

        map = scipy.optimize.fmin_ncg(f,map0,dfdX,args=(constraints,True),avextol=EPS,maxiter=MAX_ITER)
        map.resize(constraints.shape)

    elif mode=='recursion':

        eps = 1e-8
        diff = 100.0

        map = map0.copy()

        while diff>eps:

            # 4 neighbours
            map[1:-1,1:-1] = (map0[0:-2,1:-1] + map0[2:,1:-1] + map0[1:-1,0:-2] + map0[1:-1,2:])/4.0
            
            # 3 neighbours
            map[0,1:-1] = (map0[1,1:-1] + map0[0,0:-2] + map0[0,2:])/4.0
            map[-1,1:-1] = (map0[-2,1:-1] + map0[-1,0:-2] + map0[-1,2:])/4.0
            map[1:-1,0] = (map0[1:-1,1] + map0[0:-2,0] + map0[2:,0])/4.0
            map[1:-1,-1] = (map0[1:-1,-2] + map0[0:-2,-1] + map0[2:,-1])/4.0
            
            # 2 neighbours
            map[0,0] = (map0[0,1] + map0[1,0])/4.0
            map[0,-1] = (map0[0,-2] + map0[1,-1])/4.0
            map[-1,0] = (map0[-1,1] + map0[-2,0])/4.0
            map[-1,-1] = (map0[-2,-1] + map0[-1,-2])/4.0

            map = numpy.where(numpy.isfinite(constraints),constraints,map)
            diff = numpy.max(numpy.abs(map-map0))
            #print numpy.argmax(numpy.abs(map-map0))

            map0 = map.copy()
            print 'diff=%e    F=%e' %(diff,f(map0,constraints))

    return map


def mapSurface(shape,pts,init=None,mode='minimization'):
    '''Map the surface.
    '''

    # Track the points within the shape

    indices = numpy.where( (pts[:,0]>0) & (pts[:,0]<=shape[0]) & \
                           (pts[:,1]>0) & (pts[:,1]<=shape[1]))


    x = pts[indices,0].astype('int')-1
    y = pts[indices,1].astype('int')-1
    u = pts[indices,2]

    constraints = numpy.empty(shape)
    constraints[:,:] = numpy.nan
    constraints[x,y] = u

    map = inter2D(constraints,mode,init=init)

    return map


def mapSurfaceFromModel(file,indexOfObject):
    '''Map the surface.
    '''
    import modl

    f = modl.Model(file)
    f.loadFromFile()

    shape = (f.xmax,f.ymax)
    pts = f.objects[indexOfObject].points()

    return mapSurface(shape,pts)


def initRandomConstraints(n1,n2,npoints):
    '''Initialize a map of random constraints. Map size is (n1,n2).
    Variable npoints is the number of constraints.
    '''

    x = numpy.random.random_integers(n1-1,size=(npoints))
    y = numpy.random.random_integers(n2-1,size=(npoints))

    u = numpy.random.random_sample((npoints))

    constraints = numpy.empty((n1,n2))
    constraints[:,:] = numpy.nan
    constraints[x,y] = u[:]

    return constraints


def test_derivatives():
    '''Test the derivatives and hessians of the function to minimize for doing the
    interpolation scheme.
    '''

    n1,n2,npoints = 4,4,2
    constraints = initRandomConstraints(n1,n2,npoints)

    X = numpy.random.random((n1,n2))

    from tester import test_derivatives
    test_derivatives(f,dfdX,X,args=(constraints,),d2fdX2=d2fdX2,loglevel=logging.INFO)


def testRandom(n1=50,n2=50,npoints=10,mode='minimization'):
    '''Test the interpolation process for an image of size (n1,n2) with
    some contraints. Variable npoints is the number of constraints.
    '''

    constraints = initRandomConstraints(n1,n2,npoints)
    map = inter2D(constraints,mode)

    return map


if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], "dm",["test-derivatives"])
    except getopt.GetoptError, err:
        sys.exit(2)

    for option,value in opts:
        if option in ('-m'):
            map0 = mapSurfaceFromModel('/Users/sph/Electron_Tomography/txbr/data/126/sum_126.flip.volmod',0)
            map1 = mapSurfaceFromModel('/Users/sph/Electron_Tomography/txbr/data/126/sum_126.flip.volmod',1)
            plot_contour(map0.T,title='Index 0 Surface Contour for Specimen 126')
            plot_contour(map1.T,title='Index 1 Surface Contour for Specimen 126')
            pylab.show()
        elif option in ('-d','--test-derivatives'):
            test_derivatives()

    if len(opts)==0:
        npoints = 10
        plot_contour(testRandom(npoints=npoints,mode='minimization'),title='%i random scattered point-constraint' %npoints)
        pylab.show()
