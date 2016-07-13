import sys
import os.path
import numpy
import scipy.optimize
import modl
import util
import util.powerseq
import txbr.utilities

from txbr import log

MAXITER = 50000
MAXFUN = 50000

def loadXYZ(directory,basename,grid=None):

    log.info('Load Markers for project (%s,%s).' %(directory,basename))

    if grid!=None:    # The user has delimited the boundaries manually on a grid of data

        filename = os.path.join(directory, "txbr-setup", basename + ".manual.fid")

        if not os.path.exists(filename):
            log.info('File %s does not exist!' %filename)
            sys.exit()

        model = modl.Model(filename)    # y and z are flipped
        model.loadFromFile()

        model.swapAxis([1,2])

        object = model.objects[0]

        XYZ_1 = numpy.zeros((0,3))
        XYZ_2 = numpy.zeros((0,3))

        n = object.numberOfContours()

        XYZ_1 = numpy.row_stack([object.contours[i].points for i in range(0,n,2)])
        XYZ_2 = numpy.row_stack([object.contours[i].points for i in range(1,n,2)])

        for i in range(3):
            if grid[i]==None: continue
            val = numpy.array(grid[i].split(','),dtype='int')
            XYZ_1[:,i] = val[XYZ_1[:,i].astype('int')]
            XYZ_2[:,i] = val[XYZ_2[:,i].astype('int')]

        xmax = model.xmax
        ymax = 1.25*val[int(model.ymax-1)]

        return ([XYZ_1,XYZ_2],xmax,ymax)

    else:    # Read the XYZ of the fiducial marks

        filename = os.path.join(directory, "txbr-align", basename + ".mod")

        if not os.path.exists(filename):
            log.info('File %s does not exist!' %filename)
            sys.exit()

        model = modl.Model(filename)

        model.loadFromFile()

        xmax = model.xmax
        ymax = model.ymax

        XYZ = model.points()

        return (XYZ,xmax,ymax)


def eval_bounds_1(XYZ):
    """This routine calculates the best plane crossing a set of XYZ points."""

    A = numpy.ones(XYZ.shape)
    A[:,0] = 1
    A[:,1] = XYZ[:,0]
    A[:,2] = XYZ[:,1]

    B = -XYZ[:,2]

    Ainv = numpy.linalg.pinv(A)

    sol = numpy.dot(Ainv,B)

    p = numpy.ones((4))
    p[:3] = sol[:]

    return p


def eval_bounds_2( XYZ, XY_range=(100,100) ):
    """This routine evaluates the position of two boundary planes to the set
    of XYZ particles."""

    def chi_boundaries(coeff):
        chi = 0;
        for xyz in XYZ:
            f1 = coeff[0] + coeff[2]*xyz[0] + coeff[3]*xyz[1] + xyz[2]
            f2 = coeff[1] + coeff[2]*xyz[0] + coeff[3]*xyz[1] + xyz[2]
            chi += pow(f1,2)*pow(f2,2)
        return chi

    nx,ny = XY_range
    n = max(*XY_range)

    x0 = [n, -n, 0, 0]
    x0 = [0, 0, 0, 0]
    coeff = scipy.optimize.fmin(chi_boundaries, x0, xtol=1e-6, ftol=1e-6, maxiter=MAXITER, maxfun=MAXFUN,full_output=0)

    bp = [coeff[0],coeff[2],coeff[3],1.0]
    tp = [coeff[1],coeff[2],coeff[3],1.0]

    if coeff[1]<=coeff[0]:
        bp[0] = coeff[1]
        tp[0] = coeff[0]

    bp = numpy.asarray(bp)
    tp = numpy.asarray(tp)

    return ( bp, tp )


def eval_poly_bounds( XYZ, order, test=False ):
    """This routine evaluates boundary surfaces as polynomial surfaces. Order is the order
    of the polynomial"""
    
    nn = util.powerseq.numberOfTerms(order,dim=2)
    orders = numpy.array([order,order])

    pw = numpy.asarray(util.powerseq.powerOrder(order,dim=2))

    print pw.shape
    print pw
    print pw.sum(axis=0)

    global index_eval
    index_eval = 0
    
    def chi_boundaries(coeff):
        
        global index_eval

        coeff1,coeff2 = numpy.zeros(nn),numpy.zeros(nn)
        
        coeff1[0],coeff2[0] = coeff[0],coeff[1]
        coeff1[1:nn],coeff2[1:nn] = coeff[2:nn+1],coeff[2:nn+1]

        if index_eval<5000:
            coeff1 = numpy.where(pw.sum(axis=0)<=3,coeff1,0.0)
            coeff2 = numpy.where(pw.sum(axis=0)<=3,coeff2,0.0)
            coeff[1:nn+1] = numpy.where(pw.sum(axis=0)<=3,coeff[1:nn+1],0.0)
        
        f1 = util.poleval(coeff1,orders,XYZ[:,:2]) + XYZ[:,2]
        f2 = util.poleval(coeff2,orders,XYZ[:,:2]) + XYZ[:,2]
        chi = f1**2*f2**2
        chi = chi.sum(axis=0)

        index_eval = index_eval + 1

        if test and index_eval%100==0: print "index: %6i    Xi: %.3e" %(index_eval,chi)
        
        return chi

    x0 = numpy.zeros(nn+1)

    x0[0], x0[1] = 10, -10

    coeff = scipy.optimize.fmin(chi_boundaries, x0, xtol=1e-6, ftol=1e-6, maxiter=MAXITER, maxfun=MAXFUN,full_output=0)

    bp,tp = numpy.zeros(nn+1),numpy.zeros(nn+1)

    bp[0],tp[0] = coeff[0],coeff[1]
    bp[1:nn],tp[1:nn] = coeff[2:nn+1],coeff[2:nn+1]
    bp[nn],tp[nn] = 1.0,1.0
    
    bp = numpy.asarray(bp)
    tp = numpy.asarray(tp)

    return ( bp, tp )


def evalFlatteningTransformation( XYZ, order, center, doPlot=False, test=False ):
    '''Evaluate polynomial transformation (shear like along z) that
    transforms the warped volume into a flat one, and reverse.'''
    
    ( bp, tp ) = eval_poly_bounds( XYZ, order, test=test )
    
    variables = ['X','Y','Z']
    
    nn_2d = util.powerseq.numberOfTerms(order,dim=2)

    powers_2d = numpy.array(util.powerOrder(order,dim=2)).T
    powers_2d = numpy.column_stack((powers_2d,numpy.zeros(nn_2d)))
    powers_2d = numpy.row_stack((powers_2d,numpy.array([0,0,1])))
    
    nn_3d = util.powerseq.numberOfTerms(order,dim=3)
    
    powers_3d = numpy.array(util.powerOrder(order,dim=3)).T
    
    # Calculate the transformation to go from warp to flat
    
    z_tf = util.PolyM(variables, bp, powers_2d)
    
    coeffx = numpy.eye(1,nn_3d,1)
    coeffy = numpy.eye(1,nn_3d,2)
    coeffz = z_tf.extract_coefficients(powers_3d)
    
    M = numpy.row_stack((coeffx,coeffy,coeffz))
    M[:,0] = center[:]
    
    warp_to_flat_tf = util.PolynomialTransformation3D( order, M, with_center=True )
    
    # Calculate the transformation to go from flat to warp
    
    bp_inv = bp.copy()
    bp_inv[:-1] = -bp_inv[:-1]

    z_tf_inv = util.PolyM(variables, bp_inv, powers_2d)
    
    coeffx_inv = numpy.eye(1,nn_3d,1)
    coeffy_inv = numpy.eye(1,nn_3d,2)
    coeffz_inv = z_tf_inv.extract_coefficients(powers_3d)
    
    M_inv = numpy.row_stack((coeffx_inv,coeffy_inv,coeffz_inv))
    M_inv[:,0] = center[:]
    
    flat_to_warp_tf = util.PolynomialTransformation3D( order, M_inv, with_center=True )
    
    if doPlot: # Plot the beads, their flattened counterparts as well as some surfaces
        
        minimum = numpy.min(XYZ,axis=0).astype('int')
        maximum = numpy.max(XYZ,axis=0).astype('int')
        
        step = (maximum-minimum)/10.0
        step = step.astype('int')
        
        x, y = numpy.mgrid[minimum[0]:maximum[0]:step[0],minimum[1]:maximum[1]:step[1]]
        shape = x.shape
        
        scale = ( maximum[0], maximum[1], min(maximum[0],maximum[1])/3.0 )
        
        # Generate the warped surfaces
        
        powers = numpy.array(util.powerOrder(order,dim=2)).T
        f1_wrp =  util.PolyM(['X','Y'], -bp[:-1], powers)
        f2_wrp =  util.PolyM(['X','Y'], -tp[:-1], powers)
        
        z1_wrp = f1_wrp.eval(numpy.column_stack((x.ravel(),y.ravel())))
        z2_wrp = f2_wrp.eval(numpy.column_stack((x.ravel(),y.ravel())))
        
        z1_wrp.resize(shape)
        z2_wrp.resize(shape)
        
        surf1 = [ (x,y,z1_wrp), (x,y,z2_wrp) ]
        
        # Now generate the flat surfaces
                
        XYZ_f = warp_to_flat_tf.forward(XYZ)

        z1_ = f1_wrp.eval(center[:2])
        z2_ = f2_wrp.eval(center[:2])
                
        XYZ_f1 = warp_to_flat_tf.forward([center[0],center[1],z1_])
        XYZ_f2 = warp_to_flat_tf.forward([center[0],center[1],z2_])
        
        z1_flat = numpy.zeros_like(z1_wrp)
        z2_flat = numpy.zeros_like(z2_wrp)
        
        z1_flat += XYZ_f1[2]
        z2_flat += XYZ_f2[2]
        
        surf2 = [ (x,y,z1_flat), (x,y,z2_flat) ]
        
        # Do the plottings

        txbr.utilities.plot3DMarkers( XYZ[:,0], XYZ[:,1], XYZ[:,2], scale=scale, surfaces=surf1 )
        txbr.utilities.plot3DMarkers( XYZ_f[:,0], XYZ_f[:,1], XYZ_f[:,2], scale=scale, surfaces=surf2 )
        
    return ( warp_to_flat_tf, flat_to_warp_tf )
        
    


def eval_pitch(directory,basename,grid=None,doPlot=True,doSave=False):
    '''Evaluates the width of the specimen sample.'''

    if grid==None:
        log.info('Pitch calculation from fiducial markers.')
    else:
        log.info('Pitch calculation from manual tracks.')

    (XYZ,xmax,ymax) = loadXYZ(directory,basename,grid=grid)

    if grid==None:
        [p_b,p_t] = eval_bounds_2(XYZ)
    else:
        XYZ_b = XYZ[0]
        XYZ_t = XYZ[1]
        p_b = eval_bounds_1(XYZ_b)
        p_t = eval_bounds_1(XYZ_t)

    if doPlot:
        if grid==None:
            plot([p_b,p_t],pts=[XYZ],xmax=xmax,ymax=ymax)
        else:
            plot([p_b,p_t],pts=[XYZ_b,XYZ_t],xmax=xmax,ymax=ymax)

    if doSave:
        save_planes(directory,basename,[p_b,p_t])

    return [p_b,p_t]


def save_planes(directory,basename,planes):
    '''Save the boundary plane coefficients in a text file. The text file is located
    in the txbr-setup directory and called basename.pitch.txt.
    For a plane described by equation a + b*x + c*y + d*z=0, corresponds a line with
    coefficents:
    a    b    c    d'''

    bp,tp = planes

    filename = os.path.join(directory, "txbr-setup", basename + ".pitch.txt")

    f = open(filename, 'w')
    f.write("%f    %f    %f    %f\n" % (bp[0],bp[1],bp[2],bp[3]))
    f.write("%f    %f    %f    %f\n" % (tp[0],tp[1],tp[2],tp[3]))
    f.close()


def plot(planes,pts=None,xmin=1,xmax=100,ymin=1,ymax=100):
    '''Plot the tomopitch planes'''

    log.info('(xmax,ymax)=(%s,%s)' %(xmax,ymax))

    xinc = int((xmax-xmin)/10.0)
    yinc = int((ymax-ymin)/10.0)

    if xinc==0: xinc = 1
    if yinc==0: yinc = 1

    from enthought.mayavi import mlab

    fig = mlab.figure()

    for p in planes:
        def f(x,y):
            return - p[0] - p[1]*x - p[2]*y
        x, y = numpy.mgrid[xmin:xmax:xinc,ymin:ymax:yinc]
        mlab.surf(x,y,f)

    if pts!=None:
        for p in pts:
            mlab.points3d(p[:,0],p[:,1],p[:,2],scale_factor=20.0)


    mlab.show()

