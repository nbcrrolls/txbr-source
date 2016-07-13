import numpy
import numpy.linalg
import scipy.optimize
import powerseq
import transform

CONFORMAL_TRANSFORMATION = 'conformal'
AFFINE_TRANSFORMATION = 'affine'
PROJECTIVE_TRANSFORMATION = 'projective'
POLYNOMIAL_TRANSFORMATION = 'polynomial'

def polynomial_regression( x, points, orders ):
    '''This routine allows to perform polynomial regressions in multidimensional space. We
    note the dimension of the input space dim1, while the output space is dim2
    Variable x is a numpy array of shape (ndata,dim1) and represent the input parameters
    Variable points is a numpy array of shape (ndata,dim2) and contain the image of x.
    Variable orders is a numpy array of shape (dim1) and contains the polynomial order
    used for each component of x.
    Returns a polynom or an array of polynom of dimension dim2.'''
    
    x = numpy.asarray(x)
    points = numpy.asarray(points)
    orders = numpy.asarray(orders)

    (ndata,dim1) = x.shape

    if points.shape[0]!=ndata or orders.shape[0]!=dim1:
        raise ValueError

    if len(points.shape)==2:
        dim2 = points.shape[1]
    else:
        dim2 = 1
        points = numpy.resize(points,(ndata,dim2))

    nn = powerseq.numberOfTerms(orders,dim=dim1)

    coeffs_0 = numpy.ravel(numpy.eye(dim2,nn,1))
    
    coeffs = coeffs_0.copy()

    def residuals(coeffs_):
        res = numpy.zeros((dim2*ndata))
        for i in range(dim2):
           res[i*ndata:(i+1)*ndata] = points[:,i]-powerseq.poleval(coeffs_[i*nn:(i+1)*nn],orders,x)
        return res

    try:
        lsq = scipy.optimize.leastsq(residuals,coeffs)
        coeffs = lsq[0]
    except:
        coeffs = coeffs_0

    if dim2==1:
        return powerseq.Poly(coeffs,orders)
    else:
        return [powerseq.Poly(coeffs[i*nn:(i+1)*nn],orders) for i in range(dim2)]
    
    

def orthogonal_regression_2d( x, points ):
    '''This routine allows to perform a regression in a 2D space with orthogonal
    transforms.
    Variable x: a numpy array of shape (ndata,2). It represents the input parameters
    Variable points: a numpy array of shape (ndata,2). It contain the image of x.
    Returns: a tuple containing four variables (rot,theta,x1,x2). Variable rot is the 
    regression polynom; theta is the angle of the rotation and (x1,x2) represents
    the translation coefficients in the regression.'''

    ndata = x.size/2

    if points.size!=2*ndata:
        raise ValueError
    
    x.resize((ndata,2))
    points.resize((ndata,2))
    
    coeffs_0 = numpy.zeros(3)
    coeffs_0[:2] = numpy.mean(points-x, axis=0)
    
    try:
    
        R = numpy.zeros((2,2))
    
        def residuals(coeffs_):
            x1,x2,theta = coeffs_
            c,s = numpy.cos(theta),numpy.sin(theta)
            R[0,:] = [c,-s]
            R[1,:] = [s,c]
            rx = numpy.tensordot(R,x,(1,1))
            res = numpy.zeros((2*ndata))
            res[:ndata] = points[:,0] - x1 - rx[0,:]
            res[ndata:] = points[:,1] - x2 - rx[1,:]
            return res
    
        lsq = scipy.optimize.leastsq( residuals, coeffs_0 )
        
        x1,x2,theta = lsq[0]
        
    except:
        
        # In case there is a problem with the regression (e.g. not enough points)
        
        x1,x2,theta = coeffs_0
    

    c,s = numpy.cos(theta),numpy.sin(theta)
    
    rot = powerseq.Poly([x1,c,-s,x2,s,c],[1,1])

    return (rot,theta,x1,x2)



def cp2tform( input_points, base_points, transformation_type ):
    
    input_points = numpy.asarray(input_points)
    base_points = numpy.asarray(base_points)

    if input_points.shape!=base_points.shape:
        raise ValueError("Input and Base points must have the same shape.")

    if transformation_type==PROJECTIVE_TRANSFORMATION:    # Needs at least 4 pair of points
        
        n = input_points.shape[0]
        M = numpy.zeros((2*n,8))
        X = numpy.zeros(2*n)
        M[::2,0] = 1.0
        M[1::2,1] = 1.0
        M[::2,2:4] = input_points[:,:]
        M[1::2,4:6] = input_points[:,:]
        M[::2,6] = -base_points[:,0]*input_points[:,0]
        M[::2,7] = -base_points[:,0]*input_points[:,1]
        M[1::2,6] = -base_points[:,1]*input_points[:,0]
        M[1::2,7] = -base_points[:,1]*input_points[:,1]
        X[::2] = base_points[:,0]
        X[1::2] = base_points[:,1]
        M = numpy.dot(numpy.linalg.pinv(M),X)
        T = numpy.array([[M[2],M[3],M[0]],[M[4],M[5],M[1]],[M[6],M[7],1.0]])
        
        #return T
    
        T = numpy.array([[M[0],M[2],M[3]],[M[1],M[4],M[5]],[1.0,M[6],M[7]]])
        
        return transform.ProjectiveTransformation2D(T)
        
        


