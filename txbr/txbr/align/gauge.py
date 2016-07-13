import numpy
import numpy.linalg
import scipy.optimize
import txbr
import txbr.utilities
import util

log = txbr.log

def calc_rotation_axis( B, src_frame, dest_frame, ref_transform, smpl_transform, number_of_sections=5, \
                   doPlot=False, directory=None, basename=None, model=""):
    '''This routine is used to calculate the rotation axis of a series from a geometrical
    method.
    Variable B: is a numpy array (of shape (2,4,ntilt)) describing the projection map.
    Variable: ref_t_elements is a numpy array for the transfer matrix
    '''

    log.info('Calculate envelope with %i sections' %number_of_sections)
    
    ntilt = B.shape[2]
    
    src_frame = numpy.asarray(src_frame)      # Frame of the micrograph
    dest_frame = numpy.asarray(dest_frame)    # Frame of the reconstruction (sample)

    ref_transform_inv = ref_transform.inv()
    smpl_transform_inv = smpl_transform.inv()
    
    tc_r,mc_r = ref_transform.T, ref_transform.M
    tc_s,mc_s = smpl_transform.T, smpl_transform.M

    tc_r_inv,mc_r_inv = ref_transform_inv.T, ref_transform_inv.M
    tc_s_inv,mc_s_inv = smpl_transform_inv.T, smpl_transform_inv.M
    
    # Evaluate the projection map in the reference frame (as calculated by traces -- or initialized)
    
    T_r = B[:,0,:]
    P_r = B[:,1:,:]
    
    # Sample Rotation Contribution (Remove it from the general projection)
    
    T_r +=  numpy.squeeze(numpy.tensordot(P_r,tc_s_inv,axes=([1],[0])))
    P_r = numpy.swapaxes(numpy.tensordot(P_r,mc_s_inv,axes=([1],[0])),1,2)
    
    # Reference Contribution (i) to the right
    
    T_r += numpy.squeeze(numpy.tensordot(P_r,tc_r_inv,axes=([1],[0])))
    P_r = numpy.swapaxes(numpy.tensordot(P_r,mc_r_inv,axes=([1],[0])),1,2)
    
    # Reference Contribution (i) to the left
    
    T_r = numpy.tensordot(mc_r[:2,:2],T_r,axes=([1],[0]))
    T_r[0,:] += tc_r[0]
    T_r[1,:] += tc_r[1]
    P_r = numpy.swapaxes(numpy.tensordot(mc_r[:2,:2],P_r,axes=([1],[0])),1,2)

    # Source Frame - Help plotting the rays reaching the edges of the camera

    xmin,ymin = numpy.min(src_frame,axis=1)
    xmax,ymax = numpy.max(src_frame,axis=1)
    
    log.info('(xmin,xmax)=(%.2f,%.2f)' %(xmin,xmax))
    log.info('(ymin,ymax)=(%.2f,%.2f)' %(ymin,ymax))
    
    zmin = -max(xmax,ymax)/2.0
    zmax = -zmin
    
    log.info('(zmin,zmax)=(%.2f,%.2f)' %(zmin,zmax))
    
    # Destination Frame
    
    log.info("Redefine destination subvolume in the right intermediate reference frame:")
    
    dest_frame = ref_transform.forward(smpl_transform.forward(dest_frame.T)).T
  
    Xmin,Ymin,Zmin = numpy.min(dest_frame,axis=1)
    Xmax,Ymax,Zmax = numpy.max(dest_frame,axis=1)
        
    center_dest = numpy.mean(dest_frame,axis=1)
    
    log.info("(Xmin,Ymin,Zmin)=(%.1f,%.1f,%.1f)" %(Xmin,Ymin,Zmin))
    log.info("(Xmax,Ymax,Zmax)=(%.1f,%.1f,%.1f)" %(Xmax,Ymax,Zmax))
    
    half_width = (Xmax-Xmin)/2.0
    sections = numpy.linspace(Ymin,Ymax,num=number_of_sections)
    
    # define the family of tangents

    px = numpy.zeros((ntilt,number_of_sections,3))

    for isect in range(number_of_sections):
        px[:,isect,0] = T_r[0,:] + P_r[0,:,1]*sections[isect]
        px[:,isect,1] = P_r[0,:,0]
        px[:,isect,2] = P_r[0,:,2]

    tgts_x_min = px[:,:,:].copy()
    tgts_x_max = px[:,:,:].copy()

    tgts_x_min[:,:,0] = tgts_x_min[:,:,0] - xmin
    tgts_x_max[:,:,0] = tgts_x_max[:,:,0] - xmax

    tgts_x = numpy.row_stack((tgts_x_min,tgts_x_max))
    
    # Initialization values:

    a,b = half_width,half_width
    center_x,center_y = center_dest[0],center_dest[2]
    angle,eps = 0.0,1.0
    log.info('Initialization: %6.2f %6.2f  %6.2f %6.2f %6.2f %6.2f' %(a,b,center_x,center_y,angle,eps))

    xopt = numpy.array([a,b,center_x,center_y,angle,eps])
    xopt = numpy.tile(xopt,number_of_sections)

#    print ellipse_function
#    print xopt
#    print ellipse_der
#    print ellipse_hess
#    print tgts_x

    xopt = scipy.optimize.fmin_ncg(ellipse_function,xopt,ellipse_der,fhess=ellipse_hess,args=(tgts_x,),avextol=1.e-8,maxiter=50)

    xopt = numpy.resize(xopt,(number_of_sections,6))

    xopt[:,4] = numpy.where(numpy.abs(xopt[:,4])>numpy.pi/2.0,numpy.arctan(numpy.tan(xopt[:,4])),xopt[:,4])

    swapaxes = numpy.abs(xopt[:,0])<numpy.abs(xopt[:,1])
    xopt[:,4] = numpy.where(swapaxes,xopt[:,4]-numpy.sign(xopt[:,4])*numpy.pi/2.0,xopt[:,4])
    tmp = xopt[:,0].copy()
    xopt[:,0] = numpy.where(swapaxes,xopt[:,1],xopt[:,0])
    xopt[:,1] = numpy.where(swapaxes,tmp[:],xopt[:,1])
    
    ellipses = []

    log.info('  a         b         X0         Y0       angle      eps       Z0')
    
    for index,row in enumerate(xopt):
        
        a, b, x0, y0, angle, eps = row
        log.info('%6.2f    %6.2f    %6.2f    %6.2f    %6.2f    %6.2f    %6.2f' %(a,b,x0,y0,angle,eps,sections[index]))
        
        ellipse_left = util.Ellipse2D(x0,y0,a,b,angle)
        ellipse_right = util.Ellipse2D(x0,y0,a/eps,b/eps,angle)
        
        ellipses.append((ellipse_left,ellipse_right))
        
    # Do the plotting
    
    if doPlot:
        
        limits = (xmin,xmax,zmin,zmax,ymin,ymax)
        
        txbr.utilities.plotGaugeCharacterizarion( ellipses, tgts_x, limits, directory=directory, basename=basename, model=model )
 
    # Find the rotation axis

    x = xopt[:,2]
    y = xopt[:,3]
    z = numpy.asarray(sections,dtype='float')

    def residual(parameter, x, y, z):
        x0,y0,z0,u,v = parameter
        res = []
        res.append((y-y0)-(z-z0)*v)
        res.append((z-z0)*u-(x-x0))
        res.append((x-x0)*v-(y-y0)*u)
        return numpy.concatenate(res)

    params0 = [ center_dest[0], center_dest[1], center_dest[2], 0., 0. ]
    result = scipy.optimize.leastsq(residual, params0, (x,y,z))
    x0,y0,z0,u,v = result[0]

    # Rotation axis in the

    O = numpy.array([x0,z0,y0])
    N = numpy.array([u,1,v])
    N = N/numpy.sqrt(numpy.sum(N**2))

#    print 'Relative Frame'
#    print 'Axis Point: %s' %O
#    print 'Axis Vector: %s' %N

    # In the reference frame:

    OO = tc_r_inv + numpy.dot(mc_r_inv,O)
    NN = numpy.dot(mc_r_inv,N)
        
#    OO = smpl_transform.reverse(ref_transform.reverse(O))
#    NN = numpy.dot(mc_s_inv,numpy.dot(mc_r_inv,N))
    

#    print 'Absolute'
#    print 'Axis Point: %s' %OO
#    print 'Axis Vector: %s' %NN

    alpha = ((ymin+ymax)/2.0-OO[1])/NN[1]
    OO = OO + alpha*NN 

    return (OO,NN)





def ellipse_function(X,tgts_x):

    number_of_tilts = tgts_x.shape[0]/2
    number_of_sections = tgts_x.shape[1]

    phi = X[4::6]
    R = numpy.array([[numpy.cos(phi),-numpy.sin(phi)],[numpy.sin(phi),numpy.cos(phi)]])

    eps = numpy.ones((2*number_of_tilts,number_of_sections))
    eps[number_of_tilts:,:] *= X[5::6]
    eps = eps[:,:,numpy.newaxis]

    T1 = tgts_x*eps
    Y1 = numpy.row_stack((numpy.ones(number_of_sections),X[2::6],X[3::6]))
    q1 = numpy.tensordot(Y1,T1,axes=([0],[2]))
    q1 = numpy.diagonal(q1,axis1=0,axis2=2)

    T2 = numpy.tensordot(tgts_x[:,:,1:3],R,axes=([2],[0]))**2
    T2 = numpy.diagonal(T2,axis1=1,axis2=3)
    Y2 = numpy.row_stack((X[0::6]**2,X[1::6]**2))
    q2 = numpy.tensordot(Y2,T2,axes=([0],[1]))
    q2 = numpy.diagonal(q2,axis1=0,axis2=2)

    err = q1**2/q2-1

    Err= numpy.sum(err**2)

    return Err


def ellipse_der(X,tgts_x):

    der = numpy.zeros_like(X)

    number_of_tilts = tgts_x.shape[0]/2
    number_of_sections = tgts_x.shape[1]

    phi = X[4::6]

    R = numpy.array([[numpy.cos(phi),-numpy.sin(phi)],[numpy.sin(phi),numpy.cos(phi)]])
    dR = numpy.array([[-numpy.sin(phi),-numpy.cos(phi)],[numpy.cos(phi),-numpy.sin(phi)]])

    eps = numpy.ones((2*number_of_tilts,number_of_sections))
    eps[number_of_tilts:,:] *= X[5::6]
    eps = eps[:,:,numpy.newaxis]

    T1 = tgts_x*eps
    dT1 = tgts_x.copy()
    dT1[:number_of_tilts,:,:] = 0
    Y1 = numpy.row_stack((numpy.ones(number_of_sections),X[2::6],X[3::6]))
    q1 = numpy.tensordot(Y1,T1,axes=([0],[2]))
    q1 = numpy.diagonal(q1,axis1=0,axis2=2)
    dq1 = numpy.tensordot(Y1,dT1,axes=([0],[2]))
    dq1 = numpy.diagonal(dq1,axis1=0,axis2=2)

    T2 = numpy.tensordot(tgts_x[:,:,1:3],R,axes=([2],[0]))**2
    T2 = numpy.diagonal(T2,axis1=1,axis2=3)
    dT2 = 2.0*numpy.tensordot(tgts_x[:,:,1:3],dR,axes=([2],[0]))*numpy.tensordot(tgts_x[:,:,1:3],R,axes=([2],[0]))
    dT2 = numpy.diagonal(dT2,axis1=1,axis2=3)
    Y2 = numpy.row_stack((X[0::6]**2,X[1::6]**2))
    q2 = numpy.tensordot(Y2,T2,axes=([0],[1]))
    q2 = numpy.diagonal(q2,axis1=0,axis2=2)
    dq2 = numpy.tensordot(Y2,dT2,axes=([0],[1]))
    dq2 = numpy.diagonal(dq2,axis1=0,axis2=2)

    err = q1**2/q2-1;

    der_err_0 = -2.0*X[0::6]*T2[:,0,:]*q1**2/q2**2
    der_err_1 = -2.0*X[1::6]*T2[:,1,:]*q1**2/q2**2
    der_err_2 = 2.0*T1[:,:,1]*q1/q2
    der_err_3 = 2.0*T1[:,:,2]*q1/q2
    der_err_4 = -dq2*q1**2/q2**2
    der_err_5 = 2.0*dq1*q1/q2

    der[0::6] = numpy.sum(2*err*der_err_0,axis=0)
    der[1::6] = numpy.sum(2*err*der_err_1,axis=0)
    der[2::6] = numpy.sum(2*err*der_err_2,axis=0)
    der[3::6] = numpy.sum(2*err*der_err_3,axis=0)
    der[4::6] = numpy.sum(2*err*der_err_4,axis=0)
    der[5::6] = numpy.sum(2*err*der_err_5,axis=0)

    return der


def ellipse_hess(X,tgts_x):

    n = X.size
    hess = numpy.zeros((n,n))

    number_of_tilts = tgts_x.shape[0]/2
    number_of_sections = tgts_x.shape[1]

    phi = X[4::6]

    R = numpy.array([[numpy.cos(phi),-numpy.sin(phi)],[numpy.sin(phi),numpy.cos(phi)]])
    dR = numpy.array([[-numpy.sin(phi),-numpy.cos(phi)],[numpy.cos(phi),-numpy.sin(phi)]])

    eps = numpy.ones((2*number_of_tilts,number_of_sections))
    eps[number_of_tilts:,:] *= X[5::6]
    eps = eps[:,:,numpy.newaxis]

    T1 = tgts_x*eps
    dT1 = tgts_x.copy()
    dT1[:number_of_tilts,:,:] = 0
    Y1 = numpy.row_stack((numpy.ones(number_of_sections),X[2::6],X[3::6]))
    q1 = numpy.tensordot(Y1,T1,axes=([0],[2]))
    q1 = numpy.diagonal(q1,axis1=0,axis2=2)
    dq1 = numpy.tensordot(Y1,dT1,axes=([0],[2]))
    dq1 = numpy.diagonal(dq1,axis1=0,axis2=2)

    T2 = numpy.tensordot(tgts_x[:,:,1:3],R,axes=([2],[0]))**2
    T2 = numpy.diagonal(T2,axis1=1,axis2=3)
    dT2 = 2.0*numpy.tensordot(tgts_x[:,:,1:3],dR,axes=([2],[0]))*numpy.tensordot(tgts_x[:,:,1:3],R,axes=([2],[0]))
    dT2 = numpy.diagonal(dT2,axis1=1,axis2=3)
    dT2d2 = -2.0*numpy.tensordot(tgts_x[:,:,1:3],R,axes=([2],[0]))**2 + 2.0*numpy.tensordot(tgts_x[:,:,1:3],dR,axes=([2],[0]))**2
    dT2d2 = numpy.diagonal(dT2d2,axis1=1,axis2=3)
    Y2 = numpy.row_stack((X[0::6]**2,X[1::6]**2))
    q2 = numpy.tensordot(Y2,T2,axes=([0],[1]))
    q2 = numpy.diagonal(q2,axis1=0,axis2=2)
    dq2 = numpy.tensordot(Y2,dT2,axes=([0],[1]))
    dq2 = numpy.diagonal(dq2,axis1=0,axis2=2)
    dq2d2 = numpy.tensordot(Y2,dT2d2,axes=([0],[1]))
    dq2d2 = numpy.diagonal(dq2d2,axis1=0,axis2=2)

    err = q1**2/q2-1;

    der_err_0 = -2.0*X[0::6]*T2[:,0,:]*q1**2/q2**2
    der_err_1 = -2.0*X[1::6]*T2[:,1,:]*q1**2/q2**2
    der_err_2 = 2.0*T1[:,:,1]*q1/q2
    der_err_3 = 2.0*T1[:,:,2]*q1/q2
    der_err_4 = -dq2*q1**2/q2**2
    der_err_5 = 2.0*dq1*q1/q2

    hess_err_00 = der_err_0/X[0::6] + 8.0*X[0::6]**2*T2[:,0,:]**2*q1**2/q2**3
    hess_err_01 = 8.0*X[0::6]*T2[:,0,:]*X[1::6]*T2[:,1,:]*q1**2/q2**3
    hess_err_02 = -4.0*X[0::6]*T2[:,0,:]*T1[:,:,1]*q1/q2**2
    hess_err_03 = -4.0*X[0::6]*T2[:,0,:]*T1[:,:,2]*q1/q2**2
    hess_err_04 = -2.0*X[0::6]*dT2[:,0,:]*q1**2/q2**2 + 4.0*X[0::6]*T2[:,0,:]*dq2*q1**2/q2**3
    hess_err_05 = -4.0*X[0::6]*T2[:,0,:]*dq1*q1/q2**2

    hess_err_11 = der_err_1/X[1::6] + 8.0*X[1::6]**2*T2[:,1,:]**2*q1**2/q2**3
    hess_err_12 = -4.0*X[1::6]*T2[:,1,:]*T1[:,:,1]*q1/q2**2
    hess_err_13 = -4.0*X[1::6]*T2[:,1,:]*T1[:,:,2]*q1/q2**2
    hess_err_14 = -2.0*X[1::6]*dT2[:,1,:]*q1**2/q2**2 + 4.0*X[1::6]*T2[:,1,:]*dq2*q1**2/q2**3
    hess_err_15 = -4.0*X[1::6]*T2[:,1,:]*dq1*q1/q2**2

    hess_err_22 = 2.0*T1[:,:,1]**2/q2
    hess_err_23 = 2.0*T1[:,:,1]*T1[:,:,2]/q2
    hess_err_24 = -2.0*T1[:,:,1]*dq2*q1/q2**2
    hess_err_25 = 2.0*dT1[:,:,1]*q1/q2 + 2.0*T1[:,:,1]*dq1/q2

    hess_err_33 = 2.0*T1[:,:,2]*T1[:,:,2]/q2
    hess_err_34 = -2.0*T1[:,:,2]*dq2*q1/q2**2
    hess_err_35 = 2.0*dT1[:,:,2]*q1/q2 + 2.0*T1[:,:,2]*dq1/q2

    hess_err_44 = -dq2d2*q1**2/q2**2 + 2.0*dq2**2*q1**2/q2**3
    hess_err_45 = -2.0*dq2*dq1*q1/q2**2

    hess_err_55 = 2.0*dq1**2/q2

    hess[0::6,0::6] = 2*numpy.tensordot(der_err_0,der_err_0,axes=([0],[0])) + 2*numpy.sum(err*hess_err_00,axis=0)
    hess[0::6,1::6] = 2*numpy.tensordot(der_err_0,der_err_1,axes=([0],[0])) + 2*numpy.sum(err*hess_err_01,axis=0)
    hess[0::6,2::6] = 2*numpy.tensordot(der_err_0,der_err_2,axes=([0],[0])) + 2*numpy.sum(err*hess_err_02,axis=0)
    hess[0::6,3::6] = 2*numpy.tensordot(der_err_0,der_err_3,axes=([0],[0])) + 2*numpy.sum(err*hess_err_03,axis=0)
    hess[0::6,4::6] = 2*numpy.tensordot(der_err_0,der_err_4,axes=([0],[0])) + 2*numpy.sum(err*hess_err_04,axis=0)
    hess[0::6,5::6] = 2*numpy.tensordot(der_err_0,der_err_5,axes=([0],[0])) + 2*numpy.sum(err*hess_err_05,axis=0)

    hess[1::6,1::6] = 2*numpy.tensordot(der_err_1,der_err_1,axes=([0],[0])) + 2*numpy.sum(err*hess_err_11,axis=0)
    hess[1::6,2::6] = 2*numpy.tensordot(der_err_1,der_err_2,axes=([0],[0])) + 2*numpy.sum(err*hess_err_12,axis=0)
    hess[1::6,3::6] = 2*numpy.tensordot(der_err_1,der_err_3,axes=([0],[0])) + 2*numpy.sum(err*hess_err_13,axis=0)
    hess[1::6,4::6] = 2*numpy.tensordot(der_err_1,der_err_4,axes=([0],[0])) + 2*numpy.sum(err*hess_err_14,axis=0)
    hess[1::6,5::6] = 2*numpy.tensordot(der_err_1,der_err_5,axes=([0],[0])) + 2*numpy.sum(err*hess_err_15,axis=0)

    hess[2::6,2::6] = 2*numpy.tensordot(der_err_2,der_err_2,axes=([0],[0])) + 2*numpy.sum(err*hess_err_22,axis=0)
    hess[2::6,3::6] = 2*numpy.tensordot(der_err_2,der_err_3,axes=([0],[0])) + 2*numpy.sum(err*hess_err_23,axis=0)
    hess[2::6,4::6] = 2*numpy.tensordot(der_err_2,der_err_4,axes=([0],[0])) + 2*numpy.sum(err*hess_err_24,axis=0)
    hess[2::6,5::6] = 2*numpy.tensordot(der_err_2,der_err_5,axes=([0],[0])) + 2*numpy.sum(err*hess_err_25,axis=0)

    hess[3::6,3::6] = 2*numpy.tensordot(der_err_3,der_err_3,axes=([0],[0])) + 2*numpy.sum(err*hess_err_33,axis=0)
    hess[3::6,4::6] = 2*numpy.tensordot(der_err_3,der_err_4,axes=([0],[0])) + 2*numpy.sum(err*hess_err_34,axis=0)
    hess[3::6,5::6] = 2*numpy.tensordot(der_err_3,der_err_5,axes=([0],[0])) + 2*numpy.sum(err*hess_err_35,axis=0)

    hess[4::6,4::6] = 2*numpy.tensordot(der_err_4,der_err_4,axes=([0],[0])) + 2*numpy.sum(err*hess_err_44,axis=0)
    hess[4::6,5::6] = 2*numpy.tensordot(der_err_4,der_err_5,axes=([0],[0])) + 2*numpy.sum(err*hess_err_45,axis=0)

    hess[5::6,5::6] = 2*numpy.tensordot(der_err_5,der_err_5,axes=([0],[0])) + 2*numpy.sum(err*hess_err_55,axis=0)

    for i in range(number_of_sections):
        for j in range(i+1,number_of_sections):
            hess[6*i:6*i+6,6*j:6*j+6] = 0
            hess[6*j:6*j+6,6*i:6*i+6] = 0

    # Symmetryze the Hessian

    indiag = [range(n),range(n)]
    diag = numpy.diag(hess)

    hess = hess + hess.T
    hess[indiag] = diag[:]


    return hess









def test_derivatives():

    nx = 1000
    ntilt = 60
    number_of_sections = 2

    theta = numpy.linspace(-numpy.pi/3.0,numpy.pi/3.0,num=ntilt)

    tgts = numpy.zeros((ntilt,number_of_sections,3))

    for isect in range(number_of_sections):
        tgts[:,isect,0] = nx/2.0*(1-numpy.cos(theta))
        tgts[:,isect,1] = numpy.cos(theta)
        tgts[:,isect,2] = numpy.sin(theta)

    X = numpy.random.random_sample((6*number_of_sections))
    h = 0.0000001

    print 'Check First Derivatives...'

    der_app = numpy.zeros((6*number_of_sections))
    der = ellipse_der(X,tgts)

    for i in range(X.shape[0]):
        Y = X.copy()
        Y[i] += h
        der_app[i] = (ellipse_function(Y,tgts)-ellipse_function(X,tgts))/h
        print 'i=%-6i  app. der=% 10.5e  der=% 10.5e    rel. diff.=% 10.5e' %(i,der_app[i],der[i],numpy.abs((der_app[i]-der[i])/der[i]))

    print 'Check Hessians...'

    hess = ellipse_hess(X,tgts)

    for i in range(X.shape[0]):
        for j in range(i,X.shape[0]):
            Y = X.copy()
            Y[j] += h
            hess_app = (ellipse_der(Y,tgts)-ellipse_der(X,tgts))/h
            diff = numpy.where(hess_app[i]==hess[i,j],0.0,numpy.abs((hess_app[i]-hess[i,j])/hess[i,j]))
            print 'i=%-4i  j=%-4i  app. hess=% 10.5e  hess=% 10.5e    rel. diff.=% 10.5e' %(i,j,hess_app[i],hess[i,j],diff)


if __name__ == '__main__':

    test_derivatives()

