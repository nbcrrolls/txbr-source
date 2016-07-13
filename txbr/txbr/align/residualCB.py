import numpy, numpy.random, swiginac
import residual, residualSupp
import util

from residual import log, \
                     ONLY_FEATURE_DEPENDANT, \
                     ONLY_TILT_DEPENDANT, \
                     FEATURE_AND_TILT_DEPENDANT, \
                     FEATURE_AND_TILT_INDEPENDANT

class ResidualCB(residual.Residual):
    '''
    This class helps to generate the residual error during the bundle adjustment procedure
    within a constant beam approximation!
    '''

    def __init__( self, ntilt = 2, npatch = 2, n1 = 0, n2 = 1, n3 = 0, n4 = 0, \
                  fix_a = False, fix_b = True, fix_c = True, fix_d = False, \
                  prepare = True, enforceTangencyConstraint = False, \
                  structureCharacteristicLength = 100  ):

        residual.Residual.__init__( self, ntilt, npatch, n1, n2, n3, n4, 
                                    fix_a=fix_a, fix_b=fix_b, fix_c=fix_c, fix_d=fix_d, \
                                    prepare=False, enforceTangencyConstraint=enforceTangencyConstraint, \
                                    structureCharacteristicLength = structureCharacteristicLength )

        self.mode2r = ONLY_TILT_DEPENDANT    # Rotation part
        self.mode2b = FEATURE_AND_TILT_INDEPENDANT    # Beam part

        if prepare:
            self.fixVariables(fix_a,fix_b,fix_c,fix_d)
            
            
    def numberOfVariables(self, core=False):
        '''Returns the total number of variables defined in that residual.
        '''

        n = 0

        if core:
            nmode = { ONLY_FEATURE_DEPENDANT:1, \
                      ONLY_TILT_DEPENDANT:1, \
                      FEATURE_AND_TILT_DEPENDANT:1, \
                      FEATURE_AND_TILT_INDEPENDANT:1 }
        else:
            nmode = { ONLY_FEATURE_DEPENDANT:self.npatch, \
                      ONLY_TILT_DEPENDANT:self.ntilt, \
                      FEATURE_AND_TILT_DEPENDANT:self.ntilt*self.npatch, \
                      FEATURE_AND_TILT_INDEPENDANT:1 }

        n += 3*self.nn1*nmode[self.mode1]
        n += 12*nmode[self.mode2r] + 2*self.nn2*nmode[self.mode2b]
        n += 2*self.nn3*nmode[self.mode3]
        n += self.nn4*nmode[self.mode4]

        return n


    def numberOfFreeVariables(self,core=False):
        '''Returns the total number of free variables defined in that residual.
        '''

        n = 0

        if core:
            nmode = { ONLY_FEATURE_DEPENDANT:1, \
                      ONLY_TILT_DEPENDANT:1, \
                      FEATURE_AND_TILT_DEPENDANT:1, \
                      FEATURE_AND_TILT_INDEPENDANT:1 }
        else:
            nmode = { ONLY_FEATURE_DEPENDANT:self.npatch, \
                      ONLY_TILT_DEPENDANT:self.ntilt, \
                      FEATURE_AND_TILT_DEPENDANT:self.ntilt*self.npatch, \
                      FEATURE_AND_TILT_INDEPENDANT:1 }

        if not self.fix_a: n += 3*self.nn1*nmode[self.mode1]
        if not self.fix_b: n += 12*nmode[self.mode2r] + 2*self.nn2*nmode[self.mode2b]
        if not self.fix_c: n += 2*self.nn3*nmode[self.mode3]
        if not self.fix_d: n += self.nn4*nmode[self.mode4]

        return n
     
                
    def unfold_var2(self,var2):
        '''From single beam model to multiple beam model. Add some constraints
        on the fixed points.'''

        var2 = numpy.asarray(var2)

        var_l = numpy.resize(var2[:12*self.ntilt],(3,4,self.ntilt))
        var_nl = numpy.resize(var2[12*self.ntilt:],(2,self.nn2))

        # Compose the models

        var2_ = numpy.zeros((2,self.nn2,self.ntilt))

        variables = ['X','Y','Z']
        powers = numpy.array(util.powerOrder(1,dim=3)).T

        for itilt in range(self.ntilt):

            rot_x = util.PolyM(variables,var_l[0,:4,itilt],powers)
            rot_y = util.PolyM(variables,var_l[1,:4,itilt],powers)
            rot_z = util.PolyM(variables,[0.0,0.0,0.0,1.0],powers)
            rot_poly = {'X':rot_x,'Y':rot_y,'Z':rot_z}

            coefficients_x = var_nl[0,:]
            coefficients_y = var_nl[1,:]

            pp2 = numpy.array(self.pp2).T

            polyb_x = util.PolyM(variables,coefficients_x,pp2)
            polyb_x = polyb_x.compose(rot_poly)

            polyb_y = util.PolyM(variables,coefficients_y,pp2)
            polyb_y = polyb_y.compose(rot_poly)

            var2_[0,:,itilt] = polyb_x.extract_coefficients(pp2)
            var2_[1,:,itilt] = polyb_y.extract_coefficients(pp2)
            
        return var2_.ravel().tolist()


    def unfold_var2_with_scaling(self, var2, scaling):
        '''Unfold the projection coefficients (with the scaling parameters)'''

        scaling = numpy.asarray(scaling)    # Four coefficients (1,lam_X,lam_Y,lam_Z)
        var2 = numpy.asarray(var2)

        R = numpy.resize(var2[:12*self.ntilt],(3,4,self.ntilt))
        
        unfolded_scaling = numpy.tensordot(scaling[1:],R,axes=((0),(0)))    # of shape (4,ntilt)
        unfolded_scaling[0,:] += 1
        unfolded_scaling = unfolded_scaling.ravel().tolist()
        
        unfolded_var2 = self.unfold_var2(var2)

        return ( unfolded_var2, unfolded_scaling )


    def unfold_var2_a( self, var_b):     
        '''In the case we want to put a constraint on the fixed points.'''
        
        if numpy.all(~self.t_constant): return var_b[:12*self.ntilt]
        
        var_b = numpy.asarray(var_b[:12*self.ntilt])
        var_b.resize((3,4,self.ntilt))
    
        pp2 = numpy.asarray(self.pp2)
        pp2 = pp2[:,1:4]

        X0 = self.fixed_origin[0]**pp2[0,:]*self.fixed_origin[1]**pp2[1,:]*self.fixed_origin[2]**pp2[2,:]

        var2_fixed = var_b.copy()

        var2_fixed[0,0,:] = self.fixed_origin[0]
        var2_fixed[1,0,:] = self.fixed_origin[1]
        var2_fixed[2,0,:] = self.fixed_origin[2]
        
        var2_fixed[:,0,:] -= numpy.tensordot(var_b[:,1:,:],X0,axes=((1),(0)))

        indices = numpy.argwhere(~self.t_constant).ravel()    # indices on no constraint
        var2_fixed[indices] = var_b[indices]
        
        return var2_fixed.ravel().tolist()
    
    
    def unfold_var2_b( self, var_b ):
        
        var2_b = numpy.asarray(var_b)
        
        return var2_b[12*self.ntilt:]


    def fold_derivatives(self, der):
        '''Fold the derivatives to take into account possible constraints on the fixed
        points.'''

        if numpy.all(~self.t_constant): return der

        offset = [0]*5
        
        if self.fix_a: 
            l1_ = 0
        else:
            l1_ = self.l1
            
        if self.fix_b: 
            l2_ = 0
            l2a_,l2b_ = 0,0
        else:
            l2_ = self.l2
            l2a_,l2b_ = self.l2r,self.l2b
            
        if self.fix_c: 
            l3_ = 0
        else:
            l3_ = self.l3
            
        if self.fix_d: 
            l4_ = 0
        else:
            l4_ = self.l4
            
        offset[0] = l1_
        offset[1] = offset[0] + l2a_
        offset[2] = offset[1] + l2b_
        offset[3] = offset[2] + l3_
        offset[4] = offset[3] + l4_
        
        der_a = der[:offset[0]]
        der_ba = der[offset[0]:offset[1]]
        der_bb = der[offset[1]:offset[2]]
        der_c = der[offset[2]:offset[3]]
        der_d = der[offset[3]:offset[4]]
        
        indices = numpy.argwhere(self.t_constant).ravel()  # index array of the additional constraints

        pp2 = numpy.asarray(self.pp2)
        pp2 = pp2[:,1:4]

        X0 = self.fixed_origin[0]**pp2[0,:]*self.fixed_origin[1]**pp2[1,:]*self.fixed_origin[2]**pp2[2,:]

        # Positions of the patches

        der_a_ = der_a
        
        # Rotational part of the projection map
        
        der_ba_ = numpy.resize(der_ba,(3,4,self.ntilt))
                        
        der_ba_fixed = der_ba_[indices]
                
        der_ba_fixed[:,1:,:] -= numpy.rollaxis(der_ba_fixed[:,0,:,numpy.newaxis]*X0,2,1)
        der_ba_fixed[:,0,:] = 0
        
        der_ba_[indices] = der_ba_fixed
               
        der_ba_ = der_ba_.ravel()
        
        # Constant beam part of the projection map
        
        der_bb_ = der_bb
        
        # c
        
        der_c_ = der_c
        
        # d
        
        der_d_ = der_d
        
        # end

        return numpy.concatenate((der_a_, der_ba_, der_bb_, der_c_, der_d_))
    
    
    def fold_hessians(self, hess, der):
        '''From general frame to orthogonal frame'''

        if numpy.all(~self.t_constant): return hess

        offset = [0]*5
        
        if self.fix_a: 
            l1_ = 0
        else:
            l1_ = self.l1
            
        if self.fix_b: 
            l2_ = 0
            l2a_,l2b_ = 0,0
        else:
            l2_ = self.l2
            l2a_,l2b_ = self.l2r,self.l2b
            
        if self.fix_c: 
            l3_ = 0
        else:
            l3_ = self.l3
            
        if self.fix_d: 
            l4_ = 0
        else:
            l4_ = self.l4
            
        offset[0] = l1_
        offset[1] = offset[0] + l2a_
        offset[2] = offset[1] + l2b_
        offset[3] = offset[2] + l3_
        offset[4] = offset[3] + l4_
        
        der_a = der[:offset[0]]
        der_ba = der[offset[0]:offset[1]]
        der_bb = der[offset[1]:offset[2]]
        der_c = der[offset[2]:offset[3]]
        der_d = der[offset[3]:offset[4]]
        
        h_a_a = hess[:offset[0],:offset[0]]
        h_a_ba = hess[:offset[0],offset[0]:offset[1]]
        h_a_bb = hess[:offset[0],offset[1]:offset[2]]
        h_a_c = hess[:offset[0],offset[2]:offset[3]]
        h_a_d = hess[:offset[0],offset[3]:offset[4]]

        h_ba_a = hess[offset[0]:offset[1],:offset[0]]
        h_ba_ba = hess[offset[0]:offset[1],offset[0]:offset[1]]
        h_ba_bb = hess[offset[0]:offset[1],offset[1]:offset[2]]
        h_ba_c = hess[offset[0]:offset[1],offset[2]:offset[3]]
        h_ba_d = hess[offset[0]:offset[1],offset[3]:offset[4]]

        h_bb_a = hess[offset[1]:offset[2],:offset[0]]
        h_bb_ba = hess[offset[1]:offset[2],offset[0]:offset[1]]
        h_bb_bb = hess[offset[1]:offset[2],offset[1]:offset[2]]
        h_bb_c = hess[offset[1]:offset[2],offset[2]:offset[3]]
        h_bb_d = hess[offset[1]:offset[2],offset[3]:offset[4]]

        h_c_a = hess[offset[2]:offset[3],:offset[0]]
        h_c_ba = hess[offset[2]:offset[3],offset[0]:offset[1]]
        h_c_bb = hess[offset[2]:offset[3],offset[1]:offset[2]]
        h_c_c = hess[offset[2]:offset[3],offset[2]:offset[3]]
        h_c_d = hess[offset[2]:offset[3],offset[3]:offset[4]]
        
        h_d_a = hess[offset[3]:offset[4],:offset[0]]
        h_d_ba = hess[offset[3]:offset[4],offset[0]:offset[1]]
        h_d_bb = hess[offset[3]:offset[4],offset[1]:offset[2]]
        h_d_c = hess[offset[3]:offset[4],offset[2]:offset[3]]
        h_d_d = hess[offset[3]:offset[4],offset[3]:offset[4]]
        
        # Include the constraints

        indices = numpy.argwhere(self.t_constant).ravel()  # indices where there is a constraint

        pp2 = numpy.asarray(self.pp2)
        pp2 = pp2[:,1:4]

        X0 = self.fixed_origin[0]**pp2[0,:]*self.fixed_origin[1]**pp2[1,:]*self.fixed_origin[2]**pp2[2,:]
        
        # a
        
        h_a_a_ = h_a_a

        h_a_ba_ = numpy.resize(h_a_ba,(l1_,3,4,self.ntilt))
        h_a_ba_fixed = h_a_ba_[:,indices,:,:]
        
        h_a_ba_fixed[:,:,1:,:] -= numpy.rollaxis(h_a_ba_fixed[:,:,0,:,numpy.newaxis]*X0,3,2)
        h_a_ba_fixed[:,:,0,:] = 0
        h_a_ba_[:,indices,:,:] = h_a_ba_fixed
        h_a_ba_.resize(l1_,l2a_)
        
        h_a_bb_ = h_a_bb
        
        h_a_c_ = h_a_c
        
        h_a_d_ = h_a_d
        
        # b

        h_ba_a_ = h_a_ba_.T

        
        h_ba_ba = numpy.resize(h_ba_ba,(3,4,self.ntilt,3,4,self.ntilt))
        
        h_ba_ba_ = h_ba_ba.copy()
        
        i1a = numpy.ix_(indices,[0],range(self.ntilt),range(3),range(4),range(self.ntilt))
        i1b = numpy.ix_(indices,range(1,4),range(self.ntilt),range(3),range(4),range(self.ntilt))
        u1 = h_ba_ba[i1a]
        h_ba_ba_[i1b] -= numpy.rollaxis(u1[:,0,:,:,:,:,numpy.newaxis]*X0,5,1)
        
        i2a = numpy.ix_(range(3),range(4),range(self.ntilt),indices,[0],range(self.ntilt))
        i2b = numpy.ix_(range(3),range(4),range(self.ntilt),indices,range(1,4),range(self.ntilt))
        u2 = h_ba_ba[i2a]
        h_ba_ba_[i2b] -= numpy.rollaxis(u2[:,:,:,:,0,:,numpy.newaxis]*X0,5,4)
        
        i12a = numpy.ix_(indices,[0],range(self.ntilt),indices,[0],range(self.ntilt))
        i12b = numpy.ix_(indices,range(1,4),range(self.ntilt),indices,range(1,4),range(self.ntilt))
        u12a = h_ba_ba[i12a]
        h_ba_ba_[i12b] += numpy.rollaxis(numpy.rollaxis( \
                  u12a[:,0,:,:,0,:,numpy.newaxis,numpy.newaxis]*numpy.outer(X0,X0) \
                  ,5,1),5,4)

        h_ba_ba_[i1a] = 0.0
        h_ba_ba_[i2a] = 0.0

        h_ba_ba_.resize((l2a_,l2a_))
        
        
        h_ba_bb_ = numpy.resize(h_ba_bb,(3,4,self.ntilt,l2b_))
        h_ba_bb_fixed = h_ba_bb_[indices,:,:,:]
        
        h_ba_bb_fixed[:,1:,:,:] -= numpy.rollaxis(h_ba_bb_fixed[:,0,:,:,numpy.newaxis]*X0,3,1)
        h_ba_bb_fixed[:,0,:,:] = 0
        h_ba_bb_[indices,:,:,:] = h_ba_bb_fixed
        
        h_ba_bb_.resize((l2a_,l2b_))
        
        
        h_ba_c_ = numpy.resize(h_ba_c,(3,4,self.ntilt,l3_))
        h_ba_c_fixed = h_ba_c_[indices,:,:,:]
        
        h_ba_c_fixed[:,1:,:,:] -= numpy.rollaxis(h_ba_c_fixed[:,0,:,:,numpy.newaxis]*X0,3,1)
        h_ba_c_fixed[:,0,:,:] = 0
        h_ba_c_[indices,:,:,:] = h_ba_c_fixed
        h_ba_c_.resize((l2a_,l3_))
        
        
        h_ba_d_ = h_ba_d
        
        # c
        
        h_bb_a_ = h_bb_a

        h_bb_ba_ = h_ba_bb_.T
        
        h_bb_bb_ = h_bb_bb
        
        h_bb_c_ = h_bb_c
        
        h_bb_d_ = h_bb_d
        
        # c
        
        h_c_a_ = h_c_a

        h_c_ba_ = h_ba_c_.T
        
        h_c_bb_ = h_bb_c_.T
        
        h_c_c_ = h_c_c
        
        h_c_d_ = h_c_d
        
        #d
        
        h_d_a_ = h_a_d_.T
        
        h_d_ba_ = h_ba_d_.T

        h_d_bb_ = h_bb_d_.T
        
        h_d_c_ = h_c_d_.T
        
        h_d_d_ = h_d_d
       
        # Concatenate
        
        hess_ = numpy.row_stack(( \
                        numpy.column_stack((h_a_a_,h_a_ba_,h_a_bb_,h_a_c_,h_a_d_)), \
                        numpy.column_stack((h_ba_a_,h_ba_ba_,h_ba_bb_,h_ba_c_,h_ba_d_)), \
                        numpy.column_stack((h_bb_a_,h_bb_ba_,h_bb_bb_,h_bb_c_,h_bb_d_)), \
                        numpy.column_stack((h_c_a_,h_c_ba_,h_c_bb_,h_c_c_,h_c_d_)), \
                        numpy.column_stack((h_d_a_,h_d_ba_,h_d_bb_,h_d_c_,h_d_d_)), \
                                    ))

        return hess_


    def sym_error(self, type=None):

        log.info('Symbolic Calculation of the Core Cost functions, its derivatives and hessians.')

        t,u = swiginac.symbol('t'),swiginac.symbol('u')
        a = [[0]*self.nn1 for i in range(3)]
        b = [[0]*self.nn2 for i in range(2)]
        br = [[0]*4 for i in range(3)]
        c = [[0]*self.nn3 for i in range(2)]
        d = [0]*self.nn4
        mrk=[0]*2
        res=[0]*2
        XYZ, XYZ_rot, u_ = [0]*3, [0]*3, 0

        var1_,var2_,var3_,var4_=[],[],[],[]
        var2_r, var2_b = [], []

        for i in range(3):
            XYZ[i] = 0
            for j in range(self.nn1):
                a[i][j] = swiginac.symbol('a' + '_' + str(i) + '_' + str(j))
                XYZ[i] = XYZ[i] + a[i][j]*t**self.pp1[0][j]*u**self.pp1[1][j]
                var1_.append(a[i][j])

        for i in range(3):
            br[i][0] = swiginac.symbol('br' + '_' + str(i) + '_' +  str(0))
            XYZ_rot[i] = br[i][0]
            var2_r.append(br[i][0])
            for j in range(1,4):
                br[i][j] = swiginac.symbol('br' + '_' + str(i) + '_' +  str(j))
                XYZ_rot[i] += br[i][j]*XYZ[j-1]
                var2_r.append(br[i][j])

        for ii in range(2):
            for k in range(self.nn2):
                b[ii][k] = swiginac.symbol('b' + '_' + str(ii) + '_' +  str(k))
                res[ii] = res[ii] + b[ii][k]*XYZ_rot[0]**self.pp2[0][k]*XYZ_rot[1]**self.pp2[1][k]*XYZ_rot[2]**self.pp2[2][k]
                var2_b.append(b[ii][k])

        var2_ = var2_r + var2_b

        for ix in range(2):
            for i3 in range(self.nn3):
                c[ix][i3] = swiginac.symbol('c' + '_' + str(ix) + '_' +  str(i3))
                mrk[ix] = mrk[ix] + c[ix][i3]*t**self.pp3[0][i3]
                var3_.append(c[ix][i3])

        for l in range(self.nn4):
            d[l] = swiginac.symbol('d' + '_' + str(l))
            u_ = u_ + d[l]*t**self.pp4[0][l]
            var4_.append(d[l])

        log.debug('Variables (a,b,c,d) have been defined.')

        l1,l2,l3,l4 = len(var1_),len(var2_),len(var3_),len(var4_)
        l2r,l2b = len(var2_r),len(var2_b)

        log.info('Number of core variables: a=%i b=%i c=%i d=%i' %(l1,l2,l3,l4))
        log.info('Number of free variables: n=%i' %self.numberOfFreeVariables(core=True))

        res[0] = res[0].subs(u==u_)
        res[1] = res[1].subs(u==u_)

        var_ = var1_ + var2_ + var3_ + var4_ + [t]

        log.debug('Starting expansion of the x core residual...')

        res[0] = res[0].expand()
        res[1] = res[1].expand()

        P1_ = res[0]
        P2_ = res[1]

        self.Pvars = var_
        self.P1 = P1_
        self.P2 = P2_

        res[0] = res[0] - mrk[0]
        res[1] = res[1] - mrk[1]

        if type=='reprojection-only':
            # No contour alignment

            error1x = res[0]**2
            error1x = error1x.subs(t==0)
            error1x = error1x.expand()

            error1y = error1x

            substitutions = []

            for i2 in range(self.nn2): substitutions.append(b[0][i2]==b[1][i2])
            for i3 in range(self.nn3): substitutions.append(c[0][i3]==c[1][i3])

            error1y = error1y.subs(substitutions)

            self.Evars = var1_ + var2_ + var3_ + var4_
            self.E1 = error1x + error1y
            self.E2 = 0
            self.E = self.E1 + self.E2

            return

        log.debug('Residual on projection x has been calculated.')

        log.debug('res[0] = ' + str(res[0])[:70] + ' ...')

        log.debug('Starting evaluation of the x core projection error...')

        error1x = res[0]**2

        log.debug('Starting integrating the x core projection error over the feature...')

        error1x = error1x.expand()
        error1x = swiginac.integral(t,0.0,1.0,error1x).eval_integ()

        log.debug('Substitute bs and cs between expression for error in x and error in y.')

        error1y = error1x

        substitutions = []

        for i2 in range(self.nn2): substitutions.append(b[0][i2]==b[1][i2])
        for i3 in range(self.nn3): substitutions.append(c[0][i3]==c[1][i3])

        error1y = error1y.subs(substitutions)

        log.debug('Summation over both contribution x and y.')

        E1_ = error1x + error1y

        if self.enforceTangencyConstraint:

            log.debug('Starting evaluation of the core tangency error...')

            gradP1 = numpy.zeros((3),dtype=object)
            gradP2 = numpy.zeros((3),dtype=object)

            for k in range(self.nn2):

                gradP1[0] = gradP1[0] + b[0][k]*br[0][0]*self.pp2[0][k]*XYZ[0]**(self.pp2[0][k]-1)*XYZ[1]**self.pp2[1][k]*XYZ[2]**self.pp2[2][k] \
                                      + b[0][k]*br[1][0]*self.pp2[1][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**(self.pp2[1][k]-1)*XYZ[2]**self.pp2[2][k] \
                                      + b[0][k]*br[2][0]*self.pp2[2][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**self.pp2[1][k]*XYZ[2]**(self.pp2[2][k]-1)
                gradP1[0] = gradP1[0] + b[0][k]*br[0][1]*self.pp2[0][k]*XYZ[0]**(self.pp2[0][k]-1)*XYZ[1]**self.pp2[1][k]*XYZ[2]**self.pp2[2][k] \
                                      + b[0][k]*br[1][1]*self.pp2[1][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**(self.pp2[1][k]-1)*XYZ[2]**self.pp2[2][k] \
                                      + b[0][k]*br[2][1]*self.pp2[2][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**self.pp2[1][k]*XYZ[2]**(self.pp2[2][k]-1)
                gradP1[0] = gradP1[0] + b[0][k]*br[0][2]*self.pp2[0][k]*XYZ[0]**(self.pp2[0][k]-1)*XYZ[1]**self.pp2[1][k]*XYZ[2]**self.pp2[2][k] \
                                      + b[0][k]*br[1][2]*self.pp2[1][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**(self.pp2[1][k]-1)*XYZ[2]**self.pp2[2][k] \
                                      + b[0][k]*br[2][2]*self.pp2[2][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**self.pp2[1][k]*XYZ[2]**(self.pp2[2][k]-1)

                gradP1[0] = gradP1[0] + b[1][k]*br[0][0]*self.pp2[0][k]*XYZ[0]**(self.pp2[0][k]-1)*XYZ[1]**self.pp2[1][k]*XYZ[2]**self.pp2[2][k] \
                                      + b[1][k]*br[1][0]*self.pp2[1][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**(self.pp2[1][k]-1)*XYZ[2]**self.pp2[2][k] \
                                      + b[1][k]*br[2][0]*self.pp2[2][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**self.pp2[1][k]*XYZ[2]**(self.pp2[2][k]-1)
                gradP1[0] = gradP1[0] + b[1][k]*br[0][1]*self.pp2[0][k]*XYZ[0]**(self.pp2[0][k]-1)*XYZ[1]**self.pp2[1][k]*XYZ[2]**self.pp2[2][k] \
                                      + b[1][k]*br[1][1]*self.pp2[1][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**(self.pp2[1][k]-1)*XYZ[2]**self.pp2[2][k] \
                                      + b[1][k]*br[2][1]*self.pp2[2][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**self.pp2[1][k]*XYZ[2]**(self.pp2[2][k]-1)
                gradP1[0] = gradP1[0] + b[1][k]*br[0][2]*self.pp2[0][k]*XYZ[0]**(self.pp2[0][k]-1)*XYZ[1]**self.pp2[1][k]*XYZ[2]**self.pp2[2][k] \
                                      + b[1][k]*br[1][2]*self.pp2[1][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**(self.pp2[1][k]-1)*XYZ[2]**self.pp2[2][k] \
                                      + b[1][k]*br[2][2]*self.pp2[2][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**self.pp2[1][k]*XYZ[2]**(self.pp2[2][k]-1)

            Np = numpy.cross(gradP1,gradP2)

            dSdt = numpy.zeros((3),dtype=object)
            dSdu = numpy.zeros((3),dtype=object)

            for i in range(3):
                dSdt[i] = XYZ[i].diff(t)
                dSdu[i] = XYZ[i].diff(u)

            Ns = numpy.cross(dSdt,dSdu)

            for i in range(3):
                Np[i] = Np[i].expand()
                Ns[i] = Ns[i].expand()

            log.info('Np=' + str(Np)[:70] + ' ...')
            log.info('Ns=' + str(Ns)[:70] + ' ...')

            alternative = 0

            if alternative==1:
                error2_ = numpy.dot(Np,dSdt)**2 + numpy.dot(Np,dSdu)**2
                error2_ = error2_/self.L
            elif alternative==2:
                error2_ = numpy.dot(gradP1,Ns)**2 + numpy.dot(gradP2,Ns)**2
            else:
                L2 = self.L**2
                error2_ = numpy.dot(Np,Ns)
                error2_ = error2_.expand()
                error2_ = error2_*error2_/L2

            log.debug('Calculate and expand Tangency Error...')

            error2_ = error2_.expand()
            error2_ = error2_.subs(u==u_)
            error2_ = error2_.expand()

            log.debug('Starting integrating the x core projection error over the feature...')

            error2_ = swiginac.integral(t,0.0,1.0,error2_).eval_integ()
            E2_ = error2_

        else:

            E2_ = 0

        log.info('Final core error expression calculated.')

        self.Evars = var1_ + var2_ + var3_ + var4_
        self.E1 = E1_
        self.E2 = E2_
        self.E = E1_ + E2_

        E_ = E1_ + E2_

        nvars = (l1,l2r,l2b,l3,l4)    # split the contribution between single beam parameters, and sample rotation

        log.info('Number of core variables: a=%i b=%i c=%i d=%i' %(l1,l2,l3,l4))

        var_ = []
        modes = []

        if not self.fix_a:
            var_ = var_ + var1_
            modes = modes + [self.mode1]*l1

        if not self.fix_b:
            var_ = var_ + var2_
            modes = modes + [self.mode2r]*l2r
            modes = modes + [self.mode2b]*l2b

        if not self.fix_c:
            var_ = var_ + var3_
            modes = modes + [self.mode3]*l3

        if not self.fix_d:
            var_ = var_ + var4_
            modes = modes + [self.mode4]*l4

        nder = len(var_)

        log.info('error = ' + str(E_)[:70] + ' ...')

        E_ = util.asPolyM(E_,self.Evars,mode='string')

        log.info('Total number of monoms in the core cost function: %i' %numpy.size(E_.coefficients))

        log.debug('Starting to calculate Core Derivatives...')

        der_ = [0]*nder

        for i in range(nder):
            der_[i] = E_.diff(var_[i])

        log.debug('Starting to calculate Core Hessians...')

        hess_= [[]*nder for i in range(nder)]

        for i in range(nder):
            for j in range(nder):
                if j<i: hess_[i].append(hess_[j][i])
                else: hess_[i].append(der_[i].diff(var_[j]))

        log.debug('Prepare data for the C extension module')

        symbols = E_.variables

        E = (E_.coefficients,E_.powers)

        der = [0]*nder
        for i in range(nder):
            coeffs = der_[i].coefficients
            monoms = der_[i].powers
            der[i] = (coeffs,monoms)

        hess = [[0]*nder for i in range(nder)]
        for i in range(nder):
            for j in range(i,nder):
                coeffs = hess_[i][j].coefficients
                monoms = hess_[i][j].powers
                hess[i][j] = (coeffs,monoms)
                hess[j][i] = hess[i][j]

        self.core_parameters = (nvars,symbols,E,nder,modes,der,hess)

        self.var1_ = var1_
        self.var2_ = var2_
        self.var3_ = var3_
        self.var4_ = var4_

        self.E_ = E_

        self.saveCache()


    def full_sym_error(self,var1_,var2_,var3_,var4_,E_):

        log.info('Extension to Multiple tilts and patches.')

        log.info('Theoretical total number of free variables: %i.' %self.numberOfFreeVariables())

        var1,var2r,var2b,var3,var4 = [],[],[],[],[]

        for ipatch in range(self.npatch):    # for a
            for tk in var1_:
                var1.append(str(tk) + '_' + str(ipatch))

        for itilt in range(self.ntilt):    # rotations
            for tk in var2_[:12]:
                var2r.append(str(tk) + '_' + str(itilt))

        for tk in var2_[12:]:    # single beam parameter
            var2b.append(str(tk))

        for itilt in range(self.ntilt):    # for c
            for ipatch in range(self.npatch):
                for tk in var3_:
                    var3.append(str(tk) + '_' + str(itilt) + '_' + str(ipatch))

        for itilt in range(self.ntilt):    # for d
            for ipatch in range(self.npatch):
                for tk in var4_:
                    var4.append(str(tk) + '_' + str(itilt) + '_' + str(ipatch))

        l1,l2r,l2b,l3,l4 = len(var1),len(var2r),len(var2b),len(var3),len(var4)

        var2 = var2r + var2b
        l2 = len(var2)

        l = 0
        if not self.fix_a: l += l1
        if not self.fix_b: l += l2
        if not self.fix_c: l += l3
        if not self.fix_d: l += l4

        log.info('Calculated Total number of free variables (a+b+d): %i/%i. [a=%i,b=%i,c=%i,d=%i]' %(l,l1+l2+l3+l4,l1,l2,l3,l4))

        self.variables = (var1,var2,var3,var4)    # Only the general variables (see super class)

        self.nvar = (len(var1),len(var2r),len(var2b),len(var3),len(var4))
        self.varmodes = (self.mode1,self.mode2r,self.mode2b,self.mode3,self.mode4)

        self.l1 = l1
        self.l2r = l2r
        self.l2b = l2b
        self.l2 = l2
        self.l3 = l3
        self.l4 = l4


    def evaluate_error_from_vars(self, var_a=None, var_b=None, var_c=None, var_d=None):

        if var_b!=None:
            var2_ba = var_b[:12*self.ntilt]
            var2_bb = var_b[12*self.ntilt:]

        if var_a==None or self.fix_a: var_a = self.a
        if var2_ba==None or self.fix_b: var2_ba = self.b[:12*self.ntilt]
        if var2_bb==None or self.fix_b: var2_bb = self.b[12*self.ntilt:]
        if var_c==None or self.fix_c: var_c = self.c
        if var_d==None or self.fix_d: var_d = self.d

        mode_a = self.varmodes[0]
        mode_ba = self.varmodes[1]
        mode_bb = self.varmodes[2]
        mode_c = self.varmodes[3]
        mode_d = self.varmodes[4]
        
        var2_ba = self.unfold_var2_a(var2_ba)

        log.debug('Input parameter Length: [l1=%i,l2r=%i,l2b=%i,l3=%i,l4=%i]' %(len(var_a),len(var2_ba),len(var2_bb),len(var_c),len(var_d)))

        variables = ((var_a,mode_a),(var2_ba,mode_ba),(var2_bb,mode_bb),(var_c,mode_c),(var_d,mode_d))

        (E,der,hess) = residualSupp.values(self.core_parameters,self.ntilt,self.npatch,self.Emask,variables)

        log.info('Full Err = %-15.3E  Err by Tilt/Trck/Axis = %-15.4E' %(E,E/self.ntilt/self.npatch/2.0))
        
        self.var = numpy.asarray(var_a + var2_ba + var2_bb + var_c + var_d)
        self.E = E
        self.der = self.fold_derivatives( der )
        self.hess = self.fold_hessians( hess, der )

        return (E,der,hess)


    def initializeBCoefficients(self,zeros=[],skips=[]):
        """Estimation of the b coefficients from the projection map P of
        the TxBRContourAlign object. A projection at a given tilt angle can be
        forced zero if the index is included in the list zeros
        """

        nterm = self.nn2
        ntilt = self.ntilt

        b = [0]*12*ntilt + [0]*2*nterm

        for ii in range(2):
            for j in range(4):
                for itilt in range(ntilt):
                    b[ii*4*ntilt+j*ntilt+itilt] = self.P[2*itilt+ii,j]

        for itilt in range(ntilt):
            b[11*ntilt+itilt] = 1

        return b
    

    def initializeBScalingCoefficients(self):
        
        return numpy.array([1.0,0.0,0.0,0.0])


    def extractBCoefficients(self,itilt,ipatch):
        """Extract the core coefficients of the cost function for a given tilt and
        a given patch.
        """
    
        if numpy.all(self.Emask[:,ipatch]==0):
            n = 3*4 + 2*self.nn2
            return numpy.zeros(n)
        
        R = numpy.resize(self.unfold_var2_a(self.b),(3,4,self.ntilt))
        coeffs_ba = R[:,:,itilt].ravel()
        
        coeffs_bb = numpy.array(self.unfold_var2_b(self.b))

        return numpy.concatenate((coeffs_ba,coeffs_bb))
    

    def initializeBSkipping(self, projection=None, fix_tilt_range=[]):

        nobj = self.npatch
        ntilt = self.ntilt
        nn2 = self.nn2

        skip_br = numpy.zeros((12,ntilt))
        skip_bb = numpy.zeros((2,nn2))

        for itilt in range(ntilt):
            if numpy.all(self.Emask[itilt,:]==0):
                skip_br[:,itilt] = 1
                
        if len(fix_tilt_range)!=0:
            skip_br[:,fix_tilt_range] = 1

        if projection=='linear': skip_bb[:,4:] = 1
        elif projection!='non-linear': skip_bb[:,:4] = 1

        skip_b = skip_br.ravel().tolist() + skip_bb.ravel().tolist()

        return skip_b


    def relaxProjection(self, order=None):

        # For now, relax the projection only if there are only point-like markers

        if len(self.lines)!=0 or len(self.surfaces)!=0:
             return self.b

        log.info('Relax the projection map')

        # Calculate projection maps from the point-tracks position

        nn1 = self.nn1
        nn2 = self.nn2
        nn3 = self.nn3
        nobj = self.npatch
        ntilt = self.ntilt

        if order==None:
            order_min = 0
            order_max = self.n2
        else:
            order_min, order_max = order

        nn2_min = util.numberOfTerms(order_min-1,dim=3)
        nn2_max = util.numberOfTerms(order_max,dim=3)

        # Deal the 3D tracks positions

        pp2 = numpy.asarray(self.pp2)

        X = numpy.zeros((nobj,3))

        X[:,0] = self.a[:nobj]
        X[:,1] = self.a[nn1*nobj:nn1*nobj+nobj]
        X[:,2] = self.a[2*nn1*nobj:2*nn1*nobj+nobj]

        # Which one is skipped?

        skip_b = numpy.resize(self.skip_b[:12*ntilt],(3,4,ntilt))

        # Calculate the rotations

        R = numpy.resize(numpy.array(self.b[:12*ntilt]),(3,4,ntilt))
        C = numpy.resize(numpy.array(self.c),(2,nn3,ntilt,nobj))

        A = numpy.zeros((4,nobj))

        A[0,:] = 1
        A[1,:] = X[:,0]
        A[2,:] = X[:,1]
        A[3,:] = X[:,2]

        for ii in range(3):

            for itilt in range(ntilt):

                indices = numpy.where(self.Emask[itilt,:])
                nobj_ = indices[0].size
                
                if nobj_==0:
                    continue
                
                A_ = numpy.squeeze(A[:,indices])

                if ii==2:
                    C_ = A_[3,:]
                else:
                    C_ = numpy.squeeze(C[ii,0,itilt,indices])

                n_sk = numpy.sum(skip_b[ii,:,itilt]==1)
                n_nosk = numpy.sum(skip_b[ii,:,itilt]!=1)

                i_sk = numpy.where(skip_b[ii,:,itilt]==1)
                i_nosk = numpy.where(skip_b[ii,:,itilt]!=1)

                B_ = R[ii,:,itilt]

                if n_sk==0:
                    C_a = numpy.zeros_like(C_)
                else:
                    C_a = numpy.dot(B_[i_sk],A_[i_sk])
                    C_a.resize(C_.shape)

                A_nosk = numpy.resize(A_[i_nosk],(n_nosk,nobj_))

                Ainv_ = numpy.linalg.pinv(A_nosk)

                B_ = numpy.dot(C_-C_a,Ainv_)

                R[ii,i_nosk,itilt] = B_[:]


        if numpy.all(skip_b):

            R = numpy.resize(numpy.array(self.b[:12*ntilt]),(3,4,ntilt))


        # Constant Beam part

        B = numpy.resize(numpy.array(self.b[12*ntilt:],dtype=float),(2,nn2))
        C = numpy.resize(numpy.array(self.c),(2,nn3,ntilt*nobj))

        # Apply the rotation

        Xeff = numpy.zeros((nobj,3,ntilt))

        Xeff = numpy.tensordot(X,R[:,1:,:],axes=((1),(1))) + R[:,0,:]

        # Set to the different powers

        A = numpy.zeros((nn2,nobj,ntilt))

        for i in range(nn2):
            A[i,:] = Xeff[:,0,:]**pp2[0,i]*Xeff[:,1,:]**pp2[1,i]*Xeff[:,2,:]**pp2[2,i]

        A = A.swapaxes(1,2)
        A = numpy.resize(A,(nn2,ntilt*nobj))

        Ainv = numpy.linalg.pinv(A)

        #  Which one is skipped?

        skip_b = numpy.resize(self.skip_b[12*ntilt:],(2,nn2))

        skip_b[:,:nn2_min] = 1
        skip_b[:,nn2_max:] = 1

        for ii in range(2):

            indices = numpy.where(self.Emask[:,:].ravel())
            nobj_ = indices[0].size
            
            if nobj_==0:
                continue

            A_ = numpy.squeeze(A[:,indices])
            C_ = numpy.squeeze(C[ii,0,indices])

            n_sk = numpy.sum(skip_b[ii,:]==1)
            n_nosk = numpy.sum(skip_b[ii,:]!=1)

            i_sk = numpy.where(skip_b[ii,:]==1)
            i_nosk = numpy.where(skip_b[ii,:]!=1)

            B_ = B[ii,:]

            if n_sk==0:
                C_a = numpy.zeros_like(C_)
            else:
                C_a = numpy.dot(B_[i_sk],A_[i_sk])
                C_a.resize(C_.shape)

            A_nosk = numpy.resize(A_[i_nosk],(n_nosk,nobj_))

            Ainv_ = numpy.linalg.pinv(A_nosk)

            B_ = numpy.dot(C_-C_a,Ainv_)

            B[ii,i_nosk] = B_[:]

        log.info('Projection map has been relaxed!')

        return R.ravel().tolist() + B.ravel().tolist()


    def evaluateScalingParameters(self):
        '''Evaluate the scaling parameters that module the projection maps.
        At the moment, they are supposed to be related to the microscope/camera
        frame (within a constant beam model way). Calculation is performed by
        doing a regression on the linear coefficients of the projection map as well
        as the coefficients of the scaling term.
        Involving the non-linear coefficients of the system breaks the calculation...
        To be investigated more. Change of origin to the center of the camera
        '''

        # For now, relax the projection only if there are only point-like markers

        if len(self.lines)!=0 or len(self.surfaces)!=0:
             return self.b

        log.info('Evaluate the Scaling Parameters')

        # Calculate scaling parameters and non-linear projection maps coefficients from the 
        # point-tracks position

        nn1 = self.nn1
        nn2 = self.nn2
        nn3 = self.nn3
        nobj = self.npatch
        ntilt = self.ntilt

        # Deal the 3D tracks positions

        pp2 = numpy.asarray(self.pp2)

        X = numpy.zeros((nobj,3))

        X[:,0] = self.a[:nobj]
        X[:,1] = self.a[nn1*nobj:nn1*nobj+nobj]
        X[:,2] = self.a[2*nn1*nobj:2*nn1*nobj+nobj]

        # We keep the affine part of the projection map. Apply it to the track positions

        R = numpy.resize(numpy.array(self.b[:12*ntilt]),(3,4,ntilt))
        
        C = numpy.resize(numpy.array(self.c),(2,nn3,ntilt,nobj))
        
        Xeff = numpy.tensordot(X,R[:,1:,:],axes=((1),(1))) + R[:,0,:]   # of shape (nobj,3,ntilt)
        Xeff = numpy.rollaxis(Xeff,0,3)    # of shape (3,ntilt,nobj)
        
        B = numpy.concatenate(((0,0,0),self.b[12*ntilt:]))
        C = C[:,0,:,:]        # of shape (2,ntilt,nobj)
        
        # Set to the different powers

        A = numpy.zeros((3+2*nn2,2,ntilt,nobj))

        A[0,:,:,:] = -C[:,:,:]*Xeff[0,:,:]
        A[1,:,:,:] = -C[:,:,:]*Xeff[1,:,:]
        A[2,:,:,:] = -C[:,:,:]*Xeff[2,:,:]

        for i in range(nn2):
            A[3+i,0,:,:] = Xeff[0,:,:]**pp2[0,i]*Xeff[1,:,:]**pp2[1,i]*Xeff[2,:,:]**pp2[2,i]
            A[3+nn2+i,1,:,:] = A[3+i,0,:,:]


        A = numpy.resize(A,(2*nn2+3,2*ntilt*nobj))
        C = numpy.resize(C,(2*ntilt*nobj))
        
        Ainv = numpy.linalg.pinv(A)

        #  Which one is skipped?
        
        skip_b = numpy.resize(self.skip_b[12*ntilt:],(2,nn2))
        skip_b[:,4:] = 1    # skip the non-linear part of the projection map in the regression
        
        skip_b = numpy.concatenate(((0,0,0),skip_b.ravel().tolist()))
        
        mask = numpy.concatenate((self.Emask[:,:].ravel(),self.Emask[:,:].ravel()),axis=0)
        indices = numpy.where(mask)
        nobj_ = indices[0].size
        
        A_ = numpy.squeeze(A[:,indices])
        C_ = numpy.squeeze(C[indices])

        n_sk = numpy.sum(skip_b[:]==1)
        n_nosk = numpy.sum(skip_b[:]!=1)
        
        i_sk = numpy.where(skip_b[:]==1)
        i_nosk = numpy.where(skip_b[:]!=1)

        B_ = B[:]

        if n_sk==0:
            C_a = numpy.zeros_like(C_)
        else:
            C_a = numpy.dot(B_[i_sk],A_[i_sk])
            C_a.resize(C_.shape)

        A_nosk = numpy.resize(A_[i_nosk],(n_nosk,nobj_))

        Ainv_ = numpy.linalg.pinv(A_nosk)

        B_ = numpy.dot(C_-C_a,Ainv_)

        B[i_nosk] = B_[:]
            
        log.info('Projection map has been relaxed!')
        
        scaling = numpy.insert(B[:3],[0],1.0).ravel().tolist()
        b = self.b[:12*ntilt] + B[3:].ravel().tolist()
        
        return ( b, scaling )
        

    def extendProjectionMapOrder(self, n2, full_reset=True):  # to be overridden to copy self.b
        '''Extends the order of the polynomial projection map to order n2.
        Parameter full_reset allows to recalculate the whole error. It should be true
        if only a regression is used to calculate the parameters rather than doing
        the full blown minimization.'''

        # Keep track of the old parameters

        ntilt = self.ntilt

        n2_ = self.n2
        nn2_ = self.nn2

        # Set the new parameters

        residual.Residual.extendProjectionMapOrder(self, n2, full_reset)

        # Copy the old projection map

        br_ = self.b[:12*ntilt]
        bb_ = self.b[12*ntilt:]

        bb = numpy.zeros((2,self.nn2))

        for ii in range(2):
            for j in range(nn2_):
                bb[ii,j] = bb_[ii*nn2_ + j]

        self.b = br_ + bb.ravel().tolist()


    def value(self):

        if len(self.lines)==0 and len(self.surfaces)==0:

            return self.directValue()

        else:

            return residual.Residual.value(self)
        
        
    def directValue(self):
        
        # Do the direct calculation (faster) for high order projection maps 
        # because it skips polynom composition
        
        if len(self.lines)!=0 and len(self.surfaces)!=0:
            raise NotImplementedError
        
        log.info('Direct Evaluation of the residual!')
        
        ntilt = self.ntilt
        npatch = self.npatch
        nn2 = self.nn2
        nn3 = self.nn3

        pp2 = numpy.asarray(self.pp2)

        XYZ = numpy.resize(self.a,(3,self.nn1,npatch))
        XYZ = XYZ[:,0,:].swapaxes(0,1)

        R = numpy.resize(self.unfold_var2_a(self.b),(3,4,ntilt))
        Xeff = R[:,0,:] + numpy.tensordot(XYZ,R[:,1:,:],axes=((1),(1)))

        # Set to the different powers

        A = numpy.zeros((nn2,npatch,ntilt))

        for i in range(nn2):
            A[i,:] = Xeff[:,0,:]**pp2[0,i]*Xeff[:,1,:]**pp2[1,i]*Xeff[:,2,:]**pp2[2,i]

        A = A.swapaxes(1,2)

        b_c = numpy.resize(self.unfold_var2_b(self.b),(2,nn2))
        C = numpy.resize(self.c,(2,nn3,ntilt,npatch))

        res = numpy.tensordot( b_c, A, axes=((1),(0)) )
        res = res - C[:,0,:,:]
        res = res**2
        res = numpy.sum(res,axis=0)

        e1 = res*self.Emask
        e2 = 0*self.Emask

        return (e1,e2)


    def project( self, itilt, ipatch, t0=0, t1=1, ndata=20 ):
        '''Project a numpy array of points onto an image
        '''
        
        if len(self.lines)==0 and len(self.surfaces)==0:

            return self.directProject(itilt, ipatch)

        else:

            return residual.Residual.project(self,itilt,ipatch,t0=0,t1=1,ndata=20)


    def directProject( self, itilt=None, ipatch=None, t0=0, t1=1, ndata=20 ):
        
        # Do the direct calculation for high order projection maps 
        # Faster because it skips polynom composition
        
        if len(self.lines)!=0 or len(self.surfaces)!=0:
            raise NotImplementedError
        
        ntilt = self.ntilt
        npatch = self.npatch
        nn2 = self.nn2

        pp2 = numpy.asarray(self.pp2)

        XYZ = numpy.resize(self.a,(3,self.nn1,npatch))
        XYZ = numpy.squeeze(XYZ[:,0,ipatch])

        R = numpy.resize(self.b[:12*ntilt],(3,4,ntilt))

        Xeff = R[:,0,:] + numpy.tensordot(XYZ,R[:,1:],([0],[1]))

        Xeff = numpy.resize(Xeff,(npatch,3,ntilt))
        Xeff = numpy.rollaxis(Xeff,2,0)

        A = numpy.zeros((ntilt,nn2,npatch))

        for i in range(nn2):
            A[:,i,:] = Xeff[:,:,0]**pp2[0,i]*Xeff[:,:,1]**pp2[1,i]*Xeff[:,:,2]**pp2[2,i]

        b_c = numpy.resize(self.b[12*ntilt:],(2,nn2))

        proj = numpy.dot( b_c,A )
        proj = numpy.squeeze( proj[:,itilt,ipatch] )

        return ( proj[0], proj[1] )      
        

if __name__ == '__main__':

    numpy.set_printoptions(precision=4)

    ntilt = 1
    npatch = 1
    n1, n2, n3, n4 = 0, 1, 0, 0
    
    nx, ny = 500, 500
    t_constant = numpy.asarray([False, False, False])

    resid = ResidualCB(ntilt,npatch,n1,n2,n3,n4)
    
    resid.fixed_origin = numpy.array([ nx/2.0, ny/2.0, 0.0 ],dtype='float')
    resid.t_constant = numpy.asarray(t_constant)

    fix_a = False
    fix_b = False
    fix_c = False
    fix_d = False

    resid.fixVariables(fix_a=fix_a,fix_b=fix_b,fix_c=fix_c,fix_d=fix_d)

    resid.a = numpy.random.random(3*npatch*resid.nn1).tolist()
    resid.b = numpy.random.random(12*ntilt + 2*resid.nn2).tolist()
    resid.c = numpy.random.random(2*npatch*ntilt*resid.nn3).tolist()
    resid.d = numpy.random.random(npatch*ntilt*resid.nn4).tolist()

    resid.skip_a = numpy.zeros_like(resid.a).tolist()
    resid.skip_b = numpy.zeros_like(resid.b).tolist()
    resid.skip_c = numpy.zeros_like(resid.c).tolist()
    resid.skip_d = numpy.zeros_like(resid.d).tolist()

    resid.mask = numpy.zeros((ntilt,npatch))

    X = numpy.zeros((0))

    if not fix_a: X = numpy.append(X,resid.a)
    if not fix_b: X = numpy.append(X,resid.b)
    if not fix_c: X = numpy.append(X,resid.c)
    if not fix_d: X = numpy.append(X,resid.d)
    
    log.debug('len(a)=%i len(b)=%i len(c)=%i len(d)=%i' %(len(resid.a),len(resid.b),len(resid.c),len(resid.d)))

    skip = numpy.zeros_like(X)

    resid.error_struct(X,skip)
    
    with_hessians = True
    
    if with_hessians:
        d2fdX2 = resid.hess_error
    else:
        d2fdX2 = None
    
    import util.tester

    util.tester.test_derivatives(resid.error_struct,resid.der_error,X,d2fdX2=d2fdX2,args=(skip,))
    









