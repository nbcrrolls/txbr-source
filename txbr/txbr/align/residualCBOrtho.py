import numpy, numpy.random
import residual, residualCB, residualSupp
import util

from residual import log, \
                     ONLY_FEATURE_DEPENDANT, \
                     ONLY_TILT_DEPENDANT, \
                     FEATURE_AND_TILT_DEPENDANT, \
                     FEATURE_AND_TILT_INDEPENDANT

class ResidualCBOrthogonal(residualCB.ResidualCB):
    '''Orthogonal model within the single beam approximation'''
    

    def numberOfVariables(self, core=False, folded=True):
        '''Returns the total number of variables defined in that residual. 
        The orthogonal constaint is not taken into account because the symmbolic
        calculation is performed under the ResidualCB object.
        '''
        
        return residualCB.ResidualCB.numberOfVariables(self, core=core)


    def numberOfFreeVariables(self, core=False):
        '''Returns the total number of free variables defined in that residual.
        The orthogonal constaint is not taken into account because the symmbolic
        calculation is performed under the ResidualCB object.
        '''

        return residualCB.ResidualCB.numberOfFreeVariables(self, core=core)
    

    def unfold_var2(self,var2):
        '''From orthogonal model to general model'''
        
        R = self.unfold_var2_a(var2)
        R = numpy.resize(R,(3,4,self.ntilt))
        
        var2_nl = self.unfold_var2_b(var2)
        var2_nl = numpy.resize(var2_nl,(2,self.nn2))

        var2_ = numpy.zeros((2,self.nn2,self.ntilt))

        variables = ['X','Y','Z']
        powers = numpy.array(util.powerOrder(1,dim=3)).T

        for itilt in range(self.ntilt):

            rot_x = util.PolyM(variables,R[0,:,itilt],powers)
            rot_y = util.PolyM(variables,R[1,:,itilt],powers)
            rot_z = util.PolyM(variables,R[2,:,itilt],powers)
            rot_poly = {'X':rot_x,'Y':rot_y,'Z':rot_z}

            coefficients_x = var2_nl[0,:]
            coefficients_y = var2_nl[1,:]

            pp2 = numpy.array(self.pp2).T

            polyb_x = util.PolyM(variables,coefficients_x,pp2)
            polyb_x = polyb_x.compose(rot_poly)

            polyb_y = util.PolyM(variables,coefficients_y,pp2)
            polyb_y = polyb_y.compose(rot_poly)

            var2_[0,:,itilt] = polyb_x.extract_coefficients(pp2)
            var2_[1,:,itilt] = polyb_y.extract_coefficients(pp2)
            
            
        return var2_.ravel().tolist()
    
    
    def getOrthogonalParameters(self):
        
        angles = numpy.resize(self.b[3*self.ntilt:6*self.ntilt],(3,self.ntilt))
        
        return angles
    

    def unfold_var2_a( self, var2):
        '''From orthogonal model to general model'''
                
        # include the translations

        var2r = numpy.asarray(var2)

        var2_translations = var2r[:3*self.ntilt].copy()
        var2_angles = var2r[3*self.ntilt:6*self.ntilt].copy()

        var2_translations.resize((3,self.ntilt))
        var2_angles.resize((3,self.ntilt))
        
        R = numpy.zeros((3,3,self.ntilt))

        c = numpy.cos(var2_angles)
        s = numpy.sin(var2_angles)

        R[0,0,:] = c[1,:]*c[2,:]
        R[0,1,:] = c[1,:]*s[2,:]
        R[0,2,:] = s[1,:]
        R[1,0,:] = -s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:]
        R[1,1,:] = -s[0,:]*s[1,:]*s[2,:] + c[0,:]*c[2,:]
        R[1,2,:] = s[0,:]*c[1,:]
        R[2,0,:] = -c[0,:]*s[1,:]*c[2,:] + s[0,:]*s[2,:]
        R[2,1,:] = -c[0,:]*s[1,:]*s[2,:] - s[0,:]*c[2,:]
        R[2,2,:] = c[0,:]*c[1,:]
        
        # Translation part
        
        T_ = numpy.resize(var2_translations,(3,1,self.ntilt))    # if no constraint
                
        pp2 = numpy.asarray(self.pp2)
        pp2 = pp2[:,1:4]

        X0 = self.fixed_origin[0]**pp2[0,:]*self.fixed_origin[1]**pp2[1,:]*self.fixed_origin[2]**pp2[2,:]
        
        T = T_.copy()
        
        T[0,0,:] = self.fixed_origin[0]
        T[1,0,:] = self.fixed_origin[1]
        T[2,0,:] = self.fixed_origin[2]

        T[:,0,:] -= numpy.tensordot(R,X0,axes=((1),(0)))

        indices = numpy.argwhere(~self.t_constant).ravel()    # indices on no constraint
        T[indices] = T_[indices]   
        
        return numpy.concatenate((T,R),axis=1).ravel().tolist()
    
    
    def unfold_var2_b( self, var2 ):
        
        var2_b = numpy.asarray(var2)
        
        return var2_b[6*self.ntilt:]


    def fold_derivatives(self,der,angles):
        '''From general frame to orthogonal frame'''
        
        angles = numpy.asarray(angles)
        angles.resize((3,self.ntilt))

        c = numpy.cos(angles)
        s = numpy.sin(angles)
        
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
         
        J1 = numpy.zeros((3,3,3,self.ntilt))
        
        J1[0,0,0,:] = 0.0
        J1[0,0,1,:] = 0.0
        J1[0,0,2,:] = 0.0
        J1[0,1,0,:] = -c[0,:]*s[1,:]*c[2,:] + s[0,:]*s[2,:]
        J1[0,1,1,:] = -c[0,:]*s[1,:]*s[2,:] - s[0,:]*c[2,:]
        J1[0,1,2,:] = c[0,:]*c[1,:]
        J1[0,2,0,:] = s[0,:]*s[1,:]*c[2,:] + c[0,:]*s[2,:]
        J1[0,2,1,:] = s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:]
        J1[0,2,2,:] = -s[0,:]*c[1,:]

        J1[1,0,0,:] = -s[1,:]*c[2,:]
        J1[1,0,1,:] = -s[1,:]*s[2,:]
        J1[1,0,2,:] = c[1,:]
        J1[1,1,0,:] = -s[0,:]*c[1,:]*c[2,:]
        J1[1,1,1,:] = -s[0,:]*c[1,:]*s[2,:]
        J1[1,1,2,:] = -s[0,:]*s[1,:]
        J1[1,2,0,:] = -c[0,:]*c[1,:]*c[2,:]
        J1[1,2,1,:] = -c[0,:]*c[1,:]*s[2,:]
        J1[1,2,2,:] = -c[0,:]*s[1,:]

        J1[2,0,0,:] = -c[1,:]*s[2,:]
        J1[2,0,1,:] = c[1,:]*c[2,:]
        J1[2,0,2,:] = 0.0
        J1[2,1,0,:] = s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:]
        J1[2,1,1,:] = -s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:]
        J1[2,1,2,:] = 0.0
        J1[2,2,0,:] = c[0,:]*s[1,:]*s[2,:] + s[0,:]*c[2,:]
        J1[2,2,1,:] = -c[0,:]*s[1,:]*c[2,:] + s[0,:]*s[2,:]
        J1[2,2,2,:] = 0.0
        
        indices = numpy.argwhere(self.t_constant).ravel()  # indices where there is a constraint

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
        
        der_ba1_ = der_ba_[:,0,:].ravel()

        der_ba2_ = der_ba_[:,1:4,:]*J1
        der_ba2_ = der_ba2_.sum(axis=2)
        der_ba2_ = der_ba2_.sum(axis=1)    # should be of size (3,self.ntilt)
        der_ba2_ = der_ba2_.ravel()
               
        # Constant beam part of the projection map
        
        der_bb_ = der_bb
        
        # c
        
        der_c_ = der_c
        
        # d
        
        der_d_ = der_d
        
        # end

        return numpy.concatenate((der_a_, der_ba1_, der_ba2_, der_bb_, der_c_, der_d_))
        
        
    def fold_hessians(self,hess,der,angles):
        '''From general frame to orthogonal frame'''

        angles = numpy.asarray(angles)
        angles.resize((3,self.ntilt))

        c = numpy.cos(angles)
        s = numpy.sin(angles)

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
        
        J1 = numpy.zeros((3,3,3,self.ntilt))

        J1[0,0,0,:] = 0.0
        J1[0,0,1,:] = 0.0
        J1[0,0,2,:] = 0.0
        J1[0,1,0,:] = -c[0,:]*s[1,:]*c[2,:] + s[0,:]*s[2,:]
        J1[0,1,1,:] = -c[0,:]*s[1,:]*s[2,:] - s[0,:]*c[2,:]
        J1[0,1,2,:] = c[0,:]*c[1,:]
        J1[0,2,0,:] = s[0,:]*s[1,:]*c[2,:] + c[0,:]*s[2,:]
        J1[0,2,1,:] = s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:]
        J1[0,2,2,:] = -s[0,:]*c[1,:]

        J1[1,0,0,:] = -s[1,:]*c[2,:]
        J1[1,0,1,:] = -s[1,:]*s[2,:]
        J1[1,0,2,:] = c[1,:]
        J1[1,1,0,:] = -s[0,:]*c[1,:]*c[2,:]
        J1[1,1,1,:] = -s[0,:]*c[1,:]*s[2,:]
        J1[1,1,2,:] = -s[0,:]*s[1,:]
        J1[1,2,0,:] = -c[0,:]*c[1,:]*c[2,:]
        J1[1,2,1,:] = -c[0,:]*c[1,:]*s[2,:]
        J1[1,2,2,:] = -c[0,:]*s[1,:]

        J1[2,0,0,:] = -c[1,:]*s[2,:]
        J1[2,0,1,:] = c[1,:]*c[2,:]
        J1[2,0,2,:] = 0.0
        J1[2,1,0,:] = s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:]
        J1[2,1,1,:] = -s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:]
        J1[2,1,2,:] = 0.0
        J1[2,2,0,:] = c[0,:]*s[1,:]*s[2,:] + s[0,:]*c[2,:]
        J1[2,2,1,:] = -c[0,:]*s[1,:]*c[2,:] + s[0,:]*s[2,:]
        J1[2,2,2,:] = 0.0
        
        J2 = numpy.zeros((3,3,3,3,self.ntilt))

        J2[0,0,0,0,:] = 0.0
        J2[0,0,0,1,:] = 0.0
        J2[0,0,0,2,:] = 0.0
        J2[0,0,1,0,:] = s[0,:]*s[1,:]*c[2,:] + c[0,:]*s[2,:]
        J2[0,0,1,1,:] = s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:]
        J2[0,0,1,2,:] = -s[0,:]*c[1,:]
        J2[0,0,2,0,:] = c[0,:]*s[1,:]*c[2,:] - s[0,:]*s[2,:]
        J2[0,0,2,1,:] = c[0,:]*s[1,:]*s[2,:] + s[0,:]*c[2,:]
        J2[0,0,2,2,:] = -c[0,:]*c[1,:]
 
        J2[0,1,0,0,:] = 0.0
        J2[0,1,0,1,:] = 0.0
        J2[0,1,0,2,:] = 0.0
        J2[0,1,1,0,:] = -c[0,:]*c[1,:]*c[2,:]
        J2[0,1,1,1,:] = -c[0,:]*c[1,:]*s[2,:]
        J2[0,1,1,2,:] = -c[0,:]*s[1,:]
        J2[0,1,2,0,:] = s[0,:]*c[1,:]*c[2,:]
        J2[0,1,2,1,:] = s[0,:]*c[1,:]*s[2,:]
        J2[0,1,2,2,:] = s[0,:]*s[1,:]

        J2[0,2,0,0,:] = 0.0
        J2[0,2,0,1,:] = 0.0
        J2[0,2,0,2,:] = 0.0
        J2[0,2,1,0,:] = c[0,:]*s[1,:]*s[2,:] + s[0,:]*c[2,:]
        J2[0,2,1,1,:] = -c[0,:]*s[1,:]*c[2,:] + s[0,:]*s[2,:]
        J2[0,2,1,2,:] = 0.0
        J2[0,2,2,0,:] = -s[0,:]*s[1,:]*s[2,:] + c[0,:]*c[2,:]
        J2[0,2,2,1,:] = s[0,:]*s[1,:]*c[2,:] + c[0,:]*s[2,:]
        J2[0,2,2,2,:] = 0.0
        
        J2[1,0,0,0,:] = J2[0,1,0,0,:]
        J2[1,0,0,1,:] = J2[0,1,0,1,:]
        J2[1,0,0,2,:] = J2[0,1,0,2,:]
        J2[1,0,1,0,:] = J2[0,1,1,0,:]
        J2[1,0,1,1,:] = J2[0,1,1,1,:]
        J2[1,0,1,2,:] = J2[0,1,1,2,:]
        J2[1,0,2,0,:] = J2[0,1,2,0,:]
        J2[1,0,2,1,:] = J2[0,1,2,1,:]
        J2[1,0,2,2,:] = J2[0,1,2,2,:]
        
        J2[1,1,0,0,:] = -c[1,:]*c[2,:]
        J2[1,1,0,1,:] = -c[1,:]*s[2,:]
        J2[1,1,0,2,:] = -s[1,:]
        J2[1,1,1,0,:] = s[0,:]*s[1,:]*c[2,:]
        J2[1,1,1,1,:] = s[0,:]*s[1,:]*s[2,:]
        J2[1,1,1,2,:] = -s[0,:]*c[1,:]
        J2[1,1,2,0,:] = c[0,:]*s[1,:]*c[2,:]
        J2[1,1,2,1,:] = c[0,:]*s[1,:]*s[2,:]
        J2[1,1,2,2,:] = -c[0,:]*c[1,:]
        
        J2[1,2,0,0,:] = s[1,:]*s[2,:]
        J2[1,2,0,1,:] = -s[1,:]*c[2,:]
        J2[1,2,0,2,:] = 0.0
        J2[1,2,1,0,:] = s[0,:]*c[1,:]*s[2,:]
        J2[1,2,1,1,:] = -s[0,:]*c[1,:]*c[2,:]
        J2[1,2,1,2,:] = 0.0
        J2[1,2,2,0,:] = c[0,:]*c[1,:]*s[2,:]
        J2[1,2,2,1,:] = -c[0,:]*c[1,:]*c[2,:]
        J2[1,2,2,2,:] = 0
        
        J2[2,0,0,0,:] = J2[0,2,0,0,:]
        J2[2,0,0,1,:] = J2[0,2,0,1,:]
        J2[2,0,0,2,:] = J2[0,2,0,2,:]
        J2[2,0,1,0,:] = J2[0,2,1,0,:]
        J2[2,0,1,1,:] = J2[0,2,1,1,:]
        J2[2,0,1,2,:] = J2[0,2,1,2,:]
        J2[2,0,2,0,:] = J2[0,2,2,0,:]
        J2[2,0,2,1,:] = J2[0,2,2,1,:]
        J2[2,0,2,2,:] = J2[0,2,2,2,:]

        J2[2,1,0,0,:] = J2[1,2,0,0,:]
        J2[2,1,0,1,:] = J2[1,2,0,1,:]
        J2[2,1,0,2,:] = J2[1,2,0,2,:]
        J2[2,1,1,0,:] = J2[1,2,1,0,:]
        J2[2,1,1,1,:] = J2[1,2,1,1,:]
        J2[2,1,1,2,:] = J2[1,2,1,2,:]
        J2[2,1,2,0,:] = J2[1,2,2,0,:]
        J2[2,1,2,1,:] = J2[1,2,2,1,:]
        J2[2,1,2,2,:] = J2[1,2,2,2,:]

        J2[2,2,0,0,:] = -c[1,:]*c[2,:]
        J2[2,2,0,1,:] = -c[1,:]*s[2,:]
        J2[2,2,0,2,:] = 0.0
        J2[2,2,1,0,:] = s[0,:]*s[1,:]*c[2,:] + c[0,:]*s[2,:]
        J2[2,2,1,1,:] = s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:]
        J2[2,2,1,2,:] = 0.0
        J2[2,2,2,0,:] = c[0,:]*s[1,:]*c[2,:] - s[0,:]*s[2,:]
        J2[2,2,2,1,:] = c[0,:]*s[1,:]*s[2,:] + s[0,:]*c[2,:]
        J2[2,2,2,2,:] = 0.0
        
        # Include the constraints

        indices = numpy.argwhere(self.t_constant).ravel()  # indices where there is a constraint

        pp2 = numpy.asarray(self.pp2)
        pp2 = pp2[:,1:4]

        X0 = self.fixed_origin[0]**pp2[0,:]*self.fixed_origin[1]**pp2[1,:]*self.fixed_origin[2]**pp2[2,:]
        
        # Update the derivatives that will be used in ba2/ba2
        
        der_ba_ = numpy.resize(der_ba,(3,4,self.ntilt))
                        
        der_ba_fixed = der_ba_[indices]
                
        der_ba_fixed[:,1:,:] -= numpy.rollaxis(der_ba_fixed[:,0,:,numpy.newaxis]*X0,2,1)
        der_ba_fixed[:,0,:] = 0
        
        der_ba_[indices] = der_ba_fixed
        
        # a
        
        h_a_a_ = h_a_a

        h_a_ba_ = numpy.resize(h_a_ba,(l1_,3,4,self.ntilt))
        h_a_ba_fixed = h_a_ba_[:,indices,:,:]
        
        h_a_ba_fixed[:,:,1:,:] -= numpy.rollaxis(h_a_ba_fixed[:,:,0,:,numpy.newaxis]*X0,3,2)
        h_a_ba_fixed[:,:,0,:] = 0
        h_a_ba_[:,indices,:,:] = h_a_ba_fixed
        
        h_a_ba1_ = h_a_ba_[:,:,0,:]
        h_a_ba1_ = numpy.resize(h_a_ba1_,(l1_,3*self.ntilt))
        
        h_a_ba2_ = h_a_ba_[:,:,1:,:]
        h_a_ba2_ = numpy.concatenate((h_a_ba2_*J1[0],h_a_ba2_*J1[1],h_a_ba2_*J1[2]))
        h_a_ba2_ = numpy.resize(h_a_ba2_,(3,l1_,3,3,self.ntilt))
        h_a_ba2_ = numpy.sum(numpy.sum(h_a_ba2_,axis=3),axis=2)
        h_a_ba2_ = numpy.swapaxes(h_a_ba2_,0,1)
        h_a_ba2_ = numpy.resize(h_a_ba2_,(l1_,3*self.ntilt))
        
        h_a_ba_ = numpy.concatenate((h_a_ba1_,h_a_ba2_),axis=1)
        h_a_ba_ = numpy.resize(h_a_ba_,(l1_,6*self.ntilt))
        
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
        
        h_ba1_ba1_ = numpy.resize(h_ba_ba_[:,0,:,:,0,:],(3*self.ntilt,3*self.ntilt))
        
        h_ba1_ba2_ = numpy.resize(h_ba_ba_[:,0,:,:,1:,:],(3*self.ntilt,3,3,self.ntilt))
        h_ba1_ba2_ = numpy.concatenate((h_ba1_ba2_*J1[0],h_ba1_ba2_*J1[1],h_ba1_ba2_*J1[2]))
        h_ba1_ba2_ = numpy.resize(h_ba1_ba2_,(3,3*self.ntilt,3,3,self.ntilt))
        h_ba1_ba2_ = numpy.sum(numpy.sum(h_ba1_ba2_,axis=3),axis=2)
        h_ba1_ba2_ = numpy.swapaxes(h_ba1_ba2_,0,1)
        h_ba1_ba2_ = numpy.resize(h_ba1_ba2_,(3*self.ntilt,3*self.ntilt))
        
        J_ = numpy.outer(J1.ravel(),J1.ravel())
        J_ = numpy.resize(J_,(3,3,3,self.ntilt,3,3,3,self.ntilt))
        J_ = numpy.rollaxis(J_,4,1)
        
        h_ba2_ba2_ = h_ba_ba_[:,1:,:,:,1:,:]
        h_ba2_ba2_ = h_ba2_ba2_*J_
        h_ba2_ba2_ = h_ba2_ba2_.sum(axis=6)
        h_ba2_ba2_ = h_ba2_ba2_.sum(axis=5)
        h_ba2_ba2_ = h_ba2_ba2_.sum(axis=3)
        h_ba2_ba2_ = h_ba2_ba2_.sum(axis=2)
        h_ba2_ba2_ = numpy.rollaxis(h_ba2_ba2_,2,1)
                
        h_ba2_ba2_d = numpy.resize(der_ba_,(3,4,self.ntilt))
        h_ba2_ba2_d = h_ba2_ba2_d[:,1:4,:]
        h_ba2_ba2_d = h_ba2_ba2_d*J2
        h_ba2_ba2_d = h_ba2_ba2_d.sum(axis=3)
        h_ba2_ba2_d = h_ba2_ba2_d.sum(axis=2)    # should be of size (3,self.ntilt)

        for itlt in range(self.ntilt):
            h_ba2_ba2_[:,itlt,:,itlt] += h_ba2_ba2_d[:,:,itlt]

        h_ba2_ba2_ = numpy.resize(h_ba2_ba2_,(3*self.ntilt,3*self.ntilt))
        
        h_ba_ba_ = numpy.zeros((6*self.ntilt,6*self.ntilt))
        
        h_ba_ba_[:3*self.ntilt,:3*self.ntilt] = h_ba1_ba1_
        h_ba_ba_[:3*self.ntilt,3*self.ntilt:] = h_ba1_ba2_
        h_ba_ba_[3*self.ntilt:,:3*self.ntilt] = h_ba_ba_[:3*self.ntilt,3*self.ntilt:].T
        h_ba_ba_[3*self.ntilt:,3*self.ntilt:] = h_ba2_ba2_

        h_ba_bb_ = numpy.resize(h_ba_bb,(3,4,self.ntilt,l2b_))
        h_ba_bb_fixed = h_ba_bb_[indices,:,:,:]
        
        h_ba_bb_fixed[:,1:,:,:] -= numpy.rollaxis(h_ba_bb_fixed[:,0,:,:,numpy.newaxis]*X0,3,1)
        h_ba_bb_fixed[:,0,:,:] = 0
        h_ba_bb_[indices,:,:,:] = h_ba_bb_fixed
        
        #h_ba_bb_.resize((l2a_,l2b_))
        
        h_ba1_bb_ = numpy.resize(h_ba_bb_[:,0,:,:],(3*self.ntilt,l2b_))
        
        h_ba2_bb_ = h_ba_bb_[:,1:,:,:]
        h_ba2_bb_ = numpy.rollaxis(h_ba2_bb_,3,0)
        h_ba2_bb_ = numpy.concatenate((h_ba2_bb_*J1[0],h_ba2_bb_*J1[1],h_ba2_bb_*J1[2]))
        h_ba2_bb_ = numpy.resize(h_ba2_bb_,(3,2*self.nn2,3,3,self.ntilt))
        h_ba2_bb_ = numpy.sum(numpy.sum(h_ba2_bb_,axis=3),axis=2)
        h_ba2_bb_ = numpy.swapaxes(h_ba2_bb_,1,2)
        h_ba2_bb_ = numpy.resize(h_ba2_bb_,(3*self.ntilt,2*self.nn2))
        
        h_ba_bb_ = numpy.concatenate((h_ba1_bb_,h_ba2_bb_),axis=0)
        
        h_ba_c_ = numpy.resize(h_ba_c,(3,4,self.ntilt,l3_))
        h_ba_c_fixed = h_ba_c_[indices,:,:,:]
        
        h_ba_c_fixed[:,1:,:,:] -= numpy.rollaxis(h_ba_c_fixed[:,0,:,:,numpy.newaxis]*X0,3,1)
        h_ba_c_fixed[:,0,:,:] = 0
        h_ba_c_[indices,:,:,:] = h_ba_c_fixed
    
        h_ba1_c_ = numpy.resize(h_ba_c_[:,0,self.ntilt:,:],(3*self.ntilt,l3_))
        
        h_ba2_c_ = h_ba_c_[:,1:,:,:]
        h_ba2_c_ = numpy.rollaxis(h_ba2_c_,3,0)
        h_ba2_c_ = numpy.concatenate((h_ba2_c_*J1[0],h_ba2_c_*J1[1],h_ba2_c_*J1[2]))
        h_ba2_c_ = numpy.resize(h_ba2_c_,(3,l3_,3,3,self.ntilt))
        h_ba2_c_ = numpy.sum(numpy.sum(h_ba2_c_,axis=3),axis=2)
        h_ba2_c_ = numpy.swapaxes(h_ba2_c_,1,2)
        h_ba2_c_ = numpy.resize(h_ba2_c_,(3*self.ntilt,l3_))
        
        h_ba_c_ = numpy.concatenate((h_ba1_c_,h_ba2_c_),axis=0)
        
        h_ba_d_ = numpy.resize(h_ba_d[:12*self.ntilt,:],(3,4,self.ntilt,l4_))
        
        h_ba1_d_ = numpy.resize(h_ba_d_[:,0,:,:],(3*self.ntilt,l4_))
        
        h_ba2_d_ = h_ba_d_[:,1:,:,:]
        h_ba2_d_ = numpy.rollaxis(h_ba2_d_,3,0)
        h_ba2_d_ = numpy.concatenate((h_ba2_d_*J1[0],h_ba2_d_*J1[1],h_ba2_d_*J1[2]))
        h_ba2_d_ = numpy.resize(h_ba2_d_,(3,l4_,3,3,self.ntilt))
        h_ba2_d_ = numpy.sum(numpy.sum(h_ba2_d_,axis=3),axis=2)
        h_ba2_d_ = numpy.swapaxes(h_ba2_d_,1,2)
        h_ba2_d_ = numpy.resize(h_ba2_d_,(3*self.ntilt,l4_))
        
        h_ba_d_ = numpy.concatenate((h_ba1_d_,h_ba2_d_),axis=0)
        
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


    def evaluate_error_from_vars(self, var_a=None, var_b=None, var_c=None, var_d=None):

        if var_b!=None:
            var_ba = var_b[:6*self.ntilt]
            var_bb = var_b[6*self.ntilt:]

        if var_a==None or self.fix_a: var_a = self.a
        if var_ba==None or self.fix_b: var_ba = self.b[:6*self.ntilt]
        if var_bb==None or self.fix_b: var_bb = self.b[6*self.ntilt:]
        if var_c==None or self.fix_c: var_c = self.c
        if var_d==None or self.fix_d: var_d = self.d

        mode_a = self.varmodes[0]
        mode_ba = self.varmodes[1]
        mode_bb = self.varmodes[2]
        mode_c = self.varmodes[3]
        mode_d = self.varmodes[4]

        angles = var_b[3*self.ntilt:6*self.ntilt]
        var_ba = self.unfold_var2_a(var_b)
        
        log.debug('Input parameter Length: [l1=%i,l2r=%i,l2b=%i,l3=%i,l4=%i]' %(len(var_a),len(var_ba),len(var_bb),len(var_c),len(var_d)))

        variables = ((var_a,mode_a),(var_ba,mode_ba),(var_bb,mode_bb),(var_c,mode_c),(var_d,mode_d))

        (E,der,hess) = residualSupp.values(self.core_parameters,self.ntilt,self.npatch,self.Emask,variables)

        residual.log.info('Full Err = %-15.3E  Err by Tilt/Trck/Axis = %-15.4E' %(E,E/self.ntilt/self.npatch/2.0))

        self.var = numpy.asarray(var_a + var_b + var_c + var_d)
        self.E = E
        self.der = self.fold_derivatives(der,angles)
        self.hess = self.fold_hessians(hess,der,angles)

        return (E,der,hess)


    def initializeBCoefficients(self,zeros=[],skips=[]):
        """Estimation of the b coefficients from the projection map P of
        the TxBRContourAlign object. A projection at a given tilt angle can be
        forced zero if the index is included in the list zeros
        """

        nterm = self.nn2
        ntilt = self.ntilt

        b = [0]*6*ntilt + [0]*2*nterm

        for ii in range(2):
            for itilt in range(ntilt):
                b[ii*ntilt+itilt] = self.P[2*itilt+ii,0]
                
        for itilt in range(ntilt):
            b[(3+1)*ntilt+itilt] = numpy.radians(self.angles[1,itilt])

        b[6*ntilt + 1] = 1.0
        b[6*ntilt + nterm + 2] = 1.0
        
        return b
    
    
    def initializeBSkipping(self, projection=None, fix_tilt_range=[]):

        nobj = self.npatch
        ntilt = self.ntilt
        nn2 = self.nn2

        skip_br = numpy.zeros((6,ntilt))
        skip_bb = numpy.zeros((2,nn2))
        
        for itilt in range(ntilt):
            if numpy.all(self.Emask[itilt,:]==0):
                skip_br[:,itilt] = 1
        
        if len(fix_tilt_range)!=0:
            skip_br[:,fix_tilt_range] = 1
        
        skip_bb[:,:4] = 1    # Fix the Constant Beam linear part (otherwise goea to a false minimum)
                
        if projection=='linear': skip_bb[:,4:] = 1
        elif projection!='non-linear': skip_bb[:,:4] = 1
        
        skip_br[3,:] = self.phi_constant[0]
        skip_br[4,:] = self.phi_constant[1]
        skip_br[5,:] = self.phi_constant[2]

        skip_b = skip_br.ravel().tolist() + skip_bb.ravel().tolist()
        
        return skip_b
    

    def extractBCoefficients(self,itilt,ipatch):
        """Extract the core coefficients of the cost function for a given tilt and
        a given patch. Coefficients are then expressed in the frame of the 
        residualCB object.
        """
        
        return residualCB.ResidualCB.extractBCoefficients(self,itilt,ipatch)
    
    
    def relaxProjection(self, order=None):
        
        # We do not implement a relaxation of the polynomial map for an
        # orthogonal beam model
        
        return self.b
    
    
    def evaluateScalingParameters(self):
        
        # To be implemented
        log.warning("Scaling coefficients calculation not implemented for constant orthogonal beam model.")
        
        return ( self.b, self.b_scaling )


    def extendProjectionMapOrder(self, n2, full_reset=True):  # to be overridden to copy self.b
        '''Extends the order of the polynomial projection map to order n2.'''

        # Keep track of the old parameters

        n2_ = self.n2
        nn2_ = self.nn2

        ntilt = self.ntilt

        # Set the new parameters

        residual.Residual.extendProjectionMapOrder(self, n2)

        # Copy the old projection map

        b1_ = self.b[:6*ntilt]
        b2_ = self.b[6*ntilt:]

        b2_ = numpy.resize(b2_,(2,nn2_))

        b2 = numpy.zeros((2,self.nn2))
        b2[:,:nn2_] = b2_[:,:]

        self.b = b1_ + b2.ravel().tolist()


    def value(self):
        
        return residualCB.ResidualCB.value(self)
        
    
if __name__ == '__main__':

    numpy.set_printoptions(precision=4)

    ntilt = 1
    npatch = 1
    n1, n2, n3, n4 = 0, 1, 0, 0
    
    nx, ny = 500, 500
    t_constant = numpy.asarray([False, True, True])

    resid = ResidualCBOrthogonal(ntilt,npatch,n1,n2,n3,n4)
    
    resid.fixed_origin = numpy.array([ nx/2.0, ny/2.0, 0.0 ],dtype='float')
    resid.t_constant = numpy.asarray(t_constant)

    fix_a = True
    fix_b = False
    fix_c = True
    fix_d = True

    resid.fixVariables(fix_a=fix_a,fix_b=fix_b,fix_c=fix_c,fix_d=fix_d)

    resid.a = numpy.random.random(3*npatch*resid.nn1).tolist()
    resid.b = numpy.random.random(6*ntilt + 2*resid.nn2).tolist()
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
    
    with_hessians = True
    
    if with_hessians:
        d2fdX2 = resid.hess_error
    else:
        d2fdX2 = None
    
    import util.tester

    util.tester.test_derivatives(resid.error_struct,resid.der_error,X,d2fdX2=d2fdX2,args=(skip,))
