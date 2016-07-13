import numpy, numpy.random
import residual, residualMB, residualSupp

from residual import log
from residual import ONLY_FEATURE_DEPENDANT
from residual import ONLY_TILT_DEPENDANT
from residual import FEATURE_AND_TILT_DEPENDANT
from residual import FEATURE_AND_TILT_INDEPENDANT

class ResidualMBOrthogonal(residualMB.ResidualMB):
    '''Orthogonal model within the multiple beam approximation'''
    
    #SYMMETRIC_MAG = False    # Magnification is different upon the camera direction
    SYMMETRIC_MAG = True    # Magnification is the same on both camera direction

    def numberOfVariables(self, core=False):
        '''Return the number of free variables during the bundle adjustment'''

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
        n += (2*self.nn2-1)*nmode[self.mode2]    # 3 angles and two magnification scaling factors
        n += 2*self.nn3*nmode[self.mode3]
        n += self.nn4*nmode[self.mode4]

        return n


    def numberOfFreeVariables(self,core=False):
        '''Return the number of free variables during the bundle adjustment'''

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
        if not self.fix_b: n += (2*self.nn2-1)*nmode[self.mode2]    # 3 angles and two magnification scaling factors
        if not self.fix_c: n += 2*self.nn3*nmode[self.mode3]
        if not self.fix_d: n += self.nn4*nmode[self.mode4]

        return n


    def getOrthogonalParameters(self):
        
        return numpy.resize(self.b[:3*self.ntilt],(3,self.ntilt))
    
    
    def getMagnificationParameters(self):
                
        return numpy.resize(self.b[3*self.ntilt:5*self.ntilt],(2,self.ntilt))
    

    def unfold_var2( self, var2 ):
        '''Unfold the variables of an orthogonal model beam to describe a general model beam.'''
        
        var2 = numpy.asarray(var2)

        var2_angles = numpy.resize(var2[:3*self.ntilt],(3,self.ntilt))
        var2_magnifications = numpy.resize(var2[3*self.ntilt:5*self.ntilt],(2,self.ntilt))
        var2_nl = numpy.resize(var2[5*self.ntilt:],(2,self.nn2-3,self.ntilt))
        
        var2_unfold = numpy.zeros((2,self.nn2,self.ntilt))
        
        if ResidualMBOrthogonal.SYMMETRIC_MAG: 
            var2_magnifications[1,:] = var2_magnifications[0,:]

        # Linear part

        c = numpy.cos(var2_angles)
        s = numpy.sin(var2_angles)
        m = var2_magnifications

        var2_unfold[0,1,:] = m[0,:]*(c[1,:]*c[2,:])
        var2_unfold[0,2,:] = m[0,:]*(c[1,:]*s[2,:])
        var2_unfold[0,3,:] = m[0,:]*(s[1,:])
        var2_unfold[1,1,:] = m[1,:]*(-s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:])
        var2_unfold[1,2,:] = m[1,:]*(-s[0,:]*s[1,:]*s[2,:] + c[0,:]*c[2,:])
        var2_unfold[1,3,:] = m[1,:]*(s[0,:]*c[1,:])
        
        # Non-Linear part
        
        var2_unfold[:,4:,:] = var2_nl[:,1:,:]

        # Translation part
        
        var2_unfold[:,0,:] = var2_nl[:,0,:]    # if no constraint
                
        pp2 = numpy.asarray(self.pp2)

        X0 = numpy.zeros((self.nn2))
        X0 = self.fixed_origin[0]**pp2[0,:]*self.fixed_origin[1]**pp2[1,:]*self.fixed_origin[2]**pp2[2,:]
        
        var2_unfold[0,0,:] = self.fixed_origin[0]
        var2_unfold[1,0,:] = self.fixed_origin[1]
        
        var2_unfold[:,0,:] -= numpy.tensordot(var2_unfold[:,1:,:],X0[1:],axes=((1),(0)))

        if numpy.any(self.t_constant[:2]==0): 
            indices = numpy.argwhere(~self.t_constant[:2]).ravel()
            var2_unfold[indices,0,:] = var2_nl[indices,0,:]

        return var2_unfold.ravel().tolist()


    def fold_derivatives( self, der, angles, magnifications ):
        '''From general frame to orthogonal frame'''

        angles = numpy.asarray(angles)
        angles.resize((3,self.ntilt))
        m = numpy.resize(magnifications,(2,self.ntilt))
        
        if ResidualMBOrthogonal.SYMMETRIC_MAG: 
            m[1,:] = m[0,:]
            
        c = numpy.cos(angles)
        s = numpy.sin(angles)
        
        offset = [0]*4
        
        if self.fix_a: 
            l1_ = 0
        else:
            l1_ = self.l1
            
        if self.fix_b: 
            l2_ = 0
        else:
            l2_ = self.l2
            
        if self.fix_c: 
            l3_ = 0
        else:
            l3_ = self.l3
            
        if self.fix_d: 
            l4_ = 0
        else:
            l4_ = self.l4
                   
        offset[0] = l1_
        offset[1] = offset[0] + l2_
        offset[2] = offset[1] + l3_
        offset[3] = offset[2] + l4_
        
        dera = der[:offset[0]]
        derb = der[offset[0]:offset[1]]
        derc = der[offset[1]:offset[2]]
        derd = der[offset[2]:offset[3]]
        
        indices = numpy.argwhere(self.t_constant[:2]).ravel()  # index array of the additional constraints

        pp2 = numpy.asarray(self.pp2)
        X0 = self.fixed_origin[0]**pp2[0,:]*self.fixed_origin[1]**pp2[1,:]*self.fixed_origin[2]**pp2[2,:]

        J1 = numpy.zeros((5,2,3,self.ntilt))

        J1[0,0,0,:] = 0.0
        J1[0,0,1,:] = 0.0
        J1[0,0,2,:] = 0.0
        J1[0,1,0,:] = m[1,:]*(-c[0,:]*s[1,:]*c[2,:] + s[0,:]*s[2,:])
        J1[0,1,1,:] = m[1,:]*(-c[0,:]*s[1,:]*s[2,:] - s[0,:]*c[2,:])
        J1[0,1,2,:] = m[1,:]*(c[0,:]*c[1,:])

        J1[1,0,0,:] = m[0,:]*(-s[1,:]*c[2,:])
        J1[1,0,1,:] = m[0,:]*(-s[1,:]*s[2,:])
        J1[1,0,2,:] = m[0,:]*(c[1,:])
        J1[1,1,0,:] = m[1,:]*(-s[0,:]*c[1,:]*c[2,:])
        J1[1,1,1,:] = m[1,:]*(-s[0,:]*c[1,:]*s[2,:])
        J1[1,1,2,:] = m[1,:]*(-s[0,:]*s[1,:])

        J1[2,0,0,:] = m[0,:]*(-c[1,:]*s[2,:])
        J1[2,0,1,:] = m[0,:]*(c[1,:]*c[2,:])
        J1[2,0,2,:] = 0.0
        J1[2,1,0,:] = m[1,:]*(s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:])
        J1[2,1,1,:] = m[1,:]*(-s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:])
        J1[2,1,2,:] = 0.0
        
        J1[3,0,0,:] = c[1,:]*c[2,:]
        J1[3,0,1,:] = c[1,:]*s[2,:]
        J1[3,0,2,:] = s[1,:]
        J1[3,1,0,:] = 0.0
        J1[3,1,1,:] = 0.0
        J1[3,1,2,:] = 0.0
        
        J1[4,0,0,:] = 0.0
        J1[4,0,1,:] = 0.0
        J1[4,0,2,:] = 0.0
        J1[4,1,0,:] = -s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:]
        J1[4,1,1,:] = -s[0,:]*s[1,:]*s[2,:] + c[0,:]*c[2,:]
        J1[4,1,2,:] = s[0,:]*c[1,:]
        
        if ResidualMBOrthogonal.SYMMETRIC_MAG:
            
            J1[3,0,0,:] = c[1,:]*c[2,:]
            J1[3,0,1,:] = c[1,:]*s[2,:]
            J1[3,0,2,:] = s[1,:]
            J1[3,1,0,:] = -s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:]
            J1[3,1,1,:] = -s[0,:]*s[1,:]*s[2,:] + c[0,:]*c[2,:]
            J1[3,1,2,:] = s[0,:]*c[1,:]
            
            J1[4,:,:,:] = 0.0      
            
        # b

        der_proj = numpy.resize(derb,(2,self.nn2,self.ntilt))
        
        der_proj_fixed = der_proj[indices]
        der_proj_fixed[:,1:,:] -= numpy.rollaxis(der_proj_fixed[:,0,:,numpy.newaxis]*X0[1:],2,1)
        der_proj_fixed[:,0,:] = 0
        
        der_proj[indices] = der_proj_fixed
        
        der_b1 = der_proj[:,1:4,:]*J1
        der_b1 = der_b1.sum(axis=2)
        der_b1 = der_b1.sum(axis=1)    # should be of size (3,self.ntilt)
        
        der_b2 = numpy.concatenate((der_proj[:,:1,:],der_proj[:,4:,:]),axis=1)
        
        der_b1 = der_b1.ravel()
        der_b2 = der_b2.ravel()
           
        # c
        
        der_c = derc

        # concatenate
        
        der_ = numpy.concatenate((dera, der_b1, der_b2, der_c, derd))

        return der_


    def fold_hessians( self, hess, der, angles, magnifications ):
        '''From general frame to orthogonal frame'''
        
        angles = numpy.resize(angles,(3,self.ntilt))
        m = numpy.resize(magnifications,(2,self.ntilt))
        
        if ResidualMBOrthogonal.SYMMETRIC_MAG: 
            m[1,:] = m[0,:]

        c = numpy.cos(angles)
        s = numpy.sin(angles)

        offset = [0]*4
        
        if self.fix_a: 
            l1_ = 0
        else:
            l1_ = self.l1
            
        if self.fix_b: 
            l2_ = 0
        else:
            l2_ = self.l2
            
        if self.fix_c: 
            l3_ = 0
        else:
            l3_ = self.l3
            
        if self.fix_d: 
            l4_ = 0
        else:
            l4_ = self.l4
            
        offset[0] = l1_
        offset[1] = offset[0] + l2_
        offset[2] = offset[1] + l3_
        offset[3] = offset[2] + l4_
        
        dera = der[:offset[0]]
        derb = der[offset[0]:offset[1]]
        derc = der[offset[1]:offset[2]]
        derd = der[offset[2]:offset[3]]

        haa = hess[:offset[0],:offset[0]]
        hab = hess[:offset[0],offset[0]:offset[1]]
        hac = hess[:offset[0],offset[1]:offset[2]]
        had = hess[:offset[0],offset[2]:offset[3]]

        hba = hess[offset[0]:offset[1],:offset[0]]
        hbb = hess[offset[0]:offset[1],offset[0]:offset[1]]
        hbc = hess[offset[0]:offset[1],offset[1]:offset[2]]
        hbd = hess[offset[0]:offset[1],offset[2]:offset[3]]

        hca = hess[offset[1]:offset[2],:offset[0]]
        hcb = hess[offset[1]:offset[2],offset[0]:offset[1]]
        hcc = hess[offset[1]:offset[2],offset[1]:offset[2]]
        hcd = hess[offset[1]:offset[2],offset[2]:offset[3]]

        hda = hess[offset[2]:offset[3],:offset[0]]
        hdb = hess[offset[2]:offset[3],offset[0]:offset[1]]
        hdc = hess[offset[2]:offset[3],offset[1]:offset[2]]
        hdd = hess[offset[2]:offset[3],offset[2]:offset[3]]
        
        indices = numpy.argwhere(self.t_constant[:2]).ravel()  # index array of the additional constraints

        pp2 = numpy.asarray(self.pp2)
        pp2 = pp2[:,1:]

        X0 = self.fixed_origin[0]**pp2[0,:]*self.fixed_origin[1]**pp2[1,:]*self.fixed_origin[2]**pp2[2,:]

        J1 = numpy.zeros((5,2,3,self.ntilt))

        J1[0,0,0,:] = 0.0
        J1[0,0,1,:] = 0.0
        J1[0,0,2,:] = 0.0
        J1[0,1,0,:] = m[1,:]*(-c[0,:]*s[1,:]*c[2,:] + s[0,:]*s[2,:])
        J1[0,1,1,:] = m[1,:]*(-c[0,:]*s[1,:]*s[2,:] - s[0,:]*c[2,:])
        J1[0,1,2,:] = m[1,:]*(c[0,:]*c[1,:])

        J1[1,0,0,:] = m[0,:]*(-s[1,:]*c[2,:])
        J1[1,0,1,:] = m[0,:]*(-s[1,:]*s[2,:])
        J1[1,0,2,:] = m[0,:]*(c[1,:])
        J1[1,1,0,:] = m[1,:]*(-s[0,:]*c[1,:]*c[2,:])
        J1[1,1,1,:] = m[1,:]*(-s[0,:]*c[1,:]*s[2,:])
        J1[1,1,2,:] = m[1,:]*(-s[0,:]*s[1,:])

        J1[2,0,0,:] = m[0,:]*(-c[1,:]*s[2,:])
        J1[2,0,1,:] = m[0,:]*(c[1,:]*c[2,:])
        J1[2,0,2,:] = 0.0
        J1[2,1,0,:] = m[1,:]*(s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:])
        J1[2,1,1,:] = m[1,:]*(-s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:])
        J1[2,1,2,:] = 0.0
        
        J1[3,0,0,:] = c[1,:]*c[2,:]
        J1[3,0,1,:] = c[1,:]*s[2,:]
        J1[3,0,2,:] = s[1,:]
        J1[3,1,0,:] = 0.0
        J1[3,1,1,:] = 0.0
        J1[3,1,2,:] = 0.0
        
        J1[4,0,0,:] = 0.0
        J1[4,0,1,:] = 0.0
        J1[4,0,2,:] = 0.0
        J1[4,1,0,:] = -s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:]
        J1[4,1,1,:] = -s[0,:]*s[1,:]*s[2,:] + c[0,:]*c[2,:]
        J1[4,1,2,:] = s[0,:]*c[1,:]
        
        if ResidualMBOrthogonal.SYMMETRIC_MAG:
            
            J1[3,0,0,:] = c[1,:]*c[2,:]
            J1[3,0,1,:] = c[1,:]*s[2,:]
            J1[3,0,2,:] = s[1,:]
            J1[3,1,0,:] = -s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:]
            J1[3,1,1,:] = -s[0,:]*s[1,:]*s[2,:] + c[0,:]*c[2,:]
            J1[3,1,2,:] = s[0,:]*c[1,:]
            
            J1[4,:,:,:] = 0.0
                
        
        J2 = numpy.zeros((5,5,2,3,self.ntilt))

        J2[0,0,0,0,:] = 0.0
        J2[0,0,0,1,:] = 0.0
        J2[0,0,0,2,:] = 0.0
        J2[0,0,1,0,:] = m[1,:]*(s[0,:]*s[1,:]*c[2,:] + c[0,:]*s[2,:])
        J2[0,0,1,1,:] = m[1,:]*(s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:])
        J2[0,0,1,2,:] = m[1,:]*(-s[0,:]*c[1,:])

        J2[0,1,0,0,:] = 0.0
        J2[0,1,0,1,:] = 0.0
        J2[0,1,0,2,:] = 0.0
        J2[0,1,1,0,:] = m[1,:]*(-c[0,:]*c[1,:]*c[2,:])
        J2[0,1,1,1,:] = m[1,:]*(-c[0,:]*c[1,:]*s[2,:])
        J2[0,1,1,2,:] = m[1,:]*(-c[0,:]*s[1,:])

        J2[0,2,0,0,:] = 0.0
        J2[0,2,0,1,:] = 0.0
        J2[0,2,0,2,:] = 0.0
        J2[0,2,1,0,:] = m[1,:]*(c[0,:]*s[1,:]*s[2,:] + s[0,:]*c[2,:])
        J2[0,2,1,1,:] = m[1,:]*(-c[0,:]*s[1,:]*c[2,:] + s[0,:]*s[2,:])
        J2[0,2,1,2,:] = 0.0
        
        J2[0,3,0,0,:] = 0.0
        J2[0,3,0,1,:] = 0.0
        J2[0,3,0,2,:] = 0.0
        J2[0,3,1,0,:] = 0.0
        J2[0,3,1,1,:] = 0.0
        J2[0,3,1,2,:] = 0.0
        
        J2[0,4,0,0,:] = 0.0
        J2[0,4,0,1,:] = 0.0
        J2[0,4,0,2,:] = 0.0
        J2[0,4,1,0,:] = -c[0,:]*s[1,:]*c[2,:] + s[0,:]*s[2,:]
        J2[0,4,1,1,:] = -c[0,:]*s[1,:]*s[2,:] - s[0,:]*c[2,:]
        J2[0,4,1,2,:] = c[0,:]*c[1,:]
        
        J2[1,0,0,0,:] = J2[0,1,0,0,:]
        J2[1,0,0,1,:] = J2[0,1,0,1,:]
        J2[1,0,0,2,:] = J2[0,1,0,2,:]
        J2[1,0,1,0,:] = J2[0,1,1,0,:]
        J2[1,0,1,1,:] = J2[0,1,1,1,:]
        J2[1,0,1,2,:] = J2[0,1,1,2,:]

        J2[1,1,0,0,:] = m[0,:]*(-c[1,:]*c[2,:])
        J2[1,1,0,1,:] = m[0,:]*(-c[1,:]*s[2,:])
        J2[1,1,0,2,:] = m[0,:]*(-s[1,:])
        J2[1,1,1,0,:] = m[1,:]*(s[0,:]*s[1,:]*c[2,:])
        J2[1,1,1,1,:] = m[1,:]*(s[0,:]*s[1,:]*s[2,:])
        J2[1,1,1,2,:] = m[1,:]*(-s[0,:]*c[1,:])

        J2[1,2,0,0,:] = m[0,:]*(s[1,:]*s[2,:])
        J2[1,2,0,1,:] = m[0,:]*(-s[1,:]*c[2,:])
        J2[1,2,0,2,:] = 0.0
        J2[1,2,1,0,:] = m[1,:]*(s[0,:]*c[1,:]*s[2,:])
        J2[1,2,1,1,:] = m[1,:]*(-s[0,:]*c[1,:]*c[2,:])
        J2[1,2,1,2,:] = 0.0
        
        J2[1,3,0,0,:] = -s[1,:]*c[2,:]
        J2[1,3,0,1,:] = -s[1,:]*s[2,:]
        J2[1,3,0,2,:] = c[1,:]
        J2[1,3,1,0,:] = 0.0
        J2[1,3,1,1,:] = 0.0
        J2[1,3,1,2,:] = 0.0
        
        J2[1,4,0,0,:] = 0.0
        J2[1,4,0,1,:] = 0.0
        J2[1,4,0,2,:] = 0.0
        J2[1,4,1,0,:] = -s[0,:]*c[1,:]*c[2,:]
        J2[1,4,1,1,:] = -s[0,:]*c[1,:]*s[2,:]
        J2[1,4,1,2,:] = -s[0,:]*s[1,:]
        
        J2[2,0,0,0,:] = J2[0,2,0,0,:]
        J2[2,0,0,1,:] = J2[0,2,0,1,:]
        J2[2,0,0,2,:] = J2[0,2,0,2,:]
        J2[2,0,1,0,:] = J2[0,2,1,0,:]
        J2[2,0,1,1,:] = J2[0,2,1,1,:]
        J2[2,0,1,2,:] = J2[0,2,1,2,:]

        J2[2,1,0,0,:] = J2[1,2,0,0,:]
        J2[2,1,0,1,:] = J2[1,2,0,1,:]
        J2[2,1,0,2,:] = J2[1,2,0,2,:]
        J2[2,1,1,0,:] = J2[1,2,1,0,:]
        J2[2,1,1,1,:] = J2[1,2,1,1,:]
        J2[2,1,1,2,:] = J2[1,2,1,2,:]

        J2[2,2,0,0,:] = m[0,:]*(-c[1,:]*c[2,:])
        J2[2,2,0,1,:] = m[0,:]*(-c[1,:]*s[2,:])
        J2[2,2,0,2,:] = 0.0
        J2[2,2,1,0,:] = m[1,:]*(s[0,:]*s[1,:]*c[2,:] + c[0,:]*s[2,:])
        J2[2,2,1,1,:] = m[1,:]*(s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:])
        J2[2,2,1,2,:] = 0.0
        
        J2[2,3,0,0,:] = -c[1,:]*s[2,:]
        J2[2,3,0,1,:] = c[1,:]*c[2,:]
        J2[2,3,0,2,:] = 0.0
        J2[2,3,1,0,:] = 0.0
        J2[2,3,1,1,:] = 0.0
        J2[2,3,1,2,:] = 0.0
        
        J2[2,4,0,0,:] = 0.0
        J2[2,4,0,1,:] = 0.0
        J2[2,4,0,2,:] = 0.0
        J2[2,4,1,0,:] = s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:]
        J2[2,4,1,1,:] = -s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:]
        J2[2,4,1,2,:] = 0.0
        
        J2[3,0,0,0,:] = J2[0,3,0,0,:]
        J2[3,0,0,1,:] = J2[0,3,0,1,:]
        J2[3,0,0,2,:] = J2[0,3,0,2,:]
        J2[3,0,1,0,:] = J2[0,3,1,0,:]
        J2[3,0,1,1,:] = J2[0,3,1,1,:]
        J2[3,0,1,2,:] = J2[0,3,1,2,:]

        J2[3,1,0,0,:] = J2[1,3,0,0,:]
        J2[3,1,0,1,:] = J2[1,3,0,1,:]
        J2[3,1,0,2,:] = J2[1,3,0,2,:]
        J2[3,1,1,0,:] = J2[1,3,1,0,:]
        J2[3,1,1,1,:] = J2[1,3,1,1,:]
        J2[3,1,1,2,:] = J2[1,3,1,2,:]

        J2[3,2,0,0,:] = J2[2,3,0,0,:]
        J2[3,2,0,1,:] = J2[2,3,0,1,:]
        J2[3,2,0,2,:] = J2[2,3,0,2,:]
        J2[3,2,1,0,:] = J2[2,3,1,0,:]
        J2[3,2,1,1,:] = J2[2,3,1,1,:]
        J2[3,2,1,2,:] = J2[2,3,1,2,:]
        
        J2[3,3,0,0,:] = 0.0
        J2[3,3,0,1,:] = 0.0
        J2[3,3,0,2,:] = 0.0
        J2[3,3,1,0,:] = 0.0
        J2[3,3,1,1,:] = 0.0
        J2[3,3,1,2,:] = 0.0
        
        J2[3,4,0,0,:] = 0.0
        J2[3,4,0,1,:] = 0.0
        J2[3,4,0,2,:] = 0.0
        J2[3,4,1,0,:] = 0.0
        J2[3,4,1,1,:] = 0.0
        J2[3,4,1,2,:] = 0.0
        
        J2[4,0,0,0,:] = J2[0,4,0,0,:]
        J2[4,0,0,1,:] = J2[0,4,0,1,:]
        J2[4,0,0,2,:] = J2[0,4,0,2,:]
        J2[4,0,1,0,:] = J2[0,4,1,0,:]
        J2[4,0,1,1,:] = J2[0,4,1,1,:]
        J2[4,0,1,2,:] = J2[0,4,1,2,:]

        J2[4,1,0,0,:] = J2[1,4,0,0,:]
        J2[4,1,0,1,:] = J2[1,4,0,1,:]
        J2[4,1,0,2,:] = J2[1,4,0,2,:]
        J2[4,1,1,0,:] = J2[1,4,1,0,:]
        J2[4,1,1,1,:] = J2[1,4,1,1,:]
        J2[4,1,1,2,:] = J2[1,4,1,2,:]

        J2[4,2,0,0,:] = J2[2,4,0,0,:]
        J2[4,2,0,1,:] = J2[2,4,0,1,:]
        J2[4,2,0,2,:] = J2[2,4,0,2,:]
        J2[4,2,1,0,:] = J2[2,4,1,0,:]
        J2[4,2,1,1,:] = J2[2,4,1,1,:]
        J2[4,2,1,2,:] = J2[2,4,1,2,:]
        
        J2[4,3,0,0,:] = J2[3,4,0,0,:]
        J2[4,3,0,1,:] = J2[3,4,0,1,:]
        J2[4,3,0,2,:] = J2[3,4,0,2,:]
        J2[4,3,1,0,:] = J2[3,4,1,0,:]
        J2[4,3,1,1,:] = J2[3,4,1,1,:] 
        J2[4,3,1,2,:] = J2[3,4,1,2,:]
        
        J2[4,4,0,0,:] = 0.0
        J2[4,4,0,1,:] = 0.0
        J2[4,4,0,2,:] = 0.0
        J2[4,4,1,0,:] = 0.0
        J2[4,4,1,1,:] = 0.0
        J2[4,4,1,2,:] = 0.0
        
        if ResidualMBOrthogonal.SYMMETRIC_MAG:
            
            J2[0,3,0,0,:] = 0.0
            J2[0,3,0,1,:] = 0.0
            J2[0,3,0,2,:] = 0.0
            J2[0,3,1,0,:] = -c[0,:]*s[1,:]*c[2,:] + s[0,:]*s[2,:]
            J2[0,3,1,1,:] = -c[0,:]*s[1,:]*s[2,:] - s[0,:]*c[2,:]
            J2[0,3,1,2,:] = c[0,:]*c[1,:]
            
            J2[1,3,0,0,:] = -s[1,:]*c[2,:]
            J2[1,3,0,1,:] = -s[1,:]*s[2,:]
            J2[1,3,0,2,:] = c[1,:]
            J2[1,3,1,0,:] = -s[0,:]*c[1,:]*c[2,:]
            J2[1,3,1,1,:] = -s[0,:]*c[1,:]*s[2,:]
            J2[1,3,1,2,:] = -s[0,:]*s[1,:]
            
            J2[2,3,0,0,:] = -c[1,:]*s[2,:]
            J2[2,3,0,1,:] = c[1,:]*c[2,:]
            J2[2,3,0,2,:] = 0.0
            J2[2,3,1,0,:] = s[0,:]*s[1,:]*s[2,:] - c[0,:]*c[2,:]
            J2[2,3,1,1,:] = -s[0,:]*s[1,:]*c[2,:] - c[0,:]*s[2,:]
            J2[2,3,1,2,:] = 0.0
            
            J2[3,0,0,0,:] = J2[0,3,0,0,:]
            J2[3,0,0,1,:] = J2[0,3,0,1,:]
            J2[3,0,0,2,:] = J2[0,3,0,2,:]
            J2[3,0,1,0,:] = J2[0,3,1,0,:]
            J2[3,0,1,1,:] = J2[0,3,1,1,:]
            J2[3,0,1,2,:] = J2[0,3,1,2,:]

            J2[3,1,0,0,:] = J2[1,3,0,0,:]
            J2[3,1,0,1,:] = J2[1,3,0,1,:]
            J2[3,1,0,2,:] = J2[1,3,0,2,:]
            J2[3,1,1,0,:] = J2[1,3,1,0,:]
            J2[3,1,1,1,:] = J2[1,3,1,1,:]
            J2[3,1,1,2,:] = J2[1,3,1,2,:]

            J2[3,2,0,0,:] = J2[2,3,0,0,:]
            J2[3,2,0,1,:] = J2[2,3,0,1,:]
            J2[3,2,0,2,:] = J2[2,3,0,2,:]
            J2[3,2,1,0,:] = J2[2,3,1,0,:]
            J2[3,2,1,1,:] = J2[2,3,1,1,:]
            J2[3,2,1,2,:] = J2[2,3,1,2,:]
            
            J2[:,4,:,:,:] = 0.0
            J2[4,:,:,:,:] = 0.0
        
        # First row

        haa_ = haa

        hab_ = numpy.resize(hab,(l1_,2,self.nn2,self.ntilt))
        hab_fixed = hab_[:,indices,:,:]
        
        hab_fixed[:,:,1:,:] -= numpy.rollaxis(hab_fixed[:,:,0,:,numpy.newaxis]*X0,3,2)
        hab_fixed[:,:,0,:] = 0
        hab_[:,indices,:,:] = hab_fixed
        #hab_.resize(l1_,l2_)        
        
        hab1_ = hab_[:,:,1:4,:]

        hab1_ = numpy.concatenate((hab1_*J1[0],hab1_*J1[1],hab1_*J1[2],hab1_*J1[3],hab1_*J1[4]))
        hab1_ = numpy.resize(hab1_,(5,l1_,2,3,self.ntilt))
        hab1_ = numpy.sum(numpy.sum(hab1_,axis=3),axis=2)
        hab1_ = numpy.swapaxes(hab1_,0,1)
        hab1_ = numpy.resize(hab1_,(l1_,5*self.ntilt))

        hab2_ = numpy.concatenate((hab_[:,:,:1,:],hab_[:,:,4:,:]),axis=2)
        hab2_ = numpy.resize(hab2_,(l1_,2*self.nn2*self.ntilt-6*self.ntilt))

        had_ = had

        # Second row

        hb1a_ = hab1_.T
        
        der_proj = numpy.resize(derb,(2,self.nn2,self.ntilt))
        
        der_proj_fixed = der_proj[indices]
        der_proj_fixed[:,1:,:] -= numpy.rollaxis(der_proj_fixed[:,0,:,numpy.newaxis]*X0,2,1)
        der_proj_fixed[:,0,:] = 0
        
        der_proj[indices] = der_proj_fixed        
        
     
        hbb = numpy.resize(hbb,(2,self.nn2,self.ntilt,2,self.nn2,self.ntilt))
        
        hbb_ = hbb.copy()
        
        i1a = numpy.ix_(indices,[0],range(self.ntilt),range(2),range(self.nn2),range(self.ntilt))
        i1b = numpy.ix_(indices,range(1,self.nn2),range(self.ntilt),range(2),range(self.nn2),range(self.ntilt))
        u1 = hbb[i1a]
        hbb_[i1b] -= numpy.rollaxis(u1[:,0,:,:,:,:,numpy.newaxis]*X0,5,1)
        
        i2a = numpy.ix_(range(2),range(self.nn2),range(self.ntilt),indices,[0],range(self.ntilt))
        i2b = numpy.ix_(range(2),range(self.nn2),range(self.ntilt),indices,range(1,self.nn2),range(self.ntilt))
        u2 = hbb[i2a]
        hbb_[i2b] -= numpy.rollaxis(u2[:,:,:,:,0,:,numpy.newaxis]*X0,5,4)
        
        i12a = numpy.ix_(indices,[0],range(self.ntilt),indices,[0],range(self.ntilt))
        i12b = numpy.ix_(indices,range(1,self.nn2),range(self.ntilt),indices,range(1,self.nn2),range(self.ntilt))
        u12a = hbb[i12a]
        hbb_[i12b] += numpy.rollaxis(numpy.rollaxis( \
                  u12a[:,0,:,:,0,:,numpy.newaxis,numpy.newaxis]*numpy.outer(X0,X0) \
                  ,5,1),5,4)

        hbb_[i1a] = 0.0
        hbb_[i2a] = 0.0

        hbb_.resize((l2_,l2_))
        
        J_ = numpy.outer(J1.ravel(),J1.ravel())
        #print "J size: %i  ntilt: %i   size: %i " %( J_.size, self.ntilt, (3*5*2*self.ntilt)**2 )
        #J_ = numpy.resize(J_,(5,2,3,self.ntilt,5,2,3,self.ntilt))
        J_ = numpy.reshape(J_,(5,2,3,self.ntilt,5,2,3,self.ntilt))
        J_ = numpy.rollaxis(J_,4,1)

        hb1b1_ = numpy.resize(hbb_,(2,self.nn2,self.ntilt,2,self.nn2,self.ntilt))
        hb1b1_ = hb1b1_[:,1:4,:,:,1:4,:]
        hb1b1_ = hb1b1_*J_
        hb1b1_ = hb1b1_.sum(axis=6)
        hb1b1_ = hb1b1_.sum(axis=5)
        hb1b1_ = hb1b1_.sum(axis=3)
        hb1b1_ = hb1b1_.sum(axis=2)
        hb1b1_ = numpy.rollaxis(hb1b1_,2,1)

        hb1b1_d = numpy.resize(der_proj,(2,self.nn2,self.ntilt))
        hb1b1_d = hb1b1_d[:,1:4,:]
        hb1b1_d = hb1b1_d*J2
        hb1b1_d = hb1b1_d.sum(axis=3)
        hb1b1_d = hb1b1_d.sum(axis=2)    # should be of size (3,self.ntilt)

        for itlt in range(self.ntilt):
            hb1b1_[:,itlt,:,itlt] += hb1b1_d[:,:,itlt]

        hb1b1_ = numpy.resize(hb1b1_,(5*self.ntilt,5*self.ntilt))

        hb1b2_ = numpy.resize(hbb_,(2,self.nn2,self.ntilt,2,self.nn2,self.ntilt))
        hb1b2_ = numpy.concatenate((hb1b2_[:,:1,:,:,1:4,:],hb1b2_[:,4:,:,:,1:4,:]),axis=1)
        hb1b2_ = numpy.resize(hb1b2_,(2*self.nn2*self.ntilt-6*self.ntilt,2,3,self.ntilt))
        hb1b2_ = numpy.concatenate((hb1b2_*J1[0],hb1b2_*J1[1],hb1b2_*J1[2],hb1b2_*J1[3],hb1b2_*J1[4]))
        hb1b2_ = numpy.resize(hb1b2_,(5,2*self.nn2*self.ntilt-6*self.ntilt,2,3,self.ntilt))
        hb1b2_ = numpy.sum(numpy.sum(hb1b2_,axis=3),axis=2)
        hb1b2_ = numpy.swapaxes(hb1b2_,1,2)
        hb1b2_ = numpy.resize(hb1b2_,(5*self.ntilt,2*self.nn2*self.ntilt-6*self.ntilt))

        hb1d_ = numpy.zeros((5*self.ntilt,l4_))

        # Third row

        hb2a_ = hab2_.T

        hb2b1_ = hb1b2_.T

        hb2b2_ = numpy.resize(hbb_,(2,self.nn2,self.ntilt,2,self.nn2,self.ntilt))
        hb2b2_ = numpy.concatenate(( \
                                   numpy.concatenate((hb2b2_[:,:1,:,:,:1,:],hb2b2_[:,:1,:,:,4:,:]),axis=4), \
                                   numpy.concatenate((hb2b2_[:,4:,:,:,:1,:],hb2b2_[:,4:,:,:,4:,:]),axis=4) \
                                   ),axis=1)
        hb2b2_ = numpy.resize(hb2b2_,(2*self.nn2*self.ntilt-6*self.ntilt,2*self.nn2*self.ntilt-6*self.ntilt))
        
        hbc_ = numpy.resize(hbc,(2,self.nn2,self.ntilt,l3_))
        hbc_fixed = hbc_[indices,:,:,:]
        
        hbc_fixed[:,1:,:,:] -= numpy.rollaxis(hbc_fixed[:,0,:,:,numpy.newaxis]*X0,3,1)
        hbc_fixed[:,0,:,:] = 0
        hbc_[indices,:,:,:] = hbc_fixed
        hbc_.resize((l2_,l3_))

        hb2d_ = numpy.resize(hbd,(2,self.nn2,self.ntilt,l4_))
        hb2d_ = numpy.concatenate((hb2d_[:,:1,:,:],hb2d_[:,4:,:,:]),axis=1)
        hb2d_ = numpy.resize(hb2d_,(2*self.nn2*self.ntilt-6*self.ntilt,l4_))

        # Fourth row

        hda_ = hda

        hdb1_ = hb1d_.T

        hdb2_ = hb2d_.T

        hdd_ = hdd

        # Concatenate

        hess_ = numpy.row_stack(( \
                        numpy.column_stack((haa_,hab1_,hab2_,had_)), \
                        numpy.column_stack((hb1a_,hb1b1_,hb1b2_,hb1d_)), \
                        numpy.column_stack((hb2a_,hb2b1_,hb2b2_,hb2d_)), \
                        numpy.column_stack((hda_,hdb1_,hdb2_,hdd_)) \
                                    ))
        
        return hess_


    def evaluate_error_from_vars(self, var1=None, var2=None, var3=None, var4=None):

        if var1==None or self.fix_a: var1 = self.a
        if var2==None or self.fix_b: var2 = self.b
        if var3==None or self.fix_c: var3 = self.c
        if var4==None or self.fix_d: var4 = self.d

        mode1 = self.varmodes[0]
        mode2 = self.varmodes[1]
        mode3 = self.varmodes[2]
        mode4 = self.varmodes[3]
        
        if ResidualMBOrthogonal.SYMMETRIC_MAG: 
            var2[4*self.ntilt:5*self.ntilt] = var2[3*self.ntilt:4*self.ntilt]

        angles = var2[:3*self.ntilt]
        magnifications = var2[3*self.ntilt:5*self.ntilt]
        var2 = self.unfold_var2( var2 )
        
        residualMB.log.debug('Input parameter Length: [l1=%i,l2=%i,l3=%i,l4=%i]' %(len(var1),len(var2),len(var3),len(var4)))

        variables = ((var1,mode1),(var2,mode2),(var3,mode3),(var4,mode4))

        (E,der,hess) = residualSupp.values(self.core_parameters,self.ntilt,self.npatch,self.Emask,variables)

        residualMB.log.info('Full Err = %-15.3E  Err by Tilt/Trck/Axis = %-15.4E' %(E,E/self.ntilt/self.npatch/2.0))

        self.var = numpy.array(var1 + var2 + var3 + var4)
        self.E = E
        self.der = self.fold_derivatives( der, angles, magnifications )
        self.hess = self.fold_hessians( hess, der, angles, magnifications )

        return ( E, der, hess )


    def initializeBCoefficients(self,zeros=[],skips=[]):
        """Estimation of the b coefficients from the projection map P of
        the TxBRContourAlign object. A projection at a given tilt angle can be
        forced zero if the index is included in the list zeros
        """

        nn2 = self.nn2
        ntilt = self.ntilt

        b1 = numpy.zeros((5,ntilt))
        
        b1[0,:] = numpy.radians(self.angles[0])
        b1[1,:] = numpy.radians(self.angles[1])
        b1[2,:] = numpy.radians(self.angles[2])
        b1[3,:] = 1.0
        b1[4,:] = 1.0
        
        b2 = numpy.zeros((2,nn2-3,ntilt))
        
        for ii in range(2):
            for itilt in range(ntilt):
                b2[ii,0,itilt] = self.P[2*itilt+ii,0]

        return b1.ravel().tolist() + b2.ravel().tolist()


    def extractBCoefficients(self,itilt,ipatch):
        """Extract the core coefficients of the cost function for a given tilt and
        a given patch.
        """

        b = self.b

        self.b = self.unfold_var2(b)
        coeffs = residualMB.ResidualMB.extractBCoefficients(self,itilt,ipatch)

        self.b = b

        return coeffs


    def initializeBSkipping(self, projection=None, fix_tilt_range=[]):

        nobj = self.npatch
        ntilt = self.ntilt
        nn2 = self.nn2

        skip_b1 = numpy.zeros((5,ntilt))
        skip_b2 = numpy.zeros((2,nn2-3,ntilt))

        for itilt in range(ntilt):
            if numpy.all(self.Emask[itilt,:]==0):
                skip_b1[:,itilt] = 1
                skip_b2[:,:,itilt] = 1
                
        if len(fix_tilt_range)!=0:
            skip_b1[:,fix_tilt_range] = 1
            skip_b2[:,:,fix_tilt_range] = 1

        if projection=='linear': skip_b1[:,:] = 1
        elif projection!='non-linear': skip_b2[:,1:,:] = 1

        skip_b1[0,:] = self.phi_constant[0]
        skip_b1[1,:] = self.phi_constant[1]
        skip_b1[2,:] = self.phi_constant[2]
        skip_b1[3,:] = 0
        skip_b1[4,:] = 1
        
        return skip_b1.ravel().tolist() + skip_b2.ravel().tolist()


    def relaxProjection(self, order=None):
        
        # We do not implement a relaxation of the polynomial map for an
        # orthogonal beam model
        
        return residualMB.ResidualMB.relaxProjection(self)
        return self.b


    def extendProjectionMapOrder(self, n2, full_reset=True):  # to be overridden to copy self.b
        '''Extends the order of the polynomial projection map to order n2.'''

        # Keep track of the old parameters

        n2_ = self.n2
        nn2_ = self.nn2

        ntilt = self.ntilt

        # Set the new parameters

        residual.Residual.extendProjectionMapOrder(self, n2, full_reset=full_reset)

        # Copy the old projection map

        b1_ = self.b[:5*ntilt]
        b2_ = self.b[5*ntilt:]

        b2_ = numpy.resize(b2_,(2,nn2_-3,ntilt))

        b2 = numpy.zeros((2,self.nn2-3,ntilt))
        b2[:,:nn2_-3,:] = b2_[:,:,:]

        self.b = b1_ + b2.ravel().tolist()



    def value(self):

        return residualMB.ResidualMB.value( self )


if __name__ == '__main__':

    numpy.set_printoptions(precision=4)

    ntilt = 1
    npatch = 1
    n1, n2, n3, n4 = 0, 2, 0, 0
    
    nx, ny = 500, 500
    t_constant = [ False, True, False ]
    t_constant = [ False, False, False ]

    resid = ResidualMBOrthogonal(ntilt,npatch,n1,n2,n3,n4)

    resid.fixed_origin = numpy.array([ nx/2.0, ny/2.0, 0.0 ])
    resid.t_constant = numpy.asarray(t_constant)

    fix_a = False
    fix_b = False
    fix_c = True
    fix_d = True

    resid.fixVariables(fix_a=fix_a,fix_b=fix_b,fix_c=fix_c,fix_d=fix_d)

    resid.a = numpy.random.random(3*npatch*resid.nn1).tolist()
    resid.b = numpy.random.random(2*ntilt*resid.nn2-ntilt).tolist()
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

    print 'la=%i lb=%i lc=%i ld=%i' %(len(resid.a),len(resid.b),len(resid.c),len(resid.d))

    skip = numpy.zeros_like(X)

    resid.error_struct(X,skip)

    with_hessians = True
    
    if with_hessians:
        d2fdX2 = resid.hess_error
    else:
        d2fdX2 = None
    
    import util.tester
    
    util.tester.test_derivatives(resid.error_struct,resid.der_error,X,d2fdX2=d2fdX2,args=(skip,),delta_h=1.0e-6)