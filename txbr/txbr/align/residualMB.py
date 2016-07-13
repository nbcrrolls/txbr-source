import numpy
import numpy.random
import swiginac
import residual
import residualSupp
import util

from residual import log

from residual import ONLY_FEATURE_DEPENDANT
from residual import ONLY_TILT_DEPENDANT
from residual import FEATURE_AND_TILT_DEPENDANT
from residual import FEATURE_AND_TILT_INDEPENDANT


class ResidualMB(residual.Residual):
    '''Residual within the Multiple Beam Approximation'''


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
        n += 2*self.nn2*nmode[self.mode2]
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
        if not self.fix_b: n += 2*self.nn2*nmode[self.mode2]
        if not self.fix_c: n += 2*self.nn3*nmode[self.mode3]
        if not self.fix_d: n += self.nn4*nmode[self.mode4]

        return n
    

    def sym_error(self):

        log.info('Symbolic Calculation of the Core Cost functions, its derivatives and hessians.')

        t,u = swiginac.symbol('t'),swiginac.symbol('u')
        a = [[0]*self.nn1 for i in range(3)]
        b = [[0]*self.nn2 for i in range(2)]
        c = [[0]*self.nn3 for i in range(2)]
        d = [0]*self.nn4
        mrk=[0]*2
        res=[0]*2
        XYZ,u_=[0]*3,0

        var1_, var2_, var3_, var4_ = [], [], [], []

        for i in range(3):
            XYZ[i] = 0
            for j in range(self.nn1):
                a[i][j] = swiginac.symbol('a' + '_' + str(i) + '_' + str(j))
                XYZ[i] = XYZ[i] + a[i][j]*t**self.pp1[0][j]*u**self.pp1[1][j]
                var1_.append(a[i][j])

        for ii in range(2):
            for k in range(self.nn2):
                b[ii][k] = swiginac.symbol('b' + '_' + str(ii) + '_' +  str(k))
                res[ii] = res[ii] + b[ii][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**self.pp2[1][k]*XYZ[2]**self.pp2[2][k]
                var2_.append(b[ii][k])

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
                gradP1[0] = gradP1[0] + b[0][k]*self.pp2[0][k]*XYZ[0]**(self.pp2[0][k]-1)*XYZ[1]**self.pp2[1][k]*XYZ[2]**self.pp2[2][k]
                gradP1[1] = gradP1[1] + b[0][k]*self.pp2[1][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**(self.pp2[1][k]-1)*XYZ[2]**self.pp2[2][k]
                gradP1[2] = gradP1[2] + b[0][k]*self.pp2[2][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**self.pp2[1][k]*XYZ[2]**(self.pp2[2][k]-1)
                gradP2[0] = gradP2[0] + b[1][k]*self.pp2[0][k]*XYZ[0]**(self.pp2[0][k]-1)*XYZ[1]**self.pp2[1][k]*XYZ[2]**self.pp2[2][k]
                gradP2[1] = gradP2[1] + b[1][k]*self.pp2[1][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**(self.pp2[1][k]-1)*XYZ[2]**self.pp2[2][k]
                gradP2[2] = gradP2[2] + b[1][k]*self.pp2[2][k]*XYZ[0]**self.pp2[0][k]*XYZ[1]**self.pp2[1][k]*XYZ[2]**(self.pp2[2][k]-1)

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

        nvars = (l1,l2,l3,l4)

        log.info('Number of core variables: a=%i b=%i c=%i d=%i' %(l1,l2,l3,l4))

        var_ = []
        modes = []

        if not self.fix_a:
            var_ = var_ + var1_
            modes = modes + [self.mode1]*l1

        if not self.fix_b:
            var_ = var_ + var2_
            modes = modes + [self.mode2]*l2

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

        var1,var2,var3,var4 = [],[],[],[]

        for ipatch in range(self.npatch):    # for a
            for tk in var1_:
                var1.append(str(tk) + '_' + str(ipatch))

        for itilt in range(self.ntilt):    # for b
            for tk in var2_:
                var2.append(str(tk) + '_' + str(itilt))

        for itilt in range(self.ntilt):    # for c
            for ipatch in range(self.npatch):
                for tk in var3_:
                    var3.append(str(tk) + '_' + str(itilt) + '_' + str(ipatch))

        for itilt in range(self.ntilt):    # for d
            for ipatch in range(self.npatch):
                for tk in var4_:
                    var4.append(str(tk) + '_' + str(itilt) + '_' + str(ipatch))

        l1,l2,l3,l4 = len(var1),len(var2),len(var3),len(var4)

        l = 0
        if not self.fix_a: l += l1
        if not self.fix_b: l += l2
        if not self.fix_c: l += l3
        if not self.fix_d: l += l4

        log.info('Calculated Total number of free variables (a+b+d): %i/%i. [a=%i,b=%i,c=%i,d=%i]' %(l,l1+l2+l3+l4,l1,l2,l3,l4))

        self.variables = (var1,var2,var3,var4)
        self.nvar = (len(var1),len(var2),len(var3),len(var4))
        self.varmodes = (self.mode1,self.mode2,self.mode3,self.mode4)

        self.l1 = l1
        self.l2 = l2
        self.l3 = l3
        self.l4 = l4


    def unfold_var2(self, var2):
        '''In the case we want to put a constraint on the fixed points.'''

        if numpy.all(~self.t_constant[:2]): return var2

        var2 = numpy.asarray(var2)
        var2.resize((2,self.nn2,self.ntilt))

        pp2 = numpy.asarray(self.pp2)

        X0 = numpy.zeros((self.nn2))
        X0 = self.fixed_origin[0]**pp2[0,:]*self.fixed_origin[1]**pp2[1,:]*self.fixed_origin[2]**pp2[2,:]
        
        var2_unfold = var2.copy()

        var2_unfold[0,0,:] = self.fixed_origin[0]
        var2_unfold[1,0,:] = self.fixed_origin[1]

        var2_unfold[:,0,:] -= numpy.tensordot(var2[:,1:,:],X0[1:],axes=((1),(0)))

        indices = numpy.argwhere(~self.t_constant[:2]).ravel()    # indices on no constraint
        var2_unfold[indices] = var2[indices]

        return var2_unfold.ravel().tolist()


    def fold_derivatives(self, der):
        '''Fold the derivatives to take into account possible constraints on the fixed
        points.'''

        if numpy.all(~self.t_constant[:2]): return der

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
        
        der_a = der[:offset[0]]
        der_b = der[offset[0]:offset[1]]
        der_c = der[offset[1]:offset[2]]
        der_d = der[offset[2]:offset[3]]
        
        indices = numpy.argwhere(self.t_constant[:2]).ravel()  # index array of the additional constraints

        pp2 = numpy.asarray(self.pp2)
        pp2 = pp2[:,1:]

        X0 = self.fixed_origin[0]**pp2[0,:]*self.fixed_origin[1]**pp2[1,:]*self.fixed_origin[2]**pp2[2,:]
        
        # a

        der_a_ = der_a
        
        # b
        
        der_b_ = numpy.resize(der_b,(2,self.nn2,self.ntilt))
        
        ia = numpy.ix_(indices,[0],range(self.ntilt))
        ib = numpy.ix_(indices,range(1,self.nn2),range(self.ntilt))
        
        u = der_b_[ia]
        der_b_[ib] -= numpy.rollaxis(u[:,0,:,numpy.newaxis]*X0,2,1)
        der_b_[ia] = 0
        
        der_b_ = der_b_.ravel()
        
        # c
        
        der_c_ = der_c
        
        # d
        
        der_d_ = der_d
        
        # end

        return numpy.concatenate((der_a_, der_b_, der_c_, der_d_))


    def fold_hessians(self, hess, der):
        '''From general frame to orthogonal frame'''

        if numpy.all(~self.t_constant[:2]): return hess

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
        
        der_a = der[:offset[0]]
        der_b = der[offset[0]:offset[1]]
        der_c = der[offset[1]:offset[2]]
        der_d = der[offset[2]:offset[3]]

        h_a_a = hess[:offset[0],:offset[0]]
        h_a_b = hess[:offset[0],offset[0]:offset[1]]
        h_a_c = hess[:offset[0],offset[1]:offset[2]]
        h_a_d = hess[:offset[0],offset[2]:offset[3]]

        h_b_a = hess[offset[0]:offset[1],:offset[0]]
        h_b_b = hess[offset[0]:offset[1],offset[0]:offset[1]]
        h_b_c = hess[offset[0]:offset[1],offset[1]:offset[2]]
        h_b_d = hess[offset[0]:offset[1],offset[2]:offset[3]]

        h_c_a = hess[offset[1]:offset[2],:offset[0]]
        h_c_b = hess[offset[1]:offset[2],offset[0]:offset[1]]
        h_c_c = hess[offset[1]:offset[2],offset[1]:offset[2]]
        h_c_d = hess[offset[1]:offset[2],offset[2]:offset[3]]

        h_d_a = hess[offset[2]:offset[3],:offset[0]]
        h_d_b = hess[offset[2]:offset[3],offset[0]:offset[1]]
        h_d_c = hess[offset[2]:offset[3],offset[1]:offset[2]]
        h_d_d = hess[offset[2]:offset[3],offset[2]:offset[3]]
        
        # Include the constraints

        indices = numpy.argwhere(self.t_constant[:2]).ravel()  # indices where there is a constraint

        pp2 = numpy.asarray(self.pp2)
        pp2 = pp2[:,1:]

        X0 = self.fixed_origin[0]**pp2[0,:]*self.fixed_origin[1]**pp2[1,:]*self.fixed_origin[2]**pp2[2,:]

        # a
        
        h_a_a_ = h_a_a

        h_a_b_ = numpy.resize(h_a_b,(l1_,2,self.nn2,self.ntilt))
        h_a_b_fixed = h_a_b_[:,indices,:,:]
        
        h_a_b_fixed[:,:,1:,:] -= numpy.rollaxis(h_a_b_fixed[:,:,0,:,numpy.newaxis]*X0,3,2)
        h_a_b_fixed[:,:,0,:] = 0
        h_a_b_[:,indices,:,:] = h_a_b_fixed
        h_a_b_.resize(l1_,l2_)
        
        h_a_c_ = h_a_c
        
        h_a_d_ = h_a_d
        
        # b

        h_b_a_ = h_a_b_.T
        
        
        h_b_b = numpy.resize(h_b_b,(2,self.nn2,self.ntilt,2,self.nn2,self.ntilt))
        
        h_b_b_ = h_b_b.copy()
        
        i1a = numpy.ix_(indices,[0],range(self.ntilt),range(2),range(self.nn2),range(self.ntilt))
        i1b = numpy.ix_(indices,range(1,self.nn2),range(self.ntilt),range(2),range(self.nn2),range(self.ntilt))
        u1 = h_b_b[i1a]
        h_b_b_[i1b] -= numpy.rollaxis(u1[:,0,:,:,:,:,numpy.newaxis]*X0,5,1)
        
        i2a = numpy.ix_(range(2),range(self.nn2),range(self.ntilt),indices,[0],range(self.ntilt))
        i2b = numpy.ix_(range(2),range(self.nn2),range(self.ntilt),indices,range(1,self.nn2),range(self.ntilt))
        u2 = h_b_b[i2a]
        h_b_b_[i2b] -= numpy.rollaxis(u2[:,:,:,:,0,:,numpy.newaxis]*X0,5,4)
        
        i12a = numpy.ix_(indices,[0],range(self.ntilt),indices,[0],range(self.ntilt))
        i12b = numpy.ix_(indices,range(1,self.nn2),range(self.ntilt),indices,range(1,self.nn2),range(self.ntilt))
        u12a = h_b_b[i12a]
        h_b_b_[i12b] += numpy.rollaxis(numpy.rollaxis( \
                  u12a[:,0,:,:,0,:,numpy.newaxis,numpy.newaxis]*numpy.outer(X0,X0) \
                  ,5,1),5,4)

        h_b_b_[i1a] = 0.0
        h_b_b_[i2a] = 0.0

        h_b_b_.resize((l2_,l2_))
        

        h_b_c_ = numpy.resize(h_b_c,(2,self.nn2,self.ntilt,l3_))
        h_b_c_fixed = h_b_c_[indices,:,:,:]
        
        h_b_c_fixed[:,1:,:,:] -= numpy.rollaxis(h_b_c_fixed[:,0,:,:,numpy.newaxis]*X0,3,1)
        h_b_c_fixed[:,0,:,:] = 0
        h_b_c_[indices,:,:,:] = h_b_c_fixed
        h_b_c_.resize((l2_,l3_))
        
        h_b_d_ = h_b_d
        
        # c
        
        h_c_a_ = h_c_a

        h_c_b_ = h_b_c_.T
        
        h_c_c_ = h_c_c
        
        h_c_d_ = h_c_d
        
        #d
        
        h_d_a_ = h_a_d_.T

        h_d_b_ = h_b_d_.T
        
        h_d_c_ = h_c_d_.T
        
        h_d_d_ = h_d_d
       
        # Concatenate

        hess_ = numpy.row_stack(( \
                        numpy.column_stack((h_a_a_,h_a_b_,h_a_c_,h_a_d_)), \
                        numpy.column_stack((h_b_a_,h_b_b_,h_b_c_,h_b_d_)), \
                        numpy.column_stack((h_c_a_,h_c_b_,h_c_c_,h_c_d_)), \
                        numpy.column_stack((h_d_a_,h_d_b_,h_d_c_,h_d_d_)) \
                                    ))

        return hess_


    def evaluate_error_from_vars(self, var_a=None, var_b=None, var_c=None, var_d=None):

        if var_a==None or self.fix_a: var_a = self.a
        if var_b==None or self.fix_b: var_b = self.b
        if var_c==None or self.fix_c: var_c = self.c
        if var_d==None or self.fix_d: var_d = self.d

        mode_a = self.varmodes[0]
        mode_b = self.varmodes[1]
        mode_c = self.varmodes[2]
        mode_d = self.varmodes[3]

        var_b = self.unfold_var2(var_b)

        log.debug('Input parameter Length: [l1=%i,l2=%i,l3=%i,l4=%i]' %(len(var_a),len(var_b),len(var_c),len(var_d)))

        variables = ((var_a,mode_a),(var_b,mode_b),(var_c,mode_c),(var_d,mode_d))

        (E,der,hess) = residualSupp.values(self.core_parameters,self.ntilt,self.npatch,self.Emask,variables)

        log.info('Full Err = %-15.3E  Err by Tilt/Trck/Axis = %-15.4E' %(E,E/self.ntilt/self.npatch/2.0))

        self.var = numpy.asarray(var_a + var_b + var_c + var_d)
        self.E = E
        self.der = self.fold_derivatives( der )
        self.hess = self.fold_hessians( hess, der )

        return ( E, der, hess )


    def initializeBCoefficients( self, zeros=[], skips=[] ):
        """Estimation of the b coefficients from the projection map P of
        the TxBRContourAlign object. A projection at a given tilt angle can be
        forced zero if the index is included in the list zeros
        """

        nterm = self.nn2
        ntilt = self.ntilt

        b = [0]*2*nterm*ntilt

        for ii in range(2):
            for i in range(nterm):
                for itilt in range(ntilt):
                    if itilt in zeros:
                        b[ii*nterm*ntilt+i*ntilt+itilt] = 0.0
                    elif not itilt in skips:
                        b[ii*nterm*ntilt+i*ntilt+itilt] = self.P[2*itilt+ii,i]

        return b
    
    
    def initializeBScalingCoefficients(self,):
    
        scaling = numpy.zeros((4,self.ntilt))
        scaling[0,:] = 1
        
        return scaling.ravel().tolist()


    def extractBCoefficients(self,itilt,ipatch):
        """Extract the core coefficients of the cost function for a given tilt and
        a given patch.
        """

        ntilt = self.ntilt
        npatch = self.npatch

        nn2 = self.nn2

        coeffs = []

        if numpy.all(self.Emask[:,ipatch]==0):
            return [0]*2*nn2

        for ii in range(2):
            for k in range(nn2):
                coeffs.append(self.b[ii*nn2*ntilt + k*ntilt + itilt])

        return numpy.array(coeffs)


    def initializeBSkipping( self, projection=None, fix_tilt_range=[] ):

        nobj = self.npatch
        ntilt = self.ntilt

        nn2 = self.nn2

        skip_b = numpy.zeros((2,nn2,ntilt))
        
        high_order_z = numpy.array(self.pp2[2][:])>1
        #high_order_z = numpy.logical_or(high_order_z,numpy.array(self.pp2[0][:])==3)
        #high_order_z = numpy.logical_or(high_order_z,numpy.array(self.pp2[1][:])==3)
        #high_order_z = numpy.logical_or(high_order_z,numpy.array(self.pp2[0][:])==4)
        #high_order_z = numpy.logical_or(high_order_z,numpy.array(self.pp2[1][:])==4)
        
        for itilt in range(ntilt):  # Remove z terms of order larger than 1
            skip_b[:,:,itilt] = high_order_z
            
            

        for itilt in range(ntilt):
            #print "itilt=%i  %i" %(itilt,numpy.sum(self.Emask[itilt,:]))
            if numpy.all(self.Emask[itilt,:]==0):
                skip_b[:,:,itilt] = 1
              
        if len(fix_tilt_range)!=0: 
            skip_b[:,:,fix_tilt_range] = 1

        if projection=='linear': skip_b[:,4:,:] = 1
        elif projection!='non-linear': skip_b[:,:4,:] = 1
        
        skip_b = skip_b.ravel().tolist()

        return skip_b
    

    def relaxProjection( self, order=None ):
        """
        Evaluate the projection maps by inverting the correspondences between
        the 3D marker positions and their prjections.
        """
    
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

        #print "order: [%i,%i]  nn:  [%i,%i]" %(order_min,order_max,nn2_min,nn2_max)

        # Deal the 3D tracks positions

        pp2 = numpy.asarray(self.pp2)

        X = numpy.zeros((nobj,3))
        X[:,0] = self.a[:nobj]
        X[:,1] = self.a[nn1*nobj:nn1*nobj+nobj]
        X[:,2] = self.a[2*nn1*nobj:2*nn1*nobj+nobj]

        A = numpy.zeros((nn2,nobj))

        for i in range(nn2):
            A[i,:] = X[:,0]**pp2[0,i]*X[:,1]**pp2[1,i]*X[:,2]**pp2[2,i]

        Ainv = numpy.linalg.pinv(A)

        # Deal with the 2D tracks projections

        skip_b = numpy.resize(self.skip_b,(2,nn2,ntilt))

        skip_b[:,:nn2_min,:] = 1
        skip_b[:,nn2_max:,:] = 1

        if numpy.all(skip_b):
            return self.b

        B = numpy.resize(numpy.array(self.b),(2,nn2,ntilt))
        C = numpy.resize(numpy.array(self.c),(2,nn3,ntilt,nobj))

        for ii in range(2):

            for itilt in range(ntilt):

                indices = numpy.where(self.Emask[itilt,:])
                nobj_ = indices[0].size

                #print "Tilt # %i:  nobj=%i" %(itilt,nobj_)
                
                if nobj_==0:
                    log.warning("No objects on tilt #%i" %itilt)
                    continue

                A_ = numpy.squeeze(A[:,indices])
                C_ = numpy.squeeze(C[ii,0,itilt,indices])

                n_sk = numpy.sum(skip_b[ii,:,itilt]==1)
                n_nosk = numpy.sum(skip_b[ii,:,itilt]!=1)
                
                i_sk = numpy.where(skip_b[ii,:,itilt]==1)
                i_nosk = numpy.where(skip_b[ii,:,itilt]!=1)

                B_ = B[ii,:,itilt]

                if n_sk==0:
                    C_a = numpy.zeros_like(C_)
                else:
                    C_a = numpy.dot(B_[i_sk],A_[i_sk])
                    C_a.resize(C_.shape)

                A_nosk = numpy.resize(A_[i_nosk],(n_nosk,nobj_))

                # Skip the tilt if not enough information to determine adequatly the projection map
                if A_nosk.size==0 or nobj_<n_nosk: continue   
                                
                Ainv_ = numpy.linalg.pinv(A_nosk)

                B_ = numpy.dot(C_-C_a,Ainv_)

                B[ii,i_nosk,itilt] = B_[:]


        log.info('Projection map has been relaxed!')

        return B.ravel().tolist()


    def extendProjectionMapOrder(self, n2, full_reset=True):  # to be overridden to copy self.b
        '''Extends the order of the polynomial projection map to order n2.'''

        # Keep track of the old parameters

        n2_ = self.n2
        nn2_ = self.nn2

        # Set the new parameters

        residual.Residual.extendProjectionMapOrder(self, n2, full_reset=full_reset)

        # Copy the old projection map

        b = numpy.zeros((2,self.nn2,self.ntilt))

        for ii in range(2):
            for j in range(nn2_):
                for itilt in range(self.ntilt):
                    b[ii,j,itilt] = self.b[ii*nn2_*self.ntilt + j*self.ntilt + itilt]

        self.b = b.ravel().tolist()


    def value(self):

        if len(self.lines)==0 and len(self.surfaces)==0:

            return  self.directValue()

        else:

            return residual.Residual.value( self )
        
        
#    def directValue(self):
#        
#        if len(self.lines)!=0 or len(self.surfaces)!=0:
#            raise NotImplementedError
#        
#        log.info('Direct Evaluation of the residual!')
#
#        ntilt = self.ntilt
#        npatch = self.npatch
#        nn2 = self.nn2
#        nn3 = self.nn3
#
#        pp2 = numpy.asarray(self.pp2)
#
#        XYZ = numpy.resize(self.a,(3,self.nn1,npatch))
#        XYZ = XYZ[:,0,:].swapaxes(0,1)
#
#        A = numpy.zeros((nn2,npatch))
#
#        for i in range(nn2):
#            A[i,:] = XYZ[:,0]**pp2[0,i]*XYZ[:,1]**pp2[1,i]*XYZ[:,2]**pp2[2,i]
#
#        b = numpy.resize(self.b,(2,nn2,ntilt))
#        C = numpy.resize(self.c,(2,nn3,ntilt,npatch))
#
#        res = numpy.tensordot( b, A, axes=((1),(0)) )
#        res = res - C[:,0,:,:]
#        res = res**2
#        res = numpy.sum(res,axis=0)
#
#        e1 = res*self.Emask
#        e2 = 0*self.Emask
#
#        return (e1,e2)
    
    
    def project( self, itilt, ipatch, t0=0, t1=1, ndata=20 ):
        '''Project a numpy array of points onto an image
        '''
        
        if len(self.lines)==0 and len(self.surfaces)==0:

            return self.directProject(itilt=itilt, ipatch=ipatch)

        else:

            return residual.Residual.project(self,itilt,ipatch,t0=0,t1=1,ndata=20)




    
if __name__ == '__main__':

    numpy.set_printoptions(precision=4)

    ntilt = 1
    npatch = 1
    n1, n2, n3, n4 = 0, 2, 0, 0

    nx, ny = 500, 500
    t_constant = numpy.asarray([False, True, False])

    resid = ResidualMB(ntilt,npatch,n1,n2,n3,n4)
    
    resid.fixed_origin = numpy.array([ nx/2.0, ny/2.0, 0.0 ],dtype='float')
    resid.t_constant = numpy.asarray(t_constant)

    fix_a = True
    fix_b = False
    fix_c = True
    fix_d = True

    resid.fixVariables(fix_a=fix_a,fix_b=fix_b,fix_c=fix_c,fix_d=fix_d)

    resid.a = numpy.random.random(3*npatch*resid.nn1).tolist()
    resid.b = numpy.random.random(2*ntilt*resid.nn2).tolist()
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
    
    util.tester.test_derivatives(resid.error_struct,resid.der_error,X,d2fdX2=d2fdX2,args=(skip,),delta_h=1.0e-2)




