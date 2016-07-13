import os
import cPickle
import multiprocessing
import numpy
import util
import residualSupp

from txbr import log

from residualSupp import ONLY_FEATURE_DEPENDANT, \
                     ONLY_TILT_DEPENDANT, \
                     FEATURE_AND_TILT_DEPENDANT, \
                     FEATURE_AND_TILT_INDEPENDANT


def __find_obj__( iobj, mask, projection, tracks ):
    """
    Helper function to find a point object from its projection!
    """

    indices = numpy.argwhere(mask[:]==1).ravel()
    n = indices.size

    if n==0: return None

    T = numpy.resize(projection[indices,:,0],(2*n))   # translation
    P = numpy.resize(projection[indices,:,1:4],(2*n,3))   # linear projection

    mrks = tracks[indices,iobj,:,0].ravel()

    D = numpy.linalg.pinv(P)
    X = numpy.dot(D,mrks-T)

    return X


def __find_obj_star__( args ):
    """
    Wrapper to the helper function __find_obj__ used for multiprocessing purpose.
    """

    return __find_obj__(*args)


'''
A module for evaluating the residual error used in the bundle adjustment.
'''
class Residual:
    '''
    This class is useful to generate the residual error during the bundle adjustment procedure.
    '''

    global cache_directory,use_cache

    cache_directory = 'cache'
    use_cache = False

    def __init__( self, ntilt = 2, npatch = 2, n1 = 0, n2 = 1, n3 = 0, n4 = 0, \
                  fix_a = False, fix_b = True, fix_c = True, fix_d = False, \
                  prepare=True, enforceTangencyConstraint = False, \
                  structureCharacteristicLength = 100  ):

        self.ntilt = ntilt
        self.npatch = npatch
        
        self.enforceTangencyConstraint = enforceTangencyConstraint
        self.L = structureCharacteristicLength    # characteristic length

        self.n1 = n1    # for a
        self.n2 = n2    # for b
        self.n3 = n3    # for c
        self.n4 = n4    # for d

        log.info('Order: n1=%s   n2=%s   n3=%s   n4=%s',self.n1,self.n2,self.n3,self.n4)

        self.nn1 = util.numberOfTerms(self.n1,dim=2)
        self.pp1 = util.powerOrder(self.n1,dim=2)
        self.nn2 = util.numberOfTerms(self.n2,dim=3)
        self.pp2 = util.powerOrder(self.n2,dim=3)
        self.nn3 = util.numberOfTerms(self.n3,dim=1)
        self.pp3 = util.powerOrder(self.n3,dim=1)
        self.nn4 = util.numberOfTerms(self.n4,dim=1)
        self.pp4 = util.powerOrder(self.n4,dim=1)

        self.mode1 = ONLY_FEATURE_DEPENDANT
        self.mode2 = ONLY_TILT_DEPENDANT
        self.mode3 = FEATURE_AND_TILT_DEPENDANT
        self.mode4 = FEATURE_AND_TILT_DEPENDANT

        self.variables = ([],[],[],[])

        self.Emask = numpy.ones((self.ntilt,self.npatch))

        if prepare:
            self.fixVariables(fix_a,fix_b,fix_c,fix_d)


    def loadCache(self):

        key = '%i-%i-%i-%i-%s-%s-%s-%s' %(self.n1,self.n2,self.n3,self.n4,self.fix_a,self.fix_b,self.fix_c,self.fix_d)
        cache_filename = os.path.join(cache_directory,key+'.pkl')

        log.info('Load cache data from %s' %cache_filename)

        if not os.path.exists(cache_filename):
            cachedata = None
        else:
            cache_file = open(cache_filename,'rb')
            try:
                cachedata = cPickle.load(cache_file)
            except TypeError:
                cachedata = None
            cache_file.close()

        return cachedata


    def saveCache(self):

        if not use_cache: return

        if not os.path.exists(cache_directory):
            os.mkdir(cache_directory)

        key = '%i-%i-%i-%i-%s-%s-%s-%s' %(self.n1,self.n2,self.n3,self.n4,self.fix_a,self.fix_b,self.fix_c,self.fix_d)
        cache_filename = os.path.join(cache_directory,key+'.pkl')

        log.info('Cache data into %s' %cache_filename)

        cache_file = open(cache_filename,'wb')

        data = [self.Pvars,self.P1,self.P2,self.Evars,self.E1,self.E2,self.var1_,self.var2_,self.var3_,self.var4_,self.core_parameters]

        cPickle.dump(data,cache_file)

        cache_file.close()


    def dataInCache(self):

        if not use_cache: return False

        cachedata = self.loadCache()

        if cachedata!=None:
            log.info('Data in cache...')
            [self.Pvars,self.P1,self.P2,self.Evars,self.E1,self.E2,self.var1_,self.var2_,self.var3_,self.var4_,self.core_parameters] = cachedata
            log.info('Core Expressions data loaded from cache.')
            self.E_ = self.E1_ + self.E2_
            return True
        else:
            log.info('No data in Cache.')
            return False


    def numberOfVariables(self, core=False):    # to be overriden
        '''
        Return the number of variables.
        '''

        raise NotImplementedError
    
    
    def numberOfFreeVariables(self, core=False):    # to be overriden
        '''
        Return the number of free variables during the bundle adjustment.
        '''

        raise NotImplementedError


    def fixVariables(self, fix_a=False, fix_b=False, fix_c=True, fix_d=False, full_reset=False):

        try:
            if (not full_reset) and \
                fix_a==self.fix_a and fix_b==self.fix_b and fix_c==self.fix_c and fix_d==self.fix_d:

#                if not self.dataInCache():
#                    self.sym_error(type='reprojection-only')
#                self.polifyCoreFunctions()

                return

        except AttributeError:
            pass

        self.fix_a = fix_a
        self.fix_b = fix_b
        self.fix_c = fix_c
        self.fix_d = fix_d

        log.info('Fixed Variables: a=%s   b=%s   c=%s   d=%s',self.fix_a,self.fix_b,self.fix_c,self.fix_d)

        if not self.dataInCache():
            self.sym_error()

        self.polifyCoreFunctions()

        self.full_sym_error(self.var1_,self.var2_,self.var3_,self.var4_,self.E_)

        residualSupp.initialize(self.core_parameters,self.ntilt,self.npatch)

#        self.a = [0]*len(self.variables[0])
#        self.b = [0]*len(self.variables[1])
#        self.c = [0]*len(self.variables[2])
#        self.d = [0]*len(self.variables[3])
#
#        self.skip_a = [0]*len(self.variables[0])
#        self.skip_b = [0]*len(self.variables[1])
#        self.skip_c = [0]*len(self.variables[2])
#        self.skip_d = [0]*len(self.variables[3])

        self.var = numpy.zeros(0)


    def extendProjectionMapOrder(self, n2, full_reset=True):  # to be overridden to copy self.b
        '''Extends the order of the polynomial projection map to order n2.'''

        if self.n2>=n2: raise ValueError

        self.n2 = n2

        self.nn2 = util.numberOfTerms(self.n2,dim=3)
        self.pp2 = util.powerOrder(self.n2,dim=3)

        self.fixVariables(self.fix_a, self.fix_b, self.fix_c, self.fix_d, full_reset)


    def sym_error(self, type=None):  # to be overriden
        '''
        Calculate the symbolic cost function.
        '''

        self.Pvars = []
        self.P1 = 0
        self.P2 = 0

        if type=='reprojection-only':
            return

        self.Evars = var1_ + var2_ + var3_ + var4_
        self.E1 = E1_
        self.E2 = E2_
        self.E = E1_ + E2_

        self.core_parameters = None

        self.var1_ = []
        self.var2_ = []
        self.var3_ = []
        self.var4_ = []

        self.E_ = 0

        return


    def full_sym_error(self,var1_,var2_,var3_,var4_,E_):  # to be overriden
        '''
        Extend the core error to take into account all the tilts and tracks.
        '''

        var1, var2, var3, var4 = [], [], [], []
        l1, l2, l3, l4 = len(var1), len(var2), len(var3), len(var4)

        self.variables = (var1,var2,var3,var4)
        self.nvar = (len(var1),len(var2),len(var3),len(var4))
        self.varmodes = (self.mode1,self.mode2,self.mode3,self.mode4)

        self.l1 = l1
        self.l2 = l2
        self.l3 = l3
        self.l4 = l4


    def evaluate_error_from_vars( self, var1=None, var2=None, var3=None, var4=None ):  # to be overriden

        raise NoImplementedError


    def polifyCoreFunctions( self ):

        log.info('Polify the Error functions...')

        log.debug('Polify the Reprojection Error...')

        self.P1 = util.asPolyM(self.P1,self.Pvars,mode='string')
        self.P2 = util.asPolyM(self.P2,self.Pvars,mode='string')

        log.debug('Polify the Tangency Error...')

        self.E1 = util.asPolyM(self.E1,self.Evars,mode='string')
        self.E2 = util.asPolyM(self.E2,self.Evars,mode='string')


    def value(self):  # Eventually overriden

        ntilt = self.ntilt
        npatch = self.npatch

        nn = self.numberOfVariables(core=True)
        args = numpy.zeros((ntilt*npatch,nn))

        if nn!=len(self.Evars):
            log.warning('Residual object has not been reset! Cannot evaluate error accurately...')
            return

        for itilt in range(ntilt):
            for ipatch in range(npatch):
                coeffs = self.extract_coefficients(itilt,ipatch)
                args[itilt*npatch+ipatch,:] = coeffs[:]

        e1 = self.E1.eval(args)
        e2 = self.E2.eval(args)

        e1.resize((ntilt,npatch))
        e2.resize((ntilt,npatch))

        e1 = e1*self.Emask
        e2 = e2*self.Emask

        return ( e1, e2 )
    
    
    def directValue( self ):
        '''Direct residual of the error value for the point-like structures. 
        The direct calculation is done in the general frame and involve switching
        from models (from orthogonal beam to general beam for instance)'''
            
        if len(self.lines)!=0 and len(self.surfaces)!=0:
            raise NotImplementedError
        
        log.info('Direct evaluation of the reprojection error from the point markers!')

        ntilt = self.ntilt
        npatch = self.npatch
        nn2 = self.nn2
        nn3 = self.nn3

        pp2 = numpy.asarray(self.pp2)

        XYZ = numpy.resize(self.a,(3,self.nn1,npatch))
        XYZ = XYZ[:,0,:].swapaxes(0,1)

        A = numpy.zeros((nn2,npatch))

        for i in range(nn2):
            A[i,:] = XYZ[:,0]**pp2[0,i]*XYZ[:,1]**pp2[1,i]*XYZ[:,2]**pp2[2,i]
            
        b = numpy.resize(self.unfold_var2(self.b),(2,nn2,ntilt))
        C = numpy.resize(self.unfold_var3(self.c),(2,nn3,ntilt,npatch))
        
        res = numpy.tensordot( b, A, axes=((1),(0)) )
        res = res - C[:,0,:,:]
        
        res = res**2
        res = numpy.sum(res,axis=0)
        
        e1 = numpy.resize(res,(ntilt,npatch))
        
        e1 = e1*self.Emask
        e2 = 0*self.Emask
        
        return (e1,e2)   


    def project( self, itilt, ipatch, t0=0, t1=1, ndata=20 ):
        '''Project a numpy array of points onto an image
        '''

        coeffs = self.extract_coefficients(itilt,ipatch)

        nn = self.numberOfVariables(core=True)
        
        if nn!=len(self.Evars):  # This can be problematic with the orthogonal approximation
            log.warning('Residual object has not been reset! Cannot evaluate projections accurately...')
            return None

        t0, t1 = 0, 1

        if ipatch in self.points:
            ndata = 1

        data = numpy.linspace(t0,t1,num=ndata)

        args_ = numpy.column_stack((numpy.tile(coeffs,(ndata,1)),data))

        xdata = self.P1.eval(args_)
        ydata = self.P2.eval(args_)

        return (xdata, ydata)
    
    
    def directProject( self, itilt=None, ipatch=None, t0=0, t1=1, ndata=20 ):
        '''Calculate the projection of patch ipatch on tilt itilt'''

        if len(self.lines)!=0 or len(self.surfaces)!=0:
            raise NotImplementedError

        ntilt = self.ntilt
        npatch = self.npatch
        nn2 = self.nn2
        nn3 = self.nn3

        pp2 = numpy.asarray(self.pp2)

        XYZ = numpy.resize(self.a,(3,self.nn1,npatch))
        XYZ = XYZ[:,0,:]

        A = numpy.zeros((nn2,npatch))

        for i in range(nn2):
            A[i,:] = XYZ[0,:]**pp2[0,i]*XYZ[1,:]**pp2[1,i]*XYZ[2,:]**pp2[2,i]

        b = numpy.resize(self.unfold_var2(self.b),(2,nn2,ntilt))

        # Projection array of size (2,ntilt,npatch)
        # if itilt or ipatch is None, it will create a new axis in the array that will
        # be reduced by the final squeeze
        proj = numpy.tensordot( b[:,:,itilt], A[:,ipatch], axes=((1),(0)))

        return numpy.squeeze(proj)


    def evaluate_error( self, var, skip ):

        self.var = numpy.asarray(var)

        l1 = len(self.a)
        l2 = len(self.b)
        l3 = len(self.c)
        l4 = len(self.d)
        
        var = var.tolist()
        skip = skip.tolist()

        if not self.fix_a:
            var1 = var[:l1]
            skip1 = skip[:l1]
            del var[:l1]
            del skip[:l1]
        else:
             var1 = self.a
             skip1 = self.skip_a

        if not self.fix_b:
            var2 = var[:l2]
            skip2 = skip[:l2]
            del var[:l2]
            del skip[:l2]
        else:
             var2 = self.b
             skip2 = self.skip_b

        if not self.fix_c:
            var3 = var[:l3]
            skip3 = skip[:l3]
            del var[:l3]
            del skip[:l3]
        else:
             var3 = self.c
             skip3 = self.skip_c

        if not self.fix_d:
            var4 = var[:l4]
            skip4 = skip[:l4]
            del var[:l4]
            del skip[:l4]
        else:
             var4 = self.d
             skip4 = self.skip_d
             
        var1 = numpy.asarray(var1)
        var2 = numpy.asarray(var2)
        var3 = numpy.asarray(var3)
        var4 = numpy.asarray(var4)
        
        skip1 = numpy.asarray(skip1)
        skip2 = numpy.asarray(skip2)
        skip3 = numpy.asarray(skip3)
        skip4 = numpy.asarray(skip4)
        
        var1 = numpy.where(skip1,numpy.asarray(self.a),var1)
        var2 = numpy.where(skip2,numpy.asarray(self.b),var2)
        var3 = numpy.where(skip3,numpy.asarray(self.c),var3)
        var4 = numpy.where(skip4,numpy.asarray(self.d),var4)
        
        var1 = var1.tolist()
        var2 = var2.tolist()
        var3 = var3.tolist()
        var4 = var4.tolist()
             
        (E,der,hess) = self.evaluate_error_from_vars( var1, var2, var3, var4 )


    def error_struct(self,var,skip):

        self.evaluate_error(var,skip)

        return self.E


    def der_error(self,var,skip):
        
        if var.size!=self.var.size or numpy.any(var-self.var!=0):
            self.evaluate_error(var,skip)
            
        ok = numpy.where(skip==False,True,False)
        
        return self.der*ok


    def hess_error(self,var,skip):

        if var.size!=self.var.size or numpy.any(var-self.var!=0):
            self.evaluate_error(var,skip)

        ok1 = numpy.where(skip==False,True,False)
        ok2 = numpy.outer(ok1,ok1)

        return self.hess*ok2


    def unfold_var1(self, var1):
        '''Unfold the a coefficients'''

        return var1


    def unfold_var2(self, var2):
        '''Unfold the projection coefficients'''

        return var2
    
    
    def unfold_var2_with_scaling(self, var2, scaling):
        '''Unfold the projection coefficients (with the scaling parameters)'''

        return ( self.unfold_var2(var2), scaling )


    def unfold_var3(self, var3):
        '''Unfold the a coefficients'''

        return var3


    def unfold_var4(self, var4):
        '''Unfold the a coefficients'''

        return var4


    def initializeACoefficients( self, fromSingleSeries=False ):
        """
        Rough estimation of the a coefficients from the oder 0 term
        of the marker traces. Valid for ponctual markers like gold fiducials.
        For more elaborated contours, estimation is done in the initializeStructureSet()
        function
        """
                
        INITIALIZATION_FROM_SINGLE_SERIE = fromSingleSeries

        nterm = self.nn1
        nobj = self.npatch

        tilt_start = 0
        tilt_end = self.ntilt
        
        if INITIALIZATION_FROM_SINGLE_SERIE:
            tilt_end = self.tlts[0]
            
        log.info('Initialization of the (%i) track positions (from tilts #%i to %i)' %(nobj,tilt_start,tilt_end-1))

        ntilt = tilt_end - tilt_start
        
        projection = numpy.resize(self.P[2*tilt_start:2*tilt_end,:4],(ntilt,2,4))
        mask = self.Emask[tilt_start:tilt_end,:]
        
        a = [0.0]*3*nterm*nobj

        pool = multiprocessing.Pool(processes=10)

        tasks = [[iobj,mask[:,iobj],projection,self.tracks] for iobj in range(nobj)]
        it = pool.imap(__find_obj_star__,tasks)

        for iobj,X in enumerate(it):
            if X==None: continue
            a[iobj] = X[0]
            a[nterm*nobj + iobj] = X[1]
            a[2*nterm*nobj + iobj] = X[2]
            iobj += 1
            log.info('iobj=%i   X=%.2f  Y=%.2f  Z=%.2f' %(iobj,X[0],X[1],X[2]))

        pool.close()
        pool.terminate()
     
        return a
    

    def initializeBCoefficients( self, zeros=[], skips=[] ):  # to be overriden
        """Estimation of the b coefficients from the projection map P of
        the TxBRContourAlign object. A projection at a given tilt angle can be
        forced zero if the index is included in the list zeros
        """

        raise NotImplementedError
    
    
    def initializeBScalingCoefficients(self,):  # to be overriden
        """Estimation of the scaling coefficients for the projection map P of
        the TxBRContourAlign object.
        """

        raise NotImplementedError


    def initializeCCoefficients(self, zeros=[]):
        """Estimation of the c coefficients from the tracks of
        the TxBRContourAlign object. The traces at a given tilt angle can be
        forced zero if the index is included in the list zeros
        """

        nterm = self.nn3
        ntilt = self.ntilt
        npatch = self.npatch

        c = [0]*2*nterm*ntilt*npatch

        for ii in range(2):
            for k in range(nterm):
                for itilt in range(ntilt):
                    for ipatch in range(npatch):
                        if ipatch in self.points and k>0: continue
                        if itilt in zeros:
                            c[ii*nterm*ntilt*npatch + k*ntilt*npatch + itilt*npatch + ipatch] = 0.0
                        else:
                            c[ii*nterm*ntilt*npatch + k*ntilt*npatch + itilt*npatch + ipatch] = self.tracks[itilt,ipatch,ii,k]

        return c


    def initializeDCoefficients(self):
        """
        Estimation of the d coefficients.
        """

        nterm = self.nn4
        ntilt = self.ntilt
        npatch = self.npatch

        d = [0.0]*nterm*ntilt*npatch

        for itilt in range(ntilt):
            for ipatch in range(npatch):
                d[itilt*npatch + ipatch] = float(itilt)/(ntilt-1)

        return d


    def initializeASkipping(self):

        nobj = self.npatch
        ntilt = self.ntilt
        
        skip_a = [0]*len(self.a)

        nn1 = util.numberOfTerms(0,dim=2)
        nterm = self.nn1

        for iobj in range(nobj):
            if iobj in self.points:
                for iterm in range(nn1,nterm):
                    skip_a[iterm*nobj + iobj] = 1
                    skip_a[nterm*nobj + iterm*nobj + iobj] = 1
                    skip_a[2*nterm*nobj + iterm*nobj + iobj] = 1

        for iobj in range(nobj):
            if iobj in self.lines:
                for iterm in range(nterm):
                    if self.resid.pp1[1][iterm]!=0:
                        skip_a[iterm*nobj + iobj] = 1
                        skip_a[nterm*nobj + iterm*nobj + iobj] = 1
                        skip_a[2*nterm*nobj + iterm*nobj + iobj] = 1
                        self.a[iterm*nobj + iobj] = 0.0
                        self.a[nterm*nobj + iterm*nobj + iobj] = 0.0
                        self.a[2*nterm*nobj + iterm*nobj + iobj] = 0.0

        return skip_a



    def initializeBSkipping(self, projection=None):  # to be overridden

        raise NotImplementedError


    def initializeCSkipping(self):

        skip_c = [0]*len(self.c)
        
        nobj = self.npatch
        ntilt = self.ntilt

        nn3 = util.numberOfTerms(0,dim=1)
        nterm = self.nn3

        for ii in range(2):
            for k in range(nn3,nterm):
                for itilt in range(ntilt):
                    for iobj in range(nobj):
                        if iobj in self.points:
                            skip_c[ii*nterm*ntilt*nobj + k*ntilt*nobj + itilt*nobj + iobj] = 1

        return skip_c


    def initializeDSkipping(self):

        skip_d = [0]*len(self.d)

        return skip_d
    

    def extractACoefficients(self,itilt,ipatch):
        """Extract the core coefficients of the cost function for a given tilt and
        a given patch.
        """

        ntilt = self.ntilt
        npatch = self.npatch

        nn1 = self.nn1

        coeffs = []

        if numpy.all(self.Emask[:,ipatch]==0):
            return [0]*3*nn1

        for ii in range(3):
            for k in range(nn1):
                coeffs.append(self.a[ii*nn1*npatch + k*npatch + ipatch])

        return numpy.array(coeffs)


    def extractBCoefficients(self,itilt,ipatch):
        """Extract the core coefficients of the cost function for a given tilt and
        a given patch.
        """
        
        raise NotImplementedError


    def extractCCoefficients(self,itilt,ipatch):
        """Extract the core coefficients of the cost function for a given tilt and
        a given patch.
        """

        ntilt = self.ntilt
        npatch = self.npatch

        nn3 = self.nn3

        coeffs = []

        if numpy.all(self.Emask[:,ipatch]==0):
            return [0]*2*nn3

        for ii in range(2):
            for k in range(nn3):
                coeffs.append(self.c[ii*nn3*ntilt*npatch + k*ntilt*npatch + itilt*npatch + ipatch])

        return numpy.array(coeffs)


    def extractDCoefficients(self,itilt,ipatch):
        """Extract the core coefficients of the cost function for a given tilt and
        a given patch.
        """

        ntilt = self.ntilt
        npatch = self.npatch

        nn4 = self.nn4

        coeffs = []

        if numpy.all(self.Emask[:,ipatch]==0):
            return [0]*nn4

        for k in range(nn4):
            coeffs.append(self.d[k*ntilt*npatch + itilt*npatch + ipatch])

        return numpy.array(coeffs)


    def extract_coefficients(self,itilt,ipatch):
        """Extract the core coefficients of the cost function for a given tilt and
        a given patch.
        """

        a = self.extractACoefficients(itilt,ipatch)
        b = self.extractBCoefficients(itilt,ipatch)
        c = self.extractCCoefficients(itilt,ipatch)
        d = self.extractDCoefficients(itilt,ipatch)

        return numpy.concatenate((a,b,c,d))
    
    
    def extract3DPointMarkers(self):
        '''Extract the XYZ locations of the point-like structures.'''
        
        XYZ = numpy.resize(self.a,(3,self.nn1,self.npatch))
        XYZ = XYZ[:,0,self.points].transpose()
        
        return XYZ
    
    
    def getOrthogonalParameters(self):
        
        raise NotImplementedError
    
    
    def getMagnificationParameters(self):
        
        raise NotImplementedError
    

    def relaxProjection(self,order=None):  # To be overridden

        return self.b
    
    
    def evaluateScalingParameters(self):
        
        return ( self.b, self.b_scaling )


    def storeParameters(self, var, withoutA, withoutB, withoutC, withoutD):

        withA = not withoutA
        withB = not withoutB
        withC = not withoutC
        withD = not withoutD

        if withA:
            self.a = var[:len(self.a)].tolist()
            var = numpy.delete(var,numpy.arange(len(self.a)))

        if withB:
            self.b = var[:len(self.b)].tolist()
            var = numpy.delete(var,numpy.arange(len(self.b)))

        if withC:
            self.c = var[:len(self.c)].tolist()
            var = numpy.delete(var,numpy.arange(len(self.c)))

        if withD:
            self.d = var[:len(self.d)].tolist()
            var = numpy.delete(var,numpy.arange(len(self.d)))


    def joinParameters(self, withoutA, withoutB, withoutC, withoutD):

        withA = not withoutA
        withB = not withoutB
        withC = not withoutC
        withD = not withoutD

        var = []

        if withA: var += self.a
        if withB: var += self.b
        if withC: var += self.c
        if withD: var += self.d

        return numpy.array(var)


    def test(self,derivative=True,hessian=False):

        import random

        n1,n2,n3,n4 = self.nvar[0],self.nvar[1],self.nvar[2],self.nvar[3]

        var1 = [random.random() for i in range(n1)]
        var2 = [random.random() for i in range(n2)]
        var3 = [random.random() for i in range(n3)]
        var4 = [random.random() for i in range(n4)]

        var = var1 + var2 + var3 + var4
        nn = n1 + n2 + n3 + n4

        self.a = var1
        self.b = var2
        self.c = var3
        self.d = var4

        [E0,der0,hess0] = self.evaluate_error_from_vars(var1,var2,var3,var4)

        h = 0.001

        nder = self.numberOfFreeVariables()
        names = []

        if not self.fix_a:
            names += self.variables[0]
        if not self.fix_b:
            names += self.variables[1]
        if not self.fix_c:
            names += self.variables[2]
        if not self.fix_d:
            names += self.variables[3]

        nder = self.numberOfFreeVariables()

        var_p_delta = [[var[i] for i in range(nn)] for ider in range(nder)]

        ider = 0

        if not self.fix_a:
            for i1 in range(n1):
                var_p_delta[ider][i1] += h
                ider += 1
        if not self.fix_b:
            for i2 in range(n2):
                var_p_delta[ider][n1+i2] += h
                ider += 1
        if not self.fix_c:
            for i3 in range(n3):
                var_p_delta[ider][n1+n2+i3] += h
                ider += 1
        if not self.fix_d:
            for i4 in range(n4):
                var_p_delta[ider][n1+n2+n3+i4] += h
                ider += 1

        EPS = 0.000001

        # Check the derivatives

        if derivative:

            der_ = [0]*nder

            for ider in range(nder):
                v1_ = var_p_delta[ider][:n1]
                v2_ = var_p_delta[ider][n1:n1+n2]
                v3_ = var_p_delta[ider][n1+n2:n1+n2+n3]
                v4_ = var_p_delta[ider][n1+n2+n3:nn]
                (E,der,hess) = self.evaluate_error_from_vars(v1_,v2_,v3_,v4_)
                der_[ider] = (E-E0)/h

            for i in range(nder):
                log.info('%3i - %s:  der=%-6.2f  app=%-6.2f  delta=%-5.2e' %(i,names[i],der0[i],der_[i],abs((der_[i]-der0[i])/(der[i]+EPS))))

        # Check the hessians

        if hessian:

            hess_ = [[0]*nder for ider in range(nder)]

            for ider1 in range(nder):
                for ider2 in range(nder):
                    v1_ = var_p_delta[ider2][:n1]
                    v2_ = var_p_delta[ider2][n1:n1+n2]
                    v3_ = var_p_delta[ider2][n1+n2:n1+n2+n3]
                    v4_ = var_p_delta[ider2][n1+n2+n3:nn]
                    (E2,der2,hess2) = self.evaluate_error_from_vars(v1_,v2_,v3_,v4_)
                    hess_[ider1][ider2] = (der2[ider1]-der0[ider1])/h

            for ider1 in range(nder):
                for ider2 in range(nder):
                    print '%3i/%3i - %s/%s:  hess=%-6.2f  app=%-6.2f  delta=%-5.2e' %(ider1,ider2,names[ider1],names[ider2],hess0[ider1][ider2],hess_[ider1][ider2],abs((hess_[ider1][ider2]-hess0[ider1][ider2])/(hess_[ider1][ider2]+EPS)))


