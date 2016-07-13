import sys
import logging
import numpy
import numpy.linalg
import numpy.random
import optimization
import powerseq
import poly

log = logging.getLogger()

AFFINE_TRANSFORMATION = 'affine'
CONFORMAL_TRANSFORMATION = 'conformal'
PROJECTIVE_TRANSFORMATION = 'projective'
POLYNOMIAL_TRANSFORMATION = 'polynomial'

class ProjectiveTransformation2D:
    '''A class to describe a projective transformation in the 2D space.'''

    def __init__( self, M, with_center=False ):
        '''
        Constructor for the ProjectiveTransformation2D class.
        Variable M: a numpy array (or a tuple) containing the coefficents of the 2D 
        transformation. The coefficients should have the following order 
        [tx,mxx,mxy,ty,myx,myy,1,lx,ly] if with_center is false
        [cx,mxx,mxy,cy,myx,myy,1,lx,ly] if with_center is true
        Variable with_center: specifies if the translations part of the transformation
        is defined directly or from a fixed point position.
        '''
        
        M = numpy.resize(M,(3,3))   
        self.M = numpy.asarray(M[:,1:],dtype=float)
        
        self.with_center = with_center
        
        if with_center:
            self.C = numpy.asarray(M[:2,0],dtype=float)
            self.T = self.C*(1.0 + numpy.dot(M[2,1:],self.C)) - numpy.dot(self.M,self.C)
        else:
            self.C = numpy.zeros((2))
            self.T = numpy.asarray(M[:2,0],dtype=float)
            
            
    def addTranslation(self, shift):
        
        shift = numpy.asarray(shift)
        
        self.T[:2] += shift
        self.M[:2,0] += self.M[2,0]*shift
        self.M[:2,1] += self.M[2,1]*shift
        
        
    def forward(self,x):
        '''Calculates the image of x.
        Variable x: a set of 3D points.
        Returns: The image of x with this 3D transformation.
        '''
        
        x = numpy.asarray(x)
        
        if x.size==0:
            return x
        
        n = x.size/2
        
        s0 = x.shape
        s = (n,2)
        
        x.resize(s)
        
        scaling = numpy.tensordot(x,self.M[2:,:],(1,1)) + 1.0
        y = numpy.tensordot(x,self.M[:2,:],(1,1)) + self.T
        
        y = y/scaling
        
        y.resize(s0)
        
        return y
    
    
    def __repr__(self):
        
        str = "2D Projective Transf.: T=[%.1f,%.1f]  M=[[%.1f,%.1f],[%.1f,%.1f]] Scaling: [%.1f,%.1f]" \
                %tuple(numpy.row_stack((self.T,self.M)).ravel())
                    
        return str
            
            
class PolynomialTransformation2D:
    '''A class to describe a polynomial transformation in the 3D space.
    '''
    
    def __init__( self, order, M, with_center=False ):
        '''
        Constructor for the PolynomialTransformation2D class.
        Variable order: The order of the polynomial transformation
        Variable M: a numpy array (or a tuple) containing the coefficents of the 2D 
        transformation.
        Variable with_center: specifies if the translations part of the transformation
        is defined directly or from a fixed point position.
        '''        
    
        self.order = order
        self.nn = powerseq.numberOfTerms(order,dim=2)
        
        M = numpy.asarray(M)

        if M.size==1*self.nn:
            Y = numpy.zeros((self.nn))
            Y[2] = 1.0
            M = numpy.resize(M,(1,self.nn))
            M = numpy.row_stack((M,Y))
        
        if M.size!=2*self.nn:
            raise Exception("The number of elements in M does not match the order of the polynomial transformation")
       
        M = numpy.resize(M,(2,self.nn))
        
        self.M = numpy.asarray(M[:,1:],dtype=float)
        
        self.variables = ['X','Y']
        self.powers = numpy.array(powerseq.powerOrder(self.order,dim=2)).T

        self.with_center = with_center
        
        if with_center:
            self.C = numpy.asarray(M[:,0],dtype=float)
            self.resetTranslation(self.C,self.C)
        else:
            self.C = numpy.zeros((3))
            self.T = numpy.asarray(M[:,0],dtype=float)
            self.evaluateCoordinatePolynoms()
            
            
    def evaluateCoordinatePolynoms( self ):
        '''Explicitely evaluate the 1D polynomial function for each of the 
        3D coordinates'''
        
        M_ = numpy.column_stack((self.T,self.M))
                
        self.px = poly.PolyM(self.variables, M_[0,:], self.powers)
        self.py = poly.PolyM(self.variables, M_[1,:], self.powers)
       
       
    def resetTranslation( self, X, Ximg ):
        '''Recalculate the translation part of this transformation given a point and its image.
        Variable X: a 3D point
        Variable Ximg: a 3D point, image of X by this transformation'''
        
        X = numpy.asarray(X)
        Ximg = numpy.asarray(Ximg)
        
        M_ = numpy.column_stack((numpy.zeros((3)),self.M))
        
        px_ = poly.PolyM(self.variables, M_[0,:], self.powers)
        py_ = poly.PolyM(self.variables, M_[1,:], self.powers)
        
        Img = numpy.concatenate((px_.eval(X),py_.eval(X)))
        
        self.T = Ximg - Img
        
        self.evaluateCoordinatePolynoms()
        
        
    def coeffs(self):
        
        return numpy.column_stack((self.T,self.M))
        

    def forward(self,x):
        '''Calculates the image of x.
        Variable x: a set of 3D points.
        Returns: The image of x with this 3D transformation.
        '''
        
        x = numpy.asarray(x)
        
        if x.size==0:
            return x
        
        s0 = x.shape
        s = (x.size/2,2)
        
        x.resize(s)
        
        y = numpy.column_stack((self.px.eval(x),self.py.eval(x)))
        
        return numpy.resize(y,s0)


    def reverse(self,x):
        '''Calculates the inverse image of x.
        Variable x: a set of 3D points.
        Returns: The image of x with the inverse of this 3D polynomial transformation.
        '''
        
        raise NotImplementedError


    def inv(self):
        '''Returns the inverse of this 3D polynomial transformation.'''
        
        raise NotImplementedError


    def compose(self,g):
        '''Compose this 3D transformation with another one g.'''
        
        order = self.order*g.order
        nn = powerseq.numberOfTerms(order,dim=2)
        powers = numpy.array(powerseq.powerOrder(order,dim=2)).T
        
        coord_dict = {'X': g.px,'Y': g.py}
        
        px = self.px.compose(coord_dict)
        py = self.py.compose(coord_dict)
        
        M = ( px.extract_coefficients(powers),\
              py.extract_coefficients(powers) )
        
        M = numpy.row_stack(M)
        
        return PolynomialTransformation2D(order,M)


    def __repr__(self):
        
        str = "2D Polynomial Transf.: (order %i)\n%s" \
                    %(self.order,numpy.column_stack((self.T,self.M)))
                    
        return str


class AffineTransformation2D:
    '''A class to describe an affine transformation in the 2D space.'''
    

    def __init__( self, M, with_center=False ):
        '''
        Constructor for the AffineTransformation2D class.
        Variable M: a numpy array (or a tuple) containing the coefficents of the 2D 
        transformation. The coefficients should have the following order 
        [tx,mxx,mxy,ty,myx,myy] if with_center is false
        [cx,mxx,mxy,cy,myx,myy] if with_center is true
        Variable with_center: specifies if the translations part of the transformation
        is defined directly or from a fixed point position.
        '''
        
        M = numpy.resize(M,(2,3))   
        self.M = numpy.asarray(M[:,1:],dtype=float)
        
        self.with_center = with_center
        
        if with_center:
            self.C = numpy.asarray(M[:,0],dtype=float)
            self.T = self.C - numpy.dot(self.M,self.C)
        else:
            self.C = numpy.zeros((2))
            self.T = numpy.asarray(M[:,0],dtype=float)
            
        try:
            self.Minv = numpy.linalg.pinv(self.M)
        except numpy.linalg.linalg.LinAlgError:
            print "Affine Transformation is not invertible!"
            print "Singular Matrix: [[%.2f,%.2f],[%.2f,%.2f]]" %tuple(self.M.ravel())
        
    
    def forward(self,x):
        '''Calculates the image of x.
        Variable x: a set of 3D points.
        Returns: The image of x with this 3D transformation.
        '''
        
        x = numpy.asarray(x)
        
        if x.size==0:
            return x
        
        n = x.size/2
        
        s0 = x.shape
        s = (n,2)

        x.resize(s,refcheck=False)
        
        y = numpy.tensordot(x,self.M,(1,1)) + self.T

        x.resize(s0,refcheck=False)
        y.resize(s0)
        
        return y


    def reverse(self,x):
        '''Calculates the inverse image of x.
        Variable x: a set of 3D points.
        Returns: The image of x with the inverse of this 3D transformation.
        '''
        
        if self.Minv==None:
            raise Exception("Affine Transformation is not invertible!")
        
        x = numpy.asarray(x)
        
        n = x.size/2
        
        s0 = x.shape
        s = (n,2)
        
        x.resize(s)
            
        y = numpy.tensordot(x-self.T,self.Minv,(1,1))
    
        y.resize(s0)
        
        return y


    def inv(self):
        '''Returns the inverse of this 3D transformation.'''
        
        if self.Minv==None:
            raise Exception("Affine Transformation is not invertible!")

        Minv = self.Minv
        Tinv = -numpy.dot(self.Minv,self.T)
        
        return AffineTransformation2D(numpy.column_stack((Tinv,Minv)))


    def compose(self,g):
        '''Compose this 3D transformation with anotherone g.'''
        
        M = numpy.dot(self.M,g.M)
        T = numpy.dot(self.M,g.T) + self.T
        
        return AffineTransformation2D(numpy.column_stack((T,M)))
    
    
    def __repr__(self):
        
        str = "2D Affine Transf.: T=[%.1f,%.1f]  M=[[%.1f,%.1f],[%.1f,%.1f]]" \
                %tuple(numpy.row_stack((self.T,self.M)).ravel())
                    
        return str


class Translation2D(AffineTransformation2D):
    '''A class to describe a 2D translation.
    '''

    def __init__( self, tx, ty ):
        '''
        Constructor for the Translation2D class.
        Variable tx, ty: values of the translation vector
        '''

        # Create the corresponding AffineTransformation2D

        T = numpy.array([tx,ty])
        M = numpy.eye(2,2)

        AffineTransformation2D.__init__( self, numpy.column_stack((T,M)))


    def inv(self):
        '''Returns the inverse of this 3D rotation.'''

        return AffineTransformation2D( -self.T[0], -self.T[1] )


    def __repr__(self):

        return "2D Translation: T=[%.1f,%.1f]" %(self.T[0],self.T[1])


class Rotation2D(AffineTransformation2D):
    '''A class to describe a 2D rotation.
    '''
    
    def __init__( self, center, theta, in_radians=True ):
        '''
        Constructor for the Rotation2D class.
        Variable center: a numpy array (or a tuple/list) containing the 2 rotation center coefficents 
        Variable theta: angle for the rotation in radians if in_radians is True, otherwise in degrees
        '''
        
        self.center = numpy.asarray(center)
        
        if in_radians:
            self.theta = theta
        else:
            self.theta = numpy.radians(theta)
        
        # Create the corresponding AffineTransformation2D
        
        T = numpy.zeros((2))
        M = numpy.zeros((2,2))

        cos_theta = numpy.cos(self.theta)
        sin_theta = numpy.sin(self.theta)

        M[0,0] = cos_theta
        M[0,1] = -sin_theta

        M[1,0] = sin_theta
        M[1,1] = cos_theta

        T = self.center - numpy.dot( M, self.center )
        
        AffineTransformation2D.__init__( self, numpy.column_stack((T,M)))
        
        
    def inv(self):
        '''Returns the inverse of this 3D rotation.'''

        return Rotation2D( self.center, -self.theta )
        
        
    def __repr__(self):
                            
        return "2D Rotation: Center=[%.1f,%.1f]   Angle=%.2f deg" \
                %tuple(numpy.append(self.center,numpy.degrees(self.theta)))
                


class PolynomialTransformation3D:
    '''A class to describe a polynomial transformation in the 3D space.
    '''
    
    def __init__( self, order, M, with_center=False ):
        '''
        Constructor for the PolynomialTransformation3D class.
        Variable order: The order of the polynomial transformation
        Variable M: a numpy array (or a tuple) containing the coefficents of the 3D 
        transformation.
        Variable with_center: specifies if the translations part of the transformation
        is defined directly or from a fixed point position.
        '''        
    
        self.order = order
        self.nn = powerseq.numberOfTerms(order,dim=3)
        
        M = numpy.asarray(M)

        if M.size==1*self.nn:
            Y = numpy.zeros((self.nn))
            Y[2] = 1.0
            Z = numpy.zeros((self.nn))
            Z[3] = 1.0
            M = numpy.resize(M,(1,self.nn))
            M = numpy.row_stack((M,Y,Z))
        
        if M.size==2*self.nn:
            Z = numpy.zeros((self.nn))
            Z[3] = 1.0
            M = numpy.resize(M,(2,self.nn))
            M = numpy.row_stack((M,Z))
        
        if M.size!=3*self.nn:
            raise Exception("The number of elements in M does not match the order of the polynomial transformation")
       
        M = numpy.resize(M,(3,self.nn))
        
        self.M = numpy.asarray(M[:,1:],dtype=float)
        
        self.variables = ['X','Y','Z']
        self.powers = numpy.array(powerseq.powerOrder(self.order,dim=3)).T

        self.with_center = with_center
        
        if with_center:
            self.C = numpy.asarray(M[:,0],dtype=float)
            self.resetTranslation(self.C,self.C)
        else:
            self.C = numpy.zeros((3))
            self.T = numpy.asarray(M[:,0],dtype=float)
            self.evaluateCoordinatePolynoms()
            
            
    def evaluateCoordinatePolynoms( self ):
        '''Explicitely evaluate the 1D polynomial function for each of the 
        3D coordinates'''
        
        M_ = numpy.column_stack((self.T,self.M))
                
        self.px = poly.PolyM(self.variables, M_[0,:], self.powers)
        self.py = poly.PolyM(self.variables, M_[1,:], self.powers)
        self.pz = poly.PolyM(self.variables, M_[2,:], self.powers)
       
       
    def resetTranslation( self, X, Ximg ):
        '''Recalculate the translation part of this transformation given a point and its image.
        Variable X: a 3D point
        Variable Ximg: a 3D point, image of X by this transformation'''
        
        X = numpy.asarray(X)
        Ximg = numpy.asarray(Ximg)
        
        M_ = numpy.column_stack((numpy.zeros((3)),self.M))
        
        px_ = poly.PolyM(self.variables, M_[0,:], self.powers)
        py_ = poly.PolyM(self.variables, M_[1,:], self.powers)
        pz_ = poly.PolyM(self.variables, M_[2,:], self.powers)
        
        Img = numpy.concatenate((px_.eval(X),py_.eval(X),pz_.eval(X)))
        
        self.T = Ximg - Img
        
        self.evaluateCoordinatePolynoms()
        
        
    def coeffs(self):
        
        return numpy.column_stack((self.T,self.M))
        

    def forward(self,x):
        '''Calculates the image of x.
        Variable x: a set of 3D points.
        Returns: The image of x with this 3D transformation.
        '''
        
        x = numpy.asarray(x)
        
        if x.size==0:
            return x
        
        s0 = x.shape
        s = (x.size/3,3)
        
        x.resize(s)
        
        y = numpy.column_stack((self.px.eval(x),self.py.eval(x),self.pz.eval(x)))
        
        return numpy.resize(y,s0)


    def reverse(self,x):
        '''Calculates the inverse image of x.
        Variable x: a set of 3D points.
        Returns: The image of x with the inverse of this 3D polynomial transformation.
        '''
        
        raise NotImplementedError


    def inv(self):
        '''Returns the inverse of this 3D polynomial transformation.'''
        
        raise NotImplementedError


    def compose(self,g):
        '''Compose this 3D transformation with another one g.'''
        
        order = self.order*g.order
        nn = powerseq.numberOfTerms(order,dim=3)
        powers = numpy.array(powerseq.powerOrder(order,dim=3)).T
        
        coord_dict = {'X': g.px,'Y': g.py,'Z': g.pz}
        
        px = self.px.compose(coord_dict)
        py = self.py.compose(coord_dict)
        pz = self.pz.compose(coord_dict)
        
        M = ( px.extract_coefficients(powers),\
              py.extract_coefficients(powers),\
              pz.extract_coefficients(powers) )
        
        M = numpy.row_stack(M)
        
        return PolynomialTransformation3D(order,M)


    def __repr__(self):
        
        str = "3D Polynomial Transf.: (order %i)\n%s" \
                    %(self.order,numpy.column_stack((self.T,self.M)))
                    
        return str
    
    

class AffineTransformation3D(PolynomialTransformation3D):
    '''A class to describe an affine transformation in the 3D space.
    '''

    def __init__( self, M, with_center=False ):
        '''
        Constructor for the AffineTransformation3D class.
        Variable M: a numpy array (or a tuple) containing the coefficents of the 3D 
        transformation. The coefficients should have the following order 
        [tx,mxx,mxy,mxz,ty,myx,myy,myz,tz,mzx,mzy,mzz] if with_center is false
        [cx,mxx,mxy,mxz,cy,myx,myy,myz,cz,mzx,mzy,mzz] if with_center is true
        Variable with_center: specifies if the translations part of the transformation
        is defined directly or from a fixed point position.
        '''
        
        PolynomialTransformation3D.__init__(self, 1, M, with_center) # Order of this polynomial transformation is 1
            
        try:
            self.Minv = numpy.linalg.pinv(self.M)
        except numpy.linalg.linalg.LinAlgError:
            print "Affine Transformation is not invertible!"
            print "Singular Matrix: [[%.2f,%.2f,%.2f],[%.2f,%.2f,%.2f],[%.2f,%.2f,%.2f]]" %tuple(self.M.ravel())
        

    def forward(self,x):
        '''Calculates the image of x.
        Variable x: a set of 3D points.
        Returns: The image of x with this 3D transformation.
        '''
        
        x = numpy.asarray(x)
        
        if x.size==0:
            return x
        
        n = x.size/3
        
        s0 = x.shape
        s = (n,3)
        
        x = numpy.resize(x,s)
        
        y = numpy.tensordot(x,self.M,(1,1)) + self.T
        
        x = numpy.resize(x,s0)
        y.resize(s0)
        
        return y


    def reverse(self,x):
        '''Calculates the inverse image of x.
        Variable x: a set of 3D points.
        Returns: The image of x with the inverse of this 3D transformation.
        '''
        
        if self.Minv==None:
            raise Exception("Affine Transformation is not invertible!")
        
        x = numpy.asarray(x)
        
        n = x.size/3
        
        s0 = x.shape
        s = (n,3)
        
        x.resize(s)
            
        y = numpy.tensordot(x-self.T,self.Minv,(1,1))
    
        y.resize(s0)
        
        return y


    def inv(self):
        '''Returns the inverse of this 3D transformation.'''
        
        if self.Minv==None:
            raise Exception("Affine Transformation is not invertible!")

        Minv = self.Minv
        Tinv = -numpy.dot(self.Minv,self.T)
        
        return AffineTransformation3D(numpy.column_stack((Tinv,Minv)))


    def compose(self,g):
        '''Compose this 3D transformation with another one g.'''
        
        if g.order>1:
            
            return PolynomialTransformation3D.compose(self,g)
        
        elif g.order==0:
            
            M = self.M.copy()
            T = numpy.dot(self.M,g.T) + self.T
            
            return AffineTransformation3D(numpy.column_stack((T,M)))
        
        elif g.order==1:
        
            M = numpy.dot(self.M,g.M)
            T = numpy.dot(self.M,g.T) + self.T
            
            return AffineTransformation3D(numpy.column_stack((T,M)))
        
        return None
    
    
    def angles(self):
        '''Returns angle of an approximate orthogonal transformation equivalent ot
        this affine transformation'''
        
        angles = numpy.zeros((3))
        
        if self.M[0,0]==0.0:
            angles[2] = numpy.pi/2.0
        else:
            angles[2] = numpy.arctan(self.M[0,1]/self.M[0,0])
            
        angles[1] = numpy.arctan(self.M[0,2]/numpy.sqrt(self.M[0,1]**2 + self.M[0,0]**2))
        
        angles[0] = numpy.arctan(self.M[1,2]/(self.M[0,0]*self.M[1,1]+self.M[0,1]*self.M[1,0]))
        
        return angles
    
    
    def __repr__(self):
        
        str = "3D Affine Transf.: T=[%.1f,%.1f,%.1f]  M=[[%.1f,%.1f,%.1f],[%.1f,%.1f,%.1f],[%.1f,%.1f,%.1f]]" \
                    %tuple(numpy.row_stack((self.T,self.M)).ravel())
                    
        return str
    


class Rotation3D(AffineTransformation3D):
    '''A class to describe a 3D rotation.
    '''
    
    def __init__( self, center, axis, theta, in_radians=True ):
        '''
        Constructor for the Rotation3D class.
        Variable center: a numpy array (or a tuple/list) containing the 3 rotation center coefficents 
        Variable axis: a numpy array (or a tuple/list) containing the 3 axis coefficents
        Variable theta: angle for the rotation in radians if in_radians is True, otherwise in degrees
        '''
        
        self.center = numpy.asarray(center)
        
        self.axis = numpy.asarray(axis)
        self.axis = self.axis/numpy.sqrt(numpy.dot(self.axis,self.axis))    # Normalize the axis vector
        
        if in_radians:
            self.theta = theta
        else:
            self.theta = numpy.radians(theta)
        
        # Create the corresponding AffineTransformation3D
        
        T = numpy.zeros((3))
        M = numpy.zeros((3,3))

        cos_theta = numpy.cos(self.theta)
        sin_theta = numpy.sin(self.theta)

        M[0,0] = self.axis[0]*self.axis[0]*(1-cos_theta) + cos_theta
        M[0,1] = self.axis[0]*self.axis[1]*(1-cos_theta) - self.axis[2]*sin_theta
        M[0,2] = self.axis[0]*self.axis[2]*(1-cos_theta) + self.axis[1]*sin_theta

        M[1,0] = self.axis[1]*self.axis[0]*(1-cos_theta) + self.axis[2]*sin_theta
        M[1,1] = self.axis[1]*self.axis[1]*(1-cos_theta) + cos_theta
        M[1,2] = self.axis[1]*self.axis[2]*(1-cos_theta) - self.axis[0]*sin_theta

        M[2,0] = self.axis[2]*self.axis[0]*(1-cos_theta) - self.axis[1]*sin_theta
        M[2,1] = self.axis[2]*self.axis[1]*(1-cos_theta) + self.axis[0]*sin_theta
        M[2,2] = self.axis[2]*self.axis[2]*(1-cos_theta) + cos_theta

        T = self.center - numpy.dot( M, self.center )
        
        AffineTransformation3D.__init__( self, numpy.column_stack((T,M)))
        
        
    def resetCenter( self, M, Mimg):
        '''Recalculate the center of the rotation given a point and its image.
        The rotation angle is supposed to be fixed'''
        
        log.debug("Reshift rotation to have: [%.1f,%.1f,%.1f]->[%.1f,%.1f,%.1f]" \
                %tuple(numpy.concatenate((M,Mimg))))
        
        AffineTransformation3D.resetTranslation( self, M, Mimg)
        
        oldCenter = self.center
        
        self.center = numpy.dot(numpy.linalg.pinv(numpy.eye(3)-self.M),self.T)
        
        log.debug("Center: [%.1f,%.1f,%.1f]->[%.1f,%.1f,%.1f]" \
                %tuple(numpy.concatenate((oldCenter,self.center))))
        
        
    def inv(self):
        '''Returns the inverse of this 3D rotation.'''

        return Rotation3D( self.center, self.axis, -self.theta )
        
        
    def __repr__(self):
        
        u = numpy.concatenate((self.center,self.axis))
        u = numpy.append(u,numpy.degrees(self.theta))
        u = tuple(u)
                            
        return "3D Rotation 3D: Center=[%.1f,%.1f,%.1f]   Axis=[%.1f,%.1f,%.1f]   Angle=%.2f deg" %u
    
    
class Translation3D(AffineTransformation3D):
    '''A class to describe a 3D translation.
    '''

    def __init__( self, tx, ty, tz ):
        '''
        Constructor for the Translation2D class.
        Variable tx, ty: values of the translation vector
        '''

        # Create the corresponding AffineTransformation2D

        T = numpy.array([tx,ty,tz])
        M = numpy.eye(3,3)

        AffineTransformation3D.__init__( self, numpy.column_stack((T,M)))


    def inv(self):
        '''Returns the inverse of this 3D rotation.'''

        return Translation3D( -self.T[0], -self.T[1], -self.T[2] )


    def __repr__(self):

        return "3D Translation: T=[%.1f,%.1f,%.1f]" %(self.T[0],self.T[1],self.T[2])


if __name__ == '__main__':
    
    test2d = False
    test3d = False
    test_projective = False
    
    if "2d" in sys.argv: test2d = True
    if "3d" in sys.argv: test3d = True
    if "prj" in sys.argv: test_projective = True
    
    if not test2d and not test3d and not test_projective: test3d = True
        
    print "\nTesting 'transform' module...\n"
    
    if test3d:

        M3d = numpy.array([10.0,2.0,0.0,0.0,-5.0,0.0,1.0,0.0,1.0,0.0,0.0,-1.0])
        
        t3d = AffineTransformation3D(M3d)
        t3dinv = t3d.inv()
        
        print t3d
        print "Inverse: %s" %t3dinv
        print "Compose to identity (1): %s" %t3d.compose(t3dinv)
        print "Compose to identity (2): %s" %t3dinv.compose(t3d)
        print "\n"
    
        r3d = Rotation3D([0.0,0.0,0.0],[0.0,0.0,1.0],1.0)
        r3dinv = r3d.inv()
        
        print r3d
        print "Inverse: %s" %r3dinv
        print "Compose to identity (1): %s" %r3d.compose(r3dinv)
        print "Compose to identity (2): %s" %r3dinv.compose(r3d)
        print "\n"
        
        order = 2
        M = numpy.eye(2,10,1)
        M[0,3] = 2.0
        M[0,5] = 1.0
        
        p3d = PolynomialTransformation3D(order,M)
        
        X = numpy.random.random_sample((5,3))
        
        print p3d
        print "Evaluation"
        print X
        print p3d.forward(X)
        print "Composition"
        #u = p3d.compose(t3d)
        u = p3d
        print u
        C = (10.0,10.0,10.0)
        u.resetTranslation(C,C)
        print u
        print u.forward(C)
        
    if test2d:
    
        M2d = numpy.array([10.0,2.0,0.0,-5.0,0.0,1.0])
        
        t2d = AffineTransformation2D(M2d)
        t2dinv = t2d.inv()
        
        print t2d
        print "Inverse: %s" %t2dinv
        print "Compose to identity (1): %s" %t2d.compose(t2dinv)
        print "Compose to identity (2): %s" %t2dinv.compose(t2d)
        print "\n"
    
        r2d = Rotation2D([0.0,0.0],0.1)
        r2dinv = r2d.inv()
        
        print r2d
        print "Inverse: %s" %r2dinv
        print "Compose to identity (1): %s" %r2d.compose(r2dinv)
        print "Compose to identity (2): %s" %r2dinv.compose(r2d)
        print "\n"
    
    
    if test_projective:
        
#        Mproj = numpy.array([0.0,2.0,0.0,0.0,0.0,2.0,1.0,0.1,0.1])
#        
#        proj = ProjectiveTransformation2D(Mproj)
#        print proj
        
        X = numpy.array([[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0]])
        Ximg = numpy.array([[0.0,1.0],[1.1,0.0],[1.0,1.2],[0.2,1.0]])
        
        proj = optimization.cp2tform( X, Ximg, PROJECTIVE_TRANSFORMATION )
        
        print proj
        
        print 'X'
        print X
        print 'Ximg'
        print Ximg
        print
        print 'Transformation of X'
        print proj.forward(X)
    

I2D = AffineTransformation2D([[0.0,1.0,0.0],[0.0,0.0,1.0]])
I3D = AffineTransformation3D([[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]])
