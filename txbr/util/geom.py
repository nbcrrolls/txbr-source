'''Module to define some basic geometrical elements'''

import numpy
import math


def distanceBetweenTwoLines( d1, d2 ):
    '''
    d1 and d2: two numpy arrays of shape (2,4)
    '''

    K1 = d1[:,1:]
    K2 = d2[:,1:]

    M1 = numpy.dot(numpy.linalg.pinv(K1),-d1[:,0])  # A point in line 1
    M2 = numpy.dot(numpy.linalg.pinv(K2),-d2[:,0])  # A point in line 2

    a1 = numpy.cross(d1[0,1:],d1[1,1:])    # Vector director of line 1
    a2 = numpy.cross(d2[0,1:],d2[1,1:])    # Vector director of line 2

    n = numpy.cross(a1,a2)

    d = numpy.vdot((M2-M1),n)/numpy.linalg.norm(n)

    M = numpy.zeros((6,6))

    M[0,:3] = d1[0,1:]
    M[1,:3] = d1[1,1:]
    M[2,3:] = d2[0,1:]
    M[3,3:] = d2[1,1:]
    M[4,:3] = a1
    M[4,3:] = -a1
    M[5,:3] = a2
    M[5,3:] = -a2

    Minv = numpy.linalg.pinv(M)

    B = numpy.zeros(6)

    B[0] = -d1[0,0]
    B[1] = -d1[1,0]
    B[2] = -d2[0,0]
    B[3] = -d2[1,0]
    B[4] = 0.0
    B[5] = 0.0

    H = numpy.dot(Minv[:,:2],X1) + numpy.dot(Minv[:,2:4],X2)
    H1,H2 = H[:3],H[3:]
    H = H1 + H2

    return (H1+H2)/2.0,H1,H2


class Ellipse2D:
    '''A class to describe an ellipse in 2D'''
    
    def __init__(self,cx,cy,a,b,phi): # Plane Equation: ax+by+cz+d=0
        
        self.cx = cx
        self.cy = cy
        self.a = a
        self.b = b
        self.phi = phi  # Angle should be in radians
        
    def __repr__(self):
        
        return 'Ellipse2D: Center [%.2f,%.2f]  Axes semi-length (%.1f,%.1f)  Orientation: %.2f' \
                %(self.cx,self.cy,self.a,self.b,self.phi)


class Point3D:
    '''A class to describe a point in 3D'''

    def __init__(self,x,y,z):
        
        self.x = x
        self.y = y
        self.z = z

    def __repr__( self ):
        
        return "Point3D: (%.2f,%.2f,%.2f)" %( self.x, self.y, self.z )


class Vector3D:
    '''A class to describe a vector in 3D'''

    def __init__(self,x,y,z):
        
        self.x = x
        self.y = y
        self.z = z
        
    def normalize(self):
        
        n = math.sqrt(self*self)
        
        self.x /= n
        self.y /= n
        self.z /= n
        
        return self
        
    def __radd__( self, v ):
        
        return Vector3D( self.x+v.x, self.y+v.y, self.z+v.z)

    def __rsub__( self, v ):
        
        return Vector3D( self.x-v.x, self.y-v.y, self.z-v.z)

    def __rmul__( self, v ):
        
        return self.x*v.x + self.y*v.y + self.z*v.z

    def __repr__( self ):
        
        return "Vector3D: [%.2f,%.2f,%.2f]" %( self.x, self.y, self.z )


class Plane3D:
    '''A class to describe a plane in 3D'''

    def __init__(self,a,b,c,d): # Plane Equation: ax+by+cz+d=0
        
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        
        
    def normal(self):
        
        return Vector3D(self.a,self.b,self.c).normalize()
    

    def __repr__(self):
        
       return "Plane3D: %.2f x + %.2f y + %.2f z = %.2f" %( self.a, self.b, self.c, -self.d )


class Bounds3D:
    """A class to modelize a 3D parallelogram."""

    def __init__( self ):

        self.origin = Point3D(numpy.nan, numpy.nan, numpy.nan)
        self.end = Point3D(numpy.nan, numpy.nan, numpy.nan)
        self.increment = Vector3D(1.0, 1.0, 1.0) # default

        self.scale = Vector3D(numpy.nan, numpy.nan, numpy.nan)


    def isUndefined( self ):
        '''Check if the Bounds3D parameters are defined (different from numpy.nan)'''

        p = [ self.origin.x, self.origin.y, self.origin.z,
              self.end.x, self.end.y, self.end.z,
              self.increment.x, self.increment.y, self.increment.z,
              self.scale.x, self.scale.y, self.scale.z ]

        p = numpy.array(p)

        return numpy.any(numpy.isnan(p))
    

    def getFrame( self, full=True ):

        P0 = [self.origin.x,self.origin.y,self.origin.z]
        P1 = [self.end.x,self.origin.y,self.origin.z]
        P2 = [self.origin.x,self.end.y,self.origin.z]
        P3 = [self.origin.x,self.origin.y,self.end.z]
        P4 = [self.end.x,self.end.y,self.origin.z]
        P5 = [self.end.x,self.origin.y,self.end.z]
        P6 = [self.origin.x,self.end.y,self.end.z]
        P7 = [self.end.x,self.end.y,self.end.z]

        if full:
            return numpy.column_stack((P0,P1,P2,P3,P4,P5,P6,P7))
        else:
            return numpy.column_stack((P0,P7))


    def setOrigin(self,x,y,z):

        try:
            self.origin.x = float(x)
        except TypeError:
            pass
        try:
            self.origin.y = float(y)
        except TypeError:
            pass
        try:
            self.origin.z = float(z)
        except TypeError:
            pass


    def setIncrement(self,x,y,z):

        try:
            self.increment.x = float(x)
        except TypeError:
            pass
        try:
            self.increment.y = float(y)
        except TypeError:
            pass
        try:
            self.increment.z = float(z)
        except TypeError:
            pass


    def setEnd(self,x,y,z):

        try:
            self.end.x = float(x)
        except TypeError:
            pass
        try:
            self.end.y = float(y)
        except TypeError:
            pass
        try:
            self.end.z = float(z)
        except TypeError:
            pass
                
                
if __name__ == '__main__':
    
    plane = Plane3D(1.0,2.0,3.0,4.0)
    point = Point3D(1.0,2.0,3.0)
    ellipse = Ellipse2D(1.0,2.0,3.0,4.0,.2)
    
    print plane
    print 'Normal to the plane: %s' %plane.normal()
    print point
    print ellipse