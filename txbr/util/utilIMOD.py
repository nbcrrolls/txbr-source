import numpy
import sympy
import sympy.matrices


def slicerNormalAxis( phi1, phi2, phi3 ):
    '''This routine returns, given the three angles in degrees from an IMOD slicer view,
    the normal vector to the sliced plane'''

    #print "----------------------------------------------"
    #print "[phi1,phi2,phi3]=[%.2f,%.2f,%.2f]" %(phi1,phi2,phi3)

    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    z = sympy.Symbol('z')

    Rx = sympy.matrices.Matrix([[1,0,0],[0,sympy.cos(x),sympy.sin(x)],[0,-sympy.sin(x),sympy.cos(x)]])
    Ry = sympy.matrices.Matrix([[sympy.cos(y),0,sympy.sin(y)],[0,1,0],[-sympy.sin(y),0,sympy.cos(y)]])
    Rz = sympy.matrices.Matrix([[sympy.cos(z),sympy.sin(z),0],[-sympy.sin(z),sympy.cos(z),0],[0,0,1]])

    ez = sympy.matrices.Matrix([[0],[0],[1]])
    u = Rx*Ry*Rz*ez

    phi1 = -phi1    # IMOD convention
    phi2 = phi2
    phi3 = phi3

    v = u.subs(x,numpy.radians(phi1)).subs(y,numpy.radians(phi2)).subs(z,numpy.radians(phi3))
    v = [float(v[0].evalf()),float(v[1].evalf()),float(v[2].evalf())]
    
    v = numpy.asarray(v)
    
    #print "[%.3f,%.3f,%.3f]" %(v[0],v[1],v[2])

    return v


    