import numpy

def createAxes(X1,X2,X3,X4):
    '''Create a set of normalized axes from the four points X1,X2,X3 and X4.'''
    
    X1 = numpy.asarray(X1)
    X2 = numpy.asarray(X2)
    X3 = numpy.asarray(X3)
    X4 = numpy.asarray(X4)
    
    u = X2 - X1
    v = X3 - X1
    w = X4 - X1
    
    print "Norms: u: %f   v: %f   w: %f" %(numpy.sqrt(numpy.dot(u,u)),numpy.sqrt(numpy.dot(v,v)),numpy.sqrt(numpy.dot(w,w)))
    
    u = u/numpy.sqrt(numpy.dot(u,u))
    v = v/numpy.sqrt(numpy.dot(v,v))
    w = w/numpy.sqrt(numpy.dot(w,w))
    
    angles_in_degrees = ( numpy.degrees(numpy.arccos(numpy.dot(u,v))), \
                          numpy.degrees(numpy.arccos(numpy.dot(u,w))), \
                          numpy.degrees(numpy.arccos(numpy.dot(v,w))) )
    
    print "Points :   X1=%s" %X1
    print "           X2=%s" %X2
    print "           X3=%s" %X3
    print "           X4=%s" %X4
    
    print "Vectors:   u=%s" %(u)
    print "           v=%s" %(v)
    print "           w=%s" %(w)
    
    print "Angles: (u,v)=%.2f   (u,w)=%.2f   (v,w)=%.2f" %angles_in_degrees
    
    return (u,v,w)


def normalToPlane(X1,X2,X3):
    '''Given 3 points (X1,X2,X3) of a plane, this function returns a unitar vector normal
    to this plane'''
    
    X1 = numpy.asarray(X1)
    X2 = numpy.asarray(X2)
    X3 = numpy.asarray(X3)
    
    u = X2 - X1
    v = X3 - X1
    
    u = u/numpy.sqrt(numpy.dot(u,u))
    v = v/numpy.sqrt(numpy.dot(v,v))
    
    w = numpy.cross(u,v)
    w = w/numpy.sqrt(numpy.dot(w,w))
    
    return w
