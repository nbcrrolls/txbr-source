import numpy

def bcc():
    '''A method to generate a Body Centered Cubic lattice.'''
    
    a = 1.0
    
    a1 = numpy.array((-0.5,0.5,0.5))*a
    a2 = numpy.array((0.5,-0.5,0.5))*a
    a3 = numpy.array((0.5,0.5,-0.5))*a
    
    X0 = numpy.zeros((1,3))
    
    n = 3
    
    X = []

    for i1 in range(2*n):
        for i2 in range(2*n):
            for i3 in range(2*n):
                X.append(X0 + i1*a1 + i2*a2 + i3*a3)
                X.append(X0 - i1*a1 + i2*a2 + i3*a3)
                X.append(X0 + i1*a1 - i2*a2 + i3*a3)
                X.append(X0 + i1*a1 + i2*a2 - i3*a3)
                X.append(X0 - i1*a1 - i2*a2 + i3*a3)
                X.append(X0 - i1*a1 + i2*a2 - i3*a3)
                X.append(X0 + i1*a1 - i2*a2 - i3*a3)
                X.append(X0 - i1*a1 - i2*a2 - i3*a3)
                
    X = numpy.row_stack(X)
            
    indices, = numpy.where((X[:,0]>=0.0)*(X[:,0]<n)*(X[:,1]>=0.0)*(X[:,1]<n)*(X[:,2]>=0.0)*(X[:,2]<n))
    
    return X[indices]


def fcc():
    '''A method to generate a Faced Centered Cubic lattice.'''

    a = 1.0
    
    a1 = numpy.array((0.5,0.5,0.0))*a
    a2 = numpy.array((0.0,0.5,0.5))*a
    a3 = numpy.array((0.5,0.0,0.5))*a
    
    X0 = numpy.zeros((1,3))
    
    n = 2
    
    X = []

    for i1 in range(2*n):
        for i2 in range(2*n):
            for i3 in range(2*n):
                X.append(X0 + i1*a1 + i2*a2 + i3*a3)
                X.append(X0 - i1*a1 + i2*a2 + i3*a3)
                X.append(X0 + i1*a1 - i2*a2 + i3*a3)
                X.append(X0 + i1*a1 + i2*a2 - i3*a3)
                X.append(X0 - i1*a1 - i2*a2 + i3*a3)
                X.append(X0 - i1*a1 + i2*a2 - i3*a3)
                X.append(X0 + i1*a1 - i2*a2 - i3*a3)
                X.append(X0 - i1*a1 - i2*a2 - i3*a3)
                
    X = numpy.row_stack(X)
            
    indices, = numpy.where((X[:,0]>=0.0)*(X[:,0]<n)*(X[:,1]>=0.0)*(X[:,1]<n)*(X[:,2]>=0.0)*(X[:,2]<n))
    
    return X[indices]


def rhombohedral(x):
    '''A method to generate a rhombohedral lattice.'''

    alpha_in_rad = numpy.arccos(x*(3.0+2.0*x)/(2.0+(1.0+x)**2))
    alpha_in_deg = numpy.degrees(alpha_in_rad)
    print 'alpha=%f rad    %f deg' %(alpha_in_rad,alpha_in_deg)
    #peak = numpy.sqrt(((1.0+x)**2+2.0*x**2)/((1.0+2.0*x)**2+4.0*x**2))
    peak = numpy.sqrt(((1.0+x)**2+2.0*x**2))   # a1/(a1-a2)
    print 'peak: %f   %f' %(peak,1.0/peak)
    
    a = 1.0
    
    a1 = numpy.array((1.0+x,x,x))*a
    a2 = numpy.array((x,1.0+x,x))*a
    a3 = numpy.array((x,x,1.0+x))*a
    
    X0 = numpy.zeros((1,3))
    
    n = 4
    
    X = []

    for i1 in range(2*n):
        for i2 in range(2*n):
            for i3 in range(2*n):
                X.append(X0 + i1*a1 + i2*a2 + i3*a3)
                X.append(X0 - i1*a1 + i2*a2 + i3*a3)
                X.append(X0 + i1*a1 - i2*a2 + i3*a3)
                X.append(X0 + i1*a1 + i2*a2 - i3*a3)
                X.append(X0 - i1*a1 - i2*a2 + i3*a3)
                X.append(X0 - i1*a1 + i2*a2 - i3*a3)
                X.append(X0 + i1*a1 - i2*a2 - i3*a3)
                X.append(X0 - i1*a1 - i2*a2 - i3*a3)
                
    X = numpy.row_stack(X)
            
    indices, = numpy.where((X[:,0]>=0.0)*(X[:,0]<n)*(X[:,1]>=0.0)*(X[:,1]<n)*(X[:,2]>=0.0)*(X[:,2]<n))
    
    return X[indices]
