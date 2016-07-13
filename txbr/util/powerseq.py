import types
import numpy

class Poly:
    '''A multivariate polynom class'''

    def __init__(self,coeff,order,min_order=0):
        '''Initialize a polynom function.
        Variable coeff is a numpy array of shape (nn)
        '''

        self.coeff = coeff
        self.order = order
        self.min_order = min_order

        self.nn = numpy.alen(coeff)

    def eval(self,x):

        return poleval(self.coeff,self.order,x,self.min_order)

    def X(self,x):

        return polX(self.order,x,min_order=self.min_order)

    def der(self,m=None,dim=3):

        if type(m)==str and m=='der':

            polyder = []
            for i in range(self.nn):
                m_ = [0]*self.nn
                m_[i] = 1
                polyder.append(Poly(*polder(self.coeff,self.order,m_,dim)))

        elif type(m)==str and m=='hess':

            polyhess = [[] for i in range(self.order)]
            for i in range(self.nn):
                for j in range(self.nn):
                    m_ = [0]*self.nn
                    m_[i] += 1
                    m_[j] += 1
                    polyhess[i].append(Poly(*polder(self.coeff,self.order,m_,dim)))

        else:

            (coeffder,orderder) =  polder(self.coeff,self.order,m)
            return Poly(coeffder,orderder)


    def __repr__(self):
        
        return "Poly: %s \n%s" %(self.order,self.coeff)
    


def numberOfTerms(order,dim=3):

    if dim==4: return numberOfTerms4D(order)
    if dim==3: return numberOfTerms3D(order)
    if dim==2: return numberOfTerms2D(order)
    if dim==1: return numberOfTerms1D(order)

    raise Exception, "Dimension should be 1, 2, 3 or 4."


def powerOrder(order,dim=3):

    if dim==4: return powerOrder4D(order)
    if dim==3: return powerOrder3D(order)
    if dim==2: return powerOrder2D(order)
    if dim==1: return powerOrder1D(order)

    raise Exception, "Dimension should be 1, 2, 3 or 4."


def poleval(coeff,order,x,dim=3,min_order=0):

    coeff = numpy.asarray(coeff)
    x = numpy.asarray(x)

    if len(x.shape)==1:
        x.resize((x.size,1))

    dim = x.shape[1]

    if isinstance(order,numpy.ndarray):
        dim = order.shape[0]
    if isinstance(order,list):
        dim = len(order)
    if isinstance(order,int):
        order = numpy.array([order]*dim)

    if dim==4: return poleval4D(coeff,order,x,min_order)
    if dim==3: return poleval3D(coeff,order,x,min_order)
    if dim==2: return poleval2D(coeff,order,x,min_order)
    if dim==1: return poleval1D(coeff,order,x,min_order)

    raise Exception, "Dimension should be 1, 2, 3 or 4."


def polX(order,x,dim=3,min_order=0):

    x = numpy.asarray(x)

    dim = x.shape[1]

    if isinstance(order,numpy.ndarray):
        dim = order.shape[0]
    if isinstance(order,list):
        dim = len(order)
    if isinstance(order,int):
        order = numpy.array([order]*dim)

    if dim==4: return polX4D(order,x,min_order)
    if dim==3: return polX3D(order,x,min_order)
    if dim==2: return polX2D(order,x,min_order)
    if dim==1: return polX1D(order,x,min_order)

    raise Exception, "Dimension should be 1, 2, 3 or 4."


def polder(coeff,order,m,dim=3):

    coeff = numpy.asarray(coeff)
    m = numpy.asarray(m)

    if isinstance(order,numpy.ndarray):
        dim = order.shape[0]
    if isinstance(order,list):
        dim = len(order)
    if isinstance(order,int):
        order = numpy.array([order]*dim)

    if dim==4: return polder4D(coeff,order,m)
    if dim==3: return polder3D(coeff,order,m)
    if dim==2: return polder2D(coeff,order,m)
    if dim==1: return polder1D(coeff,order,m)

    raise Exception, "Dimension should be 1, 2, 3 or 4."

def powersIndexDict(order,dim=3):

    dim = order.shape[0]

    powers = numpy.array(powerOrder(order,dim))

    dict = {}
    n = powers.shape[1]

    for i in range(n): dict[tuple(powers[:,i])] = i

    return dict


def parseOrder4D(order):

    if type(order) is types.IntType:
        n1 = order
        n2 = order
        n3 = order
        n4 = order
    elif type(order) is types.TupleType or type(order) is numpy.ndarray:
        n1 = order[0]
        n2 = order[1]
        n3 = order[2]
        n4 = order[3]
    else:
        raise TypeError('Wrong argument!')

    return [n1,n2,n3,n4]

def numberOfTerms4D(order):

    [n1,n2,n3,n4] = parseOrder4D(order)

    approximationOrder = max(n1,n2,n3,n4)

    nn = 0
    for k in range(approximationOrder+1):
        for i1 in range(k,-1,-1):
            for i2 in range(k-i1,-1,-1):
                for i3 in range(k-i1-i2,-1,-1):
                    if i1>n1: continue
                    if i2>n2: continue
                    if i3>n3: continue
                    if k-i1-i2-i3>n4: continue
                    nn += 1

    return nn

def powerOrder4D(order):

    [n1,n2,n3,n4] = parseOrder4D(order)

    approximationOrder = max(n1,n2,n3,n4)

    order_1 = []
    order_2 = []
    order_3 = []
    order_4 = []

    for k in range(approximationOrder+1):
        for i1 in range(k,-1,-1):
            for i2 in range(k-i1,-1,-1):
                for i3 in range(k-i1-i2,-1,-1):
                    if i1>n1: continue
                    if i2>n2: continue
                    if i3>n3: continue
                    if k-i1-i2-i3>n4: continue
                    order_1.append(i1)
                    order_2.append(i2)
                    order_3.append(i3)
                    order_4.append(k-i1-i2-i3)

    power_order = [order_1,order_2,order_3,order_4]

    return power_order


def polX4D(order,x,minorder):

    order = numpy.asarray(order)

    max_order = numpy.max(order)

    x1_p = numpy.vander(x[:,0],N=order[0]+1)[:,::-1]
    x2_p = numpy.vander(x[:,1],N=order[1]+1)[:,::-1]
    x3_p = numpy.vander(x[:,2],N=order[2]+1)[:,::-1]
    x4_p = numpy.vander(x[:,3],N=order[3]+1)[:,::-1]

    n = x.shape[0]
    y = numpy.zeros((n))

    n1 = order[0]
    n2 = order[1]
    n3 = order[2]
    n4 = order[3]

    X = []

    for k in range(minorder,max_order+1):
        for i1 in range(k,-1,-1):
            for i2 in range(k-i1,-1,-1):
                for i3 in range(k-i1-i2,-1,-1):
                    if i1>n1: continue
                    if i2>n2: continue
                    if i3>n3: continue
                    if k-i1-i2-i3>n4: continue
                    X.append(x1_p[:,i1]*x2_p[:,i2]*x3_p[:,i3]*x4_p[:,k-i1-i2-i3])

    return numpy.row_stack(X)


def poleval4D(coeff,order,x,minorder):

    order = numpy.asarray(order)
    coeff = numpy.asarray(coeff)

    max_order = numpy.max(order)

    x1_p = numpy.vander(x[:,0],N=order[0]+1)[:,::-1]
    x2_p = numpy.vander(x[:,1],N=order[1]+1)[:,::-1]
    x3_p = numpy.vander(x[:,2],N=order[2]+1)[:,::-1]
    x4_p = numpy.vander(x[:,3],N=order[3]+1)[:,::-1]

    n = x.shape[0]
    y = numpy.zeros((n))

    n1 = order[0]
    n2 = order[1]
    n3 = order[2]
    n4 = order[3]

    index = 0
    for k in range(minorder,max_order+1):
        for i1 in range(k,-1,-1):
            for i2 in range(k-i1,-1,-1):
                for i3 in range(k-i1-i2,-1,-1):
                    if i1>n1: continue
                    if i2>n2: continue
                    if i3>n3: continue
                    if k-i1-i2-i3>n4: continue
                    y = y + coeff[index]*x1_p[:,i1]*x2_p[:,i2]*x3_p[:,i3]*x4_p[:,k-i1-i2-i3]
                    index += 1

    return y


def polder4D(coeff,order,m):
    """Return the (m[0],m[1],m[2],m[3])th derivative of the polynomial (coeff,order) of 4 variables
    """
    order = numpy.asarray(order)
    coeff = numpy.asarray(coeff)
    m =  numpy.asarray(m)

    powerIndexDict = powersIndexDict(order)
    nn = numberOfTerms4D(order)

    orderder = order-m
    powerIndexderDict = powersIndexDict(orderder)
    nnder = numberOfTerms4D(orderder)

    coeffder = numpy.ones((nnder))

    for powers,index in powerIndexderDict.iteritems():
        n1,n2,n3,n4=powers
        for i1 in range(m[0]): coeffder[index] *= n1+i1+1
        for i2 in range(m[1]): coeffder[index] *= n2+i2+1
        for i3 in range(m[2]): coeffder[index] *= n3+i3+1
        for i4 in range(m[3]): coeffder[index] *= n4+i4+1
        try:
            index0 = powerIndexDict[(n1+m[0],n2+m[1],n3+m[2],n4+m[3])]
            coeffder[index] *= coeff[index0]
        except KeyError:
            coeffder[index] *= 0.0

    return (coeffder,orderder)


def parseOrder3D(order):

    if type(order) is types.IntType:
        n1 = order
        n2 = order
        n3 = order
    elif type(order) is types.TupleType or type(order) is numpy.ndarray:
        n1 = order[0]
        n2 = order[1]
        n3 = order[2]
    else:
        raise TypeError('Wrong argument!')

    return [n1,n2,n3]


def numberOfTerms3D(order):

    [n1,n2,n3] = parseOrder3D(order)

    approximationOrder = max(n1,n2,n3)

    nn = 0
    for k in range(approximationOrder+1):
        for i1 in range(k,-1,-1):
            for i2 in range(k-i1,-1,-1):
                if i1>n1: continue
                if i2>n2: continue
                if k-i1-i2>n3: continue
                nn += 1

    return nn


def powerOrder3D(order):

    [n1,n2,n3] = parseOrder3D(order)

    approximationOrder = max(n1,n2,n3)

    order_1 = []
    order_2 = []
    order_3 = []

    for k in range(approximationOrder+1):
        for i1 in range(k,-1,-1):
            for i2 in range(k-i1,-1,-1):
                if i1>n1: continue
                if i2>n2: continue
                if k-i1-i2>n3: continue
                order_1.append(i1)
                order_2.append(i2)
                order_3.append(k-i1-i2)

    power_order = [order_1,order_2,order_3]

    return power_order


def polX3D(order,x,min_order):

    order = numpy.asarray(order)
    max_order = numpy.max(order)

    x1_p = numpy.vander(x[:,0],N=order[0]+1)[:,::-1]
    x2_p = numpy.vander(x[:,1],N=order[1]+1)[:,::-1]
    x3_p = numpy.vander(x[:,2],N=order[2]+1)[:,::-1]

    n = x.shape[0]
    y = numpy.zeros((n))

    n1 = order[0]
    n2 = order[1]
    n3 = order[2]

    X = []

    for k in range(min_order,max_order+1):
        for i1 in range(k,-1,-1):
            for i2 in range(k-i1,-1,-1):
                if i1>n1: continue
                if i2>n2: continue
                if k-i1-i2>n3: continue
                X.append(x1_p[:,i1]*x2_p[:,i2]*x3_p[:,k-i1-i2])

    return numpy.row_stack(X)


def poleval3D(coeff,order,x,min_order):

    order = numpy.asarray(order)
    coeff = numpy.asarray(coeff)

    max_order = numpy.max(order)

    x1_p = numpy.vander(x[:,0],N=order[0]+1)[:,::-1]
    x2_p = numpy.vander(x[:,1],N=order[1]+1)[:,::-1]
    x3_p = numpy.vander(x[:,2],N=order[2]+1)[:,::-1]

    n = x.shape[0]
    y = numpy.zeros((n))

    n1 = order[0]
    n2 = order[1]
    n3 = order[2]

    index = 0
    for k in range(min_order,max_order+1):
        for i1 in range(k,-1,-1):
            for i2 in range(k-i1,-1,-1):
                if i1>n1: continue
                if i2>n2: continue
                if k-i1-i2>n3: continue
                y = y + coeff[index]*x1_p[:,i1]*x2_p[:,i2]*x3_p[:,k-i1-i2]
                index += 1

    return y


def polder3D(coeff,order,m):
    """Return the (m[0],m[1],m[2])th derivative of the polynomial (coeff,order) of 3 variables
    """
    order = numpy.asarray(order)
    coeff = numpy.asarray(coeff)
    m =  numpy.asarray(m)

    powerIndexDict = powersIndexDict(order)
    nn = numberOfTerms3D(order)

    orderder = order-m
    powerIndexderDict = powersIndexDict(orderder)
    nnder = numberOfTerms3D(orderder)

    coeffder = numpy.ones((nnder))

    for powers,index in powerIndexderDict.iteritems():
        n1,n2,n3=powers
        for i1 in range(m[0]): coeffder[index] *= n1+i1+1
        for i2 in range(m[1]): coeffder[index] *= n2+i2+1
        for i3 in range(m[2]): coeffder[index] *= n3+i3+1
        try:
            index0 = powerIndexDict[(n1+m[0],n2+m[1],n3+m[2])]
            coeffder[index] *= coeff[index0]
        except KeyError:
            coeffder[index] *= 0.0

    return (coeffder,orderder)


def parseOrder2D(order):

    if type(order) is types.IntType:
        n1 = order
        n2 = order
    elif type(order) is types.TupleType or type(order) is numpy.ndarray:
        n1 = order[0]
        n2 = order[1]
    else:
        raise TypeError('Wrong argument!')

    return [n1,n2]


def numberOfTerms2D(order):

    nn = parseOrder2D(order)
    nn.sort()

    return int((nn[0]+1)*(2*nn[1]-nn[0]+2)/2.0)


def powerOrder2D(order):

    [n1,n2] = parseOrder2D(order)

    approximationOrder = max(n1,n2)

    order_1 = []
    order_2 = []

    for k in range(approximationOrder+1):
        for i1 in range(k,-1,-1):
            if i1>n1: continue
            if k-i1>n2: continue
            order_1.append(i1)
            order_2.append(k-i1)

    power_order = [order_1,order_2]

    return power_order


def polX2D(order,x,min_order):

    coeff = numpy.asarray(coeff)

    max_order = numpy.max(order)

    x1_p = numpy.vander(x[:,0],N=order[0]+1)[:,::-1]
    x2_p = numpy.vander(x[:,1],N=order[1]+1)[:,::-1]

    n = x.shape[0]
    y = numpy.zeros((n))

    n1 = order[0]
    n2 = order[1]

    X = []

    index = 0
    for k in range(min_order,max_order+1):
        for i1 in range(k,-1,-1):
            if i1>n1: continue
            if k-i1>n2: continue
            X.append(x1_p[:,i1]*x2_p[:,k-i1])

    return numpy.row_stack(X)


def poleval2D(coeff,order,x,min_order):

    order = numpy.asarray(order)
    coeff = numpy.asarray(coeff)

    max_order = numpy.max(order)

    x1_p = numpy.vander(x[:,0],N=order[0]+1)[:,::-1]
    x2_p = numpy.vander(x[:,1],N=order[1]+1)[:,::-1]

    n = x.shape[0]
    y = numpy.zeros((n))

    n1 = order[0]
    n2 = order[1]

    index = 0
    for k in range(min_order,max_order+1):
        for i1 in range(k,-1,-1):
            if i1>n1: continue
            if k-i1>n2: continue
            y = y + coeff[index]*x1_p[:,i1]*x2_p[:,k-i1]
            index += 1

    return y


def polder2D(coeff,order,m):
    """Return the (m[0],m[1])th derivative of the polynomial (coeff,order) of 2 variables
    """
    order = numpy.asarray(order)
    coeff = numpy.asarray(coeff)
    m =  numpy.asarray(m)

    powerIndexDict = powersIndexDict(order)
    nn = numberOfTerms2D(order)

    orderder = order-m
    powerIndexderDict = powersIndexDict(orderder)
    nnder = numberOfTerms2D(orderder)

    coeffder = numpy.ones((nnder))

    for powers,index in powerIndexderDict.iteritems():
        n1,n2=powers
        for i1 in range(m[0]): coeffder[index] *= n1+i1+1
        for i2 in range(m[1]): coeffder[index] *= n2+i2+1
        try:
            index0 = powerIndexDict[(n1+m[0],n2+m[1])]
            coeffder[index] *= coeff[index0]
        except KeyError:
            coeffder[index] *= 0.0

    return (coeffder,orderder)


def numberOfTerms1D(order):

    return int(order+1);


def powerOrder1D(order):

    power_order = [range(order+1)]

    return power_order


def poleval1D(coeff,order,x,min_order):

    coeff = numpy.asarray(coeff)
    order = numpy.asarray(order)

    x_p = numpy.vander(x[:,0],N=order[0]+1)[:,::-1]

    return numpy.dot(x_p[:,min_order:],coeff[min_order:])


def polX1D(order,x,min_order):

    x_p = numpy.vander(x[:,0],N=order[0]+1)[:,::-1]

    return xp[:,min_order:]



def polder1D(coeff,order,m=numpy.array((1))):
    """Return the (m[0])th derivative of the polynomial (coeff,order) of 1 variable
    """
    order = numpy.asarray(order)
    coeff = numpy.asarray(coeff)
    m =  numpy.asarray(m)

    powerIndexDict = powersIndexDict(order)
    nn = numberOfTerms1D(order)

    orderder = order-m
    powerIndexderDict = powersIndexDict(orderder)
    nnder = numberOfTerms1D(orderder)

    coeffder = numpy.ones((nnder))

    for powers,index in powerIndexderDict.iteritems():
        n1,=powers
        for i1 in range(m[0]): coeffder[index] *= n1+i1+1
        try:
            index0 = powerIndexDict[(n1+m[0])]
            coeffder[index] *= coeff[index0]
        except KeyError:
            coeffder[index] *= 0.0

    return (coeffder,orderder)



def testValue(coeff,orders,x):

    coeff = numpy.asarray(coeff)
    orders = numpy.asarray(orders)

    ndim = orders.shape[0]
    npts = x.shape[0]
    nterm = coeff.shape[0]

    power = numpy.asarray(powerOrder(orders,ndim))

    X = numpy.ones((npts,nterm))

    for i in range(npts):
        for j in range(nterm):
            for k in range(ndim):
                X[i,j] *= x[i,k]**power[k,j]

    y = numpy.dot(X,coeff)

    return y


def test(order,dim=3):

    import time

    n = numberOfTerms(order,dim=dim)
    p = powerOrder(order,dim=dim)
    print 'Order %s in %iD: %i terms' %(order,dim,n)
    for i in range(dim):
        print '  %s l=%i' %(p[i], len(p[i]))
    from random import random
    coeff = [random() for i in range(n)]
    orders = numpy.array([order]*dim)
    npts = 1000
    x = numpy.array([[random() for i in range(dim)] for j in range(npts)])

    t0 = time.time()
    y1 = poleval(coeff,orders,x)
    t1 = time.time()
    y2 = testValue(coeff,orders,x)
    t2 = time.time()

    print 'Delta t (with numpy)=%f    Delta t (without numpy)=%f' %(t1-t0,t2-t1)

    max_display = 10
    for i in range(min(max_display,npts)):
        print '    y=%e    y_=%e    delta=%e' %(y1[i],y2[i],y2[i]-y1[i])

    print

    print coeff
    m = [0]*dim
    m[0]=1
    (coeff,orders) = polder(coeff,orders,m,dim)
    print coeff


if __name__ == '__main__':
    test(2,dim=1)
    test(2,dim=2)
    test(2,dim=3)
    test(2,dim=4)
    #test((2,3,1,2),dim=4)

