import sympy

def extractCoefficient(P,x,order):
    der = P
    for i in range(1,order+1):
        der = sympy.diff(der,x)/i
    return der.subs(x,0)

def extractCoefficients(P,x):
    der,coeffs = P,[]
    i = 0
    while der!=0:
        i += 1
        coeffs.append(der.subs(x,0))
        der = sympy.diff(der,x)/i
    return coeffs

def extractAll(P):
    P = sympy.Basic.as_poly(P)
    return (P.symbols,P.coeffs,P.monoms)

def test():

    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    P = 2.7*x**3*y**2 + 2.3*x**2 + 1
    print 'P = %s' %P
    print 'Coefficients = %s' %extractCoefficients(P,x)
    print 'Coefficients at order %i = %s' %(2,extractCoefficient(P,x,2))
    print '[ Symbols=%s, Coeffs=%s, Monoms=%s ]' %extractAll(P)

if __name__ == '__main__':
    test()