import sympy

def nterms_(order1,order2,order3):

    approximationOrder = max(order1,order2,order3)

    nn = 0
    for k in range(approximationOrder+1):
        for j in range(k,-1,-1):
            for i in range(k-j,-1,-1):
                if i>order1: continue
                if j>order2: continue
                if k-i-j>order3: continue
                nn += 1

    return '(n1,n2,n3)=(%i,%i,%i) -> %i terms' %(order1,order2,order3,nn)


def nterms(order):

    n = sympy.Symbol('n')
    i,j,k = sympy.symbols('ijk')

    f1 = sympy.sum(1,(i,0,k-j))
    f2 = sympy.sum(f1,(j,0,k))
    f = sympy.sum(f2,(k,0,n))

    return 'Expression: %s: n=%i -> %i terms' %(str(f),order,f.subs(n,order))

print nterms_(0,0,0)
print nterms_(1,1,1)
print nterms_(2,2,2)
print nterms_(3,3,3)

print nterms(0)
print nterms(1)
print nterms(2)
print nterms(3)