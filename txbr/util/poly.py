import re
import numpy
import swiginac

class PolyM:
    """PolyM is a class that represents multivariate polynoms.
    It is built on * numpy* and *swiginac*, the interface to GiNac for its symbolic
    computation capabilities (composition of polynoms, and so on...)
    """

    def __init__( self, variables, coefficients, powers ):
        '''Initialize the polynom function:
        Variable variables is a list of variables (length n)
        Variable coefficients is a numpy array of shape (nn)
        Variable powers is a numpy array of integer of shape (nn,n)
        '''

        # Initialize the row data

        self.variables,self.coefficients,self.powers,self.orders = \
                        self.initialize(variables,coefficients,powers)

#        # Create the GiNac counterpart
#
#        self.dict_var_ginac,self.f_ginac = self.asGinacFunction()


    def initialize(self,variables,coefficients,powers):
        '''This routine cleans up the arrays defining the polynom. It collects
        same monoms together, and reorder the different terms.'''

        variables = numpy.asarray(variables,dtype='string')
        coefficients = numpy.asarray(coefficients,dtype='float')
        powers = numpy.asarray(powers,dtype='int')
        orders = numpy.zeros((variables.size))

        if len(coefficients)==0: return (variables,coefficients,powers,orders)

        # Reordering the terms with the cumulative powers

        acc_pow = numpy.cumsum(powers[:,::-1],axis=1)
        indices = numpy.lexsort(acc_pow.T)

        coefficients = coefficients[indices]
        powers = powers[indices]

        orders = numpy.amax(powers,axis=0)

        # Collect similar monoms together

        i=0
        while i<coefficients.size-1:
            if numpy.all(powers[i,:]==powers[i+1,:]):
                coefficients[i] += coefficients[i+1]
                coefficients = numpy.delete(coefficients,[i+1],axis=0)
                powers = numpy.delete(powers,[i+1],axis=0)
            i+=1

        return (variables,coefficients,powers,orders)


    def extract_coefficients(self,powers):
        '''Extract polynomial coefficients corresponding to the powers variable'''

        coefficients = []

        for pow in powers:
            indices = numpy.where(numpy.all(numpy.equal(pow,self.powers),axis=1)==True)
            if len(indices)!=1:
                raise ValueError
            index = indices[0]
            if index.size==0:
                coefficients.append(0.0)
            else:
                coefficients.append(self.coefficients[index[0]])

        return numpy.array(coefficients)


    def asGinacFunction(self):
        '''Return a GiNac function for this polynom.'''

        dict_var_ginac = dict([(var,swiginac.symbol(var)) for var in self.variables])

        nn = len(self.coefficients)

        if nn==0:
            return dict_var_ginac,0

        # GiNaC variables must be in the right order
        var_ = numpy.array([dict_var_ginac[var] for var in self.variables])

        X = numpy.tile(var_,(nn,1))
        X = numpy.prod(X**self.powers,axis=1)

        f = numpy.sum(numpy.dot(self.coefficients,X))
        f = f.expand()

        return dict_var_ginac,f


    def info(self):
        '''Print information on this polynom.'''

        print 'orders: %s' %self.orders
        print 'variables: %s' %self.variables
        print 'coefficients: %s' %self.coefficients
        print 'powers: %s' %self.powers.T[0]
        for i in range(1,len(self.powers.T)): print '        %s' %self.powers.T[i]
#        print 'GiNaC: %s' %self.f_ginac


    def diff(self,var):
        '''Returns the derivative of the polynom relativaly to parameter var.
        '''

        var = str(var)

        index = numpy.argwhere(self.variables==var)
        index = index[0,0]

        variables_diff = list(self.variables)
        coeffs_diff = self.coefficients.copy()
        powers_diff = self.powers.copy()

        coeffs_diff[:] = coeffs_diff[:]*powers_diff[:,index]
        powers_diff[:,index] = powers_diff[:,index]-1

        indices = numpy.where(powers_diff[:,index]>=0)

        if indices[0].size==0:
            coeffs_diff = numpy.array([])
            powers_diff = numpy.zeros((0,len(variables_diff)))
        else:
            coeffs_diff = coeffs_diff[indices]
            powers_diff = powers_diff[indices]

        return PolyM(variables_diff,coeffs_diff,powers_diff)


    def eval(self,x):
        '''Evaluate the polynom for some values x.
        Variable x: a numpy array of shape (n,nvar) where nvar is the number of
        variables involved in this polynomial function.'''
        
        x = numpy.asarray(x)

        s0 = x.shape

        nvar = len(self.variables)
        n = x.size/nvar
        s = (n,nvar)
        
        x = numpy.resize(x,s)
        
        ncoeffs = self.coefficients.size

        value = numpy.zeros((n,ncoeffs))

        for icoeff in range(ncoeffs):
            value[:,icoeff] = self.coefficients[icoeff]
            for ivar in range(nvar):
                value[:,icoeff] *= x[:,ivar]**self.powers[icoeff,ivar]

        value = numpy.sum(value,axis=1)
        
        value.resize((n))

        return value


    def compose(self,polynoms):
        '''Method to compose the polymom with other multivariate polynoms.
        Variable polynoms a dictionary that contains the (string) variable to replace
        as keys, and multivariate polynoms as values.
        This routine implemented using GiNac functionality.'''

        variables = []
#        g = self.f_ginac

        dict_var_ginac,f_ginac = self.asGinacFunction()
        g = f_ginac

        for var,p in polynoms.iteritems():
#            variables.extend(p.dict_var_ginac.values())  # keep track of the variables
#            x = self.dict_var_ginac[var]
            dict_var_ginac_,f_ginac_ = p.asGinacFunction()
            variables.extend(dict_var_ginac_.values())  # keep track of the variables
            x = dict_var_ginac[var]
            f = f_ginac_
            g = g.subs(x==f)

        g = g.expand()

        dict = {}
        for v in variables: # Create a new set of consistent swiginac variables
            dict[str(v)] = swiginac.symbol(str(v))

        for v in variables: # Replace old variables with new variables
            g = g.subs(v==dict[str(v)])

        variables = sorted(dict.keys())
        variables = [ dict[v] for v in variables ]

        return asPolyM(g,variables)


    def __add__(self,other):
        '''Add two polynoms together'''

        if numpy.any(self.variables!=other.variables):
            raise ValueError('Not the same list of variables')

        coeffs = numpy.concatenate((self.coefficients,other.coefficients))
        powers = numpy.concatenate((self.powers,other.powers))

        return PolyM(self.variables,coeffs,powers)


    def __repr__(self):

        dict_var_ginac,f_ginac = self.asGinacFunction()
        
#        return "PolyM: %s" %str(self.f_ginac)
        return "PolyM: %s" %str(f_ginac)




def asPolyM(poly,vars,mode='ginac'):
    """Return (poly,vars) as a multivariate polynom PolyM.
    """

    if mode=='ginac':
        coeffs,powers = extract_swiginac_coeff(poly,vars)
    elif mode=='string':
        coeffs,powers = extract_string_coeff(poly,vars)

    coeffs = numpy.asarray(coeffs,dtype='float')
    powers = numpy.asarray(powers,dtype='int')

    indices = numpy.where(coeffs!=0)

    variables = [str(var) for var in vars]
    coefficients = coeffs[indices]
    powers = powers[indices]

    return PolyM(variables,coefficients,powers)


def extract_swiginac_coeff(poly,variables):
    """extract_swiginac_coeff extracts the coefficients and corresponding powers from a multivariate
    swiginac polynoms given a set of variables
    """

    coefficients = [poly]
    powers = [[]]

    for ivar in range(len(variables)):
        var = variables[ivar]
        for icoeff in range(len(coefficients)-1,-1,-1):
            c = coefficients[icoeff]
            degrees = range(c.degree(var)+1)
            P = [c.coeff(var,j) for j in degrees]
            D = [[p for p in powers[icoeff]] + [j] for j in degrees]
            coefficients[icoeff:icoeff+1] = P
            powers[icoeff:icoeff+1] = D


    return (coefficients,powers)


def extract_string_coeff(poly,vars):
    """extract_string_coeff extracts the coefficients and corresponding powers from a multivariate
    string polynom expression given a set of variables
    """

#    print poly
#    print

    variables = [str(var) for var in vars ]

    # hack for the exponential numbers

    p = re.compile('([^E]|\A)([+-])')
    pstar = re.compile('\\*')

    tokens_ = p.split(str(poly))

    tokens = []
    for tk in tokens_:
        if len(tokens)==0:
            tokens.append(tk)
            continue
        lst = tokens[len(tokens)-1]
        if lst=='+' or lst=='-' or tk=='+' or tk=='-':
            tokens.append(tk)
            continue
        tokens[len(tokens)-1] +=tk


    if len(tokens)%2==1 and len(tokens[0])==0: tokens.pop(0)
    if len(tokens)%2==1 and len(tokens[0])!=0: tokens.insert(0,'+')

    n = len(tokens)/2
    nvar = len(variables)

    coefficients = numpy.zeros((n))
    powers = numpy.zeros((n,nvar))

    for i in range(n):
        sign = tokens[2*i]
        token = tokens[2*i+1]
        units = pstar.split(token)
        coeff = eval(sign + '1.0')
        while (len(units)>0):
            try:
                unit = units.pop(0)
                index = variables.index(unit)
                power = 1
                if len(units)>=2 and len(units[0])==0:
                    units.pop(0)
                    power = units.pop(0)
                powers[i,index] = int(power)
            except ValueError: # unit is not a variable. This is the coefficient
                coeff = eval(sign + '1.0*' + unit)
            except IndexError:
                pass
        coefficients[i] = coeff


    return (coefficients,powers)


