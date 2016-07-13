import os.path,logging,logging.config
import numpy,numpy.random

LOG_CONF = os.path.join(os.path.dirname(__file__),'..','txbr','log.conf')
logging.config.fileConfig(LOG_CONF)

DELTA_H = 1.0e-5

def f(X,a,b):

    return numpy.sum(X**3 + a*X)


def dfdX(X,a,b):

    return 3.0*X*X + a


def d2fdX2(X,a,b):

    shape = list(X.shape)
    shape.extend(X.shape)

    d2fdX2 = numpy.zeros(shape)

    for position1,value1 in numpy.ndenumerate(X):
        d2fdX2[position1][position1] = 6.0*value1

    return d2fdX2


def test_derivatives(f,dfdX,X,args=(),d2fdX2=None,delta_h=DELTA_H,loglevel=logging.DEBUG):
    '''Routine that test the derivatives and hessians for a function.
    '''

    #log = logging.getLogger('test')
    log = logging.getLogger()
    log.setLevel(loglevel)

    shape_der = X.shape

    shape_hess = list(shape_der)
    shape_hess.extend(shape_der)

    der = dfdX(X,*args)
    der.resize(shape_der)

    der_ = numpy.empty_like(der)

    for position,value in numpy.ndenumerate(X):
        Y = X.copy()
        Y[position] += delta_h
        der_[position] = (f(Y,*args)-f(X,*args))/delta_h

    if d2fdX2!=None:

        hess = d2fdX2(X,*args)
        hess.resize(shape_hess)
        hess_ = numpy.empty_like(hess)

        for position1,value1 in numpy.ndenumerate(X):
            for position2,value2 in numpy.ndenumerate(X):
                Y = X.copy()
                Y[position1] += delta_h
                hess_[position1] = numpy.resize((dfdX(Y,*args)-dfdX(X,*args))/delta_h,shape_der)

    # Logging statements

    diff1 = numpy.where(der==der_,0.0,numpy.abs(((der_-der)/der)))
    log.debug('First Derivative Checks')
    for position,value in numpy.ndenumerate(X):
        log.debug('i=%-5s    der=%- 10.3e  app.=%- 10.3e      Rel. Diff.=% 10.3e' %(position,der[position],der_[position],diff1[position]))

    if d2fdX2!=None:
        log.debug('Hessian Checks')
        diff2 = numpy.where(hess==hess_,0.0,numpy.abs(((hess-hess_)/hess)))
        for position1,value1 in numpy.ndenumerate(X):
            for position2,value2 in numpy.ndenumerate(X):
                log.debug('i=%-5s  j=%-5s    hess=%- 10.3e  app.=%- 10.3e      Rel. Diff.=%- 10.3e' %(position1,position2,hess[position1][position2],hess_[position1][position2],diff2[position1][position2]))

    log.info('First Derivative Error check:')
    log.info('Max Diff: %-10.2e     At: %s' %(numpy.max(diff1),numpy.argmax(diff1)))
    if d2fdX2!=None:
        log.info('Hessians Error check:')
        log.info('Max Diff: %-10.2e     At: %s' %(numpy.max(diff2),numpy.unravel_index(numpy.argmax(diff2),diff2.shape)))


def test():

    import numpy.random
    X = numpy.random.random((20,3))

    test_derivatives(f,dfdX,X,args=(2,0),d2fdX2=d2fdX2,loglevel=logging.DEBUG)


if __name__ == '__main__':

    logging.config.fileConfig(LOG_CONF)
    log = logging.getLogger('test')

    test()