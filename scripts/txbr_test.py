#!/usr/bin/python

import sys
import os.path
import time
import getopt
import numpy
import scipy
import pylab
import logging
import logging.config
import modl
import util
import txbr
import txbr.onthefly

from txbr.align import contalign,residual
from txbr import LOG_CONF
from txbr.txbrdao import loadProject

logging.config.fileConfig(LOG_CONF)

log = logging.getLogger('align')

def test_interp():

    import scipy.interpolate

    pi = numpy.pi

    x = numpy.arange(0,2*pi+pi/4,2*pi/8)
    y = numpy.sin(x)
    tck = scipy.interpolate.splrep(x,y,s=0)
    xnew = numpy.arange(-pi/4,2*pi,pi/50)
    ynew = scipy.interpolate.splev(xnew,tck,der=0)

    pylab.figure()
    pylab.scatter(x,y)
    pylab.plot(xnew,ynew)



def test_poly():
    """Test the polynomial class
    """

    import swiginac

    x,y = swiginac.symbol('x'),swiginac.symbol('y')

    variables = [x,y]

    f1 = -1 + x + 2*y + 3*x**2 - 4*x*y + 5*y**2 - 5/2*y**4
    f2 = 4*x**2 - 5*y - x*y

    f = swiginac.integral(y,0.0,1.0,f1).eval_integ()

    print 'f1 = %s' %f1
    print 'f2 = %s' %f2
    print 'f = %s' %f1

    p_ginac = util.asPolyM(f,variables,mode='ginac')
    p_str = util.asPolyM(f,variables,mode='string')

    p2 = util.asPolyM(f2,variables,mode='string')

    print 'p_ginac = %s' %p_ginac
    print 'p2 = %s' %p2
    print 'p_ginac + p2 = %s' %(p_ginac + p2)

    p2_asGiNaC = p2.asGinacFunction()
    p2_c_p2 = p2.compose({'x':p2})
    print 'p2_asGiNaC = %s' %str(p2_asGiNaC)
    print 'p2 comp p2 = %s' %p2_c_p2
    return

    X0 = numpy.random.random_sample((1,2))
    print p_ginac.eval(X0)
    print p_str.eval(X0)
    print f.subs([x==X0[0,0],y==X0[0,1]])

    stTest = True

    if stTest:
        p = p_str
    else:
        p = p_ginac

    print 'p Coefficients'
    print p.coefficients
    print 'p Powers'
    print p.powers

    return

    pdiffx = p.diff(x)
    pdiffy = p.diff(y)

    print 'p=%s' %p
    print 'dpdx=%s' %pdiffx
    print 'dpdy=%s' %pdiffy

    print 'dpdx Coefficients'
    print pdiffx.coefficients
    print 'dpdx Powers'
    print pdiffx.powers

    print 'dpdy Coefficients'
    print pdiffy.coefficients
    print 'dpdy Powers'
    print pdiffy.powers


def test_mapping():

    directory = os.path.join(os.environ['TXBR'],'examples','dataset','cb018')
    basename = 'cb018a'

    directory = os.path.join(os.environ['TXBR'],'examples','dataset','gia')
    basename = 'gia'

    project = loadProject(directory=directory,basenames=basename)


    print 'nx=%i  ny=%i' %(project.nx,project.ny)

    serie = project.series[0]
    trackModel = modl.Model(serie.nameOfFile('track'))
    trackModel.loadFromFile()

    n3 = 5
    objectMapping = trackModel.polynomialMapping(n3,u_t=[0.0,1.0,0.0])

#
#    pylab.xlim(0,project.nx)
#    pylab.ylim(0,project.ny)


def test_residual():
    """Test expression of the residual and its derivatives
    """

    log.info('Test Residual')

    tic = time.time()

    resid = residual.Residual(ntilt=1,npatch=1,n1=2,n2=1,n3=4,n4=2)
    resid.test(derivative=False,hessian=False)

    toc = time.time()

    log.info('Time %f' %(toc-tic))


def test_model():

    global directory
    global basename

    filename = os.path.join(directory,basename + '.trk')
    trackModel = modl.Model(filename)
    trackModel.loadFromFile()

    for object in trackModel.objects:
        print  'Object #%i   Z=%s' %(object.indexOfObject,object.zExtensionSet())

    for object in trackModel.objects:
        for contour in object.contours:
            print 'Object #%-5i   Contour #%-5i   Z=%s' %(object.indexOfObject,contour.indexOfContour,contour.zExtensionSet())


def test_approximation_mapping():
    """Test approximation mapping for the structure
    """

    global project

    contalign.TxBRContourAlign(project,prepare=False)


def test_occluded_contour(objects=None):

    global project

    print 'Occluded Contour'

    alignment = contalign.TxBRContourAlign(project,prepare=True)
    alignment.initializeStructureSet(objects=objects)

    project.save()


def test_alignment():

    global project

    alignment = contalign.TxBRContourAlign(project)
    alignment.process()

    project.save()


def approximate_contour(plot_contours=True,order=2):

    global directory
    global basename

    filename = os.path.join(directory,basename) + '.trk'

    trackModel = modl.Model(filename)
    trackModel.loadFromFile()

    data = []
    coeffs = numpy.zeros((0,2,order+1),dtype=float)

    object = trackModel.getObject(0)

    ncontour = object.numberOfContours()

    for icontour in range(ncontour):

        contour = trackModel.getObject(0).getContour(icontour)

        n = contour.npoints()

        data.append(contour.points)

        coeffs_ = contour.polynomialApproximationCoefficients(order)
        coeffs_.resize((1,2,order+1))

        coeffs = numpy.vstack((coeffs,coeffs_))


    if plot_contours:

        for icontour in range(ncontour):

            tdata = scipy.linspace(0,1,num=n)

            xfit = scipy.polyval(coeffs[icontour,0,::-1], tdata)
            yfit = scipy.polyval(coeffs[icontour,1,::-1], tdata)

            data_ = data[icontour]
            pylab.plot(data_[:,0],data_[:,1],'k.')
            pylab.plot(xfit,yfit, 'r-')

        pylab.show()

    print coeffs


def test():

    directory = os.path.join(os.environ['TXBR'],'examples','dataset','cb018')
    basename = 'cb018a'

    filename = os.path.join(directory,basename) + '.seg_3obj'

    trackModel = modl.Model(filename)
    trackModel.loadFromFile()

    numberOfObjects = trackModel.numberOfObjects()

    for iobject in range(numberOfObjects):

        object = trackModel.getObject(iobject)
        numberOfContours = object.numberOfContours()

        for icontour in range(numberOfContours):

            contour = object.getContour(icontour)

            print contour.points
            print contour.polynomialApproximationCoefficients(2)

def test_warp( vol_pth, trk_path, direction='Z' ):

    import txbr.utilities.warp

    txbr.utilities.warp.main(vol_pth, trk_path)


def test_detect( size=10 ):
    """Test the bead detection
    """

    log.info('Test bead detection')

    global directory
    global basename

    options = { "directory":directory, "fin":"%s.st" %basename, "size":size, "nmax":1500, "slices":[60]  }

    fin = os.path.join( options["directory"], options["fin"])

    tic = time.time()

    txbr.onthefly.findMarkersInStack( fin, nmax=options["nmax"], slices=options["slices"], size=options["size"] )

    toc = time.time()

    log.info('Time %f' %(toc-tic))


def usage():

    print "Usage: %s -b basename -d directory" % sys.argv[0]


def main():
    '''Main routine for this test module'''

    try:
        option = sys.argv[1]
        opts, args = getopt.getopt(sys.argv[2:],
            "hd:w:b:XYZ",["help","directory=","base=","size=","objects="])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    global directory
    global basename

    directory = "."
    basename = None
    objects = None
    size = None

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-d","--directory"):
            directory = value_
        elif option_ in ("-b","--basename"):
            basename = value_
        elif option_ in ("--size"):
            size = int(value_)
        elif option_ in ("-X","-Y","-Z"):
            direction = option_[1]
        elif option_ in ("--objects"):
            values = value_.split(',')
            objects = []
            for iobject in values:
                try:
                    objects.append(int(iobject))
                except ValueError:
                    pass
        else:
            assert False, "unhandled option"

    global project

    if directory!=None and basename!=None:
        project = loadProject(directory=directory,basenames=basename)


    if option=='align':

        test_alignment()

    elif option=='ginac':

        test_poly()

    elif option=='residual':

        test_residual()

    elif option=='structure':

        test_occluded_contour(objects=objects)

    elif option=='warp':

        trk_path = '/Users/sph/Electron_Tomography/txbr/dataset/g01/g01/warp/g01.surfacemod'
        vol_pth = '/Users/sph/Electron_Tomography/txbr/dataset/g01/g01/warp/g01_z_-123.0.out1'

        keywords = {}

        try:
            keywords['direction'] = direction
        except NameError:
            pass

        test_warp(vol_pth,trk_path,**keywords)

    elif option=='detect':

        keywords = {}
        if size!=None: keywords["size"]=size

        test_detect( **keywords )

main()


try:
    print 'End of TxBR align'
    pylab.show()
except:
    print sys.exc_info()


