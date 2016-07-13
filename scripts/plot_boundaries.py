#!/usr/bin/python

import sys,os.path,getopt
import numpy
import modl,util

def usage():

    print "Usage: %s file" % sys.argv[0]
    print "          file is an IMOD model file that contains points defining boundaries of a specimen"


def plotPoints( XYZ_1, XYZ_2,X,Y,Z1,Z2, scale=[1.0,1.0,10.0]):
    
    try:
        
        from enthought.mayavi import mlab
        
        XYZ_1[:,0] *= scale[0]
        XYZ_1[:,1] *= scale[1]
        XYZ_1[:,2] *= scale[2]

        XYZ_2[:,0] *= scale[0]
        XYZ_2[:,1] *= scale[1]
        XYZ_2[:,2] *= scale[2]
        
        fig = mlab.figure()

        c1 = (1.0,1.0,0.0)
        c2 = (1.0,1.0,0.0)
        
        s1 = mlab.points3d( XYZ_1[:,0], XYZ_1[:,1], XYZ_1[:,2], color=c1, scale_factor=20.0 )
        s2 = mlab.points3d( XYZ_2[:,0], XYZ_2[:,1], XYZ_2[:,2], color=c2, scale_factor=20.0 )
    
        Z1 *= scale[2]
        Z2 *= scale[2]

        mlab.surf( X, Y, Z1 )
        mlab.surf( X, Y, Z2 )

        mlab.show()

    except ImportError:

        print "Cannot complete 3d plotting:", ImportError
        

def main():

    try:
        filename = sys.argv[1]
        opts, args = getopt.getopt(sys.argv[2:],
            "h:",["help"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
        
    if filename==None:
        usage()
        raise SystemExit
    
    filename = os.path.abspath(filename)
    
    if not os.path.exists(filename):
        raise SystemExit, 'File %s does not exist!' %filename

    for option_,value in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"

    m = modl.Model(filename)
    m.loadFromFile()
    
    XYZ_1 = m.objects[0].points()
    try:
        XYZ_2 = m.objects[1].points()
    except:
        XYZ_2 = m.objects[0].points()
    
    p1 = util.polynomial_regression(XYZ_1[:,:2],XYZ_1[:,2:],[5,5])
    p2 = util.polynomial_regression(XYZ_2[:,:2],XYZ_2[:,2:],[5,5])
    
    [ X, Y ] = numpy.mgrid[ 0:910:50, 0:910:50 ]
    XY = numpy.column_stack((X.ravel(),Y.ravel()))
    
    shape = X.shape
    
    Z1 = p1.eval(XY).reshape(shape)
    Z2 = p2.eval(XY).reshape(shape)
    
    
    XYZ_1 = m.objects[0].contours[0].points
    try:
        XYZ_2 = m.objects[1].contours[0].points
    except:
        XYZ_2 = m.objects[0].contours[0].points 
    
    plotPoints( XYZ_1,XYZ_2,X,Y,Z1,Z2 )


main()
