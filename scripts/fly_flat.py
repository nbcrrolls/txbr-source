#!/usr/bin/python

import sys
import os.path
import getopt
import numpy
import modl
import txbr.setup

from enthought.mayavi import mlab


def usage():

    print
    print 'Usage: %s.py -b basename [Options]' %os.path.basename(sys.argv[0])
    print
    print "    -d directory"
    print "        Directory of the TxBR project"
    print "    --order"
    print "        Order of the polynomial map used to describe the surfaces"
    print "    --scale sx,sy,sy"
    print "        Sclaing factor in x, y and z"
    print "    --doPlot"
    print "        Do some plotting"
    print "    --test"
    print "        Some testing"


def test_fly_flattening( file, flatten_order, scale ):

    scale = numpy.asarray(scale)

    m = modl.Model(file)
    m.loadFromFile()

    XYZ = m.points()

    ss = numpy.max(XYZ,axis=0)
    ss = numpy.ones((3))
    XYZ = XYZ/ss

    print ss

    center = ( numpy.max(XYZ,axis=0) + numpy.min(XYZ,axis=0) )/2.0

    print XYZ.shape
    print center

    warp_to_flat_tf, flat_to_warp_tf = txbr.setup.evalFlatteningTransformation( XYZ, flatten_order, center, doPlot=False, test=True )

    XYZ_flat = warp_to_flat_tf.forward(XYZ)
    t = (numpy.min(XYZ_flat,axis=0)[2]+numpy.max(XYZ_flat,axis=0)[2])/2.0

    index1 = numpy.where(XYZ_flat[:,2]>=t)
    index2 = numpy.where(XYZ_flat[:,2]<t)

    # Rescale the marker locations

    XYZ = XYZ*scale
    XYZ_flat = XYZ_flat*scale

    # Partition the markers on the two different surfacess
    
    XYZ1 = XYZ[index1]
    XYZ2 = XYZ[index2]

    XYZ1_flat = XYZ_flat[index1]
    XYZ2_flat = XYZ_flat[index2]

    # Plot

    XYZ1 = XYZ1*ss
    XYZ2 = XYZ2*ss
    XYZ1_flat = XYZ1_flat*ss
    XYZ2_flat = XYZ2_flat*ss
    
    s = max(XYZ.shape)/50.0
    color1 = (1,0,0)
    color2 = (0,1,0)

    mlab.figure()

    mlab.points3d( numpy.ravel(XYZ1[:,0]), numpy.ravel(XYZ1[:,1]), numpy.ravel(XYZ1[:,2]), color=color1, scale_factor=s )
    mlab.points3d( numpy.ravel(XYZ2[:,0]), numpy.ravel(XYZ2[:,1]), numpy.ravel(XYZ2[:,2]), color=color2, scale_factor=s )

    mlab.show()

    mlab.figure()

    mlab.points3d( numpy.ravel(XYZ1_flat[:,0]), numpy.ravel(XYZ1_flat[:,1]), numpy.ravel(XYZ1_flat[:,2]), color=color1, scale_factor=s )
    mlab.points3d( numpy.ravel(XYZ2_flat[:,0]), numpy.ravel(XYZ2_flat[:,1]), numpy.ravel(XYZ2_flat[:,2]), color=color2, scale_factor=s )

    mlab.show()


def main():

    '''Main routine for this flattening module'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hb:d:",['order=','scale=','test','doPlot','help'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    directory = "."
    basename = None
    flatten_order = 4
    scale = [ 1.0, 1.0, 1.0 ]
    test = False
    doPlot = False

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-b"):
            basename = value_
        elif option_ in ("-d"):
            directory = value_
        elif option_ in ("--order"):
            flatten_order = int(value_)
        elif option_ in ("--scale"):
            scale = [ float(s) for s in value_.split(",") ]
        elif option_ in ("--test"):
            test = True
        elif option_ in ("--doPlot"):
            doPlot = True
        else:
            assert False, "unhandled option"

    if basename==None:
        sys.exit("Please provide a basename!")

    file = os.path.join("txbr-align","%s.mod"%basename)
    file = os.path.join(directory,file)

    if not os.path.exists(file):
        sys.exit("File %s does not exist!" %file)

    test_fly_flattening( file, flatten_order, scale)



main()
