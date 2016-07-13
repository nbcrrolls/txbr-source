#!/usr/bin/python

import sys
import os.path
import getopt
import txbr.stack

def usage():
    
    width_approximation_options = [txbr.stack.ISO_VOLUME_WIDTH,txbr.stack.MINIMUM_WIDTH,txbr.stack.MAXIMUM_WIDTH]
    
    print
    print 'Usage: %s.py -B boundary-model -V -inputvolume [Options]' %os.path.basename(sys.argv[0])
    print
    print "    -D direction (in [X,Y,Z], default: Z)"
    print "        Direction where to apply the flattening (normal to the slab-like specimen)"
    print "    --order orderx,ordery (default: 10,10 for polynomial mode; 2,2 for interpolation mode)"
    print "        Order of the polynomial map used to describe the surfaces"
    print "    --width-approximation w"
    print "        Hypothesis for the calculation the flattened specimen width (%s)" %",".join(width_approximation_options)
    print "    --interpolation"
    print "        Uses an interpolation scheme (rather than a polynomial regression) to evaluate the boundaries"
    print "    --doPlot"
    print "        Do some plotting"
    print "    --test"
    print "        Some testing"
    

def main():

    '''Main routine for this flattening module'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hB:D:V:",['interpolation','order=','width-approximation=','test','doPlot','help'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    vol_pth = None
    trk_path = None

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-B"):
            trk_path = value_
        elif option_ in ("-V"):
            vol_pth = value_
        elif option_ in ("-D"):
            direction = value_
        elif option_ in ("--order"):
            order = [int(s) for s in value_.split(',')]
        elif option_ in ("--interpolation"):
            mode = 'interpolation'
        elif option_ in ("--width-approximation"):
            width_type = value_
        elif option_ in ("--test"):
            test = True
        elif option_ in ("--doPlot"):
            doPlot = True
        else:
            assert False, "unhandled option"

    if vol_pth==None:
        sys.exit("No volume to flatten!")

    if not os.path.exists(vol_pth):
        sys.exit('The file %s does not exists!' %vol_pth)

    if trk_path==None:
        sys.exit("No volume to flatten!")

    if not os.path.exists(trk_path):
        sys.exit('The file %s does not exists!' %trk_path)

    keywords = {}

    try:
        keywords['direction'] = direction
    except NameError:
        pass

    try:
        keywords['mode'] = mode
    except NameError:
        pass

    try:
        keywords['test'] = test
    except NameError:
        pass

    try:
        keywords['doPlot'] = doPlot
    except NameError:
        pass

    try:
        keywords['order'] = order
    except NameError:
        pass

    try:
        keywords['width_app'] = width_type
    except NameError:
        pass

    txbr.stack.flatten(vol_pth,trk_path,**keywords)


main()
