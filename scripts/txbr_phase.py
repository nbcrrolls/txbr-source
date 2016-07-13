#!/usr/bin/python

import sys
import getopt
import txbr.phase as phase


def usage():

    print "Usage: %s -i input" % sys.argv[0]
    print "    -m model"
    print "        The model file that contains the corresponding points."
    print "    -o output"
    print "        Output MRC file containing the phase image."
    print "    --index i1,i2,ifocus"
    print "        Index of the three slices (under, over and in focus)"
    print "    --test-align"
    print "        Test the alignment of the three slices"
    print "    --help"
    print "        For usage help"


def main():

    '''Main routine for this stacking module'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],'hi:m:o:',['index=','eps=','test-align','test-diff','help'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    file = None
    indices = [0,1,2]
    model_file = None
    output = "phase-test.mrc"
    eps = 1e-4
    test_align = False
    test_diff= False

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-i"):
            file = value_
        elif option_ in ("-o"):
            output = value_
        elif option_ in ("-m"):
            model_file = value_
        elif option_ in ("--index"):
            indices = [ int(index) for index in value_.split(",") ]
        elif option_ in ("--eps"):
            eps = float(value_)
        elif option_ in ("--test-align"):
            test_align = True
        elif option_ in ("--test-diff"):
            test_diff = True
        else:
            assert False, "unhandled option"

    if file==None:
        usage()
        sys.exit(2)

    phase.getPhaseFromFile( file, index1=indices[0], index2=indices[1], indexOfFocus=indices[2], model_file=model_file, output=output, eps=eps, test_align=test_align, test_diff=test_diff )


main()
