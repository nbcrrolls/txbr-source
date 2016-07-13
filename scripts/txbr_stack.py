#!/usr/bin/python

import sys
import getopt
import txbr.stack as stack


def usage():

    print "Usage: %s -B boundary-model -V volume" % sys.argv[0]
    print "    --epsilon eps"
    print "        Epsilon condition"
    print "    --doPlot"
    print "        Do some plotting"
    print "    --test"
    print "        Some testing"
    print "    --help"
    print "        For usage help"


def main():

    '''Main routine for this stacking module'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],'hB:V:',['interface=','help','test','doPlot'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    vol_pth = None
    trk_path = None
    interface = None
    test = False
    doPlot = False

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-B"):
            trk_path = value_
        elif option_ in ("-V"):
            vol_pth = value_
        elif option_ in ("--interface"):
            interface = int(value_)
        elif option_ in ("--test"):
            test = True
        elif option_ in ("--doPlot"):
            doPlot = True
        else:
            assert False, "unhandled option"

    stack.dewarp( vol_pth, trk_path, interface=interface, doPlot=doPlot, test=test)

    if doPlot:
        import pylab
        pylab.show()

main()
