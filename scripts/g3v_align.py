#!/usr/bin/python

import sys
import getopt
import txbr.g3view


def usage():

    print "Usage: %s -f filename" % sys.argv[0]
    print "    -d directory"
    print "        Directory containing the stack to align"
    print "    --peak-threshold value [default=15.0]"
    print "        Threshold criteria delimiting good/bad image sequences. A value of 0 allows"
    print "        keeping all the slices of the volume to align."
    print "    --no-stack"
    print "        Only the correlation information is generated"
    print "    --help"
    print "        For usage help"


def main():

    '''Main routine for this stacking module'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],'hf:d:',['help','peak-threshold=','no-stack'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    filename = None

    keywords = {}
    keywords['directory'] = "."
    keywords['Lx'] = None
    keywords['Ly'] = None

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-f"):
            filename = value_
        elif option_ in ("-d"):
            keywords['directory'] = value_
        elif option_ in ("--no-stack"):
            keywords['createAlignedStack'] = False
            keywords['createBadImageStack'] = False
        elif option_ in ("--peak-threshold"):
            try:
                keywords['corr_ratio_threshold'] = float(value_)
            except:
                pass
        else:
            assert False, "unhandled option"

    txbr.g3view.alignStack( filename, **keywords )

main()
