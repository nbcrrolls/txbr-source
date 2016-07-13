#!/usr/bin/python

import sys
import getopt
import mrc


def view(filename,**kw):

    print kw

    mrc.view(filename,**kw)


def usage():

    print "Usage: %s filename -X x -Y y -Z z --help --test --doPlot" % sys.argv[0]


def main():

    '''Main routine for this unwarping module'''

    try:
        opts, args = getopt.getopt(sys.argv[2:],"hX:Y:Z:",["help"])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    filename = sys.argv[1]

    print opts

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-X"):
            X = int(value_)
        elif option_ in ("-Y"):
            Y = int(value_)
        elif option_ in ("-Z"):
            Z = int(value_)
        else:
            assert False, "unhandled option"

    keywords = {}

    try:
        keywords['x'] = X
    except NameError:
        pass

    try:
        keywords['y'] = Y
    except NameError:
        pass

    try:
        keywords['z'] = Z
    except NameError:
        pass

    view(filename,**keywords)


main()
