#!/usr/bin/python

import sys
import os.path
import getopt

import txbr.setup.txbound as tb


def usage():

    print "Usage: %s -d directory -b basename" % sys.argv[0]
    print "          -x x1[,...] -y y1[,...] -z z1[,...]"
    print "          --help --test --doPlot"


def main():
    '''Main routine for this pitch calculation module'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],'hd:b:x:y:z:',['help','doSave','doPlot'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    directory = '.'
    basename = None
    x, y, z = None, None, None
    doSave = False
    doPlot = False

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-d"):
            directory = value_
        elif option_ in ("-b"):
            basename = value_
        elif option_ in ("-x"):
            x = value_
        elif option_ in ("-y"):
            y = value_
        elif option_ in ("-z"):
            z = value_
        elif option_ in ("--doSave"):
            doSave = True
        elif option_ in ("--doPlot"):
            doPlot = True
        else:
            assert False, "unhandled option"

    if not os.path.exists(directory):
        print 'Directory %s does not exits!' %directory
        sys.exit()

    grid = None
    if x!=None or y!=None or z!=None:
        grid=[x,y,z]

    tb.eval_pitch(directory,basename,grid=grid,doPlot=doPlot,doSave=doSave)


main()
