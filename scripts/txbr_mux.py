#!/usr/bin/python

import sys
import os.path
import getopt
import txbr.txbrdao
import txbr.onthefly

def usage():

    print
    print 'Usage: %s.py -b basename [Options]' %os.path.basename(sys.argv[0])
    print
    print "    -d directory"
    print "        Directory of the TxBR project"
    print "    --wd work_directory"
    print "        Work directory of the TxBR project"
    print "    --bin n"
    print "        Bin the data in the x and y directions"
    print "    --restrict nseries"
    print "        Will only keep markers detected more than nseries series"
    print "    --dm value"
    print "        Max tolerated distance between corresponding markers."
    print "    --no-blur"
    print "        Do not blur the images before the initial correlation."
    print "    --reset"
    print "        Will recalculate the rotation angle between series."
    print "    -h or --help"
    print "        Help Information"

def main():
    '''Main routine for finding common markers between multiple tilt series.'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hb:d:",['bin=','wd=','restrict=','dm=','blur','reset','help'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    directory = "."
    work_directory = None
    basename = None
    ncommon = None
    dm = 2.5
    blur = True
    reset = False
    bin = 1

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-b"):
            basename = value_
        elif option_ in ("-d"):
            directory = value_
        elif option_ in ("--wd"):
            work_directory = value_
        elif option_ in ("--bin"):
            bin = int(value_)
        elif option_ in ("--restrict"):
            ncommon = int(value_)
        elif option_ in ("--dm"):
            dm = float(value_)
        elif option_ in ("--no-blur"):
            blur = False
        elif option_ in ("--reset"):
            reset = True
        else:
            assert False, "unhandled option"

    if basename==None:
        sys.exit("Please provide a basename!")

    if work_directory==None:
        work_directory = directory

    project = txbr.txbrdao.loadProject( directory=directory, basenames=basename, work_directory=work_directory, bin=bin, cfg_ext=".indiv" )

    txbr.onthefly.scanProjectSeries( project, ncommon=ncommon, dm=dm, blur=blur, reset=reset )


if __name__ == '__main__': main()