#!/usr/bin/python

import sys
import os.path
import getopt
import txbr.txbrdao
import txbr.onthefly
import txbr.utilities

from txbr.txbrconf import ALIGN_DIR
from txbr.txbrconf import log

def usage():

    print
    print 'Usage: %s.py -b basename [Options]' %os.path.basename(sys.argv[0])
    print
    print "    -d directory"
    print "        Directory of the TxBR project"
    print "    --wd work_directory"
    print "        Work directory of the TxBR project"
    print "    --from-bin src"
    print "        Bin value of the source"
    print "    --to-bin dest"
    print "        Bin value of the dest"
    print "    --size"
    print "        Window size use for fiducial corrections"
    print "    --test"
    print "        Some testing"


def main():

    '''Main routine for this rescaling'''

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hb:d:",['bin=','wd=','from-bin=','to-bin=','size=','help'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    directory = "."
    work_directory = None
    basename = None
    binsrc = 4
    bindest = 1
    size = 20

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
        elif option_ in ("--from-bin"):
            binsrc = int(value_)
        elif option_ in ("--to-bin"):
            bindest = int(value_)
        elif option_ in ("--size"):
            size = int(value_)
        else:
            assert False, "unhandled option"

    if basename==None:
        sys.exit("Please provide a basename!")

    if work_directory==None:
        work_directory = directory

    scale = float(bindest)/float(binsrc)

    basenames = txbr.utilities.extract_series_name_from_dir( basename, work_directory )
    base,extensions = txbr.utilities.extract_series_name(basenames)

    align_directory_src = os.path.join( work_directory, ALIGN_DIR, "bin%i" %binsrc )
    align_directory_dest = os.path.join( work_directory, ALIGN_DIR, "bin%i" %bindest )

    for ext in extensions:

        log.info("Extension: %s" %ext)

        bn = base + ext

        mrk_src = os.path.join( align_directory_src, bn + txbr.utilities.MARKER_EXTENSION )

        mrk_dest = os.path.join( align_directory_dest, bn + txbr.utilities.MARKER_EXTENSION )
        stack_dest = txbr.utilities.loadStack( directory, bn + txbr.utilities.ST_EXTENSION, work_directory, bin=bindest )

        txbr.onthefly.rescaleModel( mrk_src, mrk_dest, stack_dest, scale, size )


if __name__ == '__main__': main()