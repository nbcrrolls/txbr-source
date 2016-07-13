#!/usr/bin/python

import os.path
import subprocess
import sys
import getopt

def usage():
    """Usage for the command *txbr_quick.py*."""

    print 'Usage: %s.py -b basename[,...]  [Options]' % os.path.basename(sys.argv[0])
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files (.fid and .txbr)"
    print "    --wd work_directory (default to directory)"
    print "        Name of the directory where the output files are stored"
    print "    -h or --help"
    print "        Help Information"

def main():
    """The main routine for this iterative process"""

    try:

        flags1 = "hb:d:l:n:"
        flags2 = [ "help", "directory=", "wd=", "work_directory=", "basename=" ]

        opts, args = getopt.getopt( sys.argv[1:], flags1, flags2 )

    except getopt.GetoptError, err:

        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    keywords = {}

    keywords['directory'] = "rawdata"
    keywords['work_directory'] = "."
    keywords['basename'] = None
    keywords['bin'] = 4
    keywords['size'] = 8

    for option,value in opts:
        if option in ('-h','--help'):
            usage()
            sys.exit()
        elif option in ('-d','--directory'):
            keywords['directory'] = value
        elif option in ('--wd'):
            keywords['work_directory'] = value
        elif option in ('-b','--basename'):
            keywords['basename'] = value
        elif option in ('--bin'):
            keywords['bin'] = value
        else:
            assert False, "unhandled option"

    if keywords['basename']==None:
        usage()
        sys.exit()

    nproc = 10

    ALIGN_CMD = "txbr_auto.py -b {basename} -d {directory} --wd {work_directory} --size {size} --bin {bin} --rough".format(**keywords)
    FILTER_CMD = 'txbr_filter.py -b {basename} -d {directory} --wd {work_directory} --bin {bin} --doClean'.format(**keywords)
    BCKPRJ_CMD = 'txbr_bckprj.py -b {basename} -d {directory} --wd {work_directory} --bin {bin} --mpi'.format(**keywords)

    subprocess.call(ALIGN_CMD,shell=True)
    subprocess.call("mpiexec -n {nproc} {command}".format(nproc=nproc,command=FILTER_CMD),shell=True)
    subprocess.call("mpiexec -n {nproc} {command}".format(nproc=nproc,command=BCKPRJ_CMD),shell=True)
    

if __name__ == '__main__': main()