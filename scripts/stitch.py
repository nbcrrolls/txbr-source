#!/home/sphan/usr/local/bin/python

"""Wrapper script to imod stitching routines"""
import os.path

import os
import getopt
import sys
import numpy

def usage():
    """Usage for the command *stitch.py*."""

    print 'Usage: %s.py -b basename[,...]' % os.path.basename(sys.argv[0])
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input stack (.st) to stitch"
    print "    --wd work_directory (default to directory)"
    print "        Name of the directory where the output files are stored"
    print "    --bin n"
    print "        Bin the data in the x and y directions"


def __execute__( command ):

    print command
    os.system(command)
    

def stitch( basename, directory, work_directory, bin=1 ):
    """Stitch the montage tilt series"""

    doNormalize = True

    original_stack = "%s.st" %basename
    original_stack = os.path.join( directory, original_stack )

    normalized_series = "%s_norm.st" %basename
    piece_list_file = "%s.pl" %basename
    output_stack = "%s.bl" %basename
    normalized_series = os.path.join( work_directory, normalized_series )
    piece_list_file = os.path.join( work_directory, piece_list_file )
    output_stack = os.path.join( work_directory, output_stack )

    if bin!=1:
        command = 'newstack -bin %i %s %s' %(bin,original_stack,normalized_series)
        __execute__( command )
        
    if not os.path.lexists(): os.symlink(original_stack,normalized_series)

    if doNormalize:
        command = 'newstack -float 3 %s %s' %(normalized_series,normalized_series)
        __execute__( command )

    command = 'extractpieces %s %s' %(normalized_series,piece_list_file)
    __execute__( command )

    if bin!=1:
        pieces = numpy.loadtxt(piece_list_file)
        pieces[:,:2] /= bin
        numpy.savetxt(piece_list_file,pieces)
        
    command = 'blendmont -ImageInputFile %s -PieceListInput %s -ImageOutputFile %s -RootNameForEdges %s' %(normalized_series,piece_list_file,output_stack,basename)
    __execute__( command )


def main():
    """The main routine to stitch"""

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hb:d:",['wd=','bin=','help'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    basename = None
    directory = "."
    work_directory = "."
    bin = 1

    for option,value in opts:
        if option in ('-h','--help'):
            usage()
            sys.exit()
        elif option in ('-b'):
            basename = value
        elif option in ('-d'):
            directory = value
        elif option in ('--wd'):
            work_directory = value
        elif option in ('--bin'):
            bin = int(value)

    directory = os.path.abspath(directory)
    work_directory = os.path.abspath(work_directory)
            
    stitch( basename, directory, work_directory, bin=bin )


if __name__ == '__main__':
    main()