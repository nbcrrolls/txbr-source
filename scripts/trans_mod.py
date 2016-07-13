#!/usr/bin/python

import sys
import os.path
import getopt

from txbr.utilities import *

RESET_FILE_OPTION = "reset-file"

OPTIONS = [ RESET_FILE_OPTION ]

def usage():
    
    print
    print "Usage: %s option -i input" % os.path.basename(sys.argv[0])
    print "          option in %s" %str(OPTIONS)
    print
    print "Option \"%s\": Zero out the transformation at the mid position." %RESET_FILE_OPTION
    
    
def reset_transformation_file( directory, file ):
    """Zero out the transformation at mid position."""

    file = os.path.join( directory, file )

    txbr.utilities.zero_out_transformations( file )
        

def main():
    
    flags1 = "hd:i:"
    flags2 = [ "help", "directory=", "input=" ]

    if len(sys.argv)==1 or sys.argv[1] not in OPTIONS:
        print
        print "You need to provide an option within %s for this command!" %OPTIONS
        usage()
        sys.exit(2)

    try:
        option = sys.argv[1]
        opts, args = getopt.getopt( sys.argv[2:], flags1, flags2 )
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    directory = '.'
    file = None

    for option_,value in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-d", "--directory"):
            directory = value
        elif option_ in ("-i", "--input"):
            file = value
        else:
            assert False, "unhandled option"

    if option==RESET_FILE_OPTION:
        reset_transformation_file( directory, file )

 
main()
