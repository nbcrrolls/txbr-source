#!/usr/bin/python

import os
import sys
import getopt
import numpy
import txbr.utilities

def usage():
    
    program = os.path.basename(os.path.basename(sys.argv[0]))

    print
    print "Usage: %s -b basename1[,...] [Options]" %(program)
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files (.fid and .txbr)"
    print "    --wd work_directory (default .)"
    print


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:],
            "hd:w:m:b:",
            [ "help", "directory=", "basenames=", "wd=" ]
            )
    except getopt.GetoptError, err:    # print help information and exit:
        print str(err) # print something like "option -a not recognized"
        usage()
        sys.exit(2)

    directory = None
    basename = None
    wdirectory = None

    for option,value in opts:
        if option in ("-h", "--help"):
            usage()
            sys.exit()
        elif option in ("-d","--directory"):
            directory = value
        elif option in ("-b","--basenames"):
            basename = value
        elif option in ("--wd"):
            wdirectory = value
        else:
            assert False, "unhandled option"

    if basename==None:
        usage()
        sys.exit(2)

    if directory==None: directory = '.'

    if wdirectory==None: wdirectory = directory
        
    # First load the main project - Get the model offset
    
    project = txbr.TxBRproject(directory,basename)
    project.load()

    s = project.stackedSeries()

    ntilt = s.numberOfExposures()
    nterm = s.projection.numberOfTerms

    P =  numpy.concatenate((s.projection.x_coefficients,s.projection.y_coefficients))
    P.resize((2,ntilt,nterm))
    P = numpy.rollaxis(P,1,3)

    txbr.utilities.plot_projections( P, P, directory=wdirectory, basename=basename, save=False, show=True )
    

main()


