#!/usr/bin/python

import sys
import os
import os.path
import getopt
import txbr.onthefly

def usage():

    print 'Usage: %s -b basename[,...]  [Options]' % os.path.basename(sys.argv[0])
    print
    print "    -d directory (default .)"
    print "        Name of the directory for the \"on the fly\" reconstruction."
    print "        Images should be placed in a folder called \"src_jpegs\""
    print "    --wd work_directory (default to directory)"
    print "        Name of the directory where the output files are stored"
    print "    --scope"
    print "        Pull the configuration file of a certain electron microscope"
    print "    --status"
    print "        Gives a status on the on-fly reconstruction"
    print "    --reset"
    print "        Kill the on the fly process and clear the reconstruction files"
    print "    --extension %s (default) or %s" %( txbr.onthefly.JPG, txbr.onthefly.TIF)
    print "        Type of the tilt files"
    print "    --reset-names"
    print "        Reset the names to be processed by on-the-fly TxBR (append angles)"
    print "    --angles start,increment"
    print "        Start angle and increment angle"
    print


def main():

    try:
        
        flags1 = "hb:d:"
        flags2 = [ "help", "directory=", "wd=", "basename=", "scope=", "extension=", "reset", "reset-names", "angles=", "status" ]

        opts, args = getopt.getopt( sys.argv[1:], flags1, flags2 )

    except getopt.GetoptError, err:

        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    status = False  # Gives a status on the on-the fly reconstruction
    reset = False   # Clean up or no a reconstruction on the fly
    resetNames = False

    directory = "."
    workDirectory = "."
    basenames = None
    basename = None
    extension = txbr.onthefly.JPG
    scope = None
    
    for option,value in opts:
        if option in ('-h','--help'):
            usage()
            sys.exit()
        elif option in ('-d','--directory'):
            directory = value
        elif option in ('--wd'):
            workDirectory = value
        elif option in ('-b','--basename'):
            basenames = [str(s) for s in value.split(',')]
        elif option in ('--scope'):
            scope = value
        elif option in ('--extension'):
            extension = value
        elif option in ('--status'):
            status = True
        elif option in ('--reset'):
            reset = True
        elif option in ('--reset-names'):
            resetNames = True
        elif option in ('--angles'):
            start_angle, inc_angle = [ float(s) for s in value.split(',') ]
            
    if basenames==None and not (reset or status):
        usage()
        sys.exit()
        
    directory = os.path.abspath(directory)
    workDirectory = os.path.abspath(workDirectory)

    if basenames!=None:
        basename = ','.join(basenames)

    if not os.path.lexists(workDirectory):
        os.mkdir(workDirectory)
    
    if resetNames:
        txbr.onthefly.reset( basename, extension=extension, start_angle=start_angle, inc_angle=inc_angle )
    elif reset:
        stack = txbr.utilities.ImageStack( basename, directory, workDirectory, extension=txbr.utilities.JPG, scope="fei.titan")
        rec = txbr.onthefly.QuickReconstruction( stack, scope=scope )
     #   rec = txbr.onthefly.ReconstructionOnTheFly( directory, workDirectory, basename, scope=scope )
        rec.reset()
    elif status:
        stack = txbr.utilities.ImageStack( basename, directory, workDirectory, extension=txbr.utilities.JPG, scope="fei.titan")
        rec = txbr.onthefly.QuickReconstruction( stack, scope=scope )
        #rec = txbr.onthefly.ReconstructionOnTheFly( directory, workDirectory, basename, scope=scope )
        rec.update_status()
    else:    # On the fly reconstruction is the default option
        stack = txbr.utilities.ImageStack( basename, directory, workDirectory, extension=txbr.utilities.JPG, scope="fei.titan")
        rec = txbr.onthefly.QuickReconstruction( stack, scope="fei.titan" )
#        rec = txbr.onthefly.ReconstructionOnTheFly( directory, workDirectory, basename, extension=extension, scope=scope )
        rec.update()


main()