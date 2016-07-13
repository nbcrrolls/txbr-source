#!/usr/bin/python

import sys
import os
import os.path
import getopt
import numpy
import txbr.utilities
import txbr.onthefly

from txbr import log

global OPTIONS
global CROSS_CORRELATION_OPT,BEAD_TRACKING_OPT

CROSS_CORRELATION_OPT = "cross-correlation"
BEAD_TRACKING_OPT = "bead-tracking"

OPTIONS = [ CROSS_CORRELATION_OPT, BEAD_TRACKING_OPT ]

def usage():

    print
    print 'Usage: %s option -b basename1[,...] [Options]' %os.path.basename(sys.argv[0])
    print "          option in %s" %(", ".join(OPTIONS))
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files (st)"
    print "    --wd work_directory (default to directory)"
    print "        Name of the directory where the output files (binnned) are stored"
    print "    --bead-size size"
    print "        Size of the beads to track"
    print "    --scope"
    print "        Pull the configuration file of a certain electron microscope"
    print


def execute(cmd):
    '''Execute a system command.'''
    
    log.info(cmd)
    
    code = os.system(cmd)

    if code!=0: sys.exit("An error happened for command: %s" %cmd)


def track( basename, beadsize, scope, firstTiltAngle=-60, tiltIncrement=2 ):
    """Wrapper to etomo tracking"""


    rot_axis = txbr.utilities.loadRotationAxis( scope )

    if rot_axis[1]==0:
        angle = numpy.degrees(numpy.pi/2.0)
    else:
        angle = numpy.degrees(numpy.arctan(-rot_axis[0]/rot_axis[1]))

    print "Angle: %f" %angle

    opts = []
    opts.append('FillGaps')

    parameters = {}
    
    parameters['ImageFile'] = '%s.preali' %basename
    parameters['ImagesAreBinned'] = "1"
    parameters['InputSeedModel'] = '%s.seed' %basename
    parameters['OutputModel'] = '%s.fid' %basename
    parameters['RotationAngle'] = "%s" %angle # Angle in degrees
    parameters['FirstTiltAngle'] = firstTiltAngle
    parameters['TiltIncrement'] = tiltIncrement
    parameters['TiltDefaultGrouping'] = "7"
    parameters['MagDefaultGrouping'] = "5"
    parameters['RotDefaultGrouping'] = "1"
    parameters['BeadDiameter'] = beadsize
    parameters['MaxGapSize'] = "5"
    parameters['RoundsOfTracking'] = "2"
    parameters['LocalAreaTracking'] = "0"
    parameters['LocalAreaTargetSize'] = "1000"
    parameters['MinBeadsInArea'] = "8"
    parameters['MinOverlapBeads'] = "5"
    parameters['MinViewsForTiltalign'] = "4"
    parameters['MinTiltRangeToFindAxis'] = "10.0"
    parameters['MinTiltRangeToFindAngles'] = "20.0"
    parameters['BoxSizeXandY'] = "80,80"
    parameters['MaxBeadsToAverage'] = "4"
    parameters['PointsToFitMaxAndMin'] = "7,3"
    parameters['DensityRescueFractionAndSD'] = "0.6,1.0"
    parameters['DistanceRescueCriterion'] = 10.0
    parameters['RescueRelaxationDensityAndDistance'] = "0.7,0.9"
    parameters['PostFitRescueResidual'] = "2.5"
    parameters['DensityRelaxationPostFit'] = "0.9"
    parameters['MaxRescueDistance'] = "2.5"
    parameters['ResidualsToAnalyzeMaxAndMin'] = "9,5"
    parameters['DeletionCriterionMinAndSD'] = "0.04,2.0"

    str = []
    str.extend([ "-%s" %option for option in opts ])
    str.extend([ "-%s %s" %(key,value) for key,value in parameters.items() ])
    str = " ".join(str)

    execute( "beadtrack %s" %str )

    

def main():
    """The main function to perform the bin down"""
    
    directory = "."
    work_directory = "."
    basename = None
    beadSize = 20
    scope = None

    flags1 = "hb:d:"
    flags2 = [ "help", "directory=", "wd=", "basename=", "bead-size=","scope=" ]

    try:
        option = sys.argv[1]
        opts, args = getopt.getopt(sys.argv[2:], flags1, flags2)
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    for option_,value_ in opts:
        if option_ in ('-h','--help'):
            usage()
            sys.exit()
        elif option_ in ('-d','--directory'):
            directory = value_
        elif option_ in ('--wd'):
            work_directory = value_
        elif option_ in ('-b','--basename'):
            basename = value_
        elif option_ in ('--bead-size'):
            beadSize = int(value_)
        elif option_ in ('--scope'):
            scope = value_
        else:
            assert False, "unhandled option"
            
    if basename==None:
        usage()
        sys.exit()
        
    directory = os.path.abspath(directory)
    work_directory = os.path.abspath(work_directory)

    basenames = txbr.utilities.extract_series_name_from_dir( basename, directory, extension=txbr.utilities.ST_EXTENSION )

    if option==CROSS_CORRELATION_OPT:

        for basename in basenames:
            stack = txbr.onthefly.MRCStack( "%s%s" %(basename,txbr.utilities.ST_EXTENSION))
            reco = txbr.onthefly.QuickReconstruction( stack, scope )
            reco.createPreAlignmentFile( buildStack=True, stretch=False )

    if option==BEAD_TRACKING_OPT:

        for basename in basenames:
            stack = txbr.onthefly.MRCStack( "%s%s" %(basename,txbr.utilities.ST_EXTENSION))
            track( basename, beadSize, scope)

main()