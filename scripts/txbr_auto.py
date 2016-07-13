#!/usr/bin/python

"""
The **txbr_auto.py** script allows to run the automatic TxBR alignment
"""

import sys
import shutil
import os.path
import getopt
import numpy
import txbr.align
import txbr.onthefly
import txbr.utilities
import modl

from txbr.txbrconf import log
from txbr.txbrdao import loadProject
from txbr.txbrdao import loadMap

from txbr.utilities import TXBR_EXTENSION
from txbr.utilities import MOD_EXTENSION
from txbr.utilities import MARKER_EXTENSION

def usage():
    """Usage for this script"""

    print
    print "Usage: {0:s}.py -b basename [Options]".format(os.path.basename(sys.argv[0]))
    print
    print "    -d directory"
    print "        Directory of the TxBR project."
    print "    --wd work_directory"
    print "        Work directory used for this TxBR project"
    print "    --bin n"
    print "        Bin the aw data in the x and y directions"
    print "    --scope scope-id"
    print "        Scope selection. Default value is 'fei_titan'."
    print "    --local"
    print "        Local approximation for the auto-tomography process. By default, a global approximation is used"
    print "    --strict"
    print "        During the minimization sequence, no restriction is made on the number of beads to pick."
    print "    --reset"
    print "        The link files are recalculated at the start."
    print "    --size [default:6]"
    print "        Pixel size for the bead marker to use in the detection routine"
    print "    --lc [default:0.5]"
    print "        The minimum coverage ratio (lentgh/bumber of views) for a trace to be considered."
    print "    --nc [default:100]"
    print "        The number of projection seeds to use during the first stage of the alignment."
    print "    --nmax [default:250]"
    print "        The maximum number of projection seeds to use during the first stage of the alignment."
    print "    --test-detection"
    print "        Test the detection routine"
    print "    --slice index"
    print "        Run the detection scheme on a given slice"
    print "    --rough"
    print "        For a rough alignment"
    print "    --refine"
    print "        For refining the alignement by including more tracks"
    print "    -h"
    print "        Help usage"


def __loadMap__( workDirectory, basename, bin=None ):

    cfg_file_path = os.path.join( workDirectory, basename + TXBR_EXTENSION )

    if not os.path.lexists(cfg_file_path):
        return None

    P_ = loadMap( cfg_file_path )
    if bin!=None and bin!=1: P_.rescale( bin, bin, bin )
    projmap = {}
    for index in range(P_.numberOfExposures()):
        projmap[index] = P_.getAffineTransformation(index)

    return projmap


def __load3DPoints__( workDirectory, basename, bin=None ):

    mod_file_path = os.path.join( workDirectory, "txbr-align", "bin{0:}".format(bin), basename + ".mod" )
    
    return modl.loadAllPoints( mod_file_path )


def __build_reco__( basename, directory, workDirectory, bin, scope ):
    """To build a QuickReco object."""

    log.info( "Build reconstruction object: (directory,basename)=({0:s},{1:s})".format(directory,basename) )

    filename = "{0:s}.st".format(basename)

    if not os.path.exists(os.path.join(directory,filename)):
        sys.exit("File {0:s} does not exist!".format(file))

    stack = txbr.utilities.loadStack( directory, filename, workDirectory=workDirectory, scope=scope, bin=bin )

    if stack.isEmpty():
        print "The stack {0:s} is empty or does not exist!".format(filename)
        return

    reco = txbr.onthefly.QuickReconstruction( stack, scope=scope )

    # Load or recalculate (from rough cross-correlation) the projection map

    cfg_file_path = os.path.join( reco.workDirectory, basename + TXBR_EXTENSION )

    if os.path.lexists(cfg_file_path):
        P_ = loadMap( cfg_file_path )
        P_.rescale( bin, bin, bin)
        projmap = {}
        for index in range(P_.numberOfExposures()):
            projmap[index] = P_.getAffineTransformation(index)
        reco.prjmap = projmap
    else:
        projmap = reco.align( saveAsTxBRProject=True )

    return reco


def run_test_detection( basename, directory, workDirectory, bin, scope, size, nmin, slice, **args ):
    """Test detection scheme for the beads in the middle slice of the data"""
    
    log.info( "Test bead detection routine: (directory,basename)=({0:s},{1:s})".format(directory,basename) )

    reco = __build_reco__( basename, directory, workDirectory, bin, scope )

    if slice==None: slice = reco.stack.getReferenceIndex()

    reco.detect( size, nmin=nmin, views=[slice], overwrite=True, doShow=True )


def run_prealignment( basename, directory, workDirectory, bin, scope, **args ):
    """Run the prealignment and create the associated data stack."""

    log.info( "Test bead detection routine: (directory,basename)=({0:s},{1:s})".format(directory,basename) )

    reco = __build_reco__( basename, directory, workDirectory, bin, scope )

    reco.createStack( type="preali", showVol=True )


def run_autotomo( basename, directory=".", workDirectory=".", bin=4, scope="fei_titan", size=6, lc=0.5,
                  nc=100, nmin=0, nmax=250, local=False, closeby=False, strict=False, reset=False, rough_only=False,
                  refine_only=False, **args ):
    """Run the autotomography process."""

    log.info( "Run automo routine: (directory,basename)=({0:s},{1:s})".format(directory,basename) )

    reco = __build_reco__( basename, directory, workDirectory, bin, scope )

    # (i) Markers detection / Rough alignment

    points = reco.detect( size=size, nmin=nmin )
    projmap = reco.prjmap

    # (ii) Use an adequate TxBR

    project = loadProject(directory=directory,basenames=basename,work_directory=workDirectory,bin=bin)

    # (iii) Get the frame information

    origin,end = reco.getFrame()

    xmin,ymin,zmin = origin
    xmax,ymax,zmax = end

    log.info("Frame: x:{0:f}:{1:f}  x:{2:f}:{3:f}  x:{4:f}:{5:f}".format(xmin,xmax,ymin,ymax,zmin,zmax))

    # Perform a quick alignement based on the detected markers with the rough projection map
    # The two methods [global_align,local_align] handle the sequential steps to perform
    # The alignment on a QuickAlign object *qa*

    def global_align( **options ):  # Method handling the global alignment

        if strict: qa.nmax = None
        else: qa.nmax = nmax

        # Tolerance criteria for ( initialization, first pass, transition, final pass)
        
        dist2d_seq = (10.0,10.0,9.0,8.0)

        qa.minimumClusterSize = 2
        qa.crossDistanceThreshold = 4.0

        doShortenTracks = True

        # First pass: - Orthogonal alignment
        #             - Restrict number of markers
        #             - Restrict projection in the bundle adjustement
        # Not shortening track may help in picking sparse markers on one side...

        NUMBER_OF_PROJ_FOR_BA = 20  # Only a few projections for the first pass
        decimate = int(len(projmap)/NUMBER_OF_PROJ_FOR_BA)

        if not refine_only:

            qa.globalAlignement( nc=nc, lc=lc, dist2d=dist2d_seq[0] )

            log.info("Run standard alignment scheme (initialization)...")

            qa.projmap, _3dpts = options['initial_align_scheme']( project, decimate=decimate, init=True )

            tracks = qa.track( _3dpts, r=dist2d_seq[1] )[0]
            qa.saveTracks( tracks, lc=lc, doShortenTracks=doShortenTracks  )

            iter,niter = 0,7
            lcurrent,lnext = -1,0
            while lnext>lcurrent and iter<niter:
                log.info("Iteration #{0:}: Mean Length={1:}->{2:}".format(iter,lcurrent,lnext))
                iter += 1
                qa.saveTracks( tracks, lc=lc, doShortenTracks=doShortenTracks )
                qa.projmap, _3dpts = options['initial_align_scheme']( project, decimate=decimate, init=True )
                tracks = qa.track( _3dpts, r=dist2d_seq[2] )[0]
                L = qa.lengthDistributionInfo( tracks )
                lcurrent = lnext
                lnext = numpy.mean(L)

            # Transition step: include all projections

            qa.projmap, _3dpts = options['final_align_scheme']( project, decimate=1, init=True )

        if rough_only: return

        # final step: include all projections with a max of tracks

        qa.nmax = None
        qa.minimumClusterSize = 3
        qa.crossDistanceThreshold = 5.0

        tracks = qa.globalAlignement( nc=None, lc=lc, dist2d=dist2d_seq[2] )

        iter,niter = 0,4
        lcurrent, lnext = -1, 0
        while lnext>lcurrent and iter<niter:
            log.info("Iteration #{0:}: Mean Length={1:}->{2:}".format(iter,lcurrent,lnext))
            iter += 1
            qa.saveTracks( tracks, lc=lc, doShortenTracks=doShortenTracks )
            qa.projmap, _3dpts = options['final_align_scheme']( project, decimate=decimate, init=True )
            tracks = qa.track( _3dpts, r=dist2d_seq[3] )[0]
            L = qa.lengthDistributionInfo( tracks )
            lcurrent = lnext
            lnext = numpy.mean(L)

        qa.projmap, _3dpts = options['final_align_scheme']( project, decimate=1, init=True )
        tracks = qa.track( _3dpts, r=dist2d_seq[3] )[0]


    def local_align( **options ):     # Method handling the local alignment

        qa.localAlignement( loadLinkFiles=not options['reset'], nc=None )

        qa.projmap, _3dpts = options['final_align_scheme']( project )
        tracks = qa.track( _3dpts, r=2.5 )[0]
        qa.saveTracks( tracks, lc=0.5, doShortenTracks=True  )


    # Run the alignment scheme

    options = {}

    if local:
        method = txbr.onthefly.LOCAL_ALIGN
        align2run = local_align
        options['reset'] = reset
    else:
        method = txbr.onthefly.GLOBAL_ALIGN
        align2run = global_align
        options['initial_align_scheme'] = txbr.align.runOrthogonalAlignment
        options['final_align_scheme'] = txbr.align.runRegularAlignment
        options['rough_only'] = rough_only
        options['refine_only'] = refine_only

    qa = txbr.onthefly.QuickAlign( basename, points, projmap, xmin, xmax, ymin, ymax,
                                   directory=directory, workDirectory=workDirectory,
                                   bin=bin, method=method )
                                   
    align2run( **options )

    # Make a copy of the individual configuration file

    mod_file = os.path.join( project.align_directory, basename + MOD_EXTENSION )
    if os.path.lexists(mod_file): shutil.copyfile( mod_file , mod_file + ".indiv" )

    mrk_file = os.path.join( project.align_directory, basename + MARKER_EXTENSION )
    if os.path.lexists(mrk_file): shutil.copyfile( mrk_file , mrk_file + ".indiv" )

    txbr_file = os.path.join( project.work_directory, basename + TXBR_EXTENSION )
    if os.path.lexists(txbr_file): shutil.copyfile( txbr_file , txbr_file + ".indiv.bin{0:}".format(bin) )



def main():
    '''Main routine for running the auto-tomography module'''

    try:
        opts, args = getopt.getopt( sys.argv[1:], "hb:d:",
            [ 'bin=', 'scope=', 'wd=', 'size=', 'lc=', 'nc=',  'nmin=', 'nmax=', 'local', 'test-detection', 'slice=', 'strict', 'refine', 'rough', 'reset', 'closeby', 'help'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    method2run = run_autotomo

    keywords = {}
    keywords['basename'] = None
    keywords['directory'] = "."
    keywords['workDirectory'] = None
    keywords['bin'] = 1
    keywords['scope'] = "fei_titan"
    keywords['size'] = 6
    keywords['lc'] = 0.3
    keywords['nc'] = 150
    keywords['nmin'] = 0
    keywords['nmax'] = 250
    keywords['local'] = False
    keywords['closeby'] = False
    keywords['slice'] = None
    keywords['strict'] = False
    keywords['reset'] = False
    keywords['rough_only'] = False
    keywords['refine_only'] = False

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-b"):
            keywords['basename'] = value_
        elif option_ in ("-d"):
            keywords['directory'] = value_
        elif option_ in ("--wd"):
            keywords['workDirectory'] = value_
        elif option_ in ("--size"):
            try:
                keywords['size'] = float(value_)
            except ValueError as e:
                print "Value Error for option {0:s}: {1:s}".format(option_,e)
                pass
        elif option_ in ("--lc"):
            try:
                keywords['lc'] = float(value_)
            except ValueError as e:
                print "Value Error for option {0:s}: {1:s}".format(option_,e)
                pass
        elif option_ in ("--nc"):
            try:
                keywords['nc'] = int(value_)
            except ValueError as e:
                print "Value Error for option {0:s}: {1:s}".format(option_,e)
                pass
        elif option_ in ("--nmin"):
            try:
                keywords['nmin'] = int(value_)
            except ValueError as e:
                print "Value Error for option {0:s}: {1:s}".format(option_,e)
                pass
        elif option_ in ("--nmax"):
            try:
                keywords['nmax'] = int(value_)
            except ValueError as e:
                print "Value Error for option {0:s}: {1:s}".format(option_,e)
                pass
        elif option_ in ("--bin"):
            try:
                keywords['bin'] = int(value_)
            except ValueError as e:
                print "Value Error for option {0:s}: {1:s}".format(option_,e)
                sys.exit(1)
        elif option_ in ("--scope"):
            keywords['scope'] = value_
        elif option_ in ("--slice"):
            keywords['slice'] = int(value_)
        elif option_ in ("--local"):
            keywords['local'] = True
        elif option_ in ("--closeby"):
            keywords['closeby'] = True
        elif option_ in ("--strict"):
            keywords['strict'] = True
        elif option_ in ("--refine"):
            keywords['refine_only'] = True
        elif option_ in ("--rough"):
            keywords['rough_only'] = True
        elif option_ in ("--reset"):
            keywords['reset'] = True
        elif option_ in ("--test-detection"):
            method2run = run_test_detection
        else:
            assert False, "Unhandled option [{):s}]".format(option_)

    if keywords['basename']==None:
        sys.exit("Please provide a basename!")

    if keywords['workDirectory']==None:
        keywords['workDirectory'] = keywords['directory']

    method2run(**keywords)


if __name__ == '__main__': main()