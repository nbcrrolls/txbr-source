#!/usr/bin/python

import sys
import os
import math
import getopt
import time
import mpi4py.MPI
import mrc
import util
import txbr.bckprj
import txbr.combine
import txbr.join
import txbr.utilities

from txbr.txbrdao import loadProject

try:
    import txbr.bckprj_cu
except ImportError:
    print 'No GPU backprojection module available for this system!'

from txbr import log

SHOW_INDIVIDUAL_VOLUMES = False

def usage():

    print
    print 'Usage: %s.py -b basename1[,...] [Options]' %os.path.basename(sys.argv[0])
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files"
    print "    --wd work_directory (default to directory)"
    print "        Name of the directory where the output files are stored"
    print "    -x xstart,xstop,xinc (default from .txbr file)"
    print "        x properties of the reconstruction"
    print "    -y ystart,ystop,yinc (default from .txbr file)"
    print "        y properties of the reconstruction"
    print "    -z ztart,zstop,zinc (default from .txbr file)"
    print "        z properties of the reconstruction"
    print "    --without index1[,index2,...]"
    print "        List of indices of sections to skip within the backprojection"
    print "    -l blocksize (default from .txbr file)"
    print "        Size of z-blocks used during the backprojection"
    print "    --gpu"
    print "        Run the backprojection on a cluster of GPU boards"
    print "    --mpi"
    print "        Run the backprojection within the MPI protocol"
    print "    --mpinodes values"
    print "        "
    print "    --finalize"
    print "        Finalize the backprojection processes"
    print "    --doCombine"
    print "        Combine the different tilt series"
    print "    --doClean"
    print "        remove the temporary blocks"
    print "    --split-series"
    print "        The series are backprojected independantly [in case of a multiple tilt series]."
    print "    -h or --help"
    print "        Help Information"

    
def __execute__( cmd ):
    '''Execute a system command.

    :param cmd: the command to execute
    '''
    
    log.info(cmd)
    
    code = os.system(cmd)
    
    if code!=0: sys.exit("An error happened for command: %s" %cmd)
    
    
def __getGeneratedFiles__( directory, basenames, z, number_of_nodes, blocksize ):
    '''This routine generates a list containing the name of all the generated files after 
    the initial backprojection step is completed. Check if they exist or no.
    If there are missing files, they are inserted in a list returned as the second
    variable.

    :param directory: The directory to scan containing the volume backprojected blocks
    :param basenames: The *basenames* of the different tomographic series
    :param z: A dictionary containing the final z boundary conditions
    :param number_of_nodes: Number of nodes used if it was multiprocessed
    :param blocksize: Size of each block
    '''

    files = {}
    fails = {}
    
    for iseries,name in enumerate(basenames):

        log.debug("%s %i " %(name,iseries))

        files[name] = []
        fails[name] = []

        zstart, zstop, zinc = z[iseries]

        depth = zstop - zstart
        lpack = math.ceil(depth/number_of_nodes)

        number_of_blocks = int(math.ceil(lpack/blocksize))
        number_of_blocks = 1

        for inode in range(number_of_nodes):
            zpack = zstart + inode*lpack
            for iblock in range(number_of_blocks):
               # filename = "%s_z_%.2f-%i.out" %( name, zpack, iblock )
                filename = "%s_z_%.2f.out" %( name, zpack )
                filename = os.path.join( directory, filename )
                if not os.path.exists(filename):
                    log.warning("No file %s (%s)" %(filename,os.path.exists(filename)))
                    fails[name].append(filename)
                else:
                    log.warning("Added %s (%s)" %(filename,os.path.exists(filename)))
                    files[name].append(filename)

    for basename,missingFiles in fails.iteritems():
        if len(missingFiles)!=0:
            log.error( 'Error for series %s! The following files do not exist' %basename )
            log.error(fails)

    return files


def __deleteFiles__( inputFilesDict ):
    '''Function that allows to Delete some specified files.
    :param inputFilesDict: A dictionary containing the list of files to delete each series. The
        dictionary keys map the tilt series basename.
    '''

    for basename,series_filenames in inputFilesDict.iteritems():
        for filename in series_filenames:
            log.info("Remove file %s from the system..." %filename)
            os.remove(filename)
        
            
def __finalize_backprojection__( project, split_series, **keywords ):
    '''Collect and merge the filtered blocks.'''
    
    directory = keywords['directory']
    basename = keywords['basename']
    basenames = keywords['basenames']
    work_directory = keywords['work_directory']
    blocksize = keywords['blocksize']
    number_of_nodes = keywords['number_of_nodes']
    x = keywords['x']
    y = keywords['y']
    z = keywords['z']
    doCombine = 'doCombine' in keywords and keywords['doCombine']
    doClean = 'doClean' in keywords and keywords['doClean']

    if split_series:
        series2process = project.series
    else:
        series2process = [ project.stackedSeries() ]
        
    outputfiles = {}

    for series in series2process:

        log.info( "Finalize backprojection for %s" %basename )
        log.info( "between z=[%.1f,%.1f] with a node number of %i" %( z[0], z[1], number_of_nodes ) )

        reconstruction = series.getSubReconstruction( project.reconstruction.bottomPlane, project.reconstruction.topPlane )
        reconstruction.bin = project.bin

        reconstruction.origin.x = max(x[0],reconstruction.origin.x)
        reconstruction.end.x = min(x[1],reconstruction.end.x)
        reconstruction.increment.x = x[2]

        reconstruction.origin.y = max(y[0],reconstruction.origin.y)
        reconstruction.end.y = min(y[1],reconstruction.end.y)
        reconstruction.increment.y = y[2]

        reconstruction.origin.z = z[0]
        reconstruction.end.z = z[1]
        reconstruction.increment.z = z[2]

        output = reconstruction.getFileName()

        inputfiles = __getGeneratedFiles__( project.backprojection_directory, [series.basename], [z], number_of_nodes, blocksize )

        mrc.newstack( inputfiles[series.basename], output, bounds3D=project.reconstruction, clean=doClean )

        outputfiles[series.basename] = output

    # Eventually combine volumes from different series

    if doCombine and len(outputfiles)>1:

        xorigin = project.reconstruction.origin.x
        yorigin = project.reconstruction.origin.y
        zorigin = project.reconstruction.origin.z

        sx = project.reconstruction.scale.x
        sy = project.reconstruction.scale.y
        sz = project.reconstruction.scale.z

        input = os.path.join( outputfiles[basenames[0]] )
        output = project.reconstruction.getFileName()

        txbr.join.joinVolumes( directory, basenames, work_directory, output, zmin=z[0], zmax=z[1], bin=project.bin )
        
        mrc.update_scale( output, sx, sy, sz)
        mrc.update_origin( output, (1-xorigin)*sx, (1-yorigin)*sy, -zorigin*sz)
        mrc.copy_description( input, output )   # Copy or update some header information
        
    # Display the final volume and the fiducial model
    
    model3DFile = os.path.join( project.align_directory, '%s.mod' %project.basename )
    __execute__('imod %s %s' %( project.reconstruction.getFileName(), model3DFile ))

    if SHOW_INDIVIDUAL_VOLUMES:
        log.info("Show individual volumes.")
        for outfile in outputfiles:
            __execute__('imod %s %s' %( outfile, model3DFile ))


def run_backprojection( project, z, skip, use_fixed_segment_size, split_series, **keywords ):
    '''Routine that runs a standard backprojection.'''
      
    if split_series:
        series2process = project.series
    else:
        series2process = [ project.stackedSeries() ]

    x = keywords["x"]
    y = keywords["y"]
#    z = keywords["z"]

    for series in series2process:

        projection = series.projection

        reconstruction = series.getSubReconstruction( project.reconstruction.bottomPlane, project.reconstruction.topPlane )
        reconstruction.bin = project.bin

        reconstruction.origin.x = max(x[0],reconstruction.origin.x) or reconstruction.origin.x
        reconstruction.end.x = min(x[1],reconstruction.end.x)
        reconstruction.increment.x = x[2]

        reconstruction.origin.y = max(y[0],reconstruction.origin.y)
        reconstruction.end.y = min(y[1],reconstruction.end.y)
        reconstruction.increment.y = y[2]

        reconstruction.origin.z = z[0]
        reconstruction.end.z = z[1]
        reconstruction.increment.z = z[2]

        log.info('Run TxBR backprojection for series %s between [%f,%f]' %(series.basename,z[0],z[1]))

        txbr.bckprj.reconstruct( project.filter_directory, series.basename, project.backprojection_directory, project.work_directory, skip, use_fixed_segment_size, projection, reconstruction )

    # Finalize the backprojection

    __finalize_backprojection__( project, split_series, **keywords )


def run_MPI_backprojection( project, split_series, **keywords):
    """Run the back-projection step by using MPI to split the processes

    :param split_series: If true the different tilt series will be backprojected independently.
    """

    x = keywords['x']   # A 3-tuple containing the x boundaries and increment
    y = keywords['y']   # A 3-tuple containing the y boundaries and increment
    z = keywords['z']   # A 3-tuple containing the z boundaries and increment
    skip = keywords['skip']

    myrank = mpi4py.MPI.COMM_WORLD.Get_rank()
    nprocs = mpi4py.MPI.COMM_WORLD.Get_size()

    zstart,zstop,zinc = z

    depth = z[1]-z[0]
    block = math.ceil(depth/nprocs)

    zstart = z[0] + myrank*block
    zstop = min(zstart+block-1,z[1])
    zinc = z[2]

    if split_series:
        series2process = project.series
    else:
        series2process = [ project.stackedSeries() ]

    for series in series2process:

        projection = series.projection

        reconstruction = series.getSubReconstruction( project.reconstruction.bottomPlane, project.reconstruction.topPlane )
        reconstruction.bin = project.bin

        reconstruction.origin.x = max(x[0],reconstruction.origin.x)
        reconstruction.end.x = min(x[1],reconstruction.end.x)
        reconstruction.increment.x = x[2]

        reconstruction.origin.y = max(y[0],reconstruction.origin.y)
        reconstruction.end.y = min(y[1],reconstruction.end.y)
        reconstruction.increment.y = y[2]

        reconstruction.origin.z = zstart
        reconstruction.end.z = zstop
        reconstruction.increment.z = zinc

        log.info("Process %i/%i running between z=[%.1f,%.1f] for series %s" %( myrank, nprocs, zstart, zstop, series.basename ))

        txbr.bckprj.reconstruct( project.filter_directory, series.basename, project.backprojection_directory, project.work_directory,
                                 skip, keywords['use_fixed_segment_size'], projection, reconstruction )

    mpi4py.MPI.COMM_WORLD.Barrier()

    if mpi4py.MPI.COMM_WORLD.Get_rank()==0:
        keywords['number_of_nodes'] = nprocs
        __finalize_backprojection__( project, split_series, **keywords )

    mpi4py.MPI.COMM_WORLD.Barrier()

    
def run_GPU_backprojection( project, z, skip, use_fixed_segment_size, split_series, **keywords ):

#    if not util.isGPUEnabled():
#        return run_backprojection(keywords)
    
    myrank = mpi4py.MPI.COMM_WORLD.Get_rank()
    nprocs = mpi4py.MPI.COMM_WORLD.Get_size()
    procnm = mpi4py.MPI.Get_processor_name()

    n_GPU_tot = util.numberOfGPUBoards()
    
    zstart = float(keywords['z'][0])
    zstop = float(keywords['z'][1])  
    
    depth = (zstop-zstart)
    
    blockForGPU = math.ceil(depth/n_GPU_tot)
    
    print 'blockForGPU: %i    n_GPU_tot:%i' %(blockForGPU,n_GPU_tot)
        
    nodes, numberOfGPUBoardsForNode, infoGPUBoard = util.infoForGPUHosts()

    node, deviceID = infoGPUBoard[myrank]

    keywords['zstart'] = zstart + myrank*blockForGPU
    keywords['zstop'] = keywords['zstart'] + blockForGPU - 1
    keywords['deviceID'] = deviceID
    
    log.info('Device %i of node %s reconstructs slices for z in [%f,%f]' %(deviceID,node,keywords['zstart'],keywords['zstop']))
       
    basenames = keywords['basenames']
    directory = keywords['directory']
    work_directory = keywords['work_directory']
    z0 = keywords['zstart']
    z1 = keywords['zstop']
    zinc = 1.0
    skip = keywords['skip']
    use_fixed_segment_size = keywords['use_fixed_segment_size']
    
    deviceID = keywords['deviceID']

    if split_series:
        series2process = project.series
    else:
        series2process = [ project.stackedSeries() ]

    for series in series2process:

        log.info('Run TxBR backprojection for series %s between [%f,%f]' %(basename,z0,z1))

        txbr.bckprj_cu.reconstruct( project.filter_directory, series.basename, project.backprojection_directory, 1, z0, z1, zinc, deviceID )
    
        
def main():
    '''The main routine'''
    
    flags1 = "hd:b:wd:l:x:y:z:n:"
    flags2 = [ "help", "directory=", "wd=", "work-directory=", "basenames=", "bin=", \
               "blocksize", "split-series","without=", "mpi", "mpinodes=", "gpu", "gpuid=", "finalize", \
               "doCombine", "doClean" ]

    try:
        opts, args = getopt.getopt(sys.argv[1:], flags1, flags2)
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
        
    keywords = {}

    keywords['directory'] = "."
    keywords['basename'] = None
    keywords['bin'] = None
    keywords['deviceID'] = None
    keywords['blocksize'] = None
    keywords['number_of_nodes'] = 1
    keywords['x'] = None
    keywords['y'] = None
    keywords['z'] = None
    keywords['skip'] = []
    keywords['doCombine'] = True
    keywords['doClean'] = True
    keywords['use_fixed_segment_size'] = False
    keywords['split_series'] = False

    methods = { "regular": run_backprojection,
                "mpi": run_MPI_backprojection,
                "gpu": run_GPU_backprojection,
                "finalize": __finalize_backprojection__ }

    method2run = methods["regular"] # default backprojection

    for option,value in opts:
        if option in ("-h", "--help"):
            usage()
            sys.exit()
        elif option in ("-d","--directory"):
            keywords['directory'] = value
        elif option in ("--wd", "--work-directory"):
            keywords['work_directory'] = value
        elif option in ("-b","--basenames"):
            keywords['basename'] = value
        elif option in ("--bin"):
            keywords['bin'] = int(value)
        elif option=="--gpuid":
            keywords['deviceID'] = int(value)
        elif option in ('-l','--blocksize'):
            keywords['blocksize'] = int(value)
        elif option in ('-n','--mpinodes'):
            keywords['number_of_nodes'] = int(value)
        elif option in ("-x"):
            try:
                keywords['x'] = [ float(s) for s in value.split(',') ]
            except TypeError:
                keywords['x'] = None
        elif option in ("-y"):
            try:
                keywords['y'] = [ float(s) for s in value.split(',') ]
            except TypeError:
                keywords['y'] = None
        elif option in ("-z"):
            try:
                keywords['z'] = [ float(s) for s in value.split(',') ]
            except TypeError:
                keywords['z'] = None
        elif option in ('--without'):
            if ':' in value:
                start,stop,step = value.split(':')
                keywords['skip'] = range(int(start),int(stop),int(step))
            else:
                keywords['skip'] = [int(s) for s in value.split(',')]
        elif option in ('--doCombine'):
            keywords['doCombine'] = True
        elif option in ('--doClean'):
            keywords['doClean'] = True
        elif option in ("--gpu"):
            method2run = methods["gpu"]
        elif option in ("--mpi"):
            method2run = methods["mpi"]
        elif option in ("--finalize"):
            method2run = methods["finalize"]
        else:
            assert False, "unhandled option"

    if not 'basename' in keywords:
        usage()
        sys.exit()

    if not 'work_directory' in keywords:
        keywords['work_directory'] = keywords['directory']
        
    keywords['basenames'] = txbr.utilities.extract_series_name_from_dir( keywords['basename'], keywords['work_directory'] )

    t_start = time.time()

    project = loadProject( directory=keywords['directory'],
                           basenames=keywords['basename'],
                           work_directory=keywords['work_directory'],
                           bin=keywords['bin'] )

    if mpi4py.MPI.COMM_WORLD.Get_rank()==0: # Create the backprojection directory
        if not os.path.lexists(project.backprojection_directory):
            os.makedirs(project.backprojection_directory)

    mpi4py.MPI.COMM_WORLD.Barrier()
      
    if 'x' not in keywords or keywords['x']==None: # Global set of series might have different x values
        keywords['x'] = [ project.reconstruction.origin.x,
                          project.reconstruction.end.x,
                          project.reconstruction.increment.x ]
    else:
        project.reconstruction.origin.x,project.reconstruction.end.x,project.reconstruction.increment.x = keywords['x']
        
    if 'y' not in keywords or keywords['y']==None: # Series might have different y values
        keywords['y'] = [ project.reconstruction.origin.y,
                          project.reconstruction.end.y,
                          project.reconstruction.increment.y ]
    else:
        project.reconstruction.origin.y,project.reconstruction.end.y,project.reconstruction.increment.y = keywords['y']

    if 'z' not in keywords or keywords['z']==None: # Global set of series might have different z values
        keywords['z'] = [ project.reconstruction.origin.z,
                          project.reconstruction.end.z,
                          project.reconstruction.increment.z ]
    else:
        project.reconstruction.origin.z,project.reconstruction.end.z,project.reconstruction.increment.z = keywords['z']
        
    if keywords['blocksize']==None:
        keywords['blocksize'] = project.reconstruction.sizeOfBlock

    method2run( project, **keywords )
    
    t_stop = time.time()

    log.info('Total Time for the backprojection: %f' %(t_stop-t_start))


if __name__ == '__main__': main()

