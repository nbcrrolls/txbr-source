#!/usr/bin/python

import sys
import os.path
import math
import getopt
import logging
import numpy
import mpi4py.MPI
import mrc
import txbr.prefil
import txbr.filter

from txbr import FILTER_DIR
from txbr import FILTER_MAGNIFICATION
from txbr import LOG_CONF
from txbr.txbrdao import loadProject

logging.config.fileConfig(LOG_CONF)

log = logging.getLogger('align')

log.warning("Check the sign for the filtering angle")
ANGLE_SIGN = -1.0

def usage():
    
    print
    print 'Usage: %s.py -b basename1[,...] [Options]' %os.path.basename(sys.argv[0])
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files"
    print "    -wd work_directory (default .)"
    print "        Name of the directory where the output files are stored"
    print "    --remap"
    print "        Apply a 2D polynomial remap to the micrographs"
    print "    --filter-mag order (default 2)"
    print "        Local magnification in the remap"
    print "    --split-series"
    print "        The series are filtered independantly and stored in separate files."
    print "    --mpi"
    print "        Run the filter within the MPI protocol"
    print "    -h or --help"
    print "        Help Information"
    

def __run_mpi_filter__( project, no_remap, filter_magnification, reset=False):
    '''
    Apply the filtering scheme on each sections of the data stack.

    :param reset: Enforce recalculating the existing filtered slice even if already calculated.
    :param doClean: Clean up the individual filtered slice on the calculation is finished.
    :param fine: The filtering scheme is implemented at the full resolution, and is binned afterwards.
    '''

    myrank = mpi4py.MPI.COMM_WORLD.Get_rank()
    nprocs = mpi4py.MPI.COMM_WORLD.Get_size()

    log.info("Filter procedure: MPI-rank={0:}  MPI-size:{1:}".format(myrank,nprocs))

    wd = project.filter_directory

    for series in project.series:

        if not series.enabled: continue
        if not series.input in ['st','preali']: continue
        
        ntilts = series.numberOfExposures()
        depth = float(ntilts)
        block = math.ceil(depth/nprocs)
        
        sec_start = int(myrank*block)
        sec_stop = int(min((myrank+1)*block,ntilts)-1)
        
        if sec_start>=ntilts: break

        if no_remap:
            angle = series.getFilteringAngle()
            angle = ANGLE_SIGN*numpy.degrees(angle)
            log.info('Rank %i: filter series %s at angle %f from %s file [%i,%i].' %(myrank,series.basename, angle, series.input, sec_start, sec_stop))
        else:
            log.info('Rank %i: filter series %s with remapping procedure from %s file [%i,%i].' %(myrank,series.basename, series.input, sec_start, sec_stop))

        if series.input=='st':
            base = series.basename + '.st'
            src = os.path.join(series.stack.data_dir, base)
        elif series.input=='preali':
            base = series.basename + '.preali'
            src = os.path.join(project.align_directory, base)

        for islice in range(sec_start,sec_stop+1):
            
            input = os.path.join(wd, base + '.%i' %islice)
            output = os.path.join(wd, base + '.SL.%i' %islice)
            if os.path.lexists(output) and not reset: continue
            mrc.newstack( [src], input, secs=[[islice]] )

            if no_remap:
                txbr.prefil.rot_filter( input, output, filter_magnification, angle )
            else:
                remap_file = os.path.join(wd, serie.basename + '.remap')
                txbr.prefil.remap_filter(input, output, *txbr.filter.load2DRemap(remap_file))

            os.remove(input) # No need to keep that input file anymore

    mpi4py.MPI.COMM_WORLD.Barrier()


def __finalize_mpi_filter__( project, bin=None, doClean=False ):
    '''
    Collect and merge the filtered slices into one stack.

    :param bin: if None, the *binning* factor is specified in the variable *project*.
    :param doClean: Clean up the intermediary filtered files.
    '''

    if mpi4py.MPI.COMM_WORLD.Get_rank()!=0: return

    log.info('Finalize filters')

    wd_src = project.filter_directory
    wd_dest = project.filter_directory

    if bin==None: bin=1
    
    if bin!=1:   # The full-resolution filtered data should be binned to *bin*
        filter_directory = os.path.join( project.work_directory, FILTER_DIR, "bin{0:}".format(bin) )
        if not os.path.lexists(filter_directory): os.makedirs(filter_directory)
        wd_dest = filter_directory

    inputs = []
    
    for series in project.series:

        if not series.enabled: continue
        if not series.input in ['st','preali']: continue

        ntilts = series.numberOfExposures()
        
        if series.input=='st':
            base = series.basename + '.st'
        elif series.input=='preali':
            base = series.basename + '.preali'
            
        src = os.path.join(series.directory, base)
        wd_base_src = os.path.join(wd_src, base)
        wd_base_dest = os.path.join(wd_dest, base)

        inputs.extend([ wd_base_src + '.SL.%i' %i for i in range(ntilts) ])
                
        log.info('src: {0:s}'.format(src))
        log.info('wd_base: src={0:s}  dest={1:s}'.format(wd_base_src,wd_base_dest))


    output = os.path.join( wd_dest, project.basename + '.st.SL' )
    preali_rot_SL = os.path.join( wd_dest, project.basename + '.preali.rot.SL' )

    mrc.newstack( inputs, output, bin=bin, clean=doClean )
    mrc.copy_dimensions(src,output)
    mrc.copy_description(src,output)
    mrc.add_label(output,"TxBR: Filtered")

    if os.path.lexists(preali_rot_SL): os.unlink(preali_rot_SL)

    log.info('Link %s to %s' %(preali_rot_SL,output))

    os.symlink( os.path.basename(output), preali_rot_SL )


def main():
    '''Main routine for this filtering module'''

    try:
        opts, args = getopt.getopt(  \
                 sys.argv[1:], \
                 "hd:b:", \
                 [ "help", "directory=", "basenames=", "wd=", \
                   "remap", "filter-mag=", "bin=", "fine", "reset", "doClean" ]
                    )
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    directory = None
    work_directory = None
    basenames = None
    bin = None
    fine = False
    reset = False
    doClean = False

    no_remap = True
    filter_magnification = FILTER_MAGNIFICATION

    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-d"):
            directory = value_
        elif option_ in ("--wd"):
            work_directory = value_
        elif option_ in ("-b","--basenames"):
            basenames = value_.split(',')
        elif option_ in ("--fine"):
            fine = True
        elif option_ in ("--remap"):
            no_remap = False
        elif option_ in ("--filter-mag"):
            filter_magnification = float(value_)
        elif option_ in ("--bin"):
            bin = int(value_)
        elif option_ in ("--reset"):
            reset = True
        elif option_ in ("--doClean"):
            doClean = True
        else:
            assert False, "unhandled option"
    
    if basenames==None:
        usage()
        sys.exit()

    if directory==None: directory = '.'
    if work_directory==None: work_directory = directory
        
    basename = ",".join(basenames)
    
    directory = os.path.abspath(directory)
    work_directory = os.path.abspath(work_directory)

    # In the *fine* version, the filter is calculated at the full resolution
    
    bin_eff = fine and 1 or bin
    doClean = doClean and not fine
    
    project = loadProject( directory=directory, basenames=basename, work_directory=work_directory, bin=bin_eff )

    # Create the filter directory if it does not exist
    if mpi4py.MPI.COMM_WORLD.Get_rank()==0:
        if not os.path.lexists(project.filter_directory):
            os.makedirs(project.filter_directory)

    mpi4py.MPI.COMM_WORLD.Barrier()

    __run_mpi_filter__( project, no_remap, filter_magnification, reset=reset )
    __finalize_mpi_filter__( project, bin=fine and bin or None, doClean=doClean )

    mpi4py.MPI.COMM_WORLD.Barrier()


if __name__ == '__main__': main()