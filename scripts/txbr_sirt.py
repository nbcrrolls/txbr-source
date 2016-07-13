#!/usr/bin/python

import os
import os.path
import sys
import getopt
import glob
import multiprocessing
import shutil
import logging
import copy
import mpi4py.MPI
import numpy
import scipy.ndimage.filters
import mrc
import util
import txbr.txbrconf
import txbr.txbrdao
import txbr.bckprj

log = logging.getLogger()

def usage():
    """Usage for the command *txbr_sirt.py*."""

    print 'Usage: {0:s}.py -b basename[,...]  [Options]'.format(os.path.basename(sys.argv[0]))
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files (.fid and .txbr)"
    print "    --wd work_directory (default to directory)"
    print "        Name of the directory where the output files are stored"
    print "    -n niter [default: 25]"
    print "        Number of iteration"
    print "    --mux mix [default: 1.0]"
    print "        The mixing parameter used at each iteration."
    print "    --bin n"
    print "        Bin the data in the x and y directions"
    print "    --reset"
    print "        The iteration cycle restarts from scratch."
    print "    --no-adaptive"
    print "        In this case, the iterative method corresponds to a regular SIRT"
    print "        routine. Otherwise, the backprojection on a ray is weighted by the"
    print "        volume density."
    print "    --nproc value [default: 20 ]"
    print "        The number of CPUs to use."
    print "    --log"
    print "        Apply the logarithm function on the rawdata."
    print "    -h or --help"
    print "        Help Information"
    

def __stack__( files_in, file_out, bounds3D=None, stats=None, clean=True ):
    """
    Helper routine to stack several MRC files into a single MRC file.

    :param files_in: A list containing the paths of the MRC files to stack.
    :param file_out: The path of the output MRC file.
    :param clean: If true the input files are removed from the file system.
    """

#    if mpi4py.MPI.COMM_WORLD.Get_rank()!=0: return

    log.debug("Stack volume {0:s}!".format(file_out))

    files_in = [ f for f in files_in if os.path.lexists(f) ]
    mrc_in = [ mrc.MRCFile(fin) for fin in files_in ]

    nx,ny,nz = 0,0,0

    for fin in mrc_in:
        nx,ny = fin.nx,fin.ny
        nz += fin.nz

    if os.path.exists(file_out): os.remove(file_out)
        
    fout = mrc.MRCFile(file_out)
    fout.setHeader(nx,ny,nz)

    index = 0
    for fin in mrc_in:
        for iz in range(fin.nz):
            u = fin.getZSliceAt(iz)
            fout.setZSliceAt(index, u)
            index += 1

    fout.updateHeader()

    if bounds3D!=None:
        xstart,ystart,zstart = int(bounds3D.origin.x),int(bounds3D.origin.y),int(bounds3D.origin.z)
        xstop,ystop,zstop = int(bounds3D.end.x),int(bounds3D.end.y),int(bounds3D.end.z)
        xscale,yscale,zscale = int(bounds3D.scale.x),int(bounds3D.scale.y),int(bounds3D.scale.z)
        mrc.update_grid( file_out, xstart, ystart, zstart, xstop-xstart+1, ystop-ystart+1, zstop-zstart+1)
        mrc.update_scale( file_out, xscale, yscale, zscale)
        mrc.update_origin( file_out, (1-xstart)*xscale, (1-ystart)*yscale, -zstart*zscale)

    if stats!=None:
        mrc.update_header( file_out, stats )

    if clean:
        for fin in files_in: os.remove(fin)
            
    
def __add__( file1, file2, mux=1.0, method="a-SIRT", verbose=True ):
    """
    Helper function used in the iterative process to add the volumic correction
    to the main volumic solution.

    :param file1: The path to the volumic correction.
    :param file2: The path to the main volume. This will be modified accordingly to
    the volumic correction stored in file1.
    :param mux: The weighting used when adding the volumic correcion.
    :param method: The iterative method used can be *SIRT*, *a-SIRT*, *SMART* or *EM*.
        Default is *a-SIRT*.
    """

#    if mpi4py.MPI.COMM_WORLD.Get_rank()!=0: return

    log.debug("Add volumes {0:s} + ({2:}*) {1:s}!".format(file1,file2, mux))

    f1 = mrc.MRCFile(file1)
    f2 = mrc.MRCFile(file2)

    for slice in range(f2.nz):
        u1 = f1.getZSliceAt(slice)
        u2 = f2.getZSliceAt(slice)
        if verbose:
            mean1, std1, min1, max1 = numpy.mean(u1), numpy.std(u1), numpy.min(u1), numpy.max(u1)
            mean2, std2, min2, max2 = numpy.mean(u2), numpy.std(u2), numpy.min(u2), numpy.max(u2)
            log.info("Test: mean:{0:f} [{1:.2f}]  std: {2:f} [{3:.2f}]  min: {4:f} [{5:.2f}]  max: {6:f} [{7:.2f}]".format(mean1,100.0*mean1/mean2,std1,100.0*std1/std2,min1,100.0*min1/min2,max1,100.0*max1/max2))
        if method=="a-SIRT": # For the adaptive version
            f2.setZSliceAt(slice,mux*u1*u2+u2)  
        elif method=="EM": # Expectation maximization
            f2.setZSliceAt(slice,u1*u2)
        elif method=="SMART": # Expectation maximization
            f2.setZSliceAt(slice,u2*numpy.exp(u1))
        elif method=="SIRT":
            f2.setZSliceAt(slice,mux*u1+u2)

    f2.updateHeader()


def __susbtract__( file1, file2, alpha=1.0, skip=[], method="a-SIRT", verbose=True ):
    """
    Helper function used in the iterative process to substract the projective
    contribution of the main volume from the measured micrographs.

    :param file1: The path to MRC stack containing the original micrographs.
    :param file2: The path to MRC file containing the calculated projection of the main volume.
        This will be modified to only store the residual contibution in the projections that will
        be used to evaluate the next volumic correction.
    :param alpha: The weighting used when adding the volumic correcion. This depends on
        the backprojection normalization.
    :param method: The iterative method used can be *SIRT*, *a-SIRT*, *SMART* or *EM*.
        Default is *a-SIRT*.
    """

#    if mpi4py.MPI.COMM_WORLD.Get_rank()!=0: return

    log.debug("Substract volumes {0:s} and {1:s}!".format(file1,file2))

    apply_gaussian_filter = True
    sigma = 0.25

    f1 = mrc.MRCFile(file1)
    f2 = mrc.MRCFile(file2)

    for slice in range(f2.nz):
        if slice in skip: continue
        u1 = f1.getZSliceAt(slice)
        u2 = f2.getZSliceAt(slice)
        if method=="a-SIRT": # For the adaptive version
            u = numpy.where(u2==0,0.0,(u1-u2/alpha)/u2)  
        elif method=="EM": # Expectation maximization
            u = numpy.where(u2==0,0.0,(u1/u2))
        elif method=="SMART":
            u = numpy.where(u2==0,0.0,numpy.log(alpha*u1/u2))
        elif method=="SIRT":
            u = numpy.where(u2==0,0.0,u1-u2/alpha)
        if apply_gaussian_filter: u = scipy.ndimage.filters.gaussian_filter(u,sigma)
        if verbose:
            mean1, std1, min1, max1 = numpy.mean(u1), numpy.std(u1), numpy.min(u1), numpy.max(u1)
            mean, std, min, max = numpy.mean(u), numpy.std(u), numpy.min(u), numpy.max(u)
            log.info("Test: mean:{0:f} [{1:.2f}]  std: {2:f} [{3:.2f}]  min: {4:f} [{5:.2f}]  max: {6:f} [{7:.2f}]".format(mean,100.0*mean1/mean,std,100.0*std1/std,min,100.0*min1/min,max,100.0*max1/max))
        f2.setZSliceAt(slice,u)

    f2.updateHeader()


def __prepare_input_data__( project, original_stack, apply_flatten, apply_log, apply_inversion,
                          apply_normalization ):
    """
    Prepare the rawdata to be processed with an iterative approach.
    """

#    if mpi4py.MPI.COMM_WORLD.Get_rank()!=0: return

    log.debug("Prepare data for {0:s}!".format(project.basename))

    nx,ny = 0,0
    ntilt = 0
    
    for series in project.series:
        nx,ny = series.stack.getImageSize()
        ntilt += series.numberOfExposures()

    fout = mrc.MRCFile(original_stack)
    fout.setHeader(nx,ny,ntilt)
    index = 0
    for series in project.series:
        stack_in = series.stack
        for itilt in range(series.stack.getNumberOfImages()):
            u = stack_in.getImageAt(itilt)
            if apply_log: u = numpy.where(u>0,numpy.log(u),0.0)
            if apply_flatten: u = util.flattenImage(u, order=1)
            if apply_inversion: u = -u
            if apply_normalization: u = (u-numpy.min(u))/(numpy.max(u)-numpy.min(u))
            fout.setZSliceAt(index,u)
            index+=1
    fout.updateHeader()


def __project__( project, projection, reconstruction, tilts, skip, use_fixed_segment_size,
                 output, use_mpi=False ):
    """
    Project a volume using a multiprocessing routine. Each process handle a tilt range
    sepcified in the array *tilts*.
    """

    basename = project.basename

    if use_mpi: # Use MPI

        mpi4py.MPI.COMM_WORLD.Barrier()

        myrank = mpi4py.MPI.COMM_WORLD.Get_rank()

        if myrank<len(tilts):

            start,stop = tilts[myrank]

            proj_proc = copy.deepcopy(projection)

            proj_proc.x_coefficients = proj_proc.x_coefficients[start:stop]
            proj_proc.y_coefficients = proj_proc.y_coefficients[start:stop]
            proj_proc.scaling_coefficients = proj_proc.scaling_coefficients[start:stop]
            proj_proc.label = "{0:}-{1:}".format(start,stop)

            target = txbr.bckprj.project
            args = ( project.art_directory, basename, project.art_directory, project.work_directory, skip, use_fixed_segment_size, proj_proc, reconstruction )

            apply(target,args)

        mpi4py.MPI.COMM_WORLD.Barrier()

    else:   # Use the multiprocessing module

        processes = []

        for start,stop in tilts:

            proj_proc = copy.deepcopy(projection)

            proj_proc.x_coefficients = proj_proc.x_coefficients[start:stop]
            proj_proc.y_coefficients = proj_proc.y_coefficients[start:stop]
            proj_proc.scaling_coefficients = proj_proc.scaling_coefficients[start:stop]
            proj_proc.label = "{0:}-{1:}".format(start,stop)

            target = txbr.bckprj.project
            args = ( project.art_directory, basename, project.art_directory, project.work_directory, skip, use_fixed_segment_size, proj_proc, reconstruction )

            p = multiprocessing.Process(target=target,args=args)
            p.start()

            processes.append(p)

        for p in processes: p.join()

    # Restack the different projection MRC blocks

    if use_mpi and mpi4py.MPI.COMM_WORLD.Get_rank()==0 or not use_mpi:

        stacks = []
        for start,stop in tilts:
            filename = "{0:s}_{1:}-{2:}.project.st".format(basename,start,stop)
            stacks.append(os.path.join(project.art_directory,filename))

        __stack__( stacks, output )

    if use_mpi: mpi4py.MPI.COMM_WORLD.Barrier()


def __backproject__( project, projection, reconstruction, slices, skip, use_fixed_segment_size,
                     output, use_mpi=False ):
    """
    Backproject projection data into a 3D volume using a multiprocessing routine. Each
    process handle a section range sepcified in the array *slices*."""

    basename = project.basename

    if use_mpi: # Use MPI

        mpi4py.MPI.COMM_WORLD.Barrier()

        myrank = mpi4py.MPI.COMM_WORLD.Get_rank()

        if myrank<len(slices):

            start, stop = slices[myrank]

            reco_proc = reconstruction.copy()
            reco_proc.origin.z = start
            reco_proc.end.z = stop

            target = txbr.bckprj.reconstruct
            args = (project.art_directory, basename, project.art_directory, project.work_directory, skip, use_fixed_segment_size, projection, reco_proc)

            apply(target,args)

        mpi4py.MPI.COMM_WORLD.Barrier()

    else: # Use the multiprocessing module

        processes = []
        vols = []

        for start,stop in slices:

            reco_proc = reconstruction.copy()
            reco_proc.origin.z = start
            reco_proc.end.z = stop

            target = txbr.bckprj.reconstruct
            args = ( project.art_directory, basename, project.art_directory, project.work_directory, skip, use_fixed_segment_size, projection, reco_proc )

            p = multiprocessing.Process(target=target,args=args)
            p.start()

            processes.append(p)

        for p in processes: p.join()

    # Restack the different backprojection MRC blocks

    if use_mpi and mpi4py.MPI.COMM_WORLD.Get_rank()==0 or not use_mpi:

        vols = []
        for start,stop in slices:
            vols.append(os.path.join(project.art_directory, "{0:s}_z_{1:.2f}.out".format(basename,start)))

        __stack__( vols, output, bounds3D=reconstruction )

    mpi4py.MPI.COMM_WORLD.Barrier()


def sirt( project, niter=25, mux=0.5, method = "a-SIRT", reset=False, nproc=20, use_mpi=False,
          apply_log=True, apply_inversion=True, **keywords ):
    """
    Simultaneous iterative method for TxBR.

    :param project: The project for which the iterative method will be run.
    :param niter: The number of iterations.
    :param mux: The mixing parameter for the corrections.
    :param method: The iterative method used can be *SIRT*, *a-SIRT* or *EM*.
        Default is *a-SIRT*.
    :param reset: The mixing parameter for the corrections.
    :param nproc: The number of CPUs that will be used.
    :param apply_inversion: If true, the original data is inverted. For a-sirt, the density 
        value should correlate to the material density; a dense entity such as a gold particle
        should correspond to high pixel values.
    """

    basename = project.basename

    if not method in ["SIRT","a-SIRT","EM","SMART"]:
        method = "a-SIRT"

    # Multiprocess definitions

    is_main_process = use_mpi and mpi4py.MPI.COMM_WORLD.Get_rank()==0 or not use_mpi
    nproc = use_mpi and mpi4py.MPI.COMM_WORLD.Get_size() or nproc
    def __barrier__():
        if use_mpi: mpi4py.MPI.COMM_WORLD.Barrier()
    def __print__( message ):
        if is_main_process: log.info(message)
        
    # Some settings for the original data stack
    
    apply_flatten = True
    apply_log = apply_log
    apply_inversion = apply_inversion
    apply_normalization = True

    __print__("Start algebraic reconstruction (method:{0:s}; log:{1:}) for {2:s}".format(method,apply_log,basename))

    if is_main_process and not os.path.lexists(project.art_directory):
        os.makedirs(project.art_directory)
            
    __barrier__()

    projection = project.stackedSeries().projection
    reconstruction = project.reconstruction
    
    use_fixed_segment_size = (projection.order > 1)

    # Some global parameters

    ntilt = projection.numberOfExposures()  # The number of exposures in the rawdata
    depth = reconstruction.end.z-reconstruction.origin.z + 1     # The number of slices in the volume

    # Partition the volume and micrograph space into smaller parts

    nt = numpy.ceil(float(ntilt)/nproc).astype('int')
    ns = numpy.ceil(float(depth)/nproc).astype('int')

    slices = []
    tilts = []

    for iproc in range(nproc):
        # Divide the tilts
        tilt_start,tilt_stop = iproc*nt,min((iproc+1)*nt,ntilt)
        if tilt_stop>=tilt_start: tilts.append((tilt_start,tilt_stop))
        # Divide the slices
        start,stop = reconstruction.origin.z+iproc*ns,min(reconstruction.origin.z+(iproc+1)*ns-1,reconstruction.end.z)
        if stop>=start: slices.append((start,stop))

    __print__("Tilts: {0:s}".format(tilts))
    __print__("Slices: {0:s}".format(slices))

    original_stack = os.path.abspath(os.path.join(project.art_directory, "{0:s}.st".format(basename)))
    src = os.path.abspath(os.path.join(project.art_directory, "{0:s}.project.st".format(basename)))
    dest = os.path.abspath(os.path.join(project.art_directory, "{0:s}.preali.rot.SL".format(basename)))
    res_vol = os.path.join(project.art_directory, "{0:s}_z_{1:.1f}.res.out".format(basename, reconstruction.origin.z))
    vol = os.path.join(project.art_directory, "{0:s}_z_{1:.1f}.out".format(basename,reconstruction.origin.z))

    skip = [] # We do not want to skip any slices
    #skip = [ tilt for tilt in range(ntilt) if tilt%10!=0]
    __print__("Skip Tilts: {0:s}".format(skip))


    volume_pattern = os.path.join(project.art_directory,"{0:s}_z_{1:.1f}.out.*".format(basename,reconstruction.origin.z))
    volumes = glob.glob(volume_pattern)

    index0 = len(volume_pattern)-1
    extensions = [volume[index0:] for volume in volumes]
    iter0 = 0
    for extension in extensions:
        try:
            iter0 = max(iter0, int(extension) + 1)
        except ValueError:
            pass

    if reset or not os.path.lexists(original_stack): # Reset the iterative process
    
        iter0 = 0   # Reset iteration counter
        
        if is_main_process:
            __prepare_input_data__( project, original_stack, apply_flatten, apply_log, apply_inversion, apply_normalization )

    if iter0!=0: # Re-project the last available volume

        __barrier__()

        __project__( project, projection, reconstruction, tilts, skip, use_fixed_segment_size, src, use_mpi )

        __barrier__()

    for index in range(iter0,iter0+niter):  # Main iterative process

        if is_main_process:

            log.info("Iteration: {0:}/{1:}".format(index+1,iter0+1+niter))

            if os.path.lexists(dest): os.remove(dest)

            if index==0:
                shutil.copy(original_stack,dest)
            else:
                shutil.copy(src,dest)
#                __susbtract__( original_stack, dest, alpha=ntilt, method=method)
                __susbtract__( original_stack, dest, method=method)

        __barrier__()

        __backproject__( project, projection, reconstruction, slices, skip, use_fixed_segment_size, res_vol, use_mpi )

        __barrier__()

        if is_main_process:

            if index==0:
                shutil.copy( res_vol, vol )
            else:
                # First handle backups of the volume
                old_backup_vol = os.path.join(project.art_directory, "{0:}-vol-{1:}.mrc".format(basename,index-2))
                new_backup_vol = os.path.join(project.art_directory, "{0:}-vol-{1:}.mrc".format(basename,index-1))
                if os.path.lexists(old_backup_vol): os.remove(old_backup_vol)
                shutil.copy( vol, new_backup_vol )
                # Correct for the volume
                __add__( res_vol, vol, mux=mux, method=method )

            old_link = vol + ".{0:}".format(index-1)
            new_link = vol + ".{0:}".format(index)

            if os.path.lexists(old_link): os.remove(old_link)
            if os.path.lexists(new_link): os.remove(new_link)

            os.symlink(os.path.basename(vol),new_link)

        __barrier__()

        __project__( project, projection, reconstruction, tilts, skip, use_fixed_segment_size, src, use_mpi )

        __barrier__()


def main():
    """
    The main routine for this iterative process.
    """

    try:

        flags1 = "hb:d:l:n:"
        flags2 = [ "help", "directory=", "basename=", "wd=", "bin=", "mux=", "no-adaptive", "EM", "SMART", "log", "reset", "nproc=", "mpi" ]

        opts, args = getopt.getopt( sys.argv[1:], flags1, flags2 )

    except getopt.GetoptError, err:

        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    keywords = {}

    keywords['directory'] = "."
    keywords['work_directory'] = "."
    keywords['basenames'] = None
    keywords['bin'] = None
    keywords['reset'] = False
    keywords['method'] = "a-SIRT"
    keywords['niter'] = 50
    keywords['mux'] = 2.0
    keywords['nproc'] = 20
    keywords['apply_inversion'] = True
    keywords['apply_log'] = False
    
    for option,value in opts:
        if option in ('-h','--help'):
            usage()
            sys.exit()
        elif option in ('-d','--directory'):
            keywords['directory'] = value
        elif option in ('--wd'):
            keywords['work_directory'] = value
        elif option in ('-b','--basename'):
            keywords['basenames'] = [str(s) for s in value.split(',')]
        elif option in ('-n'):
            keywords['niter'] = int(value)
        elif option in ('--mux'):
            keywords['mux'] = float(value)
        elif option in ('--bin'):
            keywords['bin'] = int(value)
        elif option in ('--no-adaptive'):
            keywords['method'] = "SIRT"
        elif option in ('--EM'):
            keywords['method'] = "EM"
        elif option in ('--SMART'):
            keywords['method'] = "SMART"
        elif option in ('--log'):
            keywords['apply_log'] = True
        elif option in ('--reset'):
            keywords['reset'] = True
        elif option in ('--nproc'):
            keywords['nproc'] = int(value)
        elif option in ('--mpi'):
            keywords['use_mpi'] = True
        else:
            assert False, "unhandled option"

    if keywords['basenames']==None:
        usage()
        sys.exit()
        
    # Create an instance of a TxBR project

    project = txbr.txbrdao.loadProject( directory=keywords['directory'], basenames=keywords['basenames'], work_directory=keywords['work_directory'], bin=keywords['bin'] )

    sirt( project, **keywords )

if __name__ == '__main__': main()