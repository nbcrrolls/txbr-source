#!/usr/bin/python

import sys
import os
import os.path
import getopt
import util
import psutil

from txbr import log

basename = None
directory = "."
work_directory = "."
bin = 1

EXTRA_ALIGN_FLAGS = ""
EXTRA_FILTER_FLAGS = ""
EXTRA_BCKPRJ_FLAGS = ""
EXTRA_FINAL_FLAGS = ""

ALIGN_CMD = ""
FILTER_CMD = ""
BCKPRJ_CMD = ""
FINALIZE_FILTER_CMD = ""
FINALIZE_BCKPRJ_CMD = ""

skips = []

def usage():

    print
    print 'Usage: %s.py -b basename1[,...] [Options]' %os.path.basename(sys.argv[0])
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files (.fid and .txbr)"
    print "    --wd work_directory (default to directory)"
    print "        Name of the directory where the output files are stored"
    print "    --ortho"
    print "        Flag for the orthogonal model"
    print "    --multiple (default value)"
    print "        Flag for the multiple model"
    print "    --full"
    print "        Do the full minimization at high order..."
    print "    --n1 n1"
    print "        Polynomial order of the approximation patches"
    print "    --n2 n2"
    print "        Polynomial order of the projaction map"
    print "    --n3 n3"
    print "        Polynomial order of the trace approximations"
    print "    --n4 n4"
    print "        Polynomial order of the contour approximations"
    print "    --flatten-order order"
    print "        Do a on the fly flattening from the bundle adjustment"
    print "    --fix-angles 1[,2,3] (only with option --ortho)"
    print "        Fix rotation angles in the orthogonal approximation"
    print "    --fix-t 1[,2,3]"
    print "        Fix the fixed points axis"
    print "    --axis u1,u2,u3 (default 0.0,1.0,0.0)"
    print "        Estimate for the rotation axis"
    print "    -z zstart,zstop"
    print "        Limit of the reconstruction between zstart and zstop in the beam direction"
    print "    --without index1[,index2,...]"
    print "        List of slices to skip during backprojection."
    print "    -l or --blocksize lz"
    print "        Set the size of a block in the z direction"
    print "    --remap"
    print "        Apply a polynomial 2D remap on the micrographs for filtering."
    print "    --remap-order order (default 1)"
    print "        Order of the polynomial 2D remap used during filtering."
    print "    --filter-mag value (default 2)"
    print "        Local magnification in the remap"
    print "    --no-contour"
    print "        Alignment is run without the contours"
    print "    --skip align[, filter, backprojection ]"
    print "        Run just the back-projection"
    print "    --mpi"
    print "        Run on a cluster of CPUs"
    print "    --gpu"
    print "        Run on a cluster of GPUs"
    print


def execute(cmd):
    '''Execute a command on the system'''
    
    global skips

    if 'align' in skips and cmd.find('align')!=-1: return
    if 'filter' in skips and cmd.find('filter')!=-1: return
    if 'backprojection' in skips and cmd.find('bckprj')!=-1: return
    
    log.info(cmd)
    
    code = os.system(cmd)
    
    if code!=0: sys.exit("An error happened for command: %s" %cmd)
    
    
def run_GPU_TxBR():
    
    nodes, numberOfGPUBoardsForNode, infoGPUBoard = util.infoForGPUHosts()   
    numberOfGPUBoards = util.numberOfGPUBoards()
    
    n = len(nodes)
    
    if n==0: run_TxBR()
    
    log.info('Run TxBR on GPUs')
    
    global ALIGN_CMD, FILTER_CMD, FINALIZE_FILTER_CMD, BCKPRJ_CMD, FINALIZE_BCKPRJ_CMD
    
    # Create Machine File for the different nodes and for the GPU
    
    node_file_name = os.path.join( work_directory, 'nodefile.txt' )
    gpu_file_name = os.path.join( work_directory, 'gpufile.txt' )
    
    node_file = open(node_file_name,"w")
    gpu_file = open(gpu_file_name,"w")
    
    for node in nodes:
        node_file.write('%s\n' %node)
        for gpu_index in range(numberOfGPUBoardsForNode[node]):
            gpu_file.write('%s\n' %node)
            
    node_file.close()
    gpu_file.close()
    
    # Alter the different commands if needed
    
    FILTER_CMD = FILTER_CMD + ' --mpi'
    FINALIZE_FILTER_CMD = FINALIZE_FILTER_CMD + ' --finalize-mpi %i' %len(nodes)
    
    BCKPRJ_CMD = BCKPRJ_CMD + ' --gpu'
    FINALIZE_BCKPRJ_CMD = FINALIZE_BCKPRJ_CMD + ' -n %i' %numberOfGPUBoards
    
   # Do the alignment and Start spawing the processes
    
    execute(ALIGN_CMD)
            
    execute('mpirun -machinefile %s -n %i %s' %( node_file_name, len(nodes), FILTER_CMD))
    
  #  execute('mpirun -machinefile %s -n 1 %s' %( node_file_name, FINALIZE_FILTER_CMD))
    
    execute('mpirun -machinefile %s -n %i %s' %( gpu_file_name, numberOfGPUBoards, BCKPRJ_CMD))

    execute('mpirun -machinefile %s -n 1 %s' %( node_file_name, FINALIZE_BCKPRJ_CMD))
    
    
def run_MPI_TxBR():
    
    global ALIGN_CMD, FILTER_CMD, BCKPRJ_CMD
    
    nodes, numberOfCPUForNode = util.infoForCPUHosts()
    numberOfCPUs = util.numberOfCPUs()
        
    # Create a Machine File
    
    filename = os.path.join( work_directory, 'machinefile.txt' )
    
    machine_file = open(filename,"w")
    for node in nodes:
        for i in range(numberOfCPUForNode[node]):
            machine_file.write('%s\n' %node)
    machine_file.close()
    
    # Alter the different commands if needed
    
   # FILTER_CMD = FILTER_CMD + ' --mpi'
    BCKPRJ_CMD = BCKPRJ_CMD + ' --mpi'
    
   # Do the alignment and Start spawing the processes
    
    execute(ALIGN_CMD)

    if len(nodes)<=1 and (nodes[0]=='localhost' or nodes[0]=='127.0.0.1'):

#        numberOfCPUs = int(0.8*psutil.NUM_CPUS)
                
        execute('mpiexec -n %i %s' %( numberOfCPUs, FILTER_CMD))
        execute('mpiexec -n %i %s' %( numberOfCPUs, BCKPRJ_CMD) )

    else:
    
        execute('mpiexec -f %s -n %i %s' %( filename, numberOfCPUs, FILTER_CMD))    
        execute('mpiexec -f %s -n %i %s' %( filename, numberOfCPUs, BCKPRJ_CMD) )
    

def run_TxBR():
    
    global ALIGN_CMD, FILTER_CMD, BCKPRJ_CMD, FINALIZE_BCKPRJ_CMD
    
    # Run the Alignment
    
    execute(ALIGN_CMD)
    
    # Run the Filter
    
    execute(FILTER_CMD)
    
    # Run the Back-Projection
    
    execute(BCKPRJ_CMD)
    
    # Finalize

   # execute(FINALIZE_BCKPRJ_CMD)
    

def main():

    flags1 = "hb:d:x:y:z:l:"
    flags2 = [ "help", "directory=", "wd=", "basename=", "preali", "bin=", "multiple", "ortho", "full", "no-contour", \
	           "n1=", "n2=", "n3=", "n4=", "flatten-order=", 'fix-angles=', "init", "axis=", "blocksize=", \
               "without=", "skip=", "remap", "remap-order=", "filter-mag=", "mpi","gpu", "doPlot" ]

    try:
        opts, args = getopt.getopt(sys.argv[1:], flags1, flags2)
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
        
    global basename, directory, work_directory, bin
    global EXTRA_ALIGN_FLAGS, EXTRA_FILTER_FLAGS, EXTRA_BCKPRJ_FLAGS, EXTRA_FINAL_FLAGS
    global skips
    
    MPI = False
    GPU = False
    
    XYZ_FLAGS = ''

    for option,value in opts:
        if option in ('-h','--help'):
            usage()
            sys.exit()
        elif option in ('-d','--directory'):
            directory = value
        elif option in ('--wd'):
            work_directory = value
        elif option in ('-b','--basename'):
            basename = value
            basenames = [str(s) for s in basename.split(',')]
        elif option in ("--bin"):
            bin = value
        elif option in ('--multiple','--ortho','--init','--full', '--no-contour', '--doPlot'):
            EXTRA_ALIGN_FLAGS += ' ' + option
        elif option in ('--n1', '--n2', '--n3', '--n4', "--flatten-order", '--fix-angles', '--fix-t', '--axis', '--remap-order'):
            EXTRA_ALIGN_FLAGS += ' ' + option + ' ' + value
        elif option in ("--filter-mag"):
            EXTRA_ALIGN_FLAGS += ' ' + option + ' ' + value
            EXTRA_FILTER_FLAGS += ' ' + option + ' ' + value
        elif option in ("--remap"):
            EXTRA_FILTER_FLAGS += ' ' + option
        elif option in ("-x","-y","-z"):
            Z = [float(s) for s in value.split(',')]
            XYZ_FLAGS += ' ' + '%s %f,%f,1' %(option,Z[0],Z[1])
        elif option in ('-preali','-l','--blocksize'):
            EXTRA_ALIGN_FLAGS += ' -l %i' %(int(value))
        elif option in ('--without'):
            EXTRA_BCKPRJ_FLAGS += ' --without %s' %value
        elif option in ('--skip'):
            skips = [str(s) for s in value.split(',')]
        elif option in ('--mpi'):
            MPI = True
        elif option in ('--gpu'):
            GPU = True
        else:
            assert False, "unhandled option"
            
    if basenames==None:
        usage()
        sys.exit()
        
    directory = os.path.abspath(directory)
    work_directory = os.path.abspath(work_directory)

    global ALIGN_CMD, FILTER_CMD, FINALIZE_FILTER_CMD, BCKPRJ_CMD, FINALIZE_BCKPRJ_CMD
    
    ALIGN_CMD = 'txbr_align.py -b %s -d %s --wd %s --bin %s%s' %( basename, directory, work_directory, bin, EXTRA_ALIGN_FLAGS )
    FILTER_CMD = 'txbr_filter.py -b %s -d %s --wd %s --bin %s --fine --doClean%s' %( basename, directory, work_directory, bin, EXTRA_FILTER_FLAGS )
    FINALIZE_FILTER_CMD = 'txbr_filter.py -b %s -d %s --wd %s --bin %s --doClean' %( basename, directory, work_directory, bin )
    BCKPRJ_CMD = 'txbr_bckprj.py -b %s -d %s --wd %s --bin %s %s %s' %( basename, directory, work_directory, bin, XYZ_FLAGS, EXTRA_BCKPRJ_FLAGS )
    FINALIZE_BCKPRJ_CMD = 'txbr_bckprj.py -b %s -d %s --wd %s --bin %s --finalize --doClean --doCombine %s %s' %( basename, directory, work_directory, bin, XYZ_FLAGS, EXTRA_BCKPRJ_FLAGS )

    if GPU: run_GPU_TxBR()
    elif MPI: run_MPI_TxBR()
    else: run_TxBR()


main()

