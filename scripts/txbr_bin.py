#!/usr/bin/python

import sys
import os
import os.path
import shutil
import getopt
import math
import numpy
import modl
import txbr.utilities

from txbr import log

DEFAULT_BINNING = 4

global xtranslation,ytranslation

def usage():

    print
    print 'Usage: %s.py -b basename1[,...] [Options]' %os.path.basename(sys.argv[0])
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files (.st, .preali, .fid, .prexg, .xg and .rawtlt)"
    print "    --wd work_directory (default to directory)"
    print "        Name of the directory where the output files (binnned) are stored"
    print "    --bin factor (default to directory)"
    print "        Binning factor"
    print


def execute(cmd):
    '''Execute a system command.'''
    
    log.info(cmd)
    
    code = os.system(cmd)

    if code!=0: sys.exit("An error happened for command: %s" %cmd)
    
    
def evaluateBinnedParameters( n, bin ):
    '''This routine calculates the new offset and size as calculated (hopefully) 
    by the newstack utility with a binning.
    Variable n: the original length (either the width or the height of a slice images)
    Variable bin: the binning factor
    Returns: the tuple (offset, n_bin). Varibale offset corresponds a translation
    term of the origin while n_bin is the new length.
    '''
    
    n_half = n/2.0
    n_bin = int(2*math.ceil(n_half/bin) - n%2)
    
    offset = int(math.ceil((bin*n_bin-n)/2))
    
    log.info("n=%i  bin=%i :  n_half=%i  offset=%i  nbin=%i" %(n, bin, n_half, offset, n_bin))
    
    return ( offset, n_bin )


def bin_down_mrc(directory, basenames, work_directory, extension, binning):
    '''Bin down some IMOD MRC files from directory to work_directory.'''
    
    for basename in basenames:
        src = os.path.join(directory, "%s%s" %(basename,extension))
        dest = os.path.join(work_directory, "%s%s" %(basename,extension))
        if not os.path.exists(src):
            log.warning('Cannot open file %s! Skip binning MRC file!' %src)
            break
        execute('newstack -bin %i %s %s' %(binning, src, dest))
        
        
def bin_down_model(directory, basenames, work_directory, extension, binning):
    '''Bin down some IMOD model files from directory to work_directory.'''
        
    global xtranslation,ytranslation
        
    for basename in basenames:
        
        log.info("Bin down (%s) model file for %s" %(extension, basename))
        
        src = os.path.join(directory, "%s%s" %(basename,extension))
        dest = os.path.join(work_directory, "%s%s" %(basename,extension))
        
        if not os.path.exists(src):
            log.warning('Cannot open file %s! Skip binning model file!' %src)
            break
        
        model = modl.Model(src)
        model.loadFromFile()
        model.filename = dest
        
        nx,ny = model.xmax,model.ymax
        
        offsetx,nbinx = evaluateBinnedParameters( nx, binning ) # new grid parametersfor x
        offsety,nbiny = evaluateBinnedParameters( ny, binning ) # new grid parametersfor y
        
        xtranslation = offsetx*model.imgref.xscale_new + model.imgref.xtranslation_new
        xscale = binning*model.imgref.xscale_new

        ytranslation = offsety*model.imgref.yscale_new + model.imgref.ytranslation_new
        yscale = binning*model.imgref.yscale_new

        ztranslation = model.imgref.ztranslation_new
        zscale = model.imgref.zscale_new
        
        model.xmax = nbinx
        model.ymax = nbiny
        
        log.info("%s -> Offsets: (%.2f,%.2f)    Binned Size: (%i,%i)" %(basename, xtranslation, ytranslation, nbinx, nbiny))
        
        model.updateImageReference( scale=(xscale, yscale, zscale), \
                                    translation=(xtranslation, ytranslation, ztranslation) )
        
        model.save()
        
        
def copyFile(directory, basenames, work_directory, extension):
    '''Copy some files (specified with a basename and an extension) from directory 
    to work_directory.'''

    for basename in basenames:
        
        log.info("Copy (%s) transformation file for %s" %(extension,basename))
        
        src = os.path.join(directory, "%s%s" %(basename,extension))
        dest = os.path.join(work_directory, "%s%s" %(basename,extension))
        
        shutil.copyfile(src,dest)
        
        
        
def bin_down_transformation(directory, basenames, work_directory, extension, binning):
    '''Bin down some transformation files from directory to work_directory. A corresponding 
    model must have been binned first so the global offsets (xtranslation,ytranslation) are
    known.
    '''

    global xtranslation,ytranslation
    
    for basename in basenames:
        
        log.info("Bin down (%s) transformation file for %s" %(extension,basename))

        src = os.path.join(directory, "%s%s" %(basename,extension))
        dest = os.path.join(work_directory, "%s%s" %(basename,extension))
        
        if not os.path.exists(src):
            log.warning('Cannot open file %s! Skip binning transformation file!' %src)
            break

        f_src = open(src)
        transform = numpy.array([line.split() for line in f_src],dtype='float')
        f_src.close
        
        transform[:,4] /= binning
        transform[:,5] /= binning
        
        f_dest = open(dest,'w')
        for row in transform:
            f_dest.write(" %11.7f %11.7f %11.7f %11.7f %11.3f %11.3f\n" %tuple(row))
        f_dest.close
        

def main():
    """The main function to perform the bin down"""
    
    directory = "."
    work_directory = "."
    basename = None

    binning = DEFAULT_BINNING

    flags1 = "hb:d:"
    flags2 = [ "help", "directory=", "wd=", "basename=", "bin=" ]

    try:
        opts, args = getopt.getopt(sys.argv[1:], flags1, flags2)
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

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
        elif option in ('--bin'):
            binning = int(value)
        else:
            assert False, "unhandled option"
            
    if basename==None:
        usage()
        sys.exit()
        
    if directory==work_directory:
        log.warning("The source and destination directory for this bin down operation should be different!")
        sys.exit()
        
    directory = os.path.abspath(directory)
    work_directory = os.path.abspath(work_directory)

    basenames = txbr.utilities.extract_series_name_from_dir( basename, directory, extension=txbr.utilities.FID_EXTENSION )
    
    # Bin down the MRC files
        
    bin_down_mrc(directory, basenames, work_directory, txbr.utilities.PREALI_EXTENSION, binning)
    bin_down_mrc(directory, basenames, work_directory, txbr.utilities.ST_EXTENSION, binning)
    
    # Copy the rawtlt file

    copyFile(directory, basenames, work_directory, txbr.utilities.RAWTLT_EXTENSION)
    
    
    # Bin down the fiducial model
    
    bin_down_model(directory, basenames, work_directory, txbr.utilities.FID_EXTENSION, binning)
        
    # Bin down the transformation files (it is tied to the first call of bin_down_model)
    
    bin_down_transformation(directory, basenames, work_directory, txbr.utilities.PREXG_EXTENSION, binning)
    

main()