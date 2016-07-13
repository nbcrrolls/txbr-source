#!/usr/bin/python

import os
import sys
import logging
import getopt
import numpy
import txbr

from txbr.align import contalign
from txbr.txbrdao import loadProject

def usage():

    print 'Usage: %s.py -b basename[,...]  [Options]' % os.path.basename(sys.argv[0])
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files (.fid and .txbr)"
    print "    --wd work_directory (default to directory)"
    print "        Name of the directory where the output files are stored"
    print "    --ortho"
    print "        Flag for the orthogonal model"
    print "    --preali"
    print "        Forces the alignment to use the preali file"
    print "    --multiple (default value for n2<=1)"
    print "        Flag for the multiple model"
    print "    --constant (default value for n2>1)"
    print "        Flag for the constant model, Override the multiple default flag"
    print "    --flatten-order order[,order_eff]"
    print "        Do a on the fly flattening from the bundle adjustment."
    print "        Variable order corresponds to the global polynomial function order that describe the bsurface boundaries defined by the 3D markers. "
    print "        Variable order_eff corresponds to the local polynomial function order used for the local projection map. "
    print "    --fix-angles 1[,2,3] (only with option --ortho)"
    print "        Fix rotation angles in the orthogonal approximation"
    print "    --fix-t 1[,2,3]"
    print "        Fix the fixed points axis"
    print "    --fix-series extension1[,exnension2,...]"
    print "        Fix the fixed points axis"
    print "    --axis u1,u2,u3 (default 0.0,1.0,0.0)"
    print "        Estimate for the rotation axis"
    print "    --n1 n1"
    print "        Polynomial order of the approximation patches"
    print "    --n2 n2"
    print "        Polynomial order of the projaction map"
    print "    --n3 n3"
    print "        Polynomial order of the trace approximations"
    print "    --n4 n4"
    print "        Polynomial order of the contour approximations"
    print "    --flatten-order order"
    print "        Performs a flattening for the reference tracks..."
    print "    --full"
    print "        Do the full minimization at high order..."
    print "    --load"
    print "        Load last minimization value (from .ctl files)..."
    print "    --skip-ncg"
    print "        Skip the conjugated gradient minimizaton"
    print "    --init"
    print "        Load last minimization value (from .txbr)..."
    print "    --init-0"
    print "        Load last minimization value, supposetly run independently for each series (from .txbr)..."
    print "    --no-contour"
    print "        Remove Contours from bundle adjustment..."
    print "    --tg"
    print "        Enforce the tangency constraints..."
    print "    -L"
    print "        Characteristic length of the structures..."
    print "    --structure index"
    print "        Allows to monitor an estimate for a structure..."
    print "    --remap-order order (default 1)"
    print "        Order of the 2D remap"
    print "    --filter-mag order (default 2)"
    print "        Local magnification in the remap"
    print "    --3d"
    print "        Do some 3d plots"
    print "    --decimate step"
    print "        Decimate the tilt series by considering angles every *step* angle"
    print "    --mosaic"
    print "        Asumption of an image-shift mosaic when generating the alignment transformation"
    print "    --debug"
    print "        Print more debug informations"
    print "    --doPlot"
    print "        Do some plottings"


def main():

    try:
        
        flags1 = "hb:d:l:"
        flags2 = [ "help", "directory=", "basename=", "wd=", "bin=","preali", "st", \
                  "ortho", "multiple", "constant", "full", "flatten-order=", \
                  "blocksize=", "fix-angles=", "fix-t=", "fix-series=", "axis=", \
                  'n1=', 'n2=', 'n3=', 'n4=', 'iter=', 'load', 'skip-ncg','init', 'init-0', \
                  'tg', 'T=', 'structure=', \
                  'no-contour', 'remap-order=', 'filter-mag=', 'decimate=','mosaic', '3d', 'doPlot',
                  'debug']

        opts, args = getopt.getopt( sys.argv[1:], flags1, flags2 )

    except getopt.GetoptError, err:

        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    directory = "."
    work_directory = "."
    bin = None

    basenames = None

    orthogonal = False
    multiple = False
    constant = False
    fix_angles = [ False, False, False ]
    fix_t = [ False, False, False ]
    
    structureToInitialize = None

    keywords = {}

    for option,value in opts:
        if option in ('-h','--help'):
            usage()
            sys.exit()
        elif option in ('-d','--directory'):
            directory = value
        elif option in ('--wd'):
            work_directory = value
        elif option in ('-b','--basename'):
            basenames = [str(s) for s in value.split(',')]
        elif option in ('--bin'):
            bin = int(value)
        elif option in ('--preali'):
            keywords['input'] = 'preali'
        elif option in ('--st'):
            keywords['input'] = 'st'
        elif option in ('--ortho'):
            orthogonal = True
        elif option in ('--multiple'):
            multiple = True
        elif option in ('--constant'):
            constant = True
        elif option in ('--fix-angles'):
            for s in value.split(','): fix_angles[int(s)-1] = True
        elif option in ('--fix-t'):
            for s in value.split(','): fix_t[int(s)-1] = True
        elif option in ('--fix-series'):
            keywords['fix_series_projection_map'] = [ s for s in value.split(',') ]
        elif option in ('--flatten-order'):
            #keywords['flatten_order'] = int(value)
            keywords['flatten_order'] = [ int(s) for s in value.split(",") ]
        elif option in ('-l','--blocksize'):
            keywords['blocksize'] = int(value)
        elif option in ('--n1'):
            keywords['n1'] = int(value)
        elif option in ('--n2'):
            keywords['n2'] = int(value)
        elif option in ('--n3'):
            keywords['n3'] = int(value)
        elif option in ('--n4'):
            keywords['n4'] = int(value)
        elif option in ('--iter'):
            keywords['iter'] = int(value)
        elif option in ('--full'):
            keywords['shortcut'] = False
        elif option in ('--skip-ncg'):
            keywords['skip_ncg'] = True
        elif option in ('--axis'):
            u = numpy.array([float(s) for s in value.split(',')])
            keywords['axis'] = u/numpy.sqrt(numpy.sum(u*u))
        elif option in ('--load'):
            keywords['evaluateFromScratch'] = False
        elif option in ('--init'):
            keywords['init'] = True
            keywords['applyRelativeTransformAtInit'] = False
        elif option in ('--init-0'):
            keywords['init'] = True
            keywords['applyRelativeTransformAtInit'] = True
        elif option in ('--tg'):
            keywords['enforceTangencyConstraint'] = True
        elif option in ('-L'):
            keywords['structureCharacteristicLength'] = float(value)
        elif option in ('--structure'):
            try:
                structureToInitialize = [int(s) for s in value.split(',')]
            except ValueError:
                structureToInitialize = "-"
        elif option in ('--no-contour'):
            keywords['with_contour'] = False
        elif option in ('--remap-order'):
            keywords['remap_order'] = int(value)
        elif option in ('--filter-mag'):
            keywords['remap_magnification'] = float(value)
        elif option in ('--decimate'):
            keywords['decimate'] = int(value)
        elif option in ('--mosaic'):
            keywords['mosaic'] = True
        elif option in ('--debug'):
            txbr.log.setLevel(logging.DEBUG)
        elif option in ('--3d'):
            keywords['do3dPlot'] = True
        elif option in ('--doPlot'):
            keywords['doPlot'] = True
        else:
            assert False, "unhandled option"

    if basenames==None:
        usage()
        sys.exit()
        
    try:
        n2 = keywords['n2']
        if n2<=1 and not constant: # Force a multiple model if it is projective
            multiple = True
    except KeyError:
        if not constant:
            multiple = True
            
    directory = os.path.abspath(directory)
    work_directory = os.path.abspath(work_directory)
    basename = ','.join(basenames)

    keywords['model'] = contalign.__encode_model__( orthogonal, multiple, fix_angles, fix_t )

    project = loadProject(directory=directory,basenames=basename,work_directory=work_directory,bin=bin)

    alignment = contalign.TxBRContourAlign( project, **keywords )
    alignment.process( structureToInitialize )

    if 'doPlot' in keywords.keys() and keywords['doPlot']:
        
        import pylab
        pylab.show()
        
    if 'do3dPlot' in keywords.keys() and keywords['do3dPlot']:
    
        from enthought.mayavi import mlab
        mlab.show()

if __name__ == '__main__': main()
