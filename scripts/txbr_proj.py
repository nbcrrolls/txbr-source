#!/usr/bin/python

import sys
import os
import math
import getopt
import time
import numpy
import modl
import txbr.utilities

from txbr import log
from txbr import ALIGN_DIR
from txbrdao import loadProject


def usage():
    
    program = os.path.basename(os.path.basename(sys.argv[0]))

    print
    print "Usage: %s -b basename1[,...] [Options]" %(program)
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files (.fid and .txbr)"
    print "    --wd work_directory (default .)"
    print "        Name of the directory where output files will be stored"
    print "    -m model (basename.3dpkmod)"
    print "        Name of the model file that contains locations of the 3D markers"
    print "    --offset xoff,yoff,zoff"
    print "        Extra offset to add in the x,y and z directions"
    print "    --doPlot "
    print


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:],
            "hd:w:m:b:",
            [ "help", "directory=", "basenames=", "wd=", "offset=", "pts=", "doPlot" ]
            )
    except getopt.GetoptError, err:    # print help information and exit:
        print str(err) # print something like "option -a not recognized"
        usage()
        sys.exit(2)

    directory = None
    basename = None
    work_directory = None
    pts = None
    doPlot = False
    
    XYZ_filename = None
    extra_offset = [0.0,0.0,0.0]

    for option,value in opts:
        if option in ("-h", "--help"):
            usage()
            sys.exit()
        elif option in ("-d","--directory"):
            directory = value
        elif option in ("-b","--basenames"):
            basename = value
        elif option in ("--wd"):
            work_directory = value
        elif option in ("-m"):
            XYZ_filename = value
        elif option in ("--pts"):
            pts = []
            for v in value.split(","):
                rg = v.split("-")
                if len(rg)==2:
                    pts.extend(range(int(rg[0]),int(rg[1])+1))
                else:
                    pts.append(int(rg[0]))
        elif option in ("--offset"):
            extra_offset = [ float(v) for v in value.split(',') ]
        elif option in ("--doPlot"):
            doPlot = True
        else:
            assert False, "unhandled option"

    if basename==None:
        usage()
        sys.exit(2)

    if directory==None: directory = '.'

    if work_directory==None: work_directory = directory
        
    basenames = txbr.utilities.extract_series_name_from_dir( basename, directory )
    basename, extensions = txbr.utilities.extract_series_name( basenames )
    
    log.info("basename: %s" %basename)
    log.info("extensions: %s" %extensions)
        
    if XYZ_filename==None:
        XYZ_filename = os.path.join( work_directory, ALIGN_DIR, '%s.mod' %basename )
        
    if not os.path.lexists(XYZ_filename):
        log.warning('There is no 3D marker file (%s)!' %XYZ_filename)
        sys.exit(2)
        
    # First load the main project - Get the model offset
    
    project = loadProject(directory=directory,basenames=basename,work_directory=work_directory)

    rotation = project.reconstruction.getRotation()
    offset = project.reconstruction.getOffset()
    
    # Load the XYZ points
    
    log.info("Load 3D points from %s" %(XYZ_filename))
    log.info("Extra offsets %s" %(extra_offset))
    
    model = modl.Model(XYZ_filename)
    model.loadFromFile()

    XYZ = model.points()
    
    if pts!=None:
        pts = [ s-1 for s in pts if s>0 and s<=len(XYZ) ] 
        XYZ = XYZ[pts]

    XYZ[:,0] += - offset.x + extra_offset[0]
    XYZ[:,1] += - offset.y + extra_offset[1]
    XYZ[:,2] += - offset.z + extra_offset[2]

    XYZ = numpy.tensordot(XYZ,rotation,(1,0))    # Reorient the fiducials
    
    # For each series, do the projection 

    for series in project.series:

        ntilt = series.numberOfExposures()
        nx, ny = series.dimensions()
        sx, sy = series.getScales()
        sz = series.sz
        projmap = series.projection
      #  series.shiftPrealignTransform(eps=-1)  # Shift to st based pjmaps
        
        log.info('Series %s: %i exposures' %(series.basename, ntilt))

        st_name = os.path.join(directory,'%s.st' %(series.basename))
        preali_name = os.path.join(directory,'%s.preali' %(series.basename))
        template_name = os.path.join(directory,'%s.fid' %(series.basename))
        model_name = os.path.join(directory,'%s_test.mod' %(series.basename))

        log.info('Template file (%s): %s' %(template_name,os.path.lexists(template_name)))

        if os.path.lexists(template_name):
            template = modl.Model(template_name)
            template.loadFromFile()
            model = modl.Model(model_name,template=template)
        else:
            log.info('There is no template file (%s)!' %template_name)
            model = modl.Model(model_name)
            model.nx = nx
            model.ny = ny
            model.nz = ntilt
            model.xmax = nx
            model.ymax = ny
            model.zmax = ntilt
            model.imgref.xscale_new = sx
            model.imgref.yscale_new = sy
            model.imgref.zscale_new = sz/8.0

        o = model.addNewObject()

        for index,pt in enumerate(XYZ):
            contour = o.addNewContour()
            for itilt in range(ntilt):
                x = projmap.x(pt[0],pt[1],pt[2],itilt)
                y = projmap.y(pt[0],pt[1],pt[2],itilt)
                if x<=0 or x>nx or y<=0 or y>ny:
                    log.info("Point %i out of the box" %(index))
                else:
                    log.info("Point %i (%.2f,%.2f,%.2f) in series %s: (x,y)=(%.2f,%.2f,%.2f)" 
                             %(index, pt[0], pt[1], pt[2], series.basename, x, y, itilt))
                    contour.addPoint(x, y, itilt)

        model.save()
        
        # Show the volume/model
        
        if doPlot:
            os.system( 'imod %s %s' %( preali_name, model_name ) )
        

main()


