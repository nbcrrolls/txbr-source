import os.path
import re
import numpy
import logging
import logging.config
import utilities
import txbrconf

from txbr import TxBRproject
from txbr import Reconstruction
from txbr import TxBRSerie as TxBRSeries
from txbr import Map
from txbr import ProjectionMap

log = logging.getLogger('align')


def loadProject( project=None, directory=None, basenames=None, work_directory=None, input=None, bin=None, cfg_ext="" ):
    """
    Load/Reload a TxBR project.
    """

    log.info("Load project: (b='%s',d='%s',bin='%s')" %(basenames,directory,bin))

    if project==None:
        project = TxBRproject( directory, basenames, work_directory, input=input, bin=bin )

    if not os.path.lexists(project.directory) and not os.path.lexists(project.work_directory):
        return project

    # Load the info of the central TxBR reconstruction (in the main '.txbr' file' if it exists)

    loadReconstruction( project.reconstruction, bin=bin, cfg_ext=cfg_ext )

    # Load the possible series by scanning the different directories
    # A series may exist if (i) the data is available as a MRC stack (extension '.st') with header
    #                           information containing the metadata
    #                      (ii) the metadata is available ('.mdoc' file, '.txbr' file' or just '.rawtlt'
    #                           file with an asumption on the rotation axis)

    files = []

    if os.path.exists(project.directory):
        files.extend([ os.path.join(project.directory,file) for file in os.listdir(project.directory) if file.startswith(project.basename) ])
    if os.path.exists(project.work_directory) and project.directory!=project.work_directory:
        files.extend([ os.path.join(project.work_directory,file) for file in os.listdir(project.work_directory) if file.startswith(project.basename) ])

    files = [ os.path.abspath(file) for file in files ]
    
    for file in files:  # Scan for the basic core files

        st_match = re.match("%s(\S*)\\.%s$" %(os.path.abspath(os.path.join(project.directory,project.basename)),utilities.ST_EXTENSION[1:]), file)
        if st_match: ext = st_match.group(1)

        mdoc_match = re.match("%s(\S*)\\.%s$" %(os.path.abspath(os.path.join(project.directory,project.basename)),utilities.MDOC_EXTENSION[1:]), file)
        if mdoc_match: ext = mdoc_match.group(1)

        rwtlt_match = re.match("%s(\S*)\\.%s$" %(os.path.abspath(os.path.join(project.work_directory,project.basename)),utilities.RAWTLT_EXTENSION[1:]), file)
        if rwtlt_match: ext = rwtlt_match.group(1)

        txbr_match = re.match("%s(\S*)\\.%s$" %(os.path.abspath(os.path.join(project.work_directory,project.basename)),utilities.TXBR_EXTENSION[1:]), file)
        if txbr_match: ext = txbr_match.group(1)

        # A series within a project is available if its metadata is available

        if st_match:
            project.available_series[ext] = True
            continue
        if mdoc_match and not os.path.lexists(file[:-5]):   # Make sure the st was not already counted for
            project.available_series[ext] = True
            continue
        if rwtlt_match:
            project.available_series[ext] = True
            continue
        if txbr_match:
            project.available_series[ext] = True
            continue

    # Account for the possible series to add...

    series2Add = []

    if "" in project.extensions and len(project.available_series)>1:
        for extension in project.available_series.keys():
            if project.available_series[extension] and len(extension)>0:
                series2Add.append( project.basename+extension )
    else:
        for extension in project.extensions:
            if extension in project.available_series.keys() and project.available_series[extension]:
                series2Add.append( project.basename+extension )

    del project.extensions[:]

    for s in series2Add: # Series are first loaded as full resolution ones (bin 1)
        series = loadSeries( directory=directory, basename=s, work_directory=project.work_directory, bin=bin, cfg_ext=cfg_ext )
        series.rescaleForReconstruction( project.reconstruction )
        project.addSeries( series )

    # In the case the reconstruction was not loaded at first. Create it from same-scale series

    if project.reconstruction.isUndefined() and len(project.series)>0:
        # Redefine the reconstruction that contains all series
        project.reconstruction.scale.x,project.reconstruction.scale.y,project.reconstruction.scale.z \
                = numpy.max( [[series.sx,series.sy,series.sz] for series in project.series], axis=0 )
        origin = numpy.min( [[series.xmin*series.sx,series.ymin*series.sy,series.zmin*series.sz] for series in project.series], axis=0 )
        end = numpy.max( [[series.xmax*series.sx,series.ymax*series.sy,series.zmax*series.sz] for series in project.series], axis=0 )
        project.reconstruction.origin.x = origin[0]/project.reconstruction.scale.x
        project.reconstruction.origin.y = origin[1]/project.reconstruction.scale.y
        project.reconstruction.origin.z = origin[2]/project.reconstruction.scale.z
        project.reconstruction.end.x = end[0]/project.reconstruction.scale.x
        project.reconstruction.end.y = end[1]/project.reconstruction.scale.y
        project.reconstruction.end.z = end[2]/project.reconstruction.scale.z
        project.reconstruction.bin = bin

        if bin!=None and bin!=1:
            project.reconstruction.rescale( bin, bin, bin )

        # Rescale the series to the reconstruction
        for series in project.series: # isoscale the series
            series.rescaleForReconstruction( project.reconstruction )
        project.reconstruction.origin.x,project.reconstruction.origin.y,project.reconstruction.origin.z \
                = numpy.min( [[series.xmin,series.ymin,series.zmin] for series in project.series], axis=0 )
        project.reconstruction.end.x,project.reconstruction.end.y,project.reconstruction.end.z \
                = numpy.max( [[series.xmax,series.ymax,series.zmax] for series in project.series], axis=0 )


    bin = set([ s.stack.bin for s in project.series ])
   
    if len(bin)>0:
        project.setBinFactor(bin.pop())

    # Define the dimensions of the project (dimensions of the final volume)

    try:
        project.nx = int(project.reconstruction.end.x - project.reconstruction.origin.x)
        project.ny = int(project.reconstruction.end.y - project.reconstruction.origin.y)
    except ValueError:
        pass

    return project

    
def saveProject( project, directory=None, basename=None, work_directory=None ):
    """Save a TxBR project"""

    log.info("Save project")

    rep = ""

    if isinstance( project, TxBRproject ):
        rep += "basename: %s\n" %project.basename
        rep += "directory: %s\n\n" %project.directory
        rep += saveReconstruction( project.reconstruction, directory, basename, work_directory=work_directory )
        for s in project.series: rep += saveSeries( s, directory, basename, work_directory=work_directory )

    return rep
    

def loadReconstruction( reconstruction=None, directory=None, basename=None, work_directory=None, bin=None, cfg_ext=""  ):
    '''Load a reconstruction object from a configuration file'''

    if reconstruction==None or not isinstance(reconstruction,Reconstruction):
        reconstruction = Reconstruction( directory, basename, work_directory )

    wd = reconstruction.work_directory

    if bin==None: bin=1

    filename_bin = os.path.join( wd, reconstruction.basename + utilities.TXBR_EXTENSION + cfg_ext + ".bin%i" %bin)

    if not os.path.lexists(filename_bin):
        filename = filename_bin
    else:
        filename = os.path.join( wd, reconstruction.basename + utilities.TXBR_EXTENSION + cfg_ext )

    if not os.path.lexists(filename):
        return None

    log.info('Load Reconstruction from %s' %filename)

    if os.path.exists(filename):

        f = open(filename, 'r')

        for line in f:

            align_model_match = re.match('alignment: ', line)
            if align_model_match:
                index = len('alignment: ')
                reconstruction.alignment_model = line[index:].rstrip()

            input_match = re.match('input:\s*(\S+)', line)
            if input_match:
                reconstruction.input = input_match.group(1)

            xLimit_match = re.match('x->(\S*)', line)
            if xLimit_match:
                x_limits = xLimit_match.group(1).split(':')
                (reconstruction.origin.x,reconstruction.increment.x,reconstruction.end.x) = (float(x_limits[0]),float(x_limits[1]),float(x_limits[2]))

            yLimit_match = re.match('y->(\S*)', line)
            if yLimit_match:
                y_limits = yLimit_match.group(1).split(':')
                (reconstruction.origin.y,reconstruction.increment.y,reconstruction.end.y) = (float(y_limits[0]),float(y_limits[1]),float(y_limits[2]))

            zLimit_match = re.match('z->(\S*)', line)
            if zLimit_match:
                z_limits = zLimit_match.group(1).split(':')
                (reconstruction.origin.z,reconstruction.increment.z,reconstruction.end.z) = (float(z_limits[0]),float(z_limits[1]),float(z_limits[2]))

#                rotAxis_match = re.match('Rotation Axis:\s*\\[(\S*)\\]', line)
#                if rotAxis_match:
#                    rotAxis = rotAxis_match.group(1).split(',')
#                    reconstruction.rotAxis = [float(tk) for tk in rotAxis]

            sizeOfBlock_match = re.match('blocksize:\s*(\S*)', line)
            if sizeOfBlock_match:
                reconstruction.sizeOfBlock = int(sizeOfBlock_match.group(1))

            scale_match = re.match('scale:\s*(\S*)', line)
            if scale_match:
                scale_values = scale_match.group(1).split(':')
                (reconstruction.scale.x,reconstruction.scale.y,reconstruction.scale.z) = (float(scale_values[0]),float(scale_values[1]),float(scale_values[2]))

            bottomPlane_match = re.match('plane_coeffs1:\s*([\S+\s*]*)', line)
            if bottomPlane_match:
                plane = bottomPlane_match.group(1).split(' ')
                (reconstruction.bottomPlane.d,reconstruction.bottomPlane.a,reconstruction.bottomPlane.b) = (float(plane[0]),float(plane[1]),float(plane[2]))

            topPlane_match = re.match('plane_coeffs2:\s*([\S+\s*]*)', line)
            if topPlane_match:
                plane = topPlane_match.group(1).split(' ')
                (reconstruction.topPlane.d,reconstruction.topPlane.a,reconstruction.topPlane.b) = (float(plane[0]),float(plane[1]),float(plane[2]))

        f.close()

        reconstruction.update()

        log.info('Reconstruction loaded for series %s...' %filename)

    if bin!=None and bin!=1:

        reconstruction.bin = bin
        reconstruction.rescale( bin, bin, bin )

    return reconstruction



def saveReconstruction(  reconstruction,  directory=None, basename=None, work_directory=None  ):

    log.info("Save reconstruction")
    
    repr = ""

    if isinstance(reconstruction,Reconstruction):

        repr += 'x->%f:%f:%f\n' %(reconstruction.origin.x,reconstruction.increment.x,reconstruction.end.x)
        repr += 'y->%f:%f:%f\n' %(reconstruction.origin.y,reconstruction.increment.y,reconstruction.end.y)
        repr += 'z->%f:%f:%f\n' %(reconstruction.origin.z,reconstruction.increment.z,reconstruction.end.z)
        repr += '\n'
        repr += 'alignment: %s\n' %reconstruction.alignment_model
        repr += 'input: %s\n' %reconstruction.input
        repr += 'blocksize: %i\n' %(reconstruction.sizeOfBlock)
        repr += 'scale: %f:%f:%f\n' %(reconstruction.scale.x,reconstruction.scale.y,reconstruction.scale.z)
        repr += '\n'
        repr += 'plane_coeffs1: %f %f %f\n' %(reconstruction.bottomPlane.d,reconstruction.bottomPlane.a,reconstruction.bottomPlane.b)
        repr += 'plane_coeffs2: %f %f %f\n' %(reconstruction.topPlane.d,reconstruction.topPlane.a,reconstruction.topPlane.b)
        repr += '\n'

    return repr


def loadSeries( series=None, directory=None, basename=None, work_directory=None, bin=None, cfg_ext="" ):
    '''Load a series from a configuration file'''

    log.info("Load series: (b='%s',d='%s',bin='%s')" %(basename,directory,bin))

    if series==None:
        series = TxBRSeries( directory, basename, work_directory, bin=bin, load=False )

    # (i) Load the angle values [ from the header of stack, or the ".mdoc" file, or the ".rawtlt" file ]

    loadTiltAnglesForSeries( series )

    # (ii) Load the corresponding volume frame [ from a data stack (".st",".mrc",".preali") , or a model file (".fid",".mrk",".mrk")  ]

    loadFrameForSeries( series, cfg_ext )

    # (iii) Load the projection map [ from the ".txbr" file ]

    loadProjectionMapForSeries( series, cfg_ext )

    # (iv) Load the trajectory map if it exists [ from a ".traj" file ]

    loadTrajectoryMapForSeries( series, cfg_ext )

    # (v) Load the cross-correlation information (from the ".prexg" file)

    loadPrealignTransformations( series )

    # Extra post loading operations

    series.shiftPrealignTransform(eps=+1)

    series.validate()   # check conformity between the series info and its stack metadata

    return series


def saveSeries( series, directory=None, basename=None, work_directory=None ):

    log.info( "Save series" )

    rep = ""

    filename = os.path.join( work_directory, basename + utilities.TXBR_EXTENSION )

    if isinstance(series,TxBRSeries):
    
        rep +=  'Rotation Axis:[%f,%f,%f]\n\n' %(series.rotAxis[0],series.rotAxis[1],series.rotAxis[2])
        rep += saveMap(series.projection, directory, basename, work_directory=None)

    return rep


def loadTiltAnglesForSeries( series ):
    """Load the angle informations for a tomographic tilt series

    :param series: The series of interest.
    """

    if series==None: return

    log.info("Load tilt angles for series %s." %series)

    rawtlt_file = os.path.join( series.work_directory, series.basename + utilities.RAWTLT_EXTENSION )

    if series.stack!=None:
        series.tiltAngles = series.stack.getTiltAngles()
        log.info('%i tilt angles loaded for %s... from its stack' %(len(series.tiltAngles),series.basename))
    elif os.path.lexists(rawtlt_file):
        series.tiltAngles = []
        f = open(rawtlt_file, 'r')
        for line in f:  series.tiltAngles.append(float(line))
        f.close()
        log.info("%i tilt angles loaded for %s... from the '.rawtlt' file" %(len(series.tiltAngles),series.basename))
    else:
        log.warning('Unable to load tilt angles for %s!' %series.basename)


def loadFrameForSeries( series, cfg_ext="" ):
    """Load the boundary limits of the 3D volume corresponding to this series. To
    this purpose, there are three different options.
    First option, the configuration file ".txbr" is checked
    Second option, the metadata information of the stack is checked.
    Final option, the model file for the markers is checked as the last resource.

    :param series: The series of interest.
    """

    if series==None: return

    # The ".txbr" configuration file

    config_file = os.path.join( series.work_directory, series.basename + utilities.TXBR_EXTENSION + cfg_ext )

    if series.stack!=None:
        config_file_bin = os.path.join( series.work_directory, series.basename + utilities.TXBR_EXTENSION + cfg_ext + ".bin%i" %series.stack.bin )
        if os.path.lexists(config_file_bin): config_file = config_file_bin

    # The other files

    if series.stack!=None and (series.stack.bin==None or series.stack.bin==1):
        align_directory = os.path.join( series.work_directory, txbrconf.ALIGN_DIR )
    else:
        align_directory = os.path.join( series.work_directory, txbrconf.ALIGN_DIR, "bin%i" %series.stack.bin )

    # --- To fix, align directory will be only taken from the unbinned case for now --
    align_directory = os.path.join( series.work_directory, txbrconf.ALIGN_DIR )
    # --------------------------------------------------------------------------------

    marker_file = os.path.join( align_directory, series.basename + utilities.FID_EXTENSION )
    if not os.path.lexists(marker_file):
        marker_file = os.path.join( align_directory, series.basename + utilities.MARKER_EXTENSION )
    if not os.path.lexists(marker_file):
        marker_file = os.path.join( align_directory, series.basename + utilities.TRACK_EXTENSION )

    if os.path.lexists(config_file):

        f = open(config_file, 'r')

        for line in f:

            scale_match = re.match('scale:\s*(\S*)', line)
            if scale_match:
                scale_values = scale_match.group(1).split(':')
                (series.sx,series.sy,series.sz) = (float(scale_values[0]),float(scale_values[1]),float(scale_values[2]))

            xLimit_match = re.match('x->(\S*)', line)
            if xLimit_match:
                x_limits = xLimit_match.group(1).split(':')
                series.xmin = float(x_limits[0])
                series.xmax = float(x_limits[2])

            yLimit_match = re.match('y->(\S*)', line)
            if yLimit_match:
                y_limits = yLimit_match.group(1).split(':')
                series.ymin = float(y_limits[0])
                series.ymax = float(y_limits[2])

            zLimit_match = re.match('z->(\S*)', line)
            if zLimit_match:
                z_limits = zLimit_match.group(1).split(':')
                series.zmin = float(z_limits[0])
                series.zmax = float(z_limits[2])

            rotAxis_match = re.match('Rotation Axis:\s*\\[(\S*)\\]', line)
            if rotAxis_match:
                rotAxis = rotAxis_match.group(1).split(',')
                series.rotAxis = [float(tk) for tk in rotAxis]

        f.close()

    elif series.stack!=None:

        series.xmin,series.ymin = 1,1
        series.xmax,series.ymax = series.stack.getImageSize()
        series.zmax = min(series.xmax,series.ymax)/10.0
        series.zmin = - series.zmax

        series.sx,series.sy = series.stack.getImageScale()
        series.sz = max(series.sx,series.sy)

    elif os.path.lexists(marker_file):

        (series.sx,series.sy,series.sz) = modl.getScales(marker_file)
        (series.xmax,series.ymax,series.zmax) = modl.getDimensions(marker_file)

        series.xmin = 1
        series.ymin = 1
        series.zmin = -series.zmax


def loadPrealignTransformations( series ):
    """Load the prealignment information for a tomographic tilt series.

    :param series: The series of interest.
    """

    if series==None: return

    if series.stack!=None and (series.stack.bin==None or series.stack.bin==1):
        align_directory = os.path.join( series.work_directory, txbrconf.ALIGN_DIR )
    else:
        align_directory = os.path.join( series.work_directory, txbrconf.ALIGN_DIR, "bin%i" %series.stack.bin )

    # --- To fix, align directory will be only taken from the unbinned case for now --
    align_directory = os.path.join( series.work_directory, txbrconf.ALIGN_DIR )
    # --------------------------------------------------------------------------------

    preali_file = os.path.join( align_directory, series.basename + utilities.PREXG_EXTENSION )

    if os.path.exists(preali_file):
        f = open(file, 'r')
        preali = numpy.array([line.split() for line in f],dtype='float')
        series.prealignTranslations = preali[:,4:]
    else:
        series.prealignTranslations = numpy.zeros((series.numberOfExposures(),2))

    return series.prealignTranslations



def loadProjectionMapForSeries( series, cfg_ext="" ):
    """Load the projection map for a tomographic tilt series.

    :param series: The series of interest.
    """

    if series==None: return
    
    config_file = os.path.join( series.work_directory, series.basename + utilities.TXBR_EXTENSION + cfg_ext )

    if series.stack!=None:
        config_file_bin = os.path.join( series.work_directory, series.basename + utilities.TXBR_EXTENSION + cfg_ext + ".bin%i" %series.stack.bin)
        if os.path.lexists(config_file_bin): config_file = config_file_bin

    series.projection = loadMap( config_file )
    
    if series.projection==None:
        series.projection = ProjectionMap( approximationOrder=1, n=len(series.tiltAngles) )
        
    return series.projection


def loadTrajectoryMapForSeries( series, cfg_ext ):
    """Load the trajectory map involved in a tomographic tilt series.

    :param series: The series of interest.
    """

    if series==None: return

    config_file = os.path.join( series.work_directory, series.basename + utilities.TRAJ_EXTENSION + cfg_ext )

    series.trajectory = loadMap( config_file, type==utilities.TRAJ_EXTENSION )

    return series.trajectory


def loadMap( filepath, type=None ):
    '''Load a map from a configuration file.

    :param filepath: The path for the configuration file containing the projection map information.
    '''

    log.info( "Loading map from %s." %filepath )

    if not os.path.lexists(filepath):
        log.warning('Unable to load maps from %s!' %filepath)
        return None

    if type=="utilities.TRAJ_EXTENSION":
        map = TrajectoryMap()
    else:
        map = ProjectionMap()

    f = open(filepath, 'r')

    for line in f:

        scale_match = re.match('scale:\s*(\S*)', line)
        if scale_match:
            scale_values = scale_match.group(1).split(':')
            sx,sy,sz = (float(scale_values[0]),float(scale_values[1]),float(scale_values[2]))
            map.sx = sx
            map.sy = sy
        order_match = re.match('[\S,\W]*order:\W*(\S*)', line)
        if order_match:
            order = int(order_match.group(1))
            map.approximationOrder(approximationOrder=order)
        x_match = re.match('x-(\d+):\s*([\S+\s*]*)', line)
        if x_match:
            x_list = x_match.group(2).split()
            map.x_coefficients = numpy.row_stack([map.x_coefficients,numpy.array([float(s) for s in x_list])])
        y_match = re.match('y-(\d+):\s*([\S+\s*]*)', line)
        if y_match:
            y_list = y_match.group(2).split()
            map.y_coefficients = numpy.row_stack([map.y_coefficients,numpy.array([float(s) for s in y_list])])
        scaling_match = re.match('lambda-(\d+):\s*([\S+\s*]*)', line)
        if scaling_match:
            scaling_list = scaling_match.group(2).split()
            map.scaling_coefficients = numpy.row_stack([map.scaling_coefficients,numpy.array([float(s) for s in scaling_list])])

    f.close()

    log.info('Projection Map loaded for %s...' %filepath)

    return map


def saveMap( map, directory=None, basename=None, work_directory=None ):

    log.info("Save map")

    repr = ""

    if isinstance(map,Map):

        repr += 'number of tilts: %i\n\n' %map.numberOfExposures()
        repr += 'projection approximation order: %i\n\n' %map.approximationOrder()
        
        format = ''.join([' %e' for i in range(map.numberOfTerms)]) + '\n'
        for i in range(map.numberOfExposures()):
            scalingformat = 'lambda-%i:' + format
            repr += scalingformat %tuple([i] + map.scaling_coefficients[i,:].tolist())
            xformat = 'x-%i:' + format
            repr += xformat %tuple([i] + map.x_coefficients[i,:].tolist())
            yformat = 'y-%i:' + format
            repr += yformat %tuple([i] + map.y_coefficients[i,:].tolist())

    return repr
