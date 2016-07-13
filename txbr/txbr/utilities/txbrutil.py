import common
import modl
import numpy
import os
import os.path
import re
import util
import logging.config

log = logging.getLogger()

global CONFIG_FILE

CONFIG_FILE= "config.fly"

def loadRotationAxis( scope ):
    """Load Rotation Axis from a configuration file"""

    if scope==None:
        config_file = CONFIG_FILE
    else:
        config_file = "config.%s.fly" %(scope.lower())

    data_dir =  os.path.join(os.path.dirname(__file__),'..','..','..','data')
    org_config_file = os.path.join(data_dir,'config',config_file)

    print org_config_file

    if not os.path.exists(org_config_file):
        org_config_file = os.path.join(data_dir,'config',CONFIG_FILE)

    if os.path.exists(CONFIG_FILE): # Search in the current directory first
        config_file = CONFIG_FILE
    else:
        config_file = org_config_file

    try:

        f = open(org_config_file, 'r')
        for line in f:
            rotAxis_match = re.match('Rotation Axis\s*:\s*\\[(\S*)\\]', line)
            if rotAxis_match:
                rotAxis = rotAxis_match.group(1).split(',')
                rot_axis = [ float(tk) for tk in rotAxis ]
        f.close()

    except IOError:

        print "Error loading the Configuration file!"
        
    print "Scope: %s   Rotation Axis: [%.2f,%.2f,%.2f]" %( scope, rot_axis[0], rot_axis[1], rot_axis[2] )

    return rot_axis


def load_transformations(file):
    """This routine load 2D transformations from an IMOD transformation
    file (such as .xf or .xg)."""
    
    transforms = []    
    
    if file==None or not os.path.exists(file):
        log.warning('File %s does not exist!' %file)
        return transforms
        
    f = open(file)

    try:
        for line in f:
            items = line.rsplit()
            args = [float(item) for item in items]
            M = numpy.resize(args[:4],(2,2))
            T = numpy.resize(args[4:],(2,1))
            transformation = util.AffineTransformation2D(numpy.column_stack((T,M)))
            transforms.append(transformation)
            log.info(transformation) 
    finally:
        f.close()

    return transforms


def zero_out_transformations( file ):
    """This routine zero out a transformation file half way through the file."""

    transforms = numpy.loadtxt(file)
    i0 = len(transforms)/2
    transforms[:,-2:] = transforms[:,-2:] - transforms[i0,-2:]

    numpy.savetxt( file, transforms, fmt='%12.3f', delimiter='' )


def transfmod( inputPath, outputPath=None, src=common.PREALI_EXTENSION, 
               dest=common.ALI_EXTENSION, prexgPath=None, xgPath=None):
    """Transform a model from a transformed state to another transformed state. Three
    states are possible in IMOD (state='.st','.preali','.ali' or equivalently 'st','preali',
    'ali'). The variable inputPath points to the model file (.fid). One at minimum of the 
    two variables prexgPath and xgPath are needed depending on the source and destination 
    of the transformation. Model is then saved in an appropriate file.
    """
    
    if not src.startswith("."): src = "." + src
    if not dest.startswith("."): dest = "." + dest
    
    log.info("Transform model %s from %s to %s" %(inputPath,src,dest))

    directory, filename = os.path.split(inputPath)
    basename, ext = os.path.splitext(filename)

    if outputPath==None:
        outputPath = os.path.join(directory, basename + '_' + dest + ext)

    model = modl.Model(outputPath)
    model.loadFromFile(inputPath)

    prexgTransforms = load_transformations(prexgPath)
    xgTransforms = load_transformations(xgPath)

    if len(prexgTransforms)!=len(xgTransforms):
        log.warning('prexg and xg files do not have the same number of transformantions')

    n = max(len(prexgTransforms),len(xgTransforms))

    transforms = []

    if src==common.PREALI_EXTENSION and dest==common.ALI_EXTENSION:
        transforms = [ xgTransforms[index].compose(prexgTransforms[index].inv()) for index in range(n) ]
    elif src==common.PREALI_EXTENSION and dest==common.ST_EXTENSION:
        transforms = [ transform.inv() for transform in prexgTransforms ]
    elif src==common.ST_EXTENSION and dest==common.ALI_EXTENSION:
        transforms = [ transform for transform in xgTransforms ]
    elif src==common.ST_EXTENSION and dest==common.PREALI_EXTENSION:
        transforms = [ transform for transform in prexgTransforms ]
    elif src==common.ALI_EXTENSION and dest==common.ST_EXTENSION:
        transforms = [ transform.inv() for transform in xgTransforms ]
    elif src==common.ALI_EXTENSION and dest==common.PREALI_EXTENSION:
        transforms = [ prexgTransforms[index].compose(xgTransforms[index].inv()) for index in range(n) ]
        
    if transforms==None or len(transforms)==0:
        log.warning("No Transformations from %s to %s" %(src, dest))
        
    for obj in model.objects:
        for contour in obj.contours:
            for point in contour.points:
                index = int(round(point[2],0))
                point_ = transforms[index].forward(point[:2])
                point[:2] = point_[:]

    log.debug(model.info())

    model.save()

    return model


def makeMarkerModel( directory, basename, work_directory=None ):
    """Transform an IMOD fiducial file to a model where every fiducial tracks
    is an object

    :param directory: The directory containing the fiducial model.
    :param basename: The basename of the fiducial model.
    :param work_directory: The directory where the marker model will be stored in. If
        None, work_directory will be the same as *directory*.
    """

    log.info('Create a marker model for %s' %basename)
    
    if work_directory==None:
        work_directory = directory

    fidPath = os.path.join( directory, basename + '.fid' )
    mrkPath = os.path.join( work_directory, basename + '.mrk' )
    trkPath = os.path.join( work_directory, basename + '.trk' )

    if not os.path.exists(fidPath):
        log.info('File %s does not exist!' %fidPath)
        return None

    model = modl.Model(mrkPath)
    model.loadFromFile(fidPath,keepOnlyPointStructures=False) # Fiducial File for Etomo has a different structure

    object = model.objects[0]

    numberOfTracks = object.numberOfContours() # Number of tracks
    numberOfTilts = model.zmax # Number of tilts

    print "NumberOfTilts=%i" % numberOfTilts

    markers = numpy.zeros((numberOfTracks,numberOfTilts,2),dtype=float)    # markers[tilt #, marker #, X, Y]
    markers[:,:,:] = numpy.nan

    for itrack,contour in enumerate(object.contours):
        if numberOfTilts!=contour.npoints():
            log.info('Track #%i : missing %i points' %(itrack,numberOfTilts-contour.npoints()))
        log.info('track=%i    %i' %(itrack,len(contour.points)))
        for point in contour.points:
            itilt = int(round(point[2],0))
            #print " %s  %i   track: %i   tilt: %i    %f" %(basename, contour.indexOfContour, itrack,itilt,point[2])
            markers[itrack,itilt,0] = point[0]
            markers[itrack,itilt,1] = point[1]

    model.removeObject(0)

    # Save the marker model

    for itrack in range(numberOfTracks):
        if numpy.all(numpy.isnan(markers[itrack,:,:])):    # do not consider empty tracks
            # Add an empty objects
            object = model.addNewObject()
            continue
        object = model.addNewObject()
        for itilt in range(numberOfTilts):
            if numpy.any(numpy.isnan(markers[itrack,itilt,:])):
                continue
            x = markers[itrack,itilt,0]
            y = markers[itrack,itilt,1]
            z = itilt
            contour = object.addNewContour()
            contour.addPoint(x,y,z)

    model.save()

    if os.path.lexists(trkPath):
        os.unlink(trkPath)

    os.symlink(os.path.basename(mrkPath),trkPath)

    return model


def decimateModel( fidPath, n, type='contour', crush=False ):
    """Decimate the number of ojects of contours in a given model, depending
    of the value of variable type.
    The Result is stored in a new file"""
    
    if not os.path.exists(fidPath):
        log.warning("The model file %s does not exist!" %fidPath)
        return

    step_obj = 1
    step_cont = 1

    if type=='object': step_obj = n
    elif type=='contour': step_cont = n

    (directory,filename) = os.path.split(fidPath)
    (basename,ext) = os.path.splitext(filename)

    if crush:
        newPath = fidPath
    else:
        newPath = os.path.join(directory,basename + '_dec%i' %n + ext)

    model = modl.Model(newPath)
    model.loadFromFile(fidPath)

    model.objects = model.objects[::step_obj]

    for object in model.objects:
        object.contours = object.contours[::step_cont]

    model.save()

