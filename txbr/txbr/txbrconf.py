import os.path
import logging
import logging.config
import ConfigParser

DATA_DIR = 'data'
DETECT_DIR = 'txbr-detect'
ALIGN_DIR = 'txbr-align'
MULTIPLE_SERIES_DIR = 'multpl-series'
FILTER_DIR = 'txbr-filter'
SETUP_DIR = 'txbr-setup'
BACKPROJECTION_DIR = 'txbr-backprojection'
CONTOUR_DIR = 'contour'
VOL_DIR = "txbr-reco"
ART_DIR = "txbr-art"
PUB_DIR = "txbr-pub"
WEB_DIR = "web"

# Log settings

LOG_CONF = os.path.join(os.path.dirname(__file__),'log.conf')
logging.config.fileConfig(LOG_CONF)

log = logging.getLogger('align')

# Load configuration

global VERBOSE, \
       INPUT, \
       USE_MAYAVI, \
       ENFORCE_TANGENCY_CONSTRAINT, \
       STRUCTURE_CHARACTERISTIC_LENGTH, \
       MAXIMUM_NUMBER_ITERATION, \
       TOL, \
       X_PADDING, \
       Y_PADDING, \
       Z_PADDING, \
       REMAP_ORDER, \
       FILTER_MAGNIFICATION, \
       X_TAPERING_RATIO, \
       Y_TAPERING_RATIO

# Versioning

def version():
    '''Read the version of the TxBR module'''

    try:
        version_file = os.path.join('..','..','data','VERSION')
        f = open(version_file)
        return f.readline().rstrip('\n')
    except:
        return 'beta'

    return None


def loadConfig( work_directory="." ):

    global VERBOSE, INPUT, USE_MAYAVI, ENFORCE_TANGENCY_CONSTRAINT, STRUCTURE_CHARACTERISTIC_LENGTH, \
       MAXIMUM_NUMBER_ITERATION, TOL, X_PADDING, Y_PADDING, Z_PADDING, SAVE_PLOT, REMAP_ORDER, \
       FILTER_MAGNIFICATION, X_TAPERING_RATIO, Y_TAPERING_RATIO

    current_dir_file_config = os.path.join(os.path.dirname(work_directory),'txbr.cfg')

    if os.path.lexists(current_dir_file_config):
        file_config = current_dir_file_config
    else:
        file_config = os.path.join(os.path.dirname(__file__),'txbr.cfg')

    config = ConfigParser.RawConfigParser()
    config.read(file_config)

    try:
        VERBOSE = config.getboolean('Configuration', 'VERBOSE')
    except:
        log.warning("Default configuration variable: VERBOSE=False")
        VERBOSE = False

    try:
        USE_MAYAVI = config.getboolean('Configuration', 'USE-MAYAVI')
    except:
        log.warning("Default configuration variable: USE_MAYAVI = True")
        USE_MAYAVI = True

    try:
        SAVE_PLOT = config.getboolean('Configuration', 'SAVE-PLOT')
    except:
        log.warning("Default configuration variable: SAVE_PLOT = True")
        SAVE_PLOT = True

    try:
        INPUT = config.get('Align', 'INPUT')
    except:
        log.warning("Default configuration variable: INPUT = 'st'")
        INPUT = 'st'

    try:
        ENFORCE_TANGENCY_CONSTRAINT = config.getboolean('Align', 'ENFORCE-TANGENCY-CONSTRAINT')
    except:
        log.warning("Default configuration variable: ENFORCE_TANGENCY_CONSTRAINT = True")
        ENFORCE_TANGENCY_CONSTRAINT = True

    try:
        STRUCTURE_CHARACTERISTIC_LENGTH = config.getfloat('Align', 'STRUCTURE-CHARACTERISTIC-LENGTH')
    except:
        log.warning("Default configuration variable: STRUCTURE_CHARACTERISTIC_LENGTH = 100")
        STRUCTURE_CHARACTERISTIC_LENGTH = 100

    try:
        MAXIMUM_NUMBER_ITERATION = config.getint('Align', 'MAXIMUM-NUMBER-ITERATION')
    except:
        log.warning("Default configuration variable: MAXIMUM_NUMBER_ITERATION = 0")
        MAXIMUM_NUMBER_ITERATION = 0

    try:
        TOL = config.getfloat('Align', 'TOLERANCE')
    except:
        log.warning("Default configuration variable: TOLERANCE = 1e-16")
        TOL = 1e-16

    try:
        REMAP_ORDER = config.getint('Filter', 'REMAP-ORDER')
    except:
        log.warning("Default configuration variable: REMAP_ORDER = 1")
        REMAP_ORDER = 1

    try:
        FILTER_MAGNIFICATION = config.getfloat('Filter', 'FILTER-MAGNIFICATION')
    except:
        log.warning("Default configuration variable: FILTER_MAGNIFICATION = 2.0")
        FILTER_MAGNIFICATION = 2.0

    try:
        X_PADDING = config.getint('Reconstruction', 'X-PADDING')
    except:
        log.warning("Default configuration variable: X_PADDING = 0")
        X_PADDING = 0

    try:
        Y_PADDING = config.getint('Reconstruction', 'Y-PADDING')
    except:
        log.warning("Default configuration variable: Y_PADDING = 0")
        Y_PADDING = 0

    try:
        Z_PADDING = config.getint('Reconstruction', 'Z-PADDING')
    except:
        log.warning("Default configuration variable: Z_PADDING = 10")
        Z_PADDING = 10

    try:
        X_TAPERING_RATIO = config.getfloat('Merging', 'X-TAPERING-RATIO')
    except:
        log.warning("Default configuration variable: X_TAPERING_RATIO = 0.05")
        X_TAPERING_RATIO = 0.05

    try:
        Y_TAPERING_RATIO = config.getfloat('Merging', 'Y-TAPERING-RATIO')
    except:
        log.warning("Default configuration variable: Y_TAPERING_RATIO = 0.05")
        Y_TAPERING_RATIO = 0.05


loadConfig()

if VERBOSE: log.setLevel(logging.DEBUG)

log.debug('Verbose: %s' %VERBOSE)
log.debug('Use Mayavi: %s' %USE_MAYAVI)
log.debug('Save Plots: %s' %SAVE_PLOT)
log.debug('Input: %s' %INPUT)
log.debug('Enforce Tangency Constraint: %s' %ENFORCE_TANGENCY_CONSTRAINT)
log.debug('Characteristic length of the structures: %s' %STRUCTURE_CHARACTERISTIC_LENGTH)
log.debug('Maximum Number of iteration: %i' %MAXIMUM_NUMBER_ITERATION)
log.debug('Remap order: %i' %REMAP_ORDER)
log.debug('Filter Magnification: %f' %FILTER_MAGNIFICATION)
log.debug('Padding in the (X,Y,Z) direction: (%i,%i,%i)' %(X_PADDING,Y_PADDING,Z_PADDING))
log.debug('Tapering ratio in the (X,Y) direction: (%.2f,%.2f)' %(X_TAPERING_RATIO,Y_TAPERING_RATIO))