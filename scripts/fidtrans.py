#!/usr/bin/python

import sys
import os.path
import getopt

from txbr import log
from txbr.utilities import *

DEFAULT_DECIMATION = 2

CREATE_MARKER_MODEL_OPTION = "create"
DECIMATE_MODEL_OPTION = "decimate"
MODIFY_MARKER_MODEL_OPTION = "modify"
JOIN_MOSAIC_MODEL = "join-mosaic"

OPTIONS = [ CREATE_MARKER_MODEL_OPTION, MODIFY_MARKER_MODEL_OPTION, DECIMATE_MODEL_OPTION, JOIN_MOSAIC_MODEL ]

def usage():
    
    print
    print "Usage: %s option -b basename -d directory" % os.path.basename(sys.argv[0])
    print "          option in %s" %str(OPTIONS)
    print "    -n decimation (default value %i)" %DEFAULT_DECIMATION
    print "        Only keep 1 contour every n for each object of an IMOD model."
    print "    --scratch"
    print "        Ignore the .mosaic file"
    print "    --src preali[,st, ali]"
    print "        Type of the source"
    print "    --dest preali[,st, ali]"
    print "        Type of the destination"
    print "    --extension ext"
    print "        Extension of the model file to transform (.fid for a regular fiducial file)"
    print "    --common-only"
    print "        Only keep markers shared by at least two tiles in the %s option" %JOIN_MOSAIC_MODEL
    print
    print "Option \"%s\": Create a marker file from the fiducial file." %CREATE_MARKER_MODEL_OPTION
    print "   Each fiducial marker ends up being represented as one object as opposed to one"
    print "   contour (for contour alignment purpose)."
    print "Option \"%s\": Decimate the number of contours for each object of an" %DECIMATE_MODEL_OPTION
    print "    IMOD model. It is stricly applied on '.mrk' files."
    print "Option \"%s\": Allows to modify a model to match either preali, st or ali files." %MODIFY_MARKER_MODEL_OPTION
    print "    To go back and forth between preali and st, \".prexg\" file is needed."
    print "    See example below."
    print "Option \"%s\": Create marker files for a set of mosaic tiles. By default the" %JOIN_MOSAIC_MODEL
    print "    program will use a mosaic file to learn about the common marker between tiles."
    print "    With --scratch option, the program will connect contours sharing the same label in."
    print "    the fiducial files."
    print
    print "Example: - to change a fiducial file that matches a preali stack to match a st stack, type"
    print "    %s modify -b basename -[d directory] [--extension .fid] --src preali --dest st" %(os.path.basename(sys.argv[0]),)
    print "         - for the reverse if the fiducial file matching the st stack is called basename.fid [basename_st.fid]"
    print "    %s modify -b basename -[d directory] [--extension _st.fid] --src st --dest preali" %(os.path.basename(sys.argv[0]),)
    
    
def create_marker_model( directory, basenames ):
    """Create a marker model from fiducial files for a set of tomographic series."""
    
    for basename in basenames:
        
        log.info('Make marker model from (%s,%s)' %(directory,basename))
        
        makeMarkerModel(directory,basename)
        
        
def modify_marker_model( directory, basenames, src=PREALI_EXTENSION, dest=ST_EXTENSION, 
                         src_extension=FID_EXTENSION, dest_extension=FID_EXTENSION+"st" ):
    """Modify a model file to match either st, ali or preali files."""
    
    for basename in basenames:
        
        modelfile = os.path.join( directory, basename + src_extension )
        prexgPath = os.path.join( directory, basename + PREXG_EXTENSION )
        xgPath = os.path.join( directory, basename + XG_EXTENSION )
        
        log.info('Marker model from %s to be modified model from (%s,%s)' %(modelfile, src, dest))
        
        transfmod( modelfile, src=src, dest=dest, prexgPath=prexgPath, xgPath=xgPath)
        
        
def decimate_model( directory, basenames, decimation=DEFAULT_DECIMATION ):
    """Decimate an IMOD model by a factor specified in variable decimation."""

    basename, extensions = extract_series_name_on_basename_extension( basenames, directory )

    for extension in extensions:
        #model_path = os.path.join(directory, basename + extension + FID_EXTENSION)
        model_path = os.path.join(directory, basename + extension + MARKER_EXTENSION)
        log.info('Decimate %s by %i' %(model_path, decimation))
        decimateModel(model_path, decimation, type='object', crush=True)

        
def join_mosaic_model(directory, basenames, from_mosaic_file=True, keep_constraints_only=False):
    """Create marker files for a set of mosaic tiles."""
    
    if len(basenames)==0:
        return
    
    basename, extensions = extract_series_name_on_basename_extension( basenames, directory, extension=txbr.utilities.FID_EXTENSION )

    if len(extensions)==0:
        return

    if from_mosaic_file:
        mosaic_file = basename + MOSAIC_EXTENSION
        mosaic_file = os.path.join( directory, mosaic_file )
    else:
        mosaic_file = None
        
    log.info("Mosaic file: %s" %mosaic_file)
    
    mergeMosaicModels( basename, extensions, mosaic_file, keep_constraints_only=keep_constraints_only )


def main():
    
    flags1 = "hd:w:b:n:"
    flags2 = [ "help", "directory=", "base=", "extension=", "scratch", "src=", "dest=", "common-only" ]

    if len(sys.argv)==1 or sys.argv[1] not in OPTIONS:
        print
        print "You need to provide an option within %s for this command!" %OPTIONS
        usage()
        sys.exit(2)

    try:
        option = sys.argv[1]
        opts, args = getopt.getopt( sys.argv[2:], flags1, flags2 )
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    directory = '.'
    basenames = None
    n = DEFAULT_DECIMATION
    from_mosaic_file = True
    keep_constraints_only = False
    
    keywords = {}

    for option_,value in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-d", "--directory"):
            directory = value
        elif option_ in ("-b", "--basename"):
            basenames = value.split(',')
        elif option_ in ("--extension"):
            keywords['src_extension'] = value
        elif option_ in ("--src"):
            keywords['src'] = value
        elif option_ in ("--dest"):
            keywords['dest'] = value
        elif option_ in ("--scratch"):
            from_mosaic_file = False
        elif option_ in ("--common-only"):
            keep_constraints_only = True
        elif option_ in ("-n"):
            n = int(value)
        else:
            assert False, "unhandled option"

    if basenames==None:
        usage()
        raise SystemExit

    if option==CREATE_MARKER_MODEL_OPTION:
        create_marker_model( directory, basenames )

    if option==DECIMATE_MODEL_OPTION:
        decimate_model( directory, basenames, decimation=n )
        
    if option==MODIFY_MARKER_MODEL_OPTION:
        modify_marker_model( directory, basenames, **keywords )
            
    if option==JOIN_MOSAIC_MODEL:
        join_mosaic_model( directory, basenames, from_mosaic_file=from_mosaic_file, keep_constraints_only=keep_constraints_only )


main()
