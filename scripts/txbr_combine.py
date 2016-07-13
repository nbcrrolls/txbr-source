#!/usr/bin/python

import sys
import os.path
import getopt
import txbr
import txbr.join
import txbr.utilities

from txbr import log

CROSS_VALIDATION_OPTION = "crossval"
MERGE_OPTION = "combine"
JOIN_OPTION = "join"
JOIN_CHECK_OPTION = "check-join"
JOIN_MICROGRAPHS = "join-2D"
WEIGTH_VOLUME = "weight-3D"

def usage():

    print "Usage: %s option -i inputs" % sys.argv[0]
    print "    -z zstart,zstop"
    print "        Limit of the combining between zstart and zstop in the beam direction"



def execute(cmd):
    '''Execute a command on the system'''
    
    log.info(cmd)
    
    code = os.system(cmd)
    
    if code!=0: sys.exit("An error happened for command: %s" %cmd)
    

def combineMRCFiles( directory, inputs, output):

    command = os.path.join(os.environ['TXBR'],'util','icombine')
    
    args = "-i " + " 1 ".join(inputs) + " 1" + " -d " + directory + " -o " + output;

    os.system(command + ' ' + args)
    

def crossvalidate(directory,inputs,output):

    command = os.path.join(os.environ['TXBR'],'util','crossval')
    
    args = "-rm -i " + " ".join(inputs) + " -d " + directory + " -o " + output;

    os.system(command + ' ' + args)


def main():
    
    valid_options = [ CROSS_VALIDATION_OPTION, 
                      MERGE_OPTION, 
                      JOIN_OPTION, 
                      JOIN_CHECK_OPTION, 
                      JOIN_MICROGRAPHS,
                      WEIGTH_VOLUME ]
    
    flags1 = "d:b:i:o:z:h:s"
    flags2 = [ "help", "directory=", "basename=", "input=", "output=", "bin=", "test" ]

    try:
        option = sys.argv[1]
        if option in valid_options:
            n = 2
        else:
            option = JOIN_OPTION  # default option
            n = 1
        opts, args = getopt.getopt( sys.argv[n:], flags1, flags2 )
    except getopt.GetoptError, err:
        log.error(str(err)) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    except:
        usage()
        sys.exit(2)

    directory = '.'
    basename = None
    work_directory = '.'
    inputs = []
    output = None
    Z = [None,None]
    bin = 1
    test = False
    
    for option_,value_ in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-d","--directory"):
            directory = value_
        elif option_ in ("-b","--basename"):
            basename = value_
        elif option_ in ("-i","--input"):
            inputs = value_.split(',')
            inputs = [file for file in inputs if os.path.exists(file)]
        elif option_ in ("-o","--output"):
            output = value_
        elif option_ in ("-z"):
            Z = [float(s) for s in value_.split(',')]
        elif option_ in ("--bin"):
            bin = int(value_)
        elif option_ in ("--test"):
            test = True
        else:
            assert False, "unhandled option"

    if len(inputs)==0 and basename==None:
        usage()
        raise SystemExit

#    if basename!=None:
#        basenames = txbr.utilities.extract_series_name_from_dir( basename, directory )
#        basename = ",".join(basenames)

    if option==CROSS_VALIDATION_OPTION:
        crossvalidate( directory, inputs, output)
    elif option==MERGE_OPTION:
        combineMRCFiles( directory, inputs, output)
    elif option==JOIN_OPTION:
        txbr.join.joinVolumes( directory, basename, work_directory, output, zmin=Z[0], zmax=Z[1], bin=bin, test=test )
    elif option==JOIN_CHECK_OPTION:
        txbr.join.joinVolumes( directory, basename, work_directory, output, zmin=Z[0], zmax=Z[1], check_mask=True )
    elif option==JOIN_MICROGRAPHS:
        txbr.join.joinMicrographs( directory, basename, work_directory, output, test=test )
    elif option==WEIGTH_VOLUME:
        if len(inputs)!=2:
            sys.exit(2)
        volume_file = os.path.join( directory, inputs[0])
        weigth_file = os.path.join( directory, inputs[1])
        output = volume_file + '_dest'
        txbr.join.weightVolume( volume_file, weigth_file, output, test=test )

main()
