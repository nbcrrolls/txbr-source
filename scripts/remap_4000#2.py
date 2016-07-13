#!/usr/bin/python

import sys, os.path, getopt
import misc.remap.remap as remap

def usage():

    command = os.path.basename(sys.argv[0])

    print
    print "Usage: %s [Options] -b basename --model" %(command)
    print
    print "Options:"
    print "    -h"
    print "        help"
    print "    -d directory"
    print "    -n indexOfTilt"
    print "        applies the remap on tilt at indexOfTilt"
    print "    -m"
    print "        %i : remap is done with translations" %remap.TRANSLATIONS
    print "        %i : remap is done with affine transformations" %remap.AFFINE_TRANSFORMATIONS
    print "        %i : remap is done with projective transformations (version 1) (default)" %remap.PROJECTIVE_TRANSFORMATIONS_1
    print "        %i : remap is done with projective transformations (version 2)" %remap.PROJECTIVE_TRANSFORMATIONS_2
    print "    -q [0,1,2,3]"
    print "        Anchor or not the quadrants. Default: 0,1"
    print "    --interpolation"
    print "        %i : no interpolation (default)" %remap.NO_INTERPOLATION
    print "        %i : bilinear interpolation" %remap.BILINEAR_INTERPOLATION
    print "        %i : cubic interpolation" %remap.BICUBIC_INTERPOLATION


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:],
            "hb:d:n:m:q:",["help","directory=","basename=","model","interpolation="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    directory = "."
    basename = None
    mode = remap.PROJECTIVE_TRANSFORMATIONS_1
    interpolation = remap.NO_INTERPOLATION
    tilts = []
    qfix = [0,1]

    prepareStichingModel = False

    for option_,value in opts:
        if option_ in ("-h", "--help"):
            usage()
            sys.exit()
        elif option_ in ("-d","--directory"):
            directory = value
        elif option_ in ("-b","--basename"):
            basename = value
        elif option_ in ("-m"):
            mode = int(value)
        elif option_ in ("-q"):
            if value=='none':
                qfix = []
            else:
                qfix = [int(v) for v in value.split(",")]
        elif option_ in ("-n"):
            tilts = [int(v) for v in value.split(",")]
        elif option_ in ("--model"):
            prepareStichingModel = True
        elif option_ in ("--interpolation"):
            interpolation = int(value)
        else:
            assert False, "unhandled option"

    if basename==None:
        usage()
        raise SystemExit

    if prepareStichingModel:    # Create the stiching model
        remap.prepareStichingModel(directory, basename, ntilt)
    else:    # Remap the mrc files
        remap.remap(directory, basename, tilts=tilts, mode=mode,qfix=qfix,interpolation=interpolation)


main()
