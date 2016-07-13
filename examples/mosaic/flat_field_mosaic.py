import sys
import os
import os.path
import numpy
import mrc
import getopt
import txbr.utilities


def usage():

    print
    print 'Usage: %s.py -b basename1[,...] [Options]' %os.path.basename(sys.argv[0])
    print
    print "    -d directory (default .)"
    print "        Name of the directory containing the input files (.fid and .txbr)"
    print "    --wd work_directory (default to directory)"
    print "        Name of the directory where the output files are stored"
    print "    --ones"
    print "        Performs a back-projection from unity images"
    print "    --normalize"
    print "        Flat field the 3D volume"


def backProjectOnes( directory, basenames, work_directory ):
    '''Performs a back-projection of unity source images'''

    basename, extensions = txbr.utilities.extract_series_name(basenames)

    project = txbr.TxBRproject( directory, basename, work_directory=work_directory)
    project.load()

    for extension in extensions:

        print "Create the images of ones for %s%s" %(basename,extension)
        
        file_in = "%s%s.preali.rot.SL" %(basename,extension)
        file_out = "%s%s.one" %(basename,extension)

        mrc_in = mrc.MRCFile( file_in )

        nx = mrc_in.nx
        ny = mrc_in.ny
        nz = mrc_in.nz

        mrc_out = mrc.MRCFile(file_out)
        mrc_out.setHeader(nx,ny,nz)

        for iz in range(nz):
            u = numpy.ones((nx,ny))
            mrc_out.setZSliceAt(iz,u)

        mrc_out.updateHeader()

        os.unlink( file_in )
        os.symlink( file_out,file_in )

    # Run the backprojection once those files are created

    os.system( "runtxbr.py -b %s --skip align,filter" %(",".join(basenames)) )

    final_vol = "%s_z_%.1f.out" %( basename, project.reconstruction.origin.z )

    os.rename( final_vol, final_vol + ".one" )


def flatFieldVolume( directory, basenames, work_directory ):

    basename, extensions = txbr.utilities.extract_series_name(basenames)

    project = txbr.TxBRproject( directory, basename, work_directory=work_directory)
    project.load()

    final_vol = "%s_z_%.1f.out" %( basename, project.reconstruction.origin.z )

    mrc_in = mrc.MRCFile( final_vol )
    mrc_nm = mrc.MRCFile( final_vol + ".one" )

    file_out = final_vol + ".nm"
    mrc_out = mrc.MRCFile(file_out, template=final_vol)

    nz = mrc_in.nz

    for iz in range(nz):
        u = mrc_in.getZSliceAt(iz) * 61.0 / mrc_nm.getZSliceAt(iz).astype('float')
        u = numpy.where(numpy.isnan(u),0,u)
        mrc_out.setZSliceAt(iz, u)

    mrc_out.updateHeader()


def main():

    flags1 = "hb:d:z:l:"
    flags2 = [ "help", "directory=", "wd=", "basename=", "ones", "normalize" ]

    try:
        opts, args = getopt.getopt(sys.argv[1:], flags1, flags2)
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    directory = "."
    work_directory = "."
    basename = None

    bckprj_ones = False
    flat_field = False

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
        elif option in ('--ones'):
            bckprj_ones = True
        elif option in ('--normalize'):
            flat_field = True
        else:
            assert False, "unhandled option"

    if basename==None:
        usage()
        sys.exit()

    directory = os.path.abspath(directory)
    work_directory = os.path.abspath(work_directory)

    basenames = txbr.utilities.extract_series_name_from_dir( basename, directory )

    if bckprj_ones:
        backProjectOnes( directory, basenames, work_directory )
    elif flat_field:
        flatFieldVolume( directory, basenames, work_directory )


if __name__=='__main__':
    main()

