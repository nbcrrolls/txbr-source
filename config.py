import os
import os.path
import sys
import re
import getopt
import ctypes.util
import ConfigParser

'''
What is needed for TxBR:
    * numpy (tested with version 2.0)
    * scipy (tested with version 0.10.0) -- requires blans and lapack
    * matplotlib (tested with version 1.0.1)
    * psutil (tested with version 0.5.1)
    * Pyrex (tested with version 0.9.9)
    * PIL/Imaging (tested with version 1.1.7)
    * CLN (tested with version 1.3.1)
    * GiNaC (tested with version 1.5.8)
    * SWIG (tested with version 2.0.2)
    * swiginac (tested with version 1.5.1)
    * fftw (tested with version 3.2.2)
    * mpi4py (tested with version 1.2.2)
    * openCV (tested with version 2.2.0)
    * IMOD (tested with version 4.5.8)
The python interpreter was version 2.7.1 on a linux machine.
'''

def search_module( token ):

    roots = [ "/opt", "/usr/local", os.environ['HOME'] ]
    paths = []

    for root in roots:
        if len(paths)>0: break
        cmd = "find %s -name %s -print" %( root, token )
        print cmd
        for file in os.popen(cmd).readlines():
            if file.startswith("find:") or file.endswith("Permission denied"): continue
            head,tail = os.path.split(file.rstrip())
            paths.append(head)

    if len(paths)>1:
        folders = [ os.path.split(path)[1] for path in paths ]
        try:
            return paths[folders.index('lib64')]
        except:
            pass
        try:
            return paths[folders.index('lib')]
        except:
            pass
        return paths[0]
    elif len(paths)==1:
        return paths[0]
    else:
        return None


def changeAndCopy( src, dest, txbr_dir, prefix, version ):
    
    print "Modify and copy file %s to %s" %(src,dest)
    
    f = open(src)
    s = f.read()
    f.close()

    s = re.compile("\\$TXBR_DIR").sub(txbr_dir,s)
    s = re.compile("\\$INSTALL_DIR").sub(prefix,s)
    s = re.compile("\\$VERSION").sub(version,s)

    f = open(dest,"w")
    f.write(s)
    f.close()


def main():

    flags1 = ""
    flags2 = [ "gpu", "prefix=" ]

    try:
        opts, args = getopt.getopt(sys.argv[1:], flags1, flags2)
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    run_gpu = False
    prefix = os.path.join(os.environ['HOME'],"usr","local")

    for option,value in opts:
        if option in ('--prefix'):
            prefix = os.path.abspath(value)
        elif option in ('--gpu'):
            run_gpu = True

    if not os.path.exists(prefix):
        sys.exit( "Error, the directory %s does not exist!" %prefix )

    try:
        f = open('VERSION', 'r')
        version = f.readline().rstrip('\n')
    except:
        version = 'dev'
        pass

    txbr_dir = os.path.join(prefix,"TxBR-%s" %version )

    if not os.path.exists(txbr_dir):
        os.mkdir(txbr_dir)

    print "TxBR release directory: %s" %txbr_dir

    base_dict = {}
    base_dict["TXBR_BASE_DIR"] = os.path.abspath(".")
    if "MATLAB" in base_dict.keys():
        base_dict["MATLAB_BASE_DIR"] = os.environ["MATLAB"]

    # Change and Copy the configuration files

    src_txbr_bashrc = os.path.join("config","txbr_bashrc")
    dest_txbr_bashrc = os.path.join(txbr_dir,"txbr_bashrc")

    changeAndCopy( src_txbr_bashrc, dest_txbr_bashrc, txbr_dir, prefix, version)

    src_txbr_cshrc = os.path.join("config","txbr_cshrc")
    dest_txbr_cshrc = os.path.join(txbr_dir,"txbr_cshrc")
    
    changeAndCopy( src_txbr_cshrc,dest_txbr_cshrc, txbr_dir, prefix, version )
    
    # Check different libraries

    include_ext = []
    lib_ext = []

    try:
        import numpy
        include_ext.append(sys.modules['numpy'].__path__[0] + '/core/include')
    except:
        sys.exit("numpy is not properly installed")

    try:
        import scipy
    except:
        sys.exit("scipy is not properly installed")

    try:
        import matplotlib
    except:
        sys.exit("matplotlib is not properly installed.")

    try:
        import Pyrex
    except:
        sys.exit("Pyrex is not properly installed")

    try:
        import swiginac
    except:
        sys.exit("Error, GiNaC or its python bindings are not properly installed.")

    try:
        import PIL
    except:
        sys.exit("Error, PIL is not properly installed.")

    try:
        import mpi4py
    except:
        sys.exit("Error, mpi4py or its python bindings are not properly installed.")

    if ctypes.util.find_library('fft3w')==None:
        print "fft3w is not in a standard location."

    if ctypes.util.find_library('fft3w')==None:
        print "fft3wf is not in a standard location."

    if ctypes.util.find_library('imod')==None:
        print "IMOD is not in a standard location."
        libmod = search_module('libimod.so')
        if libmod!=None:
            lib_ext.append(libmod)
            head,tail = os.path.split(libmod)
            lib_ext.append(os.path.join(head,"qtlib"))
            base_dict["IMOD_BASE_DIR"] = head

    if ctypes.util.find_library('opencv_core')==None:
        print "OpenCV is not in a standard location."
        try:
            import cv
            head,tail = os.path.split(cv.__file__)
            while tail!='lib': head,tail = os.path.split(head)
            include_ext.append(os.path.join(head,'include'))
            include_ext.append(os.path.join(head,'include','opencv'))
            include_ext.append(os.path.join(head,'include','opencv2'))
            lib_ext.append(os.path.join(head,'lib'))
        except:
            sys.exit("Error, openCV is not properly installed.")

    if ctypes.util.find_library('cudart')==None:
        print "Cuda is not in a standard location."
        libmod = search_module('libcudart.so')
        if libmod!=None:
            lib_ext.append(libmod)
            head,tail = os.path.split(libmod)
            base_dict["CUDA_BASE_DIR"] = head
        incmod = search_module('cuda.h')
        if incmod!=None: include_ext.append(incmod)

    include_ext.append(os.path.join(os.path.abspath("."),"txbr","include","imod"))
    include_ext.append(os.path.join(os.path.abspath("."),"txbr","include","txbr"))
    include_ext.append(os.path.join(os.path.abspath("."),"txbr","include","pxd"))

    lib_ext.append(os.path.join(txbr_dir,"lib"))
            
    print "Include directories: %s" %(include_ext)
    print "Lib directories: %s" %(lib_ext)

    # Create the ConfigInfo.txt file

    FILE = open(os.path.join('txbr','txbr','bckprj','ConfigInfo.txt'),"w")
    for key,value in base_dict.items():
        FILE.write( "%s = %s\n" %(key,value) )
    FILE.close()

    print base_dict
   
    # Create the config object

    config = ConfigParser.RawConfigParser()

    config.add_section('install')
    config.set('install', 'install-purelib', '%s/lib' %txbr_dir)
    config.set('install', 'install-platlib', '%s/lib' %txbr_dir)
    config.set('install', 'install-scripts', '%s/scripts' %txbr_dir)
    config.set('install', 'install-data', '%s/data' %txbr_dir)

    config.add_section('build_ext')
    config.set( 'build_ext', 'include-dirs', ':'.join(include_ext) )
    config.set( 'build_ext', 'library-dirs', ':'.join(lib_ext)  )

    # Writing the setup configuration file
    configfile = open('setup.cfg', 'wb')
    config.write(configfile)

    # Extra operations for the gpu

    if run_gpu:
        txbr_dir = os.path.abspath(".")
        os.system("make -C %s" %os.path.join(txbr_dir,"txbr","txbr","bckprj","cuda"))
        os.system("python gpu_query.py")


if __name__ == '__main__': main()
