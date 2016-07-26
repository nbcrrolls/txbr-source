import os
import sys

from distutils.core import setup, Extension
from Pyrex.Distutils import build_ext


def numberOfGPUBoards():
    """Helper function to check if GPU code should be compiled"""

    txbr_dir = "."
    infofile = os.path.join(txbr_dir,'data','machinefile.gpu.txt')

    n_gpu_boards = 0

    try:
        f = open(infofile,'r')
        numberOfGPUBoardsForNode = eval(f.read().strip())
        f.close()
    except:
        return n_gpu_boards

    for value in numberOfGPUBoardsForNode.values(): n_gpu_boards += value

    return n_gpu_boards


# Check if GPU code should be built

buildGPUCode = numberOfGPUBoards()!=0

mrc = Extension( 'mrc.mrcfile',
                  libraries = ['tiff', 'cfshr', 'gomp', 'iimod', 'imod'],
                  sources = [ 'txbr/mrc/mrcfile.c',
                              'txbr/mrc/txbrstack.c',
                              'txbr/util/util_C/txbrutil.c' ]
                 )
                 
imodel = Extension( 'modl.modfile',
                    libraries = [ 'tiff', 'cfshr', 'gomp', 'iimod', 'imod' ],
                    sources = ['txbr/modl/modfile.c']
                 )

filter = Extension( 'txbr.filter.filter',
                    libraries = [ 'iimod', 'imod', 'cfshr', 'gomp', 'txbrutil' ],
                    sources = [ 'txbr/txbr/filter/3D/filter.c',
                                'txbr/txbr/filter/3D/filter_3D.c',
                                'txbr/txbr/filter/3D/txbrtraj.c' ]
                  )

contres = Extension( "txbr.align.residualSupp",
                     sources = ["txbr/txbr/align/residualSupp.pyx"],
                     include_dirs = ["./txbr/include/pxd"]
                    )

preconstruction = Extension( "txbr.prefil",
                              libraries = [ 'iimod', 'imod', 'cfshr', 'gomp',  'fftw3', 'fftw3f', 'opencv_core', 'opencv_features2d' ],
                              sources = [ "txbr/txbr/filter/run_prefil.c",
                                          "txbr/txbr/filter/filter_1D.c",
                                          "txbr/txbr/filter/utilities.c",
                                          "txbr/txbr/filter/utilities_base.c",
                                          "txbr/txbr/filter/edgetaper_2D.c" ],
                              extra_compile_args = ["-std=gnu99"]
                            )

backprojection_0 = Extension( "txbr.bckprj",
                              libraries = [ 'iimod', 'imod', 'cfshr', 'gomp' ],
                              sources = [ 'txbr/txbr/bckprj/run_bckprj.c',
                                          'txbr/txbr/bckprj/back.0/backprojection.0.c',
                                          'txbr/util/util_C/txbrutil.c',
                                          'txbr/txbr/utilities/util_C/txbrtext.c',
                                          'txbr/mrc/txbrstack.c' ],
                              extra_compile_args = ['-ffast-math','-O3']
                            )

backprojection_1 = Extension( "txbr.bckprj",
                              libraries = [ 'iimod', 'imod', 'cfshr', 'gomp' ],
                              sources = [ 'txbr/txbr/bckprj/run_bckprj.c',
                                          'txbr/txbr/bckprj/back.f/backprojection.f.c',
                                          'txbr/util/util_C/txbrutil.c',
                                          'txbr/txbr/utilities/util_C/txbrtext.c',
                                          'txbr/mrc/txbrstack.c' ],
                              extra_compile_args = ['-fno-math-errno','-ffast-math','-O3']
                            )

backprojection_2 = Extension( "txbr.bckprj",
                               libraries = [ 'iimod', 'imod', 'cfshr', 'gomp'],
                               sources = [ 'txbr/txbr/bckprj/run_bckprj.c',
                                           'txbr/txbr/bckprj/back.2/backprojection.2.c',
                                           'txbr/txbr/bckprj/back.2/calculateBox.c',
                                           'txbr/txbr/bckprj/back.2/initializeGPU.c',
                                           'txbr/txbr/bckprj/back.2/evaluateSegmentGPU.c',
                                           'txbr/txbr/bckprj/back.2/backprojectionGPUsupport.c',
                                           'txbr/util/util_C/txbrutil.c',
                                           'txbr/txbr/utilities/util_C/txbrtext.c',
                                           'txbr/mrc/txbrstack.c' ],
                            )

backprojection_cu2 = Extension( "txbr.bckprj_cu",
                                 libraries = [ 'iimod', 'imod', 'cfshr', 'gomp','cudart'],
                                 sources = [ 'txbr/txbr/bckprj/run_bckprj_cuda.c',
                                             'txbr/txbr/bckprj/back.cu.2/backprojection.2.c',
                                             'txbr/txbr/bckprj/back.cu.2/calculateBox.c',
                                             'txbr/util/util_C/txbrutil.c',
                                             'txbr/txbr/utilities/util_C/txbrtext.c',
                                              'txbr/mrc/txbrstack.c' ],
                                extra_objects = [ 'txbr/txbr/bckprj/back.cu.2/backprojectionGPUsupport.o' ],
                            )

combine = Extension( "txbr.combine",
                     libraries = [ 'iimod', 'imod', 'cfshr', 'gomp' ],
                     sources = [ 'txbr/txbr/join/combine.c',
                                 'txbr/txbr/join/merge.c',
                                 'txbr/txbr/join/crossval.c',
                                 'txbr/util/util_C/txbrutil.c' ],
                  )

remap = Extension( "misc.remap.remapSupp",
                    sources = [ "txbr/misc/remap/remapSupp.pyx" ],
                    include_dirs = ["./txbr/include/pxd"]
                    )

ginac_plug = Extension( "util.ginacplug",
                        sources = [ "txbr/util/ginacplug.cc" ],
                        libraries = [ 'cln' ]
                        )

backproj_cpu = backprojection_1
backproj_gpu = backprojection_cu2

# Setup the TxBR pakage

scripts_ = [ 'scripts/scanjobs.py', \
             'scripts/stitch.py', \
             'scripts/ccdb_publish.py', \
             'scripts/ccdb_jail_fix_mpid.py', \
             'scripts/ccdb_sync.py', \
             'scripts/txbr_version.py', \
             'scripts/runtxbr.py', \
             'scripts/txbr_align.py', \
             'scripts/txbr_auto.py', \
             'scripts/txbr_quick.py', \
             'scripts/txbr_filter.py', \
             'scripts/txbr_pitch.py', \
             'scripts/txbr_bckprj.py', \
             'scripts/txbr_phase.py', \
             'scripts/txbr_bin.py', \
             'scripts/txbr_test.py', \
             'scripts/txbr_view.py', \
             'scripts/txbr_flatten.py', \
             'scripts/txbr_onfly.py', \
             'scripts/txbr_rescale.py', \
             'scripts/txbr_sirt.py', \
             'scripts/txbr_stack.py', \
             'scripts/txbr_combine.py', \
             'scripts/txbr_mux.py', \
             'scripts/txbr_prepare.py', \
             'scripts/txbr_proj.py', \
             'scripts/txbr_plot.py', \
             'scripts/txbr_clean.py', \
             'scripts/txbr_cluster.py', \
             'scripts/txbr_drift.py', \
             'scripts/fly_flat.py', \
             'scripts/fidtrans.py', \
             'scripts/g3v_align.py', \
             'scripts/trans_mod.py', \
             'scripts/remap_4000#2.py', \
             'scripts/plot_boundaries.py', \
            ]

packages_ = [ 'txbr', \
              'txbr.align', \
              'txbr.filter', \
              'txbr.join', \
              'txbr.onthefly', \
              'txbr.phase', \
              'txbr.prj', \
              'txbr.setup', \
              'txbr.stack', \
              'txbr.utilities', \
              'txbr.g3view', \
              'modl', \
              'ccdb', \
              'mpfit', \
              'util', \
              'mrc', \
              'stats', \
              'misc', \
              'misc.remap', \
              'misc.series' ]

package_dir_ = { 'txbr':'txbr/txbr', \
                 'modl':'txbr/modl', \
                 'ccdb':'txbr/ccdb', \
                 'util':'txbr/util', \
                 'mrc':'txbr/mrc', \
                 'stats':'txbr/stats', \
                 'misc':'txbr/misc', \
                 'mpfit':'txbr/util/mpfit' }

package_data_ = { 'txbr': ['log.conf','txbr.cfg'], \
                  'ccdb': ['ccdb.cfg','machines.cfg','txbrmail.cfg'], \
                  'misc.remap': ['log.conf'] }

data_files_ = [ ('.',["VERSION"]),
                ('config',['data/machinefile.gpu.txt','data/machinefile.cpu.txt','data/config.fly','data/config.axis-0.fly','data/config.fei_spirit.fly','data/config.fei_titan.fly','data/config.fei_titan.de12.fly','data/config.j3200.fly','data/config.j3200.de12.fly','data/config.j4000_1.fly','data/config.j4000_2.fly',"VERSION"]),
                ('img',['data/txbr_logo.png','data/txbr_on_the_fly_logo.png']) ]

ext_modules_ = [ mrc, imodel, contres, backproj_cpu, combine, preconstruction, remap, ginac_plug ]

# Check if the GPU code should be built

if buildGPUCode:
    
    if len(sys.argv)>1 and sys.argv[1]=='clean':
        os.system('make --directory=./txbr/txbr/bckprj/cuda clean')
        os.system('make --directory=./txbr/txbr/bckprj/back.cu.2 clean')
        
    if len(sys.argv)>1 and ( sys.argv[1].startswith('build') or sys.argv[1]=='install' ):
        os.system('make --directory=./txbr/txbr/bckprj/cuda')
        os.system('make --directory=./txbr/txbr/bckprj/back.cu.2')
        
    ext_modules_.append(backproj_gpu)
    
# Read the version of TxBR as the first line of the file called VERSION

try:
    f = open('VERSION', 'r')
    version = f.readline().rstrip('\n')
except:
    version = 'beta'
    pass

# Do the setup

setup(
       name = 'TxBR',
       version = version,
       description = 'TxBR is an Electron Tomography package',
       author = ['Sebastien Phan','Alexander Ward Kulungowski','Raj Singh','Masako Terada','James Obayashi'],
       scripts = scripts_,
       packages = packages_,
       package_dir = package_dir_,
       package_data = package_data_,
       data_files = data_files_,
       ext_modules = ext_modules_,
       cmdclass = { 'build_ext': build_ext }
    )
