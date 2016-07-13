import os,os.path,shutil

txbr_dir = os.environ['TXBR_DIR']

lib_dir = os.path.join(txbr_dir,'lib')
scripts_dir = os.path.join(txbr_dir,'scripts')
data_dir = os.path.join(txbr_dir,'data')
txbr_cshrc = os.path.join(txbr_dir,'txbr_cshrc')
txbr_bashrc = os.path.join(txbr_dir,'txbr_bashrc')

if os.path.exists(lib_dir):
    print 'Clean up directory %s...' %lib_dir
    shutil.rmtree(lib_dir)

if os.path.exists(scripts_dir):
    print 'Clean up directory %s...' %scripts_dir
    shutil.rmtree(scripts_dir)
    
if os.path.exists(data_dir):
    print 'Clean up directory %s...' %data_dir
    shutil.rmtree(data_dir)
    
if os.path.exists(txbr_cshrc):
    print 'Clean up file %s...' %txbr_cshrc
    os.unlink(txbr_cshrc)

if os.path.exists(txbr_bashrc):
    print 'Clean up file %s...' %txbr_bashrc
    os.unlink(txbr_bashrc)

if os.path.exists('build'):
    print 'Clean up the build directory...'
    os.system('python setup.py clean -a')
