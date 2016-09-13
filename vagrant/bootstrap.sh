#!/usr/bin/env bash

echo "Installing base packages"
yum install -y cmake git epel-release tcsh
yum install -y opencv opencv-devel opencv-python
yum install -y qt qt-devel wget tcsh xauth xclock gcc-c++ mlocate time tree 
yum install -y xorg-x11-fonts-* mesa-*
yum install -y python-pip python-wheel unzip
yum install -y python-setuptools python-setuptools-devel
yum install -y cln cln-devel ginac ginac-utils ginac-devel
yum install -y libtiff libtiff-devel fftw-*
yum install -y numpy numpy-f2py scipy python-matplotlib*
yum install -y mpi4py-openmpi python-pillow* PackageKit-gtk3-module libcanberra-gtk2 libcanberra-gtk3
yum install -y swig* python-psutil sympy* Perl-Data-Dumper openmpi openmpi-devel
updatedb

pip install wheel
pip install Cycler
pip install mpi4py

# Download & install pyrex
pyrex="Pyrex-0.9.9.tar.gz"
echo "Installing pyrex"
wget http://www.cosc.canterbury.ac.nz/greg.ewing/python/Pyrex/$pyrex

tar -zxf $pyrex
cd `echo $pyrex | sed "s/\.tar.gz.*//"`
python setup.py build
python setup.py install

cd ..

# Download & attempt to install swiginac
swiginac="swiginac_1.5.1.orig.tar.gz"
wget https://launchpad.net/ubuntu/+archive/primary/+files/$swiginac
tar -zxf $swiginac

# yes whoever made the above tarball named the directory differently
# then the tar file and switched _ with -
cd `echo "swiginac-1.5.1" | sed "s/\.tar.gz//"`
python setup.py build
python setup.py install
cd ..

# imod
imodfile="imod_4.7.15_RHEL6-64_CUDA6.0.csh"

echo "Imod installation downloading $imodfile"
wget http://bio3d.colorado.edu/imod/AMD64-RHEL5/$imodfile
chmod a+x $imodfile
./$imodfile -yes

/bin/rm -f $imodfile

. /etc/profile.d/modules.sh

echo "module load mpi/openmpi-x86_64" >> /home/vagrant/.bash_profile
echo "export LD_LIBRARY_PATH=/lib64:/usr/lib64/mpich/lib:/usr/local/IMOD/lib:/usr/local/IMOD/qtlib" >> /home/vagrant/.bash_profile

echo ""
echo "Installation complete..."
echo ""
echo ""
