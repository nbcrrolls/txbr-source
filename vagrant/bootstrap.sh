#!/usr/bin/env bash

echo "Installing base packages"
yum install -y cmake git epel-release
yum install -y opencv opencv-devel mpich mpich-devel mpich-autoload qt qt-devel wget tcsh xauth xclock gcc-c++ mlocate time tree 
yum install -y xorg-x11-fonts-*
yum install -y mesa-*
yum install -y python-pip python-wheel
yum install -y python-setuptools python-setuptools-devel
yum install -y cln cln-devel ginac ginac-utils ginac-devel
yum install -y libtiff libtiff-devel fftw-*
yum install -y numpy numpy-f2py
yum install -y scipy python-matplotlib*
yum install -y mpi4py-mpich
yum install -y python-pillow*
yum install -y swig*
updatedb

pip install wheel
pip install swiginac

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
swiginac="swiginac-1.0.0.tgz"
wget https://pypi.python.org/packages/05/94/e339b91298bc06c1ab80d6572c0ffd94c0c1efcad4684f354dbaa32d0100/$swiginac
tar -zxf $swiginac
cd `echo $swiginac | sed "s/\.tgz//"`
python setup.by build
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
. /etc/profile.d/mpich-x86_64.sh

echo ""
echo "Installation complete..."
echo ""
echo ""
