#!/bin/bash

mkdir -p txbr
cp -r ../* txbr/

echo ""
echo "---- Starting VM This may take a few minutes -----"
echo ""
vagrant up
echo ""
echo "Connect to vm by invoking: vagrant ssh"
echo ""
echo "or for ssh with X redirected"
echo "vagrant ssh -- -X"
echo ""
echo "And be sure to invoke the following before building"
echo "txbr:"
echo ""
echo "module load mpi/openmpi-x86_64"
echo "export LD_LIBRARY_PATH=/lib64:/usr/lib64/mpich/lib:/usr/local/IMOD/lib:/usr/local/IMOD/qtlib"
echo ""

