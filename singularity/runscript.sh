#!/bin/bash

module load mpi/openmpi-x86_64
export LD_LIBRARY_PATH=/lib64:/usr/lib64/mpich/lib:/usr/local/IMOD/lib:/usr/local/IMOD/qtlib
. /usr/local/IMOD/IMOD-linux.sh

exec /usr/bin/runtxbr.py "$@"
