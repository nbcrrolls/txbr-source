#!/bin/bash

# dataset parameters
basename="fhv6a";
directory="$DATA_DIR/fhv6"; 
work_directory="$DATA_DIR/fhv6";

export TXBR_BCK="$TXBR/txbr/txbr/bckprj/back.cu.2"

echo "Auto-generating GPU list available on the cluster ..."

#get list of nodes
hostCount=`rocks list host | wc -l`
hostCount=`expr $hostCount - 2`

hostList=`rocks list host | tail -$hostCount | cut -f 1 -d':'`

#for each host find the number of GPUs
totalGpuCount=0
tempHostCount=0
for host in $hostList
do
	gpuCount=0
	gpuCount=`ssh -x $host $TXBR_BCK/../cuda/deviceQuery -n`
	totalGpuCount=`expr $gpuCount + $totalGpuCount`

	gpuCountOnHost[$tempHostCount]=gpuCount
	tempHostCount=`expr $tempHostCount + 1`
done

echo $totalGpuCount

# For every host spawn a process per GPU
blocksize=5
z_min=-60
z_delta=5
host_Idx=0

z_start=$z_min
z_stop=60

for host in $hostList
do
	echo "Spawning processes for $host with $gpuCountOnHost[$hostIdx] GPUs"
	gpuIdx=0
	while [ $gpuIdx -lt 4 ]
	do
		z_max=`expr $z_min + $z_delta`
		ssh -fx $host $TXBR_BCK/run_backprojection -g $gpuIdx -b $basename -d $directory -w $work_directory -zstart $z_min -zstop $z_max
		
		gpuIdx=`expr $gpuIdx + 1`

echo $z_min $z_max
		z_min=$z_max
	done

	hostIdx=`expr $hostIdx + 1`
done

# Merge individual
#./txbr_finalize_bckprj_cu.py -b $basename -d $directory -n $totalGpuCount -z $z_start,$z_stop 


#gpuCount=`$TXBR/deviceQuery -n`

# $gpuCount
