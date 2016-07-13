

/* This sample queries the properties of the CUDA devices present in the system. */

// Based on a sample app in the CUDA toolkit 

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// includes, project
#include "mycutil.h"

// forward decls
void printAllDevices();
void printUsage();

// Program will either print or return values whereever possible
// Code returns '-1' where there is no value to return (like in -printall) 
// and '-2' if there's an error or data could not be found.
int main( int argc, char** argv) 
{
	int i;
	int deviceCount = 0, devIndex = -1;
	cudaDeviceProp deviceProp;
	
	CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
	
	// parse command line args
	if(argc < 2) {
		printUsage();
		//printAllDevices();
		return(-1);
	}
	
	for(i=1; i < argc; i++)
	{
		if(strcmp(argv[i],"-h") == 0)
		{
			printUsage();
			return(-1);
		}
		
		if(strcmp(argv[i],"-printall") == 0)
		{
			printAllDevices();
			return(-1);
		}
		
		if(strcmp(argv[i],"-i") == 0)
		{
			devIndex = atoi(argv[++i]);
			if(devIndex < 0 || devIndex > (deviceCount - 1)) 
			{
				printf("\ndeviceQuery: Device number %d is out of range", devIndex);
				return (-2);
			}
			
			// get the device properties
			CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, devIndex));
			
			// some checks
			if (deviceProp.major == 9999 && deviceProp.minor == 9999)
			{
				printf("\nDevice does not support CUDA");
				return (-2);
			}
			
			continue;
		}
		
		if(strcmp(argv[i],"-n") == 0)
		{
			printf("%d", deviceCount);
			return(deviceCount);
		}
		
		if(strcmp(argv[i],"-major") == 0)
		{
			if(devIndex < 0) return (-2);
			
			printf("%d", deviceProp.major);
			return (deviceProp.major);
		}
		
		if(strcmp(argv[i],"-minor") == 0)
		{
			if(devIndex < 0) return (-2);
			
			printf("%d", deviceProp.minor);
			return (deviceProp.minor);
		}

		if(strcmp(argv[i],"-mem") == 0)
		{
			if(devIndex < 0) return (-2);
			
			printf("%d", (int)(deviceProp.totalGlobalMem / 1000000));
			return ((int)(deviceProp.totalGlobalMem / 1000000));
		}
		
		if(strcmp(argv[i],"-cores") == 0)
		{
			if(devIndex < 0) return (-2);
			
#if CUDART_VERSION >= 2000	
			printf("%d", (int)(8 * deviceProp.multiProcessorCount));
			return ((int)(8 * deviceProp.multiProcessorCount));
#endif
		}
		
		if(strcmp(argv[i],"-name") == 0)
		{
			if(devIndex < 0) return (-2);
			
			printf("%s", deviceProp.name);
			return (-1);
		}
		
		if(strcmp(argv[i],"-clock") == 0)
		{
			if(devIndex < 0) return (-2);
			
			printf("%0.2f", deviceProp.clockRate * 1e-6f);
			return (-1);
		}
	
	} // end of for (i ..
	
	// If the code reached here, we could not make sense of the command line
	printUsage();
	return(-1);
}

void printUsage()
{
	printf("\ndeviceQuery <Options>");
	printf("\n\nOptions:");
	printf("\n---------");
	printf("\n-h\t\t: Print this help");
	printf("\n-printall\t: Prints information about all CUDA devices in the system.");
	printf("\n-i <num>\t: Choose device number given by <num>. Numbered as 0,1,2..");
	printf("\n-n\t\t: Returns number of CUDA devices found");
	printf("\n-major\t\t: Returns major revision number");
	printf("\n-minor\t\t: Returns minor revision number");
	printf("\n-mem\t\t: Returns global memory on the device (in MB)");
	printf("\n-cores\t\t: Returns total number of cores on the GPU");
	printf("\n-name\t\t: Print Device name");
	printf("\n-clock\t\t: Print clock rate in Ghz");
	
	printf("\n\nCode returns '-1' where there is no value to return (like in -printall or -name)");
	printf("\nand '-2' if there's an error or data could not be found.");
	printf("\n\n");

	return;
}

void printAllDevices()
{
    int deviceCount;
    CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
    if (deviceCount == 0)
        printf("There is no device supporting CUDA\n");
    int dev;
    for (dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));
        if (dev == 0) {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999)
                printf("There is no device supporting CUDA.\n");
            else if (deviceCount == 1)
                printf("There is 1 device supporting CUDA\n");
            else
                printf("There are %d devices supporting CUDA\n", deviceCount);
        }
        printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);
        printf("  Major revision number:                         %d\n",
               deviceProp.major);
        printf("  Minor revision number:                         %d\n",
               deviceProp.minor);
        printf("  Total amount of global memory:                 %u bytes\n",
               deviceProp.totalGlobalMem);
    #if CUDART_VERSION >= 2000
        printf("  Number of multiprocessors:                     %d\n",
               deviceProp.multiProcessorCount);
        printf("  Number of cores:                               %d\n",
               8 * deviceProp.multiProcessorCount);
    #endif
        printf("  Total amount of constant memory:               %u bytes\n",
               deviceProp.totalConstMem); 
        printf("  Total amount of shared memory per block:       %u bytes\n",
               deviceProp.sharedMemPerBlock);
        printf("  Total number of registers available per block: %d\n",
               deviceProp.regsPerBlock);
        printf("  Warp size:                                     %d\n",
               deviceProp.warpSize);
        printf("  Maximum number of threads per block:           %d\n",
               deviceProp.maxThreadsPerBlock);
        printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
               deviceProp.maxThreadsDim[0],
               deviceProp.maxThreadsDim[1],
               deviceProp.maxThreadsDim[2]);
        printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
               deviceProp.maxGridSize[0],
               deviceProp.maxGridSize[1],
               deviceProp.maxGridSize[2]);
        printf("  Maximum memory pitch:                          %u bytes\n",
               deviceProp.memPitch);
        printf("  Texture alignment:                             %u bytes\n",
               deviceProp.textureAlignment);
        printf("  Clock rate:                                    %.2f GHz\n",
               deviceProp.clockRate * 1e-6f);
    #if CUDART_VERSION >= 2000
        printf("  Concurrent copy and execution:                 %s\n",
               deviceProp.deviceOverlap ? "Yes" : "No");
    #endif
    }
    
    return;

}
