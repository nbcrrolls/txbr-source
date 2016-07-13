// BEGIN: CUDA includes
// NOTE: Leave these in the .c file.
#include <cuda.h>
#include <cutil.h>
#include "utilities_CUDA.h"
// END: CUDA includes

// BEGIN: MEX includes
// NOTE: Leave this in the .c file.

#ifdef MEX
#include "mex.h"
#endif

// END: MEX includes

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Output utilities. {
////////////////////////////////////////////////////////////////////////////////

__host__ void Abort_with_info_CUDA
(
    const char *file,
    int line,
    const char *func,
    cudaError_t cudaErrorCode,
    char *format,
    ...
)
{
    const int format_buffer_max_length = 1024;
    char format_buffer[format_buffer_max_length];

    if ((format != NULL) && (strlen(format) > 0))
        snprintf(format_buffer, format_buffer_max_length, "ABORTED CUDA: %s:%i:%s() -- (%s) %s", 
            file, line, func, cudaGetErrorString(cudaErrorCode), format);
    else
        snprintf(format_buffer, format_buffer_max_length, "ABORTED CUDA: %s:%i:%s() -- (%s).",
            file, line, func, cudaGetErrorString(cudaErrorCode));

    va_list vargs;

    va_start(vargs, format);
    vprintf(format_buffer, vargs); 
    va_end(vargs);

#ifndef MEX
    //memset((void *) 1, 13, 666); // Seg. fault (useful when debugging).
    exit(0);
#else
    const int mexErrMsgTxt_buffer_max_length = 1024;
    char mexErrMsgTxt_buffer[mexErrMsgTxt_buffer_max_length];

    va_start(vargs, format);
    vsnprintf(mexErrMsgTxt_buffer, mexErrMsgTxt_buffer_max_length, format_buffer, vargs); 
    va_end(vargs);

    mexErrMsgTxt(mexErrMsgTxt_buffer);
#endif
}

////////////////////////////////////////////////////////////////////////////////
// END: Output utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Device utilities. {
////////////////////////////////////////////////////////////////////////////////

__host__ void print_device_properties_CUDA(int deviceID)
{
    cudaDeviceProp deviceProp;

    CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceProperties(&deviceProp, deviceID));    

    Print("deviceID = %d\n", deviceID);

/*
struct cudaDeviceProp {
  char name[256];
  size_t totalGlobalMem;
  size_t sharedMemPerBlock;
  int regsPerBlock;
  int warpSize;
  size_t memPitch;
  int maxThreadsPerBlock;
  int maxThreadsDim[3];
  int maxGridSize[3];
  size_t totalConstMem;
  int major;
  int minor;
  int clockRate;
  size_t textureAlignment;
  int deviceOverlap;
  int multiProcessorCount;
}
*/

    //char name[256];
    Print("deviceProp.name = \"%s\"\n", deviceProp.name);
    //size_t totalGlobalMem;
    Print("deviceProp.totalGlobalMem = %d\n", deviceProp.totalGlobalMem);
    //size_t sharedMemPerBlock;
    Print("deviceProp.sharedMemPerBlock = %d\n", deviceProp.sharedMemPerBlock);
    //int regsPerBlock;
    Print("deviceProp.regsPerBlock = %d\n", deviceProp.regsPerBlock);
    //int warpSize;
    Print("deviceProp.warpSize = %d\n", deviceProp.warpSize);
    //size_t memPitch;
    Print("deviceProp.memPitch = %d\n", deviceProp.memPitch);
    //int maxThreadsPerBlock;
    Print("deviceProp.maxThreadsPerBlock = %d\n", deviceProp.maxThreadsPerBlock);
    //int maxThreadsDim[3];
    Print("deviceProp.maxThreadsDim[0-2] = [%d %d %d]\n", deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]);
    //int maxGridSize[3];
    Print("deviceProp.maxGridSize[0-2] = [%d %d %d]\n", deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]);
    //size_t totalConstMem;
    Print("deviceProp.totalConstMem = %d\n", deviceProp.totalConstMem);
    //int major;
    Print("deviceProp.major = %d\n", deviceProp.major);
    //int minor;
    Print("deviceProp.minor = %d\n", deviceProp.minor);
    //int clockRate;
    Print("deviceProp.clockRate = %d\n", deviceProp.clockRate);
    //size_t textureAlignment;
    Print("deviceProp.textureAlignment = %d\n", deviceProp.textureAlignment);
    //int deviceOverlap;
    Print("deviceProp.deviceOverlap = %d\n", deviceProp.deviceOverlap);
    //int multiProcessorCount;
    Print("deviceProp.multiProcessorCount = %d\n", deviceProp.multiProcessorCount);
}

__host__ void check_device_properties_CUDA(int deviceID)
{
    cudaDeviceProp deviceProp;

    CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceProperties(&deviceProp, deviceID));    

    if (deviceProp.major < 1) 
        ABORT("deviceID == %d does not support CUDA.\n", deviceID);
}

__host__ void print_and_check_all_devices_properties_CUDA()
{
    int device_count = 0;
    CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceCount(&device_count));

    Print("device_count = %d\n", device_count);

    for (int i_deviceID = 0; i_deviceID < device_count; ++i_deviceID)
    {   
        print_device_properties_CUDA(i_deviceID);
        check_device_properties_CUDA(i_deviceID);
    }   
}

////////////////////////////////////////////////////////////////////////////////
// END: Device utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Memory management. {
////////////////////////////////////////////////////////////////////////////////

__host__ void mallocHost_CUDA(void **hostPtr, size_t size)
{
    cudaError_t cudaErrorCode = cudaMallocHost(hostPtr, size);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMallocHost(hostPtr, %i) failed.\n", size);
}

__host__ void freeHost_CUDA(void *hostPtr)
{
    cudaError_t cudaErrorCode = cudaFreeHost(hostPtr);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaFreeHost(hostPtr) failed.\n");
}

////////////////////////////////////////////////////////////////////////////////
// END: Memory management. }
////////////////////////////////////////////////////////////////////////////////
