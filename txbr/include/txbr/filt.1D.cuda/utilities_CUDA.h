#ifndef UTILITIES_CUDA_H
#define UTILITIES_CUDA_H

#include "utilities_base.h"

// Do NOT include CUDA headers here.

//#ifdef __cplusplus /* If this is a C++ compiler, use C linkage. */
//// NOTE: This is C99 code with CUDA extensions.
//extern "C"
//{
//#endif

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Output utilities. {
////////////////////////////////////////////////////////////////////////////////

__host__ void Abort_with_info_CUDA
(
    const char *file,
    int line,
    const char *func,
    cudaError_t cudaErrorCode, // WARNING: Undefined without CUDA headers.
    char *format,
    ...
);

#define ABORT_CUDA(cudaErrorCode_, format_, ...) \
    Abort_with_info_CUDA(__FILE__, __LINE__, __func__, (cudaErrorCode_), (format_), ## __VA_ARGS__)

////////////////////////////////////////////////////////////////////////////////
// END: Output utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Device utilities. {
////////////////////////////////////////////////////////////////////////////////

__host__ void print_device_properties_CUDA(int deviceID);
__host__ void check_device_properties_CUDA(int deviceID);
__host__ void print_and_check_all_devices_properties_CUDA();

////////////////////////////////////////////////////////////////////////////////
// END: Device utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Matrix/Matrix-like data structure access macros. {
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// END: Matrix/Matrix-like data structure access macros. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Memory management. {
////////////////////////////////////////////////////////////////////////////////

__host__ void mallocHost_CUDA(void **hostPtr, size_t size);
__host__ void freeHost_CUDA(void *hostPtr);

////////////////////////////////////////////////////////////////////////////////
// END: Memory management. }
////////////////////////////////////////////////////////////////////////////////

//#ifdef __cplusplus /* If this is a C++ compiler, use C linkage. */
//}
//#endif

#endif // UTILITIES_CUDA_H
