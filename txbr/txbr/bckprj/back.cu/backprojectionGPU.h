

#ifndef _BACKPROJECTIONCUDA_H
#define _BACKPROJECTIONCUDA_H


#include <cutil.h>
#include <cuda.h>

#define CUDABP_DEV_ID			1	// defines the GPU to be used

#define CUDABP_ENTRY_FUNC_DECL		__global__
#define CUDABP_DEV_FUNC_DECL            __device__
#define CUDABP_DEV_VAR_DECL             __device__

//
//  NOTE : For a 4Kx6K with gridnum=256 and blocknum = 64, we see a >50x speedup
//


//#define CUDABP_GRID_XDIM		1
//#define CUDABP_GRID_YDIM		1
#define CUDABP_GRID_NUM			1		

//#define CUDABP_BLOCK_XDIM               1
//#define CUDABP_BLOCK_YDIM               64
#define CUDABP_BLOCK_NUM		1	

#endif  
