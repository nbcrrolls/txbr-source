

#ifndef _BACKPROJECTIONCUDA_H
#define _BACKPROJECTIONCUDA_H

#include <cuda.h>
#include "mycutil.h"

//#define CUDABP_DEV_ID			1	// defines the GPU to be used

#define CUDABP_ENTRY_FUNC_DECL		__global__
#define CUDABP_DEV_FUNC_DECL            __device__
#define CUDABP_DEV_VAR_DECL             __device__

//
//  NOTE : For a 4Kx6K with gridnum=256 and blocknum = 64, we see a >50x speedup
//

#define CUDABP_GRID_SIZE		64 //64	
#define CUDABP_BLOCK_SIZE		64 //64

//#define MAX_ORDER       		20  	// Maximum order for the projection polynomial approximation
//#define MAX_COEFFICIENTS_NUMBER 	100 	// 1771  	// Maximum order for the projection polynomial approximation

#endif 


