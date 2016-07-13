#ifndef BCKPRJ_H

#define BCKPRJ_H

#include "mrcfiles.h"
#include "mrcslice.h"
#include "txbr.h"
#include "txbrutil.h"

#undef TEST
#define TEST 0

#undef OFFSET
#define OFFSET	1	// The offset in the projection map is usually one (IMOD viewer)

#undef EPSILON
#define EPSILON 1e-8

#undef SWAP_x_y
#define SWAP_x_y 0

#undef FULL_ROTATION
#define FULL_ROTATION 1

#undef TILT_NORMALIZATION
#define TILT_NORMALIZATION 2    // in (0,1,2)

#undef SLICE_NORMALIZATION
#define SLICE_NORMALIZATION 2    // in (0,1,2)

#undef MINIMUM_TILT_NUMBER_FOR_NORMALIZATION
#define MINIMUM_TILT_NUMBER_FOR_NORMALIZATION 0

#undef MINIMUM_SLICE_NUMBER_FOR_NORMALIZATION
#define MINIMUM_SLICE_NUMBER_FOR_NORMALIZATION 0

#undef FIX_SEGMENT_SIZE
#define FIX_SEGMENT_SIZE 10

#define slicePutVal_float(slice_data_f_, n_x_, x_, y_, v_) \
    slice_data_f_[(x_) + ((y_) * (n_x_))] = (v_)

#define sliceGetPixelMagnitude_float(slice_data_f_, n_x_, x_, y_) \
    slice_data_f_[(x_) + ((y_) * (n_x_))]

#define INTERP2(slice_data,nx,ny,i,j,wi,wj) ( \
		sliceGetPixelMagnitude_float(slice_data,nx,i,j)*(1.0-wi)*(1-wj) + \
		sliceGetPixelMagnitude_float(slice_data,nx,i,j+1)*(1-wi)*wj + \
		sliceGetPixelMagnitude_float(slice_data,nx,i+1,j)*wi*(1-wj) + \
		sliceGetPixelMagnitude_float(slice_data,nx,i+1,j+1)*wi*wj \
		)
#endif

