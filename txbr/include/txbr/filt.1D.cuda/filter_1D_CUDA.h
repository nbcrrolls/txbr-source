#ifndef FILTER_1D_CUDA_H
#define FILTER_1D_CUDA_H

// NOTE: These functions are called by code that knows nothing about CUDA.

#include "utilities_base.h"

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage. */
// NOTE: This is called by C99 code.
extern "C"
{
#endif

// BEGIN: CUDA params.

typedef struct general_CUDA_params_adt *general_CUDA_params_handle;

general_CUDA_params_handle create_general_CUDA_params_from_data_copy
(
    int dimBlock_x,
    int dimBlock_y,
    int dimBlock_z
);

int general_CUDA_params_get_dimBlock_x(general_CUDA_params_handle general_CUDA_params);
int general_CUDA_params_get_dimBlock_y(general_CUDA_params_handle general_CUDA_params);
int general_CUDA_params_get_dimBlock_z(general_CUDA_params_handle general_CUDA_params);

void general_CUDA_params_release(general_CUDA_params_handle general_CUDA_params);

// END: CUDA params.

// BEGIN: Misc.

int calc_extent_pad_CUDA(int n_extent, int dimBlock_extent);

// END: Misc.

void projection_transform_2D_CUDA
(
    general_CUDA_params_handle general_CUDA_params,
    IplImage_base *projection0_padded,
    IplImage_base *projection0,
    real_type support0[4],
    real_type rot_matrix_2x2[2][2],
    // BEGIN: remap_2D_params
    real_type local_mag,
    int map_2D_order,
    int n_coeffs,
    int power_order_n_x1,
    int power_order_n_x2,
    real_type *power_order_h, // Pointer to raw data. 
    int map_2D_n_x1,
    int map_2D_n_x2,
    //int map_2D_n_x3, // Only current map is used.
    //real_type *map_2D, // Only current map is used.
    // END: remap_2D_params
    // BEGIN: remap_2D_data
    real_type *map_2D_current_h, // Pointer to raw data.
    // END: remap_2D_data
    IplImage_base *sums_image_base, // No allocated data.
    IplImage_base *hits_image_base, // No allocated data.
    IplImage_base *projection_padded,
    IplImage_base *projection,
    real_type support[4]
);

void projection_inv_transform_2D_CUDA
(
    general_CUDA_params_handle general_CUDA_params,
    IplImage_base *projection_padded,
    IplImage_base *projection,
    real_type support[4],
    real_type rot_matrix_2x2[2][2],
    // BEGIN: remap_2D_params
    real_type local_mag,
    int map_2D_order,
    int n_coeffs,
    int power_order_n_x1,
    int power_order_n_x2,
    real_type *power_order_h, // Pointer to raw data. 
    int map_2D_n_x1,
    int map_2D_n_x2,
    //int map_2D_n_x3, // Only current map is used.
    //real_type *map_2D, // Only current map is used.
    // END: remap_2D_params
    // BEGIN: remap_2D_data
    real_type *map_2D_current_h, // Pointer to raw data.
    // END: remap_2D_data
    IplImage_base *projection0_padded,
    IplImage_base *projection0,
    real_type support0[4]
);

void projection_rotate_2D_CUDA
(
    general_CUDA_params_handle general_CUDA_params,
    IplImage_base *projection0_padded,
    IplImage_base *projection0,
    real_type support0[4],
    real_type rot_matrix_2x2[2][2],
    IplImage_base *projection_padded,
    IplImage_base *projection,
    real_type support[4]
);

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Hardware initialization. {
////////////////////////////////////////////////////////////////////////////////

void projection_series_transform_2D_CUDA_filter_1D_inv_transform_2D_CUDA_init_GPU(int deviceID);
void projection_series_transform_2D_CUDA_filter_1D_inv_transform_2D_CUDA_clear_GPU(int deviceID);

void projection_series_rotate_2D_CUDA_filter_1D_rotate_2D_CUDA_init_GPU(int deviceID);
void projection_series_rotate_2D_CUDA_filter_1D_rotate_2D_CUDA_clear_GPU(int deviceID);

////////////////////////////////////////////////////////////////////////////////
// END: Hardware initialization. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Memory management. {
////////////////////////////////////////////////////////////////////////////////

void malloc_page_locked_pixel_type_CUDA(pixel_type **data, size_t size);
void free_page_locked_pixel_type_CUDA(pixel_type *data);

////////////////////////////////////////////////////////////////////////////////
// END: Memory management. }
////////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage. */
}
#endif

#endif // FILTER_1D_CUDA_H
