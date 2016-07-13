// BEGIN: CUDA includes
// NOTE: Leave these in the .cu file.
#include <cuda.h>
#include <cutil.h>
//#include <cuda_runtime_api.h>
#include "utilities_CUDA.h"
// END: CUDA includes

/*
#ifdef __cplusplus /\* If this is a C++ compiler, use C linkage. *\/
// NOTE: This is called by C99 code.
extern "C"
{
#endif
*/

#include "filter_1D_CUDA.h"

// BEGIN: MEX includes
// NOTE: Leave this in the .cu file.

#ifdef MEX
#include "mex.h"
#endif

// END: MEX includes

#define EPSILON 1e-5

#define PROTOTYPE_COMPLIANT_INDEXING 0

//#define FLOOR_INNER_LOOP(x_) ((int_type) (x_)) // Will _not_ work with negative input.
#define FLOOR_INNER_LOOP(x_) ((int_type) floor(x_))

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Parameters. {
////////////////////////////////////////////////////////////////////////////////

// BEGIN: CUDA parameters.

typedef struct general_CUDA_params_adt // The public type is abstract.
{
    unsigned int dimBlock_x;
    unsigned int dimBlock_y;
    unsigned int dimBlock_z;
}
general_CUDA_params_type; // The private type is concrete.

general_CUDA_params_handle create_general_CUDA_params()
{
    general_CUDA_params_handle general_CUDA_params = (general_CUDA_params_handle) malloc(sizeof(general_CUDA_params_type));
    if (general_CUDA_params == NULL)
        ABORT("Cannot acquire memory for general_CUDA_params_handle general_CUDA_params.\n");

    // BEGIN: Default values.
    general_CUDA_params->dimBlock_x = 0;
    general_CUDA_params->dimBlock_y = 0;
    general_CUDA_params->dimBlock_z = 0;
    // END: Default values.

    return general_CUDA_params;
}

general_CUDA_params_handle create_general_CUDA_params_from_data_copy
(
    int dimBlock_x,
    int dimBlock_y,
    int dimBlock_z
)
{
    general_CUDA_params_handle general_CUDA_params = create_general_CUDA_params();

    if (dimBlock_x < 1)
        ABORT("dimBlock_x == %d < 1.\n", dimBlock_x);
    if (dimBlock_y < 1)
        ABORT("dimBlock_y == %d < 1.\n", dimBlock_y);
    if (dimBlock_z < 1)
        ABORT("dimBlock_z == %d < 1.\n", dimBlock_z);

    general_CUDA_params->dimBlock_x = (unsigned int) dimBlock_x;
    general_CUDA_params->dimBlock_y = (unsigned int) dimBlock_y;
    general_CUDA_params->dimBlock_z = (unsigned int) dimBlock_z;

    return general_CUDA_params;
}

int general_CUDA_params_get_dimBlock_x(general_CUDA_params_handle general_CUDA_params)
{
    return (int) general_CUDA_params->dimBlock_x;
}

int general_CUDA_params_get_dimBlock_y(general_CUDA_params_handle general_CUDA_params)
{
    return (int) general_CUDA_params->dimBlock_y;
}

int general_CUDA_params_get_dimBlock_z(general_CUDA_params_handle general_CUDA_params)
{
    return (int) general_CUDA_params->dimBlock_z;
}

void general_CUDA_params_release(general_CUDA_params_handle general_CUDA_params)
{
    //if (general_CUDA_params->member != NULL)
    //    free(general_CUDA_params->some_member);

    free(general_CUDA_params);
}

// END: CUDA parameters.

////////////////////////////////////////////////////////////////////////////////
// END: Parameters. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Misc. {
////////////////////////////////////////////////////////////////////////////////

int calc_extent_pad_CUDA(int n_extent, int dimBlock_extent)
{
    if (n_extent < 1)
        ABORT("n_extent == %d < 1.", n_extent);
    if (dimBlock_extent < 1)
        ABORT("dimBlock_extent == %d < 1.", dimBlock_extent);

    int dimGrid_extent = (int) ceil(((real_type) n_extent) / ((real_type) dimBlock_extent));
    return (dimBlock_extent * dimGrid_extent) - n_extent;
}

////////////////////////////////////////////////////////////////////////////////
// END: Misc. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Device code. {
////////////////////////////////////////////////////////////////////////////////

// WARNING: These persist for the lifetime of the application.

__device__ __constant__ static uint_address_type projection0_padded_data_min_address; // DEBUG.
__device__ __constant__ static uint_address_type projection0_padded_data_max_address; // DEBUG.
__device__ __constant__ static pixel_type *projection0_data;
__device__ __constant__ static int n0_x; // 4 bytes.
__device__ __constant__ static int n0_y; // 4 bytes.
__device__ __constant__ static int n0_x_ws; // 4 bytes.
__device__ __constant__ static real_type center0_x; // 4 or 8 bytes.
__device__ __constant__ static real_type center0_y; // 4 or 8 bytes.
__device__ __constant__ static real_type n0_x_min; // 4 or 8 bytes.
__device__ __constant__ static real_type n0_x_max; // 4 or 8 bytes.
__device__ __constant__ static real_type n0_y_min; // 4 or 8 bytes.
__device__ __constant__ static real_type n0_y_max; // 4 or 8 bytes.
__device__ __constant__ static uint_address_type projection_padded_data_min_address; // DEBUG.
__device__ __constant__ static uint_address_type projection_padded_data_max_address; // DEBUG.
__device__ __constant__ static pixel_type *projection_data;
__device__ __constant__ static int n_x; // 4 bytes.
__device__ __constant__ static int n_y; // 4 bytes.
__device__ __constant__ static int n_x_ws; // 4 bytes.
__device__ __constant__ static real_type center_x; // 4 or 8 bytes.
__device__ __constant__ static real_type center_y; // 4 or 8 bytes.
__device__ __constant__ static real_type n_x_min; // 4 or 8 bytes.
__device__ __constant__ static real_type n_x_max; // 4 or 8 bytes.
__device__ __constant__ static real_type n_y_min; // 4 or 8 bytes.
__device__ __constant__ static real_type n_y_max; // 4 or 8 bytes.
__device__ __constant__ static real_type step; // 4 or 8 bytes.
__device__ __constant__ static int n_subsamples_x; // 4 bytes.
__device__ __constant__ static int n_subsamples_y; // 4 bytes.
__device__ __constant__ static real_type sample_factor; // 4 or 8 bytes.
// BEGIN: rotate_2D only.
__device__ __constant__ static real_type rot_matrix00; // 4 or 8 bytes.
__device__ __constant__ static real_type rot_matrix01; // 4 or 8 bytes.
__device__ __constant__ static real_type rot_matrix10; // 4 or 8 bytes.
__device__ __constant__ static real_type rot_matrix11; // 4 or 8 bytes.
__device__ __constant__ static real_type x0_init; // 4 or 8 bytes.
__device__ __constant__ static real_type y0_init; // 4 or 8 bytes.
__device__ __constant__ static real_type x0_del_x; // 4 or 8 bytes.
__device__ __constant__ static real_type y0_del_x; // 4 or 8 bytes.
__device__ __constant__ static real_type x0_del_y; // 4 or 8 bytes.
__device__ __constant__ static real_type y0_del_y; // 4 or 8 bytes.
// END: rotate_2D only.
// BEGIN: transform_2D only.
__device__ __constant__ static int map_2D_order; // 4 bytes.
__device__ __constant__ static int n_coeffs; // 4 bytes.
__device__ __constant__ static int power_order_n_x1; // 4 bytes.
__device__ __constant__ static int power_order_n_x2; // 4 bytes.
__device__ __constant__ static real_type *power_order;
__device__ __constant__ static int map_2D_n_x1; // 4 bytes.
__device__ __constant__ static int map_2D_n_x2; // 4 bytes.
//__device__ __constant__ static int map_2D_n_x3; // Only current map is used.
__device__ __constant__ real_type *map_2D; // Only current map is used.
//__device__ __constant__ real_type *map_2D_current; // Only current map is used.
__device__ __constant__ pixel_type *sums_image_data; // Only needed on GPU.
__device__ __constant__ pixel_type *hits_image_data; // Only needed on GPU.
// END: transform_2D only.

////////////////////////////////////////////////////////////////////////////////
// BEGIN: transform_2D {
////////////////////////////////////////////////////////////////////////////////

#define PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY 1

#define MAX_MAP_2D_ORDER_CUDA 3

// NOTE: Consider inlining this.
// WARNING: The xy-flip looks reversed.
__device__ real_type f_x1_CUDA_d
(
#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
    real_type x1, real_type x2
#else
    real_type x2, real_type x1
#endif
)
{
    real_type x1_out = 0.0;

    for (int i = 0; i < n_coeffs; ++i)
        x1_out += 
            MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 0, i) *
            pow(x1, MATRIX_INDEX_2D(power_order, power_order_n_x2, 0, i)) *
            pow(x2, MATRIX_INDEX_2D(power_order, power_order_n_x2, 1, i));

    return x1_out;
}

// NOTE: Consider inlining this.
// WARNING: The xy-flip looks reversed.
__device__ real_type f_x2_CUDA_d
(
#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
    real_type x1, real_type x2
#else
    real_type x2, real_type x1
#endif
)
{
    real_type x1_out = 0.0;

    for (int i = 0; i < n_coeffs; ++i)
        x1_out += 
            MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 1, i) *
            pow(x1, MATRIX_INDEX_2D(power_order, power_order_n_x2, 0, i)) *
            pow(x2, MATRIX_INDEX_2D(power_order, power_order_n_x2, 1, i));

    return x1_out;
}

/*
// NOTE: Consider inlining this.
__device__ real_type f_x1_CUDA_d
(
    real_type x1, real_type x2, int n_coeffs,
    int power_order_n_x1, int power_order_n_x2, real_type *power_order,
    int map_2D_n_x1, int map_2D_n_x2, real_type *map_2D
)
{
    real_type x1_out = 0.0;

    for (int i = 0; i < n_coeffs; ++i)
        x1_out += 
            MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 0, i) *
            pow(x1, MATRIX_INDEX_2D(power_order, power_order_n_x2, 0, i)) *
            pow(x2, MATRIX_INDEX_2D(power_order, power_order_n_x2, 1, i));

    return x1_out;
}

// NOTE: Consider inlining this.
__device__ real_type f_x2_CUDA_d
(
    real_type x1, real_type x2,
    int n_coeffs,
    int power_order_n_x1, int power_order_n_x2, real_type *power_order,
    int map_2D_n_x1, int map_2D_n_x2, real_type *map_2D
)
{
    real_type x1_out = 0.0;

    for (int i = 0; i < n_coeffs; ++i)
        x1_out += 
            MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 1, i) *
            pow(x1, MATRIX_INDEX_2D(power_order, power_order_n_x2, 0, i)) *
            pow(x2, MATRIX_INDEX_2D(power_order, power_order_n_x2, 1, i));

    return x1_out;
}
*/

// NOTE: Inlined.
// WARNING: Why is the xy-flip hardcoded here?
#define projection_remap_2D_solve_ivp_CUDA_d_MACRO(i_f_x_, x_start_, y_start_, i_ivp_) \
do { \
    for(int j = 0; j < diff_2D_n_x2; ++j) \
    { \
        MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, 0, j) = \
            f_x ## i_f_x_ ## _CUDA_d \
            ( \
                y_start_, (((real_type) j) * step) + x_start_ \
            ); \
    } \
 \
   q ## i_ivp_[0] = MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, 0, 0); \
 \
   for (int i = 1; i < diff_2D_n_x1; ++i) \
   { \
       for (int j = 0; j < diff_2D_n_x2 - i; ++j) \
       { \
           MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i, j) =  \
               MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i - 1, j + 1) - \
               MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i - 1, j);\
       } \
       q ## i_ivp_[i] = MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i, 0); \
   } \
} while (0)

/*
// NOTE: Consider inlining this.
__device__ void projection_remap_2D_solve_ivp_CUDA_d
(
    int n_coeffs,
    int power_order_n_x1, int power_order_n_x2, real_type *power_order,
    int map_2D_n_x1, int map_2D_n_x2, real_type *map_2D,
    real_type (*f)(real_type, real_type, int, int, int, real_type *, int , int , real_type *),
    real_type x_start, real_type y_start, real_type step,
    int map_2D_order,
    real_type *diff_2D,
    real_type *ivp
)
{
    int diff_2D_n_x1 = map_2D_order + 1;
    int diff_2D_n_x2 = map_2D_order + 1;

    // Calculate the initial values.
    for(int j = 0; j < diff_2D_n_x2; ++j)
    {
#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
        MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, 0, j) = 
            (*f)(
                 y_start, (((real_type) j) * step) + x_start,
                 n_coeffs,
                 power_order_n_x1, power_order_n_x2, power_order, 
                 map_2D_n_x1, map_2D_n_x2, map_2D
                );
#else
        MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, 0, j) = 
            (*f)(
                 (((real_type) j) * step) + x_start, y_start,
                 n_coeffs,
                 power_order_n_x1, power_order_n_x2, power_order, 
                 map_2D_n_x1, map_2D_n_x2, map_2D
                );
#endif
    }

/\*
// BEGIN: Zero out ivp.
//    int n_ivp_bytes = (map_2D_order + 1) * sizeof(real_type);
//    memset(ivp, 0, n_ivp_bytes);
   int ivp_n_elems = map_2D_order + 1;
   for (int i = 0; i < ivp_n_elems; ++i)
       ivp[i] = 0.0;
// END: Zero out ivp.
*\/

   ivp[0] = MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, 0, 0);

   for (int i = 1; i < diff_2D_n_x1; ++i)    
   {
       for (int j = 0; j < diff_2D_n_x2 - i; ++j)    
       {
           MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i, j) =  
               MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i - 1, j + 1) - 
               MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i - 1, j);
       }
       ivp[i] = MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i, 0);
   }
}
*/

// NOTE: Inlined.
#define projection_remap_2D_scaling_CUDA_d_MACRO(x1_, x2_, lambda_pair_) \
{ \
    lambda_pair_ ## _2D[0] =  \
         MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 2, 0) + \
         ((x1_) * MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 2, 1)) + \
         ((x2_) * MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 2, 2)); \
    lambda_pair_ ## _2D[1] = \
         MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 3, 0) + \
         ((x1_) * MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 3, 1)) + \
         ((x2_) * MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 3, 2)); \
}

/*
// NOTE: Consider inlining this.
__device__ void projection_remap_2D_scaling_CUDA_d
(
    real_type x1, real_type x2,
    int map_2D_n_x1, int map_2D_n_x2, real_type *map_2D,
    real_type lambda_2D[2]
)
{
    lambda_2D[0] = 
         MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 2, 0) +
         x1 * MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 2, 1) +
         x2 * MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 2, 2);
    lambda_2D[1] =
         MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 3, 0) +
         x1 * MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 3, 1) +
         x2 * MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 3, 2);
}
*/

// Thread device code.  Uses finite difference scheme locally.
// NOTE: Multiple contiguous pixels per thread.
// NOTE: Pushes values from projection0 to projection.
// WARNING: n0_x_min == 1.0!
// WARNING: n0_y_min == 1.0!

//extern __device__ __shared__ real_type shared_memory[];

__global__ void projection_transform_2D_push_CUDA_d(int n_pixels_side_x, int n_pixels_side_y, int n_passes_x, int i_pass_x, int n_passes_y, int i_pass_y)
{
    int i0_y_pixel = (int) (((threadIdx.y + (blockIdx.y * blockDim.y)) * (n_pixels_side_y * n_passes_y)) + (n_pixels_side_y * i_pass_y));
    int i0_x_pixel = (int) (((threadIdx.x + (blockIdx.x * blockDim.x)) * (n_pixels_side_x * n_passes_x)) + (n_pixels_side_x * i_pass_x));

    //Print("i0_pixel_(%d, %d)\n", i0_x_pixel, i0_y_pixel);

    // WARNING: From _NVIDIA CUDA Programming Guide_: __syncthreads() is allowed in conditional code but only if the conditional evaluates identically across the entire thread block, otherwise the code execution is likely to hang or produce unintended side effects.

/*
    // WARNING: If thread index results in out-of-bounds memory access, return immediately.
    if ((i0_x_pixel >= n0_x) || (i0_y_pixel >= n0_y))
        return;
//    if (!((i0_x_pixel == 0) && (i0_y_pixel == 1)))
//        return;
*/

    real_type y0_pixel = (real_type) i0_y_pixel;
    real_type x0_pixel = (real_type) i0_x_pixel;

    real_type n0_y_min_local = n0_y_min + y0_pixel; // WARNING: For thread index to match pixel, n0_y_min == 1.0!
    // WARNING: From _NVIDIA CUDA Programming Guide_: __syncthreads() is allowed in conditional code but only if the conditional evaluates identically across the entire thread block, otherwise the code execution is likely to hang or produce unintended side effects.
    real_type n0_y_max_local = n0_y_min_local + ((real_type) n_pixels_side_y);
    //real_type n0_y_max_local = min(n0_y_min_local + ((real_type) n_pixels_side_y), (real_type) n0_y);

    real_type n0_x_min_local = n0_x_min + x0_pixel; // WARNING: For thread index to match pixel, n0_x_min == 1.0!
    // WARNING: From _NVIDIA CUDA Programming Guide_: __syncthreads() is allowed in conditional code but only if the conditional evaluates identically across the entire thread block, otherwise the code execution is likely to hang or produce unintended side effects.
    real_type n0_x_max_local = n0_x_min_local + ((real_type) n_pixels_side_x);
    //real_type n0_x_max_local = min(n0_x_min_local + ((real_type) n_pixels_side_x), (real_type) n0_x);

    //Print("n0_(xy)_min_local = (%.15e, %.15e)\n", n0_x_min_local, n0_x_min_local);
    //Print("n0_(xy)_max_local = (%.15e, %.15e)\n", n0_x_max_local, n0_x_max_local);

    // BEGIN: Transform 2D loop initialization. {

    // Current offset into shared_memory byte array.
    //size_t shared_memory_offset_current = 0;
    //unsigned int i_thread_block = threadIdx.x * threadIdx.y;
    //unsigned int n_threads_block = blockDim.x * blockDim.y;

    real_type lambda_2D[2];
    real_type lambda_del_x_2D[2];

    // BEGIN: q1 and q2. {

    //int ivp_n_elems = map_2D_order + 1; // NOTE: Local memory.
    //int ivp_n_elems = (map_2D_order + 1) * ((int) n_threads_block); // NOTE: Shared memory.
    ////size_t ivp_n_bytes = ivp_n_elems * sizeof(real_type); // NOTE: Shared memory.
    real_type q1[MAX_MAP_2D_ORDER_CUDA + 1]; // NOTE: Local memory.
    //real_type *q1 = (real_type *) shared_memory; // NOTE: Shared memory.
    ////real_type *q1 = (real_type *) &(shared_memory[shared_memory_offset_current]);
    ////shared_memory_offset_current += ivp_n_bytes;
    real_type q2[MAX_MAP_2D_ORDER_CUDA + 1]; // NOTE: Local memory.
    //real_type *q2 = &(q1[ivp_n_elems]); // NOTE: Shared memory.
    ////real_type *q2 = (real_type *) &(shared_memory[shared_memory_offset_current]); // NOTE: Shared memory.
    ////shared_memory_offset_current += ivp_n_bytes; // NOTE: Shared memory.

/*
    // BEGIN: DEBUG. {
    real_type *q1; // = NULL;
    cudaMalloc((void **) &q1, ivp_n_bytes);
    if (q1 == NULL)
        ABORT("Cannot acquire memory for real_type *q1.\n");
    real_type *q2; // = NULL;
    cudaMalloc((void **) &q2, ivp_n_bytes); // DEBUG.
    if (q2 == NULL)
        ABORT("Cannot acquire memory for real_type *q2.\n");
    // END: DEBUG. }
*/

    // END: q1 and q2. }

    // BEGIN: diff_2D {

    int diff_2D_n_x1 = map_2D_order + 1;
    int diff_2D_n_x2 = map_2D_order + 1;

    real_type diff_2D[(MAX_MAP_2D_ORDER_CUDA + 1) * (MAX_MAP_2D_ORDER_CUDA + 1)];

/*
    // BEGIN: DEBUG. {
    real_type *diff_2D; // = NULL;
    cudaMalloc((void **) &diff_2D, diff_2D_n_x1 * diff_2D_n_x2 * sizeof(real_type));
    if (diff_2D == NULL)
        ABORT("Cannot acquire memory for real_type *diff_2D.\n");
    // END: DEBUG. }
*/

    // END: diff_2D }

    // END: Transform 2D loop initialization. }

    //#pragma unroll 7
    for (real_type y0 = n0_y_min_local; y0 < n0_y_max_local; y0 += step) // NOTE: In ROI image coordinates.
    {

        int y0_floor = FLOOR_INNER_LOOP(y0);
        real_type y0_alpha = y0 - (real_type) y0_floor;
        --y0_floor; // NOTE: Conforms to MATLAB prototype.

/*
#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
        projection_remap_2D_scaling_CUDA_d(y0, n0_x_min_local, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling_CUDA_d(y0, n0_x_min_local + step, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#else
        projection_remap_2D_scaling_CUDA_d(n0_x_min_local, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling_CUDA_d(n0_x_min_local + step, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#endif
*/

#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
        projection_remap_2D_scaling_CUDA_d_MACRO(y0, n0_x_min_local, lambda); // QUESTION: step or n0_x_min_local?
        projection_remap_2D_scaling_CUDA_d_MACRO(y0, n0_x_min_local + step, lambda_del_x); // QUESTION: 2.0 * step or n0_x_min_local + step?
#else
        projection_remap_2D_scaling_CUDA_d_MACRO(n0_x_min_local, y0, lambda); // QUESTION: step or n0_x_min_local?
        projection_remap_2D_scaling_CUDA_d_MACRO(n0_x_min_local + step, y0, lambda_del_x); // QUESTION: 2.0 * step or n0_x_min_local + step?
#endif

        lambda_del_x_2D[0] -= lambda_2D[0];
        lambda_del_x_2D[1] -= lambda_2D[1];

        //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);
        //Print("lambda_del_x_2D[0-1] = (%.15e, %.15e)\n", lambda_del_x_2D[0], lambda_del_x_2D[1]);
        //ABORT("Abort!\n");

/*
        // BEGIN: Device function version of projection_remap_2D_solve_ivp_CUDA_d().
        projection_remap_2D_solve_ivp_CUDA_d
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x1_CUDA_d,
            n0_x_min_local, y0, step, // QUESTION: step or n0_x_min_local?
            map_2D_order,
            diff_2D,
            q1
        );

        projection_remap_2D_solve_ivp_CUDA_d
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x2_CUDA_d,
            n0_x_min_local, y0, step, // QUESTION: step or n0_x_min_local?
            map_2D_order,
            diff_2D,
            q2
        );
        // END: Device function version of projection_remap_2D_solve_ivp_CUDA_d().
*/

        // BEGIN: Macro function version of projection_remap_2D_solve_ivp_CUDA_d().
        projection_remap_2D_solve_ivp_CUDA_d_MACRO
        (
            1,
            n0_x_min_local, y0,
            1
        );

        projection_remap_2D_solve_ivp_CUDA_d_MACRO
        (
            2,
            n0_x_min_local, y0,
            2
        );
        // END: Macro function version of projection_remap_2D_solve_ivp_CUDA_d().

        //print_matrix_1D("q1", q1, ivp_n_elems);
        //print_matrix_1D("q2", q2, ivp_n_elems);

        //#pragma unroll 7
        for (real_type x0 = n0_x_min_local; x0 < n0_x_max_local; x0 += step) // NOTE: In ROI image coordinates.
        {
//            t_0 = clock();

            int x0_floor = FLOOR_INNER_LOOP(x0);
            real_type x0_alpha = x0 - (real_type) x0_floor;
            --x0_floor; // NOTE: Conforms to MATLAB prototype.

#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
            real_type y = (q1[0] / lambda_2D[0]); // NOTE: Local memory.
            real_type x = (q2[0] / lambda_2D[1]); // NOTE: Local memory.
            //real_type y = (q1[0 + i_thread_block] / lambda_2D[0]); // NOTE: Shared memory.
            //real_type x = (q2[0 + i_thread_block] / lambda_2D[1]); // NOTE: Shared memory.
#else
            real_type x = (q1[0] / lambda_2D[0]); // NOTE: Local memory.
            real_type y = (q2[0] / lambda_2D[1]); // NOTE: Local memory.
            //real_type x = (q1[0 + i_thread_block] / lambda_2D[0]); // NOTE: Shared memory.
            //real_type y = (q2[0 + i_thread_block] / lambda_2D[1]); // NOTE: Shared memory.
#endif

            //Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
            //Print("(x, y) = (%.15e, %.15e)\n", x, y);

            int x_floor = FLOOR_INNER_LOOP(x);
            int y_floor = FLOOR_INNER_LOOP(y);

            real_type x_alpha = x - (real_type) x_floor;
            real_type y_alpha = y - (real_type) y_floor;
            --x_floor;
            --y_floor;

            for (int i_power = 0; i_power < map_2D_order; ++i_power)
            {
                q1[i_power] += q1[i_power + 1];
                q2[i_power] += q2[i_power + 1];
            }

/*
            // BEGIN: Shared memory.
            unsigned int ivp_elem_offset = 0;
            unsigned int ivp_next_elem_offset = n_threads_block;
            for (int i_power = 0; i_power < map_2D_order; ++i_power)
            {
                unsigned int ivp_elem_offset_thread = ivp_elem_offset + i_thread_block;
                unsigned int ivp_next_elem_offset_thread = ivp_next_elem_offset + i_thread_block;
                q1[ivp_elem_offset_thread] += q1[ivp_next_elem_offset_thread];
                q2[ivp_elem_offset_thread] += q2[ivp_next_elem_offset_thread];
                ivp_elem_offset += n_threads_block;
                ivp_next_elem_offset += n_threads_block;
            }
            // END: Shared memory.
*/

            lambda_2D[0] += lambda_del_x_2D[0];
            lambda_2D[1] += lambda_del_x_2D[1];

            //Print("x0_floor = %i\n", x0_floor);
            //Print("y0_floor = %i\n", y0_floor);

            //Print("x0_alpha = %.15e\n", x0_alpha);
            //Print("y0_alpha = %.15e\n", y0_alpha);

            __syncthreads();

            if (
                (0 <= x0_floor) && (x0_floor < n0_x - 1) &&
                (0 <= y0_floor) && (y0_floor < n0_y - 1) &&
                (0 <= x_floor) && (x_floor < n_x - 1) &&
                (0 <= y_floor) && (y_floor < n_y - 1)
               )
            {
                ////if ((x_floor == 0) && (y_floor == 0))
                //if (y0_floor == 2)
                //{
                //    Print("=====\n");
                //    Print("(x0, y0) = (%.18e, %.18e)\n", x0, y0);
                //    Print("(x0_floor, y0_floor) = (%i, %i)\n", x0_floor, y0_floor);
                //    Print("(x, y) = (%.18e, %.18e)\n", x, y);
                //    Print("(x_floor, y_floor) = (%i, %i)\n", x_floor, y_floor);
                //    Print("=====\n");
                //}

                // BEGIN: Geometry vs. pixel intensity. {
                // NOTE: Access source global memory inside loop: real_type v0 = INTERPOLATE_2D_LL_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor, x0_alpha, y0_alpha); // NOTE: Access source global memory inside loop.
                real_type v0 = INTERPOLATE_2D_LL_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor, x0_alpha, y0_alpha); // NOTE: Access source global memory inside loop.

                real_type prod_a = (1.0 - x_alpha) * (1.0 - y_alpha);
                real_type prod_b = x_alpha * (1.0 - y_alpha);
                real_type prod_c = (1.0 - x_alpha) * y_alpha;
                real_type prod_d = x_alpha * y_alpha;

                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor, y_floor, prod_a * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor + 1, y_floor, prod_b * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor, y_floor + 1, prod_c * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor + 1, y_floor + 1, prod_d * v0); 

                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor, y_floor, prod_a); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor + 1, y_floor, prod_b); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor, y_floor + 1, prod_c); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor + 1, y_floor + 1, prod_d); 
                // END: Geometry vs. pixel intensity. }
            }

            __syncthreads();
        }
    }
}

/*
// Thread device code.  One thread per pixel in projection0 spaced from next thread by offset_(xy).  Uses finite difference scheme locally.
// NOTE: One pixel per thread.
// NOTE: Pushes values from projection0 to projection.
// WARNING: n0_x_min == 1.0!
// WARNING: n0_y_min == 1.0!

//extern __device__ __shared__ real_type shared_memory[];

__global__ void projection_transform_2D_push_CUDA_d(int n_passes_x, int n_passes_y, int i_pass_x, int i_pass_y)
{
    int i0_y_pixel = (int) ((threadIdx.y + (blockIdx.y * blockDim.y)) * n_passes_y) + i_pass_y;
    int i0_x_pixel = (int) ((threadIdx.x + (blockIdx.x * blockDim.x)) * n_passes_x) + i_pass_x;

    // WARNING: If thread index results in out-of-bounds memory access, return immediately.
    if ((i0_x_pixel >= n0_x) || (i0_y_pixel >= n0_y))
        return;
//    if (!((i0_x_pixel == 0) && (i0_y_pixel == 1)))
//        return;

    //Print("i0_pixel_(%d, %d)\n", i0_x_pixel, i0_y_pixel);

    real_type y0_pixel = (real_type) i0_y_pixel;
    real_type x0_pixel = (real_type) i0_x_pixel;

    // NOTE: Push requires only current pixel and N, NE, and E neighbors.

    real_type n0_y_min_local = n0_y_min + y0_pixel; // WARNING: For thread index to match pixel, n0_y_min == 1.0!
    real_type n0_y_max_local = n0_y_min_local + 1.0;

    real_type n0_x_min_local = n0_x_min + x0_pixel; // WARNING: For thread index to match pixel, n0_x_min == 1.0!
    real_type n0_x_max_local = n0_x_min_local + 1.0;

    //Print("n0_(xy)_min_local = (%.15e, %.15e)\n", n0_x_min_local, n0_x_min_local);
    //Print("n0_(xy)_max_local = (%.15e, %.15e)\n", n0_x_max_local, n0_x_max_local);

/\*
    // BEGIN: Access source global memory outside loop.
    int x0_floor = FLOOR_INNER_LOOP(n0_x_min_local);
    int y0_floor = FLOOR_INNER_LOOP(n0_y_min_local);
    real_type v0_00 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor);
    real_type v0_10 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor);
    real_type v0_01 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1);
    real_type v0_11 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1);
    // END: Access source global memory outside loop.
*\/

    // BEGIN: Transform 2D loop initialization. {

    // Current offset into shared_memory byte array.
    //size_t shared_memory_offset_current = 0;
    //unsigned int i_thread_block = threadIdx.x * threadIdx.y;
    //unsigned int n_threads_block = blockDim.x * blockDim.y;

    real_type lambda_2D[2];
    real_type lambda_del_x_2D[2];

    // BEGIN: q1 and q2. {

    //int ivp_n_elems = map_2D_order + 1; // NOTE: Local memory.
    //int ivp_n_elems = (map_2D_order + 1) * ((int) n_threads_block); // NOTE: Shared memory.
    ////size_t ivp_n_bytes = ivp_n_elems * sizeof(real_type); // NOTE: Shared memory.
    real_type q1[MAX_MAP_2D_ORDER_CUDA + 1]; // NOTE: Local memory.
    //real_type *q1 = (real_type *) shared_memory; // NOTE: Shared memory.
    ////real_type *q1 = (real_type *) &(shared_memory[shared_memory_offset_current]);
    ////shared_memory_offset_current += ivp_n_bytes;
    real_type q2[MAX_MAP_2D_ORDER_CUDA + 1]; // NOTE: Local memory.
    //real_type *q2 = &(q1[ivp_n_elems]); // NOTE: Shared memory.
    ////real_type *q2 = (real_type *) &(shared_memory[shared_memory_offset_current]); // NOTE: Shared memory.
    ////shared_memory_offset_current += ivp_n_bytes; // NOTE: Shared memory.

/\*
    // BEGIN: DEBUG. {
    real_type *q1; // = NULL;
    cudaMalloc((void **) &q1, ivp_n_bytes);
    if (q1 == NULL)
        ABORT("Cannot acquire memory for real_type *q1.\n");
    real_type *q2; // = NULL;
    cudaMalloc((void **) &q2, ivp_n_bytes); // DEBUG.
    if (q2 == NULL)
        ABORT("Cannot acquire memory for real_type *q2.\n");
    // END: DEBUG. }
*\/

    // END: q1 and q2. }

    // BEGIN: diff_2D {

    int diff_2D_n_x1 = map_2D_order + 1;
    int diff_2D_n_x2 = map_2D_order + 1;

    real_type diff_2D[(MAX_MAP_2D_ORDER_CUDA + 1) * (MAX_MAP_2D_ORDER_CUDA + 1)];

/\*
    // BEGIN: DEBUG. {
    real_type *diff_2D; // = NULL;
    cudaMalloc((void **) &diff_2D, diff_2D_n_x1 * diff_2D_n_x2 * sizeof(real_type));
    if (diff_2D == NULL)
        ABORT("Cannot acquire memory for real_type *diff_2D.\n");
    // END: DEBUG. }
*\/

    // END: diff_2D }

    // END: Transform 2D loop initialization. }

    //#pragma unroll 4
    for (real_type y0 = n0_y_min_local; y0 < n0_y_max_local; y0 += step) // NOTE: In ROI image coordinates.
    {

        int y0_floor = FLOOR_INNER_LOOP(y0);
        real_type y0_alpha = y0 - (real_type) y0_floor;
        --y0_floor; // NOTE: Conforms to MATLAB prototype.

/\*
#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
        projection_remap_2D_scaling_CUDA_d(y0, n0_x_min_local, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling_CUDA_d(y0, n0_x_min_local + step, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#else
        projection_remap_2D_scaling_CUDA_d(n0_x_min_local, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling_CUDA_d(n0_x_min_local + step, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#endif
*\/

#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
        projection_remap_2D_scaling_CUDA_d_MACRO(y0, n0_x_min_local, lambda); // QUESTION: step or n0_x_min_local?
        projection_remap_2D_scaling_CUDA_d_MACRO(y0, n0_x_min_local + step, lambda_del_x); // QUESTION: 2.0 * step or n0_x_min_local + step?
#else
        projection_remap_2D_scaling_CUDA_d_MACRO(n0_x_min_local, y0, lambda); // QUESTION: step or n0_x_min_local?
        projection_remap_2D_scaling_CUDA_d_MACRO(n0_x_min_local + step, y0, lambda_del_x); // QUESTION: 2.0 * step or n0_x_min_local + step?
#endif

        lambda_del_x_2D[0] -= lambda_2D[0];
        lambda_del_x_2D[1] -= lambda_2D[1];

        //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);
        //Print("lambda_del_x_2D[0-1] = (%.15e, %.15e)\n", lambda_del_x_2D[0], lambda_del_x_2D[1]);
        //ABORT("Abort!\n");

/\*
        // BEGIN: Device function version of projection_remap_2D_solve_ivp_CUDA_d().
        projection_remap_2D_solve_ivp_CUDA_d
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x1_CUDA_d,
            n0_x_min_local, y0, step, // QUESTION: step or n0_x_min_local?
            map_2D_order,
            diff_2D,
            q1
        );

        projection_remap_2D_solve_ivp_CUDA_d
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x2_CUDA_d,
            n0_x_min_local, y0, step, // QUESTION: step or n0_x_min_local?
            map_2D_order,
            diff_2D,
            q2
        );
        // END: Device function version of projection_remap_2D_solve_ivp_CUDA_d().
*\/

        // BEGIN: Macro function version of projection_remap_2D_solve_ivp_CUDA_d().
        projection_remap_2D_solve_ivp_CUDA_d_MACRO
        (
            1,
            n0_x_min_local, y0,
            1
        );

        projection_remap_2D_solve_ivp_CUDA_d_MACRO
        (
            2,
            n0_x_min_local, y0,
            2
        );
        // END: Macro function version of projection_remap_2D_solve_ivp_CUDA_d().

        //print_matrix_1D("q1", q1, ivp_n_elems);
        //print_matrix_1D("q2", q2, ivp_n_elems);

        //#pragma unroll 4
        for (real_type x0 = n0_x_min_local; x0 < n0_x_max_local; x0 += step) // NOTE: In ROI image coordinates.
        {
//            t_0 = clock();

            int x0_floor = FLOOR_INNER_LOOP(x0);
            real_type x0_alpha = x0 - (real_type) x0_floor;
            --x0_floor; // NOTE: Conforms to MATLAB prototype.

#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
            real_type y = (q1[0] / lambda_2D[0]); // NOTE: Local memory.
            real_type x = (q2[0] / lambda_2D[1]); // NOTE: Local memory.
            //real_type y = (q1[0 + i_thread_block] / lambda_2D[0]); // NOTE: Shared memory.
            //real_type x = (q2[0 + i_thread_block] / lambda_2D[1]); // NOTE: Shared memory.
#else
            real_type x = (q1[0] / lambda_2D[0]); // NOTE: Local memory.
            real_type y = (q2[0] / lambda_2D[1]); // NOTE: Local memory.
            //real_type x = (q1[0 + i_thread_block] / lambda_2D[0]); // NOTE: Shared memory.
            //real_type y = (q2[0 + i_thread_block] / lambda_2D[1]); // NOTE: Shared memory.
#endif

            //Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
            //Print("(x, y) = (%.15e, %.15e)\n", x, y);

            int x_floor = FLOOR_INNER_LOOP(x);
            int y_floor = FLOOR_INNER_LOOP(y);

            real_type x_alpha = x - (real_type) x_floor;
            real_type y_alpha = y - (real_type) y_floor;
            --x_floor;
            --y_floor;

            for (int i_power = 0; i_power < map_2D_order; ++i_power)
            {
                q1[i_power] += q1[i_power + 1];
                q2[i_power] += q2[i_power + 1];
            }

/\*
            // BEGIN: Shared memory.
            unsigned int ivp_elem_offset = 0;
            unsigned int ivp_next_elem_offset = n_threads_block;
            for (int i_power = 0; i_power < map_2D_order; ++i_power)
            {
                unsigned int ivp_elem_offset_thread = ivp_elem_offset + i_thread_block;
                unsigned int ivp_next_elem_offset_thread = ivp_next_elem_offset + i_thread_block;
                q1[ivp_elem_offset_thread] += q1[ivp_next_elem_offset_thread];
                q2[ivp_elem_offset_thread] += q2[ivp_next_elem_offset_thread];
                ivp_elem_offset += n_threads_block;
                ivp_next_elem_offset += n_threads_block;
            }
            // END: Shared memory.
*\/

            lambda_2D[0] += lambda_del_x_2D[0];
            lambda_2D[1] += lambda_del_x_2D[1];

            //Print("x0_floor = %i\n", x0_floor);
            //Print("y0_floor = %i\n", y0_floor);

            //Print("x0_alpha = %.15e\n", x0_alpha);
            //Print("y0_alpha = %.15e\n", y0_alpha);

            if (
                (0 <= x0_floor) && (x0_floor < n0_x - 1) &&
                (0 <= y0_floor) && (y0_floor < n0_y - 1) &&
                (0 <= x_floor) && (x_floor < n_x - 1) &&
                (0 <= y_floor) && (y_floor < n_y - 1)
               )
            {
                //Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
                //Print("(x0_floor, y0_floor) = (%i, %i)\n", x0_floor, y0_floor);
                //Print("(x, y) = (%.15e, %.15e)\n", x, y);
                //Print("(x_floor, y_floor) = (%i, %i)\n", x_floor, y_floor);

                // BEGIN: Geometry vs. pixel intensity. {
                real_type v0 = INTERPOLATE_2D_LL_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor, x0_alpha, y0_alpha); // NOTE: Access source global memory inside loop.
/\*
                // BEGIN: Access source global memory outside loop.
                real_type v0 = ((v0_00 * (1.0 - x0_alpha) * (1.0 - y0_alpha)) +
                                (v0_10 * x0_alpha * (1.0 - y0_alpha)) +
                                (v0_01 * (1.0 - x0_alpha) * y0_alpha) +
                                (v0_11 * x0_alpha * y0_alpha));
                // END: Access source global memory outside loop.
*\/

                real_type prod_a = (1.0 - x_alpha) * (1.0 - y_alpha);
                real_type prod_b = x_alpha * (1.0 - y_alpha);
                real_type prod_c = (1.0 - x_alpha) * y_alpha;
                real_type prod_d = x_alpha * y_alpha;

                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor, y_floor, prod_a * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor + 1, y_floor, prod_b * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor, y_floor + 1, prod_c * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor + 1, y_floor + 1, prod_d * v0); 

                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor, y_floor, prod_a); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor + 1, y_floor, prod_b); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor, y_floor + 1, prod_c); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor + 1, y_floor + 1, prod_d); 
                // END: Geometry vs. pixel intensity. }
            }
        }
    }
}
*/

__global__ void projection_transform_2D_push_fill_CUDA_d()
{
    int i_y_pixel = (int) (threadIdx.y + (blockIdx.y * blockDim.y));
    int i_x_pixel = (int) (threadIdx.x + (blockIdx.x * blockDim.x));

    //Print("i_pixel_(%d, %d)\n", i_x_pixel, i_y_pixel);

    // WARNING: If thread index results in out-of-bounds memory access, return immediately.
    if ((i_x_pixel >= n_x) || (i_y_pixel >= n_y))
        return;

    //Print("hits_image_data(%i, %i) = %.15e\n", i_x_pixel, i_y_pixel, INDEX_2D(hits_image_data, n_x, i_x_pixel, i_y_pixel));

    __syncthreads();

    // Attempt to fill in a hole in the output image covered by this pixel.
    if (INDEX_2D(hits_image_data, n_x, i_x_pixel, i_y_pixel) > 0.0)
    {
        //Print("i_(xy)_pixel = (%i, %i)\n", i_x_pixel, i_y_pixel);
        INDEX_2D(projection_data, n_x_ws, i_x_pixel, i_y_pixel) =
            INDEX_2D(sums_image_data, n_x, i_x_pixel, i_y_pixel) /
            INDEX_2D(hits_image_data, n_x, i_x_pixel, i_y_pixel);
    }
    else if ((0 < i_x_pixel) && (i_x_pixel < n_x - 1) && (0 < i_y_pixel) && (i_y_pixel < n_y - 1))
    {
        //Print("i_(xy)_pixel (fill) = (%i, %i)\n", i_x_pixel, i_y_pixel);

        pixel_type hits_local = 0.0;
        pixel_type sum_local = 0.0;
        for (int i_x_local = -1; i_x_local <= 1; ++i_x_local)
        {
            for (int i_y_local = -1; i_y_local <= 1; ++i_y_local)
            {
                hits_local += INDEX_2D(hits_image_data, n_x, i_x_pixel + i_x_local, i_y_pixel + i_y_local);
                sum_local += INDEX_2D(sums_image_data, n_x, i_x_pixel + i_x_local, i_y_pixel + i_y_local);
            }
        }
        if (hits_local > 0.0)
            PUTVAL_2D(projection_data, n_x_ws, i_x_pixel, i_y_pixel, sum_local / hits_local);
        else
            PUTVAL_2D(projection_data, n_x_ws, i_x_pixel, i_y_pixel, 0.0);
    }

    __syncthreads();
}

/*
// WARNING: Failed attempt!  Abandon this code.
// Thread device code.  One thread per pixel in projection0 spaced from next thread by offset_(xy).  Uses finite difference scheme locally.
// NOTE: One pixel per thread.
// NOTE: Pulls values from projection to projection0.
// WARNING: n0_x_min == 1.0!
// WARNING: n0_y_min == 1.0!

//extern __device__ __shared__ real_type shared_memory[];

__global__ void projection_inv_transform_2D_pull_CUDA_d(int n_passes_x, int n_passes_y, int i_pass_x, int i_pass_y)
{
    int i0_y_pixel = (int) ((threadIdx.y + (blockIdx.y * blockDim.y)) * n_passes_y) + i_pass_y;
    int i0_x_pixel = (int) ((threadIdx.x + (blockIdx.x * blockDim.x)) * n_passes_x) + i_pass_x;

    // WARNING: From _NVIDIA CUDA Programming Guide_: __syncthreads() is allowed in conditional code but only if the conditional evaluates identically across the entire thread block, otherwise the code execution is likely to hang or produce unintended side effects.

/\*
    // WARNING: If thread index results in out-of-bounds memory access, return immediately.
    if ((i0_x_pixel >= n0_x) || (i0_y_pixel >= n0_y))
        return;
//    if (!((i0_x_pixel == 0) && (i0_y_pixel == 1)))
//        return;
*\/

    //Print("i0_pixel_(%d, %d)\n", i0_x_pixel, i0_y_pixel);

    real_type y0_pixel = (real_type) i0_y_pixel;
    real_type x0_pixel = (real_type) i0_x_pixel;

    // NOTE: Push requires only current pixel and N, NE, and E neighbors.

    real_type n0_y_min_local = n0_y_min + y0_pixel; // WARNING: For thread index to match pixel, n0_y_min == 1.0!
    real_type n0_y_max_local = n0_y_min_local + 1.0;

    real_type n0_x_min_local = n0_x_min + x0_pixel; // WARNING: For thread index to match pixel, n0_x_min == 1.0!
    real_type n0_x_max_local = n0_x_min_local + 1.0;

    //Print("n0_(xy)_min_local = (%.15e, %.15e)\n", n0_x_min_local, n0_x_min_local);
    //Print("n0_(xy)_max_local = (%.15e, %.15e)\n", n0_x_max_local, n0_x_max_local);

/\*
    // BEGIN: Access source global memory outside loop.
    int x0_floor = FLOOR_INNER_LOOP(n0_x_min_local);
    int y0_floor = FLOOR_INNER_LOOP(n0_y_min_local);
    real_type v0_00 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor);
    real_type v0_10 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor);
    real_type v0_01 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1);
    real_type v0_11 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1);
    // END: Access source global memory outside loop.
*\/

    // BEGIN: Transform 2D loop initialization. {

    // Current offset into shared_memory byte array.
    //size_t shared_memory_offset_current = 0;
    //unsigned int i_thread_block = threadIdx.x * threadIdx.y;
    //unsigned int n_threads_block = blockDim.x * blockDim.y;

    real_type lambda_2D[2];
    real_type lambda_del_x_2D[2];

    // BEGIN: q1 and q2. {

    //int ivp_n_elems = map_2D_order + 1; // NOTE: Local memory.
    //int ivp_n_elems = (map_2D_order + 1) * ((int) n_threads_block); // NOTE: Shared memory.
    ////size_t ivp_n_bytes = ivp_n_elems * sizeof(real_type); // NOTE: Shared memory.
    real_type q1[MAX_MAP_2D_ORDER_CUDA + 1]; // NOTE: Local memory.
    //real_type *q1 = (real_type *) shared_memory; // NOTE: Shared memory.
    ////real_type *q1 = (real_type *) &(shared_memory[shared_memory_offset_current]);
    ////shared_memory_offset_current += ivp_n_bytes;
    real_type q2[MAX_MAP_2D_ORDER_CUDA + 1]; // NOTE: Local memory.
    //real_type *q2 = &(q1[ivp_n_elems]); // NOTE: Shared memory.
    ////real_type *q2 = (real_type *) &(shared_memory[shared_memory_offset_current]); // NOTE: Shared memory.
    ////shared_memory_offset_current += ivp_n_bytes; // NOTE: Shared memory.

/\*
    // BEGIN: DEBUG. {
    real_type *q1; // = NULL;
    cudaMalloc((void **) &q1, ivp_n_bytes);
    if (q1 == NULL)
        ABORT("Cannot acquire memory for real_type *q1.\n");
    real_type *q2; // = NULL;
    cudaMalloc((void **) &q2, ivp_n_bytes); // DEBUG.
    if (q2 == NULL)
        ABORT("Cannot acquire memory for real_type *q2.\n");
    // END: DEBUG. }
*\/

    // END: q1 and q2. }

    // BEGIN: diff_2D {

    int diff_2D_n_x1 = map_2D_order + 1;
    int diff_2D_n_x2 = map_2D_order + 1;

    real_type diff_2D[(MAX_MAP_2D_ORDER_CUDA + 1) * (MAX_MAP_2D_ORDER_CUDA + 1)];

/\*
    // BEGIN: DEBUG. {
    real_type *diff_2D; // = NULL;
    cudaMalloc((void **) &diff_2D, diff_2D_n_x1 * diff_2D_n_x2 * sizeof(real_type));
    if (diff_2D == NULL)
        ABORT("Cannot acquire memory for real_type *diff_2D.\n");
    // END: DEBUG. }
*\/

    // END: diff_2D }

    // END: Transform 2D loop initialization. }

    //#pragma unroll 4
    for (real_type y0 = n0_y_min_local; y0 < n0_y_max_local; y0 += step) // NOTE: In ROI image coordinates.
    {

        int y0_floor = FLOOR_INNER_LOOP(y0);
        real_type y0_alpha = y0 - (real_type) y0_floor;
        --y0_floor; // NOTE: Conforms to MATLAB prototype.

/\*
#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
        projection_remap_2D_scaling_CUDA_d(y0, n0_x_min_local, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling_CUDA_d(y0, n0_x_min_local + step, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#else
        projection_remap_2D_scaling_CUDA_d(n0_x_min_local, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling_CUDA_d(n0_x_min_local + step, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#endif
*\/

#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
        projection_remap_2D_scaling_CUDA_d_MACRO(y0, n0_x_min_local, lambda); // QUESTION: step or n0_x_min_local?
        projection_remap_2D_scaling_CUDA_d_MACRO(y0, n0_x_min_local + step, lambda_del_x); // QUESTION: 2.0 * step or n0_x_min_local + step?
#else
        projection_remap_2D_scaling_CUDA_d_MACRO(n0_x_min_local, y0, lambda); // QUESTION: step or n0_x_min_local?
        projection_remap_2D_scaling_CUDA_d_MACRO(n0_x_min_local + step, y0, lambda_del_x); // QUESTION: 2.0 * step or n0_x_min_local + step?
#endif

        lambda_del_x_2D[0] -= lambda_2D[0];
        lambda_del_x_2D[1] -= lambda_2D[1];

        //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);
        //Print("lambda_del_x_2D[0-1] = (%.15e, %.15e)\n", lambda_del_x_2D[0], lambda_del_x_2D[1]);
        //ABORT("Abort!\n");

/\*
        // BEGIN: Device function version of projection_remap_2D_solve_ivp_CUDA_d().
        projection_remap_2D_solve_ivp_CUDA_d
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x1_CUDA_d,
            n0_x_min_local, y0, step, // QUESTION: step or n0_x_min_local?
            map_2D_order,
            diff_2D,
            q1
        );

        projection_remap_2D_solve_ivp_CUDA_d
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x2_CUDA_d,
            n0_x_min_local, y0, step, // QUESTION: step or n0_x_min_local?
            map_2D_order,
            diff_2D,
            q2
        );
        // END: Device function version of projection_remap_2D_solve_ivp_CUDA_d().
*\/

        // BEGIN: Macro function version of projection_remap_2D_solve_ivp_CUDA_d().
        projection_remap_2D_solve_ivp_CUDA_d_MACRO
        (
            1,
            n0_x_min_local, y0,
            1
        );

        projection_remap_2D_solve_ivp_CUDA_d_MACRO
        (
            2,
            n0_x_min_local, y0,
            2
        );
        // END: Macro function version of projection_remap_2D_solve_ivp_CUDA_d().

        //print_matrix_1D("q1", q1, ivp_n_elems);
        //print_matrix_1D("q2", q2, ivp_n_elems);

        //#pragma unroll 4
        for (real_type x0 = n0_x_min_local; x0 < n0_x_max_local; x0 += step) // NOTE: In ROI image coordinates.
        {
//            t_0 = clock();

            int x0_floor = FLOOR_INNER_LOOP(x0);
            real_type x0_alpha = x0 - (real_type) x0_floor;
            --x0_floor; // NOTE: Conforms to MATLAB prototype.

#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
            real_type y = (q1[0] / lambda_2D[0]); // NOTE: Local memory.
            real_type x = (q2[0] / lambda_2D[1]); // NOTE: Local memory.
            //real_type y = (q1[0 + i_thread_block] / lambda_2D[0]); // NOTE: Shared memory.
            //real_type x = (q2[0 + i_thread_block] / lambda_2D[1]); // NOTE: Shared memory.
#else
            real_type x = (q1[0] / lambda_2D[0]); // NOTE: Local memory.
            real_type y = (q2[0] / lambda_2D[1]); // NOTE: Local memory.
            //real_type x = (q1[0 + i_thread_block] / lambda_2D[0]); // NOTE: Shared memory.
            //real_type y = (q2[0 + i_thread_block] / lambda_2D[1]); // NOTE: Shared memory.
#endif

            //Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
            //Print("(x, y) = (%.15e, %.15e)\n", x, y);

            int x_floor = FLOOR_INNER_LOOP(x);
            int y_floor = FLOOR_INNER_LOOP(y);

            real_type x_alpha = x - (real_type) x_floor;
            real_type y_alpha = y - (real_type) y_floor;
            --x_floor;
            --y_floor;

            for (int i_power = 0; i_power < map_2D_order; ++i_power)
            {
                q1[i_power] += q1[i_power + 1];
                q2[i_power] += q2[i_power + 1];
            }

/\*
            // BEGIN: Shared memory.
            unsigned int ivp_elem_offset = 0;
            unsigned int ivp_next_elem_offset = n_threads_block;
            for (int i_power = 0; i_power < map_2D_order; ++i_power)
            {
                unsigned int ivp_elem_offset_thread = ivp_elem_offset + i_thread_block;
                unsigned int ivp_next_elem_offset_thread = ivp_next_elem_offset + i_thread_block;
                q1[ivp_elem_offset_thread] += q1[ivp_next_elem_offset_thread];
                q2[ivp_elem_offset_thread] += q2[ivp_next_elem_offset_thread];
                ivp_elem_offset += n_threads_block;
                ivp_next_elem_offset += n_threads_block;
            }
            // END: Shared memory.
*\/

            lambda_2D[0] += lambda_del_x_2D[0];
            lambda_2D[1] += lambda_del_x_2D[1];

            //Print("x0_floor = %i\n", x0_floor);
            //Print("y0_floor = %i\n", y0_floor);

            //Print("x0_alpha = %.15e\n", x0_alpha);
            //Print("y0_alpha = %.15e\n", y0_alpha);

            if (
                (0 <= x0_floor) && (x0_floor < n0_x - 1) &&
                (0 <= y0_floor) && (y0_floor < n0_y - 1) &&
                (0 <= x_floor) && (x_floor < n_x - 1) &&
                (0 <= y_floor) && (y_floor < n_y - 1)
               )
            {
                //Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
                //Print("(x0_floor, y0_floor) = (%i, %i)\n", x0_floor, y0_floor);
                //Print("(x, y) = (%.15e, %.15e)\n", x, y);
                //Print("(x_floor, y_floor) = (%i, %i)\n", x_floor, y_floor);

                // BEGIN: Geometry vs. pixel intensity. {
                real_type v0 = INTERPOLATE_2D_LL_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor, x0_alpha, y0_alpha); // NOTE: Access source global memory inside loop.
/\*
                // BEGIN: Access source global memory outside loop.
                real_type v0 = ((v0_00 * (1.0 - x0_alpha) * (1.0 - y0_alpha)) +
                                (v0_10 * x0_alpha * (1.0 - y0_alpha)) +
                                (v0_01 * (1.0 - x0_alpha) * y0_alpha) +
                                (v0_11 * x0_alpha * y0_alpha));
                // END: Access source global memory outside loop.
*\/

                real_type prod_a = (1.0 - x_alpha) * (1.0 - y_alpha);
                real_type prod_b = x_alpha * (1.0 - y_alpha);
                real_type prod_c = (1.0 - x_alpha) * y_alpha;
                real_type prod_d = x_alpha * y_alpha;

                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor, y_floor, prod_a * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor + 1, y_floor, prod_b * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor, y_floor + 1, prod_c * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_floor + 1, y_floor + 1, prod_d * v0); 

                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor, y_floor, prod_a); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor + 1, y_floor, prod_b); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor, y_floor + 1, prod_c); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_floor + 1, y_floor + 1, prod_d); 
                // END: Geometry vs. pixel intensity. }
            }
        }
    }
}
*/

// Thread device code.  Uses finite difference scheme locally.
// NOTE: Multiple contiguous pixels per thread.
// NOTE: Pushes values from projection to projection0.
// WARNING: n0_x_min == 1.0!
// WARNING: n0_y_min == 1.0!

//extern __device__ __shared__ real_type shared_memory[];

__global__ void projection_inv_transform_2D_pull_CUDA_d(int n_pixels_side_x, int n_pixels_side_y, int n_passes_x, int i_pass_x, int n_passes_y, int i_pass_y)
{
    //Print("(Enter) threadIdx.(xy) = (%u, %u), blockIdx.(xy) = (%u, %u)\n", threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);

    int i0_y_pixel = (int) (((threadIdx.y + (blockIdx.y * blockDim.y)) * (n_pixels_side_y * n_passes_y)) + (n_pixels_side_y * i_pass_y));
    int i0_x_pixel = (int) (((threadIdx.x + (blockIdx.x * blockDim.x)) * (n_pixels_side_x * n_passes_x)) + (n_pixels_side_x * i_pass_x));

    // WARNING: From _NVIDIA CUDA Programming Guide_: __syncthreads() is allowed in conditional code but only if the conditional evaluates identically across the entire thread block, otherwise the code execution is likely to hang or produce unintended side effects.

/*
    // WARNING: If thread index results in out-of-bounds memory access, return immediately.
    if ((i0_x_pixel >= n0_x) || (i0_y_pixel >= n0_y))
        return;
//    if (!((i0_x_pixel == 0) && (i0_y_pixel == 1)))
//        return;
*/

    //Print("i0_pixel_(%d, %d)\n", i0_x_pixel, i0_y_pixel);

    real_type y0_pixel = (real_type) i0_y_pixel;
    real_type x0_pixel = (real_type) i0_x_pixel;

    real_type n0_y_min_local = n0_y_min + y0_pixel; // WARNING: For thread index to match pixel, n0_y_min == 1.0!
    // WARNING: From _NVIDIA CUDA Programming Guide_: __syncthreads() is allowed in conditional code but only if the conditional evaluates identically across the entire thread block, otherwise the code execution is likely to hang or produce unintended side effects.
    real_type n0_y_max_local = n0_y_min_local + ((real_type) n_pixels_side_y);
    //real_type n0_y_max_local = min(n0_y_min_local + ((real_type) n_pixels_side_y), (real_type) n0_y);

    real_type n0_x_min_local = n0_x_min + x0_pixel; // WARNING: For thread index to match pixel, n0_x_min == 1.0!
    // WARNING: From _NVIDIA CUDA Programming Guide_: __syncthreads() is allowed in conditional code but only if the conditional evaluates identically across the entire thread block, otherwise the code execution is likely to hang or produce unintended side effects.
    real_type n0_x_max_local = n0_x_min_local + ((real_type) n_pixels_side_x);
    //real_type n0_x_max_local = min(n0_x_min_local + ((real_type) n_pixels_side_x), (real_type) n0_x);

    //Print("n0_(xy)_min_local = (%.15e, %.15e)\n", n0_x_min_local, n0_x_min_local);
    //Print("n0_(xy)_max_local = (%.15e, %.15e)\n", n0_x_max_local, n0_x_max_local);

    // BEGIN: Transform 2D loop initialization. {

    // Current offset into shared_memory byte array.
    //size_t shared_memory_offset_current = 0;
    //unsigned int i_thread_block = threadIdx.x * threadIdx.y;
    //unsigned int n_threads_block = blockDim.x * blockDim.y;

    real_type lambda_2D[2];
    real_type lambda_del_x_2D[2];

    // BEGIN: q1 and q2. {

    //int ivp_n_elems = map_2D_order + 1; // NOTE: Local memory.
    //int ivp_n_elems = (map_2D_order + 1) * ((int) n_threads_block); // NOTE: Shared memory.
    ////size_t ivp_n_bytes = ivp_n_elems * sizeof(real_type); // NOTE: Shared memory.
    real_type q1[MAX_MAP_2D_ORDER_CUDA + 1]; // NOTE: Local memory.
    //real_type *q1 = (real_type *) shared_memory; // NOTE: Shared memory.
    ////real_type *q1 = (real_type *) &(shared_memory[shared_memory_offset_current]);
    ////shared_memory_offset_current += ivp_n_bytes;
    real_type q2[MAX_MAP_2D_ORDER_CUDA + 1]; // NOTE: Local memory.
    //real_type *q2 = &(q1[ivp_n_elems]); // NOTE: Shared memory.
    ////real_type *q2 = (real_type *) &(shared_memory[shared_memory_offset_current]); // NOTE: Shared memory.
    ////shared_memory_offset_current += ivp_n_bytes; // NOTE: Shared memory.

/*
    // BEGIN: DEBUG. {
    real_type *q1; // = NULL;
    cudaMalloc((void **) &q1, ivp_n_bytes);
    if (q1 == NULL)
        ABORT("Cannot acquire memory for real_type *q1.\n");
    real_type *q2; // = NULL;
    cudaMalloc((void **) &q2, ivp_n_bytes); // DEBUG.
    if (q2 == NULL)
        ABORT("Cannot acquire memory for real_type *q2.\n");
    // END: DEBUG. }
*/

    // END: q1 and q2. }

    // BEGIN: diff_2D {

    int diff_2D_n_x1 = map_2D_order + 1;
    int diff_2D_n_x2 = map_2D_order + 1;

    real_type diff_2D[(MAX_MAP_2D_ORDER_CUDA + 1) * (MAX_MAP_2D_ORDER_CUDA + 1)];

/*
    // BEGIN: DEBUG. {
    real_type *diff_2D; // = NULL;
    cudaMalloc((void **) &diff_2D, diff_2D_n_x1 * diff_2D_n_x2 * sizeof(real_type));
    if (diff_2D == NULL)
        ABORT("Cannot acquire memory for real_type *diff_2D.\n");
    // END: DEBUG. }
*/

    // END: diff_2D }

    // END: Transform 2D loop initialization. }

    //#pragma unroll 7
    for (real_type y0 = n0_y_min_local; y0 < n0_y_max_local; y0 += step) // NOTE: In ROI image coordinates.
    {

        int y0_floor = FLOOR_INNER_LOOP(y0);
        real_type y0_alpha = y0 - (real_type) y0_floor;
        --y0_floor; // NOTE: Conforms to MATLAB prototype.

/*
#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
        projection_remap_2D_scaling_CUDA_d(y0, n0_x_min_local, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling_CUDA_d(y0, n0_x_min_local + step, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#else
        projection_remap_2D_scaling_CUDA_d(n0_x_min_local, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling_CUDA_d(n0_x_min_local + step, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#endif
*/

#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
        projection_remap_2D_scaling_CUDA_d_MACRO(y0, n0_x_min_local, lambda); // QUESTION: step or n0_x_min_local?
        projection_remap_2D_scaling_CUDA_d_MACRO(y0, n0_x_min_local + step, lambda_del_x); // QUESTION: 2.0 * step or n0_x_min_local + step?
#else
        projection_remap_2D_scaling_CUDA_d_MACRO(n0_x_min_local, y0, lambda); // QUESTION: step or n0_x_min_local?
        projection_remap_2D_scaling_CUDA_d_MACRO(n0_x_min_local + step, y0, lambda_del_x); // QUESTION: 2.0 * step or n0_x_min_local + step?
#endif

        lambda_del_x_2D[0] -= lambda_2D[0];
        lambda_del_x_2D[1] -= lambda_2D[1];

        //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);
        //Print("lambda_del_x_2D[0-1] = (%.15e, %.15e)\n", lambda_del_x_2D[0], lambda_del_x_2D[1]);
        //ABORT("Abort!\n");

/*
        // BEGIN: Device function version of projection_remap_2D_solve_ivp_CUDA_d().
        projection_remap_2D_solve_ivp_CUDA_d
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x1_CUDA_d,
            n0_x_min_local, y0, step, // QUESTION: step or n0_x_min_local?
            map_2D_order,
            diff_2D,
            q1
        );

        projection_remap_2D_solve_ivp_CUDA_d
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x2_CUDA_d,
            n0_x_min_local, y0, step, // QUESTION: step or n0_x_min_local?
            map_2D_order,
            diff_2D,
            q2
        );
        // END: Device function version of projection_remap_2D_solve_ivp_CUDA_d().
*/

        // BEGIN: Macro function version of projection_remap_2D_solve_ivp_CUDA_d().
        projection_remap_2D_solve_ivp_CUDA_d_MACRO
        (
            1,
            n0_x_min_local, y0,
            1
        );

        projection_remap_2D_solve_ivp_CUDA_d_MACRO
        (
            2,
            n0_x_min_local, y0,
            2
        );
        // END: Macro function version of projection_remap_2D_solve_ivp_CUDA_d().

        //print_matrix_1D("q1", q1, ivp_n_elems);
        //print_matrix_1D("q2", q2, ivp_n_elems);

        //#pragma unroll 7
        for (real_type x0 = n0_x_min_local; x0 < n0_x_max_local; x0 += step) // NOTE: In ROI image coordinates.
        {
//            t_0 = clock();

            int x0_floor = FLOOR_INNER_LOOP(x0);
            real_type x0_alpha = x0 - (real_type) x0_floor;
            --x0_floor; // NOTE: Conforms to MATLAB prototype.

#if PROJECTION_TRANSFORM_2D_CUDA_FLIP_XY
            real_type y = (q1[0] / lambda_2D[0]); // NOTE: Local memory.
            real_type x = (q2[0] / lambda_2D[1]); // NOTE: Local memory.
            //real_type y = (q1[0 + i_thread_block] / lambda_2D[0]); // NOTE: Shared memory.
            //real_type x = (q2[0 + i_thread_block] / lambda_2D[1]); // NOTE: Shared memory.
#else
            real_type x = (q1[0] / lambda_2D[0]); // NOTE: Local memory.
            real_type y = (q2[0] / lambda_2D[1]); // NOTE: Local memory.
            //real_type x = (q1[0 + i_thread_block] / lambda_2D[0]); // NOTE: Shared memory.
            //real_type y = (q2[0 + i_thread_block] / lambda_2D[1]); // NOTE: Shared memory.
#endif

            //Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
            //Print("(x, y) = (%.15e, %.15e)\n", x, y);

            int x_floor = FLOOR_INNER_LOOP(x);
            int y_floor = FLOOR_INNER_LOOP(y);

            real_type x_alpha = x - (real_type) x_floor;
            real_type y_alpha = y - (real_type) y_floor;
            --x_floor;
            --y_floor;

            for (int i_power = 0; i_power < map_2D_order; ++i_power)
            {
                q1[i_power] += q1[i_power + 1];
                q2[i_power] += q2[i_power + 1];
            }

/*
            // BEGIN: Shared memory.
            unsigned int ivp_elem_offset = 0;
            unsigned int ivp_next_elem_offset = n_threads_block;
            for (int i_power = 0; i_power < map_2D_order; ++i_power)
            {
                unsigned int ivp_elem_offset_thread = ivp_elem_offset + i_thread_block;
                unsigned int ivp_next_elem_offset_thread = ivp_next_elem_offset + i_thread_block;
                q1[ivp_elem_offset_thread] += q1[ivp_next_elem_offset_thread];
                q2[ivp_elem_offset_thread] += q2[ivp_next_elem_offset_thread];
                ivp_elem_offset += n_threads_block;
                ivp_next_elem_offset += n_threads_block;
            }
            // END: Shared memory.
*/

            lambda_2D[0] += lambda_del_x_2D[0];
            lambda_2D[1] += lambda_del_x_2D[1];

            //Print("x0_floor = %i\n", x0_floor);
            //Print("y0_floor = %i\n", y0_floor);

            //Print("x0_alpha = %.15e\n", x0_alpha);
            //Print("y0_alpha = %.15e\n", y0_alpha);

            __syncthreads();

            if (
                (0 <= x0_floor) && (x0_floor < n0_x - 1) &&
                (0 <= y0_floor) && (y0_floor < n0_y - 1) &&
                (0 <= x_floor) && (x_floor < n_x - 1) &&
                (0 <= y_floor) && (y_floor < n_y - 1)
               )
            {
                //Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
                //Print("(x0_floor, y0_floor) = (%i, %i)\n", x0_floor, y0_floor);
                //Print("(x, y) = (%.15e, %.15e)\n", x, y);
                //Print("(x_floor, y_floor) = (%i, %i)\n", x_floor, y_floor);

                // BEGIN: Geometry vs. pixel intensity. {
                real_type v = INTERPOLATE_2D_LL_PIX(projection_data, n_x_ws, x_floor, y_floor, x_alpha, y_alpha);

                real_type prod_a = (1.0 - x0_alpha) * (1.0 - y0_alpha);
                real_type prod_b = x0_alpha * (1.0 - y0_alpha);
                real_type prod_c = (1.0 - x0_alpha) * y0_alpha;
                real_type prod_d = x0_alpha * y0_alpha;

                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor, prod_a * v); 
                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor, prod_b * v); 
                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1, prod_c * v); 
                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1, prod_d * v); 
                // END: Geometry vs. pixel intensity. }
            }

            __syncthreads();
        }
    }

    //Print("(Exit) threadIdx.(xy) = (%u, %u), blockIdx.(xy) = (%u, %u)\n", threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);
}

////////////////////////////////////////////////////////////////////////////////
// END: transform_2D }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: rotate_2D {
////////////////////////////////////////////////////////////////////////////////

// NOTE: Accessing the GPGPU's global memory is computationally more expensive than boundary checking?
#define PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS 0

// Thread device code.  One thread per pixel in projection.  Compute-intensive.  Uses finite difference scheme locally.
// WARNING: n_x_min == 0.0!
// WARNING: n_y_min == 0.0!
__global__ void projection_rotate_2D_prototype_CUDA_d()
{
    int i_y_pixel = (int) (threadIdx.y + blockIdx.y * blockDim.y);
    int i_x_pixel = (int) (threadIdx.x + blockIdx.x * blockDim.x);

    //Print("i_pixel_(%d, %d)\n", i_x_pixel, i_y_pixel);

    // WARNING: If thread index results in out-of-bounds memory access, return immediately.
    // WARNING: From _NVIDIA CUDA Programming Guide_: __syncthreads() is allowed in conditional code but only if the conditional evaluates identically across the entire thread block, otherwise the code execution is likely to hang or produce unintended side effects.
    if ((i_x_pixel >= n_x) || (i_y_pixel >= n_y))
        return;

    real_type y_pixel = (real_type) i_y_pixel;
    real_type x_pixel = (real_type) i_x_pixel;

    real_type n_y_min_local = /*n_y_min + */((y_pixel - 1.0) + step); // WARNING: For thread index to match pixel, n_y_min == 0.0!
    real_type n_y_max_local = n_y_min_local + 2.0 - step;

    real_type n_x_min_local = /*n_x_min + */((x_pixel - 1.0) + step); // WARNING: For thread index to match pixel, n_x_min == 0.0!
    real_type n_x_max_local = n_x_min_local + 2.0 - step;

    //Print("n_(xy)_min_local = (%.15e, %.15e)\n", n_x_min_local, n_x_min_local);
    //Print("n_(xy)_max_local = (%.15e, %.15e)\n", n_x_max_local, n_x_max_local);

//    real_type y0_pixel = (((x_pixel - center_x) * rot_matrix10) + ((y_pixel - center_y) * rot_matrix11)) + center0_y;

    real_type y0_init_local = (((n_x_min_local - center_x) * rot_matrix10) + ((n_y_min_local - center_y) * rot_matrix11)) + center0_y;
    real_type y0 = y0_init_local;

//    real_type x0_pixel = (((x_pixel - center_x) * rot_matrix00) + ((y_pixel - center_y) * rot_matrix01)) + center0_x;

    real_type x0_init_local = (((n_x_min_local - center_x) * rot_matrix00) + ((n_y_min_local - center_y) * rot_matrix01)) + center0_x;

    //Print("(xy)0_init_local = (%.15e, %.15e)\n", x0_init_local, y0_init_local);

    pixel_type v_pixel = 0.0;

    //#pragma unroll 7
    for (real_type y = n_y_min_local; y < n_y_max_local; y += step) // NOTE: In ROI image coordinates.
    {

        int y_floor = FLOOR_INNER_LOOP(y);
        real_type y_alpha = y - (real_type) y_floor;
        if (y >= y_pixel) 
            y_alpha = 1.0 - y_alpha;

//        if (y_floor < 0) continue; // DEBUG.

        real_type x0 = x0_init_local;

        //#pragma unroll 7
        for (real_type x = n_x_min_local; x < n_x_max_local; x += step) // NOTE: In ROI image coordinates.
        {
//            t_0 = clock();

            int x_floor = FLOOR_INNER_LOOP(x);
            real_type x_alpha = x - (real_type) x_floor;
            if (x >= x_pixel) 
                x_alpha = 1.0 - x_alpha;

//            if (x_floor < 0) continue; // DEBUG.

            //Print("x_floor = %i\n", x_floor);
            //Print("y_floor = %i\n", y_floor);

            //Print("x_alpha = %.15e\n", x_alpha);
            //Print("y_alpha = %.15e\n", y_alpha);

            // BEGIN: Geometry vs. pixel intensity. {

            // In order to avoid yet another boundary check, we have padded projection0.
//            real_type v = INTERPOLATE_2D_LL_PIX(projection_data, n_x_ws, x_floor, y_floor, x_alpha, y_alpha) * sample_factor;
//            real_type v = INDEX_2D_PIX(projection_data, n_x_ws, x_floor, y_floor) * sample_factor;
            real_type v = 0.0;

            //Print("v = %.15e\n", v);

//            Print("x0 = %.15e\n", x0);
//            Print("y0 = %.15e\n", y0);

            int x0_floor = FLOOR_INNER_LOOP(x0);
            int y0_floor = FLOOR_INNER_LOOP(y0);

            real_type x0_alpha = x0 - (real_type) x0_floor;
//            if (x0 > x0_pixel) 
//                x0_alpha = 1.0 - x0_alpha;

            real_type y0_alpha = y0 - (real_type) y0_floor;
//            if (y0 > y0_pixel) 
//                y0_alpha = 1.0 - y0_alpha;

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Pre-boundary-test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            t_0 = clock();

            real_type v0 = 0.0;

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            uint_address_type pixel0_address00 = (uint_address_type) ADDRESS_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor);
            uint_address_type pixel0_address10 = (uint_address_type) ADDRESS_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor);
            uint_address_type pixel0_address01 = (uint_address_type) ADDRESS_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1);
            uint_address_type pixel0_address11 = (uint_address_type) ADDRESS_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1); 
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
//            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            if ((projection0_padded_data_min_address <= pixel0_address00) && (pixel0_address00 < projection0_padded_data_max_address))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor) * (1.0 - x0_alpha) * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
            else
            {
                Print("(xy)_floor = (%i, %i)\n", x_floor, y_floor);
                Print("(xy)0_floor = (%i, %i)\n", x0_floor, y0_floor);
                Print("projection0_padded_data_(min, max)_address = (%llu, %llu)\n", projection0_padded_data_min_address, projection0_padded_data_max_address);
                Print("pixel0_address00 = %llu\n", pixel0_address00);
                ABORT("Source projection bounds violated (0, 0).\n"); 
            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
//            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            if ((projection0_padded_data_min_address <= pixel0_address10) && (pixel0_address10 < projection0_padded_data_max_address))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor) * x0_alpha * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
            else
            {
                Print("(xy)_floor = (%i, %i)\n", x_floor, y_floor);
                Print("(xy)0_floor = (%i, %i)\n", x0_floor + 1, y0_floor);
                Print("projection0_padded_data_(min, max)_address = (%llu, %llu).\n", projection0_padded_data_min_address, projection0_padded_data_max_address);
                Print("pixel0_address10 = %llu\n", pixel0_address10);
                ABORT("Source projection bounds violated (1, 0).\n"); 
            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
//            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            if ((projection0_padded_data_min_address <= pixel0_address01) && (pixel0_address01 < projection0_padded_data_max_address))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1) * (1.0 - x0_alpha) * y0_alpha;
                v += v0;
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor + 1, v0);
            }
            else
            {
                Print("(xy)_floor = (%i, %i)\n", x_floor, y_floor);
                Print("(xy)0_floor = (%i, %i)\n", x0_floor, y0_floor + 1);
                Print("projection0_padded_data_(min, max)_address = (%llu, %llu).\n", projection0_padded_data_min_address, projection0_padded_data_max_address);
                Print("pixel0_address01 = %llu\n", pixel0_address01);
                ABORT("Source projection bounds violated (0, 1).\n");
            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
//            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            if ((projection0_padded_data_min_address <= pixel0_address11) && (pixel0_address11 < projection0_padded_data_max_address))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1) * x0_alpha * y0_alpha;
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor + 1, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
            else
            {
                Print("(xy)_floor = (%i, %i)\n", x_floor, y_floor);
                Print("(xy)0_floor = (%i, %i)\n", x0_floor + 1, y0_floor + 1);
                Print("projection0_padded_data_(min, max)_address = (%llu, %llu).\n", projection0_padded_data_min_address, projection0_padded_data_max_address);
                Print("pixel0_address11 = %llu\n", pixel0_address11);
                ABORT("Source projection bounds violated (1, 1).\n");
            }
#endif

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Boundary test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            ADDVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v * x_alpha * y_alpha);

            v_pixel += v * x_alpha * y_alpha;

            // END: Geometry vs. pixel intensity. }

            x0 += x0_del_x;
            y0 += y0_del_x;
        }

        x0_init_local += x0_del_y;
        y0_init_local += y0_del_y;
        y0 = y0_init_local;
    }

    //// NOTE: This should not be too costly.
    //if ((i_x_pixel < n_x) && (i_y_pixel < n_y))
    //    PUTVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v_pixel * sample_factor);

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
    uint_address_type pixel_address = (uint_address_type) ADDRESS_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel);
    if (!((projection_padded_data_min_address <= pixel_address) && (pixel_address < projection_padded_data_max_address)))
    {
        Print("projection_padded_data_(min, max)_address = (%llu, %llu)\n", projection_padded_data_min_address, projection0_padded_data_max_address);
        Print("pixel_address = %llu\n", pixel_address);
        ABORT("Destination projection bounds violated (%i, %i).\n", i_x_pixel, i_y_pixel);
    }
#endif

    PUTVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v_pixel * sample_factor);
}

/*
// Thread device code.  One thread per pixel in projection.  Compute-intensive.  Uses finite difference scheme locally.
// Uses texture reference bound to 2D CUDA array to implement linear interpolation on input image.
// WARNING: n_x_min == 0.0!
// WARNING: n_y_min == 0.0!

texture<float, 2, cudaReadModeElementType> projection0_texture;

__global__ void projection_rotate_2D_prototype_CUDA_d()
{
    int i_y_pixel = (int) (threadIdx.y + blockIdx.y * blockDim.y);
    int i_x_pixel = (int) (threadIdx.x + blockIdx.x * blockDim.x);

    //Print("i_pixel_(%d, %d)\n", i_x_pixel, i_y_pixel);

    // WARNING: If thread index results in out-of-bounds memory access, return immediately.
    if ((i_x_pixel >= n_x) || (i_y_pixel >= n_y))
        return;

    real_type y_pixel = (real_type) i_y_pixel;
    real_type x_pixel = (real_type) i_x_pixel;

    real_type n_y_min_local = /\*n_y_min + *\/((y_pixel - 1.0) + step); // WARNING: For thread index to match pixel, n_y_min == 0.0!
    real_type n_y_max_local = n_y_min_local + 2.0 - step;

    real_type n_x_min_local = /\*n_x_min + *\/((x_pixel - 1.0) + step); // WARNING: For thread index to match pixel, n_x_min == 0.0!
    real_type n_x_max_local = n_x_min_local + 2.0 - step;

    //Print("n_(xy)_min_local = (%.15e, %.15e)\n", n_x_min_local, n_x_min_local);
    //Print("n_(xy)_max_local = (%.15e, %.15e)\n", n_x_max_local, n_x_max_local);

//    real_type y0_pixel = (((x_pixel - center_x) * rot_matrix10) + ((y_pixel - center_y) * rot_matrix11)) + center0_y;

    real_type y0_init_local = (((n_x_min_local - center_x) * rot_matrix10) + ((n_y_min_local - center_y) * rot_matrix11)) + center0_y;
    real_type y0 = y0_init_local;

//    real_type x0_pixel = (((x_pixel - center_x) * rot_matrix00) + ((y_pixel - center_y) * rot_matrix01)) + center0_x;

    real_type x0_init_local = (((n_x_min_local - center_x) * rot_matrix00) + ((n_y_min_local - center_y) * rot_matrix01)) + center0_x;

    //Print("(xy)0_init_local = (%.15e, %.15e)\n", x0_init_local, y0_init_local);

    pixel_type v_pixel = 0.0;

    //#pragma unroll 7
    for (real_type y = n_y_min_local; y < n_y_max_local; y += step) // NOTE: In ROI image coordinates.
    {

        int y_floor = FLOOR_INNER_LOOP(y);
        real_type y_alpha = y - (real_type) y_floor;
        if (y >= y_pixel) 
            y_alpha = 1.0 - y_alpha;

        real_type x0 = x0_init_local;

        //#pragma unroll 7
        for (real_type x = n_x_min_local; x < n_x_max_local; x += step) // NOTE: In ROI image coordinates.
        {
//            t_0 = clock();

            int x_floor = FLOOR_INNER_LOOP(x);
            real_type x_alpha = x - (real_type) x_floor;
            if (x >= x_pixel) 
                x_alpha = 1.0 - x_alpha;

            //Print("x_floor = %i\n", x_floor);
            //Print("y_floor = %i\n", y_floor);

            //Print("x_alpha = %.15e\n", x_alpha);
            //Print("y_alpha = %.15e\n", y_alpha);

            // BEGIN: Geometry vs. pixel intensity. {

            // In order to avoid yet another boundary check, we have padded projection0.
//            real_type v = INTERPOLATE_2D_LL_PIX(projection_data, n_x_ws, x_floor, y_floor, x_alpha, y_alpha) * sample_factor;
//            real_type v = INDEX_2D_PIX(projection_data, n_x_ws, x_floor, y_floor) * sample_factor;
            //real_type v = 0.0;

            //Print("v = %.15e\n", v);

//            Print("x0 = %.15e\n", x0);
//            Print("y0 = %.15e\n", y0);

/\*
            int x0_floor = FLOOR_INNER_LOOP(x0);
            int y0_floor = FLOOR_INNER_LOOP(y0);

            real_type x0_alpha = x0 - (real_type) x0_floor;
//            if (x0 > x0_pixel) 
//                x0_alpha = 1.0 - x0_alpha;

            real_type y0_alpha = y0 - (real_type) y0_floor;
//            if (y0 > y0_pixel) 
//                y0_alpha = 1.0 - y0_alpha;
*\/

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Pre-boundary-test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            t_0 = clock();

            real_type x0_norm = x0 / (real_type) n0_x;
            real_type y0_norm = y0 / (real_type) n0_y;
            real_type v0 = tex2D(projection0_texture, x0_norm, y0_norm);
            //real_type v0 = tex2D(projection0_texture, x0, y0);

/\*
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor) * (1.0 - x0_alpha) * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor, y0_floor);
//                ABORT("Source projection bounds violated (0, 0).\n"); 
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor) * x0_alpha * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor + 1, y0_floor);
//                ABORT("Source projection bounds violated (1, 0).\n"); 
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1) * (1.0 - x0_alpha) * y0_alpha;
                v += v0;
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor + 1, v0);
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor, y0_floor + 1);
//                ABORT("Source projection bounds violated (0, 1).\n");
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1) * x0_alpha * y0_alpha;
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor + 1, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor + 1, y0_floor + 1);
//                ABORT("Source projection bounds violated (1, 1).\n");
//            }
#endif
*\/

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Boundary test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            ADDVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v * x_alpha * y_alpha);

            v_pixel += v0 * x_alpha * y_alpha;

            // END: Geometry vs. pixel intensity. }

            x0 += x0_del_x;
            y0 += y0_del_x;
        }

        x0_init_local += x0_del_y;
        y0_init_local += y0_del_y;
        y0 = y0_init_local;
    }

    //// NOTE: This should not be too costly.
    //if ((i_x_pixel < n_x) && (i_y_pixel < n_y))
    //    PUTVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v_pixel * sample_factor);

    PUTVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v_pixel * sample_factor);
}
*/

/*
// Thread device code.  One thread per pixel in projection.  Compute-intensive.  Uses finite difference scheme locally.
// Uses cudaMemcpyAsync().  
// WARNING: Generates "error: calling a host function from a __device__/__global__ function is only allowed in device emulation mode".
// WARNING: n_x_min == 0.0!
// WARNING: n_y_min == 0.0!
__global__ void projection_rotate_2D_prototype_CUDA_d()
{
    int i_y_pixel = (int) (threadIdx.y + blockIdx.y * blockDim.y);
    int i_x_pixel = (int) (threadIdx.x + blockIdx.x * blockDim.x);

    //Print("i_pixel_(%d, %d)\n", i_x_pixel, i_y_pixel);

    // WARNING: If thread index results in out-of-bounds memory access, return immediately.
    if ((i_x_pixel >= n_x) || (i_y_pixel >= n_y))
        return;

    real_type y_pixel = (real_type) i_y_pixel;
    real_type x_pixel = (real_type) i_x_pixel;

    real_type n_y_min_local = /\*n_y_min + *\/((y_pixel - 1.0) + step); // WARNING: For thread index to match pixel, n_y_min == 0.0!
    real_type n_y_max_local = n_y_min_local + 2.0 - step;

    real_type n_x_min_local = /\*n_x_min + *\/((x_pixel - 1.0) + step); // WARNING: For thread index to match pixel, n_x_min == 0.0!
    real_type n_x_max_local = n_x_min_local + 2.0 - step;

    //Print("n_(xy)_min_local = (%.15e, %.15e)\n", n_x_min_local, n_x_min_local);
    //Print("n_(xy)_max_local = (%.15e, %.15e)\n", n_x_max_local, n_x_max_local);

//    real_type y0_pixel = (((x_pixel - center_x) * rot_matrix10) + ((y_pixel - center_y) * rot_matrix11)) + center0_y;

    real_type y0_init_local = (((n_x_min_local - center_x) * rot_matrix10) + ((n_y_min_local - center_y) * rot_matrix11)) + center0_y;
    real_type y0 = y0_init_local;

//    real_type x0_pixel = (((x_pixel - center_x) * rot_matrix00) + ((y_pixel - center_y) * rot_matrix01)) + center0_x;

    real_type x0_init_local = (((n_x_min_local - center_x) * rot_matrix00) + ((n_y_min_local - center_y) * rot_matrix01)) + center0_x;

    //Print("(xy)0_init_local = (%.15e, %.15e)\n", x0_init_local, y0_init_local);

    pixel_type v_pixel = 0.0;

    //#pragma unroll 7
    for (real_type y = n_y_min_local; y < n_y_max_local; y += step) // NOTE: In ROI image coordinates.
    {

        int y_floor = FLOOR_INNER_LOOP(y);
        real_type y_alpha = y - (real_type) y_floor;
        if (y >= y_pixel) 
            y_alpha = 1.0 - y_alpha;

        real_type x0 = x0_init_local;

        //#pragma unroll 7
        for (real_type x = n_x_min_local; x < n_x_max_local; x += step) // NOTE: In ROI image coordinates.
        {
//            t_0 = clock();

            int x_floor = FLOOR_INNER_LOOP(x);
            real_type x_alpha = x - (real_type) x_floor;
            if (x >= x_pixel) 
                x_alpha = 1.0 - x_alpha;

            //Print("x_floor = %i\n", x_floor);
            //Print("y_floor = %i\n", y_floor);

            //Print("x_alpha = %.15e\n", x_alpha);
            //Print("y_alpha = %.15e\n", y_alpha);

            // BEGIN: Geometry vs. pixel intensity. {

            // In order to avoid yet another boundary check, we have padded projection0.
//            real_type v = INTERPOLATE_2D_LL_PIX(projection_data, n_x_ws, x_floor, y_floor, x_alpha, y_alpha) * sample_factor;
//            real_type v = INDEX_2D_PIX(projection_data, n_x_ws, x_floor, y_floor) * sample_factor;
            real_type v = 0.0;

            //Print("v = %.15e\n", v);

//            Print("x0 = %.15e\n", x0);
//            Print("y0 = %.15e\n", y0);

            int x0_floor = FLOOR_INNER_LOOP(x0);
            int y0_floor = FLOOR_INNER_LOOP(y0);

            real_type x0_alpha = x0 - (real_type) x0_floor;
//            if (x0 > x0_pixel) 
//                x0_alpha = 1.0 - x0_alpha;

            real_type y0_alpha = y0 - (real_type) y0_floor;
//            if (y0 > y0_pixel) 
//                y0_alpha = 1.0 - y0_alpha;

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Pre-boundary-test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            t_0 = clock();

            real_type v0 = 0.0;

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor) * (1.0 - x0_alpha) * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor, y0_floor);
//                ABORT("Source projection bounds violated (0, 0).\n"); 
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor) * x0_alpha * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor + 1, y0_floor);
//                ABORT("Source projection bounds violated (1, 0).\n"); 
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1) * (1.0 - x0_alpha) * y0_alpha;
                v += v0;
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor + 1, v0);
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor, y0_floor + 1);
//                ABORT("Source projection bounds violated (0, 1).\n");
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1) * x0_alpha * y0_alpha;
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor + 1, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor + 1, y0_floor + 1);
//                ABORT("Source projection bounds violated (1, 1).\n");
//            }
#endif

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Boundary test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            ADDVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v * x_alpha * y_alpha);

            v_pixel += v * x_alpha * y_alpha;

            // END: Geometry vs. pixel intensity. }

            x0 += x0_del_x;
            y0 += y0_del_x;
        }

        x0_init_local += x0_del_y;
        y0_init_local += y0_del_y;
        y0 = y0_init_local;
    }

    //// NOTE: This should not be too costly.
    //if ((i_x_pixel < n_x) && (i_y_pixel < n_y))
    //    PUTVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v_pixel * sample_factor);

    //PUTVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v_pixel * sample_factor);

    v_pixel *= sample_factor;
    // WARNING: Generates "error: calling a host function from a __device__/__global__ function is only allowed in device emulation mode".
    // WARNING: The documentation did __not__ make it clear that this is strictly a host function.
    cudaMemcpyAsync((projection_data + i_x_pixel + (i_y_pixel * n_x_ws)), &v_pixel, sizeof(pixel_type), cudaMemcpyDeviceToHost);
}
*/

/*
#define PROJECTION_ROTATE_2D_CUDA_PROJECTION_DATA_S_MAX_ELEMS 512

// Thread device code.  One thread per pixel in projection.  Compute-intensive.  Uses finite difference scheme locally.  Attempts to use shared memory as a cache.
// WARNING: n_x_min == 0.0!
// WARNING: n_y_min == 0.0!
__global__ void projection_rotate_2D_prototype_CUDA_d()
{
    int i_y_pixel = (int) (threadIdx.y + blockIdx.y * blockDim.y);
    int i_x_pixel = (int) (threadIdx.x + blockIdx.x * blockDim.x);

    //Print("i_pixel_(%d, %d)\n", i_x_pixel, i_y_pixel);

    // WARNING: If thread index results in out-of-bounds memory access, return immediately.
    if ((i_x_pixel >= n_x) || (i_y_pixel >= n_y))
        return;

    // NOTE: This is purposefully larger than the array required by the largest possible block configuration.
    __shared__ pixel_type projection_data_s[PROJECTION_ROTATE_2D_CUDA_PROJECTION_DATA_S_MAX_ELEMS];

    real_type y_pixel = (real_type) i_y_pixel;
    real_type x_pixel = (real_type) i_x_pixel;

    real_type n_y_min_local = /\*n_y_min + *\/((y_pixel - 1.0) + step); // WARNING: For thread index to match pixel, n_y_min == 0.0!
    real_type n_y_max_local = n_y_min_local + 2.0 - step;

    real_type n_x_min_local = /\*n_x_min + *\/((x_pixel - 1.0) + step); // WARNING: For thread index to match pixel, n_x_min == 0.0!
    real_type n_x_max_local = n_x_min_local + 2.0 - step;

    //Print("n_(xy)_min_local = (%.15e, %.15e)\n", n_x_min_local, n_x_min_local);
    //Print("n_(xy)_max_local = (%.15e, %.15e)\n", n_x_max_local, n_x_max_local);

//    real_type y0_pixel = (((x_pixel - center_x) * rot_matrix10) + ((y_pixel - center_y) * rot_matrix11)) + center0_y;

    real_type y0_init_local = (((n_x_min_local - center_x) * rot_matrix10) + ((n_y_min_local - center_y) * rot_matrix11)) + center0_y;
    real_type y0 = y0_init_local;

//    real_type x0_pixel = (((x_pixel - center_x) * rot_matrix00) + ((y_pixel - center_y) * rot_matrix01)) + center0_x;

    real_type x0_init_local = (((n_x_min_local - center_x) * rot_matrix00) + ((n_y_min_local - center_y) * rot_matrix01)) + center0_x;

    //Print("(xy)0_init_local = (%.15e, %.15e)\n", x0_init_local, y0_init_local);

    pixel_type v_pixel = 0.0;

    //#pragma unroll 7
    for (real_type y = n_y_min_local; y < n_y_max_local; y += step) // NOTE: In ROI image coordinates.
    {

        int y_floor = FLOOR_INNER_LOOP(y);
        real_type y_alpha = y - (real_type) y_floor;
        if (y >= y_pixel) 
            y_alpha = 1.0 - y_alpha;

        real_type x0 = x0_init_local;

        //#pragma unroll 7
        for (real_type x = n_x_min_local; x < n_x_max_local; x += step) // NOTE: In ROI image coordinates.
        {
//            t_0 = clock();

            int x_floor = FLOOR_INNER_LOOP(x);
            real_type x_alpha = x - (real_type) x_floor;
            if (x >= x_pixel) 
                x_alpha = 1.0 - x_alpha;

            //Print("x_floor = %i\n", x_floor);
            //Print("y_floor = %i\n", y_floor);

            //Print("x_alpha = %.15e\n", x_alpha);
            //Print("y_alpha = %.15e\n", y_alpha);

            // BEGIN: Geometry vs. pixel intensity. {

            // In order to avoid yet another boundary check, we have padded projection0.
//            real_type v = INTERPOLATE_2D_LL_PIX(projection_data, n_x_ws, x_floor, y_floor, x_alpha, y_alpha) * sample_factor;
//            real_type v = INDEX_2D_PIX(projection_data, n_x_ws, x_floor, y_floor) * sample_factor;
            real_type v = 0.0;

            //Print("v = %.15e\n", v);

//            Print("x0 = %.15e\n", x0);
//            Print("y0 = %.15e\n", y0);

            int x0_floor = FLOOR_INNER_LOOP(x0);
            int y0_floor = FLOOR_INNER_LOOP(y0);

            real_type x0_alpha = x0 - (real_type) x0_floor;
//            if (x0 > x0_pixel) 
//                x0_alpha = 1.0 - x0_alpha;

            real_type y0_alpha = y0 - (real_type) y0_floor;
//            if (y0 > y0_pixel) 
//                y0_alpha = 1.0 - y0_alpha;

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Pre-boundary-test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            t_0 = clock();

            real_type v0 = 0.0;

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor) * (1.0 - x0_alpha) * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor, y0_floor);
//                ABORT("Source projection bounds violated (0, 0).\n"); 
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor) * x0_alpha * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor + 1, y0_floor);
//                ABORT("Source projection bounds violated (1, 0).\n"); 
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1) * (1.0 - x0_alpha) * y0_alpha;
                v += v0;
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor + 1, v0);
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor, y0_floor + 1);
//                ABORT("Source projection bounds violated (0, 1).\n");
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1) * x0_alpha * y0_alpha;
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor + 1, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor + 1, y0_floor + 1);
//                ABORT("Source projection bounds violated (1, 1).\n");
//            }
#endif

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Boundary test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            ADDVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v * x_alpha * y_alpha);

            v_pixel += v * x_alpha * y_alpha;

            // END: Geometry vs. pixel intensity. }

            x0 += x0_del_x;
            y0 += y0_del_x;
        }

        x0_init_local += x0_del_y;
        y0_init_local += y0_del_y;
        y0 = y0_init_local;
    }

    //// NOTE: This should not be too costly.
    //if ((i_x_pixel < n_x) && (i_y_pixel < n_y))
    //    PUTVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v_pixel * sample_factor);

    //PUTVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v_pixel * sample_factor);

    projection_data_s[threadIdx.x + threadIdx.y * blockDim.x] = v_pixel;

    __syncthreads();

    // WARNING: Even without a batched write back to global device memory, this is no faster than having each individual thread write directly to global device memory.
}
*/

/*
// Thread device code.  One thread per pixel in projection.  Compute-intensive.  Does not use finite difference scheme locally.  Slower than code using finite difference scheme locally.
// WARNING: n_x_min == 0.0!
// WARNING: n_y_min == 0.0!
__global__ void projection_rotate_2D_prototype_CUDA_d()
{
    int i_y_pixel = (int) (threadIdx.y + blockIdx.y * blockDim.y);
    int i_x_pixel = (int) (threadIdx.x + blockIdx.x * blockDim.x);
    real_type y_pixel = (real_type) i_y_pixel;
    real_type x_pixel = (real_type) i_x_pixel;

    //Print("i_pixel_(%d, %d)\n", i_x_pixel, i_y_pixel);

    real_type n_y_min_local = /\*n_y_min + *\/((y_pixel - 1.0) + step); // WARNING: For thread index to match pixel, n_y_min == 0.0!
    real_type n_y_max_local = n_y_min_local + 2.0 - step;

    real_type n_x_min_local = /\*n_x_min + *\/((x_pixel - 1.0) + step); // WARNING: For thread index to match pixel, n_x_min == 0.0!
    real_type n_x_max_local = n_x_min_local + 2.0 - step;

    //Print("n_(xy)_min_local = (%.15e, %.15e)\n", n_x_min_local, n_x_min_local);
    //Print("n_(xy)_max_local = (%.15e, %.15e)\n", n_x_max_local, n_x_max_local);

//    real_type y0_pixel = (((x_pixel - center_x) * rot_matrix10) + ((y_pixel - center_y) * rot_matrix11)) + center0_y;

//    real_type y0_init_local = (((n_x_min_local - center_x) * rot_matrix10) + ((n_y_min_local - center_y) * rot_matrix11)) + center0_y;
//    real_type y0 = y0_init_local;

//    real_type x0_pixel = (((x_pixel - center_x) * rot_matrix00) + ((y_pixel - center_y) * rot_matrix01)) + center0_x;

//    real_type x0_init_local = (((n_x_min_local - center_x) * rot_matrix00) + ((n_y_min_local - center_y) * rot_matrix01)) + center0_x;

    //Print("(xy)0_init_local = (%.15e, %.15e)\n", x0_init_local, y0_init_local);

    pixel_type v_pixel = 0.0;

    //#pragma unroll 7
    for (real_type y = n_y_min_local; y < n_y_max_local; y += step) // NOTE: In ROI image coordinates.
    {

        int y_floor = FLOOR_INNER_LOOP(y);
        real_type y_alpha = y - (real_type) y_floor;
        if (y >= y_pixel) 
            y_alpha = 1.0 - y_alpha;

        //real_type x0 = x0_init_local;
        real_type y_trans_y = y - center_y;
        real_type y_trans_y_rot_x = rot_matrix01 * y_trans_y;
        real_type y_trans_y_rot_y = rot_matrix11 * y_trans_y;
        real_type y_trans_y_rot_x_trans_x = y_trans_y_rot_x + center0_x;
        real_type y_trans_y_rot_y_trans_y = y_trans_y_rot_y + center0_y;

        //#pragma unroll 7
        for (real_type x = n_x_min_local; x < n_x_max_local; x += step) // NOTE: In ROI image coordinates.
        {
//            t_0 = clock();

            int x_floor = FLOOR_INNER_LOOP(x);
            real_type x_alpha = x - (real_type) x_floor;
            if (x >= x_pixel) 
                x_alpha = 1.0 - x_alpha;

            //Print("x_floor = %i\n", x_floor);
            //Print("y_floor = %i\n", y_floor);

            //Print("x_alpha = %.15e\n", x_alpha);
            //Print("y_alpha = %.15e\n", y_alpha);

            // BEGIN: Geometry vs. pixel intensity. {

            // In order to avoid yet another boundary check, we have padded projection0.
//            real_type v = INTERPOLATE_2D_LL_PIX(projection_data, n_x_ws, x_floor, y_floor, x_alpha, y_alpha) * sample_factor;
//            real_type v = INDEX_2D_PIX(projection_data, n_x_ws, x_floor, y_floor) * sample_factor;
            real_type v = 0.0;

            //Print("v = %.15e\n", v);

            real_type x_trans = x - center_x;

            real_type x0 = y_trans_y_rot_x_trans_x + (rot_matrix00 * x_trans);
            real_type y0 = y_trans_y_rot_y_trans_y + (rot_matrix10 * x_trans);

//            Print("x0 = %.15e\n", x0);
//            Print("y0 = %.15e\n", y0);

            int x0_floor = FLOOR_INNER_LOOP(x0);
            int y0_floor = FLOOR_INNER_LOOP(y0);

            real_type x0_alpha = x0 - (real_type) x0_floor;
//            if (x0 > x0_pixel) 
//                x0_alpha = 1.0 - x0_alpha;

            real_type y0_alpha = y0 - (real_type) y0_floor;
//            if (y0 > y0_pixel) 
//                y0_alpha = 1.0 - y0_alpha;

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Pre-boundary-test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            t_0 = clock();

            real_type v0 = 0.0;

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor) * (1.0 - x0_alpha) * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor, y0_floor);
//                ABORT("Source projection bounds violated (0, 0).\n"); 
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor) * x0_alpha * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor + 1, y0_floor);
//                ABORT("Source projection bounds violated (1, 0).\n"); 
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1) * (1.0 - x0_alpha) * y0_alpha;
                v += v0;
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor + 1, v0);
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor, y0_floor + 1);
//                ABORT("Source projection bounds violated (0, 1).\n");
//            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1) * x0_alpha * y0_alpha;
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor + 1, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
//            else
//            {
//                Print("(%i, %i)\n", x_floor, y_floor);
//                Print("(%i, %i)\n", x0_floor + 1, y0_floor + 1);
//                ABORT("Source projection bounds violated (1, 1).\n");
//            }
#endif

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Boundary test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            ADDVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v * x_alpha * y_alpha);

            v_pixel += v * x_alpha * y_alpha;

            // END: Geometry vs. pixel intensity. }

            //x0 += x0_del_x;
            //y0 += y0_del_x;
        }

        //x0_init_local += x0_del_y;
        //y0_init_local += y0_del_y;
        //y0 = y0_init_local;
    }

    // NOTE: This should not be too costly.
    if ((i_x_pixel < n_x) && (i_y_pixel < n_y))
        PUTVAL_2D_PIX(projection_data, n_x_ws, i_x_pixel, i_y_pixel, v_pixel * sample_factor);
}
*/

/*
// Thread device code.  One thread per pixel row in projection.  Attempted to optimize for memory access.  Failed.  No faster.
__global__ void projection_rotate_2D_prototype_CUDA_d(unsigned int offset_row)
{
    unsigned int i_thread = threadIdx.y + blockIdx.y * blockDim.y;
    real_type i_row_real = (real_type) ((i_thread * 2) + offset_row);

    real_type y0_outer = y0_init + (i_row_real * (y0_del_y * ((real_type) n_subsamples_y)));
    real_type x0_outer = x0_init + (i_row_real * (x0_del_y * ((real_type) n_subsamples_x)));

    // NOTE: Same order as constant list.
    real_type x0_del_x_outer = x0_del_x * ((real_type) n_subsamples_x);
    real_type y0_del_x_outer = y0_del_x * ((real_type) n_subsamples_x);
    real_type x0_del_y_outer = x0_del_y * ((real_type) n_subsamples_y);
    real_type y0_del_y_outer = y0_del_y * ((real_type) n_subsamples_y);

    // NOTE: Both of the following are constant for the whole row.
    real_type n_y_min_inner = n_y_min + i_row_real; // NOTE: n_y_min not necessarily 0.0.
    real_type n_y_max_inner = n_y_min_inner + 1.0;

    for (real_type n_x_outer = n_x_min; n_x_outer < n_x_max; n_x_outer += 1.0) // NOTE: In ROI image coordinates.
    {
        real_type y0_init_inner = y0_outer;
        real_type x0_init_inner = x0_outer;

        for (real_type y = n_y_min_inner; y < n_y_max_inner; y += step) // NOTE: In ROI image coordinates.
        {

            int y_floor = FLOOR_INNER_LOOP(y);
            real_type y_alpha = y - (real_type) y_floor;

            real_type y0 = y0_init_inner;
            real_type x0 = x0_init_inner;

            //Print("i_thread = %u, (xy)0 = (%.15e, %.15e)\n", i_thread, x0, y0);
            //Print("i_thread = %u, (xy) = (%.15e, %.15e)\n", i_thread, n_x_min, y);

            real_type n_x_min_inner = n_x_outer;
            real_type n_x_max_inner = n_x_min_inner + 1.0;
            for (real_type x = n_x_min_inner; x < n_x_max_inner; x += step) // NOTE: In ROI image coordinates.
            {
//                t_0 = clock();

                int x_floor = FLOOR_INNER_LOOP(x);
                real_type x_alpha = x - (real_type) x_floor;

//                Print("x_floor = %i\n", x_floor);
//                Print("y_floor = %i\n", y_floor);

//                Print("x_alpha = %.15e\n", x_alpha);
//                Print("y_alpha = %.15e\n", y_alpha);

                // BEGIN: Geometry vs. pixel intensity. {

                // In order to avoid yet another boundary check, we have padded projection0.
//                real_type v = INTERPOLATE_2D_LL_PIX(projection_data, n_x_ws, x_floor, y_floor, x_alpha, y_alpha) * sample_factor;
//                real_type v = INDEX_2D_PIX(projection_data, n_x_ws, x_floor, y_floor) * sample_factor;
                real_type v = 0.0;

                //Print("v = %.15e\n", v);

//                Print("x0 = %.15e\n", x0);
//                Print("y0 = %.15e\n", y0);

                int x0_floor = FLOOR_INNER_LOOP(x0);
                int y0_floor = FLOOR_INNER_LOOP(y0);

                real_type x0_alpha = x0 - (real_type) x0_floor;
                real_type y0_alpha = y0 - (real_type) y0_floor;

//                t_1 = clock();
//                Print("%i %i\n", t_1, t_0);
//                Print("Pre-boundary-test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//                t_0 = clock();

                real_type v0 = 0.0;

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//                if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
                if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
                {
#endif
                    v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor) * (1.0 - x0_alpha) * (1.0 - y0_alpha);
                    v += v0;
//                    Print("(%i, %i) = %.15e\n", x0_floor, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
                }
                else
                {
                    Print("(%i, %i)\n", x_floor, y_floor);
                    Print("(%i, %i)\n", x0_floor, y0_floor);
                    ABORT("Source projection bounds violated (0, 0).\n"); 
                }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//                if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
                if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
                {
#endif
                    v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor) * x0_alpha * (1.0 - y0_alpha);
                    v += v0;
//                    Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
                }
                else
                {
                    Print("(%i, %i)\n", x_floor, y_floor);
                    Print("(%i, %i)\n", x0_floor + 1, y0_floor);
                    ABORT("Source projection bounds violated (1, 0).\n"); 
                }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//                if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
                if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
                {
#endif
                    v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1) * (1.0 - x0_alpha) * y0_alpha;
                    v += v0;
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//                    Print("(%i, %i) = %.15e\n", x0_floor, y0_floor + 1, v0);
                }
                else
                {
                    Print("(%i, %i)\n", x_floor, y_floor);
                    Print("(%i, %i)\n", x0_floor, y0_floor + 1);
                    ABORT("Source projection bounds violated (0, 1).\n");
                }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//                if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
                if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
                {
#endif
                    v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1) * x0_alpha * y0_alpha;
                    v += v0;
//                    Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor + 1, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
                }
                else
                {
                    Print("(%i, %i)\n", x_floor, y_floor);
                    Print("(%i, %i)\n", x0_floor + 1, y0_floor + 1);
                    ABORT("Source projection bounds violated (1, 1).\n");
                }
#endif

                v *= sample_factor;

//                t_1 = clock();
//                Print("%i %i\n", t_1, t_0);
//                Print("Boundary test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//                PUTVAL_2D_PIX(projection_data, n_x, x_floor, y_floor, v * (1.0 - x_alpha) * (1.0 - y_alpha) * sample_factor);
//                PUTVAL_2D_PIX(projection_data, n_x, x_floor + 1, y_floor, v * x_alpha * (1.0 - y_alpha) * sample_factor);
//                PUTVAL_2D_PIX(projection_data, n_x, x_floor, y_floor + 1, v * (1.0 - x_alpha) * y_alpha * sample_factor);
//                PUTVAL_2D_PIX(projection_data, n_x, x_floor + 1, y_floor + 1, v * x_alpha * y_alpha * sample_factor);
                ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor, y_floor, v * (1.0 - x_alpha) * (1.0 - y_alpha));
                ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor + 1, y_floor, v * x_alpha * (1.0 - y_alpha));
                ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor, y_floor + 1, v * (1.0 - x_alpha) * y_alpha);
                ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor + 1, y_floor + 1, v * x_alpha * y_alpha);

                // END: Geometry vs. pixel intensity. }

                x0 += x0_del_x;
                y0 += y0_del_x;
            }

            y0_init_inner += y0_del_y;
            x0_init_inner += x0_del_y;
            y0 = y0_init_inner;
        }
        y0_outer += y0_del_x_outer;
        x0_outer += x0_del_x_outer;
    }
}
*/

/*
// Thread device code.  One thread per pixel row in projection.
__global__ void projection_rotate_2D_prototype_CUDA_d(unsigned int offset_row)
{
    unsigned int i_thread = threadIdx.y + blockIdx.y * blockDim.y;
    real_type i_row_real = (real_type) ((i_thread * 2) + offset_row);

    real_type y0_init_local = y0_init + (i_row_real * (y0_del_y * ((real_type) n_subsamples_y)));
    real_type y0 = y0_init_local;

    real_type x0_init_local = x0_init + (i_row_real * (x0_del_y * ((real_type) n_subsamples_y)));

    real_type n_y_min_local = n_y_min + i_row_real;
    real_type n_y_max_local = n_y_min_local + 1.0;
    for (real_type y = n_y_min_local; y < n_y_max_local; y += step) // NOTE: In ROI image coordinates.
    {

        int y_floor = FLOOR_INNER_LOOP(y);
        real_type y_alpha = y - (real_type) y_floor;

        real_type x0 = x0_init_local;

        //Print("i_thread = %u, (xy)0 = (%.15e, %.15e)\n", i_thread, x0, y0);
        //Print("i_thread = %u, (xy) = (%.15e, %.15e)\n", i_thread, n_x_min, y);

        // Unroll this loop?  Make indexing integral?
        for (real_type x = n_x_min; x <= n_x_max; x += step) // NOTE: In ROI image coordinates.
        {
//            t_0 = clock();

            int x_floor = FLOOR_INNER_LOOP(x);
            real_type x_alpha = x - (real_type) x_floor;

//            Print("x_floor = %i\n", x_floor);
//            Print("y_floor = %i\n", y_floor);

//            Print("x_alpha = %.15e\n", x_alpha);
//            Print("y_alpha = %.15e\n", y_alpha);

            // BEGIN: Geometry vs. pixel intensity. {

            // In order to avoid yet another boundary check, we have padded projection0.
//            real_type v = INTERPOLATE_2D_LL_PIX(projection_data, n_x_ws, x_floor, y_floor, x_alpha, y_alpha) * sample_factor;
//            real_type v = INDEX_2D_PIX(projection_data, n_x_ws, x_floor, y_floor) * sample_factor;
            real_type v = 0.0;

            //Print("v = %.15e\n", v);

//            Print("x0 = %.15e\n", x0);
//            Print("y0 = %.15e\n", y0);

            int x0_floor = FLOOR_INNER_LOOP(x0);
            int y0_floor = FLOOR_INNER_LOOP(y0);

            real_type x0_alpha = x0 - (real_type) x0_floor;
            real_type y0_alpha = y0 - (real_type) y0_floor;

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Pre-boundary-test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            t_0 = clock();

            real_type v0 = 0.0;

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor) * (1.0 - x0_alpha) * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
            else
            {
                Print("(%i, %i)\n", x_floor, y_floor);
                Print("(%i, %i)\n", x0_floor, y0_floor);
                ABORT("Source projection bounds violated (0, 0).\n"); 
            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor) * x0_alpha * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
            else
            {
                Print("(%i, %i)\n", x_floor, y_floor);
                Print("(%i, %i)\n", x0_floor + 1, y0_floor);
                ABORT("Source projection bounds violated (1, 0).\n"); 
            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1) * (1.0 - x0_alpha) * y0_alpha;
                v += v0;
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor + 1, v0);
            }
            else
            {
                Print("(%i, %i)\n", x_floor, y_floor);
                Print("(%i, %i)\n", x0_floor, y0_floor + 1);
                ABORT("Source projection bounds violated (0, 1).\n");
            }
#endif

#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1) * x0_alpha * y0_alpha;
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor + 1, v0);
#if PROJECTION_ROTATE_2D_CUDA_CHECK_BOUNDS
            }
            else
            {
                Print("(%i, %i)\n", x_floor, y_floor);
                Print("(%i, %i)\n", x0_floor + 1, y0_floor + 1);
                ABORT("Source projection bounds violated (1, 1).\n");
            }
#endif

            v *= sample_factor;

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Boundary test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            PUTVAL_2D_PIX(projection_data, n_x, x_floor, y_floor, v * (1.0 - x_alpha) * (1.0 - y_alpha) * sample_factor);
//            PUTVAL_2D_PIX(projection_data, n_x, x_floor + 1, y_floor, v * x_alpha * (1.0 - y_alpha) * sample_factor);
//            PUTVAL_2D_PIX(projection_data, n_x, x_floor, y_floor + 1, v * (1.0 - x_alpha) * y_alpha * sample_factor);
//            PUTVAL_2D_PIX(projection_data, n_x, x_floor + 1, y_floor + 1, v * x_alpha * y_alpha * sample_factor);
            ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor, y_floor, v * (1.0 - x_alpha) * (1.0 - y_alpha));
            ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor + 1, y_floor, v * x_alpha * (1.0 - y_alpha));
            ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor, y_floor + 1, v * (1.0 - x_alpha) * y_alpha);
            ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor + 1, y_floor + 1, v * x_alpha * y_alpha);

            // END: Geometry vs. pixel intensity. }

            x0 += x0_del_x;
            y0 += y0_del_x;
        }

        x0_init_local += x0_del_y;
        y0_init_local += y0_del_y;
        y0 = y0_init_local;
    }
}
*/

////////////////////////////////////////////////////////////////////////////////
// END: rotate_2D }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// END: Device code. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Host code. {
////////////////////////////////////////////////////////////////////////////////

// NOTE: Pushes values from projection0 to projection.
__host__ void projection_transform_2D_CUDA
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
    IplImage_base *sums_image, // No allocated data.
    IplImage_base *hits_image, // No allocated data.
    IplImage_base *projection_padded,
    IplImage_base *projection,
    real_type support[4]
)
{
// NOTE: Uses PROTOTYPE_COMPLIANT_INDEXING.
//#if !PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("!PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

    if ((projection0->depth != IPL_DEPTH_32F_BASE) && (projection0->depth != IPL_DEPTH_64F_BASE))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F_BASE) or %i (IPL_DEPTH_64F_BASE), but image->depth = %i.\n",
            IPL_DEPTH_32F_BASE, IPL_DEPTH_64F_BASE, projection0->depth);

    if ((sums_image->depth != IPL_DEPTH_32F_BASE) && (sums_image->depth != IPL_DEPTH_64F_BASE))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F_BASE) or %i (IPL_DEPTH_64F_BASE), but image->depth = %i.\n",
            IPL_DEPTH_32F_BASE, IPL_DEPTH_64F_BASE, sums_image->depth);

    if ((hits_image->depth != IPL_DEPTH_32F_BASE) && (hits_image->depth != IPL_DEPTH_64F_BASE))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F_BASE) or %i (IPL_DEPTH_64F_BASE), but image->depth = %i.\n",
            IPL_DEPTH_32F_BASE, IPL_DEPTH_64F_BASE, hits_image->depth);

    if ((projection->depth != IPL_DEPTH_32F_BASE) && (projection->depth != IPL_DEPTH_64F_BASE))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F_BASE) or %i (IPL_DEPTH_64F_BASE), but image->depth = %i.\n",
            IPL_DEPTH_32F_BASE, IPL_DEPTH_64F_BASE, projection->depth);

    //cvSetZero_base(sums_image);
    //cvSetZero_base(hits_image);
    cvSetZero_base(projection);

    cudaError_t cudaErrorCode;

    // BEGIN: projection0 {

    // BEGIN: Move projection0_padded to device, set ROI.

    // BEGIN: cudaMalloc() and cudaMallocHost(). {

    pixel_type *projection0_padded_data_h = (pixel_type *) projection0_padded->imageData;

    pixel_type *projection0_padded_data = NULL;

    cudaErrorCode = cudaMalloc((void **) &projection0_padded_data, projection0_padded->imageSize);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(projection0_padded_data, %i bytes) failed.\n", projection0_padded->imageSize);

    cudaErrorCode = cudaMemcpy(projection0_padded_data, projection0_padded_data_h, projection0_padded->imageSize, cudaMemcpyHostToDevice);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy(projection0_padded_data, projection0_padded_data, %i bytes) failed.\n", projection0_padded->imageSize);

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type).
    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / sizeof(pixel_type);
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.\n");

    // NOTE: The exterior of this image ROI is padded at the top and right with a row and column of zeros.
    uint_address_type projection0_elem_offset = calc_elem_offset_of_IplImage_base_ROI_widthStep(projection0_padded, projection0);
    pixel_type *projection0_data = projection0_padded_data + projection0_elem_offset;

    Print("projection0_padded_data = %llu\n", projection0_padded_data);
    Print("projection0_elem_offset = %llu\n", projection0_elem_offset);
    Print("projection0_data = %llu\n", projection0_data);

    // END: Move projection0_padded to device, set ROI.

    // END: cudaMalloc() and cudaMallocHost(). }

/*
    // BEGIN: cudaMallocPitch().

    // WARNING: Using data allocated by cudaMallocPitch() will require writing array indexing macros that use a widthStep in bytes and _not_ number of elements per row as is currently implemented.
    size_t projection0_padded_data_pitch;
    cudaErrorCode = cudaMallocPitch((void **) &projection0_padded_data, projection0_padded_pitch, projection0_padded->widthStep, projection0_padded->height);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMallocPitch(projection0_padded_data, (%i widthStep bytes) * (%i rows)) failed.\n", projection0_padded->widthStep, projection0_padded->height);

    // END: cudaMallocPitch().
*/

/*
    // BEGIN: CUDA array bound to texture. {

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type) and that projection0->widthStep / sizeof(pixel_type) == projection0->width.
    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / sizeof(pixel_type);
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.\n");
    //if (n0_x_ws != projection0->width)
    //    ABORT("(n0_x_ws == %i) != (projection0->width == %i)\n", n0_x_ws, projection0->width);

    // Allocate array and copy image data.
    cudaChannelFormatDesc projection0_channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

    cudaArray *projection0_array = NULL;

    cudaErrorCode = cudaMallocArray(&projection0_array, &projection0_channelDesc, projection0->width, projection0->height);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMallocArray(&projection0_array, &projection0_channelDesc, %i, %i) failed.\n", projection0->width, projection0->height);

    //cudaErrorCode = cudaMemcpyToArray(projection0_array, 0, 0, (void *) projection0->imageData, projection0->imageSize, cudaMemcpyHostToDevice);
    //if (cudaErrorCode != cudaSuccess)
    //   ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpyToArray(projection0_array, 0, 0, projection0_data, %i, %i, cudaMemcpyHostToDevice) failed.\n", projection0->width, projection0->height);

    cudaErrorCode = cudaMemcpy2DToArray(projection0_array, 0, 0, projection0->imageData, projection0->widthStep, projection0->width * sizeof(pixel_type), projection0->height, cudaMemcpyHostToDevice);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy2DToArray(projection0_array, 0, 0, projection0->imageData, %i, %i, %i, cudaMemcpyHostToDevice) failed.\n", projection0->widthStep, projection0->width * sizeof(pixel_type), projection0->height);

    // Set texture parameters.
    projection0_texture.addressMode[0] = cudaAddressModeWrap;
    projection0_texture.addressMode[1] = cudaAddressModeWrap;
//    projection0_texture.addressMode[0] = cudaAddressModeClamp;
//    projection0_texture.addressMode[1] = cudaAddressModeClamp;
    projection0_texture.filterMode = cudaFilterModeLinear;
    projection0_texture.normalized = true;

    // Bind the array to the texture.
    cudaBindTextureToArray(projection0_texture, projection0_array, projection0_channelDesc);

/\*
    // DEBUG: BEGIN: Copy CUDA array from device for later write to image.
    memset(projection0->imageData, 0, projection0->imageSize);

    //cudaErrorCode = cudaMemcpyFromArray((void *) projection0->imageData, projection0_array, 0, 0, projection0->imageSize, cudaMemcpyDeviceToHost);
    //if (cudaErrorCode != cudaSuccess)
    //   ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpyFromArray((void *) projection0->imageData, projection0_array, %i, cudaMemcpyDeviceToHost) failed.\n", projection0->imageSize);

    cudaErrorCode = cudaMemcpy2DFromArray((void *) projection0->imageData, projection0->widthStep, projection0_array, 0, 0, projection0->width * sizeof(pixel_type), projection0->height, cudaMemcpyDeviceToHost);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy2DFromArray((void *) projection0->imageData, %i, projection0_array, %i, %i, cudaMemcpyDeviceToHost) failed.\n", projection0->widthStep, projection0->width * sizeof(pixel_type), projection0->height);

    return;
    // DEBUG: END: Copy CUDA array from device for later write to image.
*\/

    // END: CUDA array bound to texture. }
*/

/*
    // BEGIN: cudaMallocHost() (WARNING: Asynchronous reads/writes are strictly host functions.  Dumb, dumb, dumb!). {

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type).
    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / sizeof(pixel_type);
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.\n");

    // Set ROI.
    pixel_type *projection0_data = (pixel_type *) projection0->imageData;

    // END: cudaMallocHost() (WARNING: Asynchronous reads/writes are strictly host functions.  Dumb, dumb, dumb!). }
*/

    int n0_x = projection0->width;
    int n0_y = projection0->height;

    real_type center0_x = support0[0];
    real_type center0_y = support0[1];

    // For support vs. projection check.
    real_type n0_x_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n0_x_max = (real_type) n0_x; // NOTE: Conforms to MATLAB prototype.
    real_type n0_y_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n0_y_max = (real_type) n0_y; // NOTE: Conforms to MATLAB prototype.

    // END: projection0 }

    // BEGIN: sums_image {

    // NOTE: This image is not padded and has the same extents as the unpadded projection.
    //pixel_type *sums_image_data_h = (pixel_type *) sums_image->imageData; // No data allocation.

    pixel_type *sums_image_data = NULL;

    cudaErrorCode = cudaMalloc((void **) &sums_image_data, sums_image->imageSize);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(sums_image_data, %i bytes) failed.\n", sums_image->imageSize);

    cudaErrorCode = cudaMemset(sums_image_data, 0, sums_image->imageSize);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemset(sums_image_data, 0, %i bytes) failed.\n", sums_image->imageSize);

    // END: sums_image }

    // BEGIN: hits_image {

    // NOTE: This image is not padded and has the same extents as the unpadded projection.
    //pixel_type *hits_image_data_h = (pixel_type *) hits_image->imageData; // No data allocation.

    pixel_type *hits_image_data = NULL;

    cudaErrorCode = cudaMalloc((void **) &hits_image_data, hits_image->imageSize);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(hits_image_data, %i bytes) failed.\n", hits_image->imageSize);

    cudaErrorCode = cudaMemset(hits_image_data, 0, hits_image->imageSize);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemset(hits_image_data, 0, %i bytes) failed.\n", hits_image->imageSize);

    // END: hits_image }

    // BEGIN: projection {

    // BEGIN: cudaMalloc(). {

    // BEGIN: Initialize projection on device, set ROI.

    pixel_type *projection_padded_data_h = (pixel_type *) projection_padded->imageData;

    pixel_type *projection_padded_data = NULL;
    cudaErrorCode = cudaMalloc((void **) &projection_padded_data, projection_padded->imageSize);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(projection_padded_data, %i bytes) failed.\n", projection_padded->imageSize);

    cudaErrorCode = cudaMemset(projection_padded_data, 0, projection_padded->imageSize);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMemset(projection_padded_data, 0, %i bytes) failed.\n", projection_padded->imageSize);

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type).
    real_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / sizeof(pixel_type);
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.\n");

    // NOTE: The exterior of this image ROI is padded at the top and right with a row and column of zeros.
    uint_address_type projection_elem_offset = calc_elem_offset_of_IplImage_base_ROI_widthStep(projection_padded, projection);
    pixel_type *projection_data = projection_padded_data + projection_elem_offset;

    Print("projection_padded_data = %llu\n", projection_padded_data);
    Print("projection_elem_offset = %llu\n", projection_elem_offset);
    Print("projection_data = %llu\n", projection_data);

    // END: Initialize projection on device, set ROI.

    // END: cudaMalloc(). }

/*
    // BEGIN: cudaMallocHost(). {

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type).
    pixel_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / sizeof(pixel_type);
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.\n");

    // Set ROI.
    pixel_type *projection_data = (pixel_type *) projection->imageData;

    // END: cudaMallocHost(). }
*/

    int n_x = projection->width;
    int n_y = projection->height;

    real_type center_x = support[0];
    real_type center_y = support[1];

    real_type n_x_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n_x_max = (real_type) n_x; // NOTE: Conforms to MATLAB prototype.
    real_type n_y_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n_y_max = (real_type) n_y; // NOTE: Conforms to MATLAB prototype.

    // END: projection }

    if ((n0_x_min < 0.0) || (n0_y_min < 0.0) || ((real_type) n0_x < n0_x_max) || ((real_type) n0_y < n0_y_max))
        ABORT("Support extends beyond projection0.\n"); 
    
    if ((n_x_min < 0.0) || (n_y_min < 0.0) || ((real_type) n_x < n_x_max) || ((real_type) n_y < n_y_max))
        ABORT("Support extends beyond projection.\n"); 

    // BEGIN: Unpack parameters and data. {

    // BEGIN: Unpack remap 2D parameters and data.

    real_type step_remap_2D = 1.0 / local_mag;

    // BEGIN: power_order {

    real_type *power_order = NULL;

    int n_power_order_bytes = power_order_n_x1 * power_order_n_x2 * sizeof(real_type);
    cudaErrorCode = cudaMalloc((void **) &power_order, n_power_order_bytes);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(power_order, %i bytes) failed.\n", n_power_order_bytes);

    cudaErrorCode = cudaMemcpy(power_order, power_order_h, n_power_order_bytes, cudaMemcpyHostToDevice);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy(power_order, power_order_h, %i bytes) failed.\n", n_power_order_bytes);

    // END: power_order }

    // BEGIN: map_2D {

    real_type *map_2D = NULL;

    // NOTE: Projection index is x3.
    int n_map_2D_bytes = map_2D_n_x1 * map_2D_n_x2 * sizeof(real_type);
    cudaErrorCode = cudaMalloc((void **) &map_2D, n_map_2D_bytes);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(map_2D, %i bytes) failed.\n", n_map_2D_bytes);

    cudaErrorCode = cudaMemcpy(map_2D, map_2D_current_h, n_map_2D_bytes, cudaMemcpyHostToDevice);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy(map_2D, map_2D_current_h, %i bytes) failed.\n", n_power_order_bytes);

    // END: map_2D }

    // END: Unpack remap 2D parameters and data.

    // BEGIN: Unpack rotation 2D parameters and data.

    real_type step_rotate_2D;
    if (test_rotation_matrix_2x2_for_subsampling(rot_matrix_2x2))
        step_rotate_2D = 0.25;
    else
        step_rotate_2D = 1.0;

    // END: Unpack rotation 2D parameters and data.

    // END: Unpack parameters and data. }

    real_type step = fmin(step_remap_2D, step_rotate_2D);
    //step = 1.0;  WARN("step is hardcoded to %.15e.\n", step);

    n0_x_min = step; // NOTE: Conforms to MATLAB prototype.
    n0_y_min = step; // NOTE: Conforms to MATLAB prototype.

    real_type n_subsamples_y_real = 1.0 / step;
    int n_subsamples_y = (int) floor(n_subsamples_y_real);
    if ((n_subsamples_y_real - ((real_type) n_subsamples_y)) != 0.0)
        ABORT("step = %.15e does not yield integral n_subsamples_y_real = %.15e.\n", step, n_subsamples_y_real);

    // NOTE: Number of subsamples in both directions is the same.
    int n_subsamples_x = (int) floor(n_subsamples_y_real);

    Print("step = %.15e\n", step);
    Print("n_subsamples_(xy) = (%d, %d)\n", n_subsamples_x, n_subsamples_y);

    //real_type sample_factor = step; // MATLAB prototype uses step for number of samples per pixel.
    real_type sample_factor = step * step;
    //real_type sample_factor = 1.0; WARN("sample_factor is hardcoded to %.15e.\n", sample_factor);

/*
    // BEGIN: rotation_2D {
 
    // NOTE: Use the inverse rotation matrix.
    real_type rot_matrix00 = rot_matrix_2x2[0][0];
    real_type rot_matrix01 = rot_matrix_2x2[1][0];
    real_type rot_matrix10 = rot_matrix_2x2[0][1];
    real_type rot_matrix11 = rot_matrix_2x2[1][1];

    real_type x0_init = 0.0;
    real_type y0_init = 0.0;

    real_type x_trans_x_init = n_x_min - center_x; // NOTE: From center of support in image coordinates, so from center of image.
    real_type x_trans_x_init_rot_x = rot_matrix00 * x_trans_x_init;
    real_type x_trans_x_init_rot_y = rot_matrix10 * x_trans_x_init;
    x0_init += x_trans_x_init_rot_x;
    y0_init += x_trans_x_init_rot_y;

    real_type y_trans_y_init = n_y_min - center_y; // NOTE: From center of support in image coordinates, so from center of image.
    real_type y_trans_y_init_rot_x = rot_matrix01 * y_trans_y_init;
    real_type y_trans_y_init_rot_y = rot_matrix11 * y_trans_y_init;
    x0_init += y_trans_y_init_rot_x;
    y0_init += y_trans_y_init_rot_y;

    x0_init += center0_x;
    y0_init += center0_y;

    real_type x0_del_x = rot_matrix00 * step;
    real_type y0_del_x = rot_matrix10 * step;

    real_type x0_del_y = rot_matrix01 * step;
    real_type y0_del_y = rot_matrix11 * step;

    // END: rotation_2D }
*/

    // BEGIN: Set constant GPU memory. {
    cudaMemcpyToSymbol("projection0_data", &projection0_data, sizeof(projection0_data)); // NOTE: cudaMalloc() and cudaMallocHost().
    cudaMemcpyToSymbol("n0_x", &n0_x, sizeof(n0_x));
    cudaMemcpyToSymbol("n0_y", &n0_y, sizeof(n0_y));
    cudaMemcpyToSymbol("n0_x_ws", &n0_x_ws, sizeof(n0_x_ws));
    cudaMemcpyToSymbol("center0_x", &center0_x, sizeof(center0_x));
    cudaMemcpyToSymbol("center0_y", &center0_y, sizeof(center0_y));
    cudaMemcpyToSymbol("n0_x_min", &n0_x_min, sizeof(n0_x_min));
    cudaMemcpyToSymbol("n0_x_max", &n0_x_max, sizeof(n0_x_max));
    cudaMemcpyToSymbol("n0_y_min", &n0_y_min, sizeof(n0_y_min));
    cudaMemcpyToSymbol("n0_y_max", &n0_y_max, sizeof(n0_y_max));
    cudaMemcpyToSymbol("projection_data", &projection_data, sizeof(projection_data));
    cudaMemcpyToSymbol("n_x", &n_x, sizeof(n_x));
    cudaMemcpyToSymbol("n_y", &n_y, sizeof(n_y));
    cudaMemcpyToSymbol("n_x_ws", &n_x_ws, sizeof(n_x_ws));
    cudaMemcpyToSymbol("center_x", &center_x, sizeof(center_x));
    cudaMemcpyToSymbol("center_y", &center_y, sizeof(center_y));
    cudaMemcpyToSymbol("n_x_min", &n_x_min, sizeof(n_x_min));
    cudaMemcpyToSymbol("n_x_max", &n_x_max, sizeof(n_x_max));
    cudaMemcpyToSymbol("n_y_min", &n_y_min, sizeof(n_y_min));
    cudaMemcpyToSymbol("n_y_max", &n_y_max, sizeof(n_y_max));
    cudaMemcpyToSymbol("step", &step, sizeof(step));
    cudaMemcpyToSymbol("n_subsamples_x", &n_subsamples_x, sizeof(n_subsamples_x));
    cudaMemcpyToSymbol("n_subsamples_y", &n_subsamples_y, sizeof(n_subsamples_y));
    cudaMemcpyToSymbol("sample_factor", &sample_factor, sizeof(sample_factor));
/*
    cudaMemcpyToSymbol("rot_matrix00", &rot_matrix00, sizeof(rot_matrix00));
    cudaMemcpyToSymbol("rot_matrix01", &rot_matrix01, sizeof(rot_matrix01));
    cudaMemcpyToSymbol("rot_matrix10", &rot_matrix10, sizeof(rot_matrix10));
    cudaMemcpyToSymbol("rot_matrix11", &rot_matrix11, sizeof(rot_matrix11));
    cudaMemcpyToSymbol("x0_init", &x0_init, sizeof(x0_init));
    cudaMemcpyToSymbol("y0_init", &y0_init, sizeof(y0_init));
    cudaMemcpyToSymbol("x0_del_x", &x0_del_x, sizeof(x0_del_x));
    cudaMemcpyToSymbol("y0_del_x", &y0_del_x, sizeof(y0_del_x));
    cudaMemcpyToSymbol("x0_del_y", &x0_del_y, sizeof(x0_del_y));
    cudaMemcpyToSymbol("y0_del_y", &y0_del_y, sizeof(y0_del_y));
*/
    cudaMemcpyToSymbol("map_2D_order", &map_2D_order, sizeof(map_2D_order));
    cudaMemcpyToSymbol("n_coeffs", &n_coeffs, sizeof(n_coeffs));
    cudaMemcpyToSymbol("power_order_n_x1", &power_order_n_x1, sizeof(power_order_n_x1));
    cudaMemcpyToSymbol("power_order_n_x2", &power_order_n_x2, sizeof(power_order_n_x2));
    cudaMemcpyToSymbol("power_order", &power_order, sizeof(power_order));
    cudaMemcpyToSymbol("map_2D_n_x1", &map_2D_n_x1, sizeof(map_2D_n_x1));
    cudaMemcpyToSymbol("map_2D_n_x2", &map_2D_n_x2, sizeof(map_2D_n_x2));
    cudaMemcpyToSymbol("map_2D", &map_2D, sizeof(map_2D));
    cudaMemcpyToSymbol("sums_image_data", &sums_image_data, sizeof(sums_image_data));
    cudaMemcpyToSymbol("hits_image_data", &hits_image_data, sizeof(hits_image_data));

    cudaErrorCode = cudaGetLastError();
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpyToSymbol() failed when setting constant GPU memory.\n");
    // END: Set constant GPU memory. }

/*
    // BEGIN: One-pass version. {

    // BEGIN: Make all of this part of future transform_2D_CUDA_params?
    unsigned int dimBlock_x = general_CUDA_params->dimBlock_x;
    unsigned int dimBlock_y = general_CUDA_params->dimBlock_y;
    dim3 dimBlock(dimBlock_x, dimBlock_y);
    unsigned int n_cols = (unsigned int) n_x;
    unsigned int n_rows = (unsigned int) n_y;
    unsigned int dimGrid_x = n_cols / dimBlock_x;
    if (n_cols % dimBlock_x != 0)
        ++dimGrid_x; 
    unsigned int dimGrid_y = n_rows / dimBlock_y;
    if (n_rows % dimBlock_y != 0)
        ++dimGrid_y; 
    dim3 dimGrid(dimGrid_x, dimGrid_y, 1u);
    // END: Make all of this part of future transform_2D_CUDA_params?

    Print("dimBlock.(xyz) = (%u, %u, %u)\n", dimBlock.x, dimBlock.y, dimBlock.z);
    Print("dimGrid.(xyz) = (%u, %u, %u)\n", dimGrid.x, dimGrid.y, dimGrid.z);
    //ABORT("ABORT!\n");

    projection_transform_2D_CUDA_d<<<dimGrid, dimBlock>>>();

    cudaErrorCode = cudaGetLastError();
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to projection_transform_2D_CUDA_d<<<dimGrid, dimBlock>>>() failed.\n");
   
    cudaErrorCode = cudaThreadSynchronize();
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaThreadSynchronize() failed.\n");

    // END: One-pass version. }
*/

    // BEGIN: N-pass version. {

    { // BEGIN: Push scope.
    
        // BEGIN: Make all of this part of future transform_2D_CUDA_params? {

/*
        // BEGIN: One pixel per thread.
        unsigned int dimBlock_x = general_CUDA_params->dimBlock_x;
        unsigned int dimBlock_y = general_CUDA_params->dimBlock_y;
        dim3 dimBlock(dimBlock_x, dimBlock_y);

        unsigned int n_passes_x = 1;
        unsigned int n_passes_y = 1;
        //unsigned int n_passes_x = 3;
        //unsigned int n_passes_y = 3;

        unsigned int n_cols_pass = ((unsigned int) n_x) / n_passes_x;
        unsigned int n_rows_pass = ((unsigned int) n_y) / n_passes_y;
    
        unsigned int dimGrid_x = n_cols_pass / dimBlock_x;
        if (n_cols_pass % dimBlock_x != 0)
            ++dimGrid_x; 
        unsigned int dimGrid_y = n_rows_pass / dimBlock_y;
        if (n_rows_pass % dimBlock_y != 0)
            ++dimGrid_y; 
        // END: One pixel per thread.
*/

        // BEGIN: Multiple pixels per thread.
        unsigned int dimBlock_x = 1;
        unsigned int dimBlock_y = 256;
        dim3 dimBlock(dimBlock_x, dimBlock_y);

        unsigned int n_passes_x = 1;
        unsigned int n_passes_y = 4;
        //unsigned int n_passes_y = 1;

        unsigned int n_cols_pass = ((unsigned int) n_x) / n_passes_x;
        if (n_cols_pass == 0)
            n_cols_pass = n_x;
        unsigned int n_rows_pass = ((unsigned int) n_y) / n_passes_y;
        if (n_rows_pass == 0)
            n_rows_pass = n_passes_y;

        unsigned int dimGrid_x = 1;
        unsigned int dimGrid_y = n_rows_pass / dimBlock_y;
        if (n_rows_pass % dimBlock_y != 0)
            ++dimGrid_y;
        // END: Multiple pixels per thread.

        // END: Make all of this part of future transform_2D_CUDA_params? }
    
        dim3 dimGrid(dimGrid_x, dimGrid_y);
    
        Print("dimBlock.(xyz) = (%u, %u, %u)\n", dimBlock.x, dimBlock.y, dimBlock.z);
        Print("dimGrid.(xyz) = (%u, %u, %u)\n", dimGrid.x, dimGrid.y, dimGrid.z);
        //ABORT("ABORT!\n");
    
        // BEGIN: Shared memory array. {
    
        unsigned int n_threads_block = dimBlock_x * dimBlock_y;
    
        size_t n_shared_memory_bytes = 0;
    
/*
        int ivp_n_elems = map_2D_order + 1; // NOTE: Local memory.
        //int ivp_n_elems = (map_2D_order + 1) * ((int) n_threads_block); // NOTE: Shared memory.
        size_t ivp_n_bytes = ivp_n_elems * sizeof(real_type);
        n_shared_memory_bytes += ivp_n_bytes * 2;
    
        int diff_2D_n_x1 = map_2D_order + 1;
        int diff_2D_n_x2 = map_2D_order + 1;
        int diff_2D_n_elems = diff_2D_n_x1 * diff_2D_n_x2; // NOTE: Local memory.
        //int diff_2D_n_elems = diff_2D_n_x1 * diff_2D_n_x2 * ((int) n_threads_block); // NOTE: Shared memory. // WARNING: Not enough shared memory.
        size_t diff_2D_n_bytes = diff_2D_n_elems * sizeof(real_type);
        n_shared_memory_bytes += diff_2D_n_bytes * n_threads_block;
*/
    
        Print("n_threads_block = %u\n", n_threads_block);
        Print("n_shared_memory_bytes = %u\n", n_shared_memory_bytes);
        //ABORT("ABORT!\n");
    
        // END: Shared memory array. }
    
        for (int i_pass_x = 0; i_pass_x < n_passes_x; ++i_pass_x)
        //for (int i_pass_x = 1; i_pass_x < 2; ++i_pass_x)
        {
            for (int i_pass_y = 0; i_pass_y < n_passes_y; ++i_pass_y)
            //for (int i_pass_y = 0; i_pass_y < 1; ++i_pass_y)
            {
                Print("i_pass_(xy) = (%i, %i)\n", i_pass_x, i_pass_y);

                // NOTE: One pixel per thread: projection_transform_2D_push_CUDA_d<<<dimGrid, dimBlock, n_shared_memory_bytes>>>(n_passes_x, n_passes_y, i_pass_x, i_pass_y); // NOTE: One pixel per thread.
                projection_transform_2D_push_CUDA_d<<<dimGrid, dimBlock, n_shared_memory_bytes>>>(n_cols_pass, 1, n_passes_x, i_pass_x, n_passes_y, i_pass_y); // NOTE: Multiple pixels per thread.
        
                cudaErrorCode = cudaGetLastError();
                if (cudaErrorCode != cudaSuccess)
                    ABORT_CUDA(cudaErrorCode, "Call to projection_transform_2D_push_CUDA_d<<<dimGrid, dimBlock, Ns>>>(%i, %i, %i, %i) failed.\n", n_passes_x, n_passes_y, i_pass_x, i_pass_y);
        
                cudaErrorCode = cudaThreadSynchronize();
                if (cudaErrorCode != cudaSuccess)
                    ABORT_CUDA(cudaErrorCode, "Call to cudaThreadSynchronize() failed.\n");
            }
        }
    
    } // END: Push scope.

    { // BEGIN: Push-fill scope.

        // BEGIN: Make all of this part of future transform_2D_CUDA_params? {
        unsigned int dimBlock_x = general_CUDA_params->dimBlock_x;
        unsigned int dimBlock_y = general_CUDA_params->dimBlock_y;
        dim3 dimBlock(dimBlock_x, dimBlock_y);
        unsigned int n_cols = (unsigned int) n_x;
        unsigned int n_rows = (unsigned int) n_y;
        unsigned int dimGrid_x = n_cols / dimBlock_x;
        if (n_cols % dimBlock_x != 0)
            ++dimGrid_x; 
        unsigned int dimGrid_y = n_rows / dimBlock_y;
        if (n_rows % dimBlock_y != 0)
            ++dimGrid_y; 
        dim3 dimGrid(dimGrid_x, dimGrid_y, 1u);
        // END: Make all of this part of future rotate_2D_CUDA_params? }

        Print("dimBlock.(xyz) = (%u, %u, %u)\n", dimBlock.x, dimBlock.y, dimBlock.z);
        Print("dimGrid.(xyz) = (%u, %u, %u)\n", dimGrid.x, dimGrid.y, dimGrid.z);
        //ABORT("ABORT!\n");

        projection_transform_2D_push_fill_CUDA_d<<<dimGrid, dimBlock>>>();
        
        cudaErrorCode = cudaGetLastError();
        if (cudaErrorCode != cudaSuccess)
            ABORT_CUDA(cudaErrorCode, "Call to projection_transform_2D_push_fill_CUDA_d<<<dimGrid, dimBlock, Ns>>>() failed.\n");
    
        cudaErrorCode = cudaThreadSynchronize();
        if (cudaErrorCode != cudaSuccess)
            ABORT_CUDA(cudaErrorCode, "Call to cudaThreadSynchronize() failed.\n");

    } // END: Push-fill scope.

    // END: N-pass version. }

    // BEGIN: cudaMalloc() and CUDA array bound to texture.
    // Copy padded projection from device to host.
    cudaErrorCode = cudaMemcpy(projection_padded_data_h, projection_padded_data, projection_padded->imageSize, cudaMemcpyDeviceToHost);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy(projection_data_h, projection_data, %i bytes) failed.\n", projection_padded->imageSize);

    cudaErrorCode = cudaFree(projection_padded_data);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaFree(projection_padded_data) failed.\n");
    // END: cudaMalloc() and CUDA array bound to texture.

    // BEGIN: cudaMalloc() and cudaMallocHost().
    cudaErrorCode = cudaFree(projection0_padded_data);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaFree(projection0_padded_data) failed.\n");
    // END: cudaMalloc() and cudaMallocHost().

/*
    // BEGIN: CUDA array bound to texture. 
    cudaErrorCode = cudaFreeArray(projection0_array);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaFreeArray(projection0_array) failed.\n");
    // END: CUDA array bound to texture.
*/

//    exit(0);
}

// NOTE: Pulls values from projection to projection0.
__host__ void projection_inv_transform_2D_CUDA
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
)
{
// NOTE: Uses PROTOTYPE_COMPLIANT_INDEXING.
//#if !PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("!PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

    if ((projection->depth != IPL_DEPTH_32F_BASE) && (projection->depth != IPL_DEPTH_64F_BASE))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F_BASE) or %i (IPL_DEPTH_64F_BASE), but image->depth = %i.\n",
            IPL_DEPTH_32F_BASE, IPL_DEPTH_64F_BASE, projection->depth);

    if ((projection0->depth != IPL_DEPTH_32F_BASE) && (projection0->depth != IPL_DEPTH_64F_BASE))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F_BASE) or %i (IPL_DEPTH_64F_BASE), but image->depth = %i.\n",
            IPL_DEPTH_32F_BASE, IPL_DEPTH_64F_BASE, projection0->depth);

    cvSetZero_base(projection0);

    cudaError_t cudaErrorCode;

    // BEGIN: projection {

    // BEGIN: cudaMalloc(). {

    // BEGIN: Initialize projection on device, set ROI.

    pixel_type *projection_padded_data_h = (pixel_type *) projection_padded->imageData;

    pixel_type *projection_padded_data = NULL;
    cudaErrorCode = cudaMalloc((void **) &projection_padded_data, projection_padded->imageSize);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(projection_padded_data, %i bytes) failed.\n", projection_padded->imageSize);

    cudaErrorCode = cudaMemcpy(projection_padded_data, projection_padded_data_h, projection_padded->imageSize, cudaMemcpyHostToDevice);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy(projection_padded_data, projection_padded_data, %i bytes) failed.\n", projection_padded->imageSize);

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type).
    real_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / sizeof(pixel_type);
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.\n");

    // NOTE: The exterior of this image ROI is padded at the top and right with a row and column of zeros.
    uint_address_type projection_elem_offset = calc_elem_offset_of_IplImage_base_ROI_widthStep(projection_padded, projection);
    pixel_type *projection_data = projection_padded_data + projection_elem_offset;

    Print("projection_padded_data = %llu\n", projection_padded_data);
    Print("projection_elem_offset = %llu\n", projection_elem_offset);
    Print("projection_data = %llu\n", projection_data);

    // END: Initialize projection on device, set ROI.

    // END: cudaMalloc(). }

    int n_x = projection->width;
    int n_y = projection->height;

    real_type center_x = support[0];
    real_type center_y = support[1];

    real_type n_x_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n_x_max = (real_type) n_x; // NOTE: Conforms to MATLAB prototype.
    real_type n_y_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n_y_max = (real_type) n_y; // NOTE: Conforms to MATLAB prototype.

    // END: projection }

    // BEGIN: projection0 {

    // BEGIN: Move projection0_padded to device, set ROI.

    // BEGIN: cudaMalloc() and cudaMallocHost(). {

    pixel_type *projection0_padded_data_h = (pixel_type *) projection0_padded->imageData;

    pixel_type *projection0_padded_data = NULL;

    cudaErrorCode = cudaMalloc((void **) &projection0_padded_data, projection0_padded->imageSize);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(projection0_padded_data, %i bytes) failed.\n", projection0_padded->imageSize);

    cudaErrorCode = cudaMemset(projection0_padded_data, 0, projection0_padded->imageSize);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMemset(projection0_padded_data, 0, %i bytes) failed.\n", projection0_padded->imageSize);

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type).
    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / sizeof(pixel_type);
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.\n");

    // NOTE: The exterior of this image ROI is padded at the top and right with a row and column of zeros.
    uint_address_type projection0_elem_offset = calc_elem_offset_of_IplImage_base_ROI_widthStep(projection0_padded, projection0);
    pixel_type *projection0_data = projection0_padded_data + projection0_elem_offset;

    Print("projection0_padded_data = %llu\n", projection0_padded_data);
    Print("projection0_elem_offset = %llu\n", projection0_elem_offset);
    Print("projection0_data = %llu\n", projection0_data);

    // END: Move projection0_padded to device, set ROI.

    // END: cudaMalloc() and cudaMallocHost(). }

    int n0_x = projection0->width;
    int n0_y = projection0->height;

    real_type center0_x = support0[0];
    real_type center0_y = support0[1];

    // For support vs. projection check.
    real_type n0_x_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n0_x_max = (real_type) n0_x; // NOTE: Conforms to MATLAB prototype.
    real_type n0_y_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n0_y_max = (real_type) n0_y; // NOTE: Conforms to MATLAB prototype.

    // END: projection0 }

    if ((n_x_min < 0.0) || (n_y_min < 0.0) || ((real_type) n_x < n_x_max) || ((real_type) n_y < n_y_max))
        ABORT("Support extends beyond projection.\n"); 

    if ((n0_x_min < 0.0) || (n0_y_min < 0.0) || ((real_type) n0_x < n0_x_max) || ((real_type) n0_y < n0_y_max))
        ABORT("Support extends beyond projection0.\n"); 
    
    // BEGIN: Unpack parameters and data. {

    // BEGIN: Unpack remap 2D parameters and data.

    real_type step_remap_2D = 1.0 / local_mag;

    // BEGIN: power_order {

    real_type *power_order = NULL;

    int n_power_order_bytes = power_order_n_x1 * power_order_n_x2 * sizeof(real_type);
    cudaErrorCode = cudaMalloc((void **) &power_order, n_power_order_bytes);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(power_order, %i bytes) failed.\n", n_power_order_bytes);

    cudaErrorCode = cudaMemcpy(power_order, power_order_h, n_power_order_bytes, cudaMemcpyHostToDevice);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy(power_order, power_order_h, %i bytes) failed.\n", n_power_order_bytes);

    // END: power_order }

    // BEGIN: map_2D {

    real_type *map_2D = NULL;

    // NOTE: Projection index is x3.
    int n_map_2D_bytes = map_2D_n_x1 * map_2D_n_x2 * sizeof(real_type);
    cudaErrorCode = cudaMalloc((void **) &map_2D, n_map_2D_bytes);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(map_2D, %i bytes) failed.\n", n_map_2D_bytes);

    cudaErrorCode = cudaMemcpy(map_2D, map_2D_current_h, n_map_2D_bytes, cudaMemcpyHostToDevice);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy(map_2D, map_2D_current_h, %i bytes) failed.\n", n_power_order_bytes);

    // END: map_2D }

    // END: Unpack remap 2D parameters and data.

    // BEGIN: Unpack rotation 2D parameters and data.

    real_type step_rotate_2D;
    if (test_rotation_matrix_2x2_for_subsampling(rot_matrix_2x2))
        step_rotate_2D = 0.25;
    else
        step_rotate_2D = 1.0;

    // END: Unpack rotation 2D parameters and data.

    // END: Unpack parameters and data. }

    real_type step = fmin(step_remap_2D, step_rotate_2D);
    //step = 1.0;  WARN("step is hardcoded to %.15e.\n", step);

    n0_x_min = step; // NOTE: Conforms to MATLAB prototype.
    n0_y_min = step; // NOTE: Conforms to MATLAB prototype.

    real_type n_subsamples_y_real = 1.0 / step;
    int n_subsamples_y = (int) floor(n_subsamples_y_real);
    if ((n_subsamples_y_real - ((real_type) n_subsamples_y)) != 0.0)
        ABORT("step = %.15e does not yield integral n_subsamples_y_real = %.15e.\n", step, n_subsamples_y_real);

    // NOTE: Number of subsamples in both directions is the same.
    int n_subsamples_x = (int) floor(n_subsamples_y_real);

    Print("step = %.15e\n", step);
    Print("n_subsamples_(xy) = (%d, %d)\n", n_subsamples_x, n_subsamples_y);

    //real_type sample_factor = step; // MATLAB prototype uses step for number of samples per pixel.
    real_type sample_factor = step * step;
    //real_type sample_factor = 1.0; WARN("sample_factor is hardcoded to %.15e.\n", sample_factor);

/*
    // BEGIN: rotation_2D {
 
    // NOTE: Use the inverse rotation matrix.
    real_type rot_matrix00 = rot_matrix_2x2[0][0];
    real_type rot_matrix01 = rot_matrix_2x2[1][0];
    real_type rot_matrix10 = rot_matrix_2x2[0][1];
    real_type rot_matrix11 = rot_matrix_2x2[1][1];

    real_type x0_init = 0.0;
    real_type y0_init = 0.0;

    real_type x_trans_x_init = n_x_min - center_x; // NOTE: From center of support in image coordinates, so from center of image.
    real_type x_trans_x_init_rot_x = rot_matrix00 * x_trans_x_init;
    real_type x_trans_x_init_rot_y = rot_matrix10 * x_trans_x_init;
    x0_init += x_trans_x_init_rot_x;
    y0_init += x_trans_x_init_rot_y;

    real_type y_trans_y_init = n_y_min - center_y; // NOTE: From center of support in image coordinates, so from center of image.
    real_type y_trans_y_init_rot_x = rot_matrix01 * y_trans_y_init;
    real_type y_trans_y_init_rot_y = rot_matrix11 * y_trans_y_init;
    x0_init += y_trans_y_init_rot_x;
    y0_init += y_trans_y_init_rot_y;

    x0_init += center0_x;
    y0_init += center0_y;

    real_type x0_del_x = rot_matrix00 * step;
    real_type y0_del_x = rot_matrix10 * step;

    real_type x0_del_y = rot_matrix01 * step;
    real_type y0_del_y = rot_matrix11 * step;

    // END: rotation_2D }
*/

    // BEGIN: Set constant GPU memory. {
    cudaMemcpyToSymbol("projection0_data", &projection0_data, sizeof(projection0_data)); // NOTE: cudaMalloc() and cudaMallocHost().
    cudaMemcpyToSymbol("n0_x", &n0_x, sizeof(n0_x));
    cudaMemcpyToSymbol("n0_y", &n0_y, sizeof(n0_y));
    cudaMemcpyToSymbol("n0_x_ws", &n0_x_ws, sizeof(n0_x_ws));
    cudaMemcpyToSymbol("center0_x", &center0_x, sizeof(center0_x));
    cudaMemcpyToSymbol("center0_y", &center0_y, sizeof(center0_y));
    cudaMemcpyToSymbol("n0_x_min", &n0_x_min, sizeof(n0_x_min));
    cudaMemcpyToSymbol("n0_x_max", &n0_x_max, sizeof(n0_x_max));
    cudaMemcpyToSymbol("n0_y_min", &n0_y_min, sizeof(n0_y_min));
    cudaMemcpyToSymbol("n0_y_max", &n0_y_max, sizeof(n0_y_max));
    cudaMemcpyToSymbol("projection_data", &projection_data, sizeof(projection_data));
    cudaMemcpyToSymbol("n_x", &n_x, sizeof(n_x));
    cudaMemcpyToSymbol("n_y", &n_y, sizeof(n_y));
    cudaMemcpyToSymbol("n_x_ws", &n_x_ws, sizeof(n_x_ws));
    cudaMemcpyToSymbol("center_x", &center_x, sizeof(center_x));
    cudaMemcpyToSymbol("center_y", &center_y, sizeof(center_y));
    cudaMemcpyToSymbol("n_x_min", &n_x_min, sizeof(n_x_min));
    cudaMemcpyToSymbol("n_x_max", &n_x_max, sizeof(n_x_max));
    cudaMemcpyToSymbol("n_y_min", &n_y_min, sizeof(n_y_min));
    cudaMemcpyToSymbol("n_y_max", &n_y_max, sizeof(n_y_max));
    cudaMemcpyToSymbol("step", &step, sizeof(step));
    cudaMemcpyToSymbol("n_subsamples_x", &n_subsamples_x, sizeof(n_subsamples_x));
    cudaMemcpyToSymbol("n_subsamples_y", &n_subsamples_y, sizeof(n_subsamples_y));
    cudaMemcpyToSymbol("sample_factor", &sample_factor, sizeof(sample_factor));
/*
    cudaMemcpyToSymbol("rot_matrix00", &rot_matrix00, sizeof(rot_matrix00));
    cudaMemcpyToSymbol("rot_matrix01", &rot_matrix01, sizeof(rot_matrix01));
    cudaMemcpyToSymbol("rot_matrix10", &rot_matrix10, sizeof(rot_matrix10));
    cudaMemcpyToSymbol("rot_matrix11", &rot_matrix11, sizeof(rot_matrix11));
    cudaMemcpyToSymbol("x0_init", &x0_init, sizeof(x0_init));
    cudaMemcpyToSymbol("y0_init", &y0_init, sizeof(y0_init));
    cudaMemcpyToSymbol("x0_del_x", &x0_del_x, sizeof(x0_del_x));
    cudaMemcpyToSymbol("y0_del_x", &y0_del_x, sizeof(y0_del_x));
    cudaMemcpyToSymbol("x0_del_y", &x0_del_y, sizeof(x0_del_y));
    cudaMemcpyToSymbol("y0_del_y", &y0_del_y, sizeof(y0_del_y));
*/
    cudaMemcpyToSymbol("map_2D_order", &map_2D_order, sizeof(map_2D_order));
    cudaMemcpyToSymbol("n_coeffs", &n_coeffs, sizeof(n_coeffs));
    cudaMemcpyToSymbol("power_order_n_x1", &power_order_n_x1, sizeof(power_order_n_x1));
    cudaMemcpyToSymbol("power_order_n_x2", &power_order_n_x2, sizeof(power_order_n_x2));
    cudaMemcpyToSymbol("power_order", &power_order, sizeof(power_order));
    cudaMemcpyToSymbol("map_2D_n_x1", &map_2D_n_x1, sizeof(map_2D_n_x1));
    cudaMemcpyToSymbol("map_2D_n_x2", &map_2D_n_x2, sizeof(map_2D_n_x2));
    cudaMemcpyToSymbol("map_2D", &map_2D, sizeof(map_2D));

    cudaErrorCode = cudaGetLastError();
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpyToSymbol() failed when setting constant GPU memory.\n");
    // END: Set constant GPU memory. }

    // BEGIN: N-pass version. {

    { // BEGIN: Pull scope.
    
        // BEGIN: Make all of this part of future transform_2D_CUDA_params? {

/*
        // BEGIN: One pixel per thread.
        unsigned int dimBlock_x = general_CUDA_params->dimBlock_x;
        unsigned int dimBlock_y = general_CUDA_params->dimBlock_y;
        dim3 dimBlock(dimBlock_x, dimBlock_y);

        unsigned int n_pixels_side = 1;
        //unsigned int n_pixels_side = 3;
        //unsigned int n_pixels_side = 5;
        //unsigned int n_pixels_side = 16;

        unsigned int n_cols_pass = ((unsigned int) n_x) / n_pixels_side;
        unsigned int n_rows_pass = ((unsigned int) n_y) / n_pixels_side;
    
        unsigned int dimGrid_x = n_cols_pass / dimBlock_x;
        if (n_cols_pass % dimBlock_x != 0)
            ++dimGrid_x; 
        unsigned int dimGrid_y = n_rows_pass / dimBlock_y;
        if (n_rows_pass % dimBlock_y != 0)
            ++dimGrid_y; 
        // END: One pixel per thread.
*/

        // BEGIN: Multiple pixels per thread.
        unsigned int dimBlock_x = 1;
        unsigned int dimBlock_y = 32;
        dim3 dimBlock(dimBlock_x, dimBlock_y);

        // NOTE: Read-before-write hazard with pull because of interpolation.
        unsigned int n_passes_x = 1;
        unsigned int n_passes_y = 2;

        unsigned int n_cols_pass = ((unsigned int) n0_x) / n_passes_x;
        unsigned int n_rows_pass = ((unsigned int) n0_y) / n_passes_y;

        unsigned int n_pixels_thread_x = n_cols_pass;
        unsigned int n_pixels_thread_y = 1;

        unsigned int dimGrid_x = 1;
//        unsigned int dimGrid_x = n_cols_pass / dimBlock_x;
        unsigned int dimGrid_y = n_rows_pass / dimBlock_y;
        if (n_rows_pass % dimBlock_y != 0)
            ++dimGrid_y;
        // END: Multiple pixels per thread.

        // END: Make all of this part of future transform_2D_CUDA_params? }
    
        dim3 dimGrid(dimGrid_x, dimGrid_y);
    
        Print("dimBlock.(xyz) = (%u, %u, %u)\n", dimBlock.x, dimBlock.y, dimBlock.z);
        Print("dimGrid.(xyz) = (%u, %u, %u)\n", dimGrid.x, dimGrid.y, dimGrid.z);
        //ABORT("ABORT!\n");
    
        // BEGIN: Shared memory array. {
    
        unsigned int n_threads_block = dimBlock_x * dimBlock_y;
    
        size_t n_shared_memory_bytes = 0;
    
/*
        int ivp_n_elems = map_2D_order + 1; // NOTE: Local memory.
        //int ivp_n_elems = (map_2D_order + 1) * ((int) n_threads_block); // NOTE: Shared memory.
        size_t ivp_n_bytes = ivp_n_elems * sizeof(real_type);
        n_shared_memory_bytes += ivp_n_bytes * 2;
    
        int diff_2D_n_x1 = map_2D_order + 1;
        int diff_2D_n_x2 = map_2D_order + 1;
        int diff_2D_n_elems = diff_2D_n_x1 * diff_2D_n_x2; // NOTE: Local memory.
        //int diff_2D_n_elems = diff_2D_n_x1 * diff_2D_n_x2 * ((int) n_threads_block); // NOTE: Shared memory. // WARNING: Not enough shared memory.
        size_t diff_2D_n_bytes = diff_2D_n_elems * sizeof(real_type);
        n_shared_memory_bytes += diff_2D_n_bytes * n_threads_block;
*/
    
        Print("n_threads_block = %u\n", n_threads_block);
        Print("n_shared_memory_bytes = %u\n", n_shared_memory_bytes);
        //ABORT("ABORT!\n");
    
        // END: Shared memory array. }
    
        for (int i_pass_x = 0; i_pass_x < n_passes_x; ++i_pass_x)
        {
            for (int i_pass_y = 0; i_pass_y < n_passes_y; ++i_pass_y)
            {
                Print("i_pass_(xy) = (%i, %i)\n", i_pass_x, i_pass_y);

                // NOTE: One pixel per thread: projection_inv_transform_2D_pull_CUDA_d<<<dimGrid, dimBlock, n_shared_memory_bytes>>>(n_pixels_side, n_pixels_side, i_pass_x, i_pass_y); // NOTE: One pixel per thread.
//                projection_inv_transform_2D_pull_CUDA_d<<<dimGrid, dimBlock, n_shared_memory_bytes>>>(n_cols_pass, 1, n_passes_x, i_pass_x, n_passes_y, i_pass_y); // NOTE: Multiple pixels per thread.
                projection_inv_transform_2D_pull_CUDA_d<<<dimGrid, dimBlock, n_shared_memory_bytes>>>(n_pixels_thread_x, n_pixels_thread_y, n_passes_x, i_pass_x, n_passes_y, i_pass_y); // NOTE: Multiple passes, multiple pixels per thread.
        
                cudaErrorCode = cudaGetLastError();
                if (cudaErrorCode != cudaSuccess)
                    ABORT_CUDA(cudaErrorCode, "Call to projection_inv_transform_2D_pull_CUDA_d<<<dimGrid, dimBlock, Ns>>>(%i, %i, %i, %i) failed.\n", n_passes_x, n_passes_y, i_pass_x, i_pass_y);
        
                cudaErrorCode = cudaThreadSynchronize();
                if (cudaErrorCode != cudaSuccess)
                    ABORT_CUDA(cudaErrorCode, "Call to cudaThreadSynchronize() failed.\n");
            }
        }
    
    } // END: Pull scope.

    // END: N-pass version. }

    // BEGIN: cudaMalloc() and CUDA array bound to texture.
    cudaErrorCode = cudaFree(projection_padded_data);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaFree(projection_padded_data) failed.\n");
    // END: cudaMalloc() and CUDA array bound to texture.

    // BEGIN: cudaMalloc() and cudaMallocHost().
    // Copy padded projection0 from device to host.
    cudaErrorCode = cudaMemcpy(projection0_padded_data_h, projection0_padded_data, projection0_padded->imageSize, cudaMemcpyDeviceToHost);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy(projection0_data_h, projection0_data, %i bytes) failed.\n", projection0_padded->imageSize);

    cudaErrorCode = cudaFree(projection0_padded_data);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaFree(projection0_padded_data) failed.\n");
    // END: cudaMalloc() and cudaMallocHost().

/*
    // BEGIN: CUDA array bound to texture. 
    cudaErrorCode = cudaFreeArray(projection0_array);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaFreeArray(projection0_array) failed.\n");
    // END: CUDA array bound to texture.
*/

//    exit(0);
}

// NOTE: This pulls pixel values from the original projection to the new one.
// NOTE: Both dimensions use finite difference scheme.
__host__ void projection_rotate_2D_prototype_CUDA
(
    general_CUDA_params_handle general_CUDA_params,
    IplImage_base *projection0_padded,
    IplImage_base *projection0,
    real_type support0[4],
    real_type rot_matrix_2x2[2][2],
    IplImage_base *projection_padded,
    IplImage_base *projection,
    real_type support[4]
)
{
// NOTE: Does NOT use PROTOTYPE_COMPLIANT_INDEXING.
//#if PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

//    clock_t t_0, t_1;
//    double secs_per_clock = 1.0 / CLOCKS_PER_SEC;

    if ((projection0->depth != IPL_DEPTH_32F_BASE) && (projection0->depth != IPL_DEPTH_64F_BASE))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F_BASE) or %i (IPL_DEPTH_64F_BASE), but image->depth = %i.\n",
            IPL_DEPTH_32F_BASE, IPL_DEPTH_64F_BASE, projection0->depth);

    if ((projection->depth != IPL_DEPTH_32F_BASE) && (projection->depth != IPL_DEPTH_64F_BASE))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F_BASE) or %i (IPL_DEPTH_64F_BASE), but image->depth = %i.\n",
            IPL_DEPTH_32F_BASE, IPL_DEPTH_64F_BASE, projection->depth);

    cvSetZero_base(projection);

    cudaError_t cudaErrorCode;

    // BEGIN: projection0 {

    // BEGIN: Move projection0_padded to device, set ROI.

    // BEGIN: cudaMalloc() and cudaMallocHost(). {

    pixel_type *projection0_padded_data_h = (pixel_type *) projection0_padded->imageData;

    pixel_type *projection0_padded_data = NULL;

    cudaErrorCode = cudaMalloc((void **) &projection0_padded_data, projection0_padded->imageSize);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(projection0_padded_data, %i bytes) failed.\n", projection0_padded->imageSize);
    WARN("Call to cudaMalloc(projection0_padded_data, %i bytes).\n", projection0_padded->imageSize);

    cudaErrorCode = cudaMemcpy(projection0_padded_data, projection0_padded_data_h, projection0_padded->imageSize, cudaMemcpyHostToDevice);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy(projection0_padded_data, projection0_padded_data_h, %i bytes) failed.\n", projection0_padded->imageSize);

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type).
    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / sizeof(pixel_type);
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.\n");

    uint_address_type projection0_elem_offset = calc_elem_offset_of_IplImage_base_ROI_widthStep(projection0_padded, projection0);
    pixel_type *projection0_data = projection0_padded_data + projection0_elem_offset;

    Print("projection0_padded_data = %llu\n", projection0_padded_data);
    Print("projection0_elem_offset = %llu\n", projection0_elem_offset);
    Print("projection0_data = %llu\n", projection0_data);
    Print("n0_x_ws = %i\n", n0_x_ws);

    // BEGIN: DEBUG.
    uint_address_type projection0_padded_data_min_address = (uint_address_type) projection0_padded_data;
    uint_address_type projection0_padded_data_max_address = projection0_padded_data_min_address + projection0_padded->imageSize - sizeof(pixel_type);
    // END: DEBUG.

    // END: Move projection0_padded to device, set ROI.

    // END: cudaMalloc() and cudaMallocHost(). }

/*
    // BEGIN: cudaMallocPitch().

    // WARNING: Using data allocated by cudaMallocPitch() will require writing array indexing macros that use a widthStep in bytes and _not_ number of elements per row as is currently implemented.
    size_t projection0_padded_data_pitch;
    cudaErrorCode = cudaMallocPitch((void **) &projection0_padded_data, projection0_padded_pitch, projection0_padded->widthStep, projection0_padded->height);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMallocPitch(projection0_padded_data, (%i widthStep bytes) * (%i rows)) failed.\n", projection0_padded->widthStep, projection0_padded->height);

    // END: cudaMallocPitch().
*/

/*
    // BEGIN: CUDA array bound to texture. {

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type) and that projection0->widthStep / sizeof(pixel_type) == projection0->width.
    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / sizeof(pixel_type);
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.\n");
    //if (n0_x_ws != projection0->width)
    //    ABORT("(n0_x_ws == %i) != (projection0->width == %i)\n", n0_x_ws, projection0->width);

    // Allocate array and copy image data.
    cudaChannelFormatDesc projection0_channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

    cudaArray *projection0_array = NULL;

    cudaErrorCode = cudaMallocArray(&projection0_array, &projection0_channelDesc, projection0->width, projection0->height);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMallocArray(&projection0_array, &projection0_channelDesc, %i, %i) failed.\n", projection0->width, projection0->height);

    //cudaErrorCode = cudaMemcpyToArray(projection0_array, 0, 0, (void *) projection0->imageData, projection0->imageSize, cudaMemcpyHostToDevice);
    //if (cudaErrorCode != cudaSuccess)
    //   ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpyToArray(projection0_array, 0, 0, projection0_data, %i, %i, cudaMemcpyHostToDevice) failed.\n", projection0->width, projection0->height);

    cudaErrorCode = cudaMemcpy2DToArray(projection0_array, 0, 0, projection0->imageData, projection0->widthStep, projection0->width * sizeof(pixel_type), projection0->height, cudaMemcpyHostToDevice);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy2DToArray(projection0_array, 0, 0, projection0->imageData, %i, %i, %i, cudaMemcpyHostToDevice) failed.\n", projection0->widthStep, projection0->width * sizeof(pixel_type), projection0->height);

    // Set texture parameters.
    projection0_texture.addressMode[0] = cudaAddressModeWrap;
    projection0_texture.addressMode[1] = cudaAddressModeWrap;
//    projection0_texture.addressMode[0] = cudaAddressModeClamp;
//    projection0_texture.addressMode[1] = cudaAddressModeClamp;
    projection0_texture.filterMode = cudaFilterModeLinear;
    projection0_texture.normalized = true;

    // Bind the array to the texture.
    cudaBindTextureToArray(projection0_texture, projection0_array, projection0_channelDesc);

/\*
    // DEBUG: BEGIN: Copy CUDA array from device for later write to image.
    memset(projection0->imageData, 0, projection0->imageSize);

    //cudaErrorCode = cudaMemcpyFromArray((void *) projection0->imageData, projection0_array, 0, 0, projection0->imageSize, cudaMemcpyDeviceToHost);
    //if (cudaErrorCode != cudaSuccess)
    //   ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpyFromArray((void *) projection0->imageData, projection0_array, %i, cudaMemcpyDeviceToHost) failed.\n", projection0->imageSize);

    cudaErrorCode = cudaMemcpy2DFromArray((void *) projection0->imageData, projection0->widthStep, projection0_array, 0, 0, projection0->width * sizeof(pixel_type), projection0->height, cudaMemcpyDeviceToHost);
    if (cudaErrorCode != cudaSuccess)
       ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy2DFromArray((void *) projection0->imageData, %i, projection0_array, %i, %i, cudaMemcpyDeviceToHost) failed.\n", projection0->widthStep, projection0->width * sizeof(pixel_type), projection0->height);

    return;
    // DEBUG: END: Copy CUDA array from device for later write to image.
*\/

    // END: CUDA array bound to texture. }
*/

/*
    // BEGIN: cudaMallocHost() (WARNING: Asynchronous reads/writes are strictly host functions.  Dumb, dumb, dumb!). {

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type).
    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / sizeof(pixel_type);
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.\n");

    // Set ROI.
    pixel_type *projection0_data = (pixel_type *) projection0->imageData;

    // END: cudaMallocHost() (WARNING: Asynchronous reads/writes are strictly host functions.  Dumb, dumb, dumb!). }
*/

    int n0_x = projection0->width;
    int n0_y = projection0->height;

    real_type center0_x = support0[0];
    real_type center0_y = support0[1];

    // For support vs. projection check.
    real_type n0_x_min = 0.0;
    real_type n0_x_max = support0[2];
    real_type n0_y_min = 0.0;
    real_type n0_y_max = support0[3];

    // END: projection0 }

    // BEGIN: projection {

    // BEGIN: cudaMalloc(). {

    // BEGIN: Initialize projection on device, set ROI.

    pixel_type *projection_padded_data_h = (pixel_type *) projection_padded->imageData;

    pixel_type *projection_padded_data = NULL;
    cudaErrorCode = cudaMalloc((void **) &projection_padded_data, projection_padded->imageSize);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMalloc(projection_padded_data, %i bytes) failed.\n", projection_padded->imageSize);
    WARN("Call to cudaMalloc(projection_padded_data, %i bytes).\n", projection_padded->imageSize);

    cudaErrorCode = cudaMemset(projection_padded_data, 0, projection_padded->imageSize);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMemset(projection_padded_data, 0, %i bytes) failed.\n", projection_padded->imageSize);

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type).
    real_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / sizeof(pixel_type);
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.\n");

    uint_address_type projection_elem_offset = calc_elem_offset_of_IplImage_base_ROI_widthStep(projection_padded, projection);
    pixel_type *projection_data = projection_padded_data + projection_elem_offset;

    Print("projection_padded_data = %llu\n", projection_padded_data);
    Print("projection_elem_offset = %llu\n", projection_elem_offset);
    Print("projection_data = %llu\n", projection_data);
    Print("n_x_ws = %i\n", n_x_ws);

    // BEGIN: DEBUG.
    uint_address_type projection_padded_data_min_address = (uint_address_type) projection_padded_data;
    uint_address_type projection_padded_data_max_address = (projection_padded_data_min_address + projection_padded->imageSize - sizeof(pixel_type));
    // END: DEBUG.

    // END: Initialize projection on device, set ROI.

    // END: cudaMalloc(). }

/*
    // BEGIN: cudaMallocHost(). {

    // NOTE: Assumes that (and checks if) widthStep of row is evenly divisible by sizeof(pixel_type).
    pixel_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / sizeof(pixel_type);
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.\n");

    // Set ROI.
    pixel_type *projection_data = (pixel_type *) projection->imageData;

    // END: cudaMallocHost(). }
*/

    int n_x = projection->width;
    int n_y = projection->height;

    real_type center_x = support[0];
    real_type center_y = support[1];

    real_type n_x_min = 0.0;
    real_type n_x_max = support[2];
    real_type n_y_min = 0.0;
    real_type n_y_max = support[3];

    // END: projection }

    if ((n0_x_min < 0.0) || (n0_y_min < 0.0) || ((real_type) n0_x < n0_x_max) || ((real_type) n0_y < n0_y_max))
        ABORT("Support extends beyond projection0.\n"); 
    
    if ((n_x_min < 0.0) || (n_y_min < 0.0) || ((real_type) n_x < n_x_max) || ((real_type) n_y < n_y_max))
        ABORT("Support extends beyond projection.\n"); 

    real_type step;
    if (test_rotation_matrix_2x2_for_subsampling(rot_matrix_2x2))
        step = 0.25;
    else
        step = 1.0;

    //step = 1.0; WARN("step is hardcoded to %.15e.\n", step);

    real_type n_subsamples_y_real = 1.0 / step;
    int n_subsamples_y = (int) floor(n_subsamples_y_real);
    if ((n_subsamples_y_real - ((real_type) n_subsamples_y)) != 0.0)
        ABORT("step = %.15e does not yield integral n_subsamples_y_real = %.15e.\n", step, n_subsamples_y_real);

    // NOTE: Number of subsamples in both directions is the same.
    int n_subsamples_x = (int) floor(n_subsamples_y_real);

    Print("step = %.15e\n", step);
    Print("n_subsamples_(xy) = (%d, %d)\n", n_subsamples_x, n_subsamples_y);

    //real_type sample_factor = step; // MATLAB prototype uses step for number of samples per pixel.
    real_type sample_factor = step * step;
    //real_type sample_factor = 1.0; WARN("sample_factor is hardcoded to %.15e.\n", sample_factor);

    // NOTE: Use the inverse rotation matrix.
    real_type rot_matrix00 = rot_matrix_2x2[0][0];
    real_type rot_matrix01 = rot_matrix_2x2[1][0];
    real_type rot_matrix10 = rot_matrix_2x2[0][1];
    real_type rot_matrix11 = rot_matrix_2x2[1][1];

    real_type x0_init = 0.0;
    real_type y0_init = 0.0;

    real_type x_trans_x_init = n_x_min - center_x; // NOTE: From center of support in image coordinates, so from center of image.
    real_type x_trans_x_init_rot_x = rot_matrix00 * x_trans_x_init;
    real_type x_trans_x_init_rot_y = rot_matrix10 * x_trans_x_init;
    x0_init += x_trans_x_init_rot_x;
    y0_init += x_trans_x_init_rot_y;

    real_type y_trans_y_init = n_y_min - center_y; // NOTE: From center of support in image coordinates, so from center of image.
    real_type y_trans_y_init_rot_x = rot_matrix01 * y_trans_y_init;
    real_type y_trans_y_init_rot_y = rot_matrix11 * y_trans_y_init;
    x0_init += y_trans_y_init_rot_x;
    y0_init += y_trans_y_init_rot_y;

    x0_init += center0_x;
    y0_init += center0_y;

    real_type x0_del_x = rot_matrix00 * step;
    real_type y0_del_x = rot_matrix10 * step;

    real_type x0_del_y = rot_matrix01 * step;
    real_type y0_del_y = rot_matrix11 * step;

    // BEGIN: Set constant GPU memory. {
    cudaMemcpyToSymbol("projection0_padded_data_min_address", &projection0_padded_data_min_address, sizeof(projection0_padded_data_min_address)); // DEBUG. // NOTE: cudaMalloc() and cudaMallocHost().
    cudaMemcpyToSymbol("projection0_padded_data_max_address", &projection0_padded_data_max_address, sizeof(projection0_padded_data_max_address)); // DEBUG. // NOTE: cudaMalloc() and cudaMallocHost().
    cudaMemcpyToSymbol("projection0_data", &projection0_data, sizeof(projection0_data)); // NOTE: cudaMalloc() and cudaMallocHost().
    cudaMemcpyToSymbol("n0_x", &n0_x, sizeof(n0_x));
    cudaMemcpyToSymbol("n0_y", &n0_y, sizeof(n0_y));
    cudaMemcpyToSymbol("n0_x_ws", &n0_x_ws, sizeof(n0_x_ws));
    cudaMemcpyToSymbol("center0_x", &center0_x, sizeof(center0_x));
    cudaMemcpyToSymbol("center0_y", &center0_y, sizeof(center0_y));
    cudaMemcpyToSymbol("n0_x_min", &n0_x_min, sizeof(n0_x_min));
    cudaMemcpyToSymbol("n0_x_max", &n0_x_max, sizeof(n0_x_max));
    cudaMemcpyToSymbol("n0_y_min", &n0_y_min, sizeof(n0_y_min));
    cudaMemcpyToSymbol("n0_y_max", &n0_y_max, sizeof(n0_y_max));
    cudaMemcpyToSymbol("projection_padded_data_min_address", &projection_padded_data_min_address, sizeof(projection_padded_data_min_address)); // DEBUG. // NOTE: cudaMalloc() _and_ cudaMallocHost()?
    cudaMemcpyToSymbol("projection_padded_data_max_address", &projection_padded_data_max_address, sizeof(projection_padded_data_max_address)); // DEBUG. // NOTE: cudaMalloc() _and_ cudaMallocHost()?
    cudaMemcpyToSymbol("projection_data", &projection_data, sizeof(projection_data)); // NOTE: cudaMalloc() _and_ cudaMallocHost()?
    cudaMemcpyToSymbol("n_x", &n_x, sizeof(n_x));
    cudaMemcpyToSymbol("n_y", &n_y, sizeof(n_y));
    cudaMemcpyToSymbol("n_x_ws", &n_x_ws, sizeof(n_x_ws));
    cudaMemcpyToSymbol("center_x", &center_x, sizeof(center_x));
    cudaMemcpyToSymbol("center_y", &center_y, sizeof(center_y));
    cudaMemcpyToSymbol("n_x_min", &n_x_min, sizeof(n_x_min));
    cudaMemcpyToSymbol("n_x_max", &n_x_max, sizeof(n_x_max));
    cudaMemcpyToSymbol("n_y_min", &n_y_min, sizeof(n_y_min));
    cudaMemcpyToSymbol("n_y_max", &n_y_max, sizeof(n_y_max));
    cudaMemcpyToSymbol("step", &step, sizeof(step));
    cudaMemcpyToSymbol("n_subsamples_x", &n_subsamples_x, sizeof(n_subsamples_x));
    cudaMemcpyToSymbol("n_subsamples_y", &n_subsamples_y, sizeof(n_subsamples_y));
    cudaMemcpyToSymbol("sample_factor", &sample_factor, sizeof(sample_factor));
    cudaMemcpyToSymbol("rot_matrix00", &rot_matrix00, sizeof(rot_matrix00));
    cudaMemcpyToSymbol("rot_matrix01", &rot_matrix01, sizeof(rot_matrix01));
    cudaMemcpyToSymbol("rot_matrix10", &rot_matrix10, sizeof(rot_matrix10));
    cudaMemcpyToSymbol("rot_matrix11", &rot_matrix11, sizeof(rot_matrix11));
    cudaMemcpyToSymbol("x0_init", &x0_init, sizeof(x0_init));
    cudaMemcpyToSymbol("y0_init", &y0_init, sizeof(y0_init));
    cudaMemcpyToSymbol("x0_del_x", &x0_del_x, sizeof(x0_del_x));
    cudaMemcpyToSymbol("y0_del_x", &y0_del_x, sizeof(y0_del_x));
    cudaMemcpyToSymbol("x0_del_y", &x0_del_y, sizeof(x0_del_y));
    cudaMemcpyToSymbol("y0_del_y", &y0_del_y, sizeof(y0_del_y));

    cudaErrorCode = cudaGetLastError();
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpyToSymbol() failed when setting constant GPU memory.\n");
    // END: Set constant GPU memory. }

    // BEGIN: One-pass version. {

    // BEGIN: Make all of this part of future rotate_2D_CUDA_params?
    unsigned int dimBlock_x = general_CUDA_params->dimBlock_x;
    unsigned int dimBlock_y = general_CUDA_params->dimBlock_y;
    dim3 dimBlock(dimBlock_x, dimBlock_y);
    unsigned int n_cols = (unsigned int) n_x;
    unsigned int n_rows = (unsigned int) n_y;
    unsigned int dimGrid_x = n_cols / dimBlock_x;
    if (n_cols % dimBlock_x != 0)
        ++dimGrid_x; 
    unsigned int dimGrid_y = n_rows / dimBlock_y;
    if (n_rows % dimBlock_y != 0)
        ++dimGrid_y; 
    dim3 dimGrid(dimGrid_x, dimGrid_y, 1u);
    // END: Make all of this part of future rotate_2D_CUDA_params?

    Print("dimBlock.(xyz) = (%u, %u, %u)\n", dimBlock.x, dimBlock.y, dimBlock.z);
    Print("dimGrid.(xyz) = (%u, %u, %u)\n", dimGrid.x, dimGrid.y, dimGrid.z);
    //ABORT("ABORT!\n");

    projection_rotate_2D_prototype_CUDA_d<<<dimGrid, dimBlock>>>();

    cudaErrorCode = cudaGetLastError();
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to projection_rotate_2D_prototype_CUDA_d<<<dimGrid, dimBlock>>>() failed.\n");
   
    cudaErrorCode = cudaThreadSynchronize();
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaThreadSynchronize() failed.\n");

    // END: One-pass version. }

/*
    // BEGIN: Two-pass version. {

    unsigned int dimBlock_y = 8u;
    dim3 dimBlock(1u, dimBlock_y);
    unsigned int n_rows_half = (unsigned int) (n_y / 2);
    unsigned int dimGrid_y = n_rows_half / dimBlock_y;
//    if (n_y % dimBlock_y != 0)
//        ++dimGrid_y; 
    dim3 dimGrid(1u, dimGrid_y);

    Print("dimBlock.(xyz) = (%u, %u, %u)\n", dimBlock.x, dimBlock.y, dimBlock.z);
    Print("dimGrid.(xyz) = (%u, %u, %u)\n", dimGrid.x, dimGrid.y, dimGrid.z);
    //ABORT("ABORT!\n");

    for (int i_offset_row = 0; i_offset_row < 2; ++i_offset_row)
    {

        projection_rotate_2D_prototype_CUDA_d<<<dimGrid, dimBlock>>>(i_offset_row);
    
        cudaErrorCode = cudaGetLastError();
        if (cudaErrorCode != cudaSuccess)
            ABORT_CUDA(cudaErrorCode, "Call to projection_rotate_2D_prototype_CUDA_d<<<dimGrid, dimBlock>>>() failed.\n");
    
        cudaErrorCode = cudaThreadSynchronize();
        if (cudaErrorCode != cudaSuccess)
            ABORT_CUDA(cudaErrorCode, "Call to cudaThreadSynchronize() failed.\n");

    }

    // END: Two-pass version. }
*/

    // BEGIN: cudaMalloc() and CUDA array bound to texture.
    // Copy padded projection from device to host.
    cudaErrorCode = cudaMemcpy(projection_padded_data_h, projection_padded_data, projection_padded->imageSize, cudaMemcpyDeviceToHost);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaMemcpy(projection_data, projection_data_d, %i bytes) failed.\n", projection_padded->imageSize);

    cudaErrorCode = cudaFree(projection_padded_data);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaFree(projection_data_d) failed.\n");
    // END: cudaMalloc() and CUDA array bound to texture.

    // BEGIN: cudaMalloc() and cudaMallocHost().
    cudaErrorCode = cudaFree(projection0_padded_data);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaFree(projection0_data_d) failed.\n");
    // END: cudaMalloc() and cudaMallocHost().

/*
    // BEGIN: CUDA array bound to texture. 
    cudaErrorCode = cudaFreeArray(projection0_array);
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaFreeArray(projection0_array) failed.\n");
    // END: CUDA array bound to texture.
*/

//    exit(0);
}

////////////////////////////////////////////////////////////////////////////////
// END: Host code. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Hardware initialization. {
////////////////////////////////////////////////////////////////////////////////

__host__ void init_GPU_CUDA(int deviceID)
{
    int device_count = 0;
    CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceCount(&device_count));

    Print("device_count = %d\n", device_count);

    if ((deviceID < 0) || (deviceID > device_count))
        ABORT("deviceID == %d.\n");

    Print("Setting device, deviceID = %i.\n", deviceID);
    CUDA_SAFE_CALL(cudaSetDevice(deviceID));

    Print("Initializing device, deviceID = %i.\n", deviceID);
    print_device_properties_CUDA(deviceID);
    check_device_properties_CUDA(deviceID);

    //print_and_check_all_devices_properties_CUDA();
    //ABORT("ABORT!\n");
}

__host__ void clear_GPU_CUDA(int deviceID)
{
    int device_count = 0;
    CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceCount(&device_count));

    Print("device_count = %d\n", device_count);

    if ((deviceID < 0) || (deviceID > device_count))
        ABORT("deviceID == %d.\n");

    Print("Setting device, deviceID = %i.\n", deviceID);
    CUDA_SAFE_CALL(cudaSetDevice(deviceID));

    Print("Clearing device, deviceID = %i.\n", deviceID);
    cudaError_t cudaErrorCode = cudaThreadExit();
    if (cudaErrorCode != cudaSuccess)
        ABORT_CUDA(cudaErrorCode, "Call to cudaThreadExit() failed.\n");
}

__host__ void projection_series_transform_2D_CUDA_filter_1D_inv_transform_2D_CUDA_init_GPU(int deviceID)
{
    init_GPU_CUDA(deviceID);
}

__host__ void projection_series_transform_2D_CUDA_filter_1D_inv_transform_2D_CUDA_clear_GPU(int deviceID)
{
    clear_GPU_CUDA(deviceID);
}

__host__ void projection_series_rotate_2D_CUDA_filter_1D_rotate_2D_CUDA_init_GPU(int deviceID)
{
    init_GPU_CUDA(deviceID);
}

__host__ void projection_series_rotate_2D_CUDA_filter_1D_rotate_2D_CUDA_clear_GPU(int deviceID)
{
    clear_GPU_CUDA(deviceID);
}

////////////////////////////////////////////////////////////////////////////////
// END: Hardware initialization. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Memory management. {
////////////////////////////////////////////////////////////////////////////////

void malloc_page_locked_pixel_type_CUDA(pixel_type **data, size_t size)
{
    // gcc produces "warning: dereferencing type-punned pointer will break strict-aliasing rules".
    mallocHost_CUDA((void **) data, size);
}

void free_page_locked_pixel_type_CUDA(pixel_type *data)
{
    freeHost_CUDA((void *) data);
}

////////////////////////////////////////////////////////////////////////////////
// END: Memory management. }
////////////////////////////////////////////////////////////////////////////////

/*
#ifdef __cplusplus /\* If this is a C++ compiler, use C linkage. *\/
}
#endif
*/
