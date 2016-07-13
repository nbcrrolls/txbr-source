#ifndef BASE_UTILITIES_H
#define BASE_UTILITIES_H

typedef unsigned long long int uint_address_type;

//typedef int int_type;
typedef long long int int_type;

//typedef float real_type;
typedef double real_type;

typedef float pixel_type;
//typedef double pixel_type;
//typedef real_type pixel_type;

#define EPSILON 1e-5

#ifdef __cplusplus

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <cctype>
#include <cstring>
#include <ctime>
#include <cmath>
#include <complex>
typedef complex<pixel_type> complex_type; // Assumed to represent pixel values.

#else

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>

// BEGIN: These are assumed to represent pixel values.
typedef float complex complex_type; // WARNING: float vs. double
//typedef double complex complex_type; // WARNING: float vs. double
// END: These are assumed to represent pixel values.

#endif

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage. */
// NOTE: This is C99 code.
extern "C"
{
#endif

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Output utilities. {
////////////////////////////////////////////////////////////////////////////////

int Print(char* format, ...);

void Warn(char* format, ...);
void Warn_with_info(const char *file, int line, const char *func, char *format, ...);

#define WARN(format, ...) \
    Warn_with_info(__FILE__, __LINE__, __func__, format, ## __VA_ARGS__)

void Abort(char* format, ...);
void Abort_with_info(const char *file, int line, const char *func, char *format, ...);

#define ABORT(format, ...) \
    Abort_with_info(__FILE__, __LINE__, __func__, format, ## __VA_ARGS__)

void print_matrix_2x2(const char *name, real_type matrix[2][2]);
void print_matrix_1D(const char *name, real_type *matrix, int n_x1);
void print_matrix_2D(const char *name, real_type *matrix, int n_x1, int n_x2);
void print_matrix_3D(const char *name, real_type *matrix, int n_x1, int n_x2, int n_x3);

////////////////////////////////////////////////////////////////////////////////
// END: Output utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Mathematical utilities. {
////////////////////////////////////////////////////////////////////////////////

#define POW2(x_) ((x_) * (x_))
#define POW3(x_) ((x_) * (x_) * (x_))

////////////////////////////////////////////////////////////////////////////////
// END: Mathematical utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Image/matrix utilities. {
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Matrix/Matrix-like data structure access macros. {
////////////////////////////////////////////////////////////////////////////////

// NOTE: 
// 
// Unless stated otherwise, the following assume (x, y) == (0, 0) is lower
// lefthand corner with x horizontal and y vertical.  In other words, x is the
// column index and y is the row index.
//
// NOTE: Use (widthStep / sizeof(elem)) to index into possible IplImage ROI.  __Not__ guaranteed to work with data padded for alignment!
//
// "MATRIX_" prefix indicates math order indexing, i.e., i-th row, j-th column, k-th plane.
// "_CM" suffix indicates column-major indexing.
// "_PIX" suffix usually indicates a typecast up to real_type on values retrieved from and down to pixel_type on values stored to the array/matrix.  Otherwise, it simply indicates that the elements are of type pixel_type.

#define PUTVAL_2D(array_2D_, n_x_, x_, y_, v_) \
    ((array_2D_)[(x_) + ((y_) * (n_x_))] = (v_))

#define PUTVAL_2D_PIX(array_2D_, n_x_, x_, y_, v_) \
    ((array_2D_)[(x_) + ((y_) * (n_x_))] = ((pixel_type) (v_)))

/*
void ADDVAL_2D(real_type array_2D_[], int n_x_, int x_, int y_, real_type v_)
{
    (array_2D_)[(x_) + ((y_) * (n_x_))] += v_;
}
*/

#define ADDVAL_2D(array_2D_, n_x_, x_, y_, v_) \
    (((array_2D_)[(x_) + ((y_) * (n_x_))]) += (v_))

#define ADDVAL_2D_PIX(array_2D_, n_x_, x_, y_, v_) \
    (((array_2D_)[(x_) + ((y_) * (n_x_))]) += ((pixel_type) (v_)))

/*
real_type INDEX_2D(real_type array_2D_[], int n_x_, int x_, int y_)
{
    return (array_2D_)[(x_) + ((y_) * (n_x_))];
}
*/

#define INDEX_2D(array_2D_, n_x_, x_, y_) \
    ((array_2D_)[(x_) + ((y_) * (n_x_))])

#define INDEX_2D_PIX(array_2D_, n_x_, x_, y_) \
    ((real_type) ((array_2D_)[(x_) + ((y_) * (n_x_))]))

/*
uint_address_type ADDRESS_2D_PIX(pixel_type array_2D_[], int n_x_, int x_, int y_)
{
    return (((uint_address_type) array_2D_) + (((uint_address_type) ((x_) + ((y_) * (n_x_)))) * (sizeof(pixel_type))));
}
*/

#define ADDRESS_2D_PIX(array_2D_, n_x_, x_, y_) \
    (((uint_address_type) array_2D_) + (((uint_address_type) (((int) x_) + (((int) y_) * ((int) n_x_)))) * (sizeof(pixel_type))))

// column-major 2D array access.
#define INDEX_2D_CM(array_2D_, n_y_, x_, y_) \
    ((array_2D_)[((x_) * (n_y_)) + (y_)])

#define MATRIX_INDEX_2D(array_2D_, n_j_, i_, j_) \
    ((array_2D_)[((i_) * (n_j_)) + (j_)])

// column-major 2D matrix access.
#define MATRIX_INDEX_2D_CM(array_2D_, n_i_, i_, j_) \
    ((array_2D_)[(i_) + ((j_) * (n_i_))])

// WARNING: Row-major, but not equivalent to C-style multidimensional array indexing.
#define MATRIX_INDEX_3D(array_3D_, n_i_, n_j_, i_, j_, k_) \
    ((array_3D_)[((i_) * (n_j_)) + (j_) + ((k_) * (n_i_) * (n_j_))])

// column-major 3D matrix array access.
#define MATRIX_INDEX_3D_CM(array_3D_, n_i_, n_j_, i_, j_, k_) \
    ((array_3D_)[(i_) + ((j_) * (n_i_)) + ((k_) * (n_i_) * (n_j_))])

// No interpolation.
#define INTERPOLATE_NONE(array_2D_, n_x_, x_, y_, x_alpha_, y_alpha_) \
    (INDEX_2D((array_2D_), (n_x_), (x_), (y_)))

#define INTERPOLATE_NONE_PIX(array_2D_, n_x_, x_, y_, x_alpha_, y_alpha_) \
    (INDEX_2D_PIX((array_2D_), (n_x_), (x_), (y_)))

// Interpolation with origin in lower left-hand corner.
#define INTERPOLATE_2D_LL(array_2D_, n_x_, x_, y_, x_alpha_, y_alpha_) \
    (INDEX_2D((array_2D_), (n_x_), (x_), (y_)) * (1.0 - (x_alpha_)) * (1.0 - (y_alpha_)) + \
    INDEX_2D((array_2D_), (n_x_), (x_) + 1, (y_)) * (x_alpha_) * (1.0 - (y_alpha_)) + \
    INDEX_2D((array_2D_), (n_x_), (x_), (y_) + 1) * (1.0 - (x_alpha_)) * (y_alpha_) + \
    INDEX_2D((array_2D_), (n_x_), (x_) + 1, (y_) + 1) * (x_alpha_) * (y_alpha_))

#define INTERPOLATE_2D_LL_PIX(array_2D_, n_x_, x_, y_, x_alpha_, y_alpha_) \
    (INDEX_2D_PIX((array_2D_), (n_x_), (x_), (y_)) * (1.0 - (x_alpha_)) * (1.0 - (y_alpha_)) + \
    INDEX_2D_PIX((array_2D_), (n_x_), (x_) + 1, (y_)) * (x_alpha_) * (1.0 - (y_alpha_)) + \
    INDEX_2D_PIX((array_2D_), (n_x_), (x_), (y_) + 1) * (1.0 - (x_alpha_)) * (y_alpha_) + \
    INDEX_2D_PIX((array_2D_), (n_x_), (x_) + 1, (y_) + 1) * (x_alpha_) * (y_alpha_))

// Interpolation with origin in upper right-hand corner.
#define INTERPOLATE_2D_UR(array_2D_, n_x_, n_y_, x_, y_, x_alpha_, y_alpha_) \
    (INDEX_2D((array_2D_), (n_x_), (n_x_) - 1 - (x_), (n_y_) - 1 - (y_)) * (1.0 - (x_alpha_)) * (1.0 - (y_alpha_)) + \
    INDEX_2D((array_2D_), (n_x_), (n_x_) - (x_), (n_y_) - 1 - (y_)) * (x_alpha_) * (1.0 - (y_alpha_)) + \
    INDEX_2D((array_2D_), (n_x_), (n_x_) - (x_) - 1, (n_y_) - (y_)) * (1.0 - (x_alpha_)) * (y_alpha_) + \
    INDEX_2D((array_2D_), (n_x_), (n_x_) - (x_), (n_y_) - (y_)) * (x_alpha_) * (y_alpha_))

#define INTERPOLATE_2D_UR_PIX(array_2D_, n_x_, n_y_, x_, y_, x_alpha_, y_alpha_) \
    (INDEX_2D_PIX((array_2D_), (n_x_), (n_x_) - 1 - (x_), (n_y_) - 1 - (y_)) * (1.0 - (x_alpha_)) * (1.0 - (y_alpha_)) + \
    INDEX_2D_PIX((array_2D_), (n_x_), (n_x_) - (x_), (n_y_) - 1 - (y_)) * (x_alpha_) * (1.0 - (y_alpha_)) + \
    INDEX_2D_PIX((array_2D_), (n_x_), (n_x_) - (x_) - 1, (n_y_) - (y_)) * (1.0 - (x_alpha_)) * (y_alpha_) + \
    INDEX_2D_PIX((array_2D_), (n_x_), (n_x_) - (x_), (n_y_) - (y_)) * (x_alpha_) * (y_alpha_))

////////////////////////////////////////////////////////////////////////////////
// END: Matrix/Matrix-like data structure access macros. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Image utilities. {
////////////////////////////////////////////////////////////////////////////////

#define IPL_DEPTH_32F_BASE 32
#define IPL_DEPTH_64F_BASE 64

typedef struct 
{
//    int nSize; // Unused.
//    int ID; // Unused.
    int nChannels; // Usually 1, but may need more at some point.
//    int alphaChannel; // Unused.
    int depth; // Pixel depth in bits.
//    char colorModel[4]; // Unused.
//    char channelSeq[4]; // Unused.
//    int dataOrder; // Unused.  Multiple channels are interleaved.
//    int origin;  // Unused.  Always lower-lefthand corner.
//    int align; // Unused.
    int width; // Number of columns.
    int height; // Number of rows.
//    struct _IplROI *roi; // Unused.
//    struct _IplImage * maskROI; // Unused.
//    void* imageId; // Unused.
//    struct _IplTileInfo* tileInfo; // Unused.
    int imageSize; // imageData size in bytes == image->height * image->widthStep
    char *imageData;
    int widthStep; // Size of aligned image row in bytes.  (Contains the number of bytes between points in the same column and successive rows.)
//    int orderMode[4]; // Unused.
//    int orderConst[4]; // Unused.
//    char* imageDataOrigin; // Unused.
}
IplImage_base;

// NOTE: Not as ridiculous as it seems.
int convert_bit_depth_to_byte_depth(int bit_depth);
// Calculate byte offset from initial byte of IplImage_base to initial byte of IplImage_base ROI.
uint_address_type calc_byte_offset_of_IplImage_base_ROI_widthStep(IplImage_base *image_base, IplImage_base *image_ROI_base);
// Calculate element offset from initial element of IplImage_base to initial element of IplImage_base ROI.
uint_address_type calc_elem_offset_of_IplImage_base_ROI_widthStep(IplImage_base *image_base, IplImage_base *image_ROI_base);

void cvSetZero_base(IplImage_base *image);

typedef struct
{
    real_type x;
    real_type y;
    real_type width;
    real_type height;
}
rect_2DR_type;

void calc_center_rect_2DR(rect_2DR_type *rect, real_type *center_x, real_type *center_y);

////////////////////////////////////////////////////////////////////////////////
// END: Image utilities. }
////////////////////////////////////////////////////////////////////////////////

void inv_rot_matrix_2x2(real_type rot_matrix_2x2[2][2], real_type inv_rot_matrix_2x2[2][2]);
void calc_rot_matrix_2x2(real_type angle, real_type rot_matrix_2x2[2][2]);
int test_rotation_matrix_2x2_for_subsampling_threshold(real_type rot_matrix_2x2[2][2]);
int test_rotation_matrix_2x2_for_subsampling(real_type rot_matrix_2x2[2][2]);

// NOTE: Used for returning multiple values.
typedef struct
{
    int center0_x;
    int center0_y;
    int center_x;
    int center_y;
    int i0_x_start;
    int i0_y_start;
    int i0_x_stop;
    int i0_y_stop;
    int n_i0_x; // NOTE: Equal to n_i_x.  Slightly redundant.
    int n_i0_y; // NOTE: Equal to n_i_y.  Slightly redundant.
    int i_x_start;
    int i_y_start;
}
center_ROI_params_type;

center_ROI_params_type calc_center_ROI_params(int n0_x, int n0_y, int n_x, int n_y);

////////////////////////////////////////////////////////////////////////////////
// END: Image/matrix utilities. }
////////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage. */
}
#endif

#endif // BASE_UTILITIES_H
