#include "utilities_base.h"

// BEGIN: MEX includes
// NOTE: Leave this in the .c file.

#ifdef MEX
#include "mex.h"
#endif

// END: MEX includes

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Output utilities. {
////////////////////////////////////////////////////////////////////////////////

int Print(char* format, ...)
{
    va_list vargs;
    int retval; 

    va_start(vargs, format);
    retval = vprintf(format, vargs); 
    va_end(vargs);

#ifdef MEX
    const int mexPrintf_buffer_max_length = 1024;
    char mexPrintf_buffer[mexPrintf_buffer_max_length];

    va_start(vargs, format);
    retval = vsnprintf(mexPrintf_buffer, mexPrintf_buffer_max_length, format, vargs); 
    va_end(vargs);

    mexPrintf(mexPrintf_buffer);
#endif

    return retval; 
}

void Warn(char* format, ...)
{
    const int format_buffer_max_length = 1024;
    char format_buffer[format_buffer_max_length];

    if ((format != NULL) && (strlen(format) > 0))
        snprintf(format_buffer, format_buffer_max_length, "WARNING: %s", format);
    else
        snprintf(format_buffer, format_buffer_max_length, "WARNING.");

    va_list vargs;
    int retval; 

    va_start(vargs, format);
    retval = vprintf(format_buffer, vargs); 
    va_end(vargs);

#ifdef MEX
    const int mexPrintf_buffer_max_length = 1024;
    char mexPrintf_buffer[mexPrintf_buffer_max_length];

    va_start(vargs, format);
    retval = vsnprintf(mexPrintf_buffer, mexPrintf_buffer_max_length, format_buffer, vargs); 
    va_end(vargs);

    mexPrintf(mexPrintf_buffer);
#endif
}

void Warn_with_info(const char *file, int line, const char *func, char *format, ...)
{
    const int format_buffer_max_length = 1024;
    char format_buffer[format_buffer_max_length];

    if ((format != NULL) && (strlen(format) > 0))
        snprintf(format_buffer, format_buffer_max_length, "WARNING: %s:%i:%s() -- %s", file, line, func, format);
    else
        snprintf(format_buffer, format_buffer_max_length, "WARNING: %s:%i:%s().", file, line, func);

    va_list vargs;
    int retval; 

    va_start(vargs, format);
    retval = vprintf(format_buffer, vargs); 
    va_end(vargs);

#ifdef MEX
    const int mexPrintf_buffer_max_length = 1024;
    char mexPrintf_buffer[mexPrintf_buffer_max_length];

    va_start(vargs, format);
    retval = vsnprintf(mexPrintf_buffer, mexPrintf_buffer_max_length, format_buffer, vargs); 
    va_end(vargs);

    mexPrintf(mexPrintf_buffer);
#endif
}

void Abort(char* format, ...)
{
    const int format_buffer_max_length = 1024;
    char format_buffer[format_buffer_max_length];

    if ((format != NULL) && (strlen(format) > 0))
        snprintf(format_buffer, format_buffer_max_length, "ABORTED: %s", format);
    else
        snprintf(format_buffer, format_buffer_max_length, "ABORTED.");

    va_list vargs;
    int retval; 

    va_start(vargs, format);
    retval = vprintf(format_buffer, vargs); 
    va_end(vargs);

#ifndef MEX
    //memset((void *) 1, 13, 666); // Seg. fault (useful when debugging).
    exit(0);
#else
    const int mexErrMsgTxt_buffer_max_length = 1024;
    char mexErrMsgTxt_buffer[mexErrMsgTxt_buffer_max_length];

    va_start(vargs, format);
    retval = vsnprintf(mexErrMsgTxt_buffer, mexErrMsgTxt_buffer_max_length, format_buffer, vargs); 
    va_end(vargs);

    mexErrMsgTxt(mexErrMsgTxt_buffer); // ORIGINAL.
    //mexPrintf(mexErrMsgTxt_buffer); // DEBUG.
#endif
}

void Abort_with_info(const char *file, int line, const char *func, char *format, ...)
{
    const int format_buffer_max_length = 1024;
    char format_buffer[format_buffer_max_length];

    if ((format != NULL) && (strlen(format) > 0))
        snprintf(format_buffer, format_buffer_max_length, "ABORTED: %s:%i:%s() -- %s", file, line, func, format);
    else
        snprintf(format_buffer, format_buffer_max_length, "ABORTED: %s:%i:%s().", file, line, func);

    va_list vargs;
    int retval; 

    va_start(vargs, format);
    retval = vprintf(format_buffer, vargs); 
    va_end(vargs);

#ifndef MEX
    //memset((void *) 1, 13, 666); // Seg. fault (useful when debugging).
    exit(0);
#else
    const int mexErrMsgTxt_buffer_max_length = 1024;
    char mexErrMsgTxt_buffer[mexErrMsgTxt_buffer_max_length];

    va_start(vargs, format);
    retval = vsnprintf(mexErrMsgTxt_buffer, mexErrMsgTxt_buffer_max_length, format_buffer, vargs); 
    va_end(vargs);

    mexErrMsgTxt(mexErrMsgTxt_buffer); // ORIGINAL.
    //mexPrintf(mexErrMsgTxt_buffer); // DEBUG.
#endif
}

// QUESTION: Why does declaring matrix const cause this warning?
// filter_1D.c:calc_rot_matrix_2x2(), line 1 and 8: warning: passing argument 2 of ‘print_matrix_2x2’ from incompatible pointer type
void print_matrix_2x2(const char *name, real_type matrix[2][2])
{
    Print("%s[0][0-1] = (%.15e, %.15e)\n", name, matrix[0][0], matrix[0][1]);
    Print("%s[1][0-1] = (%.15e, %.15e)\n", name, matrix[1][0], matrix[1][1]);
}

void print_matrix_1D(const char *name, real_type *matrix, int n_x1)
{
    Print("%s = [\n", name);
    for (int i_x1 = 0; i_x1 < n_x1; ++i_x1)
        Print("%.15e\n", matrix[i_x1]);
    Print("]\n");
}

void print_matrix_2D(const char *name, real_type *matrix, int n_x1, int n_x2)
{
    Print("%s = [\n", name);
    for (int i_x1 = 0; i_x1 < n_x1; ++i_x1) // NOTE: We want the rows to be horizontal.
    {
        for (int i_x2 = 0; i_x2 < n_x2; ++i_x2)
        {
            Print("%.15e%s", MATRIX_INDEX_2D(matrix, n_x2, i_x1, i_x2), (i_x2 < n_x2 - 1) ? " " : "");
        }
        Print("\n");
    }
    Print("]\n");
}

void print_matrix_3D(const char *name, real_type *matrix, int n_x1, int n_x2, int n_x3)
{
    Print("%s = [\n", name);
    for (int i_x3 = 0; i_x3 < n_x3; ++i_x3)
    {
        Print("plane[%i] = [\n", i_x3);
        for (int i_x1 = 0; i_x1 < n_x1; ++i_x1) // NOTE: We want the rows to be horizontal.
        {
            for (int i_x2 = 0; i_x2 < n_x2; ++i_x2)
            {
                Print("%.15e%s", MATRIX_INDEX_3D(matrix, n_x1, n_x2, i_x1, i_x2, i_x3), (i_x2 < n_x2 - 1) ? " " : "");
            }
            Print("\n");
        }
        Print("]\n");
    }
    Print("]\n");
}

////////////////////////////////////////////////////////////////////////////////
// END: Output utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Image/matrix utilities. {
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Image utilities. {
////////////////////////////////////////////////////////////////////////////////

// NOTE: Not as ridiculous as it seems.
int convert_bit_depth_to_byte_depth(int bit_depth)
{
    if (bit_depth < 8)
        ABORT("bit_depth == %d < 0\n", bit_depth);

    if (bit_depth % 8 != 0) 
        ABORT("bit_depth %% 8 == %d != 0\n", bit_depth % 8);

    return (bit_depth / 8);
}

// Calculate byte offset from initial byte of IplImage_base to initial byte of IplImage_base ROI.
uint_address_type calc_byte_offset_of_IplImage_base_ROI_widthStep(IplImage_base *image_base, IplImage_base *image_ROI_base)
{
    uint_address_type image_base_imageData_start = (uint_address_type) image_base->imageData;
    uint_address_type image_ROI_base_imageData_start = (uint_address_type) image_ROI_base->imageData;
    if (image_base_imageData_start > image_ROI_base_imageData_start)
        ABORT("image_base_imageData_start == %llu > %llu == image_ROI_base_imageData_start\n", image_base_imageData_start, image_ROI_base_imageData_start);
    uint_address_type image_ROI_base_byte_offset = image_ROI_base_imageData_start - image_base_imageData_start;
    return image_ROI_base_byte_offset;
}

// Calculate element offset from initial element of IplImage_base to initial element of IplImage_base ROI.
uint_address_type calc_elem_offset_of_IplImage_base_ROI_widthStep(IplImage_base *image_base, IplImage_base *image_ROI_base)
{
    uint_address_type image_ROI_base_byte_offset = calc_byte_offset_of_IplImage_base_ROI_widthStep(image_base, image_ROI_base);
    uint_address_type byte_depth = convert_bit_depth_to_byte_depth(image_ROI_base->depth);
    if (image_ROI_base_byte_offset % byte_depth != 0)
        ABORT("(image_ROI_base_byte_offset == %llu %% byte_depth == %llu) == %llu != 0", image_ROI_base_byte_offset, byte_depth, image_ROI_base_byte_offset % byte_depth);
    uint_address_type image_ROI_base_elem_offset = image_ROI_base_byte_offset / byte_depth;
    return image_ROI_base_elem_offset;
}

void cvSetZero_base(IplImage_base *image)
{
    memset(image->imageData, 0, image->imageSize);
}

void calc_center_rect_2DR(rect_2DR_type *rect, real_type *center_x, real_type *center_y)
{
    *center_x = (rect->x + rect->width) / 2.0;
    *center_y = (rect->y + rect->height) / 2.0;
}

////////////////////////////////////////////////////////////////////////////////
// END: Image utilities. }
////////////////////////////////////////////////////////////////////////////////

void inv_rot_matrix_2x2(real_type rot_matrix_2x2[2][2], real_type inv_rot_matrix_2x2[2][2])
{
    //print_matrix_2x2("rot_matrix_2x2", rot_matrix_2x2);

    inv_rot_matrix_2x2[0][0] = rot_matrix_2x2[0][0];
    inv_rot_matrix_2x2[1][0] = rot_matrix_2x2[0][1];
    inv_rot_matrix_2x2[0][1] = rot_matrix_2x2[1][0];
    inv_rot_matrix_2x2[1][1] = rot_matrix_2x2[1][1];

    //print_matrix_2x2("inv_rot_matrix_2x2", inv_rot_matrix_2x2);
}

void calc_rot_matrix_2x2(real_type angle, real_type rot_matrix_2x2[2][2])
{
    rot_matrix_2x2[0][0] = cos(angle); rot_matrix_2x2[0][1] = -sin(angle);
    rot_matrix_2x2[1][0] = sin(angle); rot_matrix_2x2[1][1] = cos(angle);

    //print_matrix_2x2("rot_matrix_2x2", rot_matrix_2x2);
}

// If rotation is small or close to +/-90 or +/-180, round elements to nearest
// integer and indicate the use of a smaller step.
int test_rotation_matrix_2x2_for_subsampling_threshold(real_type rot_matrix_2x2[2][2])
{
    real_type max;

    // If rotation is small.
    max = fmax(fabs(rot_matrix_2x2[0][0] - 1.0), fabs(rot_matrix_2x2[0][1]));
    max = fmax(max, fabs(rot_matrix_2x2[1][0])); max = fmax(max, fabs(rot_matrix_2x2[1][1] - 1.0));
//Print("%.15e\n", max);

    if (max < EPSILON)
    {
        rot_matrix_2x2[0][0] = 1.0; rot_matrix_2x2[0][1] = 0.0;
	rot_matrix_2x2[1][0] = 0.0; rot_matrix_2x2[1][1] = 1.0;
        return 0;
    }

    // If rotation is close to +90.
    max = fmax(fabs(rot_matrix_2x2[0][0]), fabs(rot_matrix_2x2[0][1] - 1.0));
    max = fmax(max, fabs(rot_matrix_2x2[1][0] + 1.0)); max = fmax(max, fabs(rot_matrix_2x2[1][1]));
//Print("%.15e\n", max);

    if (max < EPSILON)
    {
        rot_matrix_2x2[0][0] = 0.0; rot_matrix_2x2[0][1] = 1.0;
	rot_matrix_2x2[1][0] = -1.0; rot_matrix_2x2[1][1] = 0.0;
        return 0;
    }

    // If rotation is close to -90.
    max = fmax(fabs(rot_matrix_2x2[0][0]), fabs(rot_matrix_2x2[0][1] + 1.0));
    max = fmax(max, fabs(rot_matrix_2x2[1][0] - 1.0)); max = fmax(max, fabs(rot_matrix_2x2[1][1]));
//Print("%.15e\n", max);

    if (max < EPSILON)
    { 
        rot_matrix_2x2[0][0] = 0.0; rot_matrix_2x2[0][1] = -1.0;
	rot_matrix_2x2[1][0] = 1.0; rot_matrix_2x2[1][1] = 0.0;
        return 0;
    }

    // If rotation is close to +/-180.
    max = fmax(fabs(rot_matrix_2x2[0][0] + 1.0), fabs(rot_matrix_2x2[0][1]));
    max = fmax(max, fabs(rot_matrix_2x2[1][0])); max = fmax(max, fabs(rot_matrix_2x2[1][1] + 1.0));
//Print("%.15e\n", max);

    if (max < EPSILON)
    {
        rot_matrix_2x2[0][0] = -1.0; rot_matrix_2x2[0][1] = 0.0;
	rot_matrix_2x2[1][0] = 0.0; rot_matrix_2x2[1][1] = -1.0;
        return 0;
    }

    return 1;
}

// As above, but test only: do not threshold the matrix values.
int test_rotation_matrix_2x2_for_subsampling(real_type rot_matrix_2x2[2][2])
{
    real_type rot_matrix_copy[2][2];
    memcpy(rot_matrix_copy, rot_matrix_2x2, sizeof(real_type) * 2 * 2);

    return test_rotation_matrix_2x2_for_subsampling_threshold(rot_matrix_copy);
}

center_ROI_params_type calc_center_ROI_params(int n0_x, int n0_y, int n_x, int n_y)
{
    //Print("utilities_base.c:calc_center_ROI_params(), n0_(xy) = (%i, %i)\n", n0_x, n0_y);
    //Print("utilities_base.c:calc_center_ROI_params(), n_(xy) = (%i, %i)\n", n_x, n_y);

    //// NOTE: Assume n0_(xy) = (1, 5) and consider n_(xy) = (1, 2) and n_(xy) = (1, 3);
    //int n_x_parity = n_x % 2;
    //int n_y_parity = n_y % 2;

    center_ROI_params_type params;

    params.center0_x = n0_x / 2;
    params.center0_y = n0_y / 2;

    params.center_x = n_x / 2;
    params.center_y = n_y / 2;

    if (n0_x <= n_x)
    {
        params.i0_x_start = 0;
        params.i0_x_stop = n0_x; // NOTE: Iterate to element params.i0_x_stop - 1.

        params.i_x_start = params.center_x - params.center0_x;
    }
    else
    {
        params.i0_x_start = params.center0_x - params.center_x;
        //params.i0_x_stop = params.center0_x + params.center_x + n_x_parity; // NOTE: Iterate to element params.i0_x_stop - 1.
        params.i0_x_stop = params.i0_x_start + n_x; // NOTE: Iterate to element params.i0_x_stop - 1.

        params.i_x_start = 0;
    }

    if (n0_y <= n_y)
    {
        params.i0_y_start = 0;
        params.i0_y_stop = n0_y; // NOTE: Iterate to element i0_y_stop - 1.

        params.i_y_start = params.center_y - params.center0_y;
    }
    else
    {
        params.i0_y_start = params.center0_y - params.center_y;
        //params.i0_y_stop = params.center0_y + params.center_y + n_y_parity; // NOTE: Iterate to element i0_y_stop - 1.
        params.i0_y_stop = params.i0_y_start + n_y; // NOTE: Iterate to element i0_y_stop - 1.

        params.i_y_start = 0;
    }

    params.n_i0_x = params.i0_x_stop - params.i0_x_start;
    params.n_i0_y = params.i0_y_stop - params.i0_y_start;

//    Print("params.i0_x_start = %i\n", params.i0_x_start);
//    Print("params.i0_x_stop = %i\n", params.i0_x_stop);
//    Print("params.i_x_start = %i\n", params.i_x_start);
//    Print("params.n_i0_x = %i\n", params.n_i0_x);
//    Print("params.i0_y_start = %i\n", params.i0_y_start);
//    Print("params.i0_y_stop = %i\n", params.i0_y_stop);
//    Print("params.i_y_start = %i\n", params.i_y_start);
//    Print("params.n_i0_y = %i\n", params.n_i0_y);

    return params;
}

////////////////////////////////////////////////////////////////////////////////
// END: Image/matrix utilities. }
////////////////////////////////////////////////////////////////////////////////
