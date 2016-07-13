#ifndef UTILITIES_H
#define UTILITIES_H

#include "utilities_base.h"

/*
// BEGIN: FFTW includes

#include <fftw3.h>

// BEGIN: These are assumed to represent pixel values.

typedef fftwf_complex fftw_complex_type; // WARNING: float vs. double
//typedef fftw_complex fftw_complex_type; // WARNING: float vs. double

// END: These are assumed to represent pixel values.

// END: FFTW includes
*/

// BEGIN: IMOD includes

#include "mrcslice.h"

// END: IMOD includes

// BEGIN: OpenCV includes

#define CV_NO_BACKWARD_COMPATIBILITY
#include "cv.h"

// BEGIN: These are assumed to represent pixel values.

#define CVMAT_TYPE_REAL CV_32FC1 // WARNING: float vs. double
#define CVMAT_TYPE_COMPLEX CV_32FC2 // WARNING: float vs. double
#define IPL_DEPTH_REAL IPL_DEPTH_32F // WARNING: float vs. double

//#define CVMAT_TYPE_REAL CV_64FC1 // WARNING: float vs. double
//#define CVMAT_TYPE_COMPLEX CV_64FC2 // WARNING: float vs. double
//#define IPL_DEPTH_REAL IPL_DEPTH_64F // WARNING: float vs. double

// END: These are assumed to represent pixel values.

// END: OpenCV includes

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Parsing utilities. {
////////////////////////////////////////////////////////////////////////////////

void create_indexing_array_from_indexing_string
(
    int n_indexed_max,
    const char *indices_string,
    int *n_indexed,
    int **indices
);

void create_indexing_array_from_indexing_string_ordered_no_duplicates
(
    int n_indexed_max,
    const char *indices_string,
    int *n_indexed,
    int **indices
);

////////////////////////////////////////////////////////////////////////////////
// END: Parsing utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Image/matrix utilities. {
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: IMOD utilities. {
////////////////////////////////////////////////////////////////////////////////

void calc_mmm_of_MRC(MrcHeader *header); // Calculate mean, median, and mode and set.

void sliceWriteMRC(struct MRCheader *header, Islice *slice, int index, char axis);

Islice *sliceCreateSet(int x_size, int y_size, int mode, float val); // Create an Islice set to val.
Islice *sliceCreateZeros(int x_size, int y_size, int mode);

Islice *sliceNewMode_copy(Islice *slice, int mode);

void createNewMRCFile
(
    const char *filepath,
    int n_x, int n_y, int n_z,
    int mode
);

// Open an existing MRC file and file header.
void openMRCFile_general
(
    const char *filepath,
    MrcHeader *header,
    FILE **file
);

void get_MRC_sizes(const char *filepath, int *n_x, int *n_y, int *n_z);

void write_MRC_image_from_Islice(Islice *slice, const char *filepath);

void slice_resize_copy
(
    Islice *slice0,
    Islice *slice,
    real_type resize_factor
);

// If slice0 extent < slice extent, excess slice is filled with zeros.
// If slice0 extent > slice extent, slice is a crop of slice0.
// NOTE: Assumes SLICE_MODE_FLOAT.
void slice_recanvas_center_copy
(
    Islice *slice0,
    Islice *slice
);

// If slice0 extent < slice extent, excess slice (i.e., right and/or top) is filled with zeros.
// If slice0 extent > slice extent, slice is a crop (i.e., right and/or top) of slice0.
// NOTE: Assumes SLICE_MODE_FLOAT.
void slice_recanvas_copy
(
    Islice *slice0,
    Islice *slice
);

////////////////////////////////////////////////////////////////////////////////
// END: IMOD utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// END: OpenCV utilities. {
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: OpenCV utilities: IplImage. {
////////////////////////////////////////////////////////////////////////////////

void check_IplImage_depth_vs_size(IplImage *image, size_t size);

// BEGIN: malloc() and cvReleaseImageHeader().

/* 
NOTE:

Yes, this is confusing.  The length in bytes of a row in an OpenCV-allocated
CvMat or IplImage is not necessarily the "cols" or "width" member multiplied by
the element byte size because, for machine efficiency, matrix and image
allocation is done to the nearest four-byte boundary.
*/

// BEGIN: Use on pre-existing data; no copy.  If OpenCV allocated data, do not use free() on the actual data.
IplImage *create_IplImage_32F_from_float_data(int n_x, int n_y, float *data);
IplImage *create_IplImage_64F_from_double_data(int n_x, int n_y, double *data);
IplImage *create_IplImage_XXF_from_pixel_type_data(int n_x, int n_y, pixel_type *data);
// END: Use on pre-existing data; no copy.  If OpenCV allocated data, do not use free() on the actual data.

void free_imageData_cvReleaseImageHeader(IplImage *image); // Calls free(image->imageData) and cvReleaseImageHeader(image).

// BEGIN: Use free_image_data_cvReleaseImageHeader().
IplImage *create_IplImage_malloc(int n_x, int n_y, int ipl_depth);
IplImage *create_IplImage_from_float_data_copy_malloc(int n_x, int n_y, float *data, int ipl_depth);
IplImage *create_IplImage_from_Islice_copy_malloc(Islice *slice, int ipl_depth);
IplImage *create_IplImage_from_sliceReadMRC_malloc(MrcHeader *file_header, int i_tilt, char axis, int ipl_depth);
IplImage *create_IplImage_SetZero_malloc(int n_x, int n_y, int ipl_depth);
// END: Use free_image_data_cvReleaseImageHeader().
// END: malloc() and cvReleaseImageHeader().

// BEGIN: Use cvReleaseImage().
IplImage *create_IplImage_from_float_data_copy(int n_x, int n_y, float *data, int ipl_depth);
IplImage *create_IplImage_from_Islice_copy(Islice *slice, int ipl_depth);
IplImage *create_IplImage_from_sliceReadMRC(MrcHeader *file_header, int i_tilt, char axis, int ipl_depth);
IplImage *create_IplImage_SetZero(int n_x, int n_y, int ipl_depth);
// END: Use cvReleaseImage().

// BEGIN: IplImage_base utilities. {

// WARNING: Do NOT free/release imageData of IplImage_base.

// WARNING: Will not work with ROI widthSteps! (Not true...)
void copy_IplImage_to_IplImage_base(IplImage *image, IplImage_base *image_base);

// Initializes IplImage_base with pointer to existing ROI data inside a larger image.
// WARNING: This does not appear to be necessary.
void copy_IplImage_ROI_widthStep_to_IplImage_base_ROI_widthStep
(
    IplImage *image, // The original IplImage.
    IplImage *image_ROI, // An ROI widthStep of the original IplImage.
    IplImage_base *image_ROI_base // The output IplImage_base.
);

// END: IplImage_base utilities. }

// BEGIN: Assumes data is already allocated.

// WARNING: IPL_DEPTH_32F and IPL_DEPTH_64F only!
// WARNING: SLICE_MODE_FLOAT only!
// If image extent < slice extent, excess slice (i.e., right and/or top) is filled with zeros.
// If image extent > slice extent, slice is a crop (i.e., right and/or top) of image.
void copy_IplImage_to_Islice(IplImage *image, Islice *slice);

// WARNING: IPL_DEPTH_32F and IPL_DEPTH_64F only!
// WARNING: SLICE_MODE_FLOAT only!
// If slice extent < image extent, excess image about center is filled with zeros.
// If slice extent > image extent, image is a crop about center of slice.
void copy_IplImage_to_Islice_center(IplImage *image, Islice *slice);

// WARNING: IPL_DEPTH_32F and IPL_DEPTH_64F only!
// If slice extent < image extent, excess image (i.e., right and/or top) is filled with zeros.
// If slice extent > image extent, image is a crop (i.e., right and/or top) of slice.
void copy_Islice_to_IplImage(Islice *slice, IplImage *image);

// WARNING: IPL_DEPTH_32F and IPL_DEPTH_64F only!
// If slice extent < image extent, excess image about center is filled with zeros.
// If slice extent > image extent, image is a crop about center of slice.
void copy_Islice_to_IplImage_center(Islice *slice, IplImage *image);

// WARNING: IPL_DEPTH_32F and IPL_DEPTH_64F only!
// NOTE: Uses widthStep to define a subimage (i.e., an "ROI widthStep") pointing into image->imageData.
// WARNING: Do NOT free/release imageData of IplImage * returned by this function.
IplImage *create_IplImage_ROI_widthStep_of_IplImage(IplImage *image, CvRect rect);

// END: Assumes data is already allocated.

void write_MRC_image_from_IplImage(IplImage *image, const char *filepath);

////////////////////////////////////////////////////////////////////////////////
// END: OpenCV utilities: IplImage. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: OpenCV utilities: CvMat. {
////////////////////////////////////////////////////////////////////////////////

CvMat *create_CvMat_32F_from_float_data(int rows, int cols, float *data);
CvMat *create_CvMat_64F_from_double_data(int rows, int cols, double *data);

void copy_CvMat_to_Islice(CvMat *matrix, Islice *slice);
CvMat *create_CvMat_from_Islice_copy(Islice *slice, int type);

// WARNING: Allocates an IplImage header that must be freed by the caller with cvReleaseImageHeader().
// WARNING: IplImage's imageData points to matrix->ptr.
IplImage *create_IplImage_from_CvMat(CvMat *matrix);

void free_data_cvReleaseMat(CvMat *matrix);

void print_CvMat_packed_complex(const char *name, const CvMat *matrix);
void print_CvMat(const char *name, const CvMat *matrix);

void write_MRC_image_from_CvMat(CvMat *matrix, const char *filepath);

////////////////////////////////////////////////////////////////////////////////
// END: OpenCV utilities: CvMat. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// END: OpenCV utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// END: Image/matrix utilities. }
////////////////////////////////////////////////////////////////////////////////

#endif // UTILITIES_H
