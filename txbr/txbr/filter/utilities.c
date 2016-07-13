#include "utilities.h"

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Parsing utilities. {
////////////////////////////////////////////////////////////////////////////////

void create_indexing_array_from_indexing_string
(
    int n_indexed_max,
    const char *indices_string,
    int *n_indexed,
    int **indices
)
{
    if (n_indexed_max < 0)
        ABORT("n_indexed_max == %i < 0.\n", n_indexed_max);

    int *indices_temp = calloc(n_indexed_max, sizeof(int));
    if (indices_temp == NULL)
        ABORT("Cannot allocate memory for int *indices_temp.\n");

    *n_indexed = 0;

    int indices_string_length = 0;
    if (indices_string != NULL)
        indices_string_length = strlen(indices_string);
    int i_indices_string = 0;

    int index_buffer = 0;
    int n_chars_scanned = 0;

    int state = 0; // valid state >= 0
    while (i_indices_string < indices_string_length)
    {
        switch (state)
        {
            case 0: // Initial state.
            {
                if ((indices_string_length == 1) && (indices_string[0] == '-'))
                    --indices_string_length;
                else if (!(isdigit(indices_string[0])))
                    ABORT("Bad indexing string = \"%s\", state == %i.\n", indices_string, state);

                state = 1;

                break;
            }
            case 1: // Parse a single integer index.
            {
                int retval = sscanf(&(indices_string[i_indices_string]), "%i%n", 
                       &index_buffer, &n_chars_scanned); 

                if (retval == 0)
                    ABORT("Bad indexing string = \"%s\", state == %i.\n", indices_string, state);

                if ((index_buffer < 0) || (index_buffer >= n_indexed_max))
                    ABORT("Bad indexing string = \"%s\", state == %i.\n", indices_string, state);

                indices_temp[*n_indexed] = index_buffer;
                ++(*n_indexed);
                i_indices_string += n_chars_scanned;
                state = 2;

                break;
            }
            case 2: // Comma delimiter or hyphen delimiter.
            {
                if (i_indices_string >= indices_string_length - 1)
                    ABORT("Bad indexing string = \"%s\", state == %i.\n", indices_string, state);

                if (indices_string[i_indices_string] == ',')
                {
                    ++i_indices_string;
                    state = 1;
                }
                else if (indices_string[i_indices_string] == '-')
                {
                    ++i_indices_string;
                    state = 3;
                    break;
                }
                else
                    ABORT("Bad indexing string = \"%s\", state == %i.\n", indices_string, state);

                break;
            }
            case 3: // Parse the second integer index of a hyphenated range.
            {
                if ((*n_indexed) < 1) 
                    ABORT("Bad indexing string = \"%s\", state == %i.\n", indices_string, state);

                int min_index = index_buffer;

                int retval = sscanf(&(indices_string[i_indices_string]), "%i%n", 
                       &index_buffer, &n_chars_scanned); 
                int max_index = index_buffer;

                if (retval == 0)
                    ABORT("Bad indexing string = \"%s\", state == %i.\n", indices_string, state);

                if ((max_index >= n_indexed_max) ||
                    (max_index <= 0) ||
                    (max_index <= min_index) ||
                    ((max_index - min_index + (*n_indexed)) >= n_indexed_max))
                    ABORT("Bad indexing string = \"%s\", state == %i.\n", indices_string, state);

                i_indices_string += n_chars_scanned;

                while (++min_index <= max_index) 
                {
                    indices_temp[*n_indexed] = min_index;
                    ++(*n_indexed);
                }

                if ((i_indices_string < indices_string_length) &&
                    (indices_string[i_indices_string] != ','))
                    ABORT("Bad indexing string (indices_string[i_indices_string] != ',') = \"%s\", state == %i.\n", indices_string, state);

                state = 2;

                break;
            }
            case -1: // Failure.
            default:
                ABORT("Bad indexing string = \"%s\", invalid state == %i.\n", indices_string, state);
        }
    }

    //ABORT("*n_indexed = %i.\n", *n_indexed);

    if (indices_string_length == 0)
        *n_indexed = n_indexed_max;

    *indices = calloc(*n_indexed, sizeof(int));
    if (*indices == NULL)
        ABORT("Cannot allocate memory for int **indices.\n");

    if (indices_string_length == 0)
    {
        for (int i_index = 0; i_index < *n_indexed; ++i_index)
            (*indices)[i_index] = i_index;
    }
    else
    {
        memcpy((*indices), indices_temp, sizeof(int) * (*n_indexed));
    }

    free(indices_temp);
}

void create_indexing_array_from_indexing_string_ordered_no_duplicates
(
    int n_indexed_max,
    const char *indices_string,
    int *n_indexed,
    int **indices
)
{
    if (n_indexed_max < 0)
        ABORT("n_indexed_max == %i < 0.\n", n_indexed_max);

    int *indices_touched = calloc(n_indexed_max, sizeof(int));
    if (indices_touched == NULL)
        ABORT("Cannot allocate memory for int indices_touched[].\n");

    int n_indexed_temp;
    int *indices_temp;

    create_indexing_array_from_indexing_string
    (
        n_indexed_max,
        indices_string,
        &n_indexed_temp,
        &indices_temp
    );

    *n_indexed = 0; // NOTE: Overkill, as n_indexed_temp, if valid, is set to the correct value.

    for (int i_index = 0; i_index < n_indexed_temp; ++i_index)
    {
        int index_buffer = indices_temp[i_index];

        //Print("i_index = %i\n", i_index);
        //Print("index_buffer = %i\n", index_buffer);

        if (index_buffer >= n_indexed_max)
            ABORT("index_buffer == %i >= %i == n_indexed_max", index_buffer, n_indexed_max);

        if ((indices_touched[index_buffer])++ > 0)
            ABORT("indices_touched[%i] == %i", index_buffer, indices_touched[index_buffer]);

        ++(*n_indexed);
    }

    free(indices_temp);

    *indices = calloc(*n_indexed, sizeof(int));
    if (*indices == NULL)
        ABORT("Cannot allocate memory for int **indices.\n");

    int i_index = 0;
    for (int i_index_touched = 0; i_index_touched < n_indexed_max; ++i_index_touched)
        if (indices_touched[i_index_touched] != 0)
            (*indices)[i_index++] = i_index_touched;

    free(indices_touched);
}

////////////////////////////////////////////////////////////////////////////////
// END: Parsing utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Image/matrix utilities. {
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: IMOD utilities. {
////////////////////////////////////////////////////////////////////////////////

void calc_mmm_of_MRC(MrcHeader *header)
{
    real_type min = 0.0, max = 0.0, mean = 0.0, m2 = 0.0;

    int n_z = header->nz;
    int length = header->nx * header->ny;

    for (int i_z = 0; i_z < n_z; ++i_z) 
    {
        Islice *slice = NULL;
        slice = sliceReadMRC(header, i_z, 'z');
        if (slice == NULL)
            ABORT("Cannot read file.\n");

        if (slice->mode != SLICE_MODE_FLOAT)
            sliceNewMode(slice, SLICE_MODE_FLOAT);

        float *array = slice->data.f;

        for (int i = 0; i < length; ++i)
        {
            min = fmin(min, array[i]);
            max = fmax(max, array[i]);
            mean += array[i];
            m2 += array[i] * array[i];
        }

        sliceFree(slice);
    }
    
    int n_pixels = header->nx * header->ny * header->nz;

    mean = mean / (real_type) n_pixels;
    m2 = (m2 / (real_type) n_pixels) - (mean * mean);

//    header->amin = mean - (4.0 * sqrt(m2));
//    header->amax = mean + (4.0 * sqrt(m2));
    header->amin = min;
    header->amax = max;
    header->amean = mean;

    mrc_head_write(header->fp, header);

// PRIMARY DEBUG:     Print("header->amin = %f, header->amax = %f, header->amean = %f, m2 = %f\n", header->amin, header->amax, header->amean, m2); // PRIMARY DEBUG.
}

void sliceWriteMRC(struct MRCheader *header, Islice *slice, int index, char axis)
{
    if (header->mode != slice->mode)
    {
        sliceNewMode(slice, header->mode);
    }

    void *data = NULL;
    switch (slice->mode)
    {
        case (SLICE_MODE_FLOAT):

            data = (void *) slice->data.f;
            break;

        case (SLICE_MODE_BYTE):

            data = (void *) slice->data.b;
            break;

        case (SLICE_MODE_SHORT):

            data = (void *) slice->data.s;
            break;

        case (SLICE_MODE_USHORT):

            data = (void *) slice->data.us;
            break;

        default:

            ABORT("Unrecognized MRC data mode: %i.\n", slice->mode);
    }

// PRIMARY DEBUG:     Print("slice->xsize = %i\n", slice->xsize); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("slice->ysize = %i\n", slice->ysize); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("slice->mode = %i\n", slice->mode); // PRIMARY DEBUG.

// PRIMARY DEBUG:     Print("header->nx = %i\n", header->nx); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("header->ny = %i\n", header->ny); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("header->nz = %i\n", header->nz); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("header->mode = %i\n", header->mode); // PRIMARY DEBUG.

    if (mrc_write_slice(data, header->fp, header, index, axis) != 0)
        ABORT("Cannot write file.\n");
}

Islice *sliceCreateSet(int x_size, int y_size, int mode, float val)
{
    if (mode != SLICE_MODE_FLOAT)
        ABORT("mode == %d != %d == SLICE_MODE_FLOAT.\n", mode, SLICE_MODE_FLOAT);

    Islice *slice = sliceCreate(x_size, y_size, mode);
    if (slice == NULL)
        ABORT("Cannot acquire memory for Islice *slice.\n");
    Ival buffer;
    memset(&buffer, 0, sizeof(Ival));
    memcpy(&buffer, &val, sizeof(float));
    sliceClear(slice, buffer);
    return slice;
}

Islice *sliceCreateZeros(int x_size, int y_size, int mode)
{
    return sliceCreateSet(x_size, y_size, mode, 0.0);
}

Islice *sliceNewMode_copy(Islice *slice, int mode)
{
    Islice *slice_copy = sliceCreate(slice->xsize, slice->ysize, slice->mode);
    if (slice_copy == NULL)
        ABORT("Cannot acquire memory for Islice *slice_copy.\n");
    memcpy(slice_copy->data.b, slice->data.b, slice->xsize * slice->ysize * slice->dsize * slice->csize);
    sliceNewMode(slice_copy, mode);

    // WARNING: slice_copy->dsize does not seem to be correct!

/*
    Islice *slice_copy = sliceCreate(slice->xsize, slice->ysize, mode);
    if (slice_copy == NULL)
        ABORT("Cannot acquire memory for Islice *slice_copy.\n");
    Ival val;
    for (int i = 0; i < slice_copy->xsize; ++i)
        for (int j = 0; j < slice_copy->ysize; ++j)
        {
            sliceGetVal(slice, i, j, val);
            slicePutVal(slice_copy, i, j, val);
        }
*/

    return slice_copy;
}

void createNewMRCFile
(
    const char *filepath,
    int n_x, int n_y, int n_z,
    int mode
)
{
// PRIMARY DEBUG:     Print("Creating new MRC file: \"%s\".\n", filepath); // PRIMARY DEBUG.

    FILE* file;

    if ((file = fopen(filepath, "wb")) == NULL)
        ABORT("Cannot open file: \"%s\".\n", filepath);

    struct MRCheader header;
    mrc_head_new(&header, n_x, n_y, n_z, mode);
    mrc_head_write(file, &header);

    if (fclose(file))
        ABORT("Cannot close file: \"%s\".\n", filepath);
}

// Open an existing MRC file and file header.
void openMRCFile_general
(
    const char *filepath,
    MrcHeader *header,
    FILE **file
)
{
// PRIMARY DEBUG:     Print("openMRCFile_general(\"%s\")\n", filepath); // PRIMARY DEBUG.

    if ((*file = fopen(filepath, "r+b")) == NULL)
	if ((*file = fopen(filepath, "rb")) == NULL)
        	ABORT("Cannot open file \"%s\".\n", filepath);

    mrc_head_read(*file, header);
}

void get_MRC_sizes(const char *filepath, int *n_x, int *n_y, int *n_z)
{
    MrcHeader header;
    FILE *file;

    openMRCFile_general(filepath, &header, &file);

    *n_x = header.nx;
    *n_y = header.ny;
    *n_z = header.nz;

    if (fclose(file))
        ABORT("Cannot close file: \"%s\".\n", filepath);
}

void write_MRC_image_from_Islice(Islice *slice, const char *filepath)
{
    // BEGIN: Create and open output MRC file.
    // BEGIN: Create output MRC file.
    MrcHeader file_header;

    FILE *file = NULL;

    createNewMRCFile(filepath, slice->xsize, slice->ysize, 1, SLICE_MODE_FLOAT);
    // END: Create output MRC file.

    // BEGIN: Open output MRC file.
    openMRCFile_general(filepath, &file_header, &file);
    sliceWriteMRC(&file_header, slice, 0, 'z');
    calc_mmm_of_MRC(&file_header);
    // END: Open output MRC file.
    // END: Create and open output MRC file.

    if (fclose(file))
        ABORT("Cannot close file: \"%s\".\n", filepath);
//    exit(0);
}

// NOTE: Internal function.
void slice_resize_copy_cv
(
    Islice *slice0,
    Islice *slice,
    real_type resize_factor
)
{
    if (slice0->mode != SLICE_MODE_FLOAT)
        ABORT("slice0->mode == %d != %d == SLICE_MODE_FLOAT.\n", slice0->mode, SLICE_MODE_FLOAT);

    if (slice->mode != SLICE_MODE_FLOAT)
        ABORT("slice->mode == %d != %d == SLICE_MODE_FLOAT.\n", slice->mode, SLICE_MODE_FLOAT);

    if (slice0 == slice)
        ABORT("slice0 == slice.\n");

//    if ((slice->xsize != round((real_type) slice0->xsize * resize_factor)) || (slice->ysize != round((real_type) slice0->ysize * resize_factor)))
//        ABORT("slice0->(xsize, ysize) = (%i, %i), slice->(xsize, ysize) = (%i, %i), but resize_factor = %.15e.\n", slice0->xsize, slice0->ysize, slice->xsize, slice->ysize, resize_factor);

    IplImage *slice0_cv = create_IplImage_from_Islice_copy(slice0, IPL_DEPTH_REAL); // WARNING: float vs. double -- cvResize() requires IPL_DEPTH_32F.
    IplImage *slice_cv = create_IplImage_from_Islice_copy(slice, IPL_DEPTH_REAL); // WARNING: float vs. double -- cvResize() requires IPL_DEPTH_32F.

    cvResize(slice0_cv, slice_cv, CV_INTER_CUBIC);
    copy_IplImage_to_Islice(slice_cv, slice);

    cvReleaseImage(&slice0_cv);
    cvReleaseImage(&slice_cv);
}

void slice_resize_copy
(
    Islice *slice0,
    Islice *slice,
    real_type resize_factor
)
{
    slice_resize_copy_cv(slice0, slice, resize_factor);
}

// If slice0 extent < slice extent, excess slice is filled with zeros.
// If slice0 extent > slice extent, slice is a crop of slice0.
void slice_recanvas_center_copy
(
    Islice *slice0,
    Islice *slice
)
{
    if (slice0->mode != SLICE_MODE_FLOAT)
        ABORT("slice0->mode == %d != %d == SLICE_MODE_FLOAT.\n", slice0->mode, SLICE_MODE_FLOAT);

    if (slice->mode != SLICE_MODE_FLOAT)
        ABORT("slice->mode == %d != %d == SLICE_MODE_FLOAT.\n", slice->mode, SLICE_MODE_FLOAT);

    if (slice0 == slice)
        ABORT("slice0 == slice.\n");

    float *slice0_data_f = slice0->data.f;
    int n0_x = slice0->xsize;
    int n0_y = slice0->ysize;

    float *slice_data_f = slice->data.f;
    int n_x = slice->xsize;
    int n_y = slice->ysize;

    if ((n0_x < n_x) || (n0_y < n_y))
    {
        float v_clear = 0.0;
        sliceClear(slice, &v_clear);
    }

    center_ROI_params_type params = calc_center_ROI_params(n0_x, n0_y, n_x, n_y);

    float v = 0.0;

    int i_y = params.i_y_start;
    for (int i0_y = params.i0_y_start; i0_y < params.i0_y_stop; ++i0_y)
    {
        int i_x = params.i_x_start;
        for (int i0_x = params.i0_x_start; i0_x < params.i0_x_stop; ++i0_x)
        {
            //Print("(i0_x, i0_y) = (%i, %i)\n", i0_x, i0_y);
            //Print("(i_x, i_y) = (%i, %i)\n", i_x, i_y);

            v = INDEX_2D(slice0_data_f, n0_x, i0_x, i0_y);
            PUTVAL_2D(slice_data_f, n_x, i_x, i_y, v);
            ++i_x;
        }
        ++i_y;
    }
}

// If slice0 extent < slice extent, excess slice (i.e., right and/or top) is filled with zeros.
// If slice0 extent > slice extent, slice is a crop (i.e., right and/or top) of slice0.
void slice_recanvas_copy
(
    Islice *slice0,
    Islice *slice
)
{
    if (slice0->mode != SLICE_MODE_FLOAT)
        ABORT("slice0->mode == %d != %d == SLICE_MODE_FLOAT.\n", slice0->mode, SLICE_MODE_FLOAT);

    if (slice->mode != SLICE_MODE_FLOAT)
        ABORT("slice->mode == %d != %d == SLICE_MODE_FLOAT.\n", slice->mode, SLICE_MODE_FLOAT);

    if (slice0 == slice)
        ABORT("slice0 == slice.\n");

    float *slice0_data_f = slice0->data.f;
    int n0_x = slice0->xsize;
    int n0_y = slice0->ysize;

    float *slice_data_f = slice->data.f;
    int n_x = slice->xsize;
    int n_y = slice->ysize;

    if ((n0_x < n_x) || (n0_y < n_y))
    {
        float v_clear = 0.0;
        sliceClear(slice, &v_clear);
    }

    int i0_x_stop = fmin(n0_x, n_x);
    int i0_y_stop = fmin(n0_y, n_y);

    float v = 0.0;

    for (int i0_y = 0; i0_y < i0_y_stop; ++i0_y)
        for (int i0_x = 0; i0_x < i0_x_stop; ++i0_x)
        {
            v = INDEX_2D(slice0_data_f, n0_x, i0_x, i0_y);
            PUTVAL_2D(slice_data_f, n_x, i0_x, i0_y, v);
        }

    //for (int i_x = 0; i_x < n_x; ++i_x)
    //    Print("%.15e\n", INDEX_2D(slice_data_f, n_x, i_x, n_y - 1));
    //exit(0);
}

////////////////////////////////////////////////////////////////////////////////
// END: IMOD utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: OpenCV utilities. {
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: OpenCV utilities: IplImage. {
////////////////////////////////////////////////////////////////////////////////

void copy_padded_array_2D_float_to_padded_array_2D_float
(
     void *array_2D0, // Raw source data, float.
     void *array_2D, // Raw destination data, float.
     int height, // Assumed to be min of array_2Ds' heights.
     int width, // Assumed to be min of array_2Ds' widths.
     int step0, // Source step in bytes.
     int step // Destination step in bytes.
)
{
    for (int i_row = 0; i_row < height; ++i_row)
    {
        float *p_data0 = (float *) (array_2D0 + i_row * step0);
        float *p_data = (float *) (array_2D + i_row * step);
        for (int i_col = 0; i_col < width; ++i_col)
        {
            *p_data++ = *p_data0++;
        }
    }
}

void copy_padded_array_2D_double_to_padded_array_2D_float
(
     void *array_2D0, // Raw source data, double.
     void *array_2D, // Raw destination data, float.
     int height, // Assumed to be min of array_2Ds' heights.
     int width, // Assumed to be min of array_2Ds' widths.
     int step0, // Source step in bytes.
     int step // Destination step in bytes.
)
{
    for (int i_row = 0; i_row < height; ++i_row)
    {
        double *p_data0 = (double *) (array_2D0 + i_row * step0);
        float *p_data = (float *) (array_2D + i_row * step);
        for (int i_col = 0; i_col < width; ++i_col)
        {
            *p_data++ = (float) *p_data0++;
        }
    }
}

void copy_padded_array_2D_float_to_padded_array_2D_double
(
     void *array_2D0, // Raw source data, float.
     void *array_2D, // Raw destination data, double.
     int height, // Assumed to be min of array_2Ds' heights.
     int width, // Assumed to be min of array_2Ds' widths.
     int step0, // Source step in bytes.
     int step // Destination step in bytes.
)
{
    for (int i_row = 0; i_row < height; ++i_row)
    {
        float *p_data0 = (float *) (array_2D0 + i_row * step0);
        double *p_data = (double *) (array_2D + i_row * step);
        for (int i_col = 0; i_col < width; ++i_col)
        {
            *p_data++ = (double) *p_data0++;
        }
    }
}

void check_IplImage_widthStep_align(IplImage *image, int n_bytes_row)
{
    if (image->widthStep != n_bytes_row)
        ABORT("Expected image->widthStep = %i, but image->widthStep = %i.\n", n_bytes_row, image->widthStep);

    /* 
     * We assume data is not aligned.  Hopefully the number of the beast will
     * cause OpenCV's head to spin off if the align value is ever used.  I
     * found the following at both
     * 
     * http://opencv.willowgarage.com/wiki/CxCore
     *
     * and
     *
     * http://www.cognotics.com/opencv/docs/1.0/ref/opencvref_cxcore.htm
     *
     * :
     * 
     * =====
     * int  align;         /\* Alignment of image rows (4 or 8).
     *                         OpenCV ignores it and uses widthStep instead *\/
     *
     * "align is ignored by OpenCV, while widthStep is used to access to subsequent image rows."
     * =====
     *
     * , but the 2008 O'Reilly book neither confirms nor denies this.
     */

    image->align = 666;
}

void check_IplImage_depth_vs_size(IplImage *image, size_t size)
{
    switch (image->depth)
    {
        case IPL_DEPTH_32F:

            if (size != sizeof(float))
                ABORT("IPL_DEPTH_32F but size == %d != sizeof(float).\n", size);

            break;

        case IPL_DEPTH_64F:

            if (size != sizeof(double))
                ABORT("IPL_DEPTH_64F but size == %d != sizeof(double).\n", size);

            break;

        default:

            ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n", IPL_DEPTH_32F, IPL_DEPTH_64F, image->depth);
    }
}

IplImage *create_IplImage_32F_from_float_data(int n_x, int n_y, float *data)
{
    IplImage *image = cvCreateImageHeader(cvSize(n_x, n_y), IPL_DEPTH_32F, 1);
    image->origin = IPL_ORIGIN_BL;

    int n_bytes_row = n_x * sizeof(float);
    cvSetData(image, data, n_bytes_row);
    
    check_IplImage_widthStep_align(image, n_bytes_row);

    return image;
}

IplImage *create_IplImage_64F_from_double_data(int n_x, int n_y, double *data)
{
    IplImage *image = cvCreateImageHeader(cvSize(n_x, n_y), IPL_DEPTH_64F, 1);
    image->origin = IPL_ORIGIN_BL;

    int n_bytes_row = n_x * sizeof(double);
    cvSetData(image, data, n_bytes_row);
    
    check_IplImage_widthStep_align(image, n_bytes_row);

    return image;
}

IplImage *create_IplImage_XXF_from_pixel_type_data(int n_x, int n_y, pixel_type *data)
{
    IplImage *image = NULL;

    switch (sizeof(pixel_type))
    {
        case sizeof(float): // WARNING: Assuming float.

            image = create_IplImage_32F_from_float_data(n_x, n_y, (float *) data);
            break;

        case sizeof(double): // WARNING: Assuming double.

            image = create_IplImage_64F_from_double_data(n_x, n_y, (double *) data);
            break;

        default:
            ABORT("Expected sizeof(pixel_type) = %i (float) or %i (double), but sizeof(pixel_type) = %i.\n", sizeof(float), sizeof(double), sizeof(pixel_type));
    }

    return image;
}

void free_imageData_cvReleaseImageHeader(IplImage *image)
{
   if (image != NULL)
   {
      free(image->imageData);
      cvReleaseImageHeader(&image);
   }
}

IplImage *create_IplImage_malloc(int n_x, int n_y, int ipl_depth)
{
    IplImage *image = NULL;

    switch (ipl_depth)
    {
        case IPL_DEPTH_32F:
        {
            int n_bytes = n_x * n_y * sizeof(float);
            float *data_copy = (float *) malloc(n_bytes); 
            if (data_copy == NULL)
                ABORT("ERROR: %s -- Cannot acquire memory for float *data_copy.\n");

            image = create_IplImage_32F_from_float_data(n_x, n_y, data_copy);

            break;
        }
        case IPL_DEPTH_64F:
        {
            int n_elems = n_x * n_y;
            int n_bytes = n_elems * sizeof(double);
            double *data_copy = (double *) malloc(n_bytes); 
            if (data_copy == NULL)
                ABORT("Cannot acquire memory for double *data_copy.\n");

            image = create_IplImage_64F_from_double_data(n_x, n_y, data_copy);

            break;
        }

        default:
            ABORT("Unrecognized IPL depth, ipl_depth = %i.\n", ipl_depth);
    }

    return image;
}

IplImage *create_IplImage_from_float_data_copy_malloc(int n_x, int n_y, float *data, int ipl_depth)
{
    IplImage *image = NULL;

    switch (ipl_depth)
    {
        case IPL_DEPTH_32F:
        {
            int n_bytes = n_x * n_y * sizeof(float);
            float *data_copy = (float *) malloc(n_bytes); 
            if (data_copy == NULL)
                ABORT("Cannot acquire memory for float *data_copy.\n");

            memcpy(data_copy, data, n_bytes);

            image = create_IplImage_32F_from_float_data(n_x, n_y, data_copy);

            break;
        }
        case IPL_DEPTH_64F:
        {
            int n_elems = n_x * n_y;
            int n_bytes = n_elems * sizeof(double);
            double *data_copy = (double *) malloc(n_bytes); 
            if (data_copy == NULL)
                ABORT("Cannot acquire memory for double *data_copy.\n");

            for (int i = 0; i < n_elems; ++i)
                data_copy[i] = (double) data[i];

            image = create_IplImage_64F_from_double_data(n_x, n_y, data_copy);

            break;
        }

        default:
            ABORT("Unrecognized IPL depth, ipl_depth = %i.\n", ipl_depth);
    }

    return image;
}

IplImage *create_IplImage_from_Islice_copy_malloc(Islice *slice, int ipl_depth)
{
    Islice *slice_float = NULL;

    if (slice->mode != SLICE_MODE_FLOAT)
        slice_float = sliceNewMode_copy(slice, SLICE_MODE_FLOAT);
    else
        slice_float = slice;

    IplImage *image = create_IplImage_from_float_data_copy_malloc(slice_float->xsize, slice_float->ysize, slice_float->data.f, ipl_depth);

    if (slice->mode != SLICE_MODE_FLOAT)
        sliceFree(slice_float);

    return image;
}

IplImage *create_IplImage_from_sliceReadMRC_malloc(MrcHeader *file_header, int i_tilt, char axis, int ipl_depth)
{
    Islice *slice = NULL;
    slice = sliceReadMRC(file_header, i_tilt, axis);
    if (slice == NULL)
        ABORT("Cannot read MRC file.\n");

    IplImage *image = create_IplImage_from_Islice_copy_malloc(slice, ipl_depth);

    sliceFree(slice);

    return image;
}

IplImage *create_IplImage_SetZero_malloc(int n_x, int n_y, int ipl_depth)
{
    IplImage *image = create_IplImage_malloc(n_x, n_y, ipl_depth);
    cvSetZero(image);
    return image;
}

IplImage *create_IplImage_from_float_data_copy(int n_x, int n_y, float *data, int ipl_depth)
{
    IplImage *image0 = create_IplImage_32F_from_float_data(n_x, n_y, data);
    IplImage *image = cvCreateImage(cvSize(n_x, n_y), ipl_depth, 1);

    switch (ipl_depth)
    {
        case IPL_DEPTH_32F:
        {
            //image = cvCloneImage(image0); // NOTE: Data might not be properly aligned for efficient access.

            CvMat tmp;
            cvGetSubRect(image0, &tmp, cvRect(0, 0, image0->width, image0->height));
            cvCopy(&tmp, image, NULL);

            break;
        }
        case IPL_DEPTH_64F:
        {

            copy_padded_array_2D_float_to_padded_array_2D_double(
                 image0->imageData, image->imageData,
                 image0->height, image0->width,
                 image0->widthStep, image->widthStep);

            break;
        }

        default:
            ABORT("Unrecognized IPL depth, ipl_depth = %i.\n", ipl_depth);
    }

    cvReleaseImageHeader(&image0);

    return image;
}

IplImage *create_IplImage_from_Islice_copy(Islice *slice, int ipl_depth)
{
    Islice *slice_float = NULL;

    if (slice->mode != SLICE_MODE_FLOAT)
        slice_float = sliceNewMode_copy(slice, SLICE_MODE_FLOAT);
    else
        slice_float = slice;

    IplImage *image = create_IplImage_from_float_data_copy(slice_float->xsize, slice_float->ysize, slice_float->data.f, ipl_depth);

    if (slice->mode != SLICE_MODE_FLOAT)
        sliceFree(slice_float);

    return image;
}

IplImage *create_IplImage_from_sliceReadMRC(MrcHeader *file_header, int i_tilt, char axis, int ipl_depth)
{
    Islice *slice = NULL;
    slice = sliceReadMRC(file_header, i_tilt, axis);
    if (slice == NULL)
        ABORT("Cannot read MRC file.\n");

    IplImage *image = create_IplImage_from_Islice_copy(slice, ipl_depth);

    sliceFree(slice);

    return image;
}

IplImage *create_IplImage_SetZero(int n_x, int n_y, int ipl_depth)
{
    IplImage *image = cvCreateImage(cvSize(n_x, n_y), ipl_depth, 1);
    image->origin = IPL_ORIGIN_BL;
    cvSetZero(image);
    return image;
}

// BEGIN: IplImage_base utilities. {

// Initializes IplImage_base with pointer to existing data.
// WARNING: No integrity checks!
// WARNING: Will not work with ROI widthSteps!
void init_IplImage_base_with_IplImage(IplImage_base *image_base, IplImage *image)
{
//    int nSize; // Unused.
//    int ID; // Unused.
    image_base->nChannels = image->nChannels; // Usually 1, but may need more at some point.
//    int alphaChannel; // Unused.
    image_base->depth = image->depth; // Pixel depth in bits.
//    char colorModel[4]; // Unused.
//    char channelSeq[4]; // Unused.
//    int dataOrder; // Unused.  Multiple channels are interleaved.
//    int origin;  // Unused.  Always lower-left-hand corner.
//    int align; // Unused.
    image_base->width = image->width; // Number of columns.
    image_base->height = image->height; // Number of rows.
//    struct _IplROI* roi; // Unused.
//    struct _IplImage* maskROI; // Unused.
//    void* imageId; // Unused.
//    struct _IplTileInfo* tileInfo; // Unused.
    image_base->imageSize = image->imageSize; // imageData size in bytes == image->height * image->widthStep if and only if image is not widthStep ROI!
    image_base->imageData = image->imageData;
    image_base->widthStep = image->widthStep; // Size of aligned image row in bytes.  (Contains the number of bytes between points in the same column and successive rows.)
//    int orderMode[4]; // Unused.
//    int orderConst[4]; // Unused.
//    char* imageDataOrigin; // Unused.
}

int convert_Ipl_depth_to_Ipl_depth_base(int ipl_depth)
{
    int ipl_depth_base = 0;

    switch (ipl_depth)
    {
        case IPL_DEPTH_32F_BASE:
            ipl_depth_base = IPL_DEPTH_32F_BASE;
            break;
        case IPL_DEPTH_64F_BASE:
            ipl_depth_base = IPL_DEPTH_64F_BASE;
            break;
        default:
            ABORT("Unrecognized Ipl depth, ipl_depth = %i.\n", ipl_depth);
    };

    return ipl_depth_base;
}

uint_address_type convert_IplImage_depth_to_byte_depth(int ipl_depth)
{
    return convert_bit_depth_to_byte_depth(convert_Ipl_depth_to_Ipl_depth_base(ipl_depth));
}

// Calculate byte offset from initial byte of image to initial byte of image_ROI.
uint_address_type calc_byte_offset_of_IplImage_ROI_widthStep(IplImage *image, IplImage *image_ROI)
{
    IplImage_base image_base;
    init_IplImage_base_with_IplImage(&image_base, image);

    IplImage_base image_ROI_base;
    init_IplImage_base_with_IplImage(&image_ROI_base, image_ROI);
    
    return calc_byte_offset_of_IplImage_base_ROI_widthStep(&image_base, &image_ROI_base);
}

// Calculate element offset from initial element of image to initial element of image_ROI.
uint_address_type calc_elem_offset_of_IplImage_ROI_widthStep(IplImage *image, IplImage *image_ROI)
{
    IplImage_base image_base;
    init_IplImage_base_with_IplImage(&image_base, image);

    IplImage_base image_ROI_base;
    init_IplImage_base_with_IplImage(&image_ROI_base, image_ROI);
    
    return calc_elem_offset_of_IplImage_base_ROI_widthStep(&image_base, &image_ROI_base);
}

// WARNING: Will not work with ROI widthSteps! (Not true...)
void copy_IplImage_to_IplImage_base(IplImage *image, IplImage_base *image_base)
{
    init_IplImage_base_with_IplImage(image_base, image);

    // BEGIN: Further data integrity checks.

    image_base->depth = convert_Ipl_depth_to_Ipl_depth_base(image->depth); // Pixel depth in bits.

    if (image_base->imageSize != image_base->height * image_base->widthStep)
        ABORT("image_base->imageSize == %i != %i == (%i == image_base->height * %i == image_base->widthStep).\n",
            image_base->imageSize, image_base->height * image_base->widthStep, image_base->height, image_base->widthStep);

    // END: Further data integrity checks.
}

// WARNING: This does not appear to be necessary.
void copy_IplImage_ROI_widthStep_to_IplImage_base_ROI_widthStep
(
    IplImage *image, // The original IplImage.
    IplImage *image_ROI, // An ROI widthStep of the original IplImage.
    IplImage_base *image_ROI_base // The output IplImage_base.
)
{
    // BEGIN: Input data integrity checks.

    if (image->nChannels != image_ROI->nChannels)
        ABORT("image->nChannels == %d != %d == image_ROI->nChannels\n", image->nChannels, image_ROI->nChannels);
    if (image->depth != image_ROI->depth)
        ABORT("image->depth == %d != %d == image_ROI->depth\n", image->depth, image_ROI->depth);
    if (image->width < image_ROI->width)
        ABORT("image->width == %d < %d == image_ROI->width\n", image->width, image_ROI->width);
    if (image->height < image_ROI->height)
        ABORT("image->height == %d < %d == image_ROI->height\n", image->height, image_ROI->height);
    if (image->imageSize < image_ROI->imageSize)
        ABORT("image->imageSize == %d < %d == image_ROI->imageSize\n", image->imageSize, image_ROI->imageSize);
    calc_byte_offset_of_IplImage_ROI_widthStep(image, image_ROI); // NOTE: The return value is not used.
    if (image->widthStep != image_ROI->widthStep)
        ABORT("image->nChannels == %d != %d == image_ROI->nChannels\n", image->nChannels, image_ROI->nChannels);

    // END: Input data integrity checks.

    init_IplImage_base_with_IplImage(image_ROI_base, image_ROI);

    // BEGIN: Further data integrity checks.
    // END: Further data integrity checks.
}

// END: IplImage_base utilities. }

// WARNING: IPL_DEPTH_32F and IPL_DEPTH_64F only!
// WARNING: SLICE_MODE_FLOAT only!
// If image extent < slice extent, excess slice (i.e., right and/or top) is filled with zeros.
// If image extent > slice extent, slice is a crop (i.e., right and/or top) of image.
void copy_IplImage_to_Islice(IplImage *image, Islice *slice)
{
    if (slice->mode != SLICE_MODE_FLOAT)
        ABORT("Expected slice->mode = %i (SLICE_MODE_FLOAT), but slice->mode = %i.\n", SLICE_MODE_FLOAT, slice->mode);

    int height = fmin(image->height, slice->ysize);
    int width = fmin(image->width, slice->xsize);

    // WARNING: slice_copy->dsize does not seem to be correct!
    //int step_slice = slice->xsize * slice->dsize * slice->csize;
    int step_slice = slice->xsize * sizeof(float);

    switch (image->depth)
    {
        case IPL_DEPTH_32F:

            copy_padded_array_2D_float_to_padded_array_2D_float(
                image->imageData, slice->data.b,
                height, width,
                image->widthStep, step_slice);
/*
            for (int i_row = 0; i_row < height; ++i_row)
            {
                float *p_image = (float *) (image->imageData + i_row * image->widthStep);
                float *p_slice = (float *) (slice->data.b + i_row * step_slice);
                for (int i_col = 0; i_col < width; ++i_col)
                {
                    *p_slice++ = *p_image++;
                }
            }
*/

            // WARNING: The code below will not work with subrectangles or aligned/padded images.
            //int n_bytes = height * width * sizeof(float);
            //memcpy(p_slice, p_image, n_bytes);

            break;

        case IPL_DEPTH_64F:

            copy_padded_array_2D_double_to_padded_array_2D_float(
                image->imageData, slice->data.b,
                height, width,
                image->widthStep, step_slice);
/*
            for (int i_row = 0; i_row < height; ++i_row)
            {
                double *p_image = (double *) (image->imageData + i_row * image->widthStep);
                float *p_slice = (float *) (slice->data.b + i_row * step_slice);
                for (int i_col = 0; i_col < width; ++i_col)
                {
                    *p_slice++ = (float) *p_image++;
                }
            }
*/

            break;

        default:

            ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n", IPL_DEPTH_32F, IPL_DEPTH_64F, image->depth);
    }
}

// If slice extent < image extent, excess image about center is filled with zeros.
// If slice extent > image extent, image is a crop about center of slice.
void copy_IplImage_to_Islice_center(IplImage *image, Islice *slice)
{
    if (image->depth != IPL_DEPTH_REAL)
        ABORT("Expected image->depth = %i (IPL_DEPTH_REAL), but image->depth = %i.\n", IPL_DEPTH_REAL, image->depth);

    if (slice->mode != SLICE_MODE_FLOAT)
        ABORT("Expected slice->mode = %i (SLICE_MODE_FLOAT), but slice->mode = %i.\n", SLICE_MODE_FLOAT, slice->mode);

    cvResetImageROI(image);

    int n0_x = image->width;
    int n0_y = image->height;

    int n_x = slice->xsize;
    int n_y = slice->ysize;

    if ((n0_x < n_x) || (n0_y < n_y))
    {
        float v_clear = 0.0;
        sliceClear(slice, &v_clear);
    }

    center_ROI_params_type params = calc_center_ROI_params(n0_x, n0_y, n_x, n_y);

    int elem_size_image = sizeof(pixel_type);
    int x0_start_offset = params.i0_x_start * elem_size_image;

    int elem_size_slice = sizeof(float);
    //int step_slice = slice->xsize * slice->dsize * slice->csize; // WARNING: slice_copy->dsize does not seem to be correct!
    int step_slice = slice->xsize * elem_size_slice;
    int x_start_offset = params.i_x_start * elem_size_slice;

    int i_row = params.i_y_start;
    for (int i0_row = params.i0_y_start; i0_row < params.i0_y_stop; ++i0_row)
    {
        pixel_type *p_image = (pixel_type *) (image->imageData + i0_row * image->widthStep + x0_start_offset);
        float *p_slice = (float *) (slice->data.b + i_row * step_slice + x_start_offset);
        for (int i0_col = params.i0_x_start; i0_col < params.i0_x_stop; ++i0_col)
        {
            *p_slice++ = (float) *p_image++;
        }
        ++i_row;
    }
}

// If slice extent < image extent, excess image (i.e., right and/or top) is filled with zeros.
// If slice extent > image extent, image is a crop (i.e., right and/or top) of slice.
void copy_Islice_to_IplImage(Islice *slice, IplImage *image)
{
    IplImage *image0 = create_IplImage_from_Islice_copy(slice, image->depth);

    cvResetImageROI(image);
    cvSetZero(image);

    int height = fmin(image0->height, image->height);
    int width = fmin(image0->width, image->width);
    cvSetImageROI(image0, cvRect(0, 0, width, height));
    cvSetImageROI(image, cvRect(0, 0, width, height));
    cvCopy(image0, image, NULL);
    cvResetImageROI(image0);
    cvResetImageROI(image);

    cvReleaseImage(&image0);
}

// If slice extent < image extent, excess image about center is filled with zeros.
// If slice extent > image extent, image is a crop about center of slice.
void copy_Islice_to_IplImage_center(Islice *slice, IplImage *image)
{
    cvResetImageROI(image);

    IplImage *image0 = create_IplImage_from_Islice_copy(slice, image->depth);
    int n0_x = image0->width;
    int n0_y = image0->height;

    int n_x = image->width;
    int n_y = image->height;

    if ((n0_x < n_x) || (n0_y < n_y))
        cvSetZero(image);

    center_ROI_params_type params = calc_center_ROI_params(n0_x, n0_y, n_x, n_y);

    cvSetImageROI(image0, cvRect(params.i0_x_start, params.i0_y_start, params.n_i0_x, params.n_i0_y));
    cvSetImageROI(image, cvRect(params.i_x_start, params.i_y_start, params.n_i0_x, params.n_i0_y)); // NOTE: The extents of this ROI are the same as image0's ROI.
    cvCopy(image0, image, NULL);
    cvResetImageROI(image0);
    cvResetImageROI(image);

    cvReleaseImage(&image0);
}

// NOTE: Uses widthStep to define a subimage (i.e., an "ROI widthStep") pointing into image->imageData.
// WARNING: Do NOT free/release imageData of IplImage * returned by this function.
IplImage *create_IplImage_ROI_widthStep_of_IplImage(IplImage *image, CvRect rect)
{
    check_IplImage_depth_vs_size(image, sizeof(pixel_type));

    if ((image->height < (rect.y + rect.height)) || (image->width < (rect.x + rect.width)))
        ABORT("image extents must be >= rect extents.\n");

    IplImage *subimage = cvCreateImageHeader(
        cvSize(rect.width, rect.height),
        image->depth, image->nChannels);

    subimage->origin = image->origin;
    subimage->widthStep = image->widthStep; // NOTE: Set widthStep of subimage to widthStep of larger image.

    subimage->imageData = image->imageData + 
        rect.y * image->widthStep +
        rect.x * image->nChannels * sizeof(pixel_type);

    return subimage;
}

void write_MRC_image_from_IplImage(IplImage *image, const char *filepath)
{
    Islice *slice = sliceCreate(image->width, image->height, SLICE_MODE_FLOAT);
    if (slice == NULL)
        ABORT("Cannot acquire memory for Islice *slice.\n");
    copy_IplImage_to_Islice(image, slice);
    write_MRC_image_from_Islice(slice, filepath);
    sliceFree(slice);
}

////////////////////////////////////////////////////////////////////////////////
// END: OpenCV utilities: IplImage. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: OpenCV utilities: CvMat. {
////////////////////////////////////////////////////////////////////////////////

CvMat *create_CvMat_32F_from_float_data(int rows, int cols, float *data)
{
    CvMat *matrix = cvCreateMatHeader(rows, cols, CV_32FC1);

    int n_bytes_row = cols * sizeof(float);
    cvSetData(matrix, data, n_bytes_row);
    
    return matrix;
}

CvMat *create_CvMat_64F_from_double_data(int rows, int cols, double *data)
{
    CvMat *matrix = cvCreateMatHeader(rows, cols, CV_64FC1);

    int n_bytes_row = cols * sizeof(double);
    cvSetData(matrix, data, n_bytes_row);
    
    return matrix;
}

/*
CvMat *create_CvMat_64F_from_double_data_copy(int rows, int cols, double *data)
{
    int n_bytes = rows * cols * sizeof(double);
    double *data_copy = (double *) malloc(n_bytes); 
    if (data_copy == NULL)
        ABORT("Cannot acquire memory for double *data_copy.\n");

    memcpy(data_copy, data, n_bytes);

    return create_CvMat_64F_from_double_data(rows, cols, data_copy);
}

CvMat *create_CvMat_from_Islice_copy(Islice *slice, int type)
{
    Islice *slice_float = NULL;

    if (slice->mode != SLICE_MODE_FLOAT)
        slice_float = sliceNewMode_copy(slice, SLICE_MODE_FLOAT);
    else
        slice_float = slice;

    CvMat *matrix = create_CvMat_64F_from_float_data_copy(slice->xsize, slice->ysize, slice->data.f, type);

    if (slice->mode != SLICE_MODE_FLOAT)
        sliceFree(slice_float);

    return matrix;
}
*/

int convert_CvMat_type_to_Ipl_depth(int type)
{
    int ipl_depth = 0;

    switch (type)
    {
        case CV_32FC1:
            ipl_depth = IPL_DEPTH_32F;
            break;
        case CV_64FC1:
            ipl_depth = IPL_DEPTH_64F;
            break;
        default:
            ABORT("Unrecognized CvMat type, type = %i.\n", type);
    }

    return ipl_depth;
}

// WARNING: Allocates an IplImage header that must be freed by the caller with cvReleaseImageHeader().
// WARNING: IplImage's imageData points to matrix->ptr.
IplImage *create_IplImage_from_CvMat(CvMat *matrix)
{
    IplImage *image = cvCreateImageHeader(cvSize(matrix->cols, matrix->rows), convert_CvMat_type_to_Ipl_depth(cvGetElemType(matrix)), 1);
    IplImage *retval = cvGetImage(matrix, image);
    if (retval != image)
        ABORT("retval at address %llu != image at address %llu\n", (uint_address_type) retval, (uint_address_type) image);
    if (matrix->data.ptr != (uchar *) image->imageData)
        ABORT("matrix->data.prt == %llu != %llu == (uchar *) image->imageData\n", (uint_address_type) matrix->data.ptr, (uint_address_type) image->imageData);
    if (matrix->step != image->widthStep)
        ABORT("matrix->step == %d != %d == image->widthStep\n", matrix->step, image->widthStep);
    return retval;
}

void copy_CvMat_to_Islice(CvMat *matrix, Islice *slice)
{
    if (slice->mode != SLICE_MODE_FLOAT)
        ABORT("Expected slice->mode = %i (SLICE_MODE_FLOAT), but slice->mode = %i.\n", SLICE_MODE_FLOAT, slice->mode);

    IplImage *image = create_IplImage_from_CvMat(matrix);
    copy_IplImage_to_Islice(image, slice);
    cvReleaseImageHeader(&image);
}

CvMat *create_CvMat_from_Islice_copy(Islice *slice, int type)
{
    IplImage *image = create_IplImage_from_Islice_copy(slice, convert_CvMat_type_to_Ipl_depth(type));
    CvMat *matrix = cvCreateMatHeader(image->width, image->height, type);
    cvGetMat(image, matrix, NULL, 0);
    cvReleaseImageHeader(&image);

    return matrix;
}

void free_data_cvReleaseMat(CvMat *matrix)
{
   if (matrix != NULL)
   {
      free(matrix->data.ptr);
      cvReleaseMat(&matrix);
   }
}

//    if ((cvGetElemType(matrix) != CV_32FC1) && (cvGetElemType(matrix) != CV_64FC1))
//        ABORT("Expected cvGetElemType(matrix) = (%i, %i) (CV_32FC1, CV_64FC1), but cvGetElemType(matrix) = %i.\n", CV_32FC1, CV_64FC1, cvGetElemType(matrix));

void print_CvMat_packed_complex(const char *name, const CvMat *matrix)
{
    if (cvGetElemType(matrix) != CVMAT_TYPE_REAL)
        ABORT("Expected cvGetElemType(matrix) = %i (CVMAT_TYPE_REAL), but cvGetElemType(matrix) = %i.\n", CVMAT_TYPE_REAL, cvGetElemType(matrix));

    Print("%s(%i, %i) = [ ...\n", name, matrix->rows, matrix->cols);
    for (int i_row = 0; i_row < matrix->rows; ++i_row) 
    {
        const pixel_type *ptr = (const pixel_type *) (matrix->data.ptr + i_row * matrix->step);
        Print("%.15e + 0.0i%s", *ptr++, (matrix->cols > 1) ? " " : "");
        for (int i_col = 1; i_col < matrix->cols - 2; i_col += 2) 
        {
            Print("%.15e + %.15ei ", *ptr, *(ptr + 1));
            ptr += 2;
        }
        if ((matrix->cols > 1) && (matrix->cols % 2 == 1))
            Print("%.15e + %.15ei", *ptr, *(ptr + 1));
        else
            Print("%.15e + 0.0i", *ptr);
        Print("\n");
    }
    Print("];\n");

/*  // NOTE: Always keep this in mind, you fool!
    for (int i_row = 0; i_row < matrix->rows; ++i_row) 
    {
        for (int i_col = 0; i_col < matrix->cols; ++i_col) 
        {
            Print("%.15e ", CV_MAT_ELEM(*matrix, pixel_type, i_row, i_col));
        }
        Print("\n");
    }
    exit(0);
*/
}

void print_CvMat(const char *name, const CvMat *matrix)
{
    if (cvGetElemType(matrix) != CVMAT_TYPE_REAL)
        ABORT("Expected cvGetElemType(matrix) = %i (CVMAT_TYPE_REAL), but cvGetElemType(matrix) = %i.\n", CVMAT_TYPE_REAL, cvGetElemType(matrix));

    Print("%s(%i, %i) = [ ...\n", name, matrix->rows, matrix->cols);
    for (int i_row = 0; i_row < matrix->rows; ++i_row) 
    {
        const pixel_type *ptr = (const pixel_type *) (matrix->data.ptr + i_row * matrix->step);
        for (int i_col = 0; i_col < matrix->cols - 1; ++i_col) 
            Print("%.15e ", *ptr++);
        Print("%.15e\n", *ptr);
    }
    Print("];\n", name);
}

void write_MRC_image_from_CvMat(CvMat *matrix, const char *filepath)
{
    Islice *slice = sliceCreate(matrix->cols, matrix->rows, SLICE_MODE_FLOAT);
    if (slice == NULL)
        ABORT("Cannot acquire memory for Islice *slice.\n");
    copy_CvMat_to_Islice(matrix, slice);
    write_MRC_image_from_Islice(slice, filepath);
    sliceFree(slice);
}

////////////////////////////////////////////////////////////////////////////////
// END: OpenCV utilities: CvMat. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// END: OpenCV utilities. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// END: Image/matrix utilities. }
////////////////////////////////////////////////////////////////////////////////
