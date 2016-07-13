#include "edgetaper_2D.h"

real_type edgetaper_function_1D(real_type x)
{
    real_type y = (tanh(x) + 1.0) / 2.0;
    //real_type y = 0.0;
    return y;
}

real_type *create_edgetaper_filter_1D(int length)
{
    static const char *function_name = "create_edgetaper_filter_1D()";

    real_type *filter_1D = (real_type *) malloc(length * sizeof(real_type)); 
    if (filter_1D == NULL)
        Abort("ERROR: %s -- Memory request failed.\n", function_name);

    // Sample ((tanh(x) + 1) / 2) at length equally spaced points on the domain interval [-range, range].
    real_type max = 5;

    if (max <= 0)
        Abort("ERROR: %s -- max = %.15e <= 0.\n", function_name, max);

    real_type x_delta = (2.0 * max) / (real_type) length; 
    for (int i = 0; i < length; ++i)
        filter_1D[i] = edgetaper_function_1D(-max + (x_delta * i));

    return filter_1D;
}

void edgetaper_2D(IplImage *image, int edge_width)
{
    static const char *function_name = "edgetaper_2D()";

    check_IplImage_depth_vs_size(image, sizeof(real_type));

    if (edge_width < 1)
        Abort("ERROR %s -- edge_width = %i must be >= 1.\n", function_name, edge_width);

    if (((image->height / 2) < edge_width) || ((image->width / 2) < edge_width))
        Abort("ERROR %s -- edge_width = %i must be < 1/2 extent, but image->(height, width) = (%i, %i).\n", function_name, edge_width, image->height, image->width);

    real_type *et_filter_1D = create_edgetaper_filter_1D(edge_width);

    for (int i_row = 0; i_row < image->height; ++i_row)
    {
        real_type *p_image = (real_type *) (image->imageData + i_row * image->widthStep);
        for (int i_col = 0; i_col < edge_width; ++i_col)
        {
            p_image[i_col] *= et_filter_1D[i_col];
        }

        int i_col_start = image->width - edge_width;
        int i_et_filter_1D = edge_width;
        for (int i_col = i_col_start; i_col < image->width; ++i_col)
        {
            p_image[i_col] *= et_filter_1D[--i_et_filter_1D];
        }
    }

    for (int i_col = 0; i_col < image->width; ++i_col)
    {
        int colstep = i_col * sizeof(real_type);

        for (int i_row = 0; i_row < edge_width; ++i_row)
        {
            real_type *p_image = (real_type *) (image->imageData + i_row * image->widthStep + colstep);
            *p_image *= et_filter_1D[i_row];
        }

        int i_row_start = image->height - edge_width;
        int i_et_filter_1D = edge_width;
        for (int i_row = i_row_start; i_row < image->height; ++i_row)
        {
            real_type *p_image = (real_type *) (image->imageData + i_row * image->widthStep + colstep);
            *p_image *= et_filter_1D[--i_et_filter_1D];
        }
    }

    free(et_filter_1D);
}
