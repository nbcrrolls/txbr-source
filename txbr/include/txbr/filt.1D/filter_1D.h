#ifndef FILTER_1D_H
#define FILTER_1D_H

#include "utilities.h"
#include "edgetaper_2D.h"
#include "filter_1D_CUDA.h"

////////////////////////////////////////////////////////////////////////////////
// BEGIN: TxBR 2.0 and TxBR 3.0 code. {
////////////////////////////////////////////////////////////////////////////////

// BEGIN: filesystem params

typedef struct filesystem_params_adt *filesystem_params_handle;

filesystem_params_handle create_filesystem_params_from_data_copy
(
    const char *filepath_in,
    const char *projection_indices,
    const char *filepath_out
);

void filesystem_params_release(filesystem_params_handle filesystem_params);

// END: filesystem params

// BEGIN: remap 2D params

typedef struct remap_2D_params_adt *remap_2D_params_handle;

remap_2D_params_handle create_remap_2D_params_from_data_copy
(
    real_type local_mag,
    int map_2D_order,
    int n_coeffs,
    int power_order_n_x1,
    int power_order_n_x2,
    real_type *power_order, 
    int n_projections,
    int map_2D_n_x1,
    int map_2D_n_x2,
    int map_2D_n_x3,
    real_type *map_2D
);

remap_2D_params_handle create_remap_2D_params_from_file
(
    const char *remap_2D_params_filepath
);

void remap_2D_params_release(remap_2D_params_handle remap_2D_params);

// END: remap 2D params

// BEGIN: rotate 2D params

typedef struct rotate_2D_params_adt *rotate_2D_params_handle;

rotate_2D_params_handle create_rotate_2D_params_from_angle_list_data_copy
(
    real_type local_mag,
    int n_angles,
    real_type *angle_list,
    real_type center[2]
);

rotate_2D_params_handle create_rotate_2D_params_from_angle_data_copy
(
    real_type local_mag,
    real_type angle,
    real_type center[2]
);

rotate_2D_params_handle create_rotate_2D_params_from_rot_matrix_2x2_data_copy
(
    real_type local_mag,
    real_type rot_matrix_2x2[2][2],
    real_type center[2]
);

rotate_2D_params_handle create_rotate_2D_params_from_rot_matrix_2x2_file
(
    const char *rotate_2D_params_filepath
);

rotate_2D_params_handle create_rotate_2D_params_from_file
(
    const char *rotate_2D_params_filepath
);

void rotate_2D_params_release(rotate_2D_params_handle rotate_2D_params);

// END: rotate 2D params

// BEGIN: filter 1D params

typedef struct filter_1D_params_adt *filter_1D_params_handle;

filter_1D_params_handle create_filter_1D_params_from_data_copy
(
    const char *type_string,
    int length,
    real_type cut_off,
    real_type roll_off
);

filter_1D_params_handle create_filter_1D_params_from_strings
(
    const char *type_string,
    const char *length_string,
    const char *cut_off_string,
    const char *roll_off_string
);

int filter_1D_params_type_is_Custom_R_Weighted(filter_1D_params_handle filter_1D_params);

filter_1D_params_handle create_filter_1D_params_from_file
(
    const char *filter_1D_params_filepath
);

void filter_1D_params_release(filter_1D_params_handle filter_1D_params);

// END: filter 1D params

// BEGIN: symmetrize 2D params

typedef struct symmetrize_2D_params_adt *symmetrize_2D_params_handle;

symmetrize_2D_params_handle create_symmetrize_2D_params_from_data_copy(int symmetrize_2D_flag);

void symmetrize_2D_params_release(symmetrize_2D_params_handle symmetrize_2D_params);

// END: symmetrize 2D params

// BEGIN: edgetaper 2D params

typedef struct edgetaper_2D_params_adt *edgetaper_2D_params_handle;

edgetaper_2D_params_handle create_edgetaper_2D_params_from_data_copy(int edgetaper_width);

void edgetaper_2D_params_release(edgetaper_2D_params_handle edgetaper_2D_params);

// END: edgetaper 2D params

void projection_series_remap_2D_filter_1D_inv_remap_2D
(
    filesystem_params_handle filesystem_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    remap_2D_params_handle remap_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
);


void projection_series_rotate_2D_filter_1D_inv_rotate_2D
(
    filesystem_params_handle filesystem_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    rotate_2D_params_handle rotate_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
); 

////////////////////////////////////////////////////////////////////////////////
// END: TxBR 2.0 and TxBR 3.0 code. }
////////////////////////////////////////////////////////////////////////////////

void projection_series_filter_1D
(
    filesystem_params_handle filesystem_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    remap_2D_params_handle remap_2D_params,
    rotate_2D_params_handle rotate_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
);

#endif // FILTER_1D_H
