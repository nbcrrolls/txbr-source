#include "filter_1D.h"

// BEGIN: MEX includes
// NOTE: Leave this in the .c file.

#ifdef MEX
#include "mex.h"
#endif

// END: MEX includes

#define PROTOTYPE_COMPLIANT_INDEXING 1
#define PROTOTYPE_COMPLIANT_ROTATION 1

//#define FLOOR_INNER_LOOP(x_) ((int_type) (x_)) // Will _not_ work with negative input.
#define FLOOR_INNER_LOOP(x_) ((int_type) floor(x_))
//#define FLOOR_INNER_LOOP(x_) (cvFloor(x_))

//#define INT_FROM_REAL_IMAGE_EXTENT(x_) ((int_type) floor(x_))
//#define INT_FROM_REAL_IMAGE_EXTENT(x_) ((int_type) round(x_))
#define INT_FROM_REAL_IMAGE_EXTENT(x_) ((int_type) ceil(x_))
//#define INT_FROM_REAL_IMAGE_EXTENT(x_) (((int_type) floor(x_)) + 1) // NOTE: Least strictly greatest integer.

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Parameters. {
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: TxBR 2.0 and TxBR 3.0 code. {
////////////////////////////////////////////////////////////////////////////////

// BEGIN: File system parameters.

typedef struct filesystem_params_adt // The public type is abstract.
{
    char *filepath_in;
    int n_projections_indexed;
    int *projection_indices;
    char *filepath_out;
}
filesystem_params_type; // The private type is concrete.

filesystem_params_handle create_filesystem_params()
{
    filesystem_params_handle filesystem_params = (filesystem_params_handle) malloc(sizeof(filesystem_params_type));
    if (filesystem_params == NULL)
        ABORT("Cannot acquire memory for filesystem_params_handle filesystem_params.\n");

    // BEGIN: Default values.
    filesystem_params->filepath_in = NULL;
    filesystem_params->n_projections_indexed = 0;
    filesystem_params->projection_indices = NULL;
    filesystem_params->filepath_out= NULL;
    // END: Default values.

    return filesystem_params;
}

filesystem_params_handle create_filesystem_params_from_data_copy
(
    const char *filepath_in,
    const char *projection_indices_string,
    const char *filepath_out
)
{
    filesystem_params_handle filesystem_params = create_filesystem_params();

    if (strcmp(filepath_in, filepath_out) == 0)
        ABORT("\"%s\" == filepath_in == filepath_out == \"%s\".\n", filepath_in, filepath_out);

    filesystem_params->filepath_in = calloc(strlen(filepath_in) + 1, sizeof(char));
    if (filesystem_params->filepath_in == NULL)
        ABORT("Cannot acquire memory for filesystem_params->filepath_in.\n");
    strcpy(filesystem_params->filepath_in, filepath_in);

    // Get the input file's sizes.
    int n_x, n_y, n_projections_total;
    get_MRC_sizes(filesystem_params->filepath_in, &n_x, &n_y, &n_projections_total);

    create_indexing_array_from_indexing_string(
        n_projections_total,
        (projection_indices_string != NULL) ? projection_indices_string : "-",
        &(filesystem_params->n_projections_indexed),
        &(filesystem_params->projection_indices));

    filesystem_params->filepath_out = calloc(strlen(filepath_out) + 1, sizeof(char));
    if (filesystem_params->filepath_out == NULL)
        ABORT("Cannot acquire memory for filesystem_params->filepath_out.\n");
    strcpy(filesystem_params->filepath_out, filepath_out);

    return filesystem_params;
}

void filesystem_params_release(filesystem_params_handle filesystem_params)
{
    if (filesystem_params->filepath_in != NULL)
        free(filesystem_params->filepath_in);

    if (filesystem_params->projection_indices != NULL)
        free(filesystem_params->projection_indices);

    if (filesystem_params->filepath_in != NULL)
        free(filesystem_params->filepath_out);

    free(filesystem_params);
}

// END: File system parameters.

// BEGIN: 2D symmetrization parameters.

typedef struct symmetrize_2D_params_adt
{
    int symmetrize_2D_flag;
}
symmetrize_2D_params_type;

symmetrize_2D_params_handle create_symmetrize_2D_params()
{
    symmetrize_2D_params_handle symmetrize_2D_params = (symmetrize_2D_params_handle) malloc(sizeof(symmetrize_2D_params_type));
    if (symmetrize_2D_params == NULL)
       ABORT("Cannot aquire memory for symmetrize_2D_params_handle symmetrize_2D_params.\n");

    // BEGIN: Default values.
    symmetrize_2D_params->symmetrize_2D_flag = 0;
    // END: Default values.

    return symmetrize_2D_params;
}

symmetrize_2D_params_handle create_symmetrize_2D_params_from_data_copy(int symmetrize_2D_flag)
{
    symmetrize_2D_params_handle symmetrize_2D_params = create_symmetrize_2D_params();

    symmetrize_2D_params->symmetrize_2D_flag = symmetrize_2D_flag;

    return symmetrize_2D_params;
}

void symmetrize_2D_params_release(symmetrize_2D_params_handle symmetrize_2D_params)
{
    free(symmetrize_2D_params);
}

// END: 2D symmetrization parameters.

// BEGIN: 2D remap parameters.

typedef struct remap_2D_params_adt
{
    real_type local_mag;
    int map_2D_order;
    int n_coeffs;
    int power_order_n_x1;
    int power_order_n_x2;
    real_type *power_order;
    int map_2D_n_x1;
    int map_2D_n_x2;
    int map_2D_n_x3;
    real_type *map_2D;
}
remap_2D_params_type;

remap_2D_params_handle create_remap_2D_params()
{
    remap_2D_params_handle remap_2D_params = (remap_2D_params_handle) malloc(sizeof(remap_2D_params_type));
    if (remap_2D_params == NULL)
        ABORT("Cannot acquire memory for remap_2D_params_handle remap_2D_params.\n");

    // BEGIN: Default values.
    remap_2D_params->local_mag = 1.0;
    remap_2D_params->map_2D_order = 0;
    remap_2D_params->n_coeffs = 0;
    remap_2D_params->power_order_n_x1 = 0;
    remap_2D_params->power_order_n_x2 = 0;
    remap_2D_params->power_order = NULL;
    remap_2D_params->map_2D_n_x1 = 0;
    remap_2D_params->map_2D_n_x2 = 0;
    remap_2D_params->map_2D_n_x3 = 0;
    remap_2D_params->map_2D = NULL;
    // END: Default values.

    return remap_2D_params;
}

// NOTE: Assumes data is in row-major order.

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
)
{
    remap_2D_params_handle remap_2D_params = create_remap_2D_params();

    if (local_mag <= 0.0)
        ABORT("local_mag == %.15e <= 0.0.", local_mag);
    remap_2D_params->local_mag = local_mag;

    if (map_2D_order < 0)
        ABORT("map_2D_order == %i < 0.", map_2D_order);
    remap_2D_params->map_2D_order = map_2D_order;

    int n_coeffs_expected = ((map_2D_order + 1) * (map_2D_order + 2)) / 2; // even * odd = even
    if (n_coeffs != n_coeffs_expected)
        ABORT("n_coeffs == %i != %i == n_coeffs_expected == ((map_2D_order + 1) * (map_2D_order + 2)) / 2\n", n_coeffs, n_coeffs_expected);
    remap_2D_params->n_coeffs = n_coeffs;

    if (power_order_n_x1 != 2)
        ABORT("power_order_n_x1 == %i != 2\n", power_order_n_x1);
    remap_2D_params->power_order_n_x1 = power_order_n_x1;
    if (power_order_n_x2 != n_coeffs)
        ABORT("power_order_n_x2 == %i != %i == n_coeffs\n", power_order_n_x2, n_coeffs);
    remap_2D_params->power_order_n_x2 = power_order_n_x2;
    int n_power_order_bytes = power_order_n_x1 * power_order_n_x2 * sizeof(real_type);
    remap_2D_params->power_order = (real_type *) calloc(n_power_order_bytes, 1);
    if (remap_2D_params->power_order == NULL)
        ABORT("Cannot acquire memory for real_type *power_order.\n");
    memcpy(remap_2D_params->power_order, power_order, n_power_order_bytes);

    if (map_2D_n_x1 != 4)
        ABORT("map_2D_n_x1 == %i != 4\n", map_2D_n_x1);
    remap_2D_params->map_2D_n_x1 = map_2D_n_x1;

    if (map_2D_order == 0)
    {
        if (map_2D_n_x2 != 3)
            ABORT("map_2D_order == 0, but map_2D_n_x2 == %i != 3\n", map_2D_n_x2);
    }
    else if (map_2D_n_x2 != n_coeffs)
    {
        ABORT("map_2D_n_x2 == %i != %i == n_coeffs\n", map_2D_n_x2, n_coeffs);
    }

    remap_2D_params->map_2D_n_x2 = map_2D_n_x2;

    if (map_2D_n_x3 != n_projections)
        ABORT("map_2D_n_x3 == %i != %i == n_projections\n", map_2D_n_x3, n_projections);
    remap_2D_params->map_2D_n_x3 = map_2D_n_x3;

    int n_map_2D_bytes = map_2D_n_x1 * map_2D_n_x2 * map_2D_n_x3 * sizeof(real_type);
    remap_2D_params->map_2D = (real_type *) calloc(n_map_2D_bytes, 1);
    if (remap_2D_params->map_2D == NULL)
        ABORT("Cannot acquire memory for real_type *remap_2D_params->map_2D.\n");
/*
    for (int i_x1 = 0; i_x1 < map_2D_n_x1; ++i_x1)
        for (int i_x2 = 0; i_x2 < map_2D_n_x2; ++i_x2)
            for (int i_x3 = 0; i_x3 < map_2D_n_x3; ++i_x3)
                MATRIX_INDEX_3D(remap_2D_params->map_2D, map_2D_n_x1, map_2D_n_x2, i_x1, i_x2, i_x3) = 
                    MATRIX_INDEX_3D(map_2D, map_2D_n_x1, map_2D_n_x2, i_x1, i_x2, i_x3);
*/
    memcpy(remap_2D_params->map_2D, map_2D, n_map_2D_bytes);

// PRIMARY DEBUG:     Print("remap_2D_params->local_mag = %.15e\n", remap_2D_params->local_mag); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("remap_2D_params->map_2D_order = %i\n", remap_2D_params->map_2D_order); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("remap_2D_params->n_coeffs = %i\n", remap_2D_params->n_coeffs); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("remap_2D_params->power_order_n_x(1, 2) = (%i, %i)\n", remap_2D_params->power_order_n_x1, remap_2D_params->power_order_n_x2); // PRIMARY DEBUG.
// PRIMARY DEBUG:     print_matrix_2D("remap_2D_params->power_order", remap_2D_params->power_order, remap_2D_params->power_order_n_x1, remap_2D_params->power_order_n_x2); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("remap_2D_params->map_2D_n_x(1, 2, 3) = (%i, %i, %i)\n", remap_2D_params->map_2D_n_x1, remap_2D_params->map_2D_n_x2, remap_2D_params->map_2D_n_x3); // PRIMARY DEBUG.
// PRIMARY DEBUG:     print_matrix_3D("remap_2D_params->map_2D", remap_2D_params->map_2D, remap_2D_params->map_2D_n_x1, remap_2D_params->map_2D_n_x2, remap_2D_params->map_2D_n_x3); // PRIMARY DEBUG.

//ABORT("ABORT.\n");

    return remap_2D_params;
}

/*
// NOTE: Assumes data is in MATLAB order.

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
)
{
    remap_2D_params_handle remap_2D_params = create_remap_2D_params();

    if (local_mag <= 0.0)
        ABORT("local_mag == %.15e <= 0.0.", local_mag);
    remap_2D_params->local_mag = local_mag;

    if (map_2D_order < 0)
        ABORT("map_2D_order == %i < 0.", map_2D_order);
    remap_2D_params->map_2D_order = map_2D_order;

    int n_coeffs_expected = ((map_2D_order + 1) * (map_2D_order + 2)) / 2; // even * odd = even
    if (n_coeffs != n_coeffs_expected)
        ABORT("n_coeffs == %i != %i == n_coeffs_expected == ((map_2D_order + 1) * (map_2D_order + 2)) / 2\n", n_coeffs, n_coeffs_expected);
    remap_2D_params->n_coeffs = n_coeffs;

    if (power_order_n_x1 != 2)
        ABORT("power_order_n_x1 == %i != 2\n", power_order_n_x1);
    remap_2D_params->power_order_n_x1 = power_order_n_x1;
    if (power_order_n_x2 != n_coeffs)
        ABORT("power_order_n_x2 == %i != %i == n_coeffs\n", power_order_n_x2, n_coeffs);
    remap_2D_params->power_order_n_x2 = power_order_n_x2;
    int n_power_order_bytes = power_order_n_x1 * power_order_n_x2 * sizeof(real_type);
    remap_2D_params->power_order = (real_type *) calloc(n_power_order_bytes, 1);
    if (remap_2D_params->power_order == NULL)
        ABORT("Cannot acquire memory for real_type *power_order.\n");
    memcpy(remap_2D_params->power_order, power_order, n_power_order_bytes);

    WARN("Resizing map_2D matrix so that the projection index is dimension x3.\n");
    if (map_2D_n_x1 != n_projections)
        ABORT("map_2D_n_x1 == %i != %i == n_projections\n", map_2D_n_x1, n_projections);
    remap_2D_params->map_2D_n_x3 = map_2D_n_x1;
    if (map_2D_n_x2 != 4)
        ABORT("map_2D_n_x2 == %i != 4\n", map_2D_n_x2);
    remap_2D_params->map_2D_n_x1 = map_2D_n_x2;

    if (map_2D_order == 0)
    {
        if (map_2D_n_x3 != 3)
            ABORT("map_2D_order == 0, but map_2D_n_x3 == %i != 3\n", map_2D_n_x3);
    }
    else if (map_2D_n_x3 != n_coeffs)
    {
        ABORT("map_2D_n_x3 == %i != %i == n_coeffs\n", map_2D_n_x3, n_coeffs);
    }

    remap_2D_params->map_2D_n_x2 = map_2D_n_x3;
    int n_map_2D_bytes = map_2D_n_x1 * map_2D_n_x2 * map_2D_n_x3 * sizeof(real_type);
    remap_2D_params->map_2D = (real_type *) calloc(n_map_2D_bytes, 1);
    if (remap_2D_params->map_2D == NULL)
        ABORT("Cannot acquire memory for real_type *remap_2D_params->map_2D.\n");
    for (int i_x1 = 0; i_x1 < map_2D_n_x1; ++i_x1)
        for (int i_x2 = 0; i_x2 < map_2D_n_x2; ++i_x2)
            for (int i_x3 = 0; i_x3 < map_2D_n_x3; ++i_x3)
                MATRIX_INDEX_3D(remap_2D_params->map_2D, map_2D_n_x2, map_2D_n_x3, i_x2, i_x3, i_x1) = 
                    MATRIX_INDEX_3D(map_2D, map_2D_n_x1, map_2D_n_x2, i_x1, i_x2, i_x3);

    Print("remap_2D_params->local_mag = %.15e\n", remap_2D_params->local_mag);
    Print("remap_2D_params->map_2D_order = %i\n", remap_2D_params->map_2D_order);
    Print("remap_2D_params->n_coeffs = %i\n", remap_2D_params->n_coeffs);
    Print("remap_2D_params->power_order_n_x(1, 2) = (%i, %i)\n", remap_2D_params->power_order_n_x1, remap_2D_params->power_order_n_x2);
    print_matrix_2D("remap_2D_params->power_order", remap_2D_params->power_order, remap_2D_params->power_order_n_x1, remap_2D_params->power_order_n_x2);
    Print("remap_2D_params->map_2D_n_x(1, 2, 3) = (%i, %i, %i)\n", remap_2D_params->map_2D_n_x1, remap_2D_params->map_2D_n_x2, remap_2D_params->map_2D_n_x3);
    print_matrix_3D("remap_2D_params->map_2D", remap_2D_params->map_2D, remap_2D_params->map_2D_n_x1, remap_2D_params->map_2D_n_x2, remap_2D_params->map_2D_n_x3);

    return remap_2D_params;
}
*/

remap_2D_params_handle create_remap_2D_params_from_file
(
    const char *remap_2D_params_filepath
)
{
    FILE *remap_2D_params_file = fopen(remap_2D_params_filepath, "rb");
    if (remap_2D_params_file == NULL)
        ABORT("Cannot open remap_2D_params_filepath = \"%s\".\n", remap_2D_params_filepath);

    size_t n_elems_read_expected = 0;
    size_t n_elems_read = 0;

    real_type local_mag = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&local_mag, sizeof(real_type), 1, remap_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read real_type local_mag, n_elems_read == %u != %u == n_elems_read_expected, remap_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, remap_2D_params_filepath);

    int map_2D_order = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&map_2D_order, sizeof(int), 1, remap_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int map_2D_order, n_elems_read == %u != %u == n_elems_read_expected, remap_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, remap_2D_params_filepath);

    int n_coeffs = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&n_coeffs, sizeof(int), 1, remap_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int n_coeffs, n_elems_read == %u != %u == n_elems_read_expected, remap_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, remap_2D_params_filepath);

    int power_order_n_x1 = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&power_order_n_x1, sizeof(int), 1, remap_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int power_order_n_x1, n_elems_read == %u != %u == n_elems_read_expected, remap_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, remap_2D_params_filepath);

    int power_order_n_x2 = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&power_order_n_x2, sizeof(int), 1, remap_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int power_order_n_x2, n_elems_read == %u != %u == n_elems_read_expected, remap_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, remap_2D_params_filepath);

    int n_power_order_elems = power_order_n_x1 * power_order_n_x2;
    int n_power_order_bytes = n_power_order_elems * sizeof(real_type);
    n_elems_read_expected = n_power_order_elems;
    real_type *power_order = (real_type *) malloc(n_power_order_bytes); 
    if (power_order == NULL)
        ABORT("Cannot acquire memory for real_type *power_order.\n");

    n_elems_read = fread(power_order, sizeof(real_type), n_power_order_elems, remap_2D_params_file);
    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read real_type *power_order, n_elems_read == %u != %u == n_elems_read_expected, remap_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, remap_2D_params_filepath);

    int n_projections = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&n_projections, sizeof(int), 1, remap_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int n_projections, n_elems_read == %u != %u == n_elems_read_expected, remap_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, remap_2D_params_filepath);

    int map_2D_n_x1 = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&map_2D_n_x1, sizeof(int), 1, remap_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int map_2D_n_x1, n_elems_read == %u != %u == n_elems_read_expected, remap_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, remap_2D_params_filepath);

    int map_2D_n_x2 = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&map_2D_n_x2, sizeof(int), 1, remap_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int map_2D_n_x2, n_elems_read == %u != %u == n_elems_read_expected, remap_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, remap_2D_params_filepath);

    int map_2D_n_x3 = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&map_2D_n_x3, sizeof(int), 1, remap_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int map_2D_n_x3, n_elems_read == %u != %u == n_elems_read_expected, remap_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, remap_2D_params_filepath);

    int n_map_2D_elems = map_2D_n_x1 * map_2D_n_x2 * map_2D_n_x3;
    int n_map_2D_bytes = n_map_2D_elems * sizeof(real_type);
    n_elems_read_expected = n_map_2D_elems;
    real_type *map_2D = (real_type *) malloc(n_map_2D_bytes);
    if (map_2D == NULL)
        ABORT("Cannot acquire memory for real_type *map_2D.\n");

    n_elems_read = fread(map_2D, sizeof(real_type), n_map_2D_elems, remap_2D_params_file);
    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read real_type *map_2D, n_elems_read == %u != %u == n_elems_read_expected, remap_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, remap_2D_params_filepath);

    if (fclose(remap_2D_params_file) != 0)
        ABORT("Cannot close remap_2D_params_filepath = \"%s\".\n", remap_2D_params_filepath);

    remap_2D_params_handle remap_2D_params =
        create_remap_2D_params_from_data_copy
        (
            local_mag,
            map_2D_order,
            n_coeffs,
            power_order_n_x1,
            power_order_n_x2,
            power_order, 
            n_projections,
            map_2D_n_x1,
            map_2D_n_x2,
            map_2D_n_x3,
            map_2D
        );

    free(power_order);
    free(map_2D);

    return remap_2D_params;
}

void remap_2D_params_release(remap_2D_params_handle remap_2D_params)
{
    if (remap_2D_params->power_order != NULL)
        free(remap_2D_params->power_order);

    if (remap_2D_params->map_2D != NULL)
        free(remap_2D_params->map_2D);

    free(remap_2D_params);
}

// END: 2D remap parameters.

// BEGIN: 2D rotation parameters.

typedef struct rotate_2D_params_adt
{
    real_type local_mag;
    int n_angles;
    real_type *angle_list;
    real_type rot_matrix_2x2[2][2]; // NOTE: Temporary?  Will eventually want an array of rotation matrices?
    real_type inv_rot_matrix_2x2[2][2]; // NOTE: Temporary?  Will eventually want an array of inv. rotation matrices?
    real_type center[2]; // NOTE: Temporary?  Will eventually want an array of rotation centers?
}
rotate_2D_params_type;

rotate_2D_params_handle create_rotate_2D_params()
{
    rotate_2D_params_handle rotate_2D_params = (rotate_2D_params_handle) malloc(sizeof(rotate_2D_params_type));
    if (rotate_2D_params == NULL)
        ABORT("Cannot acquire memory for rotate_2D_params_handle rotate_2D_params.\n");

    // BEGIN: Default values.
    rotate_2D_params->local_mag = 1.0;
    rotate_2D_params->n_angles = 0;
    rotate_2D_params->angle_list = NULL;

    rotate_2D_params->rot_matrix_2x2[0][0] = 1.0; rotate_2D_params->rot_matrix_2x2[0][1] = 0.0;
    rotate_2D_params->rot_matrix_2x2[1][0] = 0.0; rotate_2D_params->rot_matrix_2x2[1][1] = 1.0;

    rotate_2D_params->inv_rot_matrix_2x2[0][0] = 1.0; rotate_2D_params->inv_rot_matrix_2x2[0][1] = 0.0;
    rotate_2D_params->inv_rot_matrix_2x2[1][0] = 0.0; rotate_2D_params->inv_rot_matrix_2x2[1][1] = 1.0;

    rotate_2D_params->center[0] = 0.0;
    rotate_2D_params->center[0] = 0.0;
    // END: Default values.

    return rotate_2D_params;
}

void print_rotate_2D_params(rotate_2D_params_handle rotate_2D_params)
{
    Print("rotate_2D_params->local_mag = %.15e\n", rotate_2D_params->local_mag);
    Print("rotate_2D_params->n_angles = %i\n", rotate_2D_params->n_angles);
    if (rotate_2D_params->n_angles > 0)
        print_matrix_1D("rotate_2D_params->angle_list", rotate_2D_params->angle_list, rotate_2D_params->n_angles);
    print_matrix_2x2("rotate_2D_params-rot_matrix_2x2", rotate_2D_params->rot_matrix_2x2);
    print_matrix_2x2("rotate_2D_params-inv_rot_matrix_2x2", rotate_2D_params->inv_rot_matrix_2x2);
    Print("rotate_2D_params->center[0-1] = (%.15e, %.15e)\n", rotate_2D_params->center[0], rotate_2D_params->center[1]);
}

rotate_2D_params_handle create_rotate_2D_params_from_angle_list_data_copy
(
    real_type local_mag,
    int n_angles,
    real_type *angle_list,
    real_type center[2]
)
{
    rotate_2D_params_handle rotate_2D_params = create_rotate_2D_params();

    if (local_mag <= 0.0)
        ABORT("local_mag == %.15e <= 0.0.", local_mag);
    rotate_2D_params->local_mag = local_mag;

    if (n_angles < 0)
        ABORT("n_angles == %i < 0.\n", n_angles);

    rotate_2D_params->n_angles = n_angles;
    int n_bytes_angle_list = n_angles * sizeof(real_type);
    rotate_2D_params->angle_list = malloc(n_bytes_angle_list);
    if (rotate_2D_params->angle_list == NULL)
        ABORT("Can't acquire memory for rotate_2D_params->angle_list.\n");
    memcpy(rotate_2D_params->angle_list, angle_list, n_bytes_angle_list);

    for (int i_angle = 0; i_angle < n_angles; ++i_angle)
    {
        // Assume degrees; convert to radians.
        real_type angle = angle_list[i_angle] * (M_PI / 180.0);
        rotate_2D_params->angle_list[i_angle] = angle;
    }

    // Initial values.
    calc_rot_matrix_2x2(angle_list[0], rotate_2D_params->rot_matrix_2x2);
    calc_rot_matrix_2x2(-angle_list[0], rotate_2D_params->inv_rot_matrix_2x2);

    print_rotate_2D_params(rotate_2D_params);

    return rotate_2D_params;
}

rotate_2D_params_handle create_rotate_2D_params_from_angle_data_copy
(
    real_type local_mag,
    real_type angle,
    real_type center[2]
)
{
    return create_rotate_2D_params_from_angle_list_data_copy(local_mag, 1, &angle, center);
}

rotate_2D_params_handle create_rotate_2D_params_from_rot_matrix_2x2_data_copy
(
    real_type local_mag,
    real_type rot_matrix_2x2[2][2],
    real_type center[2]
)
{
    rotate_2D_params_handle rotate_2D_params = create_rotate_2D_params();

    if (local_mag <= 0.0)
        ABORT("local_mag == %.15e <= 0.0.", local_mag);
    rotate_2D_params->local_mag = local_mag;

    memcpy(rotate_2D_params->rot_matrix_2x2, rot_matrix_2x2, sizeof(real_type) * 4); 
/*
    rotate_2D_params->rot_matrix_2x2[0][0] = rot_matrix_2x2[0][0]; rotate_2D_params->rot_matrix_2x2[0][1] = rot_matrix_2x2[0][1];
    rotate_2D_params->rot_matrix_2x2[1][0] = rot_matrix_2x2[1][0]; rotate_2D_params->rot_matrix_2x2[1][1] = rot_matrix_2x2[1][1];
*/

    inv_rot_matrix_2x2(rotate_2D_params->rot_matrix_2x2, rotate_2D_params->inv_rot_matrix_2x2);

    // Recovery of eqv. angle.
    real_type cos_val = rotate_2D_params->rot_matrix_2x2[0][0];
    real_type sin_val = rotate_2D_params->rot_matrix_2x2[1][0];

    real_type angle_atan2 = atan2(sin_val, cos_val);
    if (angle_atan2 < 0.0)
        angle_atan2 = (2 * M_PI) + angle_atan2;

    //ABORT("angle_atan2 = %f\n", angle_atan2);

    rotate_2D_params->n_angles = 1;
    rotate_2D_params->angle_list = malloc(sizeof(real_type));
    if (rotate_2D_params->angle_list == NULL)
        ABORT("Can't acquire memory for rotate_2D_params->angle_list.\n");
    rotate_2D_params->angle_list[0] = angle_atan2;

    rotate_2D_params->center[0] = center[0];
    rotate_2D_params->center[1] = center[1];

    print_rotate_2D_params(rotate_2D_params);

    return rotate_2D_params;
}

rotate_2D_params_handle create_rotate_2D_params_from_rot_matrix_2x2_file
(
    const char *rotate_2D_params_filepath
)
{
    FILE *rotate_2D_params_file = fopen(rotate_2D_params_filepath, "rb");
    if (rotate_2D_params_file == NULL)
        ABORT("Cannot open rotate_2D_params_filepath = \"%s\".\n", rotate_2D_params_filepath);

    size_t n_elems_read_expected = 0;
    size_t n_elems_read = 0;

    real_type local_mag = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&local_mag, sizeof(real_type), 1, rotate_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read real_type local_mag, n_elems_read == %u != %u == n_elems_read_expected, rotate_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, rotate_2D_params_filepath);

    int rot_matrix_2x2_n_x1 = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&rot_matrix_2x2_n_x1, sizeof(int), 1, rotate_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int rot_matrix_2x2_n_x1, n_elems_read == %u != %u == n_elems_read_expected, rotate_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, rotate_2D_params_filepath);

    int rot_matrix_2x2_n_x2 = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&rot_matrix_2x2_n_x2, sizeof(int), 1, rotate_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int rot_matrix_2x2_n_x2, n_elems_read == %u != %u == n_elems_read_expected, rotate_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, rotate_2D_params_filepath);

    if ((rot_matrix_2x2_n_x1 != 2) || (rot_matrix_2x2_n_x2 != 2))
        ABORT("(rot_matrix_2x2_n_x1 == %d != 2) || (rot_matrix_2x2_n_x2 == %d != 2), rotate_2D_params_filepath = \"%s\".\n", rot_matrix_2x2_n_x1, rot_matrix_2x2_n_x2, rotate_2D_params_filepath);

    int n_rot_matrix_2x2_elems = rot_matrix_2x2_n_x1 * rot_matrix_2x2_n_x2;
//    int n_rot_matrix_2x2_bytes = n_rot_matrix_2x2_elems * sizeof(real_type);
    n_elems_read_expected = n_rot_matrix_2x2_elems;
//    real_type *rot_matrix_2x2 = (real_type *) malloc(n_rot_matrix_2x2_bytes); 
    real_type rot_matrix_2x2[2][2]; 
//    if (rot_matrix_2x2 == NULL)
//        ABORT("Cannot acquire memory for real_type *rot_matrix_2x2.\n");

    n_elems_read = fread(rot_matrix_2x2, sizeof(real_type), n_rot_matrix_2x2_elems, rotate_2D_params_file);
    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read real_type *rot_matrix_2x2, n_elems_read == %u != %u == n_elems_read_expected, rotate_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, rotate_2D_params_filepath);

    real_type center[2] = {0.0, 0.0};
    n_elems_read_expected = 2;
    n_elems_read = fread(center, sizeof(real_type), 2, rotate_2D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read real_type center[2], n_elems_read == %u != %u == n_elems_read_expected, rotate_2D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, rotate_2D_params_filepath);

    rotate_2D_params_handle rotate_2D_params =
        create_rotate_2D_params_from_rot_matrix_2x2_data_copy
        (
            local_mag,
            rot_matrix_2x2,
            center
        );

//    free(rot_matrix_2x2);

    return rotate_2D_params;
}

void rotate_2D_params_release(rotate_2D_params_handle rotate_2D_params)
{
    if (rotate_2D_params != NULL)
        free(rotate_2D_params->angle_list);

    free(rotate_2D_params);
}

// END: 2D rotation parameters.

// BEGIN: 1D filtering parameters.

enum filter_1D_type_enum
{
    SHEPPLOGAN = 0,
    RAMLAK = 1,
    CUSTOMRWEIGHTED = 2
};

typedef struct filter_1D_params_adt
{
    enum filter_1D_type_enum type;
    int length_rs;
    int length_fs;
    real_type cut_off;
    real_type roll_off;
}
filter_1D_params_type;

enum filter_1D_type_enum get_filter_1D_type_from_string(const char *type_string)
{
    int type = -1;

    if (
        (strcmp(type_string, "SheppLogan") == 0) ||
        (strcmp(type_string, "Shepp-Logan") == 0)
       )
    {
        type = SHEPPLOGAN;
    }
    else if (
             (strcmp(type_string, "RamLak") == 0) ||
             (strcmp(type_string, "Ram-Lak") == 0)
            )
    {
        type = RAMLAK;
    }
    else if (
             (strcmp(type_string, "CustomRWeighted") == 0) ||
             (strcmp(type_string, "Custom R-Weighted") == 0)
            )
    {
        type = CUSTOMRWEIGHTED;
    }
    else
    {
        ABORT("Unrecognized filter type.\n");
    }

    return type;
}

filter_1D_params_handle create_filter_1D_params()
{
    filter_1D_params_handle filter_params = (filter_1D_params_handle) malloc(sizeof(filter_1D_params_type));
    if (filter_params == NULL)
       ABORT("Cannot aquire memory for filter_1D_params_handle filter_params.\n");

    // BEGIN: Default values.
    filter_params->type = SHEPPLOGAN;
    filter_params->length_rs = 251;
    filter_params->length_fs = 0;
    filter_params->cut_off = 0.0;
    filter_params->roll_off = 0.0;
    // END: Default values.

    return filter_params;
}

filter_1D_params_handle create_filter_1D_params_from_data_copy
(
    const char *type_string,
    int length,
    real_type cut_off,
    real_type roll_off
)
{
    filter_1D_params_handle filter_params = create_filter_1D_params();

    filter_params->type = get_filter_1D_type_from_string(type_string);

//    if (length < 3)
//        ABORT("Filter length cannot be less than 3: length = %i.\n", length);

    // NOTE: Only odd-length filters allowed.
    if (length % 2 == 0)
        ABORT("Only odd-length filters allowed: length = %i.\n", length);

    //if (length % 2 == 0)
    //    ++length;

    filter_params->length_rs = length;
    filter_params->length_fs = 0;

    if (filter_params->type == CUSTOMRWEIGHTED)
    {
        filter_params->cut_off = cut_off;
        filter_params->roll_off = roll_off;
    }
    else
    {
        filter_params->cut_off = 0.0;
        filter_params->roll_off = 0.0;
    }

    return filter_params;
}

// WARNING: Unsafe!  Have not written proper checks!
filter_1D_params_handle create_filter_1D_params_from_strings
(
    const char *type_string,
    const char *length_string,
    const char *cut_off_string,
    const char *roll_off_string
)
{
    int length = atoi(length_string);

    real_type cut_off = 0.0;
    real_type roll_off = 0.0;

    if (cut_off_string != NULL)
        cut_off = atof(cut_off_string);

    if (roll_off_string != NULL)
        roll_off = atof(roll_off_string);

    return create_filter_1D_params_from_data_copy(type_string, length, cut_off, roll_off);
}

int filter_1D_params_type_is_Custom_R_Weighted(filter_1D_params_handle filter_1D_params)
{
    if (filter_1D_params == NULL)
        ABORT("filter_1D_params == NULL.\n");

    if (filter_1D_params->type == CUSTOMRWEIGHTED)
        return 1;
    else
        return 0;
}

filter_1D_params_handle create_filter_1D_params_from_file
(
    const char *filter_1D_params_filepath
)
{
    FILE *filter_1D_params_file = fopen(filter_1D_params_filepath, "rb");
    if (filter_1D_params_file == NULL)
        ABORT("Cannot open filter_1D_params_filepath = \"%s\".\n", filter_1D_params_filepath);

    size_t n_elems_read_expected = 0;
    size_t n_elems_read = 0;

    int type_string_length = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&type_string_length, sizeof(int), 1, filter_1D_params_file);
    ++type_string_length; // NOTE: String in file is not NULL-terminated.

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int type_string_length, n_elems_read == %u != %u == n_elems_read_expected, filter_1D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, filter_1D_params_filepath);

    char *type_string = calloc(type_string_length, sizeof(char));
    if (type_string == NULL)
        ABORT("Cannot acquire memory for char *type_string.\n");
    if (fgets(type_string, type_string_length, filter_1D_params_file) == NULL)
        ABORT("Failed to fgets() type_string, filter_1D_params_filepath = \"%s\".\n", filter_1D_params_filepath);

    int length = 0;
    n_elems_read_expected = 1;
    n_elems_read = fread(&length, sizeof(int), 1, filter_1D_params_file);

    if (n_elems_read != n_elems_read_expected)
        ABORT("Failed to read int length, n_elems_read == %u != %u == n_elems_read_expected, filter_1D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, filter_1D_params_filepath);

    real_type cut_off = 0;
    real_type roll_off = 0;

    if (get_filter_1D_type_from_string(type_string) == CUSTOMRWEIGHTED)
    {
        n_elems_read_expected = 1;
        n_elems_read = fread(&cut_off, sizeof(real_type), 1, filter_1D_params_file);

        if (n_elems_read != n_elems_read_expected)
            ABORT("Failed to read real_type cut_off, n_elems_read == %u != %u == n_elems_read_expected, filter_1D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, filter_1D_params_filepath);

        n_elems_read_expected = 1;
        n_elems_read = fread(&roll_off, sizeof(real_type), 1, filter_1D_params_file);

        if (n_elems_read != n_elems_read_expected)
            ABORT("Failed to read real_type roll_off, n_elems_read == %u != %u == n_elems_read_expected, filter_1D_params_filepath = \"%s\".\n", n_elems_read, n_elems_read_expected, filter_1D_params_filepath);
    }

    if (fclose(filter_1D_params_file) != 0)
        ABORT("Cannot close filter_1D_params_filepath = \"%s\".\n", filter_1D_params_filepath);

    //Print("type_string_length = %d\n", type_string_length);
    //Print("type_string = \"%s\"\n", type_string);
    //Print("length = %d\n", length);
    //Print("cut_off = %f\n", cut_off);
    //Print("roll_off = %f\n", roll_off);

    //ABORT("ABORT!\n");

    filter_1D_params_handle filter_1D_params =
        create_filter_1D_params_from_data_copy
        (
            type_string,
            length,
            cut_off,
            roll_off
        );

    free(type_string);

    return filter_1D_params;
}

void filter_1D_params_release(filter_1D_params_handle filter_1D_params)
{
    free(filter_1D_params);
}

// END: 1D filtering parameters.

// BEGIN: 2D edgetapering parameters.

typedef struct edgetaper_2D_params_adt
{
    int width;
}
edgetaper_2D_params_type;

edgetaper_2D_params_handle create_edgetaper_2D_params()
{
    edgetaper_2D_params_handle edgetaper_2D_params = (edgetaper_2D_params_handle) malloc(sizeof(edgetaper_2D_params_type));
    if (edgetaper_2D_params == NULL)
       ABORT("Cannot aquire memory for edgetaper_2D_params_handle edgetaper_2D_params.\n");

    // BEGIN: Default values.
    edgetaper_2D_params->width = 0;
    // END: Default values.

    return edgetaper_2D_params;
}

edgetaper_2D_params_handle create_edgetaper_2D_params_from_data_copy(int width)
{
    edgetaper_2D_params_handle edgetaper_2D_params = create_edgetaper_2D_params();

    if (width < 0)
        ABORT("width == %i < 0.\n", width);

    edgetaper_2D_params->width = width;

    return edgetaper_2D_params;
}

void edgetaper_2D_params_release(edgetaper_2D_params_handle edgetaper_2D_params)
{
    free(edgetaper_2D_params);
}

// END: 2D edgetapering parameters.

////////////////////////////////////////////////////////////////////////////////
// END: TxBR 2.0 and TxBR 3.0 code. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// END: Parameters. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Misc. {
////////////////////////////////////////////////////////////////////////////////

// NOTE: Centers rect_2DR, rotates centered rect_2DR about origin, and then recenters.
void rect_2DR_rotate_2D(real_type rot_matrix_2x2[2][2], rect_2DR_type *rect_2DR)
{
    Print("rect_2DR = (%.15e, %.15e, %.15e, %.15e)\n", rect_2DR->x, rect_2DR->y, rect_2DR->width, rect_2DR->height);
 
    real_type p0[4][2] = // Image coordinates of corners.
        {
         {(real_type) rect_2DR->x, (real_type) rect_2DR->y},
         {(real_type) (rect_2DR->x + rect_2DR->width), (real_type) rect_2DR->y},
         {(real_type) rect_2DR->x, (real_type) (rect_2DR->y + rect_2DR->height)},
         {(real_type) (rect_2DR->x + rect_2DR->width), (real_type) (rect_2DR->y + rect_2DR->height)}
        };

    Print("p0[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n",
        p0[0][0], p0[0][1], p0[1][0], p0[1][1], p0[2][0], p0[2][1], p0[3][0], p0[3][1]);

    // BEGIN: Unpack rotate 2D parameters and data.

    real_type rot_matrix00 = rot_matrix_2x2[0][0];
    real_type rot_matrix01 = rot_matrix_2x2[0][1];
    real_type rot_matrix10 = rot_matrix_2x2[1][0];
    real_type rot_matrix11 = rot_matrix_2x2[1][1];

    // END: Unpack rotate 2D parameters and data.

/*
    // BEGIN: Calculate rotated center0. {

    real_type center0[2] = {rect_2DR->x + (rect_2DR->width / 2.0), rect_2DR->y + (rect_2DR->height / 2.0)};

    real_type center0_rotated[2] = {
                                    (rot_matrix00 * center0[0]) + (rot_matrix01 * center0[1]),
                                    (rot_matrix10 * center0[0]) + (rot_matrix11 * center0[1]),
                                   };

// PRIMARY DEBUG:     Print("center0_rotated = (%.15e, %.15e)\n", center0_rotated[0], center0_rotated[1]); // PRIMARY DEBUG.

    // END: Calculate rotated center0. }
*/

    real_type p[4][2]; // Image coordinates of remapped corners.

    for (int i_p = 0; i_p < 4; ++i_p)
    {
        //Print("p0[%i][0-1] = (%.15e, %.15e)\n", i_p0, p0[i_p][0], p0[i_p][1]);

        //p[i_p][0] -= center0[0];
        //p[i_p][1] -= center0[1];

        p[i_p][0] = (rot_matrix00 * p0[i_p][0]) + (rot_matrix01 * p0[i_p][1]);
        p[i_p][1] = (rot_matrix10 * p0[i_p][0]) + (rot_matrix11 * p0[i_p][1]);

        //Print("p[%i][0-1] = (%.15e, %.15e)\n", i_p, p[i_p][0], p[i_p][1]);
    }

    //Print("p[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n",
    //    p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]);

    // BEGIN: Build a bounding rectangle.

    real_type x_min = p[0][0];
    x_min = fmin(x_min, p[1][0]);
    x_min = fmin(x_min, p[2][0]);
    x_min = fmin(x_min, p[3][0]);

    real_type x_max = p[0][0];
    x_max = fmax(x_max, p[1][0]);
    x_max = fmax(x_max, p[2][0]);
    x_max = fmax(x_max, p[3][0]);

    real_type y_min = p[0][1];
    y_min = fmin(y_min, p[1][1]);
    y_min = fmin(y_min, p[2][1]);
    y_min = fmin(y_min, p[3][1]);

    real_type y_max = p[0][1];
    y_max = fmax(y_max, p[1][1]);
    y_max = fmax(y_max, p[2][1]);
    y_max = fmax(y_max, p[3][1]);

// PRIMARY DEBUG:     Print("(xy)_min = (%.15e, %.15e)\n", x_min, y_min); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("(xy)_max = (%.15e, %.15e)\n", x_max, y_max); // PRIMARY DEBUG.

    // END: Build a bounding rectangle.

    //rect_2DR->x = 0.0; // NOTE: The rotated rect_2DR is recentered.
    //rect_2DR->y = 0.0; // NOTE: The rotated rect_2DR is recentered.
    rect_2DR->x = x_min;
    rect_2DR->y = y_min;
    rect_2DR->width = x_max - x_min;
    rect_2DR->height = y_max - y_min;

/*
    real_type min_xy_trans[2] = { x_min - x0_min, y_min - y0_min };

    p[0][0] += min_xy_trans[0]; 
    p[0][1] += min_xy_trans[1];

    p[1][0] += min_xy_trans[0];
    p[1][1] += min_xy_trans[1];

    p[2][0] += min_xy_trans[0];
    p[2][1] += min_xy_trans[1];

    p[3][0] += min_xy_trans[0];
    p[3][1] += min_xy_trans[1];
*/

    Print("p[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n",
        p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]);

    Print("rect_2DR = (%.15e, %.15e, %.15e, %.15e)\n", rect_2DR->x, rect_2DR->y, rect_2DR->width, rect_2DR->height);

    //ABORT("Abort!\n");
}

void calc_rotated_support_corners
(
    real_type rot_matrix_2x2[2][2],
    real_type support[4],
    real_type corners[4][2] // Image coordinates of rotated corners.
)
{
// NOTE: Does NOT use PROTOTYPE_COMPLIANT_INDEXING.
//#if PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

// PRIMARY DEBUG:     Print("support[0-3] = (%.15e, %.15e, %.15e, %.15e)\n", support[0], support[1], support[2], support[3]); // PRIMARY DEBUG.

    //real_type center0_x = support[0];
    //real_type center0_y = support[1];

//#if PROTOTYPE_COMPLIANT_INDEXING
//    real_type half0_x = support[2] / 2.0; // WARNING: Expands support by 0.5 on each side.
//    real_type half0_y = support[3] / 2.0; // WARNING: Expands support by 0.5 on each side.
//#else
//    real_type half0_x = (support[2] + 1.0) / 2.0; // WARNING: Expands support by 0.5 on each side.
//    real_type half0_y = (support[3] + 1.0) / 2.0; // WARNING: Expands support by 0.5 on each side.
    real_type half0_x = support[2] / 2.0;
    real_type half0_y = support[3] / 2.0;
//#endif

    // p[1-4][2], four corners of the support in projection coordinates
    real_type p0[4][2]; 

    p0[0][0] = -half0_x; p0[0][1] = -half0_y; 
    p0[1][0] = half0_x; p0[1][1] = -half0_y;
    p0[2][0] = -half0_x; p0[2][1] = half0_y;
    p0[3][0] = half0_x; p0[3][1] = half0_y;

// PRIMARY DEBUG:     Print("p0[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n", // PRIMARY DEBUG.
// PRIMARY DEBUG:         p0[0][0], p0[0][1], p0[1][0], p0[1][1], p0[2][0], p0[2][1], p0[3][0], p0[3][1]); // PRIMARY DEBUG.

    // p[1-4][2], four rotated corners of the support in projection coordinates
    real_type p[4][2]; 

/*  
    // BEGIN: Complies with MATLAB prototype.
    p[0][0] = (rot_matrix_2x2[0][0] * p0[0][0]) + (rot_matrix_2x2[0][1] * p0[0][1]) + center0_x;
    p[0][1] = (rot_matrix_2x2[1][0] * p0[0][0]) + (rot_matrix_2x2[1][1] * p0[0][1]) + center0_y;

    p[1][0] = (rot_matrix_2x2[0][0] * p0[1][0]) + (rot_matrix_2x2[0][1] * p0[1][1]) + center0_x;
    p[1][1] = (rot_matrix_2x2[1][0] * p0[1][0]) + (rot_matrix_2x2[1][1] * p0[1][1]) + center0_y;

    p[2][0] = (rot_matrix_2x2[0][0] * p0[2][0]) + (rot_matrix_2x2[0][1] * p0[2][1]) + center0_x;
    p[2][1] = (rot_matrix_2x2[1][0] * p0[2][0]) + (rot_matrix_2x2[1][1] * p0[2][1]) + center0_y;

    p[3][0] = (rot_matrix_2x2[0][0] * p0[3][0]) + (rot_matrix_2x2[0][1] * p0[3][1]) + center0_x;
    p[3][1] = (rot_matrix_2x2[1][0] * p0[3][0]) + (rot_matrix_2x2[1][1] * p0[3][1]) + center0_y;
    // END: Complies with MATLAB prototype.
*/

    p[0][0] = (rot_matrix_2x2[0][0] * p0[0][0]) + (rot_matrix_2x2[0][1] * p0[0][1]);
    p[0][1] = (rot_matrix_2x2[1][0] * p0[0][0]) + (rot_matrix_2x2[1][1] * p0[0][1]);

    p[1][0] = (rot_matrix_2x2[0][0] * p0[1][0]) + (rot_matrix_2x2[0][1] * p0[1][1]);
    p[1][1] = (rot_matrix_2x2[1][0] * p0[1][0]) + (rot_matrix_2x2[1][1] * p0[1][1]);

    p[2][0] = (rot_matrix_2x2[0][0] * p0[2][0]) + (rot_matrix_2x2[0][1] * p0[2][1]);
    p[2][1] = (rot_matrix_2x2[1][0] * p0[2][0]) + (rot_matrix_2x2[1][1] * p0[2][1]);

    p[3][0] = (rot_matrix_2x2[0][0] * p0[3][0]) + (rot_matrix_2x2[0][1] * p0[3][1]);
    p[3][1] = (rot_matrix_2x2[1][0] * p0[3][0]) + (rot_matrix_2x2[1][1] * p0[3][1]);

// PRIMARY DEBUG:     Print("p[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n", // PRIMARY DEBUG.
// PRIMARY DEBUG:         p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]); // PRIMARY DEBUG.

    // BEGIN: Build a bounding rectangle.

    real_type x_min = p[0][0];
    x_min = fmin(x_min, p[1][0]);
    x_min = fmin(x_min, p[2][0]);
    x_min = fmin(x_min, p[3][0]);

    real_type x_max = p[0][0];
    x_max = fmax(x_max, p[1][0]);
    x_max = fmax(x_max, p[2][0]);
    x_max = fmax(x_max, p[3][0]);

    real_type y_min = p[0][1];
    y_min = fmin(y_min, p[1][1]);
    y_min = fmin(y_min, p[2][1]);
    y_min = fmin(y_min, p[3][1]);

    real_type y_max = p[0][1];
    y_max = fmax(y_max, p[1][1]);
    y_max = fmax(y_max, p[2][1]);
    y_max = fmax(y_max, p[3][1]);

    // END: Build a bounding rectangle.

//#if PROTOTYPE_COMPLIANT_INDEXING
//    support[0] = (x_max - x_min + 1.0) / 2.0;
//    support[1] = (y_max - y_min + 1.0) / 2.0;
//    support[2] = x_max - x_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
//    support[3] = y_max - y_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
//#else
    support[0] = (x_max - x_min) / 2.0;
    support[1] = (y_max - y_min) / 2.0;
    //support[0] = (x_max - x_min - 1.0) / 2.0; // NOTE: This doesn't look right, but it removes the discrepancy for the IPLIMAGE_PUSH versions.
    //support[1] = (y_max - y_min - 1.0) / 2.0; // NOTE: This doesn't look right, but it removes the discrepancy for the IPLIMAGE_PUSH versions.
    support[2] = x_max - x_min;
    support[3] = y_max - y_min;
    //support[2] = x_max - x_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
    //support[3] = y_max - y_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
//#endif

    real_type half_x = support[2] / 2.0;
    real_type half_y = support[3] / 2.0;

    real_type center_trans[2] = { half_x - half0_x, half_y - half0_y };

    p[0][0] += center_trans[0]; 
    p[0][1] += center_trans[1];

    p[1][0] += center_trans[0];
    p[1][1] += center_trans[1];

    p[2][0] += center_trans[0];
    p[2][1] += center_trans[1];

    p[3][0] += center_trans[0];
    p[3][1] += center_trans[1];

    memcpy(corners, p, sizeof(real_type) * 4 * 2);

//    Print("p[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n",
//        p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]);

// PRIMARY DEBUG:     Print("support[0-3] = (%.15e, %.15e, %.15e, %.15e)\n", support[0], support[1], support[2], support[3]); // PRIMARY DEBUG.
}

void calc_rotated_support(real_type rot_matrix_2x2[2][2], real_type support[4])
{
    real_type corners[4][2];
    calc_rotated_support_corners(rot_matrix_2x2, support, corners);
}

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Extend image values. {
////////////////////////////////////////////////////////////////////////////////

void extend_IplImage_rows(IplImage *image)
{
    if ((image->depth != IPL_DEPTH_32F) && (image->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, image->depth);

// PRIMARY DEBUG:     write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.mrc"); // PRIMARY DEBUG.

    pixel_type *image_data = (pixel_type *) image->imageData;
    int n_x = image->width;
    int n_y = image->height;

    pixel_type n_x_ws_pixel = ((pixel_type) image->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("n_x_ws_pixel - n_x_ws != 0.0.");

    for (int i_y = 0; i_y < n_y; ++i_y)
    {
        int i_x_left_stop = n_x;
        for (int i_x = 0; i_x < n_x; ++i_x)
        {
            pixel_type v = INDEX_2D(image_data, n_x_ws, i_x, i_y); // Inefficient for serial access on same row.
            if (v != 0.0)
            {
                i_x_left_stop = i_x;
                for (int i_x_left = 0; i_x_left < i_x_left_stop; ++i_x_left)
                    PUTVAL_2D(image_data, n_x_ws, i_x_left, i_y, v); // Inefficient for serial access on same row.
                break;
            }
        }
        for (int i_x = n_x - 1; i_x >= i_x_left_stop; --i_x)
        {
            pixel_type v = INDEX_2D(image_data, n_x_ws, i_x, i_y); // Inefficient for serial access on same row.
            if (v != 0.0)
            {
                for (int i_x_right = i_x; i_x_right < n_x; ++i_x_right)
                    PUTVAL_2D(image_data, n_x_ws, i_x_right, i_y, v); // Inefficient for serial access on same row.
                break;
            }
        }
    }

// PRIMARY DEBUG:     write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.extend_IplImage_rows.mrc"); // PRIMARY DEBUG.
}

// NOTE: Does not attempt to extend row left or right until discrete first derivative is nonzero for (?:n_samples|two) subsequent samples.
// NOTE: Extends a constant number of pixels into the image to pick value of extended row.
void extend_IplImage_rows_analyze(IplImage *image)
{
    if ((image->depth != IPL_DEPTH_32F) && (image->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, image->depth);

// PRIMARY DEBUG:     write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.mrc"); // PRIMARY DEBUG.

    pixel_type *image_data = (pixel_type *) image->imageData;
    int n_x = image->width;
    int n_y = image->height;

//    const int n_samples = 10;
//    const int n_pixels_in = 10;
    const int n_samples = 5;
    const int n_pixels_in = 5;
    if (n_x < n_samples + (2 * n_pixels_in))
        return;

    pixel_type n_x_ws_pixel = ((pixel_type) image->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("n_x_ws_pixel - n_x_ws != 0.0.");

    // BEGIN: n_samples samples.

    pixel_type v_samples[n_samples];
    pixel_type dv_samples[n_samples - 1];

    for (int i_y = 0; i_y < n_y; ++i_y)
    {
        int i_x_left_stop = n_x;
        for (int i_v_sample = 0; i_v_sample < n_samples - 1; ++i_v_sample)
            v_samples[i_v_sample] = INDEX_2D(image_data, n_x_ws, i_v_sample, i_y); // Inefficient for serial access on same row.
        for (int i_dv_sample = 0; i_dv_sample < n_samples - 2; ++i_dv_sample)
            dv_samples[i_dv_sample] = v_samples[i_dv_sample + 1] - v_samples[i_dv_sample];
        dv_samples[n_samples - 2] = 0.0;
        for (int i_x = n_samples - 1; i_x < n_x; ++i_x)
        {
            v_samples[n_samples - 1] = INDEX_2D(image_data, n_x_ws, i_x, i_y); // Inefficient for serial access on same row.
            dv_samples[n_samples - 2] = v_samples[n_samples - 1] - v_samples[n_samples - 2];
            int all_dv_samples_neq_zero = 1;
            for (int i_dv_sample = 0; i_dv_sample < n_samples - 1; ++i_dv_sample)
            {
                if (dv_samples[i_dv_sample] == 0.0)
                {
                    all_dv_samples_neq_zero = 0;
                    break;
                }
            }
            if (all_dv_samples_neq_zero)
            {
                i_x_left_stop = i_x + n_pixels_in;
                pixel_type v_avg = 0.0;
                for (int i_x_left = i_x + 1; i_x_left <= i_x_left_stop; ++i_x_left)
                    v_avg += INDEX_2D(image_data, n_x_ws, i_x_left, i_y); // Inefficient for serial access on same row.
                v_avg /= (pixel_type) n_pixels_in;
                i_x_left_stop = (i_x + 1) - n_samples + (int) (n_pixels_in / 2);
                for (int i_x_left = 0; i_x_left < i_x_left_stop; ++i_x_left)
                    PUTVAL_2D(image_data, n_x_ws, i_x_left, i_y, v_avg); // Inefficient for serial access on same row.
                i_x_left_stop = i_x + n_pixels_in;
                break;
            }
            for (int i_v_sample = 0; i_v_sample < n_samples - 1; ++i_v_sample)
                v_samples[i_v_sample] = v_samples[i_v_sample + 1];
            for (int i_dv_sample = 0; i_dv_sample < n_samples - 2; ++i_dv_sample)
                dv_samples[i_dv_sample] = dv_samples[i_dv_sample + 1];
        }
        for (int i_v_sample = 0; i_v_sample < n_samples - 1; ++i_v_sample)
            v_samples[i_v_sample] = INDEX_2D(image_data, n_x_ws, (n_x - 1) - i_v_sample, i_y); // Inefficient for serial access on same row.
        for (int i_dv_sample = 0; i_dv_sample < n_samples - 2; ++i_dv_sample)
            dv_samples[i_dv_sample] = v_samples[i_dv_sample + 1] - v_samples[i_dv_sample];
        dv_samples[n_samples - 2] = 0.0;
        for (int i_x = (n_x - 1) - (n_samples - 1); i_x >= i_x_left_stop; --i_x)
        {
            v_samples[n_samples - 1] = INDEX_2D(image_data, n_x_ws, i_x, i_y); // Inefficient for serial access on same row.
            dv_samples[n_samples - 2] = v_samples[n_samples - 1] - v_samples[n_samples - 2];
            int all_dv_samples_neq_zero = 1;
            for (int i_dv_sample = 0; i_dv_sample < n_samples - 1; ++i_dv_sample)
            {
                if (dv_samples[i_dv_sample] == 0.0)
                {
                    all_dv_samples_neq_zero = 0;
                    break;
                }
            }
            if (all_dv_samples_neq_zero)
            {
                int i_x_right_start = i_x - n_pixels_in;
                pixel_type v_avg = 0.0;
                for (int i_x_right = i_x_right_start; i_x_right < i_x; ++i_x_right)
                    v_avg += INDEX_2D(image_data, n_x_ws, i_x_right, i_y); // Inefficient for serial access on same row.
                v_avg /= (pixel_type) n_pixels_in;
                i_x_right_start = i_x + n_samples - (int) (n_pixels_in / 2);
                for (int i_x_right = i_x_right_start; i_x_right < n_x; ++i_x_right)
                    PUTVAL_2D(image_data, n_x_ws, i_x_right, i_y, v_avg); // Inefficient for serial access on same row.
                break;
            }
            for (int i_v_sample = 0; i_v_sample < n_samples - 1; ++i_v_sample)
                v_samples[i_v_sample] = v_samples[i_v_sample + 1];
            for (int i_dv_sample = 0; i_dv_sample < n_samples - 2; ++i_dv_sample)
                dv_samples[i_dv_sample] = dv_samples[i_dv_sample + 1];
        }
    }

    // END: n_samples samples.

/*
    // BEGIN: Three samples hardcoded.

    pixel_type v_0 = 0.0;
    pixel_type v_1 = 0.0;
    pixel_type image_dv_0_1 = 0.0;
    pixel_type image_dv_1_2 = 0.0;

    for (int i_y = 0; i_y < n_y; ++i_y)
    {
        int i_x_left_stop = n_x;
        v_0 = INDEX_2D(image_data, n_x_ws, 0, i_y); // Inefficient for serial access on same row.
        v_1 = INDEX_2D(image_data, n_x_ws, 1, i_y); // Inefficient for serial access on same row.
        image_dv_0_1 = v_1 - v_0;
        image_dv_1_2 = 0.0;
        for (int i_x = 2; i_x < n_x; ++i_x)
        {
            pixel_type v_2 = INDEX_2D(image_data, n_x_ws, i_x, i_y); // Inefficient for serial access on same row.
            image_dv_1_2 = v_2 - v_1;
            if ((image_dv_0_1 != 0.0) && (image_dv_1_2 != 0.0))
            {
                i_x_left_stop = i_x + n_pixels_in;
                pixel_type v = INDEX_2D(image_data, n_x_ws, i_x_left_stop, i_y); // Inefficient for serial access on same row.
                for (int i_x_left = 0; i_x_left < i_x_left_stop; ++i_x_left)
                    PUTVAL_2D(image_data, n_x_ws, i_x_left, i_y, v); // Inefficient for serial access on same row.
                break;
            }
            v_0 = v_1;
            v_1 = v_2;
            image_dv_0_1 = image_dv_1_2;
        }
        v_0 = INDEX_2D(image_data, n_x_ws, n_x - 1, i_y); // Inefficient for serial access on same row.
        v_1 = INDEX_2D(image_data, n_x_ws, n_x - 2, i_y); // Inefficient for serial access on same row.
        image_dv_0_1 = v_1 - v_0;
        image_dv_1_2 = 0.0;
        for (int i_x = n_x - 3; i_x >= i_x_left_stop; --i_x)
        {
            pixel_type v_2 = INDEX_2D(image_data, n_x_ws, i_x, i_y); // Inefficient for serial access on same row.
            image_dv_1_2 = v_2 - v_1;
            if ((image_dv_0_1 != 0.0) && (image_dv_1_2 != 0.0))
            {
                int i_x_right_start = i_x - n_pixels_in;
                pixel_type v = INDEX_2D(image_data, n_x_ws, i_x_right_start, i_y); // Inefficient for serial access on same row.
                for (int i_x_right = i_x_right_start + 1; i_x_right < n_x; ++i_x_right)
                    PUTVAL_2D(image_data, n_x_ws, i_x_right, i_y, v); // Inefficient for serial access on same row.
                break;
            }
            v_0 = v_1;
            v_1 = v_2;
            image_dv_0_1 = image_dv_1_2;
        }
    }

    // END: Three samples hardcoded.
*/

// PRIMARY DEBUG:     write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.extend_IplImage_rows_analyze.mrc"); // PRIMARY DEBUG.
    //ABORT("ABORT!\n");
}

void extend_IplImage_cols(IplImage *image)
{
    if ((image->depth != IPL_DEPTH_32F) && (image->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, image->depth);

// PRIMARY DEBUG:     write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.mrc"); // PRIMARY DEBUG.

    pixel_type *image_data = (pixel_type *) image->imageData;
    int n_x = image->width;
    int n_y = image->height;

    pixel_type n_x_ws_pixel = ((pixel_type) image->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("n_x_ws_pixel - n_x_ws != 0.0.");

    for (int i_x = 0; i_x < n_x; ++i_x)
    {
        int i_y_bottom_stop = n_y;
        for (int i_y = 0; i_y < n_y; ++i_y)
        {
            pixel_type v = INDEX_2D(image_data, n_x_ws, i_x, i_y); // Inefficient for serial access on same row.
            if (v != 0.0)
            {
                i_y_bottom_stop = i_y;
                for (int i_y_bottom = 0; i_y_bottom < i_y_bottom_stop; ++i_y_bottom)
                    PUTVAL_2D(image_data, n_x_ws, i_x, i_y_bottom, v); // Inefficient for serial access on same row.
                break;
            }
        }
        for (int i_y = n_y - 1; i_y >= i_y_bottom_stop; --i_y)
        {
            pixel_type v = INDEX_2D(image_data, n_x_ws, i_x, i_y); // Inefficient for serial access on same row.
            if (v != 0.0)
            {
                for (int i_y_top = i_y; i_y_top < n_y; ++i_y_top)
                    PUTVAL_2D(image_data, n_x_ws, i_x, i_y_top, v); // Inefficient for serial access on same row.
                break;
            }
        }
    }

// PRIMARY DEBUG:     write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.extend_IplImage_cols.mrc"); // PRIMARY DEBUG.
}

void extend_IplImage_values(IplImage *image)
{
    extend_IplImage_rows(image);
    extend_IplImage_cols(image);

// PRIMARY DEBUG:     write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.extend_IplImage_values.mrc"); // PRIMARY DEBUG.
    //ABORT("Abort!\n");
}

void symmetrize_IplImage_values_CvRect(IplImage *image, CvRect rect)
{
    if ((image->depth != IPL_DEPTH_32F) && (image->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, image->depth);

    //write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.mrc");

    //Print("rect = (%d, %d, %d, %d)\n", rect.x, rect.y, rect.width, rect.height);

    if ((rect.x < 0) || (rect.y < 0) || (rect.width <= 0) || (rect.height <= 0))
        ABORT("Invalid rect = (%i, %i, %i, %i).\n", rect.x, rect.y, rect.width, rect.height);

    int n_x = image->width;
    int n_y = image->height;

    CvMat matrix_src;
    CvMat matrix_dst;

/*
    // BEGIN: Extend rows left and right.
    cvGetSubRect(image, &matrix_src, cvRect(rect.x, rect.y, 1, rect.height));
    cvGetSubRect(image, &matrix_dst, cvRect(0, rect.y, rect.x, rect.height));
    cvRepeat(&matrix_src, &matrix_dst);
    cvGetSubRect(image, &matrix_src, cvRect(rect.x + rect.width - 1, rect.y, 1, rect.height));
    cvGetSubRect(image, &matrix_dst, cvRect(rect.x + rect.width, rect.y, n_x - (rect.x + rect.width), rect.height));
    cvRepeat(&matrix_src, &matrix_dst);
    // END: Extend rows left and right.

    // BEGIN: Extend cols down and up.
    cvGetSubRect(image, &matrix_src, cvRect(rect.x, rect.y, rect.width, 1));
    cvGetSubRect(image, &matrix_dst, cvRect(rect.x, 0, rect.width, rect.y));
    cvRepeat(&matrix_src, &matrix_dst);
    cvGetSubRect(image, &matrix_src, cvRect(rect.x, rect.y + rect.height - 1, rect.width, 1));
    cvGetSubRect(image, &matrix_dst, cvRect(rect.x, rect.y + rect.height, rect.width, n_y - (rect.y + rect.height)));
    cvRepeat(&matrix_src, &matrix_dst);
    // END: Extend cols down and up.
*/

    int sym_rect_width;
    int sym_rect_height;
    CvRect rect_src;
    CvRect rect_dst;

    // BEGIN: Top center.
    // NOTE: A flip about x.
    sym_rect_width = rect.width;
    sym_rect_height = (n_y - (rect.y + rect.height));
    if (sym_rect_height < 0)
        ABORT("sym_rect_height == %i.\n", sym_rect_height);
    if (sym_rect_height > 0)
    {
        rect_src = cvRect(rect.x, (rect.y + rect.height) - (n_y - (rect.y + rect.height)), sym_rect_width, sym_rect_height);
        rect_dst = cvRect(rect.x, rect.y + rect.height, sym_rect_width, sym_rect_height);
        cvGetSubRect(image, &matrix_src, rect_src);
        cvGetSubRect(image, &matrix_dst, rect_dst);
        cvFlip(&matrix_src, &matrix_dst, 0);
//        cvSet(&matrix_dst, cvRealScalar(666.0), NULL); // DEBUG.
    }
    // END: Top center.

    // BEGIN: Bottom center.
    // NOTE: A flip about x.
    sym_rect_width = rect.width;
    sym_rect_height = rect.y;
    if (sym_rect_height < 0)
        ABORT("sym_rect_height == %i.\n", sym_rect_height);
    if (sym_rect_height > 0)
    {
        rect_src = cvRect(rect.x, rect.y, sym_rect_width, sym_rect_height);
        rect_dst = cvRect(rect.x, 0, sym_rect_width, sym_rect_height);
        cvGetSubRect(image, &matrix_src, rect_src);
        cvGetSubRect(image, &matrix_dst, rect_dst);
        cvFlip(&matrix_src, &matrix_dst, 0);
//        cvSet(&matrix_dst, cvRealScalar(666.0), NULL); // DEBUG.
    }
    // END: Bottom center.

    // BEGIN: Left top, center, and bottom.
    // NOTE: A flip about y.
    sym_rect_width = rect.x;
    if (sym_rect_width < 0)
        ABORT("sym_rect_width == %i.\n", sym_rect_width);
    sym_rect_height = n_y;
    if (sym_rect_width > 0)
    {
        rect_src = cvRect(rect.x, 0, sym_rect_width, sym_rect_height);
        rect_dst = cvRect(0, 0, sym_rect_width, sym_rect_height);
        cvGetSubRect(image, &matrix_src, rect_src);
        cvGetSubRect(image, &matrix_dst, rect_dst);
        cvFlip(&matrix_src, &matrix_dst, 1);
//        cvSet(&matrix_dst, cvRealScalar(666.0), NULL); // DEBUG.
    }
    // END: Left top, center, and bottom.

    // BEGIN: Right top, center, and bottom.
    // NOTE: A flip about y.
    sym_rect_width = n_x - (rect.x + rect.width);
    if (sym_rect_width < 0)
        ABORT("sym_rect_width == %i.\n", sym_rect_width);
    sym_rect_height = n_y;
    if (sym_rect_width > 0)
    {
        rect_src = cvRect((rect.x + rect.width) - (n_x - (rect.x + rect.width)), 0, sym_rect_width, sym_rect_height);
        rect_dst = cvRect(rect.x + rect.width, 0, sym_rect_width, sym_rect_height);
        cvGetSubRect(image, &matrix_src, rect_src);
        cvGetSubRect(image, &matrix_dst, rect_dst);
        cvFlip(&matrix_src, &matrix_dst, 1);
//        cvSet(&matrix_dst, cvRealScalar(666.0), NULL); // DEBUG.
    }
    // END: Right top, center, and bottom.

// PRIMARY DEBUG:     write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.symmetrize_IplImage_values_CvRect.mrc"); // PRIMARY DEBUG.
    //ABORT("Abort!\n");
}

void cvSetZero_IplImage_values_outside_CvRect(IplImage *image, CvRect rect)
{
    if ((image->depth != IPL_DEPTH_32F) && (image->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, image->depth);

    //write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.mrc");

    int n_x = image->width;
    int n_y = image->height;

    CvMat matrix_src;
    CvMat matrix_dst;

    // BEGIN: Top center.
    cvGetSubRect(image, &matrix_src, cvRect(rect.x, (rect.y + rect.height) - (n_y - (rect.y + rect.height)), rect.width, (n_y - (rect.y + rect.height))));
    cvGetSubRect(image, &matrix_dst, cvRect(rect.x, rect.y + rect.height, rect.width, (n_y - (rect.y + rect.height))));
    cvSetZero(&matrix_dst);
    // END: Top center.

    // BEGIN: Bottom center.
    // NOTE: A flip about x.
    cvGetSubRect(image, &matrix_src, cvRect(rect.x, rect.y, rect.width, rect.y));
    cvGetSubRect(image, &matrix_dst, cvRect(rect.x, 0, rect.width, rect.y));
    cvSetZero(&matrix_dst);
    // END: Bottom center.

    // BEGIN: Left top, center, and bottom.
    // NOTE: A flip about y.
    cvGetSubRect(image, &matrix_src, cvRect(rect.x, 0, rect.x, n_y));
    cvGetSubRect(image, &matrix_dst, cvRect(0, 0, rect.x, n_y));
    cvSetZero(&matrix_dst);
    // END: Left top, center, and bottom.

    // BEGIN: Right top, center, and bottom.
    cvGetSubRect(image, &matrix_src, cvRect((rect.x + rect.width) - (n_x - (rect.x + rect.width)), 0, n_x - (rect.x + rect.width), n_y));
    cvGetSubRect(image, &matrix_dst, cvRect(rect.x + rect.width, 0, n_x - (rect.x + rect.width), n_y));
    cvSetZero(&matrix_dst);
    // END: Right top, center, and bottom.

// PRIMARY DEBUG:     write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.cvSetZero_IplImage_values_CvRect.mrc"); // PRIMARY DEBUG.
    //ABORT("Abort!\n");
}

////////////////////////////////////////////////////////////////////////////////
// END: Extend image values. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// END: Misc. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: Filtering. {
////////////////////////////////////////////////////////////////////////////////

// BEGIN: Shepp-Logan filter.

pixel_type *create_shepp_logan_filter_1D_rs(int length)
{
    if (length % 2 != 1)
        ABORT("length == %d \% 2 != 1\n", length);

    pixel_type *filter_1D = (pixel_type *) malloc(length * sizeof(pixel_type)); 
    if (filter_1D == NULL)
        ABORT("Memory request failed.\n");

    real_type n = -((real_type) length - 1.0) / 2.0;

    for (int i = 0; i < length; ++i)
    {
        //Print("n = %.15e\n", n);
        //Print("n^2 = %.15e\n", n * n);
        filter_1D[i] = (pixel_type) (-2.0 / (M_PI * ((4.0 * (n * n)) - 1.0)));
//filter_1D[i] = 1.0; // DEBUG.
//filter_1D[i] = (i != ((length - 1) / 2)) ? 0.0 : 1.0; // DEBUG.
        //Print("filter_1D[%i] = %.15e\n", i, filter_1D[i]);
         n = n + 1.0;
    }

    //for (int i = 0; i < length; ++i)
    //    Print("%.15e\n", filter_1D[i]);

    return filter_1D;
}

void create_shepp_logan_filter_1D_rsfs_cv
(
    int length_rs,
    pixel_type **filter_rs, // This is computed first.
    int length_fs, // length_rs + padding
    CvMat **filter_fs
)
{
// PRIMARY DEBUG:     Print("length_rs = %i\n", length_rs); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("length_fs = %i\n", length_fs); // PRIMARY DEBUG.

    if (length_rs % 2 != 1)
        ABORT("length_rs == %d \% 2 != 1\n", length_rs);

    if (length_fs < length_rs)
        ABORT("length_fs == %d < %d == length_rs\n", length_fs, length_rs);

    (*filter_rs) = create_shepp_logan_filter_1D_rs(length_rs);
    (*filter_fs) = cvCreateMat(1, length_fs, CVMAT_TYPE_REAL);
    cvZero(*filter_fs);
    CvMat *filter_rs_cv = create_CvMat_32F_from_float_data(1, length_rs, (*filter_rs)); // WARNING: float vs. double
    //CvMat *filter_rs_cv = create_CvMat_64F_from_double_data(1, length_rs, (*filter_rs)); // WARNING: float vs. double
    CvMat tmp;
    cvGetSubRect((*filter_fs), &tmp, cvRect(0, 0, length_rs, 1));
    cvCopy(filter_rs_cv, &tmp, NULL);
    cvReleaseMat(&filter_rs_cv); // WARNING: Do NOT free filter_rs.

    cvDFT((*filter_fs), (*filter_fs), CV_DXT_FORWARD | CV_DXT_ROWS, 1);

    //print_CvMat_packed_complex("*filter_fs", *filter_fs);
}

// END: Shepp-Logan filter.

// BEGIN: Ram-Lak filter.

pixel_type *create_ram_lak_filter_1D_rs(int length)
{
    if (length % 2 != 1)
        ABORT("length == %d \% 2 != 1\n", length);

    pixel_type *filter_1D = (pixel_type *) malloc(length * sizeof(pixel_type)); 
    if (filter_1D == NULL)
        ABORT("Cannot acquire memory for pixel_type *data.\n");

    real_type n = -((real_type) length - 1.0) / 2.0;

    for (int i = 0; i < length; ++i)
    {
        if (n == 0)
            filter_1D[i] = (pixel_type) (M_PI / 4.0);
        else if (abs(n) % 2 == 1)
            filter_1D[i] = (pixel_type) (-1 / (M_PI * (n * n)));
        else
            filter_1D[i] = 0;
        n = n + 1.0;
    }

    return filter_1D;
}

void create_ram_lak_filter_1D_rsfs_cv
(
    int length_rs,
    pixel_type **filter_rs, // This is computed first.
    int length_fs, // length_rs + padding
    CvMat **filter_fs
)
{
// PRIMARY DEBUG:     Print("length_rs = %i\n", length_rs); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("length_fs = %i\n", length_fs); // PRIMARY DEBUG.

    if (length_rs % 2 != 1)
        ABORT("length_rs == %d \% 2 != 1\n", length_rs);

    if (length_fs < length_rs)
        ABORT("length_fs == %d < %d == length_rs\n", length_fs, length_rs);

    (*filter_rs) = create_ram_lak_filter_1D_rs(length_rs);
    (*filter_fs) = cvCreateMat(1, length_fs, CVMAT_TYPE_REAL);
    cvZero(*filter_fs);
    CvMat *filter_rs_cv = create_CvMat_32F_from_float_data(1, length_rs, (*filter_rs)); // WARNING: float vs. double
    //CvMat *filter_rs_cv = create_CvMat_64F_from_double_data(1, length_rs, (*filter_rs)); // WARNING: float vs. double
    CvMat tmp;
    cvGetSubRect((*filter_fs), &tmp, cvRect(0, 0, length_rs, 1));
    cvCopy(filter_rs_cv, &tmp, NULL);
    cvReleaseMat(&filter_rs_cv); // WARNING: Do NOT free filter_rs.

    cvDFT((*filter_fs), (*filter_fs), CV_DXT_FORWARD | CV_DXT_ROWS, 1);

    //print_CvMat_packed_complex("*filter_fs", *filter_fs);
}

// END: Ram-Lak filter.

// BEGIN: Custom R-Weighted filter.

void create_custom_r_weighted_filter_1D_fsrs_cv
(
    int length_fs, // NOTE: This is the length used when calculating the real dom. filter.
    int length_fs_padded, // NOTE: This is the length used when calculating the freq. dom. filter from the real dom. filter.
    CvMat **filter_fs, // NOTE: This is computed first using length_fs, wiped out, and then recalculated using length_fs_padded.
    int length_rs,
    pixel_type **filter_rs,
    real_type cut_off,
    real_type roll_off
)
{
// PRIMARY DEBUG:     Print("length_fs = %i\n", length_fs); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("length_rs = %i\n", length_rs); // PRIMARY DEBUG.

    if (length_rs % 2 != 1)
        ABORT("length_rs == %d \% 2 != 1\n", length_rs);

    if (length_fs % 2 != 1)
        ABORT("length_fs == %d \% 2 != 1\n", length_fs);

    if (length_fs < length_rs)
        ABORT("length_fs == %d < %d == length_rs\n", length_fs, length_rs);

    // BEGIN: Use 2-channel matrix.

    *filter_fs = cvCreateMatHeader(1, length_fs, CVMAT_TYPE_COMPLEX);

    int n_bytes_row = length_fs * sizeof(complex_type);
    complex_type *filter_fs_complex_data = (complex_type *) malloc(n_bytes_row);
    if (filter_fs_complex_data == NULL)
        ABORT("Cannot acquire memory for filter_fs_data.\n");
    
    cvSetData(*filter_fs, filter_fs_complex_data, n_bytes_row);

    *filter_rs = (pixel_type *) malloc(length_rs * sizeof(pixel_type));
    if (*filter_rs == NULL)
        ABORT("Cannot acquire memory for filter_rs.\n");

    int n = (length_fs - 1) / 2;
    real_type cut_off_scaled = cut_off * (real_type) n;

// PRIMARY DEBUG:     Print("n = %i\n", n); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("cut_off_scaled = %.15e\n", cut_off_scaled); // PRIMARY DEBUG.

    complex_type *X_c = (complex_type *) (*filter_fs)->data.ptr;

    for (int i = 0; i < n + 1; ++i)
        X_c[i] = i;
    for (int i = n + 1; i < length_fs; ++i)
        X_c[i] = 2 * n - i;

//    for (int i = 0; i < length_fs; ++i)
//        Print("%.15e + %.15ei\n", creal(X_c[i]), cimag(X_c[i]));
//    Print("\n");

    for (int i = 0; i < length_fs; ++i)
        if (fabs(creal(X_c[i])) >= cut_off_scaled)
            X_c[i] = cut_off_scaled * exp(-pow(creal(X_c[i]) - cut_off_scaled, 2) / pow(roll_off, 2));

//    for (int i = 0; i < length_fs; ++i)
//        Print("%.15e + %.15ei\n", creal(X_c[i]), cimag(X_c[i]));
//    Print("\n");

    CvMat *filter_rs_cv = cvCloneMat(*filter_fs);
    complex_type *X_c_cv = (complex_type *) filter_rs_cv->data.ptr;

    //for (int i = 0; i < length_fs; ++i)
    //    Print("%.15e + %.15ei\n", creal(X_c_cv[i]), cimag(X_c_cv[i]));
    //Print("\n");

    cvDFT((*filter_fs), filter_rs_cv, CV_DXT_INV_SCALE | CV_DXT_ROWS, filter_rs_cv->rows);

//    for (int i = 0; i < length_fs; ++i)
//        Print("%.15e + %.15ei\n", creal(X_c_cv[i]), cimag(X_c_cv[i]));
//    Print("\n");

    real_type *X_r = (real_type *) malloc(length_fs * sizeof(real_type));
    if (X_r == NULL)
        ABORT("Cannot acquire memory for real_type *X_r.\n");

    // NOTE: This corresponds to MATLAB's fftshift(), not the expected ifftshift().  We have already scaled the data.
    for (int i = 0; i < n; ++i)
        X_r[i] = creal(X_c_cv[(n + 1) + i]);
    for (int i = n; i < length_fs; ++i)
        X_r[i] = creal(X_c_cv[i - n]);

    cvReleaseMat(&filter_rs_cv);

    // NOTE: We've done the ifft, so the origin (i.e., peak) is at offset n.
    // NOTE: Both length_fs and length_rs are odd.
    //int i_X_r_start = ((length_fs - 1) / 2) - ((length_rs - 1) / 2);
    //int i_X_r_stop = ((length_fs - 1) / 2) + ((length_rs - 1) / 2);
    int i_X_r_start = ((length_fs - length_rs) / 2);
    //int i_X_r_stop = (length_fs + length_rs) / 2;
    for (int i = 0; i < length_rs; ++i)
        (*filter_rs)[i] = X_r[i_X_r_start + i];

    free(X_r);

//    Print("filter_rs = [ ...\n");
//    for (int i = 0; i < length_rs; ++i) Print("%.15e\n", (*filter_rs)[i]);
//    Print("]\n");

    // BEGIN: Release 2-channel matrix, create 1-channel matrix, and take transform.
    // NOTE: I'm doing this because it produces a freq. dom. filter that matches the prototype's real dom. filter.

    free_data_cvReleaseMat(*filter_fs);
    *filter_fs = cvCreateMat(1, length_fs_padded, CVMAT_TYPE_REAL);
    cvZero(*filter_fs);

    pixel_type *filter_fs_real_data = (pixel_type *) (*filter_fs)->data.ptr;
    for (int i = 0; i < length_rs; ++i)
        filter_fs_real_data[i] = (*filter_rs)[i];

    cvDFT((*filter_fs), (*filter_fs), CV_DXT_FORWARD | CV_DXT_ROWS, (*filter_fs)->rows);

    //print_CvMat_packed_complex("*filter_fs", *filter_fs);

    // END: Release 2-channel matrix, create 1-channel matrix, and take transform.

    // END: Use 2-channel matrix.
}

// END: Custom R-Weighted filter.

typedef struct 
{
    filter_1D_params_type params;

    pixel_type *filter_rs;

// NOTE: No FFTW:     fftw_complex_type *filter_fs_fftw;

    CvMat *filter_fs_cv;

    //int workspace_row_length;
    //pixel_type *workspace; // Workspace for convolution in real space.

// NOTE: No FFTW:     int workspace_fftw_row_length;
// NOTE: No FFTW:     int workspace_fftw_n_rows;
// NOTE: No FFTW:     pixel_type *workspace_fftw; // Workspace for convolution using FFTW.

    CvMat *workspace_cv; // Workspace for convolution using OpenCV.
}
filter_1D_data_type;

filter_1D_data_type *create_filter_1D_data_n_xy(filter_1D_params_handle params, int n_x, int n_y)
{
    if (params->length_rs > n_x)
        ABORT("Filter length cannot be greater than the width of the rotated image, n_x = %i: params->length_rs = %i.\n",
            n_x, params->length_rs);

    filter_1D_data_type *filter_data = (filter_1D_data_type *) malloc(sizeof(filter_1D_data_type));
    if (filter_data == NULL)
        ABORT("Cannot acquire memory for filter_1D_data_type *filter_data.\n");
    
    filter_data->params = *params; // Shallow copy.

    filter_data->filter_rs = NULL;
// NOTE: No FFTW:     filter_data->filter_fs_fftw = NULL;
    filter_data->filter_fs_cv = NULL;
// NOTE: No real-space filtering:     filter_data->workspace_row_length = 0;
// NOTE: No real-space filtering:     filter_data->workspace = NULL;
// NOTE: No FFTW:     filter_data->workspace_fftw = NULL;
    filter_data->workspace_cv = NULL;

    int n_x_padded = filter_data->params.length_rs + n_x - 1;
    int n_x_padded_opt = cvGetOptimalDFTSize(n_x_padded);
    int n_y_opt = cvGetOptimalDFTSize(n_y);

    CvMat *filter_fs_cv;

    switch (filter_data->params.type)
    {
        case (SHEPPLOGAN) :

            filter_data->params.length_fs = n_x_padded_opt;

            create_shepp_logan_filter_1D_rsfs_cv(
                filter_data->params.length_rs,
                &(filter_data->filter_rs),
                filter_data->params.length_fs,
                &filter_fs_cv);

            break;

        case (RAMLAK) :
    
            filter_data->params.length_fs = n_x_padded_opt;

            create_ram_lak_filter_1D_rsfs_cv(
                filter_data->params.length_rs,
                &(filter_data->filter_rs),
                filter_data->params.length_fs,
                &filter_fs_cv);
    
            break;
 
        case (CUSTOMRWEIGHTED) :
        {
            filter_data->params.length_fs = (2 * n_x) + 1; // NOTE: MATLAB prototype uses this value.

            create_custom_r_weighted_filter_1D_fsrs_cv(
                filter_data->params.length_fs,
                n_x_padded_opt,
                &filter_fs_cv,
                filter_data->params.length_rs,
                &(filter_data->filter_rs),
                filter_data->params.cut_off,
                filter_data->params.roll_off);

            break;
        }
        default :

            ABORT("Unrecognized filter type.\n");
    }

// PRIMARY DEBUG:     Print("[ ...\n"); // PRIMARY DEBUG.
// PRIMARY DEBUG:     for (int i = 0; i < filter_data->params.length_rs; ++i) Print("%.15e\n", filter_data->filter_rs[i]); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("]\n"); // PRIMARY DEBUG.

// NOTE: No real-space filtering:     filter_data->workspace_row_length = n_x;
// NOTE: No real-space filtering:     int n_bytes_row = n_x * sizeof(pixel_type);
// NOTE: No real-space filtering:     filter_data->workspace = (pixel_type *) malloc(n_bytes_row);
// NOTE: No real-space filtering:     if (filter_data->workspace == NULL)
// NOTE: No real-space filtering:         ABORT("Cannot acquire memory for pixel_type *workspace.\n");
// NOTE: No real-space filtering:     memset(filter_data->workspace, 0, n_bytes_row);

    // BEGIN: OpenCV

    // NOTE: For now, OpenCV workspace contains as many filtered rows as the number of rows in the padded image.
    filter_data->workspace_cv = cvCreateMat(
        n_y_opt, 
        n_x_padded_opt,
        CVMAT_TYPE_REAL);
    cvZero(filter_data->workspace_cv);

    // Make filter the same size as the workspace and put copies of the filter in each row.
    CvMat tmp;

// PRIMARY DEBUG:     Print("filter_fs_cv->(rows, cols) = (%i, %i)\n", filter_fs_cv->rows, filter_fs_cv->cols); // PRIMARY DEBUG.
    cvGetSubRect(filter_data->workspace_cv, &tmp, cvRect(0, 0, filter_fs_cv->cols, 1));
    cvCopy(filter_fs_cv, &tmp, NULL);
    filter_data->filter_fs_cv = cvCloneMat(filter_data->workspace_cv);
    cvZero(filter_data->filter_fs_cv);
    cvGetSubRect(filter_data->workspace_cv, &tmp, cvRect(0, 0, filter_data->workspace_cv->cols, 1));
    cvRepeat(&tmp, filter_data->filter_fs_cv);

    cvReleaseMat(&filter_fs_cv);

    // END: OpenCV

    return filter_data;
}

// WARNING: This probably should not be used.
filter_1D_data_type *create_filter_1D_data(filter_1D_params_handle params, real_type support[4])
{
    WARN("Calling %s with support and without padding!\n");

#if PROTOTYPE_COMPLIANT_INDEXING
    int n_x = round(support[2]);
    int n_y = round(support[3]);
#else
    int n_x = round(support[2] + 1.0);
    int n_y = round(support[3] + 1.0);
#endif

    return create_filter_1D_data_n_xy(params, n_x, n_y);
}

void filter_1D_data_release(filter_1D_data_type *filter_data)
{
    if (filter_data == NULL)
        ABORT("Calling free() on NULL pointer.\n");

    if (filter_data->filter_rs != NULL)
        free(filter_data->filter_rs);

// NOTE: No FFTW:     if (filter_data->filter_fs_fftw != NULL)
// NOTE: No FFTW:         fftw_free(filter_data->filter_fs_fftw);

    if (filter_data->filter_fs_cv != NULL)
        cvReleaseMat(&(filter_data->filter_fs_cv));

// NOTE: No real-space filtering:     if (filter_data->workspace != NULL)
// NOTE: No real-space filtering:         free(filter_data->workspace);

// NOTE: No FFTW:     //if (filter_data->workspace_fftw != NULL)
// NOTE: No FFTW:     //    free(filter_data->workspace_fftw);

    if (filter_data->workspace_cv != NULL)
        cvReleaseMat(&(filter_data->workspace_cv));

    free(filter_data);
}

void projection_filter_1D_fs
(
    IplImage *image,
    filter_1D_data_type *filter_data,
    symmetrize_2D_params_handle symmetrize_2D_params
)
{
    int dft_rows = cvGetOptimalDFTSize(image->height); // NOTE: This may not be necessary since we are using cvDFT(CV_DXT_ROWS).
    int dft_cols = cvGetOptimalDFTSize(filter_data->params.length_rs + image->width - 1);
    //int dft_rows = image->height; // DEBUG.
    //int dft_cols = filter_data->params.length_fs; // DEBUG.

    CvMat *workspace = filter_data->workspace_cv;

    if ((dft_rows != workspace->rows) || (dft_cols != workspace->cols))
        ABORT("(dft_rows, dft_cols) = (%i, %i) != (workspace->rows, workspace->cols) = (%i, %i)\n", dft_rows, dft_cols, workspace->rows, workspace->cols);

    CvMat matrix_src;
    CvMat matrix_dst;

    int filter_offset = filter_data->params.length_rs / 2; // NOTE: Filters are assumed to be of odd length.

    if (symmetrize_2D_params->symmetrize_2D_flag)
    {
        /*
        // NOTE: Extend rows right to end of workspace.
        // NOTE: Symmetrize right only.
        cvGetSubRect(workspace, &matrix_dst, cvRect(0, 0, image->width, image->height));
        cvCopy(image, &matrix_dst, NULL);
        */

        /*
        // BEGIN: Extend rows right to end of workspace.
        cvGetSubRect(image, &matrix_src, cvRect(image->width - 1, 0, 1, image->height));
        cvRepeat(&matrix_src, &matrix_dst);
        // END: Extend rows right to end of workspace.
        */

        /*
        // BEGIN: Symmetrize right only.  Only a single mirrored copy is necessary.
        // NOTE: A flip about y.
        int image_flip_width = fmin(image->width, workspace->cols - image->width);
        cvGetSubRect(workspace, &matrix_dst, cvRect(image->width, 0, image_flip_width, image->height));
        cvGetSubRect(image, &matrix_src, cvRect(image->width - image_flip_width, 0, image_flip_width, image->height));
        cvFlip(&matrix_src, &matrix_dst, 1);
        if (image_flip_width < (workspace->cols - image->width))
        {
            cvGetSubRect(workspace, &matrix_dst, cvRect(image->width + image_flip_width, 0, workspace->cols - (image->width + image_flip_width), image->height));
            cvZero(&matrix_dst);
        }
        // END: Symmetrize right only.  Only a single mirrored copy is necessary.
        */

        // BEGIN: Symmetrize left and right.  Only a single mirrored copy on each side is necessary.
        cvGetSubRect(workspace, &matrix_dst, cvRect(filter_offset, 0, image->width, image->height));
        cvCopy(image, &matrix_dst, NULL);

        // BEGIN: A flip about y to the left.
        int image_flip_width_left = fmin(image->width, filter_offset);
        cvGetSubRect(workspace, &matrix_dst, cvRect(filter_offset - image_flip_width_left, 0, image_flip_width_left, image->height));
        cvGetSubRect(image, &matrix_src, cvRect(0, 0, image_flip_width_left, image->height));
        cvFlip(&matrix_src, &matrix_dst, 1);
        if (image->width < filter_offset)
        {
            cvGetSubRect(workspace, &matrix_dst, cvRect(0, 0, filter_offset - image->width, image->height));
            cvZero(&matrix_dst);
        }
        // END: A flip about y to the left.
        // BEGIN: A flip about y to the right.
        int image_flip_width_right = fmin(image->width, workspace->cols - (image->width + filter_offset));
        cvGetSubRect(workspace, &matrix_dst, cvRect(image->width + filter_offset, 0, image_flip_width_right, image->height));
        cvGetSubRect(image, &matrix_src, cvRect(image->width - image_flip_width_right, 0, image_flip_width_right, image->height));
        cvFlip(&matrix_src, &matrix_dst, 1);
        if (image_flip_width_right < (workspace->cols - (image->width + filter_offset)))
        {
            cvGetSubRect(workspace, &matrix_dst, cvRect(image->width + filter_offset + image_flip_width_right, 0, workspace->cols - (image->width + filter_offset + image_flip_width_right), image->height));
            cvZero(&matrix_dst);
        }
        // END: A flip about y to the right.
        // END: Symmetrize left and right.  Only a single mirrored copy on each side is necessary.
    }
    else
    {
        cvGetSubRect(workspace, &matrix_dst, cvRect(0, 0, image->width, image->height));
        cvCopy(image, &matrix_dst, NULL);

        cvGetSubRect(workspace, &matrix_dst, cvRect(image->width, 0, workspace->cols - image->width, image->height));
        cvZero(&matrix_dst);
        //cvSet(&matrix_dst, cvRealScalar(666.0), NULL);
    }

    //print_CvMat("workspace", workspace);
// PRIMARY DEBUG:     write_MRC_image_from_CvMat(workspace, "/home/akulo/filter_1D_images/workspace.mrc"); // PRIMARY DEBUG.
    //ABORT("ABORT!\n");

    cvDFT(workspace, workspace, CV_DXT_FORWARD | CV_DXT_ROWS, image->height);

    //print_CvMat_packed_complex("workspace", workspace);
    //print_CvMat_packed_complex("filter_data->filter_fs_cv", filter_data->filter_fs_cv);

    cvMulSpectrums(workspace, filter_data->filter_fs_cv, workspace, CV_DXT_ROWS);

    cvDFT(workspace, workspace, CV_DXT_INV_SCALE | CV_DXT_ROWS, image->height);

    //print_CvMat("workspace", workspace);
// PRIMARY DEBUG:     write_MRC_image_from_CvMat(workspace, "/home/akulo/filter_1D_images/workspace.filtered.mrc"); // PRIMARY DEBUG.
    //ABORT("ABORT!\n");

    if (symmetrize_2D_params->symmetrize_2D_flag)
    {
        /*
        // NOTE: Extend rows right to end of workspace.
        // NOTE: Symmetrize right only.
        cvGetSubRect(workspace, &matrix_src, cvRect(filter_offset, 0, image->width, image->height));
        cvCopy(&matrix_src, image, NULL);
        */

        // NOTE: Symmetrize left and right.
        cvGetSubRect(workspace, &matrix_src, cvRect(filter_offset * 2, 0, image->width, image->height));
        cvCopy(&matrix_src, image, NULL);
    }
    else
    {
        cvGetSubRect(workspace, &matrix_src, cvRect(filter_offset, 0, image->width, image->height));
        cvCopy(&matrix_src, image, NULL);
    }
}

void slice_filter_1D_fs_cv(Islice *slice, filter_1D_data_type *filter_data)
{
    CvMat *slice_cv = create_CvMat_from_Islice_copy(slice, CVMAT_TYPE_REAL);

    int dft_rows = cvGetOptimalDFTSize(slice_cv->rows); // NOTE: This may not be necessary since we are using cvDFT(CV_DXT_ROWS).
    int dft_cols = cvGetOptimalDFTSize(filter_data->params.length_rs + slice_cv->cols - 1);
    //int dft_rows = slice_cv->rows; // DEBUG.
    //int dft_cols = filter_data->params.length_fs; // DEBUG.

    CvMat *workspace = filter_data->workspace_cv;

    if ((dft_rows != workspace->rows) || (dft_cols != workspace->cols))
        ABORT("(dft_rows, dft_cols) = (%i, %i) != (workspace->rows, workspace->cols) = (%i, %i)\n", dft_rows, dft_cols, workspace->rows, workspace->cols);

    CvMat tmp;

    cvGetSubRect(workspace, &tmp, cvRect(0, 0, slice_cv->cols, slice_cv->rows));
    cvCopy(slice_cv, &tmp, NULL);
    cvGetSubRect(workspace, &tmp, cvRect(slice_cv->cols, 0, workspace->cols - slice_cv->cols, slice_cv->rows));
    cvZero(&tmp);
    //cvSet(&tmp, cvRealScalar(666.0), NULL);

    //print_CvMat("workspace", workspace);

    cvDFT(workspace, workspace, CV_DXT_FORWARD | CV_DXT_ROWS, slice_cv->rows);

    //print_CvMat_packed_complex("workspace", workspace);
    //print_CvMat_packed_complex("filter_data->filter_fs_cv", filter_data->filter_fs_cv);

    cvMulSpectrums(workspace, filter_data->filter_fs_cv, workspace, CV_DXT_ROWS);

    cvDFT(workspace, workspace, CV_DXT_INV_SCALE | CV_DXT_ROWS, slice_cv->rows);

    //print_CvMat("workspace", workspace);

    cvGetSubRect(workspace, &tmp, cvRect(round(filter_data->params.length_rs / 2), 0, slice_cv->cols, slice_cv->rows));
    cvCopy(&tmp, slice_cv, NULL);
    copy_CvMat_to_Islice(slice_cv, slice);

    free_data_cvReleaseMat(slice_cv);
}

void slice_filter_1D_rs(Islice *slice, pixel_type *filter_1D, int filter_length, pixel_type *row_filtered)
{
//float v_clear = 1.0;
//sliceClear(slice, &v_clear);
//return 1;

    int n_x = slice->xsize;
    int n_y = slice->ysize;
    float *slice_data_f = slice->data.f;

    int filter_half_length = filter_length / 2;
    int x_full_stop = n_x - filter_half_length;
    for (int i_y = 0; i_y < n_y; ++i_y)
    {
        int i_x = 0;
        int i_filter_start = filter_half_length;
        for (; i_x < filter_half_length; ++i_x) 
        {
            pixel_type v = 0.0;
            for (int i_filter = i_filter_start; i_filter < filter_length; ++i_filter) 
            {
                v += filter_1D[i_filter] * INDEX_2D(slice_data_f, n_x, i_x + (i_filter - filter_half_length), i_y);
            }
            row_filtered[i_x] = v;
            --i_filter_start;
        }
        //if (int i_x != filter_half_length) exit(0); // NOTE: This is what we expect.
        for (; i_x < x_full_stop; ++i_x) 
        {
            pixel_type v = 0.0;
            for (int i_filter = 0; i_filter < filter_length; ++i_filter) 
            {
                v += filter_1D[i_filter] * INDEX_2D(slice_data_f, n_x, i_x + (i_filter - filter_half_length), i_y);
            }
            row_filtered[i_x] = v;
        }
        int i_filter_stop = filter_length - 1;
        for (; i_x < n_x; ++i_x) 
        {
            pixel_type v = 0.0;
            for (int i_filter = 0; i_filter < i_filter_stop; ++i_filter) 
            {
                v += filter_1D[i_filter] * INDEX_2D(slice_data_f, n_x, i_x + (i_filter - filter_half_length), i_y);
            }
            row_filtered[i_x] = v;
            --i_filter_stop;
        }

        for (i_x = 0; i_x < n_x; ++i_x)
            PUTVAL_2D(slice_data_f, n_x, i_x, i_y, row_filtered[i_x]);
    }
}

////////////////////////////////////////////////////////////////////////////////
// END: Filtering. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: 2D remap. {
////////////////////////////////////////////////////////////////////////////////

typedef struct
{
    remap_2D_params_type params;
    real_type *map_2D_current; // The 2D map corresponding to the current projection.
}
remap_2D_data_type;

remap_2D_data_type *create_remap_2D_data(remap_2D_params_handle params)
{
    remap_2D_data_type *remap_data = (remap_2D_data_type *) malloc(sizeof(remap_2D_data_type));
    if (remap_data == NULL)
        ABORT("Cannot acquire memory for remap_2D_data_type *remap_data.\n");

    remap_data->params = *params; // Shallow copy.

    remap_data->map_2D_current = NULL; // Points into remap_data->params.map_2D.

    return remap_data;
}

void set_map_2D_current(remap_2D_data_type *remap_2D_data, int i_projection)
{
    int offset = i_projection * remap_2D_data->params.map_2D_n_x1 * remap_2D_data->params.map_2D_n_x2;
    remap_2D_data->map_2D_current = remap_2D_data->params.map_2D + offset;

    print_matrix_2D("remap_2D_data->map_2D_current", remap_2D_data->map_2D_current, remap_2D_data->params.map_2D_n_x1, remap_2D_data->params.map_2D_n_x2);
}

void remap_2D_data_release(remap_2D_data_type *remap_data)
{
    free(remap_data);
}

#define PROJECTION_REMAP_2D_FLIP_XY 0

real_type f_x1
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
            MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 0, i) *
            pow(x1, MATRIX_INDEX_2D(power_order, power_order_n_x2, 0, i)) *
            pow(x2, MATRIX_INDEX_2D(power_order, power_order_n_x2, 1, i));

    return x1_out;
}

real_type f_x2
(
    real_type x1, real_type x2,
    int n_coeffs,
    int power_order_n_x1, int power_order_n_x2, real_type *power_order,
    int map_2D_n_x1, int map_2D_n_x2, real_type *map_2D
)
{
    real_type x2_out = 0.0;

    for (int i = 0; i < n_coeffs; ++i)
        x2_out += 
            MATRIX_INDEX_2D(map_2D, map_2D_n_x2, 1, i) *
            pow(x1, MATRIX_INDEX_2D(power_order, power_order_n_x2, 0, i)) *
            pow(x2, MATRIX_INDEX_2D(power_order, power_order_n_x2, 1, i));

    return x2_out;
}

#if PROJECTION_REMAP_2D_FLIP_XY
#define projection_solve_ivp_f_x_i_MACRO(i_f_x_, x_start_, y_start_) \
do { \
    real_type x ## i_f_x_ = 0.0; \
    for (int i = 0; i < n_coeffs; ++i) \
        x ## i_f_x_ +=  \
            MATRIX_INDEX_2D(map_2D, map_2D_n_x2, i_f_x_ - 1, i) * \
            pow(y_start_, MATRIX_INDEX_2D(power_order, power_order_n_x2, 0, i)) * \
            pow((((real_type) j) * step) + x_start_, MATRIX_INDEX_2D(power_order, power_order_n_x2, 1, i)); \
    MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, 0, j) =  x ## i_f_x_; \
} while (0)
#else
#define projection_solve_ivp_f_x_i_MACRO(i_f_x_, x_start_, y_start_) \
do { \
    real_type x ## i_f_x_ = 0.0; \
    for (int i = 0; i < n_coeffs; ++i) \
        x ## i_f_x_ +=  \
            MATRIX_INDEX_2D(map_2D, map_2D_n_x2, i_f_x_ - 1, i) * \
            pow((((real_type) j) * step) + x_start_, MATRIX_INDEX_2D(power_order, power_order_n_x2, 0, i)) * \
            pow(y_start_, MATRIX_INDEX_2D(power_order, power_order_n_x2, 1, i)); \
    MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, 0, j) =  x ## i_f_x_; \
} while (0)
#endif

// NOTE: Inlined with f_x1 and f_x2 inlined.
#define projection_remap_2D_solve_ivp_MACRO(i_f_x_, x_start_, y_start_, i_ivp_) \
do { \
    for(int j = 0; j < diff_2D_n_x2; ++j) \
        projection_solve_ivp_f_x_i_MACRO(i_f_x_, x_start_, y_start_); \
 \
   q ## i_ivp_[0] = MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, 0, 0); \
 \
   for (int i = 1; i < diff_2D_n_x1; ++i)     \
   { \
       for (int j = 0; j < diff_2D_n_x2 - i; ++j)     \
       { \
           MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i, j) =  \
               MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i - 1, j + 1) - \
               MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i - 1, j);\
       } \
       q ## i_ivp_[i] = MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i, 0); \
   } \
} while (0)

/*
// NOTE: Inlined.
#define projection_remap_2D_solve_ivp_MACRO(i_f_x_, x_start_, y_start_, i_ivp_) \
do { \
    for(int j = 0; j < diff_2D_n_x2; ++j) \
    { \
        MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, 0, j) =  \
            f_x ## i_f_x_ ##_MACRO \
            ( \
                y_start_, (((real_type) j) * step) + x_start_, \ // WARNING: Why is the xy-flip hardcoded here?
                n_coeffs, \
                power_order_n_x1, power_order_n_x2, power_order, \
                map_2D_n_x1, map_2D_n_x2, map_2D \
            ); \
    } \
 \
   q ## i_ivp_[0] = MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, 0, 0); \
 \
   for (int i = 1; i < diff_2D_n_x1; ++i)     \
   { \
       for (int j = 0; j < diff_2D_n_x2 - i; ++j)     \
       { \
           MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i, j) =  \
               MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i - 1, j + 1) - \
               MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i - 1, j);\
       } \
       q ## i_ivp_[i] = MATRIX_INDEX_2D(diff_2D, diff_2D_n_x2, i, 0); \
   } \
} while (0)
*/

// NOTE: Consider inlining this.
void projection_remap_2D_solve_ivp
(
    int n_coeffs,
    int power_order_n_x1, int power_order_n_x2, real_type *power_order,
    int map_2D_n_x1, int map_2D_n_x2, real_type *map_2D,
    real_type (*f)(real_type, real_type, int, int, int, real_type *, int, int, real_type *),
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
#if PROJECTION_REMAP_2D_FLIP_XY
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

/*
// BEGIN: Zero out ivp.
//    int n_ivp_bytes = (map_2D_order + 1) * sizeof(real_type);
//    memset(ivp, 0, n_ivp_bytes);
   int ivp_n_elems = map_2D_order + 1;
   for (int i = 0; i < ivp_n_elems; ++i)
       ivp[i] = 0.0;
// END: Zero out ivp.
*/

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

// NOTE: Inlined.
#define projection_remap_2D_scaling_MACRO(x1_, x2_, lambda_pair_) \
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

void projection_remap_2D_scaling
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

////////////////////////////////////////////////////////////////////////////////
// BEGIN: TxBR 2.0 and TxBR 3.0 code. {
////////////////////////////////////////////////////////////////////////////////

// NOTE: A 1-pixel-by-1-pixel image has a support width and height of 0.
void support_remap_2D(remap_2D_data_type *remap_2D_data, rect_2DR_type *support)
{
    Print("support = (%.15e, %.15e, %.15e, %.15e)\n", support->x, support->y, support->width, support->height);
 
    real_type p0[4][2] = // Image coordinates of corners.
        {
         {(real_type) support->x, (real_type) support->y},
         {(real_type) (support->x + support->width), (real_type) support->y},
         {(real_type) support->x, (real_type) (support->y + support->height)},
         {(real_type) (support->x + support->width), (real_type) (support->y + support->height)}
        };

    Print("p0[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n",
        p0[0][0], p0[0][1], p0[1][0], p0[1][1], p0[2][0], p0[2][1], p0[3][0], p0[3][1]);

    // BEGIN: Unpack remap 2D parameters and data.

    int map_2D_order = remap_2D_data->params.map_2D_order;
    int n_coeffs = remap_2D_data->params.n_coeffs;

    int power_order_n_x1 = remap_2D_data->params.power_order_n_x1;
    int power_order_n_x2 = remap_2D_data->params.power_order_n_x2;
    real_type *power_order = remap_2D_data->params.power_order;

    int map_2D_n_x1 = remap_2D_data->params.map_2D_n_x1;
    int map_2D_n_x2 = remap_2D_data->params.map_2D_n_x2;
    real_type *map_2D = remap_2D_data->map_2D_current;

    // END: Unpack remap 2D parameters and data.

    real_type p[4][2]; // Image coordinates of remapped corners.

    real_type lambda_2D[2];

    for (int i_p = 0; i_p < 4; ++i_p)
    {
        //Print("p0[%i][0-1] = (%.15e, %.15e)\n", i_p0, p0[i_p][0], p0[i_p][1]);

        real_type x1 = f_x1
        (
#if PROJECTION_REMAP_2D_FLIP_XY
            p0[i_p][1], p0[i_p][0],
#else
            p0[i_p][0], p0[i_p][1],
#endif
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D
        );

        real_type x2 = f_x2
        (
#if PROJECTION_REMAP_2D_FLIP_XY
            p0[i_p][1], p0[i_p][0],
#else
            p0[i_p][0], p0[i_p][1],
#endif
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D
        );

        //Print("x(1, 2) = (%.15e, %.15e)\n", x1, x2);

#if PROJECTION_REMAP_2D_FLIP_XY
        projection_remap_2D_scaling(p0[i_p][1], p0[i_p][0], map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
#else
        projection_remap_2D_scaling(p0[i_p][0], p0[i_p][1], map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
#endif

        //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);

#if PROJECTION_REMAP_2D_FLIP_XY
        p[i_p][1] = (x1 / lambda_2D[0]);
        p[i_p][0] = (x2 / lambda_2D[1]);
#else
        p[i_p][0] = (x1 / lambda_2D[0]);
        p[i_p][1] = (x2 / lambda_2D[1]);
#endif

        //Print("p[%i][0-1] = (%.15e, %.15e)\n", i_p, p[i_p][0], p[i_p][1]);
    }

    //Print("p[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n",
    //    p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]);

/*
    // BEGIN: DEBUG: Calculate transformed center0. {

    real_type center0[2] = {support0->x + (support0->width / 2.0), support0-> + (support0->height / 2.0)};

    real_type x1 = f_x1
    (
#if PROJECTION_REMAP_2D_FLIP_XY
        center0[1], center0[0],
#else
        center0[0], center0[1],
#endif
        n_coeffs,
        power_order_n_x1, power_order_n_x2, power_order,
        map_2D_n_x1, map_2D_n_x2, map_2D
    );

    real_type x2 = f_x2
    (
#if PROJECTION_REMAP_2D_FLIP_XY
        center0[1], center0[0],
#else
        center0[0], center0[1],
#endif
        n_coeffs,
        power_order_n_x1, power_order_n_x2, power_order,
        map_2D_n_x1, map_2D_n_x2, map_2D
    );

    Print("x(1, 2) = (%.15e, %.15e)\n", x1, x2);

#if PROJECTION_REMAP_2D_FLIP_XY
    projection_remap_2D_scaling(center0[1], center0[0], map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
#else
    projection_remap_2D_scaling(center0[0], center0[1], map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
#endif

    //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);

    real_type center0_transformed[2];

#if PROJECTION_REMAP_2D_FLIP_XY
    center0_transformed[1] = (x1 / lambda_2D[0]);
    center0_transformed[0] = (x2 / lambda_2D[1]);
#else
    center0_transformed[0] = (x1 / lambda_2D[0]);
    center0_transformed[1] = (x2 / lambda_2D[1]);
#endif

    Print("center0_transformed = (%.15e, %.15e)\n", center0_transformed[0], center0_transformed[1]);

    // END: DEBUG: Calculate transformed center0. }
*/

    // BEGIN: Build a bounding rectangle.

    real_type x_min = p[0][0];
    x_min = fmin(x_min, p[1][0]);
    x_min = fmin(x_min, p[2][0]);
    x_min = fmin(x_min, p[3][0]);

    real_type x_max = p[0][0];
    x_max = fmax(x_max, p[1][0]);
    x_max = fmax(x_max, p[2][0]);
    x_max = fmax(x_max, p[3][0]);

    real_type y_min = p[0][1];
    y_min = fmin(y_min, p[1][1]);
    y_min = fmin(y_min, p[2][1]);
    y_min = fmin(y_min, p[3][1]);

    real_type y_max = p[0][1];
    y_max = fmax(y_max, p[1][1]);
    y_max = fmax(y_max, p[2][1]);
    y_max = fmax(y_max, p[3][1]);

// PRIMARY DEBUG:     Print("(xy)_min = (%.15e, %.15e)\n", x_min, y_min); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("(xy)_max = (%.15e, %.15e)\n", x_max, y_max); // PRIMARY DEBUG.

    // END: Build a bounding rectangle.

    support->x = x_min;
    support->y = y_min;
    support->width = x_max - x_min;
    support->height = y_max - y_min;

/*
    real_type min_xy_trans[2] = { x_min - x0_min, y_min - y0_min };

    p[0][0] += min_xy_trans[0]; 
    p[0][1] += min_xy_trans[1];

    p[1][0] += min_xy_trans[0];
    p[1][1] += min_xy_trans[1];

    p[2][0] += min_xy_trans[0];
    p[2][1] += min_xy_trans[1];

    p[3][0] += min_xy_trans[0];
    p[3][1] += min_xy_trans[1];
*/

    Print("p[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n",
        p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]);

    Print("support = (%.15e, %.15e, %.15e, %.15e)\n", support->x, support->y, support->width, support->height);

    //ABORT("Abort!\n");
}

// NOTE: Pushes values from projection0 to projection.
void projection_remap_2D
(
    IplImage *projection0,
    rect_2DR_type *support0,
    remap_2D_data_type *remap_2D_data,
    IplImage *sums_image,
    IplImage *hits_image,
    IplImage *projection,
    rect_2DR_type *support
)
{
//#if PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

    if ((projection0->depth != IPL_DEPTH_32F) && (projection0->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection0->depth);

    if ((sums_image->depth != IPL_DEPTH_32F) && (sums_image->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, sums_image->depth);

    if ((hits_image->depth != IPL_DEPTH_32F) && (hits_image->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, hits_image->depth);

    if ((projection->depth != IPL_DEPTH_32F) && (projection->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection->depth);

    cvSetZero(sums_image);
    cvSetZero(hits_image);
    cvSetZero(projection);

    // NOTE: The exterior of this image ROI is padded.
    pixel_type *projection0_data = (pixel_type *) projection0->imageData;
    int n0_x = projection0->width;
    int n0_y = projection0->height;

    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.");

    real_type x0_min = support0->x;
    real_type width0 = support0->width;
    real_type y0_min = support0->y;
    real_type height0 = support0->height;

    real_type n0_x_min = x0_min;
    real_type n0_x_max = x0_min + width0;
    real_type n0_y_min = y0_min;
    real_type n0_y_max = y0_min + height0;

// PRIMARY DEBUG:     Print("n0_(xy) = (%i, %i)\n", n0_x, n0_y); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("n0_(xy)_min = (%.15e, %.15e), n0_(xy)_max = (%.15e, %.15e)\n", n0_x_min, n0_y_min, n0_x_max, n0_y_max); // PRIMARY DEBUG.

    // NOTE: These images are not padded and have the same extents as the unpadded projection.
    pixel_type *sums_image_data = (pixel_type *) sums_image->imageData;
    pixel_type *hits_image_data = (pixel_type *) hits_image->imageData;

    // NOTE: The exterior of this image ROI is padded.
    pixel_type *projection_data = (pixel_type *) projection->imageData;
    int n_x = projection->width;
    int n_y = projection->height;

    pixel_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.");

    real_type x_min = support->x;
    real_type width = support->width;
    real_type y_min = support->y;
    real_type height = support->height;

    real_type n_x_min = x_min;
    real_type n_x_max = x_min + width;
    real_type n_y_min = y_min;
    real_type n_y_max = y_min + height;

// PRIMARY DEBUG:     Print("n_(xy) = (%i, %i)\n", n_x, n_y); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("n_(xy)_min = (%.15e, %.15e), n_(xy)_max = (%.15e, %.15e)\n", n_x_min, n_y_min, n_x_max, n_y_max); // PRIMARY DEBUG.

    // BEGIN: Unpack remap 2D parameters and data.

    real_type step = 1.0 / remap_2D_data->params.local_mag;
    //step = 1.0;  Print("DEBUG: WARNING: step is hardcoded to %.15e.\n", step);

    int map_2D_order = remap_2D_data->params.map_2D_order;
    int n_coeffs = remap_2D_data->params.n_coeffs;

    int power_order_n_x1 = remap_2D_data->params.power_order_n_x1;
    int power_order_n_x2 = remap_2D_data->params.power_order_n_x2;
    real_type *power_order = remap_2D_data->params.power_order;

    int map_2D_n_x1 = remap_2D_data->params.map_2D_n_x1;
    int map_2D_n_x2 = remap_2D_data->params.map_2D_n_x2;
    real_type *map_2D = remap_2D_data->map_2D_current;

    // END: Unpack remap 2D parameters and data.

// PRIMARY DEBUG:     Print("step = %.15e\n", step); // PRIMARY DEBUG.
//    exit(0);

    real_type sample_factor = 1.0 / POW2(remap_2D_data->params.local_mag);
    //real_type sample_factor = 1.0 / POW2((2.0 * (remap_2D_data->params.local_mag - 1.0)) + 1.0);
/*
    real_type sample_factor = 0.0;
    if (remap_2D_data->params.local_mag > 1.0)
    {
        int half_sample_extent = ((int) remap_2D_data->params.local_mag) - 1;
        real_type x_alpha = step;
        for (int i_x = -half_sample_extent; i_x <= half_sample_extent; ++i_x)
        {
            real_type y_alpha = step;
            for (int i_y = -half_sample_extent; i_y <= half_sample_extent; ++i_y)
            {
                Print("i_(xy) = (%i, %i)\n", i_x, i_y);
                Print("(xy)_alpha = (%.15e, %.15e)\n", x_alpha, y_alpha);
                sample_factor += x_alpha * y_alpha;
                y_alpha += ((i_y >= 0) ? -1 : 1) * step;
            }
            x_alpha += ((i_x >= 0) ? -1 : 1) * step;
        }
        sample_factor = 1.0 / sample_factor;
    }
    else
    {
        sample_factor = 1.0;
    }
*/
// PRIMARY DEBUG:     Print("sample_factor = %.15e\n", sample_factor); // PRIMARY DEBUG.

    // BEGIN: Remap 2D loop initialization.

    real_type lambda_2D[2];
    real_type lambda_del_x_2D[2];

    int ivp_n_elems = map_2D_order + 1;
    int ivp_n_bytes = ivp_n_elems * sizeof(real_type);
    real_type *q1 = malloc(ivp_n_bytes);
    if (q1 == NULL)
        ABORT("Cannot acquire memory for real_type *q1.\n");
    real_type *q2 = malloc(ivp_n_bytes);
    if (q2 == NULL)
        ABORT("Cannot acquire memory for real_type *q2.\n");

    int diff_2D_n_x1 = map_2D_order + 1;
    int diff_2D_n_x2 = map_2D_order + 1;

    real_type *diff_2D = malloc(diff_2D_n_x1 * diff_2D_n_x2 * sizeof(real_type));
    if (diff_2D == NULL)
        ABORT("Cannot acquire memory for real_type *diff_2D.\n");

    // END: Remap 2D loop initialization.

    // BEGIN: DEBUG: Hardcode loop parameters.
    //n0_x_min = -100.0;
    //n0_x_max = -90.0;
    //n0_y_min = 625.0;
    //n0_y_max = 625.0;
// PRIMARY DEBUG:     //WARN("Hardcoded loop parameters: n0_(xy)_min = (%.15e, %.15e), n0_(xy)_max = (%.15e, %.15e)\n", n0_x_min, n0_y_min, n0_x_max, n0_y_max); // PRIMARY DEBUG.
    // END: DEBUG: Hardcode loop parameters.

    // This pushes values from projection0 to projection.
    for (real_type y0 = n0_y_min; y0 <= n0_y_max; y0 += step)
    {
        real_type y0_trans = y0 - y0_min;
        int y0_trans_floor = FLOOR_INNER_LOOP(y0_trans);
        real_type y0_trans_alpha = y0_trans - (real_type) y0_trans_floor;

/*
        // BEGIN: Function version of projection_remap_2D_scaling().
#if PROJECTION_REMAP_2D_FLIP_XY
        projection_remap_2D_scaling(y0, n0_x_min, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
        projection_remap_2D_scaling(y0, n0_x_min + step, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D);
#else
        projection_remap_2D_scaling(n0_x_min, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
        projection_remap_2D_scaling(n0_x_min + step, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D);
#endif
        // END: Function version of projection_remap_2D_scaling().
*/

        // BEGIN: Macro version of projection_remap_2D_scaling().
#if PROJECTION_REMAP_2D_FLIP_XY
        projection_remap_2D_scaling_MACRO(y0, n0_x_min, lambda);
        projection_remap_2D_scaling_MACRO(y0, n0_x_min + step, lambda_del_x);
#else
        projection_remap_2D_scaling_MACRO(n0_x_min, y0, lambda);
        projection_remap_2D_scaling_MACRO(n0_x_min + step, y0, lambda_del_x);
#endif
        // END: Macro version of projection_remap_2D_scaling().

        lambda_del_x_2D[0] -= lambda_2D[0];
        lambda_del_x_2D[1] -= lambda_2D[1];

        //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);
        //Print("lambda_del_x_2D[0-1] = (%.15e, %.15e)\n", lambda_del_x_2D[0], lambda_del_x_2D[1]);
        //ABORT("Abort!\n");

/*
        // BEGIN: Function version of projection_remap_2D_solve_ivp().
        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x1,
            n0_x_min, y0, step,
            map_2D_order,
            diff_2D,
            q1
        );

        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x2,
            n0_x_min, y0, step,
            map_2D_order,
            diff_2D,
            q2
        );
        // END: Function version of projection_remap_2D_solve_ivp().
*/

        // BEGIN: Macro version of projection_remap_2D_solve_ivp().
        projection_remap_2D_solve_ivp_MACRO
        (
            1,
            n0_x_min, y0,
            1
        );

        projection_remap_2D_solve_ivp_MACRO
        (
            2,
            n0_x_min, y0,
            2
        );
        // END: Macro version of projection_remap_2D_solve_ivp().

        //print_matrix_1D("q1", q1, ivp_n_elems);
        //print_matrix_1D("q2", q2, ivp_n_elems);

        //for (real_type x0 = n0_x_min; x0 <= n0_x_min + 1.0; x0 += step)
        for (real_type x0 = n0_x_min; x0 <= n0_x_max; x0 += step)
        {
            real_type x0_trans = x0 - x0_min;
            int x0_trans_floor = FLOOR_INNER_LOOP(x0_trans);
            real_type x0_trans_alpha = x0_trans - (real_type) x0_trans_floor;

#if PROJECTION_REMAP_2D_FLIP_XY
            real_type y = (q1[0] / lambda_2D[0]);
            real_type x = (q2[0] / lambda_2D[1]);
#else
            real_type x = (q1[0] / lambda_2D[0]);
            real_type y = (q2[0] / lambda_2D[1]);
#endif

            real_type x_trans = x - x_min;
            real_type y_trans = y - y_min;

            int x_trans_floor = FLOOR_INNER_LOOP(x_trans);
            int y_trans_floor = FLOOR_INNER_LOOP(y_trans);

            real_type x_trans_alpha = x_trans - (real_type) x_trans_floor;
            real_type y_trans_alpha = y_trans - (real_type) y_trans_floor;

            for (int i_power = 0; i_power < map_2D_order; ++i_power)
            {
                q1[i_power] += q1[i_power + 1];
                q2[i_power] += q2[i_power + 1];
            }

            lambda_2D[0] += lambda_del_x_2D[0];
            lambda_2D[1] += lambda_del_x_2D[1];

//            Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
//            Print("(x0_trans_floor, y0_trans_floor) = (%i, %i)\n", x0_trans_floor, y0_trans_floor);
//            Print("(x, y) = (%.15e, %.15e)\n", x, y);
//            Print("(x_trans_floor, y_trans_floor) = (%i, %i)\n", x_trans_floor, y_trans_floor);

            if (
                (0 <= x0_trans_floor) && (x0_trans_floor < n0_x - 1) &&
                (0 <= y0_trans_floor) && (y0_trans_floor < n0_y - 1) &&
                (0 <= x_trans_floor) && (x_trans_floor < n_x - 1) &&
                (0 <= y_trans_floor) && (y_trans_floor < n_y - 1)
               )
            {
                //Print("Hit.\n");

                // BEGIN: Geometry vs. pixel intensity.
                real_type v0 = INTERPOLATE_2D_LL_PIX(projection0_data, n0_x_ws, x0_trans_floor, y0_trans_floor, x0_trans_alpha, y0_trans_alpha);

//                Print("v0 = %.15e\n", v0);

                real_type prod_a = (1.0 - x_trans_alpha) * (1.0 - y_trans_alpha);
                real_type prod_b = x_trans_alpha * (1.0 - y_trans_alpha);
                real_type prod_c = (1.0 - x_trans_alpha) * y_trans_alpha;
                real_type prod_d = x_trans_alpha * y_trans_alpha;

                ADDVAL_2D_PIX(sums_image_data, n_x, x_trans_floor, y_trans_floor, prod_a * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_trans_floor + 1, y_trans_floor, prod_b * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_trans_floor, y_trans_floor + 1, prod_c * v0); 
                ADDVAL_2D_PIX(sums_image_data, n_x, x_trans_floor + 1, y_trans_floor + 1, prod_d * v0); 

                ADDVAL_2D_PIX(hits_image_data, n_x, x_trans_floor, y_trans_floor, prod_a); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_trans_floor + 1, y_trans_floor, prod_b); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_trans_floor, y_trans_floor + 1, prod_c); 
                ADDVAL_2D_PIX(hits_image_data, n_x, x_trans_floor + 1, y_trans_floor + 1, prod_d); 
                // END: Geometry vs. pixel intensity.
            }
//            else // DEBUG.
//            { // DEBUG.
//                Print("Miss.\n"); // DEBUG.
//            } // DEBUG.
        }
        //ABORT("Abort!\n");
    }

    free(q2);
    free(q1);
    free(diff_2D);

    // Attempt to fill in any holes in the output image.
    for (int i_y = 0; i_y < n_y; ++i_y)
    {
        for (int i_x = 0; i_x < n_x; ++i_x)
        {
            //Print("hits_image_data(%i, %i) = %.15e\n", i_x, i_y, INDEX_2D(hits_image_data, n_x, i_x, i_y));

            if (INDEX_2D(hits_image_data, n_x, i_x, i_y) > 0.0)
            {
                //Print("i_(xy) = (%i, %i)\n", i_x, i_y);

                INDEX_2D(projection_data, n_x_ws, i_x, i_y) =
                    INDEX_2D(sums_image_data, n_x, i_x, i_y) /
                    INDEX_2D(hits_image_data, n_x, i_x, i_y);
            }
            else if ((0 < i_x) && (i_x < n_x - 1) && (0 < i_y) && (i_y < n_y - 1))
            {
                //Print("i_(xy) (fill) = (%i, %i)\n", i_x, i_y);

                pixel_type hits_local = 0.0;
                pixel_type sum_local = 0.0;
                for (int i_x_local = -1; i_x_local <= 1; ++i_x_local)
                {
                    for (int i_y_local = -1; i_y_local <= 1; ++i_y_local)
                    {
                        hits_local += INDEX_2D(hits_image_data, n_x, i_x + i_x_local, i_y + i_y_local);
                        sum_local += INDEX_2D(sums_image_data, n_x, i_x + i_x_local, i_y + i_y_local);
                    }
                }
                if (hits_local > 0.0)
                    PUTVAL_2D(projection_data, n_x_ws, i_x, i_y, sum_local / hits_local);
                else
                    PUTVAL_2D(projection_data, n_x_ws, i_x, i_y, 0.0);
            }
        }
    }

    if (sample_factor != 1.0)
        cvScale(projection, projection, sample_factor, 0);

// PRIMARY DEBUG:     write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.projection_remap_2D.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:     write_MRC_image_from_IplImage(sums_image, "/home/akulo/filter_1D_images/sums_image.projection_remap_2D.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:     write_MRC_image_from_IplImage(hits_image, "/home/akulo/filter_1D_images/hits_image.projection_remap_2D.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:     write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.projection_remap_2D.mrc"); // PRIMARY DEBUG.

    //ABORT("Abort!");
}

// NOTE: Pulls values from projection to projection0.
void projection_inv_remap_2D
(
    IplImage *projection,
    rect_2DR_type *support,
    remap_2D_data_type *remap_2D_data,
    IplImage *projection0,
    rect_2DR_type *support0
)
{
//#if PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

    if ((projection->depth != IPL_DEPTH_32F) && (projection->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection->depth);

    if ((projection0->depth != IPL_DEPTH_32F) && (projection0->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection0->depth);

    cvSetZero(projection0);

    // NOTE: The exterior of this image ROI is padded.
    pixel_type *projection_data = (pixel_type *) projection->imageData;
    int n_x = projection->width;
    int n_y = projection->height;

    pixel_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.");

    real_type x_min = support->x;
    real_type width = support->width;
    real_type y_min = support->y;
    real_type height = support->height;

    real_type n_x_min = x_min;
    real_type n_x_max = x_min + width;
    real_type n_y_min = y_min;
    real_type n_y_max = y_min + height;

// PRIMARY DEBUG:     Print("n_(xy) = (%i, %i)\n", n_x, n_y); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("n_(xy)_min = (%.15e, %.15e), n_(xy)_max = (%.15e, %.15e)\n", n_x_min, n_y_min, n_x_max, n_y_max); // PRIMARY DEBUG.

    // NOTE: The exterior of this image ROI is padded.
    pixel_type *projection0_data = (pixel_type *) projection0->imageData;
    int n0_x = projection0->width;
    int n0_y = projection0->height;

    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.");

    real_type x0_min = support0->x;
    real_type width0 = support0->width;
    real_type y0_min = support0->y;
    real_type height0 = support0->height;

    real_type n0_x_min = x0_min;
    real_type n0_x_max = x0_min + width0;
    real_type n0_y_min = y0_min;
    real_type n0_y_max = y0_min + height0;

// PRIMARY DEBUG:     Print("n0_(xy) = (%i, %i)\n", n0_x, n0_y); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("n0_(xy)_min = (%.15e, %.15e), n0_(xy)_max = (%.15e, %.15e)\n", n0_x_min, n0_y_min, n0_x_max, n0_y_max); // PRIMARY DEBUG.

    // BEGIN: Unpack remap 2D parameters and data.

    real_type step = 1.0 / remap_2D_data->params.local_mag;
    //step = 1.0;  Print("DEBUG: WARNING: step is hardcoded to %.15e.\n", step);

    int map_2D_order = remap_2D_data->params.map_2D_order;
    int n_coeffs = remap_2D_data->params.n_coeffs;

    int power_order_n_x1 = remap_2D_data->params.power_order_n_x1;
    int power_order_n_x2 = remap_2D_data->params.power_order_n_x2;
    real_type *power_order = remap_2D_data->params.power_order;

    int map_2D_n_x1 = remap_2D_data->params.map_2D_n_x1;
    int map_2D_n_x2 = remap_2D_data->params.map_2D_n_x2;
    real_type *map_2D = remap_2D_data->map_2D_current;

    // END: Unpack remap 2D parameters and data.

// PRIMARY DEBUG:     Print("step = %.15e\n", step); // PRIMARY DEBUG.
//    exit(0);

    real_type sample_factor = 1.0 / POW2(remap_2D_data->params.local_mag);
    //real_type sample_factor = 1.0 / POW2((2.0 * (remap_2D_data->params.local_mag - 1.0)) + 1.0);
/*
    real_type sample_factor = 0.0;
    if (remap_2D_data->params.local_mag > 1.0)
    {
        int half_sample_extent = ((int) remap_2D_data->params.local_mag) - 1;
        real_type x_alpha = step;
        for (int i_x = -half_sample_extent; i_x <= half_sample_extent; ++i_x)
        {
            real_type y_alpha = step;
            for (int i_y = -half_sample_extent; i_y <= half_sample_extent; ++i_y)
            {
                Print("i_(xy) = (%i, %i)\n", i_x, i_y);
                Print("(xy)_alpha = (%.15e, %.15e)\n", x_alpha, y_alpha);
                sample_factor += x_alpha * y_alpha;
                y_alpha += ((i_y >= 0) ? -1 : 1) * step;
            }
            x_alpha += ((i_x >= 0) ? -1 : 1) * step;
        }
        sample_factor = 1.0 / sample_factor;
    }
    else
    {
        sample_factor = 1.0;
    }
*/
// PRIMARY DEBUG:     Print("sample_factor = %.15e\n", sample_factor); // PRIMARY DEBUG.

    // BEGIN: Remap 2D loop initialization.

    real_type lambda_2D[2];
    real_type lambda_del_x_2D[2];

    int ivp_n_elems = map_2D_order + 1;
    int ivp_n_bytes = ivp_n_elems * sizeof(real_type);
    real_type *q1 = malloc(ivp_n_bytes);
    if (q1 == NULL)
        ABORT("Cannot acquire memory for real_type *q1.\n");
    real_type *q2 = malloc(ivp_n_bytes);
    if (q2 == NULL)
        ABORT("Cannot acquire memory for real_type *q2.\n");

    int diff_2D_n_x1 = map_2D_order + 1;
    int diff_2D_n_x2 = map_2D_order + 1;

    real_type *diff_2D = malloc(diff_2D_n_x1 * diff_2D_n_x2 * sizeof(real_type));
    if (diff_2D == NULL)
        ABORT("Cannot acquire memory for real_type *diff_2D.\n");

    // END: Remap 2D loop initialization.

    // This pulls values from projection to projection0.
//    for (real_type y0 = n0_y_min; y0 <= n0_y_max; y0 += step)
    for (real_type y0 = n0_y_min; y0 < n0_y_max; y0 += step)
    {
        real_type y0_trans = y0 - y0_min;
        int y0_trans_floor = FLOOR_INNER_LOOP(y0_trans);
        real_type y0_trans_alpha = y0_trans - (real_type) y0_trans_floor;

/*
        // BEGIN: Function version of projection_remap_2D_scaling().
#if PROJECTION_REMAP_2D_FLIP_XY
        projection_remap_2D_scaling(y0, n0_x_min, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
        projection_remap_2D_scaling(y0, n0_x_min + step, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D);
#else
        projection_remap_2D_scaling(n0_x_min, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
        projection_remap_2D_scaling(n0_x_min + step, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D);
#endif
        // END: Function version of projection_remap_2D_scaling().
*/

        // BEGIN: Macro version of projection_remap_2D_scaling().
#if PROJECTION_REMAP_2D_FLIP_XY
        projection_remap_2D_scaling_MACRO(y0, n0_x_min, lambda);
        projection_remap_2D_scaling_MACRO(y0, n0_x_min + step, lambda_del_x);
#else
        projection_remap_2D_scaling_MACRO(n0_x_min, y0, lambda);
        projection_remap_2D_scaling_MACRO(n0_x_min + step, y0, lambda_del_x);
#endif
        // END: Macro version of projection_remap_2D_scaling().

        lambda_del_x_2D[0] -= lambda_2D[0];
        lambda_del_x_2D[1] -= lambda_2D[1];

        //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);
        //Print("lambda_del_x_2D[0-1] = (%.15e, %.15e)\n", lambda_del_x_2D[0], lambda_del_x_2D[1]);
        //ABORT("Abort!\n");

/*
        // BEGIN: Function version of projection_remap_2D_solve_ivp().
        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x1,
            n0_x_min, y0, step, // QUESTION: step or n0_x_min?
            map_2D_order,
            diff_2D,
            q1
        );

        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x2,
            n0_x_min, y0, step, // QUESTION: step or n0_x_min?
            map_2D_order,
            diff_2D,
            q2
        );
        // END: Function version of projection_remap_2D_solve_ivp().
*/

        // BEGIN: Macro version of projection_remap_2D_solve_ivp().
        projection_remap_2D_solve_ivp_MACRO
        (
            1,
            n0_x_min, y0,
            1
        );

        projection_remap_2D_solve_ivp_MACRO
        (
            2,
            n0_x_min, y0,
            2
        );
        // END: Macro version of projection_remap_2D_solve_ivp().

        //print_matrix_1D("q1", q1, ivp_n_elems);
        //print_matrix_1D("q2", q2, ivp_n_elems);

        //for (real_type x0 = n0_x_min; x0 <= n0_x_max; x0 += step)
        //for (real_type x0 = n0_x_min; x0 <= n0_x_min + 1.0; x0 += step)
        for (real_type x0 = n0_x_min; x0 < n0_x_max; x0 += step)
        {
            real_type x0_trans = x0 - x0_min;
            int x0_trans_floor = FLOOR_INNER_LOOP(x0_trans);
            real_type x0_trans_alpha = x0_trans - (real_type) x0_trans_floor;

#if PROJECTION_REMAP_2D_FLIP_XY
            real_type y = (q1[0] / lambda_2D[0]);
            real_type x = (q2[0] / lambda_2D[1]);
#else
            real_type x = (q1[0] / lambda_2D[0]);
            real_type y = (q2[0] / lambda_2D[1]);
#endif

            real_type x_trans = x - n_x_min;
            real_type y_trans = y - n_y_min;

            int x_trans_floor = FLOOR_INNER_LOOP(x_trans);
            int y_trans_floor = FLOOR_INNER_LOOP(y_trans);

            real_type x_trans_alpha = x_trans - (real_type) x_trans_floor;
            real_type y_trans_alpha = y_trans - (real_type) y_trans_floor;

            for (int i_power = 0; i_power < map_2D_order; ++i_power)
            {
                q1[i_power] += q1[i_power + 1];
                q2[i_power] += q2[i_power + 1];
            }

            lambda_2D[0] += lambda_del_x_2D[0];
            lambda_2D[1] += lambda_del_x_2D[1];

//            Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
//            Print("(x0_trans_floor, y0_trans_floor) = (%i, %i)\n", x0_trans_floor, y0_trans_floor);
//            Print("(x, y) = (%.15e, %.15e)\n", x, y);
//            Print("(x_trans_floor, y_trans_floor) = (%i, %i)\n", x_trans_floor, y_trans_floor);

            if (
                (0 <= x0_trans_floor) && (x0_trans_floor < n0_x - 1) &&
                (0 <= y0_trans_floor) && (y0_trans_floor < n0_y - 1) &&
                (0 <= x_trans_floor) && (x_trans_floor < n_x - 1) &&
                (0 <= y_trans_floor) && (y_trans_floor < n_y - 1)
               )
            {
                // BEGIN: Geometry vs. pixel intensity.
                real_type v = INTERPOLATE_2D_LL_PIX(projection_data, n_x_ws, x_trans_floor, y_trans_floor, x_trans_alpha, y_trans_alpha);

//                Print("v = %.15e\n", v);

                real_type prod_a = (1.0 - x0_trans_alpha) * (1.0 - y0_trans_alpha);
                real_type prod_b = x0_trans_alpha * (1.0 - y0_trans_alpha);
                real_type prod_c = (1.0 - x0_trans_alpha) * y0_trans_alpha;
                real_type prod_d = x0_trans_alpha * y0_trans_alpha;

                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_trans_floor, y0_trans_floor, prod_a * v); 
                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_trans_floor + 1, y0_trans_floor, prod_b * v); 
                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_trans_floor, y0_trans_floor + 1, prod_c * v); 
                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_trans_floor + 1, y0_trans_floor + 1, prod_d * v); 
                // END: Geometry vs. pixel intensity.
            }
//            else // DEBUG.
//            { // DEBUG.
//                Print("Miss.\n"); // DEBUG.
//            } // DEBUG.
        }
    }

    free(q2);
    free(q1);
    free(diff_2D);

    if (sample_factor != 1.0)
        cvScale(projection0, projection0, sample_factor, 0);

    //write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.mrc");
    //write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.mrc");

    //ABORT("Abort!");
}

////////////////////////////////////////////////////////////////////////////////
// END: TxBR 2.0 and TxBR 3.0 code. }
////////////////////////////////////////////////////////////////////////////////

void projection_remap_2D_prototype
(
    remap_2D_data_type *remap_data,
    IplImage *projection0,
    IplImage *sums_image,
    IplImage *hits_image,
    IplImage *projection
)
{
//#if !PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("!PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

    if ((projection0->depth != IPL_DEPTH_32F) && (projection0->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection0->depth);

    if ((sums_image->depth != IPL_DEPTH_32F) && (sums_image->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, sums_image->depth);

    if ((hits_image->depth != IPL_DEPTH_32F) && (hits_image->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, hits_image->depth);

    if ((projection->depth != IPL_DEPTH_32F) && (projection->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection->depth);

    cvSetZero(sums_image);
    cvSetZero(hits_image);
    cvSetZero(projection);

    pixel_type *projection0_data = (pixel_type *) projection0->imageData;
    int n0_x = projection0->width;
    int n0_y = projection0->height;

    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.");

    pixel_type *sums_image_data = (pixel_type *) sums_image->imageData;
    pixel_type *hits_image_data = (pixel_type *) hits_image->imageData;

    pixel_type *projection_data = (pixel_type *) projection->imageData;

    pixel_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.");

    // BEGIN: Unpack parameters and data.

    real_type step = 1.0 / remap_data->params.local_mag;

    int map_2D_order = remap_data->params.map_2D_order;
    int n_coeffs = remap_data->params.n_coeffs;

    int power_order_n_x1 = remap_data->params.power_order_n_x1;
    int power_order_n_x2 = remap_data->params.power_order_n_x2;
    real_type *power_order = remap_data->params.power_order;

    int map_2D_n_x1 = remap_data->params.map_2D_n_x1;
    int map_2D_n_x2 = remap_data->params.map_2D_n_x2;
    real_type *map_2D = remap_data->map_2D_current;

    // END: Unpack parameters and data.

    real_type sample_factor = step * step;

    real_type lambda_2D[2];
    real_type lambda_del_x_2D[2];

    int ivp_n_elems = map_2D_order + 1;
    int ivp_n_bytes = ivp_n_elems * sizeof(real_type);
    real_type *q1 = malloc(ivp_n_bytes);
    if (q1 == NULL)
        ABORT("Cannot acquire memory for real_type *q1.\n");
    real_type *q2 = malloc(ivp_n_bytes);
    if (q2 == NULL)
        ABORT("Cannot acquire memory for real_type *q2.\n");

    int diff_2D_n_x1 = map_2D_order + 1;
    int diff_2D_n_x2 = map_2D_order + 1;

    real_type *diff_2D = malloc(diff_2D_n_x1 * diff_2D_n_x2 * sizeof(real_type));
    if (diff_2D == NULL)
        ABORT("Cannot acquire memory for real_type *diff_2D.\n");

//    real_type n0_x_min = 0.0;
//    real_type n0_x_max = (real_type) n0_x - 1.0;
    real_type n0_x_min = step; // NOTE: Conforms to MATLAB prototype.
    real_type n0_x_max = (real_type) n0_x; // NOTE: Conforms to MATLAB prototype.
//    real_type n0_y_min = 0.0;
//    real_type n0_y_max = (real_type) n0_y - 1.0;
    real_type n0_y_min = step; // NOTE: Conforms to MATLAB prototype.
    real_type n0_y_max = (real_type) n0_y; // NOTE: Conforms to MATLAB prototype.

    // This pushes values from projection0 to projection.
    for (real_type y0 = n0_y_min; y0 <= n0_y_max; y0 += step)
    {
        int y0_floor = FLOOR_INNER_LOOP(y0);
        real_type y0_alpha = y0 - (real_type) y0_floor;
        --y0_floor;

#if PROJECTION_REMAP_2D_FLIP_XY
        projection_remap_2D_scaling(y0, n0_x_min, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling(y0, n0_x_min + step, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#else
        projection_remap_2D_scaling(n0_x_min, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling(n0_x_min + step, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#endif

        lambda_del_x_2D[0] -= lambda_2D[0];
        lambda_del_x_2D[1] -= lambda_2D[1];

//        Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);
//        Print("lambda_del_x_2D[0-1] = (%.15e, %.15e)\n", lambda_del_x_2D[0], lambda_del_x_2D[1]);

        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x1,
            n0_x_min, y0, step, // QUESTION: step or n0_x_min?
            map_2D_order,
            diff_2D,
            q1
        );

        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x2,
            n0_x_min, y0, step, // QUESTION: step or n0_x_min?
            map_2D_order,
            diff_2D,
            q2
        );

//        print_matrix_1D("q1", q1, ivp_n_elems);
//        print_matrix_1D("q2", q2, ivp_n_elems);

        for (real_type x0 = n0_x_min; x0 <= n0_x_max; x0 += step)
        {
            int x0_floor = FLOOR_INNER_LOOP(x0);
            real_type x0_alpha = x0 - (real_type) x0_floor;
            --x0_floor;

#if PROJECTION_REMAP_2D_FLIP_XY
            real_type y = q1[0] / lambda_2D[0];
            real_type x = q2[0] / lambda_2D[1];
#else
            real_type x = q1[0] / lambda_2D[0];
            real_type y = q2[0] / lambda_2D[1];
#endif

//            Print("(x, y) = (%.15e, %.15e)\n", x, y);

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

            lambda_2D[0] += lambda_del_x_2D[0];
            lambda_2D[1] += lambda_del_x_2D[1];

            if (
                (0 <= x0_floor) && (x0_floor < n0_x - 1) &&
                (0 <= y0_floor) && (y0_floor < n0_y - 1) &&
                (0 <= x_floor) && (x_floor < n0_x - 1) &&
                (0 <= y_floor) && (y_floor < n0_y - 1)
               )
            {
                // BEGIN: Geometry vs. pixel intensity.
                real_type v0 = INTERPOLATE_2D_LL_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor, x0_alpha, y0_alpha);

                real_type prod_a = (1.0 - x_alpha) * (1.0 - y_alpha);
                real_type prod_b = x_alpha * (1.0 - y_alpha);
                real_type prod_c = (1.0 - x_alpha) * y_alpha;
                real_type prod_d = x_alpha * y_alpha;

                ADDVAL_2D_PIX(sums_image_data, n0_x, x_floor, y_floor, prod_a * v0); 
                ADDVAL_2D_PIX(sums_image_data, n0_x, x_floor + 1, y_floor, prod_b * v0); 
                ADDVAL_2D_PIX(sums_image_data, n0_x, x_floor, y_floor + 1, prod_c * v0); 
                ADDVAL_2D_PIX(sums_image_data, n0_x, x_floor + 1, y_floor + 1, prod_d * v0); 

                ADDVAL_2D_PIX(hits_image_data, n0_x, x_floor, y_floor, prod_a); 
                ADDVAL_2D_PIX(hits_image_data, n0_x, x_floor + 1, y_floor, prod_b); 
                ADDVAL_2D_PIX(hits_image_data, n0_x, x_floor, y_floor + 1, prod_c); 
                ADDVAL_2D_PIX(hits_image_data, n0_x, x_floor + 1, y_floor + 1, prod_d); 
                // END: Geometry vs. pixel intensity.
            }
        }
        //break; //ABORT("ABORT!");
    }

    free(q2);
    free(q1);
    free(diff_2D);

    // Attempt to fill in any holes in the output image.
    for (int i_y = 0; i_y < n0_y; ++i_y)
    {
        for (int i_x = 0; i_x < n0_x; ++i_x)
        {
            if (INDEX_2D(hits_image_data, n0_x, i_x, i_y) > 0.0)
            {
                INDEX_2D(projection_data, n_x_ws, i_x, i_y) =
                    INDEX_2D(sums_image_data, n0_x, i_x, i_y) /
                    INDEX_2D(hits_image_data, n0_x, i_x, i_y);
            }
            else if ((0 < i_x) && (i_x < n0_x - 1) && (0 < i_y) && (i_y < n0_y - 1))
            {
                pixel_type hits_local = 0.0;
                pixel_type sum_local = 0.0;
                for (int i_x_local = -1; i_x_local <= 1; ++i_x_local)
                {
                    for (int i_y_local = -1; i_y_local <= 1; ++i_y_local)
                    {
                        hits_local += INDEX_2D(hits_image_data, n0_x, i_x + i_x_local, i_y + i_y_local);
                        sum_local += INDEX_2D(sums_image_data, n0_x, i_x + i_x_local, i_y + i_y_local);
                    }
                }
                if (hits_local > 0.0)
                    PUTVAL_2D(projection_data, n_x_ws, i_x, i_y, sum_local / hits_local);
                else
                    PUTVAL_2D(projection_data, n_x_ws, i_x, i_y, 0.0);
            }
        }
    }

    //write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.mrc");
    //write_MRC_image_from_IplImage(sums_image, "/home/akulo/filter_1D_images/sums_image.mrc");
    //write_MRC_image_from_IplImage(hits_image, "/home/akulo/filter_1D_images/hits_image.mrc");
    //write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.mrc");

    //ABORT("ABORT!");
}

void projection_restore_2D_prototype
(
    remap_2D_data_type *remap_data,
    IplImage *projection,
    IplImage *projection0
)
{
//#if !PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("!PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

    if ((projection->depth != IPL_DEPTH_32F) && (projection->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection->depth);

    if ((projection0->depth != IPL_DEPTH_32F) && (projection0->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection0->depth);

    cvSetZero(projection0);

    pixel_type *projection_data = (pixel_type *) projection->imageData;

    pixel_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.");

    pixel_type *projection0_data = (pixel_type *) projection0->imageData;
    int n0_x = projection0->width;
    int n0_y = projection0->height;

    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.");

    // BEGIN: Unpack parameters and data.

    real_type step = 1.0 / remap_data->params.local_mag;

    int map_2D_order = remap_data->params.map_2D_order;
    int n_coeffs = remap_data->params.n_coeffs;

    int power_order_n_x1 = remap_data->params.power_order_n_x1;
    int power_order_n_x2 = remap_data->params.power_order_n_x2;
    real_type *power_order = remap_data->params.power_order;

    int map_2D_n_x1 = remap_data->params.map_2D_n_x1;
    int map_2D_n_x2 = remap_data->params.map_2D_n_x2;
    real_type *map_2D = remap_data->map_2D_current;

    // END: Unpack parameters and data.

    real_type sample_factor = step * step;

    real_type lambda_2D[2];
    real_type lambda_del_x_2D[2];

    int ivp_n_elems = map_2D_order + 1;
    int ivp_n_bytes = ivp_n_elems * sizeof(real_type);
    real_type *q1 = malloc(ivp_n_bytes);
    if (q1 == NULL)
        ABORT("Cannot acquire memory for real_type *q1.\n");
    real_type *q2 = malloc(ivp_n_bytes);
    if (q2 == NULL)
        ABORT("Cannot acquire memory for real_type *q2.\n");

    int diff_2D_n_x1 = map_2D_order + 1;
    int diff_2D_n_x2 = map_2D_order + 1;

    real_type *diff_2D = malloc(diff_2D_n_x1 * diff_2D_n_x2 * sizeof(real_type));
    if (diff_2D == NULL)
        ABORT("Cannot acquire memory for real_type *diff_2D.\n");

//    real_type n0_x_min = 0.0;
//    real_type n0_x_max = (real_type) n0_x - 1.0;
    real_type n0_x_min = step; // NOTE: Conforms to MATLAB prototype.
    real_type n0_x_max = (real_type) n0_x; // NOTE: Conforms to MATLAB prototype.
//    real_type n0_y_min = 0.0;
//    real_type n0_y_max = (real_type) n0_y - 1.0;
    real_type n0_y_min = step; // NOTE: Conforms to MATLAB prototype.
    real_type n0_y_max = (real_type) n0_y; // NOTE: Conforms to MATLAB prototype.

    // This pulls values from projection to projection0.
    for (real_type y0 = n0_y_min; y0 <= n0_y_max; y0 += step)
    {
        int y0_floor = FLOOR_INNER_LOOP(y0);
        //real_type y0_alpha = y0 - (real_type) y0_floor;
        --y0_floor;

#if PROJECTION_REMAP_2D_FLIP_XY
        projection_remap_2D_scaling(y0, n0_x_min, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling(y0, n0_x_min + step, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#else
        projection_remap_2D_scaling(n0_x_min, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling(n0_x_min + step, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#endif

        lambda_del_x_2D[0] -= lambda_2D[0];
        lambda_del_x_2D[1] -= lambda_2D[1];

//        Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);
//        Print("lambda_del_y_2D[0-1] = (%.15e, %.15e)\n", lambda_del_y_2D[0], lambda_del_y_2D[1]);

        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x1,
            n0_x_min, y0, step, // QUESTION: step or n0_x_min?
            map_2D_order,
            diff_2D,
            q1
        );

        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x2,
            n0_x_min, y0, step, // QUESTION: step or n0_x_min?
            map_2D_order,
            diff_2D,
            q2
        );

//        print_matrix_1D("q1", q1, ivp_n_elems);
//        print_matrix_1D("q2", q2, ivp_n_elems);

        for (real_type x0 = n0_x_min; x0 <= n0_x_max; x0 += step)
        {
            int x0_floor = FLOOR_INNER_LOOP(x0);
            //real_type x0_alpha = x0 - (real_type) x0_floor;
            --x0_floor;

#if PROJECTION_REMAP_2D_FLIP_XY
            real_type y = q1[0] / lambda_2D[0];
            real_type x = q2[0] / lambda_2D[1];
#else
            real_type x = q1[0] / lambda_2D[0];
            real_type y = q2[0] / lambda_2D[1];
#endif

//            Print("(x, y) = (%.15e, %.15e)\n", x, y);

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

            lambda_2D[0] += lambda_del_x_2D[0];
            lambda_2D[1] += lambda_del_x_2D[1];

            if (
                (0 <= x0_floor) && (x0_floor < n0_x) &&
                (0 <= y0_floor) && (y0_floor < n0_y) &&
                (0 <= x_floor) && (x_floor < n0_x - 1) &&
                (0 <= y_floor) && (y_floor < n0_y - 1)
               )
            {
                // NOTE: Geometry vs. pixel intensity.
                real_type v = INTERPOLATE_2D_LL_PIX(projection_data, n_x_ws, x_floor, y_floor, x_alpha, y_alpha);
                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor, v);
            }
        }
        //break; //ABORT("ABORT!");
    }

    free(q2);
    free(q1);
    free(diff_2D);

//    write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.mrc");
//    write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.mrc");

//    ABORT("ABORT!");
}

////////////////////////////////////////////////////////////////////////////////
// END: 2D remap. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: 2D transform. {
////////////////////////////////////////////////////////////////////////////////

void calc_transformed_support_corners
(
    remap_2D_data_type *remap_2D_data,
    real_type support[4],
    real_type corners[4][2] // Image coordinates of transformed corners.
)
{
// PRIMARY DEBUG:     Print("support[0-3] = (%.15e, %.15e, %.15e, %.15e)\n", support[0], support[1], support[2], support[3]); // PRIMARY DEBUG.

    // BEGIN: Unpack remap 2D parameters and data.

    int map_2D_order = remap_2D_data->params.map_2D_order;
    int n_coeffs = remap_2D_data->params.n_coeffs;

    int power_order_n_x1 = remap_2D_data->params.power_order_n_x1;
    int power_order_n_x2 = remap_2D_data->params.power_order_n_x2;
    real_type *power_order = remap_2D_data->params.power_order;

    int map_2D_n_x1 = remap_2D_data->params.map_2D_n_x1;
    int map_2D_n_x2 = remap_2D_data->params.map_2D_n_x2;
    real_type *map_2D = remap_2D_data->map_2D_current;

    // END: Unpack remap 2D parameters and data.

    //real_type center0_x = support[0];
    //real_type center0_y = support[1];

    // p[1-4][2], four corners of the support in projection coordinates
    real_type p0[4][2];

#if PROTOTYPE_COMPLIANT_INDEXING
    real_type x0_min = 1.0;
    real_type y0_min = 1.0;
#else
    real_type x0_min = 0.0;
    real_type y0_min = 0.0;
#endif

    p0[0][0] = x0_min; p0[0][1] = y0_min; 
    p0[1][0] = support[2]; p0[1][1] = y0_min;
    p0[2][0] = x0_min; p0[2][1] = support[3];
    p0[3][0] = support[2]; p0[3][1] = support[3];

// PRIMARY DEBUG:     Print("p0[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n", // PRIMARY DEBUG.
// PRIMARY DEBUG:         p0[0][0], p0[0][1], p0[1][0], p0[1][1], p0[2][0], p0[2][1], p0[3][0], p0[3][1]); // PRIMARY DEBUG.

    // p[1-4][2], four transformed corners of the support in projection coordinates
    real_type p[4][2]; 

    real_type lambda_2D[2];

    for (int i_p = 0; i_p < 4; ++i_p)
    {
        //Print("p0[%i][0-1] = (%.15e, %.15e)\n", i_p, p0[i_p][0], p0[i_p][1]);

        real_type x1 = f_x1
        (
#if PROJECTION_REMAP_2D_FLIP_XY
            p0[i_p][1], p0[i_p][0],
#else
            p0[i_p][0], p0[i_p][1],
#endif
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D
        );

        real_type x2 = f_x2
        (
#if PROJECTION_REMAP_2D_FLIP_XY
            p0[i_p][1], p0[i_p][0],
#else
            p0[i_p][0], p0[i_p][1],
#endif
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D
        );

        //Print("x(1, 2) = (%.15e, %.15e)\n", x1, x2);

#if PROJECTION_REMAP_2D_FLIP_XY
        projection_remap_2D_scaling(p0[i_p][1], p0[i_p][0], map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
#else
        projection_remap_2D_scaling(p0[i_p][0], p0[i_p][1], map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
#endif

        //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);

#if PROJECTION_REMAP_2D_FLIP_XY
        p[i_p][1] = (x1 / lambda_2D[0]);
        p[i_p][0] = (x2 / lambda_2D[1]);
#else
        p[i_p][0] = (x1 / lambda_2D[0]);
        p[i_p][1] = (x2 / lambda_2D[1]);
#endif

        //Print("p[%i][0-1] = (%.15e, %.15e)\n", i_p, p[i_p][0], p[i_p][1]);
    }

    //Print("p[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n",
    //    p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]);

/*
    // BEGIN: DEBUG: Calculate transformed center0. {

    real_type center0[2] = {support[0], support[1]};

    real_type x1 = f_x1
    (
#if PROJECTION_REMAP_2D_FLIP_XY
        center0[1], center0[0],
#else
        center0[0], center0[1],
#endif
        n_coeffs,
        power_order_n_x1, power_order_n_x2, power_order,
        map_2D_n_x1, map_2D_n_x2, map_2D
    );

    real_type x2 = f_x2
    (
#if PROJECTION_REMAP_2D_FLIP_XY
        center0[1], center0[0],
#else
        center0[0], center0[1],
#endif
        n_coeffs,
        power_order_n_x1, power_order_n_x2, power_order,
        map_2D_n_x1, map_2D_n_x2, map_2D
    );

    Print("x(1, 2) = (%.15e, %.15e)\n", x1, x2);

#if PROJECTION_REMAP_2D_FLIP_XY
    projection_remap_2D_scaling(center0[1], center0[0], map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
#else
    projection_remap_2D_scaling(center0[0], center0[1], map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D);
#endif

    //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);

    real_type center0_transformed[2];

#if PROJECTION_REMAP_2D_FLIP_XY
    center0_transformed[1] = (x1 / lambda_2D[0]);
    center0_transformed[0] = (x2 / lambda_2D[1]);
#else
    center0_transformed[0] = (x1 / lambda_2D[0]);
    center0_transformed[1] = (x2 / lambda_2D[1]);
#endif

    Print("center0_transformed = (%.15e, %.15e)\n", center0_transformed[0], center0_transformed[1]);

    // END: DEBUG: Calculate transformed center0. }
*/

    // BEGIN: Build a bounding rectangle.

    real_type x_min = p[0][0];
    x_min = fmin(x_min, p[1][0]);
    x_min = fmin(x_min, p[2][0]);
    x_min = fmin(x_min, p[3][0]);

    real_type x_max = p[0][0];
    x_max = fmax(x_max, p[1][0]);
    x_max = fmax(x_max, p[2][0]);
    x_max = fmax(x_max, p[3][0]);

    real_type y_min = p[0][1];
    y_min = fmin(y_min, p[1][1]);
    y_min = fmin(y_min, p[2][1]);
    y_min = fmin(y_min, p[3][1]);

    real_type y_max = p[0][1];
    y_max = fmax(y_max, p[1][1]);
    y_max = fmax(y_max, p[2][1]);
    y_max = fmax(y_max, p[3][1]);

// PRIMARY DEBUG:     Print("(xy)_min = (%.15e, %.15e)\n", x_min, y_min); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("(xy)_max = (%.15e, %.15e)\n", x_max, y_max); // PRIMARY DEBUG.

    // END: Build a bounding rectangle.

#if PROTOTYPE_COMPLIANT_INDEXING
    support[0] = (x_max - x_min + 1.0) / 2.0;
    support[1] = (y_max - y_min + 1.0) / 2.0;
    support[2] = x_max - x_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
    support[3] = y_max - y_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
#else
    support[0] = (x_max - x_min) / 2.0;
    support[1] = (y_max - y_min) / 2.0;
    //support[0] = (x_max - x_min - 1.0) / 2.0; // NOTE: This doesn't look right, but it removes the discrepancy for the IPLIMAGE_PUSH versions.
    //support[1] = (y_max - y_min - 1.0) / 2.0; // NOTE: This doesn't look right, but it removes the discrepancy for the IPLIMAGE_PUSH versions.
    support[2] = x_max - x_min;
    support[3] = y_max - y_min;
    //support[2] = x_max - x_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
    //support[3] = y_max - y_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
#endif

/*
    real_type min_xy_trans[2] = { x_min - x0_min, y_min - y0_min };

    p[0][0] += min_xy_trans[0]; 
    p[0][1] += min_xy_trans[1];

    p[1][0] += min_xy_trans[0];
    p[1][1] += min_xy_trans[1];

    p[2][0] += min_xy_trans[0];
    p[2][1] += min_xy_trans[1];

    p[3][0] += min_xy_trans[0];
    p[3][1] += min_xy_trans[1];
*/

    memcpy(corners, p, sizeof(real_type) * 4 * 2);

// PRIMARY DEBUG:     Print("p[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n", // PRIMARY DEBUG.
// PRIMARY DEBUG:         p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]); // PRIMARY DEBUG.

// PRIMARY DEBUG:     Print("support[0-3] = (%.15e, %.15e, %.15e, %.15e)\n", support[0], support[1], support[2], support[3]); // PRIMARY DEBUG.
    //ABORT("Abort!\n");
}

void calc_transformed_support(remap_2D_data_type *remap_2D_data, real_type support[4])
{
    real_type corners[4][2];
    calc_transformed_support_corners(remap_2D_data, support, corners);
}

/*
// NOTE: Uses rotation matrix instead of remap.
void calc_transformed_support_corners
(
    real_type rot_matrix_2x2[2][2],
    real_type support[4],
    real_type corners[4][2] // Image coordinates of rotated corners.
)
{
    Print("support[0-3] = (%.15e, %.15e, %.15e, %.15e)\n", support[0], support[1], support[2], support[3]);

    //real_type center0_x = support[0];
    //real_type center0_y = support[1];

#if PROTOTYPE_COMPLIANT_INDEXING
    real_type half0_x = support[2] / 2.0; // WARNING: Expands support by 0.5 on each side.
    real_type half0_y = support[3] / 2.0; // WARNING: Expands support by 0.5 on each side.
#else
//    real_type half0_x = (support[2] + 1.0) / 2.0; // WARNING: Expands support by 0.5 on each side.
//    real_type half0_y = (support[3] + 1.0) / 2.0; // WARNING: Expands support by 0.5 on each side.
    real_type half0_x = support[2] / 2.0;
    real_type half0_y = support[3] / 2.0;
#endif

    // p[1-4][2], four corners of the support in projection coordinates
    real_type p0[4][2]; 

    p0[0][0] = -half0_x; p0[0][1] = -half0_y; 
    p0[1][0] = half0_x; p0[1][1] = -half0_y;
    p0[2][0] = -half0_x; p0[2][1] = half0_y;
    p0[3][0] = half0_x; p0[3][1] = half0_y;

    Print("p0[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n",
        p0[0][0], p0[0][1], p0[1][0], p0[1][1], p0[2][0], p0[2][1], p0[3][0], p0[3][1]);

    // p[1-4][2], four rotated corners of the support in projection coordinates
    real_type p[4][2]; 

/\*  
    // BEGIN: Complies with MATLAB prototype.
    p[0][0] = (rot_matrix_2x2[0][0] * p0[0][0]) + (rot_matrix_2x2[0][1] * p0[0][1]) + center0_x;
    p[0][1] = (rot_matrix_2x2[1][0] * p0[0][0]) + (rot_matrix_2x2[1][1] * p0[0][1]) + center0_y;

    p[1][0] = (rot_matrix_2x2[0][0] * p0[1][0]) + (rot_matrix_2x2[0][1] * p0[1][1]) + center0_x;
    p[1][1] = (rot_matrix_2x2[1][0] * p0[1][0]) + (rot_matrix_2x2[1][1] * p0[1][1]) + center0_y;

    p[2][0] = (rot_matrix_2x2[0][0] * p0[2][0]) + (rot_matrix_2x2[0][1] * p0[2][1]) + center0_x;
    p[2][1] = (rot_matrix_2x2[1][0] * p0[2][0]) + (rot_matrix_2x2[1][1] * p0[2][1]) + center0_y;

    p[3][0] = (rot_matrix_2x2[0][0] * p0[3][0]) + (rot_matrix_2x2[0][1] * p0[3][1]) + center0_x;
    p[3][1] = (rot_matrix_2x2[1][0] * p0[3][0]) + (rot_matrix_2x2[1][1] * p0[3][1]) + center0_y;
    // END: Complies with MATLAB prototype.
*\/

    p[0][0] = (rot_matrix_2x2[0][0] * p0[0][0]) + (rot_matrix_2x2[0][1] * p0[0][1]);
    p[0][1] = (rot_matrix_2x2[1][0] * p0[0][0]) + (rot_matrix_2x2[1][1] * p0[0][1]);

    p[1][0] = (rot_matrix_2x2[0][0] * p0[1][0]) + (rot_matrix_2x2[0][1] * p0[1][1]);
    p[1][1] = (rot_matrix_2x2[1][0] * p0[1][0]) + (rot_matrix_2x2[1][1] * p0[1][1]);

    p[2][0] = (rot_matrix_2x2[0][0] * p0[2][0]) + (rot_matrix_2x2[0][1] * p0[2][1]);
    p[2][1] = (rot_matrix_2x2[1][0] * p0[2][0]) + (rot_matrix_2x2[1][1] * p0[2][1]);

    p[3][0] = (rot_matrix_2x2[0][0] * p0[3][0]) + (rot_matrix_2x2[0][1] * p0[3][1]);
    p[3][1] = (rot_matrix_2x2[1][0] * p0[3][0]) + (rot_matrix_2x2[1][1] * p0[3][1]);

    Print("p[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n",
        p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]);

    // BEGIN: Build a bounding rectangle.

    real_type x_min = p[0][0];
    x_min = fmin(x_min, p[1][0]);
    x_min = fmin(x_min, p[2][0]);
    x_min = fmin(x_min, p[3][0]);

    real_type x_max = p[0][0];
    x_max = fmax(x_max, p[1][0]);
    x_max = fmax(x_max, p[2][0]);
    x_max = fmax(x_max, p[3][0]);

    real_type y_min = p[0][1];
    y_min = fmin(y_min, p[1][1]);
    y_min = fmin(y_min, p[2][1]);
    y_min = fmin(y_min, p[3][1]);

    real_type y_max = p[0][1];
    y_max = fmax(y_max, p[1][1]);
    y_max = fmax(y_max, p[2][1]);
    y_max = fmax(y_max, p[3][1]);

    // END: Build a bounding rectangle.

#if PROTOTYPE_COMPLIANT_INDEXING
    support[0] = (x_max - x_min + 1.0) / 2.0;
    support[1] = (y_max - y_min + 1.0) / 2.0;
    support[2] = x_max - x_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
    support[3] = y_max - y_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
#else
    support[0] = (x_max - x_min) / 2.0;
    support[1] = (y_max - y_min) / 2.0;
    //support[0] = (x_max - x_min - 1.0) / 2.0; // NOTE: This doesn't look right, but it removes the discrepancy for the IPLIMAGE_PUSH versions.
    //support[1] = (y_max - y_min - 1.0) / 2.0; // NOTE: This doesn't look right, but it removes the discrepancy for the IPLIMAGE_PUSH versions.
    support[2] = x_max - x_min;
    support[3] = y_max - y_min;
    //support[2] = x_max - x_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
    //support[3] = y_max - y_min + 1.0; // WARNING: Expands support by an additional 0.5 on each side.
#endif

    real_type half_x = support[2] / 2.0;
    real_type half_y = support[3] / 2.0;

    real_type center_trans[2] = { half_x - half0_x, half_y - half0_y };

    p[0][0] += center_trans[0]; 
    p[0][1] += center_trans[1];

    p[1][0] += center_trans[0];
    p[1][1] += center_trans[1];

    p[2][0] += center_trans[0];
    p[2][1] += center_trans[1];

    p[3][0] += center_trans[0];
    p[3][1] += center_trans[1];

    memcpy(corners, p, sizeof(real_type) * 4 * 2);

//    Print("p[1-4][0-1] = (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e), (%.15e, %.15e)\n",
//        p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1]);

    Print("support[0-3] = (%.15e, %.15e, %.15e, %.15e)\n", support[0], support[1], support[2], support[3]);
}

void calc_transformed_support(real_type rot_matrix_2x2[2][2], real_type support[4])
{
    real_type corners[4][2];
    calc_transformed_support_corners(rot_matrix_2x2, support, corners);
}
*/

// NOTE: Pushes values from projection0 to projection.
void projection_transform_2D
(
    IplImage *projection0,
    real_type support0[4],
    remap_2D_data_type *remap_2D_data,
    IplImage *sums_image,
    IplImage *hits_image,
    IplImage *projection,
    real_type support[4]
)
{
#if !PROTOTYPE_COMPLIANT_INDEXING
    ABORT("!PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
#endif

    if ((projection0->depth != IPL_DEPTH_32F) && (projection0->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection0->depth);

    if ((sums_image->depth != IPL_DEPTH_32F) && (sums_image->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, sums_image->depth);

    if ((hits_image->depth != IPL_DEPTH_32F) && (hits_image->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, hits_image->depth);

    if ((projection->depth != IPL_DEPTH_32F) && (projection->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection->depth);

    cvSetZero(sums_image);
    cvSetZero(hits_image);
    cvSetZero(projection);

    // NOTE: The exterior of this image ROI is padded at the top and right with a row and column of zeros.
    pixel_type *projection0_data = (pixel_type *) projection0->imageData;
    int n0_x = projection0->width;
    int n0_y = projection0->height;

    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.");

    real_type center0_x = support0[0];
    real_type center0_y = support0[1];

    //real_type n0_x_min = 0.0;
    //real_type n0_x_max = support0[2];
    //real_type n0_y_min = 0.0;
    //real_type n0_y_max = support0[3];
    real_type n0_x_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n0_x_max = (real_type) n0_x; // NOTE: Conforms to MATLAB prototype.
    real_type n0_y_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n0_y_max = (real_type) n0_y; // NOTE: Conforms to MATLAB prototype.

// PRIMARY DEBUG:     Print("n0_(xy) = (%i, %i), center0_(xy) = (%.15e, %.15e)\n", n0_x, n0_y, center0_x, center0_y); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("n0_(xy)_min = (%.15e, %.15e), n0_(xy)_max = (%.15e, %.15e)\n", n0_x_min, n0_y_min, n0_x_max, n0_y_max); // PRIMARY DEBUG.

    // NOTE: These images are not padded and have the same extents as the unpadded projection.
    pixel_type *sums_image_data = (pixel_type *) sums_image->imageData;
    pixel_type *hits_image_data = (pixel_type *) hits_image->imageData;

    // NOTE: The exterior of this image ROI is padded at the top and right with a row and column of zeros.
    pixel_type *projection_data = (pixel_type *) projection->imageData;
    int n_x = projection->width;
    int n_y = projection->height;

    pixel_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.");

    real_type center_x = support[0];
    real_type center_y = support[1];

    //real_type n_x_min = 0.0;
    //real_type n_x_max = support[2];
    //real_type n_y_min = 0.0;
    //real_type n_y_max = support[3];
    real_type n_x_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n_x_max = (real_type) n_x; // NOTE: Conforms to MATLAB prototype.
    real_type n_y_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n_y_max = (real_type) n_y; // NOTE: Conforms to MATLAB prototype.

// PRIMARY DEBUG:     Print("n_(xy) = (%i, %i), center_(xy) = (%.15e, %.15e)\n", n_x, n_y, center_x, center_y); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("n_(xy)_min = (%.15e, %.15e), n_(xy)_max = (%.15e, %.15e)\n", n_x_min, n_y_min, n_x_max, n_y_max); // PRIMARY DEBUG.

    // BEGIN: Unpack remap 2D parameters and data.

    real_type step = 1.0 / remap_2D_data->params.local_mag;
    //step = 1.0;  Print("DEBUG: WARNING: step is hardcoded to %.15e.\n", step);

    int map_2D_order = remap_2D_data->params.map_2D_order;
    int n_coeffs = remap_2D_data->params.n_coeffs;

    int power_order_n_x1 = remap_2D_data->params.power_order_n_x1;
    int power_order_n_x2 = remap_2D_data->params.power_order_n_x2;
    real_type *power_order = remap_2D_data->params.power_order;

    int map_2D_n_x1 = remap_2D_data->params.map_2D_n_x1;
    int map_2D_n_x2 = remap_2D_data->params.map_2D_n_x2;
    real_type *map_2D = remap_2D_data->map_2D_current;

    // END: Unpack remap 2D parameters and data.

    n0_x_min = step; // NOTE: Conforms to MATLAB prototype.
    n0_y_min = step; // NOTE: Conforms to MATLAB prototype.

// PRIMARY DEBUG:     Print("step = %.15e\n", step); // PRIMARY DEBUG.
//    exit(0);

    real_type sample_factor = step * step;

    // BEGIN: Remap 2D loop initialization.

    real_type lambda_2D[2];
    real_type lambda_del_x_2D[2];

    int ivp_n_elems = map_2D_order + 1;
    int ivp_n_bytes = ivp_n_elems * sizeof(real_type);
    real_type *q1 = malloc(ivp_n_bytes);
    if (q1 == NULL)
        ABORT("Cannot acquire memory for real_type *q1.\n");
    real_type *q2 = malloc(ivp_n_bytes);
    if (q2 == NULL)
        ABORT("Cannot acquire memory for real_type *q2.\n");

    int diff_2D_n_x1 = map_2D_order + 1;
    int diff_2D_n_x2 = map_2D_order + 1;

    real_type *diff_2D = malloc(diff_2D_n_x1 * diff_2D_n_x2 * sizeof(real_type));
    if (diff_2D == NULL)
        ABORT("Cannot acquire memory for real_type *diff_2D.\n");

    // END: Remap 2D loop initialization.

/*
    // BEGIN: DEBUG: Hardcode loop parameters.
    n0_x_min = -100.0;
    n0_x_max = -90.0;
    n0_y_min = -100.0;
    n0_y_max = -90.0;
// PRIMARY DEBUG:     WARN("Hardcoded loop parameters: n0_(xy)_min = (%.15e, %.15e), n0_(xy)_max = (%.15e, %.15e)\n", n0_x_min, n0_y_min, n0_x_max, n0_y_max); // PRIMARY DEBUG.
    // END: DEBUG: Hardcode loop parameters.
*/

    // This pushes values from projection0 to projection.
    for (real_type y0 = n0_y_min; y0 <= n0_y_max; y0 += step)
    {
        int y0_floor = FLOOR_INNER_LOOP(y0);
        real_type y0_alpha = y0 - (real_type) y0_floor;
        --y0_floor; // NOTE: Conforms to MATLAB prototype.

/*
        // BEGIN: Function version of projection_remap_2D_scaling().
#if PROJECTION_REMAP_2D_FLIP_XY
        projection_remap_2D_scaling(y0, n0_x_min, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling(y0, n0_x_min + step, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#else
        projection_remap_2D_scaling(n0_x_min, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling(n0_x_min + step, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#endif
        // END: Function version of projection_remap_2D_scaling().
*/

        // BEGIN: Macro version of projection_remap_2D_scaling().
#if PROJECTION_REMAP_2D_FLIP_XY
        projection_remap_2D_scaling_MACRO(y0, n0_x_min, lambda); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling_MACRO(y0, n0_x_min + step, lambda_del_x); // QUESTION: 2.0 * step or n0_x_min + step?
#else
        projection_remap_2D_scaling_MACRO(n0_x_min, y0, lambda); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling_MACRO(n0_x_min + step, y0, lambda_del_x); // QUESTION: 2.0 * step or n0_x_min + step?
#endif
        // END: Macro version of projection_remap_2D_scaling().

        lambda_del_x_2D[0] -= lambda_2D[0];
        lambda_del_x_2D[1] -= lambda_2D[1];

        //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);
        //Print("lambda_del_x_2D[0-1] = (%.15e, %.15e)\n", lambda_del_x_2D[0], lambda_del_x_2D[1]);
        //ABORT("Abort!\n");

/*
        // BEGIN: Function version of projection_remap_2D_solve_ivp().
        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x1,
            n0_x_min, y0, step, // QUESTION: step or n0_x_min?
            map_2D_order,
            diff_2D,
            q1
        );

        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x2,
            n0_x_min, y0, step, // QUESTION: step or n0_x_min?
            map_2D_order,
            diff_2D,
            q2
        );
        // END: Function version of projection_remap_2D_solve_ivp().
*/

        // BEGIN: Macro version of projection_remap_2D_solve_ivp().
        projection_remap_2D_solve_ivp_MACRO
        (
            1,
            n0_x_min, y0,
            1
        );

        projection_remap_2D_solve_ivp_MACRO
        (
            2,
            n0_x_min, y0,
            2
        );
        // END: Macro version of projection_remap_2D_solve_ivp().

        //print_matrix_1D("q1", q1, ivp_n_elems);
        //print_matrix_1D("q2", q2, ivp_n_elems);

        for (real_type x0 = n0_x_min; x0 <= n0_x_max; x0 += step)
        {
            int x0_floor = FLOOR_INNER_LOOP(x0);
            real_type x0_alpha = x0 - (real_type) x0_floor;
            --x0_floor; // NOTE: Conforms to MATLAB prototype.

#if PROJECTION_REMAP_2D_FLIP_XY
            real_type y = (q1[0] / lambda_2D[0]);
            real_type x = (q2[0] / lambda_2D[1]);
#else
            real_type x = (q1[0] / lambda_2D[0]);
            real_type y = (q2[0] / lambda_2D[1]);
#endif

            //Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
            //Print("(x, y) = (%.15e, %.15e)\n", x, y);

            int x_floor = FLOOR_INNER_LOOP(x);
            int y_floor = FLOOR_INNER_LOOP(y);

            real_type x_alpha = x - (real_type) x_floor;
            real_type y_alpha = y - (real_type) y_floor;
            --x_floor; // NOTE: Conforms to MATLAB prototype.
            --y_floor; // NOTE: Conforms to MATLAB prototype.

            for (int i_power = 0; i_power < map_2D_order; ++i_power)
            {
                q1[i_power] += q1[i_power + 1];
                q2[i_power] += q2[i_power + 1];
            }

            lambda_2D[0] += lambda_del_x_2D[0];
            lambda_2D[1] += lambda_del_x_2D[1];

/*
            if (
                ((x0 == n0_x_min) && (y0 == n0_y_min)) ||
                ((x0 == n0_x_min) && (y0 == n0_y_max)) ||
                ((x0 == n0_x_max) && (y0 == n0_y_min)) ||
                ((x0 == n0_x_max) && (y0 == n0_y_max))
               )
            {
                Print("=====\n");
                Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
                Print("(x0_floor, y0_floor) = (%i, %i)\n", x0_floor, y0_floor);
                Print("(x, y) = (%.15e, %.15e)\n", x, y);
                Print("(x_floor, y_floor) = (%i, %i)\n", x_floor, y_floor);
                Print("=====\n");
            }
*/

            if (
                (0 <= x0_floor) && (x0_floor < n0_x - 1) &&
                (0 <= y0_floor) && (y0_floor < n0_y - 1) &&
                (0 <= x_floor) && (x_floor < n_x - 1) &&
                (0 <= y_floor) && (y_floor < n_y - 1)
               )
            {
/*
                if (
                    ((x0 == n0_x_min) && (y0 == n0_y_min)) ||
                    ((x0 == n0_x_min) && (y0 == n0_y_max)) ||
                    ((x0 == n0_x_max) && (y0 == n0_y_min)) ||
                    ((x0 == n0_x_max) && (y0 == n0_y_max))
                   )
                {
                    Print("Hit.\n");
                }
*/

                ////if ((x_floor == 0) && (y_floor == 0))
                //if (y0_floor == 2)
                //{
                //    Print("=====\n");
                //    Print("(x0, y0) = (%.15e, %.15e)\n", x0, y0);
                //    Print("(x0_floor, y0_floor) = (%i, %i)\n", x0_floor, y0_floor);
                //    Print("(x, y) = (%.15e, %.15e)\n", x, y);
                //    Print("(x_floor, y_floor) = (%i, %i)\n", x_floor, y_floor);
                //    Print("=====\n");
                //}

                // BEGIN: Geometry vs. pixel intensity.
                real_type v0 = INTERPOLATE_2D_LL_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor, x0_alpha, y0_alpha);

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
                // END: Geometry vs. pixel intensity.
            }
        }
        //ABORT("Abort!\n");
    }

    free(q2);
    free(q1);
    free(diff_2D);

    // Attempt to fill in any holes in the output image.
    for (int i_y = 0; i_y < n_y; ++i_y)
    {
        for (int i_x = 0; i_x < n_x; ++i_x)
        {
            //Print("hits_image_data(%i, %i) = %.15e\n", i_x, i_y, INDEX_2D(hits_image_data, n_x, i_x, i_y));

            if (INDEX_2D(hits_image_data, n_x, i_x, i_y) > 0.0)
            {
                //Print("i_(xy) = (%i, %i)\n", i_x, i_y);

                INDEX_2D(projection_data, n_x_ws, i_x, i_y) =
                    INDEX_2D(sums_image_data, n_x, i_x, i_y) /
                    INDEX_2D(hits_image_data, n_x, i_x, i_y);
            }
            else if ((0 < i_x) && (i_x < n_x - 1) && (0 < i_y) && (i_y < n_y - 1))
            {
                //Print("i_(xy) (fill) = (%i, %i)\n", i_x, i_y);

                pixel_type hits_local = 0.0;
                pixel_type sum_local = 0.0;
                for (int i_x_local = -1; i_x_local <= 1; ++i_x_local)
                {
                    for (int i_y_local = -1; i_y_local <= 1; ++i_y_local)
                    {
                        hits_local += INDEX_2D(hits_image_data, n_x, i_x + i_x_local, i_y + i_y_local);
                        sum_local += INDEX_2D(sums_image_data, n_x, i_x + i_x_local, i_y + i_y_local);
                    }
                }
                if (hits_local > 0.0)
                    PUTVAL_2D(projection_data, n_x_ws, i_x, i_y, sum_local / hits_local);
                else
                    PUTVAL_2D(projection_data, n_x_ws, i_x, i_y, 0.0);
            }
        }
    }

    if (sample_factor != 1.0)
        cvScale(projection, projection, sample_factor, 0);

    write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.mrc");
    write_MRC_image_from_IplImage(sums_image, "/home/akulo/filter_1D_images/sums_image.mrc");
    write_MRC_image_from_IplImage(hits_image, "/home/akulo/filter_1D_images/hits_image.mrc");
    write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.mrc");

    //ABORT("Abort!");
}

// NOTE: Pulls values from projection to projection0.
void projection_inv_transform_2D
(
    IplImage *projection,
    real_type support[4],
    remap_2D_data_type *remap_2D_data,
    IplImage *projection0,
    real_type support0[4]
)
{
#if !PROTOTYPE_COMPLIANT_INDEXING
    ABORT("!PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
#endif

    if ((projection->depth != IPL_DEPTH_32F) && (projection->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection->depth);

    if ((projection0->depth != IPL_DEPTH_32F) && (projection0->depth != IPL_DEPTH_64F))
        ABORT("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection0->depth);

    cvSetZero(projection0);

    // NOTE: The exterior of this image ROI is padded at the top and right with a row and column of zeros.
    pixel_type *projection_data = (pixel_type *) projection->imageData;
    int n_x = projection->width;
    int n_y = projection->height;

    pixel_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.");

    real_type center_x = support[0];
    real_type center_y = support[1];

    //real_type n_x_min = 0.0;
    //real_type n_x_max = support[2];
    //real_type n_y_min = 0.0;
    //real_type n_y_max = support[3];
    real_type n_x_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n_x_max = (real_type) n_x; // NOTE: Conforms to MATLAB prototype.
    real_type n_y_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n_y_max = (real_type) n_y; // NOTE: Conforms to MATLAB prototype.

// PRIMARY DEBUG:     Print("n_(xy) = (%i, %i), center_(xy) = (%.15e, %.15e)\n", n_x, n_y, center_x, center_y); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("n_(xy)_min = (%.15e, %.15e), n_(xy)_max = (%.15e, %.15e)\n", n_x_min, n_y_min, n_x_max, n_y_max); // PRIMARY DEBUG.

    // NOTE: The exterior of this image ROI is padded at the top and right with a row and column of zeros.
    pixel_type *projection0_data = (pixel_type *) projection0->imageData;
    int n0_x = projection0->width;
    int n0_y = projection0->height;

    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.");

    real_type center0_x = support0[0];
    real_type center0_y = support0[1];

    //real_type n0_x_min = 0.0;
    //real_type n0_x_max = support0[2];
    //real_type n0_y_min = 0.0;
    //real_type n0_y_max = support0[3];
    real_type n0_x_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n0_x_max = (real_type) n0_x; // NOTE: Conforms to MATLAB prototype.
    real_type n0_y_min = 1.0; // NOTE: Conforms to MATLAB prototype.
    real_type n0_y_max = (real_type) n0_y; // NOTE: Conforms to MATLAB prototype.

// PRIMARY DEBUG:     Print("n0_(xy) = (%i, %i), center0_(xy) = (%.15e, %.15e)\n", n0_x, n0_y, center0_x, center0_y); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("n0_(xy)_min = (%.15e, %.15e), n0_(xy)_max = (%.15e, %.15e)\n", n0_x_min, n0_y_min, n0_x_max, n0_y_max); // PRIMARY DEBUG.

    // BEGIN: Unpack remap 2D parameters and data.

    real_type step = 1.0 / remap_2D_data->params.local_mag;
    //step = 1.0;  Print("DEBUG: WARNING: step is hardcoded to %.15e.\n", step);

    int map_2D_order = remap_2D_data->params.map_2D_order;
    int n_coeffs = remap_2D_data->params.n_coeffs;

    int power_order_n_x1 = remap_2D_data->params.power_order_n_x1;
    int power_order_n_x2 = remap_2D_data->params.power_order_n_x2;
    real_type *power_order = remap_2D_data->params.power_order;

    int map_2D_n_x1 = remap_2D_data->params.map_2D_n_x1;
    int map_2D_n_x2 = remap_2D_data->params.map_2D_n_x2;
    real_type *map_2D = remap_2D_data->map_2D_current;

    // END: Unpack remap 2D parameters and data.

    n0_x_min = step; // NOTE: Conforms to MATLAB prototype.
    n0_y_min = step; // NOTE: Conforms to MATLAB prototype.

// PRIMARY DEBUG:     Print("step = %.15e\n", step); // PRIMARY DEBUG.
//    exit(0);

    real_type sample_factor = step * step;
// PRIMARY DEBUG:     Print("sample_factor = %.15e\n", sample_factor); // PRIMARY DEBUG.

    // BEGIN: Remap 2D loop initialization.

    real_type lambda_2D[2];
    real_type lambda_del_x_2D[2];

    int ivp_n_elems = map_2D_order + 1;
    int ivp_n_bytes = ivp_n_elems * sizeof(real_type);
    real_type *q1 = malloc(ivp_n_bytes);
    if (q1 == NULL)
        ABORT("Cannot acquire memory for real_type *q1.\n");
    real_type *q2 = malloc(ivp_n_bytes);
    if (q2 == NULL)
        ABORT("Cannot acquire memory for real_type *q2.\n");

    int diff_2D_n_x1 = map_2D_order + 1;
    int diff_2D_n_x2 = map_2D_order + 1;

    real_type *diff_2D = malloc(diff_2D_n_x1 * diff_2D_n_x2 * sizeof(real_type));
    if (diff_2D == NULL)
        ABORT("Cannot acquire memory for real_type *diff_2D.\n");

    // END: Remap 2D loop initialization.

    // This pulls values from projection to projection.
//    for (real_type y0 = n0_y_min; y0 <= n0_y_max; y0 += step)
    for (real_type y0 = n0_y_min; y0 < n0_y_max; y0 += step)
    {
        int y0_floor = FLOOR_INNER_LOOP(y0);
        real_type y0_alpha = y0 - (real_type) y0_floor;
        --y0_floor; // NOTE: Conforms to MATLAB prototype.

#if PROJECTION_REMAP_2D_FLIP_XY
        projection_remap_2D_scaling(y0, n0_x_min, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling(y0, n0_x_min + step, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#else
        projection_remap_2D_scaling(n0_x_min, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_2D); // QUESTION: step or n0_x_min?
        projection_remap_2D_scaling(n0_x_min + step, y0, map_2D_n_x1, map_2D_n_x2, map_2D, lambda_del_x_2D); // QUESTION: 2.0 * step or n0_x_min + step?
#endif

        lambda_del_x_2D[0] -= lambda_2D[0];
        lambda_del_x_2D[1] -= lambda_2D[1];

        //Print("lambda_2D[0-1] = (%.15e, %.15e)\n", lambda_2D[0], lambda_2D[1]);
        //Print("lambda_del_x_2D[0-1] = (%.15e, %.15e)\n", lambda_del_x_2D[0], lambda_del_x_2D[1]);
        //ABORT("Abort!\n");

        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x1,
            n0_x_min, y0, step, // QUESTION: step or n0_x_min?
            map_2D_order,
            diff_2D,
            q1
        );

        projection_remap_2D_solve_ivp
        (
            n_coeffs,
            power_order_n_x1, power_order_n_x2, power_order,
            map_2D_n_x1, map_2D_n_x2, map_2D,
            f_x2,
            n0_x_min, y0, step, // QUESTION: step or n0_x_min?
            map_2D_order,
            diff_2D,
            q2
        );

        //print_matrix_1D("q1", q1, ivp_n_elems);
        //print_matrix_1D("q2", q2, ivp_n_elems);

//        for (real_type x0 = n0_x_min; x0 <= n0_x_max; x0 += step)
        for (real_type x0 = n0_x_min; x0 < n0_x_max; x0 += step)
        {
            int x0_floor = FLOOR_INNER_LOOP(x0);
            real_type x0_alpha = x0 - (real_type) x0_floor;
            --x0_floor; // NOTE: Conforms to MATLAB prototype.

#if PROJECTION_REMAP_2D_FLIP_XY
            real_type y = (q1[0] / lambda_2D[0]);
            real_type x = (q2[0] / lambda_2D[1]);
#else
            real_type x = (q1[0] / lambda_2D[0]);
            real_type y = (q2[0] / lambda_2D[1]);
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

            lambda_2D[0] += lambda_del_x_2D[0];
            lambda_2D[1] += lambda_del_x_2D[1];

            if (
                (0 <= x0_floor) && (x0_floor < n0_x - 1) &&
                (0 <= y0_floor) && (y0_floor < n0_y - 1) &&
                (0 <= x_floor) && (x_floor < n_x - 1) &&
                (0 <= y_floor) && (y_floor < n_y - 1)
               )
            {
                // BEGIN: Geometry vs. pixel intensity.
                real_type v = INTERPOLATE_2D_LL_PIX(projection_data, n_x_ws, x_floor, y_floor, x_alpha, y_alpha);

                real_type prod_a = (1.0 - x0_alpha) * (1.0 - y0_alpha);
                real_type prod_b = x0_alpha * (1.0 - y0_alpha);
                real_type prod_c = (1.0 - x0_alpha) * y0_alpha;
                real_type prod_d = x0_alpha * y0_alpha;

                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor, prod_a * v); 
                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor, prod_b * v); 
                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1, prod_c * v); 
                ADDVAL_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1, prod_d * v); 
                // END: Geometry vs. pixel intensity.
            }
        }
    }

    free(q2);
    free(q1);
    free(diff_2D);

    if (sample_factor != 1.0)
        cvScale(projection0, projection0, sample_factor, 0);

    //write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.mrc");
    //write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.mrc");

    //ABORT("Abort!");
}

////////////////////////////////////////////////////////////////////////////////
// END: 2D transform. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: 2D rotation. {
////////////////////////////////////////////////////////////////////////////////

typedef struct
{
    rotate_2D_params_type params;
}
rotate_2D_data_type;

rotate_2D_data_type *create_rotate_2D_data(rotate_2D_params_handle params, int n_projections)
{
    if (n_projections < 1)
        ABORT("n_projections == %i < 1.\n", n_projections);

    if ((params->n_angles != 1) && (params->n_angles != n_projections))
        ABORT("params->n_angles == %i != %i == n_projections,\n", params->n_angles, n_projections);

    rotate_2D_data_type *rotate_2D_data = (rotate_2D_data_type *) malloc(sizeof(rotate_2D_data_type));
    if (rotate_2D_data == NULL)
        ABORT("Cannot acquire memory for rotate_2D_data_type *rotate_2D_data.\n");

    rotate_2D_data->params = *params; // Shallow copy.

    return rotate_2D_data;
}

void calc_rotation_matrices_current
(
    rotate_2D_data_type *rotate_2D_data,
    int i_angle,
    real_type rot_matrix_2x2[2][2],
    real_type inv_rot_matrix_2x2[2][2]
)
{
    if (i_angle < 0)
        ABORT("i_angle == %i < 0.\n", i_angle);

    if ((i_angle > 0) && (rotate_2D_data->params.n_angles > 1) && (rotate_2D_data->params.n_angles <= i_angle))
        ABORT("Bad i_angle == %i, rotate_2D_data->params.n_angles == %i.\n", i_angle, rotate_2D_data->params.n_angles);

    if (rotate_2D_data->params.n_angles == 1)
        i_angle = 0;

    calc_rot_matrix_2x2(rotate_2D_data->params.angle_list[i_angle], rot_matrix_2x2);
    calc_rot_matrix_2x2(-(rotate_2D_data->params.angle_list[i_angle]), inv_rot_matrix_2x2);
}

void rotate_2D_data_release(rotate_2D_data_type *rotate_2D_data)
{
    free(rotate_2D_data);
}

#define PROJECTION_ROTATE_2D_IPLIMAGE_PULL_2D_FINITE_DIFF
//#define PROJECTION_ROTATE_2D_IPLIMAGE_PUSH

#define PROJECTION_ROTATE_2D_CHECK_BOUNDS 1

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
IplImage *projection0_global;
CvRect *projection0_center_ROI_CvRect_global;
IplImage *projection_global;
CvRect *projection_center_ROI_CvRect_global;
#endif

#if defined PROJECTION_ROTATE_2D_IPLIMAGE_PULL_2D_FINITE_DIFF

// BEGIN: TxBR 2.0 and TxBR 3.0 code. {
// NOTE: This pulls pixel values from the original projection to the new one.
// NOTE: Both dimensions use finite difference scheme.
void projection_rotate_2D
(
    IplImage *projection0,
    rect_2DR_type *support0,
    real_type local_mag,
    real_type rot_matrix_2x2[2][2],
    IplImage *projection,
    rect_2DR_type *support
)
{
// NOTE: Does NOT use PROTOTYPE_COMPLIANT_INDEXING.
//#if PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

//    clock_t t_0, t_1;
//    double secs_per_clock = 1.0 / CLOCKS_PER_SEC;

    if ((projection0->depth != IPL_DEPTH_32F) && (projection0->depth != IPL_DEPTH_64F))
        ABORT("Expected projection0->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but projection0->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection0->depth);

    if ((projection->depth != IPL_DEPTH_32F) && (projection->depth != IPL_DEPTH_64F))
        ABORT("Expected projection->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but projection->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection->depth);

    cvSetZero(projection);

    pixel_type *projection0_data = (pixel_type *) projection0->imageData;
    int n0_x = projection0->width;
    int n0_y = projection0->height;

    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.");

    real_type center0_x = support0->x + (support0->width / 2.0);
    real_type center0_y = support0->y + (support0->height / 2.0);

    real_type n0_x_min = support0->x;
    real_type n0_x_max = support0->x + support0->width;
    real_type n0_y_min = support0->y;
    real_type n0_y_max = support0->y + support0->height;

    real_type support0_width = support0->width;
    real_type support0_height = support0->height;

// PRIMARY DEBUG:     Print("n0_(xy) = (%i, %i), center0_(xy) = (%.15e, %.15e)\n", n0_x, n0_y, center0_x, center0_y); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("n0_(xy)_min = (%.15e, %.15e), n0_(xy)_max = (%.15e, %.15e)\n", n0_x_min, n0_y_min, n0_x_max, n0_y_max); // PRIMARY DEBUG.

    pixel_type *projection_data = (pixel_type *) projection->imageData;
    int n_x = projection->width;
    int n_y = projection->height;

    pixel_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.");

    real_type center_x = support->x + (support->width / 2.0);
    real_type center_y = support->y + (support->height / 2.0);

    real_type n_x_min = support->x;
    real_type n_x_max = support->x + support->width;
    real_type n_y_min = support->y;
    real_type n_y_max = support->y + support->height;

    real_type support_width = (real_type) INT_FROM_REAL_IMAGE_EXTENT(support->width);
    real_type support_height = (real_type) INT_FROM_REAL_IMAGE_EXTENT(support->height);

// PRIMARY DEBUG:     Print("n_(xy) = (%i, %i), center_(xy) = (%.15e, %.15e)\n", n_x, n_y, center_x, center_y); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("n_(xy)_min = (%.15e, %.15e), n_(xy)_max = (%.15e, %.15e)\n", n_x_min, n_y_min, n_x_max, n_y_max); // PRIMARY DEBUG.

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
    //if (((real_type) n0_x < support0_width) || ((real_type) n0_y < support0_height)) 
    //    ABORT("support0 extends beyond projection0.\n");
    
    //if (((real_type) n_x < support_width) || ((real_type) n_y < support_height))
    //    ABORT("support extends beyond projection.\n");

    int n0_x_support_min = -projection0_center_ROI_CvRect_global->x;
    int n0_x_support_max = projection0_global->width - projection0_center_ROI_CvRect_global->x;
    int n0_y_support_min = -projection0_center_ROI_CvRect_global->y;
    int n0_y_support_max = projection0_global->height - projection0_center_ROI_CvRect_global->y;

// PRIMARY DEBUG:     Print("n0_(xy)_support_min = (%i, %i), n0_(xy)_support_max = (%i, %i)\n", n0_x_support_min, n0_y_support_min, n0_x_support_max, n0_y_support_max); // PRIMARY DEBUG.

    for (int x0 = n0_x_support_min; x0 < n0_x_support_max; ++x0)
        for (int y0 = n0_y_support_min; y0 < n0_y_support_max; ++y0)
        {
            real_type v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0, y0);
            PUTVAL_2D_PIX(projection0_data, n0_x_ws, x0, y0, v0);
            //real_type index = (real_type) (x0 + (n0_x * y0));
            //PUTVAL_2D_PIX(projection0_data, n0_x, x0, y0, index);
        }

    //int n_x_support_min = (int) floor(support->x);
    //int n_x_support_max = INT_FROM_REAL_IMAGE_EXTENT(support->x + support->width);
    //int n_y_support_min = (int) floor(support->y);
    //int n_y_support_max = INT_FROM_REAL_IMAGE_EXTENT(support->y + support->height);
#endif

    real_type step = 1.0 / local_mag;
    //step = 1.0;  Print("DEBUG: WARNING: step is hardcoded to %.15e.\n", step);

// PRIMARY DEBUG:     Print("step = %.15e\n", step); // PRIMARY DEBUG.
//    exit(0);

    real_type sample_factor = step * step;
    //real_type sample_factor = 1.0; Print("DEBUG: WARNING: sample_factor is hardcoded to %.15e.\n", sample_factor);

    // NOTE: Use the inverse rotation matrix.
    real_type rot_matrix00 = rot_matrix_2x2[0][0];
    real_type rot_matrix01 = rot_matrix_2x2[1][0];
    real_type rot_matrix10 = rot_matrix_2x2[0][1];
    real_type rot_matrix11 = rot_matrix_2x2[1][1];

    real_type x0_init = 0.0; // NOTE: Used as an accumulator.
    real_type y0_init = 0.0; // NOTE: Used as an accumulator.

    real_type x_init = n_x_min;
    real_type x_init_rot_x = rot_matrix00 * x_init;
    real_type x_init_rot_y = rot_matrix10 * x_init;
    x0_init += x_init_rot_x;
    y0_init += x_init_rot_y;

    real_type y_init = n_y_min;
    real_type y_init_rot_x = rot_matrix01 * y_init;
    real_type y_init_rot_y = rot_matrix11 * y_init;
    x0_init += y_init_rot_x;
    y0_init += y_init_rot_y;

    x0_init -= n0_x_min;
    y0_init -= n0_y_min;

    real_type x0_del_x = rot_matrix00 * step;
    real_type y0_del_x = rot_matrix10 * step;

    real_type x0_del_y = rot_matrix01 * step;
    real_type y0_del_y = rot_matrix11 * step;

    real_type y0 = y0_init;

    for (real_type y = 0.0; y < support_height; y += step) // NOTE: In image coordinates.
    {
        int y_floor = FLOOR_INNER_LOOP(y);
        real_type y_alpha = y - (real_type) y_floor;

        real_type x0 = x0_init;

        for (real_type x = 0.0; x < support_width; x += step) // NOTE: In image coordinates.
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

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
            if ((x0_floor >= n0_x_support_min) && (x0_floor < n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor < n0_y_support_max))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor) * (1.0 - x0_alpha) * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
            }
            else
            {
                Print("(%i, %i)\n", x_floor, y_floor);
                Print("(%i, %i)\n", x0_floor, y0_floor);
                ABORT("Source projection bounds violated (0, 0).\n"); 
            }
#endif

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 < n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor < n0_y_support_max))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor) * x0_alpha * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
            }
            else
            {
                Print("(%i, %i)\n", x_floor, y_floor);
                Print("(%i, %i)\n", x0_floor + 1, y0_floor);
                ABORT("Source projection bounds violated (1, 0).\n"); 
            }
#endif

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
            if ((x0_floor >= n0_x_support_min) && (x0_floor < n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 < n0_y_support_max))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1) * (1.0 - x0_alpha) * y0_alpha;
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor + 1, v0);
#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
            }
            else
            {
                Print("(%i, %i)\n", x_floor, y_floor);
                Print("(%i, %i)\n", x0_floor, y0_floor + 1);
                ABORT("Source projection bounds violated (0, 1).\n");
            }
#endif

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 < n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 < n0_y_support_max))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1) * x0_alpha * y0_alpha;
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor + 1, v0);
#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
            }
            else
            {
                Print("(%i, %i)\n", x_floor, y_floor);
                Print("(%i, %i)\n", x0_floor + 1, y0_floor + 1);
                ABORT("Source projection bounds violated (1, 1).\n");
            }
#endif

            v *= sample_factor;

//            Print("v = %.15e\n", v);

            //ABORT("Abort!\n");

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Boundary test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            PUTVAL_2D_PIX(projection_data, n_x_ws, x_floor, y_floor, v * (1.0 - x_alpha) * (1.0 - y_alpha) * sample_factor);
//            PUTVAL_2D_PIX(projection_data, n_x_ws, x_floor + 1, y_floor, v * x_alpha * (1.0 - y_alpha) * sample_factor);
//            PUTVAL_2D_PIX(projection_data, n_x_ws, x_floor, y_floor + 1, v * (1.0 - x_alpha) * y_alpha * sample_factor);
//            PUTVAL_2D_PIX(projection_data, n_x_ws, x_floor + 1, y_floor + 1, v * x_alpha * y_alpha * sample_factor);
            ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor, y_floor, v * (1.0 - x_alpha) * (1.0 - y_alpha));
            ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor + 1, y_floor, v * x_alpha * (1.0 - y_alpha));
            ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor, y_floor + 1, v * (1.0 - x_alpha) * y_alpha);
            ADDVAL_2D_PIX(projection_data, n_x_ws, x_floor + 1, y_floor + 1, v * x_alpha * y_alpha);

            // END: Geometry vs. pixel intensity. }

            x0 += x0_del_x;
            y0 += y0_del_x;
        }

        x0_init += x0_del_y;
        y0_init += y0_del_y;
        y0 = y0_init;
    }
//    exit(0);
}
// END: TxBR 2.0 and TxBR 3.0 code. }

// NOTE: This pulls pixel values from the original projection to the new one.
// NOTE: Both dimensions use finite difference scheme.
void projection_rotate_2D_prototype
(
    IplImage *projection0,
    real_type support0[4],
    real_type local_mag,
    real_type rot_matrix_2x2[2][2],
    IplImage *projection,
    real_type support[4]
)
{
// NOTE: Does NOT use PROTOTYPE_COMPLIANT_INDEXING.
//#if PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

//    clock_t t_0, t_1;
//    double secs_per_clock = 1.0 / CLOCKS_PER_SEC;

    if ((projection0->depth != IPL_DEPTH_32F) && (projection0->depth != IPL_DEPTH_64F))
        ABORT("Expected projection0->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but projection0->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection0->depth);

    if ((projection->depth != IPL_DEPTH_32F) && (projection->depth != IPL_DEPTH_64F))
        ABORT("Expected projection->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but projection->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection->depth);

    cvSetZero(projection);

    pixel_type *projection0_data = (pixel_type *) projection0->imageData;
    int n0_x = projection0->width;
    int n0_y = projection0->height;

    pixel_type n0_x_ws_pixel = ((pixel_type) projection0->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n0_x_ws = (int) floor(n0_x_ws_pixel);
    if (n0_x_ws_pixel - (pixel_type) n0_x_ws != 0.0)
        ABORT("(n0_x_ws_pixel - n0_x_ws != 0.0.");

    real_type center0_x = support0[0];
    real_type center0_y = support0[1];

    //real_type n0_x_min = 0.0;
    //real_type n0_x_max = support0[2];
    //real_type n0_y_min = 0.0;
    //real_type n0_y_max = support0[3];

/*
    // NOTE: Uncomment if not using ROI.
    real_type center0_x = ((real_type) n0_x - 1.0) / 2.0;
    real_type center0_y = ((real_type) n0_y - 1.0) / 2.0;

    real_type center0_x_support = support0[0];
    real_type center0_y_support = support0[1];
    real_type n0_x_min = center0_x - center0_x_support;
    real_type n0_x_max = center0_x + center0_x_support;
    real_type n0_y_min = center0_y - center0_y_support;
    real_type n0_y_max = center0_y + center0_y_support;

    Print("n0_(xy) = (%i, %i), center0_(xy) = (%.15e, %.15e)\n", n0_x, n0_y, center0_x, center0_y);
    Print("center0_(xy)_support = (%.15e, %.15e)\n", center0_x_support, center0_y_support);
    Print("n0_(xy)_min = (%.15e, %.15e), n0_(xy)_max = (%.15e, %.15e)\n", n0_x_min, n0_y_min, n0_x_max, n0_y_max);

    center0_x = center0_x_support;
    center0_y = center0_y_support;
*/

    pixel_type *projection_data = (pixel_type *) projection->imageData;
    //int n_x = projection->width;
    //int n_y = projection->height;

    pixel_type n_x_ws_pixel = ((pixel_type) projection->widthStep) / ((pixel_type) sizeof(pixel_type));
    int n_x_ws = (int) floor(n_x_ws_pixel);
    if (n_x_ws_pixel - (pixel_type) n_x_ws != 0.0)
        ABORT("(n_x_ws_pixel - n_x_ws != 0.0.");

    real_type center_x = support[0];
    real_type center_y = support[1];

    real_type n_x_min = 0.0;
    real_type n_x_max = support[2];
    real_type n_y_min = 0.0;
    real_type n_y_max = support[3];

/*
    // NOTE: Uncomment if not using ROI.
    real_type center_x = ((real_type) n_x - 1.0) / 2.0;
    real_type center_y = ((real_type) n_y - 1.0) / 2.0;

    real_type center_x_support = support[0];
    real_type center_y_support = support[1];
    real_type n_x_min = center_x - center_x_support;
    real_type n_x_max = center_x + center_x_support;
    real_type n_y_min = center_y - center_y_support;
    real_type n_y_max = center_y + center_y_support;

    Print("n_(xy) = (%i, %i), center_(xy) = (%.15e, %.15e)\n", n_x, n_y, center_x, center_y);
    Print("center_(xy)_support = (%.15e, %.15e)\n", center_x_support, center_y_support);
    Print("n_(xy)_min = (%.15e, %.15e), n_(xy)_max = (%.15e, %.15e)\n", n_x_min, n_y_min, n_x_max, n_y_max);

    center_x = center_x_support;
    center_y = center_y_support;
*/

/*
    if ((n0_x_min < 0.0) || (n0_y_min < 0.0) || 
        ((real_type) n0_x < n0_x_max) || ((real_type) n0_y < n0_y_max) || 
        (n_x_min < 0.0) || (n_y_min < 0.0) || 
        ((real_type) n_x < n_x_max) || ((real_type) n_y < n_y_max))
        ABORT("Support extends beyond projection.\n"); 
*/

/*
    real_type step;
    if (test_rotation_matrix_2x2_for_subsampling(rot_matrix_2x2))
        step = 0.25;
    else
        step = 1.0;
*/

    real_type step = 1.0 / local_mag;
    //step = 1.0;  Print("DEBUG: WARNING: step is hardcoded to %.15e.\n", step);

// PRIMARY DEBUG:     Print("step = %.15e\n", step); // PRIMARY DEBUG.
//    exit(0);

    //real_type sample_factor = step; // MATLAB prototype uses step for number of samples per pixel.
    real_type sample_factor = step * step;
    //real_type sample_factor = 1.0; Print("DEBUG: WARNING: sample_factor is hardcoded to %.15e.\n", sample_factor);

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

    real_type y0 = y0_init;

    for (real_type y = n_y_min; y <= n_y_max; y += step) // NOTE: In image coordinates.
    {
        int y_floor = FLOOR_INNER_LOOP(y);
        real_type y_alpha = y - (real_type) y_floor;

        real_type x0 = x0_init;

        for (real_type x = n_x_min; x <= n_x_max; x += step) // NOTE: In image coordinates.
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

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor) * (1.0 - x0_alpha) * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
            }
            else
            {
                Print("(%i, %i)\n", x_floor, y_floor);
                Print("(%i, %i)\n", x0_floor, y0_floor);
                ABORT("Source projection bounds violated (0, 0).\n"); 
            }
#endif

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor >= n0_y_support_min) && (y0_floor <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor >= 0) && (y0_floor < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor) * x0_alpha * (1.0 - y0_alpha);
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor, v0);
#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
            }
            else
            {
                Print("(%i, %i)\n", x_floor, y_floor);
                Print("(%i, %i)\n", x0_floor + 1, y0_floor);
                ABORT("Source projection bounds violated (1, 0).\n"); 
            }
#endif

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
//            if ((x0_floor >= n0_x_support_min) && (x0_floor <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor >= 0) && (x0_floor < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor, y0_floor + 1) * (1.0 - x0_alpha) * y0_alpha;
                v += v0;
#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
//                Print("(%i, %i) = %.15e\n", x0_floor, y0_floor + 1, v0);
            }
            else
            {
                Print("(%i, %i)\n", x_floor, y_floor);
                Print("(%i, %i)\n", x0_floor, y0_floor + 1);
                ABORT("Source projection bounds violated (0, 1).\n");
            }
#endif

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
//            if ((x0_floor + 1 >= n0_x_support_min) && (x0_floor + 1 <= n0_x_support_max) && (y0_floor + 1 >= n0_y_support_min) && (y0_floor + 1 <= n0_y_support_max))
            if ((x0_floor + 1 >= 0) && (x0_floor + 1 < n0_x) && (y0_floor + 1 >= 0) && (y0_floor + 1 < n0_y))
            {
#endif
                v0 = INDEX_2D_PIX(projection0_data, n0_x_ws, x0_floor + 1, y0_floor + 1) * x0_alpha * y0_alpha;
                v += v0;
//                Print("(%i, %i) = %.15e\n", x0_floor + 1, y0_floor + 1, v0);
#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
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

        x0_init += x0_del_y;
        y0_init += y0_del_y;
        y0 = y0_init;
    }
//    exit(0);
}

#elif defined PROJECTION_ROTATE_2D_IPLIMAGE_PUSH

// NOTE: This pushes pixel values from the original projection to the new one.  This is fine with rotations.
void projection_rotate_2D_prototype
(
    IplImage *projection0,
    real_type support0[4],
    real_type rot_matrix_2x2[2][2],
    IplImage *projection,
    real_type support[4]
)
{
//    clock_t t_0, t_1;
//    double secs_per_clock = 1.0 / CLOCKS_PER_SEC;

    if ((projection0->depth != IPL_DEPTH_32F) && (projection0->depth != IPL_DEPTH_64F))
        ABORT("Expected projection0->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but projection0->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection0->depth);

    if ((projection->depth != IPL_DEPTH_32F) && (projection->depth != IPL_DEPTH_64F))
        ABORT("Expected projection->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but projection->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, projection->depth);

    cvSetZero(projection);

    pixel_type *projection0_data = (pixel_type *) projection0->imageData;
    int n0_x = projection0->width; // This is required because the original projection has been padded.
    int n0_y = projection0->height; // Unnecessary.

    real_type center0_x = support0[0];
    real_type center0_y = support0[1];
#if PROTOTYPE_COMPLIANT_INDEXING
    real_type n0_x_real = round(support0[2]);
    real_type n0_y_real = round(support0[3]);
#else
    real_type n0_x_real = round(support0[2]);
    real_type n0_y_real = round(support0[3]);
    //real_type n0_x_real = round(support0[2] + 1.0);
    //real_type n0_y_real = round(support0[3] + 1.0);
#endif

// PRIMARY DEBUG:     Print("n0_x = %i, n0_y = %i, center0_x = %.15e, center0_y = %.15e, n0_x_real = %.15e, n0_y_real = %.15e\n", n0_x, n0_y, center0_x, center0_y, n0_x_real, n0_y_real); // PRIMARY DEBUG.

    pixel_type *projection_data = (pixel_type *) projection->imageData;
    int n_x = projection->width; // This is required because the transformed projection has been padded.
    int n_y = projection->height; // Unnecessary.

    real_type center_x = support[0];
    real_type center_y = support[1];
#if PROTOTYPE_COMPLIANT_INDEXING
    real_type n_x_real = round(support[2]);
    real_type n_y_real = round(support[3]);
#else
    real_type n_x_real = round(support[2] + 1.0);
    real_type n_y_real = round(support[3] + 1.0);
#endif

#if PROTOTYPE_COMPLIANT_INDEXING
    int n_x_support = (int) n_x_real - 1;
    int n_y_support = (int) n_y_real - 1;
//    // NOTE: Conforms to MATLAB prototype by not calculating topmost row and rightmost column.
//    int n_x_support = (int) n_x_real - 1 - 1;
//    int n_y_support = (int) n_y_real - 1 - 1;
#else
    // NOTE: Calculates topmost row and rightmost column.
    int n_x_support = (int) n_x_real - 1;
    int n_y_support = (int) n_y_real - 1;
#endif

// PRIMARY DEBUG:     Print("n_x = %i, n_y = %i, center_x = %.15e, center_y = %.15e, n_x_real = %.15e, n_y_real = %.15e\n", n_x, n_y, center_x, center_y, n_x_real, n_y_real); // PRIMARY DEBUG.

    real_type step;
    if (test_rotation_matrix_2x2_for_subsampling(rot_matrix_2x2))
        step = 0.25;
    else
        step = 1.0;

    //step = 0.25;  Print("DEBUG: WARNING: step is hardcoded to %.15e.\n", step);

// PRIMARY DEBUG:     Print("step = %.15e\n", step); // PRIMARY DEBUG.
//    exit(0);

    //real_type sample_factor = step; // MATLAB prototype uses step for number of samples per pixel.
    real_type sample_factor = step * step;
    //real_type sample_factor = 1.0 / (step * step); Print("DEBUG: WARNING: sample_factor is hardcoded to %.15e.\n", sample_factor);

    real_type rot_matrix00 = rot_matrix[0][0];
    real_type rot_matrix01 = rot_matrix[0][1];
    real_type rot_matrix10 = rot_matrix[1][0];
    real_type rot_matrix11 = rot_matrix[1][1];

#if PROTOTYPE_COMPLIANT_INDEXING
    for (real_type y0 = 1.0; y0 <= n0_y_real; y0 += step)
#else
    for (real_type y0 = 0.0; y0 <= n0_y_real; y0 += step)
#endif
    {
        int y0_floor = FLOOR_INNER_LOOP(y0);
        real_type y0_alpha = y0 - (real_type) y0_floor;
#if PROTOTYPE_COMPLIANT_INDEXING
        --y0_floor;
#endif
        real_type y0_trans_y = y0 - center0_y;
        real_type y0_trans_y_rot_x = rot_matrix01 * y0_trans_y;
        real_type y0_trans_y_rot_y = rot_matrix11 * y0_trans_y;
        real_type y0_trans_y_rot_x_trans_x = y0_trans_y_rot_x + center_x;
        real_type y0_trans_y_rot_y_trans_y = y0_trans_y_rot_y + center_y;

#if PROTOTYPE_COMPLIANT_INDEXING
        for (real_type x0 = 1.0; x0 <= n0_x_real; x0 += step)
#else
        for (real_type x0 = 0.0; x0 <= n0_x_real; x0 += step)
#endif
        {
//            t_0 = clock();

            int x0_floor = FLOOR_INNER_LOOP(x0);

            real_type x0_alpha = x0 - (real_type) x0_floor;

#if PROTOTYPE_COMPLIANT_INDEXING
            --x0_floor;
#endif

//            Print("x0_floor = %i\n", x0_floor);
//            Print("y0_floor = %i\n", y0_floor);

//            Print("x0_alpha = %.15e\n", x0_alpha);
//            Print("y0_alpha = %.15e\n", y0_alpha);

            // BEGIN:: Geometry vs. pixel intensity. {

            // In order to avoid yet another boundary check, we have padded projection0.
            real_type v0 = INTERPOLATE_2D_LL_PIX(projection0_data, n0_x, x0_floor, y0_floor, x0_alpha, y0_alpha) * sample_factor;
//            real_type v0 = INDEX_2D_PIX(projection0_data, n0_x, x0_floor, y0_floor) * sample_factor;

            //Print("v0 = %.15e\n", v0);

            real_type x0_trans = x0 - center0_x;

            real_type x = y0_trans_y_rot_x_trans_x + (rot_matrix00 * x0_trans);
            real_type y = y0_trans_y_rot_y_trans_y + (rot_matrix10 * x0_trans);

//            Print("x = %.15e\n", x);
//            Print("y = %.15e\n", y);

            int x_floor = FLOOR_INNER_LOOP(x);
            int y_floor = FLOOR_INNER_LOOP(y);

            real_type x_alpha = x - (real_type) x_floor;
            real_type y_alpha = y - (real_type) y_floor;

#if PROTOTYPE_COMPLIANT_INDEXING
            --x_floor;
            --y_floor;
#endif

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Pre-boundary-test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

//            t_0 = clock();

            // NOTE: We can get away with this because we know both IplImages are mode IPL_DEPTH_64F.
            pixel_type v = 0.0;

            if ((x_floor >= 0) && (x_floor <= n_x_support) && (y_floor >= 0) && (y_floor <= n_y_support))
            {
                v = INDEX_2D(projection_data, n_x, x_floor, y_floor) + v0 * (1.0 - x_alpha) * (1.0 - y_alpha);
                PUTVAL_2D(projection_data, n_x, x_floor, y_floor, v);
//                Print("(%i, %i) = %.15e\n", x_floor, y_floor, v);
            }
//            else
//            {
//                Print("(%i, %i)\n", x0_floor, y0_floor);
//                Print("(%i, %i)\n", x_floor, y_floor);
//            }

            if ((x_floor + 1 >= 0) && (x_floor + 1 <= n_x_support) && (y_floor >= 0) && (y_floor <= n_y_support))
            {
                v = INDEX_2D(projection_data, n_x, x_floor + 1, y_floor) + v0 * x_alpha * (1.0 - y_alpha);
                PUTVAL_2D(projection_data, n_x, x_floor + 1, y_floor, v);
//                Print("(%i, %i) = %.15e\n", x_floor + 1, y_floor, v);
            }
//            else
//            {
//                Print("(%i, %i)\n", x0_floor + 1, y0_floor);
//                Print("(%i, %i)\n", x_floor + 1, y_floor);
//            }

            if ((x_floor >= 0) && (x_floor <= n_x_support) && (y_floor + 1 >= 0) && (y_floor + 1 <= n_y_support))
            {
                v = INDEX_2D(projection_data, n_x, x_floor, y_floor + 1) + v0 * (1.0 - x_alpha) * y_alpha;
                PUTVAL_2D(projection_data, n_x, x_floor, y_floor + 1, v);
//                Print("(%i, %i) = %.15e\n", x_floor, y_floor + 1, v);
            }
//            else
//            {
//                Print("(%i, %i)\n", x0_floor, y0_floor + 1);
//                Print("(%i, %i)\n", x_floor, y_floor + 1);
//            }

            if ((x_floor + 1 >= 0) && (x_floor + 1 <= n_x_support) && (y_floor + 1 >= 0) && (y_floor + 1 <= n_y_support))
            {
                v = INDEX_2D(projection_data, n_x, x_floor + 1, y_floor + 1) + v0 * x_alpha * y_alpha;
                PUTVAL_2D(projection_data, n_x, x_floor + 1, y_floor + 1, v);
//                Print("(%i, %i) = %.15e\n", x_floor + 1, y_floor + 1, v);
            }
//            else
//            {
//                Print("(%i, %i)\n", x0_floor + 1, y0_floor + 1);
//                Print("(%i, %i)\n", x_floor + 1, y_floor + 1);
//            }

//            t_1 = clock();
//            Print("%i %i\n", t_1, t_0);
//            Print("Boundary test delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

            // END:: Geometry vs. pixel intensity. }
        }
    }
//    exit(0);
}

#endif

////////////////////////////////////////////////////////////////////////////////
// END: 2D rotation. }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BEGIN: projection_series_filter_1D(). {
////////////////////////////////////////////////////////////////////////////////

// BEGIN: TxBR 2.0 and TxBR 3.0 code. {

void projection_series_remap_2D_filter_1D_inv_remap_2D
(
    filesystem_params_handle filesystem_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    remap_2D_params_handle remap_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
)
{
// NOTE: Does NOT use PROTOTYPE_COMPLIANT_INDEXING.
//#if PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

    // BEGIN: Check parameter structs.
    if (filesystem_params == NULL)
        ABORT("filesystem_params == NULL.\n");

    if (symmetrize_2D_params == NULL)
        symmetrize_2D_params = create_symmetrize_2D_params();

    if (remap_2D_params == NULL)
        ABORT("remap_2D_params == NULL.\n");

    if (filter_1D_params == NULL)
        filter_1D_params = create_filter_1D_params();

    if (edgetaper_2D_params == NULL)
        edgetaper_2D_params = create_edgetaper_2D_params();
    // END: Check parameter structs.

// PRIMARY DEBUG:     Print("filesystem_params->filepath_in = \"%s\"\n", filesystem_params->filepath_in); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("filesystem_params->filepath_out = \"%s\"\n", filesystem_params->filepath_out); // PRIMARY DEBUG.

    clock_t t_0_total, t_1_total;
    clock_t t_0, t_1;

    double secs_per_clock = 1.0 / (double) CLOCKS_PER_SEC;

    t_0_total = clock();

    // BEGIN: Open input MRC file.
    MrcHeader file_header_in;

    FILE *file_in = NULL;

    openMRCFile_general(filesystem_params->filepath_in, &file_header_in, &file_in);
    // END: Open input MRC file.

    int n0_x = file_header_in.nx;
    int n0_y = file_header_in.ny;
    int n_projections = file_header_in.nz;
    //int n_projections = 1; WARN("n_projections hardcoded to %i.\n", n_projections);
    int mode = file_header_in.mode;

    // BEGIN: Create and open output MRC file.
    // BEGIN: Create output MRC file.
    if (mode != SLICE_MODE_FLOAT)
        Print("WARNING: Input mode == %i but using mode SLICE_MODE_FLOAT == %i for output!\n", mode, SLICE_MODE_FLOAT);

    MrcHeader file_header_out;

    FILE *file_out = NULL;

    createNewMRCFile(filesystem_params->filepath_out, n0_x, n0_y, filesystem_params->n_projections_indexed, /*mode*/ SLICE_MODE_FLOAT);
    // END: Create output MRC file.

    // BEGIN: Open output MRC file.
    openMRCFile_general(filesystem_params->filepath_out, &file_header_out, &file_out);
    // END: Open output MRC file.
    // END: Create and open output MRC file.

    remap_2D_data_type *remap_2D_data = create_remap_2D_data(remap_2D_params);

    //for (int i_projection = 0; i_projection < n_projections; i_projection++)
    for (int i_index = 0; i_index < filesystem_params->n_projections_indexed; ++i_index)
    {
        int i_projection = filesystem_params->projection_indices[i_index];
        Print("i_projection = %i\n", i_projection);

        if ((i_projection < 0) || (i_projection >= n_projections))
            ABORT("i_projection == %i, but n_projections == %i.\n", i_projection, n_projections);

        t_0 = clock();

        set_map_2D_current(remap_2D_data, i_projection);

        int n_neg_row_pad = 1;
        int n_pos_row_pad = 1;
        int n_neg_col_pad = 1;
        int n_pos_col_pad = 1;

        // NOTE: Treat pixels as points.
        rect_2DR_type support0 = { 0.0, 0.0, (real_type) (n0_x - 1), (real_type) (n0_y - 1) };
        rect_2DR_type support = support0;

        Print("support = support_remap_2D(support0)\n");
        support_remap_2D(remap_2D_data, &support);

        int n_x = INT_FROM_REAL_IMAGE_EXTENT(support.width) + 1;
        int n_y = INT_FROM_REAL_IMAGE_EXTENT(support.height) + 1;

        int n0_push_extend = 0;
        if (symmetrize_2D_params->symmetrize_2D_flag)
        {
            int n0_push_extend_min = 5;
            real_type diff_width = fabs(support.width - support0.width);
            real_type diff_height = fabs(support.height - support0.height);
            real_type diff_max = fmax(diff_width, diff_height);
            n0_push_extend = (int) fmax(ceil(diff_max), (real_type) n0_push_extend_min);
            n0_push_extend += (int) ceil(fmax(support0.width, support0.height) / 100.0); // NOTE: Extend by 1% of maximum extent.
// PRIMARY DEBUG:             Print("n0_push_extend = %i\n", n0_push_extend); // PRIMARY DEBUG.
        }

// PRIMARY DEBUG:         Print("n0_push_extend = %i\n", n0_push_extend); // PRIMARY DEBUG.

        int n_push_extend = n0_push_extend;

// PRIMARY DEBUG:         Print("n_push_extend = %i\n", n_push_extend); // PRIMARY DEBUG.

        int n_x_push_extended = n_x + (2 * n_push_extend);
        int n_y_push_extended = n_y + (2 * n_push_extend);

        int n_x_push_extended_padded = n_x_push_extended + n_neg_row_pad + n_pos_row_pad;
        int n_y_push_extended_padded = n_y_push_extended + n_neg_col_pad + n_pos_col_pad;

        rect_2DR_type support_push_extended = 
        { 
            support.x + ((real_type) -n_push_extend),
            support.y + ((real_type) -n_push_extend),
            (real_type) (INT_FROM_REAL_IMAGE_EXTENT(support.width) + (2 * n_push_extend)),
            (real_type) (INT_FROM_REAL_IMAGE_EXTENT(support.height) + (2 * n_push_extend))
        };

        int n0_pull_extend = 0;
        if (n0_push_extend > 0)
        {
            n0_pull_extend = 2;
            if (n0_push_extend <= n0_pull_extend)
                ABORT("n0_push_extend == %i <= %i == n0_pull_extend\n", n0_push_extend, n0_pull_extend);
        }

// PRIMARY DEBUG:         Print("n0_pull_extend = %i\n", n0_pull_extend); // PRIMARY DEBUG.

        int n0_x_push_extended = n0_x + (2 * n0_push_extend);
        int n0_y_push_extended = n0_y + (2 * n0_push_extend);

        int n0_x_push_extended_padded = n0_x_push_extended + n_neg_row_pad + n_pos_row_pad;
        int n0_y_push_extended_padded = n0_y_push_extended + n_neg_col_pad + n_neg_col_pad;

        // Ensure that the parities of the extents of the push projection0 match the parities of the extents of the source image.
        if ((n0_x % 2) != (n0_x_push_extended_padded % 2)) ++n0_x_push_extended_padded;
        if ((n0_y % 2) != (n0_y_push_extended_padded % 2)) ++n0_y_push_extended_padded;

        // NOTE: Treat pixels as points.
        rect_2DR_type support0_push_extended = { (real_type) -n0_push_extend, (real_type) -n0_push_extend, (real_type) (n0_x_push_extended - 1), (real_type) (n0_y_push_extended - 1) };

        int n0_x_pull_extended = n0_x + (2 * n0_pull_extend);
        int n0_y_pull_extended = n0_y + (2 * n0_pull_extend);

        //int n0_x_pull_extended_padded = n0_x_pull_extended + n_neg_row_pad + n_pos_row_pad;
        //int n0_y_pull_extended_padded = n0_y_pull_extended + n_neg_col_pad + n_neg_col_pad;

        // Ensure that the parities of the extents of the pull projection0 match the parities of the extents of the source image.
        //if ((n0_x % 2) != (n0_x_pull_extended_padded % 2)) ++n0_x_pull_extended_padded;
        //if ((n0_y % 2) != (n0_y_pull_extended_padded % 2)) ++n0_y_pull_extended_padded;

        // NOTE: Treat pixels as points.
        rect_2DR_type support0_pull_extended = { (real_type) -n0_pull_extend, (real_type) -n0_pull_extend, (real_type) (n0_x_pull_extended - 1), (real_type) (n0_y_pull_extended - 1) };

        IplImage *projection0 = create_IplImage_SetZero_malloc(n0_x_push_extended_padded, n0_y_push_extended_padded, IPL_DEPTH_REAL);
        IplImage *projection = create_IplImage_SetZero_malloc(n_x_push_extended_padded, n_y_push_extended_padded, IPL_DEPTH_REAL);
        IplImage *sums_image = create_IplImage_SetZero_malloc(n_x_push_extended, n_y_push_extended, IPL_DEPTH_REAL);
        IplImage *hits_image = create_IplImage_SetZero_malloc(n_x_push_extended, n_y_push_extended, IPL_DEPTH_REAL);

        center_ROI_params_type center_params;

        center_params = calc_center_ROI_params(n0_x_push_extended_padded, n0_y_push_extended_padded, n0_x, n0_y);
        CvRect projection0_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                      center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection0_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                                      projection0,
                                                      projection0_center_ROI_CvRect);

        center_params = calc_center_ROI_params(n0_x_push_extended_padded, n0_y_push_extended_padded, n0_x_push_extended, n0_y_push_extended);
        CvRect projection0_push_extended_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                                    center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection0_push_extended_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                                             projection0,
                                                             projection0_push_extended_center_ROI_CvRect);

        center_params = calc_center_ROI_params(n0_x_push_extended_padded, n0_y_push_extended_padded, n0_x_pull_extended, n0_y_pull_extended);
        CvRect projection0_pull_extended_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                                    center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection0_pull_extended_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                                             projection0,
                                                             projection0_pull_extended_center_ROI_CvRect);

        center_params = calc_center_ROI_params(n_x_push_extended_padded, n_y_push_extended_padded, n_x, n_y);
        CvRect projection_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                     center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                              projection,
                                              projection_center_ROI_CvRect);

        center_params = calc_center_ROI_params(n_x_push_extended_padded, n_y_push_extended_padded, n_x_push_extended, n_y_push_extended);
        CvRect projection_push_extended_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                                   center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection_push_extended_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                                            projection,
                                                            projection_push_extended_center_ROI_CvRect);

        // Reuse slice later to transfer data back to the MRC file.
        Islice *slice = NULL;
        slice = sliceReadMRC(&file_header_in, i_projection, 'z');
        if (slice == NULL)
            ABORT("Cannot read file: \"%s\".\n", filesystem_params->filepath_in);

        // NOTE: Output floating-point data.
        if (slice->mode != SLICE_MODE_FLOAT)
            sliceNewMode(slice, SLICE_MODE_FLOAT);

        copy_Islice_to_IplImage_center(slice, projection0_center_ROI);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0_center_ROI, "/home/akulo/filter_1D_images/projection0_center_ROI.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0_push_extended_center_ROI, "/home/akulo/filter_1D_images/projection0_push_extended_center_ROI.mrc"); // PRIMARY DEBUG.

        // Calculate filter.
        filter_1D_data_type *filter_data = create_filter_1D_data_n_xy(filter_1D_params, n_x, n_y);

        t_1 = clock();
// PRIMARY DEBUG:         Print("Initialization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        if (symmetrize_2D_params->symmetrize_2D_flag)
        {
            t_0 = clock();

            symmetrize_IplImage_values_CvRect(projection0, projection0_center_ROI_CvRect);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.symmetrized.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection0_center_ROI, "/home/akulo/filter_1D_images/projection0_center_ROI.symmetrized.mrc"); // PRIMARY DEBUG.

            t_1 = clock();
// PRIMARY DEBUG:             Print("Symmetrization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.
        }

        t_0 = clock();

        projection_remap_2D(
            projection0_push_extended_center_ROI,
            &support0_push_extended,
            remap_2D_data,
            sums_image,
            hits_image,
            projection_push_extended_center_ROI,
            &support_push_extended);

        //copy_Islice_to_IplImage(slice, projection); // DEBUG: No remap.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.remapped.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection_push_extended_center_ROI, "/home/akulo/filter_1D_images/projection_push_extended_center_ROI.remapped.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection_center_ROI, "/home/akulo/filter_1D_images/projection_center_ROI.remapped.mrc"); // PRIMARY DEBUG.

        t_0 = clock();

        // BEGIN: Extend rows left and right in image.
        //if (symmetrize_2D_params->symmetrize_2D_flag)
        //{
        //    t_0 = clock();
        //
        //    extend_IplImage_rows_analyze(projection_push_extended_center_ROI);
        //
        //    t_1 = clock();
// PRIMARY DEBUG:         //    Print("Symmetrization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.
        //
// PRIMARY DEBUG:         //    write_MRC_image_from_IplImage(projection_push_extended_center_ROI, "/home/akulo/filter_1D_images/projection_push_extended_center_ROI.remapped.symmetrized.mrc"); // PRIMARY DEBUG.
        //}
        // END: Extend rows left and right in image.

        projection_filter_1D_fs(projection_center_ROI, filter_data, symmetrize_2D_params);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.remapped.filtered.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection_push_extended_center_ROI, "/home/akulo/filter_1D_images/projection_push_extended_center_ROI.remapped.filtered.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection_center_ROI, "/home/akulo/filter_1D_images/projection_center_ROI.remapped.filtered.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Filtering delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        if (symmetrize_2D_params->symmetrize_2D_flag)
        {
            t_0 = clock();

            symmetrize_IplImage_values_CvRect(projection, projection_center_ROI_CvRect);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.remapped.filtered.symmetrized.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection_center_ROI, "/home/akulo/filter_1D_images/projection_center_ROI.remapped.filtered.symmetrized.mrc"); // PRIMARY DEBUG.

            t_1 = clock();
// PRIMARY DEBUG:             Print("Symmetrization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.
        }
        else
        {
            t_0 = clock();

            cvSetZero_IplImage_values_outside_CvRect(projection, projection_center_ROI_CvRect);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.remapped.filtered.cleared.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection_center_ROI, "/home/akulo/filter_1D_images/projection_center_ROI.remapped.filtered.cleared.mrc"); // PRIMARY DEBUG.

            t_1 = clock();
// PRIMARY DEBUG:             Print("Post-filter clearing delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.
        }

        t_0 = clock();

        cvSetZero(projection0);

        projection_inv_remap_2D(
            projection_push_extended_center_ROI,
            &support_push_extended,
            remap_2D_data,
            projection0_pull_extended_center_ROI,
            &support0_pull_extended);

        t_1 = clock();
// PRIMARY DEBUG:         Print("Restore delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        // BEGIN: Edgetapering.
        if (edgetaper_2D_params->width > 0)
        {
            ABORT("Edgetapering must occur after the inv. rotation but before the inv. transform.   Can't be done if rotation and transform are combined.\n");
            edgetaper_2D(projection0_pull_extended_center_ROI, edgetaper_2D_params->width);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection0_pull_extended_center_ROI, "/home/akulo/filter_1D_images/projection0_pull_extended_center_ROI.edgetapered.mrc"); // PRIMARY DEBUG.
        }
        // END: Edgetapering.

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.remapped.filtered.inv_remapped.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0_pull_extended_center_ROI, "/home/akulo/filter_1D_images/projection0_pull_extended_center_ROI.remapped.filtered.inv_remapped.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0_center_ROI, "/home/akulo/filter_1D_images/projection0_center_ROI.remapped.filtered.inv_remapped.mrc"); // PRIMARY DEBUG.
//        WARN("Returning...\n"); return;
        //ABORT("Abort!\n");

        copy_IplImage_to_Islice_center(projection0_center_ROI, slice);

        sliceWriteMRC(&file_header_out, slice, i_index, 'z');
        sliceFree(slice);

        // Free resources.
        cvReleaseImageHeader(&projection0_center_ROI);
        cvReleaseImageHeader(&projection0_push_extended_center_ROI);
        cvReleaseImageHeader(&projection0_pull_extended_center_ROI);
        cvReleaseImageHeader(&projection_center_ROI);
        cvReleaseImageHeader(&projection_push_extended_center_ROI);
        free_imageData_cvReleaseImageHeader(projection0);
        free_imageData_cvReleaseImageHeader(projection);
        free_imageData_cvReleaseImageHeader(sums_image);
        free_imageData_cvReleaseImageHeader(hits_image);

        filter_1D_data_release(filter_data);

        Print("Done processing projection, i_projection = %i.\n", i_projection);

#ifdef MEX
        mexEvalString("fprintf('Done processing projection.\\n');");
#endif
    }

    calc_mmm_of_MRC(&file_header_out);

    // Free resources.
    remap_2D_data_release(remap_2D_data);

    if (fclose(file_in))
        ABORT("Cannot close input file: \"%s\".\n", filesystem_params->filepath_in);

    if (fclose(file_out))
        ABORT("Cannot close output file \"%s\".\n", filesystem_params->filepath_out);

    t_1_total = clock();
 
    Print("Total duration: %.15e\n", secs_per_clock * ((double) t_1_total - (double) t_0_total));
}

/*
NOTE: Symmetrization unimplemented.
void projection_series_remap_2D_filter_1D_inv_remap_2D
(
    filesystem_params_handle filesystem_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    remap_2D_params_handle remap_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
)
{
// NOTE: Does NOT use PROTOTYPE_COMPLIANT_INDEXING.
//#if PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

    // BEGIN: Check parameter structs.
    if (filesystem_params == NULL)
        ABORT("filesystem_params == NULL.\n");

    if (symmetrize_2D_params == NULL)
        symmetrize_2D_params = create_symmetrize_2D_params();

    if (remap_2D_params == NULL)
        ABORT("remap_2D_params == NULL.\n");

    if (filter_1D_params == NULL)
        filter_1D_params = create_filter_1D_params();

    if (edgetaper_2D_params == NULL)
        edgetaper_2D_params = create_edgetaper_2D_params();
    // END: Check parameter structs.

// PRIMARY DEBUG:     Print("filesystem_params->filepath_in = \"%s\"\n", filesystem_params->filepath_in); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("filesystem_params->filepath_out = \"%s\"\n", filesystem_params->filepath_out); // PRIMARY DEBUG.

    clock_t t_0_total, t_1_total;
    clock_t t_0, t_1;

    double secs_per_clock = 1.0 / (double) CLOCKS_PER_SEC;

    t_0_total = clock();

    // BEGIN: Open input MRC file.
    MrcHeader file_header_in;

    FILE *file_in = NULL;

    openMRCFile_general(filesystem_params->filepath_in, &file_header_in, &file_in);
    // END: Open input MRC file.

    int n0_x = file_header_in.nx;
    int n0_y = file_header_in.ny;
    int n_projections = file_header_in.nz;
    //int n_projections = 1; WARN("n_projections hardcoded to %i.\n", n_projections);
    int mode = file_header_in.mode;

    // BEGIN: Create and open output MRC file.
    // BEGIN: Create output MRC file.
    if (mode != SLICE_MODE_FLOAT)
        Print("WARNING: Input mode == %i but using mode SLICE_MODE_FLOAT == %i for output!\n", mode, SLICE_MODE_FLOAT);

    MrcHeader file_header_out;

    FILE *file_out = NULL;

    createNewMRCFile(filesystem_params->filepath_out, n0_x, n0_y, filesystem_params->n_projections_indexed, /\*mode*\/ SLICE_MODE_FLOAT);
    // END: Create output MRC file.

    // BEGIN: Open output MRC file.
    openMRCFile_general(filesystem_params->filepath_out, &file_header_out, &file_out);
    // END: Open output MRC file.
    // END: Create and open output MRC file.

    remap_2D_data_type *remap_2D_data = create_remap_2D_data(remap_2D_params);

    //for (int i_projection = 0; i_projection < n_projections; i_projection++)
    for (int i_index = 0; i_index < filesystem_params->n_projections_indexed; ++i_index)
    {
        int i_projection = filesystem_params->projection_indices[i_index];
// PRIMARY DEBUG:         Print("i_projection = %i\n", i_projection); // PRIMARY DEBUG.

        if ((i_projection < 0) || (i_projection >= n_projections))
            ABORT("i_projection == %i, but n_projections == %i.\n", i_projection, n_projections);

        t_0 = clock();

        set_map_2D_current(remap_2D_data, i_projection);

        int row_pad = 1;
        int col_pad = 1;

        int n0_x_padded = n0_x + row_pad;
        int n0_y_padded = n0_y + col_pad;

        // NOTE: Treat pixels as points.
        rect_2DR_type support0 = { 0.0, 0.0, (real_type) (n0_x - 1), (real_type) (n0_y - 1) };
        rect_2DR_type support = support0;

        support_remap_2D(remap_2D_data, &support);

        int n_x = (int) floor(support.width) + 1;
        int n_y = (int) floor(support.height) + 1;

        int n_x_padded = n_x + row_pad;
        int n_y_padded = n_y + col_pad;

        // Calculate filter.
        filter_1D_data_type *filter_data = create_filter_1D_data_n_xy(filter_1D_params, n_x_padded, n_y_padded);

        IplImage *projection0 = create_IplImage_SetZero_malloc(n0_x_padded, n0_y_padded, IPL_DEPTH_REAL);
        IplImage *projection = create_IplImage_SetZero_malloc(n_x_padded, n_y_padded, IPL_DEPTH_REAL);
        IplImage *sums_image = create_IplImage_SetZero_malloc(n_x, n_y, IPL_DEPTH_REAL);
        IplImage *hits_image = create_IplImage_SetZero_malloc(n_x, n_y, IPL_DEPTH_REAL);

        IplImage *projection0_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                        projection0, cvRect(0, 0, n0_x, n0_y));

        IplImage *projection_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                       projection, cvRect(0, 0, n_x, n_y));

        // Reuse slice later to transfer data back to the MRC file.
        Islice *slice = NULL;
        slice = sliceReadMRC(&file_header_in, i_projection, 'z');
        if (slice == NULL)
            ABORT("Cannot read file: \"%s\".\n", filesystem_params->filepath_in);

        // NOTE: Output floating-point data.
        if (slice->mode != SLICE_MODE_FLOAT)
            sliceNewMode(slice, SLICE_MODE_FLOAT);

        copy_Islice_to_IplImage(slice, projection0);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0_ROI, "/home/akulo/filter_1D_images/projection0.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Initialization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        t_0 = clock();

        projection_remap_2D(
            projection0_ROI,
            &support0,
            remap_2D_data,
            sums_image,
            hits_image,
            projection_ROI,
            &support);

        //copy_Islice_to_IplImage(slice, projection); // DEBUG: No remap.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection_ROI, "/home/akulo/filter_1D_images/projection.remapped.mrc"); // PRIMARY DEBUG.

        t_0 = clock();

        // BEGIN: Extend rows left and right in image.
        if (symmetrize_2D_params->symmetrize_2D_flag)
        {
            ABORT("Symmetrization disabled in remapping code.\n");
            t_0 = clock();

            extend_IplImage_rows_analyze(projection_ROI);

            t_1 = clock();
// PRIMARY DEBUG:             Print("Symmetrization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection_ROI, "/home/akulo/filter_1D_images/projection.remapped.symmetrized.mrc"); // PRIMARY DEBUG.
        }
        // END: Extend rows left and right in image.

        projection_filter_1D_fs(projection, filter_data, symmetrize_2D_params);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.remapped.filtered.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Filtering delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        t_0 = clock();

        projection_inv_remap_2D(
            projection_ROI,
            &support,
            remap_2D_data,
            projection0_ROI,
            &support0);

        t_1 = clock();
// PRIMARY DEBUG:         Print("Restore delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        // BEGIN: Edgetapering.
        if (edgetaper_2D_params->width > 0)
        {
            ABORT("Edgetapering must occur after the inv. rotation but before the inv. transform.   Can't be done if rotation and transform are combined.\n");
            edgetaper_2D(projection0_ROI, edgetaper_2D_params->width);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.filtered.edgetapered.mrc"); // PRIMARY DEBUG.
        }
        // END: Edgetapering.

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0_ROI, "/home/akulo/filter_1D_images/projection0.remapped.filtered.inv_remapped.mrc"); // PRIMARY DEBUG.
//        WARN("Returning...\n"); return;
        //ABORT("Abort!\n");

        copy_IplImage_to_Islice(projection0, slice);

        sliceWriteMRC(&file_header_out, slice, i_index, 'z');
        sliceFree(slice);

        // Free resources.
        cvReleaseImageHeader(&projection0_ROI);
        cvReleaseImageHeader(&projection_ROI);
        free_imageData_cvReleaseImageHeader(projection0);
        free_imageData_cvReleaseImageHeader(projection);
        free_imageData_cvReleaseImageHeader(sums_image);
        free_imageData_cvReleaseImageHeader(hits_image);

        filter_1D_data_release(filter_data);

#ifdef MEX
        mexEvalString("fprintf('Done processing projection.\\n');");
#endif
    }

    calc_mmm_of_MRC(&file_header_out);

    // Free resources.
    remap_2D_data_release(remap_2D_data);

    if (fclose(file_in))
        ABORT("Cannot close input file: \"%s\".\n", filesystem_params->filepath_in);

    if (fclose(file_out))
        ABORT("Cannot close output file \"%s\".\n", filesystem_params->filepath_out);

    t_1_total = clock();
 
    Print("Total duration: %.15e\n", secs_per_clock * ((double) t_1_total - (double) t_0_total));
}
*/
// END: TxBR 2.0 and TxBR 3.0 code. }

#ifdef FILTER_1D_ENABLE_CUDA

#define PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA 1

void projection_series_transform_2D_CUDA_filter_1D_inv_transform_2D_CUDA
(
    filesystem_params_handle filesystem_params,
    general_CUDA_params_handle general_CUDA_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    remap_2D_params_handle remap_2D_params,
    rotate_2D_params_handle rotate_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
)
{
// NOTE: Uses PROTOTYPE_COMPLIANT_INDEXING.
//#if !PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("!PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

    if (PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA == 1)
        Print("PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA == 1\n");
    else
        Print("PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA == %d != 1\n", PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA);

    Print("filesystem_params->filepath_in = \"%s\"\n", filesystem_params->filepath_in);
    Print("filesystem_params->filepath_out = \"%s\"\n", filesystem_params->filepath_out);

    clock_t t_0_total, t_1_total;
    clock_t t_0, t_1;

    double secs_per_clock = 1.0 / (double) CLOCKS_PER_SEC;

    t_0_total = clock();

    // BEGIN: Open input MRC file.
    MrcHeader file_header_in;

    FILE *file_in = NULL;

    openMRCFile_general(filesystem_params->filepath_in, &file_header_in, &file_in);
    // END: Open input MRC file.

    int n0_x = file_header_in.nx;
    int n0_y = file_header_in.ny;
    int n_projections = file_header_in.nz;
    //int n_projections = 1; WARN("n_projections hardcoded to %i.\n", n_projections);
    int mode = file_header_in.mode;

    // BEGIN: Create and open output MRC file.
    // BEGIN: Create output MRC file.
    if (mode != SLICE_MODE_FLOAT)
        WARN("Input mode == %i but using mode SLICE_MODE_FLOAT == %i for output!\n", mode, SLICE_MODE_FLOAT);

    MrcHeader file_header_out;

    FILE *file_out = NULL;

    createNewMRCFile(filesystem_params->filepath_out, n0_x, n0_y, filesystem_params->n_projections_indexed, /*mode*/ SLICE_MODE_FLOAT);
    // END: Create output MRC file.

    // BEGIN: Open output MRC file.
    openMRCFile_general(filesystem_params->filepath_out, &file_header_out, &file_out);
    // END: Open output MRC file.
    // END: Create and open output MRC file.

    rotate_2D_data_type *rotate_2D_data = create_rotate_2D_data(rotate_2D_params, filesystem_params->n_projections_indexed);
    remap_2D_data_type *remap_2D_data = create_remap_2D_data(remap_2D_params);

    //for (int i_projection = 0; i_projection < n_projections; i_projection++)
    for (int i_index = 0; i_index < filesystem_params->n_projections_indexed; ++i_index)
    {
        int i_projection = filesystem_params->projection_indices[i_index];
        Print("i_projection = %i\n", i_projection);

        if ((i_projection < 0) || (i_projection >= n_projections))
            ABORT("i_projection == %i, but n_projections == %i.\n", i_projection, n_projections);

        t_0 = clock();

        real_type rot_matrix_2x2[2][2];
        real_type inv_rot_matrix_2x2[2][2];

        calc_rotation_matrices_current(rotate_2D_data, i_index, rot_matrix_2x2, inv_rot_matrix_2x2);
        set_map_2D_current(remap_2D_data, i_projection);

        // NOTE: Padding on the top and right of images.
        int row_pad = 1;
        int col_pad = 1;

        int n0_x_padded = n0_x + row_pad;
        int n0_y_padded = n0_y + col_pad;

        real_type support0[4];

//#if PROTOTYPE_COMPLIANT_INDEXING
        support0[0] = ((real_type) n0_x + 1.0) / 2.0; // y-coordinate of projection center
        support0[1] = ((real_type) n0_y + 1.0) / 2.0; // y-coordinate of projection center
        support0[2] = (real_type) n0_x; // max x-coordinate of projection
        support0[3] = (real_type) n0_y; // max y-coordinate of projection
//#else
//        support0[0] = ((real_type) n0_x - 1.0) / 2.0; // y-coordinate of projection center
//        support0[1] = ((real_type) n0_y - 1.0) / 2.0; // y-coordinate of projection center
//        support0[2] = (real_type) n0_x - 1.0; // max x-coordinate of projection
//        support0[3] = (real_type) n0_y - 1.0; // max y-coordinate of projection
//#endif

        real_type support[4];
        support[0] = support0[0]; // x-coordinate of transformed projection center
        support[1] = support0[1]; // y-coordinate of transformed projection center
        support[2] = support0[2]; // max x-coordinate of transformed projection
        support[3] = support0[3]; // max y-coordinate of transformed projection

#if PROTOTYPE_COMPLIANT_ROTATION
        // NOTE: I think the transform should be the same as that used to transform the original projection, but cf. the MATLAB prototype.
        calc_transformed_support(rotate_2D_params->rot_matrix_2x2, support);
#else
        calc_transformed_support(rotate_2D_params->inv_rot_matrix_2x2, support);
#endif

//#if PROTOTYPE_COMPLIANT_INDEXING
        int n_x = round(support[2]);
        int n_y = round(support[3]);
//#else
//        int n_x = round(support[2] + 1.0);
//        int n_y = round(support[3] + 1.0);
//#endif

        int n_x_padded = n_x + row_pad;
        int n_y_padded = n_y + col_pad;

        // Calculate filter.
        filter_1D_data_type *filter_data = create_filter_1D_data_n_xy(filter_1D_params, n_x_padded, n_y_padded);

#if PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA

        int deviceID = 1;
        projection_series_transform_2D_CUDA_filter_1D_inv_transform_2D_CUDA_init_GPU(deviceID);

        // Eventually may want to use page-locked memory.

#else // PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA

        // Eventually may want to use page-locked memory.

#endif // PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA

        IplImage *projection0 = create_IplImage_SetZero_malloc(n0_x_padded, n0_y_padded, IPL_DEPTH_REAL);
        IplImage *projection = create_IplImage_SetZero_malloc(n_x_padded, n_y_padded, IPL_DEPTH_REAL);
        IplImage *sums_image = create_IplImage_XXF_from_pixel_type_data(n_x, n_y, NULL); // NOTE: No ROI required.  No data allocation.
        IplImage *hits_image = create_IplImage_XXF_from_pixel_type_data(n_x, n_y, NULL); // NOTE: No ROI required.  No data allocation.

        IplImage *projection0_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                        projection0, cvRect(0, 0, n0_x, n0_y));

        IplImage *projection_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                       projection, cvRect(0, 0, n_x, n_y));

        // NOTE: No dynamic memory.
        IplImage_base projection0_base;
        IplImage_base projection_base;
        IplImage_base sums_image_base; // No data allocation.
        IplImage_base hits_image_base; // No data allocation.

        copy_IplImage_to_IplImage_base(projection0, &projection0_base);
        copy_IplImage_to_IplImage_base(projection, &projection_base);
        copy_IplImage_to_IplImage_base(sums_image, &sums_image_base);
        copy_IplImage_to_IplImage_base(hits_image, &hits_image_base);

        // NOTE: No dynamic memory.
        IplImage_base projection0_ROI_base;
        IplImage_base projection_ROI_base;

        copy_IplImage_ROI_widthStep_to_IplImage_base_ROI_widthStep(projection0, projection0_ROI, &projection0_ROI_base);
        copy_IplImage_ROI_widthStep_to_IplImage_base_ROI_widthStep(projection, projection_ROI, &projection_ROI_base);

        // Reuse slice later to transfer data back to the MRC file.
        Islice *slice = NULL;
        slice = sliceReadMRC(&file_header_in, i_projection, 'z');
        if (slice == NULL)
            ABORT("Cannot read file: \"%s\".\n", filesystem_params->filepath_in);

        // NOTE: Output floating-point data.
        if (slice->mode != SLICE_MODE_FLOAT)
            sliceNewMode(slice, SLICE_MODE_FLOAT);

        copy_Islice_to_IplImage(slice, projection0);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0_ROI, "/home/akulo/filter_1D_images/projection0.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
        Print("Initialization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

        t_0 = clock();

/*
#if PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA

        projection_transform_2D_CUDA(
            general_CUDA_params,
            &projection0_base,
            &projection0_ROI_base,
            support0,
    #if PROTOTYPE_COMPLIANT_ROTATION
            inv_rot_matrix_2x2, // NOTE: inv_rot_matrix_2x2 to conform to MATLAB prototype.
    #else
            rot_matrix_2x2,
    #endif
            remap_2D_data->params.local_mag,
            remap_2D_data->params.map_2D_order,
            remap_2D_data->params.n_coeffs,
            remap_2D_data->params.power_order_n_x1,
            remap_2D_data->params.power_order_n_x2,
            remap_2D_data->params.power_order,
            remap_2D_data->params.map_2D_n_x1,
            remap_2D_data->params.map_2D_n_x2,
            //remap_2D_data->params.map_2D_n_x3, // Only current map is used.
            //remap_2D_data->map_2D, // Only current map is used.
            remap_2D_data->map_2D_current,
            &sums_image_base,
            &hits_image_base,
            &projection_base,
            &projection_ROI_base,
            support);

*/
        // BEGIN: Recalculate sample factor and scale intensities.

        real_type step_remap_2D = 1.0 / remap_2D_data->params.local_mag;

        real_type step_rotate_2D;
        if (test_rotation_matrix_2x2_for_subsampling(rot_matrix_2x2))
            step_rotate_2D = 0.25;
        else 
            step_rotate_2D = 1.0; 

        real_type step = fmin(step_remap_2D, step_rotate_2D);
        //step = 1.0;  WARN("step is hardcoded to %.15e.\n", step);

        Print("step (recalculated) = %.15e\n", step);

        //real_type sample_factor = step; // MATLAB prototype uses step for number of samples per pixel.
        real_type sample_factor = step * step;
        //real_type sample_factor = 1.0; WARN("sample_factor is hardcoded to %.15e.\n", sample_factor);
        
        Print("sample_factor (recalculated) = %.15e\n", sample_factor);

        if (sample_factor != 1.0)
            cvScale(projection, projection, sample_factor, 0);

        // END: Recalculate sample factor and scale intensities.
/*

#else // PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA

        projection_transform_2D(
            projection0_ROI,
            support0,
    #if PROTOTYPE_COMPLIANT_ROTATION
            inv_rot_matrix_2x2, // NOTE: inv_rot_matrix_2x2 to conform to MATLAB prototype.
    #else
            rot_matrix_2x2,
    #endif
            remap_2D_data,
            sums_image,
            hits_image,
            projection_ROI,
            support);

#endif // PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA
*/

        copy_Islice_to_IplImage(slice, projection); // DEBUG: No remap.

        t_1 = clock();
        Print("Transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection_ROI, "/home/akulo/filter_1D_images/projection.transformed.CUDA.mrc"); // PRIMARY DEBUG.
        //ABORT("ABORT!\n");

        t_0 = clock();

        // BEGIN: Extend rows left and right in image.
        if (symmetrize_2D_params->symmetrize_2D_flag)
        {
            t_0 = clock();

            extend_IplImage_rows_analyze(projection_ROI);

            t_1 = clock();
            Print("Symmetrization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection_ROI, "/home/akulo/filter_1D_images/projection.transformed.symmetrized.mrc"); // PRIMARY DEBUG.
        }
        // END: Extend rows left and right in image.

        projection_filter_1D_fs(projection, filter_data, symmetrize_2D_params);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.transformed.filtered.CUDA.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
        Print("Filtering delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

        t_0 = clock();

#if PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA

        projection_inv_transform_2D_CUDA(
            general_CUDA_params,
            &projection_base,
            &projection_ROI_base,
            support,
    #if PROTOTYPE_COMPLIANT_ROTATION
            inv_rot_matrix_2x2, // NOTE: inv_rot_matrix_2x2 to conform to MATLAB prototype.
    #else
            rot_matrix_2x2,
    #endif
            remap_2D_data->params.local_mag,
            remap_2D_data->params.map_2D_order,
            remap_2D_data->params.n_coeffs,
            remap_2D_data->params.power_order_n_x1,
            remap_2D_data->params.power_order_n_x2,
            remap_2D_data->params.power_order,
            remap_2D_data->params.map_2D_n_x1,
            remap_2D_data->params.map_2D_n_x2,
            //remap_2D_data->params.map_2D_n_x3, // Only current map is used.
            //remap_2D_data->map_2D, // Only current map is used.
            remap_2D_data->map_2D_current,
            &projection0_base,
            &projection0_ROI_base,
            support0);

        if (sample_factor != 1.0)
            cvScale(projection0, projection0, sample_factor, 0);

#else // PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA

        projection_inv_transform_2D(
            projection_ROI,
            support,
    #if PROTOTYPE_COMPLIANT_ROTATION
            rot_matrix_2x2, // NOTE: rot_matrix_2x2 to conform to MATLAB prototype.
    #else
            inv_rot_matrix_2x2,
    #endif
            remap_2D_data,
            projection0_ROI,
            support0);

#endif // PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA

        t_1 = clock();
        Print("Restore delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

        // BEGIN: Edgetapering.
        if (edgetaper_2D_params->width > 0)
        {
            ABORT("Edgetapering must occur after the inv. rotation but before the inv. transform.   Can't be done if rotation and transform are combined.\n");
            edgetaper_2D(projection0_ROI, edgetaper_2D_params->width);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.filtered.edgetapered.mrc"); // PRIMARY DEBUG.
        }
        // END: Edgetapering.

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0_ROI, "/home/akulo/filter_1D_images/projection0.transformed.filtered.inv_transformed.CUDA.mrc"); // PRIMARY DEBUG.
        ABORT("Abort!\n");

        copy_IplImage_to_Islice(projection0, slice);

        sliceWriteMRC(&file_header_out, slice, i_index, 'z');
        sliceFree(slice);

#if PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA

        projection_series_transform_2D_CUDA_filter_1D_inv_transform_2D_CUDA_clear_GPU(deviceID);

#else // PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA

#endif // PROJECTION_SERIES_TRANSFORM_2D_CUDA_FILTER_1D_INV_TRANSFORM_2D_CUDA_ENABLE_CUDA

        // Free resources.
        cvReleaseImageHeader(&projection0_ROI);
        cvReleaseImageHeader(&projection_ROI);
        free_imageData_cvReleaseImageHeader(projection0);
        free_imageData_cvReleaseImageHeader(projection);
        free_imageData_cvReleaseImageHeader(sums_image);
        free_imageData_cvReleaseImageHeader(hits_image);

        filter_1D_data_release(filter_data);

#ifdef MEX
        mexEvalString("fprintf('Done processing projection.\\n');");
#endif
    }

    calc_mmm_of_MRC(&file_header_out);

    // Free resources.
    remap_2D_data_release(remap_2D_data);
    rotate_2D_data_release(rotate_2D_data);

    if (fclose(file_in))
        ABORT("Cannot close input file: \"%s\".\n", filesystem_params->filepath_in);

    if (fclose(file_out))
        ABORT("Cannot close output file \"%s\".\n", filesystem_params->filepath_out);

    t_1_total = clock();
 
    Print("Total duration: %.15e\n", secs_per_clock * ((double) t_1_total - (double) t_0_total));
}

#endif // FILTER_1D_ENABLE_CUDA

// NOTE: Uses support calculations from MATLAB prototype.
void projection_series_transform_2D_filter_1D_inv_transform_2D
(
    filesystem_params_handle filesystem_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    remap_2D_params_handle remap_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
)
{
#if !PROTOTYPE_COMPLIANT_INDEXING
    ABORT("!PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
#endif

// PRIMARY DEBUG:     Print("filesystem_params->filepath_in = \"%s\"\n", filesystem_params->filepath_in); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("filesystem_params->filepath_out = \"%s\"\n", filesystem_params->filepath_out); // PRIMARY DEBUG.

    clock_t t_0_total, t_1_total;
    clock_t t_0, t_1;

    double secs_per_clock = 1.0 / (double) CLOCKS_PER_SEC;

    t_0_total = clock();

    // BEGIN: Open input MRC file.
    MrcHeader file_header_in;

    FILE *file_in = NULL;

    openMRCFile_general(filesystem_params->filepath_in, &file_header_in, &file_in);
    // END: Open input MRC file.

    int n0_x = file_header_in.nx;
    int n0_y = file_header_in.ny;
    int n_projections = file_header_in.nz;
    //int n_projections = 1; WARN("n_projections hardcoded to %i.\n", n_projections);
    int mode = file_header_in.mode;

    // BEGIN: Create and open output MRC file.
    // BEGIN: Create output MRC file.
    if (mode != SLICE_MODE_FLOAT)
        Print("WARNING: Input mode == %i but using mode SLICE_MODE_FLOAT == %i for output!\n", mode, SLICE_MODE_FLOAT);

    MrcHeader file_header_out;

    FILE *file_out = NULL;

    createNewMRCFile(filesystem_params->filepath_out, n0_x, n0_y, filesystem_params->n_projections_indexed, /*mode*/ SLICE_MODE_FLOAT);
    // END: Create output MRC file.

    // BEGIN: Open output MRC file.
    openMRCFile_general(filesystem_params->filepath_out, &file_header_out, &file_out);
    // END: Open output MRC file.
    // END: Create and open output MRC file.

    remap_2D_data_type *remap_2D_data = create_remap_2D_data(remap_2D_params);

    //for (int i_projection = 0; i_projection < n_projections; i_projection++)
    for (int i_index = 0; i_index < filesystem_params->n_projections_indexed; ++i_index)
    {
        int i_projection = filesystem_params->projection_indices[i_index];
// PRIMARY DEBUG:         Print("i_projection = %i\n", i_projection); // PRIMARY DEBUG.

        if ((i_projection < 0) || (i_projection >= n_projections))
            ABORT("i_projection == %i, but n_projections == %i.\n", i_projection, n_projections);

        t_0 = clock();

        set_map_2D_current(remap_2D_data, i_projection);

        int row_pad = 1;
        int col_pad = 1;

        int n0_x_padded = n0_x + row_pad;
        int n0_y_padded = n0_y + col_pad;

        real_type support0[4];

#if PROTOTYPE_COMPLIANT_INDEXING
        support0[0] = ((real_type) n0_x + 1.0) / 2.0; // y-coordinate of projection center
        support0[1] = ((real_type) n0_y + 1.0) / 2.0; // y-coordinate of projection center
        support0[2] = (real_type) n0_x; // max x-coordinate of projection
        support0[3] = (real_type) n0_y; // max y-coordinate of projection
#else
        support0[0] = ((real_type) n0_x - 1.0) / 2.0; // y-coordinate of projection center
        support0[1] = ((real_type) n0_y - 1.0) / 2.0; // y-coordinate of projection center
        support0[2] = (real_type) n0_x - 1.0; // max x-coordinate of projection
        support0[3] = (real_type) n0_y - 1.0; // max y-coordinate of projection
#endif

        real_type support[4];
        support[0] = support0[0]; // x-coordinate of transformed projection center
        support[1] = support0[1]; // y-coordinate of transformed projection center
        support[2] = support0[2]; // max x-coordinate of transformed projection
        support[3] = support0[3]; // max y-coordinate of transformed projection

        calc_transformed_support(remap_2D_data, support);

#if PROTOTYPE_COMPLIANT_INDEXING
        int n_x = round(support[2]);
        int n_y = round(support[3]);
#else
        int n_x = round(support[2] + 1.0);
        int n_y = round(support[3] + 1.0);
#endif

        int n_x_padded = n_x + row_pad;
        int n_y_padded = n_y + col_pad;

        // Calculate filter.
        filter_1D_data_type *filter_data = create_filter_1D_data_n_xy(filter_1D_params, n_x_padded, n_y_padded);

        IplImage *projection0 = create_IplImage_SetZero_malloc(n0_x_padded, n0_y_padded, IPL_DEPTH_REAL);
        IplImage *projection = create_IplImage_SetZero_malloc(n_x_padded, n_y_padded, IPL_DEPTH_REAL);
        IplImage *sums_image = create_IplImage_SetZero_malloc(n_x, n_y, IPL_DEPTH_REAL);
        IplImage *hits_image = create_IplImage_SetZero_malloc(n_x, n_y, IPL_DEPTH_REAL);

        IplImage *projection0_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                        projection0, cvRect(0, 0, n0_x, n0_y));

        IplImage *projection_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                       projection, cvRect(0, 0, n_x, n_y));

        // Reuse slice later to transfer data back to the MRC file.
        Islice *slice = NULL;
        slice = sliceReadMRC(&file_header_in, i_projection, 'z');
        if (slice == NULL)
            ABORT("Cannot read file: \"%s\".\n", filesystem_params->filepath_in);

        // NOTE: Output floating-point data.
        if (slice->mode != SLICE_MODE_FLOAT)
            sliceNewMode(slice, SLICE_MODE_FLOAT);

        copy_Islice_to_IplImage(slice, projection0);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0_ROI, "/home/akulo/filter_1D_images/projection0.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Initialization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        t_0 = clock();

        projection_transform_2D(
            projection0_ROI,
            support0,
            remap_2D_data,
            sums_image,
            hits_image,
            projection_ROI,
            support);

        //copy_Islice_to_IplImage(slice, projection); // DEBUG: No remap.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection_ROI, "/home/akulo/filter_1D_images/projection.transformed.NOCUDA.mrc"); // PRIMARY DEBUG.

        t_0 = clock();

        // BEGIN: Extend rows left and right in image.
        if (symmetrize_2D_params->symmetrize_2D_flag)
        {
            t_0 = clock();

            extend_IplImage_rows_analyze(projection_ROI);

            t_1 = clock();
// PRIMARY DEBUG:             Print("Symmetrization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection_ROI, "/home/akulo/filter_1D_images/projection.transformed.symmetrized.mrc"); // PRIMARY DEBUG.
        }
        // END: Extend rows left and right in image.

        projection_filter_1D_fs(projection, filter_data, symmetrize_2D_params);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.transformed.filtered.NOCUDA.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Filtering delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        t_0 = clock();

        projection_inv_transform_2D(
            projection_ROI,
            support,
            remap_2D_data,
            projection0_ROI,
            support0);

        t_1 = clock();
// PRIMARY DEBUG:         Print("Restore delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        // BEGIN: Edgetapering.
        if (edgetaper_2D_params->width > 0)
        {
            ABORT("Edgetapering must occur after the inv. rotation but before the inv. transform.   Can't be done if rotation and transform are combined.\n");
            edgetaper_2D(projection0_ROI, edgetaper_2D_params->width);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.filtered.edgetapered.mrc"); // PRIMARY DEBUG.
        }
        // END: Edgetapering.

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection0_ROI, "/home/akulo/filter_1D_images/projection0.transformed.filtered.inv_transformed.NOCUDA.mrc"); // PRIMARY DEBUG.
//        WARN("Returning...\n"); return;
        //ABORT("Abort!\n");

        copy_IplImage_to_Islice(projection0, slice);

        sliceWriteMRC(&file_header_out, slice, i_index, 'z');
        sliceFree(slice);

        // Free resources.
        cvReleaseImageHeader(&projection0_ROI);
        cvReleaseImageHeader(&projection_ROI);
        free_imageData_cvReleaseImageHeader(projection0);
        free_imageData_cvReleaseImageHeader(projection);
        free_imageData_cvReleaseImageHeader(sums_image);
        free_imageData_cvReleaseImageHeader(hits_image);

        filter_1D_data_release(filter_data);

#ifdef MEX
        mexEvalString("fprintf('Done processing projection.\\n');");
#endif
    }

    calc_mmm_of_MRC(&file_header_out);

    // Free resources.
    remap_2D_data_release(remap_2D_data);

    if (fclose(file_in))
        ABORT("Cannot close input file: \"%s\".\n", filesystem_params->filepath_in);

    if (fclose(file_out))
        ABORT("Cannot close output file \"%s\".\n", filesystem_params->filepath_out);

    t_1_total = clock();
 
    Print("Total duration: %.15e\n", secs_per_clock * ((double) t_1_total - (double) t_0_total));
}

void projection_series_remap_2D_rotate_2D_filter_1D_rotate_2D_restore_2D
(
    filesystem_params_handle filesystem_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    remap_2D_params_handle remap_2D_params,
    rotate_2D_params_handle rotate_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
) 
{
#if PROTOTYPE_COMPLIANT_INDEXING
    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
#endif

    Print("filesystem_params->filepath_in = \"%s\"\n", filesystem_params->filepath_in);
    Print("filesystem_params->filepath_out = \"%s\"\n", filesystem_params->filepath_out);

    clock_t t_0_total, t_1_total;
    clock_t t_0, t_1;

    double secs_per_clock = 1.0 / (double) CLOCKS_PER_SEC;

    t_0_total = clock();

    // BEGIN: Open input MRC file.
    MrcHeader file_header_in;

    FILE *file_in = NULL;

    openMRCFile_general(filesystem_params->filepath_in, &file_header_in, &file_in);
    // END: Open input MRC file.

    int n0_x = file_header_in.nx;
    int n0_y = file_header_in.ny;
    int n_projections = file_header_in.nz;
    //int n_projections = 1; WARN("n_projections hardcoded to %i.\n", n_projections);
    int mode = file_header_in.mode;

    // BEGIN: Create and open output MRC file.
    // BEGIN: Create output MRC file.
    if (mode != SLICE_MODE_FLOAT)
        Print("WARNING: Input mode == %i but using mode SLICE_MODE_FLOAT == %i for output!\n", mode, SLICE_MODE_FLOAT);

    MrcHeader file_header_out;

    FILE *file_out = NULL;

    createNewMRCFile(filesystem_params->filepath_out, n0_x, n0_y, filesystem_params->n_projections_indexed, /*mode*/ SLICE_MODE_FLOAT);
    // END: Create output MRC file.

    // BEGIN: Open output MRC file.
    openMRCFile_general(filesystem_params->filepath_out, &file_header_out, &file_out);
    // END: Open output MRC file.
    // END: Create and open output MRC file.

    remap_2D_data_type *remap_2D_data = create_remap_2D_data(remap_2D_params);

    IplImage *sums_image = create_IplImage_SetZero_malloc(n0_x, n0_y, IPL_DEPTH_REAL);
    IplImage *hits_image = create_IplImage_SetZero_malloc(n0_x, n0_y, IPL_DEPTH_REAL);

    int row_pad = 2;
    int col_pad = 2;

/*
    //int row_pad, col_pad;
    if (test_rotation_matrix_2x2_for_subsampling_threshold(rotate_2D_params->rot_matrix_2x2))
    { 
        //row_pad = 1; col_pad = 1; 
    }
    else
    { 
        test_rotation_matrix_2x2_for_subsampling_threshold(rotate_2D_params->inv_rot_matrix_2x2);
        //row_pad = 2; col_pad = 2; 
    }
    print_matrix_2x2("rotate_2D_params->rot_matrix_2x2 (thresholded)", rotate_2D_params->rot_matrix_2x2);
    print_matrix_2x2("rotate_2D_params->inv_rot_matrix_2x2 (thresholded)", rotate_2D_params->inv_rot_matrix_2x2);
*/

    real_type support0[4];
    support0[0] = ((real_type) n0_x - 1.0) / 2.0; // y-coordinate of projection center
    support0[1] = ((real_type) n0_y - 1.0) / 2.0; // y-coordinate of projection center
    support0[2] = (real_type) n0_x - 1.0; // max x-coordinate of projection
    support0[3] = (real_type) n0_y - 1.0; // max y-coordinate of projection

    real_type support[4];
    support[0] = support0[0]; // x-coordinate of transformed projection center
    support[1] = support0[1]; // y-coordinate of transformed projection center
    support[2] = support0[2]; // max x-coordinate of transformed projection
    support[3] = support0[3]; // max y-coordinate of transformed projection

#if PROTOTYPE_COMPLIANT_ROTATION
    // NOTE: I think the transform should be the same as that used to transform the original projection, but cf. the MATLAB prototype.
    ABORT("Uses original form of calc_transformed_support().\n");
    //calc_transformed_support(rotate_2D_params->rot_matrix_2x2, support);
#else
    ABORT("Uses original form of calc_transformed_support().\n");
    //calc_transformed_support(rotate_2D_params->inv_rot_matrix_2x2, support);
#endif

    real_type support00[4];
    support00[0] = support[0]; // x-coordinate of untransformed projection center
    support00[1] = support[1]; // y-coordinate of untransformed projection center
    support00[2] = support[2]; // max x-coordinate of untransformed projection
    support00[3] = support[3]; // max y-coordinate of untransformed projection

#if PROTOTYPE_COMPLIANT_ROTATION
    // NOTE: I think the transform should be the same as that used to transform the transformed projection, but cf. the MATLAB prototype.
    ABORT("Uses original form of calc_transformed_support().\n");
    //calc_transformed_support(rotate_2D_params->inv_rot_matrix_2x2, support00);
#else
    ABORT("Uses original form of calc_transformed_support().\n");
    //calc_transformed_support(rotate_2D_params->rot_matrix_2x2, support00);
#endif

    int n_x = round(support[2] + 1.0);
    int n_y = round(support[3] + 1.0);
    int n00_x = round(support00[2] + 1.0);
    int n00_y = round(support00[3] + 1.0);

    int n00_x_padded = n00_x + row_pad;
    int n00_y_padded = n00_y + col_pad;

    // Ensure that the parities of the extents of projection00 match the parities of the extents of the source image.
    if ((n0_x % 2) != (n00_x_padded % 2)) ++n00_x_padded;
    if ((n0_y % 2) != (n00_y_padded % 2)) ++n00_y_padded;

    int n_x_padded = n_x + row_pad;
    int n_y_padded = n_y + col_pad;
//    if ((n00_x_padded % 2) != (n_x_padded % 2)) ++n_x_padded;
//    if ((n00_y_padded % 2) != (n_y_padded % 2)) ++n_y_padded;

    // Calculate filter.
    filter_1D_data_type *filter_data = create_filter_1D_data_n_xy(filter_1D_params, n_x_padded, n_y_padded);

    IplImage *projection00 = create_IplImage_SetZero_malloc(n00_x_padded, n00_y_padded, IPL_DEPTH_REAL); // NOTE: Used for both original projection and untransformed filtered transformed projection.
    IplImage *projection = create_IplImage_SetZero_malloc(n_x_padded, n_y_padded, IPL_DEPTH_REAL); // NOTE: Used for both transformed projection and filtered transformed projection.

    center_ROI_params_type center_params;

    center_params = calc_center_ROI_params(n00_x_padded, n00_y_padded, n0_x, n0_y);
    CvRect projection00_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                   center_params.n_i0_x, center_params.n_i0_y);
    IplImage *projection00_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                            projection00, 
                                            projection00_center_ROI_CvRect);

    center_params = calc_center_ROI_params(n_x_padded, n_y_padded, n0_x, n0_y);
    IplImage *projection_center_ROI_remap = create_IplImage_ROI_widthStep_of_IplImage(
                                                projection, 
                                                cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                    center_params.n_i0_x, center_params.n_i0_y));

    center_params = calc_center_ROI_params(n_x_padded, n_y_padded, n_x, n_y);
    IplImage *projection_center_ROI_transform = create_IplImage_ROI_widthStep_of_IplImage(
                                                    projection, 
                                                    cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                        center_params.n_i0_x, center_params.n_i0_y));

    ABORT("This needs to be rewritten to handle a different rotation angle for each projection.\n");
    //for (int i_projection = 0; i_projection < n_projections; i_projection++)
    for (int i_index = 0; i_index < filesystem_params->n_projections_indexed; ++i_index)
    {
        int i_projection = filesystem_params->projection_indices[i_index];
        Print("i_projection = %i\n", i_projection);

        if ((i_projection < 0) || (i_projection >= n_projections))
            ABORT("i_projection == %i, but n_projections == %i.\n", i_projection, n_projections);

        set_map_2D_current(remap_2D_data, i_projection);

        // Reuse slice later to transfer data back to the MRC file.
        Islice *slice = NULL;
        slice = sliceReadMRC(&file_header_in, i_projection, 'z');
        if (slice == NULL)
            ABORT("Cannot read file: \"%s\".\n", filesystem_params->filepath_in);

        // NOTE: Output floating-point data.
        if (slice->mode != SLICE_MODE_FLOAT)
            sliceNewMode(slice, SLICE_MODE_FLOAT);

        if (remap_2D_params->map_2D_order > 0)
        {
            copy_Islice_to_IplImage_center(slice, projection);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection_center_ROI_remap, "/home/akulo/filter_1D_images/projection.mrc"); // PRIMARY DEBUG.

            t_0 = clock();

            projection_remap_2D_prototype(
                remap_2D_data,
                projection_center_ROI_remap,
                sums_image,
                hits_image,
                projection00_center_ROI);

            t_1 = clock();
            Print("Remap 2D delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));
        }
        else
        {
            copy_Islice_to_IplImage_center(slice, projection00);
        }

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.mrc"); // PRIMARY DEBUG.

        if (symmetrize_2D_params->symmetrize_2D_flag)
        {
            t_0 = clock();

            symmetrize_IplImage_values_CvRect(projection00, projection00_center_ROI_CvRect);

            t_1 = clock();
            Print("Symmetrization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));
        }

        t_0 = clock();

        projection_rotate_2D_prototype(
            projection00_center_ROI, support0,
            rotate_2D_params->local_mag,
#if PROTOTYPE_COMPLIANT_ROTATION
            // NOTE: inv_rot_matrix_2x2 to conform to MATLAB prototype.
            rotate_2D_params->inv_rot_matrix_2x2,
#else
            rotate_2D_params->rot_matrix_2x2,
#endif
            projection_center_ROI_transform, support);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
        Print("Transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

        t_0 = clock();

        projection_filter_1D_fs(projection, filter_data, symmetrize_2D_params);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.filtered.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
        Print("Filtering delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

        t_0 = clock();

        projection_rotate_2D_prototype(
            projection_center_ROI_transform, support,
            rotate_2D_params->local_mag,
#if PROTOTYPE_COMPLIANT_ROTATION
            // NOTE: rot_matrix_2x2 to conform to MATLAB prototype.
            rotate_2D_params->rot_matrix_2x2,
#else
            rotate_2D_params->inv_rot_matrix_2x2,
#endif
            projection00_center_ROI, support0);

        t_1 = clock();
        Print("Inv. transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.filtered.mrc"); // PRIMARY DEBUG.

        // BEGIN: Edgetapering.
        if (edgetaper_2D_params->width > 0)
        {
            edgetaper_2D(projection00_center_ROI, edgetaper_2D_params->width);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.filtered.edgetapered.mrc"); // PRIMARY DEBUG.
        }
        // END: Edgetapering.

        if (remap_2D_params->map_2D_order > 0)
        {
            t_0 = clock();

            projection_restore_2D_prototype(
                remap_2D_data,
                projection00_center_ROI,
                projection_center_ROI_remap);

            t_1 = clock();
            Print("Restore 2D delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

            copy_IplImage_to_Islice_center(projection_center_ROI_remap, slice);
        }
        else
        {
            copy_IplImage_to_Islice_center(projection00_center_ROI, slice);
        }

        sliceWriteMRC(&file_header_out, slice, i_index, 'z');
        sliceFree(slice);

#ifdef MEX
        mexEvalString("fprintf('Done processing projection.\\n');");
#endif
    }

    calc_mmm_of_MRC(&file_header_out);

    // Free resources.

    cvReleaseImageHeader(&projection_center_ROI_remap);
    cvReleaseImageHeader(&projection_center_ROI_transform);
    cvReleaseImageHeader(&projection00_center_ROI);
    free_imageData_cvReleaseImageHeader(projection00);
    free_imageData_cvReleaseImageHeader(projection);
    free_imageData_cvReleaseImageHeader(sums_image);
    free_imageData_cvReleaseImageHeader(hits_image);

    filter_1D_data_release(filter_data);

    if (fclose(file_in))
        ABORT("Cannot close input file: \"%s\".\n", filesystem_params->filepath_in);

    if (fclose(file_out))
        ABORT("Cannot close output file \"%s\".\n", filesystem_params->filepath_out);

    t_1_total = clock();
 
    Print("Total duration: %.15e\n", secs_per_clock * ((double) t_1_total - (double) t_0_total));
}

#ifdef FILTER_1D_ENABLE_CUDA

#define PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA 1

// NOTE: Uses IplImage for workspace.  Pulls w/out index check.
void projection_series_rotate_2D_CUDA_filter_1D_rotate_2D_CUDA
(
    filesystem_params_handle filesystem_params,
    general_CUDA_params_handle general_CUDA_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    rotate_2D_params_handle rotate_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
) 
{
// NOTE: Does NOT use PROTOTYPE_COMPLIANT_INDEXING.
//#if PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

    if (PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA == 1)
        Print("PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA == 1\n");
    else
        Print("PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA == %d != 1\n", PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA);

    Print("filesystem_params->filepath_in = \"%s\"\n", filesystem_params->filepath_in);
    Print("filesystem_params->filepath_out = \"%s\"\n", filesystem_params->filepath_out);

    clock_t t_0_total, t_1_total;
    clock_t t_0, t_1;

    double secs_per_clock = 1.0 / (double) CLOCKS_PER_SEC;

    t_0_total = clock();

    // BEGIN: Open input MRC file.
    MrcHeader file_header_in;

    FILE *file_in = NULL;

    openMRCFile_general(filesystem_params->filepath_in, &file_header_in, &file_in);
    // END: Open input MRC file.

    int n0_x = file_header_in.nx;
    int n0_y = file_header_in.ny;
    int n_projections = file_header_in.nz;
    //int n_projections = 1; WARN("n_projections hardcoded to %i.\n", n_projections);
    int mode = file_header_in.mode;

    // BEGIN: Create and open output MRC file.
    // BEGIN: Create output MRC file.
    if (mode != SLICE_MODE_FLOAT)
        WARN("Input mode == %i but using mode SLICE_MODE_FLOAT == %i for output!\n", mode, SLICE_MODE_FLOAT);

    MrcHeader file_header_out;

    FILE *file_out = NULL;

    createNewMRCFile(filesystem_params->filepath_out, n0_x, n0_y, filesystem_params->n_projections_indexed, /*mode*/ SLICE_MODE_FLOAT);
    // END: Create output MRC file.

    // BEGIN: Open output MRC file.
    openMRCFile_general(filesystem_params->filepath_out, &file_header_out, &file_out);
    // END: Open output MRC file.
    // END: Create and open output MRC file.

    rotate_2D_data_type *rotate_2D_data = create_rotate_2D_data(rotate_2D_params, filesystem_params->n_projections_indexed);

    //for (int i_projection = 0; i_projection < n_projections; i_projection++)
    for (int i_index = 0; i_index < filesystem_params->n_projections_indexed; ++i_index)
    {
        int i_projection = filesystem_params->projection_indices[i_index];
        Print("i_projection = %i\n", i_projection);

        if ((i_projection < 0) || (i_projection >= n_projections))
            ABORT("i_projection == %i, but n_projections == %i.\n", i_projection, n_projections);

        t_0 = clock();

        real_type rot_matrix_2x2[2][2];
        real_type inv_rot_matrix_2x2[2][2];

        calc_rotation_matrices_current(rotate_2D_data, i_index, rot_matrix_2x2, inv_rot_matrix_2x2);

        // NOTE: Padding on all four sides of the image.
        //int row_pad = 2;
        //int col_pad = 2;
        int row_pad = 4;
        int col_pad = 4;

/*
        //int row_pad, col_pad;
        if (test_rotation_matrix_2x2_for_subsampling_threshold(rot_matrix_2x2))
        { 
            //row_pad = 1; col_pad = 1; 
        }
        else
        { 
            test_rotation_matrix_2x2_for_subsampling_threshold(inv_rot_matrix_2x2);
            //row_pad = 2; col_pad = 2; 
        }
        print_matrix_2x2("rot_matrix_2x2 (thresholded)", rot_matrix_2x2);
        print_matrix_2x2("inv_rot_matrix_2x2 (thresholded)", inv_rot_matrix_2x2);
*/

        real_type support0[4];
        support0[0] = ((real_type) n0_x - 1.0) / 2.0; // y-coordinate of projection center
        support0[1] = ((real_type) n0_y - 1.0) / 2.0; // y-coordinate of projection center
        support0[2] = (real_type) n0_x - 1.0; // max x-coordinate of projection
        support0[3] = (real_type) n0_y - 1.0; // max y-coordinate of projection

        real_type support[4];
        support[0] = support0[0]; // x-coordinate of transformed projection center
        support[1] = support0[1]; // y-coordinate of transformed projection center
        support[2] = support0[2]; // max x-coordinate of transformed projection
        support[3] = support0[3]; // max y-coordinate of transformed projection

#if PROTOTYPE_COMPLIANT_ROTATION
        // NOTE: I think the transform should be the same as that used to transform the original projection, but cf. the MATLAB prototype.
        calc_rotated_support(rot_matrix_2x2, support);
#else
        calc_rotated_support(inv_rot_matrix_2x2, support);
#endif

        real_type support00[4];
        support00[0] = support[0]; // x-coordinate of untransformed projection center
        support00[1] = support[1]; // y-coordinate of untransformed projection center
        support00[2] = support[2]; // max x-coordinate of untransformed projection
        support00[3] = support[3]; // max y-coordinate of untransformed projection

#if PROTOTYPE_COMPLIANT_ROTATION
        // NOTE: I think the transform should be the same as that used to transform the transformed projection, but cf. the MATLAB prototype.
        calc_rotated_support(inv_rot_matrix_2x2, support00);
#else
        calc_rotated_support(rot_matrix_2x2, support00);
#endif

        int n_x = round(support[2] + 1.0);
        int n_y = round(support[3] + 1.0);
        int n00_x = round(support00[2] + 1.0);
        int n00_y = round(support00[3] + 1.0);

        int n00_x_padded_centered = n00_x + col_pad;
        int n00_y_padded_centered = n00_y + row_pad;

        //int dimBlock_x = general_CUDA_params_get_dimBlock_x(general_CUDA_params);
        //int dimBlock_y = general_CUDA_params_get_dimBlock_y(general_CUDA_params);

        // NOTE: Padding on strictly the top and right of the image.
        //int n00_col_pad_CUDA = calc_extent_pad_CUDA(n00_x, dimBlock_x);
        //int n00_row_pad_CUDA = calc_extent_pad_CUDA(n00_y, dimBlock_y);

        //int n00_x_padded = n00_x_padded_centered + n00_col_pad_CUDA;
        //int n00_y_padded = n00_y_padded_centered + n00_row_pad_CUDA;
        int n00_x_padded = n00_x_padded_centered;
        int n00_y_padded = n00_y_padded_centered;

        // Ensure that the parities of the extents of projection00 match the parities of the extents of the source image.
        if ((n0_x % 2) != (n00_x_padded % 2)) { ++n00_x_padded; /*++n00_col_pad_CUDA;*/ }
        if ((n0_y % 2) != (n00_y_padded % 2)) { ++n00_y_padded; /*++n00_row_pad_CUDA;*/ }

        Print("n00_(xy)_padded_centered = (%d, %d)\n", n00_x_padded_centered, n00_y_padded_centered);
        //Print("n00_(col, row)_pad_CUDA = (%d, %d)\n", n00_col_pad_CUDA, n00_row_pad_CUDA);
        Print("n00_(xy)_padded = (%d, %d)\n", n00_x_padded, n00_y_padded);

        int n_x_padded_centered = n_x + row_pad;
        int n_y_padded_centered = n_y + col_pad;
        //int n_x_padded = n_x + 1;
        //int n_y_padded = n_y + 1;

        // NOTE: Padding strictly on the top and right of the image.
        //int n_col_pad_CUDA = calc_extent_pad_CUDA(n_x, dimBlock_x);
        //int n_row_pad_CUDA = calc_extent_pad_CUDA(n_y, dimBlock_y);

        //int n_x_padded = n_x_padded_centered + n_col_pad_CUDA;
        //int n_y_padded = n_y_padded_centered + n_row_pad_CUDA;
        int n_x_padded = n_x_padded_centered;
        int n_y_padded = n_y_padded_centered;

        // Ensure that the parities of the extents of projection match the parities of the extents of projection00.
        if ((n00_x_padded % 2) != (n_x_padded % 2)) { ++n_x_padded; /*++n_col_pad_CUDA;*/ }
        if ((n00_y_padded % 2) != (n_y_padded % 2)) { ++n_y_padded; /*++n_row_pad_CUDA;*/ }

        Print("n_(xy)_padded_centered = (%d, %d)\n", n_x_padded_centered, n_y_padded_centered);
        //Print("n_(col, row)_pad_CUDA = (%d, %d)\n", n_col_pad_CUDA, n_row_pad_CUDA);
        Print("n_(xy)_padded = (%d, %d)\n", n_x_padded, n_y_padded);

        // Calculate filter.
        filter_1D_data_type *filter_data = create_filter_1D_data_n_xy(filter_1D_params, n_x, n_y);

#if PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        int deviceID = 1;
        projection_series_rotate_2D_CUDA_filter_1D_rotate_2D_CUDA_init_GPU(deviceID);

        pixel_type *projection00_data = NULL;
        malloc_page_locked_pixel_type_CUDA(&projection00_data, n00_x_padded * n00_y_padded * sizeof(pixel_type));
        IplImage *projection00 = create_IplImage_XXF_from_pixel_type_data(n00_x_padded, n00_y_padded, projection00_data); // NOTE: Used for both original projection and untransformed filtered transformed projection.
        cvSetZero(projection00);

        pixel_type *projection_data = NULL;
        malloc_page_locked_pixel_type_CUDA(&projection_data, n_x_padded * n_y_padded * sizeof(pixel_type));
        IplImage *projection = create_IplImage_XXF_from_pixel_type_data(n_x_padded, n_y_padded, projection_data); // NOTE: Used for both transformed projection and filtered transformed projection.
        cvSetZero(projection);

#else // PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        //IplImage *projection0 = create_IplImage_SetZero_malloc(n0_x, n0_y, IPL_DEPTH_REAL); // DEBUG: CUDA texturing. NOTE: Used for the original projection only.
        IplImage *projection00 = create_IplImage_SetZero_malloc(n00_x_padded, n00_y_padded, IPL_DEPTH_REAL); // NOTE: Used for both original projection and untransformed filtered transformed projection.
        IplImage *projection = create_IplImage_SetZero_malloc(n_x_padded, n_y_padded, IPL_DEPTH_REAL); // NOTE: Used for both transformed projection and filtered transformed projection.

#endif // PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        center_ROI_params_type center_params;

        center_params = calc_center_ROI_params(n00_x_padded/*_centered*/, n00_y_padded/*_centered*/, n0_x, n0_y);
        CvRect projection00_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                       center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection00_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                                projection00, 
                                                projection00_center_ROI_CvRect);

        center_params = calc_center_ROI_params(n_x_padded/*_centered*/, n_y_padded/*_centered*/, n_x, n_y);
        CvRect projection_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                     center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                              projection, 
                                              projection_center_ROI_CvRect);

//        IplImage *projection_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
//                                              projection, 
//                                              cvRect(0, 0, n_x, n_y));

        // NOTE: No dynamic memory.
        //IplImage_base projection0_base; // DEBUG: CUDA texturing.
        IplImage_base projection00_base;
        IplImage_base projection_base;

        //copy_IplImage_to_IplImage_base(projection0, &projection0_base); // DEBUG: CUDA texturing.
        copy_IplImage_to_IplImage_base(projection00, &projection00_base);
        copy_IplImage_to_IplImage_base(projection, &projection_base);

        // NOTE: No dynamic memory.
        IplImage_base projection00_center_ROI_base;
        IplImage_base projection_center_ROI_base;

        copy_IplImage_ROI_widthStep_to_IplImage_base_ROI_widthStep(projection00, projection00_center_ROI, &projection00_center_ROI_base);
        copy_IplImage_ROI_widthStep_to_IplImage_base_ROI_widthStep(projection, projection_center_ROI, &projection_center_ROI_base);

        // Reuse slice later to transfer data back to the MRC file.
        Islice *slice = NULL;
        slice = sliceReadMRC(&file_header_in, i_projection, 'z');
        if (slice == NULL)
            ABORT("Cannot read file: \"%s\".\n", filesystem_params->filepath_in);

        // NOTE: Output floating-point data.
        if (slice->mode != SLICE_MODE_FLOAT)
            sliceNewMode(slice, SLICE_MODE_FLOAT);

        //copy_Islice_to_IplImage_center(slice, projection0); // DEBUG: CUDA texturing.
        copy_Islice_to_IplImage_center(slice, projection00_center_ROI);

        //t_0 = clock();

        //write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection0.mrc"); // DEBUG: CUDA texturing.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.mrc"); // PRIMARY DEBUG.
        //copy_Islice_to_IplImage_center(slice, projection_center_ROI); // DEBUG: Copy slice to projection.
        //write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.mrc");
        //ABORT("ABORT!\n");

        //t_1 = clock();
        //Print("write_MRC_image_from_IplImage() t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));
        //ABORT("Abort!\n");

        t_1 = clock();
        Print("Initialization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

        if (symmetrize_2D_params->symmetrize_2D_flag)
        {
            t_0 = clock();

            symmetrize_IplImage_values_CvRect(projection00, projection00_center_ROI_CvRect);

            t_1 = clock();
            Print("Symmetrization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));
        }

        t_0 = clock();

#if PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        projection_rotate_2D_prototype_CUDA(
            general_CUDA_params,
            &projection00_base, &projection00_center_ROI_base, support0, // ORIGINAL.
            //&projection0_base, &projection0_base, support0, // DEBUG: CUDA texturing.
    #if PROTOTYPE_COMPLIANT_ROTATION
            inv_rot_matrix_2x2, // NOTE: inv_rot_matrix_2x2 to conform to MATLAB prototype.
    #else
            rot_matrix_2x2,
    #endif
            &projection_base, &projection_center_ROI_base, support);

#else // PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        projection_rotate_2D_prototype(
            projection00_center_ROI, support0,
    #if PROTOTYPE_COMPLIANT_ROTATION
            inv_rot_matrix_2x2, // NOTE: inv_rot_matrix_2x2 to conform to MATLAB prototype.
    #else
            rot_matrix_2x2,
    #endif
            projection_center_ROI, support);

#endif // PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        t_1 = clock();
        Print("Transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

        //write_MRC_image_from_IplImage(projection0, "/home/akulo/filter_1D_images/projection0.mrc"); // DEBUG: CUDA texturing.
        //write_MRC_image_from_IplImage(projection_center_ROI, "/home/akulo/filter_1D_images/projection_center_ROI.mrc");
        //ABORT("ABORT!\n");

        t_0 = clock();

        projection_filter_1D_fs(projection_center_ROI, filter_data, symmetrize_2D_params);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.filtered.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
        Print("Filtering delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

        t_0 = clock();

#if PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        projection_rotate_2D_prototype_CUDA(
            general_CUDA_params,
            &projection_base, &projection_center_ROI_base, support,
    #if PROTOTYPE_COMPLIANT_ROTATION
            rot_matrix_2x2, // NOTE: rot_matrix_2x2 to conform to MATLAB prototype.
    #else
            inv_rot_matrix_2x2,
    #endif
            &projection00_base, &projection00_center_ROI_base, support0);

#else // PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        projection_rotate_2D_prototype(
            projection_center_ROI, support,
    #if PROTOTYPE_COMPLIANT_ROTATION
            rot_matrix_2x2, // NOTE: rot_matrix_2x2 to conform to MATLAB prototype.
    #else
            inv_rot_matrix_2x2,
    #endif
            projection00_center_ROI, support0);

#endif // PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        t_1 = clock();
        Print("Inv. transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0));

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.filtered.mrc"); // PRIMARY DEBUG.

        // BEGIN: Edgetapering.
        if (edgetaper_2D_params->width > 0)
        {
            edgetaper_2D(projection00_center_ROI, edgetaper_2D_params->width);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.filtered.edgetapered.mrc"); // PRIMARY DEBUG.
        }
        // END: Edgetapering.

        copy_IplImage_to_Islice_center(projection00_center_ROI, slice);

        sliceWriteMRC(&file_header_out, slice, i_index, 'z');
        sliceFree(slice);

        // BEGIN: Free resources.

        cvReleaseImageHeader(&projection_center_ROI);
        cvReleaseImageHeader(&projection00_center_ROI);

#if PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        free_page_locked_pixel_type_CUDA(projection00_data);
        cvReleaseImageHeader(&projection00);

        free_page_locked_pixel_type_CUDA(projection_data);
        cvReleaseImageHeader(&projection);

        projection_series_rotate_2D_CUDA_filter_1D_rotate_2D_CUDA_clear_GPU(deviceID);

#else // PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        free_imageData_cvReleaseImageHeader(projection00);
        free_imageData_cvReleaseImageHeader(projection);

#endif // PROJECTION_SERIES_ROTATE_2D_CUDA_FILTER_1D_ROTATE_2D_CUDA_ENABLE_CUDA

        filter_1D_data_release(filter_data);

        // END: Free resources.

#ifdef MEX
        mexEvalString("fprintf('Done processing projection.\\n');");
#endif
    }

    calc_mmm_of_MRC(&file_header_out);

    // Free resources.
    rotate_2D_data_release(rotate_2D_data);

    if (fclose(file_in))
        ABORT("Cannot close input file: \"%s\".\n", filesystem_params->filepath_in);

    if (fclose(file_out))
        ABORT("Cannot close output file \"%s\".\n", filesystem_params->filepath_out);

    t_1_total = clock();
 
    Print("Total duration: %.15e\n", secs_per_clock * ((double) t_1_total - (double) t_0_total));
}

#endif // FILTER_1D_ENABLE_CUDA

// BEGIN: TxBR 2.0 and TxBR 3.0 code. {
// NOTE: Uses IplImage for workspace.  Pulls w/out index check.
void projection_series_rotate_2D_filter_1D_rotate_2D
(
    filesystem_params_handle filesystem_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    rotate_2D_params_handle rotate_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
) 
{
// NOTE: Does NOT use PROTOTYPE_COMPLIANT_INDEXING.
//#if PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif
// NOTE: Does NOT use PROTOTYPE_COMPLIANT_ROTATION.
//#if PROTOTYPE_COMPLIANT_ROTATION
//    ABORT("PROTOTYPE_COMPLIANT_ROTATION not implemented!\n");
//#endif

    // BEGIN: Check parameter structs.
    if (filesystem_params == NULL)
        ABORT("filesystem_params == NULL.\n");

    if (symmetrize_2D_params == NULL)
        symmetrize_2D_params = create_symmetrize_2D_params();

    if (rotate_2D_params == NULL)
        ABORT("rotate_2D_params == NULL.\n");

    if (filter_1D_params == NULL)
        filter_1D_params = create_filter_1D_params();

    if (edgetaper_2D_params == NULL)
        edgetaper_2D_params = create_edgetaper_2D_params();
    // END: Check parameter structs.

// PRIMARY DEBUG:     Print("filesystem_params->filepath_in = \"%s\"\n", filesystem_params->filepath_in); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("filesystem_params->filepath_out = \"%s\"\n", filesystem_params->filepath_out); // PRIMARY DEBUG.

    clock_t t_0_total, t_1_total;
    clock_t t_0, t_1;

    double secs_per_clock = 1.0 / (double) CLOCKS_PER_SEC;

    t_0_total = clock();

    // BEGIN: Open input MRC file.
    MrcHeader file_header_in;

    FILE *file_in = NULL;

    openMRCFile_general(filesystem_params->filepath_in, &file_header_in, &file_in);
    // END: Open input MRC file.

    int n0_x = file_header_in.nx;
    int n0_y = file_header_in.ny;
    int n_projections = file_header_in.nz;
    //int n_projections = 1; WARN("n_projections hardcoded to %i.\n", n_projections);
    int mode = file_header_in.mode;

    // BEGIN: Create and open output MRC file.
    // BEGIN: Create output MRC file.
    if (mode != SLICE_MODE_FLOAT)
        WARN("Input mode == %i but using mode SLICE_MODE_FLOAT == %i for output!\n", mode, SLICE_MODE_FLOAT);

    MrcHeader file_header_out;

    FILE *file_out = NULL;

    createNewMRCFile(filesystem_params->filepath_out, n0_x, n0_y, filesystem_params->n_projections_indexed, /*mode*/ SLICE_MODE_FLOAT);
    // END: Create output MRC file.

    // BEGIN: Open output MRC file.
    openMRCFile_general(filesystem_params->filepath_out, &file_header_out, &file_out);
    // END: Open output MRC file.
    // END: Create and open output MRC file.

    rotate_2D_data_type *rotate_2D_data = create_rotate_2D_data(rotate_2D_params, filesystem_params->n_projections_indexed);

    //for (int i_projection = 0; i_projection < n_projections; i_projection++)
    for (int i_index = 0; i_index < filesystem_params->n_projections_indexed; ++i_index)
    {
        int i_projection = filesystem_params->projection_indices[i_index];
        Print("i_projection = %i\n", i_projection);

        if ((i_projection < 0) || (i_projection >= n_projections))
            ABORT("i_projection == %i, but n_projections == %i.\n", i_projection, n_projections);

        t_0 = clock();

        real_type rot_matrix_2x2[2][2];
        real_type inv_rot_matrix_2x2[2][2];

        calc_rotation_matrices_current(rotate_2D_data, i_index, rot_matrix_2x2, inv_rot_matrix_2x2);

        int n_neg_row_pad = 2;
        int n_pos_row_pad = 2;
        int n_neg_col_pad = 2;
        int n_pos_col_pad = 2;

        int n_neg_row_extend = 2;
        int n_pos_row_extend = 2;
        int n_neg_col_extend = 2;
        int n_pos_col_extend = 2;

        // NOTE: Treat pixels as points.
        rect_2DR_type support0 = { 0.0, 0.0, (real_type) (n0_x - 1), (real_type) (n0_y - 1) };
        rect_2DR_type support0_extended =
        {
            support0.x - (real_type) n_neg_row_extend,
            support0.y - (real_type) n_neg_col_extend,
            (real_type) (INT_FROM_REAL_IMAGE_EXTENT(support0.width) + n_neg_row_extend + n_pos_row_extend),
            (real_type) (INT_FROM_REAL_IMAGE_EXTENT(support0.height) + n_neg_col_extend + n_pos_col_extend)
        };

        rect_2DR_type support = support0;

        Print("support = rect_2DR_rotate_2D(rot_matrix_2x2, support0)\n");
        rect_2DR_rotate_2D(rot_matrix_2x2, &support);

        rect_2DR_type support_extended = 
        {
            support.x - (real_type) n_neg_row_extend,
            support.y - (real_type) n_neg_col_extend,
            (real_type) (INT_FROM_REAL_IMAGE_EXTENT(support.width) + n_neg_row_extend + n_pos_row_extend),
            (real_type) (INT_FROM_REAL_IMAGE_EXTENT(support.height) + n_neg_col_extend + n_pos_col_extend)
        };

        rect_2DR_type support00 = support;
        rect_2DR_type support00_extended = support_extended;

        Print("support00 = rect_2DR_rotate_2D(inv_rot_matrix_2x2, support)\n");
        rect_2DR_rotate_2D(inv_rot_matrix_2x2, &support00);
        Print("support00_extended = rect_2DR_rotate_2D(inv_rot_matrix_2x2, support_extended)\n");
        rect_2DR_rotate_2D(inv_rot_matrix_2x2, &support00_extended);

        int n_x = INT_FROM_REAL_IMAGE_EXTENT(support.width) + 1;
        int n_y = INT_FROM_REAL_IMAGE_EXTENT(support.height) + 1;
        //int n00_x = INT_FROM_REAL_IMAGE_EXTENT(support00.width) + 1;
        //int n00_y = INT_FROM_REAL_IMAGE_EXTENT(support00.height) + 1;

        //int n0_x_extended = INT_FROM_REAL_IMAGE_EXTENT(support0_extended.width) + 1;
        ////int n0_y_extended = INT_FROM_REAL_IMAGE_EXTENT(support0_extended.height) + 1;
        int n0_x_extended = n0_x + (n_neg_row_extend + n_pos_row_extend);
        int n0_y_extended = n0_y + (n_neg_col_extend + n_pos_col_extend);
        int n_x_extended = INT_FROM_REAL_IMAGE_EXTENT(support_extended.width) + 1;
        int n_y_extended = INT_FROM_REAL_IMAGE_EXTENT(support_extended.height) + 1;
        int n00_x_extended = INT_FROM_REAL_IMAGE_EXTENT(support00_extended.width) + 1;
        int n00_y_extended = INT_FROM_REAL_IMAGE_EXTENT(support00_extended.height) + 1;

        int n00_x_extended_padded = n00_x_extended + n_neg_row_pad + n_pos_row_pad;
        int n00_y_extended_padded = n00_y_extended + n_neg_col_pad + n_pos_col_pad;

        // Ensure that the parities of the extents of projection00 match the parities of the extents of the source image.
        if ((n0_x % 2) != (n00_x_extended_padded % 2)) ++n00_x_extended_padded;
        if ((n0_y % 2) != (n00_y_extended_padded % 2)) ++n00_y_extended_padded;

        int n_x_extended_padded = n_x_extended + n_neg_row_pad + n_pos_row_pad;
        int n_y_extended_padded = n_y_extended + n_neg_col_pad + n_pos_col_pad;

        // Ensure that the parities of the extents of projection match the parities of the extents of projection00.
        if ((n00_x_extended_padded % 2) != (n_x_extended_padded % 2)) ++n_x_extended_padded;
        if ((n00_y_extended_padded % 2) != (n_y_extended_padded % 2)) ++n_y_extended_padded;

        // Calculate filter.
        filter_1D_data_type *filter_data = create_filter_1D_data_n_xy(filter_1D_params, n_x, n_y);

        IplImage *projection00 = create_IplImage_SetZero_malloc(n00_x_extended_padded, n00_y_extended_padded, IPL_DEPTH_REAL); // NOTE: Used for both original projection and untransformed filtered transformed projection.
        IplImage *projection = create_IplImage_SetZero_malloc(n_x_extended_padded, n_y_extended_padded, IPL_DEPTH_REAL); // NOTE: Used for both transformed projection and filtered transformed projection.

        center_ROI_params_type center_params;

        center_params = calc_center_ROI_params(n00_x_extended_padded, n00_y_extended_padded, n0_x, n0_y);
        CvRect projection00_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                       center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection00_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                                projection00, 
                                                projection00_center_ROI_CvRect);

        center_params = calc_center_ROI_params(n00_x_extended_padded, n00_y_extended_padded, n0_x_extended, n0_y_extended);
        CvRect projection00_extended_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                                center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection00_extended_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                                         projection00, 
                                                         projection00_extended_center_ROI_CvRect);

        center_params = calc_center_ROI_params(n_x_extended_padded, n_y_extended_padded, n_x, n_y);
        CvRect projection_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                     center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                              projection, 
                                              projection_center_ROI_CvRect);

        center_params = calc_center_ROI_params(n_x_extended_padded, n_y_extended_padded, n_x_extended, n_y_extended);
        CvRect projection_extended_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                              center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection_extended_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                                       projection, 
                                                       projection_extended_center_ROI_CvRect);

//        IplImage *projection_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
//                                              projection, 
//                                              cvRect(0, 0, n_x, n_y));

        // Reuse slice later to transfer data back to the MRC file.
        Islice *slice = NULL;
        slice = sliceReadMRC(&file_header_in, i_projection, 'z');
        if (slice == NULL)
            ABORT("Cannot read file: \"%s\".\n", filesystem_params->filepath_in);

        // NOTE: Output floating-point data.
        if (slice->mode != SLICE_MODE_FLOAT)
            sliceNewMode(slice, SLICE_MODE_FLOAT);

        copy_Islice_to_IplImage_center(slice, projection00_center_ROI);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection00_center_ROI, "/home/akulo/filter_1D_images/projection00_center_ROI.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Initialization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        if (symmetrize_2D_params->symmetrize_2D_flag)
        {
            t_0 = clock();

            symmetrize_IplImage_values_CvRect(projection00, projection00_center_ROI_CvRect);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.symmetrized.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection00_center_ROI, "/home/akulo/filter_1D_images/projection00_center_ROI.symmetrized.mrc"); // PRIMARY DEBUG.

            t_1 = clock();
// PRIMARY DEBUG:             Print("Symmetrization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.
        }

        t_0 = clock();

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
        projection0_global = projection00;
        projection0_center_ROI_CvRect_global = &projection00_center_ROI_CvRect;
        projection_global = projection;
        projection_center_ROI_CvRect_global = &projection_extended_center_ROI_CvRect;
#endif

        projection_rotate_2D(
            projection00_center_ROI, &support0,
            rotate_2D_params->local_mag,
            rot_matrix_2x2,
            projection_extended_center_ROI, &support_extended);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection_center_ROI, "/home/akulo/filter_1D_images/projection_center_ROI.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        t_0 = clock();

        projection_filter_1D_fs(projection_center_ROI, filter_data, symmetrize_2D_params);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.filtered.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection_center_ROI, "/home/akulo/filter_1D_images/projection_center_ROI.filtered.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Filtering delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        if (symmetrize_2D_params->symmetrize_2D_flag)
        {
            t_0 = clock();

            symmetrize_IplImage_values_CvRect(projection, projection_center_ROI_CvRect);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.symmetrized.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection_center_ROI, "/home/akulo/filter_1D_images/projection_center_ROI.symmetrized.mrc"); // PRIMARY DEBUG.

            t_1 = clock();
// PRIMARY DEBUG:             Print("Symmetrization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.
        }
        else
        {
            t_0 = clock();

            cvSetZero_IplImage_values_outside_CvRect(projection, projection_center_ROI_CvRect);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.cleared.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection_center_ROI, "/home/akulo/filter_1D_images/projection_center_ROI.cleared.mrc"); // PRIMARY DEBUG.

            t_1 = clock();
// PRIMARY DEBUG:             Print("Post-filter clearing delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.
        }

        t_0 = clock();

#if PROJECTION_ROTATE_2D_CHECK_BOUNDS
        projection0_global = projection;
        projection0_center_ROI_CvRect_global = &projection_center_ROI_CvRect;
        projection_global = projection00;
        projection_center_ROI_CvRect_global = &projection00_extended_center_ROI_CvRect;
#endif

        projection_rotate_2D(
            projection_center_ROI, &support,
            rotate_2D_params->local_mag,
            inv_rot_matrix_2x2,
            projection00_extended_center_ROI, &support0_extended);

        t_1 = clock();
// PRIMARY DEBUG:         Print("Inv. transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.filtered.mrc"); // PRIMARY DEBUG.

        // BEGIN: Edgetapering.
        if (edgetaper_2D_params->width > 0)
        {
            edgetaper_2D(projection00_center_ROI, edgetaper_2D_params->width);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.filtered.edgetapered.mrc"); // PRIMARY DEBUG.
        }
        // END: Edgetapering.

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection00_center_ROI, "/home/akulo/filter_1D_images/projection00_center_ROI.filtered.mrc"); // PRIMARY DEBUG.

        copy_IplImage_to_Islice_center(projection00_center_ROI, slice);

        sliceWriteMRC(&file_header_out, slice, i_index, 'z');
        sliceFree(slice);

        // Free resources.
        cvReleaseImageHeader(&projection00_center_ROI);
        cvReleaseImageHeader(&projection00_extended_center_ROI);
        cvReleaseImageHeader(&projection_center_ROI);
        cvReleaseImageHeader(&projection_extended_center_ROI);
        free_imageData_cvReleaseImageHeader(projection00);
        free_imageData_cvReleaseImageHeader(projection);

        filter_1D_data_release(filter_data);

        Print("Done processing projection, i_projection = %i.\n", i_projection);

#ifdef MEX
        mexEvalString("fprintf('Done processing projection.\\n');");
#endif
    }

    calc_mmm_of_MRC(&file_header_out);

    // Free resources.
    rotate_2D_data_release(rotate_2D_data);

    if (fclose(file_in))
        ABORT("Cannot close input file: \"%s\".\n", filesystem_params->filepath_in);

    if (fclose(file_out))
        ABORT("Cannot close output file \"%s\".\n", filesystem_params->filepath_out);

    t_1_total = clock();
 
    Print("Total duration: %.15e\n", secs_per_clock * ((double) t_1_total - (double) t_0_total));
}

// NOTE: Uses IplImage for workspace.  Pulls w/out index check.
void projection_series_rotate_2D_prototype_filter_1D_rotate_2D_prototype
(
    filesystem_params_handle filesystem_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    rotate_2D_params_handle rotate_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
) 
{
// NOTE: Does NOT use PROTOTYPE_COMPLIANT_INDEXING.
//#if PROTOTYPE_COMPLIANT_INDEXING
//    ABORT("PROTOTYPE_COMPLIANT_INDEXING not implemented!\n");
//#endif

    // BEGIN: Check parameter structs.
    if (filesystem_params == NULL)
        ABORT("filesystem_params == NULL.\n");

    if (symmetrize_2D_params == NULL)
        symmetrize_2D_params = create_symmetrize_2D_params();

    if (rotate_2D_params == NULL)
        ABORT("rotate_2D_params == NULL.\n");

    if (filter_1D_params == NULL)
        filter_1D_params = create_filter_1D_params();

    if (edgetaper_2D_params == NULL)
        edgetaper_2D_params = create_edgetaper_2D_params();
    // END: Check parameter structs.

// PRIMARY DEBUG:     Print("filesystem_params->filepath_in = \"%s\"\n", filesystem_params->filepath_in); // PRIMARY DEBUG.
// PRIMARY DEBUG:     Print("filesystem_params->filepath_out = \"%s\"\n", filesystem_params->filepath_out); // PRIMARY DEBUG.

    clock_t t_0_total, t_1_total;
    clock_t t_0, t_1;

    double secs_per_clock = 1.0 / (double) CLOCKS_PER_SEC;

    t_0_total = clock();

    // BEGIN: Open input MRC file.
    MrcHeader file_header_in;

    FILE *file_in = NULL;

    openMRCFile_general(filesystem_params->filepath_in, &file_header_in, &file_in);
    // END: Open input MRC file.

    int n0_x = file_header_in.nx;
    int n0_y = file_header_in.ny;
    int n_projections = file_header_in.nz;
    //int n_projections = 1; WARN("n_projections hardcoded to %i.\n", n_projections);
    int mode = file_header_in.mode;

    // BEGIN: Create and open output MRC file.
    // BEGIN: Create output MRC file.
    if (mode != SLICE_MODE_FLOAT)
        WARN("Input mode == %i but using mode SLICE_MODE_FLOAT == %i for output!\n", mode, SLICE_MODE_FLOAT);

    MrcHeader file_header_out;

    FILE *file_out = NULL;

    createNewMRCFile(filesystem_params->filepath_out, n0_x, n0_y, filesystem_params->n_projections_indexed, /*mode*/ SLICE_MODE_FLOAT);
    // END: Create output MRC file.

    // BEGIN: Open output MRC file.
    openMRCFile_general(filesystem_params->filepath_out, &file_header_out, &file_out);
    // END: Open output MRC file.
    // END: Create and open output MRC file.

    rotate_2D_data_type *rotate_2D_data = create_rotate_2D_data(rotate_2D_params, filesystem_params->n_projections_indexed);

    //for (int i_projection = 0; i_projection < n_projections; i_projection++)
    for (int i_index = 0; i_index < filesystem_params->n_projections_indexed; ++i_index)
    {
        int i_projection = filesystem_params->projection_indices[i_index];
// PRIMARY DEBUG:         Print("i_projection = %i\n", i_projection); // PRIMARY DEBUG.

        if ((i_projection < 0) || (i_projection >= n_projections))
            ABORT("i_projection == %i, but n_projections == %i.\n", i_projection, n_projections);

        t_0 = clock();

        real_type rot_matrix_2x2[2][2];
        real_type inv_rot_matrix_2x2[2][2];

        calc_rotation_matrices_current(rotate_2D_data, i_index, rot_matrix_2x2, inv_rot_matrix_2x2);

        int row_pad = 2;
        int col_pad = 2;

/*
        //int row_pad, col_pad;
        if (test_rotation_matrix_2x2_for_subsampling_threshold(rot_matrix_2x2))
        { 
            //row_pad = 1; col_pad = 1; 
        }
        else
        { 
            test_rotation_matrix_2x2_for_subsampling_threshold(inv_rot_matrix_2x2);
            //row_pad = 2; col_pad = 2; 
        }
        print_matrix_2x2("rot_matrix_2x2 (thresholded)", rot_matrix_2x2);
        print_matrix_2x2("inv_rot_matrix_2x2 (thresholded)", inv_rot_matrix_2x2);
*/

        real_type support0[4];
        support0[0] = ((real_type) n0_x - 1.0) / 2.0; // y-coordinate of projection center
        support0[1] = ((real_type) n0_y - 1.0) / 2.0; // y-coordinate of projection center
        support0[2] = (real_type) n0_x - 1.0; // max x-coordinate of projection
        support0[3] = (real_type) n0_y - 1.0; // max y-coordinate of projection

        real_type support[4];
        support[0] = support0[0]; // x-coordinate of transformed projection center
        support[1] = support0[1]; // y-coordinate of transformed projection center
        support[2] = support0[2]; // max x-coordinate of transformed projection
        support[3] = support0[3]; // max y-coordinate of transformed projection

#if PROTOTYPE_COMPLIANT_ROTATION
        // NOTE: I think the transform should be the same as that used to transform the original projection, but cf. the MATLAB prototype.
        calc_rotated_support(rot_matrix_2x2, support);
#else
        calc_rotated_support(inv_rot_matrix_2x2, support);
#endif

        real_type support00[4];
        support00[0] = support[0]; // x-coordinate of untransformed projection center
        support00[1] = support[1]; // y-coordinate of untransformed projection center
        support00[2] = support[2]; // max x-coordinate of untransformed projection
        support00[3] = support[3]; // max y-coordinate of untransformed projection

#if PROTOTYPE_COMPLIANT_ROTATION
        // NOTE: I think the transform should be the same as that used to transform the transformed projection, but cf. the MATLAB prototype.
        calc_rotated_support(inv_rot_matrix_2x2, support00);
#else
        calc_rotated_support(rot_matrix_2x2, support00);
#endif

        int n_x = round(support[2] + 1.0);
        int n_y = round(support[3] + 1.0);
        int n00_x = round(support00[2] + 1.0);
        int n00_y = round(support00[3] + 1.0);

        int n00_x_padded = n00_x + row_pad;
        int n00_y_padded = n00_y + col_pad;

        // Ensure that the parities of the extents of projection00 match the parities of the extents of the source image.
        if ((n0_x % 2) != (n00_x_padded % 2)) ++n00_x_padded;
        if ((n0_y % 2) != (n00_y_padded % 2)) ++n00_y_padded;

        int n_x_padded = n_x + row_pad;
        int n_y_padded = n_y + col_pad;
        //int n_x_padded = n_x + 1;
        //int n_y_padded = n_y + 1;

        // Ensure that the parities of the extents of projection match the parities of the extents of projection00.
        if ((n00_x_padded % 2) != (n_x_padded % 2)) ++n_x_padded;
        if ((n00_y_padded % 2) != (n_y_padded % 2)) ++n_y_padded;

        // Calculate filter.
        filter_1D_data_type *filter_data = create_filter_1D_data_n_xy(filter_1D_params, n_x_padded, n_y_padded);

        IplImage *projection00 = create_IplImage_SetZero_malloc(n00_x_padded, n00_y_padded, IPL_DEPTH_REAL); // NOTE: Used for both original projection and untransformed filtered transformed projection.
        IplImage *projection = create_IplImage_SetZero_malloc(n_x_padded, n_y_padded, IPL_DEPTH_REAL); // NOTE: Used for both transformed projection and filtered transformed projection.

        center_ROI_params_type center_params;

        center_params = calc_center_ROI_params(n00_x_padded, n00_y_padded, n0_x, n0_y);
        CvRect projection00_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                       center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection00_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                                projection00, 
                                                projection00_center_ROI_CvRect);

        center_params = calc_center_ROI_params(n_x_padded, n_y_padded, n_x, n_y);
        CvRect projection_center_ROI_CvRect = cvRect(center_params.i0_x_start, center_params.i0_y_start,
                                                     center_params.n_i0_x, center_params.n_i0_y);
        IplImage *projection_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
                                              projection, 
                                              projection_center_ROI_CvRect);

//        IplImage *projection_center_ROI = create_IplImage_ROI_widthStep_of_IplImage(
//                                              projection, 
//                                              cvRect(0, 0, n_x, n_y));

        // Reuse slice later to transfer data back to the MRC file.
        Islice *slice = NULL;
        slice = sliceReadMRC(&file_header_in, i_projection, 'z');
        if (slice == NULL)
            ABORT("Cannot read file: \"%s\".\n", filesystem_params->filepath_in);

        // NOTE: Output floating-point data.
        if (slice->mode != SLICE_MODE_FLOAT)
            sliceNewMode(slice, SLICE_MODE_FLOAT);

        copy_Islice_to_IplImage_center(slice, projection00);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection00_center_ROI, "/home/akulo/filter_1D_images/projection00_center_ROI.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Initialization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        if (symmetrize_2D_params->symmetrize_2D_flag)
        {
            t_0 = clock();

            symmetrize_IplImage_values_CvRect(projection00, projection00_center_ROI_CvRect);

            t_1 = clock();
// PRIMARY DEBUG:             Print("Symmetrization delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.
        }

        t_0 = clock();

        projection_rotate_2D_prototype(
            projection00_center_ROI, support0,
            rotate_2D_params->local_mag,
#if PROTOTYPE_COMPLIANT_ROTATION
            inv_rot_matrix_2x2, // NOTE: inv_rot_matrix_2x2 to conform to MATLAB prototype.
#else
            rot_matrix_2x2,
#endif
            projection_center_ROI, support);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection_center_ROI, "/home/akulo/filter_1D_images/projection_center_ROI.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        t_0 = clock();

        projection_filter_1D_fs(projection_center_ROI, filter_data, symmetrize_2D_params);

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection, "/home/akulo/filter_1D_images/projection.filtered.mrc"); // PRIMARY DEBUG.
// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection_center_ROI, "/home/akulo/filter_1D_images/projection_center_ROI.filtered.mrc"); // PRIMARY DEBUG.

        t_1 = clock();
// PRIMARY DEBUG:         Print("Filtering delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

        t_0 = clock();

        projection_rotate_2D_prototype(
            projection_center_ROI, support,
            rotate_2D_params->local_mag,
#if PROTOTYPE_COMPLIANT_ROTATION
            rot_matrix_2x2, // NOTE: rot_matrix_2x2 to conform to MATLAB prototype.
#else
            inv_rot_matrix_2x2,
#endif
            projection00_center_ROI, support0);

        t_1 = clock();
// PRIMARY DEBUG:         Print("Inv. transform delta t:%.15e\n", secs_per_clock * ((double) t_1 - (double) t_0)); // PRIMARY DEBUG.

// PRIMARY DEBUG:         write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.filtered.mrc"); // PRIMARY DEBUG.

        // BEGIN: Edgetapering.
        if (edgetaper_2D_params->width > 0)
        {
            edgetaper_2D(projection00_center_ROI, edgetaper_2D_params->width);

// PRIMARY DEBUG:             write_MRC_image_from_IplImage(projection00, "/home/akulo/filter_1D_images/projection00.filtered.edgetapered.mrc"); // PRIMARY DEBUG.
        }
        // END: Edgetapering.

        copy_IplImage_to_Islice_center(projection00_center_ROI, slice);

        sliceWriteMRC(&file_header_out, slice, i_index, 'z');
        sliceFree(slice);

        // Free resources.
        cvReleaseImageHeader(&projection_center_ROI);
        cvReleaseImageHeader(&projection00_center_ROI);
        free_imageData_cvReleaseImageHeader(projection00);
        free_imageData_cvReleaseImageHeader(projection);

        filter_1D_data_release(filter_data);

#ifdef MEX
        mexEvalString("fprintf('Done processing projection.\\n');");
#endif
    }

    calc_mmm_of_MRC(&file_header_out);

    // Free resources.
    rotate_2D_data_release(rotate_2D_data);

    if (fclose(file_in))
        ABORT("Cannot close input file: \"%s\".\n", filesystem_params->filepath_in);

    if (fclose(file_out))
        ABORT("Cannot close output file \"%s\".\n", filesystem_params->filepath_out);

    t_1_total = clock();
 
    Print("Total duration: %.15e\n", secs_per_clock * ((double) t_1_total - (double) t_0_total));
}

// NOTE: A wrapper to make name "consistent" with projection_series_remap_2D_filter_1D_inv_remap_2D.
// NOTE: Original name is already consistent, because implementation of remap_2D is different from inv_remap_2D.
void projection_series_rotate_2D_filter_1D_inv_rotate_2D
(
    filesystem_params_handle filesystem_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    rotate_2D_params_handle rotate_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
) 
{
    projection_series_rotate_2D_filter_1D_rotate_2D(
    //projection_series_rotate_2D_prototype_filter_1D_rotate_2D_prototype(
        filesystem_params,
        symmetrize_2D_params,
        rotate_2D_params,
        filter_1D_params,
        edgetaper_2D_params);
}
// END: TxBR 2.0 and TxBR 3.0 code. }

void projection_series_filter_1D
(
    filesystem_params_handle filesystem_params,
    symmetrize_2D_params_handle symmetrize_2D_params,
    remap_2D_params_handle remap_2D_params,
    rotate_2D_params_handle rotate_2D_params,
    filter_1D_params_handle filter_1D_params,
    edgetaper_2D_params_handle edgetaper_2D_params
) 
{
    // BEGIN: Check parameter structs.
    if (filesystem_params == NULL)
        ABORT("filesystem_params == NULL.\n");

    if (symmetrize_2D_params == NULL)
        symmetrize_2D_params = create_symmetrize_2D_params();

    if (remap_2D_params == NULL)
        remap_2D_params = create_remap_2D_params();

    if (rotate_2D_params == NULL)
        rotate_2D_params = create_rotate_2D_params();

    if (filter_1D_params == NULL)
        filter_1D_params = create_filter_1D_params();

    if (edgetaper_2D_params == NULL)
        edgetaper_2D_params = create_edgetaper_2D_params();
    // END: Check parameter structs.

// PRIMARY DEBUG:     if (PROTOTYPE_COMPLIANT_INDEXING == 1) // PRIMARY DEBUG.
// PRIMARY DEBUG:         Print("PROTOTYPE_COMPLIANT_INDEXING == 1\n"); // PRIMARY DEBUG.
// PRIMARY DEBUG:     else // PRIMARY DEBUG.
// PRIMARY DEBUG:         Print("PROTOTYPE_COMPLIANT_INDEXING == %d != 1\n", PROTOTYPE_COMPLIANT_INDEXING); // PRIMARY DEBUG.

    if (remap_2D_params->map_2D_order == 0)
    {
        // NOTE: Pure rotation, filter, and inv. rotation without boundary check.

#ifdef FILTER_1D_ENABLE_CUDA

        general_CUDA_params_handle general_CUDA_params = create_general_CUDA_params_from_data_copy(8, 8, 1);

        projection_series_rotate_2D_CUDA_filter_1D_rotate_2D_CUDA(
            filesystem_params,
            general_CUDA_params,
            symmetrize_2D_params,
            rotate_2D_params,
            filter_1D_params,
            edgetaper_2D_params);

#else // FILTER_1D_ENABLE_CUDA

// BEGIN: TxBR 2.0 and TxBR 3.0 code. {
        projection_series_rotate_2D_filter_1D_inv_rotate_2D(
            filesystem_params,
            symmetrize_2D_params,
            rotate_2D_params,
            filter_1D_params,
            edgetaper_2D_params);
// END: TxBR 2.0 and TxBR 3.0 code. }

#endif // FILTER_1D_ENABLE_CUDA

    }
    else
    {
        // BEGIN: Single remap function.
#ifdef FILTER_1D_ENABLE_CUDA

        // BEGIN: DEBUG.
        projection_series_transform_2D_filter_1D_inv_transform_2D(
            filesystem_params,
            symmetrize_2D_params,
            remap_2D_params,
            rotate_2D_params,
            filter_1D_params,
            edgetaper_2D_params);
        // END: DEBUG.

        general_CUDA_params_handle general_CUDA_params = create_general_CUDA_params_from_data_copy(16, 16, 1);

        projection_series_transform_2D_CUDA_filter_1D_inv_transform_2D_CUDA(
            filesystem_params,
            general_CUDA_params,
            symmetrize_2D_params,
            remap_2D_params,
            rotate_2D_params,
            filter_1D_params,
            edgetaper_2D_params);

#else // FILTER_1D_ENABLE_CUDA

/*
        projection_series_transform_2D_filter_1D_inv_transform_2D(
            filesystem_params,
            symmetrize_2D_params,
            remap_2D_params,
            filter_1D_params,
            edgetaper_2D_params);
*/

// BEGIN: TxBR 2.0 and TxBR 3.0 code. {
        projection_series_remap_2D_filter_1D_inv_remap_2D(
            filesystem_params,
            symmetrize_2D_params,
            remap_2D_params,
            filter_1D_params,
            edgetaper_2D_params);
// END: TxBR 2.0 and TxBR 3.0 code. }

#endif // FILTER_1D_ENABLE_CUDA
        // END: Single remap function.
    
/*
        // BEGIN: Separate remap and rotate functions.
        projection_series_remap_2D_rotate_2D_filter_1D_rotate_2D_restore_2D(
            filesystem_params,
            symmetrize_2D_params,
            remap_2D_params,
            rotate_2D_params,
            filter_1D_params,
            edgetaper_2D_params); 
        // END: Separate remap and rotate functions.
*/
    }
}

////////////////////////////////////////////////////////////////////////////////
// END: projection_series_filter_1D(). }
////////////////////////////////////////////////////////////////////////////////
