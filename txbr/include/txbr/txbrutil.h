#ifndef TXBRUTIL_H

#define TXBRUTIL_H

#include "txbr.h"
#include "mrcfiles.h"

#define MAX_TXBR_ERROR_STRING  512

#define NO_ERROR 0
#define FILE_ERROR 1
#define MEMORY_ERROR 2
#define PARAMETER_ERROR 3

#define MAX_NUMBER_OF_FILES_TO_MERGE 500

#define TXBR_MIN(a,b) ((a) < (b) ? (a) : (b))
#define TXBR_MAX(a,b) ((a) > (b) ? (a) : (b))

void txbr_set_store_error(int value);
void txbr_error(FILE *out, char *format, ...);
char *txbr_get_error(void);

int openMRCFile(char* directory, char* basename, MrcHeader* header, FILE** fp);

int number_of_polynomial_terms(int order);

void init_power_order(int order, int* order_in_X, int* order_in_Y, int* order_in_Z);
void print_power_order(int order, int* order_in_Axis);

int save_block(char* directory, char* basename, double z0, int indexOfBlock, int nx, int ny, int blocksize, float* block);
int merge_blocks(char* directory, char* basename, double z0, int numberOfBlocks);
int recompose_stack(char* directory, char* filename, int n_1, int n_2, double alpha_1, double alpha_2);

void set_path_projection_coefficients(PathProjection *path, int tilt, int order, int axis, double *coefficients);

void print_TxBR_setup(TxBRsetup *setup);
void print_path_projection(PathProjection *path);

void load_configuration(TxBRsetup* setup, PathProjection* projection, char* configuration_file_name);
int write_configuration(TxBRsetup* setup, PathProjection* projection, char* configuration_file_name);

PathProjection* load_path(char* directory, char* basename, int order);
TxBRsetup* load_TxBR_setup(char* directory, char* basename);

#endif
