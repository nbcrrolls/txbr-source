#ifndef TXBR_H
#define TXBR_H

#define PROJECTIVE_ORDER 1
#define QUADRATIC_ORDER 2
#define CUBIC_ORDER 3

#define X_MAX	8192	/*  Maximum size in the X direction	*/
#define Y_MAX	8192	/*  Maximum size in the Y direction	*/
#define BLOCK_SIZE_MAX	250  /*	Maximum Number of blocks	*/
#define MAX_TILT	2000  /*	Maximum Number of tilts	 - Enforced*/
#define MAX_BLOCK_VOL	500000000  /*   Maximum Volume In Pixels For A Block - Enforced */

#define MAX_ORDER	20  /*	Maximum order for the projection polynomial approximation   */
#define MAX_COEFFICIENTS_NUMBER	1771  /*    Maximum order for the projection polynomial approximation   */
#define FILENAME_LEN	256	/*  Maximum length for a file name  */

typedef struct {
    char basename[FILENAME_LEN];
    char directory[FILENAME_LEN];
    int blocksize;
    int numberOfBlocks;
    int use_fixed_segment_size;   /* Use fix segments in the polynomial evaluation trick. */
    int error;
    double x_0, y_0, z_0; /*	x_0 and y_0 should be more than 1	*/
    double x_1, y_1, z_1; /*	x_1 and y_1 should be less than nx and ny	*/
    double x_inc, y_inc, z_inc;
    double sx, sy, sz;
    double z_start;
    double plane_coeffs1[3], plane_coeffs2[3];
} TxBRsetup;

typedef struct {
    TxBRsetup *setup;
    int indexOfBlock_x, indexOfBlock_y, indexOfBlock_z;
    double x_start, y_start, z_start;
    double x_stop, y_stop, z_stop;
} Volume;

typedef struct {
    char label[FILENAME_LEN];
    int order;
    int number_of_terms;
    int numberOfTilts;
    int skipView[MAX_TILT];
    int itilt_start, itilt_stop;
    int nx, ny;
    double lambda[MAX_TILT][4];
    double coefficients_1[MAX_TILT][MAX_COEFFICIENTS_NUMBER];
    double coefficients_2[MAX_TILT][MAX_COEFFICIENTS_NUMBER];
    double coefficients_3[MAX_TILT][MAX_COEFFICIENTS_NUMBER];
    int error;
} PathProjection;

typedef struct {
    int order;
    int number_of_terms;
    int numberOfTilts;
    int skipView[MAX_TILT];
    double lambda[MAX_TILT][4];
    double coefficients_1[MAX_TILT][MAX_COEFFICIENTS_NUMBER];
    double coefficients_2[MAX_TILT][MAX_COEFFICIENTS_NUMBER];
    int error;
} TrajectoryMap;

typedef struct {
    float origin;
    float* data;
    int tiltsize;
    int rsize;
} Ifan;

#endif
