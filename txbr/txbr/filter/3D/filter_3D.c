#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "txbrutil.h"
#include "mrcslice.h"

Ifan *fanCreate(int tiltsize, int rsize);
int fanInit(Ifan *f, int tiltsize, int rsize, void *data);
int fanFree(Ifan *f);
void fanClear(Ifan *s, float val);
int fanGetTiltSize(Ifan *fan);
int fanGetRadialSize(Ifan *fan);
int fanGetVal(Ifan *s, int x, int y, float* val);
int fanPutVal(Ifan *s, int x, int y, float val);
Ifan *fanReadMRC(MrcHeader *hin, int x, int y, int z, int tiltsize, int rsize, TrajectoryMap *trajectoryMap);
float *mrc_mread_fan(FILE *fin, MrcHeader *hdata, int x, int y, int z, int tiltsize, int rsize, TrajectoryMap *trajectoryMap);
int mrc_read_fan(float *data, FILE *fin, MrcHeader *hdata, int x, int y, int z, int tiltsize, int rsize, TrajectoryMap *trajectoryMap);

void load_trajectory(char* directory, char* basename, TrajectoryMap* trajectory);
void print_trajectory( TrajectoryMap *trajectory );

int number_of_terms_traj;
int order_in_X_traj[MAX_COEFFICIENTS_NUMBER];
int order_in_Y_traj[MAX_COEFFICIENTS_NUMBER];
int order_in_Z_traj[MAX_COEFFICIENTS_NUMBER];

/*
 *
 */
double lambda( TrajectoryMap *trajectoryMap, int index_of_tilt, double x, double y, double z ) {

	int i;

	double value = 0.0;

	for (i=0; i<number_of_terms_traj; i++) {
		value += trajectoryMap->lambda[index_of_tilt][i]*
			pow(x,order_in_X_traj[i])*pow(y,order_in_Y_traj[i])*pow(z,order_in_Z_traj[i]);
	}

	return value;

}

/*
 *
 */
double X( TrajectoryMap *trajectoryMap, int index_of_tilt, double x, double y, double z ) {

	int i;

	double value = 0.0;

	for (i=0; i<number_of_terms_traj; i++) {
		value += trajectoryMap->coefficients_1[index_of_tilt][i]*
			pow(x,order_in_X_traj[i])*pow(y,order_in_Y_traj[i])*pow(z,order_in_Z_traj[i]);
	}

	return value;

}

/*
 *
 */
double Y( TrajectoryMap *trajectoryMap, int index_of_tilt, double x, double y, double z ) {

	int i;

	double value = 0.0;

	for (i=0; i<number_of_terms_traj; i++) {
		value += trajectoryMap->coefficients_1[index_of_tilt][i]*
			pow(x,order_in_X_traj[i])*pow(y,order_in_Y_traj[i])*pow(z,order_in_Z_traj[i]);
	}

	return value;

}

/*
 * Creates an @@Ifan structure@. Allocates a data array of float appropriate for the
 * size (tiltsize,rsize).  Returns a pointer to the fan or the NULL for error
 */
Ifan *fanCreate(int tiltsize, int rsize) {

	Ifan *f;
	size_t fullsize = (size_t)tiltsize * (size_t)rsize;

	if (fullsize/tiltsize!=rsize) {
		txbr_error(stderr, "ERROR: fanCreate - fan is too large.\n");
		return NULL;
	}

	f = (Ifan *)malloc(sizeof(Ifan));
	if (!f) return NULL;

	f->tiltsize = tiltsize;
	f->rsize = rsize;
	f->data = (float *)malloc(tiltsize*rsize*sizeof(float));

	if (!f->data){
		txbr_error(stderr, "ERROR: fanCreate - not enough memory.\n");
	    free(f);
	    return NULL;
	}

	return f;

}

/*
 * Initialize elements in the @@Ifan structure@ [s] with the given size
 * [xsize], [ysize] and mode [mode], and sets the data member to [data].
 * Returns -1 for an undefined mode.
 */
int fanInit(Ifan *f, int tiltsize, int rsize, void *data) {

	f->tiltsize = tiltsize;
	f->rsize = rsize;
	f->data = data;

	return NO_ERROR;
}

/*!
 * Frees the @@Ifan structure@ fan and its data array if any.
 */
int fanFree(Ifan *f) {

	if (!f) return (-1);
	if (f->data) free(f->data);
	free(f);

	return NO_ERROR;

}

/*!
 * Sets the entire fan [s] to the value given by the float [val]
 */
void fanClear(Ifan *f, float val) {

	int i,j;

	f->origin = val;

	for(i=0; i<f->tiltsize; i++)
		for (j=0; j<f->rsize; j++)
			fanPutVal(f,i,j,val);

	return;

}

/*! Returns the Tilt size of [fan] */
int fanGetTiltSize(Ifan *fan) {

	return fan->tiltsize;

}

/*! Returns the Radial size of [fan] */

int fanGetRadialSize(Ifan *fan) {

	return fan->rsize;

}

/*!
 * Gets the value of point at [tiltsize], [rsize] in fan [f] and paases it
 * into the pointer [val].  Returns -1 for a point out of bounds.
 */
int fanGetVal(Ifan *f, int t, int r, float* val) {

	int index = t + (r * f->tiltsize);

	if ( (t < 0) || (r < 0) || (t >= f->tiltsize) || (r >= f->rsize)) {
		return(-1);
	}

	val[0] = f->data[index];

	return NO_ERROR;

}

/*!
 * Puts the value [val] into the pixel at [x], [y] in fan [s].
 * Returns -1 if the point is out of bounds.
 */
int fanPutVal(Ifan *f, int t, int r, float val) {

	if ( (r < 0) || (t < 0) || (r >= f->tiltsize) || (t >= f->rsize)) {
		return -1;
	} else if (r==0) {
		int t,index;
		for (t=0; t<f->tiltsize; t++) {
			index = t * f->rsize;
			f->data[index] = val;
		}
	} else {
		int index = r + (t * f->rsize);
		f->data[index] = val;
	}

	return NO_ERROR;

}

/*!
 * Returns a fan of data at point [x][y][z] from the file described by the
 * @@MrcHeader structure@ [hin]. The file pointer in [hin] is used.
 * Returns NULL for errors.
 */
Ifan *fanReadMRC(MrcHeader *hin, int x, int y, int z, int tiltsize, int rsize, TrajectoryMap *trajectoryMap) {

	Ifan *fan = (Ifan *)malloc(sizeof(Ifan));

	if (!fan){
		txbr_error(stderr, "ERROR: fanReadMRC - couldn't get memory.\n");
		return(NULL);
	}

	float *buf = (float *)mrc_mread_fan(hin->fp, hin, x, y, z, tiltsize, rsize, trajectoryMap);

	if (!buf) {
		free(fan);
		return(NULL);
	}

 	if (fanInit(fan, tiltsize, rsize, buf)) {
 		free(fan);
 		return NULL;
	}

 	return fan;

}

/*!
 * Allocates and returns one fan of data at the coordinate given by [x][y][z].
 * Reads from the file with pointer [fin] according to the header in [hdata].
 * MRC file should be of type float. Returns NULL for errors.
 */
float *mrc_mread_fan(FILE *fin, MrcHeader *hdata, int x, int y, int z, int tiltsize, int rsize, TrajectoryMap *trajectoryMap) {

	float *buf = NULL;

	buf = (float *)malloc(tiltsize*rsize*sizeof(float));

	if (!buf){
		txbr_error(stderr, "ERROR: mrc_mread_fan - couldn't get memory.\n");
		return NULL;
	}

	if (!mrc_read_fan(buf, fin, hdata, x, y, z, tiltsize, rsize, trajectoryMap)) return buf;

 	free(buf);

 	return NULL;

}

/*!
 * Reads a fan of data into the buffer [data] at the coordinate given by
 * [x][y][z]. Reads from the file with pointer [fin] according to the header in
 * [hdata]. Returns -1 for errors.
 */
int mrc_read_fan(float *data, FILE *fin, MrcHeader *hdata, int x, int y, int z, int tiltsize, int rsize, TrajectoryMap *trajectoryMap) {

	int itilt,j,index,nx,ny,nz;
	float value;

	nx = hdata->nx;
	ny = hdata->ny;
	nz = hdata->nz;

	x = x-1;
	y = y-1;
	z = z-1;

	if ( x<0 || x>=nx || y<0 || y>=ny || z<0 || z>=nz ) return -1;

	printf("%i\n",rsize);
	printf("%i\n",tiltsize);

	//float lambda,X,Y;

	index = 0;

	for (itilt=0; itilt<tiltsize; itilt++) {

		for(j=0; j<rsize; j++) {

			x = itilt;
			y = j;
			value = (float)mrc_read_point( fin, hdata, x, y, z );
			data[index++] = value;

		}

	}

	return NO_ERROR;

}

































/*
 * Filter!
 */
int filter( char* directory, char* basename ) {

	/*	Load electron trajectories	*/

	TrajectoryMap *trajectoryMap = (TrajectoryMap *)malloc(sizeof(TrajectoryMap));

	if (!trajectoryMap) {
		txbr_error(stderr, "ERROR: main - couldn't get memory.\n");
		return MEMORY_ERROR;
	}

	char configuration_file_name[FILENAME_LEN];
	sprintf(&configuration_file_name[0],"%s/%s.traj",directory,basename);

	load_trajectory( directory, basename, trajectoryMap );

	print_trajectory( trajectoryMap );

	/*	Do some initialization	*/

	int n_ = trajectoryMap->order;
	number_of_terms_traj = number_of_polynomial_terms(n_);

	init_power_order(n_, &order_in_X_traj[0], &order_in_Y_traj[0], &order_in_Z_traj[0]);

	/*	Do the filtering	*/

	char file_in[FILENAME_LEN], file_out[FILENAME_LEN];

	sprintf(&file_in[0],"%s/%s.preali",&directory[0],&basename[0]);
	sprintf(&file_out[0],"%s/%s.filter3D",&directory[0],&basename[0]);

	printf("File: %s\n",&file_in[0] );

	MrcHeader *header_in  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header_in) {
		txbr_error(stderr,"ERROR: filter - getting memory.\n");
		return MEMORY_ERROR;
	}

	FILE *fin, *fout;

	fin = fopen(&file_in[0],"r");

	if (fin==NULL) {
		txbr_error(stderr, "ERROR: filter - cannot open file %s.\n", &file_in[0]);
		return FILE_ERROR;
	}

	mrc_head_read(fin,header_in);

	int nx = header_in->nx;
	int ny = header_in->ny;
	int numberOfSlices = header_in->nz;
// TEST  TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
numberOfSlices = 1;

	fout = fopen(&file_out[0],"w");

	if (fout==NULL) {
		txbr_error(stderr, "ERROR: filter - cannot open file %s.\n", &file_out[0]);
		return FILE_ERROR;
	}

	/* Creates the MRC header for the main output file	*/

	struct MRCheader header_out;

	mrc_head_new(&header_out,nx,ny,numberOfSlices,SLICE_MODE_FLOAT);

	mrc_head_write(fout, &header_out);

	/* Do the filtering	*/

	Islice *slice = NULL;
	Ifan *fan = NULL;
	int islice;

	for (islice=0;islice<numberOfSlices;islice++) {

		printf("Reading Slice #%i\n",islice);

		fan = fanReadMRC(header_in, 1, 1, islice+1, nx, ny, trajectoryMap);

//		slice = sliceReadMRC(header_in,islice,'z');

//		printf("min %f max %f\n",slice->min,slice->max);

		mrc_write_slice( fan->data, fout, &header_out, islice, 'z' );

		if (slice) sliceFree(slice);

	}

	if (fclose(fin)) {
		txbr_error(stderr, "ERROR: merge_blocks - closing file %s.\n", &file_in[0]);
		return FILE_ERROR;
	}

	if (fclose(fout)) {
		txbr_error(stderr, "ERROR: merge_blocks - closing file %s.\n", &file_out[0]);
		return FILE_ERROR;
	}

	/* Free stuff	*/

	free(header_in);

	return NO_ERROR;

}

/*

int main(int argc, char *argv[]) {

	char *directory = "/Users/sph/Electron_Tomography/txbr-data/stacked_grids";
	char *basename = "x";

	filter( directory, basename );

}
*/