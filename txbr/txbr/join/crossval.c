#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mrcfiles.h"
#include "mrcslice.h"
#include "txbrutil.h"

#define MAX_NUMBER_OF_FILES 50

#define NO_CROSS_VALIDATION	100
#define REMOVE_POSITIVE_ARTEFACTS	101
#define REMOVE_NEGATIVE_ARTEFACTS	102

int calculate_moments(char* directory, char* filein, float* min, float* max, float* moments, int n);
int cross_validate_slices(int nx, int ny, int nc, float* data_2D, float* offset, int mode, float* data_c);
int cross_validate_stacks(int nc, char* directory, char** filenames, float* offset, float* deviation, int mode, char* fileout);

/*
 * Calculate Moments
 */
int calculate_moments(char* directory, char* filein, float* min, float* max, float* moments, int n) {

	int error = NO_ERROR;

	/* Open the input MRC file	*/

	MrcHeader input_file_header, output_file_header;
	FILE *input_file;
	Islice* slice = NULL;

	char imodfile_in[FILENAME_LEN];

	sprintf(&imodfile_in[0], "%s/%s", directory, filein);

	input_file = fopen(&imodfile_in[0],"rbe");

	if (input_file==NULL) {
		txbr_error(stderr, "ERROR: calculate_moments - open file %s.\n", &imodfile_in[0]);
		return FILE_ERROR;
	} else {
		printf("File %s was successfully opened.\n", &imodfile_in[0]);
	}

	mrc_head_read(input_file, &input_file_header);

	int nx = input_file_header.nx;
	int ny = input_file_header.ny;
	int nz = input_file_header.nz;

	int itilt, ix, iy, index, imoment;
	float moment1 = 0;

	float *data;

	min[0] = FLT_MAX;
	max[0] = -FLT_MAX;

	for (itilt=0; itilt<nz; itilt++) {

		slice = sliceReadMRC(&input_file_header,itilt,'z');

		data = &slice->data.f[0];

		index=0;

		for (iy=0; iy<ny; iy++) {

			for (ix=0; ix<nx; ix++) {

				min[0] = (min[0] <data[index]) ? min[0] : data[index];
				max[0] = (max[0]<data[index]) ? data[index] : max[0];
				moment1 = moment1 + data[index];

				index++;

			}

		}

		if (slice) sliceFree(slice);

	}

	moment1 = moment1/nx/ny/nz;

	printf("Average Value for %s: %f.\n", &filein[0], moment1);
	printf("Minimal Value for %s: %f.\n", &filein[0], min[0]);
	printf("Maximal Value for %s: %f.\n", &filein[0], max[0]);


	for (imoment=2; imoment<n; imoment++) {

		moments[imoment] = 0.0;

	}


	for (itilt=0; itilt<nz; itilt++) {

		slice = sliceReadMRC(&input_file_header,itilt,'z');

		data = &slice->data.f[0];

		for (imoment=2; imoment<n; imoment++) {

			index=0;

			for (iy=0; iy<ny; iy++) {

				for (ix=0; ix<nx; ix++) {

					moments[imoment] = moments[imoment] + pow(data[index]-moment1,imoment);

					index++;

				}

			}

		}

		if (slice) sliceFree(slice);

	}

	moments[0] = 1;
	moments[1] = moment1;

	for (imoment=2; imoment<n; imoment++) {

		moments[imoment] = moments[imoment]/nx/ny/nz;

	}

	if (n>=2) {
		printf("2nd Moment for %s: %f.\n", &filein[0], moments[2]);
	}

	/* Free the resources	*/

	if (fclose(input_file)) {
		txbr_error(stderr, "ERROR: calculate_moments - closing file %s.\n",input_file);
		return FILE_ERROR;
	}

	return error;

}

/*
 * Combine 2D arrays of data in one array with a given weight
 */
int cross_validate_slices(int nx, int ny, int nc, float* data_2D, float* offset, int mode, float* data_c) {

	int ix, iy, ic, ic1, ic2, index;
	float value, value1, value2, value12, value12_p, value12_m;
	float cross_val = 0;

	switch (mode) {

		case NO_CROSS_VALIDATION :
			cross_val = 0;
			break;

		case REMOVE_POSITIVE_ARTEFACTS :
			cross_val = 1;
			break;

		case REMOVE_NEGATIVE_ARTEFACTS :
			cross_val = -1;
			break;
	}

	index = 0;

	for (iy=0; iy<ny; iy++) {

		for (ix=0; ix<nx; ix++) {

			data_c[index] = 0.0;

			for (ic=0; ic<nc; ic++) {

				value = offset[0] + data_2D[ic*nx*ny+index];
				data_c[index] = data_c[index] + value;

			}

			for (ic=0; ic<nc; ic++) {

				ic1 = ic;
				ic2 = (ic+1)%nc;

				value1 = offset[ic1] + data_2D[ic1*nx*ny+index];
				value2 = offset[ic2] + data_2D[ic2*nx*ny+index];

				value12 = value2 - value1;
				value12_m = (value12>=0) ? 0.0 : value12;
				value12_p = (value12>=0) ? value12 : 0.0;

				data_c[index] = data_c[index] - cross_val*(value12_p-value12_m);

			}

			index++;

		}

	}

	return;

}

/*
 * Combine nc MRC stacks located in a given directory with different possible weights.
 * The output is stored in directory/fileout.
 */
int cross_validate_stacks(int nc, char* directory, char** filenames, float* offset, float* deviation, int mode, char* fileout) {

	int error = NO_ERROR;

	/* Open the input MRC file	*/

	MrcHeader input_file_header[MAX_NUMBER_OF_FILES], output_file_header;
	FILE *input_file[MAX_NUMBER_OF_FILES], *output_file;
	Islice* slice = NULL;

	char imodfile_in[FILENAME_LEN];
	char imodfile_out[FILENAME_LEN];

	int ic;

	for (ic=0; ic<nc; ic++) {

		sprintf(&imodfile_in[0], "%s/%s", directory, filenames[ic]);

		input_file[ic] = fopen(&imodfile_in[0],"rbe");

		if (input_file[ic]==NULL) {
			txbr_error(stderr, "ERROR: cross_validate_stacks - open file %s.\n", &imodfile_in[0]);
			return FILE_ERROR;
		} else {
			printf("File %s was successfully opened.\n", &imodfile_in[0]);
		}

		mrc_head_read(input_file[ic], &input_file_header[ic]);

	}

	sprintf(&imodfile_out[0], "%s/%s", directory, fileout);

	output_file = fopen(&imodfile_out[0],"w");

	if (output_file==NULL) {
		txbr_error(stderr, "ERROR: cross_validate_stacks - open file %s.\n", &imodfile_out[0]);
		return FILE_ERROR;
	} else {
		printf("File %s was successfully opened.\n", &imodfile_out[0]);
	}

	int nx = input_file_header[0].nx;
	int ny = input_file_header[0].ny;
	int nz = input_file_header[0].nz;

	float amin = input_file_header[0].amin;
	float amax = input_file_header[0].amax;
	float amean = input_file_header[0].amean;

	for (ic=0; ic<nc; ic++) {
		amin = amin + input_file_header[ic].amin;
		amax = amax + input_file_header[ic].amax;
		amean = amean + input_file_header[ic].amean;
	}

	printf("nx=%i,ny=%i,nz=%i\n",nx,ny,nz);
	printf("min=%f,max=%f,mean=%f\n",amin,amax,amean);

	mrc_head_new(&output_file_header,nx,ny,nz,SLICE_MODE_FLOAT);

	float min = -deviation[0];
	float max = deviation[0];

	output_file_header.amin = min;
	output_file_header.amax = max;
	output_file_header.amean = 0.0;

	mrc_head_write(output_file, &output_file_header);

	int i=0, j=0, index = 0;

	float *data_in = malloc(nx*ny*nc*sizeof(float));
	float *data = malloc(nx*ny*sizeof(float));

	if (!data) {
		txbr_error(stderr, "ERROR: cross_validate_stacks - getting memory.\n");
		return MEMORY_ERROR;
	}

	int itilt;

	for (itilt=0; itilt<input_file_header[1].nz; itilt++) {

		for (ic=0; ic<nc; ic++) {

			slice = sliceReadMRC(&input_file_header[ic],itilt,'z');

			memcpy(&data_in[ic*nx*ny], &slice->data.f[0], nx*ny*sizeof(float));

			if (slice) sliceFree(slice);

		}

		cross_validate_slices(nx, ny, nc, data_in, offset, mode, data);

		mrc_write_slice(data, output_file, &output_file_header, itilt, 'z');

	}

	/* Free the resources	*/

	for (ic=0; ic<nc; ic++) {
		if (fclose(input_file[ic])) {
		txbr_error(stderr, "ERROR: cross_validate_stacks - closing file %s.\n", input_file[ic]);
		return FILE_ERROR;
		}
	}

	if (fclose(output_file)) {
		txbr_error(stderr, "ERROR: cross_validate_stacks - closing file %s.\n", output_file);
		return FILE_ERROR;
	}

	free(data_in);
	free(data);

	return error;

}

/*
 * Do a cross-validation
 */
int cross_validate(int nc, char* directory, char** filenames, int mode, char* fileout) {

	int error = NO_ERROR;

	float *offset = malloc(nc*sizeof(float));

	if (!offset) {
		txbr_error(stderr, "ERROR: crossvalidate - getting memory.\n");
		return MEMORY_ERROR;
	}

	float *deviation = malloc(nc*sizeof(float));

	if (!deviation) {
		txbr_error(stderr, "ERROR: crossvalidate - getting memory.\n");
		return MEMORY_ERROR;
	}

	int i, n = 3;
	float min, max, moments[3];

	for (i=0; i<nc; i++) {

		error = calculate_moments(directory, &filenames[i][0], &min, &max, &moments[0], n);

		offset[i] = -moments[1];
		deviation[i] = sqrt(moments[2]);
		deviation[i] = max-min;

	}

	error = cross_validate_stacks(nc, directory, filenames, offset, deviation, mode, fileout);

	free(deviation);
	free(offset);

	return error;

}


int main2(int argc, char *argv[]) {

	int error, i, nc;
	int mode = NO_CROSS_VALIDATION;

	char *directory, *fileout;
	char **filenames = malloc(MAX_NUMBER_OF_FILES*FILENAME_LEN*sizeof(char));

	if (!filenames) {
		txbr_error(stderr, "ERROR: crossval - getting memory.\n");
		return MEMORY_ERROR;
	}

	for (i=1; i<argc; i++){

		if (argv[i][0]=='-') {

			switch (argv[i][1]) {

				case 'd' :
					directory = argv[++i];
					break;

				case 'i' :

					nc = 0;
					while (argv[i+1][0]!='-') {
						filenames[nc] = argv[++i];
						nc++;
					}
					break;

				case 'o' :

					fileout = argv[++i];
					break;

				case 'r' :

					if (argv[i][2]=='p') mode = REMOVE_POSITIVE_ARTEFACTS;
					if (argv[i][2]=='m') mode = REMOVE_NEGATIVE_ARTEFACTS;
					break;


			}

		}

	}

	printf("directory: %s\n", &directory[0]);
	printf("fileout: %s\n", &fileout[0]);
	printf("Number of Files to Cross-Validate: %i\n", nc);
	for (i=0; i<nc; i++) {
		printf("filename: %s\n", &filenames[i][0]);
	}
	printf("Mode: %i\n", mode);

	error = cross_validate(nc, directory, filenames, mode, fileout);

	printf("Done\n");

}
