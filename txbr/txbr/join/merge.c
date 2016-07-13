#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mrcfiles.h"
#include "mrcslice.h"
#include "txbrutil.h"

#define MAX_NUMBER_OF_FILES 50

int combine_slices(int nx, int ny, int nc, float* data_2D, float* weight, float* data_c);
int combine_stacks(int nc, char** filenames, float* weight, char* fileout);

/*
 * Combine 2D arrays of data in one array with a given weight
 */
int combine_slices(int nx, int ny, int nc, float* data_2D, float* weight, float* data_c) {

	int ix, iy, ic, index;

	for (index=0; index<nx*ny; index++) {
		data_c[index] = 0;
	}

	for (ic=0; ic<nc; ic++) {
		index = 0;
		for (iy=0; iy<ny; iy++) {
			for (ix=0; ix<nx; ix++) {
				data_c[index] = data_c[index] + weight[ic]*data_2D[ic*nx*ny+index];
				index++;
			}
		}
	}

	return;

}


/*
 * Combine nc MRC stacks, whose locations in the file system is specified by
 * filenames, with different possible weights.
 * The output is stored in fileout.
 * Variables filenames: a pointer of pointer specifying location of the input 
 * files.
 * Variables filenames: a pointer specifying the first element of an array
 * of nc elements to weight the MRC file sommation.
 */
int combine_stacks(int nc, char** filenames, float* weight, char* fileout) {

	int error = NO_ERROR;

	/* Open the input MRC file	*/

	MrcHeader input_file_header[MAX_NUMBER_OF_FILES], output_file_header;
	FILE *input_file[MAX_NUMBER_OF_FILES], *output_file;
	Islice* slice = NULL;

	char imodfile_in[FILENAME_LEN];
	char imodfile_out[FILENAME_LEN];

	int ic;

	for (ic=0; ic<nc; ic++) {

		printf("%s to combine\n",filenames[ic]);

		sprintf(&imodfile_in[0], "%s", filenames[ic]);

		input_file[ic] = fopen(&imodfile_in[0],"rbe");

		if (input_file[ic]==NULL) {
			txbr_error(stderr, "ERROR: combine_stacks - open file %s.\n", &imodfile_in[0]);
			return FILE_ERROR;
		} else {
			printf("File %s was successfully opened.\n", &imodfile_in[0]);
		}

		mrc_head_read(input_file[ic], &input_file_header[ic]);

	}

	sprintf(&imodfile_out[0], "%s", fileout);

	output_file = fopen(&imodfile_out[0],"w");

	if (output_file==NULL) {
		txbr_error(stderr, "ERROR: rotate_stack - open file %s.\n", &imodfile_out[0]);
		return FILE_ERROR;
	} else {
		printf("File %s was successfully opened.\n", &imodfile_out[0]);
	}

	int nx = input_file_header[0].nx;
	int ny = input_file_header[0].ny;
	int nz = input_file_header[0].nz;

	float amin = weight[0]*input_file_header[0].amin;
	float amax = weight[0]*input_file_header[0].amax;
	float amean = weight[0]*input_file_header[0].amean;

	for (ic=0; ic<nc; ic++) {
		amin = amin + input_file_header[ic].amin;
		amax = amax + input_file_header[ic].amax;
		amean = amean + weight[ic]*input_file_header[ic].amean;
	}

	printf("nx=%i,ny=%i,nz=%i\n",nx,ny,nz);
	printf("min=%f,max=%f,mean=%f\n",amin,amax,amean);

	mrc_head_new(&output_file_header,nx,ny,nz,SLICE_MODE_FLOAT);

	output_file_header.amin = amin;
	output_file_header.amax = amax;
	output_file_header.amean = amean;

	mrc_head_write(output_file, &output_file_header);

	int i=0, j=0, index = 0;

	float *data_in = malloc(nx*ny*nc*sizeof(float));
	float *data = malloc(nx*ny*sizeof(float));

	if (!data) {
		txbr_error(stderr, "ERROR: rotate_stack - getting memory.\n");
		return MEMORY_ERROR;
	}

	int itilt;

	for (itilt=0; itilt<input_file_header[1].nz; itilt++) {

		for (ic=0; ic<nc; ic++) {

			slice = sliceReadMRC(&input_file_header[ic],itilt,'z');

			memcpy(&data_in[ic*nx*ny], &slice->data.f[0], nx*ny*sizeof(float));

			if (slice) sliceFree(slice);

		}

		combine_slices(nx, ny, nc, data_in, weight, data);

		mrc_write_slice(data, output_file, &output_file_header, itilt, 'z');

	}

	/* Free the resources	*/

	for (ic=0; ic<nc; ic++) {
		if (fclose(input_file[ic])) {
		txbr_error(stderr, "ERROR: do_block_reconstruction - closing file %s.\n", input_file[ic]);
		return FILE_ERROR;
		}
	}

	if (fclose(output_file)) {
		txbr_error(stderr, "ERROR: do_block_reconstruction - closing file %s.\n", output_file);
		return FILE_ERROR;
	}

	free(data_in);
	free(data);

	return error;

}


int main1(int argc, char *argv[]) {

	int i, nc;

	char *directory, *fileout;

	char **filenames = malloc(MAX_NUMBER_OF_FILES*FILENAME_LEN*sizeof(char));

	if (!filenames) {
		txbr_error(stderr, "ERROR: icombine - getting memory.\n");
		return MEMORY_ERROR;
	}

	float *weight = malloc(MAX_NUMBER_OF_FILES*sizeof(float));

	if (!weight) {
		txbr_error(stderr, "ERROR: icombine - getting memory.\n");
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
						weight[nc] = (float)atof(argv[++i]);
						nc++;
					}
					break;

				case 'o' :

					fileout = argv[++i];
					break;


			}

		}

	}
	
	char *imodfile_in[FILENAME_LEN];
	char imodfile_out[FILENAME_LEN];

	printf("directory: %s\n", &directory[0]);
	printf("fileout: %s\n", &fileout[0]);
	printf("Number of Files to combine: %i\n", nc);
	for (i=0; i<nc; i++) {
		sprintf(&imodfile_in[0], "%s/%s", directory, filenames[i]);
		printf("filename: %s/%s\n", directory, &filenames[i][0]);
		printf("weight: %f\n", weight[i]);
	}
	sprintf(&imodfile_out[0], "%s/%s", directory, fileout);

	combine_stacks(nc, &imodfile_in, weight, &imodfile_out);

	printf("Done\n");

}
