#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "txbrutil.h"
#include "mrcslice.h"

#define PROJECTION 1
#define BACKPROJECTION 2

int mrcReadFloatSlice(b3dFloat *buf, MrcHeader *hdata, int slice);

void print_axis_path_projection_for_tilt(PathProjection *path, int indexOfTilt, int axis);

static int storeError = 0;
static char errorMessage[MAX_TXBR_ERROR_STRING] = "";

/*
 * Set to not print messages
 */
void txbr_set_store_error(int value) {

    storeError = value;

}

/*
 * Store an error, possibly print it
 */
void txbr_error(FILE *out, char *format, ...) {

  va_list args;
  va_start(args, format);

  vsprintf(errorMessage, format, args);

  if (out && !storeError) fprintf(out, errorMessage);

}

/*
 * Return the error string
 */
char* txbr_get_error() {

  return &errorMessage[0];

}

/*
 * Initialize an array of short int to 0.
 */
void init_short_int_array(short int *values, int n) {

	int i;
	for (i=0; i<n; i++) values[i] = 0;

}

/*
 * Initialize an array of int to 0.
 */
void init_int_array(int *values, int n) {

	int i;
	for (i=0; i<n; i++) values[i] = 0;

}

/*
 * Initialize an array of float to 0.
 */
void init_float_array(float *values, int n) {

	int i;
	for (i=0; i<n; i++) values[i] = 0.0;

}

/*
 * Initialize an array of double to 0.
 */
void init_double_array(double *values, int n) {

	int i;
	for (i=0; i<n; i++) values[i] = 0.0;

}

/*
 * Print a 1D array of int.
 */
void print_int_array(int* values, int n) {

	int i=0;

	printf("[");
	for (i=0; i<n; i++) printf("%i%s",values[i],i!=n-1?",":"");
	printf("]\n");

	return;

}

/*
 * Print a 1D array of float.
 */
void print_float_array(float* values, int n) {

	int i=0;

	printf("[");
	for (i=0; i<n; i++) printf("%.4f%s",values[i],i!=n-1?",":"");
	printf("]\n");

	return;

}

/*
 * Print a 1D array of double.
 */
void print_double_array(double* values, int n) {

	int i=0;

	printf("[");
	for (i=0; i<n; i++) printf("%.3e%s",values[i],i!=n-1?",":"");
	printf("]\n");

	return;

}

/*
 * For a given order, this function returns the total number of terms
 * for the polynomial description of the electron path projection.
 */
int number_of_polynomial_terms(int order) {

	return order<0 ? 0 : (order+1)*(((order+1)+3)*(order+1)+2.0)/6.0;

}


/*
 * Initialize the power order mapping for the 3 axis X, Y and Z.
 */
void init_power_order(int order, int* order_in_X, int* order_in_Y, int* order_in_Z) {

	int i,j,k,index = 0;
	int n = number_of_polynomial_terms(order);

	/*	First set all elements to zero	*/
	for (i=0; i<n+1; i++) {
		order_in_X[i] = 0;
		order_in_Y[i] = 0;
		order_in_Z[i] = 0;
	}

	/*	Set the non-zero coefficients	*/
	for (k=0; k<order+1; k++) {
		for (i=k; i>=0; i--) {
			for (j=k-i; j>=0; j--) {
				order_in_X[index] = i;
				order_in_Y[index] = j;
				order_in_Z[index] = k-i-j;
				index++;
			}
		}
	}

}


/*
 * Print the power order for each terms in the polynom for X, Y and Z
 */
void print_power_order(int order, int* order_in_Axis) {

	int i,k,index,n;

	index=0;

	for (k=0; k<order+1; k++) {
		n = (k==0) ? 1 : (k+1)*(k+2)/2;
		for (i=0; i<n; i++) {
			printf("%3i",order_in_Axis[index++]);
		}
		printf("\t");
	}
	printf("\n");

}


/*
 * Set the path correction.
 */
void set_path_projection_coefficients(PathProjection *path, int tilt, int order, int axis, double *coefficients) {

	/* for a transformation whose order is order, there is order(order+1) terms */
	/* order = 0 corresponds to the translation component */

	int i;
	int index0 = number_of_polynomial_terms(order-1);
	int n = number_of_polynomial_terms(order);

	switch (axis) {
		case 1 :
			for (i=0; i<n; i++) {
				path->coefficients_1[tilt][index0+i] = coefficients[i];
			}
			break;
		case 2 :
			for (i=0; i<n; i++) {
				path->coefficients_2[tilt][index0+i] = coefficients[i];
			}
			break;
	}


}

/*
 * Print informations about the electron paths projection coefficients.
 */
void print_axis_path_projection_for_tilt(PathProjection *path, int indexOfTilt, int axis) {

	int i1,i2,index0,n;
	double coeffs[MAX_COEFFICIENTS_NUMBER];
	char* title = "";

	switch (axis) {
		case 0 : memcpy(coeffs, &path->lambda[indexOfTilt][0], 4*sizeof(double)); title="lambda"; break;
		case 1 : memcpy(coeffs, &path->coefficients_1[indexOfTilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double)); title="     x"; break;
		case 2 : memcpy(coeffs, &path->coefficients_2[indexOfTilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double)); title="     y"; break;
	}

	printf("%s:\t",title);

	int effective_order = axis==0 ? 1 : path->order;

	for (i1=0; i1<effective_order+1; i1++) {
		n = i1==0 ? 1 : (i1+1)*(i1+2)/2;
		index0 = number_of_polynomial_terms(i1-1);
		for (i2=0; i2<n; i2++) {
			printf("%.1E ",coeffs[index0 + i2]);
		}
		printf(" ");
	}

	printf("\n");

}

/*
 * Print informations about the electron path projection.
 */
void print_path_projection(PathProjection *path) {

	int numberOfTilts = path->numberOfTilts;
	int i;

	printf("-------------------------- Projection Coefficients --------------------------------------------\n");
	printf("Order           : %i\n",path->order);
	printf("Number of Tilts : %i\n",path->numberOfTilts);

        printf("Skip            :");
        for (i=0; i<numberOfTilts; i++) {
            if (path->skipView[i]) printf(" %i",i);
        }
        printf("\n");

	for (i=0; i<numberOfTilts; i++) {
		if (i%20!=0) continue; /* We only print the info every 20 tilts	*/
		printf("%i/%i\n",i,numberOfTilts);
		print_axis_path_projection_for_tilt(path,i,0);
		print_axis_path_projection_for_tilt(path,i,1);
		print_axis_path_projection_for_tilt(path,i,2);
	}
	printf("-----------------------------------------------------------------------------------------------\n");

}

/*
 * Print informations about the reconstruction setup
 */
void print_TxBR_setup(TxBRsetup *setup) {

	printf("------------------ Backprojection Setup --------------------\n");
	printf("Directory               : %s\n",&setup->directory[0]);
	printf("Basename                : %s\n",&setup->basename[0]);
	printf("Size of a Block         : %i\n",setup->blocksize);
	printf("Number Of Blocks        : %i\n",setup->numberOfBlocks);
	printf("Coefficients of Plane 1 : %10.4E %10.4E %10.4E\n",setup->plane_coeffs1[0],setup->plane_coeffs1[1],setup->plane_coeffs1[2]);
	printf("Coefficients of Plane 2 : %f %f %f\n",setup->plane_coeffs2[0],setup->plane_coeffs2[1],setup->plane_coeffs2[2]);
	printf("X                       : %f:%f:%f\n",setup->x_0,setup->x_inc,setup->x_1);
	printf("Y                       : %f:%f:%f\n",setup->y_0,setup->y_inc,setup->y_1);
	printf("Z                       : %f:%f:%f\n",setup->z_0,setup->z_inc,setup->z_1);
	printf("Scale                   : %f:%f:%f\n",setup->sx,setup->sy,setup->sz);
	//printf("Z Start                 : %f\n",setup->z_start);
	printf("-----------------------------------------------------------\n");

}

/*
 * Open an MRC File Header;
 */
int openMRCFile(char* directory, char* basename, MrcHeader* header, FILE** fp) {

	char imodfile_in[FILENAME_LEN];

	sprintf(imodfile_in, "%s/%s.preali.rot.SL", directory, basename);
        
	printf("File to process: %s.\n",imodfile_in);

	printf("Reading the File %s...\n",imodfile_in);

	if((*fp = fopen(imodfile_in,"rbe"))==NULL) {

		txbr_error(stderr, "ERROR: openMRCFile - cannot open file %s.\n", imodfile_in);

		return FILE_ERROR;

	}

	printf("Read the MRC header...\n");

	mrc_head_read(*fp,header);

	printf("Finish Read the MRC header...\n");

	printf("\n");

	return NO_ERROR;

}

/*
 * Load a slice of MRC stack
 */
int load_projected_slice( char* directory, char* basename, int index, float* slice ) {

    char proj_src[FILENAME_LEN];
    MrcHeader *input_file_header;
    FILE *input_file;

    sprintf(proj_src, "%s/%s.st", directory, basename);

    printf("%s\n",proj_src);

    input_file_header  = (MrcHeader *)malloc(sizeof(MrcHeader));

    if (!input_file_header) {
        txbr_error(stderr, "ERROR: load_projected_slice - getting memory.\n");
        return MEMORY_ERROR;
    }

    if ((input_file = fopen(proj_src, "rbe")) == NULL) { /* OPen file */
        txbr_error(stderr, "ERROR: load_projected_slice - cannot open file %s.\n", proj_src);
        return FILE_ERROR;
    }

    mrc_head_read( proj_src, input_file_header );

    mrcReadFloatSlice( slice, input_file_header, index );

    if (fclose(input_file)) { /* Close file */
        txbr_error(stderr, "ERROR: load_projected_slice - closing file %s.\n", input_file);
        return FILE_ERROR;
    }

    return NO_ERROR;

}

int save_projection_block( char* directory, char* basename, int indexOfBlock,
                 int nx, int ny, int blocksize, float* block ) {

    double x0 = 0.0;
    double y0 = 0.0;
    double z0 = 0.0;
    double sx = 1.0;
    double sy = 1.0;
    double sz = 1.0;

    return save_block_( directory, basename, indexOfBlock, x0, y0, z0, nx, ny, blocksize, sx, sy, sz,
                        block, PROJECTION );

}

int save_backprojection_block( char* directory, char* basename, int indexOfBlock, double x0,
                               double y0, double z0, int nx, int ny, int blocksize, double sx,
                               double sy, double sz, float* block ) {

    return save_block_( directory, basename, indexOfBlock, x0, y0, z0, nx, ny, blocksize, sx, sy, sz,
                        block, BACKPROJECTION );

}


/*
 * Save a block of data!
 */
int save_block_( char* directory, char* basename, int indexOfBlock, double x0,
                 double y0, double z0, int nx, int ny, int blocksize, double sx,
                 double sy, double sz, float* block, int type ) {

    char imodfile_out[FILENAME_LEN];
    struct MRCheader hdata;
    FILE* output;
    char *mode;

    int offset = 0;
    int nn = nx*ny;
    int nz = blocksize;
    int isect;

    double min=0.0, max=0.0, mean=0.0, m2=0.0;

    if (type==PROJECTION) {
        sprintf(imodfile_out, "%s/%s.project.st", directory, basename);
        printf("Saving block #%i. %i slices in file %s.\n", indexOfBlock + 1, blocksize, imodfile_out);
    } else if (type==BACKPROJECTION) {
        sprintf(imodfile_out, "%s/%s_z_%.2f.out", directory, basename, z0);
        printf("Saving block #%i. %i slices in file %s.\n", indexOfBlock + 1, blocksize, imodfile_out);
    }

    if (indexOfBlock==0) {
        mode = "w+";
    } else {
        mode = "rw+";
    }

    if ((output = fopen(&imodfile_out[0], mode)) == NULL) {
        txbr_error(stderr, "ERROR: save_block - cannot open file %s.\n", &imodfile_out[0]);
        return FILE_ERROR;
    }

    if (indexOfBlock==0) {
        mrc_head_new(&hdata, nx, ny, blocksize, SLICE_MODE_FLOAT);
        hdata.xorg = -x0*sx+sx;
        hdata.yorg = -y0*sy+sy;
        hdata.zorg = -z0*sz;
        hdata.nxstart = (int)x0;
	hdata.nystart = (int)y0;
	hdata.nzstart = (int)z0;
    } else {
        mrc_head_read(output, &hdata);
        offset = hdata.nz;
        nz = offset + blocksize;
        hdata.nz = nz;
        hdata.mz = nz;
        mean = nn * offset * hdata.amean;
        m2 = (hdata.amax - hdata.amin) / 8.0;
        m2 = nn * offset * (m2 * m2 + hdata.amean * hdata.amean);
    }

    for (isect = 0; isect < nn * blocksize; isect++) {
        min = TXBR_MIN(min,block[isect]);
        max = TXBR_MAX(max,block[isect]);
        mean += block[isect];
        m2 += block[isect]*block[isect];
    }

    mean = mean/nn/nz;
    m2 = m2/nn/nz - mean*mean;

    hdata.amin = mean - 4.0 * sqrt(m2);
    hdata.amax = mean + 4.0 * sqrt(m2);
    hdata.amean = mean;

    int digits = 2;
    double factor = pow(10.0,digits-ceil(log10(fabs(10*sqrt(m2)))));

    hdata.amin = floor(hdata.amin*factor)/factor;
    hdata.amax = ceil(hdata.amax*factor)/factor;
    
    mrc_set_scale( &hdata, sx, sy, sz );

    mrc_head_write(output, &hdata);

    printf("min: %f, max %f, mean: %f, m2: %f\n", min, max, mean, m2);

    for (isect = 0; isect < blocksize; isect++) {
        mrc_write_slice(&block[isect * nn], output, &hdata, isect+ offset, 'z');
    }

    if (fclose(output)) {
        txbr_error(stderr, "ERROR: save_block - closing file %s.\n", &imodfile_out[0]);
        return FILE_ERROR;
    }

    return NO_ERROR;

}










/*
 * Save a block of data!
 */
int save_block(char* directory, char* basename, double z0, int indexOfBlock, int nx, int ny, int blocksize, float* block) {

	char imodfile_out[FILENAME_LEN];

	sprintf(imodfile_out, "%s/%s_z_%.2f-%i.out", directory, basename, z0, indexOfBlock);

	printf("Saving block #%i. %i slices in file %s.\n",indexOfBlock+1,blocksize,imodfile_out);
//	printf("nx=%i  ny=%i\n",nx,ny);

	FILE* output;

	if((output=fopen(&imodfile_out[0],"w"))==NULL) {

		txbr_error(stderr, "ERROR: save_block - cannot open file %s.\n",&imodfile_out[0]);

		return FILE_ERROR;

	}

	struct MRCheader hdata;

	mrc_head_new(&hdata, nx, ny, blocksize, SLICE_MODE_FLOAT);

	int nn = nx*ny;
	int isect;

	double min=0, max = 0, mean = 0, m2 = 0;

	for (isect=0; isect<nn*blocksize; isect++) {
		min = TXBR_MIN(min,block[isect]);
		max = TXBR_MAX(min,block[isect]);
		mean += block[isect];
		m2 += block[isect]*block[isect];
	}

	mean = mean/nn/blocksize;
	m2 = m2/nn/blocksize-mean*mean;

	hdata.amin = mean - 4.0*sqrt(m2);
	hdata.amax = mean + 4.0*sqrt(m2);
	hdata.amean = mean;

        hdata.zorg = z0 + indexOfBlock*blocksize;


	mrc_head_write(output, &hdata);

	printf("min: %f, max %f, mean: %f, m2: %f\n",min,max,mean,m2);

	for (isect=0; isect<blocksize; isect++) {

		mrc_write_slice(&block[isect*nn], output, &hdata, isect, 'z');

	}

	if(fclose(output)) {

		txbr_error(stderr, "ERROR: save_block - closing file %s.\n",&imodfile_out[0]);

		return FILE_ERROR;

	}

	return NO_ERROR;

}

///*
// * Merge blocks together!
// */
//int merge_blocks(char* directory, char* basename, double z0, int numberOfBlocks) {
//
//	char imodfile_out[FILENAME_LEN];
//
//	int nx = 0;
//	int ny = 0;
//
//	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));
//
//	if (!header) {
//
//		txbr_error(stderr, "ERROR: merge_blocks - getting memory.\n");
//
//		return MEMORY_ERROR;
//
//	}
//
//	FILE* fp;
//	int iblock;
//
//	double min = 0.0;
//	double max = 0.0;
//	double mean = 0.0;
//
//	char files[MAX_NUMBER_OF_FILES_TO_MERGE][FILENAME_LEN];
//
//	int number_of_files_to_merge = 0;
//
//	int totalNumberOfSlices = 0;
//
//	/*	Determine first the right parameters (nx,ny,numberOfSlices)	by reading the MRC header.*/
//
//	for (iblock=0; iblock<numberOfBlocks; iblock++) {
//
//		sprintf(&files[number_of_files_to_merge][0], "%s/%s_z_%.2f-%i.out", directory, basename, z0, iblock);
//
//		if((fp=fopen(files[number_of_files_to_merge],"r"))==NULL) {
//
//			txbr_error(stderr, "ERROR: merge_blocks - cannot open file %s.\n", &files[number_of_files_to_merge][0]);
//
//			return FILE_ERROR;
//
//		}
//
//		printf("Cheak header: %s\n",files[number_of_files_to_merge]);
//
//		mrc_head_read(fp,header);
//
//		if (nx==0 && ny==0) {
//			nx = header->nx;
//			ny = header->ny;
//		}
//
//		if (header->nx==nx && header->ny==ny) {
//
//			number_of_files_to_merge++;
//
//			min = TXBR_MIN(min,header->amin);
//			max = TXBR_MAX(max,header->amax);
//			mean += header->amean*header->nz;
//
//			totalNumberOfSlices += header->nz;
//
//		}
//
//		if(fclose(fp)) {
//
//			txbr_error(stderr, "ERROR: merge_blocks - closing file %s.\n", &files[number_of_files_to_merge][0]);
//
//			return FILE_ERROR;
//
//		}
//
//	}
//
//	/*	Right parameters are now known	*/
//
//	sprintf(imodfile_out, "%s/%s_z_%.2f.out", directory, basename, z0);
//
//	FILE* output;
//
//	if((output=fopen(&imodfile_out[0],"w"))==NULL) {
//
//		txbr_error(stderr, "ERROR: merge_blocks - cannot open file %s.\n", &imodfile_out[0]);
//
//		return FILE_ERROR;
//
//	}
//
//	/* MRC header for the main output file	*/
//
//	struct MRCheader hdata;
//
//	mrc_head_new(&hdata, nx, ny, totalNumberOfSlices, SLICE_MODE_FLOAT);
//
//	hdata.amin = min;
//	hdata.amax = max;
//	hdata.amean = mean/number_of_files_to_merge/totalNumberOfSlices;
//
//	mrc_head_write(output, &hdata);
//
//	/* MRC header for current file	*/
//
//	printf("Merging #%i blocks together\n",number_of_files_to_merge);
//
//	Islice* slice = NULL;
//	int numberOfSlices = 0;	/*	normally should be blocksize	*/
//
//	int islice=0, index=0;
//
//	for (iblock=0; iblock<number_of_files_to_merge; iblock++) {
//
//		printf("Processing block #%i... File %s.\n",iblock,files[iblock]);
//
//		if((fp=fopen(files[iblock],"r"))==NULL) {
//
//			txbr_error(stderr, "ERROR: merge_blocks - cannot open file %s.\n", &files[iblock][0]);
//
//			return FILE_ERROR;
//
//		}
//
//		mrc_head_read(fp,header);
//
//		numberOfSlices = header->nz;
//
//		for (islice=0; islice<numberOfSlices; islice++) {
//
//			printf("Reading Slice #%i of block %i\n",islice,iblock);
//
//			slice = sliceReadMRC(header,islice,'z');
//
//			mrc_write_slice(slice->data.f, output, &hdata, index++, 'z');
//
//			if (slice) sliceFree(slice);
//
//		}
//
//		if(fclose(fp)) {
//
//			txbr_error(stderr, "ERROR: merge_blocks - closing file %s.\n", &files[iblock][0]);
//
//			return FILE_ERROR;
//
//		}
//
//	}
//
//	if(fclose(output)) {
//
//		txbr_error(stderr, "ERROR: merge_blocks - closing file %s.\n", &imodfile_out[0]);
//
//		return FILE_ERROR;
//
//	}
//
//	/* Free stuff	*/
//
//	free(header);
//
//	return NO_ERROR;
//
//}

/*
 * This routine takes the fist n_1 slice of a stack (specified by directory/filename), multiply them by alpha_1.
 * Then takes the n_2 next ones and multiply them by alpha_2.
 */
int recompose_stack(char* directory, char* filename, int n_1, int n_2, double alpha_1, double alpha_2) {

	char file_in[FILENAME_LEN];
	sprintf(file_in, "%s/%s", directory, filename);

	char imodfile_out[FILENAME_LEN];
	sprintf(imodfile_out, "%s/%s-recomp_%i_%i_%f_%f", directory, filename, n_1, n_2, alpha_1, alpha_2);

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) {

		txbr_error(stderr, "ERROR: recompose_stack - getting memory.\n");

		return MEMORY_ERROR;

	}

	FILE* fp;
	int iblock;

	if((fp=fopen(&file_in[0],"r"))==NULL) {

		txbr_error(stderr, "ERROR: recompose_stack - cannot open file %s.\n", &file_in[0]);

		return FILE_ERROR;

	}

	printf("Cheak header: %s\n",&file_in[0]);

	mrc_head_read(fp,header);

	int nx = header->nx;
	int ny = header->ny;
	int numberOfSlices = header->nz;

	if (numberOfSlices!=n_1+n_2) {
		printf("The number of slices %i is different from %i + %i\n", numberOfSlices, n_1, n_2);
	}

	float* block = (float*)malloc(nx*ny*sizeof(float));

	/* Handle for the output file	*/

	FILE* output;

	if((output=fopen(&imodfile_out[0],"w"))==NULL) {

		txbr_error(stderr, "ERROR: recompose_stack - cannot open file %s.\n", &imodfile_out[0]);

		return FILE_ERROR;

	}

	/* MRC header for the main output file	*/

	struct MRCheader hdata;

	mrc_head_new(&hdata, nx, ny, numberOfSlices, SLICE_MODE_FLOAT);

	double min=0, max = 0, mean = 0, m2 = 0;

	Islice* slice = NULL;
	int islice,index;
	double alpha;

	for (islice=0; islice<numberOfSlices; islice++) {

		printf("%s %i\0",islice==0 ? "Reading Slice #" : ",",islice);
		fflush(stdout);

		alpha = (islice<n_1) ? alpha_1 : (islice<n_1+n_2) ? alpha_2 : 1;

		slice = sliceReadMRC(header,islice,'z');

		index = 0;

		for (index=0; index<nx*ny; index++) {

			block[index] = alpha*(slice->data.f[index]);

			min = TXBR_MIN(min,block[index]);
			max = TXBR_MAX(min,block[index]);
			mean += block[index];
			m2 += block[index]*block[index];

		}

		mrc_write_slice(&block[0], output, &hdata, islice, 'z');

		if (slice) sliceFree(slice);

	}

	printf("\n");

	mean = mean/nx/ny/numberOfSlices;
	m2 = m2/nx/ny/numberOfSlices-mean*mean;

	hdata.amin = mean - 4.0*sqrt(m2);
	hdata.amax = mean + 4.0*sqrt(m2);
	hdata.amean = mean;

	mrc_head_write(output, &hdata);

	if(fclose(fp)) {

		txbr_error(stderr, "ERROR: recompose_stack - closing file %s.\n", &file_in[0]);

		return FILE_ERROR;

	}

	if(fclose(output)) {

		txbr_error(stderr, "ERROR: recompose_stack - closing file %s.\n", &imodfile_out[0]);

		return FILE_ERROR;

	}

	/* Free stuff	*/

	free(block);

	free(header);

	return NO_ERROR;

}


