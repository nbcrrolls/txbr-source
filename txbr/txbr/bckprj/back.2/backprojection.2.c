#include <string.h>
#include <math.h>
#include <time.h>

#include "bckprj.h"

MrcHeader* input_file_header;
FILE* input_file;

int x_inc = 1, y_inc = 1;

int x0_dest_f, y0_dest_f;
int nx_dest_f, ny_dest_f;

float *block_dest;

/*	 Define some eventual test parameters	*/

#if TEST>=1

int ix_1, iy_1, iz_1, itilt_1;
int ix_2, iy_2, iz_2, itilt_2;
double ll=0, ff1 = 0, ff2 = 0;
double max_absolute_error1 = 0, max_absolute_error2 = 0;

#endif

#if SWAP_x_y==1

double tmp;

#endif

/*	 Functions	*/

void calculate_box( int x_start, int y_start, int z_start, int x_stop, int y_stop, int z_stop,
					int nx, int ny, int *x0_src, int *y0_src, int *nx_src, int *ny_src);

void initializeGPU( int n__, int number_of_terms_, int *order_in_X_, int *order_in_Y_, int *order_in_Z_,
					double *lambda__, double *a1__, double *a2__, double *bottom_plane__, double *rot__, int x_inc);

int calculate_segment_size(int order, int x0, int y0, int z0, int x_inc, int x1, double *coefficients, double eps);

int evaluateBlock( float *block, int x0, int y0, int z0, int x_inc, int y_inc, int z_inc, int x1, int y1, int z1,
				   float* slice_data_f, int x0_src, int y0_src, int nx_src, int ny_src, int nsegment, int segment_size);

int mrcReadSectionFloat(MrcHeader *hdata, IloadInfo *li, b3dFloat *buf, int z);

void sliceWriteBlock( float *block_dest, int x0_dest, int y0_dest, int nx_dest, int ny_dest, int nz_dest,
					  float *block_src, int x0_src, int y0_src, int nx_src, int ny_src, int nz_src);

/*
 * Initialize the rotation coefficients from affine coefficents
 * of the bottom/median plane of the sample.
 */
void initialize_rotation_coefficients(double* plane_coefficients, double *_rot_) {

	double d1_ = plane_coefficients[1]*plane_coefficients[1] + plane_coefficients[2]*plane_coefficients[2];
	double d1_sqr_root_ = sqrt(d1_);

	double d2_ = 1.0 + d1_;
	double d2_sqr_root_ = sqrt(d2_);

	/* theta defines the angle of the rotation	*/
	/* phi + pi/2 is the angle between the rotation axis and the x axis in the xy plane	*/
	double cos_phi, sin_phi, cos_theta, sin_theta;

	cos_phi = plane_coefficients[1]/d1_sqr_root_;
	sin_phi = plane_coefficients[2]/d1_sqr_root_;

	if (d1_==0) {
		cos_phi = 1.0;
		sin_phi = 0.0;
	}

	cos_theta = 1.0/d2_sqr_root_;
	sin_theta = d1_sqr_root_/d2_sqr_root_;

	_rot_[0] = cos_phi*cos_phi*(cos_theta-1.0) + 1.0;	// _rot_[0][0]
	_rot_[1] = cos_phi*sin_phi*(cos_theta-1.0);	// _rot_[0][1]
	_rot_[2] = -cos_phi*sin_theta;	// _rot_[0][2]

	_rot_[3] = cos_phi*sin_phi*(cos_theta-1.0);	// _rot_[1][0]
	_rot_[4] = sin_phi*sin_phi*(cos_theta-1.0) + 1.0;	// _rot_[1][1]
	_rot_[5] = -sin_phi*sin_theta;	// _rot_[1][2]

	_rot_[6] = cos_phi*sin_theta;	// _rot_[2][0]
	_rot_[7] = sin_phi*sin_theta;	// _rot_[2][1]
	_rot_[8] = cos_theta ;	// _rot_[2][2]

}



/*****************************************************************************
* Function:     do_reconstruction
*
* Argument:     directory		The directory for this project
* 				basename		The basename for this reconstruction project
* 				work_directory	The work directory
* 				vol				The coarse volume that will be reconstructed
* 				proj_map		The projection Map
*
* Returns:      void
*
* Description:  This routine reconstructs a tomographic volume from a set of
*  				projection maps.
*
*****************************************************************************/

int do_reconstruction( char* directory, char* basename, char* work_directory,
					   Volume *vol, PathProjection* proj_map) {

	TxBRsetup *setup = vol->setup;

	int n_;
	int number_of_terms;
	int order_in_X[MAX_COEFFICIENTS_NUMBER];
	int order_in_Y[MAX_COEFFICIENTS_NUMBER];
	int order_in_Z[MAX_COEFFICIENTS_NUMBER];

	double lambda_[4];
	double a1_[MAX_COEFFICIENTS_NUMBER];
	double a2_[MAX_COEFFICIENTS_NUMBER];
	double bottom_plane_[3];
	double rot_[3][3];

	// Boundaries of the Volume to reconstruct

	int x_start = (int)vol->x_start;
	int y_start = (int)vol->y_start;
	int z_start = (int)vol->z_start;
	int x_stop = (int)vol->x_stop;
	int y_stop = (int)vol->y_stop;
	int z_stop = (int)vol->z_stop;

	// others

	int x_inc = (int)(setup->x_inc);
	int y_inc = (int)(setup->y_inc);
	int z_inc = (int)(setup->z_inc);

	int blocksize_x = (x_stop-x_start+1)/x_inc;
	int blocksize_y = (y_stop-y_start+1)/y_inc;
	int blocksize_z = (z_stop-z_start+1)/z_inc;

	printf("(blocksize_x,blocksize_y,blocksize_z)->(%i,%i,%i)\n",blocksize_x,blocksize_y,blocksize_z);

	if (blocksize_z>BLOCK_SIZE_MAX) {
		printf("Block size %i is higher than recommanded %i.\n", blocksize_z, BLOCK_SIZE_MAX);
		return (PARAMETER_ERROR);
	}

	int error = NO_ERROR;

	clock_t t_0,t_1;
	float ratio = 1./CLOCKS_PER_SEC;

	t_0 = clock();

	memcpy(bottom_plane_, setup->plane_coeffs2, sizeof(setup->plane_coeffs2));

	initialize_rotation_coefficients(bottom_plane_,rot_);

	printf("rot[0][0]=%4.5f\trot[0][1]=%4.5f\trot[0][2]=%4.5f\t\n",rot_[0][0],rot_[0][1],rot_[0][2]);
	printf("rot[0][0]=%4.5f\trot[0][1]=%4.5f\trot[0][2]=%4.5f\t\n",rot_[1][0],rot_[1][1],rot_[1][2]);
	printf("rot[0][0]=%4.5f\trot[0][1]=%4.5f\trot[0][2]=%4.5f\t\n",rot_[2][0],rot_[2][1],rot_[2][2]);

	n_ = proj_map->order;
	number_of_terms = number_of_polynomial_terms(n_);

	init_power_order(n_, &order_in_X[0], &order_in_Y[0], &order_in_Z[0]);

	int nx = input_file_header->nx;
	int ny = input_file_header->ny;
	int numtlts = input_file_header->nz;

	printf("Reconstruction for %s from z=%f to z=%f\n", basename, vol->z_start, vol->z_stop);
	printf("(nx,ny,ntilts)->(%i,%i,%i)\n",nx,ny,numtlts);

	if (nx>X_MAX || ny>Y_MAX) {
		printf("Image size (%i,%i) is higher than the recommanded (%i,%i).\n", nx, ny, X_MAX,Y_MAX);
	}

	// Image Source

	int x0_src = x_start;
	int y0_src = y_start;

	int nx_src = x_stop - x_start + 1;	// should be independent of the increment
	int ny_src = y_stop - y_start + 1;

	// Load Info

	IloadInfo li;

    rewind(input_file_header->fp);
    fseek(input_file_header->fp, input_file_header->headerSize, SEEK_SET);

	// Destination block

	int nx_eff = blocksize_x;
	int ny_eff = blocksize_y;

	printf("(blocksize_x,blocksize_y,blocksize_z)->(%i,%i,%i)\n",blocksize_x,blocksize_y,blocksize_z);
	printf("(nx_eff,ny_eff,ntilts)->(%i,%i,%i)\n",nx_eff,ny_eff,numtlts);

	float* block = (float*)malloc(nx_eff*ny_eff*blocksize_z*sizeof(float));

	if (!block) {
		txbr_error(stderr, "ERROR: backprojection - getting memory for block.\n");
		return MEMORY_ERROR;
	}

	init_float_array(&block[0],nx_eff*ny_eff*blocksize_z);	/*	Initialization of the Output Image Buffer	*/

	// Misc

	int itilt;

	memcpy(lambda_, &proj_map->lambda[0], 4*sizeof(double));
	memcpy(a1_, &proj_map->coefficients_1[0], MAX_COEFFICIENTS_NUMBER*sizeof(double));
	memcpy(a2_, &proj_map->coefficients_2[0], MAX_COEFFICIENTS_NUMBER*sizeof(double));

	initializeGPU(n_, number_of_terms, &order_in_X[0], &order_in_Y[0], &order_in_Z[0], lambda_, a1_, a2_, bottom_plane_, rot_, setup->x_inc);

	int seg_size = nx_eff;

	seg_size = calculate_segment_size( proj_map->order, setup->x_0, setup->y_0, setup->z_0,
									   setup->x_inc, setup->x_1, &proj_map->coefficients_1[0], EPSILON);

	int nseg = floor(nx_eff/seg_size) + 1;

	#if TEST>=1

	printf("Reconstruction routine:\n");
	printf("	(x0,y0,z0)=(%i,%i,%i)\n",x_start,y_start,z_start);
	printf("	(x_inc,y_inc,z_inc)=(%i,%i,%i)\n",x_inc,y_inc,z_inc);
	printf("	(x1,y1,z1)=(%i,%i,%i)\n",x_stop,y_stop,z_stop);
	printf("	(nx,ny)=(%i,%i)\n",nx,ny);
	printf("	(nx_eff,ny_eff)=(%i,%i)\n",nx_eff,ny_eff);
	printf("	(nseg,seg_size)=(%i,%i)\n",nseg,seg_size);

	#endif

	int data_f_size = 2*nx_src*ny_src;
	float* data_f = (float*)malloc(data_f_size*sizeof(float));

	if (!data_f) {
		txbr_error(stderr, "ERROR: backprojection - getting memory for data_f.\n");
		return MEMORY_ERROR;
	}

   	for (itilt=0; itilt<numtlts; itilt++) {

		if (proj_map->skipView[itilt]==1) {
			printf("Exposure %i has been skipped!",itilt);
			continue;
		}

		//	 Copy the bundle adjustment parameters for a given tilt

		memcpy(lambda_, &proj_map->lambda[itilt][0], 4*sizeof(double));
		memcpy(a1_, &proj_map->coefficients_1[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));
		memcpy(a2_, &proj_map->coefficients_2[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));

		 //	 Read input image

		printf("%s %i\0",itilt==0 ? "Reading Image #" : ",",itilt);
		fflush(stdout);

		calculate_box(x_start, y_start, z_start, x_stop, y_stop, z_stop, nx, ny, &x0_src, &y0_src, &nx_src, &ny_src);

		if (data_f_size<nx_src*ny_src) {

			free(data_f);

			data_f_size = 2*nx_src*ny_src;
			data_f = (float*)malloc(data_f_size*sizeof(float));

			if (!data_f) {
				txbr_error(stderr, "ERROR: backprojection - getting memory for data_f.\n");
				return MEMORY_ERROR;
			}

		}

	    li.xmin = x0_src - OFFSET;
	    li.xmax = x0_src + nx_src - 1 - OFFSET;
	    li.ymin = y0_src - OFFSET;
	    li.ymax = y0_src + ny_src - 1 - OFFSET;

	    mrcReadSectionFloat(input_file_header, &li, data_f, itilt);

		#if TEST>=1

		printf("\n");
		printf("order in x: "); print_int_array(order_in_X,number_of_terms);
		printf("order in y: "); print_int_array(order_in_Y,number_of_terms);
		printf("order in z: "); print_int_array(order_in_Z,number_of_terms);
		printf("Lambda: "); print_double_array(lambda_,4);
		printf("x: "); print_double_array(a1_,number_of_terms);
		printf("y: "); print_double_array(a2_,number_of_terms);

		printf("LoadInfo->(%i,%i,%i,%i)\n",li.xmin,li.xmax,li.ymin,li.ymax);

		#endif

		initializeGPU( n_, number_of_terms, order_in_X, order_in_Y, order_in_Z,
					   lambda_, a1_, a2_, bottom_plane_, rot_, setup->x_inc);

		evaluateBlock( block, x_start, y_start, z_start, x_inc, y_inc, z_inc, x_stop, y_stop, z_stop,
					   (float *)data_f, x0_src, y0_src, nx_src, ny_src, nseg, seg_size);

		#if TEST>=1

		printf("Max Error in x for tilt %i:%f for (ix,iy,isect)=(%i,%i,%i)\n",max_absolute_error1,itilt,ix_1, iy_1, iz_1);
		printf("Max Error in y for tilt %i:%f for (ix,iy,isect)=(%i,%i,%i)\n",max_absolute_error2,itilt,ix_2, iy_2, iz_2);

		#endif

	}

	printf("\n");

	#if TEST>=1	/*	Test Results	*/

	printf("Max Error in x:%f for (ix,iy,isect,itilt)=(%i,%i,%i,%i)\n",max_absolute_error1,ix_1, iy_1, iz_1,itilt_1);
	printf("Max Error in y:%f for (six,iy,isect,itilt)=(%i,%i,%i,%i)\n",max_absolute_error2,ix_2, iy_2, iz_2,itilt_2);

	#endif	/* End of the Test Results	*/

	/*	Save the reconstruction data block */

	#if TEST>=1

	printf("Embed Volume:\n");
	printf("	x0=%i	y0=%i\n", x_start, y_start);
	printf("	nx_eff=%i	ny_eff=%i\n", nx_eff, ny_eff);
	printf("in Final Destination Block:\n");
	printf("	x0_dest_f=%i	y0_dest_f=%i\n", x0_dest_f, y0_dest_f);
	printf("	nx_dest_f=%i	ny_dest_f=%i\n", nx_dest_f, ny_dest_f);

	# endif

	embedBlock( block_dest, x0_dest_f, y0_dest_f, nx_dest_f, ny_dest_f,
				block, x_start, y_start, nx_eff, ny_eff, x_inc, y_inc, blocksize_z);

	t_1 = clock();

	printf("Delta t:%E\n",ratio*(long)t_1-ratio*(long)t_0);

	free(data_f);

	free(block);

	return NO_ERROR;

}


/*
 * This routine reconstructs the 3D view between two z-values.
 */
int do_full_reconstruction(char* directory, char* basename, char* work_directory, TxBRsetup *setup, PathProjection* proj_map) {

	// Some initialization

	int error = NO_ERROR;

	double x_start = setup->x_0;
	double x_stop = setup->x_1;
	double y_start = setup->y_0;
	double y_stop = setup->y_1;
	double z_start = setup->z_0;
	double z_stop = setup->z_1;

	double x_inc = setup->x_inc;
	double y_inc = setup->y_inc;
	double z_inc = setup->z_inc;

	int blocksize_x = (x_stop-x_start+1)/x_inc/1;
	int blocksize_y = (y_stop-y_start+1)/y_inc/1;
	int blocksize_z = setup->blocksize;

	int numberOfXBlocks = ceil((x_stop-x_start)/x_inc/blocksize_x);
	int numberOfYBlocks = ceil((y_stop-y_start)/y_inc/blocksize_y);
	int numberOfZBlocks = ceil((z_stop-z_start)/z_inc/blocksize_z);

	printf("Size of blocks\n");
	printf("blocksize_x: %i\n",blocksize_x);
	printf("blocksize_y: %i\n",blocksize_y);
	printf("blocksize_z: %i\n",blocksize_z);

	// Final block

	x0_dest_f = setup->x_0;
	y0_dest_f = setup->y_0;
	nx_dest_f = (setup->x_1-setup->x_0 + 1)/x_inc;
	ny_dest_f = (setup->y_1-setup->y_0 + 1)/y_inc;

	block_dest = (float*)malloc(nx_dest_f*ny_dest_f*blocksize_z*sizeof(float));

	if (!block_dest) {
		txbr_error(stderr, "ERROR: backprojection - getting memory.\n");
		return MEMORY_ERROR;
	}

	init_float_array(&block_dest[0],nx_dest_f*ny_dest_f*blocksize_z);	/*	Initialization of the Output Image Buffer	*/

	// Define the coarse volume that will be reconstructed

	Volume *vol  = (Volume *)malloc(sizeof(Volume));

	if (!vol) {
		txbr_error(stderr, "ERROR: do_full_reconstruction - getting memory.\n");
		return MEMORY_ERROR;
	}

	vol->setup = setup;

	vol->x_start = setup->x_0;
	vol->x_stop = setup->x_1;
	vol->y_start = setup->y_0;
	vol->y_stop = setup->y_1;
	vol->z_start = setup->z_0;
	vol->z_stop = setup->z_1;

	/* Open the input file	*/

	input_file_header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!input_file_header) {
		txbr_error(stderr, "ERROR: do_full_reconstruction - getting memory.\n");
		return MEMORY_ERROR;
	}

	error = openMRCFile(directory, basename, input_file_header, &input_file);

	if (error!=NO_ERROR) return error;

	int numtlts = input_file_header->nz;

	if (numtlts>MAX_TILT) {
		txbr_error(stderr, "ERROR: do_full_reconstruction - the maximum authorized number of tilts is %i.\n", MAX_TILT);
		return (PARAMETER_ERROR);
	}

	int iblock_x,iblock_y,iblock_z;

	/*	 Start the main loop. Output blocks	*/

	for (iblock_z=0; iblock_z<numberOfZBlocks; iblock_z++) {  // block in the z direction

		vol->indexOfBlock_z = iblock_z;

		vol->z_start = z_start + z_inc*iblock_z*blocksize_z;
		vol->z_stop = TXBR_MIN(vol->z_start + z_inc*(blocksize_z-1), z_stop);

		printf("Z block %i/%i. Starts at z=%f.\n", iblock_z+1, numberOfZBlocks, vol->z_start);

		for (iblock_y=0; iblock_y<numberOfYBlocks; iblock_y++) {  // block in the y direction

			vol->indexOfBlock_y = iblock_y;

			vol->y_start = y_start + y_inc*iblock_y*blocksize_y;
			vol->y_stop = TXBR_MIN(vol->y_start + y_inc*blocksize_y-1, y_stop);

			printf("Y block %i/%i. Starts at y=%f.\n", iblock_y+1, numberOfYBlocks, vol->y_start);

			for (iblock_x=0; iblock_x<numberOfXBlocks; iblock_x++) {  // block in the x direction

				vol->indexOfBlock_x = iblock_x;

				vol->x_start = x_start + x_inc*iblock_x*blocksize_x;
				vol->x_stop = TXBR_MIN(vol->x_start + x_inc*blocksize_x-1, x_stop);

				printf("X block %i/%i. Starts at x=%f.\n", iblock_x+1, numberOfXBlocks, vol->x_start);

				if ((error = do_reconstruction(directory, basename, work_directory, vol, proj_map)) != NO_ERROR) {
					return error;
				};

			}

		}

		save_block(work_directory, basename, z_start, iblock_z, nx_dest_f, ny_dest_f, blocksize_z, &block_dest[0]);

	}

	/* Free the resources	*/

	if (fclose(input_file)) {
		txbr_error(stderr, "ERROR: do_full_reconstruction - closing file %s.\n",input_file);
		return FILE_ERROR;\
	}

	free(input_file_header);

	free(vol);

	free(block_dest);

	return error;

}

/*
 * Embed a small volume within a larger one.
 */
void embedBlock( float *block_dest, int x0_dest, int y0_dest, int nx_dest, int ny_dest,
			     float *block_src, int x0_src, int y0_src, int nx_src, int ny_src,
			     int x_inc, int y_inc, int nz) {

	// Source block (the smaller) needs to inserted in the destination block
	// Eventually do it directly to the hard drive

	int nz_dest = nz;
	int nz_src = nz;

	int delta_x = (x0_src - x0_dest)/x_inc;
	int delta_y = (y0_src - y0_dest)/y_inc;

	int i_,j_,k_;
	int x__,y__,z__;

	for (i_=0; i_<nx_src; i_++) {
		x__ = delta_x + i_;
		if (x__>=x0_dest+nx_dest) continue;
		for (j_=0; j_<ny_src; j_++) {
			y__ = delta_y + j_;
			if (y__>=y0_dest+ny_dest) continue;
			for (k_=0; k_<nz_src; k_++) {
				z__ = k_;
				block_dest[ny_dest*nx_dest*z__+nx_dest*y__+x__] = block_src[ny_src*nx_src*k_+nx_src*j_+i_];
			}
		}
	}

}