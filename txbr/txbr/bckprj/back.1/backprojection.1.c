#include <string.h>
#include <math.h>
#include <time.h>
#include <bckprj.h>

MrcHeader* input_file_header;
FILE* input_file;

double q1_[MAX_ORDER],q2_[MAX_ORDER]; /*	 For a quick evaluation of the polynoms	along 1 direction - y*/
double diff2_[MAX_ORDER][MAX_ORDER];
double bottom_plane_[3];

int n_ = 2; /*	 Order of the reconstruction	*/
int number_of_terms;

int order_in_X[MAX_COEFFICIENTS_NUMBER];
int order_in_Y[MAX_COEFFICIENTS_NUMBER];
int order_in_Z[MAX_COEFFICIENTS_NUMBER];

double lambda_[4];
double a1_[MAX_COEFFICIENTS_NUMBER];
double a2_[MAX_COEFFICIENTS_NUMBER];

double rot_[3][3];
double XYZ_[3], XYZ_rot[3];

int x_inc=1, y_inc=1, z_inc=1;

/*	 Define some test parameters	*/

#if TEST>=1

int ix_1, iy_1, iz_1, itilt_1;
int ix_2, iy_2, iz_2, itilt_2;
double ll=0, ff1 = 0, ff2 = 0;
double max_absolute_error1 = 0, max_absolute_error2 = 0;

#endif

#if SWAP_x_y==1

double tmp;

#endif

/*	 No more test parameters	*/

/*
 *
 */
double lambda(double x, double y, double z) {

//	int i;
//
//	double value = 0.0;
//
//	for (i=0; i<4; i++) {
//		value += lambda_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
//	}
//
//	return value;

	return lambda_[0] + lambda_[1]*x + lambda_[2]*y + lambda_[3]*z;

}

/*
 *
 */
double f1(double x, double y, double z) {

	int i;

//	double value = 0.0;
//
//	for (i=0; i<number_of_terms; i++) {
//		value += a1_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
//	}

	double value = a1_[0] + a1_[1]*x + a1_[2]*y + a1_[3]*z;

	for (i=4; i<number_of_terms; i++) {
		value += a1_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
	}

	return value;

}

/*
 *
 */
double f2(double x, double y, double z) {

	int i;

//	double value = 0.0;
//
//	for (i=0; i<number_of_terms; i++) {
//		value += a2_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
//	}

	double value = a2_[0] + a2_[1]*x + a2_[2]*y + a2_[3]*z;

	for (i=4; i<number_of_terms; i++) {
		value += a2_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
	}

	return value;

}

/*
 * Rotation function
 */
double* Rot(double x, double y, double z) {

	XYZ_rot[0] = rot_[0][0]*x + rot_[1][0]*y + rot_[2][0]*z;
	XYZ_rot[1] = rot_[0][1]*x + rot_[1][1]*y + rot_[2][1]*z;
	XYZ_rot[2] = rot_[0][2]*x + rot_[1][2]*y + rot_[2][2]*z;

	return &XYZ_rot[0];

}

/*
 * Inverse rotation function.
 */
double* Rot_inverse(double x, double y, double z) {

	/* The matrix is transposed - equivalent to theta in -theta	*/

	XYZ_rot[0] = rot_[0][0]*x + rot_[0][1]*y + rot_[0][2]*z;
	XYZ_rot[1] = rot_[1][0]*x + rot_[1][1]*y + rot_[1][2]*z;
	XYZ_rot[2] = rot_[2][0]*x + rot_[2][1]*y + rot_[2][2]*z;

	return &XYZ_rot[0];

}

/*
 * Initialize the rotation coefficients from affine coefficents
 * of the bottom/median plane of the sample.
 */
double* initialize_rotation_coefficients(double* plane_coefficients) {

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

	rot_[0][0] = cos_phi*cos_phi*(cos_theta-1.0) + 1.0;
	rot_[0][1] = cos_phi*sin_phi*(cos_theta-1.0);
	rot_[0][2] = -cos_phi*sin_theta;

	rot_[1][0] = cos_phi*sin_phi*(cos_theta-1.0);
	rot_[1][1] = sin_phi*sin_phi*(cos_theta-1.0) + 1.0;
	rot_[1][2] = -sin_phi*sin_theta;

	rot_[2][0] = cos_phi*sin_theta;
	rot_[2][1] = sin_phi*sin_theta;
	rot_[2][2] = cos_theta;

	return &rot_[0][0];

}

/*
 * This function swap the x, y and follows the bottom plane in z
 */
double* XYZ_t(double x, double y, double z) {

	#if SWAP_x_y==1

	tmp = x;
	x = y;
	y = tmp;

	#endif

	#if FULL_ROTATION==1

	Rot_inverse(x,y,z);

	XYZ_[0] = XYZ_rot[0];
	XYZ_[1] = XYZ_rot[1];
	XYZ_[2] = XYZ_rot[2];

	#else

	/*	Original pseudo rotation case	*/
	XYZ_[0] = x;
	XYZ_[1] = y;
	XYZ_[2] = bottom_plane_[1]*x + bottom_plane_[2]*y + z; /* Planes are defined in the wrong reference, reversed.	*/

	#endif

	return &XYZ_[0];

}

/*
 *
 */
double lambda_implicit(double x, double y, double z) {

	return lambda(x,y,z);

//	XYZ_t(x,y,z);
//
//	return lambda(XYZ_[0], XYZ_[1], XYZ_[2]);

}

/*
 *
 */
double f1_implicit(double x, double y, double z) {

	return f1(x,y,z);

//	XYZ_t(x,y,z);
//
//	return f1(XYZ_[0], XYZ_[1], XYZ_[2]);

}

/*
 *
 */
double f2_implicit(double x, double y, double z) {

	return f2(x,y,z);

//	XYZ_t(x,y,z);
//
//	return f2(XYZ_[0], XYZ_[1], XYZ_[2]);

}


/*
 * Solve the initial value problem for the quick polynomial evaluation scheme
 */
void solve_ivp(double x0, double y0, double z0, double (*f_)(double,double,double), double *ivp) {

	int i,j;

	 /*	 Initialize the array diff2_ to zeros	*/
	for (i=0; i<MAX_ORDER; i++) for (j=0; j<MAX_ORDER; j++) diff2_[i][j] = 0;

	 /*	 Calculate the first elements	*/
	for (j=0; j<n_+1; j++) {
		diff2_[0][j] = f_(x0-x_inc*j,y0,z0);
	}
	ivp[0] = diff2_[0][0];

	for (i=1; i<n_+1; i++) {
		for (j=0; j<n_-i+1; j++) {
			diff2_[i][j] = diff2_[i-1][j+1] - diff2_[i-1][j];
		}
		ivp[i] = diff2_[i][0];
		ivp[i] = (i%2==0) ? ivp[i] : -ivp[i];
	}

	return;

}

/*
 * Caculate the segment size for the propagation error to remain less than eps
 */
int calculate_segment_size( TxBRsetup *setup, PathProjection* path, int z0, int itilt, double eps) {

	int order = path->order;
	int ix_0 = setup->x_0-1;
	int ix_1 = setup->x_1-1;
	int x_inc = setup->x_inc;
	int nx_eff = (ix_1-ix_0)/x_inc + 1;
	int nx_lim = ix_0 + x_inc*nx_eff;

	/*	Make sure f1_implicit is defined for itilt	*/

	memcpy(a1_, &path->coefficients_1[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));

	int m, segmentSize=0;
	int x0=setup->x_0, y0=setup->y_0;
	double value, value_app, diff=0.0;

	solve_ivp(x0,y0,z0,&f1_implicit,&q1_[0]);

	while (diff<eps && segmentSize<nx_lim) {

		value_app = q1_[0];
		value = f1_implicit(x0,y0,z0);

		diff = fabs(value_app-value);

		#if TEST>=1
		printf("%f	%f	%e\n", value_app, value, diff);
		#endif

		for (m=order-1;m>=0;m--) {
			q1_[m] += q1_[m+1];
		}

		x0 = x0 + x_inc;

		segmentSize++;

	}

	printf("Segment Size used for this recursion algorithm: %i\n", segmentSize);

	return segmentSize;

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
					Volume* vol, PathProjection* proj_map) {

	printf("Reconstruction for %s from z=%f to z=%f\n", basename, vol->z_start, vol->z_stop);

	// Boundaries of the Volume to reconstruct

	int x_start = (int)vol->x_start;
	int y_start = (int)vol->y_start;
	int z_start = (int)vol->z_start;
	int x_stop = (int)vol->x_stop;
	int y_stop = (int)vol->y_stop;
	int z_stop = (int)vol->z_stop;

	TxBRsetup *setup = vol->setup;

	x_inc = setup->x_inc;
	y_inc = (int)setup->y_inc;
	z_inc = (int)setup->z_inc;

	int x0 = setup->x_0;
	int y0 = setup->y_0;
	int x1 = setup->x_1;
	int y1 = setup->y_1;

	int blocksize = (z_stop-z_start+1)/z_inc;
	int indexOfBlock = vol->indexOfBlock_z;

	printf("x0=%i	y0=%i\n",x0,y0);
	printf("x1=%i	y1=%i\n",x1,y1);

	/* The z_start variable purpose is mainly for file storage	*/

	int error = NO_ERROR;

	clock_t t_0,t_1;
	float ratio = 1./CLOCKS_PER_SEC;

	t_0 = clock();

	if (blocksize>BLOCK_SIZE_MAX) {
		printf("Block size %i is higher than recommanded %i.\n", blocksize, BLOCK_SIZE_MAX);
		return (PARAMETER_ERROR);
	}

	memcpy(bottom_plane_, setup->plane_coeffs2, sizeof(setup->plane_coeffs2));

	initialize_rotation_coefficients(&bottom_plane_[0]);

	printf("rot[0][0]=%4.5f\trot[0][1]=%4.5f\trot[0][2]=%4.5f\t\n",rot_[0][0],rot_[0][1],rot_[0][2]);
	printf("rot[0][0]=%4.5f\trot[0][1]=%4.5f\trot[0][2]=%4.5f\t\n",rot_[1][0],rot_[1][1],rot_[1][2]);
	printf("rot[0][0]=%4.5f\trot[0][1]=%4.5f\trot[0][2]=%4.5f\t\n",rot_[2][0],rot_[2][1],rot_[2][2]);

	n_ = proj_map->order;
	number_of_terms = number_of_polynomial_terms(n_);

	init_power_order(n_, &order_in_X[0], &order_in_Y[0], &order_in_Z[0]);

	int nx = input_file_header->nx;
	int ny = input_file_header->ny;
	int numtlts = input_file_header->nz;

	printf("(nx,ny,ntilts)->(%i,%i,%i)\n",nx,ny,numtlts);

    rewind(input_file_header->fp);
    fseek(input_file_header->fp, input_file_header->headerSize, SEEK_SET);

	if (nx>X_MAX || ny>Y_MAX) {

		printf("Image size (%i,%i) is higher than the recommanded (%i,%i).\n", nx, ny, X_MAX,Y_MAX);

	}

	int nx_eff = (x_stop-x_start+1)/x_inc;
	int ny_eff = (y_stop-y_start+1)/y_inc;

	float* block = (float*)malloc(nx_eff*ny_eff*blocksize*sizeof(float));

	if (!block) {
		txbr_error(stderr, "ERROR: backprojection - getting memory.\n");
		return MEMORY_ERROR;
	}

	init_float_array(&block[0],nx_eff*ny_eff*blocksize);	/*	Initialization of the Output Image Buffer	*/

	printf("(nx_eff,ny_eff,ntilts)->(%i,%i,%i)\n",nx_eff,ny_eff,numtlts);

//	Islice* slice = NULL;
//	float* slice_data_f = NULL;

	float* slice_data_f = (float*)malloc(nx*ny*sizeof(float));

	if (!slice_data_f) {
		txbr_error(stderr, "ERROR: backprojection - getting memory for data_f.\n");
		return MEMORY_ERROR;
	}

	int itilt,isect,ix,iy,index,m;
	double z_section;

	int seg_size = nx_eff;

	seg_size = calculate_segment_size( setup, proj_map, setup->z_0, 0, EPSILON );	//	itilt=0 here

	int nseg = floor(nx_eff/seg_size) + 1;

	double lambda, dely_lambda;
	int isegment, istart, istop;
	double XYZ[3],xproj,yproj,xdiff,ydiff;
	long xv,yv;

	for (itilt=0; itilt<numtlts; itilt++) {

		if (proj_map->skipView[itilt]==1) {
			printf("Exposure %i has been skipped!",itilt);
			continue;
		}

		 /*	 Read input image	*/

		printf("%s %i\0",itilt==0 ? "Reading Image #" : ",",itilt);
		fflush(stdout);

//		slice = sliceReadMRC(input_file_header,itilt,'z');
//		slice_data_f = slice->data.f;

	    mrcReadFloatSlice(slice_data_f, input_file_header, itilt);

		/*	 Copy the bundle adjustment parameters	*/

		memcpy(lambda_, &proj_map->lambda[itilt][0], 4*sizeof(double));
		memcpy(a1_, &proj_map->coefficients_1[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));
		memcpy(a2_, &proj_map->coefficients_2[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));

		#if TEST>=1

		printf("\n");
		printf("order in x: "); print_int_array(order_in_X,number_of_terms);
		printf("order in y: "); print_int_array(order_in_Y,number_of_terms);
		printf("order in z: "); print_int_array(order_in_Z,number_of_terms);
		printf("Lambda: "); print_double_array(lambda_,4);
		printf("x: "); print_double_array(a1_,number_of_terms);
		printf("y: "); print_double_array(a2_,number_of_terms);

		#endif

		for (isect=0; isect<blocksize; isect++) {	/*	 Step through the output block	*/

			z_section = z_start + z_inc*isect;	/*	 z-section coordinate of the output slice */

			index = isect*nx_eff*ny_eff;	/*	counter for the block	*/

			/*	 Loop along y image coordinate */

			 for (iy=y0; iy<y_stop+1; iy=iy+y_inc) {

				 /*	 Now increment along x direction. Iterative evaluation echeme for each segment.	*/

				 for (isegment=0; isegment<nseg; isegment++) {

					 istart = x0 + isegment*seg_size;
					 istop = TXBR_MIN( istart+seg_size, x_stop );

					 #if OFFSET==1

					 XYZ[0] = istart;
					 XYZ[1] = iy;
					 XYZ[2] = z_section;

					 #elif

					 XYZ[0] = istart - 1;
					 XYZ[1] = iy - 1;
					 XYZ[2] = z_section;

					 #endif

					 lambda = lambda_implicit(XYZ[0],XYZ[1],XYZ[2]);
					 dely_lambda = lambda_implicit(XYZ[0]+x_inc,XYZ[1],XYZ[2])-lambda;

//					 xproj = f1_implicit(XYZ[0],XYZ[1],XYZ[2]);
//					 yproj = f2_implicit(XYZ[0],XYZ[1],XYZ[2]);
//
//					 solve_ivp(XYZ[0],XYZ[1],XYZ[2],&f1_implicit,&q1_[0]);
//					 solve_ivp(XYZ[0],XYZ[1],XYZ[2],&f2_implicit,&q2_[0]);

					 xproj = f1(XYZ[0],XYZ[1],XYZ[2]);
					 yproj = f2(XYZ[0],XYZ[1],XYZ[2]);

					 solve_ivp(XYZ[0],XYZ[1],XYZ[2],&f1,&q1_[0]);
					 solve_ivp(XYZ[0],XYZ[1],XYZ[2],&f2,&q2_[0]);

					 for (ix=istart; ix<istop+1; ix=ix+x_inc) {

						xproj = q1_[0];
						yproj = q2_[0];

						#if TEST==10	/*	Test the accuracy	*/

						XYZ[0] = ix;
						XYZ[1] = iy;
						XYZ[2] = z_section;

						ll = lambda_implicit(XYZ[0],XYZ[1],XYZ[2]);
						ff1 = f1_implicit(XYZ[0],XYZ[1],XYZ[2]);
						ff2 = f2_implicit(XYZ[0],XYZ[1],XYZ[2]);

						if (fabs((xproj-ff1))>max_absolute_error1) {
							ix_1 = ix;
							iy_1 = iy;
							iz_1 = isect;
							itilt_1 = itilt;
						}

						if (fabs((yproj-ff2))>max_absolute_error2) {
							ix_2 = ix;
							iy_2 = iy;
							iz_2 = isect;
							itilt_2 = itilt;
						}

						max_absolute_error1 = TXBR_MAX(fabs((xproj-ff1)),max_absolute_error1);
						max_absolute_error2 = TXBR_MAX(fabs((yproj-ff2)),max_absolute_error2);

						if (ix%40==0 & iy==30) printf("%i %i %i -> %f %f %f\n",ix,iy,isect,xproj,ff1,fabs((xproj-ff1)));
						if (ix%40==0 & iy==30) printf("%i %i %i -> %f %f %f\n",ix,iy,isect,yproj,ff2,fabs((yproj-ff2)));
						if (ix%40==0 & iy==30) printf("%i %i %i -> %f %f %f\n",ix,iy,isect,lambda,ll,fabs((lambda-ll)));

						#endif	/* end of test	*/

						xproj /= lambda;
						yproj /= lambda;

						/* Eventually reswap xproj and yproj like in the original code	*/

						#if SWAP_x_y==1

						tmp = xproj;
						xproj = yproj;
						yproj = tmp;

						#endif

						#if OFFSET==0

						xproj -= 1.0;
						yproj -= 1.0;

						#endif


//						xv = floor(xproj);
//						yv = floor(yproj);

						xv = (long)(xproj+1.0) - 1;
						yv = (long)(yproj+1.0) - 1;

						xdiff = xproj - xv;	/*	 Weight for Linear Interpolation	*/
						ydiff = yproj - yv;

						/*	 Calculate linear interpolation error, making sure that we are within the coordinate limits	*/

						if (xv>=0 && xv<nx-1 && yv>=0 && yv<ny-1) {
							block[index] = block[index] + INTERP2(slice_data_f, nx, ny, xv, yv, xdiff, ydiff);
						}

						for (m=n_-1;m>=0;m--) {
							q1_[m] += q1_[m+1];
							q2_[m] += q2_[m+1];
						}

						lambda += dely_lambda;

						index++;

					}

				 }

			}

		}

		#if TEST>=1

		printf("Max Error in x for tilt %i:%f for (ix,iy,isect)=(%i,%i,%i)\n",max_absolute_error1,itilt,ix_1, iy_1, iz_1);
		printf("Max Error in y for tilt %i:%f for (ix,iy,isect)=(%i,%i,%i)\n",max_absolute_error2,itilt,ix_2, iy_2, iz_2);

		#endif

//		slice_data_f = NULL;
//		if (slice) sliceFree(slice);	/* do not move	*/

	}

	free(slice_data_f);

	printf("\n");

	#if TEST>=1	/*	Test Results	*/

	printf("Max Error in x:%f for (ix,iy,isect,itilt)=(%i,%i,%i,%i)\n",max_absolute_error1,ix_1, iy_1, iz_1,itilt_1);
	printf("Max Error in y:%f for (ix,iy,isect,itilt)=(%i,%i,%i,%i)\n",max_absolute_error2,ix_2, iy_2, iz_2,itilt_2);

	#endif	/* End of the Test Results	*/

	/*	Save the reconstruction data	*/

	save_block(work_directory, basename, setup->z_0, indexOfBlock, nx_eff, ny_eff, blocksize, &block[0]);

	t_1 = clock();

	printf("Delta t:%E\n",ratio*(long)t_1-ratio*(long)t_0);

	free(block);

	return NO_ERROR;

}


/*
 * This routine reconstructs the 3D view between two z-values.
 */
int do_full_reconstruction(char* directory, char* basename, char* work_directory, TxBRsetup *setup, PathProjection* proj_map) {

	double z_start = setup->z_0;
	double z_stop = setup->z_1;

	/*	If z_start or z_stop are float instead of double, they are put to 0	*/

	int blocksize_z = setup->blocksize;
	double z_inc = setup->z_inc;

	/* Make sure the number of blocks is correct*/

	setup->numberOfBlocks = ceil((z_stop - z_start + 1)/z_inc/blocksize_z);

	int numberOfBlocks = setup->numberOfBlocks;

	double depth = (z_stop - z_start);	 /*	 Distance in increment unit between output planes	*/

	printf("Start the reconstruction...\n");

	printf("z start: %f\n",z_start);
	printf("z stop: %f\n",z_stop);
	printf("Block Size: %i\n",blocksize_z);
	printf("Number of Blocks: %i\n",numberOfBlocks);
	printf("Depth in z: %f\n",depth);
	printf("\n");

	int error = NO_ERROR;

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

	int iblock_z;

	/*	 Start the main loop. Output blocks	*/

	for (iblock_z=0; iblock_z<numberOfBlocks; iblock_z++) {

		vol->indexOfBlock_z = iblock_z;

		vol->z_start = z_start + z_inc*iblock_z*blocksize_z;
		vol->z_stop = TXBR_MIN(vol->z_start + z_inc*(blocksize_z-1), z_stop);

		printf("Block %i/%i. Starts at z=%f.\n", iblock_z+1, numberOfBlocks, vol->z_start);

		if ((error = do_reconstruction(directory, basename, work_directory, vol, proj_map))!=NO_ERROR) {

			return error;

		};

	}

	/* Free the resources	*/

	if(fclose(input_file)) {

		txbr_error(stderr, "ERROR: do_full_reconstruction - closing file %s.\n",input_file);

		return FILE_ERROR;

	}

	free(input_file_header);

	/*	Merge the files together	*/

	// merge_blocks(work_directory, basename, z_start, numberOfBlocks);

	return error;

}

