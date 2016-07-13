#include <string.h>
#include <math.h>
#include <time.h>
#include <bckprj.h>

MrcHeader* input_file_header;
MrcHeader* output_file_header;
FILE* input_file;
FILE* output_file;

double q1_[MAX_ORDER],q2_[MAX_ORDER]; /*    For a quick evaluation of the polynoms along 1 direction (y) */
double diff2_[MAX_ORDER][MAX_ORDER];
double bottom_plane_[3];

int order_in_X[MAX_COEFFICIENTS_NUMBER];
int order_in_Y[MAX_COEFFICIENTS_NUMBER];
int order_in_Z[MAX_COEFFICIENTS_NUMBER];

double lambda_[4];
double a1_[MAX_COEFFICIENTS_NUMBER];
double a2_[MAX_COEFFICIENTS_NUMBER];

double rot_[3][3];

float *slice_data_f, *block;

#if TILT_NORMALIZATION==1 | TILT_NORMALIZATION==2 | SLICE_NORMALIZATION==1 | SLICE_NORMALIZATION==2
float *norm;
#endif

void init_short_int_array(short int*, int );
void init_float_array(float*, int );

int mrcReadFloatSlice(b3dFloat *buf, MrcHeader *hdata, int slice);

#define FLOOR_PLUS( x )	(long)(x+1.0) - 1

#define ROT( x, y, z)																		\
	x_ = x, y_=y, z_ = z;																	\
	x = rot_[0][0]*x_ + rot_[1][0]*y_ + rot_[2][0]*z_;										\
	y = rot_[0][1]*x_ + rot_[1][1]*y_ + rot_[2][1]*z_;										\
	z = rot_[0][2]*x_ + rot_[1][2]*y_ + rot_[2][2]*z_;										\


#define ROT_INV( x, y, z)																	\
	x_ = x, y_=y, z_ = z;																	\
	x = rot_[0][0]*x_ + rot_[0][1]*y_ + rot_[0][2]*z_;										\
	y = rot_[1][0]*x_ + rot_[1][1]*y_ + rot_[1][2]*z_;										\
	z = rot_[2][0]*x_ + rot_[2][1]*y_ + rot_[2][2]*z_;										\

#define SHEAR( x, y, z)																		\
	z = bottom_plane_[1]*x + bottom_plane_[2]*y + z;

#define NO_TRANSFORMATION_MODE 0
#define XYSWAP_ROTATION_MODE 1
#define XYSWAP_SHEAR_MODE 2
#define ROTATION_MODE 3
#define SHEAR_MODE 4


#define TRANSF_0( x, y, z)										\
	x_ = x; y_ = y; z_ = z;                                                                         \


#define TRANSF_1( x, y, z)										\
	x_ = y; y_ = x; z_ = z;										\
	ROT_INV( x_, y_, z_)                                                                            \


#define TRANSF_2( x, y, z)										\
	x_ = y; y_ = x; z_ = z;                                                                 	\
	SHEAR( x_, y_, z_)                                                                              \


#define TRANSF_3( x, y, z)										\
	ROT_INV( x, y, z)


#define TRANSF_4( x, y, z)										\
	SHEAR( x, y, z)


#define LAMBDA( lambda, x, y, z, mode)                                                                  \
	TRANSF_##mode( x, y, z)										\
	lambda = lambda_[0] + lambda_[1]*x_ + lambda_[2]*y_ + lambda_[3]*z_;                            \


#define F1( f1, x, y, z, mode)										\
	TRANSF_##mode( x, y, z)										\
	f1 = a1_[0] + a1_[1]*x_ + a1_[2]*y_ + a1_[3]*z_;						\
	for (i=4; i<number_of_terms; i++) {								\
		f1 += a1_[i]*pow(x_,order_in_X[i])*pow(y_,order_in_Y[i])*pow(z_,order_in_Z[i]);		\
	}                                                                                               \


#define F2( f2, x, y, z, mode)                                                                      	\
	TRANSF_##mode( x, y, z)										\
	f2 = a2_[0] + a2_[1]*x_ + a2_[2]*y_ + a2_[3]*z_;                                                \
	for (i=4; i<number_of_terms; i++) {								\
		f2 += a2_[i]*pow(x_,order_in_X[i])*pow(y_,order_in_Y[i])*pow(z_,order_in_Z[i]);		\
	}


#define DEL_LAMBDA( x, y, z, dx, mode)														\
	LAMBDA( lambda, x, y, z, mode)									\
	LAMBDA( dely_lambda, x + dx, y, z, mode)							\
	dely_lambda -= lambda;                                                                          \

#define DEL_F1( x, y, z, dx, mode)                                                      \
	F1( xproj, x, y, z, mode)																\
	F1( dely_xproj, x + dx, y, z, mode)						\
	dely_xproj -= xproj;


#define DEL_F2( x, y, z, dx, mode)                      				\
	F2( yproj, x, y, z, mode)							\
	F2( dely_yproj, x + dx, y, z, mode)                     			\
	dely_yproj -= yproj;


#define IVP( x0, y0, z0, axis, ivp, mode)													\
																							\
	 /*	 Initialize the array diff2_ to zeros	*/											\
	for (i=0; i<MAX_ORDER; i++) for (j=0; j<MAX_ORDER; j++) diff2_[i][j] = 0;				\
																							\
	 /*	 Calculate the first elements	*/													\
	for (j=0; j<n_+1; j++) {																\
		F##axis(diff2_[0][j],x0-x_inc*j,y0,z0, mode);										\
	}																						\
	ivp[0] = diff2_[0][0];																	\
																							\
	for (i=1; i<n_+1; i++) {																\
		for (j=0; j<n_-i+1; j++) {															\
			diff2_[i][j] = diff2_[i-1][j+1] - diff2_[i-1][j];								\
		}																					\
		ivp[i] = diff2_[i][0];																\
		ivp[i] = (i%2==0) ? ivp[i] : -ivp[i];												\
	}																						\


#define INIT_PRJ()                                                                          \
                                                                                            \
    /* Boundaries of the Volume to reconstruct	*/                                          \
                                                                                            \
    int x_start = (int)vol->x_start;                                                        \
    int y_start = (int)vol->y_start;                                                        \
    int z_start = (int)vol->z_start;                                                        \
    int x_stop = (int)vol->x_stop;                                                          \
    int y_stop = (int)vol->y_stop;                                                          \
    int z_stop = (int)vol->z_stop;                                                          \
                                                                                            \
    TxBRsetup *setup = vol->setup;                                                          \
                                                                                            \
    int x_inc = (int)setup->x_inc;                                                          \
    int y_inc = (int)setup->y_inc;                                                          \
    int z_inc = (int)setup->z_inc;                                                          \
                                                                                            \
    int x0 = setup->x_0;                                                                    \
    int y0 = setup->y_0;                                                                    \
                                                                                            \
    int itilt_start = proj_map->itilt_start;                                                \
    int itilt_stop = proj_map->itilt_stop;                                                  \
    int blocksize = itilt_stop-itilt_start;                                                 \
    int indexOfBlock = itilt_start/blocksize;                                               \
                                                                                            \
    /* The z_start variable purpose is mainly for file storage	*/                          \
                                                                                            \
    int error = NO_ERROR;                                                                   \
                                                                                            \
    clock_t t_0,t_1;                                                                        \
    float ratio = 1./CLOCKS_PER_SEC;                                                        \
                                                                                            \
    t_0 = clock();                                                                          \
                                                                                            \
    int number_of_terms = proj_map->number_of_terms;                                        \
    int nslices = input_file_header->nz;                                                    \
                                                                                            \
    rewind(input_file_header->fp);                                                          \
    fseek(input_file_header->fp, input_file_header->headerSize, SEEK_SET);                  \
                                                                                            \
    /* Destination image characteristics    */                                              \
                                                                                            \
    int nx = proj_map->nx;                                                                  \
    int ny = proj_map->ny;                                                                  \
    int numtlts = proj_map->numberOfTilts;                                                  \
                                                                                            \
    /* Final block characteristics	*/                                                  \
                                                                                            \
    int nx_eff = (x_stop-x_start+1)/x_inc;                                                  \
    int ny_eff = (y_stop-y_start+1)/y_inc;                                                  \
                                                                                            \
    /*	Initialization of the Output Image Buffer	*/                                  \
    init_float_array( &block[0], nx*ny*blocksize );                                         \
    /* Segment size calculated for tilt_ref=numtlts/2	*/                                  \
                                                                                            \
    double z_med = (setup->z_0+setup->z_1)/2.0;                                             \
    int tilt_ref = numtlts/2;                                                               \
    int seg_size =  use_fixed_segment_size ? FIX_SEGMENT_SIZE:                              \
            calculate_segment_size( setup, proj_map, z_med, tilt_ref, EPSILON );            \
    int nseg = FLOOR_PLUS(nx_eff/seg_size) + 1;                                             \


#define INIT_BCKPRJ()                                                                 		\
												\
	/* Boundaries of the Volume to reconstruct	*/					\
												\
	int x_start = (int)vol->x_start;							\
	int y_start = (int)vol->y_start;							\
	int z_start = (int)vol->z_start;							\
	int x_stop = (int)vol->x_stop;								\
	int y_stop = (int)vol->y_stop;								\
	int z_stop = (int)vol->z_stop;								\
												\
	TxBRsetup *setup = vol->setup;								\
												\
	int x_inc = (int)setup->x_inc;								\
	int y_inc = (int)setup->y_inc;								\
	int z_inc = (int)setup->z_inc;								\
												\
	int x0 = setup->x_0;									\
	int y0 = setup->y_0;									\
												\
	int blocksize = (z_stop-z_start+1)/z_inc;						\
	int indexOfBlock = vol->indexOfBlock_z;							\
												\
	/* The z_start variable purpose is mainly for file storage	*/			\
												\
	int error = NO_ERROR;									\
												\
	clock_t t_0,t_1;									\
	float ratio = 1./CLOCKS_PER_SEC;							\
												\
	t_0 = clock();										\
												\
	int number_of_terms = proj_map->number_of_terms;					\
												\
	/* Image source characteristics	*/							\
												\
	int nx = input_file_header->nx;								\
	int ny = input_file_header->ny;								\
	int numtlts = input_file_header->nz;							\
												\
	rewind(input_file_header->fp);								\
	fseek(input_file_header->fp, input_file_header->headerSize, SEEK_SET);                  \
												\
	/* Final block characteristics	*/							\
												\
	int nx_eff = (x_stop-x_start+1)/x_inc;                                                  \
	int ny_eff = (y_stop-y_start+1)/y_inc;                                                 	\
												\
	/*	Initialization of the Output Image Buffer	*/                              \
        init_float_array( &block[0], nx_eff*ny_eff*blocksize );                                 \
	/* Segment size calculated for tilt_ref=numtlts/2	*/                              \
												\
	double z_med = (setup->z_0+setup->z_1)/2.0;                                             \
	int tilt_ref = numtlts/2;								\
        int seg_size =  use_fixed_segment_size ? FIX_SEGMENT_SIZE:                              \
                   calculate_segment_size( setup, proj_map, z_med, tilt_ref, EPSILON );         \
                                                                                                \
	int nseg = FLOOR_PLUS(nx_eff/seg_size) + 1;												\


// ---------------------------------------------------------------------------------------------//


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

	printf("rot[0][0]=%4.5f\trot[0][1]=%4.5f\trot[0][2]=%4.5f\t\n",rot_[0][0],rot_[0][1],rot_[0][2]);
	printf("rot[0][0]=%4.5f\trot[0][1]=%4.5f\trot[0][2]=%4.5f\t\n",rot_[1][0],rot_[1][1],rot_[1][2]);
	printf("rot[0][0]=%4.5f\trot[0][1]=%4.5f\trot[0][2]=%4.5f\t\n",rot_[2][0],rot_[2][1],rot_[2][2]);

	return &rot_[0][0];

}


/*
 * Caculate the segment size for the propagation error to remain less than eps
 */
int calculate_segment_size( TxBRsetup *setup, PathProjection* proj_map, int z0, int itilt, double eps) {

	int order = proj_map->order;

	int ix_0 = setup->x_0-1;
	int ix_1 = setup->x_1-1;
	int x_inc = setup->x_inc;

	int nx_eff = (ix_1-ix_0)/x_inc + 1;
	//int nx_lim = ix_0 + x_inc*nx_eff;

	int n_ = proj_map->order;
	int number_of_terms = proj_map->number_of_terms;

	int i,j;
	double x_,y_,z_;

	/*	Make sure f1_implicit is defined for itilt	*/

	memcpy(a1_, &proj_map->coefficients_1[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));

	int m, segmentSize=0;
	int x0=setup->x_0, y0=setup->y_0;
	double value, value_app, diff=0.0;

	IVP(x0,y0,z0,1,q1_,0)

	//while (diff<eps && segmentSize<nx_lim) {
        while (diff<eps && x0<=setup->x_1) {

		value_app = q1_[0];
		F1(value,x0,y0,z0,0)

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

/*
        segmentSize = 20;
*/ 
	printf("Segment Size used for this recursion algorithm: %i (nx_eff=%i)\n", segmentSize, nx_eff);

	return segmentSize;

}

/****************************************************************************
* Function:     do_project;
*
* Argument:     directory	The directory for this project
* 		basename	The basename for this reconstruction project
* 		work_directory	The work directory
* 		vol		The coarse volume that will be reconstructed
* 		proj_map	The projection Map
*
* Returns:      void
*
* Description:  This routine reconstructs a tomographic volume from a set of
*  				projection maps.
*
*****************************************************************************/

int do_project( char* directory, char* basename, char* work_directory, int use_fixed_segment_size,
		Volume* vol, PathProjection* proj_map) {
    
    printf("Project the volume\n");

    int n_ = proj_map->order;

    INIT_PRJ()

    #if SLICE_NORMALIZATION==1 | SLICE_NORMALIZATION==2
    init_float_array( &norm[0], nx*ny*blocksize );
    #endif

    int itilt, isect, isegment, index;
    long ix, iy, istart, istop;
    long xv, yv;
    double lambda, dely_lambda;
    double xproj, yproj;
    double dely_xproj, dely_yproj;
    double XYZ[3], xdiff, ydiff, z_section;

    int i;
    double x_,y_,z_;
    float value_;

    printf("tilt: %i->%i  y: %i->%i  y: %i->%i\n",itilt_start,itilt_stop,x0,x_stop,y0,y_stop);
    printf("nx: %i  ny: %i\n",nx,ny);
    printf("nx_eff:%i  ny_eff:%i\n",nx_eff,ny_eff);
    printf("nslices: %i\n",nslices); 

    for ( isect=0; isect<nslices; isect++) {  /*  Read input volume	*/

        printf("%s %i", isect==0 ? "Reading Slice #" : ",", isect);
        fflush(stdout);

        z_section = z_start + z_inc*isect; /*	 z-section coordinate of the output slice */

        mrcReadFloatSlice(slice_data_f, input_file_header, isect);

        for (itilt=itilt_start; itilt<itilt_stop; itilt++) { /*	 Step through the output block	*/

            if (proj_map->skipView[itilt]==1) {
                continue;
            }

            /*	 Copy the bundle adjustment parameters	*/

            memcpy(lambda_, &proj_map->lambda[itilt][0], 4 * sizeof (double));
            memcpy(a1_, &proj_map->coefficients_1[itilt][0], MAX_COEFFICIENTS_NUMBER * sizeof (double));
            memcpy(a2_, &proj_map->coefficients_2[itilt][0], MAX_COEFFICIENTS_NUMBER * sizeof (double));

            /*	 Loop along Y image coordinate */

            for (iy = y0; iy < y_stop + 1; iy = iy + y_inc) {

                /*	 Now increment along x direction.	\
                        Iterative evaluation echeme for each segment.	*/

                for (isegment = 0; isegment < nseg; isegment++) {

                    istart = x0 + isegment*seg_size;
                    istop = TXBR_MIN(istart + seg_size - 1, x_stop);

                    XYZ[0] = istart;
                    XYZ[1] = iy;
                    XYZ[2] = z_section;

                    DEL_LAMBDA(XYZ[0], XYZ[1], XYZ[2], x_inc, 0)
                    DEL_F1(XYZ[0], XYZ[1], XYZ[2], x_inc, 0)
                    DEL_F2(XYZ[0], XYZ[1], XYZ[2], x_inc, 0)

                    xproj -= 1; /*  Correct for the offset   */
                    yproj -= 1;

                    for (ix=istart; ix<istop+1; ix=ix+x_inc) {

                        /* Eventually reswap xproj and yproj like in the original code	*/

                        #if SWAP_x_y==1

                        tmp = xproj;
                        xproj = yproj;
                        yproj = tmp;

                        #endif

                        xv = FLOOR_PLUS(xproj);
                        yv = FLOOR_PLUS(yproj);

                        xdiff = xproj - xv; /*	 Weight for Linear Interpolation	*/
                        ydiff = yproj - yv;

                        /*	 Calculate linear interpolation error, making sure that we are within the coordinate limits	*/

                        // Split the value
                        if (xv >= 0 && xv < nx - 1 && yv >= 0 && yv < ny - 1 ) {

                            index = (itilt-itilt_start)*nx*ny + yv*nx + xv;
                            value_ = sliceGetPixelMagnitude_float(slice_data_f,nx_eff,ix-x0,iy-y0);

                            block[index] = block[index] + (1.0-xdiff)*(1.0-ydiff)*value_;
                            block[index+1] = block[index+1] + xdiff*(1.0-ydiff)*value_;
                            block[index+nx] = block[index+nx] + (1.0-xdiff)*ydiff*value_;
                            block[index+nx+1] = block[index+nx+1] + xdiff*ydiff*value_;

                            #if SLICE_NORMALIZATION==1 | SLICE_NORMALIZATION==2
                            norm[index] = norm[index] + (1.0-xdiff)*(1.0-ydiff);
                            norm[index+1] = norm[index+1] + xdiff*(1.0-ydiff);
                            norm[index+nx] = norm[index+nx] + (1.0-xdiff)*ydiff;
                            norm[index+nx+1] = norm[index+nx+1] + xdiff*ydiff;
                            #endif

                        }

                        lambda += dely_lambda;
                        xproj += dely_xproj;
                        yproj += dely_yproj;

                    }

                }

            }

        }

    }

    printf("\n");

    #if SLICE_NORMALIZATION==1
    for (index=0; index<nx*ny*blocksize; index++) {
        block[index] = norm[index]<=MINIMUM_SLICE_NUMBER_FOR_NORMALIZATION ? 0.0 : block[index] / norm[index];
    }
    #elif SLICE_NORMALIZATION==2
    for (index=0; index<nx*ny*blocksize; index++) {
        block[index] = norm[index]!=0 ? block[index]/norm[index] : 0.0;
    }
    #else
    for (index=0; index < nx*ny*blocksize; index++) {
        block[index] / nslices;
    }
    #endif

    /*	Save the reconstruction data	*/

    //save_projection_block( work_directory, basename, indexOfBlock, nx, ny, blocksize, &block[0] );

    char bs[FILENAME_LEN];
    sprintf(bs, "%s_%s", basename, &(proj_map->label));
    
    save_projection_block( work_directory, &bs, indexOfBlock, nx, ny, blocksize, &block[0] );

    t_1 = clock();

    printf("Delta t:%E\n", ratio * (long) t_1 - ratio * (long) t_0);

    return error;

    
}



/****************************************************************************
* Function:     do_reconstruction; do_reconstruction_order1;
* 		do_reconstruction_order1_noscaling
*
* Argument:     directory	The directory for this project
* 		basename	The basename for this reconstruction project
* 		work_directory	The work directory
* 		vol		The coarse volume that will be reconstructed
* 		proj_map	The projection Map
*
* Returns:      void
*
* Description:  This routine reconstructs a tomographic volume from a set of
*  				projection maps.
*
*****************************************************************************/


int do_reconstruction( char* directory, char* basename, char* work_directory, int use_fixed_segment_size,
					Volume* vol, PathProjection* proj_map) {

	int n_ = proj_map->order;

	INIT_BCKPRJ()

        #if TILT_NORMALIZATION==1 | TILT_NORMALIZATION==2
        init_float_array( &norm[0], nx_eff*ny_eff*blocksize );
        #endif

	int itilt, isect, isegment, index;
	int m;
	long ix, iy, istart, istop;
	long xv, yv;
	double lambda, dely_lambda;
	double xproj, yproj;
	double XYZ[3], xdiff, ydiff, z_section;

	int i,j;
	double x_,y_,z_;

	for (itilt=0; itilt<numtlts; itilt++) {

                if (proj_map->skipView[itilt]==0) {
                    printf( "%s%i" , itilt==0 ? "Reading Image #" : ", ", itilt);
                } else {
                    printf( "%s%i]" , itilt==0 ? "Reading Image # [" : ", [", itilt);
                }
		fflush(stdout);

                if (proj_map->skipView[itilt]==1) {
                    continue;
		}

		 /*	 Read input image	*/


                mrcReadFloatSlice(slice_data_f, input_file_header, itilt);

		/*	 Copy the bundle adjustment parameters	*/

		memcpy(lambda_, &proj_map->lambda[itilt][0], 4*sizeof(double));
		memcpy(a1_, &proj_map->coefficients_1[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));
		memcpy(a2_, &proj_map->coefficients_2[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));

		for (isect=0; isect<blocksize; isect++) {	/*	 Step through the output block	*/

			z_section = z_start + z_inc*isect;	/*	 z-section coordinate of the output slice */

			index = isect*nx_eff*ny_eff;	/*	counter for the block	*/

			/*	 Loop along y image coordinate */

			for (iy=y0; iy<y_stop+1; iy=iy+y_inc) {

				/*	 Now increment along x direction.	\
				 	Iterative evaluation echeme for each segment.	*/

				for (isegment=0; isegment<nseg; isegment++) {

					istart = x0 + isegment*seg_size;
					istop = TXBR_MIN( istart+seg_size-1, x_stop );

					XYZ[0] = istart;
					XYZ[1] = iy;
					XYZ[2] = z_section;

					DEL_LAMBDA( XYZ[0], XYZ[1], XYZ[2], x_inc, 0)

					F1( xproj, XYZ[0],XYZ[1],XYZ[2],0)
					F2( yproj, XYZ[0],XYZ[1],XYZ[2],0)

					IVP(XYZ[0],XYZ[1],XYZ[2],1,q1_,0)
					IVP(XYZ[0],XYZ[1],XYZ[2],2,q2_,0)

					for (ix=istart; ix<istop+1; ix=ix+x_inc) {

						xproj = q1_[0];
						yproj = q2_[0];

						xproj /= lambda;
						yproj /= lambda;

						xproj -= 1;	/* Correct for the offset	*/
						yproj -= 1;

						/* Eventually reswap xproj and yproj like in the original code	*/

						#if SWAP_x_y==1

						tmp = xproj;
						xproj = yproj;
						yproj = tmp;

						#endif

						xv = FLOOR_PLUS(xproj);
						yv = FLOOR_PLUS(yproj);

						xdiff = xproj - xv;	/*	 Weight for Linear Interpolation	*/
						ydiff = yproj - yv;

						/*	 Calculate linear interpolation error, making sure that we are within the coordinate limits	*/

						if (xv>=0 && xv<nx-1 && yv>=0 && yv<ny-1) {
							block[index] = block[index] + INTERP2(slice_data_f, nx, ny, xv, yv, xdiff, ydiff);
                                                        #if TILT_NORMALIZATION==1 | TILT_NORMALIZATION==2
                                                        norm[index] = norm[index] + 1;
                                                        #endif
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

	}

	printf("\n");

        #if TILT_NORMALIZATION==1
        for (index=0;index<nx_eff*ny_eff*blocksize;index++) {
            block[index] = norm[index]<=MINIMUM_TILT_NUMBER_FOR_NORMALIZATION ? 0.0 : block[index]*numtlts/norm[index];
        }
        #elif TILT_NORMALIZATION==2
        for (index=0;index<nx_eff*ny_eff*blocksize;index++) {
            block[index] = norm[index]!=0 ? block[index]/norm[index] : 0.0;
        }
        #endif

	/*	Save the reconstruction data	*/

        save_backprojection_block( work_directory, basename, indexOfBlock, setup->x_0,
                                   setup->y_0, setup->z_0, nx_eff, ny_eff, blocksize, setup->sx,
                                   setup->sy, setup->sz, &block[0] );

	t_1 = clock();

	printf("Delta t:%E\n",ratio*(long)t_1-ratio*(long)t_0);

	return error;

}


int do_reconstruction_order1( char* directory, char* basename, char* work_directory, int use_fixed_segment_size,
					Volume* vol, PathProjection* proj_map) {

	INIT_BCKPRJ()

        #if TILT_NORMALIZATION==1 | TILT_NORMALIZATION==2
        init_float_array( &norm[0], nx_eff*ny_eff*blocksize );
        #endif

	int itilt, isect, isegment, index;
	long ix, iy, istart, istop;
	long xv, yv;
	double lambda, dely_lambda;
	double xproj, yproj;
	double dely_xproj, dely_yproj;
	double XYZ[3], xdiff, ydiff, z_section;

	int i;
	double x_,y_,z_;

	for (itilt=0; itilt<numtlts; itilt++) {

		if (proj_map->skipView[itilt]==1) {
			printf("-%i-",itilt);
			continue;
		}

		 /*	 Read input image	*/

		printf("%s %i",itilt==0 ? "Reading Image #" : ",",itilt);
		fflush(stdout);

                mrcReadFloatSlice(slice_data_f, input_file_header, itilt);

		/*	 Copy the bundle adjustment parameters	*/

		memcpy(lambda_, &proj_map->lambda[itilt][0], 4*sizeof(double));
		memcpy(a1_, &proj_map->coefficients_1[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));
		memcpy(a2_, &proj_map->coefficients_2[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));

		for (isect=0; isect<blocksize; isect++) {	/*	 Step through the output block	*/

			z_section = z_start + z_inc*isect;	/*	 z-section coordinate of the output slice */

			index = isect*nx_eff*ny_eff;	/*	counter for the block	*/

			/*	 Loop along y image coordinate */

			for (iy=y0; iy<y_stop+1; iy=iy+y_inc) {

				/*	 Now increment along x direction.	\
				 	Iterative evaluation echeme for each segment.	*/

				for (isegment=0; isegment<nseg; isegment++) {

					istart = x0 + isegment*seg_size;
					istop = TXBR_MIN( istart+seg_size-1, x_stop );

					XYZ[0] = istart;
					XYZ[1] = iy;
					XYZ[2] = z_section;

					DEL_LAMBDA( XYZ[0], XYZ[1], XYZ[2], x_inc, 0)
					DEL_F1( XYZ[0], XYZ[1], XYZ[2], x_inc, 0)
					DEL_F2( XYZ[0], XYZ[1], XYZ[2], x_inc, 0)

					for (ix=istart; ix<istop+1; ix=ix+x_inc) {

						xproj /= lambda;
						yproj /= lambda;

						/* Eventually reswap xproj and yproj like in the original code	*/

						#if SWAP_x_y==1

						tmp = xproj;
						xproj = yproj;
						yproj = tmp;

						#endif

						xv = FLOOR_PLUS(xproj);
						yv = FLOOR_PLUS(yproj);

						xdiff = xproj - xv;	/*	 Weight for Linear Interpolation	*/
						ydiff = yproj - yv;

						xv -= 1;	/* Correct for the offset	*/
						yv -= 1;

						/*	 Calculate linear interpolation error, making sure that we are within the coordinate limits	*/

						if (xv>=0 && xv<nx-1 && yv>=0 && yv<ny-1) {
							block[index] = block[index] + INTERP2(slice_data_f, nx, ny, xv, yv, xdiff, ydiff);
                                                        #if TILT_NORMALIZATION==1 | TILT_NORMALIZATION==2
                                                        norm[index] = norm[index] + 1;
                                                        #endif
						}

						lambda += dely_lambda;
						xproj += dely_xproj;
						yproj += dely_yproj;

						index++;

					}

				 }

			}

		}

	}

	printf("\n");

        #if TILT_NORMALIZATION==1
        for (index=0;index<nx_eff*ny_eff*blocksize;index++) {
            block[index] = norm[index]<=MINIMUM_TILT_NUMBER_FOR_NORMALIZATION ? 0.0 : block[index]*numtlts/norm[index];
        }
        #elif TILT_NORMALIZATION==2
        for (index=0;index<nx_eff*ny_eff*blocksize;index++) {
            block[index] = norm[index]!=0.0 ? block[index]/norm[index] : 0.0;
        }
        #endif

	/*	Save the reconstruction data	*/

        save_backprojection_block( work_directory, basename, indexOfBlock, setup->x_0,
                                   setup->y_0, setup->z_0, nx_eff, ny_eff, blocksize, setup->sx,
                                   setup->sy, setup->sz, &block[0] );

	t_1 = clock();

	printf("Delta t:%E\n",ratio*(long)t_1-ratio*(long)t_0);

	return error;

}


int do_reconstruction_order1_noscaling( char* directory, char* basename, char* work_directory, int use_fixed_segment_size,
					Volume* vol, PathProjection* proj_map) {

	INIT_BCKPRJ()

        #if TILT_NORMALIZATION==1 | TILT_NORMALIZATION==2
            init_float_array( &norm[0], nx_eff*ny_eff*blocksize );
        #endif

	int itilt, isect, isegment, index;
	long ix, iy, istart, istop;
	long xv, yv;
	double lambda, dely_lambda;
	double xproj, yproj;
	double dely_xproj, dely_yproj;
	double XYZ[3], xdiff, ydiff, z_section;

	int i;
	double x_,y_,z_;

	for (itilt=0; itilt<numtlts; itilt++) {

		if (proj_map->skipView[itilt]==1) {
			printf("-%i-",itilt);
			continue;
		}

		 /*	 Read input image	*/

		printf("%s %i",itilt==0 ? "Reading Image #" : ",",itilt);
		fflush(stdout);

                mrcReadFloatSlice(slice_data_f, input_file_header, itilt);

		/*	 Copy the bundle adjustment parameters	*/

		memcpy(lambda_, &proj_map->lambda[itilt][0], 4*sizeof(double));
		memcpy(a1_, &proj_map->coefficients_1[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));
		memcpy(a2_, &proj_map->coefficients_2[itilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double));

		for (isect=0; isect<blocksize; isect++) {	/*	 Step through the output block	*/

			z_section = z_start + z_inc*isect;	/*	 z-section coordinate of the output slice */

			index = isect*nx_eff*ny_eff;	/*	counter for the block	*/

			/*	 Loop along y image coordinate */

			for (iy=y0; iy<y_stop+1; iy=iy+y_inc) {

				/*	 Now increment along x direction.	\
				 	Iterative evaluation echeme for each segment.	*/

				for (isegment=0; isegment<nseg; isegment++) {

					istart = x0 + isegment*seg_size;
					istop = TXBR_MIN( istart+seg_size-1, x_stop );

					XYZ[0] = istart;
					XYZ[1] = iy;
					XYZ[2] = z_section;

					DEL_LAMBDA( XYZ[0], XYZ[1], XYZ[2], x_inc, 0)
					DEL_F1( XYZ[0], XYZ[1], XYZ[2], x_inc, 0)
					DEL_F2( XYZ[0], XYZ[1], XYZ[2], x_inc, 0)

					xproj -= 1;	/* Correct for the offset	*/
					yproj -= 1;

					for (ix=istart; ix<istop+1; ix=ix+x_inc) {

						/* Eventually reswap xproj and yproj like in the original code	*/

						#if SWAP_x_y==1

						tmp = xproj;
						xproj = yproj;
						yproj = tmp;

						#endif

						xv = FLOOR_PLUS(xproj);
						yv = FLOOR_PLUS(yproj);

						xdiff = xproj - xv;	/*	 Weight for Linear Interpolation	*/
						ydiff = yproj - yv;


						/*	 Calculate linear interpolation error, making sure that we are within the coordinate limits	*/

						if (xv>=0 && xv<nx-1 && yv>=0 && yv<ny-1) {
							block[index] = block[index] + INTERP2(slice_data_f, nx, ny, xv, yv, xdiff, ydiff);
                                                        #if TILT_NORMALIZATION==1 | TILT_NORMALIZATION==2
                                                        norm[index] = norm[index] + 1.0;
                                                        #endif
						}

						lambda += dely_lambda;
						xproj += dely_xproj;
						yproj += dely_yproj;

						index++;

					}

				 }

			}

		}

	}
 
	printf("\n");

        #if TILT_NORMALIZATION==1
        for (index=0;index<nx_eff*ny_eff*blocksize;index++) {
            block[index] = norm[index]<=MINIMUM_TILT_NUMBER_FOR_NORMALIZATION ? 0.0 : block[index]*numtlts/norm[index];
        }
        #elif TILT_NORMALIZATION==2
        for (index=0;index<nx_eff*ny_eff*blocksize;index++) {
            block[index] = norm[index]!=0.0 ? block[index]/norm[index] : 0.0;
        }
        #else
        for (index=0;index<nx_eff*ny_eff*blocksize;index++) {
            block[index] /= numtlts;
        }
        #endif

	/*	Save the reconstruction data	*/

        save_backprojection_block( work_directory, basename, indexOfBlock, setup->x_0,
                                   setup->y_0, setup->z_0, nx_eff, ny_eff, blocksize, setup->sx,
                                   setup->sy, setup->sz, &block[0] );

	t_1 = clock();

	printf("Delta t:%E\n",ratio*(long)t_1-ratio*(long)t_0);

	return error;

}

// --------------------------------------------------------------------------------------------- //
/*
 * This routine reconstructs the 3D view between two z-values.
 */
int do_full_reconstruction( char* directory, char* basename, char* work_directory, TxBRsetup *setup,
        PathProjection* proj_map ) {

    double z_start = setup->z_0;    /*  If z_start or z_stop are float instead of double, they are put to 0 */
    double z_stop = setup->z_1;    
    double z_inc = setup->z_inc;
    double depth = (z_stop-z_start); /*   Distance in increment unit between output planes    */

    int blocksize = setup->blocksize;
    int numberOfBlocks = ceil((z_stop-z_start+1)/z_inc/blocksize);
    int error = NO_ERROR;
    int use_fixed_segment_size = setup->use_fixed_segment_size;
    
    setup->numberOfBlocks = numberOfBlocks;
    
    int n_ = proj_map->order;   /*	Power Order */
    init_power_order(n_, &order_in_X[0], &order_in_Y[0], &order_in_Z[0]);

    memcpy(bottom_plane_, setup->plane_coeffs2, sizeof (setup->plane_coeffs2)); /*  Boundary Planes */
    initialize_rotation_coefficients(&bottom_plane_[0]);

    if (blocksize > BLOCK_SIZE_MAX) {
        printf("Block size %i is higher than recommanded %i.\n", blocksize, BLOCK_SIZE_MAX);
        return (PARAMETER_ERROR);
    }

    printf("Start the reconstruction...\n");

    printf("z start: %f\n", z_start);
    printf("z stop: %f\n", z_stop);
    printf("Block Size: %i\n", blocksize);
    printf("Number of Blocks: %i\n", numberOfBlocks);
    printf("Depth in z: %f\n", depth);
    printf("\n");

    // Define the coarse volume that will be reconstructed

    Volume *vol = (Volume *) malloc(sizeof (Volume));

    if (!vol) {
        txbr_error(stderr, "ERROR: do_full_reconstruction - getting memory.\n");
        return MEMORY_ERROR;
    }

    vol->setup = setup;

    vol->x_start = setup->x_0;
    vol->y_start = setup->y_0;
    vol->z_start = setup->z_0;
    
    vol->x_stop = setup->x_1;
    vol->y_stop = setup->y_1;
    vol->z_stop = setup->z_1;

    /* Open the input file	*/

    input_file_header = (MrcHeader *) malloc(sizeof (MrcHeader));

    if (!input_file_header) {
        txbr_error(stderr, "ERROR: do_full_reconstruction - getting memory.\n");
        return MEMORY_ERROR;

    }

    error = openMRCFile(directory, basename, input_file_header, &input_file);

    if (error!=NO_ERROR) return error;

    int numtlts = input_file_header->nz;

    if (numtlts > MAX_TILT) {
        txbr_error(stderr, "ERROR: do_full_reconstruction - the maximum authorized number of tilts is %i.\n", MAX_TILT);
        return (PARAMETER_ERROR);
    }

    /*	Scaling term or not	*/

    int scaling = 0;
    int itilt;

    for (itilt = 0; itilt < numtlts; itilt++) {
        if (proj_map->lambda[itilt][0] != 1.0 || proj_map->lambda[itilt][1] != 0.0 ||
                proj_map->lambda[itilt][1] != 0.0 || proj_map->lambda[itilt][1] != 0.0) scaling = 1;
    }

    /* Image source	*/

    int nx = input_file_header->nx;
    int ny = input_file_header->ny;

    if (nx > X_MAX || ny > Y_MAX) {
        printf("Image size (%i,%i) is higher than the recommanded (%i,%i).\n", nx, ny, X_MAX, Y_MAX);
    }

    slice_data_f = (float*) malloc(nx * ny * sizeof (float));

    if (!slice_data_f) {
        txbr_error(stderr, "ERROR: backprojection - getting memory for data_f.\n");
        return MEMORY_ERROR;
    }

    int nx_eff = (vol->x_stop - vol->x_start + 1) / setup->x_inc;
    int ny_eff = (vol->y_stop - vol->y_start + 1) / setup->y_inc;

    block = (float*) malloc(nx_eff * ny_eff * blocksize * sizeof (float));

    if (!block) {
        txbr_error(stderr, "ERROR: backprojection - getting memory.\n");
        return MEMORY_ERROR;
    }

    #if TILT_NORMALIZATION==1 | TILT_NORMALIZATION==2
    norm = (float*) malloc(nx_eff * ny_eff * blocksize * sizeof (float));

    if (!norm) {
        txbr_error(stderr, "ERROR: backprojection - getting memory.\n");
        return MEMORY_ERROR;
    }
    #endif

    int iblock_z;

    /*	 Start the main loop. Output blocks	*/

    for (iblock_z = 0; iblock_z < numberOfBlocks; iblock_z++) {

        vol->indexOfBlock_z = iblock_z;

        vol->z_start = z_start + z_inc * iblock_z*blocksize;
        vol->z_stop = TXBR_MIN(vol->z_start + z_inc * (blocksize - 1), z_stop);

        printf("Block %i/%i. From z=%f to z=%f.\n", iblock_z + 1, numberOfBlocks, vol->z_start, vol->z_stop);

        if (n_ == 1 && scaling == 0) error = do_reconstruction_order1_noscaling(directory, basename, work_directory, use_fixed_segment_size, vol, proj_map);
        else if (n_ == 1 && scaling != 0) error = do_reconstruction_order1(directory, basename, work_directory, use_fixed_segment_size, vol, proj_map);
        else error = do_reconstruction(directory, basename, work_directory, use_fixed_segment_size, vol, proj_map);

        if (error!=NO_ERROR) return error;

    }

    /* Close file and free the resources	*/

    if (fclose(input_file)) {
        txbr_error(stderr, "ERROR: do_full_reconstruction - closing file %s.\n", input_file);
        return FILE_ERROR;
    }

    free(input_file_header);
    free(slice_data_f);
    free(block);

    #if TILT_NORMALIZATION==1 | TILT_NORMALIZATION==2
    free(norm);
    #endif

    return error;

}


/*
 * This routine projects a volume on a tilt series.
 */
int do_full_projection( char* directory, char* basename, char* work_directory,
                        TxBRsetup *setup, PathProjection* proj_map ) {

    printf("Full projection\n");

    int error = NO_ERROR;
    char imodfile_in[FILENAME_LEN];
    char imodfile_out[FILENAME_LEN];

    int iblock_tilt;
    int numtlts = proj_map->numberOfTilts;
    int blocksize = setup->blocksize;
    int numberOfBlocks = 1 + (numtlts-1)/setup->blocksize;
    int use_fixed_segment_size = setup->use_fixed_segment_size;

    setup->numberOfBlocks = numberOfBlocks;
    proj_map->numberOfTilts = numtlts;

    int n_ = proj_map->order; /*  Power Order */
    init_power_order(n_, &order_in_X[0], &order_in_Y[0], &order_in_Z[0]);

    memcpy(bottom_plane_, setup->plane_coeffs2, sizeof (setup->plane_coeffs2)); /*	Boundary Planes	*/
    initialize_rotation_coefficients(&bottom_plane_[0]);

    if (numtlts>MAX_TILT) {
        txbr_error(stderr, "ERROR: do_full_projection - the maximum authorized number of tilts is %i.\n", MAX_TILT);
        return (PARAMETER_ERROR);
    }

    if (blocksize > BLOCK_SIZE_MAX) {
        printf("Block size %i is higher than recommanded %i.\n", blocksize, BLOCK_SIZE_MAX);
        return (PARAMETER_ERROR);
    }

    printf("Number of tilts: %i\n", numtlts);
    printf("Block Size: %i\n", blocksize);
    printf("Number of Blocks: %i\n", numberOfBlocks);
    printf("\n");

    // Define the volume to be reprojected

    Volume *vol = (Volume *) malloc(sizeof (Volume));

    if (!vol) {
        txbr_error(stderr, "ERROR: do_full_projection - getting memory.\n");
        return MEMORY_ERROR;
    }

    vol->setup = setup;

    vol->x_start = setup->x_0;
    vol->x_stop = setup->x_1;
    vol->y_start = setup->y_0;
    vol->y_stop = setup->y_1;
    vol->z_start = setup->z_0;
    vol->z_stop = setup->z_1;

    /*  Input file - Data source  */

    sprintf(imodfile_in, "%s/%s_z_%.1f.out", directory, basename, vol->z_start);
    input_file_header = (MrcHeader *) malloc(sizeof (MrcHeader));

    if (!input_file_header) {
        txbr_error(stderr, "ERROR: do_full_projection - getting memory.\n");
        return MEMORY_ERROR;
    }

    printf("Reading the File %s...\n", imodfile_in);

    if ((input_file = fopen(imodfile_in, "rbe")) == NULL) {
        txbr_error(stderr, "ERROR: do_full_projection - cannot open file %s.\n", imodfile_in);
        return FILE_ERROR;
    }

    printf("Read the MRC header...\n");

    mrc_head_read( input_file, input_file_header );

    printf("Finish Read the MRC header...\n");

    printf("\n");

    int nx_eff = (vol->x_stop - vol->x_start + 1) / setup->x_inc;
    int ny_eff = (vol->y_stop - vol->y_start + 1) / setup->y_inc;

    slice_data_f = (float*) malloc(nx_eff * ny_eff * sizeof (float));

    if (!slice_data_f) {
        txbr_error(stderr, "ERROR: do_full_projection - getting memory for data_f.\n");
        return MEMORY_ERROR;
    }

    /*  Ouput file - Data destination  */

    output_file_header = (MrcHeader *) malloc(sizeof (MrcHeader));

    if (!output_file_header) {
        txbr_error(stderr, "ERROR: do_full_projection - getting memory.\n");
        return MEMORY_ERROR;
    }

    sprintf(imodfile_out, "%s/%s.st", work_directory, basename);
//    sprintf(imodfile_out, "%s/%s_z_%.1f.out", directory, basename, vol->z_start);

    if ((output_file = fopen(imodfile_out, "rbe"))==NULL) {
        txbr_error(stderr, "ERROR: do_full_projection - cannot open file %s.\n", imodfile_out);
        return FILE_ERROR;
    }

    mrc_head_read(output_file, output_file_header);

    int nx = output_file_header->nx;
    int ny = output_file_header->ny;

    if (nx >X_MAX || ny>Y_MAX) {
        printf("Image size (%i,%i) is higher than the recommanded (%i,%i).\n", nx, ny, X_MAX, Y_MAX);
    }

    if (fclose(output_file)) {
        txbr_error(stderr, "ERROR: do_full_reconstruction - closing file %s.\n", output_file);
        return FILE_ERROR;
    }

    proj_map->nx =nx;
    proj_map->ny =ny;

    free(output_file_header);

    block = (float*) malloc(nx*ny*blocksize*sizeof(float));

    if (!block) {
        txbr_error(stderr, "ERROR: do_full_projection - getting memory.\n");
        return MEMORY_ERROR;
    }

    #if SLICE_NORMALIZATION==1 | SLICE_NORMALIZATION==2
    norm  = (float *)malloc(nx*ny*blocksize*sizeof(float));
    
    if (!norm) {
        txbr_error(stderr, "ERROR: do_full_projection - getting memory.\n");
        return MEMORY_ERROR;
    }
    #endif

    /*  Start the main loop. Output blocks	*/

    for (iblock_tilt=0; iblock_tilt<numberOfBlocks; iblock_tilt++) {

        proj_map->itilt_start = iblock_tilt*blocksize;
        proj_map->itilt_stop = TXBR_MIN((iblock_tilt + 1) * blocksize, numtlts);

        printf("Block %i/%i. From tilt=%i to tilt=%i.\n", iblock_tilt+1, numberOfBlocks, proj_map->itilt_start, proj_map->itilt_stop);

        do_project(directory, basename, work_directory, use_fixed_segment_size, vol, proj_map);

        if (error!=NO_ERROR) return error;

    }

    /* Close input file and free the resources  */

    if (fclose(input_file)) {
        txbr_error(stderr, "ERROR: do_full_reconstruction - closing file %s.\n", input_file);
        return FILE_ERROR;
    }

    free(input_file_header);
    free(slice_data_f);
    free(block);

    #if SLICE_NORMALIZATION==1 | SLICE_NORMALIZATION==2
    free(norm);
    #endif

    return error;

}

