#include <math.h>
#include <string.h>
#include <stdio.h>

#include <bckprj.h>

// CPU version of the backprojection functions
double f1_implicit_cb(double x, double y, double z, int number_of_terms, 
			int* order_in_X, int* order_in_Y, int* order_in_Z, 
			double rot_[3][3], double *bottom_plane_, double *a1_);
double f2_implicit_cb(double x, double y, double z, int number_of_terms, 
			int* order_in_X, int* order_in_Y, int* order_in_Z, 
			double rot_[3][3], double *bottom_plane_, double *a2_);

void XYZ_t_cb(double *XYZ, double  rot_[3][3], double *bottom_plane_);
void Rot_cb(double *XYZ, double rot_[3][3]);
void Rot_inverse_cb(double *XYZ, double rot_[3][3]);



/*
 * Caculate the Box to use!
 */
void calculate_box( int x_start, int y_start, int z_start, int x_stop, int y_stop, int z_stop,
		int nx, int ny, int *x0_src, int *y0_src, int *nx_src, int *ny_src, 
		int number_of_terms, int *order_in_X, int *order_in_Y, int *order_in_Z,
		double rot_[3][3], double* bottom_plane_, double *a1_, double *a2_) {

	// Very crude and naive approximation to find the boundary box
	// To change at some point
	int padding = (int)(0.01*TXBR_MAX(nx,ny));

	float xmin = 0.0, xmax = 0.0;
	float ymin = 0.0, ymax = 0.0;

	float x[8],y[8];

	#if OFFSET==0

	x_start -= 1;
	y_start -= 1;
	z_start -= 1;
	x_stop -= 1;
	y_stop -= 1;
	z_stop -= 1;

	#endif

	x[0] = f1_implicit_cb(x_start,y_start,z_start, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a1_);
	x[1] = f1_implicit_cb(x_start,y_start,z_stop, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a1_);
	x[2] = f1_implicit_cb(x_start,y_stop,z_start, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a1_);
	x[3] = f1_implicit_cb(x_start,y_stop,z_stop, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a1_);
	x[4] = f1_implicit_cb(x_stop,y_start,z_start, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a1_);
	x[5] = f1_implicit_cb(x_stop,y_start,z_stop, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a1_);
	x[6] = f1_implicit_cb(x_stop,y_stop,z_start, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a1_);
	x[7] = f1_implicit_cb(x_stop,y_stop,z_stop, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a1_);

	y[0] = f2_implicit_cb(x_start,y_start,z_start, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a2_);
	y[1] = f2_implicit_cb(x_start,y_start,z_stop, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a2_);
	y[2] = f2_implicit_cb(x_start,y_stop,z_start, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a2_);
	y[3] = f2_implicit_cb(x_start,y_stop,z_stop, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a2_);
	y[4] = f2_implicit_cb(x_stop,y_start,z_start, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a2_);
	y[5] = f2_implicit_cb(x_stop,y_start,z_stop, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a2_);
	y[6] = f2_implicit_cb(x_stop,y_stop,z_start, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a2_);
	y[7] = f2_implicit_cb(x_stop,y_stop,z_stop, number_of_terms, order_in_X, order_in_Y, order_in_Z, rot_, bottom_plane_, a2_);

	xmin = (float)x[0];
	xmax = (float)x[0];

	int i=0;

	for (i=1;i<8;i++) {
		xmin = (xmin<x[i]) ? xmin : x[i];
		xmax = (xmax>x[i]) ? xmax : x[i];
	}

	xmin -= padding;
	xmax += padding;

	xmin = TXBR_MAX(xmin,OFFSET);
	xmax = TXBR_MIN(xmax,nx-1+OFFSET);

	ymin = (float)y[0];
	ymax = (float)y[0];

	for (i=1;i<8;i++) {
		ymin = (ymin<y[i]) ? ymin : y[i];
		ymax = (ymax>y[i]) ? ymax : y[i];
	}

	ymin -= padding;
	ymax += padding;

	ymin = TXBR_MAX(ymin,OFFSET);
	ymax = TXBR_MIN(ymax,ny-1+OFFSET);

	#if TEST>0

	printf("X boundaries:\n");
	printf("(x,y,z)=(%i,%i,%i) -> %f\n", x_start, y_start, z_start, x[0]);
	printf("(x,y,z)=(%i,%i,%i) -> %f\n", x_start, y_start, z_stop, x[1]);
	printf("(x,y,z)=(%i,%i,%i) -> %f\n", x_start, y_stop, z_start, x[2]);
	printf("(x,y,z)=(%i,%i,%i) -> %f\n", x_start, y_stop, z_stop, x[3]);
	printf("(x,y,z)=(%i,%i,%i) -> %f\n", x_stop, y_start, z_start, x[4]);
	printf("(x,y,z)=(%i,%i,%i) -> %f\n", x_stop, y_start, z_stop, x[5]);
	printf("(x,y,z)=(%i,%i,%i) -> %f\n", x_stop, y_stop, z_start, x[6]);
	printf("(x,y,z)=(%i,%i,%i) -> %f\n", x_stop, y_stop, z_stop, x[7]);

	#endif

	*x0_src = (int)xmin;
	*y0_src = (int)ymin;

	*nx_src = (int)(xmax - xmin + 1);
	*ny_src = (int)(ymax - ymin + 1);

	#if TEST>0

	printf("(xmin,xmax)=(%f,%f)\n", xmin, xmax);
	printf("(ymin,ymax)=(%f,%f)\n", ymin, ymax);

	printf("(x0_src,y0_src)=(%i,%i)\n",*x0_src,*y0_src);
	printf("(nx_src,ny_src)=(%i,%i)\n",*nx_src,*ny_src);

	#endif

}


/*
 *
 */
double f1_implicit_cb(double x, double y, double z, 
	int number_of_terms, int* order_in_X, int* order_in_Y, int* order_in_Z, 
	double rot_[3][3], double *bottom_plane_, double *a1_) {

	double XYZ[3];
	
	XYZ[0] = x;
	XYZ[1] = y;
	XYZ[2] = z;
	
	XYZ_t_cb(XYZ, rot_, bottom_plane_);
	// Old codel; is unrolled -> return f1(XYZ_[0], XYZ_[1], XYZ_[2]);
	
	int i;
	double value = 0.0;
	
	for (i=0; i<number_of_terms; i++) {
		value += a1_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
	}
	
	return value;

}

/*
 *
 */
double f2_implicit_cb(double x, double y, double z, 
	int number_of_terms, int* order_in_X, int* order_in_Y, int* order_in_Z, 
	double rot_[3][3], double *bottom_plane_, double *a2_) {

	double XYZ[3];
		
	XYZ[0] = x;
	XYZ[1] = y;
	XYZ[2] = z;
	
	XYZ_t_cb(XYZ, rot_, bottom_plane_);

	// Old codel; is unrolled -> return f2(XYZ_[0], XYZ_[1], XYZ_[2]);
	int i;
	double value = 0.0;
	
	for (i=0; i<number_of_terms; i++) {
		value += a2_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
	}
	
	return value;

}

/*
 * This function swap the x, y and follows the bottom plane in z
 */
void XYZ_t_cb(double *XYZ, double  rot_[3][3], double *bottom_plane_) {

	double x = XYZ[0];
	double y = XYZ[1];
	double z = XYZ[2];
	
	#if SWAP_x_y==1

	tmp = x;
	x = y;
	y = tmp;

	#endif

	#if FULL_ROTATION==1

	Rot_inverse_cb(XYZ, rot_);

//	XYZ_[0] = XYZ_rot[0];
//	XYZ_[1] = XYZ_rot[1];
//	XYZ_[2] = XYZ_rot[2];

	#else

	/*	Original pseudo rotation case	*/
	XYZ[0] = x;
	XYZ[1] = y;
	XYZ[2] = bottom_plane_[1]*x + bottom_plane_[2]*y + z; /* Planes are defined in the wrong reference, reversed.	*/

	#endif

	//return &XYZ_[0];

}



/*
 * Rotation function
 */
void Rot_cb(double *XYZ, double rot_[3][3]) {

	double x = XYZ[0];
	double y = XYZ[1];
	double z = XYZ[2];
	
	XYZ[0] = rot_[0][0]*x + rot_[1][0]*y + rot_[2][0]*z;
	XYZ[1] = rot_[0][1]*x + rot_[1][1]*y + rot_[2][1]*z;
	XYZ[2] = rot_[0][2]*x + rot_[1][2]*y + rot_[2][2]*z;

	//return &XYZ_rot[0];
}

/*
 * Inverse rotation function.
 */
void Rot_inverse_cb(double *XYZ, double rot_[3][3]) {

	/* The matrix is transposed - equivalent to theta in -theta	*/

	double x = XYZ[0];
	double y = XYZ[1];
	double z = XYZ[2];
	
	XYZ[0] = rot_[0][0]*x + rot_[0][1]*y + rot_[0][2]*z;
	XYZ[1] = rot_[1][0]*x + rot_[1][1]*y + rot_[1][2]*z;
	XYZ[2] = rot_[2][0]*x + rot_[2][1]*y + rot_[2][2]*z;

}

