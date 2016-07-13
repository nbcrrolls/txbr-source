#include <math.h>
#include "bckprj.h"

double q1_[MAX_ORDER];
double q2_[MAX_ORDER];
double diff2_[MAX_ORDER][MAX_ORDER];

double XYZ_[3], XYZ_rot[3];

double lambda(double x, double y, double z);
double f1(double x, double y, double z);
double f2(double x, double y, double z);
double* Rot(double x, double y, double z);
double* Rot_inverse(double x, double y, double z);
double* XYZ_t(double x, double y, double z);
double lambda_implicit(double x, double y, double z);
double f1_implicit(double x, double y, double z);
double f2_implicit(double x, double y, double z);
void solve_ivp1(double x0, double y0, double z0, double *ivp1);
void solve_ivp2(double x0, double y0, double z0, double *ivp2);

// EXTERN variables
extern int number_of_terms, n_;
extern long order_in_X[MAX_COEFFICIENTS_NUMBER], order_in_Y[MAX_COEFFICIENTS_NUMBER], order_in_Z[MAX_COEFFICIENTS_NUMBER];
extern double a1_[MAX_COEFFICIENTS_NUMBER], a2_[MAX_COEFFICIENTS_NUMBER];
extern double lambda_[4];
extern double rot_[3][3];
extern double bottom_plane_[3];
extern int x_inc, y_inc;

/*
 *
 */
double lambda(double x, double y, double z) {

	int i;
	double value = 0.0;

	value = lambda_[0];
	for (i=1; i<4; i++) {
		value += lambda_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
	}

	return value;

}


double f1(double x, double y, double z) {

	int i;
	double value = 0.0;

	for (i=0; i<number_of_terms; i++) {
		value += a1_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
	}

	return value;

}


double f2(double x, double y, double z) {

	int i;

	double value = 0.0;

	for (i=0; i<number_of_terms; i++) {
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

	XYZ_t(x,y,z);

	return lambda(XYZ_[0], XYZ_[1], XYZ_[2]);

}

/*
 *
 */
double f1_implicit(double x, double y, double z) {

	XYZ_t(x,y,z);

	return f1(XYZ_[0], XYZ_[1], XYZ_[2]);

}

/*
 *
 */
double f2_implicit(double x, double y, double z) {

	XYZ_t(x,y,z);

	return f2(XYZ_[0], XYZ_[1], XYZ_[2]);

}


/*
 * Solve the initial value problem for the quick polynomial evaluation scheme
 */
void solve_ivp1(double x0, double y0, double z0, double *ivp) {

	int i,j;

	 /*	 Initialize the array diff2_ to zeros	*/
	for (i=0; i<MAX_ORDER; i++) for (j=0; j<MAX_ORDER; j++) diff2_[i][j] = 0;

	 /*	 Calculate the first elements	*/
	for (j=0; j<n_+1; j++) {
		diff2_[0][j] = f1_implicit(x0-j*x_inc,y0,z0);
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
 * Solve the initial value problem for the quick polynomial evaluation scheme
 */
void solve_ivp2(double x0, double y0, double z0, double *ivp) {

	int i,j;

	 /*	 Initialize the array diff2_ to zeros	*/
	for (i=0; i<MAX_ORDER; i++) for (j=0; j<MAX_ORDER; j++) diff2_[i][j] = 0;

	 /*	 Calculate the first elements	*/
	for (j=0; j<n_+1; j++) {
		diff2_[0][j] = f2_implicit(x0-j*x_inc,y0,z0);
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


void evaluateBlock( float *block, int x0, int y0, int z0, int x_inc, int y_inc, int z_inc, int x1, int y1, int z1,
				   float* slice_data_f, int x0_src, int y0_src, int nx_src, int ny_src, int nseg, int seg_size) {

	int nx_dest = floor((x1-x0+1)/x_inc);
	int ny_dest = floor((y1-y0 + 1)/y_inc);
	int blocksize = floor((z1-z0+1)/z_inc);

	#if TEST>=1

	printf("\nEvaluate block routine ->\n");
	printf("	(x0,y0,z0)=(%i,%i,%i)\n",x0,y0,z0);
	printf("	(x_inc,y_inc,z_inc)=(%i,%i,%i)\n",x_inc,y_inc,z_inc);
	printf("	(x1,y1,z1)=(%i,%i,%i)\n",x1,y1,z1);
	printf("	(x0_src,y0_src)=(%i,%i)\n",x0_src,y0_src);
	printf("	(nx_src,ny_src)=(%i,%i)\n",nx_src,ny_src);
	printf("	(nx_dest,ny_dest)=(%i,%i)\n",nx_dest,ny_dest);
	printf("	nseg=%i\n",nseg);
	printf("	seg_size=%i\n",seg_size);

	#endif

	int isegment;
	int ix, iy, isect, istart, istop;
	double XYZ[3],z_section;
	double lam, dely_lam;

	double xproj, yproj, xdiff, ydiff;
	long xv = 0, yv = 0;
	int index, m;

	for (isect=0; isect<blocksize; isect++) {	/*	 Step through the output block	*/

		z_section = z0 + z_inc*isect;	/*	 z-section coordinate of the output slice */

		index = isect*nx_dest*ny_dest;	/*	counter for the block	*/

		/*	 Loop along y image coordinate */

		for (iy=y0; iy<y1+1; iy=iy+y_inc) {
			
		/*	 Now increment along x direction. Iterative evaluation scheme for each segment.	*/

			for (isegment=0; isegment<nseg; isegment++) {

				istart = x0 + isegment*seg_size;
				istop = TXBR_MIN( istart + seg_size, x1+1 );

				#if OFFSET==1

				XYZ[0] = istart;
				XYZ[1] = iy;
				XYZ[2] = z_section;

				#else

				XYZ[0] = istart - 1;
				XYZ[1] = iy - 1;
				XYZ[2] = z_section;

				#endif

				xproj = f1_implicit(XYZ[0],XYZ[1],XYZ[2]);
				yproj = f2_implicit(XYZ[0],XYZ[1],XYZ[2]);

				lam = lambda_implicit(XYZ[0],XYZ[1],XYZ[2]);
				dely_lam = lambda_implicit(XYZ[0]+x_inc,XYZ[1],XYZ[2])-lam;

				solve_ivp1(XYZ[0],XYZ[1],XYZ[2],&q1_[0]);
				solve_ivp2(XYZ[0],XYZ[1],XYZ[2],&q2_[0]);

				//for (ix=istart; ix<istop+1; ix=ix+x_inc) {
				for (ix=istart; ix<istop; ix=ix+x_inc) {
					
					xproj = q1_[0];
					yproj = q2_[0];

					xproj /= lam;
					yproj /= lam;

					/* Eventually reswap xproj and yproj like in the original code	*/

					#if OFFSET==0

					xproj -= 1.0;
					yproj -= 1.0;

					#endif

					xproj -= x0_src;
					yproj -= y0_src;

//					xv = floor(xproj);	// floor causes problem ??
//					yv = floor(yproj);

					xv = (long)xproj;
					yv = (long)yproj;

					xdiff = xproj - xv;	/*	 Weight for Linear Interpolation	*/
					ydiff = yproj - yv;

					#if TEST==1

					if ((ix-istart)==0.5*(istop-istart) && (iy-y0+1)==0.5*(y1-y0+1)) {

						printf( "ix=%i iy=%i  isect=%i  index=%i   xproj=%f  yproj=%f  xv=%i  yv=%i\n",
								q1_[0],q2_[0],f1_implicit(ix,iy,z_section), f2_implicit(ix,iy,z_section));

						printf( "ix=%i iy=%i  isect=%i  index=%i   xproj=%f  yproj=%f  xv=%i  yv=%i\n",
								 ix,iy,isect,index,xproj,yproj,xv,yv);

					}

					#endif

					/*	 Calculate linear interpolation error, making sure that we are within the coordinate limits	*/

					if (xv>=0 && xv<nx_src-1 && yv>=0 && yv<ny_src-1) {
						block[index] = block[index] + INTERP2(slice_data_f, nx_src, ny_src, xv, yv, xdiff, ydiff);
					}

					for (m=n_-1;m>=0;m--) {
						q1_[m] += q1_[m+1];
						q2_[m] += q2_[m+1];
					}

					lam += dely_lam;

					index++;

				}

			}

		}

	}

}

