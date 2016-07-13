#include <math.h>
#include <string.h>
#include <stdio.h>

#include <bckprj.h>

double f1_implicit(double x, double y, double z);
double f2_implicit(double x, double y, double z);


/*
 * Caculate the Box to use!
 */
void calculate_box( int x_start, int y_start, int z_start, int x_stop, int y_stop, int z_stop,
					int nx, int ny, int *x0_src, int *y0_src, int *nx_src, int *ny_src) {

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

	x[0] = f1_implicit(x_start,y_start,z_start);
	x[1] = f1_implicit(x_start,y_start,z_stop);
	x[2] = f1_implicit(x_start,y_stop,z_start);
	x[3] = f1_implicit(x_start,y_stop,z_stop);
	x[4] = f1_implicit(x_stop,y_start,z_start);
	x[5] = f1_implicit(x_stop,y_start,z_stop);
	x[6] = f1_implicit(x_stop,y_stop,z_start);
	x[7] = f1_implicit(x_stop,y_stop,z_stop);

	y[0] = f2_implicit(x_start,y_start,z_start);
	y[1] = f2_implicit(x_start,y_start,z_stop);
	y[2] = f2_implicit(x_start,y_stop,z_start);
	y[3] = f2_implicit(x_start,y_stop,z_stop);
	y[4] = f2_implicit(x_stop,y_start,z_start);
	y[5] = f2_implicit(x_stop,y_start,z_stop);
	y[6] = f2_implicit(x_stop,y_stop,z_start);
	y[7] = f2_implicit(x_stop,y_stop,z_stop);

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
