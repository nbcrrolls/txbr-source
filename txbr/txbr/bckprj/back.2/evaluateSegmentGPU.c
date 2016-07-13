#include <math.h>
#include <string.h>
#include <stdio.h>

#include <bckprj.h>

double f1_implicit(double x, double y, double z);
double f2_implicit(double x, double y, double z);

void solve_ivp1(double x0, double y0, double z0, double *ivp1);
void solve_ivp2(double x0, double y0, double z0, double *ivp1);

// EXTERN variables
extern int x_inc;
extern double a1_[MAX_COEFFICIENTS_NUMBER];
extern double q1_[MAX_ORDER];

/*
 * Caculate the segment size for the propagation error to remain less than eps
 */
int calculate_segment_size(int order, int x0, int y0, int z0, int x_inc, int x1, double *coefficients, double eps) {

	int ix_0 = x0-1;
	int ix_1 = x1-1;

	int blocksize_x = (ix_1-ix_0+1)/x_inc;

	memcpy(a1_, coefficients, MAX_COEFFICIENTS_NUMBER*sizeof(double));

	int m, segmentSize=0;

	double value, value_app, diff=0.0;

	solve_ivp1(x0, y0, z0, &q1_[0]);

	while (diff<eps && segmentSize<blocksize_x) {

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

	segmentSize = segmentSize*x_inc;

	printf("Segment Size used for this recursion algorithm: %i\n", segmentSize);

	return segmentSize;

}

