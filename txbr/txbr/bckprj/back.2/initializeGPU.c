#include <string.h>

#define MAX_ORDER       20  /*  Maximum order for the projection polynomial approximation*/
#define MAX_COEFFICIENTS_NUMBER 1771  /*        Maximum order for the projection polynomial approximation*/

long order_in_X[MAX_COEFFICIENTS_NUMBER];
long order_in_Y[MAX_COEFFICIENTS_NUMBER];
long order_in_Z[MAX_COEFFICIENTS_NUMBER];

int n_ = 2; /*	 Order of the reconstruction	*/
int number_of_terms;

double a1_[MAX_COEFFICIENTS_NUMBER];
double a2_[MAX_COEFFICIENTS_NUMBER];
double lambda_[4];
double rot_[3][3];
double bottom_plane_[3];

// EXTERN variable
extern int x_inc, y_inc;

void initializeGPU(int n__, int number_of_terms_, int *order_in_X_, int *order_in_Y_, int *order_in_Z_, double *lambda__,
		 	double *a1__, double *a2__, double *bottom_plane__, double *rot__, int x_inc_, int y_inc_) {

	n_ = n__;
	number_of_terms = number_of_terms_;

	x_inc = x_inc_;
	y_inc = y_inc_;

	//	memcpy(&order_in_X[0], &order_in_X_[0], MAX_COEFFICIENTS_NUMBER*sizeof(int));
	//	memcpy(&order_in_Y[0], &order_in_Y_[0], MAX_COEFFICIENTS_NUMBER*sizeof(int));
	//	memcpy(&order_in_Z[0], &order_in_Z_[0], MAX_COEFFICIENTS_NUMBER*sizeof(int));

	int i;

	for (i=0;i<MAX_COEFFICIENTS_NUMBER;i++) { // change int in long
		order_in_X[i] = order_in_X_[i];
		order_in_Y[i] = order_in_Y_[i];
		order_in_Z[i] = order_in_Z_[i];
	}

 	memcpy(lambda_, lambda__, 4*sizeof(double));
	memcpy(a1_, a1__, MAX_COEFFICIENTS_NUMBER*sizeof(double));
	memcpy(a2_, a2__, MAX_COEFFICIENTS_NUMBER*sizeof(double));

	/* Make sure it is 9 elements array */
	memcpy((void*)bottom_plane_, (void *)bottom_plane__, 3 * sizeof(double));
	memcpy((void *)rot_, (void *)rot__, 9*sizeof(double));


}


