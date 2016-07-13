#include <stdio.h>
#include <math.h>

#define	MAXN	30

float a2_[MAXN][MAXN];
float a3_[MAXN][MAXN][MAXN];
int nn_[MAXN];

int bc[MAXN][MAXN];	/* table for binomial coefficient */
int fac[MAXN];	/* table for factorial */
float alpha[MAXN][MAXN];	/* table for alpha */

float diff2[MAXN][MAXN];

/*****************************************************************************
* Function:     max2
*
* Argument:		n1 an integer
*				n2 an integer
*
* Returns:      The maximum between n1 and n2.
*
* Description:  This routine calculates the maximum of two integers.
*****************************************************************************/
int max2(int n1, int n2) {
	
	int max = n1<n2 ? n2 : n1;
	
	return max;
	
}

/*****************************************************************************
* Function:     max3
*
* Argument:		n1 an integer
*				n2 an integer
*				n3 an integer
*
* Returns:      The maximum between n1, n2 and n3.
*
* Description:  This routine calculates the maximum of three integers.
*****************************************************************************/
int max3(int n1, int n2, int n3) {
	
	int max = n2<n3 ? n3 : n2;
	max = n1<max ? max : n1;
	
	return max;
	
}

/*****************************************************************************
* Function:     binomial_coefficient
*
* Argument:		n	maximum order.
*
* Returns:      int*	a pointer to a table containing binomial coefficient
*						values
*
* Description:  This routine calculates the binomial coefficients.
*****************************************************************************/
int *binomial_coefficients(int n) {
	
	int i,j;
	
	// Set every elements to 0
	for (i=0; i<=MAXN; i++) for (j=0; j<=MAXN; j++) bc[i][j] = 0;

	// Calculate the first elements
	for (i=0; i<=n; i++) bc[i][0] = 1;
	for (j=0; j<=n; j++) bc[j][j] = 1;

	for (i=1; i<=n; i++)
		for (j=1; j<i; j++)
			bc[i][j] = bc[i-1][j-1] + bc[i-1][j];
	
	return &bc[0][0];

}

/*****************************************************************************
* Function:     factorial_coefficients
*
* Argument:		n	maximum order.
*
* Returns:      int*	a pointer to a table containing binomial coefficient
*						values
*
* Description:  This routine calculates factorials up to order n.
*****************************************************************************/
int *factorial_coefficients(int n) {
	
	int i;
	
	// Set every elements to 0
	for (i=0; i<=MAXN; i++) fac[i] = 0;
	
	// Calculate the first n elements
	fac[0] = 1;
	for (i=1; i<=n; i++) fac[i] = i*fac[i-1];
	
	return &fac[0];
	
}

/*****************************************************************************
* Function:     Calculate what we call the alpha coefficient.
*
* Argument:		n	maximum order.
*
* Returns:      double*	a pointer to a table containing alpha coefficient
*						values
*
* Description:  This routine calculates alpha.
*****************************************************************************/
float *alpha_coefficients(int n) {

	float sum;
	int i1,i2,index;
	
	/*
	
	// Direct calculation
	
	binomial_coefficients(n);

	for (i1=0; i1<n+1; i1++) {
		for (i2=0; i2<n+1; i2++) {
			sum = 0.0;
			for (index=1;index<i1+1;index++) {
				sum += pow(-1,index+i2)*pow(index,i2)*bc[i1][index];
			}
			alpha[i1][i2] = sum;
		}
	}
	
	*/
	
	// Recursive calculation. Improves accuracy.
	
	for (i1=1; i1<n+1; i1++) alpha[i1][0] = -1;
	for (i2=0; i2<n+1; i2++) alpha[0][i2] = 0;
	
	for (i2=1; i2<n+1; i2++) {
		for (i1=1; i1<n+1; i1++) {
			alpha[i1][i2] = -i1*alpha[i1][i2-1]+i1*alpha[i1-1][i2-1];
		}
	}
	
	/*
	for (i1=0; i1<n+1; i1++) {
		for (i2=0; i2<n+1; i2++) {
			printf("(%i,%i) %f\n",i1,i2,alpha[i1][i2]);
		}
	}
	*/

	return &alpha[0][0];
	
}

/*****************************************************************************
* Function:     f2
*
* Argument:		x1	first coordinate
*
* Returns:      
*
* Description:  This routine calculates alpha.
*****************************************************************************/
float f2(float x1, float x2) {
	
	float x = 0.0;
	int k1,k2;
	
	for (k1=0; k1<nn_[0]+1; k1++) {
		for (k2=0; k2<nn_[1]+1; k2++) {
			x += a2_[k1][k2]*pow(x1,k1)*pow(x2,k2);
		}
	}
			
	return x;
	
}

/*****************************************************************************
* Function:     f3
*
* Argument:		x1	first coordinate
*
* Returns:      
*
* Description:  This routine calculates alpha.
*****************************************************************************/
float f3(float x1, float x2, float x3) {
	
	float x = 0.0;
	int k1,k2,k3;
	
	for (k1=0; k1<nn_[0]+1; k1++) {
		for (k2=0; k2<nn_[1]+1; k2++) {
			for (k3=0; k2<nn_[2]+1; k3++) {
				x += a3_[k1][k2][k3]*pow(x1,k1)*pow(x2,k2)*pow(x3,x3);
			}
		}
	}
			
	return x;
	
}
	

/*****************************************************************************
* Function:     differences_1D
*
* Argument:		n	dimension of the polynomial coefficient
*				*f_1	pointer to the polynomial function
*				*ivp
*
*
* Description:  This routine calculates alpha.
*****************************************************************************/
void differences_1D(int n, float (*f_1)(float), float *ivp) {
	
	int i,j;
	
	// Set every elements to 0
	for (i=0; i<=MAXN; i++) for (j=0; j<=MAXN; j++) diff2[i][j] = 0;

	// Calculate the first elements
	for (j=0; j<=n; j++) {
		diff2[0][j] = f_1(-j); // printf("%.2e  ",diff2[0][j]);
	}
	ivp[0] = diff2[0][0]; printf("%.2e",ivp[0]);

	printf("\n");

	for (i=1; i<n+1; i++) {
		for (j=0; j<n-i+1; j++) {
			diff2[i][j] = diff2[i-1][j+1] - diff2[i-1][j];  // printf("%.2e  ",diff2[i][j]);
		}
		ivp[i] = diff2[i][0];
		ivp[i] = (i%2==0) ? ivp[i] : -ivp[i]; printf("%.2e",ivp[i]);
		printf("\n");
	}
	printf("\n");
	
}

