
/*****************************************************************************
*
*	Author(s)	:	Sebastien Phan, Raj Singh
*	Emails		:	sph@ncmir.ucsd.edu, raj@ncmir.ucsd.edu
*
*	Status		:	Experimental
*
*	Description	:	This is the CUDA kernel portion for GPU accelerated
*				backprojection code in TxBR. The original code 
*				was taken mostly from backprojection.c 
*
*****************************************************************************/



/* Old Headers .. all data structures used by MRC or IMOD have been unrolled to
   simpler data types.
#include <math.h>
#include "mrcslice.h"
#include "txbrutil.h"
*/

#include <stdio.h>
#include "backprojectionGPU.h"

#define TXBR_MIN(a,b) ((a) < (b) ? (a) : (b))
#define TXBR_MAX(a,b) ((a) > (b) ? (a) : (b))

#define INTERPOLE2(data, xsize, ysize, mean, mode, csize,i,j,wi,wj) (mrc_slice_getmagnitude_gpu(data, xsize, ysize, mean, mode, csize,i,j)*(1.0-wi)*(1-wj) + \ 
								mrc_slice_getmagnitude_gpu(data, xsize, ysize, mean, mode, csize,i+1,j)*wi*(1-wj) + \
								mrc_slice_getmagnitude_gpu(data, xsize, ysize, mean, mode, csize,i,j+1)*(1-wi)*wj + \
								mrc_slice_getmagnitude_gpu(data, xsize, ysize, mean, mode, csize,i+1,j+1)*wi*wj)

#define MAX_ORDER       		20  	// Maximum order for the projection polynomial approximation
#define MAX_COEFFICIENTS_NUMBER 	100 // 1771  	// Maximum order for the projection polynomial approximation


/*
// The following variables are temporary storage for their faster shared memory counterparts
CUDABP_DEV_VAR_DECL int order_in_X_temp[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL int order_in_Y_temp[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL int order_in_Z_temp[MAX_COEFFICIENTS_NUMBER];

CUDABP_DEV_VAR_DECL double bottom_plane__temp[3];
CUDABP_DEV_VAR_DECL double lambda__temp[4];
CUDABP_DEV_VAR_DECL double a1__temp[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL double a2__temp[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL double rot__temp[3][3];
CUDABP_DEV_VAR_DECL int n__temp; //	 Order of the reconstruction
CUDABP_DEV_VAR_DECL int x_inc_temp;
CUDABP_DEV_VAR_DECL int number_of_terms_temp;
*/

// Global variables

CUDABP_DEV_VAR_DECL int order_in_X[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL int order_in_Y[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL int order_in_Z[MAX_COEFFICIENTS_NUMBER];

CUDABP_DEV_VAR_DECL int n_; /*	 Order of the reconstruction	*/
CUDABP_DEV_VAR_DECL int x_inc;
CUDABP_DEV_VAR_DECL int number_of_terms;

CUDABP_DEV_VAR_DECL double bottom_plane_[3];
CUDABP_DEV_VAR_DECL double lambda_[4];
CUDABP_DEV_VAR_DECL double a1_[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL double a2_[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL double rot_[3][3];


CUDABP_DEV_VAR_DECL int segment_size_d = 0;

void *slice_data_d_ptr = NULL;
void *block_d_ptr = NULL;

cudaError_t cErr;



// other forward decls
CUDABP_DEV_FUNC_DECL float mrc_slice_getmagnitude_gpu(void *data, int xsize, int ysize, 
		float mean, int mode, int csize, 
		int x, int y);
CUDABP_DEV_FUNC_DECL int sliceGetVal_gpu(void *data, int xsize, int ysize, 
		float mean, int mode, 
		int x, int y, float * val);

CUDABP_DEV_FUNC_DECL double lambda(double x, double y, double z);
CUDABP_DEV_FUNC_DECL double f1(double x, double y, double z);
CUDABP_DEV_FUNC_DECL double f2(double x, double y, double z);

//CUDABP_DEV_FUNC_DECL double* Rot(double *XYZ);
//CUDABP_DEV_FUNC_DECL double* Rot_inverse(double *XYZ);
//CUDABP_DEV_FUNC_DECL double* XYZ_t(double *XYZ);
CUDABP_DEV_FUNC_DECL void Rot(double *XYZ);
CUDABP_DEV_FUNC_DECL void Rot_inverse(double *XYZ);
CUDABP_DEV_FUNC_DECL void XYZ_t(double *XYZ);

CUDABP_DEV_FUNC_DECL double lambda_implicit(double x, double y, double z);
CUDABP_DEV_FUNC_DECL double f1_implicit(double x, double y, double z);
CUDABP_DEV_FUNC_DECL double f2_implicit(double x, double y, double z);
CUDABP_DEV_FUNC_DECL void solve_ivp1(double x0, double y0, double z0, double *ivp);
CUDABP_DEV_FUNC_DECL void solve_ivp2(double x0, double y0, double z0, double *ivp);


// The Entry functions which are defined as __global__
CUDABP_ENTRY_FUNC_DECL void initializeGPU_d(int n__, int number_of_terms_, int x_inc_);

CUDABP_ENTRY_FUNC_DECL void calculate_segment_size_d(int order, int x0, int y0, int z0, int x1, 
		double eps);

CUDABP_ENTRY_FUNC_DECL void evaluateBlock_d(int blocksize, int x0, int y0, int z0, 
		int x_inc, int y_inc,int z_inc, 
		int nx, int ny, int nx_eff, int ny_eff, int nx_lim, int ny_lim, 
		int nsegment, int segment_size, void * slice_data, int slice_xsize, 
		int slice_ysize, float slice_mean, int slice_mode, int slice_csize,
		float *block);


// This is the host version of the corresponding initializeGPU_d() function invoked on
// the GPU which is called with unpacked data
extern "C" 
int initializeGPU(int n__, int number_of_terms_, 
			int *order_in_X_, int *order_in_Y_, int *order_in_Z_, 
			double *lambda__, double *a1__, double *a2__, 
			double *bottom_plane__, double *rot__, int x_inc_, 
			int _argc, char **_argv)
{
	static int hwInitFlag = 0;
	dim3 grids(1, 1);
    	dim3 block(1, 1);
	void *symPtr;
	int dev;
	
	if(!hwInitFlag)
	{
		// Init the GPUs
		printf("\ninitializeGPU(): Initializing GPUs ... \n");
		CUT_DEVICE_INIT(_argc, _argv);
		
		// since the command line parsing for txbr is incompatible with 
		// the cutil command parsing, we will set the device to be used 
		// explicitly here. The id for the device should come from a config
		// file or the command line
		
		dev = CUDABP_DEV_ID;
		printf("\n\n====================================================");
		printf("\ninitializeGPU(): Resetting device usage for device %d : ", dev);
		cudaDeviceProp deviceProp;                                               
		CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceProperties(&deviceProp, dev));       
		if (deviceProp.major < 1) {                                              
			fprintf(stderr, "\ninitializeGPU()(): Error: device does not support CUDA.\n");
			return -1;
		}
		
		printf("%s\n", deviceProp.name);
		printf("\nShared memory per block = %d", deviceProp.sharedMemPerBlock);
		printf("\nWarp Size = %d", deviceProp.warpSize);
		printf("\nMax threads per block = %d", deviceProp.maxThreadsPerBlock);
		printf("\nNumber of multiprocessors = %d", deviceProp.multiProcessorCount);
		
		printf("\n====================================================\n");
		
		CUDA_SAFE_CALL(cudaSetDevice(dev));
		
		hwInitFlag = 1;
	}
	
	// since we cannot pass arrays to the device through the kernel invocation
	// we will do memcpys to  symbols on the device
	cudaGetSymbolAddress(&symPtr, "order_in_X");
	cErr = cudaMemcpy(symPtr, (void *)order_in_X_, MAX_COEFFICIENTS_NUMBER * sizeof(int), 
		cudaMemcpyHostToDevice);
	
	cudaGetSymbolAddress(&symPtr, "order_in_Y");
	cErr = cudaMemcpy(symPtr, (void *)order_in_Y_, MAX_COEFFICIENTS_NUMBER * sizeof(int), 
		cudaMemcpyHostToDevice);
	
	cudaGetSymbolAddress(&symPtr, "order_in_Z");
	cErr = cudaMemcpy(symPtr, (void *)order_in_Z_, MAX_COEFFICIENTS_NUMBER * sizeof(int), 
		cudaMemcpyHostToDevice);

	// check the error returned occasionally
	if(cErr != cudaSuccess)
	{
		printf("\ninitializeGPU(): %s", cudaGetErrorString(cErr));
		return -1;
	}
	
	cudaGetSymbolAddress(&symPtr, "lambda_");
	cErr = cudaMemcpy(symPtr, (void *)lambda__, 4 * sizeof(double), 
		cudaMemcpyHostToDevice);
	
	cudaGetSymbolAddress(&symPtr, "a1_");
	cErr = cudaMemcpy(symPtr, (void *)a1__, MAX_COEFFICIENTS_NUMBER * sizeof(double), 
		cudaMemcpyHostToDevice);
	
	cudaGetSymbolAddress(&symPtr, "a2_");
	cErr = cudaMemcpy(symPtr, (void *)a2__, MAX_COEFFICIENTS_NUMBER * sizeof(double), 
		cudaMemcpyHostToDevice);

	cudaGetSymbolAddress(&symPtr, "bottom_plane_");
	cErr = cudaMemcpy(symPtr, (void *)bottom_plane__, 3 * sizeof(double), 
		cudaMemcpyHostToDevice);

	cudaGetSymbolAddress(&symPtr, "rot_");
	cErr = cudaMemcpy(symPtr, (void *)rot__, 9 * sizeof(double), 
		cudaMemcpyHostToDevice);
	
	// check the error returned occasionally
	if(cErr != cudaSuccess)
	{
		printf("\ninitializeGPU(): %s", cudaGetErrorString(cErr));
		return -1;
	}
	
	// initialize simple variables on the device
	initializeGPU_d<<< grids, block >>>(n__, number_of_terms_, x_inc_);
	
	
	// everything seems ok
	return 1;
}


CUDABP_ENTRY_FUNC_DECL void initializeGPU_d(int n__, int number_of_terms_, int x_inc_) 
{
	// nothing much happening here except for simple local copies
	n_ = n__;
	x_inc = x_inc_;
	number_of_terms = number_of_terms_;

/*	
	// copy the slow global memory variables to faster shared memory
	int i, j;
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) order_in_X[i] = order_in_X_temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) order_in_Y[i] = order_in_Y_temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) order_in_Z[i] = order_in_Z_temp[i];

	for(i=0; i< 3; i++) bottom_plane_[i] = bottom_plane__temp[i];
	for(i=0; i< 4; i++) lambda_[i] = lambda__temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) a1_[i] = a1__temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) a2_[i] = a2__temp[i];
	for(i = 0; i<3; i++){
		for(j = 0;  j<3 ; j++)	rot_[i][j] = rot__temp[i][j];
	}
	n_ = n__temp;
	x_inc = x_inc_temp;
	number_of_terms = number_of_terms_temp;
*/	
	
	// I guess you need to do this to make sure the modified shared mem variables
	// are visible to all threads.
	__syncthreads();
	
	return;
}


//
// Caculate the segment size for the propagation error to remain less than eps
//
// This is the non kernel version of the corresponding calculate_segment_size_d() 
// function which is called in turn with unpacked data
 extern "C"
 int calculate_segment_size(int order, int x0, int y0, int z0, int x1, 
	 double *coefficients, double eps) {

 	void *symPtr;
	int l_segment_size = 0;
	dim3 grids(1, 1);
    	dim3 block(1, 1);
 
	// Need to copy the coefficient array to a1_
	cudaGetSymbolAddress(&symPtr, "a1_");
	cErr = cudaMemcpy(symPtr, (void *)coefficients, MAX_COEFFICIENTS_NUMBER * sizeof(double), 
		cudaMemcpyHostToDevice);
	
	// invoke the kernel on the GPU to calculate the segment size but just on
	// one core
	calculate_segment_size_d <<< grids, block >>> (order, x0, y0, z0, x1, eps);
	
	// get the symbol pointer to the result on the GPU and copy the result back to host
	cudaGetSymbolAddress((void **)&symPtr, "segment_size_d");
	cErr = cudaMemcpy((void *)&l_segment_size, symPtr, sizeof(int), cudaMemcpyDeviceToHost);
	
 	return l_segment_size;
}

 
CUDABP_ENTRY_FUNC_DECL void calculate_segment_size_d(int order, int x0, int y0, int z0, 
	int x1, double eps) {

	double q1_[MAX_ORDER];
	int ix_0 = x0-1;
	int ix_1 = x1-1;

	int nx_eff = (ix_1-ix_0)/x_inc + 1;
	int nx_lim = ix_0 + x_inc*nx_eff;	

	//memcpy(a1_, coefficients_, MAX_COEFFICIENTS_NUMBER*sizeof(double));

	int m, segmentSize=0;
	double value, value_app, diff=0.0;
	
// DEBUG
/*
	// copy the slow global memory variables to faster shared memory
	int i, j;
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) order_in_X[i] = order_in_X_temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) order_in_Y[i] = order_in_Y_temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) order_in_Z[i] = order_in_Z_temp[i];

	for(i=0; i< 3; i++) bottom_plane_[i] = bottom_plane__temp[i];
	for(i=0; i< 4; i++) lambda_[i] = lambda__temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) a1_[i] = a1__temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) a2_[i] = a2__temp[i];
	for(i = 0; i<3; i++){
		for(j = 0;  j<3 ; j++)	rot_[i][j] = rot__temp[i][j];
	}
	n_ = n__temp;
	x_inc = x_inc_temp;
	number_of_terms = number_of_terms_temp;
*/
// DEBUG

	solve_ivp1(x0, y0, z0, &q1_[0]);
	
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
	
	//printf("Segment Size used for this recursion algorithm: %i\n", segmentSize);
	
	segment_size_d = segmentSize;
	
	return;
}


/*
 *
 */
CUDABP_DEV_FUNC_DECL double lambda(double x, double y, double z) {
	int i;
	double value = 0.0;
	
	for (i=0; i<4; i++) {
		value += lambda_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
	}
	
	return value;
}


CUDABP_DEV_FUNC_DECL double f1(double x, double y, double z) {
	int i;
	double value = 0.0;
	
	for (i=0; i<number_of_terms; i++) {
		value += a1_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
	}
	
	return value;	
}


CUDABP_DEV_FUNC_DECL double f2(double x, double y, double z) {
	
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
//CUDABP_DEV_FUNC_DECL double* Rot(double *XYZ) {
CUDABP_DEV_FUNC_DECL void Rot(double *XYZ) {
	
	double x = XYZ[0];
	double y = XYZ[1];
	double z = XYZ[2];
	
	XYZ[0] = rot_[0][0]*x + rot_[1][0]*y + rot_[2][0]*z;
	XYZ[1] = rot_[0][1]*x + rot_[1][1]*y + rot_[2][1]*z;
	XYZ[2] = rot_[0][2]*x + rot_[1][2]*y + rot_[2][2]*z;	
	
	//return XYZ;	
}

/*
 * Inverse rotation function.
 */
//CUDABP_DEV_FUNC_DECL double* Rot_inverse(double *XYZ) {
CUDABP_DEV_FUNC_DECL void Rot_inverse(double *XYZ) {
	
	/* The matrix is transposed - equivalent to theta in -theta	*/
	
	double x = XYZ[0];
	double y = XYZ[1];
	double z = XYZ[2];
	
	XYZ[0] = rot_[0][0]*x + rot_[0][1]*y + rot_[0][2]*z;
	XYZ[1] = rot_[1][0]*x + rot_[1][1]*y + rot_[1][2]*z;
	XYZ[2] = rot_[2][0]*x + rot_[2][1]*y + rot_[2][2]*z;
	
}

/*
 * This function swap the x, y and follows the bottom plane in z
 */
// CUDABP_DEV_FUNC_DECL double* XYZ_t(double *XYZ) {
CUDABP_DEV_FUNC_DECL void XYZ_t(double *XYZ) {
	 
	double x = XYZ[0];
	double y = XYZ[1];
	double z = XYZ[2];
	
	x = x_inc*x;
	
	#if SWAP_x_y==1
	
	tmp = x;
	x = y;
	y = tmp;
	
	#endif
	
	XYZ[0] = x;
	XYZ[1] = y;
	XYZ[2] = z;
	
	#if FULL_ROTATION==1
	
	XYZ = Rot_inverse(XYZ);	
	
	/*
	XYZ[0] = XYZ_rot[0];
	XYZ_[1] = XYZ_rot[1];
	XYZ_[2] = XYZ_rot[2];
	*/
	
	#else
	
	/*	Original pseudo rotation case	*/
	XYZ[0] = x;
	XYZ[1] = y;
	XYZ[2] = bottom_plane_[1]*x + bottom_plane_[2]*y + z; /* Planes are defined in the wrong reference, reversed.	*/
	
	#endif
	
	//return XYZ;
	
}

/*
 *
 */
CUDABP_DEV_FUNC_DECL double lambda_implicit(double x, double y, double z) {
	
	double XYZ[3];

/* OC	
	XYZ[0] = x;
	XYZ[1] = y;
	XYZ[2] = z;
	
	XYZ_t(XYZ);
*/
	XYZ[0] = x* x_inc;
	XYZ[1] = y;
	XYZ[2] = bottom_plane_[1]*x + bottom_plane_[2]*y + z;	

//OC	return lambda(XYZ[0], XYZ[1], XYZ[2]);
	int i;
	double value = 0.0;
	
	for (i=0; i<4; i++) {
		value += lambda_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
	}
	return value;
	
}

/*
 *
 */
CUDABP_DEV_FUNC_DECL double f1_implicit(double x, double y, double z) {
		
	double XYZ[3];
/* OC	
	XYZ[0] = x;
	XYZ[1] = y;
	XYZ[2] = z;
	
	XYZ_t(XYZ);
*/	
	XYZ[0] = x* x_inc;
	XYZ[1] = y;
	XYZ[2] = bottom_plane_[1]*x + bottom_plane_[2]*y + z;
	
//OC	return f1(XYZ[0], XYZ[1], XYZ[2]);

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
CUDABP_DEV_FUNC_DECL double f2_implicit(double x, double y, double z) {
	
	double XYZ[3];
	
/* OC	
	XYZ[0] = x;
	XYZ[1] = y;
	XYZ[2] = z;
	
	XYZ_t(XYZ);
*/

	XYZ[0] = x* x_inc;
	XYZ[1] = y;
	XYZ[2] = bottom_plane_[1]*x + bottom_plane_[2]*y + z;

//OC	return f2(XYZ[0], XYZ[1], XYZ[2]);

	int i;
	double value = 0.0;
	
	for (i=0; i<number_of_terms; i++) {
		value += a2_[i]*pow(x,order_in_X[i])*pow(y,order_in_Y[i])*pow(z,order_in_Z[i]);
	}
	
	return value;
	
}


/*
 * Solve the initial value problem for the quick polynomial evaluation scheme
 */
CUDABP_DEV_FUNC_DECL void solve_ivp1(double x0, double y0, double z0, double *ivp) {

	int i,j;
	double diff2_[MAX_ORDER][MAX_ORDER];

	 //	 Initialize the array diff2_ to zeros
	for (i=0; i<MAX_ORDER; i++) for (j=0; j<MAX_ORDER; j++) diff2_[i][j] = 0;
	
	 //	 Calculate the first elements
	for (j=0; j<n_+1; j++) {
		diff2_[0][j] = f1_implicit(x0-j,y0,z0);
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
CUDABP_DEV_FUNC_DECL void solve_ivp2(double x0, double y0, double z0, double *ivp) {


	int i,j;
	double diff2_[MAX_ORDER][MAX_ORDER];
	
	 //	 Initialize the array diff2_ to zeros
	for (i=0; i<MAX_ORDER; i++) for (j=0; j<MAX_ORDER; j++) diff2_[i][j] = 0;
	
	 //	 Calculate the first elements
	for (j=0; j<n_+1; j++) {
		diff2_[0][j] = f2_implicit(x0-j,y0,z0);
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


//
// This is the host interface to the evaluateBlock_d() kernel which is used for the bulk of 
// computing. We can hide the actual kernel invocation here. Not sure yet but for 
// performance tweaks, I might have to so multiple invokations with smaller tasks
extern "C"
int evaluateBlock(int blocksize, int x0, int y0, int z0, int x_inc, int y_inc,int z_inc, 
		int nx, int ny, int nx_eff, int ny_eff, int nx_lim, int ny_lim, 
		int nsegment, int segment_size, void *slice_data, int slice_xsize, 
		int slice_ysize, float slice_mean, int slice_mode, int slice_csize, 
		int slice_dsize, float *block) {
	
	int memSize;

//	dim3 gpuGrids(CUDABP_GRID_XDIM, CUDABP_GRID_YDIM);
//	dim3 gpuBlock(CUDABP_BLOCK_XDIM, CUDABP_BLOCK_YDIM);

	// If no memory is allocated for in or out buffers, do that here
	if(slice_data_d_ptr == NULL || block_d_ptr == NULL)
	{
		
		// Allocate memory for slice data
		memSize = slice_xsize * slice_ysize * slice_csize * slice_dsize;
		if(CUDA_SAFE_CALL((cErr = cudaMalloc((void **)&slice_data_d_ptr, memSize))
			== cudaErrorMemoryAllocation))
		{
		    printf("evaluateBlock(): Error allocating %d bytes for slice on GPU", 
						memSize);
		    return -1;
		}
		
		// Allocate memory for output block
		memSize = nx_eff*ny_eff*blocksize*sizeof(float);
		if(CUDA_SAFE_CALL((cErr = cudaMalloc((void **)&block_d_ptr, memSize))
			== cudaErrorMemoryAllocation))
		{
		    printf("evaluateBlock(): Error allocating %d bytes for block on GPU", 
						memSize);
		    return -1;
		}
		
		// set the block values to zero
		cudaMemset(block_d_ptr, 0, memSize);
	}
	
	
	// Copy the slice data to the GPU memory
	memSize = slice_xsize * slice_ysize * slice_csize * slice_dsize;
	cErr = cudaMemcpy(slice_data_d_ptr, (void *)slice_data, memSize, 
		cudaMemcpyHostToDevice);
	
	if(cErr != cudaSuccess)
	{
		printf("\nevaluateBlock(): %s", cudaGetErrorString(cErr));
		return -1;
	}
	
	evaluateBlock_d <<< CUDABP_GRID_NUM, CUDABP_BLOCK_NUM >>>(blocksize, x0, y0, z0, x_inc, y_inc, z_inc, 
		nx, ny, nx_eff, ny_eff, nx_lim, ny_lim, 
		nsegment, segment_size, slice_data_d_ptr, slice_xsize, 
		slice_ysize, slice_mean, slice_mode, slice_csize,
		(float *)block_d_ptr);
	
	CUT_CHECK_ERROR("evaluateBlock(): Kernel execution failed");

	// copy back the computed block array. Assuming things went well
	memSize = nx_eff*ny_eff*blocksize*sizeof(float);
	cErr = cudaMemcpy(block, block_d_ptr, memSize, 
			cudaMemcpyDeviceToHost);
	
	if(cErr != cudaSuccess)
	{
		printf("\nevaluateBlock(): %s", cudaGetErrorString(cErr));
		return -1;
	}
	

	return 1;
}


CUDABP_ENTRY_FUNC_DECL void evaluateBlock_d(int blocksize, int x0, int y0, int z0, int x_inc, int y_inc,int z_inc, 
		int nx, int ny, int nx_eff, int ny_eff, int nx_lim, int ny_lim, 
		int nsegment, int segment_size, void *slice_data, int slice_xsize, 
		int slice_ysize, float slice_mean, int slice_mode, int slice_csize,
		float *block) 
{

	int isect;
	double z_section;
	int iy, ix, isegment, istart, istop;
	double XYZ[3];
	double lambda, dely_lambda;

	double xproj, yproj;

	int xv, yv, m;
	double xdiff, ydiff;
	
	double q1_[MAX_ORDER];
	double q2_[MAX_ORDER];

	// first thing find out the index of the core
        // LWC -> unsigned int thIndex = threadIdx.x + threadIdx.y * blockDim.y;
	// LWC -> unsigned int totalNoOfGpus = blockDim.x * blockDim.y;
	//unsigned int thIndex = threadIdx.x + blockIdx.x * gridDim.x;
	unsigned int thIndex = threadIdx.x + blockIdx.x * blockDim.x;
        unsigned int totalNoOfThreads = gridDim.x * blockDim.x;

//DEBUG
/*
if(thIndex == 0)
{
	// copy the slow global memory variables to faster shared memory
	int i,j;
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) order_in_X[i] = order_in_X_temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) order_in_Y[i] = order_in_Y_temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) order_in_Z[i] = order_in_Z_temp[i];

	for(i=0; i< 3; i++) bottom_plane_[i] = bottom_plane__temp[i];
	for(i=0; i< 4; i++) lambda_[i] = lambda__temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) a1_[i] = a1__temp[i];
	for(i = 0; i<MAX_COEFFICIENTS_NUMBER; i++) a2_[i] = a2__temp[i];
	for(i = 0; i<3; i++){
		for(j = 0;  j<3 ; j++)	rot_[i][j] = rot__temp[i][j];
	}
	n_ = n__temp;
	x_inc = x_inc_temp;
	number_of_terms = number_of_terms_temp;
}
__syncthreads();
*/
	
	for (isect=0; isect<blocksize; isect++) {	/*	 Step through the output block	*/
			 
		z_section = z0 + z_inc*isect;	/*	 z-section coordinate of the output slice */

		/*	 Loop along y image coordinate */
		//LWC -> for (iy=y0-1 + threadIdx.y; iy<ny_lim; iy += y_inc + blockDim.y - 1) {
		for (iy=y0-1 + thIndex; iy<ny_lim; iy += y_inc + totalNoOfThreads - 1) {
		// OC 1 -> for (iy=y0-1; iy<ny_lim; iy=iy+y_inc) {
			
			// Now increment along x direction. Iterative evaluation scheme for each segment
			for (isegment=0; isegment<nsegment; isegment++) {
			
				istart = x0 + isegment*segment_size - 1;
				istop = TXBR_MIN( istart+segment_size, nx_lim );				

				XYZ[0] = istart + 1;
				XYZ[1] = iy + 1;
				XYZ[2] = z_section;
	
				lambda = lambda_implicit(XYZ[0],XYZ[1],XYZ[2]);
				dely_lambda = lambda_implicit(XYZ[0]+x_inc,XYZ[1],XYZ[2])-lambda;
			
				xproj = f1_implicit(XYZ[0],XYZ[1],XYZ[2]);
				yproj = f2_implicit(XYZ[0],XYZ[1],XYZ[2]);
			
//				solve_ivp1(XYZ[0],XYZ[1],XYZ[2],q1_);
//				solve_ivp2(XYZ[0],XYZ[1],XYZ[2],q2_);

//  **************************************   DEBUG CODE
			
				// solve_ivp1 unrolled
				int r,s;
				double diff2_[MAX_ORDER][MAX_ORDER];
			
				 //	 Initialize the array diff2_ to zeros
				for (r=0; r<MAX_ORDER; r++) for (s=0; s<MAX_ORDER; s++) diff2_[r][s] = 0;
				
				 //	 Calculate the first elements
				for (s=0; s<n_+1; s++) {
					diff2_[0][s] = f1_implicit(XYZ[0]-s,XYZ[1],XYZ[2]);
				}
				q1_[0] = diff2_[0][0];
				
				for (r=1; r<n_+1; r++) {
					for (s=0; s<n_-r+1; s++) {
						diff2_[r][s] = diff2_[r-1][s+1] - diff2_[r-1][s];
					}
					q1_[r] = diff2_[r][0];
					q1_[r] = (r%2==0) ? q1_[r] : - q1_[r];
				}
				
				
				// solve_ivp2 unrolled
				//	 Initialize the array diff2_ to zeros
				for (r=0; r<MAX_ORDER; r++) for (s=0; s<MAX_ORDER; s++) diff2_[r][s] = 0;
				
				 //	 Calculate the first elements
				for (s=0; s<n_+1; s++) {
					diff2_[0][s] = f2_implicit(XYZ[0]-s,XYZ[1],XYZ[2]);
				}
				q2_[0] = diff2_[0][0];
				
				for (r=1; r<n_+1; r++) {
					for (s=0; s<n_-r+1; s++) {
						diff2_[r][s] = diff2_[r-1][s+1] - diff2_[r-1][s];
					}
					q2_[r] = diff2_[r][0];
					q2_[r] = (r%2==0) ? q2_[r] : - q2_[r];
				}

//  **************************************   DEBUG CODE

				//for(ix = istart + thIndex; ix < istop; ix += x_inc + (totalNoOfThreads - 1)) {
				for (ix=istart; ix<istop; ix=ix+x_inc) {
				// LWC->for(ix = istart + threadIdx.y; ix < istop; ix += x_inc + (blockDim.x - 1)) {
						
					xproj = q1_[0]; 
					yproj = q2_[0];
					
					xproj /= lambda;
					yproj /= lambda;
										
					// Eventually reswap xproj and yproj like in the original code
		
					xproj-=1.0;
					yproj-=1.0;
					
					xv = floor(xproj);
					yv = floor(yproj);
                    
					xdiff = xproj - xv;	//	 Weight for Linear Interpolation
					ydiff = yproj - yv;
					                
					//	 Calculate linear interpolation error, making sure that we are within the coordinate limits
					 
					if (xv>=0 && xv<nx-1 && yv>=0 && yv<ny-1) {						
						block[isect*nx_eff*ny_eff + iy*nx_eff+ix] += INTERPOLE2(slice_data, slice_xsize, slice_ysize, 
												slice_mean, slice_mode, slice_csize, 
												xv, yv, xdiff, ydiff);
						
					}
					
					//block[isect*nx*ny + iy*nx+ix] = mrc_slice_getmagnitude_gpu(slice_data, slice_xsize, slice_ysize, slice_mean, slice_mode, slice_csize,ix,iy);
						
					for (m=n_-1;m>=0;m--) {
						q1_[m] += q1_[m+1];
						q2_[m] += q2_[m+1];
					}
					
					lambda += dely_lambda;
						
				} // for (ix ..				
			} // for (isegment ..			
		} // for(iy ..		
	} // for (isect ..
	
	// sync threads here before returning
	__syncthreads();
	
	return;

}



//
//      Following are some dependecies from IMOD that have been rolled into the code here to provide
//      independence from libraries


CUDABP_DEV_FUNC_DECL
float mrc_slice_getmagnitude_gpu(void *data, int xsize, int ysize, float mean, int mode, int csize, 
	int x, int y)
{

  float val[3];
  float m = 0.0;

  sliceGetVal_gpu(data, xsize, ysize, mean, mode, x, y, val);

  if (csize == 1)
    return(val[0]);

  if (csize == 2){
    m = (val[0] * val[0]) + (val[1] * val[1]);
    return((float)sqrt(m));
  }

  //m = val[0] * 0.3f;
  //m += val[1] * 0.59f;
  //m += val[2] * 0.11f;
  
  m = val[0] * 0.3f + val[1] * 0.59f + val[2] * 0.11f ;
    
  return(m);
}


CUDABP_DEV_FUNC_DECL
int sliceGetVal_gpu(void *data, int xsize, int ysize, float mean, int mode, 
		int x, int y, float * val)
{
  int index = x + (y * xsize);

  if ( (x < 0) || (y < 0) || (x >= xsize) || (y >= ysize)){
    val[0] = mean;
    return(-1);
  }

  val[1] = 0;
  switch (mode)
    {
    case 0: // MRC_MODE_BYTE
      //val[0] = (float)s->data.b[index];
      val[0] = (float)(((unsigned char *)data)[index]);
      break;
    case 1: //MRC_MODE_SHORT
      //val[0] = (float)s->data.s[index];
      val[0] = (float)(((short int *)data)[index]);
      break;
    case 2: // MRC_MODE_FLOAT
      //val[0] = s->data.f[index];
      val[0] = ((float *)data)[index];
      break;
    case 3: // MRC_MODE_COMPLEX_SHORT
      index *= 2;
      //val[0] = (float)s->data.s[index];
      //val[1] = (float)s->data.s[index + 1];
      val[0] = (float)(((short int *)data)[index]);
      val[1] = (float)(((short int *)data)[index + 1]);
      break;
    case 4:// MRC_MODE_COMPLEX_FLOAT
      index *= 2;
      //val[0] = s->data.f[index];
      //val[1] = s->data.f[index + 1];
      val[0] = ((float *)data)[index];
      val[1] = ((float *)data)[index + 1];
      break;
    case 16: // MRC_MODE_RGB
      index *= 3;
      //val[0] = (float)s->data.b[index++];
      //val[1] = (float)s->data.b[index++];
      //val[2] = (float)s->data.b[index];
      val[0] = (float)(((unsigned char *)data)[index++]);
      val[1] = (float)(((unsigned char *)data)[index++]);
      val[2] = (float)(((unsigned char *)data)[index]);
      
      break;
    case 99: // SLICE_MODE_MAX
      index *= 3;
      //val[0] = s->data.f[index++];
      //val[1] = s->data.f[index++];
      //val[2] = s->data.f[index];
      val[0] = ((float *)data)[index++];
      val[1] = ((float *)data)[index++];
      val[2] = ((float *)data)[index];
      break;
    default:
      //b3dError(stderr, "sliceGetVal: unknown mode.\n");
      return(-1);
    }
  return(0);
}



extern "C"
void releaseGPU()
{
	// clean up memory allocated	
	if(slice_data_d_ptr)
	{
		CUDA_SAFE_CALL(cudaFree(slice_data_d_ptr));
		slice_data_d_ptr = NULL;
	}
	if(block_d_ptr)
	{
		CUDA_SAFE_CALL(cudaFree(block_d_ptr));
		block_d_ptr = NULL;
	}
}






