
/*****************************************************************************
*
*	Author(s)	:	Sebastien Phan, Raj Singh
*	Emails		:	sph@ncmir.ucsd.edu, raj@ncmir.ucsd.edu
*
*	Status		:	Experimental
*
*	Description	:	This is the CUDA kernel portion for GPU accelerated
*				backprojection code in TxBR. 
*
*****************************************************************************/


#include "bckprjGPU.h"
#include "bckprj.h"
#include <stdio.h>
#include "mycutil.h"

CUDABP_DEV_VAR_DECL int order_in_X[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL int order_in_Y[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL int order_in_Z[MAX_COEFFICIENTS_NUMBER];

CUDABP_DEV_VAR_DECL int n_;
CUDABP_DEV_VAR_DECL int number_of_terms;
CUDABP_DEV_VAR_DECL int x_inc, y_inc;

CUDABP_DEV_VAR_DECL double bottom_plane_[3];
CUDABP_DEV_VAR_DECL double lambda_[4];
CUDABP_DEV_VAR_DECL double a1_[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL double a2_[MAX_COEFFICIENTS_NUMBER];
CUDABP_DEV_VAR_DECL double rot_[3][3];

/* NOT SURE ABOUT WHY THESE ARE GLOBAL 
double q1_[MAX_ORDER];
double q2_[MAX_ORDER];
double diff2_[MAX_ORDER][MAX_ORDER];
double XYZ_[3], XYZ_rot[3];
*/

CUDABP_DEV_VAR_DECL int segment_size_d = 0;

float *slice_data_d_ptr = NULL;
float *block_d_ptr = NULL;

cudaError_t cErr;

// Forward declarations

CUDABP_DEV_FUNC_DECL double lambda(double x, double y, double z);
CUDABP_DEV_FUNC_DECL double f1(double x, double y, double z);
CUDABP_DEV_FUNC_DECL double f2(double x, double y, double z);

//double* Rot(double x, double y, double z);
//double* Rot_inverse(double x, double y, double z);
//double* XYZ_t(double x, double y, double z);

CUDABP_DEV_FUNC_DECL void Rot(double *XYZ);
CUDABP_DEV_FUNC_DECL void Rot_inverse(double *XYZ);
CUDABP_DEV_FUNC_DECL void XYZ_t(double *XYZ);

CUDABP_DEV_FUNC_DECL double lambda_implicit(double x, double y, double z);
CUDABP_DEV_FUNC_DECL double f1_implicit(double x, double y, double z);
CUDABP_DEV_FUNC_DECL double f2_implicit(double x, double y, double z);
CUDABP_DEV_FUNC_DECL void solve_ivp1(double x0, double y0, double z0, double *ivp1);
CUDABP_DEV_FUNC_DECL void solve_ivp2(double x0, double y0, double z0, double *ivp2);

CUDABP_ENTRY_FUNC_DECL void evaluateBlock_d( float *block, int x0, int y0, int z0, int x_inc, int y_inc, int z_inc, 
		int x1, int y1, int z1, float* slice_data_f, int x0_src, int y0_src, 
		int nx_src, int ny_src, int nseg, int seg_size);
CUDABP_ENTRY_FUNC_DECL void initializeGPU_d(int n__, int number_of_terms_, int x_inc_);
CUDABP_ENTRY_FUNC_DECL void calculate_segment_size_d(int order, int x0, int y0, int z0, int x_inc,
		int x1, double eps);


/*
 *
 */
CUDABP_DEV_FUNC_DECL double lambda(double x, double y, double z) {

	int i;
	double value = 0.0;

	value = lambda_[0];
	for (i=1; i<4; i++) {
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


CUDABP_DEV_FUNC_DECL  double f2(double x, double y, double z) {

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
//double* Rot(double x, double y, double z) {
CUDABP_DEV_FUNC_DECL void Rot(double *XYZ) {

	double x = XYZ[0];
	double y = XYZ[1];
	double z = XYZ[2];
	
	XYZ[0] = rot_[0][0]*x + rot_[1][0]*y + rot_[2][0]*z;
	XYZ[1] = rot_[0][1]*x + rot_[1][1]*y + rot_[2][1]*z;
	XYZ[2] = rot_[0][2]*x + rot_[1][2]*y + rot_[2][2]*z;

}

/*
 * Inverse rotation function.
 */
//double* Rot_inverse(double x, double y, double z) {
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
//double* XYZ_t(double x, double y, double z) {
CUDABP_DEV_FUNC_DECL void XYZ_t(double *XYZ) {

	double x = XYZ[0];
	double y = XYZ[1];
	double z = XYZ[2];
	
	#if SWAP_x_y==1

	tmp = x;
	x = y;
	y = tmp;

	#endif

	#if FULL_ROTATION==1

	Rot_inverse(XYZ);

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
 *
 */
CUDABP_DEV_FUNC_DECL double lambda_implicit(double x, double y, double z) {

	double XYZ[3];
	
	XYZ[0] = x;
	XYZ[1] = y;
	XYZ[2] = z;
	
	XYZ_t(XYZ);

	// Orig Code; is unrolled -> return lambda(XYZ_[0], XYZ_[1], XYZ_[2]);
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
	
	XYZ[0] = x;
	XYZ[1] = y;
	XYZ[2] = z;
	
	XYZ_t(XYZ);

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
CUDABP_DEV_FUNC_DECL double f2_implicit(double x, double y, double z) {

	double XYZ[3];
		
	XYZ[0] = x;
	XYZ[1] = y;
	XYZ[2] = z;
	
	XYZ_t(XYZ);

	// Old codel; is unrolled -> return f2(XYZ_[0], XYZ_[1], XYZ_[2]);
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

	 /*	 Initialize the array diff2_ to zeros	*/
	for (i=0; i<MAX_ORDER; i++) for (j=0; j<MAX_ORDER; j++) diff2_[i][j] = 0.0f;

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
CUDABP_DEV_FUNC_DECL void solve_ivp2(double x0, double y0, double z0, double *ivp) {

	int i,j;
	double diff2_[MAX_ORDER][MAX_ORDER];

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


extern "C" 
int initializeGPU(int n__, int number_of_terms_, int *order_in_X_, int *order_in_Y_, int *order_in_Z_, double *lambda__,
		 	double *a1__, double *a2__, double *bottom_plane__, double *rot__, int x_inc_, int y_inc_,
			int _argc, char **_argv, int devId) {

	static int hwInitFlag = 0;
	dim3 grids(1, 1);
    	dim3 block(1, 1);
	void *symPtr;
	int i;
	
	if(!hwInitFlag)
	{			
		// Init the GPUs
		printf("\ninitializeGPU(): Initializing GPU # %d ... \n", devId);
		//CUT_DEVICE_INIT(_argc, _argv);
		
		// since the command line parsing for txbr is incompatible with 
		// the cutil command parsing, we will set the device to be used 
		// explicitly here. The id for the device should come from a config
		// file or the command line
	
		printf("\n\n====================================================");
		printf("\ninitializeGPU(): Resetting device usage for device %d : ", devId);
		cudaDeviceProp deviceProp;                                               
		cudaGetDeviceProperties(&deviceProp, devId);       
		if (deviceProp.major < 1) {                                              
			fprintf(stderr, "\ninitializeGPU()(): Error: device does not support CUDA.\n");
			return -1;
		}
		
		printf("%s\n", deviceProp.name);
		printf("\nTotal Global Memory = %d", deviceProp.totalGlobalMem);
		printf("\nShared memory per block = %d", deviceProp.sharedMemPerBlock);
		printf("\nWarp Size = %d", deviceProp.warpSize);
		printf("\nMax threads per block = %d", deviceProp.maxThreadsPerBlock);
		printf("\nNumber of multiprocessors = %d", deviceProp.multiProcessorCount);
		
		printf("\n====================================================\n");
		
		cudaSetDevice(devId);
		
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
	
	// initialize simple variables on the device
	initializeGPU_d<<< grids, block >>>(n__, number_of_terms_, x_inc_);
	
	return 1;
}

CUDABP_ENTRY_FUNC_DECL void initializeGPU_d(int n__, int number_of_terms_, int x_inc_) 
{
	// nothing much happening here except for simple local copies
	n_ = n__;
	x_inc = x_inc_;
	number_of_terms = number_of_terms_;

	// I guess you need to do this to make sure the modified shared mem variables
	// are visible to all threads.
	__syncthreads();
	
	return;
}

/*
 * Caculate the segment size for the propagation error to remain less than eps
 */
extern "C"
int calculate_segment_size(int order, int x0, int y0, int z0, int x_inc, int x1, double *coefficients, double eps) {

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
	calculate_segment_size_d <<< grids, block >>> (order, x0, y0, z0, x_inc, x1, eps);
	
	// get the symbol pointer to the result on the GPU and copy the result back to host
	cudaGetSymbolAddress((void **)&symPtr, "segment_size_d");
	cErr = cudaMemcpy((void *)&l_segment_size, symPtr, sizeof(int), cudaMemcpyDeviceToHost);
	
 	return l_segment_size;

}


CUDABP_ENTRY_FUNC_DECL void calculate_segment_size_d(int order, int x0, int y0, int z0, int x_inc,
	int x1, double eps) {

	double q1_[MAX_ORDER];
	int ix_0 = x0-1;
	int ix_1 = x1-1;

	int blocksize_x = (ix_1-ix_0+1)/x_inc;
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
	
	segment_size_d = segmentSize*x_inc;
	
	return;
}


//
// This is the host interface to the evaluateBlock_d() kernel which is used for the bulk of 
// computing. We can hide the actual kernel invocation here. Not sure yet but for 
// performance tweaks, I might have to so multiple invokations with smaller tasks
extern "C"
int evaluateBlock( float *block, int x0, int y0, int z0, int x_inc, int y_inc, int z_inc, 
		int x1, int y1, int z1, float* slice_data_f, int x0_src, int y0_src, 
		int nx_src, int ny_src, int nseg, int seg_size, unsigned int sliceMemSize, 
		unsigned int blockMemSize, int itilt) {

	static unsigned int l_sliceMemSize = 0, l_blockMemSize = 0;
	int i;
	
	// If no memory needs to be allocated, do it here 
//	if(slice_data_d_ptr == NULL || block_d_ptr == NULL)
	if(sliceMemSize > l_sliceMemSize)
	{		
		if(slice_data_d_ptr) free(slice_data_d_ptr);
		
		// Allocate memory for slice data
		if(CUDA_SAFE_CALL((cErr = cudaMalloc((void **)&slice_data_d_ptr, sliceMemSize))
			== cudaErrorMemoryAllocation))
		{
		    printf("evaluateBlock(): Error allocating %d bytes for slice on GPU", 
						sliceMemSize);
		    return -1;
		}
		
		l_sliceMemSize = sliceMemSize;		
	}
	
	if(blockMemSize > l_blockMemSize)
	{
	
		if(block_d_ptr) free(block_d_ptr);
		
		// Allocate memory for output block
		if(CUDA_SAFE_CALL((cErr = cudaMalloc((void **)&block_d_ptr, blockMemSize))
			== cudaErrorMemoryAllocation))
		{
		    printf("evaluateBlock(): Error allocating %d bytes for block on GPU", 
						blockMemSize);
		    return -1;
		}
		
		l_blockMemSize = blockMemSize;		
	}
	
	// set the block values to zero
	if(itilt == 0)
		cudaMemset(block_d_ptr, 0, blockMemSize);
	
	// Copy the slice data to the GPU memory
	cErr = cudaMemcpy(slice_data_d_ptr, (void *)slice_data_f, sliceMemSize, 
		cudaMemcpyHostToDevice);

	if(cErr != cudaSuccess)
	{
		printf("\nevaluateBlock(): %s", cudaGetErrorString(cErr));
		return -1;
	}
		
	evaluateBlock_d <<< CUDABP_GRID_SIZE, CUDABP_BLOCK_SIZE >>>( (float *)block_d_ptr, 
		x0, y0, z0, x_inc, y_inc, z_inc, 
		x1, y1, z1, (float *)slice_data_d_ptr, x0_src, y0_src, 
		nx_src, ny_src, nseg, seg_size);

	CUT_CHECK_ERROR("evaluateBlock(): Kernel execution failed");

	// copy back the computed block array. Assuming things went well
	cErr = cudaMemcpy(block, block_d_ptr, blockMemSize, 
			cudaMemcpyDeviceToHost);
	
	
	if(cErr != cudaSuccess)
	{
		printf("\nevaluateBlock(): %s", cudaGetErrorString(cErr));
		return -1;
	}
	
	return 1;
}


CUDABP_ENTRY_FUNC_DECL void evaluateBlock_d( float *block, int x0, int y0, int z0, int x_inc, int y_inc, int z_inc, 
		int x1, int y1, int z1, float* slice_data_f, int x0_src, int y0_src, 
		int nx_src, int ny_src, int nseg, int seg_size)
{
	
	int nx_dest = (int)floor((double)(x1-x0+1)/(double)x_inc);
	int ny_dest = (int)floor((double)(y1-y0 + 1)/(double)y_inc);
	int blocksize = (int)floor((double)(z1-z0+1)/(double)z_inc);
	
	int nx_eff = (x1-x0+1)/x_inc;
	int ny_eff = (y1-y0+1)/y_inc;
//	int nx_eff = nx_dest;
//	int ny_eff = ny_dest;
		
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
	int xv = 0, yv = 0;
	//int index;
	int m;
	
	double q1_[MAX_ORDER];
	double q2_[MAX_ORDER];
	
	// first thing find out the index of the core
	unsigned int thIndex = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int totalNoOfThreads = gridDim.x * blockDim.x;

	for (isect=0; isect<blocksize; isect++) {	/*	 Step through the output block	*/

		z_section = z0 + z_inc*isect;	/*	 z-section coordinate of the output slice */	

		// index = isect*nx_dest*ny_dest;	/*	counter for the block	*/

		/*	 Loop along y image coordinate */
		
		for (iy=y0 + thIndex; iy<y1+1; iy += y_inc + totalNoOfThreads - 1) {
		//for (iy=y0-1 + thIndex; iy<y1; iy += y_inc + totalNoOfThreads - 1) {
		//OC -> for (iy=y0-1; iy<y1; iy+=y_inc) {
		/* Now increment along x direction. Iterative evaluation scheme for each segment.	*/

			for (isegment=0; isegment<nseg; isegment++) {

				//istart = x0 + isegment*seg_size - 1;
				istart = x0 + isegment*seg_size;
				istop = TXBR_MIN( istart + seg_size, x1+1 );

				#if OFFSET==1

				//XYZ[0] = istart + 1;
				//XYZ[1] = iy + 1;
				XYZ[0] = istart;
				XYZ[1] = iy;
				XYZ[2] = z_section;

				#else

				//XYZ[0] = istart;
				//XYZ[1] = iy;
				XYZ[0] = istart-1;
				XYZ[1] = iy-1;
				XYZ[2] = z_section;

				#endif

				xproj = f1_implicit(XYZ[0],XYZ[1],XYZ[2]);
				yproj = f2_implicit(XYZ[0],XYZ[1],XYZ[2]);

				lam = lambda_implicit(XYZ[0],XYZ[1],XYZ[2]);
				dely_lam = lambda_implicit(XYZ[0]+x_inc,XYZ[1],XYZ[2])-lam;

				// In the last version of the code the solve_ivp functions
				// were unrolled
				solve_ivp1(XYZ[0],XYZ[1],XYZ[2],q1_);
				solve_ivp2(XYZ[0],XYZ[1],XYZ[2],q2_);

				//for (ix=istart; ix<istop+1; ix=ix+x_inc) {
				for (ix=istart; ix<istop; ix=ix+x_inc) {

					xproj = q1_[0];
					yproj = q2_[0];

					xproj /= lam;
					yproj /= lam;

					/* Eventually reswap xproj and yproj like in the original code	*/

					//#if OFFSET==1
					#if OFFSET==0

					xproj -= 1.0;
					yproj -= 1.0;

					#endif

					xproj -= x0_src;
					yproj -= y0_src;

//					xv = floor(xproj);	// floor causes problem ??
//					yv = floor(yproj);

					xv = (int)xproj;
					yv = (int)yproj;

					xdiff = xproj - xv;	/*	 Weight for Linear Interpolation	*/
					ydiff = yproj - yv;

					#if TEST==1

					if ((ix-istart)==0.5*(istop-istart) && (iy-y0+1)==0.5*(y1-y0+1)) {

						printf( "ix=%i iy=%i  isect=%i  index=%i   xproj=%f  yproj=%f  xv=%i  yv=%i\n",
								q1[0],q2[0],f1_implicit(ix,iy,z_section), f2_implicit(ix,iy,z_section));

						printf( "ix=%i iy=%i  isect=%i  index=%i   xproj=%f  yproj=%f  xv=%i  yv=%i\n",
								 ix,iy,isect,index,xproj,yproj,xv,yv);

					}

					#endif

					/*	 Calculate linear interpolation error, making sure that we are within the coordinate limits	*/

				//	if (xv>=0 && xv<nx_src-1 && yv>=0 && yv<ny_src-1) {
				//	if (xv>=0 && xv<nx_src && yv>=0 && yv<ny_src) {
					if (xv>=0 && xv<nx_src-1 && yv>=0 && yv<ny_src-1) {					
//							block[isect*nx_eff*ny_eff + (iy - y0 + 1)*nx_eff+(ix - istart)] += 
//							INTERP2(slice_data_f, nx_src, ny_src, xv, yv, xdiff, ydiff);
							//block[isect*nx_eff*ny_eff + (iy - y0 + 1)*nx_eff+(ix - x0 + 1)] += 
							//INTERP2(slice_data_f, nx_src, ny_src, xv, yv, xdiff, ydiff);
							block[isect*nx_eff*ny_eff + (iy - y0)*nx_eff+(ix - x0)] += 
							INTERP2(slice_data_f, nx_src, ny_src, xv, yv, xdiff, ydiff);							
					}

					for (m=n_-1;m>=0;m--) {
						q1_[m] += q1_[m+1];
						q2_[m] += q2_[m+1];
					}

					lam += dely_lam;

					//index++;

				} // for (ix ..
			} // for (isegment ..
		} // for(iy ..
	} // for (isect ..
	
	// sync threads here before returning
	__syncthreads();


	return;
	
} // End of evaluateBlock()




extern "C"
void releaseGPUMemory()
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



