#include <Python.h>
#include <numpy/arrayobject.h>

#include "txbrutil.h"
#include "filter_1D.h"


/*
 * Create a 2D remap structure object from python (double) arrays
 */
remap_2D_params_handle create_remap_2D_params_from_PY_array ( 
				float local_magnification,
				int order,
				PyArrayObject *power_order,
				PyArrayObject *map2D )
{

	npy_intp *pow_order_dimensions, *map2D_dimensions;
	
	pow_order_dimensions = power_order->dimensions;
	map2D_dimensions = map2D->dimensions;
	
    int i, j, k, n1, n2, n3, index;
    char *data;
    
	n1 = (int)pow_order_dimensions[0];
 	n2 = (int)pow_order_dimensions[1];
    
 	data = (char *)PyArray_DATA(power_order);

   	real_type *power_order_ = (real_type *)malloc(n1*n2*sizeof(real_type));
   	
   	index = 0;
    
	for (i=0; i<n1; i++)
	{
		for (j=0; j<n2; j++)
		{
			power_order_[index++] = (real_type)*(double *)(data + i*power_order->strides[0] + j*power_order->strides[1]);
		}
	}
    
    n1 = (int)map2D_dimensions[0];
    n2 = (int)map2D_dimensions[1];
    n3 = (int)map2D_dimensions[2];
    
 	real_type *map2D_ = (real_type *)malloc(n1*n2*n3*sizeof(real_type));
   
    index = 0;
    
	data = (char *)PyArray_DATA(map2D);
	
	for (i=0; i<n1; i++)
	{
		for (j=0; j<n2; j++)
		{
			for (k=0; k<n3; k++)
			{
				map2D_[index++] = (real_type)*(double *)(data + i*map2D->strides[0] + j*map2D->strides[1] + k*map2D->strides[2]);
			}
		}
	}
	
	remap_2D_params_handle remap_2D_params = create_remap_2D_params_from_data_copy
	(
	    local_magnification,
	    order,
	    power_order->dimensions[1],
	    power_order->dimensions[0],
	    power_order->dimensions[1],
	    power_order_, 
	    map2D_dimensions[0],
	    map2D_dimensions[1],
	    map2D_dimensions[2],
	    map2D_dimensions[0],
	    map2D_
	);
	
	//displayRemap2DParameters(remap_2D_params);
	
	return remap_2D_params;
	
}


/*
 * Filter the micrographs
 */
static PyObject *rot_filter(PyObject *self, PyObject *args) {

	char *filepath_in, *filepath_out;
	float local_magnification = 1.0;
	float angle;
	
	if (!PyArg_ParseTuple(args, "ssff", &filepath_in, &filepath_out, 
								&local_magnification, &angle)) return NULL;
	
	printf("Filtering with no remap with angle %f!\n",angle);	
	printf("Local Magnification: %f\n",local_magnification);
	
    filesystem_params_handle filesystem_params = NULL;
    symmetrize_2D_params_handle symmetrize_2D_params = NULL;
    rotate_2D_params_handle rotate_2D_params = NULL;
    filter_1D_params_handle filter_1D_params = NULL;
    edgetaper_2D_params_handle edgetaper_2D_params = NULL;

    filesystem_params = create_filesystem_params_from_data_copy(filepath_in, "-", filepath_out);

    symmetrize_2D_params = create_symmetrize_2D_params_from_data_copy(1);
    
	filter_1D_params = create_filter_1D_params_from_strings(
		"SheppLogan",
    	"251",
    	NULL,
    	NULL );
        	
        	
 	real_type rotate_2D_center[2] = {0.0, 0.0};
    rotate_2D_params = create_rotate_2D_params_from_angle_data_copy(local_magnification,angle,rotate_2D_center);

    edgetaper_2D_params = create_edgetaper_2D_params(1);

	/* Run the filtering routine */
    
	projection_series_rotate_2D_filter_1D_rotate_2D (
        filesystem_params,
        symmetrize_2D_params,
        rotate_2D_params,
        filter_1D_params,
        edgetaper_2D_params
 	);

	return Py_BuildValue( "i", 1 );

}


/*
 * Remap and Filter the micrographs
 */
static PyObject *remap_filter(PyObject *self, PyObject *args) {

	char *filepath_in, *filepath_out;
	float local_magnification = 1.0;
	int order = 1.0;
	PyArrayObject *power_order, *map2D;
	
	if (!PyArg_ParseTuple(args, "ssfiOO", &filepath_in, &filepath_out, 
						&local_magnification, &order, &power_order, &map2D)) return NULL;
	
	printf("Filtering with 2D remap (order %i)!\n",order);
	printf("Local Magnification: %f\n",local_magnification);
	
    filesystem_params_handle filesystem_params = NULL;
    symmetrize_2D_params_handle symmetrize_2D_params = NULL;
    remap_2D_params_handle remap_2D_params = NULL;
    filter_1D_params_handle filter_1D_params = NULL;
    edgetaper_2D_params_handle edgetaper_2D_params = NULL;

    filesystem_params = create_filesystem_params_from_data_copy(filepath_in, NULL, filepath_out);

    //symmetrize_2D_params = create_symmetrize_2D_params();
    symmetrize_2D_params = create_symmetrize_2D_params_from_data_copy(1);
    
	filter_1D_params = create_filter_1D_params_from_strings(
		"SheppLogan",
    	"251",
    	NULL,
    	NULL );
    
    remap_2D_params = create_remap_2D_params_from_PY_array( local_magnification, order, power_order, map2D );
   
    edgetaper_2D_params = create_edgetaper_2D_params();

	/* Run the filtering routine */
	
	projection_series_remap_2D_filter_1D_inv_remap_2D (
        filesystem_params,
        symmetrize_2D_params,
        remap_2D_params,
        filter_1D_params,
        edgetaper_2D_params
	);

	return Py_BuildValue( "i", 1 );

}


static PyMethodDef FilterMethods[] = {
	
	{ "rot_filter", rot_filter, METH_VARARGS, "Filter Images prior to the backprojection." },
	{ "remap_filter", remap_filter, METH_VARARGS, "Remap and Filter Images prior to the backprojection." },
	{NULL, NULL, 0, NULL} /* Sentinel */
	
};


PyMODINIT_FUNC initprefil(void) {

	(void)Py_InitModule("prefil", FilterMethods);

}

int main(int argc, char *argv[]) {

	Py_SetProgramName(argv[0]);	/* Pass argv[0] to the Python interpreter */

	Py_Initialize(); /* Initialize the Python interpreter.	*/

	initprefil();	/* Add a static module */

	return 1;

}
