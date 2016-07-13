#include <Python.h>
#include <numpy/arrayobject.h>

#include "txbrutil.h"


static PyObject *mapper(PyObject *self, PyObject *args, int direction) {

    char *directory, *basename, *work_directory, *cfg_directory;
    int i,j, itilt;
    int use_fixed_segment_size;
    char *data1, *data2, *lambda;

    PyObject *skip = NULL;
    PyObject *projmap = NULL;
    PyObject *reconstruction = NULL;

    if (!PyArg_ParseTuple(args, "ssssOiOO", &directory, &basename, &work_directory, &cfg_directory,
            &skip, &use_fixed_segment_size, &projmap, &reconstruction)) return NULL;

    /* Build the projection map*/

    PathProjection *prj = (PathProjection*)malloc(sizeof(PathProjection));

    if (!prj) {
        txbr_error( stderr, "ERROR: reconstruct() projmap - getting memory.\n" );
        return;
    }

    /* Register the skipping slices in the projection map*/

    for (itilt=0; itilt<MAX_TILT; itilt++) {
        prj->skipView[itilt] = 0;
    }

    if (!PyList_Check(skip)) return NULL;
    Py_ssize_t nskip = PyList_Size(skip);

    for (i=0; i<nskip; i++) {
        itilt = (int)PyInt_AsLong(PyList_GetItem(skip, i));
        prj->skipView[itilt] = 1;
    }

    PyArrayObject *x_coeff = (PyArrayObject *)PyObject_GetAttrString(projmap, "x_coefficients");
    PyArrayObject *y_coeff = (PyArrayObject *)PyObject_GetAttrString(projmap, "y_coefficients");
    PyArrayObject *scaling_coefficients = (PyArrayObject *)PyObject_GetAttrString(projmap, "scaling_coefficients");

    npy_intp *dimensions = x_coeff->dimensions;

    sprintf( prj->label, "%s", PyString_AsString(PyObject_GetAttrString(projmap, "label")));

    prj->order = (int)PyInt_AsLong(PyObject_GetAttrString(projmap, "order"));
    prj->number_of_terms = (int)PyInt_AsLong(PyObject_GetAttrString(projmap, "numberOfTerms"));
    prj->numberOfTilts = (int)dimensions[0];
    
    data1 = (char *)PyArray_DATA(x_coeff);
    data2 = (char *)PyArray_DATA(y_coeff);
    lambda = (char *)PyArray_DATA(scaling_coefficients);

    for (i = 0; i < prj->numberOfTilts; i++) {
        for (j = 0; j < 4; j++) {
            prj->lambda[i][j] = *(double *) (lambda + i * scaling_coefficients->strides[0] + j * scaling_coefficients->strides[1]);
        }
        for (j = 0; j < prj->number_of_terms; j++) {
            prj->coefficients_1[i][j] = *(double *) (data1 + i * x_coeff->strides[0] + j * x_coeff->strides[1]);
            prj->coefficients_2[i][j] = *(double *) (data2 + i * y_coeff->strides[0] + j * y_coeff->strides[1]);
        }
    }

    Py_DECREF(x_coeff);
    Py_DECREF(y_coeff);
    Py_DECREF(scaling_coefficients);

    /* Build the setup of the reconstruction */

    TxBRsetup *setup;

    setup = (TxBRsetup *)malloc(sizeof(TxBRsetup));

    if (!setup) {
        txbr_error(stderr, "ERROR: reconstruct() setup - getting memory.\n");
        return;
    }

    setup->blocksize = (int)PyInt_AsLong(PyObject_GetAttrString(reconstruction, "sizeOfBlock"));
    setup->use_fixed_segment_size = use_fixed_segment_size;

    PyObject *origin = PyObject_GetAttrString(reconstruction, "origin");
    PyObject *end = PyObject_GetAttrString(reconstruction, "end");
    PyObject *increment = PyObject_GetAttrString(reconstruction, "increment");
    PyObject *bottomPlane = PyObject_GetAttrString(reconstruction, "bottomPlane");
    PyObject *topPlane = PyObject_GetAttrString(reconstruction, "topPlane");
    PyObject *scale = PyObject_GetAttrString(reconstruction, "scale");

    setup->x_0 = (double)PyFloat_AsDouble(PyObject_GetAttrString(origin, "x"));
    setup->y_0 = (double)PyFloat_AsDouble(PyObject_GetAttrString(origin, "y"));
    setup->z_0 = (double)PyFloat_AsDouble(PyObject_GetAttrString(origin, "z"));

    setup->x_1 = (double)PyFloat_AsDouble(PyObject_GetAttrString(end, "x"));
    setup->y_1 = (double)PyFloat_AsDouble(PyObject_GetAttrString(end, "y"));
    setup->z_1 = (double)PyFloat_AsDouble(PyObject_GetAttrString(end, "z"));

    setup->x_inc = (double)PyFloat_AsDouble(PyObject_GetAttrString(increment, "x"));
    setup->y_inc = (double)PyFloat_AsDouble(PyObject_GetAttrString(increment, "x"));
    setup->z_inc = (double)PyFloat_AsDouble(PyObject_GetAttrString(increment, "x"));

    setup->plane_coeffs1[0] = (double)PyFloat_AsDouble(PyObject_GetAttrString(bottomPlane, "d"));
    setup->plane_coeffs1[1] = (double)PyFloat_AsDouble(PyObject_GetAttrString(bottomPlane, "a"));
    setup->plane_coeffs1[2] = (double)PyFloat_AsDouble(PyObject_GetAttrString(bottomPlane, "b"));

    setup->plane_coeffs2[0] = (double)PyFloat_AsDouble(PyObject_GetAttrString(topPlane, "d"));
    setup->plane_coeffs2[1] = (double)PyFloat_AsDouble(PyObject_GetAttrString(topPlane, "a"));
    setup->plane_coeffs2[2] = (double)PyFloat_AsDouble(PyObject_GetAttrString(topPlane, "b"));

    setup->sx = (double)PyFloat_AsDouble(PyObject_GetAttrString(scale, "x"));
    setup->sy = (double)PyFloat_AsDouble(PyObject_GetAttrString(scale, "y"));
    setup->sz = (double)PyFloat_AsDouble(PyObject_GetAttrString(scale, "z"));

    setup->z_start = setup->z_0;

    Py_DECREF(origin);
    Py_DECREF(end);
    Py_DECREF(increment);
    Py_DECREF(bottomPlane);
    Py_DECREF(topPlane);
    Py_DECREF(scale);

    sprintf( setup->basename, "%s", PyString_AsString(PyObject_GetAttrString(reconstruction, "basename")));
    sprintf( setup->directory, "%s", PyString_AsString(PyObject_GetAttrString(reconstruction, "directory")));

    /* Print out informations	*/

    print_TxBR_setup(setup);
    printf("\n");

    print_path_projection(prj);
    printf("\n");

    /* Run the reconstruction routine */

    if (direction==1) {
        do_full_reconstruction( directory, basename, work_directory, setup, prj);
    } else if (direction==0) {
        do_full_projection( directory, basename, work_directory, setup, prj);
    }

    free(prj);

    free(setup);

    return Py_BuildValue("i", 1);

}


static PyObject *reconstruct(PyObject *self, PyObject *args) {

    return mapper( self, args, 1 );

}


static PyObject *project(PyObject *self, PyObject *args) {

    return mapper( self, args, 0 );

}


static PyMethodDef BackprojectionMethods[] = {
	{ "reconstruct", reconstruct, METH_VARARGS, "Reconstruct a volume from a tilt series." },
	{ "project", project, METH_VARARGS, "Project a volume on a tilt series." },
	{NULL, NULL, 0, NULL} /* Sentinel */
};


PyMODINIT_FUNC initbckprj(void) {

	(void)Py_InitModule("bckprj", BackprojectionMethods);

}

int main(int argc, char *argv[]) {

	Py_SetProgramName(argv[0]);	/* Pass argv[0] to the Python interpreter */

	Py_Initialize(); /* Initialize the Python interpreter.	*/

	initbckprj();	/* Add a static module */

	return 1;

}
