#include <Python.h>

static PyObject *filter_3D(PyObject *self, PyObject *args) {

    const char *directory, *basename;

	if (!PyArg_ParseTuple(args, "ss", &directory, &basename)) return NULL;

    int sts;

    sts = filter( directory, basename);

    return NULL;

}

static PyMethodDef TxBRMethods[] = {

    {"filter_3D",  filter_3D, METH_VARARGS, "Execute a 3D filtering."},
	{NULL, NULL, 0, NULL} /* Sentinel */

};

PyMODINIT_FUNC initfilter(void) {

	(void)Py_InitModule("filter", TxBRMethods);

}

void main(int argc, char *argv[]){

	Py_SetProgramName(argv[0]);	/* Pass argv[0] to the Python interpreter */

	Py_Initialize(); /* Initialize the Python interpreter.	*/

	initfilter();	/* Add a static module */

}
