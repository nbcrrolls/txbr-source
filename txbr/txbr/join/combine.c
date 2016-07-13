#include <Python.h>
#include "txbrutil.h"

#define MAX_NUMBER_OF_FILES 50

/*
 *	Sum series together.
 */
static PyObject *average(PyObject *self, PyObject *args) 
{

	char *fileout;
	PyObject *filesin;

	if (!PyArg_ParseTuple(args, "Os", &filesin, &fileout)) return NULL;

	if (!PyList_Check(filesin)) return NULL;

	char **filenames = malloc(MAX_NUMBER_OF_FILES*FILENAME_LEN*sizeof(char));

	if (!filenames) {
		txbr_error(stderr, "ERROR: icombine - getting memory.\n");
		return MEMORY_ERROR;
	}

	float *weight = malloc(MAX_NUMBER_OF_FILES*sizeof(float));

	if (!weight) {
		txbr_error(stderr, "ERROR: icombine - getting memory.\n");
		return MEMORY_ERROR;
	}
	
	Py_ssize_t n = PyList_Size(filesin);
	int i;
	
	for (i=0; i<n; i++)
	{
		filenames[i] = PyString_AsString(PyList_GetItem(filesin,i));
		weight[i] = 1.0;
	}
	
	combine_stacks( n, filenames, weight, fileout );

	free(filenames);
	free(weight);

	return Py_BuildValue( "i", 1 );

}


/*
 *	Combine series stacks together.
 */
static PyObject *merge(PyObject *self, PyObject *args) 
{

	char *fileout;
	PyObject *fileDict;

	if (!PyArg_ParseTuple(args, "Os", &fileDict, &fileout)) return NULL;

	if (!PyDict_Check(fileDict)) return NULL;

	char **filenames = malloc(MAX_NUMBER_OF_FILES*FILENAME_LEN*sizeof(char));

	if (!filenames) {
		txbr_error(stderr, "ERROR: icombine - getting memory.\n");
		return MEMORY_ERROR;
	}

	float *weight = malloc(MAX_NUMBER_OF_FILES*sizeof(float));

	if (!weight) {
		txbr_error(stderr, "ERROR: icombine - getting memory.\n");
		return MEMORY_ERROR;
	}

	int count = 0;

	PyObject *key, *value;
	Py_ssize_t pos = 0;

	while (PyDict_Next(fileDict, &pos, &key, &value)) {
		filenames[count] = PyString_AsString(key);
		weight[count] = (float)PyFloat_AsDouble(value);
		count++;
	}

	combine_stacks(count, filenames, weight, fileout);

	free(filenames);
	free(weight);

	return Py_BuildValue( "i", 1 );

}

/*
 *	Cross Validate series stacks.
 */
static PyObject *crossval(PyObject *self, PyObject *args) 
{

	char *fileout, *directory;
	PyObject *fileList;

	if (!PyArg_ParseTuple(args, "sOs", &directory, &fileList, &fileout)) return NULL;

	if (!PyList_Check(fileList)) return NULL;

	char **filenames = malloc(MAX_NUMBER_OF_FILES*FILENAME_LEN*sizeof(char));

	if (!filenames) {
		txbr_error(stderr, "ERROR: icombine - getting memory.\n");
		return MEMORY_ERROR;
	}

	Py_ssize_t count = PyList_Size(fileList);
	int i;

	for (i=0; i<count; i++) {

		filenames[count] = PyString_AsString(PyList_GetItem(fileList,i));

	}

	cross_validate((int)count, directory, filenames, "p", fileout);

	free(filenames);

	return Py_BuildValue( "i", 1 );

}


static PyMethodDef CombineMethods[] = {
	{ "average", average, METH_VARARGS, "Average tomographic series together." },
	{ "merge", merge, METH_VARARGS, "Merge tomographic series together." },
	{ "crossval", merge, METH_VARARGS, "Cross Validate series." },
	{NULL, NULL, 0, NULL} /* Sentinel */
};


PyMODINIT_FUNC initcombine(void) {

	(void)Py_InitModule("combine", CombineMethods);

}


int main(int argc, char *argv[]) {

	Py_SetProgramName(argv[0]);	/* Pass argv[0] to the Python interpreter */

	Py_Initialize(); /* Initialize the Python interpreter.	*/

	initcombine();	/* Add a static module */

	return 1;

}
