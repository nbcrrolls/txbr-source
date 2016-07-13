#include <Python.h>
#include <numpy/arrayobject.h>
#include <iostream>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

static PyObject *coefficients(PyObject *self, PyObject *args) {

    PyObject *f;
    PyObject *variables;

	if (!PyArg_ParseTuple(args, "OO", &f, &variables)) return NULL;

	cout << f << endl;



	symbol x("x"), y("y");
	ex poly;

	for (int i=0; i<3; ++i) {
		poly += factorial(i+16)*pow(x,i)*pow(y,2-i);
     	cout << poly << endl;
	}

	for (int i=poly.ldegree(x); i<=poly.degree(x); ++i) {
		cout << "The x^" << i << "-coefficient is "
         	 << poly.coeff(x,i) << endl;
	}

	return Py_BuildValue( "i", 1 );

}

static PyMethodDef GINACPLUGMethods[] = {
	{ "coefficients", coefficients, METH_VARARGS, "Return the GiNaC coefficients of a polynom." },
	{NULL, NULL, 0, NULL} /* Sentinel */
};


PyMODINIT_FUNC initginacplug(void) {

	(void)Py_InitModule("ginacplug", GINACPLUGMethods);

	import_array();

}

int main(int argc, char *argv[]){

	Py_SetProgramName(argv[0]);	/* Pass argv[0] to the Python interpreter */

	Py_Initialize(); /* Initialize the Python interpreter.	*/

	initginacplug();	/* Add a static module */

	return 1;

}