#include <Python.h>

#include "txbrutil.h"

static PyObject *reconstruct(PyObject *self, PyObject *args) {

	char *directory, *basename, *work_directory, *cfg_directory;
	int override_z;
	int deviceId;
	float zstart, zinc, zstop;

	if (!PyArg_ParseTuple(args, "sssifffi", &directory, &basename, &work_directory, &cfg_directory, &override_z, &zstart, &zstop, &zinc, &deviceId)) return NULL;

	printf("%s\n",&directory[0]);
	printf("%s\n",&basename[0]);

	TxBRsetup *setup;
	PathProjection *projection;

	setup  = (TxBRsetup *)malloc(sizeof(TxBRsetup));

	if (!setup) {
		txbr_error(stderr, "ERROR: reconstruct() setup - getting memory.\n");
		return;
	}

	projection  = (PathProjection*)malloc(sizeof(PathProjection));

	if (!projection) {
		txbr_error(stderr, "ERROR: reconstruct() projection - getting memory.\n");
		return;
	}

	char configuration_file_name[FILENAME_LEN];
	//sprintf(configuration_file_name, "%s/%s.txbr", directory, basename);
	sprintf(configuration_file_name, "%s/%s.txbr", cfg_directory, basename);

	load_configuration(setup, projection, configuration_file_name);

	/* Give the boundaries condition	*/

	if (override_z==1) {
		setup->z_0 = zstart;
		setup->z_1 = zstop;
	}

	/* Print out informations	*/

	print_TxBR_setup(setup);
	printf("\n");

	// projection->skipView[21] = 1;

	print_path_projection(projection);
	printf("\n");

	/* Run the reconstruction routine */

	do_full_reconstruction(directory, basename, work_directory, setup, projection, deviceId);

	free(projection);

	free(setup);

	return Py_BuildValue( "i", 1 );

}


static PyMethodDef BackprojectionMethods[] = {
	{ "reconstruct", reconstruct, METH_VARARGS, "Reconstruct a volume." },
	{NULL, NULL, 0, NULL} /* Sentinel */
};


PyMODINIT_FUNC initbckprj_cu(void) {

	(void)Py_InitModule("bckprj_cu", BackprojectionMethods);

}

int main(int argc, char *argv[]) {

	Py_SetProgramName(argv[0]);	/* Pass argv[0] to the Python interpreter */

	Py_Initialize(); /* Initialize the Python interpreter.	*/

	initbckprj_cu();	/* Add a static module */

	return 1;

}
