#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>


#include "txbrutil.h"
#include "mrcslice.h"

int header_info(char* mrcFileName);
int get_mode_of(char* mrcFileName, int *mode);
int get_size_of(char* mrcFileName, int *nx, int *ny, int *nz);
int get_tilt_angles(char* mrcFileName, int nangle, float* angles);
int get_scale_of(char* mrcFileName, float *nx, float *ny, float *nz);
int update_scale_of(char* mrcFileName, float sx, float sy, float sz);
int get_origin_of(char* mrcFileName, float *xorg, float *yorg, float *zorg);
int update_origin_of(const char* mrcFileName, const float xorg, const float yorg, const float zorg);
int get_length_of(const char* mrcFileName, float* xlen, float* ylen, float* zlen);
int update_length_of(const char* mrcFileName, const float xlen, const float ylen, const float zlen);
int get_grid_origin_of(const char* mrcFileName, int* nxstart, int* nystart, int* nzstart);
int update_grid_origin_of(const char* mrcFileName, const int nxstart, const int nystart, const int nzstart);
int get_grid_size_of(const char* mrcFileName, int* mx, int* my, int* mz);
int update_grid_size_of(const char* mrcFileName, const int mx, const int my, const int mz);
int get_description_of(char* mrcFileName, int *nlabels, char** labels);
int update_description_of(const char* mrcFileName, const int nlabl, const char** labels);
int get_statistics_of(char* mrcFileName, float *m0, float *m1, float *m2, float *min, float *max);


static PyObject *mrc_header(PyObject *self, PyObject *args) {

    char *fileName;

	if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    header_info(fileName);

    return Py_None;

}


static PyObject *mrc_mode(PyObject *self, PyObject *args) {

    char *fileName;

	if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    int sts, mode=0;

    sts = get_mode_of(fileName, &mode);

	if (sts!=NO_ERROR) return NULL;

    return Py_BuildValue( "i", mode );

}


static PyObject *mrc_size(PyObject *self, PyObject *args) {

    char *fileName;

	if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    int sts, nx=0, ny=0, nz=0;

    sts = get_size_of(fileName, &nx, &ny, &nz);

    if (sts!=NO_ERROR) return NULL;

    return Py_BuildValue( "(i,i,i)", nx, ny, nz );

}


static PyObject *mrc_scale(PyObject *self, PyObject *args) {

    char *fileName;

	if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    int sts;
    float sx=0.0, sy=0.0, sz=0.0;

    sts = get_scale_of(fileName, &sx, &sy, &sz);

    if (sts!=NO_ERROR) return NULL;

    return Py_BuildValue( "(f,f,f)", sx, sy, sz );

}

/*
 * Update the scale.
 */
static PyObject *update_scale(PyObject *self, PyObject *args) {

    char *fileName;
 	float sx, sy, sz;

	if (!PyArg_ParseTuple(args, "sfff", &fileName, &sx, &sy, &sz)) return NULL;

	int sts;

	sts = update_scale_of( fileName, sx, sy, sz );
	
	if (sts!=NO_ERROR) return NULL;
	
    return Py_None;

}

static PyObject *mrc_grid_parameters(PyObject *self, PyObject *args) {

    char *fileName;

	if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    int sts;
    int nxstart=0, nystart=0, nzstart=0;
    int mx=0, my=0, mz=0;

    sts = get_grid_origin_of(fileName, &nxstart, &nystart, &nzstart);

	if (sts!=NO_ERROR) return NULL;
	
    sts = get_grid_size_of(fileName, &mx, &my, &mz);

	if (sts!=NO_ERROR) return NULL;

    return Py_BuildValue( "(i,i,i,i,i,i)", nxstart, nystart, nzstart, mx, my, mz );

}


/*
 * Update the grid parameters.
 */
static PyObject *update_grid_parameters(PyObject *self, PyObject *args) {

    char *fileName;
 	int nxstart, nystart, nzstart, mx, my, mz;

	if (!PyArg_ParseTuple(args, "siiiiii", &fileName, &nxstart, &nystart, &nzstart, &mx, &my, &mz)) return NULL;

	int sts;

	sts = update_grid_origin_of( fileName, nxstart, nystart, nzstart );
	
	if (sts!=NO_ERROR) return NULL;
	
	sts = update_grid_size_of( fileName, mx, my, mz );
	
	if (sts!=NO_ERROR) return NULL;
	
    return Py_None;

}


static PyObject *mrc_coords(PyObject *self, PyObject *args) {

    char *fileName;

	if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    int sts;
    float xorg=0.0, yorg=0.0, zorg=0.0;
    float xlen=0.0, ylen=0.0, zlen=0.0;

    sts = get_origin_of(fileName, &xorg, &yorg, &zorg);

	if (sts!=NO_ERROR) return NULL;
	
    sts = get_length_of(fileName, &xlen, &ylen, &zlen);

	if (sts!=NO_ERROR) return NULL;

	return Py_BuildValue( "(f,f,f,f,f,f)", xorg, yorg, zorg, xlen, ylen, zlen );
	
}

/*
 * Update the coordinates of the MRC File.
 */
static PyObject *update_coords(PyObject *self, PyObject *args) {

    char *fileName;
    float xorg, yorg, zorg, xlen, ylen, zlen;

	if (!PyArg_ParseTuple(args, "sffffff", &fileName, &xorg, &yorg, &zorg, &xlen, &ylen, &zlen)) return NULL;

	int sts;

	sts = update_origin_of( fileName, xorg, yorg, zorg );
	
	if (sts!=NO_ERROR) return NULL;
	
	sts = update_length_of( fileName, xlen, ylen, zlen );
	
	if (sts!=NO_ERROR) return NULL;
	
    return Py_None;

}

/*
 * Get the origin.
 */
static PyObject *mrc_origin(PyObject *self, PyObject *args) {

    char *fileName;

	if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    int sts;
    float xorg=0.0, yorg=0.0, zorg=0.0;

    sts = get_origin_of(fileName, &xorg, &yorg, &zorg);

	if (sts!=NO_ERROR) return NULL;

    return Py_BuildValue( "(f,f,f)", xorg, yorg, zorg );

}


/*
 * Update the origin.
 */
static PyObject *update_origin(PyObject *self, PyObject *args) {

    char *fileName;
 	float xorg, yorg, zorg;

	if (!PyArg_ParseTuple(args, "sfff", &fileName, &xorg, &yorg, &zorg)) return NULL;

	int sts;

	sts = update_origin_of( fileName, xorg, yorg, zorg );
	
	if (sts!=NO_ERROR) return NULL;
	
    return Py_None;

}


static PyObject *mrc_length(PyObject *self, PyObject *args) {

    char *fileName;

	if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    int sts;
    float xlen=0.0, ylen=0.0, zlen=0.0;

    sts = get_length_of(fileName, &xlen, &ylen, &zlen);

	if (sts!=NO_ERROR) return NULL;

    return Py_BuildValue( "(f,f,f)", xlen, ylen, zlen );

}

/*
 * Update the length.
 */
static PyObject *update_length(PyObject *self, PyObject *args) {

    char *fileName;
 	float xlen, ylen, zlen;

	if (!PyArg_ParseTuple(args, "sfff", &fileName, &xlen, &ylen, &zlen)) return NULL;

	int sts;

	sts = update_length_of( fileName, xlen, ylen, zlen );
	
	if (sts!=NO_ERROR) return NULL;
	
    return Py_None;

}



static PyObject *mrc_description(PyObject *self, PyObject *args) {

    char *fileName;

    if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    int sts, i, nlabels;
    char *labels[MRC_NLABELS];

    for (i = 0; i < MRC_NLABELS; i++) {
        labels[i] = (char *) malloc((MRC_LABEL_SIZE + 1) * sizeof (char));
    }

    sts = get_description_of(fileName, &nlabels, &labels[0]);

    PyObject *description = PyList_New(nlabels);

    for (i = 0; i < nlabels; i++) {
        PyList_SetItem(description, i, PyString_FromString(labels[i]));
        free(labels[i]);
    }

    if (sts!=NO_ERROR) return NULL;

    return description;

}


static PyObject *update_description(PyObject *self, PyObject *args) {

    char *fileName;
    PyObject *description = NULL;

	if (!PyArg_ParseTuple(args, "sO", &fileName, &description)) return NULL;

	int sts,i,nlabels;
	char *labels[MRC_NLABELS];

	nlabels = PyList_GET_SIZE(description);
	
	for (i=0;i<nlabels;i++) 
	{
		labels[i] = PyString_AsString(PyList_GetItem(description,i));
	}	
	
	sts = update_description_of(fileName, nlabels, &labels[0]);
	
	if (sts!=NO_ERROR) return NULL;
	
    return Py_None;

}


static PyObject *copy_description(PyObject *self, PyObject *args) {

    char *input_fileName, *output_fileName;

	if (!PyArg_ParseTuple(args, "ss", &input_fileName, &output_fileName)) return NULL;

	int sts,i,nlabels;
	char *labels[MRC_NLABELS];

	for (i=0;i<MRC_NLABELS;i++) {
		labels[i]  = (char *)malloc((MRC_LABEL_SIZE+1)*sizeof(char));
	}

    sts = get_description_of(input_fileName, &nlabels, &labels[0]);

	if (sts!=NO_ERROR) return NULL;
	
	sts = update_description_of(output_fileName, nlabels, &labels[0]);
	
    for (i=0; i<nlabels; i++) 
    {
		free(labels[i]);
	}	
	
	if (sts!=NO_ERROR) return NULL;
	
    return Py_None;

}

/*
 *  Add a label to the description.
 */
static PyObject *add_label_description(PyObject *self, PyObject *args) {

    char *fileName, *description;

	if (!PyArg_ParseTuple(args, "ss", &fileName, &description)) return NULL;

	int sts = add_label(fileName, description);
	
	if (sts!=NO_ERROR) return NULL;
	
    return Py_None;

}


/*
 * Copy dimensions informations (mx,xlen,..) contained in a header from
 * one MRC file to another. Only informations on the two directions X and
 * Y are copied by default.
 */
static PyObject *copy_dimensions(PyObject *self, PyObject *args) {

    char *input_fileName, *output_fileName;
    PyObject *description = NULL;

	if (!PyArg_ParseTuple(args, "ss", &input_fileName, &output_fileName)) return NULL;

	int sts;
	// Sample dimensions
	float xorg, yorg, zorg;
	float xlen, ylen, zlen;
	// Grid to load
	int nxstart, nystart, nzstart;
	int mx, my, mz;
	
	/* Copy the sample dimensions */
	
    sts = get_origin_of(input_fileName, &xorg, &yorg, &zorg);

	if (sts!=NO_ERROR) return NULL;
	
	sts = update_origin_of(output_fileName, xorg, yorg, zorg);
	
	if (sts!=NO_ERROR) return NULL;

    sts = get_length_of(input_fileName, &xlen, &ylen, &zlen);

	if (sts!=NO_ERROR) return NULL;
	
	sts = update_length_of(output_fileName, xlen, ylen, zlen);
	
	if (sts!=NO_ERROR) return NULL;
	
	/* Copy the grid-to-load parameters */
	
    sts = get_grid_origin_of(input_fileName, &nxstart, &nystart, &nzstart);

	if (sts!=NO_ERROR) return NULL;
	
	sts = update_grid_origin_of(output_fileName, nxstart, nystart, nzstart);
	
	if (sts!=NO_ERROR) return NULL;
	
	sts = get_grid_size_of(input_fileName, &mx, &my, &mz);

	if (sts!=NO_ERROR) return NULL;
	
	sts = update_grid_size_of(output_fileName, mx, my, mz);
	
	if (sts!=NO_ERROR) return NULL;
	
    return Py_None;

}


static PyObject *mrc_stats(PyObject *self, PyObject *args) {

    char *fileName;

	if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    int sts;
    float m0=0, m1=0, m2=0, min=0, max=0;

    sts = get_statistics_of(fileName, &m0, &m1, &m2, &min, &max);

    if (sts!=NO_ERROR) return NULL;

    return Py_BuildValue( "(f,f,f,f,f)", m0, m1, m2, min, max );

}

static PyObject *mrc_tilt_angles(PyObject *self, PyObject *args) {

    char *fileName;
    int nx, ny, nangles;
    float *angles;

    if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;
    
    int sts = get_size_of( fileName, &nx, &ny, &nangles );

    if (sts!=NO_ERROR) return NULL;

    angles = (float *)malloc( nangles*sizeof(float) );

    sts = get_tilt_angles( fileName, nangles, angles);

    if (sts!=NO_ERROR) return NULL;

    if (nangles==0) return Py_None;

    PyObject *description = PyList_New(nangles);

    int i;

    for (i = 0; i < nangles; i++) {
        PyList_SetItem(description, i, PyFloat_FromDouble(angles[i]));
    }
    
    free(angles);

    return description;

}

static PyObject *update_header(PyObject *self, PyObject *args) {

    char *fileName;
    PyObject *stats;

    if (!PyArg_ParseTuple(args, "sO", &fileName, &stats)) return NULL;

    MrcHeader *input_file_header  = (MrcHeader *)malloc(sizeof(MrcHeader));
    FILE* input_file;

    if (!input_file_header) {
            txbr_error(stderr, "ERROR: read_slice - getting memory.\n");
            return NULL;
    }

    if ((input_file = fopen(fileName,"r+"))==NULL) {
            txbr_error(stderr, "ERROR: read_slice - cannot open file %s.\n", fileName);
            return NULL;
    }

    mrc_head_read(input_file,input_file_header);

    int sts=NO_ERROR;

    float m0=0, m1=0, m2=0, min=0, max=0;

    if (PyDict_Check(stats)) {
            m0 =  (float)PyFloat_AsDouble(PyDict_GetItemString(stats,"m0"));
            m1 =  (float)PyFloat_AsDouble(PyDict_GetItemString(stats,"mean"));
            m2 =  (float)PyFloat_AsDouble(PyDict_GetItemString(stats,"m2"));
            min =  (float)PyFloat_AsDouble(PyDict_GetItemString(stats,"min"));
            max =  (float)PyFloat_AsDouble(PyDict_GetItemString(stats,"max"));
    } else {
            sts = get_statistics_of(fileName, &m0, &m1, &m2, &min, &max);
    }

    int digits = 2;
    double factor = pow(10.0,digits-ceil(log10(fabs(max-min))));

    input_file_header->amean = m1;
    input_file_header->amin = floor(min*factor)/factor;
    input_file_header->amax = ceil(max*factor)/factor;

    mrc_head_write(input_file,input_file_header);

    if (fclose(input_file)) {
            txbr_error(stderr, "ERROR: write_slice - closing file %s.\n", &input_file[0]);
            return NULL;
    }

    free(input_file_header);

    if (sts!=NO_ERROR) return NULL;

    return Py_BuildValue( "(f,f,f)", min, m0, max );

}


static PyObject *create_header(PyObject *self, PyObject *args) {

    char *fileName;
    int nx=0, ny=0, nz=0;
    int mode = MRC_MODE_FLOAT;
    float xorg=0.0, yorg=0.0, zorg=0.0;

	if (!PyArg_ParseTuple(args, "siiifffi", &fileName, &nx, &ny, &nz, &xorg, &yorg, &zorg, &mode)) return NULL;

	FILE* input_file;

	if ((input_file = fopen(fileName,"w+"))==NULL) {
		txbr_error(stderr, "ERROR: create_header - cannot open file %s.\n", fileName);
		return NULL;
	}

	MrcHeader *input_file_header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!input_file_header) {
		txbr_error(stderr, "ERROR: create_header - getting memory.\n");
		return NULL;
	}

	mrc_head_new(input_file_header, nx, ny, nz, mode);
	
	input_file_header->xorg = xorg;
	input_file_header->yorg = yorg;
	input_file_header->zorg = zorg;

	mrc_head_write(input_file, input_file_header);

	if (fclose(input_file)) {
		txbr_error(stderr, "ERROR: create_header - closing file %s.\n", &input_file[0]);
		return NULL;
	}

	free(input_file_header);

	return Py_BuildValue( "i", 1 );

}


static PyObject *create_header_from(PyObject *self, PyObject *args) {

    char *fileName, *tplFileName;

	if (!PyArg_ParseTuple(args, "ss", &fileName, &tplFileName)) return NULL;

	MrcHeader *tpl_file_header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!tpl_file_header) {
		txbr_error(stderr, "ERROR: create_header_from - getting memory.\n");
		return NULL;
	}

	FILE* template_file;

	if ((template_file = fopen(tplFileName,"r"))==NULL) {
		txbr_error(stderr, "ERROR: create_header_from - cannot open template file %s.\n", tplFileName);
		return NULL;
	}

	mrc_head_read(template_file, tpl_file_header);

        if (fclose(template_file)) {
		txbr_error(stderr, "ERROR: create_header_from - closing template file %s.\n", tplFileName);
		return NULL;
	}

	FILE* file;

	if ((file = fopen(fileName,"w+"))==NULL) {
		txbr_error(stderr, "ERROR: create_header_from - cannot open new file %s.\n", fileName);
		return NULL;
	}

	mrc_head_write(file, tpl_file_header);
	mrc_data_new(file, tpl_file_header);

	if (fclose(file)) {
		txbr_error(stderr, "ERROR: create_header_from - closing new file %s.\n", fileName);
		return NULL;
	}

	free(tpl_file_header);

	return Py_BuildValue( "i", 1 );

}


static PyObject *read_slice(PyObject *self, PyObject *args) {
	
	typedef double ARRAY_TYPE;

    char *fileName, *axis;
    int slice = 0;

    if (!PyArg_ParseTuple(args, "ssi", &fileName, &axis, &slice)) return NULL;

    MrcHeader *input_file_header  = (MrcHeader *)malloc(sizeof(MrcHeader));
    FILE* input_file;

    if (!input_file_header) {
            txbr_error(stderr, "ERROR: read_slice - getting memory.\n");
            return NULL;
    }

    if ((input_file = fopen(fileName,"rbe"))==NULL) {
            txbr_error(stderr, "ERROR: read_slice - cannot open file %s.\n", fileName);
            return NULL;
    }

    mrc_head_read(input_file,input_file_header);

    int nx = input_file_header->nx;
    int ny = input_file_header->ny;
    int nz = input_file_header->nz;

    int n1,n2;

    switch (axis[0]) {

            case 'x':
            case 'X':
                    n1 = ny;
                    n2 = nz;
                    if (slice >= nx) return NULL;
            break;

            case 'y':
            case 'Y':
                    n1 = nx;
                    n2 = nz;
                    if (slice >= ny) return NULL;
            break;

            case 'z':
            case 'Z':
                    n1 = nx;
                    n2 = ny;
                    if (slice >= nz) return NULL;
            break;

            default:
                    txbr_error(stderr, "ERROR: mrc_slice - axis error.\n");
                    return NULL;

    }

    PyArrayObject *array = NULL;

    int i,j,k,index = 0;

//	char* c_array_buffer = NULL;
    unsigned char* uc_array_buffer = NULL;
    short* s_array_buffer = NULL;
    float* f_array_buffer = NULL;
    unsigned short* us_array_buffer = NULL;

    if (input_file_header->mode==SLICE_MODE_BYTE) {

            npy_intp dimensions[2] = { n1, n2 };
            array = PyArray_SimpleNew(2, dimensions, PyArray_DOUBLE);
            char *data = (char *)PyArray_DATA(array);

            uc_array_buffer = (unsigned char*)mrc_mread_slice(input_file, input_file_header, slice, axis[0]);

            for (j=0; j<dimensions[1]; j++) {
                    for (i=0; i<dimensions[0]; i++) {
                            *(double *)(data + i*array->strides[0] + j*array->strides[1]) = uc_array_buffer[index++];
                    }
            }

    }

    if (input_file_header->mode==SLICE_MODE_SHORT) {

            npy_intp dimensions[2] = { n1, n2 };
            array = PyArray_SimpleNew(2, dimensions, PyArray_DOUBLE);
            char *data = (char *)PyArray_DATA(array);

            s_array_buffer = (short*)mrc_mread_slice(input_file, input_file_header, slice, axis[0]);

            for (j=0; j<dimensions[1]; j++) {
                    for (i=0; i<dimensions[0]; i++) {
                            *(double *)(data + i*array->strides[0] + j*array->strides[1]) = s_array_buffer[index++];
                    }
            }

    }

    if (input_file_header->mode==SLICE_MODE_FLOAT) {

            npy_intp dimensions[2] = { n1, n2 };
            array = PyArray_SimpleNew(2, dimensions, PyArray_DOUBLE);
            char *data = (char *)PyArray_DATA(array);

            f_array_buffer = (float*)mrc_mread_slice(input_file, input_file_header, slice, axis[0]);

            for (j=0; j<dimensions[1]; j++) {
                    for (i=0; i<dimensions[0]; i++) {
                            *(double *)(data + i*array->strides[0] + j*array->strides[1]) = f_array_buffer[index++];
                    }
            }
            /*
            npy_intp dimensions[2] = { n1, n2 };
            array = PyArray_SimpleNew(2, dimensions, PyArray_FLOAT);
            char *data = (char *)PyArray_DATA(array);

            f_array_buffer = (float*)mrc_mread_slice(input_file, input_file_header, slice, axis[0]);

            for (j=0; j<dimensions[1]; j++) {
                    for (i=0; i<dimensions[0]; i++) {
                            *(float *)(data + i*array->strides[0] + j*array->strides[1]) = f_array_buffer[index++];
                    }
            }
    */
    }

    if (input_file_header->mode==SLICE_MODE_COMPLEX_SHORT) {

            npy_intp dimensions[2] = { n1, n2 };
            array = PyArray_SimpleNew(2, dimensions, PyArray_CDOUBLE);
            char *data = (char *)PyArray_DATA(array);

            s_array_buffer = (short*)mrc_mread_slice(input_file, input_file_header, slice, axis[0]);

            for (j=0; j<dimensions[1]; j++) {
                    for (i=0; i<dimensions[0]; i++) {
                            *(double *)(data + i*array->strides[0] + j*array->strides[1]) = s_array_buffer[index++];
                            *(double *)(data + i*array->strides[0] + j*array->strides[1] + sizeof(double)) = s_array_buffer[index++];
                    }
            }

    }

    if (input_file_header->mode==SLICE_MODE_COMPLEX_FLOAT) {

            npy_intp dimensions[2] = { n1, n2 };
            array = PyArray_SimpleNew(2, dimensions, PyArray_CDOUBLE);
            char *data = (char *)PyArray_DATA(array);

            f_array_buffer = (float*)mrc_mread_slice(input_file, input_file_header, slice, axis[0]);

            for (j=0; j<dimensions[1]; j++) {
                    for (i=0; i<dimensions[0]; i++) {
                            *(double *)(data + i*array->strides[0] + j*array->strides[1]) = f_array_buffer[index++];
                            *(double *)(data + i*array->strides[0] + j*array->strides[1] + sizeof(double)) = f_array_buffer[index++];
                    }
            }

    }

    if (input_file_header->mode==SLICE_MODE_USHORT) {

            npy_intp dimensions[2] = { n1, n2 };
            array = PyArray_SimpleNew(2, dimensions, PyArray_DOUBLE);
            char *data = (char *)PyArray_DATA(array);

            us_array_buffer = (unsigned short*)mrc_mread_slice(input_file, input_file_header, slice, axis[0]);

            for (j=0; j<dimensions[1]; j++) {
                    for (i=0; i<dimensions[0]; i++) {
                            *(double *)(data + i*array->strides[0] + j*array->strides[1]) = us_array_buffer[index++];
                    }
            }

    }

    if (input_file_header->mode==SLICE_MODE_RGB) {

            npy_intp dimensions[3] = { n1, n2, 3 };
            array = PyArray_SimpleNew(3, dimensions, PyArray_UBYTE);
            char *data = (char *)PyArray_DATA(array);

            uc_array_buffer = (unsigned char*)mrc_mread_slice(input_file, input_file_header, slice, axis[0]);

            for (j=0; j<dimensions[1]; j++) {
                    for (i=0; i<dimensions[0]; i++) {
                            for (k=0; k<dimensions[2]; k++) {
                                    *(unsigned char *)(data + i*array->strides[0] + j*array->strides[1] + k*array->strides[2]) = uc_array_buffer[index++];
                            }
                    }
            }

    }

    if (fclose(input_file)) {
            txbr_error(stderr, "ERROR: read_slice - closing file %s.\n", &input_file[0]);
            return NULL;
    }

    if (input_file_header->mode==SLICE_MODE_BYTE) free(uc_array_buffer);
    if (input_file_header->mode==SLICE_MODE_SHORT) free(s_array_buffer);
    if (input_file_header->mode==SLICE_MODE_FLOAT) free(f_array_buffer);
    if (input_file_header->mode==SLICE_MODE_COMPLEX_SHORT) free(f_array_buffer);
    if (input_file_header->mode==SLICE_MODE_COMPLEX_FLOAT) free(f_array_buffer);
    if (input_file_header->mode==SLICE_MODE_USHORT) free(us_array_buffer);
    if (input_file_header->mode==SLICE_MODE_RGB) free(uc_array_buffer);

    free(input_file_header);

    if (array==NULL) return NULL;

    return PyArray_Return(array);

}


static PyObject *write_slice(PyObject *self, PyObject *args) {

    char *fileName, *axis;
    int slice = 0;

    PyArrayObject *array;

    if (!PyArg_ParseTuple(args, "ssiO", &fileName, &axis, &slice, &array)) return NULL;

//    if (array->nd!=2) {
//    	PyErr_SetString(PyExc_ValueError, "Array must be two-dimensional!");
//    	return NULL;
//    }
//
//    if ( array->descr->type_num!=PyArray_FLOAT && array->descr->type_num!=PyArray_DOUBLE) {
//    	PyErr_SetString(PyExc_ValueError, "Array must be of type float");
//    	return NULL;
//    }

	MrcHeader *input_file_header  = (MrcHeader *)malloc(sizeof(MrcHeader));
	FILE* input_file;

	if (!input_file_header) {
		txbr_error(stderr, "ERROR: write_slice - getting memory.\n");
		return NULL;
	}

	if ((input_file = fopen(fileName,"r+"))==NULL) {
		txbr_error(stderr, "ERROR: openMRCFile - cannot open file %s.\n", fileName);
		return NULL;
	}

	mrc_head_read(input_file,input_file_header);

	int n1,n2,length = 0;
	int nx = input_file_header->nx;
	int ny = input_file_header->ny;
	int nz = input_file_header->nz;
	int mode = input_file_header->mode;

	switch (axis[0]) {

		case 'x':
		case 'X':
			n1 = ny;
			n2 = nz;
			if (slice >= nx) return NULL;
		break;

		case 'y':
		case 'Y':
			n1 = nx;
			n2 = nz;
			if (slice >= ny) return NULL;
		break;

		case 'z':
		case 'Z':
			n1 = nx;
			n2 = ny;
			if (slice >= nz) return NULL;
		break;

		default:
			txbr_error(stderr, "ERROR: write_slice - axis error.\n");
			return NULL;

	}

	length = n1*n2;

	int i,j,k;

	unsigned char* uc_array_buffer = NULL;
	char* c_array_buffer = NULL;
	short* s_array_buffer = NULL;
	float* f_array_buffer = NULL;
	unsigned short* us_array_buffer = NULL;

	if (mode==SLICE_MODE_BYTE) {

		uc_array_buffer  = (float *)malloc(length*sizeof(short));

		if (!uc_array_buffer) {
			txbr_error(stderr, "ERROR: write_slice - getting memory.\n");
			return NULL;
		}

		if ( array->descr->type_num==PyArray_FLOAT ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					uc_array_buffer[index++] = *(unsigned char *)(array->data + i*array->strides[0] + j*array->strides[1]);
				}
			}
		}

		if ( array->descr->type_num==PyArray_DOUBLE ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					uc_array_buffer[index++] = (unsigned char)*(double *)(array->data + i*array->strides[0] + j*array->strides[1]);
				}
			}
		}

		mrc_write_slice(&uc_array_buffer[0], input_file, input_file_header, slice, axis[0]);

	}


	if (mode==SLICE_MODE_SHORT) {

		s_array_buffer  = (short *)malloc(length*sizeof(short));

		if (!s_array_buffer) {
			txbr_error(stderr, "ERROR: write_slice - getting memory.\n");
			return NULL;
		}

		if ( array->descr->type_num==PyArray_FLOAT ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					s_array_buffer[index++] = *(short *)(array->data + i*array->strides[0] + j*array->strides[1]);
				}
			}
		}

		if ( array->descr->type_num==PyArray_DOUBLE ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					s_array_buffer[index++] = (short)*(double *)(array->data + i*array->strides[0] + j*array->strides[1]);
				}
			}
		}

		mrc_write_slice(&s_array_buffer[0], input_file, input_file_header, slice, axis[0]);

	}

	if (mode==SLICE_MODE_FLOAT) {

		f_array_buffer  = (float *)malloc(length*sizeof(float));

		if (!f_array_buffer) {
			txbr_error(stderr, "ERROR: write_slice - getting memory.\n");
			return NULL;
		}

		if ( array->descr->type_num==PyArray_FLOAT ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					f_array_buffer[index++] = *(float *)(array->data + i*array->strides[0] + j*array->strides[1]);
				}
			}
		}

		if ( array->descr->type_num==PyArray_DOUBLE ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					f_array_buffer[index++] = (float)*(double *)(array->data + i*array->strides[0] + j*array->strides[1]);
				}
			}
		}

		mrc_write_slice(&f_array_buffer[0], input_file, input_file_header, slice, axis[0]);

	}

	if (mode==SLICE_MODE_COMPLEX_SHORT) {

		s_array_buffer  = (float *)malloc(length*sizeof(short));

		if (!s_array_buffer) {
			txbr_error(stderr, "ERROR: write_slice - getting memory.\n");
			return NULL;
		}

		if ( array->descr->type_num==PyArray_CFLOAT ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					s_array_buffer[index++] = *(short *)(array->data + i*array->strides[0] + j*array->strides[1]);
					s_array_buffer[index++] = *(short *)(array->data + i*array->strides[0] + j*array->strides[1] + sizeof(float));
				}
			}
		}

		if ( array->descr->type_num==PyArray_CDOUBLE ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					s_array_buffer[index++] = (short)*(double *)(array->data + i*array->strides[0] + j*array->strides[1]);
					s_array_buffer[index++] = (short)*(double *)(array->data + i*array->strides[0] + j*array->strides[1] + sizeof(double));
				}
			}
		}

		mrc_write_slice(&s_array_buffer[0], input_file, input_file_header, slice, axis[0]);

	}

	if (mode==SLICE_MODE_COMPLEX_FLOAT) {

		f_array_buffer  = (float *)malloc(2*length*sizeof(float));

		if (!f_array_buffer) {
			txbr_error(stderr, "ERROR: write_slice - getting memory.\n");
			return NULL;
		}


		if ( array->descr->type_num==PyArray_CFLOAT ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					f_array_buffer[index++] = *(float *)(array->data + i*array->strides[0] + j*array->strides[1]);
					f_array_buffer[index++] = *(float *)(array->data + i*array->strides[0] + j*array->strides[1] + sizeof(float));
				}
			}
		}

		if ( array->descr->type_num==PyArray_CDOUBLE ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					f_array_buffer[index++] = (float)*(double *)(array->data + i*array->strides[0] + j*array->strides[1]);
					f_array_buffer[index++] = (float)*(double *)(array->data + i*array->strides[0] + j*array->strides[1] + sizeof(double));
				}
			}
		}

		mrc_write_slice(&f_array_buffer[0], input_file, input_file_header, slice, axis[0]);

	}

	if (mode==SLICE_MODE_USHORT) {

		us_array_buffer = (unsigned short*)mrc_mread_slice(input_file, input_file_header, slice, axis[0]);

		if (!us_array_buffer) {
			txbr_error(stderr, "ERROR: write_slice - getting memory.\n");
			return NULL;
		}

		if ( array->descr->type_num==PyArray_FLOAT ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					us_array_buffer[index++] = *(unsigned short *)(array->data + i*array->strides[0] + j*array->strides[1]);
				}
			}
		}

		if ( array->descr->type_num==PyArray_DOUBLE ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					us_array_buffer[index++] = (unsigned char)*(double *)(array->data + i*array->strides[0] + j*array->strides[1]);
				}
			}
		}

		mrc_write_slice(&us_array_buffer[0], input_file, input_file_header, slice, axis[0]);

	}

	if (mode==SLICE_MODE_RGB) {

		uc_array_buffer  = (unsigned char *)malloc(3*length*sizeof(unsigned char));

		if (!uc_array_buffer) {
			txbr_error(stderr, "ERROR: write_slice - getting memory.\n");
			return NULL;
		}

		if ( array->descr->type_num==PyArray_UBYTE ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					for (k=0; k<array->dimensions[2]; k++) {
						uc_array_buffer[index++] = *(unsigned char *)(array->data + i*array->strides[0] + j*array->strides[1] + k*array->strides[2]);
					}
				}
			}
		}

		if ( array->descr->type_num==PyArray_FLOAT ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					for (k=0; k<array->dimensions[2]; k++) {
						uc_array_buffer[index++] = (unsigned char)*(float *)(array->data + i*array->strides[0] + j*array->strides[1] + k*array->strides[2]);
					}
				}
			}
		}

		if ( array->descr->type_num==PyArray_DOUBLE ) {
			int index = 0;
			for (j=0; j<array->dimensions[1]; j++) {
				for (i=0; i<array->dimensions[0]; i++) {
					for (k=0; k<array->dimensions[2]; k++) {
						uc_array_buffer[index++] = (unsigned char)*(double *)(array->data + i*array->strides[0] + j*array->strides[1] + k*array->strides[2]);
					}
				}
			}
		}

		mrc_write_slice(&uc_array_buffer[0], input_file, input_file_header, slice, axis[0]);

	}


	if (fclose(input_file)) {
		txbr_error(stderr, "ERROR: write_slice - closing file %s.\n", &input_file[0]);
		return NULL;
	}

	if (uc_array_buffer) free(uc_array_buffer);
	if (c_array_buffer) free(c_array_buffer);
	if (s_array_buffer) free(s_array_buffer);
	if (f_array_buffer) free(f_array_buffer);
	if (us_array_buffer) free(us_array_buffer);

	free(input_file_header);

	return Py_BuildValue( "i", 1 );

}


static PyMethodDef MRCMethods[] = {
	{ "header", mrc_header, METH_VARARGS, "Display the header info." },
	{ "mode", mrc_mode, METH_VARARGS, "Get the MRC File mode." },
	{ "size", mrc_size, METH_VARARGS, "Give the 3D size of an MRC file." },
	{ "scale", mrc_scale, METH_VARARGS, "Give the 3D scale of an MRC file." },
	{ "update_scale", update_scale, METH_VARARGS, "Update the 3D scale of an MRC file." },
	{ "grid", mrc_grid_parameters, METH_VARARGS, "Get the grid parameters from the area to load in an MRC file." },
	{ "update_grid", update_grid_parameters, METH_VARARGS, "Update the grid parameters from the area to load in an MRC file." },
	{ "coords", mrc_coords, METH_VARARGS, "Give coordinates corresponding to the grid of an MRC file." },
	{ "update_coords", update_coords, METH_VARARGS, "Update the coordinates corresponding to the grid of an MRC file." },
	{ "origin", mrc_origin, METH_VARARGS, "Give coordinates of the origin point in an MRC file." },
	{ "update_origin", update_origin, METH_VARARGS, "Update the coordinates of the origin point in an MRC file." },
	{ "length", mrc_length, METH_VARARGS, "Give the lengths for this MRC file." },
	{ "update_length", update_length, METH_VARARGS, "Update the lengths for this MRC file." },
	{ "description", mrc_description, METH_VARARGS, "Get the MRC File description." },
	{ "update_description", update_description, METH_VARARGS, "Update the MRC File description." },
        { "tilts", mrc_tilt_angles, METH_VARARGS, "Load the tilt values from the ader." },
	{ "copy_description", copy_description, METH_VARARGS, "Copy a MRC File description to another file." },
	{ "add_label", add_label_description, METH_VARARGS, "Add a label to a MRC File description." },
	{ "copy_dimensions", copy_dimensions, METH_VARARGS, "Copy a MRC File dimensions to another file." },
	{ "stats", mrc_stats, METH_VARARGS, "Calculate Statistics on an MRC File." },
	{ "update_header", update_header, METH_VARARGS, "Write the an MRC File header." },
	{ "write_header", create_header, METH_VARARGS, "Write the an MRC File header." },
	{ "create_header_from", create_header_from, METH_VARARGS, "Create an MRC File header from another MRCFile." },
	{ "read_slice", read_slice, METH_VARARGS, "Get a slice from an MRC File." },
	{ "write_slice", write_slice, METH_VARARGS, "Write a slice to an MRC File." },
	{NULL, NULL, 0, NULL} /* Sentinel */
};

PyMODINIT_FUNC initmrcfile(void) {

    PyObject *m;

    m = Py_InitModule("mrcfile", MRCMethods);

    if (m != NULL) {

        PyObject *mode_byte;
        PyObject *mode_short;
        PyObject *mode_float;
        PyObject *mode_complex_short;
        PyObject *mode_complex_float;
        PyObject *mode_ushort;
        PyObject *mode_rgb;
        
        mode_byte = PyInt_FromSize_t(MRC_MODE_BYTE);
        Py_INCREF(mode_byte);
        PyModule_AddObject(m, "MODE_BYTE", mode_byte);

        mode_short = PyInt_FromSize_t(MRC_MODE_SHORT);
        Py_INCREF(mode_short);
        PyModule_AddObject(m, "MODE_SHORT", mode_short);

        mode_float = PyInt_FromSize_t(MRC_MODE_FLOAT);
        Py_INCREF(mode_float);
        PyModule_AddObject(m, "MODE_FLOAT", mode_float);

        mode_complex_short = PyInt_FromSize_t(MRC_MODE_COMPLEX_SHORT);
        Py_INCREF(mode_complex_short);
        PyModule_AddObject(m, "MODE_COMPLEX_SHORT", mode_complex_short);

        mode_complex_float = PyInt_FromSize_t(MRC_MODE_COMPLEX_FLOAT);
        Py_INCREF(mode_complex_float);
        PyModule_AddObject(m, "MODE_COMPLEX_FLOAT", mode_complex_float);

        mode_ushort = PyInt_FromSize_t(MRC_MODE_USHORT);
        Py_INCREF(mode_ushort);
        PyModule_AddObject(m, "MODE_USHORT", mode_ushort);

        mode_rgb = PyInt_FromSize_t(MRC_MODE_RGB);
        Py_INCREF(mode_rgb);
        PyModule_AddObject(m, "MODE_RGB", mode_rgb);

    }

    import_array();


}


void main(int argc, char *argv[]){

	Py_SetProgramName(argv[0]);	/* Pass argv[0] to the Python interpreter */

	Py_Initialize(); /* Initialize the Python interpreter.	*/

	initmrcfile();	/* Add a static module */

}
