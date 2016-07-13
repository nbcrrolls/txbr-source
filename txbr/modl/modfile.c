#include <Python.h>

#include <numpy/arrayobject.h>
#include "imodel.h"


static PyObject *getXYZMax(PyObject *self, PyObject *args) {

    char *fileName;

    if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    Imod *imodel = imodRead(fileName);

    int xmax = imodel->xmax;
    int ymax = imodel->ymax;
    int zmax = imodel->zmax;

    imodDelete(imodel);

    return Py_BuildValue("(i,i,i)",xmax,ymax,zmax);

}

static PyObject *getDimensions(PyObject *self, PyObject *args) {

    char *fileName;

    if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    Imod *imodel = imodRead(fileName);

    int xlen = (int)(imodel->xmax - imodel->xoffset);
    int ylen = (int)(imodel->ymax - imodel->yoffset);
    int zlen = (int)(imodel->zmax - imodel->zoffset);

    float sx = imodel->refImage->cscale.x;
    float sy = imodel->refImage->cscale.y;
    float sz = imodel->refImage->cscale.z;

    int nx = (int)rint(xlen/sx);
    int ny = (int)rint(ylen/sy);
    int nz = (int)rint(zlen/sz);

    imodDelete(imodel);

    return Py_BuildValue("(i,i,i)",nx,ny,nz);

}

static PyObject *getScales(PyObject *self, PyObject *args) {

    char *fileName;

    if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

    Imod *imodel = imodRead(fileName);

    float sx = imodel->refImage->cscale.x;
    float sy = imodel->refImage->cscale.y;
    float sz = imodel->refImage->cscale.z;

    imodDelete(imodel);

    return Py_BuildValue("(f,f,f)",sx,sy,sz);

}


static PyObject *getImgRefParameters(PyObject *self, PyObject *args) {

	const char *filename;

	if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;

	Imod *imod_model = imodRead(filename);

	double tx=0.0, ty=0.0, tz=0.0, sx=1.0, sy=1.0, sz=1.0, rx=0.0, ry=0.0, rz=0.0;

	if (imod_model->refImage) {

		tx = (double)imod_model->refImage->ctrans.x;
		ty = (double)imod_model->refImage->ctrans.y;
		tz = (double)imod_model->refImage->ctrans.z;

		sx = (double)imod_model->refImage->cscale.x;
		sy = (double)imod_model->refImage->cscale.y;
		sz = (double)imod_model->refImage->cscale.z;
		
		rx = (double)imod_model->refImage->crot.x;
		ry = (double)imod_model->refImage->crot.y;
		rz = (double)imod_model->refImage->crot.z;

	}

	return Py_BuildValue("(f,f,f,f,f,f,f,f,f)",tx,ty,tz,sx,sy,sz,rx,ry,rz);

}


static PyObject *numberOfObjects(PyObject *self, PyObject *args) {

    char *fileName;

	if (!PyArg_ParseTuple(args, "s", &fileName)) return NULL;

	Imod *imodel = imodRead(fileName);

	int number_of_objects = imodGetMaxObject(imodel);

	imodDelete(imodel);

    return Py_BuildValue( "i", number_of_objects );

}


static PyObject *numberOfContours(PyObject *self, PyObject *args) {

    char *fileName;
    int indexOfObject = 0;

	if (!PyArg_ParseTuple(args, "si", &fileName, &indexOfObject)) return NULL;

	Imod *imodel = imodRead(fileName);

	imodSetIndex(imodel, indexOfObject, -1, -1);

	Iobj *object = imodObjectGet(imodel);

	int numberOfContours =  imodObjectGetMaxContour(object);

	imodDelete(imodel);

    return Py_BuildValue( "i", numberOfContours );

}


static PyObject *numberOfPoints(PyObject *self, PyObject *args) {

    char *fileName;
    int indexOfObject = -1, indexOfContour=-1;

	if (!PyArg_ParseTuple(args, "sii", &fileName, &indexOfObject, &indexOfContour)) return NULL;

	Imod *imodel = imodRead(fileName);

	imodSetIndex(imodel, indexOfObject, indexOfContour, -1);

	Icont *contour = imodContourGet(imodel);

	int numberOfPoints = imodContourGetMaxPoint(contour);

	imodDelete(imodel);

    return Py_BuildValue( "i", numberOfPoints );

}


static PyObject *loadPoints(PyObject *self, PyObject *args) {

    char *fileName;
    int indexOfObject = -1, indexOfContour=-1;

    if (!PyArg_ParseTuple(args, "sii", &fileName, &indexOfObject, &indexOfContour)) return NULL;

	Imod *imodel = imodRead(fileName);

	imodSetIndex(imodel, indexOfObject, indexOfContour, -1);

	Icont *contour = imodContourGet(imodel);

	int numberOfPoints = imodContourGetMaxPoint(contour);
	Ipoint *points = imodContourGetPoints(contour);

	npy_intp dimensions[2] = { numberOfPoints, 3 };

	PyArrayObject *array = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, PyArray_DOUBLE);

	char *data = (char *)PyArray_DATA(array);

	int i,index = 0;

	for (i=0; i<numberOfPoints; i++) {
		*(double *)(data + i*array->strides[0] + 0*array->strides[1]) = points[index].x;
		*(double *)(data + i*array->strides[0] + 1*array->strides[1]) = points[index].y;
		*(double *)(data + i*array->strides[0] + 2*array->strides[1]) = points[index].z;
		index++;
	}

	imodDelete(imodel);

    return PyArray_Return(array);

}


static PyObject *loadModel(PyObject *self, PyObject *args) {

	PyObject *model;
	const char *filename;

	if (!PyArg_ParseTuple(args,"Os",&model,&filename)) return NULL;

	PyObject *object = NULL, *contour = NULL, *name = NULL;

	Imod *imod_model = imodRead(filename);
	Iobj *imod_obj;
	Icont *imod_cont;
	Ipoint *pts;

	int number_of_objects, numberOfContours, numberOfPoints;
	int iobj,icont;

	number_of_objects = imodGetMaxObject(imod_model);

	PyObject *xmax = PyInt_FromLong((long)imod_model->xmax);
	PyObject_SetAttrString(model, "xmax", xmax);
	Py_DECREF(xmax);

	PyObject *ymax = PyInt_FromLong((long)imod_model->ymax);
	PyObject_SetAttrString(model, "ymax", ymax);
	Py_DECREF(ymax);

	PyObject *zmax = PyInt_FromLong((long)imod_model->zmax);
	PyObject_SetAttrString(model, "zmax", zmax);
	Py_DECREF(zmax);

	PyObject *xoffset = PyInt_FromLong((long)imod_model->xoffset);
	PyObject_SetAttrString(model, "xoffset", xoffset);
	Py_DECREF(xoffset);

	PyObject *yoffset = PyInt_FromLong((long)imod_model->yoffset);
	PyObject_SetAttrString(model, "yoffset", yoffset);
	Py_DECREF(yoffset);

	PyObject *zoffset = PyInt_FromLong((long)imod_model->zoffset);
	PyObject_SetAttrString(model, "zoffset", zoffset);
	Py_DECREF(zoffset);


	// Image reference values

	PyObject *imgref = PyObject_GetAttrString(model, "imgref");

	if (imod_model->refImage) {

		PyObject *number;

		number = PyFloat_FromDouble((double)imod_model->refImage->oscale.x);
		PyObject_SetAttrString(imgref, "xscale_old", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->oscale.y);
		PyObject_SetAttrString(imgref, "yscale_old", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->oscale.z);
		PyObject_SetAttrString(imgref, "zscale_old", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->cscale.x);
		PyObject_SetAttrString(imgref, "xscale_new", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->cscale.y);
		PyObject_SetAttrString(imgref, "yscale_new", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->cscale.z);
		PyObject_SetAttrString(imgref, "zscale_new", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->otrans.x);
		PyObject_SetAttrString(imgref, "xtranslation_old", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->otrans.y);
		PyObject_SetAttrString(imgref, "ytranslation_old", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->otrans.z);
		PyObject_SetAttrString(imgref, "ztranslation_old", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->ctrans.x);
		PyObject_SetAttrString(imgref, "xtranslation_new", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->ctrans.y);
		PyObject_SetAttrString(imgref, "ytranslation_new", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->ctrans.z);
		PyObject_SetAttrString(imgref, "ztranslation_new", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->orot.x);
		PyObject_SetAttrString(imgref, "xrotation_old", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->orot.y);
		PyObject_SetAttrString(imgref, "yrotation_old", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->orot.z);
		PyObject_SetAttrString(imgref, "zrotation_old", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->crot.x);
		PyObject_SetAttrString(imgref, "xrotation_new", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->crot.y);
		PyObject_SetAttrString(imgref, "yrotation_new", number);
		Py_DECREF(number);

		number = PyFloat_FromDouble((double)imod_model->refImage->crot.z);
		PyObject_SetAttrString(imgref, "zrotation_new", number);
		Py_DECREF(number);

	}

	Py_DECREF(imgref);
	imgref = NULL;

	for (iobj=0; iobj<number_of_objects; iobj++) {

		imodSetIndex(imod_model, iobj, -1, -1);

		object = PyObject_CallMethod(model,"addNewObject",NULL);
		if (object==NULL) return NULL;

		imod_obj = imodObjectGet(imod_model);

		// Do not sort. Can lead to problems of mismatch when dealing with multiple model file.
		//imodObjectSort(imod_obj);

		numberOfContours =  imodObjectGetMaxContour(imod_obj);

		name = PyString_FromString((char*)imod_obj->name);
		PyObject_SetAttrString(object, "name", name);
		Py_DECREF(name);

		for (icont=0; icont<numberOfContours; icont++) {

			//imodSetIndex(imod_model, iobj, icont, -1);

			contour = PyObject_CallMethod(object,"addNewContour",NULL);

			if (contour==NULL) {
				Py_DECREF(object);
				return NULL;
			}

			//imod_cont = imodContourGet(imod_model);
			imod_cont = imodContourGetNext(imod_model);
			numberOfPoints = imodContourGetMaxPoint(imod_cont);

			//printf("contour: %i n=%i\n",imod_cont->surf,numberOfPoints);
			
			if (imod_cont->label) {
				name = PyString_FromString((char*)imod_cont->label[0].name);
				PyObject_SetAttrString(contour, "label", name);
				Py_DECREF(name);
			}

			npy_intp dimensions[2] = { numberOfPoints, 3 };

			PyArrayObject *array = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, PyArray_DOUBLE);

			if (array==NULL) {
				Py_DECREF(contour);
				Py_DECREF(object);
				return NULL;
			}

			PyObject_SetAttrString(contour,"points",(PyObject *)array);

			pts = imodContourGetPoints(imod_cont);

			char *data = (char *)PyArray_DATA(array);

			int i,index = 0;
			//printf("%i	%i	%i	%f\n",icont,imod_cont->surf,imod_model->cindex.contour,pts[index].x);

			for (i=0; i<numberOfPoints; i++) {
				*(double *)(data + i*array->strides[0] + 0*array->strides[1]) = pts[index].x;
				*(double *)(data + i*array->strides[0] + 1*array->strides[1]) = pts[index].y;
				*(double *)(data + i*array->strides[0] + 2*array->strides[1]) = pts[index].z;
				index++;
			}

			Py_DECREF(contour);
			Py_DECREF(array);

		}

		Py_DECREF(object);

	}

	imodDelete(imod_model);

	return Py_BuildValue( "i", 1 );

}


static PyObject *saveModel(PyObject *self, PyObject *args) {

	PyObject *model;

	if (!PyArg_ParseTuple(args,"O",&model)) return NULL;

	PyObject *name = NULL, *number = NULL;

	name = PyObject_GetAttrString(model, "filename");
	const char *filename = PyString_AsString(name);
	Py_DECREF(name);

	number = PyObject_GetAttrString(model, "xmax");
	float xmax =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(model, "ymax");
	float ymax =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(model, "zmax");
	float zmax =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(model, "xoffset");
	float xoffset =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(model, "yoffset");
	float yoffset =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(model, "zoffset");
	float zoffset =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	// The image reference

	PyObject *imgref = PyObject_GetAttrString(model, "imgref");

	number = PyObject_GetAttrString(imgref, "xscale_old");
	float xscale_old =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "yscale_old");
	float yscale_old =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "zscale_old");
	float zscale_old =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "xscale_new");
	float xscale_new =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "yscale_new");
	float yscale_new =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "zscale_new");
	float zscale_new =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "xtranslation_old");
	float xtranslation_old =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "ytranslation_old");
	float ytranslation_old =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "ztranslation_old");
	float ztranslation_old =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "xtranslation_new");
	float xtranslation_new =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "ytranslation_new");
	float ytranslation_new =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "ztranslation_new");
	float ztranslation_new =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "xrotation_old");
	float xrotation_old =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "yrotation_old");
	float yrotation_old =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "zrotation_old");
	float zrotation_old =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "xrotation_new");
	float xrotation_new =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "yrotation_new");
	float yrotation_new =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	number = PyObject_GetAttrString(imgref, "zrotation_new");
	float zrotation_new =  (float)PyFloat_AsDouble(number);
	Py_DECREF(number);

	Py_DECREF(imgref);

//	printf("Save model in file %s\n",filename);

	PyObject *objects = NULL, *contours = NULL;
	PyArrayObject *points;
	PyObject *object = NULL, *contour = NULL;

	Imod *imod_model = imodNew();

	imod_model->xoffset = xoffset;
	imod_model->yoffset = yoffset;
	imod_model->zoffset = zoffset;

	imod_model->xmax = xmax;
	imod_model->ymax = ymax;
	imod_model->zmax = zmax;

	IrefImage *iref = (IrefImage *)malloc(sizeof(IrefImage));
	if (!iref) return NULL;

	imod_model->refImage = iref;

	imod_model->refImage->oscale.x = xscale_old;
	imod_model->refImage->oscale.y = yscale_old;
	imod_model->refImage->oscale.z = zscale_old;

	imod_model->refImage->cscale.x = xscale_new;
	imod_model->refImage->cscale.y = yscale_new;
	imod_model->refImage->cscale.z = zscale_new;

	imod_model->refImage->otrans.x = xtranslation_old;
	imod_model->refImage->otrans.y = ytranslation_old;
	imod_model->refImage->otrans.z = ztranslation_old;

	imod_model->refImage->ctrans.x = xtranslation_new;
	imod_model->refImage->ctrans.y = ytranslation_new;
	imod_model->refImage->ctrans.z = ztranslation_new;

	imod_model->refImage->orot.x = xrotation_old;
	imod_model->refImage->orot.y = yrotation_old;
	imod_model->refImage->orot.z = zrotation_old;

	imod_model->refImage->crot.x = xrotation_new;
	imod_model->refImage->crot.y = yrotation_new;
	imod_model->refImage->crot.z = zrotation_new;

	imod_model->pixsize = 1;
	imod_model->units = IMOD_UNIT_PIXEL;

	Py_ssize_t number_of_objects,number_of_contours,number_of_points;
	Py_ssize_t iobj,icont,k;

	objects = PyObject_GetAttrString(model, "objects");

	if (!PyList_Check(objects)) {
		PyErr_SetString(PyExc_StandardError, "Attribute objects from the model should be a list!");
		return NULL;
	}

	number_of_objects = PyList_Size(objects);

//	printf("Number of objects: %i\n",(int)number_of_objects);

	Iobj *imod_object;

	npy_intp *dimensions;
	char *data;

	for (iobj=0; iobj<number_of_objects; iobj++) {

		object = PyList_GetItem(objects,iobj);

		contours = PyObject_GetAttrString(object, "contours");

		if (!PyList_Check(contours)) {
			PyErr_SetString(PyExc_StandardError, "Attribute contours from an object should be a list!");
			return NULL;
		}

		number_of_contours = PyList_Size(contours);

		imodNewObject(imod_model);

		imod_object = imodObjectGet(imod_model);

//		To fix
//		name = PyObject_GetAttrString(object, "name");
//		char *nn = PyString_AsString(name);
//		imodObjectSetName(imod_object, nn);
//		Py_DECREF(name);

		imod_object->symbol = IOBJ_SYM_CIRCLE;	// to see the symbols
		imod_object->symsize = 5;	// size of the markers (for the volume)
                imod_object->pdrawsize = 5;	// sphere radius for point (for the model)
		imod_object->flags |= IMOD_OBJFLAG_SCAT;  // open objects
		imod_object->flags |= IMOD_OBJFLAG_PNT_ON_SEC;  // show on points of sections

		for (icont=0; icont<number_of_contours; icont++) {

			contour = PyList_GetItem(contours,icont);

			points = (PyArrayObject *)PyObject_GetAttrString(contour, "points");

			if (points==NULL) {
				PyErr_SetString(PyExc_StandardError, "Attribute points from an contour should be a numeric array!");
				return NULL;
			}

			dimensions = points->dimensions;
			number_of_points = dimensions[0];

			data = (char *)PyArray_DATA(points);

			imodNewContour(imod_model);

			for (k=0; k<number_of_points; k++) {
				
				Ipoint point;

				point.x = (float)*(double *)(data + k*points->strides[0] + 0*points->strides[1]);
				point.y = (float)*(double *)(data + k*points->strides[0] + 1*points->strides[1]);
				point.z = (float)*(double *)(data + k*points->strides[0] + 2*points->strides[1]);

				imodNewPoint(imod_model, &point);

			}

			Py_DECREF(points);

		}

		Py_DECREF(contours);

	}
	
	Py_DECREF(objects);

	// Save the model in a file

	FILE* fout;

	if ((fout = fopen(filename,"w+"))==NULL) {
		printf("Cannot open file %s",filename);
		return NULL;
	}
	
	imodWrite(imod_model, fout);

	if (fclose(fout)) {
		return NULL;
	}

	imodDelete(imod_model);

    return Py_BuildValue( "i", 1 );

}


static PyMethodDef ModelMethods[] = {
	{ "XYZMax", getXYZMax, METH_VARARGS, "Get the Maximum Size in X Y and Z." },
	{ "getDimensions", getDimensions, METH_VARARGS, "Get the dimensions in X Y and Z." },
	{ "getScales", getScales, METH_VARARGS, "Get the scales in X Y and Z." },
	{ "getImgRefParameters", getImgRefParameters, METH_VARARGS, "Get the image reference parameters for an IMOD model." },
	{ "numberOfObjects", numberOfObjects, METH_VARARGS, "Get the number of Objects for a given model." },
	{ "numberOfContours", numberOfContours, METH_VARARGS, "Get the number of Contours for a given Object." },
	{ "numberOfPoints", numberOfPoints, METH_VARARGS, "Get the number of Points for a given Contour." },
	{ "loadPoints", loadPoints, METH_VARARGS, "Load the Points for a given Contour." },
	{ "loadModel", loadModel, METH_VARARGS, "Load a Model." },
	{ "saveModel", saveModel, METH_VARARGS, "Save a Model." },
	{NULL, NULL, 0, NULL} /* Sentinel */
};


PyMODINIT_FUNC initmodfile(void) {

	(void)Py_InitModule("modfile", ModelMethods);

	import_array();

}

int main(int argc, char *argv[]) {

	Py_SetProgramName(argv[0]);	/* Pass argv[0] to the Python interpreter */

	Py_Initialize(); /* Initialize the Python interpreter.	*/

	initmodfile();	/* Add a static module */

	return 1;

}
