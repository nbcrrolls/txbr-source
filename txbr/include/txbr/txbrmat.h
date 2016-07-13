#ifndef TXBRMAT_H
#define TXBRMAT_H 

#include "mat.h"
#include "matrix.h"

void print_mxArray(mxArray* mx, char* array_name);
double* getDataPtr(MATFile* mat, char** path, int* index, int n, int length);

#endif
