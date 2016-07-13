#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "txbrutil.h"
#include "mrcfiles.h"
#include "mrcslice.h"

#define TILT_FLAG 1

int mrc_data_new(FILE *fout, MrcHeader *hdata)
{
  int dsize, i;
  unsigned char cdata = 0;
  short sdata = 0;
  float fdata = 0;

  /* 6/26/01: change from error message to swapping fdata */
  if (hdata->swapped) {
    mrc_swap_floats(&fdata, 1);
  }

  dsize = hdata->nx * hdata->ny * hdata->nz;
  rewind(fout);

  fseek(fout, hdata->headerSize, SEEK_CUR);

  switch(hdata->mode)
    {
    case MRC_MODE_BYTE:
      for (i = 0; i < dsize; i++)
        if (!fwrite(&cdata, 1, 1, fout))
          return(-1);
      break;

    case MRC_MODE_SHORT:
      for (i = 0; i < dsize; i++)
        if (!fwrite(&sdata, 2, 1, fout))
          return(-1);
      break;

    case MRC_MODE_FLOAT:
      for (i = 0; i < dsize; i++)
        if (!fwrite(&fdata, 4, 1, fout))
          return(-1);
      break;

    default:
      return(-1);

    }
  return(0);
}

/*
 * Display information contained in an MRC file
 */
int header_info(char* mrcFileName) {
	
	FILE* fp = fopen( mrcFileName, "r");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: header_info - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: header_info - getting memory.\n");
		return MEMORY_ERROR;
	}

	float xs,ys,zs;

	mrc_head_read(fp, header);
	mrc_get_scale(header, &xs, &ys, &zs);
	
	printf("\n");
	
	printf("Size   : nx,ny,nz = %i,%i,%i\n", header->nx, header->ny, header->nz);
	printf("Grid   : mx,my,mz = %i,%i,%i\n", header->mx, header->my, header->mz);
	printf("Length : xlen,ylen,zlen = %.1f,%.1f,%.1f\n", header->xlen, header->ylen, header->zlen);
	printf("Scale  : xs,ys,zs = %.2f,%.2f,%.2f\n",xs,ys,zs);

	int i=0;

        if (header->nreal & TILT_FLAG)
        {
            short angle = 0;
            fseek(fp, 1024, SEEK_SET);
            printf("Tilt angles: ");
            for (i = 0; i < header->nz; i++)
            {
                fread(&angle, 2, 1, fp);
                printf("%d->%.1f", i+1, angle/100.0);
                if (i==header->nz-1) printf("\n");
                else printf(", ");
                fseek(fp, header->nint - 2, SEEK_CUR);

            }
        }

	printf("\n");

        printf("Description:\n");

	for (i=0;i<header->nlabl;i++)
	{
		printf("%s\n",header->labels[i]);
	}

	if (fclose(fp))
	{
		txbr_error(stderr, "ERROR: header_info - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


int get_tilt_angles(char* mrcFileName, int nangles, float *angles) {

    FILE* fp = fopen( mrcFileName, "r");

    if (fp == NULL) {
        txbr_error(stderr, "ERROR: header_info - cannot open file %s.\n", &mrcFileName[0]);
        return FILE_ERROR;
    }

    MrcHeader* header = (MrcHeader *) malloc(sizeof (MrcHeader));

    if (!header) {
        txbr_error(stderr, "ERROR: header_info - getting memory.\n");
        return MEMORY_ERROR;
    }

    mrc_head_read(fp, header);

    if (header->nreal & TILT_FLAG) {

        short angle = 0;
        fseek(fp, 1024, SEEK_SET);

        int i;

        if (nangles!=(header->nz)) return PARAMETER_ERROR;
        
        for (i = 0; i < header->nz; i++) {
            fread(&angle, 2, 1, fp);
            angles[i] = (float)angle/100.0;
            fseek(fp, header->nint - 2, SEEK_CUR);
        }

    }

    if (fclose(fp)) {
        txbr_error(stderr, "ERROR: header_info - closing file %s.\n", &mrcFileName[0]);
        return FILE_ERROR;
    }

    free(header);

    return NO_ERROR;

}


/*
 * Get the spatial dimensions (columns,rows,sections) of an MRC image stack stored in file 
 * called mrcFileName. Dimensions are stored in 3 integers whom pointers (nx,ny,nz) are passed 
 * as variables.
 */
int get_size_of(const char* mrcFileName, int *nx, int *ny, int *nz) {

	FILE* fp = fopen( mrcFileName, "r");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_size_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_size_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read(fp, header);

	*nx = header->nx;
	*ny = header->ny;
	*nz = header->nz;

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_size_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


/*
 * Get the scale of an MRC image stack stored in file called mrcFileName. 
 */
int get_scale_of(const char* mrcFileName, float *sx, float *sy, float *sz) {

	FILE* fp = fopen( mrcFileName, "r");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_scale_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_scale_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read(fp, header);
	mrc_get_scale( header, sx, sy, sz);

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_scale_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


int update_scale_of( char* mrcFileName, float sx, float sy, float sz )
{

	FILE* fp = fopen( mrcFileName, "r+");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: update_scale_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: update_scale_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read( fp, header );
	
	mrc_set_scale( header, sx, sy, sz );
	
	mrc_head_write( fp, header );

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: update_scale_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


/*
 * Get the type mode of an MRC image stack stored in file called mrcFileName.
 */
int get_mode_of(const char* mrcFileName, int *mode) {

	FILE* fp = fopen( mrcFileName, "r");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_size_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_size_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read( fp, header );

	*mode = header->mode;

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_size_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}

/*
 * Get the description of an MRC image stack stored in file called mrcFileName.
 */
int get_description_of(const char* mrcFileName, int* nlabl, char** labels) {

	FILE* fp = fopen( mrcFileName, "r");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_size_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_size_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read(fp, header);

	*nlabl = header->nlabl;

	int i=0;

	for ( i=0; i<header->nlabl; i++ ) 
	{
		memcpy(labels[i], header->labels[i], (MRC_LABEL_SIZE+1)*sizeof(char));
	}

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_size_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


/*
 * Update the description of an MRC image stack stored in file called mrcFileName.
 */
int update_description_of(const char* mrcFileName, const int nlabl, const char** labels) {

	FILE* fp = fopen( mrcFileName, "r+");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_size_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_size_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read( fp, header );
	
	header->nlabl = TXBR_MIN(nlabl,MRC_NLABELS);	
	
	int i;
		
	for (i=0; i<header->nlabl; i++) 
	{
		strncpy(header->labels[i], labels[i], MRC_LABEL_SIZE);
	}
	
	mrc_head_write( fp, header );
	
	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_size_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}
	
	free(header);

	return NO_ERROR;

}

/*
 * Add label to the description
 */
int add_label(const char* mrcFileName, const char* label) {

	FILE* fp = fopen( mrcFileName, "r+");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_size_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_size_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read( fp, header );
	
	mrc_head_label(header,label);
	
	mrc_head_write( fp, header );
	
	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_size_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}
	
	free(header);

	return NO_ERROR;

}


/*
 * Get the length of a MRC volume in the three directions.
 */
int get_length_of(const char* mrcFileName, float* xlen, float* ylen, float* zlen) {

	FILE* fp = fopen( mrcFileName, "r");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_length_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_length_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read(fp, header);
	
	*xlen = header->xlen;
	*ylen = header->ylen;
	*zlen = header->zlen;

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_length_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


/*
 * Update the length dimensions of a MRC volume in the three directions. If a value in one
 * direction is negative, nothing is done for this axis.
 */
int update_length_of(const char* mrcFileName, const float xlen, const float ylen, const float zlen) {

	FILE* fp = fopen( mrcFileName, "r+");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_size_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_size_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read(fp, header);

	header->xlen = xlen;
	header->ylen = ylen;
	header->zlen = zlen;
	
	mrc_head_write( fp, header );

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_size_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


/*
 * Get the point origin of an MRC image stack stored in file called mrcFileName. 
 * Coordinates are stored in 3 float whom pointers (xorg,yorg,zorg) are passed 
 * as variables.
 */
int get_origin_of(const char* mrcFileName, float *xorg, float *yorg, float *zorg) {

	FILE* fp = fopen( mrcFileName, "r");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_origin_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_origin_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read(fp, header);

	*xorg = header->xorg;
	*yorg = header->yorg;
	*zorg = header->zorg;

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_origin_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


/*
 * Update the origin point of an MRC volume in the three directions.
 */
int update_origin_of(const char* mrcFileName, const float xorg, const float yorg, const float zorg) {

	FILE* fp = fopen( mrcFileName, "r+");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: update_origin_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: update_origin_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read(fp, header);

	header->xorg = xorg;
	header->yorg = yorg;
	header->zorg = zorg;
	
	mrc_head_write( fp, header );

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: update_origin_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


/*
 * Get the grid origin of a MRC volume in the three directions.
 */
int get_grid_origin_of(const char* mrcFileName, int* nxstart, int* nystart, int* nzstart) {

	FILE* fp = fopen( mrcFileName, "r");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_grid_origin_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_grid_origin_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read(fp, header);

	*nxstart = header->nxstart;
	*nystart = header->nystart;
	*nystart = header->nystart;

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_grid_origin_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


/*
 * Update the grid size of a MRC volume in the three directions. If a value in one
 * direction is negative, nothing is done for this axis.
 */
int update_grid_origin_of(const char* mrcFileName, const int nxstart, const int nystart, const int nzstart) {

	FILE* fp = fopen( mrcFileName, "r+");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: update_grid_origin_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: update_grid_origin_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read(fp, header);

	header->nxstart = nxstart;
	header->nystart = nystart;
	header->nzstart = nzstart;
	
	mrc_head_write( fp, header );

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: update_grid_origin_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


/*
 * Get the grid size of a MRC volume in the three directions.
 */
int get_grid_size_of(const char* mrcFileName, int* mx, int* my, int* mz) {

	FILE* fp = fopen( mrcFileName, "r");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_size_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_size_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read(fp, header);

	*mx = header->mx;
	*my = header->my;
	*mz = header->mz;

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_size_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


/*
 * Update the grid size of a MRC volume in the three directions. If a value in one
 * direction is negative, nothing is done for this axis.
 */
int update_grid_size_of(const char* mrcFileName, const int mx, const int my, const int mz) {

	FILE* fp = fopen( mrcFileName, "r+");

	if (fp==NULL) 
	{
		txbr_error(stderr, "ERROR: get_size_of - cannot open file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) 
	{
		txbr_error(stderr, "ERROR: get_size_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	mrc_head_read(fp, header);

	header->mx = mx;
	header->my = my;
	header->mz = mz;
	
	mrc_head_write( fp, header );

	if (fclose(fp)) 
	{
		txbr_error(stderr, "ERROR: get_size_of - closing file %s.\n", &mrcFileName[0]);
		return FILE_ERROR;
	}

	free(header);

	return NO_ERROR;

}


/*
 * Extract data (of ndata items) from a MRC slice of a given mode, and copy
 * it into an array of float (whose size is supposed to be already allocated)
 * The routine also provides some statistics about the data.
 */
void extract_float_array( Islice* slice, int ndata, int mode, float* data) {

	if (mode==SLICE_MODE_FLOAT) {
		memcpy(&data[0], &slice->data.f[0], ndata*sizeof(float));
	}

	if (mode==SLICE_MODE_BYTE) {
		int i;
		for (i=0; i<ndata; i++) {
			data[i] = (float)(slice->data.b[i]);
		}
	}

	if (mode==SLICE_MODE_SHORT) {
		int i;
		for (i=0; i<ndata; i++) {
			data[i] = (float)(slice->data.s[i]);
		}
	}

	if (mode==SLICE_MODE_USHORT) {
		int i;
		for (i=0; i<ndata; i++) {
			data[i] = (float)(slice->data.us[i]);
		}
	}

}


/*
 * Get some statistics on an MRC file mrcFileName.This routine calculates:
 * m0 is the number of elements.
 * m1 is the mean value.
 * m2 is the standard deviation
 * min is the minimal value
 * max is the maximal value
 */
int get_statistics_of(const char* mrcFileName, float *m0, float *m1, float *m2, float *min, float *max) {

	/* First get the mode of the MRC stack file - Everything will be converted to float	*/

	int mode = SLICE_MODE_FLOAT;
	int error = get_mode_of(mrcFileName, &mode);

	if (error!=NO_ERROR) return error;

	/* Now calculate the statistics	*/

	FILE* fp = fopen( mrcFileName, "r");

	if (fp==NULL) {

		txbr_error(stderr, "ERROR: get_statistics_of - cannot open file %s.\n", &mrcFileName[0]);

		return FILE_ERROR;

	}

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) {

		txbr_error(stderr, "ERROR: get_statistics_of - getting memory.\n");

		return MEMORY_ERROR;

	}

	mrc_head_read( fp, header );

	int islice=0, i=0;
	int numberOfSlices = header->nz;
        numberOfSlices = 1;
	int ndata = (header->nx)*(header->ny);
	Islice* slice = NULL;

	float *data = (float *)malloc(ndata*sizeof(float));

	if (!data) {
		txbr_error(stderr, "ERROR: get_statistics_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	*m0 = ndata*numberOfSlices;
	*m1 = 0.0;
	*m2 = 0.0;
	*min = FLT_MAX;
	*max = -FLT_MAX;

	for (islice=0; islice<numberOfSlices; islice++) {

		slice = sliceReadMRC(header,islice,'z');

		extract_float_array( slice, ndata, mode, data );

		for (i=0; i<ndata; i++) {

			*m1 = *m1 + data[i]/ndata/numberOfSlices;
			*m2 = *m2 + pow(data[i],2)/ndata/numberOfSlices;
			*min = TXBR_MIN(*min,data[i]);
			*max = TXBR_MAX(*max,data[i]);

		}

		if (slice) sliceFree(slice);

	}

	if (fclose(fp)) {

		txbr_error(stderr, "ERROR: get_statistics_of - closing file %s.\n", &mrcFileName[0]);

		return FILE_ERROR;

	}

	/*	free stuff	*/

	free(data);
	free(header);

	return NO_ERROR;

}

/*
 * Calculate the total number of slices for a serie of MRC Files! All MRC files
 * should have the same xy dimensions - if not returns a PARAMETER_ERROR. If a
 * MRC file cannot be opened , the routine returns a The pointer
 * number_of_slices points to an integer showing the total number of slices.
 * Incidentally the function also returns the size of the stack in x and y.
 */
int get_dimensions(char** mrcFileNames, int number_of_files, int *nx, int *ny, int *number_of_slices ) {

	int error=NO_ERROR;

	*nx=-1;
	*ny=-1;

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) {
		txbr_error(stderr, "ERROR: get_total_number_of_slices - getting memory.\n");
		return MEMORY_ERROR;
	}

	int index;
	FILE* fp;

	*number_of_slices = 0;

	for (index=0; index<number_of_files; index++) {

		fp = fopen(&mrcFileNames[index][0],"r");

		if (fp==NULL) {	/* We skip the file	*/
			txbr_error(stderr, "ERROR: get_total_number_of_slices - cannot open file %s.\n",
				&mrcFileNames[index][0]);
			error = FILE_ERROR;
			continue;
		}

		mrc_head_read( fp, header );

		if (*nx==-1 && *ny==-1) {
			*nx = header->nx;
			*ny = header->ny;
		}

		if (*nx!=header->nx && *ny!=header->ny) error = PARAMETER_ERROR;

		*number_of_slices += *number_of_slices + 1;

		if( fclose(fp) ) {
			txbr_error(stderr, "ERROR: get_total_number_of_slices - closing file %s.\n", &mrcFileNames[index][0]);
			error = FILE_ERROR;
		}

	}

	free(header);

	return error;

}

/*
 * Stack together a series of MRC Files! xy Slices should be of the same
 * size! Final stack is of float type and stored in the file outputFileName.
 */
int stack_together(char** mrcFileNames, int number_of_files, char* outputFileName) {

	/*	First evaluate the dimensions of the new stack	*/

	int nx, ny, number_of_slices;

	int error = get_dimensions( mrcFileNames, number_of_files, &nx, &ny, &number_of_slices );

	if (error!=NO_ERROR) return error;

	FILE* output = fopen(&outputFileName[0],"w");

	if(output==NULL) {
		txbr_error(stderr, "ERROR: stack_together - cannot open file %s.\n", &outputFileName[0]);
		return FILE_ERROR;
	}

	/* Creates the MRC header for the main output file	*/

	struct MRCheader hdata;

	mrc_head_new(&hdata, nx, ny, number_of_slices, SLICE_MODE_FLOAT);

	mrc_head_write(output, &hdata);

	/* Allocate a MRC header for current file	*/

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) {
		txbr_error(stderr, "ERROR: get_total_number_of_slices - getting memory.\n");
		return MEMORY_ERROR;
	}

	printf("Merging #%i blocks together\n",number_of_files);

	FILE* fp;
	Islice* slice = NULL;

	int iblock, islice, n, ndata, mode, i, index=0;
	float mean=0, dev=0, min=FLT_MAX, max=-FLT_MAX;

	float *data = (float *)malloc(ndata*sizeof(float));

	if (!data) {
		txbr_error(stderr, "ERROR: get_statistics_of - getting memory.\n");
		return MEMORY_ERROR;
	}

	for (iblock=0; iblock<number_of_files; iblock++) {

		printf( "Processing block #%i... File %s.\n", iblock, mrcFileNames[iblock]);

		fp = fopen(mrcFileNames[iblock],"r");

		if (fp==NULL) {
			txbr_error(stderr, "ERROR: merge_blocks - cannot open file %s.\n", &mrcFileNames[iblock][0]);
			return FILE_ERROR;
		}

		mrc_head_read( fp, header );

		n = header->nz;
		ndata = (header->nx)*(header->ny);
		mode = header->mode;

		for (islice=0; islice<n; islice++) {

			printf( "Reading Slice #%i of block %i\n", islice, iblock);

			slice = sliceReadMRC(header,islice,'z');

			extract_float_array( slice, ndata, mode, data);

			for (i=0; i<ndata; i++) {
				mean = mean + data[i]/ndata/number_of_slices;
				dev = dev + pow(data[i],2)/ndata/number_of_files;
				min = TXBR_MIN(min,data[i]);
				max = TXBR_MAX(max,data[i]);
			}

			mrc_write_slice(slice->data.f, output, &hdata, index++, 'z');

			if (slice) sliceFree(slice);

		}

		if (fclose(fp)) {
			txbr_error(stderr, "ERROR: merge_blocks - closing file %s.\n", &mrcFileNames[iblock][0]);
			return FILE_ERROR;
		}

	}

	dev = dev - pow(mean,2);

	hdata.amin = min;
	hdata.amax = max;
	hdata.amean = mean;

	mrc_head_write(output, &hdata);

	if(fclose(output)) {
		txbr_error(stderr, "ERROR: merge_blocks - closing file %s.\n", &outputFileName[0]);
		return FILE_ERROR;
	}

	/* Free stuff	*/

	free(data);
	free(header);

	return NO_ERROR;

}









/*
 * Merge blocks together!
 */
int merge_blocks(char* directory, char* basename, double z0, int numberOfBlocks) {

	char imodfile_out[FILENAME_LEN];

	int nx = 0;
	int ny = 0;

	MrcHeader* header  = (MrcHeader *)malloc(sizeof(MrcHeader));

	if (!header) {

		txbr_error(stderr, "ERROR: merge_blocks - getting memory.\n");

		return MEMORY_ERROR;

	}

	FILE* fp;
	int iblock;

	double min = 0.0;
	double max = 0.0;
	double mean = 0.0;

	char files[MAX_NUMBER_OF_FILES_TO_MERGE][FILENAME_LEN];

	int number_of_files_to_merge = 0;

	int totalNumberOfSlices = 0;

	/*	Determine first the right parameters (nx,ny,numberOfSlices)	by reading the MRC header.*/

	for (iblock=0; iblock<numberOfBlocks; iblock++) {

		sprintf(&files[number_of_files_to_merge][0], "%s/%s_z_%.2f-%i.out", directory, basename, z0, iblock);

		if ((fp=fopen(files[number_of_files_to_merge],"r"))==NULL) {

			txbr_error(stderr, "ERROR: merge_blocks - cannot open file %s.\n", &files[number_of_files_to_merge][0]);

			return FILE_ERROR;

		}

		printf("Cheak header: %s\n",files[number_of_files_to_merge]);

		mrc_head_read(fp,header);

		if (nx==0 && ny==0) {
			nx = header->nx;
			ny = header->ny;
		}

		if (header->nx==nx && header->ny==ny) {

			number_of_files_to_merge++;

			min = TXBR_MIN(min,header->amin);
			max = TXBR_MAX(max,header->amax);
			mean += header->amean*header->nz;

			totalNumberOfSlices += header->nz;

		}

		if(fclose(fp)) {

			txbr_error(stderr, "ERROR: merge_blocks - closing file %s.\n", &files[number_of_files_to_merge][0]);

			return FILE_ERROR;

		}

	}

	/*	Right parameters are now known	*/

	sprintf(imodfile_out, "%s/%s_z_%.2f.out", directory, basename, z0);

	FILE* output;

	if((output=fopen(&imodfile_out[0],"w"))==NULL) {

		txbr_error(stderr, "ERROR: merge_blocks - cannot open file %s.\n", &imodfile_out[0]);

		return FILE_ERROR;

	}

	/* MRC header for the main output file	*/

	struct MRCheader hdata;

	mrc_head_new(&hdata, nx, ny, totalNumberOfSlices, SLICE_MODE_FLOAT);

	hdata.amin = min;
	hdata.amax = max;
	hdata.amean = mean/number_of_files_to_merge/totalNumberOfSlices;

	mrc_head_write(output, &hdata);

	/* MRC header for current file	*/

	printf("Merging #%i blocks together\n",number_of_files_to_merge);

	Islice* slice = NULL;
	int numberOfSlices = 0;	/*	normally should be blocksize	*/

	int islice=0, index=0;

	for (iblock=0; iblock<number_of_files_to_merge; iblock++) {

		printf("Processing block #%i... File %s.\n",iblock,files[iblock]);

		if((fp=fopen(files[iblock],"r"))==NULL) {

			txbr_error(stderr, "ERROR: merge_blocks - cannot open file %s.\n", &files[iblock][0]);

			return FILE_ERROR;

		}

		mrc_head_read(fp,header);

		numberOfSlices = header->nz;

		for (islice=0; islice<numberOfSlices; islice++) {

			printf("Reading Slice #%i of block %i\n",islice,iblock);

			slice = sliceReadMRC(header,islice,'z');

			mrc_write_slice(slice->data.f, output, &hdata, index++, 'z');

			if (slice) sliceFree(slice);

		}

		if(fclose(fp)) {

			txbr_error(stderr, "ERROR: merge_blocks - closing file %s.\n", &files[iblock][0]);

			return FILE_ERROR;

		}

	}

	if(fclose(output)) {

		txbr_error(stderr, "ERROR: merge_blocks - closing file %s.\n", &imodfile_out[0]);

		return FILE_ERROR;

	}

	/* Free stuff	*/

	free(header);

	return NO_ERROR;

}


