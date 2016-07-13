#include <string.h>
#include "txbrutil.h"

#define CR 13            /* Decimal code of Carriage Return char */
#define LF 10            /* Decimal code of Line Feed char */
#define EOF_MARKER 26    /* Decimal code of DOS end-of-file marker */
#define MAX_REC_LEN 2046 /* Maximum size of input buffer */

char basename_[FILENAME_LEN];
char directory_[FILENAME_LEN];
int order_, number_of_tilts_;
int skipView[MAX_TILT];
float lambda_[MAX_TILT][4];
float coefficients_1_[MAX_TILT][MAX_COEFFICIENTS_NUMBER];
float coefficients_2_[MAX_TILT][MAX_COEFFICIENTS_NUMBER];
float coefficients_3_[MAX_TILT][MAX_COEFFICIENTS_NUMBER];

int read_trajectory_file(char* configuration_file_name) {

	FILE *input = fopen(configuration_file_name,"rt");

	long file_length, index;
	char *file_buffer, *pointer;
	char line[MAX_REC_LEN];

	fseek(input, 0L, SEEK_END);
	file_length = ftell(input);
	rewind(input);

	file_buffer = calloc(file_length + 1, sizeof(char));

	if(!file_buffer) {

		txbr_error(stderr, "ERROR: read_trajectory_file - getting memory.\n");

		return MEMORY_ERROR;
	}

	fread(file_buffer, file_length, 1, input); /* Read the entire file into file_buffer */

	pointer = file_buffer;

	while (*pointer) {

		index = 0;

		while (*pointer) {

			if (*pointer == CR || *pointer == LF) {

				*pointer++;
				break;

			}

			line[index++] = *pointer++;

		}

		line[index] = '\0';

		if (index!=0) read_parameters(&line[0]);

  }

  free(file_buffer);

  fclose(input);

  return NO_ERROR;

}

int read_parameters(char* line) {

	int index_of_tilt, index;
	float t;

	char format0[MAX_REC_LEN] = "%*s";
	char format[MAX_REC_LEN] = "%*s %%f";

	/* sscanf works with float	*/

	sscanf(line,"basename: %s",&basename_[0]);
	sscanf(line,"directory: %s",&directory_[0]);

	sscanf(line,"trajectory approximation order: %i",&order_);
	sscanf(line,"number of tilts: %i",&number_of_tilts_);

	if (sscanf(line,"lambda-%i:",&index_of_tilt)>=1) {
		index = 0;
		sprintf(&format[0],"%s %%f",&format0[0]);
		while (sscanf(line,&format[0],&t)>=1) {
			sprintf(&format0[0],"%s %%*f",&format0[0]);
			sprintf(&format[0],"%s %%f",&format0[0]);
			lambda_[index_of_tilt][index] = t;
			index++;
		}
	}

	if (sscanf(line,"x-%i:",&index_of_tilt)>=1) {
		index = 0;
		sprintf(&format[0],"%s %%f",&format0[0]);
		while (sscanf(line,&format[0],&t)>=1) {
			sprintf(&format0[0],"%s %%*f",&format0[0]);
			sprintf(&format[0],"%s %%f",&format0[0]);
			coefficients_1_[index_of_tilt][index] = t;
			index++;
		}
	}

	if (sscanf(line,"y-%i:",&index_of_tilt)>=1) {
		index = 0;
		sprintf(&format[0],"%s %%f",&format0[0]);
		while (sscanf(line,&format[0],&t)>=1) {
			sprintf(&format0[0],"%s %%*f",&format0[0]);
			sprintf(&format[0],"%s %%f",&format0[0]);
			coefficients_2_[index_of_tilt][index] = t;
			index++;
		}
	}

}

/*
 * Load TxBRsetup and PathProjection structures from the configuration file
 */
void load_trajectory(char* directory, char* basename, TrajectoryMap* trajectory) {

	char configuration_file_name[FILENAME_LEN];
	sprintf(&configuration_file_name[0],"%s/%s.traj",directory,basename);

	printf("Load trajectory from %s\n", &configuration_file_name[0]);

	int i,j;

	init_int_array(&skipView[0],MAX_TILT);
	init_float_array(&lambda_[0][0],4*MAX_TILT);
	init_float_array(&coefficients_1_[0][0],MAX_COEFFICIENTS_NUMBER*MAX_TILT);
	init_float_array(&coefficients_2_[0][0],MAX_COEFFICIENTS_NUMBER*MAX_TILT);
	init_float_array(&coefficients_3_[0][0],MAX_COEFFICIENTS_NUMBER*MAX_TILT);

	read_trajectory_file(configuration_file_name);

	trajectory->order = order_;
	trajectory->numberOfTilts = number_of_tilts_;
	trajectory->number_of_terms = number_of_polynomial_terms(order_);

	int n = number_of_polynomial_terms(order_);

	for (i=0; i<number_of_tilts_; i++) {
		trajectory->skipView[i] = skipView[i];
		for (j=0; j<4; j++) {
			trajectory->lambda[i][j] = lambda_[i][j];
		}
		for (j=0; j<n; j++) {
			trajectory->coefficients_1[i][j] = coefficients_1_[i][j];
			trajectory->coefficients_2[i][j] = coefficients_2_[i][j];
		}
	}

}

/*
 * Print informations about the electron trajectory coefficients.
 */
void print_axis_path_projection_for_tilt( TrajectoryMap *trajectory, int indexOfTilt, int axis ) {

	int i1,i2,index0,n;
	double coeffs[MAX_COEFFICIENTS_NUMBER];
	char* title = "";

	switch (axis) {
		case 0 : memcpy(coeffs, &trajectory->lambda[indexOfTilt][0], 4*sizeof(double)); title="lambda"; break;
		case 1 : memcpy(coeffs, &trajectory->coefficients_1[indexOfTilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double)); title="     x"; break;
		case 2 : memcpy(coeffs, &trajectory->coefficients_2[indexOfTilt][0], MAX_COEFFICIENTS_NUMBER*sizeof(double)); title="     y"; break;
	}

	printf("%s:\t",title);

	int effective_order = axis==0 ? 1 : trajectory->order;

	for (i1=0; i1<effective_order+1; i1++) {
		n = i1==0 ? 1 : (i1+1)*(i1+2)/2;
		index0 = number_of_polynomial_terms(i1-1);
		for (i2=0; i2<n; i2++) {
			printf("%.1E ",coeffs[index0 + i2]);
		}
		printf(" ");
	}

	printf("\n");

}

/*
 * Print informations about the electron trajectories.
 */
void print_trajectory( TrajectoryMap *trajectory ) {

	int numberOfTilts = trajectory->numberOfTilts;
	int i;

	printf("-------------------------- Projection Coefficients --------------------------------------------\n");
	printf("Order           : %i\n",trajectory->order);
	printf("Number of Tilts : %i\n",trajectory->numberOfTilts);

	for (i=0; i<numberOfTilts; i++) {
		if (i%20!=0) continue; /* We only print the info every 20 tilts	*/
		printf("%i/%i\n",i,numberOfTilts);
		print_axis_path_projection_for_tilt(trajectory,i,0);
		print_axis_path_projection_for_tilt(trajectory,i,1);
		print_axis_path_projection_for_tilt(trajectory,i,2);
	}
	printf("-----------------------------------------------------------------------------------------------\n");

}

int write_trajectory_file(TxBRsetup* setup, TrajectoryMap* trajectory, char* configuration_file_name) {

	FILE *output = fopen(configuration_file_name,"wt");

	fprintf(output,"basename: %s\n",setup->basename);
	fprintf(output,"directory: %s\n",setup->directory);

	fprintf(output,"\n");

	fprintf(output,"trajectory approximation order: %i\n",trajectory->order);

	fprintf(output,"\n");

	fprintf(output,"number of tilts: %i\n",trajectory->numberOfTilts);

	fprintf(output,"\n");

	int i,j;

	int n = number_of_polynomial_terms(trajectory->order);

	for (i=0; i<trajectory->numberOfTilts; i++) {

		fprintf(output,"lambda-%i: ",i);
		for (j=0; j<4; j++) fprintf(output,"%e ",trajectory->lambda[i][j]);
		fprintf(output,"\n");

		fprintf(output,"x-%i: ",i);
		for (j=0; j<n; j++) fprintf(output,"%e ",trajectory->coefficients_1[i][j]);
		fprintf(output,"\n");

		fprintf(output,"y-%i: ",i);
		for (j=0; j<n; j++) fprintf(output,"%e ",trajectory->coefficients_2[i][j]);
		fprintf(output,"\n");

	}

	fclose(output);

	return NO_ERROR;

}

