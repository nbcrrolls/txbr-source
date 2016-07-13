#include <string.h>
#include "txbrutil.h"

#define CR 13            /* Decimal code of Carriage Return char */
#define LF 10            /* Decimal code of Line Feed char */
#define EOF_MARKER 26    /* Decimal code of DOS end-of-file marker */
#define MAX_REC_LEN 2046 /* Maximum size of input buffer */

char basename_[FILENAME_LEN];
char directory_[FILENAME_LEN];
int blocksize_;
float x_0_, y_0_, z_0_, x_1_, y_1_, z_1_;
float x_inc_=1.0, y_inc_=1.0, z_inc_=1.0;
float sx_, sy_, sz_;
float plane_coeffs1_[3],plane_coeffs2_[3];
int order_, number_of_tilts_;
int skipView[MAX_TILT];
float coefficients_lambda_[MAX_TILT][4];
float coefficients_1_[MAX_TILT][MAX_COEFFICIENTS_NUMBER];
float coefficients_2_[MAX_TILT][MAX_COEFFICIENTS_NUMBER];
float coefficients_3_[MAX_TILT][MAX_COEFFICIENTS_NUMBER];

int read_configuration_file(char* configuration_file_name) {

	FILE *input = fopen(configuration_file_name,"rt");

	long file_length, index;
	char *file_buffer, *pointer;
	char line[MAX_REC_LEN];

	fseek(input, 0L, SEEK_END);
	file_length = ftell(input);
	rewind(input);

	file_buffer = calloc(file_length + 1, sizeof(char));

	if(!file_buffer) {

		txbr_error(stderr, "ERROR: save_block - getting memory.\n");

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
	char format1[MAX_REC_LEN] = "%*s %f";
        char format2[MAX_REC_LEN] = "%*s %i";

	/* sscanf works with float	*/

	sscanf(line,"basename: %s",&basename_[0]);
	sscanf(line,"directory: %s",&directory_[0]);
	sscanf(line,"x->%f:%f:%f",&x_0_,&x_inc_,&x_1_);
	sscanf(line,"y->%f:%f:%f",&y_0_,&y_inc_,&y_1_);
	sscanf(line,"z->%f:%f:%f",&z_0_,&z_inc_,&z_1_);
	sscanf(line,"blocksize: %i",&blocksize_);
        sscanf(line,"scale->%f:%f:%f",&sx_,&sy_,&sz_);
	sscanf(line,"plane_coeffs1: %f %f %f",&plane_coeffs1_[0],&plane_coeffs1_[1],&plane_coeffs1_[2]);
	sscanf(line,"plane_coeffs2: %f %f %f",&plane_coeffs2_[0],&plane_coeffs2_[1],&plane_coeffs2_[2]);
	sscanf(line,"projection approximation order: %i",&order_);
	sscanf(line,"number of tilts: %i",&number_of_tilts_);

	if (strncmp(line, "skip: ", 6)==0) {
            while (sscanf(line, &format2[0], &index) >= 1) {
                strcat(format0, " %*i");
                strcpy(format2, format0);
                strcat(format2, " %i");
                skipView[index] = 1;
            }
        }

	if (sscanf(line, "lambda-%i:", &index_of_tilt) >= 1) {
            index = 0;
            while (sscanf(line, &format1[0], &t) >= 1) {
                strcat(format0, " %*f");
                strcpy(format1, format0);
                strcat(format1, " %f");
                coefficients_lambda_[index_of_tilt][index] = t;
                index++;
            }
        }

	if (sscanf(line,"x-%i:",&index_of_tilt)>=1) {
		index = 0;
		while (sscanf(line,&format1[0],&t)>=1) {
                    strcat(format0," %*f");
                    strcpy(format1,format0);
                    strcat(format1," %f");
                    coefficients_1_[index_of_tilt][index] = t;
                    index++;
		}
	}

	if (sscanf(line,"y-%i:",&index_of_tilt)>=1) {
		index = 0;
		while (sscanf(line,&format1[0],&t)>=1) {
                    strcat(format0," %*f");
                    strcpy(format1,format0);
                    strcat(format1," %f");
                    coefficients_2_[index_of_tilt][index] = t;
                    index++;
		}
	}

}

/*
 * Load TxBRsetup and PathProjection structures from the configuration file
 */
void load_configuration(TxBRsetup* setup, PathProjection* projection, char* configuration_file_name) {

	printf("Load configuration from %s\n",configuration_file_name);

	int i,j;

	init_int_array(&skipView[0],MAX_TILT);
	init_float_array(&coefficients_lambda_[0][0],4*MAX_TILT);
	init_float_array(&coefficients_1_[0][0],MAX_COEFFICIENTS_NUMBER*MAX_TILT);
	init_float_array(&coefficients_2_[0][0],MAX_COEFFICIENTS_NUMBER*MAX_TILT);
	init_float_array(&coefficients_3_[0][0],MAX_COEFFICIENTS_NUMBER*MAX_TILT);

	read_configuration_file(configuration_file_name);

	memcpy(&setup->basename[0],basename_,FILENAME_LEN);
	memcpy(&setup->directory[0],directory_,FILENAME_LEN);

	setup->blocksize = TXBR_MIN(blocksize_,BLOCK_SIZE_MAX);

	setup->x_0 = x_0_;
	setup->x_inc = x_inc_;
	setup->x_1 = x_1_;

	setup->y_0 = y_0_;
	setup->y_inc = y_inc_;
	setup->y_1 = y_1_;

	setup->z_0 = z_0_;
	setup->z_inc = z_inc_;
	setup->z_1 = z_1_;

        setup->sx = sx_;
	setup->sy = sy_;
	setup->sz = sz_;
        
	for (i=0; i<3; i++) {
		setup->plane_coeffs1[i] = plane_coeffs1_[i];
		setup->plane_coeffs2[i] = plane_coeffs2_[i];
	}

	projection->order = order_;
	projection->number_of_terms = number_of_polynomial_terms(order_);
	projection->numberOfTilts = number_of_tilts_;

	int n = number_of_polynomial_terms(order_);

	for (i=0; i<number_of_tilts_; i++) {
		projection->skipView[i] = skipView[i];
		for (j=0; j<4; j++) {
			projection->lambda[i][j] = coefficients_lambda_[i][j];
		}
		for (j=0; j<n; j++) {
			projection->coefficients_1[i][j] = coefficients_1_[i][j];
			projection->coefficients_2[i][j] = coefficients_2_[i][j];
			projection->coefficients_3[i][j] = coefficients_3_[i][j];
		}
	}

}

int write_configuration(TxBRsetup* setup, PathProjection* projection, char* configuration_file_name) {

	FILE *output = fopen(configuration_file_name,"wt");

	fprintf(output,"basename: %s\n",setup->basename);
	fprintf(output,"directory: %s\n",setup->directory);

	fprintf(output,"\n");

	fprintf(output,"x->%f:%f:%f",setup->x_0,setup->x_inc,setup->x_1);
	fprintf(output,"y->%f:%f:%f",setup->y_0,setup->y_inc,setup->y_1);
	fprintf(output,"z->%f:%f:%f",setup->z_0,setup->z_inc,setup->z_1);

	fprintf(output,"\n");

	fprintf(output,"blocksize: %i\n",setup->blocksize);

	fprintf(output,"\n");

	fprintf(output,"plane_coeffs1: %f %f %f\n",setup->plane_coeffs1[0],setup->plane_coeffs1[1],setup->plane_coeffs1[2]);
	fprintf(output,"plane_coeffs2: %f %f %f\n",setup->plane_coeffs2[0],setup->plane_coeffs2[1],setup->plane_coeffs2[2]);

	fprintf(output,"\n");

	fprintf(output,"projection approximation order: %i\n",projection->order);

	fprintf(output,"\n");

	fprintf(output,"number of tilts: %i\n",projection->numberOfTilts);

	fprintf(output,"\n");

	int i,j;

	int n = number_of_polynomial_terms(projection->order);

	for (i=0; i<projection->numberOfTilts; i++) {

		fprintf(output,"lambda-%i: ",i);
		for (j=0; j<4; j++) fprintf(output,"%e ",projection->lambda[i][j]);
		fprintf(output,"\n");

		fprintf(output,"x-%i: ",i);
		for (j=0; j<n; j++) fprintf(output,"%e ",projection->coefficients_1[i][j]);
		fprintf(output,"\n");

		fprintf(output,"y-%i: ",i);
		for (j=0; j<n; j++) fprintf(output,"%e ",projection->coefficients_2[i][j]);
		fprintf(output,"\n");

	}

	fclose(output);

}

