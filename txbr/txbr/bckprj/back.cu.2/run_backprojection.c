#include <stdio.h>
#include <time.h>
#include "txbr.h"
#include "txbrutil.h"

// global variables
extern int g_argc;
extern char **g_argv;

/* To reconstruct a volume */
void reconstruct(char* directory, char* basename, char* work_directory, int order, int config, 
	double zstart, double zinc, double zstop, int override_z, int gpuDevId) {

	TxBRsetup *setup;
	PathProjection *projection;

	printf("Read configuration from a text file\n");

	setup  = (TxBRsetup *)malloc(sizeof(TxBRsetup));

	if (!setup) {

		txbr_error(stderr, "ERROR: setup - getting memory.\n");

		return;

	}

	projection  = (PathProjection*)malloc(sizeof(PathProjection));

	if (!projection) {

		txbr_error(stderr, "ERROR: projection - getting memory.\n");

		return;

	}

	char configuration_file_name[FILENAME_LEN];
	sprintf(configuration_file_name, "%s/%s.txbr", directory, basename);

	load_configuration(setup, projection, configuration_file_name);

	/* Print out informations	*/

	print_TxBR_setup(setup);
	printf("\n");

	print_path_projection(projection);
	printf("\n");

	/* Give the boundaries condition	*/

	if (override_z==1) {
		setup->z_0 = zstart;
		setup->z_1 = zstop;
	}

	/* Run the reconstruction routine */

	do_full_reconstruction(directory, basename, work_directory, setup, projection,
		gpuDevId);

	free(projection);

	free(setup);

}


int main(int argc, char *argv[]) {

	int i;

	char *directory, *basename, *work_directory=NULL;
	double zstart = 0, zstop = 0, zinc = 1;
	int override_z = 0;
	int order = 1, config = 1, blocksize = 5;
	int gpuDevId = 0;
	
	g_argc = argc;
	g_argv = argv;

	for (i=1; i<argc; i++){

		if (argv[i][0]=='-') {

			switch (argv[i][1]) {

				case 'b' :

					if (strcmp(&argv[i][0],"-b")==0) basename = argv[++i];
					if (strcmp(&argv[i][0],"-block")==0) blocksize = atoi(argv[++i]);
					break;

				case 'c' :

					config = atoi(argv[++i]);
					break;

				case 'd' :
					directory = argv[++i];
					break;

				case 'o' :

					order = atoi(argv[++i]);
					break;

				case 'w' :

					work_directory = argv[++i];
					break;
				
				case 'g' :
					gpuDevId = atoi(argv[++i]);
					break;
					
				case 'z' :

					override_z = 1;
					if (strcmp(&argv[i][0],"-zstart")==0) zstart = atof(argv[++i]);
					if (strcmp(&argv[i][0],"-zinc")==0) zinc = atof(argv[++i]);
					if (strcmp(&argv[i][0],"-zstop")==0) zstop = atof(argv[++i]);
					break;

			}
		}

	}

	if (!work_directory) {

		work_directory = (char*)malloc(FILENAME_LEN*sizeof(char));

		if (!work_directory) {
			txbr_error(stderr, "ERROR: main - getting memory.\n");
			return MEMORY_ERROR;
		}

		sprintf(&work_directory[0],".");

	}

	printf("directory: %s\n", &directory[0]);
	printf("basename: %s\n", &basename[0]);
	printf("work_directory: %s\n", &work_directory[0]);
	printf("alignment order: %i\n", order);

	clock_t t_0,t_1;
	float ratio = 1./CLOCKS_PER_SEC;

	t_0 = clock();

	reconstruct(directory, basename, work_directory, order, config, 
		zstart, zinc, zstop, override_z, gpuDevId);

	t_1 = clock();

	printf("Total Durations:%E\n",ratio*(long)t_1-ratio*(long)t_0);

}
