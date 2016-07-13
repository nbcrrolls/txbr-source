#include "iostream.h"
#include "soapH.h"
#include "txbr.nsmap"

static const int RUN_MODE = 0;
static const int RESET_MODE = 1;

int main(int argc, char** argv) {

    int return_value;
    int mode = RUN_MODE;

    // Initialize the soap structure for invoking calls.
    struct soap soap;
    soap_init(&soap);

    if (argc>1 && strcmp(argv[1],"reset")==0) {
        mode = RESET_MODE;
    }

    printf("About to make the soap call\n");

    char* host = "http://jane.crbs.ucsd.edu:7654";
    char* src_dir = "/ncmirdata5/sphan/work-fiducialless/on_the_fly/onthefly1";
    char* work_dir = "/ncmirdata5/scopes/onthefly1";
    char* basename = "zap";
    enum FILE_TYPE type = JPG;
    enum SCOPE scope = FEI_titan;
    
    switch (mode) {

        case RUN_MODE:
            if ( soap_call_txbr__runOFFTxbr(&soap, host, "", src_dir, work_dir, basename, type, scope, return_value) == SOAP_OK )
                printf("The service call returned");
            else
                printf("FAILED...");
            break;

        case RESET_MODE:
            if ( soap_call_txbr__resetOFFTxbr(&soap, host, "", work_dir, return_value) == SOAP_OK )
                printf("The service call returned");
            else
                printf("FAILED...");
            break;
            
    }

    printf("\ndone with call\n");

    // clean up the soap environment
    //most likely this bit of code would be in a destructor of a C++ class....
    //soap_destroy(&soap);
    soap_end(&soap);
    soap_done(&soap);

    return 0;

}

