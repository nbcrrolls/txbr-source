#include "soapH.h"
#include "txbr.nsmap"

//#define HOST "jane.crbs.ucsd.edu"
//#define PORT 7654
#define HOST "gelfand.crbs.ucsd.edu"
#define PORT 7654

int main()
{
   struct soap soap;
   int m, s; // master and slave sockets
   soap_init(&soap);
   m = soap_bind(&soap, HOST, PORT, 100);
   if (m < 0)
      soap_print_fault(&soap, stderr);
   else
   {
      fprintf(stderr, "Socket connection successful: master socket = %d\n", m);
      for (int i = 1; ; i++)
      {
         s = soap_accept(&soap);
         if (s < 0)
         {
            soap_print_fault(&soap, stderr);
            break;
         }
         fprintf(stderr, "%d: accepted connection from IP= %d.%d.%d.%d socket= %d", i,
             (soap.ip >> 24)&0xFF, (soap.ip >> 16)&0xFF, (soap.ip >> 8)&0xFF, soap.ip&0xFF, s);
         if (soap_serve(&soap) != SOAP_OK) // process RPC request
            soap_print_fault(&soap, stderr); // print error
         fprintf(stderr, "request served\n");
         soap_destroy(&soap); // clean up class instances
         soap_end(&soap); // clean up everything and close socket
      }
   }
   soap_done(&soap); // close master socket and detach context
} 


int txbr__runOFFTxbr( struct soap *soap, char* d, char* wd, char* basename, enum FILE_TYPE type, enum SCOPE scope, int& return_value)
{

    fprintf(stderr, "\nin the correct method...\n");
    fprintf(stderr, "\ngot this from client for path: %s and for basename %s\n", d, basename);

    std::string txbr_command = "txbr_onfly.py -d " + std::string(d) + 
                                            " -b " + std::string(basename) +
                                          " --wd " + std::string(wd);
    
    if ( type==JPG) {
        txbr_command +=  " --extension jpg";
    } else if ( type==TIF) {
        txbr_command +=  " --extension tif";
    }

    if ( scope==J3200) {
        txbr_command +=  " --scope J3200";
    } else if ( scope==J4000_1) {
        txbr_command +=  " --scope J4000_1";
    } else if ( scope==J4000_2) {
        txbr_command +=  " --scope J4000_2";
    } else if ( scope==FEI_titan) {
        txbr_command +=  " --scope FEI_titan";
    } else if ( scope==FEI_spirit) {
        txbr_command +=  " --scope FEI_spirit";
    }

    const char *command = txbr_command.c_str();

    fprintf(stderr, "\n%s\n", command);
    std::cout << command << "\n\n";

    return_value = system( command );

    return SOAP_OK;

}


int txbr__resetOFFTxbr( struct soap *soap, char* wd, int &return_value )
{

    fprintf( stderr, "\nin the correct method...\n" );
    fprintf( stderr, "\ngot this from client for path: %s\n", wd );

    std::string txbr_command = "txbr_onfly.py --wd " + std::string(wd) + " --reset ";

    const char *command = txbr_command.c_str();

    fprintf(stderr, "\n%s\n", command);
    std::cout << command << "\n\n";

    return_value = system( command );

    return SOAP_OK;

}


int txbr__statOFFTxbr( struct soap *soap, char* wd, int &return_value )
{

    fprintf( stderr, "\nin the correct method...\n" );
    fprintf( stderr, "\ngot this from client for path: %s\n", wd );

    std::string txbr_command = "txbr_onfly.py --wd " + std::string(wd) + " --status ";

    const char *command = txbr_command.c_str();

    fprintf(stderr, "\n%s\n", command);
    std::cout << command << "\n\n";

    return_value = system( command );

    return SOAP_OK;

}

