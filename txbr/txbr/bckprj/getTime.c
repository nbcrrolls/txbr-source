
/********************************************************************************************

	Description		:	Only for Win32. Emulates the UNIX gettimeofday() function


	Parameters		:	same as the unix gettimeofday()
	
	Returns			:	same as the unix gettimeofday()
********************************************************************************************/

#include "getTime.h"

#ifdef WIN32

#include <sys\timeb.h>
#include <math.h>
#include <winsock2.h>

// Note that tzp is going to be ignored in this function and only milliseconds will
// be reported instead of microseconds.
int gettimeofday(struct timeval *tp, void *tzp)
{
    struct _timeb t;

    _ftime(&t);
	tp->tv_sec = t.time;
	tp->tv_usec = t.millitm * 1000;

    return 0;
}

#endif



/********************************************************************************************

	Description		:	The function returns the current time in microseconds. Its used for
					the various performance measurement parts of the code.


	Parameters		:	None.
	
	Returns			:	The time in millisecs
********************************************************************************************/

double getTime()
{
	struct timeval t;
	double timeStamp;
	gettimeofday(&t, NULL);	

	timeStamp = (double)t.tv_sec * 1000000.0 + (double)t.tv_usec;

	return ( timeStamp );
}


