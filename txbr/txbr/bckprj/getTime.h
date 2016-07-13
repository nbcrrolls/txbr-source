
 /****************************************************************************************
 *
 *	Author	:		Raj Singh
 *	email	:		rajvikrams@gmail.com
 *
 *	Status	;		Experimental
 *	Version :	
 *
 *	Description:	Cross platform time function 
 *
 *	Notes	: 	
 *			
 * 
******************************************************************************************/

#ifndef _PTGETTIME_H
#define _PTGETTIME_H

#include <stdlib.h>

#ifdef WIN32
	#include <time.h>
	
	int gettimeofday(struct timeval *tp, void *tzp);

#else
	#include <sys/time.h>
	#include <unistd.h>
	#include <sys/types.h>
	
#endif

double getTime();

#endif 
