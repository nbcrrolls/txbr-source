// dllexport.h : Include this file to define DLL export or import for Windows
// define EXPORT_FOR_DLL to do export instead of import
/*  $Author: mast $

$Date: 2004/12/22 15:24:18 $

$Revision: 4ddbc87b5613 $

$Log$
Revision 1.1  2004/06/04 02:55:43  mast
Creation for making libdiaqt a DLL

*/
#include "imodconfig.h"

#ifndef DLLEXPORT_H
#define DLLEXPORT_H

// Note that DLL_EX_IM is used for exporting functions from 3dmod so we use
// DLL_IM_EX for library exports to avoid conflict
#ifdef _WIN32
#ifdef EXPORT_FOR_DLL
#define DLL_IM_EX _declspec(dllexport)
#else
#define DLL_IM_EX _declspec(dllimport)
#endif
#else
#define DLL_IM_EX
#endif

#endif
