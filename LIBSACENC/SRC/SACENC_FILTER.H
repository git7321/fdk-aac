/*********************** MPEG surround encoder library *************************

   Author(s):   M. Multrus

   Description: Encoder Library Interface
                Filter functions

*******************************************************************************/

#ifndef SACENC_FILTER_H
#define SACENC_FILTER_H

/* Includes ******************************************************************/
#include "common_fix.h"
#include "sacenc_lib.h"
#include "FDK_matrixCalloc.h"

/* Defines *******************************************************************/
#define DC_FILTER_SF 1

/* Data Types ****************************************************************/
typedef struct T_DC_FILTER *HANDLE_DC_FILTER;

FDK_SACENC_ERROR fdk_sacenc_createDCFilter(HANDLE_DC_FILTER *hDCFilter);

FDK_SACENC_ERROR fdk_sacenc_initDCFilter(HANDLE_DC_FILTER hDCFilter,
                                         const UINT sampleRate);

FDK_SACENC_ERROR fdk_sacenc_destroyDCFilter(HANDLE_DC_FILTER *hDCFilter);

FDK_SACENC_ERROR fdk_sacenc_applyDCFilter(HANDLE_DC_FILTER hDCFilter,
                                          const INT_PCM *const signalIn,
                                          INT_PCM *const signalOut,
                                          const INT signalLength);

#endif /* SACENC_FILTER_H */
