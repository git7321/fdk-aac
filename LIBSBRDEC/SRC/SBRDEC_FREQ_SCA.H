/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef SBRDEC_FREQ_SCA_H
#define SBRDEC_FREQ_SCA_H

#include "sbrdecoder.h"
#include "env_extr.h"

typedef enum { DUAL, QUAD } SBR_RATE;

SBR_ERROR
sbrdecUpdateFreqScale(UCHAR *v_k_master, UCHAR *numMaster, UINT fs,
                      HANDLE_SBR_HEADER_DATA headerData, UINT flags);

void sbrdecDownSampleLoRes(UCHAR *v_result, UCHAR num_result,
                           UCHAR *freqBandTableRef, UCHAR num_Ref);

void shellsort(UCHAR *in, UCHAR n);

SBR_ERROR
resetFreqBandTables(HANDLE_SBR_HEADER_DATA hHeaderData, const UINT flags);

#endif
