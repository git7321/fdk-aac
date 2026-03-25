/*********************** MPEG surround encoder library *************************

   Author(s):   Max Neuendorf

   Description: Encoder Library Interface
                Determine TES flags

*******************************************************************************/

#ifndef SACENC_ONSETDETECT_H
#define SACENC_ONSETDETECT_H

/**************************************************************************/ /**
   \file
   Description of file contents
 ******************************************************************************/

/* Includes ******************************************************************/
#include "common_fix.h"
#include "FDK_matrixCalloc.h"
#include "sacenc_lib.h"
#include "sacenc_bitstream.h"

#define MAX_NUM_TRANS (MAX_NUM_PARAMS / 2)

typedef struct T_ONSET_DETECT_CONFIG {
  INT maxTimeSlots;

  INT lowerBoundOnsetDetection;
  INT upperBoundOnsetDetection;

} ONSET_DETECT_CONFIG;

typedef struct ONSET_DETECT *HANDLE_ONSET_DETECT;

FDK_SACENC_ERROR fdk_sacenc_onsetDetect_Open(HANDLE_ONSET_DETECT *phOnset,
                                             const UINT maxTimeSlots);

FDK_SACENC_ERROR fdk_sacenc_onsetDetect_Init(
    HANDLE_ONSET_DETECT hOnset,
    const ONSET_DETECT_CONFIG *const pOnsetDetectConfig, const UINT initFlags);

FDK_SACENC_ERROR fdk_sacenc_onsetDetect_Close(HANDLE_ONSET_DETECT *phOnset);

FDK_SACENC_ERROR fdk_sacenc_onsetDetect_Update(HANDLE_ONSET_DETECT hOnset,
                                               const INT timeSlots);

FDK_SACENC_ERROR fdk_sacenc_onsetDetect_Apply(
    HANDLE_ONSET_DETECT hOnset, const INT nTimeSlots, const INT nHybridBands,
    FIXP_DPK *const *const ppHybridData__FDK, const INT hybridDataScale,
    const INT prevPos, INT pTransientPos[MAX_NUM_TRANS]);

#endif /* SACENC_ONSETDETECT_H */
