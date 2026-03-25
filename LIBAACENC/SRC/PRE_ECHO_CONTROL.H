/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Pre echo control

*******************************************************************************/

#ifndef PRE_ECHO_CONTROL_H
#define PRE_ECHO_CONTROL_H

#include "common_fix.h"

void FDKaacEnc_InitPreEchoControl(FIXP_DBL *pbThresholdnm1, INT *calcPreEcho,
                                  INT numPb, FIXP_DBL *sfbPcmQuantThreshold,
                                  INT *mdctScalenm1);

void FDKaacEnc_PreEchoControl(FIXP_DBL *pbThresholdNm1, INT calcPreEcho,
                              INT numPb, INT maxAllowedIncreaseFactor,
                              FIXP_SGL minRemainingThresholdFactor,
                              FIXP_DBL *pbThreshold, INT mdctScale,
                              INT *mdctScalenm1);

#endif
