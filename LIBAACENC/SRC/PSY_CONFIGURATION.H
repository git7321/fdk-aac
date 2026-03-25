/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Psychoaccoustic configuration

*******************************************************************************/

#ifndef PSY_CONFIGURATION_H
#define PSY_CONFIGURATION_H

#include "aacenc.h"
#include "common_fix.h"

#include "psy_const.h"
#include "aacenc_tns.h"
#include "aacenc_pns.h"

#define THR_SHIFTBITS 4
#define PCM_QUANT_THR_SCALE 16
#define BITS_PER_LINE_SHIFT 3

#define C_RATIO \
  (FIXP_DBL)0x02940a10

typedef struct {
  INT sfbCnt;
  INT sfbActive;
  INT sfbActiveLFE;
  INT sfbOffset[MAX_SFB + 1];

  INT filterbank;

  FIXP_DBL sfbPcmQuantThreshold[MAX_SFB];

  INT maxAllowedIncreaseFactor;
  FIXP_SGL minRemainingThresholdFactor;

  INT lowpassLine;
  INT lowpassLineLFE;
  FIXP_DBL clipEnergy;

  FIXP_DBL sfbMaskLowFactor[MAX_SFB];
  FIXP_DBL sfbMaskHighFactor[MAX_SFB];

  FIXP_DBL sfbMaskLowFactorSprEn[MAX_SFB];
  FIXP_DBL sfbMaskHighFactorSprEn[MAX_SFB];

  FIXP_DBL sfbMinSnrLdData[MAX_SFB];

  TNS_CONFIG tnsConf;
  PNS_CONFIG pnsConf;

  INT granuleLength;
  INT allowIS;
  INT allowMS;
} PSY_CONFIGURATION;

typedef struct {
  UCHAR sfbCnt;
  UCHAR sfbWidth[MAX_SFB_LONG];
} SFB_PARAM_LONG;

typedef struct {
  UCHAR sfbCnt;
  UCHAR
  sfbWidth[MAX_SFB_SHORT];
} SFB_PARAM_SHORT;

AAC_ENCODER_ERROR FDKaacEnc_InitPsyConfiguration(INT bitrate, INT samplerate,
                                                 INT bandwidth, INT blocktype,
                                                 INT granuleLength, INT useIS,
                                                 INT useMS,
                                                 PSY_CONFIGURATION *psyConf,
                                                 FB_TYPE filterbank);

#endif /* PSY_CONFIGURATION_H */
