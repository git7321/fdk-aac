/**************************** AAC encoder library ******************************

   Author(s):   M. Schug / A. Groeschel

   Description: threshold calculations

*******************************************************************************/

#ifndef ADJ_THR_DATA_H
#define ADJ_THR_DATA_H

#include "psy_const.h"

typedef enum {
  AACENC_BD_MODE_INTER_ELEMENT = 0,
  AACENC_BD_MODE_INTRA_ELEMENT = 1
} AACENC_BIT_DISTRIBUTION_MODE;

typedef enum {
  AACENC_BR_MODE_FULL = 0,
  AACENC_BR_MODE_REDUCED = 1,
  AACENC_BR_MODE_DISABLED = 2
} AACENC_BITRES_MODE;

typedef struct {
  FIXP_DBL clipSaveLow, clipSaveHigh;
  FIXP_DBL minBitSave, maxBitSave;
  FIXP_DBL clipSpendLow, clipSpendHigh;
  FIXP_DBL minBitSpend, maxBitSpend;
} BRES_PARAM;

typedef struct {
  INT modifyMinSnr;
  INT startSfbL, startSfbS;
} AH_PARAM;

typedef struct {
  FIXP_DBL maxRed;
  FIXP_DBL startRatio;
  FIXP_DBL maxRatio;
  FIXP_DBL redRatioFac;
  FIXP_DBL redOffs;
} MINSNR_ADAPT_PARAM;

typedef struct {
  INT peMin, peMax;
  INT peOffset;
  FIXP_DBL bits2PeFactor_m;
  INT bits2PeFactor_e;
  AH_PARAM ahParam;
  MINSNR_ADAPT_PARAM minSnrAdaptParam;

  INT peLast;
  INT dynBitsLast;
  FIXP_DBL peCorrectionFactor_m;
  INT peCorrectionFactor_e;

  FIXP_DBL vbrQualFactor;
  FIXP_DBL chaosMeasureOld;

  FIXP_DBL chaosMeasureEnFac[(2)];
  INT lastEnFacPatch[(2)];

} ATS_ELEMENT;

typedef struct {
  BRES_PARAM bresParamLong, bresParamShort;
  ATS_ELEMENT* adjThrStateElem[((8))];
  AACENC_BIT_DISTRIBUTION_MODE bitDistributionMode;
  INT maxIter2ndGuess;
} ADJ_THR_STATE;

#endif
