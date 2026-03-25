/**************************** AAC encoder library ******************************

   Author(s):   Alex Groeschel

   Description: Temporal noise shaping

*******************************************************************************/

#ifndef AACENC_TNS_H
#define AACENC_TNS_H

#include "common_fix.h"
#include "psy_const.h"

#ifndef PI
#define PI 3.1415926535897931f
#endif

#define TNS_ENABLE_MASK 0xf
#define TNS_MAX_ORDER 12
#define MAX_NUM_OF_FILTERS 2

#define HIFILT 0
#define LOFILT 1

typedef struct {
  INT filterEnabled[MAX_NUM_OF_FILTERS];
  INT threshOn[MAX_NUM_OF_FILTERS];
  INT filterStartFreq[MAX_NUM_OF_FILTERS];
  INT tnsLimitOrder[MAX_NUM_OF_FILTERS];
  INT tnsFilterDirection[MAX_NUM_OF_FILTERS];
  INT acfSplit[MAX_NUM_OF_FILTERS];
  FIXP_DBL tnsTimeResolution[MAX_NUM_OF_FILTERS];
  INT seperateFiltersAllowed;
} TNS_PARAMETER_TABULATED;

typedef struct {
  TNS_PARAMETER_TABULATED confTab;
  INT isLowDelay;
  INT tnsActive;
  INT maxOrder;
  INT coefRes;
  FIXP_DBL acfWindow[MAX_NUM_OF_FILTERS][TNS_MAX_ORDER + 3 + 1];
  INT lpcStartBand[MAX_NUM_OF_FILTERS];
  INT lpcStartLine[MAX_NUM_OF_FILTERS];
  INT lpcStopBand;
  INT lpcStopLine;

} TNS_CONFIG;

typedef struct {
  INT tnsActive[MAX_NUM_OF_FILTERS];
  INT predictionGain[MAX_NUM_OF_FILTERS];
} TNS_SUBBLOCK_INFO;

typedef struct {
  TNS_SUBBLOCK_INFO subBlockInfo[TRANS_FAC];
  FIXP_DBL ratioMultTable[TRANS_FAC][MAX_SFB_SHORT];
} TNS_DATA_SHORT;

typedef struct {
  TNS_SUBBLOCK_INFO subBlockInfo;
  FIXP_DBL ratioMultTable[MAX_SFB_LONG];
} TNS_DATA_LONG;

typedef shouldBeUnion {
  TNS_DATA_LONG Long;
  TNS_DATA_SHORT Short;
}
TNS_DATA_RAW;

typedef struct {
  INT numOfSubblocks;
  TNS_DATA_RAW dataRaw;
  INT tnsMaxScaleSpec;
  INT filtersMerged;
} TNS_DATA;

typedef struct {
  INT numOfFilters[TRANS_FAC];
  INT coefRes[TRANS_FAC];
  INT length[TRANS_FAC][MAX_NUM_OF_FILTERS];
  INT order[TRANS_FAC][MAX_NUM_OF_FILTERS];
  INT direction[TRANS_FAC][MAX_NUM_OF_FILTERS];
  INT coefCompress[TRANS_FAC][MAX_NUM_OF_FILTERS];
  INT coef[TRANS_FAC][MAX_NUM_OF_FILTERS][TNS_MAX_ORDER];
} TNS_INFO;

INT FDKaacEnc_FreqToBandWidthRounding(const INT freq, const INT fs,
                                      const INT numOfBands,
                                      const INT *bandStartOffset);

#endif /* AACENC_TNS_H */
