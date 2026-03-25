/*********************** MPEG surround encoder library *************************

   Author(s):   M. Multrus

   Description: Parameter Extraction

*******************************************************************************/

#ifndef SACENC_PARAMEXTRACT_H
#define SACENC_PARAMEXTRACT_H

/* Includes ******************************************************************/
#include "common_fix.h"
#include "sacenc_lib.h"
#include "sacenc_const.h"
#include "sacenc_bitstream.h"

/* Defines *******************************************************************/
#define MAX_CLD_QUANT_FINE (31)
#define MAX_CLD_QUANT_COARSE (15)
#define OFFSET_CLD_QUANT_COARSE (7)
#define OFFSET_CLD_QUANT_FINE (15)

#define MAX_ICC_QUANT_COARSE (4)
#define MAX_ICC_QUANT_FINE (8)
#define OFFSET_ICC_QUANT_COARSE (0)
#define OFFSET_ICC_QUANT_FINE (0)

#define MAX_NUM_PARAM_BANDS (28)

#define NUM_MAPPED_HYBRID_BANDS (16)

/* Data Types ****************************************************************/
typedef struct T_TTO_BOX *HANDLE_TTO_BOX;

typedef enum {
  BOX_SUBBANDS_INVALID = 0,
  BOX_SUBBANDS_4 = 4,
  BOX_SUBBANDS_5 = 5,
  BOX_SUBBANDS_7 = 7,
  BOX_SUBBANDS_9 = 9,
  BOX_SUBBANDS_12 = 12,
  BOX_SUBBANDS_15 = 15,
  BOX_SUBBANDS_23 = 23

} BOX_SUBBAND_CONFIG;

typedef enum {
  BOX_QUANTMODE_INVALID = -1,
  BOX_QUANTMODE_FINE = 0,
  BOX_QUANTMODE_EBQ1 = 1,
  BOX_QUANTMODE_EBQ2 = 2,
  BOX_QUANTMODE_RESERVED3 = 3,
  BOX_QUANTMODE_RESERVED4 = 4,
  BOX_QUANTMODE_RESERVED5 = 5,
  BOX_QUANTMODE_RESERVED6 = 6,
  BOX_QUANTMODE_RESERVED7 = 7

} BOX_QUANTMODE;

typedef struct T_TTO_BOX_CONFIG {
  UCHAR bUseCoarseQuantCld;
  UCHAR bUseCoarseQuantIcc;
  UCHAR bUseCoherenceIccOnly;

  BOX_SUBBAND_CONFIG subbandConfig;
  BOX_QUANTMODE boxQuantMode;

  UCHAR nHybridBandsMax;

  UCHAR bFrameKeep;

} TTO_BOX_CONFIG;

FDK_SACENC_ERROR fdk_sacenc_createTtoBox(HANDLE_TTO_BOX *hTtoBox);

FDK_SACENC_ERROR fdk_sacenc_initTtoBox(HANDLE_TTO_BOX hTtoBox,
                                       const TTO_BOX_CONFIG *const ttoBoxConfig,
                                       UCHAR *pParameterBand2HybridBandOffset);

FDK_SACENC_ERROR fdk_sacenc_destroyTtoBox(HANDLE_TTO_BOX *hTtoBox);

FDK_SACENC_ERROR fdk_sacenc_applyTtoBox(
    HANDLE_TTO_BOX hTtoBox, const INT nTimeSlots, const INT startTimeSlot,
    const INT nHybridBands, const FIXP_DPK *const *const ppHybridData1__FDK,
    const FIXP_DPK *const *const ppHybridData2__FDK, SCHAR *const pIccIdx,
    UCHAR *const pbIccQuantCoarse, SCHAR *const pCldIdx,
    UCHAR *const pbCldQuantCoarse, const INT bUseBBCues, INT *scaleCh0,
    INT *scaleCh1);

INT fdk_sacenc_subband2ParamBand(const BOX_SUBBAND_CONFIG boxSubbandConfig,
                                 const INT nSubband);

const INT *fdk_sacenc_getSubbandImagSign();

void fdk_sacenc_calcParameterBand2HybridBandOffset(
    const BOX_SUBBAND_CONFIG subbandConfig, const INT nHybridBands,
    UCHAR *pParameterBand2HybridBandOffset);

/* Function / Class Definition ***********************************************/
static inline UCHAR fdk_sacenc_getCldQuantOffset(const INT bUseCoarseQuant) {
  return ((bUseCoarseQuant) ? OFFSET_CLD_QUANT_COARSE : OFFSET_CLD_QUANT_FINE);
}
static inline UCHAR fdk_sacenc_getIccQuantOffset(const INT bUseCoarseQuant) {
  return ((bUseCoarseQuant) ? OFFSET_ICC_QUANT_COARSE : OFFSET_ICC_QUANT_FINE);
}

static inline UCHAR fdk_sacenc_getNumberCldQuantLevels(
    const INT bUseCoarseQuant) {
  return ((bUseCoarseQuant) ? MAX_CLD_QUANT_COARSE : MAX_CLD_QUANT_FINE);
}
static inline UCHAR fdk_sacenc_getNumberIccQuantLevels(
    const INT bUseCoarseQuant) {
  return ((bUseCoarseQuant) ? MAX_ICC_QUANT_COARSE : MAX_ICC_QUANT_FINE);
}

#endif /* SACENC_PARAMEXTRACT_H */
