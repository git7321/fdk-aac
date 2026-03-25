/*********************** MPEG surround encoder library *************************

   Author(s):   M. Multrus

   Description: Parameter Extraction

*******************************************************************************/

#include "sacenc_paramextract.h"
#include "sacenc_tree.h"
#include "sacenc_vectorfunctions.h"

#define LOG10_2_10 (3.01029995664f)
#define SCALE_CLDE_SF (7)
#define SCALE_CLDD_SF (8)

typedef struct T_TTO_BOX {
  FIXP_DBL pCld__FDK[MAX_NUM_PARAM_BANDS];
  FIXP_DBL pIcc__FDK[MAX_NUM_PARAM_BANDS];
  FIXP_DBL pCldQuant__FDK[MAX_NUM_PARAM_BANDS];

  const FIXP_DBL *pIccQuantTable__FDK;
  const FIXP_DBL *pCldQuantTableDec__FDK;
  const FIXP_DBL *pCldQuantTableEnc__FDK;

  SCHAR pCldEbQIdx[MAX_NUM_PARAM_BANDS];
  SCHAR pIccDownmixIdx[MAX_NUM_PARAM_BANDS];

  UCHAR *pParameterBand2HybridBandOffset;
  const INT *pSubbandImagSign;
  UCHAR nHybridBandsMax;
  UCHAR nParameterBands;
  UCHAR bFrameKeep;

  UCHAR iccCorrelationCoherenceBorder;
  BOX_QUANTMODE boxQuantMode;

  UCHAR nIccQuantSteps;
  UCHAR nIccQuantOffset;

  UCHAR nCldQuantSteps;
  UCHAR nCldQuantOffset;

  UCHAR bUseCoarseQuantCld;
  UCHAR bUseCoarseQuantIcc;

} TTO_BOX;

struct BOX_SUBBAND_SETUP {
  BOX_SUBBAND_CONFIG subbandConfig;
  UCHAR nParameterBands;
  const UCHAR *pSubband2ParameterIndexLd;
  UCHAR iccCorrelationCoherenceBorder;
};

static const UCHAR subband2Parameter4_Ld[NUM_QMF_BANDS] = {
    0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

static const UCHAR subband2Parameter5_Ld[NUM_QMF_BANDS] = {
    0, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

static const UCHAR subband2Parameter7_Ld[NUM_QMF_BANDS] = {
    0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6};

static const UCHAR subband2Parameter9_Ld[NUM_QMF_BANDS] = {
    0, 1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};

static const UCHAR subband2Parameter12_Ld[NUM_QMF_BANDS] = {
    0,  1,  2,  3,  4,  4,  5,  5,  6,  6,  6,  7,  7,  7,  8,  8,
    8,  8,  9,  9,  9,  9,  9,  10, 10, 10, 10, 10, 10, 10, 10, 10,
    10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
    11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11};

static const UCHAR subband2Parameter15_Ld[NUM_QMF_BANDS] = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  9,  10, 10, 10, 11, 11,
    11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
    14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14};

static const UCHAR subband2Parameter23_Ld[NUM_QMF_BANDS] = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 12, 13, 13,
    14, 14, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 18, 19, 19,
    19, 19, 19, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21,
    22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22};

static const INT subbandImagSign_Ld[NUM_QMF_BANDS] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
};

#define SCALE_CLDE(a) (FL2FXCONST_DBL(a / (float)(1 << SCALE_CLDE_SF)))
static const FIXP_DBL cldQuantTableFineEnc__FDK[MAX_CLD_QUANT_FINE] = {
    SCALE_CLDE(-50.0), SCALE_CLDE(-45.0), SCALE_CLDE(-40.0), SCALE_CLDE(-35.0),
    SCALE_CLDE(-30.0), SCALE_CLDE(-25.0), SCALE_CLDE(-22.0), SCALE_CLDE(-19.0),
    SCALE_CLDE(-16.0), SCALE_CLDE(-13.0), SCALE_CLDE(-10.0), SCALE_CLDE(-8.0),
    SCALE_CLDE(-6.0),  SCALE_CLDE(-4.0),  SCALE_CLDE(-2.0),  SCALE_CLDE(0.0),
    SCALE_CLDE(2.0),   SCALE_CLDE(4.0),   SCALE_CLDE(6.0),   SCALE_CLDE(8.0),
    SCALE_CLDE(10.0),  SCALE_CLDE(13.0),  SCALE_CLDE(16.0),  SCALE_CLDE(19.0),
    SCALE_CLDE(22.0),  SCALE_CLDE(25.0),  SCALE_CLDE(30.0),  SCALE_CLDE(35.0),
    SCALE_CLDE(40.0),  SCALE_CLDE(45.0),  SCALE_CLDE(50.0)};

static const FIXP_DBL cldQuantTableCoarseEnc__FDK[MAX_CLD_QUANT_COARSE] = {
    SCALE_CLDE(-50.0), SCALE_CLDE(-35.0), SCALE_CLDE(-25.0), SCALE_CLDE(-19.0),
    SCALE_CLDE(-13.0), SCALE_CLDE(-8.0),  SCALE_CLDE(-4.0),  SCALE_CLDE(0.0),
    SCALE_CLDE(4.0),   SCALE_CLDE(8.0),   SCALE_CLDE(13.0),  SCALE_CLDE(19.0),
    SCALE_CLDE(25.0),  SCALE_CLDE(35.0),  SCALE_CLDE(50.0)};

#define SCALE_CLDD(a) (FL2FXCONST_DBL(a / (float)(1 << SCALE_CLDD_SF)))
static const FIXP_DBL cldQuantTableFineDec__FDK[MAX_CLD_QUANT_FINE] = {
    SCALE_CLDD(-150.0), SCALE_CLDD(-45.0), SCALE_CLDD(-40.0), SCALE_CLDD(-35.0),
    SCALE_CLDD(-30.0),  SCALE_CLDD(-25.0), SCALE_CLDD(-22.0), SCALE_CLDD(-19.0),
    SCALE_CLDD(-16.0),  SCALE_CLDD(-13.0), SCALE_CLDD(-10.0), SCALE_CLDD(-8.0),
    SCALE_CLDD(-6.0),   SCALE_CLDD(-4.0),  SCALE_CLDD(-2.0),  SCALE_CLDD(0.0),
    SCALE_CLDD(2.0),    SCALE_CLDD(4.0),   SCALE_CLDD(6.0),   SCALE_CLDD(8.0),
    SCALE_CLDD(10.0),   SCALE_CLDD(13.0),  SCALE_CLDD(16.0),  SCALE_CLDD(19.0),
    SCALE_CLDD(22.0),   SCALE_CLDD(25.0),  SCALE_CLDD(30.0),  SCALE_CLDD(35.0),
    SCALE_CLDD(40.0),   SCALE_CLDD(45.0),  SCALE_CLDD(150.0)};

static const FIXP_DBL cldQuantTableCoarseDec__FDK[MAX_CLD_QUANT_COARSE] = {
    SCALE_CLDD(-150.0), SCALE_CLDD(-35.0), SCALE_CLDD(-25.0), SCALE_CLDD(-19.0),
    SCALE_CLDD(-13.0),  SCALE_CLDD(-8.0),  SCALE_CLDD(-4.0),  SCALE_CLDD(0.0),
    SCALE_CLDD(4.0),    SCALE_CLDD(8.0),   SCALE_CLDD(13.0),  SCALE_CLDD(19.0),
    SCALE_CLDD(25.0),   SCALE_CLDD(35.0),  SCALE_CLDD(150.0)};

#define SCALE_ICC(a) (FL2FXCONST_DBL(a))
static const FIXP_DBL iccQuantTableFine__FDK[MAX_ICC_QUANT_FINE] = {
    SCALE_ICC(0.99999999953), SCALE_ICC(0.937f),   SCALE_ICC(0.84118f),
    SCALE_ICC(0.60092f),      SCALE_ICC(0.36764f), SCALE_ICC(0.0f),
    SCALE_ICC(-0.589f),       SCALE_ICC(-0.99f)};

static const FIXP_DBL iccQuantTableCoarse__FDK[MAX_ICC_QUANT_COARSE] = {
    SCALE_ICC(0.99999999953), SCALE_ICC(0.84118f), SCALE_ICC(0.36764f),
    SCALE_ICC(-0.5890f)};

static const BOX_SUBBAND_SETUP boxSubbandSetup[] = {
    {BOX_SUBBANDS_4, 4, subband2Parameter4_Ld, 1},
    {BOX_SUBBANDS_5, 5, subband2Parameter5_Ld, 2},
    {BOX_SUBBANDS_7, 7, subband2Parameter7_Ld, 3},
    {BOX_SUBBANDS_9, 9, subband2Parameter9_Ld, 4},
    {BOX_SUBBANDS_12, 12, subband2Parameter12_Ld, 4},
    {BOX_SUBBANDS_15, 15, subband2Parameter15_Ld, 5},
    {BOX_SUBBANDS_23, 23, subband2Parameter23_Ld, 8}};

static const BOX_SUBBAND_SETUP *getBoxSubbandSetup(
    const BOX_SUBBAND_CONFIG subbandConfig) {
  int i;
  const BOX_SUBBAND_SETUP *setup = NULL;

  for (i = 0; i < (int)(sizeof(boxSubbandSetup) / sizeof(BOX_SUBBAND_SETUP));
       i++) {
    if (boxSubbandSetup[i].subbandConfig == subbandConfig) {
      setup = &boxSubbandSetup[i];
      break;
    }
  }
  return setup;
}

static inline void ApplyBBCuesFDK(FIXP_DBL *const pData,
                                  const INT nParamBands) {
  int i, s;
  FIXP_DBL tmp, invParamBands;

  invParamBands = fDivNormHighPrec((FIXP_DBL)1, (FIXP_DBL)nParamBands, &s);
  s = -s;

  tmp = fMult(pData[0], invParamBands) >> s;
  for (i = 1; i < nParamBands; i++) {
    tmp += fMult(pData[i], invParamBands) >> s;
  }

  for (i = 0; i < nParamBands; i++) {
    pData[i] = tmp;
  }
}

static INT getNumberParameterBands(const BOX_SUBBAND_CONFIG subbandConfig) {
  const BOX_SUBBAND_SETUP *setup = getBoxSubbandSetup(subbandConfig);
  return ((setup == NULL) ? 0 : setup->nParameterBands);
}

static const UCHAR *getSubband2ParameterIndex(
    const BOX_SUBBAND_CONFIG subbandConfig) {
  const BOX_SUBBAND_SETUP *setup = getBoxSubbandSetup(subbandConfig);

  return ((setup == NULL) ? NULL : (setup->pSubband2ParameterIndexLd));
}

void fdk_sacenc_calcParameterBand2HybridBandOffset(
    const BOX_SUBBAND_CONFIG subbandConfig, const INT nHybridBands,
    UCHAR *pParameterBand2HybridBandOffset) {
  const BOX_SUBBAND_SETUP *setup = getBoxSubbandSetup(subbandConfig);
  const UCHAR *pSubband2ParameterIndex;

  int i, pb;

  pSubband2ParameterIndex = setup->pSubband2ParameterIndexLd;

  for (pb = 0, i = 0; i < nHybridBands - 1; i++) {
    if (pSubband2ParameterIndex[i + 1] - pSubband2ParameterIndex[i]) {
      pParameterBand2HybridBandOffset[pb++] = (i + 1);
    }
  }
  pParameterBand2HybridBandOffset[pb++] = (i + 1);
}

const INT *fdk_sacenc_getSubbandImagSign() {
  const INT *pImagSign = NULL;

  pImagSign = subbandImagSign_Ld;

  return (pImagSign);
}

static INT getIccCorrelationCoherenceBorder(
    const BOX_SUBBAND_CONFIG subbandConfig, const INT bUseCoherenceOnly) {
  const BOX_SUBBAND_SETUP *setup = getBoxSubbandSetup(subbandConfig);
  return (
      (setup == NULL)
          ? 0
          : ((bUseCoherenceOnly) ? 0 : setup->iccCorrelationCoherenceBorder));
}

FDK_SACENC_ERROR fdk_sacenc_createTtoBox(HANDLE_TTO_BOX *hTtoBox) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == hTtoBox) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDK_ALLOCATE_MEMORY_1D(*hTtoBox, 1, TTO_BOX);
  }
  return error;

bail:
  fdk_sacenc_destroyTtoBox(hTtoBox);
  return ((SACENC_OK == error) ? SACENC_MEMORY_ERROR : error);
}

FDK_SACENC_ERROR fdk_sacenc_initTtoBox(HANDLE_TTO_BOX hTtoBox,
                                       const TTO_BOX_CONFIG *const ttoBoxConfig,
                                       UCHAR *pParameterBand2HybridBandOffset) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((hTtoBox == NULL) || (ttoBoxConfig == NULL) ||
      (pParameterBand2HybridBandOffset == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDKmemclear(hTtoBox, sizeof(TTO_BOX));

    hTtoBox->bUseCoarseQuantCld = ttoBoxConfig->bUseCoarseQuantCld;
    hTtoBox->bUseCoarseQuantIcc = ttoBoxConfig->bUseCoarseQuantIcc;
    hTtoBox->boxQuantMode = ttoBoxConfig->boxQuantMode;
    hTtoBox->iccCorrelationCoherenceBorder = getIccCorrelationCoherenceBorder(
        ttoBoxConfig->subbandConfig, ttoBoxConfig->bUseCoherenceIccOnly);
    hTtoBox->nHybridBandsMax = ttoBoxConfig->nHybridBandsMax;
    hTtoBox->nParameterBands =
        getNumberParameterBands(ttoBoxConfig->subbandConfig);
    hTtoBox->bFrameKeep = ttoBoxConfig->bFrameKeep;

    hTtoBox->nIccQuantSteps =
        fdk_sacenc_getNumberIccQuantLevels(hTtoBox->bUseCoarseQuantIcc);
    hTtoBox->nIccQuantOffset =
        fdk_sacenc_getIccQuantOffset(hTtoBox->bUseCoarseQuantIcc);

    hTtoBox->pIccQuantTable__FDK = hTtoBox->bUseCoarseQuantIcc
                                       ? iccQuantTableCoarse__FDK
                                       : iccQuantTableFine__FDK;
    hTtoBox->pCldQuantTableDec__FDK = hTtoBox->bUseCoarseQuantCld
                                          ? cldQuantTableCoarseDec__FDK
                                          : cldQuantTableFineDec__FDK;
    hTtoBox->pCldQuantTableEnc__FDK = hTtoBox->bUseCoarseQuantCld
                                          ? cldQuantTableCoarseEnc__FDK
                                          : cldQuantTableFineEnc__FDK;

    hTtoBox->nCldQuantSteps =
        fdk_sacenc_getNumberCldQuantLevels(hTtoBox->bUseCoarseQuantCld);
    hTtoBox->nCldQuantOffset =
        fdk_sacenc_getCldQuantOffset(hTtoBox->bUseCoarseQuantCld);

    if (NULL == (hTtoBox->pParameterBand2HybridBandOffset =
                     pParameterBand2HybridBandOffset)) {
      error = SACENC_INIT_ERROR;
      goto bail;
    }

    if (NULL == (hTtoBox->pSubbandImagSign = fdk_sacenc_getSubbandImagSign())) {
      error = SACENC_INIT_ERROR;
    }

    if ((hTtoBox->boxQuantMode != BOX_QUANTMODE_FINE) &&
        (hTtoBox->boxQuantMode != BOX_QUANTMODE_EBQ1) &&
        (hTtoBox->boxQuantMode != BOX_QUANTMODE_EBQ2)) {
      error = SACENC_INIT_ERROR;
      goto bail;
    }
  }
bail:
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_destroyTtoBox(HANDLE_TTO_BOX *hTtoBox) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (*hTtoBox != NULL) {
    FDKfree(*hTtoBox);
    *hTtoBox = NULL;
  }

  return error;
}

static FDK_SACENC_ERROR calculateIccFDK(const INT nParamBand,
                                        const INT correlationCoherenceBorder,
                                        const FIXP_DBL *const pPwr1,
                                        const FIXP_DBL *const pPwr2,
                                        const FIXP_DBL *const pProdReal,
                                        FIXP_DBL const *const pProdImag,
                                        FIXP_DBL *const pIcc) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((pPwr1 == NULL) || (pPwr2 == NULL) || (pProdReal == NULL) ||
      (pProdImag == NULL) || (pIcc == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    if (correlationCoherenceBorder > nParamBand) {
      error = SACENC_INVALID_CONFIG;
    } else {
      FDKcalcCorrelationVec(pIcc, pProdReal, pPwr1, pPwr2,
                            correlationCoherenceBorder);

      calcCoherenceVec(&pIcc[correlationCoherenceBorder],
                       &pProdReal[correlationCoherenceBorder],
                       &pProdImag[correlationCoherenceBorder],
                       &pPwr1[correlationCoherenceBorder],
                       &pPwr2[correlationCoherenceBorder], 0, 0,
                       nParamBand - correlationCoherenceBorder);
    }
  }

  return error;
}

static void QuantizeCoefFDK(const FIXP_DBL *const input, const INT nBands,
                            const FIXP_DBL *const quantTable,
                            const INT idxOffset, const INT nQuantSteps,
                            SCHAR *const quantOut) {
  int band;
  const int reverse = (quantTable[0] > quantTable[1]);

  for (band = 0; band < nBands; band++) {
    FIXP_DBL qVal;
    FIXP_DBL curVal = input[band];

    int lower = 0;
    int upper = nQuantSteps - 1;

    if (reverse) {
      while (upper - lower > 1) {
        int idx = (lower + upper) >> 1;
        qVal = quantTable[idx];
        if (curVal >= qVal) {
          upper = idx;
        } else {
          lower = idx;
        }
      }

      if ((curVal - quantTable[lower]) >= (quantTable[upper] - curVal)) {
        quantOut[band] = lower - idxOffset;
      } else {
        quantOut[band] = upper - idxOffset;
      }
    }
    else {
      while (upper - lower > 1) {
        int idx = (lower + upper) >> 1;
        qVal = quantTable[idx];
        if (curVal <= qVal) {
          upper = idx;
        } else {
          lower = idx;
        }
      }

      if ((curVal - quantTable[lower]) <= (quantTable[upper] - curVal)) {
        quantOut[band] = lower - idxOffset;
      } else {
        quantOut[band] = upper - idxOffset;
      }
    }
  }
}

static void deQuantizeCoefFDK(const SCHAR *const input, const INT nBands,
                              const FIXP_DBL *const quantTable,
                              const INT idxOffset, FIXP_DBL *const dequantOut) {
  int band;

  for (band = 0; band < nBands; band++) {
    dequantOut[band] = quantTable[input[band] + idxOffset];
  }
}

static void CalculateCldFDK(FIXP_DBL *const pCld, const FIXP_DBL *const pPwr1,
                            const FIXP_DBL *const pPwr2, const INT scaleCh1,
                            const INT *const pbScaleCh1, const INT scaleCh2,
                            const INT *const pbScaleCh2, const int nParamBand) {
  INT i;
  FIXP_DBL ldPwr1, ldPwr2, cld;
  FIXP_DBL maxPwr = FL2FXCONST_DBL(
      30.0f /
      (1 << (LD_DATA_SHIFT +
             1)));

  for (i = 0; i < nParamBand; i++) {
    ldPwr1 =
        (CalcLdData(pPwr1[i]) >> 1) + ((FIXP_DBL)(scaleCh1 + pbScaleCh1[i])
                                       << (DFRACT_BITS - 1 - LD_DATA_SHIFT));
    ldPwr2 =
        (CalcLdData(pPwr2[i]) >> 1) + ((FIXP_DBL)(scaleCh2 + pbScaleCh2[i])
                                       << (DFRACT_BITS - 1 - LD_DATA_SHIFT));

    ldPwr1 = fixMax(fixMin(ldPwr1, maxPwr), -maxPwr);
    ldPwr2 = fixMax(fixMin(ldPwr2, maxPwr), -maxPwr);

    cld = fMultDiv2(FL2FXCONST_DBL(LOG10_2_10 / (1 << SCALE_CLDE_SF)),
                    ldPwr1 - ldPwr2);

    cld =
        fixMin(cld, (FIXP_DBL)(((FIXP_DBL)MAXVAL_DBL) >> (LD_DATA_SHIFT + 2)));
    cld =
        fixMax(cld, (FIXP_DBL)(((FIXP_DBL)MINVAL_DBL) >> (LD_DATA_SHIFT + 2)));
    pCld[i] = cld << (LD_DATA_SHIFT + 2);
  }
}

FDK_SACENC_ERROR fdk_sacenc_applyTtoBox(
    HANDLE_TTO_BOX hTtoBox, const INT nTimeSlots, const INT startTimeSlot,
    const INT nHybridBands, const FIXP_DPK *const *const ppHybridData1__FDK,
    const FIXP_DPK *const *const ppHybridData2__FDK, SCHAR *const pIccIdx,
    UCHAR *const pbIccQuantCoarse, SCHAR *const pCldIdx,
    UCHAR *const pbCldQuantCoarse, const INT bUseBBCues, INT *scaleCh1,
    INT *scaleCh2) {
  FDK_SACENC_ERROR error = SACENC_OK;

  C_ALLOC_SCRATCH_START(powerHybridData1__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_START(powerHybridData2__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_START(prodHybridDataReal__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_START(prodHybridDataImag__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)

  C_ALLOC_SCRATCH_START(IccDownmix__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_START(IccDownmixQuant__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_START(pbScaleCh1, INT, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_START(pbScaleCh2, INT, MAX_NUM_PARAM_BANDS)

  if ((hTtoBox == NULL) || (pCldIdx == NULL) || (pbCldQuantCoarse == NULL) ||
      (ppHybridData1__FDK == NULL) || (ppHybridData2__FDK == NULL) ||
      (pIccIdx == NULL) || (pbIccQuantCoarse == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    int j, pb;
    const int nParamBands = hTtoBox->nParameterBands;
    const int bUseEbQ = (hTtoBox->boxQuantMode == BOX_QUANTMODE_EBQ1) ||
                        (hTtoBox->boxQuantMode == BOX_QUANTMODE_EBQ2);

    if ((nHybridBands < 0) || (nHybridBands > hTtoBox->nHybridBandsMax)) {
      error = SACENC_INVALID_CONFIG;
      goto bail;
    }

    int outScale;
    int inScale = 5;

    FDKcalcPbScaleFactor(ppHybridData1__FDK,
                         hTtoBox->pParameterBand2HybridBandOffset, pbScaleCh1,
                         startTimeSlot, nTimeSlots, nParamBands);
    FDKcalcPbScaleFactor(ppHybridData2__FDK,
                         hTtoBox->pParameterBand2HybridBandOffset, pbScaleCh2,
                         startTimeSlot, nTimeSlots, nParamBands);

    for (j = 0, pb = 0; pb < nParamBands; pb++) {
      FIXP_DBL data1, data2;
      data1 = data2 = (FIXP_DBL)0;
      for (; j < hTtoBox->pParameterBand2HybridBandOffset[pb]; j++) {
        data1 += sumUpCplxPow2Dim2(ppHybridData1__FDK, SUM_UP_STATIC_SCALE,
                                   inScale + pbScaleCh1[pb], &outScale,
                                   startTimeSlot, nTimeSlots, j, j + 1);
        data2 += sumUpCplxPow2Dim2(ppHybridData2__FDK, SUM_UP_STATIC_SCALE,
                                   inScale + pbScaleCh2[pb], &outScale,
                                   startTimeSlot, nTimeSlots, j, j + 1);
      }
      powerHybridData1__FDK[pb] = data1;
      powerHybridData2__FDK[pb] = data2;
    }

    {
      for (j = 0, pb = 0; pb < nParamBands; pb++) {
        FIXP_DBL dataReal, dataImag;
        dataReal = dataImag = (FIXP_DBL)0;
        for (; j < hTtoBox->pParameterBand2HybridBandOffset[pb]; j++) {
          FIXP_DPK scalarProd;
          cplx_cplxScalarProduct(&scalarProd, ppHybridData1__FDK,
                                 ppHybridData2__FDK, inScale + pbScaleCh1[pb],
                                 inScale + pbScaleCh2[pb], &outScale,
                                 startTimeSlot, nTimeSlots, j, j + 1);
          dataReal += scalarProd.v.re;
          if (hTtoBox->pSubbandImagSign[j] < 0) {
            dataImag -= scalarProd.v.im;
          } else {
            dataImag += scalarProd.v.im;
          }
        }
        prodHybridDataReal__FDK[pb] = dataReal;
        prodHybridDataImag__FDK[pb] = dataImag;
      }

      if (SACENC_OK != (error = calculateIccFDK(
                            nParamBands, hTtoBox->iccCorrelationCoherenceBorder,
                            powerHybridData1__FDK, powerHybridData2__FDK,
                            prodHybridDataReal__FDK, prodHybridDataImag__FDK,
                            hTtoBox->pIcc__FDK))) {
        goto bail;
      }

      if (SACENC_OK != (error = calculateIccFDK(
                            nParamBands, nParamBands, powerHybridData1__FDK,
                            powerHybridData2__FDK, prodHybridDataReal__FDK,
                            prodHybridDataImag__FDK, IccDownmix__FDK))) {
        goto bail;
      }
    }

    if (!bUseEbQ) {
      CalculateCldFDK(hTtoBox->pCld__FDK, powerHybridData1__FDK,
                      powerHybridData2__FDK, *scaleCh1 + inScale + 1,
                      pbScaleCh1, *scaleCh2 + inScale + 1, pbScaleCh2,
                      nParamBands);
    }

    if (bUseBBCues) {
      ApplyBBCuesFDK(&hTtoBox->pCld__FDK[0], nParamBands);

      { ApplyBBCuesFDK(&hTtoBox->pIcc__FDK[0], nParamBands); }
    }

    {
      QuantizeCoefFDK(hTtoBox->pIcc__FDK, nParamBands,
                      hTtoBox->pIccQuantTable__FDK, hTtoBox->nIccQuantOffset,
                      hTtoBox->nIccQuantSteps, pIccIdx);
      QuantizeCoefFDK(IccDownmix__FDK, nParamBands,
                      hTtoBox->pIccQuantTable__FDK, hTtoBox->nIccQuantOffset,
                      hTtoBox->nIccQuantSteps, hTtoBox->pIccDownmixIdx);
      deQuantizeCoefFDK(hTtoBox->pIccDownmixIdx, nParamBands,
                        hTtoBox->pIccQuantTable__FDK, hTtoBox->nIccQuantOffset,
                        IccDownmixQuant__FDK);

      *pbIccQuantCoarse = hTtoBox->bUseCoarseQuantIcc;
    }

    if (!bUseEbQ) {
      QuantizeCoefFDK(hTtoBox->pCld__FDK, nParamBands,
                      hTtoBox->pCldQuantTableEnc__FDK, hTtoBox->nCldQuantOffset,
                      hTtoBox->nCldQuantSteps, pCldIdx);
      deQuantizeCoefFDK(pCldIdx, nParamBands, hTtoBox->pCldQuantTableDec__FDK,
                        hTtoBox->nCldQuantOffset, hTtoBox->pCldQuant__FDK);
    } else {
      FDKmemcpy(pCldIdx, hTtoBox->pCldEbQIdx, nParamBands * sizeof(SCHAR));
    }
    *pbCldQuantCoarse = hTtoBox->bUseCoarseQuantCld;
  }

bail:
  C_ALLOC_SCRATCH_END(pbScaleCh2, INT, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_END(pbScaleCh1, INT, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_END(IccDownmixQuant__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_END(IccDownmix__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)

  C_ALLOC_SCRATCH_END(prodHybridDataImag__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_END(prodHybridDataReal__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_END(powerHybridData2__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)
  C_ALLOC_SCRATCH_END(powerHybridData1__FDK, FIXP_DBL, MAX_NUM_PARAM_BANDS)

  return error;
}

INT fdk_sacenc_subband2ParamBand(const BOX_SUBBAND_CONFIG boxSubbandConfig,
                                 const INT nSubband) {
  INT nParamBand = -1;
  const UCHAR *pSubband2ParameterIndex =
      getSubband2ParameterIndex(boxSubbandConfig);

  if (pSubband2ParameterIndex != NULL) {
    const int hybrid_resolution = 64;

    if ((nSubband > -1) && (nSubband < hybrid_resolution)) {
      nParamBand = pSubband2ParameterIndex[nSubband];
    }
  }

  return nParamBand;
}
