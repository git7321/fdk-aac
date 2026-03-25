/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: independent channel concealment

*******************************************************************************/

#include "conceal.h"

#include "aac_rom.h"
#include "genericStds.h"
#include "aacdec_pns.h"
#include "block.h"

#define CONCEAL_DFLT_COMF_NOISE_LEVEL (0x100000)

#define CONCEAL_NOT_DEFINED ((UCHAR)-1)

#define CONCEAL_DFLT_FADEOUT_FRAMES (6)
#define CONCEAL_DFLT_FADEIN_FRAMES (5)
#define CONCEAL_DFLT_MUTE_RELEASE_FRAMES (0)

#define CONCEAL_DFLT_FADE_FACTOR (0.707106781186548f)

#define FIXP_ZERO FL2FXCONST_DBL(0.0f)
#define FIXP_ONE FL2FXCONST_DBL(1.0f)
#define FIXP_FL_CORRECTION FL2FXCONST_DBL(0.53333333333333333f)

#define CONCEAL_PARAMETER_BITS (8)
#define CONCEAL_MAX_QUANT_FACTOR ((1 << CONCEAL_PARAMETER_BITS) - 1)
#define CONCEAL_MIN_ATTENUATION_FACTOR_025_LD \
  FL2FXCONST_DBL(-0.041524101186092029596853445212299)
#define CONCEAL_MIN_ATTENUATION_FACTOR_050_LD \
  FL2FXCONST_DBL(-0.083048202372184059253597008145293)

typedef enum {
  CConcealment_NoExpand,
  CConcealment_Expand,
  CConcealment_Compress
} CConcealmentExpandType;

static const FIXP_SGL facMod4Table[4] = {
    FL2FXCONST_SGL(0.500000000f),
    FL2FXCONST_SGL(0.594603558f),
    FL2FXCONST_SGL(0.707106781f),
    FL2FXCONST_SGL(0.840896415f)
};

static void CConcealment_CalcBandEnergy(
    FIXP_DBL *spectrum, const SamplingRateInfo *pSamplingRateInfo,
    const int blockType, CConcealmentExpandType ex, int *sfbEnergy);

static void CConcealment_InterpolateBuffer(FIXP_DBL *spectrum,
                                           SHORT *pSpecScalePrev,
                                           SHORT *pSpecScaleAct,
                                           SHORT *pSpecScaleOut, int *enPrv,
                                           int *enAct, int sfbCnt,
                                           const SHORT *pSfbOffset);

static int CConcealment_ApplyInter(
    CConcealmentInfo *pConcealmentInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const int samplesPerFrame,
    const int improveTonal, const int frameOk, const int mute_release_active);

static int CConcealment_ApplyNoise(
    CConcealmentInfo *pConcealmentInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const int samplesPerFrame,
    const UINT flags);

static void CConcealment_UpdateState(
    CConcealmentInfo *pConcealmentInfo, int frameOk,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const int samplesPerFrame, CAacDecoderChannelInfo *pAacDecoderChannelInfo);

static void CConcealment_ApplyRandomSign(int iRandomPhase, FIXP_DBL *spec,
                                         int samplesPerFrame);

static void CConcealment_TDFadePcmAtt(int start, int len, FIXP_DBL fadeStart,
                                      FIXP_DBL fadeStop, PCM_DEC *pcmdata);
static void CConcealment_TDFadeFillFadingStations(FIXP_DBL *fadingStations,
                                                  int *fadingSteps,
                                                  FIXP_DBL fadeStop,
                                                  FIXP_DBL fadeStart,
                                                  TDfadingType fadingType);
static void CConcealment_TDFading_doLinearFadingSteps(int *fadingSteps);

static int CConcealment_ApplyFadeOut(
    int mode, CConcealmentInfo *pConcealmentInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const int samplesPerFrame, CAacDecoderChannelInfo *pAacDecoderChannelInfo);

static int CConcealment_TDNoise_Random(ULONG *seed);
static void CConcealment_TDNoise_Apply(CConcealmentInfo *const pConcealmentInfo,
                                       const int len,
                                       const INT aacOutDataHeadroom,
                                       PCM_DEC *const pcmdata);

static BLOCK_TYPE CConcealment_GetWinSeq(int prevWinSeq) {
  BLOCK_TYPE newWinSeq = BLOCK_LONG;

  if (prevWinSeq == BLOCK_START || prevWinSeq == BLOCK_SHORT) {
    newWinSeq = BLOCK_STOP;
  }

  return (newWinSeq);
}

void CConcealment_InitCommonData(CConcealParams *pConcealCommonData) {
  if (pConcealCommonData != NULL) {
    int i;

    pConcealCommonData->method = ConcealMethodInter;

    pConcealCommonData->numFadeOutFrames = CONCEAL_DFLT_FADEOUT_FRAMES;
    pConcealCommonData->numFadeInFrames = CONCEAL_DFLT_FADEIN_FRAMES;
    pConcealCommonData->numMuteReleaseFrames = CONCEAL_DFLT_MUTE_RELEASE_FRAMES;

    pConcealCommonData->comfortNoiseLevel =
        (FIXP_DBL)CONCEAL_DFLT_COMF_NOISE_LEVEL;

    pConcealCommonData->fadeOutFactor[0] =
        FL2FXCONST_SGL(CONCEAL_DFLT_FADE_FACTOR);
    pConcealCommonData->fadeInFactor[0] = pConcealCommonData->fadeOutFactor[0];

    for (i = 1; i < CONCEAL_MAX_NUM_FADE_FACTORS; i++) {
      pConcealCommonData->fadeOutFactor[i] =
          FX_DBL2FX_SGL(fMult(pConcealCommonData->fadeOutFactor[i - 1],
                              FL2FXCONST_SGL(CONCEAL_DFLT_FADE_FACTOR)));
      pConcealCommonData->fadeInFactor[i] =
          pConcealCommonData->fadeOutFactor[i];
    }
  }
}

CConcealmentMethod CConcealment_GetMethod(CConcealParams *pConcealCommonData) {
  CConcealmentMethod method = ConcealMethodNone;

  if (pConcealCommonData != NULL) {
    method = pConcealCommonData->method;
  }

  return (method);
}

void CConcealment_InitChannelData(CConcealmentInfo *pConcealChannelInfo,
                                  CConcealParams *pConcealCommonData,
                                  AACDEC_RENDER_MODE initRenderMode,
                                  int samplesPerFrame) {
  int i;
  pConcealChannelInfo->TDNoiseSeed = 0;
  FDKmemclear(pConcealChannelInfo->TDNoiseStates,
              sizeof(pConcealChannelInfo->TDNoiseStates));
  pConcealChannelInfo->TDNoiseCoef[0] = FL2FXCONST_SGL(0.05f);
  pConcealChannelInfo->TDNoiseCoef[1] = FL2FXCONST_SGL(0.5f);
  pConcealChannelInfo->TDNoiseCoef[2] = FL2FXCONST_SGL(0.45f);

  pConcealChannelInfo->pConcealParams = pConcealCommonData;

  pConcealChannelInfo->lastRenderMode = initRenderMode;

  pConcealChannelInfo->windowShape = CONCEAL_NOT_DEFINED;
  pConcealChannelInfo->windowSequence = BLOCK_LONG;
  pConcealChannelInfo->lastWinGrpLen = 1;

  pConcealChannelInfo->concealState = ConcealState_Ok;

  FDKmemclear(pConcealChannelInfo->spectralCoefficient,
              1024 * sizeof(FIXP_CNCL));

  for (i = 0; i < 8; i++) {
    pConcealChannelInfo->specScale[i] = 0;
  }

  pConcealChannelInfo->iRandomPhase = 0;

  pConcealChannelInfo->prevFrameOk[0] = 1;
  pConcealChannelInfo->prevFrameOk[1] = 1;

  pConcealChannelInfo->cntFadeFrames = 0;
  pConcealChannelInfo->cntValidFrames = 0;
  pConcealChannelInfo->fade_old = (FIXP_DBL)MAXVAL_DBL;
  pConcealChannelInfo->winGrpOffset[0] = 0;
  pConcealChannelInfo->winGrpOffset[1] = 0;
  pConcealChannelInfo->attGrpOffset[0] = 0;
  pConcealChannelInfo->attGrpOffset[1] = 0;
}

AAC_DECODER_ERROR
CConcealment_SetParams(CConcealParams *concealParams, int method,
                       int fadeOutSlope, int fadeInSlope, int muteRelease,
                       FIXP_DBL comfNoiseLevel) {
  if (method != AACDEC_CONCEAL_PARAM_NOT_SPECIFIED) {
    switch ((CConcealmentMethod)method) {
      case ConcealMethodMute:
      case ConcealMethodNoise:
      case ConcealMethodInter:
        if (concealParams == NULL) {
          return AAC_DEC_INVALID_HANDLE;
        } else {
          concealParams->method = (CConcealmentMethod)method;
        }
        break;

      default:
        return AAC_DEC_SET_PARAM_FAIL;
    }
  }

  if (fadeOutSlope != AACDEC_CONCEAL_PARAM_NOT_SPECIFIED) {
    if ((fadeOutSlope < CONCEAL_MAX_NUM_FADE_FACTORS) && (fadeOutSlope >= 0)) {
      if (concealParams == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      } else {
        concealParams->numFadeOutFrames = fadeOutSlope;
      }
    } else {
      return AAC_DEC_SET_PARAM_FAIL;
    }
  }

  if (fadeInSlope != AACDEC_CONCEAL_PARAM_NOT_SPECIFIED) {
    if ((fadeInSlope < CONCEAL_MAX_NUM_FADE_FACTORS) && (fadeInSlope >= 0)) {
      if (concealParams == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      } else {
        concealParams->numFadeInFrames = fadeInSlope;
      }
    } else {
      return AAC_DEC_SET_PARAM_FAIL;
    }
  }

  if (muteRelease != AACDEC_CONCEAL_PARAM_NOT_SPECIFIED) {
    if ((muteRelease < (CONCEAL_MAX_NUM_FADE_FACTORS << 1)) &&
        (muteRelease >= 0)) {
      if (concealParams == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      } else {
        concealParams->numMuteReleaseFrames = muteRelease;
      }
    } else {
      return AAC_DEC_SET_PARAM_FAIL;
    }
  }

  if (comfNoiseLevel != (FIXP_DBL)AACDEC_CONCEAL_PARAM_NOT_SPECIFIED) {
    if ((comfNoiseLevel < (FIXP_DBL)0) ||
        (comfNoiseLevel > (FIXP_DBL)MAXVAL_DBL)) {
      return AAC_DEC_SET_PARAM_FAIL;
    }
    if (concealParams == NULL) {
      return AAC_DEC_INVALID_HANDLE;
    } else {
      concealParams->comfortNoiseLevel = (FIXP_DBL)comfNoiseLevel;
    }
  }

  return (AAC_DEC_OK);
}

AAC_DECODER_ERROR
CConcealment_SetAttenuation(CConcealParams *concealParams,
                            const SHORT *fadeOutAttenuationVector,
                            const SHORT *fadeInAttenuationVector) {
  if ((fadeOutAttenuationVector == NULL) && (fadeInAttenuationVector == NULL)) {
    return AAC_DEC_SET_PARAM_FAIL;
  }

  if (fadeOutAttenuationVector != NULL) {
    int i;

    for (i = 0; i < CONCEAL_MAX_NUM_FADE_FACTORS; i++) {
      if ((fadeOutAttenuationVector[i] < 0) ||
          (fadeOutAttenuationVector[i] > CONCEAL_MAX_QUANT_FACTOR)) {
        return AAC_DEC_SET_PARAM_FAIL;
      }
    }
    if (concealParams == NULL) {
      return AAC_DEC_INVALID_HANDLE;
    }

    for (i = 0; i < CONCEAL_MAX_NUM_FADE_FACTORS; i++) {
      concealParams->fadeOutFactor[i] =
          FX_DBL2FX_SGL(fLdPow(CONCEAL_MIN_ATTENUATION_FACTOR_025_LD, 0,
                               (FIXP_DBL)((INT)(FL2FXCONST_DBL(1.0 / 2.0) >>
                                                (CONCEAL_PARAMETER_BITS - 1)) *
                                          (INT)fadeOutAttenuationVector[i]),
                               CONCEAL_PARAMETER_BITS));
    }
  }

  if (fadeInAttenuationVector != NULL) {
    int i;

    for (i = 0; i < CONCEAL_MAX_NUM_FADE_FACTORS; i++) {
      if ((fadeInAttenuationVector[i] < 0) ||
          (fadeInAttenuationVector[i] > CONCEAL_MAX_QUANT_FACTOR)) {
        return AAC_DEC_SET_PARAM_FAIL;
      }
    }
    if (concealParams == NULL) {
      return AAC_DEC_INVALID_HANDLE;
    }

    for (i = 0; i < CONCEAL_MAX_NUM_FADE_FACTORS; i++) {
      concealParams->fadeInFactor[i] = FX_DBL2FX_SGL(
          fLdPow(CONCEAL_MIN_ATTENUATION_FACTOR_025_LD, 0,
                 (FIXP_DBL)((INT)(FIXP_ONE >> CONCEAL_PARAMETER_BITS) *
                            (INT)fadeInAttenuationVector[i]),
                 CONCEAL_PARAMETER_BITS));
    }
  }

  return (AAC_DEC_OK);
}

CConcealmentState CConcealment_GetState(CConcealmentInfo *pConcealChannelInfo) {
  CConcealmentState state = ConcealState_Ok;

  if (pConcealChannelInfo != NULL) {
    state = pConcealChannelInfo->concealState;
  }

  return (state);
}

void CConcealment_Store(
    CConcealmentInfo *hConcealmentInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo) {
  UCHAR nbDiv = NB_DIV;

  if (!(pAacDecoderChannelInfo->renderMode == AACDEC_RENDER_LPD &&
        pAacDecoderChannelInfo->data.usac.mod[nbDiv - 1] == 0))

  {
    FIXP_DBL *pSpectralCoefficient =
        SPEC_LONG(pAacDecoderChannelInfo->pSpectralCoefficient);
    SHORT *pSpecScale = pAacDecoderChannelInfo->specScale;
    CIcsInfo *pIcsInfo = &pAacDecoderChannelInfo->icsInfo;

    SHORT tSpecScale[8];
    UCHAR tWindowShape;
    BLOCK_TYPE tWindowSequence;

    tWindowSequence = hConcealmentInfo->windowSequence;
    tWindowShape = hConcealmentInfo->windowShape;

    FDKmemcpy(tSpecScale, hConcealmentInfo->specScale, 8 * sizeof(SHORT));

    hConcealmentInfo->windowSequence = GetWindowSequence(pIcsInfo);
    hConcealmentInfo->windowShape = GetWindowShape(pIcsInfo);
    hConcealmentInfo->lastWinGrpLen =
        *(GetWindowGroupLengthTable(pIcsInfo) + GetWindowGroups(pIcsInfo) - 1);

    FDKmemcpy(hConcealmentInfo->specScale, pSpecScale, 8 * sizeof(SHORT));

    if (hConcealmentInfo->pConcealParams->method < ConcealMethodInter) {
#if (CNCL_FRACT_BITS == DFRACT_BITS)
      FDKmemcpy(hConcealmentInfo->spectralCoefficient, pSpectralCoefficient,
                1024 * sizeof(FIXP_CNCL));
#else
      FIXP_CNCL *RESTRICT pCncl =
          &hConcealmentInfo->spectralCoefficient[1024 - 1];
      FIXP_DBL *RESTRICT pSpec = &pSpectralCoefficient[1024 - 1];
      int i;
      for (i = 1024; i != 0; i--) {
        *pCncl-- = FX_DBL2FX_CNCL(*pSpec--);
      }
#endif
    } else {
#if (FIXP_CNCL == FIXP_DBL)
      C_ALLOC_SCRATCH_START(pSpecTmp, FIXP_DBL, 1024);
      FDKmemcpy(pSpecTmp, pSpectralCoefficient, 1024 * sizeof(FIXP_DBL));
      FDKmemcpy(pSpectralCoefficient, hConcealmentInfo->spectralCoefficient,
                1024 * sizeof(FIXP_DBL));
      FDKmemcpy(hConcealmentInfo->spectralCoefficient, pSpecTmp,
                1024 * sizeof(FIXP_DBL));
      C_ALLOC_SCRATCH_END(pSpecTmp, FIXP_DBL, 1024);
#else
      FIXP_CNCL *RESTRICT pCncl =
          &hConcealmentInfo->spectralCoefficient[1024 - 1];
      FIXP_DBL *RESTRICT pSpec = &pSpectralCoefficient[1024 - 1];
      FIXP_DBL tSpec;

      for (int i = 1024; i != 0; i--) {
        tSpec = *pSpec;
        *pSpec-- = FX_CNCL2FX_DBL(*pCncl);
        *pCncl-- = FX_DBL2FX_CNCL(tSpec);
      }
#endif

      pIcsInfo->WindowSequence = tWindowSequence;
      pIcsInfo->WindowShape = tWindowShape;

      FDKmemcpy(pSpecScale, tSpecScale, 8 * sizeof(SHORT));
    }
  }

  if (pAacDecoderChannelInfo->renderMode == AACDEC_RENDER_LPD) {
    FDKmemcpy(hConcealmentInfo->lsf4, pAacDecoderStaticChannelInfo->lpc4_lsf,
              sizeof(hConcealmentInfo->lsf4));
    hConcealmentInfo->last_tcx_gain =
        pAacDecoderStaticChannelInfo->last_tcx_gain;
    hConcealmentInfo->last_tcx_gain_e =
        pAacDecoderStaticChannelInfo->last_tcx_gain_e;
  }
}

int CConcealment_Apply(
    CConcealmentInfo *hConcealmentInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const int samplesPerFrame,
    const UCHAR lastLpdMode, const int frameOk, const UINT flags) {
  int appliedProcessing = 0;
  const int mute_release_active =
      frameOk && (hConcealmentInfo->concealState >= ConcealState_Mute) &&
      (hConcealmentInfo->cntValidFrames + 1 <=
       hConcealmentInfo->pConcealParams->numMuteReleaseFrames);

  if (hConcealmentInfo->windowShape == CONCEAL_NOT_DEFINED) {
    hConcealmentInfo->windowShape = pAacDecoderChannelInfo->icsInfo.WindowShape;
  }

  if (frameOk && !mute_release_active) {
    hConcealmentInfo->lastRenderMode =
        (SCHAR)pAacDecoderChannelInfo->renderMode;

    CConcealment_Store(hConcealmentInfo, pAacDecoderChannelInfo,
                       pAacDecoderStaticChannelInfo);
    hConcealmentInfo->iRandomPhase = 0;
  } else {
    if (hConcealmentInfo->lastRenderMode == AACDEC_RENDER_INVALID) {
      hConcealmentInfo->lastRenderMode = AACDEC_RENDER_IMDCT;
    }
    pAacDecoderChannelInfo->renderMode =
        (AACDEC_RENDER_MODE)hConcealmentInfo->lastRenderMode;
  }

  CConcealment_UpdateState(hConcealmentInfo, frameOk,
                           pAacDecoderStaticChannelInfo, samplesPerFrame,
                           pAacDecoderChannelInfo);

  {
    if (!frameOk && pAacDecoderChannelInfo->renderMode == AACDEC_RENDER_IMDCT) {
      CLpc_Conceal(pAacDecoderChannelInfo->data.usac.lsp_coeff,
                   pAacDecoderStaticChannelInfo->lpc4_lsf,
                   pAacDecoderStaticChannelInfo->lsf_adaptive_mean,
                   hConcealmentInfo->lastRenderMode == AACDEC_RENDER_IMDCT);
      FDKmemcpy(hConcealmentInfo->lsf4, pAacDecoderStaticChannelInfo->lpc4_lsf,
                sizeof(pAacDecoderStaticChannelInfo->lpc4_lsf));
    }

    if ((!frameOk || mute_release_active) &&
        (pAacDecoderChannelInfo->renderMode == AACDEC_RENDER_LPD)) {
      FDKmemcpy(pAacDecoderStaticChannelInfo->lpc4_lsf, hConcealmentInfo->lsf4,
                sizeof(pAacDecoderStaticChannelInfo->lpc4_lsf));
      pAacDecoderStaticChannelInfo->last_tcx_gain =
          hConcealmentInfo->last_tcx_gain;
      pAacDecoderStaticChannelInfo->last_tcx_gain_e =
          hConcealmentInfo->last_tcx_gain_e;
    }

    if (!(pAacDecoderChannelInfo->renderMode == AACDEC_RENDER_LPD &&
          pAacDecoderStaticChannelInfo->last_lpd_mode == 0)) {
      switch (hConcealmentInfo->pConcealParams->method) {
        default:
        case ConcealMethodMute:
          if (!frameOk) {
            FDKmemclear(pAacDecoderChannelInfo->pSpectralCoefficient,
                        samplesPerFrame * sizeof(FIXP_DBL));
            pAacDecoderChannelInfo->icsInfo.WindowShape =
                hConcealmentInfo->windowShape;
            appliedProcessing = 1;
          }
          break;

        case ConcealMethodNoise:
          appliedProcessing = CConcealment_ApplyNoise(
              hConcealmentInfo, pAacDecoderChannelInfo,
              pAacDecoderStaticChannelInfo, pSamplingRateInfo, samplesPerFrame,
              flags);
          break;

        case ConcealMethodInter:
          appliedProcessing = CConcealment_ApplyInter(
              hConcealmentInfo, pAacDecoderChannelInfo, pSamplingRateInfo,
              samplesPerFrame, 0,
              frameOk, mute_release_active);
          break;
      }
    } else if (!frameOk || mute_release_active) {
      FIXP_DBL *pSpectralCoefficient =
          SPEC_LONG(pAacDecoderChannelInfo->pSpectralCoefficient);
      SHORT *pSpecScale = pAacDecoderChannelInfo->specScale;
      CIcsInfo *pIcsInfo = &pAacDecoderChannelInfo->icsInfo;
#if (CNCL_FRACT_BITS != DFRACT_BITS)
      FIXP_CNCL *RESTRICT pCncl =
          &hConcealmentInfo->spectralCoefficient[1024 - 1];
      FIXP_DBL *RESTRICT pSpec = &pSpectralCoefficient[1024 - 1];
      int i;
#endif

      pIcsInfo->WindowSequence = hConcealmentInfo->windowSequence;
      pIcsInfo->WindowShape = hConcealmentInfo->windowShape;

      if (hConcealmentInfo->concealState != ConcealState_Mute) {
        FDKmemcpy(pSpecScale, hConcealmentInfo->specScale, 8 * sizeof(SHORT));

#if (CNCL_FRACT_BITS == DFRACT_BITS)
        FDKmemcpy(pSpectralCoefficient, hConcealmentInfo->spectralCoefficient,
                  1024 * sizeof(FIXP_DBL));
#else
        for (i = 1024; i != 0; i--) {
          *pSpec-- = FX_CNCL2FX_DBL(*pCncl--);
        }
#endif
      } else {
        FDKmemclear(pSpecScale, 8 * sizeof(SHORT));

        FDKmemclear(pSpectralCoefficient, 1024 * sizeof(FIXP_CNCL));
      }
    }
  }
  hConcealmentInfo->prevFrameOk[0] = hConcealmentInfo->prevFrameOk[1];
  hConcealmentInfo->prevFrameOk[1] = frameOk;

  return mute_release_active ? -1 : appliedProcessing;
}

static int CConcealment_ApplyNoise(
    CConcealmentInfo *pConcealmentInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const int samplesPerFrame,
    const UINT flags) {
  FIXP_DBL *pSpectralCoefficient =
      SPEC_LONG(pAacDecoderChannelInfo->pSpectralCoefficient);
  CIcsInfo *pIcsInfo = &pAacDecoderChannelInfo->icsInfo;

  int appliedProcessing = 0;

  FDK_ASSERT(pConcealmentInfo != NULL);
  FDK_ASSERT((samplesPerFrame >= 120) && (samplesPerFrame <= 1024));

  switch (pConcealmentInfo->concealState) {
    case ConcealState_Ok:
      break;

    case ConcealState_Single:
    case ConcealState_FadeOut:
      appliedProcessing = CConcealment_ApplyFadeOut(
          1, pConcealmentInfo, pAacDecoderStaticChannelInfo,
          samplesPerFrame, pAacDecoderChannelInfo);
      break;

    case ConcealState_Mute: {
      pIcsInfo->Valid = 0;
      pIcsInfo->WindowShape =
          pConcealmentInfo->windowShape;
      pIcsInfo->WindowSequence =
          CConcealment_GetWinSeq(pConcealmentInfo->windowSequence);
      pConcealmentInfo->windowSequence =
          pIcsInfo->WindowSequence;

      FDKmemclear(pSpectralCoefficient, samplesPerFrame * sizeof(FIXP_DBL));
      FDKmemclear(pConcealmentInfo->spectralCoefficient,
                  samplesPerFrame * sizeof(FIXP_DBL));

      appliedProcessing = 1;
    } break;

    case ConcealState_FadeIn: {

      appliedProcessing = 1;
    } break;

    default:
      FDK_ASSERT(0);
      break;
  }

  return appliedProcessing;
}

static int CConcealment_ApplyInter(
    CConcealmentInfo *pConcealmentInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const int samplesPerFrame,
    const int improveTonal, const int frameOk, const int mute_release_active) {
#if defined(FDK_ASSERT_ENABLE)
  CConcealParams *pConcealCommonData = pConcealmentInfo->pConcealParams;
#endif

  FIXP_DBL *pSpectralCoefficient =
      SPEC_LONG(pAacDecoderChannelInfo->pSpectralCoefficient);
  CIcsInfo *pIcsInfo = &pAacDecoderChannelInfo->icsInfo;
  SHORT *pSpecScale = pAacDecoderChannelInfo->specScale;

  int sfbEnergyPrev[64];
  int sfbEnergyAct[64];

  int i, appliedProcessing = 0;

  FDKmemclear(sfbEnergyPrev, 64 * sizeof(int));
  FDKmemclear(sfbEnergyAct, 64 * sizeof(int));

  if (!frameOk || mute_release_active) {
    pIcsInfo->WindowShape = pConcealmentInfo->windowShape;
    pIcsInfo->WindowSequence = pConcealmentInfo->windowSequence;

    for (i = 0; i < samplesPerFrame; i++) {
      pSpectralCoefficient[i] =
          FX_CNCL2FX_DBL(pConcealmentInfo->spectralCoefficient[i]);
    }

    FDKmemcpy(pSpecScale, pConcealmentInfo->specScale, 8 * sizeof(SHORT));
  }

  if (!pConcealmentInfo->prevFrameOk[1] || mute_release_active) {
    if (frameOk && pConcealmentInfo->prevFrameOk[0] && !mute_release_active) {
      appliedProcessing = 1;

      if (pIcsInfo->WindowSequence == BLOCK_SHORT) {
        int wnd;

        if (pConcealmentInfo->windowSequence ==
            BLOCK_SHORT) {

          int scaleFactorBandsTotal =
              pSamplingRateInfo->NumberOfScaleFactorBands_Short;
          const SHORT *pSfbOffset = pSamplingRateInfo->ScaleFactorBands_Short;
          pIcsInfo->WindowShape = (samplesPerFrame <= 512) ? 2 : 1;
          pIcsInfo->WindowSequence = BLOCK_SHORT;

          for (wnd = 0; wnd < 8; wnd++) {
            CConcealment_CalcBandEnergy(
                &pSpectralCoefficient[wnd *
                                      (samplesPerFrame / 8)],
                pSamplingRateInfo, BLOCK_SHORT, CConcealment_NoExpand,
                sfbEnergyPrev);

            CConcealment_CalcBandEnergy(
                &pConcealmentInfo->spectralCoefficient[wnd * (samplesPerFrame /
                                                              8)],
                pSamplingRateInfo, BLOCK_SHORT, CConcealment_NoExpand,
                sfbEnergyAct);

            CConcealment_InterpolateBuffer(
                &pSpectralCoefficient[wnd *
                                      (samplesPerFrame / 8)],
                &pSpecScale[wnd], &pConcealmentInfo->specScale[wnd],
                &pSpecScale[wnd], sfbEnergyPrev, sfbEnergyAct,
                scaleFactorBandsTotal, pSfbOffset);
          }
        } else {

          int scaleFactorBandsTotal =
              pSamplingRateInfo->NumberOfScaleFactorBands_Long;
          const SHORT *pSfbOffset = pSamplingRateInfo->ScaleFactorBands_Long;
          SHORT specScaleOut;

          CConcealment_CalcBandEnergy(
              &pSpectralCoefficient[samplesPerFrame -
                                    (samplesPerFrame /
                                     8)],
              pSamplingRateInfo, BLOCK_SHORT, CConcealment_Expand,
              sfbEnergyAct);

          CConcealment_CalcBandEnergy(
              pConcealmentInfo->spectralCoefficient,
              pSamplingRateInfo, BLOCK_LONG, CConcealment_NoExpand,
              sfbEnergyPrev);

          pIcsInfo->WindowShape = 0;
          pIcsInfo->WindowSequence = BLOCK_STOP;

          for (i = 0; i < samplesPerFrame; i++) {
            pSpectralCoefficient[i] =
                pConcealmentInfo->spectralCoefficient[i];
          }

          for (i = 0; i < 8; i++) {
            if (pSpecScale[i] > pSpecScale[0]) {
              pSpecScale[0] = pSpecScale[i];
            }
          }

          CConcealment_InterpolateBuffer(
              pSpectralCoefficient,
              &pConcealmentInfo->specScale[0], &pSpecScale[0], &specScaleOut,
              sfbEnergyPrev, sfbEnergyAct, scaleFactorBandsTotal, pSfbOffset);

          pSpecScale[0] = specScaleOut;
        }
      } else {

        int scaleFactorBandsTotal =
            pSamplingRateInfo->NumberOfScaleFactorBands_Long;
        const SHORT *pSfbOffset = pSamplingRateInfo->ScaleFactorBands_Long;
        SHORT specScaleAct = pConcealmentInfo->specScale[0];

        CConcealment_CalcBandEnergy(pSpectralCoefficient,
                                    pSamplingRateInfo, BLOCK_LONG,
                                    CConcealment_NoExpand, sfbEnergyPrev);

        if (pConcealmentInfo->windowSequence ==
            BLOCK_SHORT) {

          pIcsInfo->WindowShape = (samplesPerFrame <= 512) ? 2 : 1;
          pIcsInfo->WindowSequence = BLOCK_START;

          for (i = 1; i < 8; i++) {
            if (pConcealmentInfo->specScale[i] > specScaleAct) {
              specScaleAct = pConcealmentInfo->specScale[i];
            }
          }

          CConcealment_CalcBandEnergy(
              pConcealmentInfo->spectralCoefficient,
              pSamplingRateInfo, BLOCK_SHORT, CConcealment_Expand,
              sfbEnergyAct);
        } else {

          pIcsInfo->WindowShape = 0;
          pIcsInfo->WindowSequence = BLOCK_LONG;

          CConcealment_CalcBandEnergy(
              pConcealmentInfo->spectralCoefficient,
              pSamplingRateInfo, BLOCK_LONG, CConcealment_NoExpand,
              sfbEnergyAct);
        }

        CConcealment_InterpolateBuffer(
            pSpectralCoefficient,
            &pSpecScale[0], &specScaleAct, &pSpecScale[0], sfbEnergyPrev,
            sfbEnergyAct, scaleFactorBandsTotal, pSfbOffset);
      }
    }

    CConcealment_ApplyRandomSign(pConcealmentInfo->iRandomPhase,
                                 pSpectralCoefficient, samplesPerFrame);
    pConcealmentInfo->iRandomPhase =
        (pConcealmentInfo->iRandomPhase + 1) & (AAC_NF_NO_RANDOM_VAL - 1);
  }

  switch (pConcealmentInfo->concealState) {
    case ConcealState_Single:
      appliedProcessing = 1;
      break;

    case ConcealState_FadeOut: {
      FDK_ASSERT(pConcealmentInfo->cntFadeFrames >= 0);
      FDK_ASSERT(pConcealmentInfo->cntFadeFrames <
                 CONCEAL_MAX_NUM_FADE_FACTORS);
      FDK_ASSERT(pConcealmentInfo->cntFadeFrames <
                 pConcealCommonData->numFadeOutFrames);

      appliedProcessing = 1;
    } break;

    case ConcealState_FadeIn: {
      FDK_ASSERT(pConcealmentInfo->cntFadeFrames >= 0);
      FDK_ASSERT(pConcealmentInfo->cntFadeFrames <
                 CONCEAL_MAX_NUM_FADE_FACTORS);
      FDK_ASSERT(pConcealmentInfo->cntFadeFrames <
                 pConcealCommonData->numFadeInFrames);

      appliedProcessing = 1;
    } break;

    case ConcealState_Mute: {
      pIcsInfo->Valid = 0;
      pIcsInfo->WindowShape =
          pConcealmentInfo->windowShape;
      pIcsInfo->WindowSequence =
          CConcealment_GetWinSeq(pConcealmentInfo->windowSequence);
      pConcealmentInfo->windowSequence =
          pIcsInfo->WindowSequence;

      FDKmemclear(pSpectralCoefficient, samplesPerFrame * sizeof(FIXP_DBL));

      appliedProcessing = 1;
    } break;

    default:
      break;
  }

  return appliedProcessing;
}

static void CConcealment_CalcBandEnergy(
    FIXP_DBL *spectrum, const SamplingRateInfo *pSamplingRateInfo,
    const int blockType, CConcealmentExpandType expandType, int *sfbEnergy) {
  const SHORT *pSfbOffset;
  int line, sfb, scaleFactorBandsTotal = 0;

  line = 0;

  switch (blockType) {
    case BLOCK_LONG:
    case BLOCK_START:
    case BLOCK_STOP:

      if (expandType == CConcealment_NoExpand) {
        scaleFactorBandsTotal =
            pSamplingRateInfo->NumberOfScaleFactorBands_Long;
        pSfbOffset = pSamplingRateInfo->ScaleFactorBands_Long;

        for (sfb = 0; sfb < scaleFactorBandsTotal; sfb++) {
          FIXP_DBL enAccu = (FIXP_DBL)(LONG)1;
          int sfbScale =
              (sizeof(LONG) << 3) -
              CntLeadingZeros(pSfbOffset[sfb + 1] - pSfbOffset[sfb]) - 1;
          for (; line < pSfbOffset[sfb + 1]; line++) {
            enAccu += fPow2Div2(*(spectrum + line)) >> sfbScale;
          }
          *(sfbEnergy + sfb) = CntLeadingZeros(enAccu) - 1;
        }
      } else {
        scaleFactorBandsTotal =
            pSamplingRateInfo->NumberOfScaleFactorBands_Short;
        pSfbOffset = pSamplingRateInfo->ScaleFactorBands_Short;

        for (sfb = 0; sfb < scaleFactorBandsTotal; sfb++) {
          FIXP_DBL enAccu = (FIXP_DBL)(LONG)1;
          int sfbScale =
              (sizeof(LONG) << 3) -
              CntLeadingZeros(pSfbOffset[sfb + 1] - pSfbOffset[sfb]) - 1;
          for (; line < pSfbOffset[sfb + 1] << 3; line++) {
            enAccu +=
                (enAccu + (fPow2Div2(*(spectrum + line)) >> sfbScale)) >> 3;
          }
          *(sfbEnergy + sfb) = CntLeadingZeros(enAccu) - 1;
        }
      }
      break;

    case BLOCK_SHORT:

      if (expandType == CConcealment_NoExpand) {
        scaleFactorBandsTotal =
            pSamplingRateInfo->NumberOfScaleFactorBands_Short;
        pSfbOffset = pSamplingRateInfo->ScaleFactorBands_Short;

        for (sfb = 0; sfb < scaleFactorBandsTotal; sfb++) {
          FIXP_DBL enAccu = (FIXP_DBL)(LONG)1;
          int sfbScale =
              (sizeof(LONG) << 3) -
              CntLeadingZeros(pSfbOffset[sfb + 1] - pSfbOffset[sfb]) - 1;
          for (; line < pSfbOffset[sfb + 1]; line++) {
            enAccu += fPow2Div2(*(spectrum + line)) >> sfbScale;
          }
          *(sfbEnergy + sfb) = CntLeadingZeros(enAccu) - 1;
        }
      } else {
        scaleFactorBandsTotal =
            pSamplingRateInfo->NumberOfScaleFactorBands_Long;
        pSfbOffset = pSamplingRateInfo->ScaleFactorBands_Long;

        for (sfb = 0; sfb < scaleFactorBandsTotal; sfb++) {
          FIXP_DBL enAccu = (FIXP_DBL)(LONG)1;
          int sfbScale =
              (sizeof(LONG) << 3) -
              CntLeadingZeros(pSfbOffset[sfb + 1] - pSfbOffset[sfb]) - 1;
          for (; line < pSfbOffset[sfb + 1]; line++) {
            enAccu += fPow2Div2(*(spectrum + (line >> 3))) >> sfbScale;
          }
          *(sfbEnergy + sfb) = CntLeadingZeros(enAccu) - 1;
        }
      }
      break;
  }
}

static void CConcealment_InterpolateBuffer(FIXP_DBL *spectrum,
                                           SHORT *pSpecScalePrv,
                                           SHORT *pSpecScaleAct,
                                           SHORT *pSpecScaleOut, int *enPrv,
                                           int *enAct, int sfbCnt,
                                           const SHORT *pSfbOffset) {
  int sfb, line = 0;
  int fac_shift;
  int fac_mod;

  for (sfb = 0; sfb < sfbCnt; sfb++) {
    fac_shift =
        enPrv[sfb] - enAct[sfb] + ((*pSpecScaleAct - *pSpecScalePrv) << 1);
    fac_mod = fac_shift & 3;
    fac_shift = (fac_shift >> 2) + 1;
    fac_shift += *pSpecScalePrv - fixMax(*pSpecScalePrv, *pSpecScaleAct);
    fac_shift = fMax(fMin(fac_shift, DFRACT_BITS - 1), -(DFRACT_BITS - 1));

    for (; line < pSfbOffset[sfb + 1]; line++) {
      FIXP_DBL accu = fMult(*(spectrum + line), facMod4Table[fac_mod]);
      *(spectrum + line) = scaleValue(accu, fac_shift);
    }
  }
  *pSpecScaleOut = fixMax(*pSpecScalePrv, *pSpecScaleAct);
}

static INT findEquiFadeFrame(CConcealParams *pConcealCommonData,
                             INT actFadeIndex, int direction) {
  FIXP_SGL *pFactor;
  FIXP_SGL referenceVal;
  FIXP_SGL minDiff = (FIXP_SGL)MAXVAL_SGL;

  INT nextFadeIndex = 0;

  int i;

  if (direction == 0) {
    if (actFadeIndex < 0) {
      referenceVal = (FIXP_SGL)MAXVAL_SGL;
    } else {
      referenceVal = pConcealCommonData->fadeOutFactor[actFadeIndex] >> 1;
    }
    pFactor = pConcealCommonData->fadeInFactor;
  } else {
    if (actFadeIndex < 0) {
      referenceVal = (FIXP_SGL)MAXVAL_SGL;
    } else {
      referenceVal = pConcealCommonData->fadeInFactor[actFadeIndex] >> 1;
    }
    pFactor = pConcealCommonData->fadeOutFactor;
  }

  for (i = 0; i < CONCEAL_MAX_NUM_FADE_FACTORS; i++) {
    FIXP_SGL diff = fixp_abs((pFactor[i] >> 1) - referenceVal);
    if (diff < minDiff) {
      minDiff = diff;
      nextFadeIndex = i;
    }
  }

  if (direction == 0) {
    if (nextFadeIndex > pConcealCommonData->numFadeInFrames) {
      nextFadeIndex = fMax(pConcealCommonData->numFadeInFrames - 1, 0);
    }
    if (((pFactor[nextFadeIndex] >> 1) <= referenceVal) &&
        (nextFadeIndex > 0)) {
      nextFadeIndex -= 1;
    }
  } else {
    if (((pFactor[nextFadeIndex] >> 1) >= referenceVal) &&
        (nextFadeIndex < CONCEAL_MAX_NUM_FADE_FACTORS - 1)) {
      nextFadeIndex += 1;
    }
  }

  return (nextFadeIndex);
}

static void CConcealment_UpdateState(
    CConcealmentInfo *pConcealmentInfo, int frameOk,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const int samplesPerFrame, CAacDecoderChannelInfo *pAacDecoderChannelInfo) {
  CConcealParams *pConcealCommonData = pConcealmentInfo->pConcealParams;

  switch (pConcealCommonData->method) {
    case ConcealMethodNoise: {
      if (pConcealmentInfo->concealState != ConcealState_Ok) {
        if (frameOk) {
          pConcealmentInfo->cntValidFrames += 1;
        } else {
          pConcealmentInfo->cntValidFrames = 0;
        }
      }

      switch (pConcealmentInfo->concealState) {
        case ConcealState_Ok:
          if (!frameOk) {
            pConcealmentInfo->cntFadeFrames = 0;
            pConcealmentInfo->cntValidFrames = 0;
            pConcealmentInfo->attGrpOffset[0] = 0;
            pConcealmentInfo->attGrpOffset[1] = 0;
            pConcealmentInfo->winGrpOffset[0] = 0;
            pConcealmentInfo->winGrpOffset[1] = 0;
            if (pConcealCommonData->numFadeOutFrames > 0) {
              pConcealmentInfo->concealState = ConcealState_Single;
              CConcealment_ApplyFadeOut(
                  0, pConcealmentInfo, pAacDecoderStaticChannelInfo,
                  samplesPerFrame, pAacDecoderChannelInfo);

            } else {
              pConcealmentInfo->concealState = ConcealState_Mute;
            }
          }
          break;

        case ConcealState_Single:
          if (frameOk) {
            pConcealmentInfo->concealState = ConcealState_Ok;
          } else {
            if (pConcealmentInfo->cntFadeFrames >=
                pConcealCommonData->numFadeOutFrames) {
              pConcealmentInfo->concealState = ConcealState_Mute;
            } else {
              pConcealmentInfo->concealState = ConcealState_FadeOut;
              CConcealment_ApplyFadeOut(
                  0, pConcealmentInfo, pAacDecoderStaticChannelInfo,
                  samplesPerFrame, pAacDecoderChannelInfo);
            }
          }
          break;

        case ConcealState_FadeOut:
          if (pConcealmentInfo->cntValidFrames >
              pConcealCommonData->numMuteReleaseFrames) {
            if (pConcealCommonData->numFadeInFrames > 0) {
              pConcealmentInfo->concealState = ConcealState_FadeIn;
              pConcealmentInfo->cntFadeFrames = findEquiFadeFrame(
                  pConcealCommonData, pConcealmentInfo->cntFadeFrames,
                  0);
            } else {
              pConcealmentInfo->concealState = ConcealState_Ok;
            }
          } else {
            if (frameOk) {
              pConcealmentInfo->winGrpOffset[0] = 0;
              pConcealmentInfo->winGrpOffset[1] = 0;
              pConcealmentInfo->attGrpOffset[0] = 0;
              pConcealmentInfo->attGrpOffset[1] = 0;
            }
            if (pConcealmentInfo->cntFadeFrames >=
                pConcealCommonData->numFadeOutFrames) {
              pConcealmentInfo->concealState = ConcealState_Mute;
            } else
            {
              CConcealment_ApplyFadeOut(
                  0, pConcealmentInfo, pAacDecoderStaticChannelInfo,
                  samplesPerFrame, pAacDecoderChannelInfo);
            }
          }
          break;

        case ConcealState_Mute:
          if (pConcealmentInfo->cntValidFrames >
              pConcealCommonData->numMuteReleaseFrames) {
            if (pConcealCommonData->numFadeInFrames > 0) {
              pConcealmentInfo->concealState = ConcealState_FadeIn;
              pConcealmentInfo->cntFadeFrames =
                  pConcealCommonData->numFadeInFrames - 1;
            } else {
              pConcealmentInfo->concealState = ConcealState_Ok;
            }
          } else {
            if (frameOk) {
              pConcealmentInfo->winGrpOffset[0] = 0;
              pConcealmentInfo->winGrpOffset[1] = 0;
              pConcealmentInfo->attGrpOffset[0] = 0;
              pConcealmentInfo->attGrpOffset[1] = 0;
            }
          }
          break;

        case ConcealState_FadeIn:
          pConcealmentInfo->cntFadeFrames -= 1;
          if (frameOk) {
            if (pConcealmentInfo->cntFadeFrames < 0) {
              pConcealmentInfo->concealState = ConcealState_Ok;
            }
          } else {
            if (pConcealCommonData->numFadeOutFrames > 0) {
              pConcealmentInfo->concealState = ConcealState_FadeOut;
              pConcealmentInfo->cntFadeFrames = findEquiFadeFrame(
                  pConcealCommonData, pConcealmentInfo->cntFadeFrames + 1,
                  1);
              pConcealmentInfo->winGrpOffset[0] = 0;
              pConcealmentInfo->winGrpOffset[1] = 0;
              pConcealmentInfo->attGrpOffset[0] = 0;
              pConcealmentInfo->attGrpOffset[1] = 0;

              pConcealmentInfo
                  ->cntFadeFrames--;
              CConcealment_ApplyFadeOut(
                  0, pConcealmentInfo, pAacDecoderStaticChannelInfo,
                  samplesPerFrame, pAacDecoderChannelInfo);
            } else {
              pConcealmentInfo->concealState = ConcealState_Mute;
            }
          }
          break;

        default:
          FDK_ASSERT(0);
          break;
      }
    } break;

    case ConcealMethodInter:
    case ConcealMethodTonal: {
      if (pConcealmentInfo->concealState != ConcealState_Ok) {
        if (pConcealmentInfo->prevFrameOk[1] ||
            (pConcealmentInfo->prevFrameOk[0] &&
             !pConcealmentInfo->prevFrameOk[1] && frameOk)) {
          pConcealmentInfo->cntValidFrames += 1;
        } else {
          pConcealmentInfo->cntValidFrames = 0;
        }
      }

      switch (pConcealmentInfo->concealState) {
        case ConcealState_Ok:
          if (!(pConcealmentInfo->prevFrameOk[1] ||
                (pConcealmentInfo->prevFrameOk[0] &&
                 !pConcealmentInfo->prevFrameOk[1] && frameOk))) {
            if (pConcealCommonData->numFadeOutFrames > 0) {
              pConcealmentInfo->concealState = ConcealState_FadeOut;
            } else {
              pConcealmentInfo->concealState = ConcealState_Mute;
            }
            pConcealmentInfo->cntFadeFrames = 0;
            pConcealmentInfo->cntValidFrames = 0;
          }
          break;

        case ConcealState_Single:
          pConcealmentInfo->concealState = ConcealState_Ok;
          break;

        case ConcealState_FadeOut:
          pConcealmentInfo->cntFadeFrames += 1;

          if (pConcealmentInfo->cntValidFrames >
              pConcealCommonData->numMuteReleaseFrames) {
            if (pConcealCommonData->numFadeInFrames > 0) {
              pConcealmentInfo->concealState = ConcealState_FadeIn;
              pConcealmentInfo->cntFadeFrames = findEquiFadeFrame(
                  pConcealCommonData, pConcealmentInfo->cntFadeFrames - 1,
                  0);
            } else {
              pConcealmentInfo->concealState = ConcealState_Ok;
            }
          } else {
            if (pConcealmentInfo->cntFadeFrames >=
                pConcealCommonData->numFadeOutFrames) {
              pConcealmentInfo->concealState = ConcealState_Mute;
            }
          }
          break;

        case ConcealState_Mute:
          if (pConcealmentInfo->cntValidFrames >
              pConcealCommonData->numMuteReleaseFrames) {
            if (pConcealCommonData->numFadeInFrames > 0) {
              pConcealmentInfo->concealState = ConcealState_FadeIn;
              pConcealmentInfo->cntFadeFrames =
                  pConcealCommonData->numFadeInFrames - 1;
            } else {
              pConcealmentInfo->concealState = ConcealState_Ok;
            }
          }
          break;

        case ConcealState_FadeIn:
          pConcealmentInfo->cntFadeFrames -=
              1;

          if (frameOk || pConcealmentInfo->prevFrameOk[1]) {
            if (pConcealmentInfo->cntFadeFrames < 0) {
              pConcealmentInfo->concealState = ConcealState_Ok;
            }
          } else {
            if (pConcealCommonData->numFadeOutFrames > 0) {
              pConcealmentInfo->concealState = ConcealState_FadeOut;
              pConcealmentInfo->cntFadeFrames = findEquiFadeFrame(
                  pConcealCommonData, pConcealmentInfo->cntFadeFrames + 1,
                  1);
            } else {
              pConcealmentInfo->concealState = ConcealState_Mute;
            }
          }
          break;
      }
    } break;

    default:
      break;
  }
}

static void CConcealment_ApplyRandomSign(int randomPhase, FIXP_DBL *spec,
                                         int samplesPerFrame) {
  int i;
  USHORT packedSign = 0;

  packedSign = AacDec_randomSign[randomPhase >> 4];
  packedSign >>= (randomPhase & 0xf);

  for (i = 0; i < samplesPerFrame; i++) {
    if ((randomPhase & 0xf) == 0) {
      packedSign = AacDec_randomSign[randomPhase >> 4];
    }

    if (packedSign & 0x1) {
      spec[i] = -fMax(spec[i], (FIXP_DBL)(MINVAL_DBL + 1));
    }
    packedSign >>= 1;

    randomPhase = (randomPhase + 1) & (AAC_NF_NO_RANDOM_VAL - 1);
  }
}

int CConcealment_GetLastFrameOk(CConcealmentInfo *hConcealmentInfo,
                                const int fBeforeApply) {
  int prevFrameOk = 1;

  if (hConcealmentInfo != NULL) {
    prevFrameOk = hConcealmentInfo->prevFrameOk[fBeforeApply & 0x1];
  }

  return prevFrameOk;
}

UINT CConcealment_GetDelay(CConcealParams *pConcealCommonData) {
  UINT frameDelay = 0;

  if (pConcealCommonData != NULL) {
    switch (pConcealCommonData->method) {
      case ConcealMethodTonal:
      case ConcealMethodInter:
        frameDelay = 1;
        break;
      default:
        break;
    }
  }

  return frameDelay;
}

static int CConcealment_ApplyFadeOut(
    int mode, CConcealmentInfo *pConcealmentInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const int samplesPerFrame, CAacDecoderChannelInfo *pAacDecoderChannelInfo) {
  int srcWin, dstWin, numWindows = 1;
  int windowLen = samplesPerFrame;
  int srcGrpStart = 0;
  int winIdxStride = 1;
  int numWinGrpPerFac, attIdx, attIdxStride;
  int i;
  int appliedProcessing = 0;

  CIcsInfo *pIcsInfo = &pAacDecoderChannelInfo->icsInfo;
  FIXP_DBL *pSpectralCoefficient =
      SPEC_LONG(pAacDecoderChannelInfo->pSpectralCoefficient);
  SHORT *pSpecScale = pAacDecoderChannelInfo->specScale;

  if (pConcealmentInfo->lastRenderMode == AACDEC_RENDER_LPD) {
    switch (pAacDecoderStaticChannelInfo->last_lpd_mode) {
      case 1:
        numWindows = 4;
        srcGrpStart = 3;
        windowLen = samplesPerFrame >> 2;
        break;
      case 2:
        numWindows = 2;
        srcGrpStart = 1;
        windowLen = samplesPerFrame >> 1;
        winIdxStride = 2;
        break;
      case 3:
        numWindows = 1;
        srcGrpStart = 0;
        windowLen = samplesPerFrame;
        winIdxStride = 4;
        break;
    }
    pConcealmentInfo->lastWinGrpLen = 1;
  } else {
    pIcsInfo->WindowShape = pConcealmentInfo->windowShape;
    pIcsInfo->WindowSequence = pConcealmentInfo->windowSequence;

    if (pConcealmentInfo->windowSequence == BLOCK_SHORT) {
      numWindows = 8;
      windowLen = samplesPerFrame >> 3;
      srcGrpStart = numWindows - pConcealmentInfo->lastWinGrpLen;
    }
  }

  attIdxStride =
      fMax(1, (int)(numWindows / (pConcealmentInfo->lastWinGrpLen + 1)));

  attIdx = pConcealmentInfo->cntFadeFrames;
  numWinGrpPerFac = pConcealmentInfo->attGrpOffset[mode];
  srcWin = srcGrpStart + pConcealmentInfo->winGrpOffset[mode];

  FDK_ASSERT((srcGrpStart * windowLen + windowLen) <= samplesPerFrame);
  FDK_ASSERT((srcWin * windowLen + windowLen) <= 1024);

  for (dstWin = 0; dstWin < numWindows; dstWin += 1) {
    FIXP_CNCL *pCncl =
        pConcealmentInfo->spectralCoefficient + (srcWin * windowLen);
    FIXP_DBL *pOut = pSpectralCoefficient + (dstWin * windowLen);

    if (mode == 1) {
      if (attIdx > pConcealmentInfo->pConcealParams->numFadeOutFrames) {
        FDKmemclear(pCncl, sizeof(FIXP_DBL) * windowLen);
      }

      for (i = 0; i < windowLen; i++) {
        pOut[i] = pCncl[i];
      }

      CConcealment_ApplyRandomSign(pConcealmentInfo->iRandomPhase, pOut,
                                   windowLen);

      pConcealmentInfo->iRandomPhase =
          (pConcealmentInfo->iRandomPhase + 1) & (AAC_NF_NO_RANDOM_VAL - 1);

      pSpecScale[dstWin * winIdxStride] =
          pConcealmentInfo->specScale[srcWin * winIdxStride];
    }

    srcWin += 1;

    if (srcWin >= numWindows) {
      srcWin = srcGrpStart;
      numWinGrpPerFac += 1;
      if (numWinGrpPerFac >= attIdxStride) {
        numWinGrpPerFac = 0;
        attIdx += 1;
      }
    }
  }

  pConcealmentInfo->winGrpOffset[mode] = srcWin - srcGrpStart;
  FDK_ASSERT((pConcealmentInfo->winGrpOffset[mode] >= 0) &&
             (pConcealmentInfo->winGrpOffset[mode] < 8));
  pConcealmentInfo->attGrpOffset[mode] = numWinGrpPerFac;
  FDK_ASSERT((pConcealmentInfo->attGrpOffset[mode] >= 0) &&
             (pConcealmentInfo->attGrpOffset[mode] < attIdxStride));

  if (mode == 0) {
    pConcealmentInfo->cntFadeFrames = attIdx;
  }

  appliedProcessing = 1;

  return appliedProcessing;
}

INT CConcealment_TDFading(
    int len, CAacDecoderStaticChannelInfo **ppAacDecoderStaticChannelInfo,
    const INT aacOutDataHeadroom, PCM_DEC *pcmdata, PCM_DEC *pcmdata_1) {
  FIXP_DBL fadeStop, attMute = (FIXP_DBL)0;
  int idx = 0, ii;
  CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo =
      *ppAacDecoderStaticChannelInfo;
  CConcealmentInfo *pConcealmentInfo =
      &pAacDecoderStaticChannelInfo->concealmentInfo;
  CConcealParams *pConcealParams = pConcealmentInfo->pConcealParams;
  const CConcealmentState concealState = pConcealmentInfo->concealState;
  TDfadingType fadingType;
  FIXP_DBL fadingStations[9] = {0};
  int fadingSteps[8] = {0};
  const FIXP_DBL fadeStart =
      pConcealmentInfo
          ->fade_old;
  FIXP_SGL *fadeFactor = pConcealParams->fadeOutFactor;
  const INT cntFadeFrames = pConcealmentInfo->cntFadeFrames;
  int TDFadeOutStopBeforeMute = 1;
  int TDFadeInStopBeforeFullLevel = 1;

  switch (concealState) {
    case ConcealState_Single:
    case ConcealState_Mute:
    case ConcealState_FadeOut:
      idx = (pConcealParams->method == ConcealMethodNoise) ? cntFadeFrames - 1
                                                           : cntFadeFrames;
      fadingType = FADE_TIMEDOMAIN;

      if (concealState == ConcealState_Mute ||
          (cntFadeFrames + TDFadeOutStopBeforeMute) >
              pConcealmentInfo->pConcealParams->numFadeOutFrames) {
        fadingType = FADE_TIMEDOMAIN_TOSPECTRALMUTE;
      }

      break;
    case ConcealState_FadeIn:
      idx = cntFadeFrames;
      idx -= TDFadeInStopBeforeFullLevel;
      FDK_FALLTHROUGH;
    case ConcealState_Ok:
      fadeFactor = pConcealParams->fadeInFactor;
      idx = (concealState == ConcealState_Ok) ? -1 : idx;
      fadingType = (pConcealmentInfo->concealState_old == ConcealState_Mute)
                       ? FADE_TIMEDOMAIN_FROMSPECTRALMUTE
                       : FADE_TIMEDOMAIN;
      break;
    default:
      FDK_ASSERT(0);
      fadingType = FADE_TIMEDOMAIN_TOSPECTRALMUTE;
      break;
  }

  switch (fadingType) {
    case FADE_TIMEDOMAIN_FROMSPECTRALMUTE:
      fadeStop =
          (idx < 0) ? (FIXP_DBL)MAXVAL_DBL : FX_SGL2FX_DBL(fadeFactor[idx]);
      if (pConcealmentInfo->pConcealParams->numFadeInFrames == 0) {
        fadingSteps[0] = 1;
        break;
      }
      CConcealment_TDFading_doLinearFadingSteps(&fadingSteps[0]);
      break;
    case FADE_TIMEDOMAIN:
      fadeStop =
          (idx < 0) ? (FIXP_DBL)MAXVAL_DBL : FX_SGL2FX_DBL(fadeFactor[idx]);
      CConcealment_TDFading_doLinearFadingSteps(&fadingSteps[0]);
      break;
    case FADE_TIMEDOMAIN_TOSPECTRALMUTE:
      fadeStop = attMute;
      if (pConcealmentInfo->pConcealParams->numFadeOutFrames == 0) {
        fadingSteps[0] = 1;
        break;
      }
      CConcealment_TDFading_doLinearFadingSteps(&fadingSteps[0]);
      break;
  }

  len >>= 3;
  CConcealment_TDFadeFillFadingStations(fadingStations, fadingSteps, fadeStop,
                                        fadeStart, fadingType);

  if ((fadingStations[8] != (FIXP_DBL)MAXVAL_DBL) ||
      (fadingStations[7] != (FIXP_DBL)MAXVAL_DBL) ||
      (fadingStations[6] != (FIXP_DBL)MAXVAL_DBL) ||
      (fadingStations[5] != (FIXP_DBL)MAXVAL_DBL) ||
      (fadingStations[4] != (FIXP_DBL)MAXVAL_DBL) ||
      (fadingStations[3] != (FIXP_DBL)MAXVAL_DBL) ||
      (fadingStations[2] != (FIXP_DBL)MAXVAL_DBL) ||
      (fadingStations[1] != (FIXP_DBL)MAXVAL_DBL) ||
      (fadingStations[0] !=
       (FIXP_DBL)MAXVAL_DBL))
  {
    int start = 0;
    for (ii = 0; ii < 8; ii++) {
      CConcealment_TDFadePcmAtt(start, len, fadingStations[ii],
                                fadingStations[ii + 1], pcmdata);
      start += len;
    }
  }
  CConcealment_TDNoise_Apply(pConcealmentInfo, len, aacOutDataHeadroom,
                             pcmdata);

  pConcealmentInfo->lastFadingType = fadingType;
  pConcealmentInfo->fade_old = fadeStop;
  pConcealmentInfo->concealState_old = concealState;

  return 1;
}

static void CConcealment_TDFadePcmAtt(int start, int len, FIXP_DBL fadeStart,
                                      FIXP_DBL fadeStop, PCM_DEC *pcmdata) {
  int i;
  FIXP_DBL dStep;
  FIXP_DBL dGain;
  FIXP_DBL dGain_apply;

  dGain = fadeStart;
  dStep = (FIXP_DBL)((int)((fadeStart >> 1) - (fadeStop >> 1)) / len) << 1;

  for (i = start; i < (start + len); i++) {
    dGain -= dStep;
    dGain_apply = fMax((FIXP_DBL)0, dGain);
    pcmdata[i] = FIXP_DBL2PCM_DEC(fMult(pcmdata[i], dGain_apply));
  }
}

static void CConcealment_TDFadeFillFadingStations(FIXP_DBL *fadingStations,
                                                  int *fadingSteps,
                                                  FIXP_DBL fadeStop,
                                                  FIXP_DBL fadeStart,
                                                  TDfadingType fadingType) {
  int i;
  INT fadingSteps_sum = 0;
  INT fadeDiff;

  fadingSteps_sum = fadingSteps[0] + fadingSteps[1] + fadingSteps[2] +
                    fadingSteps[3] + fadingSteps[4] + fadingSteps[5] +
                    fadingSteps[6] + fadingSteps[7];
  fadeDiff = ((INT)(fadeStop - fadeStart) / fMax(fadingSteps_sum, (INT)1));
  fadingStations[0] = fadeStart;
  for (i = 1; i < 8; i++) {
    fadingStations[i] =
        fadingStations[i - 1] + (FIXP_DBL)(fadeDiff * fadingSteps[i - 1]);
  }
  fadingStations[8] = fadeStop;
}

static void CConcealment_TDFading_doLinearFadingSteps(int *fadingSteps) {
  fadingSteps[0] = fadingSteps[1] = fadingSteps[2] = fadingSteps[3] =
      fadingSteps[4] = fadingSteps[5] = fadingSteps[6] = fadingSteps[7] = 1;
}

static int CConcealment_TDNoise_Random(ULONG *seed) {
  *seed = (ULONG)(((UINT64)(*seed) * 69069) + 5);
  return (int)(*seed);
}

static void CConcealment_TDNoise_Apply(CConcealmentInfo *const pConcealmentInfo,
                                       const int len,
                                       const INT aacOutDataHeadroom,
                                       PCM_DEC *const pcmdata) {
  PCM_DEC *states = pConcealmentInfo->TDNoiseStates;
  PCM_DEC noiseVal;
  FIXP_DBL noiseValLong;
  FIXP_SGL *coef = pConcealmentInfo->TDNoiseCoef;
  FIXP_DBL TDNoiseAtt;
  ULONG seed = pConcealmentInfo->TDNoiseSeed =
      (ULONG)CConcealment_TDNoise_Random(&pConcealmentInfo->TDNoiseSeed) + 1;

  TDNoiseAtt = pConcealmentInfo->pConcealParams->comfortNoiseLevel;

  int ii;

  if ((pConcealmentInfo->concealState != ConcealState_Ok ||
       pConcealmentInfo->concealState_old != ConcealState_Ok) &&
      TDNoiseAtt != (FIXP_DBL)0) {
    for (ii = 0; ii < (len << 3); ii++) {
      states[2] = states[1];
      states[1] = states[0];
      states[0] =
          FIXP_DBL2PCM_DEC((FIXP_DBL)CConcealment_TDNoise_Random(&seed));
      noiseValLong = fMult(states[0], coef[0]) + fMult(states[1], coef[1]) +
                     fMult(states[2], coef[2]);
      noiseVal = FIXP_DBL2PCM_DEC(fMult(noiseValLong, TDNoiseAtt) >>
                                  aacOutDataHeadroom);

      if (noiseVal > (PCM_DEC)0 &&
          pcmdata[ii] > (PCM_DEC)MAXVAL_PCM_DEC - noiseVal) {
        noiseVal = noiseVal * (PCM_DEC)-1;
      } else if (noiseVal < (PCM_DEC)0 &&
                 pcmdata[ii] < (PCM_DEC)MINVAL_PCM_DEC - noiseVal) {
        noiseVal = noiseVal * (PCM_DEC)-1;
      }

      pcmdata[ii] += noiseVal;
    }
  }
}
