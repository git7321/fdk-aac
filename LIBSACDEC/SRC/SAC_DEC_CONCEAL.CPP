/*********************** MPEG surround decoder library *************************

   Author(s):   Christian Ertel, Christian Griebel

   Description: SAC Dec error concealment

*******************************************************************************/

#include "sac_dec_conceal.h"

void SpatialDecConcealment_Init(SpatialDecConcealmentInfo *info,
                                const UINT resetFlags) {
  FDK_ASSERT(info != NULL);

  if (resetFlags & MPEGS_CONCEAL_RESET_STATE) {
    info->concealState = SpatialDecConcealState_Init;
  }

  if (resetFlags & MPEGS_CONCEAL_RESET_PARAMETER) {
    info->concealParams.method = MPEGS_CONCEAL_DEFAULT_METHOD;
    info->concealParams.numKeepFrames = MPEGS_CONCEAL_DEFAULT_NUM_KEEP_FRAMES;
    info->concealParams.numFadeOutFrames =
        MPEGS_CONCEAL_DEFAULT_FADE_OUT_SLOPE_LENGTH;
    info->concealParams.numFadeInFrames =
        MPEGS_CONCEAL_DEFAULT_FADE_IN_SLOPE_LENGTH;
    info->concealParams.numReleaseFrames =
        MPEGS_CONCEAL_DEFAULT_NUM_RELEASE_FRAMES;
  }

  return;
}

int SpatialDecConcealment_Apply(
    SpatialDecConcealmentInfo *info,
    const SCHAR (*cmpIdxData)[MAX_PARAMETER_BANDS], SCHAR **diffIdxData,
    SCHAR *
        idxPrev,
    SCHAR *bsXXXDataMode, const int startBand, const int stopBand,
    const SCHAR defaultValue, const int paramType, const int numParamSets) {
  int appliedProcessing = 0;
  int band, dataMode = -1;

  FDK_ASSERT(info != NULL);
  FDK_ASSERT(cmpIdxData != NULL);
  FDK_ASSERT(idxPrev != NULL);
  FDK_ASSERT(bsXXXDataMode != NULL);

  switch (info->concealState) {
    case SpatialDecConcealState_Init:
      dataMode = 0;
      break;

    case SpatialDecConcealState_Ok:
      break;

    case SpatialDecConcealState_Keep:
      dataMode = 1;
      break;

    case SpatialDecConcealState_FadeToDefault: {
      FIXP_DBL fac = fDivNorm(info->cntStateFrames + 1,
                              info->concealParams.numFadeOutFrames + 1);

      for (band = startBand; band < stopBand; band += 1) {
        idxPrev[band] =
            fMultI(fac, defaultValue - idxPrev[band]) + idxPrev[band];
      }
      dataMode = 1;
      appliedProcessing = 1;
    } break;

    case SpatialDecConcealState_Default:
      for (band = startBand; band < stopBand; band += 1) {
        idxPrev[band] = defaultValue;
      }
      dataMode = 1;
      appliedProcessing = 1;
      break;

    case SpatialDecConcealState_FadeFromDefault: {
      FIXP_DBL fac = fDivNorm(info->cntValidFrames + 1,
                              info->concealParams.numFadeInFrames + 1);

      for (band = startBand; band < stopBand; band += 1) {
        idxPrev[band] =
            fMultI(fac, cmpIdxData[numParamSets - 1][band] - defaultValue) +
            defaultValue;
      }
      dataMode = 1;
      appliedProcessing = 1;
    } break;

    default:
      FDK_ASSERT(0);
      break;
  }

  if (dataMode >= 0) {
    int i;
    for (i = 0; i < numParamSets; i += 1) {
      bsXXXDataMode[i] = dataMode;
      if (diffIdxData != NULL) {
        for (band = startBand; band < stopBand; band += 1) {
          diffIdxData[i][band] = 0;
        }
      }
    }
  }

  return appliedProcessing;
}

void SpatialDecConcealment_UpdateState(SpatialDecConcealmentInfo *info,
                                       const int frameOk) {
  FDK_ASSERT(info != NULL);

  if (frameOk) {
    info->cntValidFrames += 1;
  } else {
    info->cntValidFrames = 0;
  }

  switch (info->concealState) {
    case SpatialDecConcealState_Init:
      if (frameOk) {
        info->concealState = SpatialDecConcealState_Ok;
        info->cntStateFrames = 0;
      }
      break;

    case SpatialDecConcealState_Ok:
      if (!frameOk) {
        info->concealState = SpatialDecConcealState_Keep;
        info->cntStateFrames = 0;
      }
      break;

    case SpatialDecConcealState_Keep:
      info->cntStateFrames += 1;
      if (frameOk) {
        info->concealState = SpatialDecConcealState_Ok;
      } else {
        if (info->cntStateFrames >= info->concealParams.numKeepFrames) {
          if (info->concealParams.numFadeOutFrames == 0) {
            info->concealState = SpatialDecConcealState_Default;
          } else {
            info->concealState = SpatialDecConcealState_FadeToDefault;
            info->cntStateFrames = 0;
          }
        }
      }
      break;

    case SpatialDecConcealState_FadeToDefault:
      info->cntStateFrames += 1;
      if (info->cntValidFrames > 0) {
        info->concealState = SpatialDecConcealState_FadeFromDefault;
        info->cntStateFrames = 0;
      } else {
        if (info->cntStateFrames >= info->concealParams.numFadeOutFrames) {
          info->concealState = SpatialDecConcealState_Default;
        }
      }
      break;

    case SpatialDecConcealState_Default:
      if (info->cntValidFrames > 0) {
        if (info->concealParams.numFadeInFrames == 0) {
          info->concealState = SpatialDecConcealState_Ok;
        } else {
          info->concealState = SpatialDecConcealState_FadeFromDefault;
          info->cntValidFrames = 0;
        }
      }
      break;

    case SpatialDecConcealState_FadeFromDefault:
      info->cntValidFrames += 1;
      if (frameOk) {
        if (info->cntValidFrames >= info->concealParams.numFadeInFrames) {
          info->concealState = SpatialDecConcealState_Ok;
        }
      } else {
        info->concealState = SpatialDecConcealState_FadeToDefault;
        info->cntStateFrames = 0;
      }
      break;

    default:
      FDK_ASSERT(0);
      break;
  }
}

SACDEC_ERROR SpatialDecConcealment_SetParam(SpatialDecConcealmentInfo *self,
                                            const SAC_DEC_CONCEAL_PARAM param,
                                            const INT value) {
  SACDEC_ERROR err = MPS_OK;

  switch (param) {
    case SAC_DEC_CONCEAL_METHOD:
      switch ((SpatialDecConcealmentMethod)value) {
        case SAC_DEC_CONCEAL_WITH_ZERO_VALUED_OUTPUT:
        case SAC_DEC_CONCEAL_BY_FADING_PARAMETERS:
          break;
        default:
          err = MPS_INVALID_PARAMETER;
          goto bail;
      }
      if (self != NULL) {
        self->concealParams.method = (SpatialDecConcealmentMethod)value;
      } else {
        err = MPS_INVALID_HANDLE;
        goto bail;
      }
      break;
    case SAC_DEC_CONCEAL_NUM_KEEP_FRAMES:
      if (value < 0) {
        err = MPS_INVALID_PARAMETER;
        goto bail;
      }
      if (self != NULL) {
        self->concealParams.numKeepFrames = (UINT)value;
      } else {
        err = MPS_INVALID_HANDLE;
        goto bail;
      }
      break;
    case SAC_DEC_CONCEAL_FADE_OUT_SLOPE_LENGTH:
      if (value < 0) {
        err = MPS_INVALID_PARAMETER;
        goto bail;
      }
      if (self != NULL) {
        self->concealParams.numFadeOutFrames = (UINT)value;
      } else {
        err = MPS_INVALID_HANDLE;
        goto bail;
      }
      break;
    case SAC_DEC_CONCEAL_FADE_IN_SLOPE_LENGTH:
      if (value < 0) {
        err = MPS_INVALID_PARAMETER;
        goto bail;
      }
      if (self != NULL) {
        self->concealParams.numFadeInFrames = (UINT)value;
      } else {
        err = MPS_INVALID_HANDLE;
        goto bail;
      }
      break;
    case SAC_DEC_CONCEAL_NUM_RELEASE_FRAMES:
      if (value < 0) {
        err = MPS_INVALID_PARAMETER;
        goto bail;
      }
      if (self != NULL) {
        self->concealParams.numReleaseFrames = (UINT)value;
      } else {
        err = MPS_INVALID_HANDLE;
        goto bail;
      }
      break;
    default:
      err = MPS_INVALID_PARAMETER;
      goto bail;
  }

bail:
  return err;
}
