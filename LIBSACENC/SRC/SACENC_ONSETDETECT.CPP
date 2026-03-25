/*********************** MPEG surround encoder library *************************

   Author(s):   Max Neuendorf

   Description: Encoder Library Interface
                Detect Onset in current frame

*******************************************************************************/

#include "sacenc_onsetdetect.h"
#include "genericStds.h"
#include "sacenc_vectorfunctions.h"

#define SPACE_ONSET_THRESHOLD (3.0)
#define SPACE_ONSET_THRESHOLD_SF (3)
#define SPACE_ONSET_THRESHOLD_SQUARE                                        \
  (FL2FXCONST_DBL((1.0 / (SPACE_ONSET_THRESHOLD * SPACE_ONSET_THRESHOLD)) * \
                  (float)(1 << SPACE_ONSET_THRESHOLD_SF)))

struct ONSET_DETECT {
  INT maxTimeSlots;
  INT minTransientDistance;
  INT avgEnergyDistance;
  INT lowerBoundOnsetDetection;
  INT upperBoundOnsetDetection;
  FIXP_DBL *pEnergyHist__FDK;
  SCHAR *pEnergyHistScale;
  SCHAR avgEnergyDistanceScale;
};

FDK_SACENC_ERROR fdk_sacenc_onsetDetect_Open(HANDLE_ONSET_DETECT *phOnset,
                                             const UINT maxTimeSlots) {
  FDK_SACENC_ERROR error = SACENC_OK;
  HANDLE_ONSET_DETECT hOnset = NULL;

  if (NULL == phOnset) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDK_ALLOCATE_MEMORY_1D(hOnset, 1, struct ONSET_DETECT);
    FDK_ALLOCATE_MEMORY_1D(hOnset->pEnergyHist__FDK, 16 + maxTimeSlots,
                           FIXP_DBL);
    FDK_ALLOCATE_MEMORY_1D(hOnset->pEnergyHistScale, 16 + maxTimeSlots, SCHAR);

    hOnset->maxTimeSlots = maxTimeSlots;
    hOnset->minTransientDistance =
        8;
    hOnset->avgEnergyDistance = 16;

    hOnset->avgEnergyDistanceScale = 4;
    *phOnset = hOnset;
  }
  return error;

bail:
  fdk_sacenc_onsetDetect_Close(&hOnset);
  return ((SACENC_OK == error) ? SACENC_MEMORY_ERROR : error);
}

FDK_SACENC_ERROR fdk_sacenc_onsetDetect_Init(
    HANDLE_ONSET_DETECT hOnset,
    const ONSET_DETECT_CONFIG *const pOnsetDetectConfig, const UINT initFlags) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((NULL == hOnset) || (pOnsetDetectConfig == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    if ((pOnsetDetectConfig->maxTimeSlots > hOnset->maxTimeSlots) ||
        (pOnsetDetectConfig->upperBoundOnsetDetection <
         hOnset->lowerBoundOnsetDetection)) {
      error = SACENC_INVALID_CONFIG;
      goto bail;
    }

    hOnset->maxTimeSlots = pOnsetDetectConfig->maxTimeSlots;
    hOnset->lowerBoundOnsetDetection =
        pOnsetDetectConfig->lowerBoundOnsetDetection;
    hOnset->upperBoundOnsetDetection =
        pOnsetDetectConfig->upperBoundOnsetDetection;

    hOnset->minTransientDistance =
        8;
    hOnset->avgEnergyDistance = 16;

    hOnset->avgEnergyDistanceScale = 4;

    if (initFlags) {
      int i;
      for (i = 0; i < hOnset->avgEnergyDistance + hOnset->maxTimeSlots; i++)
        hOnset->pEnergyHistScale[i] = -(DFRACT_BITS - 3);

      FDKmemset_flex(
          hOnset->pEnergyHist__FDK,
          FL2FXCONST_DBL(SACENC_FLOAT_EPSILON * (1 << (DFRACT_BITS - 3))),
          hOnset->avgEnergyDistance + hOnset->maxTimeSlots);
    }
  }

bail:
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_onsetDetect_Close(HANDLE_ONSET_DETECT *phOnset) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((NULL != phOnset) && (NULL != *phOnset)) {
    if (NULL != (*phOnset)->pEnergyHist__FDK) {
      FDKfree((*phOnset)->pEnergyHist__FDK);
    }
    (*phOnset)->pEnergyHist__FDK = NULL;

    if (NULL != (*phOnset)->pEnergyHistScale) {
      FDKfree((*phOnset)->pEnergyHistScale);
    }
    (*phOnset)->pEnergyHistScale = NULL;
    FDKfree(*phOnset);
    *phOnset = NULL;
  }
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_onsetDetect_Update(HANDLE_ONSET_DETECT hOnset,
                                               const INT timeSlots) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == hOnset) {
    error = SACENC_INVALID_HANDLE;
  } else {
    if (timeSlots > hOnset->maxTimeSlots) {
      error = SACENC_INVALID_CONFIG;
    } else {
      int i;
      for (i = 0; i < hOnset->avgEnergyDistance; i++) {
        hOnset->pEnergyHist__FDK[i] = hOnset->pEnergyHist__FDK[i + timeSlots];
        hOnset->pEnergyHistScale[i] = hOnset->pEnergyHistScale[i + timeSlots];
      }

      FDKmemset_flex(&hOnset->pEnergyHist__FDK[hOnset->avgEnergyDistance],
                     FL2FXCONST_DBL(SACENC_FLOAT_EPSILON), timeSlots);
    }
  }
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_onsetDetect_Apply(
    HANDLE_ONSET_DETECT hOnset, const INT nTimeSlots, const INT nHybridBands,
    FIXP_DPK *const *const ppHybridData__FDK, const INT hybridDataScale,
    const INT prevPos, INT pTransientPos[MAX_NUM_TRANS]) {
  FDK_SACENC_ERROR error = SACENC_OK;

  C_ALLOC_SCRATCH_START(envs, FIXP_DBL, (16 + MAX_TIME_SLOTS))
  FDKmemclear(envs, (16 + MAX_TIME_SLOTS) * sizeof(FIXP_DBL));

  if ((hOnset == NULL) || (pTransientPos == NULL) ||
      (ppHybridData__FDK == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    int i, ts, trCnt, currPos;

    if ((nTimeSlots < 0) || (nTimeSlots > hOnset->maxTimeSlots) ||
        (hOnset->lowerBoundOnsetDetection < -1) ||
        (hOnset->upperBoundOnsetDetection > nHybridBands)) {
      error = SACENC_INVALID_CONFIG;
      goto bail;
    }

    const int lowerBoundOnsetDetection = hOnset->lowerBoundOnsetDetection;
    const int upperBoundOnsetDetection = hOnset->upperBoundOnsetDetection;
    const int M = hOnset->avgEnergyDistance;

    {
      SCHAR *envScale = hOnset->pEnergyHistScale;
      FIXP_DBL *env = hOnset->pEnergyHist__FDK;
      const FIXP_DBL threshold_square = SPACE_ONSET_THRESHOLD_SQUARE;

      trCnt = 0;

      FDKmemset_flex(pTransientPos, -1, MAX_NUM_TRANS);

      if (prevPos > 0) {
        currPos = FDKmax(nTimeSlots,
                         prevPos - nTimeSlots + hOnset->minTransientDistance);
      } else {
        currPos = nTimeSlots;
      }

      int outScale;
      int inScale = 3;
      for (ts = 0; ts < nTimeSlots; ts++) {
        env[M + ts] = sumUpCplxPow2(
            &ppHybridData__FDK[ts][lowerBoundOnsetDetection + 1],
            SUM_UP_DYNAMIC_SCALE, inScale, &outScale,
            upperBoundOnsetDetection - lowerBoundOnsetDetection - 1);
        envScale[M + ts] = outScale + (hybridDataScale << 1);
      }

      SCHAR maxScale = -(DFRACT_BITS - 1);
      for (i = 0; i < (nTimeSlots + M); i++) {
        maxScale = fixMax(maxScale, envScale[i]);
      }

      for (i = 0; i < (nTimeSlots + M); i++) {
        envs[i] = env[i] >> fixMin((maxScale - envScale[i]), (DFRACT_BITS - 1));
      }

      FIXP_DBL maxVal = FL2FXCONST_DBL(0.0f);
      for (i = 0; i < (nTimeSlots + M); i++) {
        maxVal |= fAbs(envs[i]);
      }

      int s = fixMax(0, CntLeadingZeros(maxVal) - 1);

      for (i = 0; i < (nTimeSlots + M); i++) {
        envs[i] = envs[i] << s;
      }

      int currPosPrev = currPos;
      FIXP_DBL p1, p2;
      p2 = FL2FXCONST_DBL(0.0f);
      for (; (currPos < (nTimeSlots << 1)) && (trCnt < MAX_NUM_TRANS);
           currPos++) {
        p1 = fMultDiv2(envs[currPos - nTimeSlots + M], threshold_square) >>
             (SPACE_ONSET_THRESHOLD_SF - 1);

        if (currPosPrev == (currPos - 1)) {
          p2 -= (envs[currPosPrev - nTimeSlots] >>
                 (int)hOnset->avgEnergyDistanceScale);
          p2 += (envs[currPos - nTimeSlots + M - 1] >>
                 (int)hOnset->avgEnergyDistanceScale);
        } else {
          p2 = FL2FXCONST_DBL(0.0f);
          for (ts = 0; ts < M; ts++) {
            p2 += (envs[currPos - nTimeSlots + ts] >>
                   (int)hOnset->avgEnergyDistanceScale);
          }
        }
        currPosPrev = currPos;

        {
          if (p1 > p2) {
            pTransientPos[trCnt++] = currPos;
            currPos += hOnset->minTransientDistance;
          }
        }
      }
    }
  }
bail:

  C_ALLOC_SCRATCH_END(envs, FIXP_DBL, (16 + MAX_TIME_SLOTS))

  return error;
}
