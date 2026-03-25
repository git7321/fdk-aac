/**************************** PCM utility library ******************************

   Author(s):   Matthias Neusinger

   Description: Hard limiter for clipping prevention

*******************************************************************************/

#include "limiter.h"
#include "FDK_core.h"
#include "version.h"

#define TDLIMIT_LIB_TITLE "TD Limiter Lib"

TDLimiterPtr pcmLimiter_Create(unsigned int maxAttackMs, unsigned int releaseMs,
                               FIXP_DBL threshold, unsigned int maxChannels,
                               UINT maxSampleRate) {
  TDLimiterPtr limiter = NULL;
  unsigned int attack, release;
  FIXP_DBL attackConst, releaseConst, exponent;
  INT e_ans;

  attack = (unsigned int)(maxAttackMs * maxSampleRate / 1000);
  release = (unsigned int)(releaseMs * maxSampleRate / 1000);

  limiter = (TDLimiterPtr)FDKcalloc(1, sizeof(struct TDLimiter));
  if (!limiter) return NULL;

  limiter->maxBuf = (FIXP_DBL*)FDKcalloc(attack + 1, sizeof(FIXP_DBL));
  limiter->delayBuf =
      (FIXP_DBL*)FDKcalloc(attack * maxChannels, sizeof(FIXP_DBL));

  if (!limiter->maxBuf || !limiter->delayBuf) {
    pcmLimiter_Destroy(limiter);
    return NULL;
  }

  exponent = invFixp(attack + 1);
  attackConst = fPow(FL2FXCONST_DBL(0.1f), 0, exponent, 0, &e_ans);
  attackConst = scaleValue(attackConst, e_ans);

  exponent = invFixp(release + 1);
  releaseConst = fPow(FL2FXCONST_DBL(0.1f), 0, exponent, 0, &e_ans);
  releaseConst = scaleValue(releaseConst, e_ans);

  limiter->attackMs = maxAttackMs;
  limiter->maxAttackMs = maxAttackMs;
  limiter->releaseMs = releaseMs;
  limiter->attack = attack;
  limiter->attackConst = attackConst;
  limiter->releaseConst = releaseConst;
  limiter->threshold = threshold;
  limiter->channels = maxChannels;
  limiter->maxChannels = maxChannels;
  limiter->sampleRate = maxSampleRate;
  limiter->maxSampleRate = maxSampleRate;

  pcmLimiter_Reset(limiter);

  return limiter;
}

TDLIMITER_ERROR pcmLimiter_Apply(TDLimiterPtr limiter, PCM_LIM* samplesIn,
                                 INT_PCM* samplesOut, FIXP_DBL* pGainPerSample,
                                 const INT scaling, const UINT nSamples) {
  unsigned int i, j;
  FIXP_DBL tmp2;
  FIXP_DBL tmp, old, gain, additionalGain = 0;
  FIXP_DBL minGain = FL2FXCONST_DBL(1.0f / (1 << 1));
  UINT additionalGainAvailable = 1;

  if (limiter == NULL) return TDLIMIT_INVALID_HANDLE;

  {
    unsigned int channels = limiter->channels;
    unsigned int attack = limiter->attack;
    FIXP_DBL attackConst = limiter->attackConst;
    FIXP_DBL releaseConst = limiter->releaseConst;
    FIXP_DBL threshold = limiter->threshold >> scaling;

    FIXP_DBL max = limiter->max;
    FIXP_DBL* maxBuf = limiter->maxBuf;
    unsigned int maxBufIdx = limiter->maxBufIdx;
    FIXP_DBL cor = limiter->cor;
    FIXP_DBL* delayBuf = limiter->delayBuf;
    unsigned int delayBufIdx = limiter->delayBufIdx;

    FIXP_DBL smoothState0 = limiter->smoothState0;

    if (limiter->scaling != scaling) {
      scaleValuesSaturate(delayBuf, attack * channels,
                          limiter->scaling - scaling);
      scaleValuesSaturate(maxBuf, attack + 1, limiter->scaling - scaling);
      max = scaleValueSaturate(max, limiter->scaling - scaling);
      limiter->scaling = scaling;
    }

    if (pGainPerSample == NULL) {
      additionalGainAvailable = 0;
    }

    for (i = 0; i < nSamples; i++) {
      tmp = (FIXP_DBL)0;
      for (j = 0; j < channels; j++) {
        tmp2 = PCM_LIM2FIXP_DBL(samplesIn[j]);
        tmp2 =
            (tmp2 == (FIXP_DBL)MINVAL_DBL) ? (FIXP_DBL)MAXVAL_DBL : fAbs(tmp2);
        tmp = fMax(tmp, tmp2);
      }

      if (additionalGainAvailable) {
        additionalGain = pGainPerSample[i];
        tmp = fMult(tmp, additionalGain);
      }

      tmp = fMax(tmp, threshold);

      old = maxBuf[maxBufIdx];
      maxBuf[maxBufIdx] = tmp;

      if (tmp >= max) {
        max = tmp;
      } else if (old < max) {
      } else {
        max = maxBuf[0];
        for (j = 1; j <= attack; j++) {
          max = fMax(max, maxBuf[j]);
        }
      }
      maxBufIdx++;
      if (maxBufIdx >= attack + 1) maxBufIdx = 0;

      if (max > threshold) {
        gain = fDivNorm(threshold, max) >> 1;
      } else {
        gain = FL2FXCONST_DBL(1.0f / (1 << 1));
      }

      if (gain < smoothState0) {
        cor = fMin(cor,
                   fMultDiv2((gain - fMultDiv2(FL2FXCONST_SGL(0.1f * (1 << 1)),
                                               smoothState0)),
                             FL2FXCONST_SGL(1.11111111f / (1 << 1)))
                       << 2);
      } else {
        cor = gain;
      }

      if (cor < smoothState0) {
        smoothState0 =
            fMult(attackConst, (smoothState0 - cor)) + cor;
        smoothState0 = fMax(smoothState0, gain);
      } else {
        smoothState0 =
            -fMult(releaseConst, -(smoothState0 - cor)) + cor;
      }

      gain = smoothState0;

      FIXP_DBL* p_delayBuf = &delayBuf[delayBufIdx * channels + 0];
      if (gain < FL2FXCONST_DBL(1.0f / (1 << 1))) {
        gain <<= 1;
        for (j = 0; j < channels; j++) {
          tmp = p_delayBuf[j];

          if (additionalGainAvailable) {
            p_delayBuf[j] = fMult((FIXP_PCM_LIM)samplesIn[j], additionalGain);
          } else {
            p_delayBuf[j] = PCM_LIM2FIXP_DBL(samplesIn[j]);
          }

          tmp = fMultDiv2(tmp, gain);
#if (SAMPLE_BITS == DFRACT_BITS)
          samplesOut[j] = (INT_PCM)FX_DBL2FX_PCM(
              (FIXP_DBL)SATURATE_LEFT_SHIFT(tmp, scaling + 1, DFRACT_BITS));
#else
          samplesOut[j] = (INT_PCM)FX_DBL2FX_PCM((FIXP_DBL)SATURATE_LEFT_SHIFT(
              tmp + ((FIXP_DBL)0x8000 >> (scaling + 1)), scaling + 1,
              DFRACT_BITS));
#endif
        }
        gain >>= 1;
      } else {
        for (j = 0; j < channels; j++) {
          tmp = p_delayBuf[j];
          if (additionalGainAvailable) {
            p_delayBuf[j] = fMult((FIXP_PCM_LIM)samplesIn[j], additionalGain);
          } else {
            p_delayBuf[j] = PCM_LIM2FIXP_DBL(samplesIn[j]);
          }

#if (SAMPLE_BITS == DFRACT_BITS)
          samplesOut[j] = (INT_PCM)FX_DBL2FX_PCM(
              (FIXP_DBL)SATURATE_LEFT_SHIFT(tmp, scaling, DFRACT_BITS));
#else
          samplesOut[j] = (INT_PCM)FX_DBL2FX_PCM((FIXP_DBL)SATURATE_LEFT_SHIFT(
              (tmp >> 1) + ((FIXP_DBL)0x8000 >> (scaling + 1)), scaling + 1,
              DFRACT_BITS));
#endif
        }
      }

      delayBufIdx++;
      if (delayBufIdx >= attack) {
        delayBufIdx = 0;
      }

      if (gain < minGain) {
        minGain = gain;
      }

      samplesIn += channels;
      samplesOut += channels;
    }

    limiter->max = max;
    limiter->maxBufIdx = maxBufIdx;
    limiter->cor = cor;
    limiter->delayBufIdx = delayBufIdx;

    limiter->smoothState0 = smoothState0;

    limiter->minGain = minGain;

    return TDLIMIT_OK;
  }
}

TDLIMITER_ERROR pcmLimiter_SetThreshold(TDLimiterPtr limiter,
                                        FIXP_DBL threshold) {
  if (limiter == NULL) return TDLIMIT_INVALID_HANDLE;

  limiter->threshold = threshold;

  return TDLIMIT_OK;
}

TDLIMITER_ERROR pcmLimiter_Reset(TDLimiterPtr limiter) {
  if (limiter != NULL) {
    limiter->maxBufIdx = 0;
    limiter->delayBufIdx = 0;
    limiter->max = (FIXP_DBL)0;
    limiter->cor = FL2FXCONST_DBL(1.0f / (1 << 1));
    limiter->smoothState0 = FL2FXCONST_DBL(1.0f / (1 << 1));
    limiter->minGain = FL2FXCONST_DBL(1.0f / (1 << 1));
    limiter->scaling = 0;

    FDKmemset(limiter->maxBuf, 0, (limiter->attack + 1) * sizeof(FIXP_DBL));
    FDKmemset(limiter->delayBuf, 0,
              limiter->attack * limiter->channels * sizeof(FIXP_DBL));
  } else {
    return TDLIMIT_INVALID_HANDLE;
  }

  return TDLIMIT_OK;
}

TDLIMITER_ERROR pcmLimiter_Destroy(TDLimiterPtr limiter) {
  if (limiter != NULL) {
    FDKfree(limiter->maxBuf);
    FDKfree(limiter->delayBuf);

    FDKfree(limiter);
  } else {
    return TDLIMIT_INVALID_HANDLE;
  }
  return TDLIMIT_OK;
}

unsigned int pcmLimiter_GetDelay(TDLimiterPtr limiter) {
  FDK_ASSERT(limiter != NULL);
  return limiter->attack;
}

INT pcmLimiter_GetMaxGainReduction(TDLimiterPtr limiter) {
  int e_ans;
  FIXP_DBL loggain, maxGainReduction;

  FDK_ASSERT(limiter != NULL);

  loggain = fLog2(limiter->minGain, 1, &e_ans);

  maxGainReduction = fMult(loggain, FL2FXCONST_DBL(-6.0206f / (1 << 3)));

  return fixp_roundToInt(maxGainReduction, (e_ans + 3));
}

TDLIMITER_ERROR pcmLimiter_SetNChannels(TDLimiterPtr limiter,
                                        unsigned int nChannels) {
  if (limiter == NULL) return TDLIMIT_INVALID_HANDLE;

  if (nChannels > limiter->maxChannels) return TDLIMIT_INVALID_PARAMETER;

  limiter->channels = nChannels;

  return TDLIMIT_OK;
}

TDLIMITER_ERROR pcmLimiter_SetSampleRate(TDLimiterPtr limiter,
                                         UINT sampleRate) {
  unsigned int attack, release;
  FIXP_DBL attackConst, releaseConst, exponent;
  INT e_ans;

  if (limiter == NULL) return TDLIMIT_INVALID_HANDLE;

  if (sampleRate > limiter->maxSampleRate) return TDLIMIT_INVALID_PARAMETER;

  attack = (unsigned int)(limiter->attackMs * sampleRate / 1000);
  release = (unsigned int)(limiter->releaseMs * sampleRate / 1000);

  exponent = invFixp(attack + 1);
  attackConst = fPow(FL2FXCONST_DBL(0.1f), 0, exponent, 0, &e_ans);
  attackConst = scaleValue(attackConst, e_ans);

  exponent = invFixp(release + 1);
  releaseConst = fPow(FL2FXCONST_DBL(0.1f), 0, exponent, 0, &e_ans);
  releaseConst = scaleValue(releaseConst, e_ans);

  limiter->attack = attack;
  limiter->attackConst = attackConst;
  limiter->releaseConst = releaseConst;
  limiter->sampleRate = sampleRate;

  return TDLIMIT_OK;
}

TDLIMITER_ERROR pcmLimiter_SetAttack(TDLimiterPtr limiter,
                                     unsigned int attackMs) {
  unsigned int attack;
  FIXP_DBL attackConst, exponent;
  INT e_ans;

  if (limiter == NULL) return TDLIMIT_INVALID_HANDLE;

  if (attackMs > limiter->maxAttackMs) return TDLIMIT_INVALID_PARAMETER;

  attack = (unsigned int)(attackMs * limiter->sampleRate / 1000);

  exponent = invFixp(attack + 1);
  attackConst = fPow(FL2FXCONST_DBL(0.1f), 0, exponent, 0, &e_ans);
  attackConst = scaleValue(attackConst, e_ans);

  limiter->attack = attack;
  limiter->attackConst = attackConst;
  limiter->attackMs = attackMs;

  return TDLIMIT_OK;
}

TDLIMITER_ERROR pcmLimiter_SetRelease(TDLimiterPtr limiter,
                                      unsigned int releaseMs) {
  unsigned int release;
  FIXP_DBL releaseConst, exponent;
  INT e_ans;

  if (limiter == NULL) return TDLIMIT_INVALID_HANDLE;

  release = (unsigned int)(releaseMs * limiter->sampleRate / 1000);

  exponent = invFixp(release + 1);
  releaseConst = fPow(FL2FXCONST_DBL(0.1f), 0, exponent, 0, &e_ans);
  releaseConst = scaleValue(releaseConst, e_ans);

  limiter->releaseConst = releaseConst;
  limiter->releaseMs = releaseMs;

  return TDLIMIT_OK;
}

TDLIMITER_ERROR pcmLimiter_GetLibInfo(LIB_INFO* info) {
  int i;

  if (info == NULL) {
    return TDLIMIT_INVALID_PARAMETER;
  }

  for (i = 0; i < FDK_MODULE_LAST; i++) {
    if (info[i].module_id == FDK_NONE) break;
  }
  if (i == FDK_MODULE_LAST) {
    return TDLIMIT_UNKNOWN;
  }

  info[i].module_id = FDK_TDLIMIT;
  info[i].version =
      LIB_VERSION(PCMUTIL_LIB_VL0, PCMUTIL_LIB_VL1, PCMUTIL_LIB_VL2);
  LIB_VERSION_STRING(info + i);
  info[i].build_date = PCMUTIL_LIB_BUILD_DATE;
  info[i].build_time = PCMUTIL_LIB_BUILD_TIME;
  info[i].title = TDLIMIT_LIB_TITLE;

  info[i].flags = CAPF_LIMITER;

  FDK_toolsGetLibInfo(info);

  return TDLIMIT_OK;
}
