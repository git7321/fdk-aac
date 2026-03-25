/************************* MPEG-D DRC decoder library **************************

   Author(s):

   Description:

*******************************************************************************/

#include "drcDec_types.h"
#include "drcDec_tools.h"
#include "fixpoint_math.h"
#include "drcDecoder.h"

int getDeltaTmin(const int sampleRate) {
  int half_ms = (sampleRate + 1000) / 2000;
  int deltaTmin = 1;
  if (sampleRate < 1000) {
    return DE_NOT_OK;
  }
  while (deltaTmin <= half_ms) {
    deltaTmin = deltaTmin << 1;
  }
  return deltaTmin;
}

DRC_COEFFICIENTS_UNI_DRC* selectDrcCoefficients(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, const int location) {
  int n;
  int c = -1;
  for (n = 0; n < hUniDrcConfig->drcCoefficientsUniDrcCount; n++) {
    if (hUniDrcConfig->drcCoefficientsUniDrc[n].drcLocation == location) {
      c = n;
    }
  }
  if (c >= 0) {
    return &(hUniDrcConfig->drcCoefficientsUniDrc[c]);
  }
  return NULL;
}

DRC_INSTRUCTIONS_UNI_DRC* selectDrcInstructions(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, const int drcSetId) {
  int i;
  for (i = 0; i < hUniDrcConfig->drcInstructionsCountInclVirtual; i++) {
    if (hUniDrcConfig->drcInstructionsUniDrc[i].drcSetId == drcSetId) {
      return &(hUniDrcConfig->drcInstructionsUniDrc[i]);
    }
  }
  return NULL;
}

DOWNMIX_INSTRUCTIONS* selectDownmixInstructions(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, const int downmixId) {
  int i;
  for (i = 0; i < hUniDrcConfig->downmixInstructionsCount; i++) {
    if (hUniDrcConfig->downmixInstructions[i].downmixId == downmixId) {
      return &(hUniDrcConfig->downmixInstructions[i]);
    }
  }
  return NULL;
}

DRC_ERROR
deriveDrcChannelGroups(
    const int drcSetEffect,
    const int channelCount,
    const SCHAR* gainSetIndex,
    const DUCKING_MODIFICATION* duckingModificationForChannel,
    UCHAR* nDrcChannelGroups,
    SCHAR* uniqueIndex,
    SCHAR* groupForChannel,
    DUCKING_MODIFICATION* duckingModificationForChannelGroup)
{
  int duckingSequence = -1;
  int c, n, g, match, idx;
  FIXP_SGL factor;
  FIXP_SGL uniqueScaling[8];

  for (g = 0; g < 8; g++) {
    uniqueIndex[g] = -10;
    uniqueScaling[g] = FIXP_SGL(-1.0f);
  }

  g = 0;

  if (drcSetEffect & EB_DUCK_OTHER) {
    for (c = 0; c < channelCount; c++) {
      match = 0;
      if (c >= 8) return DE_MEMORY_ERROR;
      idx = gainSetIndex[c];
      factor = duckingModificationForChannel[c].duckingScaling;
      if (idx < 0) {
        for (n = 0; n < g; n++) {
          if (uniqueScaling[n] == factor) {
            match = 1;
            groupForChannel[c] = n;
            break;
          }
        }
        if (match == 0) {
          if (g >= 8) return DE_MEMORY_ERROR;
          uniqueIndex[g] = idx;
          uniqueScaling[g] = factor;
          groupForChannel[c] = g;
          g++;
        }
      } else {
        if ((duckingSequence > 0) && (duckingSequence != idx)) {
          return DE_NOT_OK;
        }
        duckingSequence = idx;
        groupForChannel[c] = -1;
      }
    }
    if (duckingSequence == -1) {
      return DE_NOT_OK;
    }
  } else if (drcSetEffect & EB_DUCK_SELF) {
    for (c = 0; c < channelCount; c++) {
      match = 0;
      if (c >= 8) return DE_MEMORY_ERROR;
      idx = gainSetIndex[c];
      factor = duckingModificationForChannel[c].duckingScaling;
      if (idx >= 0) {
        for (n = 0; n < g; n++) {
          if ((uniqueIndex[n] == idx) && (uniqueScaling[n] == factor)) {
            match = 1;
            groupForChannel[c] = n;
            break;
          }
        }
        if (match == 0) {
          if (g >= 8) return DE_MEMORY_ERROR;
          uniqueIndex[g] = idx;
          uniqueScaling[g] = factor;
          groupForChannel[c] = g;
          g++;
        }
      } else {
        groupForChannel[c] = -1;
      }
    }
  } else {
    for (c = 0; c < channelCount; c++) {
      if (c >= 8) return DE_MEMORY_ERROR;
      idx = gainSetIndex[c];
      match = 0;
      if (idx >= 0) {
        for (n = 0; n < g; n++) {
          if (uniqueIndex[n] == idx) {
            match = 1;
            groupForChannel[c] = n;
            break;
          }
        }
        if (match == 0) {
          if (g >= 8) return DE_MEMORY_ERROR;
          uniqueIndex[g] = idx;
          groupForChannel[c] = g;
          g++;
        }
      } else {
        groupForChannel[c] = -1;
      }
    }
  }
  *nDrcChannelGroups = g;

  if (drcSetEffect & (EB_DUCK_OTHER | EB_DUCK_SELF)) {
    for (g = 0; g < *nDrcChannelGroups; g++) {
      if (drcSetEffect & EB_DUCK_OTHER) {
        uniqueIndex[g] = duckingSequence;
      }
      duckingModificationForChannelGroup[g].duckingScaling = uniqueScaling[g];
      if (uniqueScaling[g] != FL2FXCONST_SGL(1.0f / (float)(1 << 2))) {
        duckingModificationForChannelGroup[g].duckingScalingPresent = 1;
      } else {
        duckingModificationForChannelGroup[g].duckingScalingPresent = 0;
      }
    }
  }

  return DE_OK;
}

FIXP_DBL
dB2lin(const FIXP_DBL dB_m, const int dB_e, int* pLin_e) {
  FIXP_DBL lin_m =
      f2Pow(fMult(dB_m, FL2FXCONST_DBL(0.1660964f * (float)(1 << 2))), dB_e - 2,
            pLin_e);

  return lin_m;
}

FIXP_DBL
lin2dB(const FIXP_DBL lin_m, const int lin_e, int* pDb_e) {
  FIXP_DBL dB_m;

  if (lin_m == (FIXP_DBL)0) {
    dB_m = (FIXP_DBL)MINVAL_DBL;
    *pDb_e = DFRACT_BITS - 1;
  } else {
    dB_m = fMultDiv2(FL2FXCONST_DBL(6.02059991f / (float)(1 << 3)),
                     fLog2(lin_m, lin_e, pDb_e));
    *pDb_e += 3 + 1;
  }

  return dB_m;
}

FIXP_DBL
approxDb2lin(const FIXP_DBL dB_m, const int dB_e, int* pLin_e) {
  FIXP_DBL lin_m =
      f2Pow(fMult(dB_m, FL2FXCONST_DBL(0.1666667f * (float)(1 << 2))), dB_e - 2,
            pLin_e);

  return lin_m;
}

int bitstreamContainsMultibandDrc(HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
                                  const int downmixId) {
  int i, g, d, seq;
  DRC_INSTRUCTIONS_UNI_DRC* pInst;
  DRC_COEFFICIENTS_UNI_DRC* pCoef = NULL;
  int isMultiband = 0;

  pCoef = selectDrcCoefficients(hUniDrcConfig, LOCATION_SELECTED);
  if (pCoef == NULL) return 0;

  for (i = 0; i < hUniDrcConfig->drcInstructionsUniDrcCount; i++) {
    pInst = &(hUniDrcConfig->drcInstructionsUniDrc[i]);
    for (d = 0; d < pInst->downmixIdCount; d++) {
      if (downmixId == pInst->downmixId[d]) {
        for (g = 0; g < pInst->nDrcChannelGroups; g++) {
          seq = pInst->gainSetIndexForChannelGroup[g];
          if (pCoef->gainSet[seq].bandCount > 1) {
            isMultiband = 1;
          }
        }
      }
    }
  }

  return isMultiband;
}

FIXP_DBL getDownmixOffset(DOWNMIX_INSTRUCTIONS* pDown, int baseChannelCount) {
  FIXP_DBL downmixOffset = FL2FXCONST_DBL(1.0f / (1 << 1));
  if ((pDown->bsDownmixOffset == 1) || (pDown->bsDownmixOffset == 2)) {
    int e_a, e_downmixOffset;
    FIXP_DBL a, q;
    if (baseChannelCount <= pDown->targetChannelCount) return downmixOffset;

    q = fDivNorm((FIXP_DBL)pDown->targetChannelCount,
                 (FIXP_DBL)baseChannelCount);
    a = lin2dB(q, 0, &e_a);
    if (pDown->bsDownmixOffset == 2) {
      e_a += 1;
    }
    a = fixp_round(a, e_a) >> 1;
    downmixOffset = dB2lin(a, e_a, &e_downmixOffset);
    downmixOffset = scaleValue(downmixOffset, e_downmixOffset - 1);
  }
  return downmixOffset;
}
