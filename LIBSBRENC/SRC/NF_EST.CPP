/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "nf_est.h"
#include "sbr_misc.h"
#include "genericStds.h"

static const FIXP_DBL smoothFilter[4] = {0x077f813d, 0x19999995, 0x2bb3b1f5,
                                         0x33333335};

static const FIXP_DBL QuantOffset = (INT)0xfc000000;

#ifndef min
#define min(a, b) (a < b ? a : b)
#endif

#ifndef max
#define max(a, b) (a > b ? a : b)
#endif

#define NOISE_FLOOR_OFFSET_SCALING (4)

static void smoothingOfNoiseLevels(
    FIXP_DBL *NoiseLevels,
    INT nEnvelopes,
    INT noNoiseBands,
    FIXP_DBL prevNoiseLevels[NF_SMOOTHING_LENGTH]
                            [MAX_NUM_NOISE_VALUES],
    const FIXP_DBL *
        pSmoothFilter,
    INT transientFlag)
{
  INT i, band, env;
  FIXP_DBL accu;

  for (env = 0; env < nEnvelopes; env++) {
    if (transientFlag) {
      for (i = 0; i < NF_SMOOTHING_LENGTH; i++) {
        FDKmemcpy(prevNoiseLevels[i], NoiseLevels + env * noNoiseBands,
                  noNoiseBands * sizeof(FIXP_DBL));
      }
    } else {
      for (i = 1; i < NF_SMOOTHING_LENGTH; i++) {
        FDKmemcpy(prevNoiseLevels[i - 1], prevNoiseLevels[i],
                  noNoiseBands * sizeof(FIXP_DBL));
      }
      FDKmemcpy(prevNoiseLevels[NF_SMOOTHING_LENGTH - 1],
                NoiseLevels + env * noNoiseBands,
                noNoiseBands * sizeof(FIXP_DBL));
    }

    for (band = 0; band < noNoiseBands; band++) {
      accu = FL2FXCONST_DBL(0.0f);
      for (i = 0; i < NF_SMOOTHING_LENGTH; i++) {
        accu += fMultDiv2(pSmoothFilter[i], prevNoiseLevels[i][band]);
      }
      FDK_ASSERT((band + env * noNoiseBands) < MAX_NUM_NOISE_VALUES);
      NoiseLevels[band + env * noNoiseBands] = accu << 1;
    }
  }
}

static void qmfBasedNoiseFloorDetection(
    FIXP_DBL *noiseLevel,
    FIXP_DBL **quotaMatrixOrig,
    SCHAR *indexVector,
    INT startIndex,
    INT stopIndex,
    INT startChannel,
    INT stopChannel,
    FIXP_DBL ana_max_level,
    FIXP_DBL noiseFloorOffset,
    INT missingHarmonicFlag,
    FIXP_DBL weightFac,
    INVF_MODE diffThres,
    INVF_MODE inverseFilteringLevel)
{
  INT scale, l, k;
  FIXP_DBL meanOrig = FL2FXCONST_DBL(0.0f), meanSbr = FL2FXCONST_DBL(0.0f),
           diff;
  FIXP_DBL invIndex = GetInvInt(stopIndex - startIndex);
  FIXP_DBL invChannel = GetInvInt(stopChannel - startChannel);
  FIXP_DBL accu;

  if (missingHarmonicFlag == 1) {
    for (l = startChannel; l < stopChannel; l++) {
      accu = FL2FXCONST_DBL(0.0f);
      for (k = startIndex; k < stopIndex; k++) {
        accu += fMultDiv2(quotaMatrixOrig[k][l], invIndex);
      }
      meanOrig = fixMax(meanOrig, (accu << 1));

      accu = FL2FXCONST_DBL(0.0f);
      for (k = startIndex; k < stopIndex; k++) {
        accu += fMultDiv2(quotaMatrixOrig[k][indexVector[l]], invIndex);
      }
      meanSbr = fixMax(meanSbr, (accu << 1));
    }
  } else {
    for (l = startChannel; l < stopChannel; l++) {
      accu = FL2FXCONST_DBL(0.0f);
      for (k = startIndex; k < stopIndex; k++) {
        accu += fMultDiv2(quotaMatrixOrig[k][l], invIndex);
      }
      meanOrig += fMult((accu << 1), invChannel);

      accu = FL2FXCONST_DBL(0.0f);
      for (k = startIndex; k < stopIndex; k++) {
        accu += fMultDiv2(quotaMatrixOrig[k][indexVector[l]], invIndex);
      }
      meanSbr += fMult((accu << 1), invChannel);
    }
  }

  if (meanOrig <= FL2FXCONST_DBL(0.000976562f * RELAXATION_FLOAT) &&
      meanSbr <= FL2FXCONST_DBL(0.000976562f * RELAXATION_FLOAT)) {
    meanOrig = FL2FXCONST_DBL(101.5936673f * RELAXATION_FLOAT);
    meanSbr = FL2FXCONST_DBL(101.5936673f * RELAXATION_FLOAT);
  }

  meanOrig = fixMax(meanOrig, RELAXATION);
  meanSbr = fixMax(meanSbr, RELAXATION);

  if (missingHarmonicFlag == 1 || inverseFilteringLevel == INVF_MID_LEVEL ||
      inverseFilteringLevel == INVF_LOW_LEVEL ||
      inverseFilteringLevel == INVF_OFF || inverseFilteringLevel <= diffThres) {
    diff = RELAXATION;
  } else {
    accu = fDivNorm(meanSbr, meanOrig, &scale);

    diff = fixMax(RELAXATION, fMult(RELAXATION_FRACT, fMult(weightFac, accu)) >>
                                  (RELAXATION_SHIFT - scale));
  }

  accu = fDivNorm(diff, meanOrig, &scale);
  scale -= 2;

  if ((scale > 0) && (accu > ((FIXP_DBL)MAXVAL_DBL) >> scale)) {
    *noiseLevel = (FIXP_DBL)MAXVAL_DBL;
  } else {
    *noiseLevel = scaleValue(accu, scale);
  }

  if (!missingHarmonicFlag) {
    *noiseLevel = fixMin(fMult(*noiseLevel, noiseFloorOffset),
                         (FIXP_DBL)MAXVAL_DBL >> NOISE_FLOOR_OFFSET_SCALING)
                  << NOISE_FLOOR_OFFSET_SCALING;
  }

  *noiseLevel =
      fixMin(*noiseLevel,
             ana_max_level);
}

void FDKsbrEnc_sbrNoiseFloorEstimateQmf(
    HANDLE_SBR_NOISE_FLOOR_ESTIMATE
        h_sbrNoiseFloorEstimate,
    const SBR_FRAME_INFO
        *frame_info,
    FIXP_DBL
        *noiseLevels,
    FIXP_DBL **quotaMatrixOrig,
    SCHAR *indexVector,
    INT missingHarmonicsFlag,
    INT startIndex,
    UINT numberOfEstimatesPerFrame,
    int transientFrame,
    INVF_MODE *pInvFiltLevels,
    UINT sbrSyntaxFlags)

{
  INT nNoiseEnvelopes, startPos[2], stopPos[2], env, band;

  INT noNoiseBands = h_sbrNoiseFloorEstimate->noNoiseBands;
  INT *freqBandTable = h_sbrNoiseFloorEstimate->freqBandTableQmf;

  nNoiseEnvelopes = frame_info->nNoiseEnvelopes;

  startPos[0] = startIndex;

  if (nNoiseEnvelopes == 1) {
    stopPos[0] = startIndex + min(numberOfEstimatesPerFrame, 2);
  } else {
    stopPos[0] = startIndex + 1;
    startPos[1] = startIndex + 1;
    stopPos[1] = startIndex + min(numberOfEstimatesPerFrame, 2);
  }

  for (env = 0; env < nNoiseEnvelopes; env++) {
    for (band = 0; band < noNoiseBands; band++) {
      FDK_ASSERT((band + env * noNoiseBands) < MAX_NUM_NOISE_VALUES);
      qmfBasedNoiseFloorDetection(
          &noiseLevels[band + env * noNoiseBands], quotaMatrixOrig, indexVector,
          startPos[env], stopPos[env], freqBandTable[band],
          freqBandTable[band + 1], h_sbrNoiseFloorEstimate->ana_max_level,
          h_sbrNoiseFloorEstimate->noiseFloorOffset[band], missingHarmonicsFlag,
          h_sbrNoiseFloorEstimate->weightFac,
          h_sbrNoiseFloorEstimate->diffThres, pInvFiltLevels[band]);
    }
  }

  smoothingOfNoiseLevels(noiseLevels, nNoiseEnvelopes,
                         h_sbrNoiseFloorEstimate->noNoiseBands,
                         h_sbrNoiseFloorEstimate->prevNoiseLevels,
                         h_sbrNoiseFloorEstimate->smoothFilter, transientFrame);

  for (env = 0; env < nNoiseEnvelopes; env++) {
    for (band = 0; band < noNoiseBands; band++) {
      FDK_ASSERT((band + env * noNoiseBands) < MAX_NUM_NOISE_VALUES);
      noiseLevels[band + env * noNoiseBands] =
          (FIXP_DBL)NOISE_FLOOR_OFFSET_64 -
          (FIXP_DBL)CalcLdData(noiseLevels[band + env * noNoiseBands] +
                               (FIXP_DBL)1) +
          QuantOffset;
    }
  }
}

static INT downSampleLoRes(INT *v_result,
                           INT num_result,
                           const UCHAR *freqBandTableRef,
                           INT num_Ref)
{
  INT step;
  INT i, j;
  INT org_length, result_length;
  INT v_index[MAX_FREQ_COEFFS / 2];

  org_length = num_Ref;
  result_length = num_result;

  v_index[0] = 0;
  i = 0;
  while (org_length > 0)
  {
    i++;
    step = org_length / result_length;
    org_length = org_length - step;
    result_length--;
    v_index[i] = v_index[i - 1] + step;
  }

  if (i != num_result)
    return (1);

  for (j = 0; j <= i;
       j++)
  {
    v_result[j] = freqBandTableRef[v_index[j]];
  }

  return (0);
}

INT FDKsbrEnc_InitSbrNoiseFloorEstimate(
    HANDLE_SBR_NOISE_FLOOR_ESTIMATE
        h_sbrNoiseFloorEstimate,
    INT ana_max_level,
    const UCHAR *freqBandTable,
    INT nSfb,
    INT noiseBands,
    INT noiseFloorOffset,
    INT timeSlots,
    UINT useSpeechConfig
) {
  INT i, qexp, qtmp;
  FIXP_DBL tmp, exp;

  FDKmemclear(h_sbrNoiseFloorEstimate, sizeof(SBR_NOISE_FLOOR_ESTIMATE));

  h_sbrNoiseFloorEstimate->smoothFilter = smoothFilter;
  if (useSpeechConfig) {
    h_sbrNoiseFloorEstimate->weightFac = (FIXP_DBL)MAXVAL_DBL;
    h_sbrNoiseFloorEstimate->diffThres = INVF_LOW_LEVEL;
  } else {
    h_sbrNoiseFloorEstimate->weightFac = FL2FXCONST_DBL(0.25f);
    h_sbrNoiseFloorEstimate->diffThres = INVF_MID_LEVEL;
  }

  h_sbrNoiseFloorEstimate->timeSlots = timeSlots;
  h_sbrNoiseFloorEstimate->noiseBands = noiseBands;

  switch (ana_max_level) {
    case 6:
      h_sbrNoiseFloorEstimate->ana_max_level = (FIXP_DBL)MAXVAL_DBL;
      break;
    case 3:
      h_sbrNoiseFloorEstimate->ana_max_level = FL2FXCONST_DBL(0.5);
      break;
    case -3:
      h_sbrNoiseFloorEstimate->ana_max_level = FL2FXCONST_DBL(0.125);
      break;
    default:
      h_sbrNoiseFloorEstimate->ana_max_level = (FIXP_DBL)MAXVAL_DBL;
      break;
  }

  if (FDKsbrEnc_resetSbrNoiseFloorEstimate(h_sbrNoiseFloorEstimate,
                                           freqBandTable, nSfb))
    return (1);

  if (noiseFloorOffset == 0) {
    tmp = ((FIXP_DBL)MAXVAL_DBL) >> NOISE_FLOOR_OFFSET_SCALING;
  } else {
    FDK_ASSERT(noiseFloorOffset < 12);

    exp = fDivNorm((FIXP_DBL)noiseFloorOffset, 3, &qexp);
    tmp = fPow(2, DFRACT_BITS - 1, exp, qexp, &qtmp);
    tmp = scaleValue(tmp, qtmp - NOISE_FLOOR_OFFSET_SCALING);
  }

  for (i = 0; i < h_sbrNoiseFloorEstimate->noNoiseBands; i++) {
    h_sbrNoiseFloorEstimate->noiseFloorOffset[i] = tmp;
  }

  return (0);
}

INT FDKsbrEnc_resetSbrNoiseFloorEstimate(
    HANDLE_SBR_NOISE_FLOOR_ESTIMATE
        h_sbrNoiseFloorEstimate,
    const UCHAR *freqBandTable,
    INT nSfb
) {
  INT k2, kx;

  k2 = freqBandTable[nSfb];
  kx = freqBandTable[0];
  if (h_sbrNoiseFloorEstimate->noiseBands == 0) {
    h_sbrNoiseFloorEstimate->noNoiseBands = 1;
  } else {
    FIXP_DBL tmp, ratio, lg2;
    INT ratio_e, qlg2, nNoiseBands;

    ratio = fDivNorm(k2, kx, &ratio_e);
    lg2 = fLog2(ratio, ratio_e, &qlg2);
    tmp = fMult((FIXP_DBL)(h_sbrNoiseFloorEstimate->noiseBands << 24), lg2);
    tmp = scaleValue(tmp, qlg2 - 23);

    nNoiseBands = (INT)((tmp + (FIXP_DBL)1) >> 1);

    if (nNoiseBands > MAX_NUM_NOISE_COEFFS) {
      nNoiseBands = MAX_NUM_NOISE_COEFFS;
    }

    if (nNoiseBands == 0) {
      nNoiseBands = 1;
    }

    h_sbrNoiseFloorEstimate->noNoiseBands = nNoiseBands;
  }

  return (downSampleLoRes(h_sbrNoiseFloorEstimate->freqBandTableQmf,
                          h_sbrNoiseFloorEstimate->noNoiseBands, freqBandTable,
                          nSfb));
}

void FDKsbrEnc_deleteSbrNoiseFloorEstimate(
    HANDLE_SBR_NOISE_FLOOR_ESTIMATE
        h_sbrNoiseFloorEstimate)
{
  if (h_sbrNoiseFloorEstimate) {
  }
}
