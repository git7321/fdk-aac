/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "env_dec.h"
#include "env_extr.h"
#include "transcendent.h"
#include "genericStds.h"

static void decodeEnvelope(HANDLE_SBR_HEADER_DATA hHeaderData,
                           HANDLE_SBR_FRAME_DATA h_sbr_data,
                           HANDLE_SBR_PREV_FRAME_DATA h_prev_data,
                           HANDLE_SBR_PREV_FRAME_DATA h_prev_data_otherChannel);
static void sbr_envelope_unmapping(HANDLE_SBR_HEADER_DATA hHeaderData,
                                   HANDLE_SBR_FRAME_DATA h_data_left,
                                   HANDLE_SBR_FRAME_DATA h_data_right);
static void requantizeEnvelopeData(HANDLE_SBR_FRAME_DATA h_sbr_data,
                                   int ampResolution);
static void deltaToLinearPcmEnvelopeDecoding(
    HANDLE_SBR_HEADER_DATA hHeaderData, HANDLE_SBR_FRAME_DATA h_sbr_data,
    HANDLE_SBR_PREV_FRAME_DATA h_prev_data);
static void decodeNoiseFloorlevels(HANDLE_SBR_HEADER_DATA hHeaderData,
                                   HANDLE_SBR_FRAME_DATA h_sbr_data,
                                   HANDLE_SBR_PREV_FRAME_DATA h_prev_data);
static void timeCompensateFirstEnvelope(HANDLE_SBR_HEADER_DATA hHeaderData,
                                        HANDLE_SBR_FRAME_DATA h_sbr_data,
                                        HANDLE_SBR_PREV_FRAME_DATA h_prev_data);
static int checkEnvelopeData(HANDLE_SBR_HEADER_DATA hHeaderData,
                             HANDLE_SBR_FRAME_DATA h_sbr_data,
                             HANDLE_SBR_PREV_FRAME_DATA h_prev_data);

#define SBR_ENERGY_PAN_OFFSET (12 << ENV_EXP_FRACT)
#define SBR_MAX_ENERGY (35 << ENV_EXP_FRACT)

#define DECAY (1 << ENV_EXP_FRACT)

#if ENV_EXP_FRACT
#define DECAY_COUPLING \
  (1 << (ENV_EXP_FRACT - 1))
#else
#define DECAY_COUPLING \
  1
#endif

static int indexLow2High(int offset,
                         int index,
                         int res)
{
  if (res == 0) {
    if (offset >= 0) {
      if (index < offset)
        return (index);
      else
        return (2 * index - offset);
    } else {
      offset = -offset;
      if (index < offset)
        return (2 * index + index);
      else
        return (2 * index + offset);
    }
  } else
    return (index);
}

static void mapLowResEnergyVal(
    FIXP_SGL currVal,
    FIXP_SGL *prevData,
    int offset,
    int index,
    int res)
{
  if (res == 0) {
    if (offset >= 0) {
      if (index < offset)
        prevData[index] = currVal;
      else {
        prevData[2 * index - offset] = currVal;
        prevData[2 * index + 1 - offset] = currVal;
      }
    } else {
      offset = -offset;
      if (index < offset) {
        prevData[3 * index] = currVal;
        prevData[3 * index + 1] = currVal;
        prevData[3 * index + 2] = currVal;
      } else {
        prevData[2 * index + offset] = currVal;
        prevData[2 * index + 1 + offset] = currVal;
      }
    }
  } else
    prevData[index] = currVal;
}

void decodeSbrData(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA
        h_data_left,
    HANDLE_SBR_PREV_FRAME_DATA
        h_prev_data_left,
    HANDLE_SBR_FRAME_DATA
        h_data_right,
    HANDLE_SBR_PREV_FRAME_DATA
        h_prev_data_right)
{
  FIXP_SGL tempSfbNrgPrev[MAX_FREQ_COEFFS];
  int errLeft;

  FDKmemcpy(tempSfbNrgPrev, h_prev_data_left->sfb_nrg_prev,
            MAX_FREQ_COEFFS * sizeof(FIXP_SGL));

  if (hHeaderData->frameErrorFlag || hHeaderData->bs_info.pvc_mode == 0) {
    decodeEnvelope(hHeaderData, h_data_left, h_prev_data_left,
                   h_prev_data_right);
  } else {
    FDK_ASSERT(h_data_right == NULL);
  }
  decodeNoiseFloorlevels(hHeaderData, h_data_left, h_prev_data_left);

  if (h_data_right != NULL) {
    errLeft = hHeaderData->frameErrorFlag;
    decodeEnvelope(hHeaderData, h_data_right, h_prev_data_right,
                   h_prev_data_left);
    decodeNoiseFloorlevels(hHeaderData, h_data_right, h_prev_data_right);

    if (!errLeft && hHeaderData->frameErrorFlag) {
      FDKmemcpy(h_prev_data_left->sfb_nrg_prev, tempSfbNrgPrev,
                MAX_FREQ_COEFFS * sizeof(FIXP_SGL));
      decodeEnvelope(hHeaderData, h_data_left, h_prev_data_left,
                     h_prev_data_right);
    }

    if (h_data_left->coupling) {
      sbr_envelope_unmapping(hHeaderData, h_data_left, h_data_right);
    }
  }
}

static void sbr_envelope_unmapping(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_data_left,
    HANDLE_SBR_FRAME_DATA h_data_right)
{
  int i;
  FIXP_SGL tempL_m, tempR_m, tempRplus1_m, newL_m, newR_m;
  SCHAR tempL_e, tempR_e, tempRplus1_e, newL_e, newR_e;

  for (i = 0; i < h_data_left->nScaleFactors; i++) {
    tempR_m = (FIXP_SGL)((LONG)h_data_right->iEnvelope[i] & MASK_M);
    tempR_e = (SCHAR)((LONG)h_data_right->iEnvelope[i] & MASK_E);

    tempR_e -= (18 + NRG_EXP_OFFSET);
    tempL_m = (FIXP_SGL)((LONG)h_data_left->iEnvelope[i] & MASK_M);
    tempL_e = (SCHAR)((LONG)h_data_left->iEnvelope[i] & MASK_E);

    tempL_e -= NRG_EXP_OFFSET;

    FDK_add_MantExp(tempR_m, tempR_e, FL2FXCONST_SGL(0.5f), 1,
                    &tempRplus1_m, &tempRplus1_e);

    FDK_divide_MantExp(tempL_m, tempL_e + 1,
                       tempRplus1_m, tempRplus1_e, &newR_m, &newR_e);

    if (newR_m >= ((FIXP_SGL)MAXVAL_SGL - ROUNDING)) {
      newR_m >>= 1;
      newR_e += 1;
    }

    newL_m = FX_DBL2FX_SGL(fMult(tempR_m, newR_m));
    newL_e = tempR_e + newR_e;

    h_data_right->iEnvelope[i] =
        ((FIXP_SGL)((SHORT)(FIXP_SGL)(newR_m + ROUNDING) & MASK_M)) +
        (FIXP_SGL)((SHORT)(FIXP_SGL)(newR_e + NRG_EXP_OFFSET) & MASK_E);
    h_data_left->iEnvelope[i] =
        ((FIXP_SGL)((SHORT)(FIXP_SGL)(newL_m + ROUNDING) & MASK_M)) +
        (FIXP_SGL)((SHORT)(FIXP_SGL)(newL_e + NRG_EXP_OFFSET) & MASK_E);
  }

  for (i = 0; i < hHeaderData->freqBandData.nNfb *
                      h_data_left->frameInfo.nNoiseEnvelopes;
       i++) {
    tempL_e = (SCHAR)(6 - (LONG)h_data_left->sbrNoiseFloorLevel[i]);
    tempR_e = (SCHAR)((LONG)h_data_right->sbrNoiseFloorLevel[i] -
                      12);

    FDK_add_MantExp(FL2FXCONST_SGL(0.5f), 1 + tempR_e,
                    FL2FXCONST_SGL(0.5f), 1,
                    &tempRplus1_m, &tempRplus1_e);

    FDK_divide_MantExp(FL2FXCONST_SGL(0.5f), tempL_e + 2,
                       tempRplus1_m, tempRplus1_e, &newR_m, &newR_e);

    newL_m = newR_m;
    newL_e = newR_e + tempR_e;
    h_data_right->sbrNoiseFloorLevel[i] =
        ((FIXP_SGL)((SHORT)(FIXP_SGL)(newR_m + ROUNDING) & MASK_M)) +
        (FIXP_SGL)((SHORT)(FIXP_SGL)(newR_e + NOISE_EXP_OFFSET) & MASK_E);
    h_data_left->sbrNoiseFloorLevel[i] =
        ((FIXP_SGL)((SHORT)(FIXP_SGL)(newL_m + ROUNDING) & MASK_M)) +
        (FIXP_SGL)((SHORT)(FIXP_SGL)(newL_e + NOISE_EXP_OFFSET) & MASK_E);
  }
}

static void leanSbrConcealment(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_sbr_data,
    HANDLE_SBR_PREV_FRAME_DATA h_prev_data
) {
  FIXP_SGL target;
  FIXP_SGL step;
  int i;

  int currentStartPos =
      fMax(0, h_prev_data->stopPos - hHeaderData->numberTimeSlots);
  int currentStopPos = hHeaderData->numberTimeSlots;

  h_sbr_data->ampResolutionCurrentFrame = h_prev_data->ampRes;
  h_sbr_data->coupling = h_prev_data->coupling;
  for (i = 0; i < MAX_INVF_BANDS; i++)
    h_sbr_data->sbr_invf_mode[i] = h_prev_data->sbr_invf_mode[i];

  h_sbr_data->frameInfo.nEnvelopes = 1;
  h_sbr_data->frameInfo.borders[0] = currentStartPos;
  h_sbr_data->frameInfo.borders[1] = currentStopPos;
  h_sbr_data->frameInfo.freqRes[0] = 1;
  h_sbr_data->frameInfo.tranEnv = -1;
  h_sbr_data->frameInfo.nNoiseEnvelopes = 1;
  h_sbr_data->frameInfo.bordersNoise[0] = currentStartPos;
  h_sbr_data->frameInfo.bordersNoise[1] = currentStopPos;

  h_sbr_data->nScaleFactors = hHeaderData->freqBandData.nSfb[1];

  h_sbr_data->domain_vec[0] = 1;

  if (h_sbr_data->coupling == COUPLING_BAL) {
    target = (FIXP_SGL)SBR_ENERGY_PAN_OFFSET;
    step = (FIXP_SGL)DECAY_COUPLING;
  } else {
    target = FL2FXCONST_SGL(0.0f);
    step = (FIXP_SGL)DECAY;
  }
  if (hHeaderData->bs_info.ampResolution == 0) {
    target <<= 1;
    step <<= 1;
  }

  for (i = 0; i < h_sbr_data->nScaleFactors; i++) {
    if (h_prev_data->sfb_nrg_prev[i] > target)
      h_sbr_data->iEnvelope[i] = -step;
    else
      h_sbr_data->iEnvelope[i] = step;
  }

  h_sbr_data->domain_vec_noise[0] = 1;
  FDKmemclear(h_sbr_data->sbrNoiseFloorLevel,
              sizeof(h_sbr_data->sbrNoiseFloorLevel));

  FDKmemclear(h_sbr_data->addHarmonics,
              sizeof(ULONG) * ADD_HARMONICS_FLAGS_SIZE);
}

static void decodeEnvelope(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_sbr_data,
    HANDLE_SBR_PREV_FRAME_DATA
        h_prev_data,
    HANDLE_SBR_PREV_FRAME_DATA
        otherChannel
) {
  int i;
  int fFrameError = hHeaderData->frameErrorFlag;
  FIXP_SGL tempSfbNrgPrev[MAX_FREQ_COEFFS];

  if (!fFrameError) {
    if (h_prev_data->frameErrorFlag) {
      if (h_sbr_data->domain_vec[0] != 0) {
        fFrameError = 1;
      }
    } else {
      if (h_sbr_data->frameInfo.borders[0] !=
          h_prev_data->stopPos - hHeaderData->numberTimeSlots) {
        if (h_sbr_data->domain_vec[0] == 1) {
          fFrameError = 1;
        } else {
          fFrameError = 1;
        }
      }
    }
  }

  if (fFrameError)
  {
    leanSbrConcealment(hHeaderData, h_sbr_data, h_prev_data);

    deltaToLinearPcmEnvelopeDecoding(hHeaderData, h_sbr_data, h_prev_data);
  } else
  {
    if (h_prev_data->frameErrorFlag) {
      timeCompensateFirstEnvelope(hHeaderData, h_sbr_data, h_prev_data);
      if (h_sbr_data->coupling != h_prev_data->coupling) {
        for (i = 0; i < hHeaderData->freqBandData.nSfb[1]; i++) {
          if (h_prev_data->coupling == COUPLING_BAL) {
            h_prev_data->sfb_nrg_prev[i] =
                (otherChannel != NULL) ? otherChannel->sfb_nrg_prev[i]
                                       : (FIXP_SGL)SBR_ENERGY_PAN_OFFSET;
          }
          else if (h_sbr_data->coupling == COUPLING_LEVEL &&
                   otherChannel != NULL) {
            h_prev_data->sfb_nrg_prev[i] = (h_prev_data->sfb_nrg_prev[i] +
                                            otherChannel->sfb_nrg_prev[i]) >>
                                           1;
          } else if (h_sbr_data->coupling == COUPLING_BAL) {
            h_prev_data->sfb_nrg_prev[i] = (FIXP_SGL)SBR_ENERGY_PAN_OFFSET;
          }
        }
      }
    }
    FDKmemcpy(tempSfbNrgPrev, h_prev_data->sfb_nrg_prev,
              MAX_FREQ_COEFFS * sizeof(FIXP_SGL));

    deltaToLinearPcmEnvelopeDecoding(hHeaderData, h_sbr_data, h_prev_data);

    fFrameError = checkEnvelopeData(hHeaderData, h_sbr_data, h_prev_data);

    if (fFrameError) {
      hHeaderData->frameErrorFlag = 1;
      FDKmemcpy(h_prev_data->sfb_nrg_prev, tempSfbNrgPrev,
                MAX_FREQ_COEFFS * sizeof(FIXP_SGL));
      decodeEnvelope(hHeaderData, h_sbr_data, h_prev_data, otherChannel);
      return;
    }
  }

  requantizeEnvelopeData(h_sbr_data, h_sbr_data->ampResolutionCurrentFrame);

  hHeaderData->frameErrorFlag = fFrameError;
}

static int checkEnvelopeData(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_sbr_data,
    HANDLE_SBR_PREV_FRAME_DATA h_prev_data
) {
  FIXP_SGL *iEnvelope = h_sbr_data->iEnvelope;
  FIXP_SGL *sfb_nrg_prev = h_prev_data->sfb_nrg_prev;
  int i = 0, errorFlag = 0;
  FIXP_SGL sbr_max_energy = (h_sbr_data->ampResolutionCurrentFrame == 1)
                                ? SBR_MAX_ENERGY
                                : (SBR_MAX_ENERGY << 1);

  for (i = 0; i < h_sbr_data->nScaleFactors; i++) {
    if (iEnvelope[i] > sbr_max_energy) {
      errorFlag = 1;
    }
    if (iEnvelope[i] < FL2FXCONST_SGL(0.0f)) {
      errorFlag = 1;
    }
  }

  for (i = 0; i < hHeaderData->freqBandData.nSfb[1]; i++) {
    sfb_nrg_prev[i] = fixMax(sfb_nrg_prev[i], FL2FXCONST_SGL(0.0f));
    sfb_nrg_prev[i] = fixMin(sfb_nrg_prev[i], sbr_max_energy);
  }

  return (errorFlag);
}

static void limitNoiseLevels(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_sbr_data)
{
  int i;
  int nNfb = hHeaderData->freqBandData.nNfb;

#define lowerLimit \
  ((FIXP_SGL)0)
#define upperLimit \
  ((FIXP_SGL)35)

  for (i = 0; i < h_sbr_data->frameInfo.nNoiseEnvelopes * nNfb; i++) {
    h_sbr_data->sbrNoiseFloorLevel[i] =
        fixMin(h_sbr_data->sbrNoiseFloorLevel[i], upperLimit);
    h_sbr_data->sbrNoiseFloorLevel[i] =
        fixMax(h_sbr_data->sbrNoiseFloorLevel[i], lowerLimit);
  }
}

static void timeCompensateFirstEnvelope(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_sbr_data,
    HANDLE_SBR_PREV_FRAME_DATA
        h_prev_data)
{
  int i, nScalefactors;
  FRAME_INFO *pFrameInfo = &h_sbr_data->frameInfo;
  UCHAR *nSfb = hHeaderData->freqBandData.nSfb;
  int estimatedStartPos =
      fMax(0, h_prev_data->stopPos - hHeaderData->numberTimeSlots);
  int refLen, newLen, shift;
  FIXP_SGL deltaExp;

  refLen = pFrameInfo->borders[1] - pFrameInfo->borders[0];
  newLen = pFrameInfo->borders[1] - estimatedStartPos;

  if (newLen <= 0) {
    newLen = refLen;
    estimatedStartPos = pFrameInfo->borders[0];
  }

  deltaExp = FDK_getNumOctavesDiv8(newLen, refLen);

  shift = (FRACT_BITS - 1 - ENV_EXP_FRACT - 1 +
           h_sbr_data->ampResolutionCurrentFrame - 3);
  deltaExp = deltaExp >> shift;
  pFrameInfo->borders[0] = estimatedStartPos;
  pFrameInfo->bordersNoise[0] = estimatedStartPos;

  if (h_sbr_data->coupling != COUPLING_BAL) {
    nScalefactors = (pFrameInfo->freqRes[0]) ? nSfb[1] : nSfb[0];

    for (i = 0; i < nScalefactors; i++)
      h_sbr_data->iEnvelope[i] = h_sbr_data->iEnvelope[i] + deltaExp;
  }
}

static void requantizeEnvelopeData(HANDLE_SBR_FRAME_DATA h_sbr_data,
                                   int ampResolution) {
  int i;
  FIXP_SGL mantissa;
  int ampShift = 1 - ampResolution;
  int exponent;

#if ENV_EXP_FRACT
  static const FIXP_SGL pow2[ENV_EXP_FRACT] = {
      FL2FXCONST_SGL(0.5f * pow(2.0f, pow(0.5f, 1))),
      FL2FXCONST_SGL(0.5f * pow(2.0f, pow(0.5f, 2))),
      FL2FXCONST_SGL(0.5f * pow(2.0f, pow(0.5f, 3))),
      FL2FXCONST_SGL(0.5f * pow(2.0f, pow(0.5f, 4))),
      FL2FXCONST_SGL(0.5f * pow(2.0f, pow(0.5f, 5))),
      FL2FXCONST_SGL(0.5f * pow(2.0f, pow(0.5f, 6))),
      FL2FXCONST_SGL(0.5f * pow(2.0f, pow(0.5f, 7))),
      FL2FXCONST_SGL(0.5f * pow(2.0f, pow(0.5f, 8)))
  };

  int bit, mask;
#endif

  for (i = 0; i < h_sbr_data->nScaleFactors; i++) {
    exponent = (LONG)h_sbr_data->iEnvelope[i];

#if ENV_EXP_FRACT

    exponent = exponent >> ampShift;
    mantissa = 0.5f;

    mask = 1;

    for (bit = ENV_EXP_FRACT - 1; bit >= 0; bit--) {
      if (exponent & mask) {
        mantissa = (FIXP_SGL)((mantissa * pow2[bit]) << 1);
      }
      mask = mask << 1;
    }

    exponent = exponent >> ENV_EXP_FRACT;

#else
    mantissa = (exponent & ampShift) ? FL2FXCONST_SGL(0.707106781186548f)
                                     : FL2FXCONST_SGL(0.5f);
    exponent = exponent >> ampShift;
#endif

    exponent += 7 + NRG_EXP_OFFSET;

    h_sbr_data->iEnvelope[i] =
        ((FIXP_SGL)((SHORT)(FIXP_SGL)mantissa & MASK_M)) +
        (FIXP_SGL)((SHORT)(FIXP_SGL)exponent & MASK_E);
  }
}

static void deltaToLinearPcmEnvelopeDecoding(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_sbr_data,
    HANDLE_SBR_PREV_FRAME_DATA h_prev_data)
{
  int i, domain, no_of_bands, band, freqRes;

  FIXP_SGL *sfb_nrg_prev = h_prev_data->sfb_nrg_prev;
  FIXP_SGL *ptr_nrg = h_sbr_data->iEnvelope;

  int offset =
      2 * hHeaderData->freqBandData.nSfb[0] - hHeaderData->freqBandData.nSfb[1];

  for (i = 0; i < h_sbr_data->frameInfo.nEnvelopes; i++) {
    domain = h_sbr_data->domain_vec[i];
    freqRes = h_sbr_data->frameInfo.freqRes[i];

    FDK_ASSERT(freqRes >= 0 && freqRes <= 1);

    no_of_bands = hHeaderData->freqBandData.nSfb[freqRes];

    FDK_ASSERT(no_of_bands < (64));

    if (domain == 0) {
      mapLowResEnergyVal(*ptr_nrg, sfb_nrg_prev, offset, 0, freqRes);
      ptr_nrg++;
      for (band = 1; band < no_of_bands; band++) {
        *ptr_nrg = *ptr_nrg + *(ptr_nrg - 1);
        mapLowResEnergyVal(*ptr_nrg, sfb_nrg_prev, offset, band, freqRes);
        ptr_nrg++;
      }
    } else {
      for (band = 0; band < no_of_bands; band++) {
        *ptr_nrg =
            *ptr_nrg + sfb_nrg_prev[indexLow2High(offset, band, freqRes)];
        mapLowResEnergyVal(*ptr_nrg, sfb_nrg_prev, offset, band, freqRes);
        ptr_nrg++;
      }
    }
  }
}

static void decodeNoiseFloorlevels(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_sbr_data,
    HANDLE_SBR_PREV_FRAME_DATA h_prev_data)
{
  int i;
  int nNfb = hHeaderData->freqBandData.nNfb;
  int nNoiseFloorEnvelopes = h_sbr_data->frameInfo.nNoiseEnvelopes;

  if (h_sbr_data->domain_vec_noise[0] == 0) {
    FIXP_SGL noiseLevel = h_sbr_data->sbrNoiseFloorLevel[0];
    for (i = 1; i < nNfb; i++) {
      noiseLevel += h_sbr_data->sbrNoiseFloorLevel[i];
      h_sbr_data->sbrNoiseFloorLevel[i] = noiseLevel;
    }
  } else {
    for (i = 0; i < nNfb; i++) {
      h_sbr_data->sbrNoiseFloorLevel[i] += h_prev_data->prevNoiseLevel[i];
    }
  }

  if (nNoiseFloorEnvelopes > 1) {
    if (h_sbr_data->domain_vec_noise[1] == 0) {
      FIXP_SGL noiseLevel = h_sbr_data->sbrNoiseFloorLevel[nNfb];
      for (i = nNfb + 1; i < 2 * nNfb; i++) {
        noiseLevel += h_sbr_data->sbrNoiseFloorLevel[i];
        h_sbr_data->sbrNoiseFloorLevel[i] = noiseLevel;
      }
    } else {
      for (i = 0; i < nNfb; i++) {
        h_sbr_data->sbrNoiseFloorLevel[i + nNfb] +=
            h_sbr_data->sbrNoiseFloorLevel[i];
      }
    }
  }

  limitNoiseLevels(hHeaderData, h_sbr_data);

  for (i = 0; i < nNfb; i++)
    h_prev_data->prevNoiseLevel[i] =
        h_sbr_data->sbrNoiseFloorLevel[i + nNfb * (nNoiseFloorEnvelopes - 1)];

  if (!h_sbr_data->coupling) {
    int nf_e;

    for (i = 0; i < nNoiseFloorEnvelopes * nNfb; i++) {
      nf_e = 6 - (LONG)h_sbr_data->sbrNoiseFloorLevel[i] + 1 + NOISE_EXP_OFFSET;

      h_sbr_data->sbrNoiseFloorLevel[i] =
          (FIXP_SGL)(((LONG)FL2FXCONST_SGL(0.5f)) +
                     (nf_e & MASK_E));
    }
  }
}
