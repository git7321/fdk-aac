/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "env_calc.h"
#include "sbrdec_freq_sca.h"
#include "env_extr.h"
#include "transcendent.h"
#include "sbr_ram.h"
#include "sbr_rom.h"
#include "genericStds.h"

#define MAX_SFB_NRG_HEADROOM (1)
#define MAX_VAL_NRG_HEADROOM ((((FIXP_DBL)MAXVAL_DBL) >> MAX_SFB_NRG_HEADROOM))

typedef struct {
  FIXP_DBL nrgRef[MAX_FREQ_COEFFS];
  FIXP_DBL nrgEst[MAX_FREQ_COEFFS];
  FIXP_DBL nrgGain[MAX_FREQ_COEFFS];
  FIXP_DBL noiseLevel[MAX_FREQ_COEFFS];
  FIXP_DBL nrgSine[MAX_FREQ_COEFFS];

  SCHAR nrgRef_e[MAX_FREQ_COEFFS];
  SCHAR nrgEst_e[MAX_FREQ_COEFFS];
  SCHAR nrgGain_e[MAX_FREQ_COEFFS];
  SCHAR noiseLevel_e[MAX_FREQ_COEFFS];
  SCHAR nrgSine_e[MAX_FREQ_COEFFS];
  SCHAR exponent[2];
} ENV_CALC_NRGS;

static void equalizeFiltBufferExp(FIXP_DBL *filtBuffer, SCHAR *filtBuffer_e,
                                  FIXP_DBL *NrgGain, SCHAR *NrgGain_e,
                                  int subbands);

static void calcNrgPerSubband(FIXP_DBL **analysBufferReal,
                              FIXP_DBL **analysBufferImag, int lowSubband,
                              int highSubband, int start_pos, int next_pos,
                              SCHAR frameExp, FIXP_DBL *nrgEst,
                              SCHAR *nrgEst_e);

static void calcNrgPerSfb(FIXP_DBL **analysBufferReal,
                          FIXP_DBL **analysBufferImag, int nSfb,
                          UCHAR *freqBandTable, int start_pos, int next_pos,
                          SCHAR input_e, FIXP_DBL *nrg_est, SCHAR *nrg_est_e);

static void calcSubbandGain(FIXP_DBL nrgRef, SCHAR nrgRef_e,
                            ENV_CALC_NRGS *nrgs, int c, FIXP_DBL tmpNoise,
                            SCHAR tmpNoise_e, UCHAR sinePresentFlag,
                            UCHAR sineMapped, int noNoiseFlag);

static void calcAvgGain(ENV_CALC_NRGS *nrgs, int lowSubband, int highSubband,
                        FIXP_DBL *sumRef_m, SCHAR *sumRef_e,
                        FIXP_DBL *ptrAvgGain_m, SCHAR *ptrAvgGain_e);

static void adjustTimeSlot_EldGrid(FIXP_DBL *ptrReal, ENV_CALC_NRGS *nrgs,
                                   UCHAR *ptrHarmIndex, int lowSubbands,
                                   int noSubbands, int scale_change,
                                   int noNoiseFlag, int *ptrPhaseIndex,
                                   int scale_diff_low);

static void adjustTimeSlotLC(FIXP_DBL *ptrReal, ENV_CALC_NRGS *nrgs,
                             UCHAR *ptrHarmIndex, int lowSubbands,
                             int noSubbands, int scale_change, int noNoiseFlag,
                             int *ptrPhaseIndex);

static void adjustTimeSlotHQ_GainAndNoise(
    FIXP_DBL *ptrReal, FIXP_DBL *ptrImag,
    HANDLE_SBR_CALCULATE_ENVELOPE h_sbr_cal_env, ENV_CALC_NRGS *nrgs,
    int lowSubbands, int noSubbands, int scale_change, FIXP_SGL smooth_ratio,
    int noNoiseFlag, int filtBufferNoiseShift);
static void adjustTimeSlotHQ_AddHarmonics(
    FIXP_DBL *ptrReal, FIXP_DBL *ptrImag,
    HANDLE_SBR_CALCULATE_ENVELOPE h_sbr_cal_env, ENV_CALC_NRGS *nrgs,
    int lowSubbands, int noSubbands, int scale_change);

static void adjustTimeSlotHQ(FIXP_DBL *ptrReal, FIXP_DBL *ptrImag,
                             HANDLE_SBR_CALCULATE_ENVELOPE h_sbr_cal_env,
                             ENV_CALC_NRGS *nrgs, int lowSubbands,
                             int noSubbands, int scale_change,
                             FIXP_SGL smooth_ratio, int noNoiseFlag,
                             int filtBufferNoiseShift);

static void mapSineFlags(
    UCHAR *freqBandTable,
    int nSfb,
    ULONG *addHarmonics,
    ULONG *harmFlagsPrev,
    ULONG *harmFlagsPrevActive,
    int tranEnv,
    SCHAR *sineMapped)

{
  int i;
  int bitcount = 31;
  ULONG harmFlagsQmfBands[ADD_HARMONICS_FLAGS_SIZE] = {0};
  ULONG *curFlags = addHarmonics;

  FDKmemset(sineMapped, 32,
            MAX_FREQ_COEFFS * sizeof(SCHAR));
  FDKmemclear(harmFlagsPrevActive, ADD_HARMONICS_FLAGS_SIZE * sizeof(ULONG));
  for (i = 0; i < nSfb; i++) {
    ULONG maskSfb =
        1 << bitcount;

    if (*curFlags & maskSfb) {
      const int lsb = freqBandTable[0];
      const int qmfBand = (freqBandTable[i] + freqBandTable[i + 1]) >> 1;
      const int qmfBandDiv32 = qmfBand >> 5;
      const int maskQmfBand =
          1 << (qmfBand &
                31);

      harmFlagsQmfBands[qmfBandDiv32] |= maskQmfBand;

      sineMapped[qmfBand - lsb] =
          (harmFlagsPrev[qmfBandDiv32] & maskQmfBand) ? 0 : tranEnv;
      if (sineMapped[qmfBand - lsb] < PVC_NTIMESLOT) {
        harmFlagsPrevActive[qmfBandDiv32] |= maskQmfBand;
      }
    }

    if (bitcount-- == 0) {
      bitcount = 31;
      curFlags++;
    }
  }
  FDKmemcpy(harmFlagsPrev, harmFlagsQmfBands,
            sizeof(ULONG) * ADD_HARMONICS_FLAGS_SIZE);
}

void mapSineFlagsPvc(
    UCHAR *freqBandTable,
    int nSfb,
    ULONG *harmFlagsPrev,
    ULONG *harmFlagsPrevActive,
    SCHAR *sineMapped,
    int sinusoidalPos,
    SCHAR *sinusoidalPosPrev,
    int trailingSbrFrame)
{
  FDKmemset(sineMapped, 32, MAX_FREQ_COEFFS);

  if (trailingSbrFrame) {
    int i;
    const int lsb = freqBandTable[0];
    const int usb = freqBandTable[nSfb];
    for (i = lsb; i < usb; i++) {
      const int qmfBandDiv32 = i >> 5;
      const int maskQmfBand =
          1 << (i & 31);

      if (harmFlagsPrevActive[qmfBandDiv32] & maskQmfBand) {
        sineMapped[i - lsb] = 0;
      } else if (harmFlagsPrev[qmfBandDiv32] & maskQmfBand) {
        sineMapped[i - lsb] =
            *sinusoidalPosPrev - PVC_NTIMESLOT;
      }
    }
  }
  *sinusoidalPosPrev = sinusoidalPos;
}

void aliasingReduction(
    FIXP_DBL *degreeAlias,
    ENV_CALC_NRGS *nrgs,
    UCHAR *useAliasReduction,
    int noSubbands)
{
  FIXP_DBL *nrgGain = nrgs->nrgGain;
  SCHAR *nrgGain_e =
      nrgs->nrgGain_e;
  FIXP_DBL *nrgEst = nrgs->nrgEst;
  SCHAR *nrgEst_e =
      nrgs->nrgEst_e;
  int grouping = 0, index = 0, noGroups, k;
  int groupVector[MAX_FREQ_COEFFS];

  for (k = 0; k < noSubbands - 1; k++) {
    if ((degreeAlias[k + 1] != FL2FXCONST_DBL(0.0f)) && useAliasReduction[k]) {
      if (grouping == 0) {
        groupVector[index++] = k;
        grouping = 1;
      } else {
        if (groupVector[index - 1] + 3 == k) {
          groupVector[index++] = k + 1;
          grouping = 0;
        }
      }
    } else {
      if (grouping) {
        if (useAliasReduction[k])
          groupVector[index++] = k + 1;
        else
          groupVector[index++] = k;
        grouping = 0;
      }
    }
  }

  if (grouping) {
    groupVector[index++] = noSubbands;
  }
  noGroups = index >> 1;

  for (int group = 0; group < noGroups; group++) {
    FIXP_DBL nrgOrig = FL2FXCONST_DBL(
        0.0f);
    SCHAR nrgOrig_e = 0;
    FIXP_DBL nrgAmp = FL2FXCONST_DBL(
        0.0f);
    SCHAR nrgAmp_e = 0;
    FIXP_DBL nrgMod = FL2FXCONST_DBL(
        0.0f);
    SCHAR nrgMod_e = 0;
    FIXP_DBL groupGain;
    SCHAR groupGain_e;
    FIXP_DBL compensation;
    SCHAR compensation_e;

    int startGroup = groupVector[2 * group];
    int stopGroup = groupVector[2 * group + 1];

    for (k = startGroup; k < stopGroup; k++) {
      FIXP_DBL tmp = nrgEst[k];
      SCHAR tmp_e = nrgEst_e[k];

      FDK_add_MantExp(tmp, tmp_e, nrgOrig, nrgOrig_e, &nrgOrig, &nrgOrig_e);

      tmp = fMult(tmp, nrgGain[k]);
      tmp_e = tmp_e + nrgGain_e[k];

      FDK_add_MantExp(tmp, tmp_e, nrgAmp, nrgAmp_e, &nrgAmp, &nrgAmp_e);
    }

    FDK_divide_MantExp(nrgAmp, nrgAmp_e, nrgOrig, nrgOrig_e, &groupGain,
                       &groupGain_e);

    for (k = startGroup; k < stopGroup; k++) {
      FIXP_DBL tmp;
      SCHAR tmp_e;

      FIXP_DBL alpha = degreeAlias[k];
      if (k < noSubbands - 1) {
        if (degreeAlias[k + 1] > alpha) alpha = degreeAlias[k + 1];
      }

      FDK_add_MantExp(
          fMult(alpha, groupGain), groupGain_e,
          fMult((FIXP_DBL)MAXVAL_DBL - alpha,
                nrgGain[k]),
          nrgGain_e[k], &nrgGain[k], &nrgGain_e[k]);

      tmp = fMult(nrgGain[k], nrgEst[k]);
      tmp_e = nrgGain_e[k] + nrgEst_e[k];

      FDK_add_MantExp(tmp, tmp_e, nrgMod, nrgMod_e, &nrgMod, &nrgMod_e);
    }

    FDK_divide_MantExp(nrgAmp, nrgAmp_e, nrgMod, nrgMod_e, &compensation,
                       &compensation_e);

    for (k = startGroup; k < stopGroup; k++) {
      nrgGain[k] = fMult(nrgGain[k], compensation);
      nrgGain_e[k] = nrgGain_e[k] + compensation_e;
    }
  }
}

#define INTER_TES_SF_CHANGE 4

typedef struct {
  FIXP_DBL subsample_power_low[(((1024) / (32) * (4) / 2) + (3 * (4)))];
  FIXP_DBL subsample_power_high[(((1024) / (32) * (4) / 2) + (3 * (4)))];
  FIXP_DBL gain[(((1024) / (32) * (4) / 2) + (3 * (4)))];
  SCHAR subsample_power_low_sf[(((1024) / (32) * (4) / 2) + (3 * (4)))];
  SCHAR subsample_power_high_sf[(((1024) / (32) * (4) / 2) + (3 * (4)))];
} ITES_TEMP;

static void apply_inter_tes(FIXP_DBL **qmfReal, FIXP_DBL **qmfImag,
                            const QMF_SCALE_FACTOR *sbrScaleFactor,
                            const SCHAR exp[2], const int RATE,
                            const int startPos, const int stopPos,
                            const int lowSubband, const int nbSubband,
                            const UCHAR gamma_idx) {
  int highSubband = lowSubband + nbSubband;
  FIXP_DBL *subsample_power_high, *subsample_power_low;
  SCHAR *subsample_power_high_sf, *subsample_power_low_sf;
  FIXP_DBL total_power_high = (FIXP_DBL)0;
  FIXP_DBL total_power_low = (FIXP_DBL)0;
  FIXP_DBL *gain;
  int gain_sf[(((1024) / (32) * (4) / 2) + (3 * (4)))];

  int gamma_sf =
      (int)gamma_idx - 1;

  int nbSubsample = stopPos - startPos;
  int i, j;

  C_ALLOC_SCRATCH_START(pTmp, ITES_TEMP, 1);
  subsample_power_high = pTmp->subsample_power_high;
  subsample_power_low = pTmp->subsample_power_low;
  subsample_power_high_sf = pTmp->subsample_power_high_sf;
  subsample_power_low_sf = pTmp->subsample_power_low_sf;
  gain = pTmp->gain;

  if (gamma_idx > 0) {
    int preShift2 = 32 - fNormz((FIXP_DBL)nbSubsample);
    int total_power_low_sf = 1 - DFRACT_BITS;
    int total_power_high_sf = 1 - DFRACT_BITS;

    for (i = 0; i < nbSubsample; ++i) {
      FIXP_DBL bufferReal[(((1024) / (32) * (4) / 2) + (3 * (4)))];
      FIXP_DBL bufferImag[(((1024) / (32) * (4) / 2) + (3 * (4)))];
      FIXP_DBL maxVal = (FIXP_DBL)0;

      int ts = startPos + i;

      int low_sf = (ts < 3 * RATE) ? sbrScaleFactor->ov_lb_scale
                                   : sbrScaleFactor->lb_scale;
      low_sf = 15 - low_sf;

      for (j = 0; j < lowSubband; ++j) {
        bufferImag[j] = qmfImag[startPos + i][j];
        maxVal |= (FIXP_DBL)((LONG)(bufferImag[j]) ^
                             ((LONG)bufferImag[j] >> (DFRACT_BITS - 1)));
        bufferReal[j] = qmfReal[startPos + i][j];
        maxVal |= (FIXP_DBL)((LONG)(bufferReal[j]) ^
                             ((LONG)bufferReal[j] >> (DFRACT_BITS - 1)));
      }

      subsample_power_low[i] = (FIXP_DBL)0;
      subsample_power_low_sf[i] = 0;

      if (maxVal != FL2FXCONST_DBL(0.f)) {
        int preShift = 1 - CntLeadingZeros(maxVal);
        int postShift = 32 - fNormz((FIXP_DBL)lowSubband);

        if (preShift != 0) preShift++;

        subsample_power_low_sf[i] += (low_sf + preShift) * 2 + postShift + 1;

        scaleValues(bufferReal, lowSubband, -preShift);
        scaleValues(bufferImag, lowSubband, -preShift);
        for (j = 0; j < lowSubband; ++j) {
          FIXP_DBL addme;
          addme = fPow2Div2(bufferReal[j]);
          subsample_power_low[i] += addme >> postShift;
          addme = fPow2Div2(bufferImag[j]);
          subsample_power_low[i] += addme >> postShift;
        }
      }

      maxVal = (FIXP_DBL)0;

      int high_sf = exp[(ts < 16 * RATE) ? 0 : 1];

      for (j = lowSubband; j < highSubband; ++j) {
        bufferImag[j] = qmfImag[startPos + i][j];
        maxVal |= (FIXP_DBL)((LONG)(bufferImag[j]) ^
                             ((LONG)bufferImag[j] >> (DFRACT_BITS - 1)));
        bufferReal[j] = qmfReal[startPos + i][j];
        maxVal |= (FIXP_DBL)((LONG)(bufferReal[j]) ^
                             ((LONG)bufferReal[j] >> (DFRACT_BITS - 1)));
      }

      subsample_power_high[i] = (FIXP_DBL)0;
      subsample_power_high_sf[i] = 0;

      if (maxVal != FL2FXCONST_DBL(0.f)) {
        int preShift = 1 - CntLeadingZeros(maxVal);
        if (preShift != 0) preShift++;

        int postShift = 32 - fNormz((FIXP_DBL)(highSubband - lowSubband));
        subsample_power_high_sf[i] += (high_sf + preShift) * 2 + postShift + 1;

        scaleValues(&bufferReal[lowSubband], highSubband - lowSubband,
                    -preShift);
        scaleValues(&bufferImag[lowSubband], highSubband - lowSubband,
                    -preShift);
        for (j = lowSubband; j < highSubband; j++) {
          subsample_power_high[i] += fPow2Div2(bufferReal[j]) >> postShift;
          subsample_power_high[i] += fPow2Div2(bufferImag[j]) >> postShift;
        }
      }

      FIXP_DBL new_summand = subsample_power_low[i];
      int new_summand_sf = subsample_power_low_sf[i];

      if (new_summand_sf > total_power_low_sf) {
        int diff = fMin(DFRACT_BITS - 1, new_summand_sf - total_power_low_sf);
        total_power_low >>= diff;
        total_power_low_sf = new_summand_sf;
      } else if (new_summand_sf < total_power_low_sf) {
        new_summand >>=
            fMin(DFRACT_BITS - 1, total_power_low_sf - new_summand_sf);
      }

      total_power_low += (new_summand >> preShift2);

      new_summand = subsample_power_high[i];
      new_summand_sf = subsample_power_high_sf[i];
      if (new_summand_sf > total_power_high_sf) {
        total_power_high >>=
            fMin(DFRACT_BITS - 1, new_summand_sf - total_power_high_sf);
        total_power_high_sf = new_summand_sf;
      } else if (new_summand_sf < total_power_high_sf) {
        new_summand >>=
            fMin(DFRACT_BITS - 1, total_power_high_sf - new_summand_sf);
      }

      total_power_high += (new_summand >> preShift2);
    }

    total_power_low_sf += preShift2;
    total_power_high_sf += preShift2;

    for (i = 0; i < nbSubsample; ++i) {
      int sf2;
      FIXP_DBL mult =
          fMultNorm(subsample_power_low[i], (FIXP_DBL)nbSubsample, &sf2);
      int mult_sf = subsample_power_low_sf[i] + DFRACT_BITS - 1 + sf2;

      if (total_power_low != FIXP_DBL(0)) {
        gain[i] = fDivNorm(mult, total_power_low, &sf2);
        gain_sf[i] = mult_sf - total_power_low_sf + sf2;
        gain[i] = sqrtFixp_lookup(gain[i], &gain_sf[i]);
        if (gain_sf[i] < 0) {
          gain[i] >>= fMin(DFRACT_BITS - 1, -gain_sf[i]);
          gain_sf[i] = 0;
        }
      } else {
        if (mult == FIXP_DBL(0)) {
          gain[i] = FIXP_DBL(0);
          gain_sf[i] = 0;
        } else {
          gain[i] = (FIXP_DBL)MAXVAL_DBL;
          gain_sf[i] = 0;
        }
      }
    }

    FIXP_DBL total_power_high_after = (FIXP_DBL)0;
    int total_power_high_after_sf = 1 - DFRACT_BITS;

    for (i = 0; i < nbSubsample; ++i) {
      FIXP_DBL one = (FIXP_DBL)MAXVAL_DBL >>
                     gain_sf[i];

      FIXP_DBL mult = (gain[i] - one) >> 1;
      int mult_sf = gain_sf[i] + gamma_sf;

      one = FL2FXCONST_DBL(0.5f) >> mult_sf;
      gain[i] = one + mult;
      gain_sf[i] += gamma_sf + 1;

      FIXP_DBL gain_pow2;
      int gain_pow2_sf;

      if (fIsLessThan(gain[i], gain_sf[i], FL2FXCONST_DBL(0.2f), 0)) {
        gain[i] = FL2FXCONST_DBL(0.8f);
        gain_sf[i] = -2;
        gain_pow2 = FL2FXCONST_DBL(0.64f);
        gain_pow2_sf = -4;
      } else {
        int r = CountLeadingBits(gain[i]);
        gain[i] <<= r;
        gain_sf[i] -= r;

        gain_pow2 = fPow2(gain[i]);
        gain_pow2_sf = gain_sf[i] << 1;
      }

      int room;
      subsample_power_high[i] =
          fMultNorm(subsample_power_high[i], gain_pow2, &room);
      subsample_power_high_sf[i] =
          subsample_power_high_sf[i] + gain_pow2_sf + room;

      int new_summand_sf = subsample_power_high_sf[i];
      if (new_summand_sf > total_power_high_after_sf) {
        total_power_high_after >>=
            fMin(DFRACT_BITS - 1, new_summand_sf - total_power_high_after_sf);
        total_power_high_after_sf = new_summand_sf;
      } else if (new_summand_sf < total_power_high_after_sf) {
        subsample_power_high[i] >>=
            fMin(DFRACT_BITS - 1, total_power_high_after_sf - new_summand_sf);
      }
      total_power_high_after += subsample_power_high[i] >> preShift2;
    }

    total_power_high_after_sf += preShift2;

    int sf2 = 0;
    FIXP_DBL gain_adj_2 = FL2FX_DBL(0.5f);
    int gain_adj_2_sf = 1;

    if ((total_power_high != (FIXP_DBL)0) &&
        (total_power_high_after != (FIXP_DBL)0)) {
      gain_adj_2 = fDivNorm(total_power_high, total_power_high_after, &sf2);
      gain_adj_2_sf = total_power_high_sf - total_power_high_after_sf + sf2;
    }

    FIXP_DBL gain_adj = sqrtFixp_lookup(gain_adj_2, &gain_adj_2_sf);
    int gain_adj_sf = gain_adj_2_sf;

    for (i = 0; i < nbSubsample; ++i) {
      int gain_e = fMax(
          fMin(gain_sf[i] + gain_adj_sf - INTER_TES_SF_CHANGE, DFRACT_BITS - 1),
          -(DFRACT_BITS - 1));
      FIXP_DBL gain_final = fMult(gain[i], gain_adj);
      gain_final = scaleValueSaturate(gain_final, gain_e);

      for (j = lowSubband; j < highSubband; j++) {
        qmfReal[startPos + i][j] = fMult(qmfReal[startPos + i][j], gain_final);
        qmfImag[startPos + i][j] = fMult(qmfImag[startPos + i][j], gain_final);
      }
    }
  } else {
    for (i = 0; i < nbSubsample; ++i) {
      for (j = lowSubband; j < highSubband; j++) {
        qmfReal[startPos + i][j] >>= INTER_TES_SF_CHANGE;
        qmfImag[startPos + i][j] >>= INTER_TES_SF_CHANGE;
      }
    }
  }
  C_ALLOC_SCRATCH_END(pTmp, ITES_TEMP, 1);
}

void calculateSbrEnvelope(
    QMF_SCALE_FACTOR *sbrScaleFactor,
    HANDLE_SBR_CALCULATE_ENVELOPE
        h_sbr_cal_env,
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA hFrameData,
    PVC_DYNAMIC_DATA *pPvcDynamicData,
    FIXP_DBL *
        *analysBufferReal,
    FIXP_DBL *
        *analysBufferImag,
    const int useLP,
    FIXP_DBL *degreeAlias,
    const UINT flags, const int frameErrorFlag) {
  int c, i, i_stop, j, envNoise = 0;
  UCHAR *borders = hFrameData->frameInfo.borders;
  UCHAR *bordersPvc = hFrameData->frameInfo.pvcBorders;
  int pvc_mode = pPvcDynamicData->pvc_mode;
  int first_start =
      ((pvc_mode > 0) ? bordersPvc[0] : borders[0]) * hHeaderData->timeStep;
  FIXP_SGL *noiseLevels = hFrameData->sbrNoiseFloorLevel;
  HANDLE_FREQ_BAND_DATA hFreq = &hHeaderData->freqBandData;
  UCHAR **pFreqBandTable = hFreq->freqBandTable;
  UCHAR *pFreqBandTableNoise = hFreq->freqBandTableNoise;

  int lowSubband = hFreq->lowSubband;
  int highSubband = hFreq->highSubband;
  int noSubbands = highSubband - lowSubband;

  int ov_highSubband = hFreq->highSubband;

  int noNoiseBands = hFreq->nNfb;
  UCHAR *noSubFrameBands = hFreq->nSfb;
  int no_cols = hHeaderData->numberTimeSlots * hHeaderData->timeStep;

  SCHAR sineMapped[MAX_FREQ_COEFFS];
  SCHAR ov_adj_e = SCALE2EXP(sbrScaleFactor->ov_hb_scale);
  SCHAR adj_e = 0;
  SCHAR output_e;
  SCHAR final_e = 0;
  const int iTES_enable = hFrameData->iTESactive;
  const int iTES_scale_change = (iTES_enable) ? INTER_TES_SF_CHANGE : 0;
  SCHAR maxGainLimit_e = (frameErrorFlag) ? MAX_GAIN_CONCEAL_EXP : MAX_GAIN_EXP;

  UCHAR smooth_length = 0;

  FIXP_SGL *pIenv = hFrameData->iEnvelope;

  C_ALLOC_SCRATCH_START(useAliasReduction, UCHAR, 64)

  if (hFreq->highSubband < hFreq->ov_highSubband) {
    ov_highSubband = hFreq->ov_highSubband;
  }

  if (pvc_mode > 0) {
    if (hFrameData->frameInfo.bordersNoise[0] > bordersPvc[0]) {
      envNoise = -1;
      noiseLevels = h_sbr_cal_env->prevSbrNoiseFloorLevel;
      noNoiseBands = h_sbr_cal_env->prevNNfb;
      noSubFrameBands = h_sbr_cal_env->prevNSfb;
      lowSubband = h_sbr_cal_env->prevLoSubband;
      highSubband = h_sbr_cal_env->prevHiSubband;

      noSubbands = highSubband - lowSubband;
      ov_highSubband = highSubband;
      if (highSubband < h_sbr_cal_env->prev_ov_highSubband) {
        ov_highSubband = h_sbr_cal_env->prev_ov_highSubband;
      }

      pFreqBandTable[0] = h_sbr_cal_env->prevFreqBandTableLo;
      pFreqBandTable[1] = h_sbr_cal_env->prevFreqBandTableHi;
      pFreqBandTableNoise = h_sbr_cal_env->prevFreqBandTableNoise;
    }

    mapSineFlagsPvc(pFreqBandTable[1], noSubFrameBands[1],
                    h_sbr_cal_env->harmFlagsPrev,
                    h_sbr_cal_env->harmFlagsPrevActive, sineMapped,
                    hFrameData->sinusoidal_position,
                    &h_sbr_cal_env->sinusoidal_positionPrev,
                    (borders[0] > bordersPvc[0]) ? 1 : 0);
  } else {
    mapSineFlags(pFreqBandTable[1], noSubFrameBands[1],
                 hFrameData->addHarmonics, h_sbr_cal_env->harmFlagsPrev,
                 h_sbr_cal_env->harmFlagsPrevActive,
                 hFrameData->frameInfo.tranEnv, sineMapped);
  }

  if (!useLP)
    adj_e = h_sbr_cal_env->filtBufferNoise_e -
            getScalefactor(h_sbr_cal_env->filtBufferNoise, noSubbands) +
            (INT)MAX_SFB_NRG_HEADROOM;

  if (pvc_mode > 0) {
    INT maxSfbNrg_e = pPvcDynamicData->predEsg_expMax;

    maxSfbNrg_e =
        (maxSfbNrg_e + 1) >> 1;

    maxSfbNrg_e += (6 + MAX_SFB_NRG_HEADROOM);

    adj_e = maxSfbNrg_e;
  } else {
    for (i = 0; i < hFrameData->frameInfo.nEnvelopes; i++) {
      INT maxSfbNrg_e =
          -FRACT_BITS + NRG_EXP_OFFSET;

      for (j = noSubFrameBands[hFrameData->frameInfo.freqRes[i]]; j != 0; j--) {
        maxSfbNrg_e = fixMax(maxSfbNrg_e, (INT)((LONG)(*pIenv++) & MASK_E));
      }
      maxSfbNrg_e -= NRG_EXP_OFFSET;

      maxSfbNrg_e =
          (maxSfbNrg_e + 1) >> 1;

      maxSfbNrg_e += (6 + MAX_SFB_NRG_HEADROOM);

      if (borders[i] < hHeaderData->numberTimeSlots)
        adj_e = fMax(maxSfbNrg_e, adj_e);

      if (borders[i + 1] > hHeaderData->numberTimeSlots)
        final_e = fMax(maxSfbNrg_e, final_e);
    }
  }
  pIenv = hFrameData->iEnvelope;

  if (pvc_mode > 0) {
    i = bordersPvc[0];
    i_stop = PVC_NTIMESLOT;
    FDK_ASSERT(bordersPvc[hFrameData->frameInfo.nEnvelopes] == PVC_NTIMESLOT);
  } else {
    i = 0;
    i_stop = hFrameData->frameInfo.nEnvelopes;
  }
  for (; i < i_stop; i++) {
    int k, noNoiseFlag;
    SCHAR noise_e, input_e = SCALE2EXP(sbrScaleFactor->hb_scale);
    C_ALLOC_SCRATCH_START(pNrgs, ENV_CALC_NRGS, 1);

    int start_pos, stop_pos, freq_res;
    if (pvc_mode > 0) {
      start_pos =
          hHeaderData->timeStep *
          i;
      stop_pos = hHeaderData->timeStep * (i + 1);
      freq_res =
          hFrameData->frameInfo
              .freqRes[0];
      FDK_ASSERT(
          freq_res ==
          hFrameData->frameInfo.freqRes[hFrameData->frameInfo.nEnvelopes - 1]);
    } else {
      start_pos = hHeaderData->timeStep *
                  borders[i];
      stop_pos = hHeaderData->timeStep *
                 borders[i + 1];
      freq_res =
          hFrameData->frameInfo
              .freqRes[i];
    }

    FDKmemclear(pNrgs, sizeof(ENV_CALC_NRGS));

    if (pvc_mode > 0) {
      expandPredEsg(pPvcDynamicData, i, (int)MAX_FREQ_COEFFS, pNrgs->nrgRef,
                    pNrgs->nrgRef_e);

      if (i == borders[0]) {
        mapSineFlags(pFreqBandTable[1], noSubFrameBands[1],
                     hFrameData->addHarmonics, h_sbr_cal_env->harmFlagsPrev,
                     h_sbr_cal_env->harmFlagsPrevActive,
                     hFrameData->sinusoidal_position, sineMapped);
      }

      if (i >= hFrameData->frameInfo.bordersNoise[envNoise + 1]) {
        if (envNoise >= 0) {
          noiseLevels += noNoiseBands;
        } else {
          noNoiseBands = hFreq->nNfb;
          noSubFrameBands = hFreq->nSfb;
          noiseLevels = hFrameData->sbrNoiseFloorLevel;

          lowSubband = hFreq->lowSubband;
          highSubband = hFreq->highSubband;

          noSubbands = highSubband - lowSubband;
          ov_highSubband = highSubband;
          if (highSubband < hFreq->ov_highSubband) {
            ov_highSubband = hFreq->ov_highSubband;
          }

          pFreqBandTable[0] = hFreq->freqBandTableLo;
          pFreqBandTable[1] = hFreq->freqBandTableHi;
          pFreqBandTableNoise = hFreq->freqBandTableNoise;
        }
        envNoise++;
      }
    } else {
      if (borders[i] == hFrameData->frameInfo.bordersNoise[envNoise + 1]) {
        noiseLevels += noNoiseBands;
        envNoise++;
      }
    }
    if (i == hFrameData->frameInfo.tranEnv ||
        i == h_sbr_cal_env->prevTranEnv)
    {
      noNoiseFlag = 1;
      if (!useLP) smooth_length = 0;
    } else {
      noNoiseFlag = 0;
      if (!useLP)
        smooth_length = (1 - hHeaderData->bs_data.smoothingLength)
                        << 2;
    }

    if (hHeaderData->bs_data.interpolFreq)
      calcNrgPerSubband(analysBufferReal, (useLP) ? NULL : analysBufferImag,
                        lowSubband, highSubband, start_pos, stop_pos, input_e,
                        pNrgs->nrgEst, pNrgs->nrgEst_e);
    else
      calcNrgPerSfb(analysBufferReal, (useLP) ? NULL : analysBufferImag,
                    noSubFrameBands[freq_res], pFreqBandTable[freq_res],
                    start_pos, stop_pos, input_e, pNrgs->nrgEst,
                    pNrgs->nrgEst_e);

    {
      UCHAR *table = pFreqBandTable[freq_res];
      UCHAR *pUiNoise =
          &pFreqBandTableNoise[1];

      FIXP_SGL *pNoiseLevels = noiseLevels;

      FIXP_DBL tmpNoise =
          FX_SGL2FX_DBL((FIXP_SGL)((LONG)(*pNoiseLevels) & MASK_M));
      SCHAR tmpNoise_e =
          (UCHAR)((LONG)(*pNoiseLevels++) & MASK_E) - NOISE_EXP_OFFSET;

      int cc = 0;
      c = 0;
      if (pvc_mode > 0) {
        for (j = 0; j < noSubFrameBands[freq_res]; j++) {
          UCHAR sinePresentFlag = 0;
          int li = table[j];
          int ui = table[j + 1];

          for (k = li; k < ui; k++) {
            sinePresentFlag |= (i >= sineMapped[cc]);
            cc++;
          }

          for (k = li; k < ui; k++) {
            FIXP_DBL refNrg = pNrgs->nrgRef[k - lowSubband];
            SCHAR refNrg_e = pNrgs->nrgRef_e[k - lowSubband];

            if (k >= *pUiNoise) {
              tmpNoise =
                  FX_SGL2FX_DBL((FIXP_SGL)((LONG)(*pNoiseLevels) & MASK_M));
              tmpNoise_e =
                  (SCHAR)((LONG)(*pNoiseLevels++) & MASK_E) - NOISE_EXP_OFFSET;

              pUiNoise++;
            }

            FDK_ASSERT(k >= lowSubband);

            if (useLP) useAliasReduction[k - lowSubband] = !sinePresentFlag;

            pNrgs->nrgSine[c] = FL2FXCONST_DBL(0.0f);
            pNrgs->nrgSine_e[c] = 0;

            calcSubbandGain(refNrg, refNrg_e, pNrgs, c, tmpNoise, tmpNoise_e,
                            sinePresentFlag, i >= sineMapped[c], noNoiseFlag);

            c++;
          }
        }
      } else {
        for (j = 0; j < noSubFrameBands[freq_res]; j++) {
          FIXP_DBL refNrg = FX_SGL2FX_DBL((FIXP_SGL)((LONG)(*pIenv) & MASK_M));
          SCHAR refNrg_e = (SCHAR)((LONG)(*pIenv) & MASK_E) - NRG_EXP_OFFSET;

          UCHAR sinePresentFlag = 0;
          int li = table[j];
          int ui = table[j + 1];

          for (k = li; k < ui; k++) {
            sinePresentFlag |= (i >= sineMapped[cc]);
            cc++;
          }

          for (k = li; k < ui; k++) {
            if (k >= *pUiNoise) {
              tmpNoise =
                  FX_SGL2FX_DBL((FIXP_SGL)((LONG)(*pNoiseLevels) & MASK_M));
              tmpNoise_e =
                  (SCHAR)((LONG)(*pNoiseLevels++) & MASK_E) - NOISE_EXP_OFFSET;

              pUiNoise++;
            }

            FDK_ASSERT(k >= lowSubband);

            if (useLP) useAliasReduction[k - lowSubband] = !sinePresentFlag;

            pNrgs->nrgSine[c] = FL2FXCONST_DBL(0.0f);
            pNrgs->nrgSine_e[c] = 0;

            calcSubbandGain(refNrg, refNrg_e, pNrgs, c, tmpNoise, tmpNoise_e,
                            sinePresentFlag, i >= sineMapped[c], noNoiseFlag);

            pNrgs->nrgRef[c] = refNrg;
            pNrgs->nrgRef_e[c] = refNrg_e;

            c++;
          }
          pIenv++;
        }
      }
    }

    for (c = 0; c < hFreq->noLimiterBands; c++) {
      FIXP_DBL sumRef, boostGain, maxGain;
      FIXP_DBL accu = FL2FXCONST_DBL(0.0f);
      SCHAR sumRef_e, boostGain_e, maxGain_e, accu_e = 0;
      int maxGainLimGainSum_e = 0;

      calcAvgGain(pNrgs, hFreq->limiterBandTable[c],
                  hFreq->limiterBandTable[c + 1], &sumRef, &sumRef_e, &maxGain,
                  &maxGain_e);

      maxGain = fMult(
          maxGain,
          FDK_sbrDecoder_sbr_limGains_m[hHeaderData->bs_data.limiterGains]);
      maxGainLimGainSum_e =
          maxGain_e +
          FDK_sbrDecoder_sbr_limGains_e[hHeaderData->bs_data.limiterGains];
      maxGain_e =
          (maxGainLimGainSum_e > 127) ? (SCHAR)127 : (SCHAR)maxGainLimGainSum_e;

      if (maxGain == FL2FXCONST_DBL(0.0f))
        maxGain_e = -FRACT_BITS;
      else {
        SCHAR charTemp = CountLeadingBits(maxGain);
        maxGain_e -= charTemp;
        maxGain <<= (int)charTemp;
      }

      if (maxGain_e >= maxGainLimit_e) {
        maxGain = FL2FXCONST_DBL(0.5f);
        maxGain_e = maxGainLimit_e;
      }

      for (k = hFreq->limiterBandTable[c]; k < hFreq->limiterBandTable[c + 1];
           k++) {
        if ((pNrgs->nrgGain_e[k] > maxGain_e) ||
            (pNrgs->nrgGain_e[k] == maxGain_e && pNrgs->nrgGain[k] > maxGain)) {
          FIXP_DBL noiseAmp;
          SCHAR noiseAmp_e;

          FDK_divide_MantExp(maxGain, maxGain_e, pNrgs->nrgGain[k],
                             pNrgs->nrgGain_e[k], &noiseAmp, &noiseAmp_e);
          pNrgs->noiseLevel[k] = fMult(pNrgs->noiseLevel[k], noiseAmp);
          pNrgs->noiseLevel_e[k] += noiseAmp_e;
          pNrgs->nrgGain[k] = maxGain;
          pNrgs->nrgGain_e[k] = maxGain_e;
        }
      }

      for (k = hFreq->limiterBandTable[c]; k < hFreq->limiterBandTable[c + 1];
           k++) {
        FIXP_DBL tmp = fMult(pNrgs->nrgGain[k], pNrgs->nrgEst[k]);
        SCHAR tmp_e = pNrgs->nrgGain_e[k] + pNrgs->nrgEst_e[k];
        FDK_add_MantExp(tmp, tmp_e, accu, accu_e, &accu, &accu_e);

        if (pNrgs->nrgSine[k] != FL2FXCONST_DBL(0.0f)) {
          FDK_add_MantExp(pNrgs->nrgSine[k], pNrgs->nrgSine_e[k], accu, accu_e,
                          &accu, &accu_e);
        } else {
          if (noNoiseFlag == 0) {
            FDK_add_MantExp(pNrgs->noiseLevel[k], pNrgs->noiseLevel_e[k], accu,
                            accu_e, &accu, &accu_e);
          }
        }
      }

      if (accu == (FIXP_DBL)0) {
        boostGain = FL2FXCONST_DBL(0.6279716f);
        boostGain_e = 2;
      } else {
        INT div_e;
        boostGain = fDivNorm(sumRef, accu, &div_e);
        boostGain_e = sumRef_e - accu_e + div_e;
      }

      if ((boostGain_e > 3) ||
          (boostGain_e == 2 && boostGain > FL2FXCONST_DBL(0.6279716f)) ||
          (boostGain_e == 3 && boostGain > FL2FXCONST_DBL(0.3139858f))) {
        boostGain = FL2FXCONST_DBL(0.6279716f);
        boostGain_e = 2;
      }
      for (k = hFreq->limiterBandTable[c]; k < hFreq->limiterBandTable[c + 1];
           k++) {
        pNrgs->nrgGain[k] = fMultDiv2(pNrgs->nrgGain[k], boostGain);
        pNrgs->nrgGain_e[k] = pNrgs->nrgGain_e[k] + boostGain_e + 1;

        pNrgs->nrgSine[k] = fMultDiv2(pNrgs->nrgSine[k], boostGain);
        pNrgs->nrgSine_e[k] = pNrgs->nrgSine_e[k] + boostGain_e + 1;

        pNrgs->noiseLevel[k] = fMultDiv2(pNrgs->noiseLevel[k], boostGain);
        pNrgs->noiseLevel_e[k] = pNrgs->noiseLevel_e[k] + boostGain_e + 1;
      }
    }

    if (useLP)
      aliasingReduction(degreeAlias + lowSubband, pNrgs, useAliasReduction,
                        noSubbands);

    noise_e = (start_pos < no_cols) ? adj_e : final_e;

    if (start_pos >= no_cols) {
      int diff = h_sbr_cal_env->filtBufferNoise_e - noise_e;
      if (diff > 0) {
        int s = getScalefactor(h_sbr_cal_env->filtBufferNoise, noSubbands);
        if (diff > s) {
          final_e += diff - s;
          noise_e = final_e;
        }
      }
    }

    for (k = 0; k < noSubbands; k++) {
      FDK_sqrt_MantExp(&pNrgs->nrgSine[k], &pNrgs->nrgSine_e[k], &noise_e);
      FDK_sqrt_MantExp(&pNrgs->nrgGain[k], &pNrgs->nrgGain_e[k],
                       &pNrgs->nrgGain_e[k]);
      FDK_sqrt_MantExp(&pNrgs->noiseLevel[k], &pNrgs->noiseLevel_e[k],
                       &noise_e);
    }

    {
      int scale_change, sc_change;
      FIXP_SGL smooth_ratio;
      int filtBufferNoiseShift = 0;

      if (h_sbr_cal_env->startUp) {
        if (!useLP) {
          h_sbr_cal_env->filtBufferNoise_e = noise_e;

          FDKmemcpy(h_sbr_cal_env->filtBuffer_e, pNrgs->nrgGain_e,
                    noSubbands * sizeof(SCHAR));
          FDKmemcpy(h_sbr_cal_env->filtBufferNoise, pNrgs->noiseLevel,
                    noSubbands * sizeof(FIXP_DBL));
          FDKmemcpy(h_sbr_cal_env->filtBuffer, pNrgs->nrgGain,
                    noSubbands * sizeof(FIXP_DBL));
        }
        h_sbr_cal_env->startUp = 0;
      }

      if (!useLP) {
        equalizeFiltBufferExp(h_sbr_cal_env->filtBuffer,
                              h_sbr_cal_env->filtBuffer_e,
                              pNrgs->nrgGain,
                              pNrgs->nrgGain_e,
                              noSubbands);

        if ((h_sbr_cal_env->filtBufferNoise_e - noise_e) >= 0) {
          int shift = fixMin(DFRACT_BITS - 1,
                             (int)(h_sbr_cal_env->filtBufferNoise_e - noise_e));
          for (k = 0; k < noSubbands; k++)
            h_sbr_cal_env->filtBufferNoise[k] <<= shift;
        } else {
          int shift =
              fixMin(DFRACT_BITS - 1,
                     -(int)(h_sbr_cal_env->filtBufferNoise_e - noise_e));
          for (k = 0; k < noSubbands; k++)
            h_sbr_cal_env->filtBufferNoise[k] >>= shift;
        }

        h_sbr_cal_env->filtBufferNoise_e = noise_e;
      }

      scale_change = -(DFRACT_BITS - 1);
      for (k = 0; k < noSubbands; k++) {
        scale_change = fixMax(scale_change, (int)pNrgs->nrgGain_e[k]);
      }
      sc_change = (start_pos < no_cols) ? adj_e - input_e : final_e - input_e;

      if ((scale_change - sc_change + 1) < 0)
        scale_change -= (scale_change - sc_change + 1);

      scale_change = (scale_change - sc_change) + 1;

      for (k = 0; k < noSubbands; k++) {
        int sc = scale_change - pNrgs->nrgGain_e[k] + (sc_change - 1);
        pNrgs->nrgGain[k] >>= fixMin(sc, DFRACT_BITS - 1);
        pNrgs->nrgGain_e[k] += sc;
      }

      if (!useLP) {
        for (k = 0; k < noSubbands; k++) {
          int sc =
              scale_change - h_sbr_cal_env->filtBuffer_e[k] + (sc_change - 1);
          h_sbr_cal_env->filtBuffer[k] >>= fixMin(sc, DFRACT_BITS - 1);
        }
      }

      for (j = start_pos; j < stop_pos; j++) {
        if ((j == no_cols) && (start_pos < no_cols)) {
          int shift = (int)(noise_e - final_e);
          if (!useLP)
            filtBufferNoiseShift = shift;
          if (shift >= 0) {
            shift = fixMin(DFRACT_BITS - 1, shift);
            for (k = 0; k < noSubbands; k++) {
              pNrgs->nrgSine[k] <<= shift;
              pNrgs->noiseLevel[k] <<= shift;
            }
          } else {
            shift = fixMin(DFRACT_BITS - 1, -shift);
            for (k = 0; k < noSubbands; k++) {
              pNrgs->nrgSine[k] >>= shift;
              pNrgs->noiseLevel[k] >>= shift;
            }
          }

          noise_e = final_e;
          if (!useLP)
            h_sbr_cal_env->filtBufferNoise_e =
                noise_e;

          sc_change -= (final_e - input_e);

          if (sc_change < 0) {
            for (k = 0; k < noSubbands; k++) {
              pNrgs->nrgGain[k] >>= -sc_change;
              pNrgs->nrgGain_e[k] += -sc_change;
            }
            if (!useLP) {
              for (k = 0; k < noSubbands; k++) {
                h_sbr_cal_env->filtBuffer[k] >>= -sc_change;
              }
            }
          } else {
            scale_change += sc_change;
          }
        }

        if (!useLP) {
          if (j - start_pos < smooth_length)
            smooth_ratio = FDK_sbrDecoder_sbr_smoothFilter[j - start_pos];
          else
            smooth_ratio = FL2FXCONST_SGL(0.0f);

          if (iTES_enable) {
            adjustTimeSlotHQ_GainAndNoise(
                &analysBufferReal[j][lowSubband],
                &analysBufferImag[j][lowSubband], h_sbr_cal_env, pNrgs,
                lowSubband, noSubbands, fMin(scale_change, DFRACT_BITS - 1),
                smooth_ratio, noNoiseFlag, filtBufferNoiseShift);
          } else {
            adjustTimeSlotHQ(&analysBufferReal[j][lowSubband],
                             &analysBufferImag[j][lowSubband], h_sbr_cal_env,
                             pNrgs, lowSubband, noSubbands,
                             fMin(scale_change, DFRACT_BITS - 1), smooth_ratio,
                             noNoiseFlag, filtBufferNoiseShift);
          }
        } else {
          FDK_ASSERT(!iTES_enable);
          if (flags & SBRDEC_ELD_GRID) {
            adjustTimeSlot_EldGrid(
                &analysBufferReal[j][lowSubband], pNrgs,
                &h_sbr_cal_env->harmIndex, lowSubband, noSubbands,
                fMin(scale_change, DFRACT_BITS - 1), noNoiseFlag,
                &h_sbr_cal_env->phaseIndex,
                fMax(EXP2SCALE(adj_e) - sbrScaleFactor->lb_scale,
                     -(DFRACT_BITS - 1)));
          } else {
            adjustTimeSlotLC(&analysBufferReal[j][lowSubband], pNrgs,
                             &h_sbr_cal_env->harmIndex, lowSubband, noSubbands,
                             fMin(scale_change, DFRACT_BITS - 1), noNoiseFlag,
                             &h_sbr_cal_env->phaseIndex);
          }
        }
        pNrgs->exponent[(j < no_cols) ? 0 : 1] =
            (SCHAR)((15 - sbrScaleFactor->hb_scale) + pNrgs->nrgGain_e[0] + 1 -
                    scale_change);
      }

      if (iTES_enable) {
        apply_inter_tes(
            analysBufferReal,
            analysBufferImag,
            sbrScaleFactor, pNrgs->exponent, hHeaderData->timeStep, start_pos,
            stop_pos, lowSubband, noSubbands,
            hFrameData
                ->interTempShapeMode[i]
        );

        for (j = start_pos; j < stop_pos; j++) {
          scale_change = 0;

          if ((start_pos <= no_cols) && (stop_pos > no_cols)) {
            scale_change = pNrgs->exponent[1] - pNrgs->exponent[0];
          }

          adjustTimeSlotHQ_AddHarmonics(
              &analysBufferReal[j][lowSubband],
              &analysBufferImag[j][lowSubband], h_sbr_cal_env, pNrgs,
              lowSubband, noSubbands,
              -iTES_scale_change + ((j < no_cols) ? scale_change : 0));
        }
      }

      if (!useLP) {
        FDKmemcpy(h_sbr_cal_env->filtBuffer, pNrgs->nrgGain,
                  noSubbands * sizeof(FIXP_DBL));
        FDKmemcpy(h_sbr_cal_env->filtBuffer_e, pNrgs->nrgGain_e,
                  noSubbands * sizeof(SCHAR));
        FDKmemcpy(h_sbr_cal_env->filtBufferNoise, pNrgs->noiseLevel,
                  noSubbands * sizeof(FIXP_DBL));
      }
    }
    C_ALLOC_SCRATCH_END(pNrgs, ENV_CALC_NRGS, 1);
  }

  adj_e += iTES_scale_change;

  {
    FIXP_DBL maxVal;
    int ov_reserve, reserve;

    maxVal =
        maxSubbandSample(analysBufferReal, (useLP) ? NULL : analysBufferImag,
                         lowSubband, ov_highSubband, 0, first_start);

    ov_reserve = fNorm(maxVal);

    maxVal =
        maxSubbandSample(analysBufferReal, (useLP) ? NULL : analysBufferImag,
                         lowSubband, highSubband, first_start, no_cols);

    reserve = fNorm(maxVal);

    output_e = fMax(ov_adj_e - ov_reserve, adj_e - reserve);

    rescaleSubbandSamples(analysBufferReal, (useLP) ? NULL : analysBufferImag,
                          lowSubband, ov_highSubband, 0, first_start,
                          ov_adj_e - output_e);

    rescaleSubbandSamples(analysBufferReal, (useLP) ? NULL : analysBufferImag,
                          lowSubband, highSubband, first_start, no_cols,
                          adj_e - output_e);
  }

  sbrScaleFactor->hb_scale = EXP2SCALE(output_e);

  sbrScaleFactor->ov_hb_scale = EXP2SCALE(final_e + iTES_scale_change);

  if (hFrameData->frameInfo.tranEnv == hFrameData->frameInfo.nEnvelopes)
    h_sbr_cal_env->prevTranEnv = 0;
  else
    h_sbr_cal_env->prevTranEnv = -1;

  if (pvc_mode > 0) {
    FDK_ASSERT(hFrameData->frameInfo
                   .bordersNoise[hFrameData->frameInfo.nNoiseEnvelopes - 1] <
               PVC_NTIMESLOT);
    if (hFrameData->frameInfo
            .bordersNoise[hFrameData->frameInfo.nNoiseEnvelopes] >
        PVC_NTIMESLOT) {
      FDK_ASSERT(noiseLevels ==
                 (hFrameData->sbrNoiseFloorLevel +
                  (hFrameData->frameInfo.nNoiseEnvelopes - 1) * noNoiseBands));
      h_sbr_cal_env->prevNNfb = noNoiseBands;

      h_sbr_cal_env->prevNSfb[0] = noSubFrameBands[0];
      h_sbr_cal_env->prevNSfb[1] = noSubFrameBands[1];

      h_sbr_cal_env->prevLoSubband = lowSubband;
      h_sbr_cal_env->prevHiSubband = highSubband;
      h_sbr_cal_env->prev_ov_highSubband = ov_highSubband;

      FDKmemcpy(h_sbr_cal_env->prevFreqBandTableLo, pFreqBandTable[0],
                noSubFrameBands[0] + 1);
      FDKmemcpy(h_sbr_cal_env->prevFreqBandTableHi, pFreqBandTable[1],
                noSubFrameBands[1] + 1);
      FDKmemcpy(h_sbr_cal_env->prevFreqBandTableNoise,
                hFreq->freqBandTableNoise, sizeof(hFreq->freqBandTableNoise));

      FDKmemcpy(h_sbr_cal_env->prevSbrNoiseFloorLevel, noiseLevels,
                MAX_NOISE_COEFFS * sizeof(FIXP_SGL));
    }
  }

  C_ALLOC_SCRATCH_END(useAliasReduction, UCHAR, 64)
}

SBR_ERROR
createSbrEnvelopeCalc(
    HANDLE_SBR_CALCULATE_ENVELOPE hs,
    HANDLE_SBR_HEADER_DATA
        hHeaderData,
    const int chan,
    const UINT flags) {
  SBR_ERROR err = SBRDEC_OK;
  int i;

  for (i = 0; i < ADD_HARMONICS_FLAGS_SIZE; i++) {
    hs->harmFlagsPrev[i] = 0;
    hs->harmFlagsPrevActive[i] = 0;
  }
  hs->harmIndex = 0;

  FDKmemclear(hs->prevSbrNoiseFloorLevel, sizeof(hs->prevSbrNoiseFloorLevel));
  hs->prevNNfb = 0;
  FDKmemclear(hs->prevFreqBandTableNoise, sizeof(hs->prevFreqBandTableNoise));
  hs->sinusoidal_positionPrev = 0;

  hs->prevTranEnv = -1;

  resetSbrEnvelopeCalc(hs);

  if (chan == 0) {
    err = resetFreqBandTables(hHeaderData, flags);
  }

  return err;
}

int deleteSbrEnvelopeCalc(HANDLE_SBR_CALCULATE_ENVELOPE hs) { return 0; }

void resetSbrEnvelopeCalc(
    HANDLE_SBR_CALCULATE_ENVELOPE hCalEnv)
{
  hCalEnv->phaseIndex = 0;

  hCalEnv->filtBufferNoise_e = 0;

  hCalEnv->startUp = 1;
}

static void equalizeFiltBufferExp(
    FIXP_DBL *filtBuffer,
    SCHAR *filtBuffer_e,
    FIXP_DBL *nrgGain,
    SCHAR *nrgGain_e,
    int subbands)
{
  int band;
  int diff;

  for (band = 0; band < subbands; band++) {
    diff = (int)(nrgGain_e[band] - filtBuffer_e[band]);
    if (diff > 0) {
      filtBuffer[band] >>=
          fMin(diff, DFRACT_BITS - 1);
      filtBuffer_e[band] += diff;
    } else if (diff < 0) {
      int reserve = CntLeadingZeros(fixp_abs(filtBuffer[band])) - 1;

      if ((-diff) <= reserve) {
        filtBuffer[band] <<= (-diff);
        filtBuffer_e[band] += diff;
      } else {
        filtBuffer[band] <<=
            reserve;
        filtBuffer_e[band] -= reserve;

        diff = -(reserve + diff);
        nrgGain[band] >>= fMin(diff, DFRACT_BITS - 1);
        nrgGain_e[band] += diff;
      }
    }
  }
}

void rescaleSubbandSamples(
    FIXP_DBL **re,
    FIXP_DBL **im,
    int lowSubband,
    int highSubband,
    int start_pos,
    int next_pos,
    int shift)
{
  int width = highSubband - lowSubband;

  if ((width > 0) && (shift != 0)) {
    if (im != NULL) {
      for (int l = start_pos; l < next_pos; l++) {
        scaleValues(&re[l][lowSubband], width, shift);
        scaleValues(&im[l][lowSubband], width, shift);
      }
    } else {
      for (int l = start_pos; l < next_pos; l++) {
        scaleValues(&re[l][lowSubband], width, shift);
      }
    }
  }
}

static inline FIXP_DBL FDK_get_maxval_real(FIXP_DBL maxVal, FIXP_DBL *reTmp,
                                           INT width) {
  maxVal = (FIXP_DBL)0;
  while (width-- != 0) {
    FIXP_DBL tmp = *(reTmp++);
    maxVal |= (FIXP_DBL)((LONG)(tmp) ^ ((LONG)tmp >> (DFRACT_BITS - 1)));
  }

  return maxVal;
}

FIXP_DBL maxSubbandSample(
    FIXP_DBL **re,
    FIXP_DBL **im,
    int lowSubband,
    int highSubband,
    int start_pos,
    int next_pos
) {
  FIXP_DBL maxVal = FL2FX_DBL(0.0f);
  unsigned int width = highSubband - lowSubband;

  FDK_ASSERT(width <= (64));

  if (width > 0) {
    if (im != NULL) {
      for (int l = start_pos; l < next_pos; l++) {
        int k = width;
        FIXP_DBL *reTmp = &re[l][lowSubband];
        FIXP_DBL *imTmp = &im[l][lowSubband];
        do {
          FIXP_DBL tmp1 = *(reTmp++);
          FIXP_DBL tmp2 = *(imTmp++);
          maxVal |=
              (FIXP_DBL)((LONG)(tmp1) ^ ((LONG)tmp1 >> (DFRACT_BITS - 1)));
          maxVal |=
              (FIXP_DBL)((LONG)(tmp2) ^ ((LONG)tmp2 >> (DFRACT_BITS - 1)));
        } while (--k != 0);
      }
    } else {
      for (int l = start_pos; l < next_pos; l++) {
        maxVal |= FDK_get_maxval_real(maxVal, &re[l][lowSubband], width);
      }
    }
  }

  if (maxVal > (FIXP_DBL)0) {
    FIXP_DBL lowerPow2 =
        (FIXP_DBL)(1 << (DFRACT_BITS - 1 - CntLeadingZeros(maxVal)));
    if (maxVal == lowerPow2) maxVal += (FIXP_DBL)1;
  }

  return (maxVal);
}

#define SHIFT_BEFORE_SQUARE (4)

static void calcNrgPerSubband(
    FIXP_DBL **analysBufferReal,
    FIXP_DBL **analysBufferImag,
    int lowSubband,
    int highSubband,
    int start_pos,
    int next_pos,
    SCHAR frameExp,
    FIXP_DBL *nrgEst,
    SCHAR *nrgEst_e)
{
  FIXP_SGL invWidth;
  SCHAR preShift;
  SCHAR shift;
  FIXP_DBL sum;
  int k;

  invWidth = FX_DBL2FX_SGL(GetInvInt(next_pos - start_pos));
  frameExp = frameExp << 1;

  for (k = lowSubband; k < highSubband; k++) {
    FIXP_DBL bufferReal[(((1024) / (32) * (4) / 2) + (3 * (4)))];
    FIXP_DBL bufferImag[(((1024) / (32) * (4) / 2) + (3 * (4)))];
    FIXP_DBL maxVal;

    if (analysBufferImag != NULL) {
      int l;
      maxVal = FL2FX_DBL(0.0f);
      for (l = start_pos; l < next_pos; l++) {
        bufferImag[l] = analysBufferImag[l][k];
        maxVal |= (FIXP_DBL)((LONG)(bufferImag[l]) ^
                             ((LONG)bufferImag[l] >> (DFRACT_BITS - 1)));
        bufferReal[l] = analysBufferReal[l][k];
        maxVal |= (FIXP_DBL)((LONG)(bufferReal[l]) ^
                             ((LONG)bufferReal[l] >> (DFRACT_BITS - 1)));
      }
    } else {
      int l;
      maxVal = FL2FX_DBL(0.0f);
      for (l = start_pos; l < next_pos; l++) {
        bufferReal[l] = analysBufferReal[l][k];
        maxVal |= (FIXP_DBL)((LONG)(bufferReal[l]) ^
                             ((LONG)bufferReal[l] >> (DFRACT_BITS - 1)));
      }
    }

    if (maxVal != FL2FXCONST_DBL(0.f)) {
      FIXP_DBL accu;
      preShift = CntLeadingZeros(maxVal) - 1;
      preShift -= SHIFT_BEFORE_SQUARE;

      preShift = fMin(preShift, (SCHAR)25);

      accu = FL2FXCONST_DBL(0.0f);
      if (preShift >= 0) {
        int l;
        if (analysBufferImag != NULL) {
          for (l = start_pos; l < next_pos; l++) {
            FIXP_DBL temp1 = bufferReal[l] << (int)preShift;
            FIXP_DBL temp2 = bufferImag[l] << (int)preShift;
            accu = fPow2AddDiv2(accu, temp1);
            accu = fPow2AddDiv2(accu, temp2);
          }
        } else {
          for (l = start_pos; l < next_pos; l++) {
            FIXP_DBL temp = bufferReal[l] << (int)preShift;
            accu = fPow2AddDiv2(accu, temp);
          }
        }
      } else {
        int l;
        int negpreShift = -preShift;
        if (analysBufferImag != NULL) {
          for (l = start_pos; l < next_pos; l++) {
            FIXP_DBL temp1 = bufferReal[l] >> (int)negpreShift;
            FIXP_DBL temp2 = bufferImag[l] >> (int)negpreShift;
            accu = fPow2AddDiv2(accu, temp1);
            accu = fPow2AddDiv2(accu, temp2);
          }
        } else {
          for (l = start_pos; l < next_pos; l++) {
            FIXP_DBL temp = bufferReal[l] >> (int)negpreShift;
            accu = fPow2AddDiv2(accu, temp);
          }
        }
      }
      accu <<= 1;

      shift = fNorm(accu);
      sum = accu << (int)shift;

      *nrgEst++ = fMult(sum, invWidth);
      shift += 2 * preShift;
      if (analysBufferImag != NULL)
        *nrgEst_e++ = frameExp - shift;
      else
        *nrgEst_e++ = frameExp - shift + 1;
    }
    else {
      *nrgEst++ = FL2FXCONST_DBL(0.0f);
      *nrgEst_e++ = 0;
    }
  }
}

static void calcNrgPerSfb(
    FIXP_DBL **analysBufferReal,
    FIXP_DBL **analysBufferImag,
    int nSfb,
    UCHAR *freqBandTable,
    int start_pos,
    int next_pos,
    SCHAR input_e,
    FIXP_DBL *nrgEst,
    SCHAR *nrgEst_e)
{
  FIXP_SGL invWidth;
  FIXP_DBL temp;
  SCHAR preShift;
  SCHAR shift, sum_e;
  FIXP_DBL sum;

  int j, k, l, li, ui;
  FIXP_DBL sumAll, sumLine;

  invWidth = FX_DBL2FX_SGL(GetInvInt(next_pos - start_pos));
  input_e = input_e << 1;

  for (j = 0; j < nSfb; j++) {
    li = freqBandTable[j];
    ui = freqBandTable[j + 1];

    FIXP_DBL maxVal = maxSubbandSample(analysBufferReal, analysBufferImag, li,
                                       ui, start_pos, next_pos);

    if (maxVal != FL2FXCONST_DBL(0.f)) {
      preShift = CntLeadingZeros(maxVal) - 1;

      preShift -= SHIFT_BEFORE_SQUARE;

      sumAll = FL2FXCONST_DBL(0.0f);

      for (k = li; k < ui; k++) {
        sumLine = FL2FXCONST_DBL(0.0f);

        if (analysBufferImag != NULL) {
          if (preShift >= 0) {
            for (l = start_pos; l < next_pos; l++) {
              temp = analysBufferReal[l][k] << (int)preShift;
              sumLine += fPow2Div2(temp);
              temp = analysBufferImag[l][k] << (int)preShift;
              sumLine += fPow2Div2(temp);
            }
          } else {
            for (l = start_pos; l < next_pos; l++) {
              temp = analysBufferReal[l][k] >> -(int)preShift;
              sumLine += fPow2Div2(temp);
              temp = analysBufferImag[l][k] >> -(int)preShift;
              sumLine += fPow2Div2(temp);
            }
          }
        } else {
          if (preShift >= 0) {
            for (l = start_pos; l < next_pos; l++) {
              temp = analysBufferReal[l][k] << (int)preShift;
              sumLine += fPow2Div2(temp);
            }
          } else {
            for (l = start_pos; l < next_pos; l++) {
              temp = analysBufferReal[l][k] >> -(int)preShift;
              sumLine += fPow2Div2(temp);
            }
          }
        }

        sumLine = sumLine >> (4 - 1);
        sumAll += sumLine;
      }

      shift = fNorm(sumAll);
      sum = sumAll << (int)shift;

      sum = fMult(sum, invWidth);

      sum = fMult(sum, FX_DBL2FX_SGL(GetInvInt(ui - li)));

      if (analysBufferImag != NULL)
        sum_e = input_e + 4 - shift;
      else
        sum_e = input_e + 4 + 1 -
                shift;

      sum_e -= 2 * preShift;
    }
    else {
      sum = FL2FXCONST_DBL(0.0f);
      sum_e = 0;
    }

    for (k = li; k < ui; k++) {
      *nrgEst++ = sum;
      *nrgEst_e++ = sum_e;
    }
  }
}

static void calcSubbandGain(
    FIXP_DBL nrgRef,
    SCHAR
        nrgRef_e,
    ENV_CALC_NRGS *nrgs, int i, FIXP_DBL tmpNoise,
    SCHAR tmpNoise_e,
    UCHAR sinePresentFlag,
    UCHAR sineMapped,
    int noNoiseFlag)
{
  FIXP_DBL nrgEst = nrgs->nrgEst[i];
  SCHAR nrgEst_e =
      nrgs->nrgEst_e[i];
  FIXP_DBL *ptrNrgGain = &nrgs->nrgGain[i];
  SCHAR *ptrNrgGain_e =
      &nrgs->nrgGain_e[i];
  FIXP_DBL *ptrNoiseLevel =
      &nrgs->noiseLevel[i];
  SCHAR *ptrNoiseLevel_e =
      &nrgs->noiseLevel_e[i];
  FIXP_DBL *ptrNrgSine = &nrgs->nrgSine[i];
  SCHAR *ptrNrgSine_e =
      &nrgs->nrgSine_e[i];

  FIXP_DBL a, b, c;
  SCHAR a_e, b_e, c_e;

  b_e = (int)(nrgEst_e - 1);
  if (b_e >= 0) {
    nrgEst = (FL2FXCONST_DBL(0.5f) >> (INT)fixMin(b_e + 1, DFRACT_BITS - 1)) +
             (nrgEst >> 1);
    nrgEst_e += 1;

  } else {
    nrgEst = (nrgEst >> (INT)(fixMin(-b_e + 1, DFRACT_BITS - 1))) +
             (FL2FXCONST_DBL(0.5f) >> 1);
    nrgEst_e = 2;
  }

  a = fMult(nrgRef, tmpNoise);
  a_e = nrgRef_e + tmpNoise_e;

  b_e = (int)(tmpNoise_e - 1);
  if (b_e >= 0) {
    b = (FL2FXCONST_DBL(0.5f) >> (INT)fixMin(b_e + 1, DFRACT_BITS - 1)) +
        (tmpNoise >> 1);
    b_e = tmpNoise_e + 1;
  } else {
    b = (tmpNoise >> (INT)(fixMin(-b_e + 1, DFRACT_BITS - 1))) +
        (FL2FXCONST_DBL(0.5f) >> 1);
    b_e = 2;
  }

  FDK_divide_MantExp(a, a_e, b, b_e, ptrNoiseLevel, ptrNoiseLevel_e);

  if (sinePresentFlag) {
    c = fMult(b, nrgEst);
    c_e = b_e + nrgEst_e;

    FDK_divide_MantExp(a, a_e, c, c_e, ptrNrgGain, ptrNrgGain_e);

    if (sineMapped) {
      FDK_divide_MantExp(nrgRef, nrgRef_e, b, b_e, ptrNrgSine, ptrNrgSine_e);
    }
  } else {
    if (noNoiseFlag) {
      b = nrgEst;
      b_e = nrgEst_e;
    } else {
      b = fMult(b, nrgEst);
      b_e = b_e + nrgEst_e;
    }

    INT result_exp = 0;
    *ptrNrgGain = fDivNorm(nrgRef, b, &result_exp);
    *ptrNrgGain_e = (SCHAR)result_exp + (nrgRef_e - b_e);

    int headroom = CountLeadingBits(*ptrNrgGain);
    *ptrNrgGain <<= headroom;
    *ptrNrgGain_e -= headroom;
  }
}

static void calcAvgGain(
    ENV_CALC_NRGS *nrgs, int lowSubband,
    int highSubband,
    FIXP_DBL *ptrSumRef, SCHAR *ptrSumRef_e,
    FIXP_DBL *ptrAvgGain,
    SCHAR *ptrAvgGain_e)
{
  FIXP_DBL *nrgRef =
      nrgs->nrgRef;
  SCHAR *nrgRef_e =
      nrgs->nrgRef_e;
  FIXP_DBL *nrgEst = nrgs->nrgEst;
  SCHAR *nrgEst_e =
      nrgs->nrgEst_e;

  FIXP_DBL sumRef = 1;
  FIXP_DBL sumEst = 1;
  SCHAR sumRef_e = -FRACT_BITS;
  SCHAR sumEst_e = -FRACT_BITS;
  int k;

  for (k = lowSubband; k < highSubband; k++) {
    FDK_add_MantExp(sumRef, sumRef_e, nrgRef[k], nrgRef_e[k], &sumRef,
                    &sumRef_e);

    FDK_add_MantExp(sumEst, sumEst_e, nrgEst[k], nrgEst_e[k], &sumEst,
                    &sumEst_e);
  }

  FDK_divide_MantExp(sumRef, sumRef_e, sumEst, sumEst_e, ptrAvgGain,
                     ptrAvgGain_e);

  *ptrSumRef = sumRef;
  *ptrSumRef_e = sumRef_e;
}

static void adjustTimeSlot_EldGrid(
    FIXP_DBL *RESTRICT
        ptrReal,
    ENV_CALC_NRGS *nrgs, UCHAR *ptrHarmIndex,
    int lowSubband,
    int noSubbands,
    int scale_change,
    int noNoiseFlag,
    int *ptrPhaseIndex,
    int scale_diff_low)

{
  int k;
  FIXP_DBL signalReal, sbNoise;
  int tone_count = 0;

  FIXP_DBL *pGain = nrgs->nrgGain;
  FIXP_DBL *RESTRICT pNoiseLevel =
      nrgs->noiseLevel;
  FIXP_DBL *RESTRICT pSineLevel = nrgs->nrgSine;

  int phaseIndex = *ptrPhaseIndex;
  UCHAR harmIndex = *ptrHarmIndex;

  static const INT harmonicPhase[4][2] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};

  static const FIXP_DBL harmonicPhaseX[4][2] = {
      {FL2FXCONST_DBL(2.0 * 1.245183154539139e-001),
       FL2FXCONST_DBL(2.0 * 1.245183154539139e-001)},
      {FL2FXCONST_DBL(2.0 * -1.123767859325028e-001),
       FL2FXCONST_DBL(2.0 * 1.123767859325028e-001)},
      {FL2FXCONST_DBL(2.0 * -1.245183154539139e-001),
       FL2FXCONST_DBL(2.0 * -1.245183154539139e-001)},
      {FL2FXCONST_DBL(2.0 * 1.123767859325028e-001),
       FL2FXCONST_DBL(2.0 * -1.123767859325028e-001)}};

  const FIXP_DBL *p_harmonicPhaseX = &harmonicPhaseX[harmIndex][0];
  const INT *p_harmonicPhase = &harmonicPhase[harmIndex][0];

  const FIXP_DBL max_val = MAX_VAL_NRG_HEADROOM >> scale_change;
  const FIXP_DBL min_val = -max_val;

  *(ptrReal - 1) = fAddSaturate(
      *(ptrReal - 1),
      SATURATE_SHIFT(fMultDiv2(p_harmonicPhaseX[lowSubband & 1], pSineLevel[0]),
                     scale_diff_low, DFRACT_BITS));
  FIXP_DBL pSineLevel_prev = (FIXP_DBL)0;

  int idx_k = lowSubband & 1;

  for (k = 0; k < noSubbands; k++) {
    FIXP_DBL sineLevel_curr = *pSineLevel++;
    phaseIndex = (phaseIndex + 1) & (SBR_NF_NO_RANDOM_VAL - 1);

    signalReal = fMax(fMin(fMultDiv2(*ptrReal, *pGain++), max_val), min_val)
                 << scale_change;
    sbNoise = *pNoiseLevel++;
    if (((INT)sineLevel_curr | noNoiseFlag) == 0) {
      signalReal +=
          fMult(FDK_sbrDecoder_sbr_randomPhase[phaseIndex][0], sbNoise);
    }
    signalReal += sineLevel_curr * p_harmonicPhase[0];
    signalReal =
        fMultAddDiv2(signalReal, pSineLevel_prev, p_harmonicPhaseX[idx_k]);
    pSineLevel_prev = sineLevel_curr;
    idx_k = !idx_k;
    if (k < noSubbands - 1) {
      signalReal =
          fMultAddDiv2(signalReal, pSineLevel[0], p_harmonicPhaseX[idx_k]);
    } else
    {
      if (k + lowSubband + 1 < 63) {
        *(ptrReal + 1) += fMultDiv2(pSineLevel_prev, p_harmonicPhaseX[idx_k]);
      }
    }
    *ptrReal++ = signalReal;

    if (pSineLevel_prev != FL2FXCONST_DBL(0.0f)) {
      if (++tone_count == 16) {
        k++;
        break;
      }
    }
  }
  for (; k < noSubbands; k++) {
    FIXP_DBL sineLevel_curr = *pSineLevel++;
    phaseIndex = (phaseIndex + 1) & (SBR_NF_NO_RANDOM_VAL - 1);

    signalReal = fMax(fMin(fMultDiv2(*ptrReal, *pGain++), max_val), min_val)
                 << scale_change;
    sbNoise = *pNoiseLevel++;
    if (((INT)sineLevel_curr | noNoiseFlag) == 0) {
      signalReal +=
          fMult(FDK_sbrDecoder_sbr_randomPhase[phaseIndex][0], sbNoise);
    }
    signalReal += sineLevel_curr * p_harmonicPhase[0];
    *ptrReal++ = signalReal;
  }

  *ptrHarmIndex = (harmIndex + 1) & 3;
  *ptrPhaseIndex = phaseIndex & (SBR_NF_NO_RANDOM_VAL - 1);
}

static void adjustTimeSlotLC(
    FIXP_DBL *ptrReal,
    ENV_CALC_NRGS *nrgs, UCHAR *ptrHarmIndex,
    int lowSubband,
    int noSubbands,
    int scale_change,
    int noNoiseFlag,
    int *ptrPhaseIndex)
{
  FIXP_DBL *pGain = nrgs->nrgGain;
  FIXP_DBL *pNoiseLevel =
      nrgs->noiseLevel;
  FIXP_DBL *pSineLevel = nrgs->nrgSine;

  int k;
  int index = *ptrPhaseIndex;
  UCHAR harmIndex = *ptrHarmIndex;
  UCHAR freqInvFlag = (lowSubband & 1);
  FIXP_DBL signalReal, sineLevel, sineLevelNext, sineLevelPrev;
  int tone_count = 0;
  int sineSign = 1;
  const FIXP_DBL max_val = MAX_VAL_NRG_HEADROOM >> scale_change;
  const FIXP_DBL min_val = -max_val;

#define C1 ((FIXP_SGL)FL2FXCONST_SGL(2.f * 0.00815f))
#define C1_CLDFB ((FIXP_SGL)FL2FXCONST_SGL(2.f * 0.16773f))

  index = (index + 1) & (SBR_NF_NO_RANDOM_VAL - 1);

  signalReal = fMax(fMin(fMultDiv2(*ptrReal, *pGain++), max_val), min_val)
               << scale_change;
  sineLevel = *pSineLevel++;
  sineLevelNext = (noSubbands > 1) ? pSineLevel[0] : FL2FXCONST_DBL(0.0f);

  if (sineLevel != FL2FXCONST_DBL(0.0f))
    tone_count++;
  else if (!noNoiseFlag)
    signalReal +=
        fMult(FDK_sbrDecoder_sbr_randomPhase[index][0], pNoiseLevel[0]);

  {
    if (!(harmIndex & 0x1)) {
      signalReal += (harmIndex & 0x2) ? -sineLevel : sineLevel;
      *ptrReal++ = signalReal;
    } else {
      int shift = (int)(scale_change + 1);
      shift = (shift >= 0) ? fixMin(DFRACT_BITS - 1, shift)
                           : fixMax(-(DFRACT_BITS - 1), shift);

      FIXP_DBL tmp1 = (shift >= 0) ? (fMultDiv2(C1, sineLevel) >> shift)
                                   : (fMultDiv2(C1, sineLevel) << (-shift));
      FIXP_DBL tmp2 = fMultDiv2(C1, sineLevelNext);

      if (((harmIndex >> 1) & 0x1) ^ freqInvFlag) {
        *(ptrReal - 1) = fAddSaturate(*(ptrReal - 1), tmp1);
        signalReal -= tmp2;
      } else {
        *(ptrReal - 1) = fAddSaturate(*(ptrReal - 1), -tmp1);
        signalReal += tmp2;
      }
      *ptrReal++ = signalReal;
      freqInvFlag = !freqInvFlag;
    }
  }

  pNoiseLevel++;

  if (noSubbands > 2) {
    if (!(harmIndex & 0x1)) {
      if (!harmIndex) {
        sineSign = 0;
      }

      for (k = noSubbands - 2; k != 0; k--) {
        FIXP_DBL sinelevel = *pSineLevel++;
        index++;
        if (((signalReal = (sineSign ? -sinelevel : sinelevel)) ==
             FL2FXCONST_DBL(0.0f)) &&
            !noNoiseFlag) {
          index &= (SBR_NF_NO_RANDOM_VAL - 1);
          signalReal +=
              fMult(FDK_sbrDecoder_sbr_randomPhase[index][0], pNoiseLevel[0]);
        }

        signalReal +=
            fMax(fMin(fMultDiv2(*ptrReal, *pGain++), max_val), min_val)
            << scale_change;

        pNoiseLevel++;
        *ptrReal++ = signalReal;
      }
    } else {
      if (harmIndex == 1) freqInvFlag = !freqInvFlag;

      for (k = noSubbands - 2; k != 0; k--) {
        index++;
        signalReal = fMax(fMin(fMultDiv2(*ptrReal, *pGain++), max_val), min_val)
                     << scale_change;

        if (*pSineLevel++ != FL2FXCONST_DBL(0.0f))
          tone_count++;
        else if (!noNoiseFlag) {
          index &= (SBR_NF_NO_RANDOM_VAL - 1);
          signalReal +=
              fMult(FDK_sbrDecoder_sbr_randomPhase[index][0], pNoiseLevel[0]);
        }

        pNoiseLevel++;

        if (tone_count <= 16) {
          FIXP_DBL addSine = fMultDiv2((pSineLevel[-2] - pSineLevel[0]), C1);
          signalReal += (freqInvFlag) ? (-addSine) : (addSine);
        }

        *ptrReal++ = signalReal;
        freqInvFlag = !freqInvFlag;
      }
    }
  }

  if (noSubbands > -1) {
    index++;
    signalReal = fMax(fMin(fMultDiv2(*ptrReal, *pGain), max_val), min_val)
                 << scale_change;
    sineLevelPrev = fMultDiv2(pSineLevel[-1], FL2FX_SGL(0.0163f));
    sineLevel = pSineLevel[0];

    if (pSineLevel[0] != FL2FXCONST_DBL(0.0f))
      tone_count++;
    else if (!noNoiseFlag) {
      index &= (SBR_NF_NO_RANDOM_VAL - 1);
      signalReal = signalReal + fMult(FDK_sbrDecoder_sbr_randomPhase[index][0],
                                      pNoiseLevel[0]);
    }

    if (!(harmIndex & 0x1)) {
      *ptrReal = signalReal + ((sineSign) ? -sineLevel : sineLevel);
    } else {
      if (tone_count <= 16) {
        if (freqInvFlag) {
          *ptrReal++ = signalReal - sineLevelPrev;
          if (noSubbands + lowSubband < 63)
            *ptrReal = *ptrReal + fMultDiv2(C1, sineLevel);
        } else {
          *ptrReal++ = signalReal + sineLevelPrev;
          if (noSubbands + lowSubband < 63)
            *ptrReal = *ptrReal - fMultDiv2(C1, sineLevel);
        }
      } else
        *ptrReal = signalReal;
    }
  }
  *ptrHarmIndex = (harmIndex + 1) & 3;
  *ptrPhaseIndex = index & (SBR_NF_NO_RANDOM_VAL - 1);
}

static void adjustTimeSlotHQ_GainAndNoise(
    FIXP_DBL *RESTRICT
        ptrReal,
    FIXP_DBL *RESTRICT
        ptrImag,
    HANDLE_SBR_CALCULATE_ENVELOPE h_sbr_cal_env, ENV_CALC_NRGS *nrgs,
    int lowSubband,
    int noSubbands,
    int scale_change,
    FIXP_SGL smooth_ratio,
    int noNoiseFlag,
    int filtBufferNoiseShift)
{
  FIXP_DBL *RESTRICT gain = nrgs->nrgGain;
  FIXP_DBL *RESTRICT noiseLevel =
      nrgs->noiseLevel;
  FIXP_DBL *RESTRICT pSineLevel = nrgs->nrgSine;

  FIXP_DBL *RESTRICT filtBuffer =
      h_sbr_cal_env->filtBuffer;
  FIXP_DBL *RESTRICT filtBufferNoise =
      h_sbr_cal_env->filtBufferNoise;
  int *RESTRICT ptrPhaseIndex =
      &h_sbr_cal_env->phaseIndex;

  int k;
  FIXP_DBL signalReal, signalImag;
  FIXP_DBL noiseReal, noiseImag;
  FIXP_DBL smoothedGain, smoothedNoise;
  FIXP_SGL direct_ratio =
      (FIXP_SGL)MAXVAL_SGL - smooth_ratio;
  int index = *ptrPhaseIndex;
  int shift;
  FIXP_DBL max_val_noise = 0, min_val_noise = 0;
  const FIXP_DBL max_val = MAX_VAL_NRG_HEADROOM >> scale_change;
  const FIXP_DBL min_val = -max_val;

  *ptrPhaseIndex = (index + noSubbands) & (SBR_NF_NO_RANDOM_VAL - 1);

  filtBufferNoiseShift +=
      1;
  if (filtBufferNoiseShift < 0) {
    shift = fixMin(DFRACT_BITS - 1, -filtBufferNoiseShift);
  } else {
    shift = fixMin(DFRACT_BITS - 1, filtBufferNoiseShift);
    max_val_noise = MAX_VAL_NRG_HEADROOM >> shift;
    min_val_noise = -max_val_noise;
  }

  if (smooth_ratio > FL2FXCONST_SGL(0.0f)) {
    for (k = 0; k < noSubbands; k++) {
      smoothedGain =
          fMult(smooth_ratio, filtBuffer[k]) + fMult(direct_ratio, gain[k]);

      if (filtBufferNoiseShift < 0) {
        smoothedNoise = (fMultDiv2(smooth_ratio, filtBufferNoise[k]) >> shift) +
                        fMult(direct_ratio, noiseLevel[k]);
      } else {
        smoothedNoise = fMultDiv2(smooth_ratio, filtBufferNoise[k]);
        smoothedNoise =
            (fMax(fMin(smoothedNoise, max_val_noise), min_val_noise) << shift) +
            fMult(direct_ratio, noiseLevel[k]);
      }

      smoothedNoise = fMax(fMin(smoothedNoise, (FIXP_DBL)(MAXVAL_DBL / 2)),
                           (FIXP_DBL)(MINVAL_DBL / 2));

      signalReal =
          fMax(fMin(fMultDiv2(*ptrReal, smoothedGain), max_val), min_val)
          << scale_change;
      signalImag =
          fMax(fMin(fMultDiv2(*ptrImag, smoothedGain), max_val), min_val)
          << scale_change;

      index++;

      if ((pSineLevel[k] != FL2FXCONST_DBL(0.0f)) || noNoiseFlag) {
        *ptrReal++ = signalReal;
        *ptrImag++ = signalImag;
      } else {
        index &= (SBR_NF_NO_RANDOM_VAL - 1);
        noiseReal =
            fMult(FDK_sbrDecoder_sbr_randomPhase[index][0], smoothedNoise);
        noiseImag =
            fMult(FDK_sbrDecoder_sbr_randomPhase[index][1], smoothedNoise);
        *ptrReal++ = (signalReal + noiseReal);
        *ptrImag++ = (signalImag + noiseImag);
      }
    }
  } else {
    for (k = 0; k < noSubbands; k++) {
      smoothedGain = gain[k];
      signalReal =
          fMax(fMin(fMultDiv2(*ptrReal, smoothedGain), max_val), min_val)
          << scale_change;
      signalImag =
          fMax(fMin(fMultDiv2(*ptrImag, smoothedGain), max_val), min_val)
          << scale_change;

      index++;

      if ((pSineLevel[k] == FL2FXCONST_DBL(0.0f)) && (noNoiseFlag == 0)) {
        smoothedNoise = noiseLevel[k];
        index &= (SBR_NF_NO_RANDOM_VAL - 1);
        noiseReal =
            fMult(FDK_sbrDecoder_sbr_randomPhase[index][0], smoothedNoise);
        noiseImag =
            fMult(FDK_sbrDecoder_sbr_randomPhase[index][1], smoothedNoise);

        signalReal += noiseReal;
        signalImag += noiseImag;
      }
      *ptrReal++ = signalReal;
      *ptrImag++ = signalImag;
    }
  }
}

static void adjustTimeSlotHQ_AddHarmonics(
    FIXP_DBL *RESTRICT
        ptrReal,
    FIXP_DBL *RESTRICT
        ptrImag,
    HANDLE_SBR_CALCULATE_ENVELOPE h_sbr_cal_env, ENV_CALC_NRGS *nrgs,
    int lowSubband,
    int noSubbands,
    int scale_change
) {
  FIXP_DBL *RESTRICT pSineLevel = nrgs->nrgSine;
  UCHAR *RESTRICT ptrHarmIndex =
      &h_sbr_cal_env->harmIndex;

  int k;
  FIXP_DBL signalReal, signalImag;
  UCHAR harmIndex = *ptrHarmIndex;
  int freqInvFlag = (lowSubband & 1);
  FIXP_DBL sineLevel;

  *ptrHarmIndex = (harmIndex + 1) & 3;

  for (k = 0; k < noSubbands; k++) {
    sineLevel = pSineLevel[k];
    freqInvFlag ^= 1;
    if (sineLevel != FL2FXCONST_DBL(0.f)) {
      signalReal = ptrReal[k];
      signalImag = ptrImag[k];
      sineLevel = scaleValue(sineLevel, scale_change);
      if (harmIndex & 2) {
        sineLevel = -sineLevel;
      }
      if (!(harmIndex & 1)) {
        ptrReal[k] = signalReal + sineLevel;
      } else {
        if (!freqInvFlag) sineLevel = -sineLevel;
        ptrImag[k] = signalImag + sineLevel;
      }
    }
  }
}

static void adjustTimeSlotHQ(
    FIXP_DBL *RESTRICT
        ptrReal,
    FIXP_DBL *RESTRICT
        ptrImag,
    HANDLE_SBR_CALCULATE_ENVELOPE h_sbr_cal_env, ENV_CALC_NRGS *nrgs,
    int lowSubband,
    int noSubbands,
    int scale_change,
    FIXP_SGL smooth_ratio,
    int noNoiseFlag,
    int filtBufferNoiseShift)
{
  FIXP_DBL *RESTRICT gain = nrgs->nrgGain;
  FIXP_DBL *RESTRICT noiseLevel =
      nrgs->noiseLevel;
  FIXP_DBL *RESTRICT pSineLevel = nrgs->nrgSine;

  FIXP_DBL *RESTRICT filtBuffer =
      h_sbr_cal_env->filtBuffer;
  FIXP_DBL *RESTRICT filtBufferNoise =
      h_sbr_cal_env->filtBufferNoise;
  UCHAR *RESTRICT ptrHarmIndex =
      &h_sbr_cal_env->harmIndex;
  int *RESTRICT ptrPhaseIndex =
      &h_sbr_cal_env->phaseIndex;

  int k;
  FIXP_DBL signalReal, signalImag;
  FIXP_DBL noiseReal, noiseImag;
  FIXP_DBL smoothedGain, smoothedNoise;
  FIXP_SGL direct_ratio =
      (FIXP_SGL)MAXVAL_SGL - smooth_ratio;
  int index = *ptrPhaseIndex;
  UCHAR harmIndex = *ptrHarmIndex;
  int freqInvFlag = (lowSubband & 1);
  FIXP_DBL sineLevel;
  int shift;
  FIXP_DBL max_val_noise = 0, min_val_noise = 0;
  const FIXP_DBL max_val = MAX_VAL_NRG_HEADROOM >> scale_change;
  const FIXP_DBL min_val = -max_val;

  *ptrPhaseIndex = (index + noSubbands) & (SBR_NF_NO_RANDOM_VAL - 1);
  *ptrHarmIndex = (harmIndex + 1) & 3;

  filtBufferNoiseShift +=
      1;
  if (filtBufferNoiseShift < 0) {
    shift = fixMin(DFRACT_BITS - 1, -filtBufferNoiseShift);
  } else {
    shift = fixMin(DFRACT_BITS - 1, filtBufferNoiseShift);
    max_val_noise = MAX_VAL_NRG_HEADROOM >> shift;
    min_val_noise = -max_val_noise;
  }

  if (smooth_ratio > FL2FXCONST_SGL(0.0f)) {
    for (k = 0; k < noSubbands; k++) {

      smoothedGain =
          fMult(smooth_ratio, filtBuffer[k]) + fMult(direct_ratio, gain[k]);

      if (filtBufferNoiseShift < 0) {
        smoothedNoise = (fMultDiv2(smooth_ratio, filtBufferNoise[k]) >> shift) +
                        fMult(direct_ratio, noiseLevel[k]);
      } else {
        smoothedNoise = fMultDiv2(smooth_ratio, filtBufferNoise[k]);
        smoothedNoise =
            (fMax(fMin(smoothedNoise, max_val_noise), min_val_noise) << shift) +
            fMult(direct_ratio, noiseLevel[k]);
      }

      smoothedNoise = fMax(fMin(smoothedNoise, (FIXP_DBL)(MAXVAL_DBL / 2)),
                           (FIXP_DBL)(MINVAL_DBL / 2));

      signalReal =
          fMax(fMin(fMultDiv2(*ptrReal, smoothedGain), max_val), min_val)
          << scale_change;
      signalImag =
          fMax(fMin(fMultDiv2(*ptrImag, smoothedGain), max_val), min_val)
          << scale_change;

      index++;

      if (pSineLevel[k] != FL2FXCONST_DBL(0.0f)) {
        sineLevel = pSineLevel[k];

        switch (harmIndex) {
          case 0:
            *ptrReal++ = (signalReal + sineLevel);
            *ptrImag++ = (signalImag);
            break;
          case 2:
            *ptrReal++ = (signalReal - sineLevel);
            *ptrImag++ = (signalImag);
            break;
          case 1:
            *ptrReal++ = (signalReal);
            if (freqInvFlag)
              *ptrImag++ = (signalImag - sineLevel);
            else
              *ptrImag++ = (signalImag + sineLevel);
            break;
          case 3:
            *ptrReal++ = signalReal;
            if (freqInvFlag)
              *ptrImag++ = (signalImag + sineLevel);
            else
              *ptrImag++ = (signalImag - sineLevel);
            break;
        }
      } else {
        if (noNoiseFlag) {
          *ptrReal++ = (signalReal);
          *ptrImag++ = (signalImag);
        } else {
          index &= (SBR_NF_NO_RANDOM_VAL - 1);
          noiseReal =
              fMult(FDK_sbrDecoder_sbr_randomPhase[index][0], smoothedNoise);
          noiseImag =
              fMult(FDK_sbrDecoder_sbr_randomPhase[index][1], smoothedNoise);
          *ptrReal++ = (signalReal + noiseReal);
          *ptrImag++ = (signalImag + noiseImag);
        }
      }
      freqInvFlag ^= 1;
    }

  } else {
    for (k = 0; k < noSubbands; k++) {
      smoothedGain = gain[k];
      signalReal =
          fMax(fMin(fMultDiv2(*ptrReal, smoothedGain), max_val), min_val)
          << scale_change;
      signalImag =
          fMax(fMin(fMultDiv2(*ptrImag, smoothedGain), max_val), min_val)
          << scale_change;

      index++;

      if ((sineLevel = pSineLevel[k]) != FL2FXCONST_DBL(0.0f)) {
        switch (harmIndex) {
          case 0:
            signalReal += sineLevel;
            break;
          case 1:
            if (freqInvFlag)
              signalImag -= sineLevel;
            else
              signalImag += sineLevel;
            break;
          case 2:
            signalReal -= sineLevel;
            break;
          case 3:
            if (freqInvFlag)
              signalImag += sineLevel;
            else
              signalImag -= sineLevel;
            break;
        }
      } else {
        if (noNoiseFlag == 0) {
          smoothedNoise = noiseLevel[k];
          index &= (SBR_NF_NO_RANDOM_VAL - 1);
          noiseReal =
              fMult(FDK_sbrDecoder_sbr_randomPhase[index][0], smoothedNoise);
          noiseImag =
              fMult(FDK_sbrDecoder_sbr_randomPhase[index][1], smoothedNoise);

          signalReal += noiseReal;
          signalImag += noiseImag;
        }
      }
      *ptrReal++ = signalReal;
      *ptrImag++ = signalImag;

      freqInvFlag ^= 1;
    }
  }
}

SBR_ERROR
ResetLimiterBands(
    UCHAR *limiterBandTable,
    UCHAR *noLimiterBands,
    UCHAR *freqBandTable,
    int noFreqBands,
    const PATCH_PARAM *patchParam,
    int noPatches,
    int limiterBands,
    UCHAR sbrPatchingMode, int xOverQmf[MAX_NUM_PATCHES], int b41Sbr) {
  int i, k, isPatchBorder[2], loLimIndex, hiLimIndex, tempNoLim, nBands;
  UCHAR workLimiterBandTable[MAX_FREQ_COEFFS / 2 + MAX_NUM_PATCHES + 1];
  int patchBorders[MAX_NUM_PATCHES + 1];
  int kx, k2;

  int lowSubband = freqBandTable[0];
  int highSubband = freqBandTable[noFreqBands];

  if (limiterBands == 0) {
    limiterBandTable[0] = 0;
    limiterBandTable[1] = highSubband - lowSubband;
    nBands = 1;
  } else {
    if (!sbrPatchingMode && xOverQmf != NULL) {
      noPatches = 0;

      if (b41Sbr == 1) {
        for (i = 1; i < MAX_NUM_PATCHES_HBE; i++)
          if (xOverQmf[i] != 0) noPatches++;
      } else {
        for (i = 1; i < MAX_STRETCH_HBE; i++)
          if (xOverQmf[i] != 0) noPatches++;
      }
      for (i = 0; i < noPatches; i++) {
        patchBorders[i] = xOverQmf[i] - lowSubband;
      }
    } else {
      for (i = 0; i < noPatches; i++) {
        patchBorders[i] = patchParam[i].guardStartBand - lowSubband;
      }
    }
    patchBorders[i] = highSubband - lowSubband;

    for (k = 0; k <= noFreqBands; k++) {
      workLimiterBandTable[k] = freqBandTable[k] - lowSubband;
    }
    for (k = 1; k < noPatches; k++) {
      workLimiterBandTable[noFreqBands + k] = patchBorders[k];
    }

    tempNoLim = nBands = noFreqBands + noPatches - 1;
    shellsort(workLimiterBandTable, tempNoLim + 1);

    loLimIndex = 0;
    hiLimIndex = 1;

    while (hiLimIndex <= tempNoLim) {
      FIXP_DBL div_m, oct_m, temp;
      INT div_e = 0, oct_e = 0, temp_e = 0;

      k2 = workLimiterBandTable[hiLimIndex] + lowSubband;
      kx = workLimiterBandTable[loLimIndex] + lowSubband;

      div_m = fDivNorm(k2, kx, &div_e);

      oct_m = fLog2(div_m, div_e, &oct_e);

      temp = fMultNorm(
          oct_m, FDK_sbrDecoder_sbr_limiterBandsPerOctaveDiv4_DBL[limiterBands],
          &temp_e);

      temp_e += oct_e + 2;

      if (temp >> (5 - temp_e) < FL2FXCONST_DBL(0.49f) >> 5) {
        if (workLimiterBandTable[hiLimIndex] ==
            workLimiterBandTable[loLimIndex]) {
          workLimiterBandTable[hiLimIndex] = highSubband;
          nBands--;
          hiLimIndex++;
          continue;
        }
        isPatchBorder[0] = isPatchBorder[1] = 0;
        for (k = 0; k <= noPatches; k++) {
          if (workLimiterBandTable[hiLimIndex] == patchBorders[k]) {
            isPatchBorder[1] = 1;
            break;
          }
        }
        if (!isPatchBorder[1]) {
          workLimiterBandTable[hiLimIndex] = highSubband;
          nBands--;
          hiLimIndex++;
          continue;
        }
        for (k = 0; k <= noPatches; k++) {
          if (workLimiterBandTable[loLimIndex] == patchBorders[k]) {
            isPatchBorder[0] = 1;
            break;
          }
        }
        if (!isPatchBorder[0]) {
          workLimiterBandTable[loLimIndex] = highSubband;
          nBands--;
        }
      }
      loLimIndex = hiLimIndex;
      hiLimIndex++;
    }
    shellsort(workLimiterBandTable, tempNoLim + 1);

    if (nBands > MAX_NUM_LIMITERS || nBands <= 0) {
      return SBRDEC_UNSUPPORTED_CONFIG;
    }

    if (workLimiterBandTable[tempNoLim] > highSubband) {
      return SBRDEC_UNSUPPORTED_CONFIG;
    }

    for (k = 0; k <= nBands; k++) {
      limiterBandTable[k] = workLimiterBandTable[k];
    }
  }
  *noLimiterBands = nBands;

  return SBRDEC_OK;
}
