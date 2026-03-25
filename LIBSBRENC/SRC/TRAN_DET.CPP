/**************************** SBR encoder library ******************************

   Author(s):   Tobias Chalupka

   Description: SBR encoder transient detector

*******************************************************************************/

#include "tran_det.h"
#include "fram_gen.h"
#include "sbrenc_ram.h"
#include "sbr_misc.h"
#include "genericStds.h"

#define NORM_QMF_ENERGY 9.31322574615479E-10

#define ABS_THRES ((FIXP_DBL)16)

#define NRG_SHIFT 3

static FIXP_DBL spectralChange(
    FIXP_DBL Energies[NUMBER_TIME_SLOTS_2304][MAX_FREQ_COEFFS],
    INT *scaleEnergies, FIXP_DBL EnergyTotal, INT nSfb, INT start, INT border,
    INT YBufferWriteOffset, INT stop, INT *result_e) {
  INT i, j;
  INT len1, len2;
  SCHAR energies_e_diff[NUMBER_TIME_SLOTS_2304], energies_e, energyTotal_e = 19,
                                                             energies_e_add;
  SCHAR prevEnergies_e_diff, newEnergies_e_diff;
  FIXP_DBL tmp0, tmp1;
  FIXP_DBL delta, delta_sum;
  INT accu_e, tmp_e;

  delta_sum = FL2FXCONST_DBL(0.0f);
  *result_e = 0;

  len1 = border - start;
  len2 = stop - border;

  FIXP_DBL pos_weight;
  pos_weight = FL2FXCONST_DBL(0.5f) - (len1 * GetInvInt(len1 + len2));
  pos_weight = /*FL2FXCONST_DBL(1.0)*/ (FIXP_DBL)MAXVAL_DBL -
               (fMult(pos_weight, pos_weight) << 2);

  /*** Calc scaling for energies ***/
  FDK_ASSERT(scaleEnergies[0] >= 0);
  FDK_ASSERT(scaleEnergies[1] >= 0);

  energies_e = 19 - fMin(scaleEnergies[0], scaleEnergies[1]);

  /* limit shift for energy accumulation, energies_e can be -10 min. */
  if (energies_e < -10) {
    energies_e_add = -10 - energies_e;
    energies_e = -10;
  } else if (energies_e > 17) {
    energies_e_add = energies_e - 17;
    energies_e = 17;
  } else {
    energies_e_add = 0;
  }

  /* compensate scaling differences between scaleEnergies[0] and
   * scaleEnergies[1]  */
  prevEnergies_e_diff = scaleEnergies[0] -
                        fMin(scaleEnergies[0], scaleEnergies[1]) +
                        energies_e_add + NRG_SHIFT;
  newEnergies_e_diff = scaleEnergies[1] -
                       fMin(scaleEnergies[0], scaleEnergies[1]) +
                       energies_e_add + NRG_SHIFT;

  prevEnergies_e_diff = fMin(prevEnergies_e_diff, DFRACT_BITS - 1);
  newEnergies_e_diff = fMin(newEnergies_e_diff, DFRACT_BITS - 1);

  for (i = start; i < YBufferWriteOffset; i++) {
    energies_e_diff[i] = prevEnergies_e_diff;
  }
  for (i = YBufferWriteOffset; i < stop; i++) {
    energies_e_diff[i] = newEnergies_e_diff;
  }

  /* Sum up energies of all QMF-timeslots for both halfs */
  FDK_ASSERT(len1 <= 8); /* otherwise an overflow is possible */
  FDK_ASSERT(len2 <= 8); /* otherwise an overflow is possible */

  for (j = 0; j < nSfb; j++) {
    FIXP_DBL accu1 = FL2FXCONST_DBL(0.f);
    FIXP_DBL accu2 = FL2FXCONST_DBL(0.f);
    accu_e = energies_e + 3;

    /* Sum up energies in first half */
    for (i = start; i < border; i++) {
      accu1 += scaleValue(Energies[i][j], -energies_e_diff[i]);
    }

    /* Sum up energies in second half */
    for (i = border; i < stop; i++) {
      accu2 += scaleValue(Energies[i][j], -energies_e_diff[i]);
    }

    /* Ensure certain energy to prevent division by zero and to prevent
     * splitting for very low levels */
    accu1 = fMax(accu1, (FIXP_DBL)len1);
    accu2 = fMax(accu2, (FIXP_DBL)len2);

/* Energy change in current band */
#define LN2 FL2FXCONST_DBL(0.6931471806f) /* ln(2) */
    tmp0 = fLog2(accu2, accu_e) - fLog2(accu1, accu_e);
    tmp1 = fLog2((FIXP_DBL)len1, 31) - fLog2((FIXP_DBL)len2, 31);
    delta = fMult(LN2, (tmp0 + tmp1));
    delta = (FIXP_DBL)fAbs(delta);

    /* Weighting with amplitude ratio of this band */
    accu_e++; /* scale at least one bit due to (accu1+accu2) */
    accu1 >>= 1;
    accu2 >>= 1;

    if (accu_e & 1) {
      accu_e++; /* for a defined square result exponent, the exponent has to be
                   even */
      accu1 >>= 1;
      accu2 >>= 1;
    }

    delta_sum += fMult(sqrtFixp(accu1 + accu2), delta);
    *result_e = ((accu_e >> 1) + LD_DATA_SHIFT);
  }

  if (energyTotal_e & 1) {
    energyTotal_e += 1; /* for a defined square result exponent, the exponent
                           has to be even */
    EnergyTotal >>= 1;
  }

  delta_sum = fMult(delta_sum, invSqrtNorm2(EnergyTotal, &tmp_e));
  *result_e = *result_e + (tmp_e - (energyTotal_e >> 1));

  return fMult(delta_sum, pos_weight);
}

static FIXP_DBL addLowbandEnergies(FIXP_DBL **Energies, int *scaleEnergies,
                                   int YBufferWriteOffset, int nrgSzShift,
                                   int tran_off, UCHAR *freqBandTable,
                                   int slots) {
  INT nrgTotal_e;
  FIXP_DBL nrgTotal_m;
  FIXP_DBL accu1 = FL2FXCONST_DBL(0.0f);
  FIXP_DBL accu2 = FL2FXCONST_DBL(0.0f);
  int tran_offdiv2 = tran_off >> nrgSzShift;
  const int sc1 =
      DFRACT_BITS -
      fNormz((FIXP_DBL)fMax(
          1, (freqBandTable[0] * (YBufferWriteOffset - tran_offdiv2) - 1)));
  const int sc2 =
      DFRACT_BITS -
      fNormz((FIXP_DBL)fMax(
          1, (freqBandTable[0] *
                  (tran_offdiv2 + (slots >> nrgSzShift) - YBufferWriteOffset) -
              1)));
  int ts, k;

  for (ts = tran_offdiv2; ts < YBufferWriteOffset; ts++) {
    for (k = 0; k < freqBandTable[0]; k++) {
      accu1 += Energies[ts][k] >> sc1;
    }
  }
  for (; ts < tran_offdiv2 + (slots >> nrgSzShift); ts++) {
    for (k = 0; k < freqBandTable[0]; k++) {
      accu2 += Energies[ts][k] >> sc2;
    }
  }

  nrgTotal_m = fAddNorm(accu1, (sc1 - 5) - scaleEnergies[0], accu2,
                        (sc2 - 5) - scaleEnergies[1], &nrgTotal_e);
  nrgTotal_m = scaleValueSaturate(nrgTotal_m, nrgTotal_e);

  return (nrgTotal_m);
}

static FIXP_DBL addHighbandEnergies(
    FIXP_DBL **RESTRICT Energies, /*!< input */
    INT *scaleEnergies, INT YBufferWriteOffset,
    FIXP_DBL EnergiesM[NUMBER_TIME_SLOTS_2304]
                      [MAX_FREQ_COEFFS], /*!< Combined output */
    UCHAR *RESTRICT freqBandTable, INT nSfb, INT sbrSlots, INT timeStep) {
  INT i, j, k, slotIn, slotOut, scale[2];
  INT li, ui;
  FIXP_DBL nrgTotal;
  FIXP_DBL accu = FL2FXCONST_DBL(0.0f);

  for (slotOut = 0; slotOut < sbrSlots; slotOut++) {
    slotIn = slotOut;

    for (j = 0; j < nSfb; j++) {
      accu = FL2FXCONST_DBL(0.0f);

      li = freqBandTable[j];
      ui = freqBandTable[j + 1];

      for (k = li; k < ui; k++) {
        for (i = 0; i < timeStep; i++) {
          accu += Energies[slotIn][k] >> 5;
        }
      }
      EnergiesM[slotOut][j] = accu;
    }
  }

  /* scale energies down before add up */
  scale[0] = fixMin(8, scaleEnergies[0]);
  scale[1] = fixMin(8, scaleEnergies[1]);

  if ((scaleEnergies[0] - scale[0]) > (DFRACT_BITS - 1) ||
      (scaleEnergies[1] - scale[1]) > (DFRACT_BITS - 1))
    nrgTotal = FL2FXCONST_DBL(0.0f);
  else {
    /* Now add all energies */
    accu = FL2FXCONST_DBL(0.0f);

    for (slotOut = 0; slotOut < YBufferWriteOffset; slotOut++) {
      for (j = 0; j < nSfb; j++) {
        accu += (EnergiesM[slotOut][j] >> scale[0]);
      }
    }
    nrgTotal = accu >> (scaleEnergies[0] - scale[0]);

    for (slotOut = YBufferWriteOffset; slotOut < sbrSlots; slotOut++) {
      for (j = 0; j < nSfb; j++) {
        accu += (EnergiesM[slotOut][j] >> scale[0]);
      }
    }
    nrgTotal = fAddSaturate(nrgTotal, accu >> (scaleEnergies[1] - scale[1]));
  }

  return (nrgTotal);
}

void FDKsbrEnc_frameSplitter(
    FIXP_DBL **Energies, INT *scaleEnergies,
    HANDLE_SBR_TRANSIENT_DETECTOR h_sbrTransientDetector, UCHAR *freqBandTable,
    UCHAR *tran_vector, int YBufferWriteOffset, int YBufferSzShift, int nSfb,
    int timeStep, int no_cols, FIXP_DBL *tonality) {
  if (tran_vector[1] == 0) /* no transient was detected */
  {
    FIXP_DBL delta;
    INT delta_e;
    FIXP_DBL(*EnergiesM)[MAX_FREQ_COEFFS];
    FIXP_DBL EnergyTotal, newLowbandEnergy, newHighbandEnergy;
    INT border;
    INT sbrSlots = fMultI(GetInvInt(timeStep), no_cols);
    C_ALLOC_SCRATCH_START(_EnergiesM, FIXP_DBL,
                          NUMBER_TIME_SLOTS_2304 * MAX_FREQ_COEFFS)

    FDK_ASSERT(sbrSlots * timeStep == no_cols);

    EnergiesM = (FIXP_DBL(*)[MAX_FREQ_COEFFS])_EnergiesM;

    newLowbandEnergy = addLowbandEnergies(
        Energies, scaleEnergies, YBufferWriteOffset, YBufferSzShift,
        h_sbrTransientDetector->tran_off, freqBandTable, no_cols);

    newHighbandEnergy =
        addHighbandEnergies(Energies, scaleEnergies, YBufferWriteOffset,
                            EnergiesM, freqBandTable, nSfb, sbrSlots, timeStep);

    {
      EnergyTotal = (newLowbandEnergy >> 1) +
                    (h_sbrTransientDetector->prevLowBandEnergy >>
                     1); /* mean of new and prev LB NRG */
      EnergyTotal =
          fAddSaturate(EnergyTotal, newHighbandEnergy); /* Add HB NRG */
      border = (sbrSlots + 1) >> 1;

      if ((INT)EnergyTotal & 0xffffffe0 &&
          (scaleEnergies[0] < 32 || scaleEnergies[1] < 32)) /* i.e. > 31 */ {
        delta = spectralChange(EnergiesM, scaleEnergies, EnergyTotal, nSfb, 0,
                               border, YBufferWriteOffset, sbrSlots, &delta_e);
      } else {
        delta = FL2FXCONST_DBL(0.0f);
        delta_e = 0;

        *tonality = FL2FXCONST_DBL(0.0f);
      }

      if (fIsLessThan(h_sbrTransientDetector->split_thr_m,
                      h_sbrTransientDetector->split_thr_e, delta, delta_e)) {
        tran_vector[0] = 1; /* Set flag for splitting */
      } else {
        tran_vector[0] = 0;
      }
    }

    h_sbrTransientDetector->prevLowBandEnergy = newLowbandEnergy;
    h_sbrTransientDetector->prevHighBandEnergy = newHighbandEnergy;
    C_ALLOC_SCRATCH_END(_EnergiesM, FIXP_DBL,
                        NUMBER_TIME_SLOTS_2304 * MAX_FREQ_COEFFS)
  }
}

static void calculateThresholds(FIXP_DBL **RESTRICT Energies,
                                INT *RESTRICT scaleEnergies,
                                FIXP_DBL *RESTRICT thresholds,
                                int YBufferWriteOffset, int YBufferSzShift,
                                int noCols, int noRows, int tran_off) {
  FIXP_DBL mean_val, std_val, temp;
  FIXP_DBL i_noCols;
  FIXP_DBL i_noCols1;
  FIXP_DBL accu, accu0, accu1;
  int scaleFactor0, scaleFactor1, commonScale;
  int i, j;

  i_noCols = GetInvInt(noCols + tran_off) << YBufferSzShift;
  i_noCols1 = GetInvInt(noCols + tran_off - 1) << YBufferSzShift;

  commonScale = fixMin(scaleEnergies[0], scaleEnergies[1]);

  scaleFactor0 = fixMin((scaleEnergies[0] - commonScale), (DFRACT_BITS - 1));
  scaleFactor1 = fixMin((scaleEnergies[1] - commonScale), (DFRACT_BITS - 1));

  FDK_ASSERT((scaleFactor0 >= 0) && (scaleFactor1 >= 0));

  for (i = 0; i < noRows; i++) {
    int startEnergy = (tran_off >> YBufferSzShift);
    int endEnergy = ((noCols >> YBufferSzShift) + tran_off);
    int shift;

    accu0 = accu1 = FL2FXCONST_DBL(0.0f);

    for (j = startEnergy; j < YBufferWriteOffset; j++)
      accu0 = fMultAddDiv2(accu0, Energies[j][i], i_noCols);
    for (; j < endEnergy; j++)
      accu1 = fMultAddDiv2(accu1, Energies[j][i], i_noCols);

    mean_val = ((accu0 << 1) >> scaleFactor0) +
               ((accu1 << 1) >> scaleFactor1);
    shift = fixMax(
        0, CountLeadingBits(mean_val) -
               6);

    accu = FL2FXCONST_DBL(0.0f);

    for (j = startEnergy; j < YBufferWriteOffset; j++) {
      temp = ((FIXP_DBL)mean_val - ((FIXP_DBL)Energies[j][i] >> scaleFactor0))
             << shift;
      temp = fPow2Div2(temp);
      accu = fMultAddDiv2(accu, temp, i_noCols1);
    }
    for (; j < endEnergy; j++) {
      temp = ((FIXP_DBL)mean_val - ((FIXP_DBL)Energies[j][i] >> scaleFactor1))
             << shift;
      temp = fPow2Div2(temp);
      accu = fMultAddDiv2(accu, temp, i_noCols1);
    }
    accu <<= 2;
    std_val = sqrtFixp(accu) >> shift;

    temp = (commonScale <= (DFRACT_BITS - 1))
               ? fMult(FL2FXCONST_DBL(0.66f), thresholds[i]) +
                     (fMult(FL2FXCONST_DBL(0.34f), std_val) >> commonScale)
               : (FIXP_DBL)0;

    thresholds[i] = fixMax(ABS_THRES, temp);

    FDK_ASSERT(commonScale >= 0);
  }
}

static void extractTransientCandidates(
    FIXP_DBL **RESTRICT Energies, INT *RESTRICT scaleEnergies,
    FIXP_DBL *RESTRICT thresholds, FIXP_DBL *RESTRICT transients,
    int YBufferWriteOffset, int YBufferSzShift, int noCols, int start_band,
    int stop_band, int tran_off, int addPrevSamples) {
  FIXP_DBL i_thres;
  C_ALLOC_SCRATCH_START(EnergiesTemp, FIXP_DBL, 2 * 32)
  int tmpScaleEnergies0, tmpScaleEnergies1;
  int endCond;
  int startEnerg, endEnerg;
  int i, j, jIndex, jpBM;

  tmpScaleEnergies0 = scaleEnergies[0];
  tmpScaleEnergies1 = scaleEnergies[1];

  tmpScaleEnergies0 = fixMin(tmpScaleEnergies0, MAX_SHIFT_DBL);
  tmpScaleEnergies1 = fixMin(tmpScaleEnergies1, MAX_SHIFT_DBL);

  FDK_ASSERT((tmpScaleEnergies0 >= 0) && (tmpScaleEnergies1 >= 0));

  FDKmemmove(transients, transients + noCols - addPrevSamples,
             (tran_off + addPrevSamples) * sizeof(FIXP_DBL));
  FDKmemclear(transients + tran_off + addPrevSamples,
              noCols * sizeof(FIXP_DBL));

  endCond = noCols;
  startEnerg = (tran_off - 3) >> YBufferSzShift;
  endEnerg =
      ((noCols + (YBufferWriteOffset << YBufferSzShift)) - 1) >>
      YBufferSzShift;

  for (i = start_band; i < stop_band; i++) {
    FIXP_DBL thres = thresholds[i];

    if ((LONG)thresholds[i] >= 256)
      i_thres = (LONG)((LONG)MAXVAL_DBL / ((((LONG)thresholds[i])) + 1))
                << (32 - 24);
    else
      i_thres = (LONG)MAXVAL_DBL;

    if (YBufferSzShift == 1) {
      for (j = startEnerg; j < YBufferWriteOffset; j++) {
        FIXP_DBL tmp = Energies[j][i];
        EnergiesTemp[(j << 1) + 1] = EnergiesTemp[j << 1] =
            tmp >> tmpScaleEnergies0;
      }
      for (; j <= endEnerg; j++) {
        FIXP_DBL tmp = Energies[j][i];
        EnergiesTemp[(j << 1) + 1] = EnergiesTemp[j << 1] =
            tmp >> tmpScaleEnergies1;
      }
    } else {
      for (j = startEnerg; j < YBufferWriteOffset; j++) {
        FIXP_DBL tmp = Energies[j][i];
        EnergiesTemp[j] = tmp >> tmpScaleEnergies0;
      }
      for (; j <= endEnerg; j++) {
        FIXP_DBL tmp = Energies[j][i];
        EnergiesTemp[j] = tmp >> tmpScaleEnergies1;
      }
    }

    jIndex = tran_off;
    jpBM = jIndex + addPrevSamples;

    for (j = endCond; j--; jIndex++, jpBM++) {
      FIXP_DBL delta, tran;
      int d;

      delta = (FIXP_DBL)0;
      tran = (FIXP_DBL)0;

      for (d = 1; d < 4; d++) {
        delta += EnergiesTemp[jIndex + d];
        delta -= EnergiesTemp[jIndex - d];
        delta -= thres;

        if (delta > (FIXP_DBL)0) {
          tran = fMultAddDiv2(tran, i_thres, delta);
        }
      }
      transients[jpBM] += (tran << 1);
    }
  }
  C_ALLOC_SCRATCH_END(EnergiesTemp, FIXP_DBL, 2 * 32)
}

void FDKsbrEnc_transientDetect(HANDLE_SBR_TRANSIENT_DETECTOR h_sbrTran,
                               FIXP_DBL **Energies, INT *scaleEnergies,
                               UCHAR *transient_info, int YBufferWriteOffset,
                               int YBufferSzShift, int timeStep,
                               int frameMiddleBorder) {
  int no_cols = h_sbrTran->no_cols;
  int qmfStartSample;
  int addPrevSamples;
  int timeStepShift = 0;
  int i, cond;

  qmfStartSample = timeStep * frameMiddleBorder;
  addPrevSamples = (qmfStartSample > 0) ? 0 : 1;

  switch (timeStep) {
    case 1:
      timeStepShift = 0;
      break;
    case 2:
      timeStepShift = 1;
      break;
    case 4:
      timeStepShift = 2;
      break;
  }

  calculateThresholds(Energies, scaleEnergies, h_sbrTran->thresholds,
                      YBufferWriteOffset, YBufferSzShift, h_sbrTran->no_cols,
                      h_sbrTran->no_rows, h_sbrTran->tran_off);

  extractTransientCandidates(
      Energies, scaleEnergies, h_sbrTran->thresholds, h_sbrTran->transients,
      YBufferWriteOffset, YBufferSzShift, h_sbrTran->no_cols, 0,
      h_sbrTran->no_rows, h_sbrTran->tran_off, addPrevSamples);

  transient_info[0] = 0;
  transient_info[1] = 0;
  transient_info[2] = 0;

  qmfStartSample += addPrevSamples;

  for (i = qmfStartSample; i < qmfStartSample + no_cols; i++) {
    cond = (h_sbrTran->transients[i] <
            fMult(FL2FXCONST_DBL(0.9f), h_sbrTran->transients[i - 1])) &&
           (h_sbrTran->transients[i - 1] > h_sbrTran->tran_thr);

    if (cond) {
      transient_info[0] = (i - qmfStartSample) >> timeStepShift;
      transient_info[1] = 1;
      break;
    }
  }

  if (h_sbrTran->frameShift != 0) {
    for (i = qmfStartSample + no_cols;
         i < qmfStartSample + no_cols + h_sbrTran->frameShift; i++) {
      cond = (h_sbrTran->transients[i] <
              fMult(FL2FXCONST_DBL(0.9f), h_sbrTran->transients[i - 1])) &&
             (h_sbrTran->transients[i - 1] > h_sbrTran->tran_thr);

      if (cond) {
        int pos = (int)((i - qmfStartSample - no_cols) >> timeStepShift);
        if ((pos < 3) && (transient_info[1] == 0)) {
          transient_info[2] = 1;
        }
        break;
      }
    }
  }
}

int FDKsbrEnc_InitSbrTransientDetector(
    HANDLE_SBR_TRANSIENT_DETECTOR h_sbrTransientDetector,
    UINT sbrSyntaxFlags,
    INT frameSize, INT sampleFreq, sbrConfigurationPtr params, int tran_fc,
    int no_cols, int no_rows, int YBufferWriteOffset, int YBufferSzShift,
    int frameShift, int tran_off) {
  INT totalBitrate =
      params->codecSettings.standardBitrate * params->codecSettings.nChannels;
  INT codecBitrate = params->codecSettings.bitRate;
  FIXP_DBL bitrateFactor_m, framedur_fix;
  INT bitrateFactor_e, tmp_e;

  FDKmemclear(h_sbrTransientDetector, sizeof(SBR_TRANSIENT_DETECTOR));

  h_sbrTransientDetector->frameShift = frameShift;
  h_sbrTransientDetector->tran_off = tran_off;

  if (codecBitrate) {
    bitrateFactor_m = fDivNorm((FIXP_DBL)totalBitrate,
                               (FIXP_DBL)(codecBitrate << 2), &bitrateFactor_e);
    bitrateFactor_e += 2;
  } else {
    bitrateFactor_m = FL2FXCONST_DBL(1.0 / 4.0);
    bitrateFactor_e = 2;
  }

  framedur_fix = fDivNorm(frameSize, sampleFreq);

  FIXP_DBL tmp = framedur_fix - FL2FXCONST_DBL(0.010);

  tmp = fixMax(tmp, FL2FXCONST_DBL(0.0001));
  tmp = fDivNorm(FL2FXCONST_DBL(0.000075), fPow2(tmp), &tmp_e);

  bitrateFactor_e = (tmp_e + bitrateFactor_e);

  if (sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY) {
    bitrateFactor_e--;
  }

  FDK_ASSERT(no_cols <= 32);
  FDK_ASSERT(no_rows <= 64);

  h_sbrTransientDetector->no_cols = no_cols;
  h_sbrTransientDetector->tran_thr =
      (FIXP_DBL)((params->tran_thr << (32 - 24 - 1)) / no_rows);
  h_sbrTransientDetector->tran_fc = tran_fc;
  h_sbrTransientDetector->split_thr_m = fMult(tmp, bitrateFactor_m);
  h_sbrTransientDetector->split_thr_e = bitrateFactor_e;
  h_sbrTransientDetector->no_rows = no_rows;
  h_sbrTransientDetector->mode = params->tran_det_mode;
  h_sbrTransientDetector->prevLowBandEnergy = FL2FXCONST_DBL(0.0f);

  return (0);
}

#define ENERGY_SCALING_SIZE 32

INT FDKsbrEnc_InitSbrFastTransientDetector(
    HANDLE_FAST_TRAN_DET h_sbrFastTransientDetector,
    const INT time_slots_per_frame, const INT bandwidth_qmf_slot,
    const INT no_qmf_channels, const INT sbr_qmf_1st_band) {
  int i;
  int buff_size;
  FIXP_DBL myExp;
  FIXP_DBL myExpSlot;

  h_sbrFastTransientDetector->lookahead = TRAN_DET_LOOKAHEAD;
  h_sbrFastTransientDetector->nTimeSlots = time_slots_per_frame;

  buff_size = h_sbrFastTransientDetector->nTimeSlots +
              h_sbrFastTransientDetector->lookahead;

  for (i = 0; i < buff_size; i++) {
    h_sbrFastTransientDetector->delta_energy[i] = FL2FXCONST_DBL(0.0f);
    h_sbrFastTransientDetector->energy_timeSlots[i] = FL2FXCONST_DBL(0.0f);
    h_sbrFastTransientDetector->lowpass_energy[i] = FL2FXCONST_DBL(0.0f);
    h_sbrFastTransientDetector->transientCandidates[i] = 0;
  }

  FDK_ASSERT(bandwidth_qmf_slot > 0.f);
  h_sbrFastTransientDetector->stopBand =
      fMin(TRAN_DET_STOP_FREQ / bandwidth_qmf_slot, no_qmf_channels);
  h_sbrFastTransientDetector->startBand =
      fMin(sbr_qmf_1st_band,
           h_sbrFastTransientDetector->stopBand - TRAN_DET_MIN_QMFBANDS);

  FDK_ASSERT(h_sbrFastTransientDetector->startBand < no_qmf_channels);
  FDK_ASSERT(h_sbrFastTransientDetector->startBand <
             h_sbrFastTransientDetector->stopBand);
  FDK_ASSERT(h_sbrFastTransientDetector->startBand > 1);
  FDK_ASSERT(h_sbrFastTransientDetector->stopBand > 1);

  FDK_ASSERT(h_sbrFastTransientDetector->stopBand -
                 h_sbrFastTransientDetector->startBand <=
             64);

#define EXP_E 7
  myExp = fMultNorm(QMF_HP_dBd_SLOPE_FIX, 0, (FIXP_DBL)bandwidth_qmf_slot,
                    DFRACT_BITS - 1, EXP_E);
  myExpSlot = myExp;

  for (i = 0; i < 64; i++) {
    FIXP_DBL dBf_m;
    INT dBf_e;
    INT tmp;

    INT dBf_int;
    FIXP_DBL dBf_fract;

    FIXP_DBL myExp_int =
        (myExpSlot & (FIXP_DBL)0xfe000000) + (FIXP_DBL)0x02000000;

    FIXP_DBL myExp_fract = myExp_int - myExpSlot;

    dBf_int = CalcInvLdData(myExp_int);

    if (dBf_int <=
        46340) {
      dBf_int *= dBf_int;
      dBf_fract = CalcInvLdData(-myExp_fract);
      dBf_fract = fMultNorm(dBf_fract, dBf_fract, &tmp);
      dBf_e = (DFRACT_BITS - 1 - tmp) - CountLeadingBits(dBf_int);
      dBf_m = fMultNorm(dBf_int, DFRACT_BITS - 1, dBf_fract, tmp, dBf_e);

      myExpSlot += myExp;
    } else {
      dBf_m = (FIXP_DBL)0;
      dBf_e = 0;
    }

    h_sbrFastTransientDetector->dBf_m[i] = dBf_m;
    h_sbrFastTransientDetector->dBf_e[i] = dBf_e;
  }

  return 0;
}

void FDKsbrEnc_fastTransientDetect(
    const HANDLE_FAST_TRAN_DET h_sbrFastTransientDetector,
    const FIXP_DBL *const *Energies, const int *const scaleEnergies,
    const INT YBufferWriteOffset, UCHAR *const tran_vector) {
  int timeSlot, band;

  FIXP_DBL max_delta_energy;
  int max_delta_energy_scale;
  int ind_max = 0;
  int isTransientInFrame = 0;

  const int nTimeSlots = h_sbrFastTransientDetector->nTimeSlots;
  const int lookahead = h_sbrFastTransientDetector->lookahead;
  const int startBand = h_sbrFastTransientDetector->startBand;
  const int stopBand = h_sbrFastTransientDetector->stopBand;

  int *transientCandidates = h_sbrFastTransientDetector->transientCandidates;

  FIXP_DBL *energy_timeSlots = h_sbrFastTransientDetector->energy_timeSlots;
  int *energy_timeSlots_scale =
      h_sbrFastTransientDetector->energy_timeSlots_scale;

  FIXP_DBL *delta_energy = h_sbrFastTransientDetector->delta_energy;
  int *delta_energy_scale = h_sbrFastTransientDetector->delta_energy_scale;

  const FIXP_DBL thr = TRAN_DET_THRSHLD;
  const INT thr_scale = TRAN_DET_THRSHLD_SCALE;

  tran_vector[2] = 0;

  FDKmemclear(transientCandidates + lookahead, nTimeSlots * sizeof(int));

  for (timeSlot = lookahead; timeSlot < nTimeSlots + lookahead; timeSlot++) {
    int i, norm;
    FIXP_DBL tmpE = FL2FXCONST_DBL(0.0f);
    int headroomEnSlot = DFRACT_BITS - 1;

    FIXP_DBL smallNRG = FL2FXCONST_DBL(1e-2f);
    FIXP_DBL denominator;
    INT denominator_scale;

    for (band = startBand; band < stopBand; band++) {
      int tmp_headroom = fNormz(Energies[timeSlot][band]) - 1;
      if (tmp_headroom < headroomEnSlot) {
        headroomEnSlot = tmp_headroom;
      }
    }

    for (i = 0, band = startBand; band < stopBand; band++, i++) {
      FIXP_DBL weightedEnergy =
          fMult(Energies[timeSlot][band] << headroomEnSlot,
                h_sbrFastTransientDetector->dBf_m[i]);

      tmpE +=
          weightedEnergy >> (6 + (10 - h_sbrFastTransientDetector->dBf_e[i]));
    }

    energy_timeSlots[timeSlot] = tmpE;

    if (timeSlot < YBufferWriteOffset) {
      energy_timeSlots_scale[timeSlot] =
          (-scaleEnergies[0] + 2 * QMF_SCALE_OFFSET) + (10 + 6) -
          headroomEnSlot;
    } else {
      energy_timeSlots_scale[timeSlot] =
          (-scaleEnergies[1] + 2 * QMF_SCALE_OFFSET) + (10 + 6) -
          headroomEnSlot;
    }

    if (-energy_timeSlots_scale[timeSlot - 1] + 1 > 5) {
      denominator = smallNRG;
      denominator_scale = 0;
    } else {
      smallNRG =
          scaleValue(smallNRG, -(energy_timeSlots_scale[timeSlot - 1] + 1));
      denominator = (energy_timeSlots[timeSlot - 1] >> 1) + smallNRG;
      denominator_scale = energy_timeSlots_scale[timeSlot - 1] + 1;
    }

    delta_energy[timeSlot] =
        fDivNorm(energy_timeSlots[timeSlot], denominator, &norm);
    delta_energy_scale[timeSlot] =
        energy_timeSlots_scale[timeSlot] - denominator_scale + norm;
  }

  FDK_ASSERT(lookahead >= 2);
  for (timeSlot = lookahead; timeSlot < nTimeSlots + lookahead; timeSlot++) {
    FIXP_DBL energy_cur_slot_weighted =
        fMult(energy_timeSlots[timeSlot], FL2FXCONST_DBL(1.0f / 1.4f));
    if (!fIsLessThan(delta_energy[timeSlot], delta_energy_scale[timeSlot], thr,
                     thr_scale) &&
        (((transientCandidates[timeSlot - 2] == 0) &&
          (transientCandidates[timeSlot - 1] == 0)) ||
         !fIsLessThan(energy_cur_slot_weighted,
                      energy_timeSlots_scale[timeSlot],
                      energy_timeSlots[timeSlot - 1],
                      energy_timeSlots_scale[timeSlot - 1]) ||
         !fIsLessThan(energy_cur_slot_weighted,
                      energy_timeSlots_scale[timeSlot],
                      energy_timeSlots[timeSlot - 2],
                      energy_timeSlots_scale[timeSlot - 2]))) {
      transientCandidates[timeSlot] = 1;
    }
  }

  max_delta_energy = FL2FXCONST_DBL(0.0f);
  max_delta_energy_scale = 0;
  ind_max = 0;
  isTransientInFrame = 0;
  for (timeSlot = 0; timeSlot < nTimeSlots; timeSlot++) {
    int scale = fMax(delta_energy_scale[timeSlot], max_delta_energy_scale);
    if (transientCandidates[timeSlot] &&
        ((delta_energy[timeSlot] >> (scale - delta_energy_scale[timeSlot])) >
         (max_delta_energy >> (scale - max_delta_energy_scale)))) {
      max_delta_energy = delta_energy[timeSlot];
      max_delta_energy_scale = scale;
      ind_max = timeSlot;
      isTransientInFrame = 1;
    }
  }

  if (isTransientInFrame) {
    tran_vector[0] = ind_max;
    tran_vector[1] = 1;
  } else {
    tran_vector[0] = tran_vector[1] = 0;
  }

  for (timeSlot = nTimeSlots; timeSlot < nTimeSlots + lookahead; timeSlot++) {
    if (transientCandidates[timeSlot]) {
      tran_vector[2] = 1;
    }
  }

  for (timeSlot = 0; timeSlot < lookahead; timeSlot++) {
    transientCandidates[timeSlot] = transientCandidates[nTimeSlots + timeSlot];

    energy_timeSlots[timeSlot] = energy_timeSlots[nTimeSlots + timeSlot];
    energy_timeSlots_scale[timeSlot] =
        energy_timeSlots_scale[nTimeSlots + timeSlot];

    delta_energy[timeSlot] = delta_energy[nTimeSlots + timeSlot];
    delta_energy_scale[timeSlot] = delta_energy_scale[nTimeSlots + timeSlot];
  }
}
