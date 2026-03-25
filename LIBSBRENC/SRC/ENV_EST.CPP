/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "env_est.h"
#include "tran_det.h"
#include "qmf.h"
#include "fram_gen.h"
#include "bit_sbr.h"
#include "cmondata.h"
#include "sbrenc_ram.h"
#include "genericStds.h"

#define QUANT_ERROR_THRES 200
#define Y_NRG_SCALE 5
#define MAX_NRG_SLOTS_LD 16

static const UCHAR panTable[2][10] = {{0, 2, 4, 6, 8, 12, 16, 20, 24},
                                      {0, 2, 4, 8, 12, 0, 0, 0, 0}};
static const UCHAR maxIndex[2] = {9, 5};

static FIXP_DBL FDKsbrEnc_GetTonality(const FIXP_DBL *const *quotaMatrix,
                                      const INT noEstPerFrame,
                                      const INT startIndex,
                                      const FIXP_DBL *const *Energies,
                                      const UCHAR startBand, const INT stopBand,
                                      const INT numberCols) {
  UCHAR b, e, k;
  INT no_enMaxBand[SBR_MAX_ENERGY_VALUES] = {-1, -1, -1, -1, -1};
  FIXP_DBL energyMax[SBR_MAX_ENERGY_VALUES] = {
      FL2FXCONST_DBL(0.0f), FL2FXCONST_DBL(0.0f), FL2FXCONST_DBL(0.0f),
      FL2FXCONST_DBL(0.0f), FL2FXCONST_DBL(0.0f)};
  FIXP_DBL energyMaxMin = MAXVAL_DBL;
  UCHAR posEnergyMaxMin = 0;
  FIXP_DBL tonalityBand[SBR_MAX_ENERGY_VALUES] = {
      FL2FXCONST_DBL(0.0f), FL2FXCONST_DBL(0.0f), FL2FXCONST_DBL(0.0f),
      FL2FXCONST_DBL(0.0f), FL2FXCONST_DBL(0.0f)};
  FIXP_DBL globalTonality = FL2FXCONST_DBL(0.0f);
  FIXP_DBL energyBand[64];
  INT maxNEnergyValues;

  FDK_ASSERT(numberCols == 15 || numberCols == 16);
  if (numberCols == 15) {
    for (b = startBand; b < stopBand; b++) {
      energyBand[b] = FL2FXCONST_DBL(0.0f);
    }
  } else {
    for (b = startBand; b < stopBand; b++) {
      energyBand[b] = Energies[15][b] >> 4;
    }
  }

  for (k = 0; k < 15; k++) {
    for (b = startBand; b < stopBand; b++) {
      energyBand[b] += Energies[k][b] >> 4;
    }
  }

  maxNEnergyValues = fMin(SBR_MAX_ENERGY_VALUES, stopBand - startBand);

  energyMaxMin = energyMax[0] = energyBand[startBand];
  no_enMaxBand[0] = startBand;
  posEnergyMaxMin = 0;
  for (k = 1; k < maxNEnergyValues; k++) {
    energyMax[k] = energyBand[startBand + k];
    no_enMaxBand[k] = startBand + k;
    if (energyMaxMin > energyMax[k]) {
      energyMaxMin = energyMax[k];
      posEnergyMaxMin = k;
    }
  }

  for (b = startBand + maxNEnergyValues; b < stopBand; b++) {
    if (energyBand[b] > energyMaxMin) {
      energyMax[posEnergyMaxMin] = energyBand[b];
      no_enMaxBand[posEnergyMaxMin] = b;

      energyMaxMin = energyMax[0];
      posEnergyMaxMin = 0;
      for (k = 1; k < maxNEnergyValues; k++) {
        if (energyMaxMin > energyMax[k]) {
          energyMaxMin = energyMax[k];
          posEnergyMaxMin = k;
        }
      }
    }
  }
  for (e = 0; e < maxNEnergyValues; e++) {
    tonalityBand[e] = FL2FXCONST_DBL(0.0f);
    for (k = 0; k < noEstPerFrame; k++) {
      tonalityBand[e] += quotaMatrix[startIndex + k][no_enMaxBand[e]] >> 1;
    }
    globalTonality +=
        tonalityBand[e] >> 2;
  }

  return globalTonality;
}

LNK_SECTION_CODE_L1
static void FDKsbrEnc_getEnergyFromCplxQmfData(
    FIXP_DBL **RESTRICT energyValues,
    FIXP_DBL **RESTRICT realValues,
    FIXP_DBL **RESTRICT
        imagValues,
    INT numberBands,
    INT numberCols,
    INT *qmfScale,
    INT *energyScale)
{
  int j, k;
  int scale;
  FIXP_DBL max_val = FL2FXCONST_DBL(0.0f);

  C_ALLOC_SCRATCH_START(tmpNrg, FIXP_DBL, 32 * 64 / 2)

  scale = DFRACT_BITS;
  for (k = 0; k < numberCols; k++) {
    scale = fixMin(scale, fixMin(getScalefactor(realValues[k], numberBands),
                                 getScalefactor(imagValues[k], numberBands)));
  }

  if (scale >= DFRACT_BITS - 1) {
    scale = (FRACT_BITS - 1 - *qmfScale);
  }
  scale = fixMax(0, scale - 1);

  *qmfScale += scale;

  {
    FIXP_DBL *nrgValues = tmpNrg;
    for (k = 0; k < numberCols; k += 2) {
      FIXP_DBL *RESTRICT r0 = realValues[k];
      FIXP_DBL *RESTRICT i0 = imagValues[k];
      FIXP_DBL *RESTRICT r1 = realValues[k + 1];
      FIXP_DBL *RESTRICT i1 = imagValues[k + 1];
      for (j = 0; j < numberBands; j++) {
        FIXP_DBL energy;
        FIXP_DBL tr0, tr1, ti0, ti1;

        tr0 = r0[j];
        tr1 = r1[j];
        ti0 = i0[j];
        ti1 = i1[j];

        tr0 <<= scale;
        ti0 <<= scale;
        energy = fPow2AddDiv2(fPow2Div2(tr0), ti0) >> 1;

        tr1 <<= scale;
        ti1 <<= scale;
        energy += fPow2AddDiv2(fPow2Div2(tr1), ti1) >> 1;

        *nrgValues++ = energy;
        max_val = fixMax(max_val, energy);

        r0[j] = tr0;
        r1[j] = tr1;
        i0[j] = ti0;
        i1[j] = ti1;
      }
    }
  }
  *energyScale =
      2 * (*qmfScale) -
      1;

  scale = CountLeadingBits(max_val);
  {
    FIXP_DBL *nrgValues = tmpNrg;
    for (k = 0; k<numberCols>> 1; k++) {
      scaleValues(energyValues[k], nrgValues, numberBands, scale);
      nrgValues += numberBands;
    }
    *energyScale += scale;
  }

  C_ALLOC_SCRATCH_END(tmpNrg, FIXP_DBL, 32 * 64 / 2)
}

LNK_SECTION_CODE_L1
static void FDKsbrEnc_getEnergyFromCplxQmfDataFull(
    FIXP_DBL **RESTRICT energyValues,
    FIXP_DBL **RESTRICT realValues,
    FIXP_DBL **RESTRICT
        imagValues,
    int numberBands,
    int numberCols,
    int *qmfScale,
    int *energyScale)
{
  int j, k;
  int scale;
  FIXP_DBL max_val = FL2FXCONST_DBL(0.0f);

  C_ALLOC_SCRATCH_START(tmpNrg, FIXP_DBL, MAX_NRG_SLOTS_LD * 64)

  FDK_ASSERT(numberCols <= MAX_NRG_SLOTS_LD);
  FDK_ASSERT(numberBands <= 64);

  scale = DFRACT_BITS;
  for (k = 0; k < numberCols; k++) {
    scale = fixMin(scale, fixMin(getScalefactor(realValues[k], numberBands),
                                 getScalefactor(imagValues[k], numberBands)));
  }

  if (scale >= DFRACT_BITS - 1) {
    scale = (FRACT_BITS - 1 - *qmfScale);
  }
  scale = fixMax(0, scale - 1);

  *qmfScale += scale;

  {
    FIXP_DBL *nrgValues = tmpNrg;
    for (k = 0; k < numberCols; k++) {
      FIXP_DBL *RESTRICT r0 = realValues[k];
      FIXP_DBL *RESTRICT i0 = imagValues[k];
      for (j = 0; j < numberBands; j++) {
        FIXP_DBL energy;
        FIXP_DBL tr0, ti0;

        tr0 = r0[j];
        ti0 = i0[j];

        tr0 <<= scale;
        ti0 <<= scale;
        energy = fPow2AddDiv2(fPow2Div2(tr0), ti0);
        *nrgValues++ = energy;

        max_val = fixMax(max_val, energy);

        r0[j] = tr0;
        i0[j] = ti0;
      }
    }
  }
  *energyScale =
      2 * (*qmfScale) -
      1;

  scale = CountLeadingBits(max_val);
  {
    FIXP_DBL *nrgValues = tmpNrg;
    for (k = 0; k < numberCols; k++) {
      scaleValues(energyValues[k], nrgValues, numberBands, scale);
      nrgValues += numberBands;
    }
    *energyScale += scale;
  }

  C_ALLOC_SCRATCH_END(tmpNrg, FIXP_DBL, MAX_NRG_SLOTS_LD * 64)
}

static INT mapPanorama(INT nrgVal,
                       INT ampRes,
                       INT *quantError
) {
  int i;
  INT min_val, val;
  UCHAR panIndex;
  INT sign;

  sign = nrgVal > 0 ? 1 : -1;

  nrgVal *= sign;

  min_val = FDK_INT_MAX;
  panIndex = 0;
  for (i = 0; i < maxIndex[ampRes]; i++) {
    val = fixp_abs((nrgVal - (INT)panTable[ampRes][i]));

    if (val < min_val) {
      min_val = val;
      panIndex = i;
    }
  }

  *quantError = min_val;

  return panTable[ampRes][maxIndex[ampRes] - 1] +
         sign * panTable[ampRes][panIndex];
}

static void sbrNoiseFloorLevelsQuantisation(
    SCHAR *RESTRICT iNoiseLevels,
    FIXP_DBL *RESTRICT
        NoiseLevels,
    INT coupling
) {
  INT i;
  INT tmp, dummy;

  for (i = 0; i < MAX_NUM_NOISE_VALUES; i++) {
    if ((FIXP_DBL)NoiseLevels[i] > FL2FXCONST_DBL(0.46875f)) {
      tmp = 30;
    } else {
      tmp = ((((INT)NoiseLevels[i]) >>
              (DFRACT_BITS - 1 - LD_DATA_SHIFT)));
      if (tmp != 0) tmp += 1;
    }

    if (coupling) {
      tmp = tmp < -30 ? -30 : tmp;
      tmp = mapPanorama(tmp, 1, &dummy);
    }
    iNoiseLevels[i] = tmp;
  }
}

static void coupleNoiseFloor(
    FIXP_DBL *RESTRICT noise_level_left,
    FIXP_DBL *RESTRICT noise_level_right
) {
  FIXP_DBL cmpValLeft, cmpValRight;
  INT i;
  FIXP_DBL temp1, temp2;

  for (i = 0; i < MAX_NUM_NOISE_VALUES; i++) {
    cmpValLeft = NOISE_FLOOR_OFFSET_64 - noise_level_left[i];
    cmpValRight = NOISE_FLOOR_OFFSET_64 - noise_level_right[i];

    if (cmpValRight < FL2FXCONST_DBL(0.0f)) {
      temp1 = CalcInvLdData(NOISE_FLOOR_OFFSET_64 - noise_level_right[i]);
    } else {
      temp1 = CalcInvLdData(NOISE_FLOOR_OFFSET_64 - noise_level_right[i]);
      temp1 = temp1 << (DFRACT_BITS - 1 - LD_DATA_SHIFT -
                        1);
    }

    if (cmpValLeft < FL2FXCONST_DBL(0.0f)) {
      temp2 = CalcInvLdData(NOISE_FLOOR_OFFSET_64 - noise_level_left[i]);
    } else {
      temp2 = CalcInvLdData(NOISE_FLOOR_OFFSET_64 - noise_level_left[i]);
      temp2 = temp2 << (DFRACT_BITS - 1 - LD_DATA_SHIFT -
                        1);
    }

    if ((cmpValLeft < FL2FXCONST_DBL(0.0f)) &&
        (cmpValRight < FL2FXCONST_DBL(0.0f))) {
      noise_level_left[i] =
          NOISE_FLOOR_OFFSET_64 -
          (CalcLdData(
              ((temp1 >> 1) +
               (temp2 >> 1))));
      noise_level_right[i] = CalcLdData(temp2) - CalcLdData(temp1);
    }

    if ((cmpValLeft >= FL2FXCONST_DBL(0.0f)) &&
        (cmpValRight >= FL2FXCONST_DBL(0.0f))) {
      noise_level_left[i] = NOISE_FLOOR_OFFSET_64 -
                            (CalcLdData(((temp1 >> 1) + (temp2 >> 1))) +
                             FL2FXCONST_DBL(0.109375f));
      noise_level_right[i] = CalcLdData(temp2) - CalcLdData(temp1);
    }

    if ((cmpValLeft >= FL2FXCONST_DBL(0.0f)) &&
        (cmpValRight < FL2FXCONST_DBL(0.0f))) {
      noise_level_left[i] = NOISE_FLOOR_OFFSET_64 -
                            (CalcLdData(((temp1 >> (7 + 1)) + (temp2 >> 1))) +
                             FL2FXCONST_DBL(0.109375f));
      noise_level_right[i] =
          (CalcLdData(temp2) + FL2FXCONST_DBL(0.109375f)) - CalcLdData(temp1);
    }

    if ((cmpValLeft < FL2FXCONST_DBL(0.0f)) &&
        (cmpValRight >= FL2FXCONST_DBL(0.0f))) {
      noise_level_left[i] = NOISE_FLOOR_OFFSET_64 -
                            (CalcLdData(((temp1 >> 1) + (temp2 >> (7 + 1)))) +
                             FL2FXCONST_DBL(0.109375f));
      noise_level_right[i] = CalcLdData(temp2) -
                             (CalcLdData(temp1) +
                              FL2FXCONST_DBL(0.109375f));
    }
  }
}

static FIXP_DBL getEnvSfbEnergy(
    INT li,
    INT ui,
    INT start_pos,
    INT stop_pos,
    INT border_pos,
    FIXP_DBL **YBuffer,
    INT YBufferSzShift,
    INT scaleNrg0,
    INT scaleNrg1)
{
  int sc0, sc1, k, l;

  FIXP_DBL nrgSum, nrg1, nrg2, accu1, accu2;
  INT dynScale, dynScale1, dynScale2;
  if (ui - li == 0)
    dynScale = DFRACT_BITS - 1;
  else
    dynScale = CalcLdInt(ui - li) >> (DFRACT_BITS - 1 - LD_DATA_SHIFT);

  sc0 = fixMin(scaleNrg0, Y_NRG_SCALE);
  sc1 = fixMin(scaleNrg1, Y_NRG_SCALE);
  dynScale1 = fixMin((scaleNrg0 - sc0), dynScale);
  dynScale2 = fixMin((scaleNrg1 - sc1), dynScale);
  nrgSum = accu1 = accu2 = (FIXP_DBL)0;

  for (k = li; k < ui; k++) {
    nrg1 = nrg2 = (FIXP_DBL)0;
    for (l = start_pos; l < border_pos; l++) {
      nrg1 += YBuffer[l >> YBufferSzShift][k] >> sc0;
    }
    for (; l < stop_pos; l++) {
      nrg2 += YBuffer[l >> YBufferSzShift][k] >> sc1;
    }
    accu1 = fAddSaturate(accu1, (nrg1 >> dynScale1));
    accu2 = fAddSaturate(accu2, (nrg2 >> dynScale2));
  }
  nrgSum +=
      (accu1 >> fixMin((scaleNrg0 - sc0 - dynScale1), (DFRACT_BITS - 1))) +
      (accu2 >> fixMin((scaleNrg1 - sc1 - dynScale2), (DFRACT_BITS - 1)));

  return nrgSum;
}

static FIXP_DBL mhLoweringEnergy(FIXP_DBL nrg, INT M) {
  if (M > 2) {
    INT tmpScale;
    tmpScale = CountLeadingBits(nrg);
    nrg <<= tmpScale;
    nrg = fMult(nrg, FL2FXCONST_DBL(0.398107267f));
    nrg >>= tmpScale;
  } else {
    if (M > 1) {
      nrg >>= 1;
    }
  }

  return nrg;
}

static FIXP_DBL nmhLoweringEnergy(FIXP_DBL nrg, const FIXP_DBL nrgSum,
                                  const INT nrgSum_scale, const INT M) {
  if (nrg > FL2FXCONST_DBL(0)) {
    int sc = 0;
    FIXP_DBL gain = fMult(fDivNorm(nrgSum, nrg, &sc), GetInvInt(M + 1));
    sc += nrgSum_scale;

    if (!((sc >= 0) && (gain > ((FIXP_DBL)MAXVAL_DBL >> sc)))) {
      nrg = fMult(scaleValue(gain, sc), nrg);
    }
  }
  return nrg;
}

static void calculateSbrEnvelope(
    FIXP_DBL **RESTRICT YBufferLeft,
    FIXP_DBL **RESTRICT YBufferRight,
    int *RESTRICT YBufferScaleLeft,
    int *RESTRICT YBufferScaleRight,
    const SBR_FRAME_INFO *frame_info,
    SCHAR *RESTRICT sfb_nrgLeft,
    SCHAR *RESTRICT sfb_nrgRight,
    HANDLE_SBR_CONFIG_DATA h_con,
    HANDLE_ENV_CHANNEL h_sbr,
    SBR_STEREO_MODE stereoMode,
    INT *maxQuantError,
    int YBufferSzShift)

{
  int env, j, m = 0;
  INT no_of_bands, start_pos, stop_pos, li, ui;
  FREQ_RES freq_res;

  INT ca = 2 - h_sbr->encEnvData.init_sbr_amp_res;
  INT oneBitLess = 0;
  if (ca == 2)
    oneBitLess =
        1;

  INT quantError;
  INT nEnvelopes = frame_info->nEnvelopes;
  INT short_env = frame_info->shortEnv - 1;
  INT timeStep = h_sbr->sbrExtractEnvelope.time_step;
  INT commonScale, scaleLeft0, scaleLeft1;
  INT scaleRight0 = 0, scaleRight1 = 0;

  commonScale = fixMin(YBufferScaleLeft[0], YBufferScaleLeft[1]);

  if (stereoMode == SBR_COUPLING) {
    commonScale = fixMin(commonScale, YBufferScaleRight[0]);
    commonScale = fixMin(commonScale, YBufferScaleRight[1]);
  }

  commonScale = commonScale - 7;

  scaleLeft0 = YBufferScaleLeft[0] - commonScale;
  scaleLeft1 = YBufferScaleLeft[1] - commonScale;
  FDK_ASSERT((scaleLeft0 >= 0) && (scaleLeft1 >= 0));

  if (stereoMode == SBR_COUPLING) {
    scaleRight0 = YBufferScaleRight[0] - commonScale;
    scaleRight1 = YBufferScaleRight[1] - commonScale;
    FDK_ASSERT((scaleRight0 >= 0) && (scaleRight1 >= 0));
    *maxQuantError = 0;
  }

  for (env = 0; env < nEnvelopes; env++) {
    FIXP_DBL pNrgLeft[32];
    FIXP_DBL pNrgRight[32];
    int envNrg_scale;
    FIXP_DBL envNrgLeft = FL2FXCONST_DBL(0.0f);
    FIXP_DBL envNrgRight = FL2FXCONST_DBL(0.0f);
    int missingHarmonic[32];
    int count[32];

    start_pos = timeStep * frame_info->borders[env];
    stop_pos = timeStep * frame_info->borders[env + 1];
    freq_res = frame_info->freqRes[env];
    no_of_bands = h_con->nSfb[freq_res];
    envNrg_scale = DFRACT_BITS - fNormz((FIXP_DBL)no_of_bands);
    if (env == short_env) {
      j = fMax(2, timeStep);
      if ((stop_pos - start_pos - j) > 0) {
        stop_pos = stop_pos - j;
      }
    }
    for (j = 0; j < no_of_bands; j++) {
      FIXP_DBL nrgLeft = FL2FXCONST_DBL(0.0f);
      FIXP_DBL nrgRight = FL2FXCONST_DBL(0.0f);

      li = h_con->freqBandTable[freq_res][j];
      ui = h_con->freqBandTable[freq_res][j + 1];

      if (freq_res == FREQ_RES_HIGH) {
        if (j == 0 && ui - li > 1) {
          li++;
        }
      } else {
        if (j == 0 && ui - li > 2) {
          li++;
        }
      }

      missingHarmonic[j] = 0;

      if (h_sbr->encEnvData.addHarmonicFlag) {
        if (freq_res == FREQ_RES_HIGH) {
          if (h_sbr->encEnvData
                  .addHarmonic[j]) {
            missingHarmonic[j] = 1;
          }
        } else {
          INT i;
          INT startBandHigh = 0;
          INT stopBandHigh = 0;

          while (h_con->freqBandTable[FREQ_RES_HIGH][startBandHigh] <
                 h_con->freqBandTable[FREQ_RES_LOW][j])
            startBandHigh++;
          while (h_con->freqBandTable[FREQ_RES_HIGH][stopBandHigh] <
                 h_con->freqBandTable[FREQ_RES_LOW][j + 1])
            stopBandHigh++;

          for (i = startBandHigh; i < stopBandHigh; i++) {
            if (h_sbr->encEnvData.addHarmonic[i]) {
              missingHarmonic[j] = 1;
            }
          }
        }
      }

      int border_pos =
          fixMin(stop_pos, h_sbr->sbrExtractEnvelope.YBufferWriteOffset
                               << YBufferSzShift);

      if (missingHarmonic[j]) {
        int k;
        count[j] = stop_pos - start_pos;
        nrgLeft = FL2FXCONST_DBL(0.0f);

        for (k = li; k < ui; k++) {
          FIXP_DBL tmpNrg;
          tmpNrg = getEnvSfbEnergy(k, k + 1, start_pos, stop_pos, border_pos,
                                   YBufferLeft, YBufferSzShift, scaleLeft0,
                                   scaleLeft1);

          nrgLeft = fixMax(nrgLeft, tmpNrg);
        }

        nrgLeft = mhLoweringEnergy(nrgLeft, ui - li);

        if (stereoMode == SBR_COUPLING) {
          nrgRight = FL2FXCONST_DBL(0.0f);

          for (k = li; k < ui; k++) {
            FIXP_DBL tmpNrg;
            tmpNrg = getEnvSfbEnergy(k, k + 1, start_pos, stop_pos, border_pos,
                                     YBufferRight, YBufferSzShift, scaleRight0,
                                     scaleRight1);

            nrgRight = fixMax(nrgRight, tmpNrg);
          }

          nrgRight = mhLoweringEnergy(nrgRight, ui - li);
        }
      }
      else {
        count[j] = (stop_pos - start_pos) * (ui - li);

        nrgLeft = getEnvSfbEnergy(li, ui, start_pos, stop_pos, border_pos,
                                  YBufferLeft, YBufferSzShift, scaleLeft0,
                                  scaleLeft1);

        if (stereoMode == SBR_COUPLING) {
          nrgRight = getEnvSfbEnergy(li, ui, start_pos, stop_pos, border_pos,
                                     YBufferRight, YBufferSzShift, scaleRight0,
                                     scaleRight1);
        }
      }

      pNrgLeft[j] = nrgLeft;
      pNrgRight[j] = nrgRight;
      envNrgLeft += (nrgLeft >> envNrg_scale);
      envNrgRight += (nrgRight >> envNrg_scale);
    }

    for (j = 0; j < no_of_bands; j++) {
      FIXP_DBL nrgLeft2 = FL2FXCONST_DBL(0.0f);
      FIXP_DBL nrgLeft = pNrgLeft[j];
      FIXP_DBL nrgRight = pNrgRight[j];

      if (!missingHarmonic[j] && h_sbr->fLevelProtect) {
        nrgLeft =
            nmhLoweringEnergy(nrgLeft, envNrgLeft, envNrg_scale, no_of_bands);
        if (stereoMode == SBR_COUPLING) {
          nrgRight = nmhLoweringEnergy(nrgRight, envNrgRight, envNrg_scale,
                                       no_of_bands);
        }
      }

      if (stereoMode == SBR_COUPLING) {
        nrgLeft2 = nrgLeft;
        nrgLeft = (nrgRight + nrgLeft) >> 1;
      }

      if (nrgLeft > FL2FXCONST_DBL(0.0f)) {
        FIXP_DBL tmp0, tmp1, tmp2, tmp3;
        INT tmpScale;

        tmpScale = CountLeadingBits(nrgLeft);
        nrgLeft = nrgLeft << tmpScale;

        tmp0 = CalcLdData(nrgLeft);
        tmp1 = ((FIXP_DBL)(commonScale + tmpScale))
               << (DFRACT_BITS - 1 - LD_DATA_SHIFT - 1);
        tmp2 = ((FIXP_DBL)(count[j] * 64)) << (DFRACT_BITS - 1 - 14 - 1);
        tmp2 = CalcLdData(tmp2);
        tmp3 = FL2FXCONST_DBL(0.6875f - 0.21875f - 0.015625f) >>
               1;

        nrgLeft = ((tmp0 - tmp2) >> 1) + (tmp3 - tmp1);
      } else {
        nrgLeft = FL2FXCONST_DBL(-1.0f);
      }

      nrgLeft = fixMin(fixMax(nrgLeft, FL2FXCONST_DBL(0.0f)),
                       (FL2FXCONST_DBL(0.5f) >> oneBitLess));
      nrgLeft = (FIXP_DBL)(LONG)nrgLeft >>
                (DFRACT_BITS - 1 - LD_DATA_SHIFT - 1 - oneBitLess - 1);
      sfb_nrgLeft[m] = ((INT)nrgLeft + 1) >> 1;

      if (stereoMode == SBR_COUPLING) {
        FIXP_DBL scaleFract;
        int sc0, sc1;

        nrgLeft2 = fixMax((FIXP_DBL)0x1, nrgLeft2);
        nrgRight = fixMax((FIXP_DBL)0x1, nrgRight);

        sc0 = CountLeadingBits(nrgLeft2);
        sc1 = CountLeadingBits(nrgRight);

        scaleFract =
            ((FIXP_DBL)(sc0 - sc1))
            << (DFRACT_BITS - 1 -
                LD_DATA_SHIFT);
        nrgRight = CalcLdData(nrgLeft2 << sc0) - CalcLdData(nrgRight << sc1) -
                   scaleFract;

        nrgRight = (FIXP_DBL)(LONG)(nrgRight) >>
                   (DFRACT_BITS - 1 - LD_DATA_SHIFT - 1 - oneBitLess);
        nrgRight = (nrgRight + (FIXP_DBL)1) >> 1;

        sfb_nrgRight[m] = mapPanorama(
            nrgRight, h_sbr->encEnvData.init_sbr_amp_res, &quantError);

        *maxQuantError = fixMax(quantError, *maxQuantError);
      }

      m++;
    }

    if (h_con->useParametricCoding) {
      m -= no_of_bands;
      for (j = 0; j < no_of_bands; j++) {
        if (freq_res == FREQ_RES_HIGH &&
            h_sbr->sbrExtractEnvelope.envelopeCompensation[j]) {
          sfb_nrgLeft[m] -=
              (ca *
               fixp_abs(
                   (INT)h_sbr->sbrExtractEnvelope.envelopeCompensation[j]));
        }
        sfb_nrgLeft[m] = fixMax(0, sfb_nrgLeft[m]);
        m++;
      }
    }
  }
}

LNK_SECTION_CODE_L1
void FDKsbrEnc_extractSbrEnvelope1(
    HANDLE_SBR_CONFIG_DATA h_con,
    HANDLE_SBR_HEADER_DATA sbrHeaderData,
    HANDLE_SBR_BITSTREAM_DATA sbrBitstreamData, HANDLE_ENV_CHANNEL hEnvChan,
    HANDLE_COMMON_DATA hCmonData, SBR_ENV_TEMP_DATA *eData,
    SBR_FRAME_TEMP_DATA *fData) {
  HANDLE_SBR_EXTRACT_ENVELOPE sbrExtrEnv = &hEnvChan->sbrExtractEnvelope;

  if (sbrExtrEnv->YBufferSzShift == 0)
    FDKsbrEnc_getEnergyFromCplxQmfDataFull(
        &sbrExtrEnv->YBuffer[sbrExtrEnv->YBufferWriteOffset],
        sbrExtrEnv->rBuffer + sbrExtrEnv->rBufferReadOffset,
        sbrExtrEnv->iBuffer + sbrExtrEnv->rBufferReadOffset, h_con->noQmfBands,
        sbrExtrEnv->no_cols, &hEnvChan->qmfScale, &sbrExtrEnv->YBufferScale[1]);
  else
    FDKsbrEnc_getEnergyFromCplxQmfData(
        &sbrExtrEnv->YBuffer[sbrExtrEnv->YBufferWriteOffset],
        sbrExtrEnv->rBuffer + sbrExtrEnv->rBufferReadOffset,
        sbrExtrEnv->iBuffer + sbrExtrEnv->rBufferReadOffset, h_con->noQmfBands,
        sbrExtrEnv->no_cols, &hEnvChan->qmfScale, &sbrExtrEnv->YBufferScale[1]);

  FDKsbrEnc_CalculateTonalityQuotas(
      &hEnvChan->TonCorr, sbrExtrEnv->rBuffer, sbrExtrEnv->iBuffer,
      h_con->freqBandTable[HI][h_con->nSfb[HI]], hEnvChan->qmfScale);

  if (h_con->sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY) {
    FIXP_DBL tonality = FDKsbrEnc_GetTonality(
        hEnvChan->TonCorr.quotaMatrix,
        hEnvChan->TonCorr.numberOfEstimatesPerFrame,
        hEnvChan->TonCorr.startIndexMatrix,
        sbrExtrEnv->YBuffer + sbrExtrEnv->YBufferWriteOffset,
        h_con->freqBandTable[HI][0] + 1, h_con->noQmfBands,
        sbrExtrEnv->no_cols);

    hEnvChan->encEnvData.ton_HF[1] = hEnvChan->encEnvData.ton_HF[0];
    hEnvChan->encEnvData.ton_HF[0] = tonality;

    hEnvChan->encEnvData.global_tonality =
        (hEnvChan->encEnvData.ton_HF[0] >> 1) +
        (hEnvChan->encEnvData.ton_HF[1] >> 1);
  }

  if (h_con->sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY) {
    FDKsbrEnc_fastTransientDetect(&hEnvChan->sbrFastTransientDetector,
                                  sbrExtrEnv->YBuffer, sbrExtrEnv->YBufferScale,
                                  sbrExtrEnv->YBufferWriteOffset,
                                  eData->transient_info);

  } else {
    FDKsbrEnc_transientDetect(
        &hEnvChan->sbrTransientDetector, sbrExtrEnv->YBuffer,
        sbrExtrEnv->YBufferScale, eData->transient_info,
        sbrExtrEnv->YBufferWriteOffset, sbrExtrEnv->YBufferSzShift,
        sbrExtrEnv->time_step, hEnvChan->SbrEnvFrame.frameMiddleSlot);
  }

  FDKsbrEnc_frameSplitter(
      sbrExtrEnv->YBuffer, sbrExtrEnv->YBufferScale,
      &hEnvChan->sbrTransientDetector, h_con->freqBandTable[1],
      eData->transient_info, sbrExtrEnv->YBufferWriteOffset,
      sbrExtrEnv->YBufferSzShift, h_con->nSfb[1], sbrExtrEnv->time_step,
      sbrExtrEnv->no_cols, &hEnvChan->encEnvData.global_tonality);
}

LNK_SECTION_CODE_L1
void FDKsbrEnc_extractSbrEnvelope2(
    HANDLE_SBR_CONFIG_DATA h_con,
    HANDLE_SBR_HEADER_DATA sbrHeaderData,
    HANDLE_PARAMETRIC_STEREO hParametricStereo,
    HANDLE_SBR_BITSTREAM_DATA sbrBitstreamData, HANDLE_ENV_CHANNEL h_envChan0,
    HANDLE_ENV_CHANNEL h_envChan1, HANDLE_COMMON_DATA hCmonData,
    SBR_ENV_TEMP_DATA *eData, SBR_FRAME_TEMP_DATA *fData, int clearOutput) {
  HANDLE_ENV_CHANNEL h_envChan[MAX_NUM_CHANNELS] = {h_envChan0, h_envChan1};
  int ch, i, j, c, YSzShift = h_envChan[0]->sbrExtractEnvelope.YBufferSzShift;

  SBR_STEREO_MODE stereoMode = h_con->stereoMode;
  int nChannels = h_con->nChannels;
  const int *v_tuning;
  static const int v_tuningHEAAC[6] = {0, 2, 4, 0, 0, 0};

  static const int v_tuningELD[6] = {0, 2, 3, 0, 0, 0};

  if (h_con->sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY)
    v_tuning = v_tuningELD;
  else
    v_tuning = v_tuningHEAAC;

  if (stereoMode == SBR_COUPLING) {
    if (eData[0].transient_info[1] && eData[1].transient_info[1]) {
      eData[0].transient_info[0] =
          fixMin(eData[1].transient_info[0], eData[0].transient_info[0]);
      eData[1].transient_info[0] = eData[0].transient_info[0];
    } else {
      if (eData[0].transient_info[1] && !eData[1].transient_info[1]) {
        eData[1].transient_info[0] = eData[0].transient_info[0];
      } else {
        if (!eData[0].transient_info[1] && eData[1].transient_info[1])
          eData[0].transient_info[0] = eData[1].transient_info[0];
        else {
          eData[0].transient_info[0] =
              fixMax(eData[1].transient_info[0], eData[0].transient_info[0]);
          eData[1].transient_info[0] = eData[0].transient_info[0];
        }
      }
    }
  }

  eData[0].frame_info = FDKsbrEnc_frameInfoGenerator(
      &h_envChan[0]->SbrEnvFrame, eData[0].transient_info,
      sbrBitstreamData->rightBorderFIX,
      h_envChan[0]->sbrExtractEnvelope.pre_transient_info,
      h_envChan[0]->encEnvData.ldGrid, v_tuning);

  h_envChan[0]->encEnvData.hSbrBSGrid = &h_envChan[0]->SbrEnvFrame.SbrGrid;

  if (h_envChan[0]->encEnvData.ldGrid && eData[0].transient_info[2]) {
    h_envChan[0]->SbrEnvFrame.SbrFrameInfo.shortEnv =
        h_envChan[0]->SbrEnvFrame.SbrFrameInfo.nEnvelopes;
  }

  switch (stereoMode) {
    case SBR_LEFT_RIGHT:
    case SBR_SWITCH_LRC:
      eData[1].frame_info = FDKsbrEnc_frameInfoGenerator(
          &h_envChan[1]->SbrEnvFrame, eData[1].transient_info,
          sbrBitstreamData->rightBorderFIX,
          h_envChan[1]->sbrExtractEnvelope.pre_transient_info,
          h_envChan[1]->encEnvData.ldGrid, v_tuning);

      h_envChan[1]->encEnvData.hSbrBSGrid = &h_envChan[1]->SbrEnvFrame.SbrGrid;

      if (h_envChan[1]->encEnvData.ldGrid && eData[1].transient_info[2]) {
        h_envChan[1]->SbrEnvFrame.SbrFrameInfo.shortEnv =
            h_envChan[1]->SbrEnvFrame.SbrFrameInfo.nEnvelopes;
      }

      if (eData[0].frame_info->nEnvelopes != eData[1].frame_info->nEnvelopes) {
        stereoMode = SBR_LEFT_RIGHT;
      } else {
        for (i = 0; i < eData[0].frame_info->nEnvelopes + 1; i++) {
          if (eData[0].frame_info->borders[i] !=
              eData[1].frame_info->borders[i]) {
            stereoMode = SBR_LEFT_RIGHT;
            break;
          }
        }
        for (i = 0; i < eData[0].frame_info->nEnvelopes; i++) {
          if (eData[0].frame_info->freqRes[i] !=
              eData[1].frame_info->freqRes[i]) {
            stereoMode = SBR_LEFT_RIGHT;
            break;
          }
        }
        if (eData[0].frame_info->shortEnv != eData[1].frame_info->shortEnv) {
          stereoMode = SBR_LEFT_RIGHT;
        }
      }
      break;
    case SBR_COUPLING:
      eData[1].frame_info = eData[0].frame_info;
      h_envChan[1]->encEnvData.hSbrBSGrid = &h_envChan[0]->SbrEnvFrame.SbrGrid;
      break;
    case SBR_MONO:
      break;
    default:
      FDK_ASSERT(0);
  }

  for (ch = 0; ch < nChannels; ch++) {
    HANDLE_ENV_CHANNEL hEnvChan = h_envChan[ch];
    HANDLE_SBR_EXTRACT_ENVELOPE sbrExtrEnv = &hEnvChan->sbrExtractEnvelope;
    SBR_ENV_TEMP_DATA *ed = &eData[ch];

    sbrExtrEnv->pre_transient_info[0] = ed->transient_info[0];
    sbrExtrEnv->pre_transient_info[1] = ed->transient_info[1];
    hEnvChan->encEnvData.noOfEnvelopes = ed->nEnvelopes =
        ed->frame_info->nEnvelopes;
    hEnvChan->encEnvData.currentAmpResFF = (AMP_RES)h_con->initAmpResFF;

    if ((hEnvChan->encEnvData.hSbrBSGrid->frameClass == FIXFIX) &&
        (ed->nEnvelopes == 1)) {
      AMP_RES currentAmpResFF = SBR_AMP_RES_1_5;
      if (h_con->sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY) {
        if (fIsLessThan(
                        h_con->thresholdAmpResFF_m, h_con->thresholdAmpResFF_e,
                        hEnvChan->encEnvData.global_tonality,
                        RELAXATION_SHIFT + 2)) {
          hEnvChan->encEnvData.currentAmpResFF = SBR_AMP_RES_1_5;
        } else {
          hEnvChan->encEnvData.currentAmpResFF = SBR_AMP_RES_3_0;
        }
        currentAmpResFF = hEnvChan->encEnvData.currentAmpResFF;
      }

      if (currentAmpResFF != hEnvChan->encEnvData.init_sbr_amp_res) {
        FDKsbrEnc_InitSbrHuffmanTables(
            &hEnvChan->encEnvData, &hEnvChan->sbrCodeEnvelope,
            &hEnvChan->sbrCodeNoiseFloor, currentAmpResFF);
      }
    } else {
      if (sbrHeaderData->sbr_amp_res != hEnvChan->encEnvData.init_sbr_amp_res) {
        FDKsbrEnc_InitSbrHuffmanTables(
            &hEnvChan->encEnvData, &hEnvChan->sbrCodeEnvelope,
            &hEnvChan->sbrCodeNoiseFloor, sbrHeaderData->sbr_amp_res);
      }
    }

    if (!clearOutput) {
      FDKsbrEnc_TonCorrParamExtr(
          &hEnvChan->TonCorr, hEnvChan->encEnvData.sbr_invf_mode_vec,
          ed->noiseFloor, &hEnvChan->encEnvData.addHarmonicFlag,
          hEnvChan->encEnvData.addHarmonic, sbrExtrEnv->envelopeCompensation,
          ed->frame_info, ed->transient_info, h_con->freqBandTable[HI],
          h_con->nSfb[HI], hEnvChan->encEnvData.sbr_xpos_mode,
          h_con->sbrSyntaxFlags);
    }

    if (hEnvChan->sbrTransientDetector.prevLowBandEnergy <
            hEnvChan->sbrTransientDetector.prevHighBandEnergy &&
        hEnvChan->sbrTransientDetector.prevHighBandEnergy > FL2FX_DBL(0.03)
        && !(h_con->sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY)) {
      hEnvChan->fLevelProtect = 1;

      for (i = 0; i < MAX_NUM_NOISE_VALUES; i++)
        hEnvChan->encEnvData.sbr_invf_mode_vec[i] = INVF_HIGH_LEVEL;
    } else {
      hEnvChan->fLevelProtect = 0;
    }

    hEnvChan->encEnvData.sbr_invf_mode =
        hEnvChan->encEnvData.sbr_invf_mode_vec[0];

    hEnvChan->encEnvData.noOfnoisebands =
        hEnvChan->TonCorr.sbrNoiseFloorEstimate.noNoiseBands;
  }

  for (ch = 0; ch < nChannels; ch++) {
    for (i = 0; i < eData[ch].nEnvelopes; i++) {
      h_envChan[ch]->encEnvData.noScfBands[i] =
          (eData[ch].frame_info->freqRes[i] == FREQ_RES_HIGH
               ? h_con->nSfb[FREQ_RES_HIGH]
               : h_con->nSfb[FREQ_RES_LOW]);
    }
  }

  if (h_con->sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY &&
      stereoMode == SBR_SWITCH_LRC &&
      h_envChan[0]->encEnvData.currentAmpResFF !=
          h_envChan[1]->encEnvData.currentAmpResFF) {
    stereoMode = SBR_LEFT_RIGHT;
  }

  switch (stereoMode) {
    case SBR_MONO:
      calculateSbrEnvelope(h_envChan[0]->sbrExtractEnvelope.YBuffer, NULL,
                           h_envChan[0]->sbrExtractEnvelope.YBufferScale, NULL,
                           eData[0].frame_info, eData[0].sfb_nrg, NULL, h_con,
                           h_envChan[0], SBR_MONO, NULL, YSzShift);
      break;
    case SBR_LEFT_RIGHT:
      calculateSbrEnvelope(h_envChan[0]->sbrExtractEnvelope.YBuffer, NULL,
                           h_envChan[0]->sbrExtractEnvelope.YBufferScale, NULL,
                           eData[0].frame_info, eData[0].sfb_nrg, NULL, h_con,
                           h_envChan[0], SBR_MONO, NULL, YSzShift);
      calculateSbrEnvelope(h_envChan[1]->sbrExtractEnvelope.YBuffer, NULL,
                           h_envChan[1]->sbrExtractEnvelope.YBufferScale, NULL,
                           eData[1].frame_info, eData[1].sfb_nrg, NULL, h_con,
                           h_envChan[1], SBR_MONO, NULL, YSzShift);
      break;
    case SBR_COUPLING:
      calculateSbrEnvelope(h_envChan[0]->sbrExtractEnvelope.YBuffer,
                           h_envChan[1]->sbrExtractEnvelope.YBuffer,
                           h_envChan[0]->sbrExtractEnvelope.YBufferScale,
                           h_envChan[1]->sbrExtractEnvelope.YBufferScale,
                           eData[0].frame_info, eData[0].sfb_nrg,
                           eData[1].sfb_nrg, h_con, h_envChan[0], SBR_COUPLING,
                           &fData->maxQuantError, YSzShift);
      break;
    case SBR_SWITCH_LRC:
      calculateSbrEnvelope(h_envChan[0]->sbrExtractEnvelope.YBuffer, NULL,
                           h_envChan[0]->sbrExtractEnvelope.YBufferScale, NULL,
                           eData[0].frame_info, eData[0].sfb_nrg, NULL, h_con,
                           h_envChan[0], SBR_MONO, NULL, YSzShift);
      calculateSbrEnvelope(h_envChan[1]->sbrExtractEnvelope.YBuffer, NULL,
                           h_envChan[1]->sbrExtractEnvelope.YBufferScale, NULL,
                           eData[1].frame_info, eData[1].sfb_nrg, NULL, h_con,
                           h_envChan[1], SBR_MONO, NULL, YSzShift);
      calculateSbrEnvelope(h_envChan[0]->sbrExtractEnvelope.YBuffer,
                           h_envChan[1]->sbrExtractEnvelope.YBuffer,
                           h_envChan[0]->sbrExtractEnvelope.YBufferScale,
                           h_envChan[1]->sbrExtractEnvelope.YBufferScale,
                           eData[0].frame_info, eData[0].sfb_nrg_coupling,
                           eData[1].sfb_nrg_coupling, h_con, h_envChan[0],
                           SBR_COUPLING, &fData->maxQuantError, YSzShift);
      break;
  }

  switch (stereoMode) {
    case SBR_MONO:
      sbrNoiseFloorLevelsQuantisation(eData[0].noise_level, eData[0].noiseFloor,
                                      0);

      FDKsbrEnc_codeEnvelope(eData[0].noise_level, fData->res,
                             &h_envChan[0]->sbrCodeNoiseFloor,
                             h_envChan[0]->encEnvData.domain_vec_noise, 0,
                             (eData[0].frame_info->nEnvelopes > 1 ? 2 : 1), 0,
                             sbrBitstreamData->HeaderActive);

      break;
    case SBR_LEFT_RIGHT:
      sbrNoiseFloorLevelsQuantisation(eData[0].noise_level, eData[0].noiseFloor,
                                      0);

      FDKsbrEnc_codeEnvelope(eData[0].noise_level, fData->res,
                             &h_envChan[0]->sbrCodeNoiseFloor,
                             h_envChan[0]->encEnvData.domain_vec_noise, 0,
                             (eData[0].frame_info->nEnvelopes > 1 ? 2 : 1), 0,
                             sbrBitstreamData->HeaderActive);

      sbrNoiseFloorLevelsQuantisation(eData[1].noise_level, eData[1].noiseFloor,
                                      0);

      FDKsbrEnc_codeEnvelope(eData[1].noise_level, fData->res,
                             &h_envChan[1]->sbrCodeNoiseFloor,
                             h_envChan[1]->encEnvData.domain_vec_noise, 0,
                             (eData[1].frame_info->nEnvelopes > 1 ? 2 : 1), 0,
                             sbrBitstreamData->HeaderActive);

      break;

    case SBR_COUPLING:
      coupleNoiseFloor(eData[0].noiseFloor, eData[1].noiseFloor);

      sbrNoiseFloorLevelsQuantisation(eData[0].noise_level, eData[0].noiseFloor,
                                      0);

      FDKsbrEnc_codeEnvelope(eData[0].noise_level, fData->res,
                             &h_envChan[0]->sbrCodeNoiseFloor,
                             h_envChan[0]->encEnvData.domain_vec_noise, 1,
                             (eData[0].frame_info->nEnvelopes > 1 ? 2 : 1), 0,
                             sbrBitstreamData->HeaderActive);

      sbrNoiseFloorLevelsQuantisation(eData[1].noise_level, eData[1].noiseFloor,
                                      1);

      FDKsbrEnc_codeEnvelope(eData[1].noise_level, fData->res,
                             &h_envChan[1]->sbrCodeNoiseFloor,
                             h_envChan[1]->encEnvData.domain_vec_noise, 1,
                             (eData[1].frame_info->nEnvelopes > 1 ? 2 : 1), 1,
                             sbrBitstreamData->HeaderActive);

      break;
    case SBR_SWITCH_LRC:
      sbrNoiseFloorLevelsQuantisation(eData[0].noise_level, eData[0].noiseFloor,
                                      0);
      sbrNoiseFloorLevelsQuantisation(eData[1].noise_level, eData[1].noiseFloor,
                                      0);
      coupleNoiseFloor(eData[0].noiseFloor, eData[1].noiseFloor);
      sbrNoiseFloorLevelsQuantisation(eData[0].noise_level_coupling,
                                      eData[0].noiseFloor, 0);
      sbrNoiseFloorLevelsQuantisation(eData[1].noise_level_coupling,
                                      eData[1].noiseFloor, 1);
      break;
  }

  switch (stereoMode) {
    case SBR_MONO:
      sbrHeaderData->coupling = 0;
      h_envChan[0]->encEnvData.balance = 0;
      FDKsbrEnc_codeEnvelope(
          eData[0].sfb_nrg, eData[0].frame_info->freqRes,
          &h_envChan[0]->sbrCodeEnvelope, h_envChan[0]->encEnvData.domain_vec,
          sbrHeaderData->coupling, eData[0].frame_info->nEnvelopes, 0,
          sbrBitstreamData->HeaderActive);
      break;
    case SBR_LEFT_RIGHT:
      sbrHeaderData->coupling = 0;

      h_envChan[0]->encEnvData.balance = 0;
      h_envChan[1]->encEnvData.balance = 0;

      FDKsbrEnc_codeEnvelope(
          eData[0].sfb_nrg, eData[0].frame_info->freqRes,
          &h_envChan[0]->sbrCodeEnvelope, h_envChan[0]->encEnvData.domain_vec,
          sbrHeaderData->coupling, eData[0].frame_info->nEnvelopes, 0,
          sbrBitstreamData->HeaderActive);
      FDKsbrEnc_codeEnvelope(
          eData[1].sfb_nrg, eData[1].frame_info->freqRes,
          &h_envChan[1]->sbrCodeEnvelope, h_envChan[1]->encEnvData.domain_vec,
          sbrHeaderData->coupling, eData[1].frame_info->nEnvelopes, 0,
          sbrBitstreamData->HeaderActive);
      break;
    case SBR_COUPLING:
      sbrHeaderData->coupling = 1;
      h_envChan[0]->encEnvData.balance = 0;
      h_envChan[1]->encEnvData.balance = 1;

      FDKsbrEnc_codeEnvelope(
          eData[0].sfb_nrg, eData[0].frame_info->freqRes,
          &h_envChan[0]->sbrCodeEnvelope, h_envChan[0]->encEnvData.domain_vec,
          sbrHeaderData->coupling, eData[0].frame_info->nEnvelopes, 0,
          sbrBitstreamData->HeaderActive);
      FDKsbrEnc_codeEnvelope(
          eData[1].sfb_nrg, eData[1].frame_info->freqRes,
          &h_envChan[1]->sbrCodeEnvelope, h_envChan[1]->encEnvData.domain_vec,
          sbrHeaderData->coupling, eData[1].frame_info->nEnvelopes, 1,
          sbrBitstreamData->HeaderActive);
      break;
    case SBR_SWITCH_LRC: {
      INT payloadbitsLR;
      INT payloadbitsCOUPLING;

      SCHAR sfbNrgPrevTemp[MAX_NUM_CHANNELS][MAX_FREQ_COEFFS];
      SCHAR noisePrevTemp[MAX_NUM_CHANNELS][MAX_NUM_NOISE_COEFFS];
      INT upDateNrgTemp[MAX_NUM_CHANNELS];
      INT upDateNoiseTemp[MAX_NUM_CHANNELS];
      INT domainVecTemp[MAX_NUM_CHANNELS][MAX_ENVELOPES];
      INT domainVecNoiseTemp[MAX_NUM_CHANNELS][MAX_ENVELOPES];

      INT tempFlagRight = 0;
      INT tempFlagLeft = 0;

      for (ch = 0; ch < nChannels; ch++) {
        FDKmemcpy(sfbNrgPrevTemp[ch],
                  h_envChan[ch]->sbrCodeEnvelope.sfb_nrg_prev,
                  MAX_FREQ_COEFFS * sizeof(SCHAR));

        FDKmemcpy(noisePrevTemp[ch],
                  h_envChan[ch]->sbrCodeNoiseFloor.sfb_nrg_prev,
                  MAX_NUM_NOISE_COEFFS * sizeof(SCHAR));

        upDateNrgTemp[ch] = h_envChan[ch]->sbrCodeEnvelope.upDate;
        upDateNoiseTemp[ch] = h_envChan[ch]->sbrCodeNoiseFloor.upDate;

        if (sbrHeaderData->prev_coupling) {
          h_envChan[ch]->sbrCodeEnvelope.upDate = 0;
          h_envChan[ch]->sbrCodeNoiseFloor.upDate = 0;
        }
      }

      FDKsbrEnc_codeEnvelope(eData[0].sfb_nrg, eData[0].frame_info->freqRes,
                             &h_envChan[0]->sbrCodeEnvelope,
                             h_envChan[0]->encEnvData.domain_vec, 0,
                             eData[0].frame_info->nEnvelopes, 0,
                             sbrBitstreamData->HeaderActive);
      FDKsbrEnc_codeEnvelope(eData[1].sfb_nrg, eData[1].frame_info->freqRes,
                             &h_envChan[1]->sbrCodeEnvelope,
                             h_envChan[1]->encEnvData.domain_vec, 0,
                             eData[1].frame_info->nEnvelopes, 0,
                             sbrBitstreamData->HeaderActive);

      c = 0;
      for (i = 0; i < eData[0].nEnvelopes; i++) {
        for (j = 0; j < h_envChan[0]->encEnvData.noScfBands[i]; j++) {
          h_envChan[0]->encEnvData.ienvelope[i][j] = eData[0].sfb_nrg[c];
          h_envChan[1]->encEnvData.ienvelope[i][j] = eData[1].sfb_nrg[c];
          c++;
        }
      }

      FDKsbrEnc_codeEnvelope(eData[0].noise_level, fData->res,
                             &h_envChan[0]->sbrCodeNoiseFloor,
                             h_envChan[0]->encEnvData.domain_vec_noise, 0,
                             (eData[0].frame_info->nEnvelopes > 1 ? 2 : 1), 0,
                             sbrBitstreamData->HeaderActive);

      for (i = 0; i < MAX_NUM_NOISE_VALUES; i++)
        h_envChan[0]->encEnvData.sbr_noise_levels[i] = eData[0].noise_level[i];

      FDKsbrEnc_codeEnvelope(eData[1].noise_level, fData->res,
                             &h_envChan[1]->sbrCodeNoiseFloor,
                             h_envChan[1]->encEnvData.domain_vec_noise, 0,
                             (eData[1].frame_info->nEnvelopes > 1 ? 2 : 1), 0,
                             sbrBitstreamData->HeaderActive);

      for (i = 0; i < MAX_NUM_NOISE_VALUES; i++)
        h_envChan[1]->encEnvData.sbr_noise_levels[i] = eData[1].noise_level[i];

      sbrHeaderData->coupling = 0;
      h_envChan[0]->encEnvData.balance = 0;
      h_envChan[1]->encEnvData.balance = 0;

      payloadbitsLR = FDKsbrEnc_CountSbrChannelPairElement(
          sbrHeaderData, hParametricStereo, sbrBitstreamData,
          &h_envChan[0]->encEnvData, &h_envChan[1]->encEnvData, hCmonData,
          h_con->sbrSyntaxFlags);

      for (ch = 0; ch < nChannels; ch++) {
        INT itmp;
        for (i = 0; i < MAX_FREQ_COEFFS; i++) {
          itmp = h_envChan[ch]->sbrCodeEnvelope.sfb_nrg_prev[i];
          h_envChan[ch]->sbrCodeEnvelope.sfb_nrg_prev[i] =
              sfbNrgPrevTemp[ch][i];
          sfbNrgPrevTemp[ch][i] = itmp;
        }
        for (i = 0; i < MAX_NUM_NOISE_COEFFS; i++) {
          itmp = h_envChan[ch]->sbrCodeNoiseFloor.sfb_nrg_prev[i];
          h_envChan[ch]->sbrCodeNoiseFloor.sfb_nrg_prev[i] =
              noisePrevTemp[ch][i];
          noisePrevTemp[ch][i] = itmp;
        }
        itmp = h_envChan[ch]->sbrCodeEnvelope.upDate;
        h_envChan[ch]->sbrCodeEnvelope.upDate = upDateNrgTemp[ch];
        upDateNrgTemp[ch] = itmp;

        itmp = h_envChan[ch]->sbrCodeNoiseFloor.upDate;
        h_envChan[ch]->sbrCodeNoiseFloor.upDate = upDateNoiseTemp[ch];
        upDateNoiseTemp[ch] = itmp;

        FDKmemcpy(domainVecTemp[ch], h_envChan[ch]->encEnvData.domain_vec,
                  sizeof(INT) * MAX_ENVELOPES);
        FDKmemcpy(domainVecNoiseTemp[ch],
                  h_envChan[ch]->encEnvData.domain_vec_noise,
                  sizeof(INT) * MAX_ENVELOPES);

        if (!sbrHeaderData->prev_coupling) {
          h_envChan[ch]->sbrCodeEnvelope.upDate = 0;
          h_envChan[ch]->sbrCodeNoiseFloor.upDate = 0;
        }
      }

      FDKsbrEnc_codeEnvelope(
          eData[0].sfb_nrg_coupling, eData[0].frame_info->freqRes,
          &h_envChan[0]->sbrCodeEnvelope, h_envChan[0]->encEnvData.domain_vec,
          1, eData[0].frame_info->nEnvelopes, 0,
          sbrBitstreamData->HeaderActive);

      FDKsbrEnc_codeEnvelope(
          eData[1].sfb_nrg_coupling, eData[1].frame_info->freqRes,
          &h_envChan[1]->sbrCodeEnvelope, h_envChan[1]->encEnvData.domain_vec,
          1, eData[1].frame_info->nEnvelopes, 1,
          sbrBitstreamData->HeaderActive);

      c = 0;
      for (i = 0; i < eData[0].nEnvelopes; i++) {
        for (j = 0; j < h_envChan[0]->encEnvData.noScfBands[i]; j++) {
          h_envChan[0]->encEnvData.ienvelope[i][j] =
              eData[0].sfb_nrg_coupling[c];
          h_envChan[1]->encEnvData.ienvelope[i][j] =
              eData[1].sfb_nrg_coupling[c];
          c++;
        }
      }

      FDKsbrEnc_codeEnvelope(eData[0].noise_level_coupling, fData->res,
                             &h_envChan[0]->sbrCodeNoiseFloor,
                             h_envChan[0]->encEnvData.domain_vec_noise, 1,
                             (eData[0].frame_info->nEnvelopes > 1 ? 2 : 1), 0,
                             sbrBitstreamData->HeaderActive);

      for (i = 0; i < MAX_NUM_NOISE_VALUES; i++)
        h_envChan[0]->encEnvData.sbr_noise_levels[i] =
            eData[0].noise_level_coupling[i];

      FDKsbrEnc_codeEnvelope(eData[1].noise_level_coupling, fData->res,
                             &h_envChan[1]->sbrCodeNoiseFloor,
                             h_envChan[1]->encEnvData.domain_vec_noise, 1,
                             (eData[1].frame_info->nEnvelopes > 1 ? 2 : 1), 1,
                             sbrBitstreamData->HeaderActive);

      for (i = 0; i < MAX_NUM_NOISE_VALUES; i++)
        h_envChan[1]->encEnvData.sbr_noise_levels[i] =
            eData[1].noise_level_coupling[i];

      sbrHeaderData->coupling = 1;

      h_envChan[0]->encEnvData.balance = 0;
      h_envChan[1]->encEnvData.balance = 1;

      tempFlagLeft = h_envChan[0]->encEnvData.addHarmonicFlag;
      tempFlagRight = h_envChan[1]->encEnvData.addHarmonicFlag;

      payloadbitsCOUPLING = FDKsbrEnc_CountSbrChannelPairElement(
          sbrHeaderData, hParametricStereo, sbrBitstreamData,
          &h_envChan[0]->encEnvData, &h_envChan[1]->encEnvData, hCmonData,
          h_con->sbrSyntaxFlags);

      h_envChan[0]->encEnvData.addHarmonicFlag = tempFlagLeft;
      h_envChan[1]->encEnvData.addHarmonicFlag = tempFlagRight;

      if (payloadbitsCOUPLING < payloadbitsLR) {
        for (ch = 0; ch < nChannels; ch++) {
          SBR_ENV_TEMP_DATA *ed = &eData[ch];
          FDKmemcpy(ed->sfb_nrg, ed->sfb_nrg_coupling,
                    MAX_NUM_ENVELOPE_VALUES * sizeof(SCHAR));
          FDKmemcpy(ed->noise_level, ed->noise_level_coupling,
                    MAX_NUM_NOISE_VALUES * sizeof(SCHAR));
        }

        sbrHeaderData->coupling = 1;
        h_envChan[0]->encEnvData.balance = 0;
        h_envChan[1]->encEnvData.balance = 1;
      } else {
        for (ch = 0; ch < nChannels; ch++) {
          FDKmemcpy(h_envChan[ch]->sbrCodeEnvelope.sfb_nrg_prev,
                    sfbNrgPrevTemp[ch], MAX_FREQ_COEFFS * sizeof(SCHAR));

          h_envChan[ch]->sbrCodeEnvelope.upDate = upDateNrgTemp[ch];

          FDKmemcpy(h_envChan[ch]->sbrCodeNoiseFloor.sfb_nrg_prev,
                    noisePrevTemp[ch], MAX_NUM_NOISE_COEFFS * sizeof(SCHAR));

          FDKmemcpy(h_envChan[ch]->encEnvData.domain_vec, domainVecTemp[ch],
                    sizeof(INT) * MAX_ENVELOPES);
          FDKmemcpy(h_envChan[ch]->encEnvData.domain_vec_noise,
                    domainVecNoiseTemp[ch], sizeof(INT) * MAX_ENVELOPES);

          h_envChan[ch]->sbrCodeNoiseFloor.upDate = upDateNoiseTemp[ch];
        }

        sbrHeaderData->coupling = 0;
        h_envChan[0]->encEnvData.balance = 0;
        h_envChan[1]->encEnvData.balance = 0;
      }
    } break;
  }

  if (stereoMode == SBR_MONO) {
    if (h_envChan[0]->encEnvData.domain_vec[0] == TIME)
      h_envChan[0]->sbrCodeEnvelope.dF_edge_incr_fac++;
    else
      h_envChan[0]->sbrCodeEnvelope.dF_edge_incr_fac = 0;
  } else {
    if (h_envChan[0]->encEnvData.domain_vec[0] == TIME ||
        h_envChan[1]->encEnvData.domain_vec[0] == TIME) {
      h_envChan[0]->sbrCodeEnvelope.dF_edge_incr_fac++;
      h_envChan[1]->sbrCodeEnvelope.dF_edge_incr_fac++;
    } else {
      h_envChan[0]->sbrCodeEnvelope.dF_edge_incr_fac = 0;
      h_envChan[1]->sbrCodeEnvelope.dF_edge_incr_fac = 0;
    }
  }

  for (ch = 0; ch < nChannels; ch++) {
    SBR_ENV_TEMP_DATA *ed = &eData[ch];
    c = 0;
    for (i = 0; i < ed->nEnvelopes; i++) {
      for (j = 0; j < h_envChan[ch]->encEnvData.noScfBands[i]; j++) {
        h_envChan[ch]->encEnvData.ienvelope[i][j] = ed->sfb_nrg[c];

        c++;
      }
    }
    for (i = 0; i < MAX_NUM_NOISE_VALUES; i++) {
      h_envChan[ch]->encEnvData.sbr_noise_levels[i] = ed->noise_level[i];
    }
  }

  if (nChannels == 2) {
    FDKsbrEnc_WriteEnvChannelPairElement(
        sbrHeaderData, hParametricStereo, sbrBitstreamData,
        &h_envChan[0]->encEnvData, &h_envChan[1]->encEnvData, hCmonData,
        h_con->sbrSyntaxFlags);
  } else {
    FDKsbrEnc_WriteEnvSingleChannelElement(
        sbrHeaderData, hParametricStereo, sbrBitstreamData,
        &h_envChan[0]->encEnvData, hCmonData, h_con->sbrSyntaxFlags);
  }

  for (ch = 0; ch < nChannels; ch++) {
    int YBufferLength = h_envChan[ch]->sbrExtractEnvelope.no_cols >>
                        h_envChan[ch]->sbrExtractEnvelope.YBufferSzShift;
    for (i = 0; i < h_envChan[ch]->sbrExtractEnvelope.YBufferWriteOffset; i++) {
      FDKmemcpy(h_envChan[ch]->sbrExtractEnvelope.YBuffer[i],
                h_envChan[ch]->sbrExtractEnvelope.YBuffer[i + YBufferLength],
                sizeof(FIXP_DBL) * 64);
    }
    h_envChan[ch]->sbrExtractEnvelope.YBufferScale[0] =
        h_envChan[ch]->sbrExtractEnvelope.YBufferScale[1];
  }

  sbrHeaderData->prev_coupling = sbrHeaderData->coupling;
}

INT FDKsbrEnc_CreateExtractSbrEnvelope(HANDLE_SBR_EXTRACT_ENVELOPE hSbrCut,
                                       INT channel, INT chInEl,
                                       UCHAR *dynamic_RAM) {
  INT i;
  FIXP_DBL *rBuffer, *iBuffer;
  INT n;
  FIXP_DBL *YBufferDyn;

  FDKmemclear(hSbrCut, sizeof(SBR_EXTRACT_ENVELOPE));

  if (NULL == (hSbrCut->p_YBuffer = GetRam_Sbr_envYBuffer(channel))) {
    goto bail;
  }

  for (i = 0; i < (32 >> 1); i++) {
    hSbrCut->YBuffer[i] = hSbrCut->p_YBuffer + (i * 64);
  }
  YBufferDyn = GetRam_Sbr_envYBuffer(chInEl, dynamic_RAM);
  for (n = 0; i < 32; i++, n++) {
    hSbrCut->YBuffer[i] = YBufferDyn + (n * 64);
  }

  rBuffer = GetRam_Sbr_envRBuffer(0, dynamic_RAM);
  iBuffer = GetRam_Sbr_envIBuffer(0, dynamic_RAM);

  for (i = 0; i < 32; i++) {
    hSbrCut->rBuffer[i] = rBuffer + (i * 64);
    hSbrCut->iBuffer[i] = iBuffer + (i * 64);
  }

  return 0;

bail:
  FDKsbrEnc_deleteExtractSbrEnvelope(hSbrCut);

  return -1;
}

INT FDKsbrEnc_InitExtractSbrEnvelope(HANDLE_SBR_EXTRACT_ENVELOPE hSbrCut,
                                     int no_cols, int no_rows, int start_index,
                                     int time_slots, int time_step,
                                     int tran_off, ULONG statesInitFlag,
                                     int chInEl, UCHAR *dynamic_RAM,
                                     UINT sbrSyntaxFlags) {
  int YBufferLength, rBufferLength;
  int i;

  if (sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY) {
    int off = TRANSIENT_OFFSET_LD;
    hSbrCut->YBufferWriteOffset = (no_cols >> 1) + off * time_step;
  } else {
    hSbrCut->YBufferWriteOffset = tran_off * time_step;
  }
  hSbrCut->rBufferReadOffset = 0;

  YBufferLength = hSbrCut->YBufferWriteOffset + no_cols;
  rBufferLength = no_cols;

  hSbrCut->pre_transient_info[0] = 0;
  hSbrCut->pre_transient_info[1] = 0;

  hSbrCut->no_cols = no_cols;
  hSbrCut->no_rows = no_rows;
  hSbrCut->start_index = start_index;

  hSbrCut->time_slots = time_slots;
  hSbrCut->time_step = time_step;

  FDK_ASSERT(no_rows <= 64);

  if (time_step >= 2)
    hSbrCut->YBufferSzShift = 1;
  else
    hSbrCut->YBufferSzShift = 0;

  YBufferLength >>= hSbrCut->YBufferSzShift;
  hSbrCut->YBufferWriteOffset >>= hSbrCut->YBufferSzShift;

  FDK_ASSERT(YBufferLength <= 32);

  FIXP_DBL *YBufferDyn = GetRam_Sbr_envYBuffer(chInEl, dynamic_RAM);
  INT n = 0;
  for (i = (32 >> 1); i < 32; i++, n++) {
    hSbrCut->YBuffer[i] = YBufferDyn + (n * 64);
  }

  if (statesInitFlag) {
    for (i = 0; i < YBufferLength; i++) {
      FDKmemclear(hSbrCut->YBuffer[i], 64 * sizeof(FIXP_DBL));
    }
  }

  for (i = 0; i < rBufferLength; i++) {
    FDKmemclear(hSbrCut->rBuffer[i], 64 * sizeof(FIXP_DBL));
    FDKmemclear(hSbrCut->iBuffer[i], 64 * sizeof(FIXP_DBL));
  }

  FDKmemclear(hSbrCut->envelopeCompensation, sizeof(UCHAR) * MAX_FREQ_COEFFS);

  if (statesInitFlag) {
    hSbrCut->YBufferScale[0] = hSbrCut->YBufferScale[1] = FRACT_BITS - 1;
  }

  return (0);
}

void FDKsbrEnc_deleteExtractSbrEnvelope(HANDLE_SBR_EXTRACT_ENVELOPE hSbrCut) {
  if (hSbrCut) {
    FreeRam_Sbr_envYBuffer(&hSbrCut->p_YBuffer);
  }
}

INT FDKsbrEnc_GetEnvEstDelay(HANDLE_SBR_EXTRACT_ENVELOPE hSbr) {
  return hSbr->no_rows *
         ((hSbr->YBufferWriteOffset) *
              2
          - hSbr->rBufferReadOffset);
}
