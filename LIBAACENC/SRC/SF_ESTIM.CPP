/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Scale factor estimation

*******************************************************************************/

#include "sf_estim.h"
#include "aacEnc_rom.h"
#include "quantize.h"
#include "bit_cnt.h"

#ifdef __arm__
#endif

#define UPCOUNT_LIMIT 1
#define AS_PE_FAC_SHIFT 7
#define DIST_FAC_SHIFT 3
#define AS_PE_FAC_FLOAT (float)(1 << AS_PE_FAC_SHIFT)
static const INT MAX_SCF_DELTA = 60;

static const FIXP_DBL PE_C1 = FL2FXCONST_DBL(
    3.0f / AS_PE_FAC_FLOAT);
static const FIXP_DBL PE_C2 = FL2FXCONST_DBL(
    1.3219281f / AS_PE_FAC_FLOAT);
static const FIXP_DBL PE_C3 = FL2FXCONST_DBL(0.5593573f);

static void FDKaacEnc_FDKaacEnc_CalcFormFactorChannel(
    FIXP_DBL *RESTRICT sfbFormFactorLdData,
    PSY_OUT_CHANNEL *RESTRICT psyOutChan) {
  INT j, sfb, sfbGrp;
  FIXP_DBL formFactor;

  int tmp0 = psyOutChan->sfbCnt;
  int tmp1 = psyOutChan->maxSfbPerGroup;
  int step = psyOutChan->sfbPerGroup;
  for (sfbGrp = 0; sfbGrp < tmp0; sfbGrp += step) {
    for (sfb = 0; sfb < tmp1; sfb++) {
      formFactor = FL2FXCONST_DBL(0.0f);
      for (j = psyOutChan->sfbOffsets[sfbGrp + sfb];
           j < psyOutChan->sfbOffsets[sfbGrp + sfb + 1]; j++) {
        formFactor +=
            sqrtFixp(fixp_abs(psyOutChan->mdctSpectrum[j])) >> FORM_FAC_SHIFT;
      }
      sfbFormFactorLdData[sfbGrp + sfb] = CalcLdData(formFactor);
    }
    for (; sfb < psyOutChan->sfbPerGroup; sfb++) {
      sfbFormFactorLdData[sfbGrp + sfb] = FL2FXCONST_DBL(-1.0f);
    }
  }
}

void FDKaacEnc_CalcFormFactor(QC_OUT_CHANNEL *qcOutChannel[(2)],
                              PSY_OUT_CHANNEL *psyOutChannel[(2)],
                              const INT nChannels) {
  INT j;
  for (j = 0; j < nChannels; j++) {
    FDKaacEnc_FDKaacEnc_CalcFormFactorChannel(
        qcOutChannel[j]->sfbFormFactorLdData, psyOutChannel[j]);
  }
}

static void FDKaacEnc_calcSfbRelevantLines(
    const FIXP_DBL *const sfbFormFactorLdData,
    const FIXP_DBL *const sfbEnergyLdData,
    const FIXP_DBL *const sfbThresholdLdData, const INT *const sfbOffsets,
    const INT sfbCnt, const INT sfbPerGroup, const INT maxSfbPerGroup,
    FIXP_DBL *sfbNRelevantLines) {
  INT sfbOffs, sfb;
  FIXP_DBL sfbWidthLdData;
  FIXP_DBL asPeFacLdData =
      FL2FXCONST_DBL(0.109375);
  FIXP_DBL accu;

  FDKmemclear(sfbNRelevantLines, sfbCnt * sizeof(FIXP_DBL));

  for (sfbOffs = 0; sfbOffs < sfbCnt; sfbOffs += sfbPerGroup) {
    for (sfb = 0; sfb < maxSfbPerGroup; sfb++) {
      if ((FIXP_DBL)sfbEnergyLdData[sfbOffs + sfb] >
          (FIXP_DBL)sfbThresholdLdData[sfbOffs + sfb]) {
        INT sfbWidth =
            sfbOffsets[sfbOffs + sfb + 1] - sfbOffsets[sfbOffs + sfb];

        sfbWidthLdData =
            (FIXP_DBL)(sfbWidth << (DFRACT_BITS - 1 - AS_PE_FAC_SHIFT));
        sfbWidthLdData = CalcLdData(sfbWidthLdData);

        accu = sfbEnergyLdData[sfbOffs + sfb] - sfbWidthLdData - asPeFacLdData;
        accu = sfbFormFactorLdData[sfbOffs + sfb] - (accu >> 2);

        sfbNRelevantLines[sfbOffs + sfb] = CalcInvLdData(accu) >> 1;
      }
    }
  }
}

static FIXP_DBL FDKaacEnc_countSingleScfBits(INT scf, INT scfLeft,
                                             INT scfRight) {
  FIXP_DBL scfBitsFract;

  scfBitsFract = (FIXP_DBL)(FDKaacEnc_bitCountScalefactorDelta(scfLeft - scf) +
                            FDKaacEnc_bitCountScalefactorDelta(scf - scfRight));

  scfBitsFract = scfBitsFract << (DFRACT_BITS - 1 - (2 * AS_PE_FAC_SHIFT));

  return scfBitsFract;
}

static FIXP_DBL FDKaacEnc_calcSingleSpecPe(INT scf, FIXP_DBL sfbConstPePart,
                                           FIXP_DBL nLines) {
  FIXP_DBL specPe = FL2FXCONST_DBL(0.0f);
  FIXP_DBL ldRatio;
  FIXP_DBL scfFract;

  scfFract = (FIXP_DBL)(scf << (DFRACT_BITS - 1 - AS_PE_FAC_SHIFT));

  ldRatio = sfbConstPePart - fMult(FL2FXCONST_DBL(0.375f), scfFract);

  if (ldRatio >= PE_C1) {
    specPe = fMult(FL2FXCONST_DBL(0.7f), fMult(nLines, ldRatio));
  } else {
    specPe = fMult(FL2FXCONST_DBL(0.7f),
                   fMult(nLines, (PE_C2 + fMult(PE_C3, ldRatio))));
  }

  return specPe;
}

static FIXP_DBL FDKaacEnc_countScfBitsDiff(INT *scfOld, INT *scfNew, INT sfbCnt,
                                           INT startSfb, INT stopSfb) {
  FIXP_DBL scfBitsFract;
  INT scfBitsDiff = 0;
  INT sfb = 0, sfbLast;
  INT sfbPrev, sfbNext;

  sfbLast = startSfb;
  while ((sfbLast < stopSfb) && (scfOld[sfbLast] == FDK_INT_MIN)) sfbLast++;
  sfbPrev = startSfb - 1;
  while ((sfbPrev >= 0) && (scfOld[sfbPrev] == FDK_INT_MIN)) sfbPrev--;
  if (sfbPrev >= 0)
    scfBitsDiff +=
        FDKaacEnc_bitCountScalefactorDelta(scfNew[sfbPrev] - scfNew[sfbLast]) -
        FDKaacEnc_bitCountScalefactorDelta(scfOld[sfbPrev] - scfOld[sfbLast]);
  for (sfb = sfbLast + 1; sfb < stopSfb; sfb++) {
    if (scfOld[sfb] != FDK_INT_MIN) {
      scfBitsDiff +=
          FDKaacEnc_bitCountScalefactorDelta(scfNew[sfbLast] - scfNew[sfb]) -
          FDKaacEnc_bitCountScalefactorDelta(scfOld[sfbLast] - scfOld[sfb]);
      sfbLast = sfb;
    }
  }
  sfbNext = stopSfb;
  while ((sfbNext < sfbCnt) && (scfOld[sfbNext] == FDK_INT_MIN)) sfbNext++;
  if (sfbNext < sfbCnt)
    scfBitsDiff +=
        FDKaacEnc_bitCountScalefactorDelta(scfNew[sfbLast] - scfNew[sfbNext]) -
        FDKaacEnc_bitCountScalefactorDelta(scfOld[sfbLast] - scfOld[sfbNext]);

  scfBitsFract =
      (FIXP_DBL)(scfBitsDiff << (DFRACT_BITS - 1 - (2 * AS_PE_FAC_SHIFT)));

  return scfBitsFract;
}

static FIXP_DBL FDKaacEnc_calcSpecPeDiff(
    PSY_OUT_CHANNEL *psyOutChan, QC_OUT_CHANNEL *qcOutChannel, INT *scfOld,
    INT *scfNew, FIXP_DBL *sfbConstPePart, FIXP_DBL *sfbFormFactorLdData,
    FIXP_DBL *sfbNRelevantLines, INT startSfb, INT stopSfb) {
  FIXP_DBL specPeDiff = FL2FXCONST_DBL(0.0f);
  FIXP_DBL scfFract = FL2FXCONST_DBL(0.0f);
  INT sfb;

  for (sfb = startSfb; sfb < stopSfb; sfb++) {
    if (scfOld[sfb] != FDK_INT_MIN) {
      FIXP_DBL ldRatioOld, ldRatioNew, pOld, pNew;

      if (sfbConstPePart[sfb] == (FIXP_DBL)FDK_INT_MIN)
        sfbConstPePart[sfb] =
            ((psyOutChan->sfbEnergyLdData[sfb] - sfbFormFactorLdData[sfb] -
              FL2FXCONST_DBL(0.09375f)) >>
             1) +
            FL2FXCONST_DBL(0.02152255861f);

      scfFract = (FIXP_DBL)(scfOld[sfb] << (DFRACT_BITS - 1 - AS_PE_FAC_SHIFT));
      ldRatioOld =
          sfbConstPePart[sfb] - fMult(FL2FXCONST_DBL(0.375f), scfFract);

      scfFract = (FIXP_DBL)(scfNew[sfb] << (DFRACT_BITS - 1 - AS_PE_FAC_SHIFT));
      ldRatioNew =
          sfbConstPePart[sfb] - fMult(FL2FXCONST_DBL(0.375f), scfFract);

      if (ldRatioOld >= PE_C1)
        pOld = ldRatioOld;
      else
        pOld = PE_C2 + fMult(PE_C3, ldRatioOld);

      if (ldRatioNew >= PE_C1)
        pNew = ldRatioNew;
      else
        pNew = PE_C2 + fMult(PE_C3, ldRatioNew);

      specPeDiff += fMult(FL2FXCONST_DBL(0.7f),
                          fMult(sfbNRelevantLines[sfb], (pNew - pOld)));
    }
  }

  return specPeDiff;
}

static INT FDKaacEnc_improveScf(const FIXP_DBL *spec, SHORT *quantSpec,
                                SHORT *quantSpecTmp, INT sfbWidth,
                                FIXP_DBL threshLdData, INT scf, INT minScf,
                                FIXP_DBL *distLdData, INT *minScfCalculated,
                                INT dZoneQuantEnable) {
  FIXP_DBL sfbDistLdData;
  INT scfBest = scf;
  INT k;
  FIXP_DBL distFactorLdData = FL2FXCONST_DBL(-0.0050301265);

  sfbDistLdData =
      FDKaacEnc_calcSfbDist(spec, quantSpec, sfbWidth, scf, dZoneQuantEnable);
  *minScfCalculated = scf;
  if (sfbDistLdData > (threshLdData - distFactorLdData)) {
    INT scfEstimated = scf;
    FIXP_DBL sfbDistBestLdData = sfbDistLdData;
    INT cnt;
    cnt = 0;

    while ((sfbDistLdData > (threshLdData - distFactorLdData)) &&
           (cnt++ < UPCOUNT_LIMIT)) {
      scf++;
      sfbDistLdData = FDKaacEnc_calcSfbDist(spec, quantSpecTmp, sfbWidth, scf,
                                            dZoneQuantEnable);

      if (sfbDistLdData < sfbDistBestLdData) {
        scfBest = scf;
        sfbDistBestLdData = sfbDistLdData;
        for (k = 0; k < sfbWidth; k++) quantSpec[k] = quantSpecTmp[k];
      }
    }
    cnt = 0;
    scf = scfEstimated;
    sfbDistLdData = sfbDistBestLdData;
    while ((sfbDistLdData > (threshLdData - distFactorLdData)) && (cnt++ < 1) &&
           (scf > minScf)) {
      scf--;
      sfbDistLdData = FDKaacEnc_calcSfbDist(spec, quantSpecTmp, sfbWidth, scf,
                                            dZoneQuantEnable);

      if (sfbDistLdData < sfbDistBestLdData) {
        scfBest = scf;
        sfbDistBestLdData = sfbDistLdData;
        for (k = 0; k < sfbWidth; k++) quantSpec[k] = quantSpecTmp[k];
      }
      *minScfCalculated = scf;
    }
    *distLdData = sfbDistBestLdData;
  } else {
    FIXP_DBL sfbDistBestLdData = sfbDistLdData;
    FIXP_DBL sfbDistAllowedLdData =
        fixMin(sfbDistLdData - distFactorLdData, threshLdData);
    int cnt;
    for (cnt = 0; cnt < UPCOUNT_LIMIT; cnt++) {
      scf++;
      sfbDistLdData = FDKaacEnc_calcSfbDist(spec, quantSpecTmp, sfbWidth, scf,
                                            dZoneQuantEnable);

      if (sfbDistLdData < sfbDistAllowedLdData) {
        *minScfCalculated = scfBest + 1;
        scfBest = scf;
        sfbDistBestLdData = sfbDistLdData;
        for (k = 0; k < sfbWidth; k++) quantSpec[k] = quantSpecTmp[k];
      }
    }
    *distLdData = sfbDistBestLdData;
  }

  return scfBest;
}

static void FDKaacEnc_assimilateSingleScf(
    const PSY_OUT_CHANNEL *psyOutChan, const QC_OUT_CHANNEL *qcOutChannel,
    SHORT *quantSpec, SHORT *quantSpecTmp, INT dZoneQuantEnable, INT *scf,
    const INT *minScf, FIXP_DBL *sfbDist, FIXP_DBL *sfbConstPePart,
    const FIXP_DBL *sfbFormFactorLdData, const FIXP_DBL *sfbNRelevantLines,
    INT *minScfCalculated, INT restartOnSuccess) {
  INT sfbLast, sfbAct, sfbNext;
  INT scfAct, *scfLast, *scfNext, scfMin, scfMax;
  INT sfbWidth, sfbOffs;
  FIXP_DBL enLdData;
  FIXP_DBL sfbPeOld, sfbPeNew;
  FIXP_DBL sfbDistNew;
  INT i, k;
  INT success = 0;
  FIXP_DBL deltaPe = FL2FXCONST_DBL(0.0f);
  FIXP_DBL deltaPeNew, deltaPeTmp;
  INT prevScfLast[MAX_GROUPED_SFB], prevScfNext[MAX_GROUPED_SFB];
  FIXP_DBL deltaPeLast[MAX_GROUPED_SFB];
  INT updateMinScfCalculated;

  for (i = 0; i < psyOutChan->sfbCnt; i++) {
    prevScfLast[i] = FDK_INT_MAX;
    prevScfNext[i] = FDK_INT_MAX;
    deltaPeLast[i] = (FIXP_DBL)FDK_INT_MAX;
  }

  sfbLast = -1;
  sfbAct = -1;
  sfbNext = -1;
  scfLast = 0;
  scfNext = 0;
  scfMin = FDK_INT_MAX;
  scfMax = FDK_INT_MAX;
  do {
    sfbNext++;
    while ((sfbNext < psyOutChan->sfbCnt) && (scf[sfbNext] == FDK_INT_MIN))
      sfbNext++;
    if ((sfbLast >= 0) && (sfbAct >= 0) && (sfbNext < psyOutChan->sfbCnt)) {
      scfAct = scf[sfbAct];
      scfLast = scf + sfbLast;
      scfNext = scf + sfbNext;
      scfMin = fixMin(*scfLast, *scfNext);
      scfMax = fixMax(*scfLast, *scfNext);
    } else if ((sfbLast == -1) && (sfbAct >= 0) &&
               (sfbNext < psyOutChan->sfbCnt)) {
      scfAct = scf[sfbAct];
      scfLast = &scfAct;
      scfNext = scf + sfbNext;
      scfMin = *scfNext;
      scfMax = *scfNext;
    } else if ((sfbLast >= 0) && (sfbAct >= 0) &&
               (sfbNext == psyOutChan->sfbCnt)) {
      scfAct = scf[sfbAct];
      scfLast = scf + sfbLast;
      scfNext = &scfAct;
      scfMin = *scfLast;
      scfMax = *scfLast;
    }
    if (sfbAct >= 0) scfMin = fixMax(scfMin, minScf[sfbAct]);

    if ((sfbAct >= 0) && (sfbLast >= 0 || sfbNext < psyOutChan->sfbCnt) &&
        (scfAct > scfMin) && (scfAct <= scfMin + MAX_SCF_DELTA) &&
        (scfAct >= scfMax - MAX_SCF_DELTA) &&
        (scfAct <=
         fixMin(scfMin, fixMin(*scfLast, *scfNext)) + MAX_SCF_DELTA) &&
        (*scfLast != prevScfLast[sfbAct] || *scfNext != prevScfNext[sfbAct] ||
         deltaPe < deltaPeLast[sfbAct])) {
      success = 0;

      sfbWidth =
          psyOutChan->sfbOffsets[sfbAct + 1] - psyOutChan->sfbOffsets[sfbAct];
      sfbOffs = psyOutChan->sfbOffsets[sfbAct];

      enLdData = qcOutChannel->sfbEnergyLdData[sfbAct];

      if (sfbConstPePart[sfbAct] == (FIXP_DBL)FDK_INT_MIN) {
        sfbConstPePart[sfbAct] = ((enLdData - sfbFormFactorLdData[sfbAct] -
                                   FL2FXCONST_DBL(0.09375f)) >>
                                  1) +
                                 FL2FXCONST_DBL(0.02152255861f);
      }

      sfbPeOld = FDKaacEnc_calcSingleSpecPe(scfAct, sfbConstPePart[sfbAct],
                                            sfbNRelevantLines[sfbAct]) +
                 FDKaacEnc_countSingleScfBits(scfAct, *scfLast, *scfNext);

      deltaPeNew = deltaPe;
      updateMinScfCalculated = 1;

      do {
        scfAct--;
        if (scfAct < minScfCalculated[sfbAct] &&
            scfAct >= scfMax - MAX_SCF_DELTA) {
          sfbPeNew = FDKaacEnc_calcSingleSpecPe(scfAct, sfbConstPePart[sfbAct],
                                                sfbNRelevantLines[sfbAct]) +
                     FDKaacEnc_countSingleScfBits(scfAct, *scfLast, *scfNext);

          deltaPeTmp = deltaPe + sfbPeNew - sfbPeOld;
          if (deltaPeTmp < FL2FXCONST_DBL(0.0006103515625f)) {
            sfbDistNew = FDKaacEnc_calcSfbDist(
                qcOutChannel->mdctSpectrum + sfbOffs, quantSpecTmp + sfbOffs,
                sfbWidth, scfAct, dZoneQuantEnable);

            if (sfbDistNew < sfbDist[sfbAct]) {
              scf[sfbAct] = scfAct;
              sfbDist[sfbAct] = sfbDistNew;

              for (k = 0; k < sfbWidth; k++)
                quantSpec[sfbOffs + k] = quantSpecTmp[sfbOffs + k];

              deltaPeNew = deltaPeTmp;
              success = 1;
            }
            if (updateMinScfCalculated) minScfCalculated[sfbAct] = scfAct;
          } else {
            updateMinScfCalculated = 0;
          }
        }
      } while (scfAct > scfMin);

      deltaPe = deltaPeNew;

      prevScfLast[sfbAct] = *scfLast;
      prevScfNext[sfbAct] = *scfNext;
      deltaPeLast[sfbAct] = deltaPe;
    }

    if (success && restartOnSuccess) {
      sfbLast = -1;
      sfbAct = -1;
      sfbNext = -1;
      scfLast = 0;
      scfNext = 0;
      scfMin = FDK_INT_MAX;
      scfMax = FDK_INT_MAX;
      success = 0;
    } else {
      sfbLast = sfbAct;
      sfbAct = sfbNext;
    }
  } while (sfbNext < psyOutChan->sfbCnt);
}

static void FDKaacEnc_assimilateMultipleScf(
    PSY_OUT_CHANNEL *psyOutChan, QC_OUT_CHANNEL *qcOutChannel, SHORT *quantSpec,
    SHORT *quantSpecTmp, INT dZoneQuantEnable, INT *scf, const INT *minScf,
    FIXP_DBL *sfbDist, FIXP_DBL *sfbConstPePart, FIXP_DBL *sfbFormFactorLdData,
    FIXP_DBL *sfbNRelevantLines) {
  INT sfb, startSfb, stopSfb;
  INT scfTmp[MAX_GROUPED_SFB], scfMin, scfMax, scfAct;
  INT possibleRegionFound;
  INT sfbWidth, sfbOffs, i, k;
  FIXP_DBL sfbDistNew[MAX_GROUPED_SFB], distOldSum, distNewSum;
  INT deltaScfBits;
  FIXP_DBL deltaSpecPe;
  FIXP_DBL deltaPe = FL2FXCONST_DBL(0.0f);
  FIXP_DBL deltaPeNew;
  INT sfbCnt = psyOutChan->sfbCnt;

  scfMin = FDK_INT_MAX;
  scfMax = FDK_INT_MIN;
  for (sfb = 0; sfb < sfbCnt; sfb++) {
    if (scf[sfb] != FDK_INT_MIN) {
      scfMin = fixMin(scfMin, scf[sfb]);
      scfMax = fixMax(scfMax, scf[sfb]);
    }
  }

  if (scfMax != FDK_INT_MIN && scfMax <= scfMin + MAX_SCF_DELTA) {
    scfAct = scfMax;

    do {
      scfAct--;
      for (i = 0; i < MAX_GROUPED_SFB; i++) scfTmp[i] = scf[i];
      stopSfb = 0;
      do {
        sfb = stopSfb;
        while (sfb < sfbCnt && (scf[sfb] == FDK_INT_MIN || scf[sfb] <= scfAct))
          sfb++;
        startSfb = sfb;
        sfb++;
        while (sfb < sfbCnt && (scf[sfb] == FDK_INT_MIN || scf[sfb] > scfAct))
          sfb++;
        stopSfb = sfb;

        possibleRegionFound = 0;
        if (startSfb < sfbCnt) {
          possibleRegionFound = 1;
          for (sfb = startSfb; sfb < stopSfb; sfb++) {
            if (scf[sfb] != FDK_INT_MIN)
              if (scfAct < minScf[sfb]) {
                possibleRegionFound = 0;
                break;
              }
          }
        }

        if (possibleRegionFound) {

          for (sfb = startSfb; sfb < stopSfb; sfb++) {
            if (scfTmp[sfb] != FDK_INT_MIN) scfTmp[sfb] = scfAct;
          }

          deltaScfBits = FDKaacEnc_countScfBitsDiff(scf, scfTmp, sfbCnt,
                                                    startSfb, stopSfb);

          deltaSpecPe = FDKaacEnc_calcSpecPeDiff(
              psyOutChan, qcOutChannel, scf, scfTmp, sfbConstPePart,
              sfbFormFactorLdData, sfbNRelevantLines, startSfb, stopSfb);

          deltaPeNew = deltaPe + (FIXP_DBL)deltaScfBits + deltaSpecPe;

          if (deltaPeNew < FL2FXCONST_DBL(0.0006103515625f)) {
            distOldSum = distNewSum = FL2FXCONST_DBL(0.0f);
            for (sfb = startSfb; sfb < stopSfb; sfb++) {
              if (scfTmp[sfb] != FDK_INT_MIN) {
                distOldSum += CalcInvLdData(sfbDist[sfb]) >> DIST_FAC_SHIFT;

                sfbWidth = psyOutChan->sfbOffsets[sfb + 1] -
                           psyOutChan->sfbOffsets[sfb];
                sfbOffs = psyOutChan->sfbOffsets[sfb];

                sfbDistNew[sfb] = FDKaacEnc_calcSfbDist(
                    qcOutChannel->mdctSpectrum + sfbOffs,
                    quantSpecTmp + sfbOffs, sfbWidth, scfAct, dZoneQuantEnable);

                if (sfbDistNew[sfb] > qcOutChannel->sfbThresholdLdData[sfb]) {
                  distNewSum = distOldSum << 1;
                  break;
                }
                distNewSum += CalcInvLdData(sfbDistNew[sfb]) >> DIST_FAC_SHIFT;
              }
            }
            if (distNewSum < distOldSum) {
              deltaPe = deltaPeNew;
              for (sfb = startSfb; sfb < stopSfb; sfb++) {
                if (scf[sfb] != FDK_INT_MIN) {
                  sfbWidth = psyOutChan->sfbOffsets[sfb + 1] -
                             psyOutChan->sfbOffsets[sfb];
                  sfbOffs = psyOutChan->sfbOffsets[sfb];
                  scf[sfb] = scfAct;
                  sfbDist[sfb] = sfbDistNew[sfb];

                  for (k = 0; k < sfbWidth; k++)
                    quantSpec[sfbOffs + k] = quantSpecTmp[sfbOffs + k];
                }
              }
            }
          }
        }

      } while (stopSfb <= sfbCnt);

    } while (scfAct > scfMin);
  }
}

static void FDKaacEnc_FDKaacEnc_assimilateMultipleScf2(
    PSY_OUT_CHANNEL *psyOutChan, QC_OUT_CHANNEL *qcOutChannel, SHORT *quantSpec,
    SHORT *quantSpecTmp, INT dZoneQuantEnable, INT *scf, const INT *minScf,
    FIXP_DBL *sfbDist, FIXP_DBL *sfbConstPePart, FIXP_DBL *sfbFormFactorLdData,
    FIXP_DBL *sfbNRelevantLines) {
  INT sfb, startSfb, stopSfb;
  INT scfTmp[MAX_GROUPED_SFB], scfAct, scfNew;
  INT scfPrev, scfNext, scfPrevNextMin, scfPrevNextMax, scfLo, scfHi;
  INT scfMin, scfMax;
  INT *sfbOffs = psyOutChan->sfbOffsets;
  FIXP_DBL sfbDistNew[MAX_GROUPED_SFB], sfbDistMax[MAX_GROUPED_SFB];
  FIXP_DBL distOldSum, distNewSum;
  INT deltaScfBits;
  FIXP_DBL deltaSpecPe;
  FIXP_DBL deltaPe = FL2FXCONST_DBL(0.0f);
  FIXP_DBL deltaPeNew = FL2FXCONST_DBL(0.0f);
  INT sfbCnt = psyOutChan->sfbCnt;
  INT bSuccess, bCheckScf;
  INT i, k;

  scfMin = FDK_INT_MAX;
  scfMax = FDK_INT_MIN;
  for (sfb = 0; sfb < sfbCnt; sfb++) {
    if (scf[sfb] != FDK_INT_MIN) {
      scfMin = fixMin(scfMin, scf[sfb]);
      scfMax = fixMax(scfMax, scf[sfb]);
    }
  }

  stopSfb = 0;
  scfAct = FDK_INT_MIN;
  do {
    scfPrev = scfAct;

    sfb = stopSfb;
    while (sfb < sfbCnt && (scf[sfb] == FDK_INT_MIN)) sfb++;
    startSfb = sfb;
    scfAct = scf[startSfb];
    sfb++;
    while (sfb < sfbCnt &&
           ((scf[sfb] == FDK_INT_MIN) || (scf[sfb] == scf[startSfb])))
      sfb++;
    stopSfb = sfb;

    if (stopSfb < sfbCnt)
      scfNext = scf[stopSfb];
    else
      scfNext = scfAct;

    if (scfPrev == FDK_INT_MIN) scfPrev = scfAct;

    scfPrevNextMax = fixMax(scfPrev, scfNext);
    scfPrevNextMin = fixMin(scfPrev, scfNext);

    scfHi = fixMax(scfPrevNextMax, scfAct);
    if (scfPrevNextMax >= scfAct)
      scfLo = fixMin(scfAct, scfPrevNextMin);
    else
      scfLo = scfPrevNextMax;

    if (startSfb < sfbCnt &&
        scfHi - scfLo <= MAX_SCF_DELTA) {
      if (scfHi > scf[startSfb]) {
        for (sfb = startSfb; sfb < stopSfb; sfb++) {
          if (scf[sfb] != FDK_INT_MIN) {
            sfbDistMax[sfb] = fMult(FL2FXCONST_DBL(1.0f / 3.0f),
                                    qcOutChannel->sfbThresholdLdData[sfb]) +
                              fMult(FL2FXCONST_DBL(1.0f / 3.0f), sfbDist[sfb]) +
                              fMult(FL2FXCONST_DBL(1.0f / 3.0f), sfbDist[sfb]);
            sfbDistMax[sfb] =
                fixMax(sfbDistMax[sfb], qcOutChannel->sfbEnergyLdData[sfb] -
                                            FL2FXCONST_DBL(0.15571537944));
            sfbDistMax[sfb] =
                fixMin(sfbDistMax[sfb], qcOutChannel->sfbThresholdLdData[sfb]);
          }
        }

        bCheckScf = 1;
        for (scfNew = scf[startSfb] + 1; scfNew <= scfHi; scfNew++) {
          for (k = 0; k < MAX_GROUPED_SFB; k++) scfTmp[k] = scf[k];

          for (sfb = startSfb; sfb < stopSfb; sfb++) {
            if (scfTmp[sfb] != FDK_INT_MIN) scfTmp[sfb] = scfNew;
          }

          deltaScfBits = FDKaacEnc_countScfBitsDiff(scf, scfTmp, sfbCnt,
                                                    startSfb, stopSfb);

          deltaSpecPe = FDKaacEnc_calcSpecPeDiff(
              psyOutChan, qcOutChannel, scf, scfTmp, sfbConstPePart,
              sfbFormFactorLdData, sfbNRelevantLines, startSfb, stopSfb);

          deltaPeNew = deltaPe + (FIXP_DBL)deltaScfBits + deltaSpecPe;

          if (deltaPeNew < FL2FXCONST_DBL(0.0f)) {
            bSuccess = 1;

            for (sfb = startSfb; sfb < stopSfb; sfb++) {
              if (scfTmp[sfb] != FDK_INT_MIN) {
                sfbDistNew[sfb] = FDKaacEnc_calcSfbDist(
                    qcOutChannel->mdctSpectrum + sfbOffs[sfb],
                    quantSpecTmp + sfbOffs[sfb],
                    sfbOffs[sfb + 1] - sfbOffs[sfb], scfNew, dZoneQuantEnable);

                if (sfbDistNew[sfb] > sfbDistMax[sfb]) {
                  bSuccess = 0;
                  if (sfbDistNew[sfb] == qcOutChannel->sfbEnergyLdData[sfb]) {
                    bCheckScf = 0;
                  }
                  break;
                }
              }
            }
            if (bCheckScf == 0)
              break;
            if (bSuccess) {
              deltaPe = deltaPeNew;
              for (sfb = startSfb; sfb < stopSfb; sfb++) {
                if (scf[sfb] != FDK_INT_MIN) {
                  scf[sfb] = scfNew;
                  sfbDist[sfb] = sfbDistNew[sfb];

                  for (k = 0; k < sfbOffs[sfb + 1] - sfbOffs[sfb]; k++)
                    quantSpec[sfbOffs[sfb] + k] =
                        quantSpecTmp[sfbOffs[sfb] + k];
                }
              }
            }
          }
        }
      }

      if (scfAct == scf[startSfb] && scfLo < scfAct &&
          scfMax - scfMin <= MAX_SCF_DELTA) {
        int bminScfViolation = 0;

        for (k = 0; k < MAX_GROUPED_SFB; k++) scfTmp[k] = scf[k];

        scfNew = scfLo;

        for (sfb = startSfb; sfb < stopSfb; sfb++) {
          if (scfTmp[sfb] != FDK_INT_MIN) {
            scfTmp[sfb] = scfNew;
            if (scfNew < minScf[sfb]) bminScfViolation = 1;
          }
        }

        if (!bminScfViolation) {
          deltaScfBits = FDKaacEnc_countScfBitsDiff(scf, scfTmp, sfbCnt,
                                                    startSfb, stopSfb);

          deltaSpecPe = FDKaacEnc_calcSpecPeDiff(
              psyOutChan, qcOutChannel, scf, scfTmp, sfbConstPePart,
              sfbFormFactorLdData, sfbNRelevantLines, startSfb, stopSfb);

          deltaPeNew = deltaPe + (FIXP_DBL)deltaScfBits + deltaSpecPe;
        }

        if (!bminScfViolation && deltaPeNew < FL2FXCONST_DBL(0.0f)) {
          distOldSum = distNewSum = FL2FXCONST_DBL(0.0f);
          for (sfb = startSfb; sfb < stopSfb; sfb++) {
            if (scfTmp[sfb] != FDK_INT_MIN) {
              distOldSum += CalcInvLdData(sfbDist[sfb]) >> DIST_FAC_SHIFT;

              sfbDistNew[sfb] = FDKaacEnc_calcSfbDist(
                  qcOutChannel->mdctSpectrum + sfbOffs[sfb],
                  quantSpecTmp + sfbOffs[sfb], sfbOffs[sfb + 1] - sfbOffs[sfb],
                  scfNew, dZoneQuantEnable);

              if (sfbDistNew[sfb] > qcOutChannel->sfbThresholdLdData[sfb]) {
                distNewSum = distOldSum << 1;
                break;
              }
              distNewSum += CalcInvLdData(sfbDistNew[sfb]) >> DIST_FAC_SHIFT;
            }
          }
          if (distNewSum < fMult(FL2FXCONST_DBL(0.8f), distOldSum)) {
            deltaPe = deltaPeNew;
            for (sfb = startSfb; sfb < stopSfb; sfb++) {
              if (scf[sfb] != FDK_INT_MIN) {
                scf[sfb] = scfNew;
                sfbDist[sfb] = sfbDistNew[sfb];

                for (k = 0; k < sfbOffs[sfb + 1] - sfbOffs[sfb]; k++)
                  quantSpec[sfbOffs[sfb] + k] = quantSpecTmp[sfbOffs[sfb] + k];
              }
            }
          }
        }
      }

      if (scfMax - scfMin <=
          MAX_SCF_DELTA - 3) {

        for (k = 0; k < sfbCnt; k++) scfTmp[k] = scf[k];

        for (i = 0; i < 3; i++) {
          scfNew = scfTmp[startSfb] - 1;
          for (sfb = startSfb; sfb < stopSfb; sfb++) {
            if (scfTmp[sfb] != FDK_INT_MIN) scfTmp[sfb] = scfNew;
          }
          deltaScfBits = FDKaacEnc_countScfBitsDiff(scf, scfTmp, sfbCnt,
                                                    startSfb, stopSfb);
          deltaPeNew = deltaPe + (FIXP_DBL)deltaScfBits;
          if (deltaPeNew <= FL2FXCONST_DBL(0.0f)) {
            bSuccess = 1;
            distOldSum = distNewSum = FL2FXCONST_DBL(0.0f);
            for (sfb = startSfb; sfb < stopSfb; sfb++) {
              if (scfTmp[sfb] != FDK_INT_MIN) {
                FIXP_DBL sfbEnQ;
                FDKaacEnc_calcSfbQuantEnergyAndDist(
                    qcOutChannel->mdctSpectrum + sfbOffs[sfb],
                    quantSpec + sfbOffs[sfb], sfbOffs[sfb + 1] - sfbOffs[sfb],
                    scfNew, &sfbEnQ, &sfbDistNew[sfb]);

                distOldSum += CalcInvLdData(sfbDist[sfb]) >> DIST_FAC_SHIFT;
                distNewSum += CalcInvLdData(sfbDistNew[sfb]) >> DIST_FAC_SHIFT;

                if ((sfbDistNew[sfb] >
                     (sfbDist[sfb] + FL2FXCONST_DBL(0.00259488556167f))) ||
                    (sfbEnQ < (qcOutChannel->sfbEnergyLdData[sfb] -
                               FL2FXCONST_DBL(0.00778722686652f)))) {
                  bSuccess = 0;
                  break;
                }
              }
            }
            if (distNewSum < distOldSum && bSuccess) {
              deltaPe = deltaPeNew;
              for (sfb = startSfb; sfb < stopSfb; sfb++) {
                if (scf[sfb] != FDK_INT_MIN) {
                  scf[sfb] = scfNew;
                  sfbDist[sfb] = sfbDistNew[sfb];
                }
              }
            }
          }
        }
      }
    }
  } while (stopSfb <= sfbCnt);
}

static void FDKaacEnc_EstimateScaleFactorsChannel(
    QC_OUT_CHANNEL *qcOutChannel, PSY_OUT_CHANNEL *psyOutChannel,
    INT *RESTRICT scf, INT *RESTRICT globalGain,
    FIXP_DBL *RESTRICT sfbFormFactorLdData, const INT invQuant,
    SHORT *RESTRICT quantSpec, const INT dZoneQuantEnable) {
  INT i, j, sfb, sfbOffs;
  INT scfInt;
  INT maxSf;
  INT minSf;
  FIXP_DBL threshLdData;
  FIXP_DBL energyLdData;
  FIXP_DBL energyPartLdData;
  FIXP_DBL thresholdPartLdData;
  FIXP_DBL scfFract;
  FIXP_DBL maxSpec;
  INT minScfCalculated[MAX_GROUPED_SFB];
  FIXP_DBL sfbDistLdData[MAX_GROUPED_SFB];
  C_ALLOC_SCRATCH_START(quantSpecTmp, SHORT, (1024))
  INT minSfMaxQuant[MAX_GROUPED_SFB];

  FIXP_DBL threshConstLdData =
      FL2FXCONST_DBL(0.04304511722f);
  FIXP_DBL convConst = FL2FXCONST_DBL(0.30102999566f);
  FIXP_DBL c1Const =
      FL2FXCONST_DBL(-0.27083183594f);

  if (invQuant > 0) {
    FDKmemclear(quantSpec, (1024) * sizeof(SHORT));
  }

  for (i = 0; i < psyOutChannel->sfbCnt; i++) {
    scf[i] = FDK_INT_MIN;
  }

  for (i = 0; i < MAX_GROUPED_SFB; i++) {
    minSfMaxQuant[i] = FDK_INT_MIN;
  }

  for (sfbOffs = 0; sfbOffs < psyOutChannel->sfbCnt;
       sfbOffs += psyOutChannel->sfbPerGroup) {
    for (sfb = 0; sfb < psyOutChannel->maxSfbPerGroup; sfb++) {
      threshLdData = qcOutChannel->sfbThresholdLdData[sfbOffs + sfb];
      energyLdData = qcOutChannel->sfbEnergyLdData[sfbOffs + sfb];

      sfbDistLdData[sfbOffs + sfb] = energyLdData;

      if (energyLdData > threshLdData) {
        FIXP_DBL tmp;

        energyPartLdData =
            sfbFormFactorLdData[sfbOffs + sfb] + FL2FXCONST_DBL(0.09375f);

        thresholdPartLdData = threshConstLdData + threshLdData;

        scfFract = thresholdPartLdData - energyPartLdData;
        scfFract = fMult(convConst, scfFract);
        scfFract = scfFract + fMult(FL2FXCONST_DBL(0.8585f), scfFract >> 3);

        scfInt =
            (INT)(scfFract >>
                  ((DFRACT_BITS - 1) - 3 -
                   LD_DATA_SHIFT));

        maxSpec = FL2FXCONST_DBL(0.0f);

        DWORD_ALIGNED(qcOutChannel->mdctSpectrum);
        for (j = psyOutChannel->sfbOffsets[sfbOffs + sfb];
             j < psyOutChannel->sfbOffsets[sfbOffs + sfb + 1]; j += 4) {
          maxSpec = fMax(maxSpec,
                         fMax(fMax(fAbs(qcOutChannel->mdctSpectrum[j + 0]),
                                   fAbs(qcOutChannel->mdctSpectrum[j + 1])),
                              fMax(fAbs(qcOutChannel->mdctSpectrum[j + 2]),
                                   fAbs(qcOutChannel->mdctSpectrum[j + 3]))));
        }
        tmp = CalcLdData(maxSpec);
        if (c1Const > FL2FXCONST_DBL(-1.f) - tmp) {
          minSfMaxQuant[sfbOffs + sfb] =
              ((INT)((c1Const + tmp) >> ((DFRACT_BITS - 1) - 8))) + 1;
        } else {
          minSfMaxQuant[sfbOffs + sfb] =
              ((INT)(FL2FXCONST_DBL(-1.f) >> ((DFRACT_BITS - 1) - 8))) + 1;
        }

        scfInt = fixMax(scfInt, minSfMaxQuant[sfbOffs + sfb]);

        if (invQuant > 0) {
          scfInt = FDKaacEnc_improveScf(
              qcOutChannel->mdctSpectrum +
                  psyOutChannel->sfbOffsets[sfbOffs + sfb],
              quantSpec + psyOutChannel->sfbOffsets[sfbOffs + sfb],
              quantSpecTmp + psyOutChannel->sfbOffsets[sfbOffs + sfb],
              psyOutChannel->sfbOffsets[sfbOffs + sfb + 1] -
                  psyOutChannel->sfbOffsets[sfbOffs + sfb],
              threshLdData, scfInt, minSfMaxQuant[sfbOffs + sfb],
              &sfbDistLdData[sfbOffs + sfb], &minScfCalculated[sfbOffs + sfb],
              dZoneQuantEnable);
        }
        scf[sfbOffs + sfb] = scfInt;
      }
    }
  }

  if (invQuant > 0) {
    FIXP_DBL sfbConstPePart[MAX_GROUPED_SFB];
    FIXP_DBL sfbNRelevantLines[MAX_GROUPED_SFB];

    for (i = 0; i < psyOutChannel->sfbCnt; i++)
      sfbConstPePart[i] = (FIXP_DBL)FDK_INT_MIN;

    FDKaacEnc_calcSfbRelevantLines(
        sfbFormFactorLdData, qcOutChannel->sfbEnergyLdData,
        qcOutChannel->sfbThresholdLdData, psyOutChannel->sfbOffsets,
        psyOutChannel->sfbCnt, psyOutChannel->sfbPerGroup,
        psyOutChannel->maxSfbPerGroup, sfbNRelevantLines);

    FDKaacEnc_assimilateSingleScf(
        psyOutChannel, qcOutChannel, quantSpec, quantSpecTmp, dZoneQuantEnable,
        scf, minSfMaxQuant, sfbDistLdData, sfbConstPePart, sfbFormFactorLdData,
        sfbNRelevantLines, minScfCalculated, 1);

    if (invQuant > 1) {
      FDKaacEnc_assimilateMultipleScf(
          psyOutChannel, qcOutChannel, quantSpec, quantSpecTmp,
          dZoneQuantEnable, scf, minSfMaxQuant, sfbDistLdData, sfbConstPePart,
          sfbFormFactorLdData, sfbNRelevantLines);

      FDKaacEnc_FDKaacEnc_assimilateMultipleScf2(
          psyOutChannel, qcOutChannel, quantSpec, quantSpecTmp,
          dZoneQuantEnable, scf, minSfMaxQuant, sfbDistLdData, sfbConstPePart,
          sfbFormFactorLdData, sfbNRelevantLines);
    }
  }

  minSf = FDK_INT_MAX;
  for (sfbOffs = 0; sfbOffs < psyOutChannel->sfbCnt;
       sfbOffs += psyOutChannel->sfbPerGroup) {
    for (sfb = 0; sfb < psyOutChannel->maxSfbPerGroup; sfb++) {
      if (scf[sfbOffs + sfb] != FDK_INT_MIN)
        minSf = fixMin(minSf, scf[sfbOffs + sfb]);
    }
  }

  for (sfbOffs = 0; sfbOffs < psyOutChannel->sfbCnt;
       sfbOffs += psyOutChannel->sfbPerGroup) {
    for (sfb = 0; sfb < psyOutChannel->maxSfbPerGroup; sfb++) {
      if ((scf[sfbOffs + sfb] != FDK_INT_MIN) &&
          (minSf + MAX_SCF_DELTA) < scf[sfbOffs + sfb]) {
        scf[sfbOffs + sfb] = minSf + MAX_SCF_DELTA;
        if (invQuant > 0) {
          sfbDistLdData[sfbOffs + sfb] = FDKaacEnc_calcSfbDist(
              qcOutChannel->mdctSpectrum +
                  psyOutChannel->sfbOffsets[sfbOffs + sfb],
              quantSpec + psyOutChannel->sfbOffsets[sfbOffs + sfb],
              psyOutChannel->sfbOffsets[sfbOffs + sfb + 1] -
                  psyOutChannel->sfbOffsets[sfbOffs + sfb],
              scf[sfbOffs + sfb], dZoneQuantEnable);
        }
      }
    }
  }

  maxSf = FDK_INT_MIN;
  for (sfbOffs = 0; sfbOffs < psyOutChannel->sfbCnt;
       sfbOffs += psyOutChannel->sfbPerGroup) {
    for (sfb = 0; sfb < psyOutChannel->maxSfbPerGroup; sfb++) {
      maxSf = fixMax(maxSf, scf[sfbOffs + sfb]);
    }
  }

  if (maxSf > FDK_INT_MIN) {
    *globalGain = maxSf;
    for (sfbOffs = 0; sfbOffs < psyOutChannel->sfbCnt;
         sfbOffs += psyOutChannel->sfbPerGroup) {
      for (sfb = 0; sfb < psyOutChannel->maxSfbPerGroup; sfb++) {
        if (scf[sfbOffs + sfb] == FDK_INT_MIN) {
          scf[sfbOffs + sfb] = 0;
          for (j = psyOutChannel->sfbOffsets[sfbOffs + sfb];
               j < psyOutChannel->sfbOffsets[sfbOffs + sfb + 1]; j++) {
            qcOutChannel->mdctSpectrum[j] = FL2FXCONST_DBL(0.0f);
          }
        } else {
          scf[sfbOffs + sfb] = maxSf - scf[sfbOffs + sfb];
        }
      }
    }
  } else {
    *globalGain = 0;
    for (sfbOffs = 0; sfbOffs < psyOutChannel->sfbCnt;
         sfbOffs += psyOutChannel->sfbPerGroup) {
      for (sfb = 0; sfb < psyOutChannel->maxSfbPerGroup; sfb++) {
        scf[sfbOffs + sfb] = 0;
        for (j = psyOutChannel->sfbOffsets[sfbOffs + sfb];
             j < psyOutChannel->sfbOffsets[sfbOffs + sfb + 1]; j++) {
          qcOutChannel->mdctSpectrum[j] = FL2FXCONST_DBL(0.0f);
        }
      }
    }
  }

  C_ALLOC_SCRATCH_END(quantSpecTmp, SHORT, (1024))
}

void FDKaacEnc_EstimateScaleFactors(PSY_OUT_CHANNEL *psyOutChannel[],
                                    QC_OUT_CHANNEL *qcOutChannel[],
                                    const INT invQuant,
                                    const INT dZoneQuantEnable,
                                    const INT nChannels) {
  int ch;

  for (ch = 0; ch < nChannels; ch++) {
    FDKaacEnc_EstimateScaleFactorsChannel(
        qcOutChannel[ch], psyOutChannel[ch], qcOutChannel[ch]->scf,
        &qcOutChannel[ch]->globalGain, qcOutChannel[ch]->sfbFormFactorLdData,
        invQuant, qcOutChannel[ch]->quantSpec, dZoneQuantEnable);
  }
}
