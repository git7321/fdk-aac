/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "mh_det.h"
#include "sbrenc_ram.h"
#include "sbr_misc.h"
#include "genericStds.h"

#define SFM_SHIFT 2
#define SFM_SCALE (MAXVAL_DBL >> SFM_SHIFT)

static const DETECTOR_PARAMETERS_MH paramsAac = {
    9,
    {
        FL2FXCONST_DBL(20.0f * RELAXATION_FLOAT),
        FL2FXCONST_DBL(1.26f * RELAXATION_FLOAT),
        FL2FXCONST_DBL(15.0f * RELAXATION_FLOAT),
        FL2FXCONST_DBL((1.0f / 15.0f) *
                       RELAXATION_FLOAT),
        FL2FXCONST_DBL(1.26f * RELAXATION_FLOAT),
        FL2FXCONST_DBL(0.3f) >> SFM_SHIFT,
        FL2FXCONST_DBL(0.1f) >> SFM_SHIFT,
        FL2FXCONST_DBL(0.3f),
        FL2FXCONST_DBL(0.5f),
        FL2FXCONST_DBL(-0.000112993269),
        FL2FXCONST_DBL(-0.000112993269),
        FL2FXCONST_DBL(
            -0.005030126483f)
    },
    50
};

static const DETECTOR_PARAMETERS_MH paramsAacLd = {
    16,
    {
        FL2FXCONST_DBL(25.0f * RELAXATION_FLOAT),
        FL2FXCONST_DBL(1.26f * RELAXATION_FLOAT),
        FL2FXCONST_DBL(15.0f * RELAXATION_FLOAT),
        FL2FXCONST_DBL((1.0f / 15.0f) *
                       RELAXATION_FLOAT),
        FL2FXCONST_DBL(1.26f * RELAXATION_FLOAT),
        FL2FXCONST_DBL(0.3f) >> SFM_SHIFT,
        FL2FXCONST_DBL(0.1f) >> SFM_SHIFT,
        FL2FXCONST_DBL(0.3f),
        FL2FXCONST_DBL(0.2f),
        FL2FXCONST_DBL(-0.000112993269),
        FL2FXCONST_DBL(-0.000112993269),
        FL2FXCONST_DBL(
            -0.005030126483f)
    },
    50
};

static void diff(FIXP_DBL *RESTRICT pTonalityOrig, FIXP_DBL *pDiffMapped2Scfb,
                 const UCHAR *RESTRICT pFreqBandTable, INT nScfb,
                 SCHAR *indexVector) {
  UCHAR i, ll, lu, k;
  FIXP_DBL maxValOrig, maxValSbr, tmp;
  INT scale;

  for (i = 0; i < nScfb; i++) {
    ll = pFreqBandTable[i];
    lu = pFreqBandTable[i + 1];

    maxValOrig = FL2FXCONST_DBL(0.0f);
    maxValSbr = FL2FXCONST_DBL(0.0f);

    for (k = ll; k < lu; k++) {
      maxValOrig = fixMax(maxValOrig, pTonalityOrig[k]);
      maxValSbr = fixMax(maxValSbr, pTonalityOrig[indexVector[k]]);
    }

    if ((maxValSbr >= RELAXATION)) {
      tmp = fDivNorm(maxValOrig, maxValSbr, &scale);
      pDiffMapped2Scfb[i] =
          scaleValue(fMult(tmp, RELAXATION_FRACT),
                     fixMax(-(DFRACT_BITS - 1), (scale - RELAXATION_SHIFT)));
    } else {
      pDiffMapped2Scfb[i] = maxValOrig;
    }
  }
}

static void calculateFlatnessMeasure(FIXP_DBL *pQuotaBuffer, SCHAR *indexVector,
                                     FIXP_DBL *pSfmOrigVec,
                                     FIXP_DBL *pSfmSbrVec,
                                     const UCHAR *pFreqBandTable, INT nSfb) {
  INT i, j;
  FIXP_DBL invBands, tmp1, tmp2;
  INT shiftFac0, shiftFacSum0;
  INT shiftFac1, shiftFacSum1;
  FIXP_DBL accu;

  for (i = 0; i < nSfb; i++) {
    INT ll = pFreqBandTable[i];
    INT lu = pFreqBandTable[i + 1];
    pSfmOrigVec[i] = (FIXP_DBL)(MAXVAL_DBL >> 2);
    pSfmSbrVec[i] = (FIXP_DBL)(MAXVAL_DBL >> 2);

    if (lu - ll > 1) {
      FIXP_DBL amOrig, amTransp, gmOrig, gmTransp, sfmOrig, sfmTransp;
      invBands = GetInvInt(lu - ll);
      shiftFacSum0 = 0;
      shiftFacSum1 = 0;
      amOrig = amTransp = FL2FXCONST_DBL(0.0f);
      gmOrig = gmTransp = (FIXP_DBL)MAXVAL_DBL;

      for (j = ll; j < lu; j++) {
        sfmOrig = pQuotaBuffer[j];
        sfmTransp = pQuotaBuffer[indexVector[j]];

        amOrig += fMult(sfmOrig, invBands);
        amTransp += fMult(sfmTransp, invBands);

        shiftFac0 = CountLeadingBits(sfmOrig);
        shiftFac1 = CountLeadingBits(sfmTransp);

        gmOrig = fMult(gmOrig, sfmOrig << shiftFac0);
        gmTransp = fMult(gmTransp, sfmTransp << shiftFac1);

        shiftFacSum0 += shiftFac0;
        shiftFacSum1 += shiftFac1;
      }

      if (gmOrig > FL2FXCONST_DBL(0.0f)) {
        tmp1 = CalcLdData(gmOrig);
        tmp1 = fMult(invBands, tmp1);

        accu = (FIXP_DBL)-shiftFacSum0 << (DFRACT_BITS - 1 - 8);
        tmp2 = fMultDiv2(invBands, accu) << (2 + 1);

        tmp2 = tmp1 + tmp2;
        gmOrig = CalcInvLdData(tmp2);
      } else {
        gmOrig = FL2FXCONST_DBL(0.0f);
      }

      if (gmTransp > FL2FXCONST_DBL(0.0f)) {
        tmp1 = CalcLdData(gmTransp);
        tmp1 = fMult(invBands, tmp1);

        accu = (FIXP_DBL)-shiftFacSum1 << (DFRACT_BITS - 1 - 8);
        tmp2 = fMultDiv2(invBands, accu) << (2 + 1);

        tmp2 = tmp1 + tmp2;
        gmTransp = CalcInvLdData(tmp2);
      } else {
        gmTransp = FL2FXCONST_DBL(0.0f);
      }
      if (amOrig != FL2FXCONST_DBL(0.0f))
        pSfmOrigVec[i] =
            FDKsbrEnc_LSI_divide_scale_fract(gmOrig, amOrig, SFM_SCALE);

      if (amTransp != FL2FXCONST_DBL(0.0f))
        pSfmSbrVec[i] =
            FDKsbrEnc_LSI_divide_scale_fract(gmTransp, amTransp, SFM_SCALE);
    }
  }
}

static void calculateDetectorInput(
    FIXP_DBL **RESTRICT pQuotaBuffer,
    SCHAR *RESTRICT indexVector, FIXP_DBL **RESTRICT tonalityDiff,
    FIXP_DBL **RESTRICT pSfmOrig, FIXP_DBL **RESTRICT pSfmSbr,
    const UCHAR *freqBandTable, INT nSfb, INT noEstPerFrame, INT move) {
  INT est;

  for (est = 0; est < noEstPerFrame; est++) {
    diff(pQuotaBuffer[est + move], tonalityDiff[est + move], freqBandTable,
         nSfb, indexVector);

    calculateFlatnessMeasure(pQuotaBuffer[est + move], indexVector,
                             pSfmOrig[est + move], pSfmSbr[est + move],
                             freqBandTable, nSfb);
  }
}

static void removeLowPassDetection(UCHAR *RESTRICT pAddHarmSfb,
                                   UCHAR **RESTRICT pDetectionVectors,
                                   INT start, INT stop, INT nSfb,
                                   const UCHAR *RESTRICT pFreqBandTable,
                                   FIXP_DBL *RESTRICT pNrgVector,
                                   THRES_HOLDS mhThresh)

{
  INT i, est;
  INT maxDerivPos = pFreqBandTable[nSfb];
  INT numBands = pFreqBandTable[nSfb];
  FIXP_DBL nrgLow, nrgHigh;
  FIXP_DBL nrgLD64, nrgLowLD64, nrgHighLD64, nrgDiffLD64;
  FIXP_DBL valLD64, maxValLD64, maxValAboveLD64;
  INT bLPsignal = 0;

  maxValLD64 = FL2FXCONST_DBL(-1.0f);
  for (i = numBands - 1 - 2; i > pFreqBandTable[0]; i--) {
    nrgLow = pNrgVector[i];
    nrgHigh = pNrgVector[i + 2];

    if (nrgLow != FL2FXCONST_DBL(0.0f) && nrgLow > nrgHigh) {
      nrgLowLD64 = CalcLdData(nrgLow >> 1);
      nrgDiffLD64 = CalcLdData((nrgLow >> 1) - (nrgHigh >> 1));
      valLD64 = nrgDiffLD64 - nrgLowLD64;
      if (valLD64 > maxValLD64) {
        maxDerivPos = i;
        maxValLD64 = valLD64;
      }
      if (maxValLD64 > mhThresh.derivThresMaxLD64) {
        break;
      }
    }
  }

  maxValAboveLD64 = FL2FXCONST_DBL(-1.0f);
  for (i = numBands - 1 - 2; i > maxDerivPos + 2; i--) {
    nrgLow = pNrgVector[i];
    nrgHigh = pNrgVector[i + 2];

    if (nrgLow != FL2FXCONST_DBL(0.0f) && nrgLow > nrgHigh) {
      nrgLowLD64 = CalcLdData(nrgLow >> 1);
      nrgDiffLD64 = CalcLdData((nrgLow >> 1) - (nrgHigh >> 1));
      valLD64 = nrgDiffLD64 - nrgLowLD64;
      if (valLD64 > maxValAboveLD64) {
        maxValAboveLD64 = valLD64;
      }
    } else {
      if (nrgHigh != FL2FXCONST_DBL(0.0f) && nrgHigh > nrgLow) {
        nrgHighLD64 = CalcLdData(nrgHigh >> 1);
        nrgDiffLD64 = CalcLdData((nrgHigh >> 1) - (nrgLow >> 1));
        valLD64 = nrgDiffLD64 - nrgHighLD64;
        if (valLD64 > maxValAboveLD64) {
          maxValAboveLD64 = valLD64;
        }
      }
    }
  }

  if (maxValLD64 > mhThresh.derivThresMaxLD64 &&
      maxValAboveLD64 < mhThresh.derivThresAboveLD64) {
    bLPsignal = 1;

    for (i = maxDerivPos - 1; i > maxDerivPos - 5 && i >= 0; i--) {
      if (pNrgVector[i] != FL2FXCONST_DBL(0.0f) &&
          pNrgVector[i] > pNrgVector[maxDerivPos + 2]) {
        nrgDiffLD64 = CalcLdData((pNrgVector[i] >> 1) -
                                 (pNrgVector[maxDerivPos + 2] >> 1));
        nrgLD64 = CalcLdData(pNrgVector[i] >> 1);
        valLD64 = nrgDiffLD64 - nrgLD64;
        if (valLD64 < mhThresh.derivThresBelowLD64) {
          bLPsignal = 0;
          break;
        }
      } else {
        bLPsignal = 0;
        break;
      }
    }
  }

  if (bLPsignal) {
    for (i = 0; i < nSfb; i++) {
      if (maxDerivPos >= pFreqBandTable[i] &&
          maxDerivPos < pFreqBandTable[i + 1])
        break;
    }

    if (pAddHarmSfb[i]) {
      pAddHarmSfb[i] = 0;
      for (est = start; est < stop; est++) {
        pDetectionVectors[est][i] = 0;
      }
    }
  }
}

static INT isDetectionOfNewToneAllowed(
    const SBR_FRAME_INFO *pFrameInfo, INT *pDetectionStartPos,
    INT noEstPerFrame, INT prevTransientFrame, INT prevTransientPos,
    INT prevTransientFlag, INT transientPosOffset, INT transientFlag,
    INT transientPos, INT deltaTime,
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR h_sbrMissingHarmonicsDetector) {
  INT transientFrame, newDetectionAllowed;

  transientFrame = 0;
  if (transientFlag) {
    if (transientPos + transientPosOffset <
        pFrameInfo->borders[pFrameInfo->nEnvelopes]) {
      transientFrame = 1;
      if (noEstPerFrame > 1) {
        if (transientPos + transientPosOffset >
            h_sbrMissingHarmonicsDetector->timeSlots >> 1) {
          *pDetectionStartPos = noEstPerFrame;
        } else {
          *pDetectionStartPos = noEstPerFrame >> 1;
        }

      } else {
        *pDetectionStartPos = noEstPerFrame;
      }
    }
  } else {
    if (prevTransientFlag && !prevTransientFrame) {
      transientFrame = 1;
      *pDetectionStartPos = 0;
    }
  }

  newDetectionAllowed = 0;
  if (transientFrame) {
    newDetectionAllowed = 1;
  } else {
    if (prevTransientFrame &&
        fixp_abs(pFrameInfo->borders[0] -
                 (prevTransientPos + transientPosOffset -
                  h_sbrMissingHarmonicsDetector->timeSlots)) < deltaTime) {
      newDetectionAllowed = 1;
      *pDetectionStartPos = 0;
    }
  }

  h_sbrMissingHarmonicsDetector->previousTransientFlag = transientFlag;
  h_sbrMissingHarmonicsDetector->previousTransientFrame = transientFrame;
  h_sbrMissingHarmonicsDetector->previousTransientPos = transientPos;

  return (newDetectionAllowed);
}

static void transientCleanUp(FIXP_DBL **quotaBuffer, INT nSfb,
                             UCHAR **detectionVectors, UCHAR *pAddHarmSfb,
                             UCHAR *pPrevAddHarmSfb, INT **signBuffer,
                             const UCHAR *pFreqBandTable, INT start, INT stop,
                             INT newDetectionAllowed, FIXP_DBL *pNrgVector,
                             THRES_HOLDS mhThresh) {
  INT i, j, est;

  for (est = start; est < stop; est++) {
    for (i = 0; i < nSfb; i++) {
      pAddHarmSfb[i] = pAddHarmSfb[i] || detectionVectors[est][i];
    }
  }

  if (newDetectionAllowed == 1) {
    for (i = 0; i < nSfb - 1; i++) {
      if (pAddHarmSfb[i] && pAddHarmSfb[i + 1]) {
        FIXP_DBL maxVal1, maxVal2;
        INT maxPos1, maxPos2, maxPosTime1, maxPosTime2;

        INT li = pFreqBandTable[i];
        INT ui = pFreqBandTable[i + 1];

        maxPosTime1 = start;
        maxPos1 = li;
        maxVal1 = quotaBuffer[start][li];
        for (est = start; est < stop; est++) {
          for (j = li; j < ui; j++) {
            if (quotaBuffer[est][j] > maxVal1) {
              maxVal1 = quotaBuffer[est][j];
              maxPos1 = j;
              maxPosTime1 = est;
            }
          }
        }

        li = pFreqBandTable[i + 1];
        ui = pFreqBandTable[i + 2];

        maxPosTime2 = start;
        maxPos2 = li;
        maxVal2 = quotaBuffer[start][li];
        for (est = start; est < stop; est++) {
          for (j = li; j < ui; j++) {
            if (quotaBuffer[est][j] > maxVal2) {
              maxVal2 = quotaBuffer[est][j];
              maxPos2 = j;
              maxPosTime2 = est;
            }
          }
        }

        if (maxPos2 - maxPos1 < 2) {
          if (pPrevAddHarmSfb[i] == 1 && pPrevAddHarmSfb[i + 1] == 0) {
            pAddHarmSfb[i + 1] = 0;
            for (est = start; est < stop; est++) {
              detectionVectors[est][i + 1] = 0;
            }
          } else {
            if (pPrevAddHarmSfb[i] == 0 && pPrevAddHarmSfb[i + 1] == 1) {
              pAddHarmSfb[i] = 0;
              for (est = start; est < stop; est++) {
                detectionVectors[est][i] = 0;
              }
            } else {
              if (maxVal1 > maxVal2) {
                if (signBuffer[maxPosTime1][maxPos2] < 0 &&
                    signBuffer[maxPosTime1][maxPos1] > 0) {
                  pAddHarmSfb[i + 1] = 0;
                  for (est = start; est < stop; est++) {
                    detectionVectors[est][i + 1] = 0;
                  }
                }
              } else {
                if (signBuffer[maxPosTime2][maxPos2] < 0 &&
                    signBuffer[maxPosTime2][maxPos1] > 0) {
                  pAddHarmSfb[i] = 0;
                  for (est = start; est < stop; est++) {
                    detectionVectors[est][i] = 0;
                  }
                }
              }
            }
          }
        }
      }
    }

    removeLowPassDetection(pAddHarmSfb, detectionVectors, start, stop, nSfb,
                           pFreqBandTable, pNrgVector, mhThresh);
  } else {
    for (i = 0; i < nSfb; i++) {
      if (pAddHarmSfb[i] - pPrevAddHarmSfb[i] > 0) pAddHarmSfb[i] = 0;
    }
  }
}

static void detection(FIXP_DBL *quotaBuffer, FIXP_DBL *pDiffVecScfb, INT nSfb,
                      UCHAR *pHarmVec, const UCHAR *pFreqBandTable,
                      FIXP_DBL *sfmOrig, FIXP_DBL *sfmSbr,
                      GUIDE_VECTORS guideVectors, GUIDE_VECTORS newGuideVectors,
                      THRES_HOLDS mhThresh) {
  INT i, j, ll, lu;
  FIXP_DBL thresTemp, thresOrig;

  for (i = 0; i < nSfb; i++) {
    thresTemp = (guideVectors.guideVectorDiff[i] != FL2FXCONST_DBL(0.0f))
                    ? fMax(fMult(mhThresh.decayGuideDiff,
                                 guideVectors.guideVectorDiff[i]),
                           mhThresh.thresHoldDiffGuide)
                    : mhThresh.thresHoldDiff;

    thresTemp = fMin(thresTemp, mhThresh.thresHoldDiff);

    if (pDiffVecScfb[i] > thresTemp) {
      pHarmVec[i] = 1;
      newGuideVectors.guideVectorDiff[i] = pDiffVecScfb[i];
    } else {
      if (guideVectors.guideVectorDiff[i] != FL2FXCONST_DBL(0.0f)) {
        guideVectors.guideVectorOrig[i] = mhThresh.thresHoldToneGuide;
      }
    }
  }

  for (i = 0; i < nSfb; i++) {
    ll = pFreqBandTable[i];
    lu = pFreqBandTable[i + 1];

    thresOrig =
        fixMax(fMult(guideVectors.guideVectorOrig[i], mhThresh.decayGuideOrig),
               mhThresh.thresHoldToneGuide);
    thresOrig = fixMin(thresOrig, mhThresh.thresHoldTone);

    if (guideVectors.guideVectorOrig[i] != FL2FXCONST_DBL(0.0f)) {
      for (j = ll; j < lu; j++) {
        if (quotaBuffer[j] > thresOrig) {
          pHarmVec[i] = 1;
          newGuideVectors.guideVectorOrig[i] = quotaBuffer[j];
        }
      }
    }
  }

  thresOrig = mhThresh.thresHoldTone;

  for (i = 0; i < nSfb; i++) {
    ll = pFreqBandTable[i];
    lu = pFreqBandTable[i + 1];

    if (pHarmVec[i] == 0) {
      if (lu - ll > 1) {
        for (j = ll; j < lu; j++) {
          if (quotaBuffer[j] > thresOrig &&
              (sfmSbr[i] > mhThresh.sfmThresSbr &&
               sfmOrig[i] < mhThresh.sfmThresOrig)) {
            pHarmVec[i] = 1;
            newGuideVectors.guideVectorOrig[i] = quotaBuffer[j];
          }
        }
      } else {
        if (i < nSfb - 1) {
          ll = pFreqBandTable[i];

          if (i > 0) {
            if (quotaBuffer[ll] > mhThresh.thresHoldTone &&
                (pDiffVecScfb[i + 1] < mhThresh.invThresHoldTone ||
                 pDiffVecScfb[i - 1] < mhThresh.invThresHoldTone)) {
              pHarmVec[i] = 1;
              newGuideVectors.guideVectorOrig[i] = quotaBuffer[ll];
            }
          } else {
            if (quotaBuffer[ll] > mhThresh.thresHoldTone &&
                pDiffVecScfb[i + 1] < mhThresh.invThresHoldTone) {
              pHarmVec[i] = 1;
              newGuideVectors.guideVectorOrig[i] = quotaBuffer[ll];
            }
          }
        }
      }
    }
  }
}

static void detectionWithPrediction(
    FIXP_DBL **quotaBuffer, FIXP_DBL **pDiffVecScfb, INT **signBuffer, INT nSfb,
    const UCHAR *pFreqBandTable, FIXP_DBL **sfmOrig, FIXP_DBL **sfmSbr,
    UCHAR **detectionVectors, UCHAR *pPrevAddHarmSfb,
    GUIDE_VECTORS *guideVectors, INT noEstPerFrame, INT detectionStart,
    INT totNoEst, INT newDetectionAllowed, INT *pAddHarmFlag,
    UCHAR *pAddHarmSfb, FIXP_DBL *pNrgVector,
    const DETECTOR_PARAMETERS_MH *mhParams) {
  INT est = 0, i;
  INT start;

  FDKmemclear(pAddHarmSfb, nSfb * sizeof(UCHAR));

  if (newDetectionAllowed) {
    if (totNoEst > 1) {
      start = detectionStart + 1;

      if (start != 0) {
        FDKmemcpy(guideVectors[start].guideVectorDiff,
                  guideVectors[0].guideVectorDiff, nSfb * sizeof(FIXP_DBL));
        FDKmemcpy(guideVectors[start].guideVectorOrig,
                  guideVectors[0].guideVectorOrig, nSfb * sizeof(FIXP_DBL));
        FDKmemclear(guideVectors[start - 1].guideVectorDetected,
                    nSfb * sizeof(UCHAR));
      }
    } else {
      start = 0;
    }
  } else {
    start = 0;
  }

  for (est = start; est < totNoEst; est++) {
    if (est > 0) {
      FDKmemcpy(guideVectors[est].guideVectorDetected,
                detectionVectors[est - 1], nSfb * sizeof(UCHAR));
    }

    FDKmemclear(detectionVectors[est], nSfb * sizeof(UCHAR));

    if (est < totNoEst - 1) {
      FDKmemclear(guideVectors[est + 1].guideVectorDiff,
                  nSfb * sizeof(FIXP_DBL));
      FDKmemclear(guideVectors[est + 1].guideVectorOrig,
                  nSfb * sizeof(FIXP_DBL));
      FDKmemclear(guideVectors[est + 1].guideVectorDetected,
                  nSfb * sizeof(UCHAR));

      detection(quotaBuffer[est], pDiffVecScfb[est], nSfb,
                detectionVectors[est], pFreqBandTable, sfmOrig[est],
                sfmSbr[est], guideVectors[est], guideVectors[est + 1],
                mhParams->thresHolds);
    } else {
      FDKmemclear(guideVectors[est].guideVectorDiff, nSfb * sizeof(FIXP_DBL));
      FDKmemclear(guideVectors[est].guideVectorOrig, nSfb * sizeof(FIXP_DBL));
      FDKmemclear(guideVectors[est].guideVectorDetected, nSfb * sizeof(UCHAR));

      detection(quotaBuffer[est], pDiffVecScfb[est], nSfb,
                detectionVectors[est], pFreqBandTable, sfmOrig[est],
                sfmSbr[est], guideVectors[est], guideVectors[est],
                mhParams->thresHolds);
    }
  }

  transientCleanUp(quotaBuffer, nSfb, detectionVectors, pAddHarmSfb,
                   pPrevAddHarmSfb, signBuffer, pFreqBandTable, start, totNoEst,
                   newDetectionAllowed, pNrgVector, mhParams->thresHolds);

  *pAddHarmFlag = 0;
  for (i = 0; i < nSfb; i++) {
    if (pAddHarmSfb[i]) {
      *pAddHarmFlag = 1;
      break;
    }
  }

  FDKmemcpy(pPrevAddHarmSfb, pAddHarmSfb, nSfb * sizeof(UCHAR));
  FDKmemcpy(guideVectors[0].guideVectorDetected, pAddHarmSfb,
            nSfb * sizeof(INT));

  for (i = 0; i < nSfb; i++) {
    guideVectors[0].guideVectorDiff[i] = FL2FXCONST_DBL(0.0f);
    guideVectors[0].guideVectorOrig[i] = FL2FXCONST_DBL(0.0f);

    if (pAddHarmSfb[i] == 1) {
      for (est = start; est < totNoEst; est++) {
        if (guideVectors[est].guideVectorDiff[i] != FL2FXCONST_DBL(0.0f)) {
          guideVectors[0].guideVectorDiff[i] =
              guideVectors[est].guideVectorDiff[i];
        }
        if (guideVectors[est].guideVectorOrig[i] != FL2FXCONST_DBL(0.0f)) {
          guideVectors[0].guideVectorOrig[i] =
              guideVectors[est].guideVectorOrig[i];
        }
      }
    }
  }
}

static void calculateCompVector(UCHAR *pAddHarmSfb, FIXP_DBL **pTonalityMatrix,
                                INT **pSignMatrix, UCHAR *pEnvComp, INT nSfb,
                                const UCHAR *freqBandTable, INT totNoEst,
                                INT maxComp, UCHAR *pPrevEnvComp,
                                INT newDetectionAllowed) {
  INT scfBand, est, l, ll, lu, maxPosF, maxPosT;
  FIXP_DBL maxVal;
  INT compValue;
  FIXP_DBL tmp;

  FDKmemclear(pEnvComp, nSfb * sizeof(UCHAR));

  for (scfBand = 0; scfBand < nSfb; scfBand++) {
    if (pAddHarmSfb[scfBand]) {
      ll = freqBandTable[scfBand];
      lu = freqBandTable[scfBand + 1];

      maxPosF = 0;
      maxPosT = 0;
      maxVal = FL2FXCONST_DBL(0.0f);

      for (est = 0; est < totNoEst; est++) {
        for (l = ll; l < lu; l++) {
          if (pTonalityMatrix[est][l] > maxVal) {
            maxVal = pTonalityMatrix[est][l];
            maxPosF = l;
            maxPosT = est;
          }
        }
      }

      if (maxPosF == ll && scfBand) {
        if (!pAddHarmSfb[scfBand - 1]) {
          if (pSignMatrix[maxPosT][maxPosF - 1] > 0 &&
              pSignMatrix[maxPosT][maxPosF] < 0) {
            tmp = fixp_abs(
                (FIXP_DBL)CalcLdData(pTonalityMatrix[maxPosT][maxPosF - 1]) +
                RELAXATION_LD64);
            tmp = (tmp >> (DFRACT_BITS - 1 - LD_DATA_SHIFT - 1)) +
                  (FIXP_DBL)1;
            compValue = ((INT)(LONG)tmp) >> 1;

            if (compValue > maxComp) compValue = maxComp;

            pEnvComp[scfBand - 1] = compValue;
          }
        }
      }

      if (maxPosF == lu - 1 && scfBand + 1 < nSfb) {
        if (!pAddHarmSfb[scfBand + 1]) {
          if (pSignMatrix[maxPosT][maxPosF] > 0 &&
              pSignMatrix[maxPosT][maxPosF + 1] < 0) {
            tmp = fixp_abs(
                (FIXP_DBL)CalcLdData(pTonalityMatrix[maxPosT][maxPosF + 1]) +
                RELAXATION_LD64);
            tmp = (tmp >> (DFRACT_BITS - 1 - LD_DATA_SHIFT - 1)) +
                  (FIXP_DBL)1;
            compValue = ((INT)(LONG)tmp) >> 1;

            if (compValue > maxComp) compValue = maxComp;

            pEnvComp[scfBand + 1] = compValue;
          }
        }
      }
    }
  }

  if (newDetectionAllowed == 0) {
    for (scfBand = 0; scfBand < nSfb; scfBand++) {
      if (pEnvComp[scfBand] != 0 && pPrevEnvComp[scfBand] == 0)
        pEnvComp[scfBand] = 0;
    }
  }

  FDKmemcpy(pPrevEnvComp, pEnvComp, nSfb * sizeof(UCHAR));
}

void FDKsbrEnc_SbrMissingHarmonicsDetectorQmf(
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR h_sbrMHDet, FIXP_DBL **pQuotaBuffer,
    INT **pSignBuffer, SCHAR *indexVector, const SBR_FRAME_INFO *pFrameInfo,
    const UCHAR *pTranInfo, INT *pAddHarmonicsFlag,
    UCHAR *pAddHarmonicsScaleFactorBands, const UCHAR *freqBandTable, INT nSfb,
    UCHAR *envelopeCompensation, FIXP_DBL *pNrgVector) {
  INT transientFlag = pTranInfo[1];
  INT transientPos = pTranInfo[0];
  INT newDetectionAllowed;
  INT transientDetStart = 0;

  UCHAR **detectionVectors = h_sbrMHDet->detectionVectors;
  INT move = h_sbrMHDet->move;
  INT noEstPerFrame = h_sbrMHDet->noEstPerFrame;
  INT totNoEst = h_sbrMHDet->totNoEst;
  INT prevTransientFlag = h_sbrMHDet->previousTransientFlag;
  INT prevTransientFrame = h_sbrMHDet->previousTransientFrame;
  INT transientPosOffset = h_sbrMHDet->transientPosOffset;
  INT prevTransientPos = h_sbrMHDet->previousTransientPos;
  GUIDE_VECTORS *guideVectors = h_sbrMHDet->guideVectors;
  INT deltaTime = h_sbrMHDet->mhParams->deltaTime;
  INT maxComp = h_sbrMHDet->mhParams->maxComp;

  int est;

  FDK_ASSERT(move <= (MAX_NO_OF_ESTIMATES >> 1));
  FDK_ASSERT(noEstPerFrame <= (MAX_NO_OF_ESTIMATES >> 1));

  FIXP_DBL *sfmSbr[MAX_NO_OF_ESTIMATES];
  FIXP_DBL *sfmOrig[MAX_NO_OF_ESTIMATES];
  FIXP_DBL *tonalityDiff[MAX_NO_OF_ESTIMATES];

  for (est = 0; est < MAX_NO_OF_ESTIMATES / 2; est++) {
    sfmSbr[est] = h_sbrMHDet->sfmSbr[est];
    sfmOrig[est] = h_sbrMHDet->sfmOrig[est];
    tonalityDiff[est] = h_sbrMHDet->tonalityDiff[est];
  }

  C_ALLOC_SCRATCH_START(_scratch, FIXP_DBL,
                        3 * MAX_NO_OF_ESTIMATES / 2 * MAX_FREQ_COEFFS)
  FIXP_DBL *scratch = _scratch;
  for (; est < MAX_NO_OF_ESTIMATES; est++) {
    sfmSbr[est] = scratch;
    scratch += MAX_FREQ_COEFFS;
    sfmOrig[est] = scratch;
    scratch += MAX_FREQ_COEFFS;
    tonalityDiff[est] = scratch;
    scratch += MAX_FREQ_COEFFS;
  }

  newDetectionAllowed = isDetectionOfNewToneAllowed(
      pFrameInfo, &transientDetStart, noEstPerFrame, prevTransientFrame,
      prevTransientPos, prevTransientFlag, transientPosOffset, transientFlag,
      transientPos, deltaTime, h_sbrMHDet);

  calculateDetectorInput(pQuotaBuffer, indexVector, tonalityDiff, sfmOrig,
                         sfmSbr, freqBandTable, nSfb, noEstPerFrame, move);

  detectionWithPrediction(pQuotaBuffer, tonalityDiff, pSignBuffer, nSfb,
                          freqBandTable, sfmOrig, sfmSbr, detectionVectors,
                          h_sbrMHDet->guideScfb, guideVectors, noEstPerFrame,
                          transientDetStart, totNoEst, newDetectionAllowed,
                          pAddHarmonicsFlag, pAddHarmonicsScaleFactorBands,
                          pNrgVector, h_sbrMHDet->mhParams);

  calculateCompVector(pAddHarmonicsScaleFactorBands, pQuotaBuffer, pSignBuffer,
                      envelopeCompensation, nSfb, freqBandTable, totNoEst,
                      maxComp, h_sbrMHDet->prevEnvelopeCompensation,
                      newDetectionAllowed);

  for (est = 0; est < move; est++) {
    FDKmemcpy(tonalityDiff[est], tonalityDiff[est + noEstPerFrame],
              sizeof(FIXP_DBL) * MAX_FREQ_COEFFS);
    FDKmemcpy(sfmOrig[est], sfmOrig[est + noEstPerFrame],
              sizeof(FIXP_DBL) * MAX_FREQ_COEFFS);
    FDKmemcpy(sfmSbr[est], sfmSbr[est + noEstPerFrame],
              sizeof(FIXP_DBL) * MAX_FREQ_COEFFS);
  }
  C_ALLOC_SCRATCH_END(_scratch, FIXP_DBL,
                      3 * MAX_NO_OF_ESTIMATES / 2 * MAX_FREQ_COEFFS)
}

INT FDKsbrEnc_CreateSbrMissingHarmonicsDetector(
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR hSbrMHDet, INT chan) {
  HANDLE_SBR_MISSING_HARMONICS_DETECTOR hs = hSbrMHDet;
  INT i;

  UCHAR *detectionVectors = GetRam_Sbr_detectionVectors(chan);
  UCHAR *guideVectorDetected = GetRam_Sbr_guideVectorDetected(chan);
  FIXP_DBL *guideVectorDiff = GetRam_Sbr_guideVectorDiff(chan);
  FIXP_DBL *guideVectorOrig = GetRam_Sbr_guideVectorOrig(chan);

  FDKmemclear(hs, sizeof(SBR_MISSING_HARMONICS_DETECTOR));

  hs->prevEnvelopeCompensation = GetRam_Sbr_prevEnvelopeCompensation(chan);
  hs->guideScfb = GetRam_Sbr_guideScfb(chan);

  if ((NULL == detectionVectors) || (NULL == guideVectorDetected) ||
      (NULL == guideVectorDiff) || (NULL == guideVectorOrig) ||
      (NULL == hs->prevEnvelopeCompensation) || (NULL == hs->guideScfb)) {
    goto bail;
  }

  for (i = 0; i < MAX_NO_OF_ESTIMATES; i++) {
    hs->guideVectors[i].guideVectorDiff =
        guideVectorDiff + (i * MAX_FREQ_COEFFS);
    hs->guideVectors[i].guideVectorOrig =
        guideVectorOrig + (i * MAX_FREQ_COEFFS);
    hs->detectionVectors[i] = detectionVectors + (i * MAX_FREQ_COEFFS);
    hs->guideVectors[i].guideVectorDetected =
        guideVectorDetected + (i * MAX_FREQ_COEFFS);
  }

  return 0;

bail:
  hs->guideVectors[0].guideVectorDiff = guideVectorDiff;
  hs->guideVectors[0].guideVectorOrig = guideVectorOrig;
  hs->detectionVectors[0] = detectionVectors;
  hs->guideVectors[0].guideVectorDetected = guideVectorDetected;

  FDKsbrEnc_DeleteSbrMissingHarmonicsDetector(hs);
  return -1;
}

INT FDKsbrEnc_InitSbrMissingHarmonicsDetector(
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR hSbrMHDet, INT sampleFreq,
    INT frameSize, INT nSfb, INT qmfNoChannels, INT totNoEst, INT move,
    INT noEstPerFrame, UINT sbrSyntaxFlags) {
  HANDLE_SBR_MISSING_HARMONICS_DETECTOR hs = hSbrMHDet;
  int i;

  FDK_ASSERT(totNoEst <= MAX_NO_OF_ESTIMATES);

  if (sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY) {
    switch (frameSize) {
      case 1024:
      case 512:
        hs->transientPosOffset = FRAME_MIDDLE_SLOT_512LD;
        hs->timeSlots = 16;
        break;
      case 960:
      case 480:
        hs->transientPosOffset = FRAME_MIDDLE_SLOT_512LD;
        hs->timeSlots = 15;
        break;
      default:
        return -1;
    }
  } else {
    switch (frameSize) {
      case 2048:
      case 1024:
        hs->transientPosOffset = FRAME_MIDDLE_SLOT_2048;
        hs->timeSlots = NUMBER_TIME_SLOTS_2048;
        break;
      case 1920:
      case 960:
        hs->transientPosOffset = FRAME_MIDDLE_SLOT_1920;
        hs->timeSlots = NUMBER_TIME_SLOTS_1920;
        break;
      default:
        return -1;
    }
  }

  if (sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY) {
    hs->mhParams = &paramsAacLd;
  } else
    hs->mhParams = &paramsAac;

  hs->qmfNoChannels = qmfNoChannels;
  hs->sampleFreq = sampleFreq;
  hs->nSfb = nSfb;

  hs->totNoEst = totNoEst;
  hs->move = move;
  hs->noEstPerFrame = noEstPerFrame;

  for (i = 0; i < totNoEst; i++) {
    FDKmemclear(hs->guideVectors[i].guideVectorDiff,
                sizeof(FIXP_DBL) * MAX_FREQ_COEFFS);
    FDKmemclear(hs->guideVectors[i].guideVectorOrig,
                sizeof(FIXP_DBL) * MAX_FREQ_COEFFS);
    FDKmemclear(hs->detectionVectors[i], sizeof(UCHAR) * MAX_FREQ_COEFFS);
    FDKmemclear(hs->guideVectors[i].guideVectorDetected,
                sizeof(UCHAR) * MAX_FREQ_COEFFS);
  }

  for (i = 0; i < MAX_NO_OF_ESTIMATES / 2; i++) {
    FDKmemclear(hs->tonalityDiff[i], sizeof(FIXP_DBL) * MAX_FREQ_COEFFS);
    FDKmemclear(hs->sfmOrig[i], sizeof(FIXP_DBL) * MAX_FREQ_COEFFS);
    FDKmemclear(hs->sfmSbr[i], sizeof(FIXP_DBL) * MAX_FREQ_COEFFS);
  }

  FDKmemclear(hs->prevEnvelopeCompensation, sizeof(UCHAR) * MAX_FREQ_COEFFS);
  FDKmemclear(hs->guideScfb, sizeof(UCHAR) * MAX_FREQ_COEFFS);

  hs->previousTransientFlag = 0;
  hs->previousTransientFrame = 0;
  hs->previousTransientPos = 0;

  return (0);
}

void FDKsbrEnc_DeleteSbrMissingHarmonicsDetector(
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR hSbrMHDet) {
  if (hSbrMHDet) {
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR hs = hSbrMHDet;

    FreeRam_Sbr_detectionVectors(&hs->detectionVectors[0]);
    FreeRam_Sbr_guideVectorDetected(&hs->guideVectors[0].guideVectorDetected);
    FreeRam_Sbr_guideVectorDiff(&hs->guideVectors[0].guideVectorDiff);
    FreeRam_Sbr_guideVectorOrig(&hs->guideVectors[0].guideVectorOrig);
    FreeRam_Sbr_prevEnvelopeCompensation(&hs->prevEnvelopeCompensation);
    FreeRam_Sbr_guideScfb(&hs->guideScfb);
  }
}

INT FDKsbrEnc_ResetSbrMissingHarmonicsDetector(
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR hSbrMissingHarmonicsDetector,
    INT nSfb) {
  int i;
  FIXP_DBL tempGuide[MAX_FREQ_COEFFS];
  UCHAR tempGuideInt[MAX_FREQ_COEFFS];
  INT nSfbPrev;

  nSfbPrev = hSbrMissingHarmonicsDetector->nSfb;
  hSbrMissingHarmonicsDetector->nSfb = nSfb;

  FDKmemcpy(tempGuideInt, hSbrMissingHarmonicsDetector->guideScfb,
            nSfbPrev * sizeof(UCHAR));

  if (nSfb > nSfbPrev) {
    for (i = 0; i < (nSfb - nSfbPrev); i++) {
      hSbrMissingHarmonicsDetector->guideScfb[i] = 0;
    }

    for (i = 0; i < nSfbPrev; i++) {
      hSbrMissingHarmonicsDetector->guideScfb[i + (nSfb - nSfbPrev)] =
          tempGuideInt[i];
    }
  } else {
    for (i = 0; i < nSfb; i++) {
      hSbrMissingHarmonicsDetector->guideScfb[i] =
          tempGuideInt[i + (nSfbPrev - nSfb)];
    }
  }

  FDKmemcpy(tempGuide,
            hSbrMissingHarmonicsDetector->guideVectors[0].guideVectorDiff,
            nSfbPrev * sizeof(FIXP_DBL));

  if (nSfb > nSfbPrev) {
    for (i = 0; i < (nSfb - nSfbPrev); i++) {
      hSbrMissingHarmonicsDetector->guideVectors[0].guideVectorDiff[i] =
          FL2FXCONST_DBL(0.0f);
    }

    for (i = 0; i < nSfbPrev; i++) {
      hSbrMissingHarmonicsDetector->guideVectors[0]
          .guideVectorDiff[i + (nSfb - nSfbPrev)] = tempGuide[i];
    }
  } else {
    for (i = 0; i < nSfb; i++) {
      hSbrMissingHarmonicsDetector->guideVectors[0].guideVectorDiff[i] =
          tempGuide[i + (nSfbPrev - nSfb)];
    }
  }

  FDKmemcpy(tempGuide,
            hSbrMissingHarmonicsDetector->guideVectors[0].guideVectorOrig,
            nSfbPrev * sizeof(FIXP_DBL));

  if (nSfb > nSfbPrev) {
    for (i = 0; i < (nSfb - nSfbPrev); i++) {
      hSbrMissingHarmonicsDetector->guideVectors[0].guideVectorOrig[i] =
          FL2FXCONST_DBL(0.0f);
    }

    for (i = 0; i < nSfbPrev; i++) {
      hSbrMissingHarmonicsDetector->guideVectors[0]
          .guideVectorOrig[i + (nSfb - nSfbPrev)] = tempGuide[i];
    }
  } else {
    for (i = 0; i < nSfb; i++) {
      hSbrMissingHarmonicsDetector->guideVectors[0].guideVectorOrig[i] =
          tempGuide[i + (nSfbPrev - nSfb)];
    }
  }

  FDKmemcpy(tempGuideInt,
            hSbrMissingHarmonicsDetector->guideVectors[0].guideVectorDetected,
            nSfbPrev * sizeof(UCHAR));

  if (nSfb > nSfbPrev) {
    for (i = 0; i < (nSfb - nSfbPrev); i++) {
      hSbrMissingHarmonicsDetector->guideVectors[0].guideVectorDetected[i] = 0;
    }

    for (i = 0; i < nSfbPrev; i++) {
      hSbrMissingHarmonicsDetector->guideVectors[0]
          .guideVectorDetected[i + (nSfb - nSfbPrev)] = tempGuideInt[i];
    }
  } else {
    for (i = 0; i < nSfb; i++) {
      hSbrMissingHarmonicsDetector->guideVectors[0].guideVectorDetected[i] =
          tempGuideInt[i + (nSfbPrev - nSfb)];
    }
  }

  FDKmemcpy(tempGuideInt,
            hSbrMissingHarmonicsDetector->prevEnvelopeCompensation,
            nSfbPrev * sizeof(UCHAR));

  if (nSfb > nSfbPrev) {
    for (i = 0; i < (nSfb - nSfbPrev); i++) {
      hSbrMissingHarmonicsDetector->prevEnvelopeCompensation[i] = 0;
    }

    for (i = 0; i < nSfbPrev; i++) {
      hSbrMissingHarmonicsDetector
          ->prevEnvelopeCompensation[i + (nSfb - nSfbPrev)] = tempGuideInt[i];
    }
  } else {
    for (i = 0; i < nSfb; i++) {
      hSbrMissingHarmonicsDetector->prevEnvelopeCompensation[i] =
          tempGuideInt[i + (nSfbPrev - nSfb)];
    }
  }

  return 0;
}
