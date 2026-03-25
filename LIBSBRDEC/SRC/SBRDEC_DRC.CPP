/**************************** SBR decoder library ******************************

   Author(s):   Christian Griebel

   Description: Dynamic range control (DRC) decoder tool for SBR

*******************************************************************************/

#include "sbrdec_drc.h"

static const UCHAR winBorderToColMappingTab[2][16] = {
    {0, 0, 4, 8, 12, 16, 20, 24, 28, 32, 32, 32, 32, 32, 32,
     32},
    {0, 0, 4, 8, 11, 15, 19, 23, 26, 30, 30, 30, 30, 30, 30,
     30}
};

void sbrDecoder_drcInitChannel(HANDLE_SBR_DRC_CHANNEL hDrcData) {
  int band;
  if (hDrcData == NULL) {
    return;
  }

  for (band = 0; band < (64); band++) {
    hDrcData->prevFact_mag[band] = FL2FXCONST_DBL(0.5f);
  }
  for (band = 0; band < SBRDEC_MAX_DRC_BANDS; band++) {
    hDrcData->currFact_mag[band] = FL2FXCONST_DBL(0.5f);
    hDrcData->nextFact_mag[band] = FL2FXCONST_DBL(0.5f);
  }

  hDrcData->prevFact_exp = 1;
  hDrcData->currFact_exp = 1;
  hDrcData->nextFact_exp = 1;
  hDrcData->numBandsCurr = 1;
  hDrcData->numBandsNext = 1;

  hDrcData->winSequenceCurr = 0;
  hDrcData->winSequenceNext = 0;

  hDrcData->drcInterpolationSchemeCurr = 0;
  hDrcData->drcInterpolationSchemeNext = 0;

  hDrcData->enable = 0;
}

void sbrDecoder_drcUpdateChannel(HANDLE_SBR_DRC_CHANNEL hDrcData) {
  if (hDrcData == NULL) {
    return;
  }
  if (hDrcData->enable != 1) {
    return;
  }

  FDKmemcpy(hDrcData->currFact_mag, hDrcData->nextFact_mag,
            SBRDEC_MAX_DRC_BANDS * sizeof(FIXP_DBL));

  hDrcData->currFact_exp = hDrcData->nextFact_exp;
  hDrcData->numBandsCurr = hDrcData->numBandsNext;

  FDKmemcpy(hDrcData->bandTopCurr, hDrcData->bandTopNext,
            SBRDEC_MAX_DRC_BANDS * sizeof(USHORT));

  hDrcData->drcInterpolationSchemeCurr = hDrcData->drcInterpolationSchemeNext;
  hDrcData->winSequenceCurr = hDrcData->winSequenceNext;
}

void sbrDecoder_drcApplySlot(HANDLE_SBR_DRC_CHANNEL hDrcData,
                             FIXP_DBL *qmfRealSlot, FIXP_DBL *qmfImagSlot,
                             int col, int numQmfSubSamples, int maxShift) {
  const UCHAR *winBorderToColMap;

  int band, bottomMdct, topMdct, bin, useLP;
  int indx = numQmfSubSamples - (numQmfSubSamples >> 1) - 10;
  int frameLenFlag = (numQmfSubSamples == 30) ? 1 : 0;
  int frameSize = (frameLenFlag == 1) ? 960 : 1024;

  const FIXP_DBL *fact_mag = NULL;
  INT fact_exp = 0;
  UINT numBands = 0;
  USHORT *bandTop = NULL;
  int shortDrc = 0;

  FIXP_DBL alphaValue = FL2FXCONST_DBL(0.0f);

  if (hDrcData == NULL) {
    return;
  }
  if (hDrcData->enable != 1) {
    return;
  }

  winBorderToColMap = winBorderToColMappingTab[frameLenFlag];

  useLP = (qmfImagSlot == NULL) ? 1 : 0;

  col += indx;
  bottomMdct = 0;

  if (col < (numQmfSubSamples >> 1)) {
    if (hDrcData->winSequenceCurr != 2) {
      int j = col + (numQmfSubSamples >> 1);
      if (j < winBorderToColMap[15]) {
        if (hDrcData->drcInterpolationSchemeCurr == 0) {
          INT k = (frameLenFlag) ? 0x4444445 : 0x4000000;
          alphaValue = (FIXP_DBL)(j * k);
        } else {
          if (j >=
              (int)winBorderToColMap[hDrcData->drcInterpolationSchemeCurr]) {
            alphaValue = (FIXP_DBL)MAXVAL_DBL;
          }
        }
      } else {
        alphaValue = (FIXP_DBL)MAXVAL_DBL;
      }
    } else {
      shortDrc = 1;
    }

    fact_mag = hDrcData->currFact_mag;
    fact_exp = hDrcData->currFact_exp;
    numBands = hDrcData->numBandsCurr;
    bandTop = hDrcData->bandTopCurr;
  } else if (col < numQmfSubSamples) {
    if (hDrcData->winSequenceNext != 2) {
      int j = col - (numQmfSubSamples >> 1);

      if (j < winBorderToColMap[15]) {
        if (hDrcData->drcInterpolationSchemeNext == 0) {
          INT k = (frameLenFlag) ? 0x4444445 : 0x4000000;
          alphaValue = (FIXP_DBL)(j * k);
        } else {
          if (j >=
              (int)winBorderToColMap[hDrcData->drcInterpolationSchemeNext]) {
            alphaValue = (FIXP_DBL)MAXVAL_DBL;
          }
        }
      } else {
        alphaValue = (FIXP_DBL)MAXVAL_DBL;
      }

      fact_mag = hDrcData->nextFact_mag;
      fact_exp = hDrcData->nextFact_exp;
      numBands = hDrcData->numBandsNext;
      bandTop = hDrcData->bandTopNext;
    } else {
      if (hDrcData->winSequenceCurr != 2) {
        alphaValue = (FIXP_DBL)0;
        fact_mag = hDrcData->nextFact_mag;
        fact_exp = hDrcData->nextFact_exp;
        numBands = hDrcData->numBandsNext;
        bandTop = hDrcData->bandTopNext;
      } else {
        shortDrc = 1;
        fact_mag = hDrcData->currFact_mag;
        fact_exp = hDrcData->currFact_exp;
        numBands = hDrcData->numBandsCurr;
        bandTop = hDrcData->bandTopCurr;
      }
    }
  } else {
    if (hDrcData->winSequenceNext != 2) {
      int j = col - (numQmfSubSamples >> 1);

      if (j < winBorderToColMap[15]) {
        if (hDrcData->drcInterpolationSchemeNext == 0) {
          INT k = (frameLenFlag) ? 0x4444445 : 0x4000000;
          alphaValue = (FIXP_DBL)(j * k);
        } else {
          if (j >=
              (int)winBorderToColMap[hDrcData->drcInterpolationSchemeNext]) {
            alphaValue = (FIXP_DBL)MAXVAL_DBL;
          }
        }
      } else {
        alphaValue = (FIXP_DBL)MAXVAL_DBL;
      }
    } else {
      shortDrc = 1;
    }

    fact_mag = hDrcData->nextFact_mag;
    fact_exp = hDrcData->nextFact_exp;
    numBands = hDrcData->numBandsNext;
    bandTop = hDrcData->bandTopNext;

    col -= numQmfSubSamples;
  }

  for (band = 0; band < (int)numBands; band++) {
    int bottomQmf, topQmf;

    FIXP_DBL drcFact_mag = (FIXP_DBL)MAXVAL_DBL;

    topMdct = (bandTop[band] + 1) << 2;

    if (!shortDrc) {
      if (frameLenFlag) {
        bottomQmf = fMultIfloor((FIXP_DBL)0x4444445, bottomMdct);
        topQmf = fMultIfloor((FIXP_DBL)0x4444445, topMdct);

        topMdct = 30 * topQmf;
      } else {
        topMdct &= ~0x1f;
        bottomQmf = bottomMdct >> 5;
        topQmf = topMdct >> 5;
      }

      if (band == ((int)numBands - 1)) {
        topQmf = (64);
      }

      for (bin = bottomQmf; bin < topQmf; bin++) {
        FIXP_DBL drcFact1_mag = hDrcData->prevFact_mag[bin];
        FIXP_DBL drcFact2_mag = fact_mag[band];

        if (hDrcData->prevFact_exp < maxShift) {
          drcFact1_mag >>= maxShift - hDrcData->prevFact_exp;
        }
        if (fact_exp < maxShift) {
          drcFact2_mag >>= maxShift - fact_exp;
        }
        if (alphaValue == (FIXP_DBL)0) {
          drcFact_mag = drcFact1_mag;
        } else if (alphaValue == (FIXP_DBL)MAXVAL_DBL) {
          drcFact_mag = drcFact2_mag;
        } else {
          drcFact_mag =
              fMult(alphaValue, drcFact2_mag) +
              fMult(((FIXP_DBL)MAXVAL_DBL - alphaValue), drcFact1_mag);
        }

        qmfRealSlot[bin] = fMult(qmfRealSlot[bin], drcFact_mag);
        if (!useLP) {
          qmfImagSlot[bin] = fMult(qmfImagSlot[bin], drcFact_mag);
        }
        if (col == (numQmfSubSamples >> 1) - 1) {
          hDrcData->prevFact_mag[bin] = fact_mag[band];
        }
      }
    } else {
      unsigned startWinIdx, stopWinIdx;
      int startCol, stopCol;
      FIXP_DBL invFrameSizeDiv8 =
          (frameLenFlag) ? (FIXP_DBL)0x1111112 : (FIXP_DBL)0x1000000;

      if (topMdct < 0) {
        topMdct = 0;
      }
      if (topMdct >= frameSize) {
        topMdct = frameSize - 1;
      }
      if (frameLenFlag) {
        topMdct = fMultIfloor((FIXP_DBL)0x78000000,
                              fMultIfloor((FIXP_DBL)0x22222223, topMdct) << 2);

        startWinIdx = fMultIfloor(invFrameSizeDiv8, bottomMdct) +
                      1;
        stopWinIdx = fMultIceil(invFrameSizeDiv8 - (FIXP_DBL)1, topMdct) + 1;
      } else {
        topMdct &= ~0x03;

        startWinIdx = fMultIfloor(invFrameSizeDiv8, bottomMdct) + 1;
        stopWinIdx = fMultIceil(invFrameSizeDiv8, topMdct) + 1;
      }

      startCol = (int)winBorderToColMap[startWinIdx];
      stopCol = (int)winBorderToColMap[stopWinIdx];

      bottomQmf = fMultIfloor(invFrameSizeDiv8,
                              ((bottomMdct % (numQmfSubSamples << 2)) << 5));
      topQmf = fMultIfloor(invFrameSizeDiv8,
                           ((topMdct % (numQmfSubSamples << 2)) << 5));

      if (band == ((int)numBands - 1)) {
        topQmf = (64);
        stopCol = numQmfSubSamples;
        stopWinIdx = 10;
      }

      if (topQmf == 0) {
        if (frameLenFlag) {
          FIXP_DBL rem = fMult(invFrameSizeDiv8,
                               (FIXP_DBL)(topMdct << (DFRACT_BITS - 12)));
          if ((LONG)rem & (LONG)0x1F) {
            stopWinIdx -= 1;
            stopCol = (int)winBorderToColMap[stopWinIdx];
          }
        }
        topQmf = (64);
      }

      if (stopCol == numQmfSubSamples) {
        int tmpBottom = bottomQmf;
        if ((int)winBorderToColMap[8] > startCol) {
          tmpBottom = 0;
        }

        for (bin = tmpBottom; bin < topQmf; bin++) {
          hDrcData->prevFact_mag[bin] = fact_mag[band];
        }
      }

      if ((col >= startCol) && (col < stopCol)) {
        if (col >= (int)winBorderToColMap[startWinIdx + 1]) {
          bottomQmf = 0;
        }
        if (col < (int)winBorderToColMap[stopWinIdx - 1]) {
          topQmf = (64);
        }

        drcFact_mag = fact_mag[band];

        if (fact_exp < maxShift) {
          drcFact_mag >>= maxShift - fact_exp;
        }

        for (bin = bottomQmf; bin < topQmf; bin++) {
          qmfRealSlot[bin] = fMult(qmfRealSlot[bin], drcFact_mag);
          if (!useLP) {
            qmfImagSlot[bin] = fMult(qmfImagSlot[bin], drcFact_mag);
          }
        }
      }
    }

    bottomMdct = topMdct;
  }

  if (col == (numQmfSubSamples >> 1) - 1) {
    hDrcData->prevFact_exp = fact_exp;
  }
}

void sbrDecoder_drcApply(HANDLE_SBR_DRC_CHANNEL hDrcData,
                         FIXP_DBL **QmfBufferReal, FIXP_DBL **QmfBufferImag,
                         int numQmfSubSamples, int *scaleFactor) {
  int col;
  int maxShift = 0;

  if (hDrcData == NULL) {
    return;
  }
  if (hDrcData->enable == 0) {
    return;
  }
  if (hDrcData->prevFact_exp > maxShift) {
    maxShift = hDrcData->prevFact_exp;
  }
  if (hDrcData->currFact_exp > maxShift) {
    maxShift = hDrcData->currFact_exp;
  }
  if (hDrcData->nextFact_exp > maxShift) {
    maxShift = hDrcData->nextFact_exp;
  }

  for (col = 0; col < numQmfSubSamples; col++) {
    FIXP_DBL *qmfSlotReal = QmfBufferReal[col];
    FIXP_DBL *qmfSlotImag = (QmfBufferImag == NULL) ? NULL : QmfBufferImag[col];

    sbrDecoder_drcApplySlot(hDrcData, qmfSlotReal, qmfSlotImag, col,
                            numQmfSubSamples, maxShift);
  }

  *scaleFactor += maxShift;
}
