/**************************** AAC decoder library ******************************

   Author(s):   Robert Weidner (DSP Solutions)

   Description: HCR Decoder: HCR initialization, preprocess HCR sideinfo,
                decode priority codewords (PCWs)

*******************************************************************************/

#include "aacdec_hcr.h"
#include "aacdec_hcr_types.h"
#include "aacdec_hcr_bit.h"
#include "aacdec_hcrs.h"
#include "aac_ram.h"
#include "aac_rom.h"
#include "channel.h"
#include "block.h"
#include "aacdecoder.h"
#include "FDK_bitstream.h"

extern int mlFileChCurr;

static void errDetectorInHcrSideinfoShrt(SCHAR cb, SHORT numLine,
                                         UINT *errorWord);

static void errDetectorInHcrLengths(SCHAR lengthOfLongestCodeword,
                                    SHORT lengthOfReorderedSpectralData,
                                    UINT *errorWord);

static void HcrCalcNumCodeword(H_HCR_INFO pHcr);
static void HcrSortCodebookAndNumCodewordInSection(H_HCR_INFO pHcr);
static void HcrPrepareSegmentationGrid(H_HCR_INFO pHcr);
static void HcrExtendedSectionInfo(H_HCR_INFO pHcr);

static void DeriveNumberOfExtendedSortedSectionsInSets(
    UINT numSegment, USHORT *pNumExtendedSortedCodewordInSection,
    int numExtendedSortedCodewordInSectionIdx,
    USHORT *pNumExtendedSortedSectionsInSets,
    int numExtendedSortedSectionsInSetsIdx);

static INT DecodeEscapeSequence(HANDLE_FDK_BITSTREAM bs, const INT bsAnchor,
                                INT quantSpecCoef, INT *pLeftStartOfSegment,
                                SCHAR *pRemainingBitsInSegment,
                                int *pNumDecodedBits, UINT *errorWord);

static int DecodePCW_Sign(HANDLE_FDK_BITSTREAM bs, const INT bsAnchor,
                          UINT codebookDim, const SCHAR *pQuantVal,
                          FIXP_DBL *pQuantSpecCoef, int *quantSpecCoefIdx,
                          INT *pLeftStartOfSegment,
                          SCHAR *pRemainingBitsInSegment, int *pNumDecodedBits);

static const SCHAR *DecodePCW_Body(HANDLE_FDK_BITSTREAM bs, const INT bsAnchor,
                                   const UINT *pCurrentTree,
                                   const SCHAR *pQuantValBase,
                                   INT *pLeftStartOfSegment,
                                   SCHAR *pRemainingBitsInSegment,
                                   int *pNumDecodedBits);

static void DecodePCWs(HANDLE_FDK_BITSTREAM bs, H_HCR_INFO pHcr);

static void HcrReorderQuantizedSpectralCoefficients(
    H_HCR_INFO pHcr, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo);

static UCHAR errDetectPcwSegmentation(SCHAR remainingBitsInSegment,
                                      H_HCR_INFO pHcr, PCW_TYPE kind,
                                      FIXP_DBL *qsc_base_of_cw,
                                      UCHAR dimension);

static void errDetectWithinSegmentationFinal(H_HCR_INFO pHcr);

static void errDetectorInHcrSideinfoShrt(SCHAR cb, SHORT numLine,
                                         UINT *errorWord) {
  if (cb < ZERO_HCB || cb >= MAX_CB_CHECK || cb == BOOKSCL) {
    *errorWord |= CB_OUT_OF_RANGE_SHORT_BLOCK;
  }
  if (numLine < 0 || numLine > 1024) {
    *errorWord |= LINE_IN_SECT_OUT_OF_RANGE_SHORT_BLOCK;
  }
}

static void errDetectorInHcrLengths(SCHAR lengthOfLongestCodeword,
                                    SHORT lengthOfReorderedSpectralData,
                                    UINT *errorWord) {
  if (lengthOfReorderedSpectralData < lengthOfLongestCodeword) {
    *errorWord |= HCR_SI_LENGTHS_FAILURE;
  }
}

void CHcr_Read(HANDLE_FDK_BITSTREAM bs,
               CAacDecoderChannelInfo *pAacDecoderChannelInfo,
               const MP4_ELEMENT_ID globalHcrType) {
  SHORT lengOfReorderedSpectralData;
  SCHAR lengOfLongestCodeword;

  pAacDecoderChannelInfo->pDynData->specificTo.aac.lenOfReorderedSpectralData =
      0;
  pAacDecoderChannelInfo->pDynData->specificTo.aac.lenOfLongestCodeword = 0;

  lengOfReorderedSpectralData = FDKreadBits(bs, 14) + ERROR_LORSD;
  if (globalHcrType == ID_CPE) {
    if ((lengOfReorderedSpectralData >= 0) &&
        (lengOfReorderedSpectralData <= CPE_TOP_LENGTH)) {
      pAacDecoderChannelInfo->pDynData->specificTo.aac
          .lenOfReorderedSpectralData =
          lengOfReorderedSpectralData;
    } else {
      if (lengOfReorderedSpectralData > CPE_TOP_LENGTH) {
        pAacDecoderChannelInfo->pDynData->specificTo.aac
            .lenOfReorderedSpectralData =
            CPE_TOP_LENGTH;
      }
    }
  } else if (globalHcrType == ID_SCE || globalHcrType == ID_LFE ||
             globalHcrType == ID_CCE) {
    if ((lengOfReorderedSpectralData >= 0) &&
        (lengOfReorderedSpectralData <= SCE_TOP_LENGTH)) {
      pAacDecoderChannelInfo->pDynData->specificTo.aac
          .lenOfReorderedSpectralData =
          lengOfReorderedSpectralData;
    } else {
      if (lengOfReorderedSpectralData > SCE_TOP_LENGTH) {
        pAacDecoderChannelInfo->pDynData->specificTo.aac
            .lenOfReorderedSpectralData =
            SCE_TOP_LENGTH;
      }
    }
  }

  lengOfLongestCodeword = FDKreadBits(bs, 6) + ERROR_LOLC;
  if ((lengOfLongestCodeword >= 0) &&
      (lengOfLongestCodeword <= LEN_OF_LONGEST_CW_TOP_LENGTH)) {
    pAacDecoderChannelInfo->pDynData->specificTo.aac.lenOfLongestCodeword =
        lengOfLongestCodeword;
  } else {
    if (lengOfLongestCodeword > LEN_OF_LONGEST_CW_TOP_LENGTH) {
      pAacDecoderChannelInfo->pDynData->specificTo.aac.lenOfLongestCodeword =
          LEN_OF_LONGEST_CW_TOP_LENGTH;
    }
  }
}

UINT HcrInit(H_HCR_INFO pHcr, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
             const SamplingRateInfo *pSamplingRateInfo,
             HANDLE_FDK_BITSTREAM bs) {
  CIcsInfo *pIcsInfo = &pAacDecoderChannelInfo->icsInfo;
  SHORT *pNumLinesInSec;
  UCHAR *pCodeBk;
  SHORT numSection;
  SCHAR cb;
  int numLine;
  int i;

  pHcr->decInOut.lengthOfReorderedSpectralData =
      pAacDecoderChannelInfo->pDynData->specificTo.aac
          .lenOfReorderedSpectralData;
  pHcr->decInOut.lengthOfLongestCodeword =
      pAacDecoderChannelInfo->pDynData->specificTo.aac.lenOfLongestCodeword;
  pHcr->decInOut.pQuantizedSpectralCoefficientsBase =
      pAacDecoderChannelInfo->pSpectralCoefficient;
  pHcr->decInOut.quantizedSpectralCoefficientsIdx = 0;
  pHcr->decInOut.pCodebook =
      pAacDecoderChannelInfo->pDynData->specificTo.aac.aCodeBooks4Hcr;
  pHcr->decInOut.pNumLineInSect =
      pAacDecoderChannelInfo->pDynData->specificTo.aac.aNumLineInSec4Hcr;
  pHcr->decInOut.numSection =
      pAacDecoderChannelInfo->pDynData->specificTo.aac.numberSection;
  pHcr->decInOut.errorLog = 0;
  pHcr->nonPcwSideinfo.pResultBase =
      SPEC_LONG(pAacDecoderChannelInfo->pSpectralCoefficient);

  FDKsyncCache(bs);
  pHcr->decInOut.bitstreamAnchor = (INT)FDKgetValidBits(bs);

  if (!IsLongBlock(&pAacDecoderChannelInfo->icsInfo))
  {
    SHORT band;
    SHORT maxBand;
    SCHAR group;
    SCHAR winGroupLen;
    SCHAR window;
    SCHAR numUnitInBand;
    SCHAR cntUnitInBand;
    SCHAR groupWin;
    SCHAR cb_prev;

    UCHAR *pCodeBook;
    const SHORT *BandOffsets;
    SCHAR numOfGroups;

    pCodeBook = pAacDecoderChannelInfo->pDynData->aCodeBook;
    pNumLinesInSec = pHcr->decInOut.pNumLineInSect;
    pCodeBk = pHcr->decInOut.pCodebook;
    BandOffsets =
        GetScaleFactorBandOffsets(pIcsInfo, pSamplingRateInfo);
    numOfGroups = GetWindowGroups(pIcsInfo);

    numLine = 0;
    numSection = 0;
    cb = pCodeBook[0];
    cb_prev = pCodeBook[0];

    *pCodeBk++ = cb_prev;

    maxBand = GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo->icsInfo);
    for (band = 0; band < maxBand;
         band++) {
      numUnitInBand =
          ((BandOffsets[band + 1] - BandOffsets[band]) >>
           FOUR_LOG_DIV_TWO_LOG);
      for (cntUnitInBand = numUnitInBand; cntUnitInBand != 0;
           cntUnitInBand--) {
        for (window = 0, group = 0; group < numOfGroups; group++) {
          winGroupLen = (SCHAR)GetWindowGroupLength(
              &pAacDecoderChannelInfo->icsInfo, group);
          for (groupWin = winGroupLen; groupWin != 0; groupWin--, window++) {
            cb = pCodeBook[group * 16 + band];
            if (cb != cb_prev) {
              errDetectorInHcrSideinfoShrt(cb, numLine,
                                           &pHcr->decInOut.errorLog);
              if (pHcr->decInOut.errorLog != 0) {
                return (pHcr->decInOut.errorLog);
              }
              *pCodeBk++ = cb;
              *pNumLinesInSec++ = numLine;
              numSection++;

              cb_prev = cb;
              numLine = LINES_PER_UNIT;
            } else {
              numLine += LINES_PER_UNIT;
            }
          }
        }
      }
    }

    numSection++;

    errDetectorInHcrSideinfoShrt(cb, numLine, &pHcr->decInOut.errorLog);
    errDetectorInHcrLengths(pHcr->decInOut.lengthOfLongestCodeword,
                            pHcr->decInOut.lengthOfReorderedSpectralData,
                            &pHcr->decInOut.errorLog);
    if (pHcr->decInOut.errorLog != 0) {
      return (pHcr->decInOut.errorLog);
    }

    *pCodeBk = cb;
    *pNumLinesInSec = numLine;
    pHcr->decInOut.numSection = numSection;

  } else
  {
    errDetectorInHcrLengths(pHcr->decInOut.lengthOfLongestCodeword,
                            pHcr->decInOut.lengthOfReorderedSpectralData,
                            &pHcr->decInOut.errorLog);
    numSection = pHcr->decInOut.numSection;
    pNumLinesInSec = pHcr->decInOut.pNumLineInSect;
    pCodeBk = pHcr->decInOut.pCodebook;
    if (numSection <= 0 || numSection > 64) {
      numSection = 0;
    }

    for (i = numSection; i != 0; i--) {
      cb = *pCodeBk++;

      numLine = *pNumLinesInSec++;
    }
    if (pHcr->decInOut.errorLog != 0) {
      return (pHcr->decInOut.errorLog);
    }
  }

  pCodeBk = pHcr->decInOut.pCodebook;
  for (i = 0; i < numSection; i++) {
    if ((*pCodeBk == NOISE_HCB) || (*pCodeBk == INTENSITY_HCB2) ||
        (*pCodeBk == INTENSITY_HCB)) {
      *pCodeBk = 0;
    }
    pCodeBk++;
  }

  return (pHcr->decInOut.errorLog);
}

UINT HcrDecoder(H_HCR_INFO pHcr, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                const SamplingRateInfo *pSamplingRateInfo,
                HANDLE_FDK_BITSTREAM bs) {
  int pTmp1, pTmp2, pTmp3, pTmp4;
  int pTmp5;

  INT bitCntOffst;
  INT saveBitCnt = (INT)FDKgetValidBits(bs);

  HcrCalcNumCodeword(pHcr);

  HcrSortCodebookAndNumCodewordInSection(pHcr);

  HcrPrepareSegmentationGrid(pHcr);

  HcrExtendedSectionInfo(pHcr);

  if ((pHcr->decInOut.errorLog & HCR_FATAL_PCW_ERROR_MASK) != 0) {
    return (pHcr->decInOut.errorLog);
  }

  DeriveNumberOfExtendedSortedSectionsInSets(
      pHcr->segmentInfo.numSegment,
      pHcr->sectionInfo.pNumExtendedSortedCodewordInSection,
      pHcr->sectionInfo.numExtendedSortedCodewordInSectionIdx,
      pHcr->sectionInfo.pNumExtendedSortedSectionsInSets,
      pHcr->sectionInfo.numExtendedSortedSectionsInSetsIdx);

  pTmp1 = pHcr->sectionInfo.numExtendedSortedCodewordInSectionIdx;
  pTmp2 = pHcr->sectionInfo.extendedSortedCodebookIdx;
  pTmp3 = pHcr->sectionInfo.numExtendedSortedSectionsInSetsIdx;
  pTmp4 = pHcr->decInOut.quantizedSpectralCoefficientsIdx;
  pTmp5 = pHcr->sectionInfo.maxLenOfCbInExtSrtSecIdx;

  DecodePCWs(bs, pHcr);

  if ((pHcr->decInOut.errorLog & HCR_FATAL_PCW_ERROR_MASK) == 0) {
    DecodeNonPCWs(bs, pHcr);
  }

  errDetectWithinSegmentationFinal(pHcr);

  pHcr->sectionInfo.numExtendedSortedCodewordInSectionIdx = pTmp1;
  pHcr->sectionInfo.extendedSortedCodebookIdx = pTmp2;
  pHcr->sectionInfo.numExtendedSortedSectionsInSetsIdx = pTmp3;
  pHcr->decInOut.quantizedSpectralCoefficientsIdx = pTmp4;
  pHcr->sectionInfo.maxLenOfCbInExtSrtSecIdx = pTmp5;

  HcrReorderQuantizedSpectralCoefficients(pHcr, pAacDecoderChannelInfo,
                                          pSamplingRateInfo);

  bitCntOffst = (INT)FDKgetValidBits(bs) - saveBitCnt;
  if (bitCntOffst) {
    FDKpushBiDirectional(bs, bitCntOffst);
  }

  return (pHcr->decInOut.errorLog);
}

static void HcrReorderQuantizedSpectralCoefficients(
    H_HCR_INFO pHcr, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo) {
  INT qsc;
  UINT abs_qsc;
  UINT i, j;
  USHORT numSpectralValuesInSection;
  FIXP_DBL *pTeVa;
  USHORT lavErrorCnt = 0;

  UINT numSection = pHcr->decInOut.numSection;
  SPECTRAL_PTR pQuantizedSpectralCoefficientsBase =
      pHcr->decInOut.pQuantizedSpectralCoefficientsBase;
  FIXP_DBL *pQuantizedSpectralCoefficients =
      SPEC_LONG(pHcr->decInOut.pQuantizedSpectralCoefficientsBase);
  const UCHAR *pCbDimShift = aDimCbShift;
  const USHORT *pLargestAbsVal = aLargestAbsoluteValue;
  UCHAR *pSortedCodebook = pHcr->sectionInfo.pSortedCodebook;
  USHORT *pNumSortedCodewordInSection =
      pHcr->sectionInfo.pNumSortedCodewordInSection;
  USHORT *pReorderOffset = pHcr->sectionInfo.pReorderOffset;
  FIXP_DBL pTempValues[1024];
  FIXP_DBL *pBak = pTempValues;

  FDKmemclear(pTempValues, 1024 * sizeof(FIXP_DBL));

  for (i = numSection; i != 0; i--) {
    numSpectralValuesInSection = *pNumSortedCodewordInSection++
                                 << pCbDimShift[*pSortedCodebook];
    pTeVa = &pTempValues[*pReorderOffset++];
    for (j = numSpectralValuesInSection; j != 0; j--) {
      qsc = *pQuantizedSpectralCoefficients++;
      abs_qsc = fAbs(qsc);
      if (abs_qsc <= pLargestAbsVal[*pSortedCodebook]) {
        *pTeVa++ = (FIXP_DBL)qsc;
      } else {
        if (abs_qsc ==
            Q_VALUE_INVALID) {
          *pTeVa++ = (FIXP_DBL)qsc;
        } else {
          *pTeVa++ = (FIXP_DBL)Q_VALUE_INVALID;
          lavErrorCnt += 1;
        }
      }
    }
    pSortedCodebook++;
  }

  if (!IsLongBlock(&pAacDecoderChannelInfo->icsInfo)) {
    FIXP_DBL *pOut;
    FIXP_DBL locMax;
    FIXP_DBL tmp;
    SCHAR groupoffset;
    SCHAR group;
    SCHAR band;
    SCHAR groupwin;
    SCHAR window;
    SCHAR numWinGroup;
    SHORT interm;
    SCHAR numSfbTransm;
    SCHAR winGroupLen;
    SHORT index;
    INT msb;
    INT lsb;

    SHORT *pScaleFacHcr = pAacDecoderChannelInfo->pDynData->aScaleFactor;
    SHORT *pSfbSclHcr = pAacDecoderChannelInfo->pDynData->aSfbScale;
    const SHORT *BandOffsets = GetScaleFactorBandOffsets(
        &pAacDecoderChannelInfo->icsInfo, pSamplingRateInfo);

    pBak = pTempValues;
    for (window = 0; window < (8); window++) {
      pOut = SPEC(pQuantizedSpectralCoefficientsBase, window,
                  pAacDecoderChannelInfo->granuleLength);
      for (i = 0; i < (LINES_PER_UNIT_GROUP); i++) {
        pTeVa = pBak + (window << FOUR_LOG_DIV_TWO_LOG) +
                i * 32;
        for (j = (LINES_PER_UNIT); j != 0; j--) {
          *pOut++ = *pTeVa++;
        }
      }
    }

    groupoffset = 0;

    numWinGroup = GetWindowGroups(&pAacDecoderChannelInfo->icsInfo);
    numSfbTransm =
        GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo->icsInfo);

    for (group = 0; group < numWinGroup; group++) {
      winGroupLen =
          GetWindowGroupLength(&pAacDecoderChannelInfo->icsInfo, group);
      for (band = 0; band < numSfbTransm; band++) {
        interm = group * 16 + band;
        msb = pScaleFacHcr[interm] >> 2;
        lsb = pScaleFacHcr[interm] & 3;
        for (groupwin = 0; groupwin < winGroupLen; groupwin++) {
          window = groupoffset + groupwin;
          pBak = SPEC(pQuantizedSpectralCoefficientsBase, window,
                      pAacDecoderChannelInfo->granuleLength);
          locMax = FL2FXCONST_DBL(0.0f);
          for (index = BandOffsets[band]; index < BandOffsets[band + 1];
               index += LINES_PER_UNIT) {
            pTeVa = &pBak[index];
            for (i = LINES_PER_UNIT; i != 0; i--) {
              tmp = (*pTeVa < FL2FXCONST_DBL(0.0f)) ? -*pTeVa++ : *pTeVa++;
              locMax = fixMax(tmp, locMax);
            }
          }
          if (fixp_abs(locMax) > (FIXP_DBL)MAX_QUANTIZED_VALUE) {
            locMax = (FIXP_DBL)MAX_QUANTIZED_VALUE;
          }
          pSfbSclHcr[window * 16 + band] =
              msb - GetScaleFromValue(
                        locMax, lsb);
        }
      }
      groupoffset +=
          GetWindowGroupLength(&pAacDecoderChannelInfo->icsInfo, group);
    }
  } else {
    pQuantizedSpectralCoefficients =
        SPEC_LONG(pQuantizedSpectralCoefficientsBase);
    for (i = 1024; i != 0; i--) {
      *pQuantizedSpectralCoefficients++ = *pBak++;
    }
  }
}

static void HcrCalcNumCodeword(H_HCR_INFO pHcr) {
  int hcrSection;
  UINT numCodeword;

  UINT numSection = pHcr->decInOut.numSection;
  UCHAR *pCodebook = pHcr->decInOut.pCodebook;
  SHORT *pNumLineInSection = pHcr->decInOut.pNumLineInSect;
  const UCHAR *pCbDimShift = aDimCbShift;

  USHORT *pNumCodewordInSection = pHcr->sectionInfo.pNumCodewordInSection;

  numCodeword = 0;
  for (hcrSection = numSection; hcrSection != 0; hcrSection--) {
    *pNumCodewordInSection = *pNumLineInSection++ >> pCbDimShift[*pCodebook];
    if (*pCodebook != 0) {
      numCodeword += *pNumCodewordInSection;
    }
    pNumCodewordInSection++;
    pCodebook++;
  }
  pHcr->sectionInfo.numCodeword = numCodeword;
}

static void HcrSortCodebookAndNumCodewordInSection(H_HCR_INFO pHcr) {
  UINT i, j, k;
  UCHAR temp;
  UINT counter;
  UINT startOffset;
  UINT numZeroSection;
  UCHAR *pDest;
  UINT numSectionDec;

  UINT numSection = pHcr->decInOut.numSection;
  UCHAR *pCodebook = pHcr->decInOut.pCodebook;
  UCHAR *pSortedCodebook = pHcr->sectionInfo.pSortedCodebook;
  USHORT *pNumCodewordInSection = pHcr->sectionInfo.pNumCodewordInSection;
  USHORT *pNumSortedCodewordInSection =
      pHcr->sectionInfo.pNumSortedCodewordInSection;
  UCHAR *pCodebookSwitch = pHcr->sectionInfo.pCodebookSwitch;
  USHORT *pReorderOffset = pHcr->sectionInfo.pReorderOffset;
  const UCHAR *pCbPriority = aCbPriority;
  const UCHAR *pMinOfCbPair = aMinOfCbPair;
  const UCHAR *pMaxOfCbPair = aMaxOfCbPair;
  const UCHAR *pCbDimShift = aDimCbShift;

  UINT searchStart = 0;

  pDest = pSortedCodebook;
  numZeroSection = 0;
  for (i = numSection; i != 0; i--) {
    if (pCbPriority[*pCodebook] == 0) {
      numZeroSection += 1;
    }
    *pDest++ = pCbPriority[*pCodebook++];
  }
  pHcr->sectionInfo.numSortedSection =
      numSection - numZeroSection;
  pCodebook = pHcr->decInOut.pCodebook;

  numSectionDec = numSection - 1;
  if (numSectionDec > 0) {
    counter = numSectionDec;
    for (j = numSectionDec; j != 0; j--) {
      for (i = 0; i < counter; i++) {
        if (pSortedCodebook[i + 1] > pSortedCodebook[i]) {
          temp = pSortedCodebook[i];
          pSortedCodebook[i] = pSortedCodebook[i + 1];
          pSortedCodebook[i + 1] = temp;
        }
      }
      counter -= 1;
    }
  }

  for (i = numSection; i != 0; i--) {
    *pCodebookSwitch++ = 0;
  }
  pCodebookSwitch = pHcr->sectionInfo.pCodebookSwitch;

  for (j = 0; j < numSection; j++) {
    for (i = searchStart; i < numSection; i++) {
      if (pCodebookSwitch[i] == 0 &&
          (pMinOfCbPair[pSortedCodebook[j]] == pCodebook[i] ||
           pMaxOfCbPair[pSortedCodebook[j]] == pCodebook[i])) {
        pCodebookSwitch[i] = 1;
        pSortedCodebook[j] = pCodebook[i];
        pNumSortedCodewordInSection[j] =
            pNumCodewordInSection[i];

        startOffset = 0;
        for (k = 0; k < i; k++) {
          startOffset += pNumCodewordInSection[k] << pCbDimShift[pCodebook[k]];
        }
        pReorderOffset[j] =
            startOffset;

        if (i == searchStart) {
          k = i;
          while (pCodebookSwitch[k++] == 1) searchStart++;
        }
        break;
      }
    }
  }
}

static void HcrPrepareSegmentationGrid(H_HCR_INFO pHcr) {
  USHORT i, j;
  USHORT numSegment = 0;
  INT segmentStart = 0;
  UCHAR segmentWidth;
  UCHAR lastSegmentWidth;
  UCHAR sortedCodebook;
  UCHAR endFlag = 0;
  INT intermediateResult;

  SCHAR lengthOfLongestCodeword = pHcr->decInOut.lengthOfLongestCodeword;
  SHORT lengthOfReorderedSpectralData =
      pHcr->decInOut.lengthOfReorderedSpectralData;
  UINT numSortedSection = pHcr->sectionInfo.numSortedSection;
  UCHAR *pSortedCodebook = pHcr->sectionInfo.pSortedCodebook;
  USHORT *pNumSortedCodewordInSection =
      pHcr->sectionInfo.pNumSortedCodewordInSection;
  INT *pLeftStartOfSegment = pHcr->segmentInfo.pLeftStartOfSegment;
  INT *pRightStartOfSegment = pHcr->segmentInfo.pRightStartOfSegment;
  SCHAR *pRemainingBitsInSegment = pHcr->segmentInfo.pRemainingBitsInSegment;
  const UCHAR *pMaxCwLength = aMaxCwLen;

  for (i = numSortedSection; i != 0; i--) {
    sortedCodebook = *pSortedCodebook++;
    segmentWidth =
        fMin((INT)pMaxCwLength[sortedCodebook], (INT)lengthOfLongestCodeword);

    for (j = *pNumSortedCodewordInSection; j != 0; j--) {
      intermediateResult = segmentStart;
      if ((segmentStart + segmentWidth) <= lengthOfReorderedSpectralData) {
        *pLeftStartOfSegment++ = intermediateResult;
        *pRightStartOfSegment++ = intermediateResult + segmentWidth - 1;
        *pRemainingBitsInSegment++ = segmentWidth;
        segmentStart += segmentWidth;
        numSegment += 1;
      }
      else {
        pLeftStartOfSegment--;
        pRightStartOfSegment--;
        pRemainingBitsInSegment--;
        segmentStart = *pLeftStartOfSegment;

        lastSegmentWidth = lengthOfReorderedSpectralData - segmentStart;
        *pRemainingBitsInSegment = lastSegmentWidth;
        *pRightStartOfSegment = segmentStart + lastSegmentWidth - 1;
        endFlag = 1;
        break;
      }
    }
    pNumSortedCodewordInSection++;
    if (endFlag != 0) {
      break;
    }
  }
  pHcr->segmentInfo.numSegment = numSegment;
}

static void HcrExtendedSectionInfo(H_HCR_INFO pHcr) {
  UINT srtSecCnt = 0;
  UINT xSrtScCnt = 0;
  UINT remainNumCwInSortSec;
  UINT inSegmentRemainNumCW;

  UINT numSortedSection = pHcr->sectionInfo.numSortedSection;
  UCHAR *pSortedCodebook = pHcr->sectionInfo.pSortedCodebook;
  USHORT *pNumSortedCodewordInSection =
      pHcr->sectionInfo.pNumSortedCodewordInSection;
  UCHAR *pExtendedSortedCoBo = pHcr->sectionInfo.pExtendedSortedCodebook;
  USHORT *pNumExtSortCwInSect =
      pHcr->sectionInfo.pNumExtendedSortedCodewordInSection;
  UINT numSegment = pHcr->segmentInfo.numSegment;
  UCHAR *pMaxLenOfCbInExtSrtSec = pHcr->sectionInfo.pMaxLenOfCbInExtSrtSec;
  SCHAR lengthOfLongestCodeword = pHcr->decInOut.lengthOfLongestCodeword;
  const UCHAR *pMaxCwLength = aMaxCwLen;

  remainNumCwInSortSec = pNumSortedCodewordInSection[srtSecCnt];
  inSegmentRemainNumCW = numSegment;

  while (srtSecCnt < numSortedSection) {
    if (inSegmentRemainNumCW < remainNumCwInSortSec) {
      pNumExtSortCwInSect[xSrtScCnt] = inSegmentRemainNumCW;
      pExtendedSortedCoBo[xSrtScCnt] = pSortedCodebook[srtSecCnt];

      remainNumCwInSortSec -= inSegmentRemainNumCW;
      inSegmentRemainNumCW = numSegment;
    } else if (inSegmentRemainNumCW == remainNumCwInSortSec) {
      pNumExtSortCwInSect[xSrtScCnt] = inSegmentRemainNumCW;
      pExtendedSortedCoBo[xSrtScCnt] = pSortedCodebook[srtSecCnt];

      srtSecCnt++;
      remainNumCwInSortSec = pNumSortedCodewordInSection[srtSecCnt];
      inSegmentRemainNumCW = numSegment;
    } else {
      pNumExtSortCwInSect[xSrtScCnt] = remainNumCwInSortSec;
      pExtendedSortedCoBo[xSrtScCnt] = pSortedCodebook[srtSecCnt];

      inSegmentRemainNumCW -= remainNumCwInSortSec;
      srtSecCnt++;
      remainNumCwInSortSec = pNumSortedCodewordInSection[srtSecCnt];
    }
    pMaxLenOfCbInExtSrtSec[xSrtScCnt] =
        fMin((INT)pMaxCwLength[pExtendedSortedCoBo[xSrtScCnt]],
             (INT)lengthOfLongestCodeword);

    xSrtScCnt += 1;

    if (xSrtScCnt >= (MAX_SFB_HCR + MAX_HCR_SETS)) {
      return;
    }
  }
  pNumExtSortCwInSect[xSrtScCnt] = 0;
}

static void DeriveNumberOfExtendedSortedSectionsInSets(
    UINT numSegment, USHORT *pNumExtendedSortedCodewordInSection,
    int numExtendedSortedCodewordInSectionIdx,
    USHORT *pNumExtendedSortedSectionsInSets,
    int numExtendedSortedSectionsInSetsIdx) {
  USHORT counter = 0;
  UINT cwSum = 0;
  USHORT *pNumExSortCwInSec = pNumExtendedSortedCodewordInSection;
  USHORT *pNumExSortSecInSets = pNumExtendedSortedSectionsInSets;

  while (pNumExSortCwInSec[numExtendedSortedCodewordInSectionIdx] != 0) {
    cwSum += pNumExSortCwInSec[numExtendedSortedCodewordInSectionIdx];
    numExtendedSortedCodewordInSectionIdx++;
    if (numExtendedSortedCodewordInSectionIdx >= (MAX_SFB_HCR + MAX_HCR_SETS)) {
      return;
    }
    if (cwSum > numSegment) {
      return;
    }
    counter++;
    if (counter > 1024 / 4) {
      return;
    }
    if (cwSum == numSegment) {
      pNumExSortSecInSets[numExtendedSortedSectionsInSetsIdx] = counter;
      numExtendedSortedSectionsInSetsIdx++;
      if (numExtendedSortedSectionsInSetsIdx >= MAX_HCR_SETS) {
        return;
      }
      counter = 0;
      cwSum = 0;
    }
  }
  pNumExSortSecInSets[numExtendedSortedSectionsInSetsIdx] =
      counter;
}

static void DecodePCWs(HANDLE_FDK_BITSTREAM bs, H_HCR_INFO pHcr) {
  UINT i;
  USHORT extSortSec;
  USHORT curExtSortCwInSec;
  UCHAR codebook;
  UCHAR dimension;
  const UINT *pCurrentTree;
  const SCHAR *pQuantValBase;
  const SCHAR *pQuantVal;

  USHORT *pNumExtendedSortedCodewordInSection =
      pHcr->sectionInfo.pNumExtendedSortedCodewordInSection;
  int numExtendedSortedCodewordInSectionIdx =
      pHcr->sectionInfo.numExtendedSortedCodewordInSectionIdx;
  UCHAR *pExtendedSortedCodebook = pHcr->sectionInfo.pExtendedSortedCodebook;
  int extendedSortedCodebookIdx = pHcr->sectionInfo.extendedSortedCodebookIdx;
  USHORT *pNumExtendedSortedSectionsInSets =
      pHcr->sectionInfo.pNumExtendedSortedSectionsInSets;
  int numExtendedSortedSectionsInSetsIdx =
      pHcr->sectionInfo.numExtendedSortedSectionsInSetsIdx;
  FIXP_DBL *pQuantizedSpectralCoefficients =
      SPEC_LONG(pHcr->decInOut.pQuantizedSpectralCoefficientsBase);
  int quantizedSpectralCoefficientsIdx =
      pHcr->decInOut.quantizedSpectralCoefficientsIdx;
  INT *pLeftStartOfSegment = pHcr->segmentInfo.pLeftStartOfSegment;
  SCHAR *pRemainingBitsInSegment = pHcr->segmentInfo.pRemainingBitsInSegment;
  UCHAR *pMaxLenOfCbInExtSrtSec = pHcr->sectionInfo.pMaxLenOfCbInExtSrtSec;
  int maxLenOfCbInExtSrtSecIdx = pHcr->sectionInfo.maxLenOfCbInExtSrtSecIdx;
  UCHAR maxAllowedCwLen;
  int numDecodedBits;
  const UCHAR *pCbDimension = aDimCb;
  const UCHAR *pCbSign = aSignCb;

  FDKmemclear(pQuantizedSpectralCoefficients + quantizedSpectralCoefficientsIdx,
              1024 * sizeof(FIXP_DBL));

  for (extSortSec =
           pNumExtendedSortedSectionsInSets[numExtendedSortedSectionsInSetsIdx];
       extSortSec != 0; extSortSec--) {
    codebook =
        pExtendedSortedCodebook[extendedSortedCodebookIdx];
    extendedSortedCodebookIdx++;
    if (extendedSortedCodebookIdx >= (MAX_SFB_HCR + MAX_HCR_SETS)) {
      return;
    }
    dimension = pCbDimension[codebook];
    pCurrentTree =
        aHuffTable[codebook];
    pQuantValBase =
        aQuantTable[codebook];
    maxAllowedCwLen = pMaxLenOfCbInExtSrtSec[maxLenOfCbInExtSrtSecIdx];
    maxLenOfCbInExtSrtSecIdx++;
    if (maxLenOfCbInExtSrtSecIdx >= (MAX_SFB_HCR + MAX_HCR_SETS)) {
      return;
    }

    if (pCbSign[codebook] ==
        0) {

      for (curExtSortCwInSec = pNumExtendedSortedCodewordInSection
               [numExtendedSortedCodewordInSectionIdx];
           curExtSortCwInSec != 0; curExtSortCwInSec--) {
        numDecodedBits = 0;

        pQuantVal = DecodePCW_Body(
            bs, pHcr->decInOut.bitstreamAnchor, pCurrentTree, pQuantValBase,
            pLeftStartOfSegment, pRemainingBitsInSegment, &numDecodedBits);

        for (i = dimension; i != 0; i--) {
          pQuantizedSpectralCoefficients[quantizedSpectralCoefficientsIdx] =
              (FIXP_DBL)*pQuantVal++;
          quantizedSpectralCoefficientsIdx++;
          if (quantizedSpectralCoefficientsIdx >= 1024) {
            return;
          }
        }

        if (1 == errDetectPcwSegmentation(
                     *pRemainingBitsInSegment - ERROR_PCW_BODY, pHcr, PCW_BODY,
                     pQuantizedSpectralCoefficients +
                         quantizedSpectralCoefficientsIdx - dimension,
                     dimension)) {
          return;
        }
        pLeftStartOfSegment++;
        pRemainingBitsInSegment++;
      }
    } else if ((codebook < 11) && (pCbSign[codebook] ==
                                   1)) {

      for (curExtSortCwInSec = pNumExtendedSortedCodewordInSection
               [numExtendedSortedCodewordInSectionIdx];
           curExtSortCwInSec != 0; curExtSortCwInSec--) {
        int err;
        numDecodedBits = 0;

        pQuantVal = DecodePCW_Body(
            bs, pHcr->decInOut.bitstreamAnchor, pCurrentTree, pQuantValBase,
            pLeftStartOfSegment, pRemainingBitsInSegment, &numDecodedBits);

        err = DecodePCW_Sign(
            bs, pHcr->decInOut.bitstreamAnchor, dimension, pQuantVal,
            pQuantizedSpectralCoefficients, &quantizedSpectralCoefficientsIdx,
            pLeftStartOfSegment, pRemainingBitsInSegment, &numDecodedBits);
        if (err != 0) {
          return;
        }

        if (1 == errDetectPcwSegmentation(
                     *pRemainingBitsInSegment - ERROR_PCW_BODY_SIGN, pHcr,
                     PCW_BODY_SIGN,
                     pQuantizedSpectralCoefficients +
                         quantizedSpectralCoefficientsIdx - dimension,
                     dimension)) {
          return;
        }
        pLeftStartOfSegment++;
        pRemainingBitsInSegment++;
      }
    } else if ((pCbSign[codebook] == 1) &&
               (codebook >= 11)) {

      for (curExtSortCwInSec = pNumExtendedSortedCodewordInSection
               [numExtendedSortedCodewordInSectionIdx];
           curExtSortCwInSec != 0; curExtSortCwInSec--) {
        int err;
        numDecodedBits = 0;

        pQuantVal = DecodePCW_Body(
            bs, pHcr->decInOut.bitstreamAnchor, pCurrentTree, pQuantValBase,
            pLeftStartOfSegment, pRemainingBitsInSegment, &numDecodedBits);

        err = DecodePCW_Sign(
            bs, pHcr->decInOut.bitstreamAnchor, dimension, pQuantVal,
            pQuantizedSpectralCoefficients, &quantizedSpectralCoefficientsIdx,
            pLeftStartOfSegment, pRemainingBitsInSegment, &numDecodedBits);
        if (err != 0) {
          return;
        }

        quantizedSpectralCoefficientsIdx -= DIMENSION_OF_ESCAPE_CODEBOOK;

        if (fixp_abs(pQuantizedSpectralCoefficients
                         [quantizedSpectralCoefficientsIdx]) ==
            (FIXP_DBL)ESCAPE_VALUE) {
          pQuantizedSpectralCoefficients[quantizedSpectralCoefficientsIdx] =
              (FIXP_DBL)DecodeEscapeSequence(
                  bs, pHcr->decInOut.bitstreamAnchor,
                  pQuantizedSpectralCoefficients
                      [quantizedSpectralCoefficientsIdx],
                  pLeftStartOfSegment, pRemainingBitsInSegment, &numDecodedBits,
                  &pHcr->decInOut.errorLog);
        }
        quantizedSpectralCoefficientsIdx++;
        if (quantizedSpectralCoefficientsIdx >= 1024) {
          return;
        }

        if (fixp_abs(pQuantizedSpectralCoefficients
                         [quantizedSpectralCoefficientsIdx]) ==
            (FIXP_DBL)ESCAPE_VALUE) {
          pQuantizedSpectralCoefficients[quantizedSpectralCoefficientsIdx] =
              (FIXP_DBL)DecodeEscapeSequence(
                  bs, pHcr->decInOut.bitstreamAnchor,
                  pQuantizedSpectralCoefficients
                      [quantizedSpectralCoefficientsIdx],
                  pLeftStartOfSegment, pRemainingBitsInSegment, &numDecodedBits,
                  &pHcr->decInOut.errorLog);
        }
        quantizedSpectralCoefficientsIdx++;
        if (quantizedSpectralCoefficientsIdx >= 1024) {
          return;
        }

        if (1 == errDetectPcwSegmentation(
                     *pRemainingBitsInSegment - ERROR_PCW_BODY_SIGN_ESC, pHcr,
                     PCW_BODY_SIGN_ESC,
                     pQuantizedSpectralCoefficients +
                         quantizedSpectralCoefficientsIdx -
                         DIMENSION_OF_ESCAPE_CODEBOOK,
                     DIMENSION_OF_ESCAPE_CODEBOOK)) {
          return;
        }
        pLeftStartOfSegment++;
        pRemainingBitsInSegment++;
      }
    }

    numExtendedSortedCodewordInSectionIdx++;
    if (numExtendedSortedCodewordInSectionIdx >= MAX_SFB_HCR + MAX_HCR_SETS) {
      return;
    }
  }

  numExtendedSortedSectionsInSetsIdx++;
  if (numExtendedSortedSectionsInSetsIdx >= MAX_HCR_SETS) {
    return;
  }

  pHcr->sectionInfo.numExtendedSortedCodewordInSectionIdx =
      numExtendedSortedCodewordInSectionIdx;
  pHcr->sectionInfo.extendedSortedCodebookIdx = extendedSortedCodebookIdx;
  pHcr->sectionInfo.numExtendedSortedSectionsInSetsIdx =
      numExtendedSortedSectionsInSetsIdx;
  pHcr->decInOut.quantizedSpectralCoefficientsIdx =
      quantizedSpectralCoefficientsIdx;
  pHcr->sectionInfo.maxLenOfCbInExtSrtSecIdx = maxLenOfCbInExtSrtSecIdx;
}

static UCHAR errDetectPcwSegmentation(SCHAR remainingBitsInSegment,
                                      H_HCR_INFO pHcr, PCW_TYPE kind,
                                      FIXP_DBL *qsc_base_of_cw,
                                      UCHAR dimension) {
  SCHAR i;
  if (remainingBitsInSegment < 0) {
    switch (kind) {
      case PCW_BODY:
        pHcr->decInOut.errorLog |= SEGMENT_OVERRIDE_ERR_PCW_BODY;
        break;
      case PCW_BODY_SIGN:
        pHcr->decInOut.errorLog |= SEGMENT_OVERRIDE_ERR_PCW_BODY_SIGN;
        break;
      case PCW_BODY_SIGN_ESC:
        pHcr->decInOut.errorLog |= SEGMENT_OVERRIDE_ERR_PCW_BODY_SIGN_ESC;
        break;
    }
    for (i = dimension; i != 0; i--) {
      *qsc_base_of_cw++ = (FIXP_DBL)Q_VALUE_INVALID;
    }
    return 1;
  }
  return 0;
}

static void errDetectWithinSegmentationFinal(H_HCR_INFO pHcr) {
  UCHAR segmentationErrorFlag = 0;
  USHORT i;
  SCHAR *pRemainingBitsInSegment = pHcr->segmentInfo.pRemainingBitsInSegment;
  UINT numSegment = pHcr->segmentInfo.numSegment;

  for (i = numSegment; i != 0; i--) {
    if (*pRemainingBitsInSegment++ != 0) {
      segmentationErrorFlag = 1;
    }
  }
}

void CarryBitToBranchValue(UCHAR carryBit, UINT treeNode, UINT *branchValue,
                           UINT *branchNode) {
  if (carryBit == 0) {
    *branchNode =
        (treeNode & MASK_LEFT) >> LEFT_OFFSET;
  } else {
    *branchNode = treeNode & MASK_RIGHT;
  }

  *branchValue = *branchNode & CLR_BIT_10;
}

static const SCHAR *DecodePCW_Body(HANDLE_FDK_BITSTREAM bs, const INT bsAnchor,
                                   const UINT *pCurrentTree,
                                   const SCHAR *pQuantValBase,
                                   INT *pLeftStartOfSegment,
                                   SCHAR *pRemainingBitsInSegment,
                                   int *pNumDecodedBits) {
  UCHAR carryBit;
  UINT branchNode;
  UINT treeNode;
  UINT branchValue;
  const SCHAR *pQuantVal;

  treeNode = *pCurrentTree;

  while (1) {
    carryBit = HcrGetABitFromBitstream(bs, bsAnchor, pLeftStartOfSegment,
                                       pLeftStartOfSegment,
                                       FROM_LEFT_TO_RIGHT);
    *pRemainingBitsInSegment -= 1;
    *pNumDecodedBits += 1;

    CarryBitToBranchValue(carryBit, treeNode, &branchValue, &branchNode);

    if ((branchNode & TEST_BIT_10) ==
        TEST_BIT_10) {
      break;
    } else {
      treeNode = *(
          pCurrentTree +
          branchValue);
    }
  }

  pQuantVal =
      pQuantValBase + branchValue;

  return pQuantVal;
}

static INT DecodeEscapeSequence(HANDLE_FDK_BITSTREAM bs, const INT bsAnchor,
                                INT quantSpecCoef, INT *pLeftStartOfSegment,
                                SCHAR *pRemainingBitsInSegment,
                                int *pNumDecodedBits, UINT *errorWord) {
  UINT i;
  INT sign;
  UINT escapeOnesCounter = 0;
  UINT carryBit;
  INT escape_word = 0;

  while (1) {
    carryBit = HcrGetABitFromBitstream(bs, bsAnchor, pLeftStartOfSegment,
                                       pLeftStartOfSegment,
                                       FROM_LEFT_TO_RIGHT);
    *pRemainingBitsInSegment -= 1;
    *pNumDecodedBits += 1;
    if (*pRemainingBitsInSegment < 0) {
      return Q_VALUE_INVALID;
    }

    if (carryBit != 0) {
      escapeOnesCounter += 1;
    } else {
      escapeOnesCounter += 4;
      break;
    }
  }

  for (i = escapeOnesCounter; i != 0; i--) {
    carryBit = HcrGetABitFromBitstream(bs, bsAnchor, pLeftStartOfSegment,
                                       pLeftStartOfSegment,
                                       FROM_LEFT_TO_RIGHT);
    *pRemainingBitsInSegment -= 1;
    *pNumDecodedBits += 1;
    if (*pRemainingBitsInSegment < 0) {
      return Q_VALUE_INVALID;
    }

    escape_word <<= 1;
    escape_word = escape_word | carryBit;
  }

  sign = (quantSpecCoef >= 0) ? 1 : -1;

  if (escapeOnesCounter < 13) {
    quantSpecCoef = sign * (((INT)1 << escapeOnesCounter) + escape_word);
  } else {
    *errorWord |= TOO_MANY_PCW_BODY_SIGN_ESC_BITS_DECODED;
    quantSpecCoef = Q_VALUE_INVALID;
  }
  return quantSpecCoef;
}

static int DecodePCW_Sign(HANDLE_FDK_BITSTREAM bs, const INT bsAnchor,
                          UINT codebookDim, const SCHAR *pQuantVal,
                          FIXP_DBL *pQuantSpecCoef, int *quantSpecCoefIdx,
                          INT *pLeftStartOfSegment,
                          SCHAR *pRemainingBitsInSegment,
                          int *pNumDecodedBits) {
  UINT i;
  UINT carryBit;
  INT quantSpecCoef;

  for (i = codebookDim; i != 0; i--) {
    quantSpecCoef = *pQuantVal++;
    if (quantSpecCoef != 0) {
      carryBit = HcrGetABitFromBitstream(bs, bsAnchor, pLeftStartOfSegment,
                                         pLeftStartOfSegment,
                                         FROM_LEFT_TO_RIGHT);
      *pRemainingBitsInSegment -= 1;
      *pNumDecodedBits += 1;
      if (*pRemainingBitsInSegment < 0 || *pNumDecodedBits >= (1024 >> 1)) {
        return -1;
      }

      if (carryBit != 0) {
        pQuantSpecCoef[*quantSpecCoefIdx] = -(FIXP_DBL)quantSpecCoef;
      } else {
        pQuantSpecCoef[*quantSpecCoefIdx] = (FIXP_DBL)quantSpecCoef;
      }
    } else {
      pQuantSpecCoef[*quantSpecCoefIdx] = FL2FXCONST_DBL(0.0f);
    }
    *quantSpecCoefIdx += 1;
    if (*quantSpecCoefIdx >= 1024) {
      return -1;
    }
  }
  return 0;
}

void HcrMuteErroneousLines(H_HCR_INFO hHcr) {
  int c;
  FIXP_DBL *RESTRICT pLong =
      SPEC_LONG(hHcr->decInOut.pQuantizedSpectralCoefficientsBase);

  for (c = 0; c < 1024; c++) {
    if (pLong[c] == (FIXP_DBL)Q_VALUE_INVALID) {
      pLong[c] = FL2FXCONST_DBL(0.0f);
    }
  }
}
