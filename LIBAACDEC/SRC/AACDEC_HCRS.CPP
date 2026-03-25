/**************************** AAC decoder library ******************************

   Author(s):   Robert Weidner (DSP Solutions)

   Description: HCR Decoder: Prepare decoding of non-PCWs, segmentation- and
                bitfield-handling, HCR-Statemachine

*******************************************************************************/

#include "aacdec_hcrs.h"
#include "aacdec_hcr.h"
#include "aacdec_hcr_bit.h"
#include "aac_rom.h"
#include "aac_ram.h"

static UINT InitSegmentBitfield(UINT *pNumSegment,
                                SCHAR *pRemainingBitsInSegment,
                                UINT *pSegmentBitfield,
                                UCHAR *pNumWordForBitfield,
                                USHORT *pNumBitValidInLastWord);

static void InitNonPCWSideInformationForCurrentSet(H_HCR_INFO pHcr);

static INT ModuloValue(INT input, INT bufferlength);

static void ClearBitFromBitfield(STATEFUNC *ptrState, UINT offset,
                                 UINT *pBitfield);

void DecodeNonPCWs(HANDLE_FDK_BITSTREAM bs, H_HCR_INFO pHcr) {
  UINT numValidSegment;
  INT segmentOffset;
  INT codewordOffsetBase;
  INT codewordOffset;
  UINT trial;

  UINT *pNumSegment;
  SCHAR *pRemainingBitsInSegment;
  UINT *pSegmentBitfield;
  UCHAR *pNumWordForBitfield;
  USHORT *pNumBitValidInLastWord;
  UINT *pCodewordBitfield;
  INT bitfieldWord;
  INT bitInWord;
  UINT tempWord;
  UINT interMediateWord;
  INT tempBit;
  INT carry;

  UINT numCodeword;
  UCHAR numSet;
  UCHAR currentSet;
  UINT codewordInSet;
  UINT remainingCodewordsInSet;
  SCHAR *pSta;
  UINT ret;

  pNumSegment = &(pHcr->segmentInfo.numSegment);
  pRemainingBitsInSegment = pHcr->segmentInfo.pRemainingBitsInSegment;
  pSegmentBitfield = pHcr->segmentInfo.pSegmentBitfield;
  pNumWordForBitfield = &(pHcr->segmentInfo.numWordForBitfield);
  pNumBitValidInLastWord = &(pHcr->segmentInfo.pNumBitValidInLastWord);
  pSta = pHcr->nonPcwSideinfo.pSta;

  numValidSegment = InitSegmentBitfield(pNumSegment, pRemainingBitsInSegment,
                                        pSegmentBitfield, pNumWordForBitfield,
                                        pNumBitValidInLastWord);

  if (numValidSegment != 0) {
    numCodeword = pHcr->sectionInfo.numCodeword;
    numSet = ((numCodeword - 1) / *pNumSegment) + 1;

    pHcr->segmentInfo.readDirection = FROM_RIGHT_TO_LEFT;

    numSet = fMin(numSet, (UCHAR)MAX_HCR_SETS);
    for (currentSet = 1; currentSet < numSet; currentSet++) {

      numCodeword -=
          *pNumSegment;
      if (numCodeword < *pNumSegment) {
        codewordInSet = numCodeword;
      } else {
        codewordInSet = *pNumSegment;
      }

      tempWord = 0xFFFFFFFF;
      pCodewordBitfield = pHcr->segmentInfo.pCodewordBitfield;

      for (bitfieldWord = *pNumWordForBitfield; bitfieldWord != 0;
           bitfieldWord--) {
        if (codewordInSet > NUMBER_OF_BIT_IN_WORD) {
          *pCodewordBitfield++ = tempWord;
          codewordInSet -= NUMBER_OF_BIT_IN_WORD;
        } else {
          for (remainingCodewordsInSet = codewordInSet;
               remainingCodewordsInSet < NUMBER_OF_BIT_IN_WORD;
               remainingCodewordsInSet++) {
            tempWord =
                tempWord &
                ~(1
                  << (NUMBER_OF_BIT_IN_WORD - 1 -
                      remainingCodewordsInSet));
          }
          *pCodewordBitfield++ = tempWord;
          tempWord = 0x00000000;
        }
      }
      pCodewordBitfield = pHcr->segmentInfo.pCodewordBitfield;

      InitNonPCWSideInformationForCurrentSet(pHcr);

      codewordOffsetBase = 0;
      for (trial = *pNumSegment; trial > 0; trial--) {
        segmentOffset = 0;
        pHcr->segmentInfo.segmentOffset =
            segmentOffset;
        codewordOffset = codewordOffsetBase;
        pHcr->nonPcwSideinfo.codewordOffset =
            codewordOffset;

        for (bitfieldWord = 0; bitfieldWord < *pNumWordForBitfield;
             bitfieldWord++) {
          tempWord =
              pSegmentBitfield[bitfieldWord] & pCodewordBitfield[bitfieldWord];

          if (tempWord != 0) {
            for (bitInWord = NUMBER_OF_BIT_IN_WORD; bitInWord > 0;
                 bitInWord--) {
              interMediateWord = ((UINT)1 << (bitInWord - 1));
              if ((tempWord & interMediateWord) == interMediateWord) {
                pHcr->nonPcwSideinfo.pState =
                    aStateConstant2State[pSta[codewordOffset]];

                while (pHcr->nonPcwSideinfo.pState) {
                  ret = ((STATEFUNC)pHcr->nonPcwSideinfo.pState)(bs, pHcr);
                  if (ret != 0) {
                    return;
                  }
                }
              }

              segmentOffset += 1;
              pHcr->segmentInfo.segmentOffset = segmentOffset;
              codewordOffset += 1;
              codewordOffset =
                  ModuloValue(codewordOffset,
                              *pNumSegment);
              pHcr->nonPcwSideinfo.codewordOffset = codewordOffset;
            }
          } else {
            segmentOffset +=
                NUMBER_OF_BIT_IN_WORD;
            pHcr->segmentInfo.segmentOffset = segmentOffset;
            codewordOffset +=
                NUMBER_OF_BIT_IN_WORD;
            codewordOffset = ModuloValue(
                codewordOffset,
                *pNumSegment);
            pHcr->nonPcwSideinfo.codewordOffset = codewordOffset;
          }
        }

        codewordOffsetBase -= 1;
        codewordOffsetBase =
            ModuloValue(codewordOffsetBase, *pNumSegment);

        tempBit = pCodewordBitfield[*pNumWordForBitfield - 1] &
                  (1 << (NUMBER_OF_BIT_IN_WORD - *pNumBitValidInLastWord));
        tempBit = tempBit >> (NUMBER_OF_BIT_IN_WORD - *pNumBitValidInLastWord);

        pCodewordBitfield[*pNumWordForBitfield - 1] =
            pCodewordBitfield[*pNumWordForBitfield - 1] &
            ~(1 << (NUMBER_OF_BIT_IN_WORD - *pNumBitValidInLastWord));

        pCodewordBitfield[*pNumWordForBitfield - 1] =
            pCodewordBitfield[*pNumWordForBitfield - 1] >> 1;

        for (bitfieldWord = *pNumWordForBitfield - 2; bitfieldWord > -1;
             bitfieldWord--) {
          carry = pCodewordBitfield[bitfieldWord] & 1;

          pCodewordBitfield[bitfieldWord + 1] =
              pCodewordBitfield[bitfieldWord + 1] |
              (carry << (NUMBER_OF_BIT_IN_WORD - 1));

          pCodewordBitfield[bitfieldWord] =
              pCodewordBitfield[bitfieldWord] >> 1;
        }

        pCodewordBitfield[0] =
            pCodewordBitfield[0] | (tempBit << (NUMBER_OF_BIT_IN_WORD - 1));

      }

      pHcr->segmentInfo.readDirection =
          ToggleReadDirection(pHcr->segmentInfo.readDirection);
    }
  }
}

static UINT InitSegmentBitfield(UINT *pNumSegment,
                                SCHAR *pRemainingBitsInSegment,
                                UINT *pSegmentBitfield,
                                UCHAR *pNumWordForBitfield,
                                USHORT *pNumBitValidInLastWord) {
  SHORT i;
  USHORT r;
  UCHAR bitfieldWord;
  UINT tempWord;
  USHORT numValidSegment;

  *pNumWordForBitfield =
      (*pNumSegment == 0)
          ? 0
          : ((*pNumSegment - 1) >> THIRTYTWO_LOG_DIV_TWO_LOG) + 1;

  numValidSegment = 0;
  *pNumBitValidInLastWord = *pNumSegment;

  for (bitfieldWord = 0; bitfieldWord < *pNumWordForBitfield - 1;
       bitfieldWord++) {
    tempWord = 0xFFFFFFFF;
    r = bitfieldWord << THIRTYTWO_LOG_DIV_TWO_LOG;
    for (i = 0; i < NUMBER_OF_BIT_IN_WORD; i++) {
      if (pRemainingBitsInSegment[r + i] == 0) {
        tempWord = tempWord & ~(1 << (NUMBER_OF_BIT_IN_WORD - 1 -
                                      i));
      } else {
        numValidSegment += 1;
      }
    }
    pSegmentBitfield[bitfieldWord] = tempWord;
    *pNumBitValidInLastWord -= NUMBER_OF_BIT_IN_WORD;
  }

  tempWord = 0xFFFFFFFF;
  for (i = 0; i < (NUMBER_OF_BIT_IN_WORD - *pNumBitValidInLastWord); i++) {
    tempWord = tempWord & ~(1 << i);
  }

  r = bitfieldWord << THIRTYTWO_LOG_DIV_TWO_LOG;
  for (i = 0; i < *pNumBitValidInLastWord; i++) {
    if (pRemainingBitsInSegment[r + i] == 0) {
      tempWord = tempWord & ~(1 << (NUMBER_OF_BIT_IN_WORD - 1 -
                                    i));
    } else {
      numValidSegment += 1;
    }
  }
  pSegmentBitfield[bitfieldWord] = tempWord;

  return numValidSegment;
}

static void InitNonPCWSideInformationForCurrentSet(H_HCR_INFO pHcr) {
  USHORT i, k;
  UCHAR codebookDim;
  UINT startNode;

  UCHAR *pCodebook = pHcr->nonPcwSideinfo.pCodebook;
  UINT *iNode = pHcr->nonPcwSideinfo.iNode;
  UCHAR *pCntSign = pHcr->nonPcwSideinfo.pCntSign;
  USHORT *iResultPointer = pHcr->nonPcwSideinfo.iResultPointer;
  UINT *pEscapeSequenceInfo = pHcr->nonPcwSideinfo.pEscapeSequenceInfo;
  SCHAR *pSta = pHcr->nonPcwSideinfo.pSta;
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
  int quantizedSpectralCoefficientsIdx =
      pHcr->decInOut.quantizedSpectralCoefficientsIdx;
  const UCHAR *pCbDimension = aDimCb;
  int iterationCounter = 0;

  for (i = pNumExtendedSortedSectionsInSets[numExtendedSortedSectionsInSetsIdx];
       i != 0; i--) {
    codebookDim =
        pCbDimension[pExtendedSortedCodebook[extendedSortedCodebookIdx]];
    startNode = *aHuffTable[pExtendedSortedCodebook[extendedSortedCodebookIdx]];

    for (k = pNumExtendedSortedCodewordInSection
             [numExtendedSortedCodewordInSectionIdx];
         k != 0; k--) {
      iterationCounter++;
      if (iterationCounter > (1024 >> 2)) {
        return;
      }
      *pSta++ = aCodebook2StartInt
          [pExtendedSortedCodebook[extendedSortedCodebookIdx]];
      *pCodebook++ = pExtendedSortedCodebook[extendedSortedCodebookIdx];
      *iNode++ = startNode;
      *pCntSign++ = 0;
      *iResultPointer++ = quantizedSpectralCoefficientsIdx;
      *pEscapeSequenceInfo++ = 0;
      quantizedSpectralCoefficientsIdx +=
          codebookDim;
      if (quantizedSpectralCoefficientsIdx >= 1024) {
        return;
      }
    }
    numExtendedSortedCodewordInSectionIdx++;
    extendedSortedCodebookIdx++;
    if (numExtendedSortedCodewordInSectionIdx >= (MAX_SFB_HCR + MAX_HCR_SETS) ||
        extendedSortedCodebookIdx >= (MAX_SFB_HCR + MAX_HCR_SETS)) {
      return;
    }
  }
  numExtendedSortedSectionsInSetsIdx++;
  if (numExtendedSortedCodewordInSectionIdx >= (MAX_SFB_HCR + MAX_HCR_SETS)) {
    return;
  }

  pHcr->sectionInfo.numExtendedSortedCodewordInSectionIdx =
      numExtendedSortedCodewordInSectionIdx;
  pHcr->sectionInfo.extendedSortedCodebookIdx = extendedSortedCodebookIdx;
  pHcr->sectionInfo.numExtendedSortedSectionsInSetsIdx =
      numExtendedSortedSectionsInSetsIdx;
  pHcr->sectionInfo.numExtendedSortedCodewordInSectionIdx =
      numExtendedSortedCodewordInSectionIdx;
  pHcr->decInOut.quantizedSpectralCoefficientsIdx =
      quantizedSpectralCoefficientsIdx;
}

static INT ModuloValue(INT input, INT bufferlength) {
  if (input > (bufferlength - 1)) {
    return (input - bufferlength);
  }
  if (input < 0) {
    return (input + bufferlength);
  }
  return input;
}

static void ClearBitFromBitfield(STATEFUNC *ptrState, UINT offset,
                                 UINT *pBitfield) {
  UINT numBitfieldWord;
  UINT numBitfieldBit;

  numBitfieldWord = offset >> THIRTYTWO_LOG_DIV_TWO_LOG;
  numBitfieldBit = offset - (numBitfieldWord
                             << THIRTYTWO_LOG_DIV_TWO_LOG);

  pBitfield[numBitfieldWord] =
      pBitfield[numBitfieldWord] &
      ~(1 << (NUMBER_OF_BIT_IN_WORD - 1 - numBitfieldBit));

  *ptrState = NULL;
}

UINT Hcr_State_BODY_ONLY(HANDLE_FDK_BITSTREAM bs, void *ptr) {
  H_HCR_INFO pHcr = (H_HCR_INFO)ptr;
  UINT *pSegmentBitfield;
  UINT *pCodewordBitfield;
  UINT segmentOffset;
  FIXP_DBL *pResultBase;
  UINT *iNode;
  USHORT *iResultPointer;
  UINT codewordOffset;
  UINT branchNode;
  UINT branchValue;
  UINT iQSC;
  UINT treeNode;
  UCHAR carryBit;
  INT *pLeftStartOfSegment;
  INT *pRightStartOfSegment;
  SCHAR *pRemainingBitsInSegment;
  UCHAR readDirection;
  UCHAR *pCodebook;
  UCHAR dimCntr;
  const UINT *pCurrentTree;
  const UCHAR *pCbDimension;
  const SCHAR *pQuantVal;
  const SCHAR *pQuantValBase;

  pRemainingBitsInSegment = pHcr->segmentInfo.pRemainingBitsInSegment;
  pLeftStartOfSegment = pHcr->segmentInfo.pLeftStartOfSegment;
  pRightStartOfSegment = pHcr->segmentInfo.pRightStartOfSegment;
  readDirection = pHcr->segmentInfo.readDirection;
  pSegmentBitfield = pHcr->segmentInfo.pSegmentBitfield;
  pCodewordBitfield = pHcr->segmentInfo.pCodewordBitfield;
  segmentOffset = pHcr->segmentInfo.segmentOffset;

  pCodebook = pHcr->nonPcwSideinfo.pCodebook;
  iNode = pHcr->nonPcwSideinfo.iNode;
  pResultBase = pHcr->nonPcwSideinfo.pResultBase;
  iResultPointer = pHcr->nonPcwSideinfo.iResultPointer;
  codewordOffset = pHcr->nonPcwSideinfo.codewordOffset;

  pCbDimension = aDimCb;

  treeNode = iNode[codewordOffset];
  pCurrentTree = aHuffTable[pCodebook[codewordOffset]];

  for (; pRemainingBitsInSegment[segmentOffset] > 0;
       pRemainingBitsInSegment[segmentOffset] -= 1) {
    carryBit = HcrGetABitFromBitstream(
        bs, pHcr->decInOut.bitstreamAnchor, &pLeftStartOfSegment[segmentOffset],
        &pRightStartOfSegment[segmentOffset], readDirection);

    CarryBitToBranchValue(carryBit,
                          treeNode, &branchValue, &branchNode);

    if ((branchNode & TEST_BIT_10) ==
        TEST_BIT_10) {
      pQuantValBase = aQuantTable[pCodebook[codewordOffset]];
      pQuantVal = pQuantValBase + branchValue;

      iQSC = iResultPointer[codewordOffset];

      for (dimCntr = pCbDimension[pCodebook[codewordOffset]]; dimCntr != 0;
           dimCntr--) {
        pResultBase[iQSC++] =
            (FIXP_DBL)*pQuantVal++;
      }

      ClearBitFromBitfield(&(pHcr->nonPcwSideinfo.pState), segmentOffset,
                           pCodewordBitfield);
      pRemainingBitsInSegment[segmentOffset] -= 1;
      break;
    } else {
      treeNode = *(
          pCurrentTree +
          branchValue);
    }
  }
  iNode[codewordOffset] = treeNode;

  if (pRemainingBitsInSegment[segmentOffset] <= 0) {
    ClearBitFromBitfield(&(pHcr->nonPcwSideinfo.pState), segmentOffset,
                         pSegmentBitfield);

    if (pRemainingBitsInSegment[segmentOffset] < 0) {
      return BODY_ONLY;
    }
  }

  return STOP_THIS_STATE;
}

UINT Hcr_State_BODY_SIGN__BODY(HANDLE_FDK_BITSTREAM bs, void *ptr) {
  H_HCR_INFO pHcr = (H_HCR_INFO)ptr;
  SCHAR *pRemainingBitsInSegment;
  INT *pLeftStartOfSegment;
  INT *pRightStartOfSegment;
  UCHAR readDirection;
  UINT *pSegmentBitfield;
  UINT *pCodewordBitfield;
  UINT segmentOffset;

  UCHAR *pCodebook;
  UINT *iNode;
  UCHAR *pCntSign;
  FIXP_DBL *pResultBase;
  USHORT *iResultPointer;
  UINT codewordOffset;

  UINT iQSC;
  UINT cntSign;
  UCHAR dimCntr;
  UCHAR carryBit;
  SCHAR *pSta;
  UINT treeNode;
  UINT branchValue;
  UINT branchNode;
  const UCHAR *pCbDimension;
  const UINT *pCurrentTree;
  const SCHAR *pQuantValBase;
  const SCHAR *pQuantVal;

  pRemainingBitsInSegment = pHcr->segmentInfo.pRemainingBitsInSegment;
  pLeftStartOfSegment = pHcr->segmentInfo.pLeftStartOfSegment;
  pRightStartOfSegment = pHcr->segmentInfo.pRightStartOfSegment;
  readDirection = pHcr->segmentInfo.readDirection;
  pSegmentBitfield = pHcr->segmentInfo.pSegmentBitfield;
  pCodewordBitfield = pHcr->segmentInfo.pCodewordBitfield;
  segmentOffset = pHcr->segmentInfo.segmentOffset;

  pCodebook = pHcr->nonPcwSideinfo.pCodebook;
  iNode = pHcr->nonPcwSideinfo.iNode;
  pCntSign = pHcr->nonPcwSideinfo.pCntSign;
  pResultBase = pHcr->nonPcwSideinfo.pResultBase;
  iResultPointer = pHcr->nonPcwSideinfo.iResultPointer;
  codewordOffset = pHcr->nonPcwSideinfo.codewordOffset;
  pSta = pHcr->nonPcwSideinfo.pSta;

  pCbDimension = aDimCb;

  treeNode = iNode[codewordOffset];
  pCurrentTree = aHuffTable[pCodebook[codewordOffset]];

  for (; pRemainingBitsInSegment[segmentOffset] > 0;
       pRemainingBitsInSegment[segmentOffset] -= 1) {
    carryBit = HcrGetABitFromBitstream(
        bs, pHcr->decInOut.bitstreamAnchor, &pLeftStartOfSegment[segmentOffset],
        &pRightStartOfSegment[segmentOffset], readDirection);

    CarryBitToBranchValue(carryBit,
                          treeNode, &branchValue, &branchNode);

    if ((branchNode & TEST_BIT_10) ==
        TEST_BIT_10) {
      pQuantValBase = aQuantTable[pCodebook[codewordOffset]];
      pQuantVal = pQuantValBase + branchValue;

      iQSC = iResultPointer[codewordOffset];

      cntSign = 0;
      for (dimCntr = pCbDimension[pCodebook[codewordOffset]]; dimCntr != 0;
           dimCntr--) {
        pResultBase[iQSC++] =
            (FIXP_DBL)*pQuantVal;
        if (*pQuantVal++ != 0) {
          cntSign += 1;
        }
      }

      if (cntSign == 0) {
        ClearBitFromBitfield(
            &(pHcr->nonPcwSideinfo.pState), segmentOffset,
            pCodewordBitfield);
      } else {
        pCntSign[codewordOffset] = cntSign;
        pSta[codewordOffset] = BODY_SIGN__SIGN;
        pHcr->nonPcwSideinfo.pState =
            aStateConstant2State[pSta[codewordOffset]];
      }
      pRemainingBitsInSegment[segmentOffset] -= 1;
      break;
    } else {
      treeNode = *(
          pCurrentTree +
          branchValue);
    }
  }
  iNode[codewordOffset] = treeNode;

  if (pRemainingBitsInSegment[segmentOffset] <= 0) {
    ClearBitFromBitfield(&(pHcr->nonPcwSideinfo.pState), segmentOffset,
                         pSegmentBitfield);

    if (pRemainingBitsInSegment[segmentOffset] < 0) {
      return BODY_SIGN__BODY;
    }
  }

  return STOP_THIS_STATE;
}

UINT Hcr_State_BODY_SIGN__SIGN(HANDLE_FDK_BITSTREAM bs, void *ptr) {
  H_HCR_INFO pHcr = (H_HCR_INFO)ptr;
  SCHAR *pRemainingBitsInSegment;
  INT *pLeftStartOfSegment;
  INT *pRightStartOfSegment;
  UCHAR readDirection;
  UINT *pSegmentBitfield;
  UINT *pCodewordBitfield;
  UINT segmentOffset;

  UCHAR *pCntSign;
  FIXP_DBL *pResultBase;
  USHORT *iResultPointer;
  UINT codewordOffset;

  UCHAR carryBit;
  UINT iQSC;
  UCHAR cntSign;

  pRemainingBitsInSegment = pHcr->segmentInfo.pRemainingBitsInSegment;
  pLeftStartOfSegment = pHcr->segmentInfo.pLeftStartOfSegment;
  pRightStartOfSegment = pHcr->segmentInfo.pRightStartOfSegment;
  readDirection = pHcr->segmentInfo.readDirection;
  pSegmentBitfield = pHcr->segmentInfo.pSegmentBitfield;
  pCodewordBitfield = pHcr->segmentInfo.pCodewordBitfield;
  segmentOffset = pHcr->segmentInfo.segmentOffset;

  pCntSign = pHcr->nonPcwSideinfo.pCntSign;
  pResultBase = pHcr->nonPcwSideinfo.pResultBase;
  iResultPointer = pHcr->nonPcwSideinfo.iResultPointer;
  codewordOffset = pHcr->nonPcwSideinfo.codewordOffset;

  iQSC = iResultPointer[codewordOffset];
  cntSign = pCntSign[codewordOffset];

  for (; pRemainingBitsInSegment[segmentOffset] > 0;
       pRemainingBitsInSegment[segmentOffset] -= 1) {
    carryBit = HcrGetABitFromBitstream(
        bs, pHcr->decInOut.bitstreamAnchor, &pLeftStartOfSegment[segmentOffset],
        &pRightStartOfSegment[segmentOffset], readDirection);
    cntSign -=
        1;

    while (pResultBase[iQSC] == (FIXP_DBL)0) {
      if (++iQSC >= 1024) {
        return BODY_SIGN__SIGN;
      }
    }

    if (carryBit != 0) {
      pResultBase[iQSC] = -pResultBase[iQSC];
    }

    iQSC++;

    if (cntSign == 0) {
      ClearBitFromBitfield(&(pHcr->nonPcwSideinfo.pState), segmentOffset,
                           pCodewordBitfield);
      pRemainingBitsInSegment[segmentOffset] -= 1;
      break;
    }
  }
  pCntSign[codewordOffset] = cntSign;
  iResultPointer[codewordOffset] = iQSC;

  if (pRemainingBitsInSegment[segmentOffset] <= 0) {
    ClearBitFromBitfield(&(pHcr->nonPcwSideinfo.pState), segmentOffset,
                         pSegmentBitfield);

    if (pRemainingBitsInSegment[segmentOffset] < 0) {
      return BODY_SIGN__SIGN;
    }
  }

  return STOP_THIS_STATE;
}

UINT Hcr_State_BODY_SIGN_ESC__BODY(HANDLE_FDK_BITSTREAM bs, void *ptr) {
  H_HCR_INFO pHcr = (H_HCR_INFO)ptr;
  SCHAR *pRemainingBitsInSegment;
  INT *pLeftStartOfSegment;
  INT *pRightStartOfSegment;
  UCHAR readDirection;
  UINT *pSegmentBitfield;
  UINT *pCodewordBitfield;
  UINT segmentOffset;

  UINT *iNode;
  UCHAR *pCntSign;
  FIXP_DBL *pResultBase;
  USHORT *iResultPointer;
  UINT codewordOffset;

  UCHAR carryBit;
  UINT iQSC;
  UINT cntSign;
  UINT dimCntr;
  UINT treeNode;
  SCHAR *pSta;
  UINT branchNode;
  UINT branchValue;
  const UINT *pCurrentTree;
  const SCHAR *pQuantValBase;
  const SCHAR *pQuantVal;

  pRemainingBitsInSegment = pHcr->segmentInfo.pRemainingBitsInSegment;
  pLeftStartOfSegment = pHcr->segmentInfo.pLeftStartOfSegment;
  pRightStartOfSegment = pHcr->segmentInfo.pRightStartOfSegment;
  readDirection = pHcr->segmentInfo.readDirection;
  pSegmentBitfield = pHcr->segmentInfo.pSegmentBitfield;
  pCodewordBitfield = pHcr->segmentInfo.pCodewordBitfield;
  segmentOffset = pHcr->segmentInfo.segmentOffset;

  iNode = pHcr->nonPcwSideinfo.iNode;
  pCntSign = pHcr->nonPcwSideinfo.pCntSign;
  pResultBase = pHcr->nonPcwSideinfo.pResultBase;
  iResultPointer = pHcr->nonPcwSideinfo.iResultPointer;
  codewordOffset = pHcr->nonPcwSideinfo.codewordOffset;
  pSta = pHcr->nonPcwSideinfo.pSta;

  treeNode = iNode[codewordOffset];
  pCurrentTree = aHuffTable[ESCAPE_CODEBOOK];

  for (; pRemainingBitsInSegment[segmentOffset] > 0;
       pRemainingBitsInSegment[segmentOffset] -= 1) {
    carryBit = HcrGetABitFromBitstream(
        bs, pHcr->decInOut.bitstreamAnchor, &pLeftStartOfSegment[segmentOffset],
        &pRightStartOfSegment[segmentOffset], readDirection);

    CarryBitToBranchValue(carryBit, treeNode, &branchValue, &branchNode);

    if ((branchNode & TEST_BIT_10) ==
        TEST_BIT_10) {

      pQuantValBase = aQuantTable[ESCAPE_CODEBOOK];
      pQuantVal = pQuantValBase + branchValue;

      iNode[codewordOffset] = iResultPointer[codewordOffset];

      iQSC = iResultPointer[codewordOffset];

      cntSign = 0;

      for (dimCntr = DIMENSION_OF_ESCAPE_CODEBOOK; dimCntr != 0; dimCntr--) {
        pResultBase[iQSC++] =
            (FIXP_DBL)*pQuantVal;
        if (*pQuantVal++ != 0) {
          cntSign += 1;
        }
      }

      if (cntSign == 0) {
        ClearBitFromBitfield(
            &(pHcr->nonPcwSideinfo.pState), segmentOffset,
            pCodewordBitfield);
      } else {
        pCntSign[codewordOffset] = cntSign;
        pSta[codewordOffset] = BODY_SIGN_ESC__SIGN;
        pHcr->nonPcwSideinfo.pState =
            aStateConstant2State[pSta[codewordOffset]];
      }
      pRemainingBitsInSegment[segmentOffset] -= 1;
      break;
    } else {
      treeNode = *(pCurrentTree + branchValue);
      iNode[codewordOffset] = treeNode;
    }
  }

  if (pRemainingBitsInSegment[segmentOffset] <= 0) {
    ClearBitFromBitfield(&(pHcr->nonPcwSideinfo.pState), segmentOffset,
                         pSegmentBitfield);

    if (pRemainingBitsInSegment[segmentOffset] < 0) {
      return BODY_SIGN_ESC__BODY;
    }
  }

  return STOP_THIS_STATE;
}

UINT Hcr_State_BODY_SIGN_ESC__SIGN(HANDLE_FDK_BITSTREAM bs, void *ptr) {
  H_HCR_INFO pHcr = (H_HCR_INFO)ptr;
  SCHAR *pRemainingBitsInSegment;
  INT *pLeftStartOfSegment;
  INT *pRightStartOfSegment;
  UCHAR readDirection;
  UINT *pSegmentBitfield;
  UINT *pCodewordBitfield;
  UINT segmentOffset;

  UINT *iNode;
  UCHAR *pCntSign;
  FIXP_DBL *pResultBase;
  USHORT *iResultPointer;
  UINT *pEscapeSequenceInfo;
  UINT codewordOffset;

  UINT iQSC;
  UCHAR cntSign;
  UINT flagA;
  UINT flagB;
  UINT flags;
  UCHAR carryBit;
  SCHAR *pSta;

  pRemainingBitsInSegment = pHcr->segmentInfo.pRemainingBitsInSegment;
  pLeftStartOfSegment = pHcr->segmentInfo.pLeftStartOfSegment;
  pRightStartOfSegment = pHcr->segmentInfo.pRightStartOfSegment;
  readDirection = pHcr->segmentInfo.readDirection;
  pSegmentBitfield = pHcr->segmentInfo.pSegmentBitfield;
  pCodewordBitfield = pHcr->segmentInfo.pCodewordBitfield;
  segmentOffset = pHcr->segmentInfo.segmentOffset;

  iNode = pHcr->nonPcwSideinfo.iNode;
  pCntSign = pHcr->nonPcwSideinfo.pCntSign;
  pResultBase = pHcr->nonPcwSideinfo.pResultBase;
  iResultPointer = pHcr->nonPcwSideinfo.iResultPointer;
  pEscapeSequenceInfo = pHcr->nonPcwSideinfo.pEscapeSequenceInfo;
  codewordOffset = pHcr->nonPcwSideinfo.codewordOffset;
  pSta = pHcr->nonPcwSideinfo.pSta;

  iQSC = iResultPointer[codewordOffset];
  cntSign = pCntSign[codewordOffset];

  for (; pRemainingBitsInSegment[segmentOffset] > 0;
       pRemainingBitsInSegment[segmentOffset] -= 1) {
    carryBit = HcrGetABitFromBitstream(
        bs, pHcr->decInOut.bitstreamAnchor, &pLeftStartOfSegment[segmentOffset],
        &pRightStartOfSegment[segmentOffset], readDirection);

    cntSign -= 1;
    pCntSign[codewordOffset] = cntSign;

    while (pResultBase[iQSC] == (FIXP_DBL)0) {
      if (++iQSC >= 1024) {
        return BODY_SIGN_ESC__SIGN;
      }
    }
    iResultPointer[codewordOffset] = iQSC;

    if (carryBit != 0) {
      pResultBase[iQSC] = -pResultBase[iQSC];
    }
    iQSC++;
    iResultPointer[codewordOffset] = iQSC;

    if (cntSign == 0) {
      pRemainingBitsInSegment[segmentOffset] -= 1;

      iQSC = iNode[codewordOffset];

      flagA = 0;
      if (fixp_abs(pResultBase[iQSC++]) == (FIXP_DBL)ESCAPE_VALUE) {
        flagA = 1;
      }

      flagB = 0;
      if (fixp_abs(pResultBase[iQSC]) == (FIXP_DBL)ESCAPE_VALUE) {
        flagB = 1;
      }

      if (!flagA && !flagB) {
        ClearBitFromBitfield(
            &(pHcr->nonPcwSideinfo.pState), segmentOffset,
            pCodewordBitfield);
      } else {
        flags = flagA << POSITION_OF_FLAG_A;
        flags |= (flagB << POSITION_OF_FLAG_B);
        pEscapeSequenceInfo[codewordOffset] = flags;

        pSta[codewordOffset] = BODY_SIGN_ESC__ESC_PREFIX;
        pHcr->nonPcwSideinfo.pState =
            aStateConstant2State[pSta[codewordOffset]];

        iResultPointer[codewordOffset] = iNode[codewordOffset];

        if (!flagA && flagB) {
          iQSC = iResultPointer[codewordOffset];
          iQSC++;
          iResultPointer[codewordOffset] = iQSC;
        }
      }
      break;

    }
  }

  if (pRemainingBitsInSegment[segmentOffset] <= 0) {
    ClearBitFromBitfield(&(pHcr->nonPcwSideinfo.pState), segmentOffset,
                         pSegmentBitfield);

    if (pRemainingBitsInSegment[segmentOffset] < 0) {
      return BODY_SIGN_ESC__SIGN;
    }
  }
  return STOP_THIS_STATE;
}

UINT Hcr_State_BODY_SIGN_ESC__ESC_PREFIX(HANDLE_FDK_BITSTREAM bs, void *ptr) {
  H_HCR_INFO pHcr = (H_HCR_INFO)ptr;
  SCHAR *pRemainingBitsInSegment;
  INT *pLeftStartOfSegment;
  INT *pRightStartOfSegment;
  UCHAR readDirection;
  UINT *pSegmentBitfield;
  UINT segmentOffset;
  UINT *pEscapeSequenceInfo;
  UINT codewordOffset;
  UCHAR carryBit;
  UINT escapePrefixUp;
  SCHAR *pSta;

  pRemainingBitsInSegment = pHcr->segmentInfo.pRemainingBitsInSegment;
  pLeftStartOfSegment = pHcr->segmentInfo.pLeftStartOfSegment;
  pRightStartOfSegment = pHcr->segmentInfo.pRightStartOfSegment;
  readDirection = pHcr->segmentInfo.readDirection;
  pSegmentBitfield = pHcr->segmentInfo.pSegmentBitfield;
  segmentOffset = pHcr->segmentInfo.segmentOffset;
  pEscapeSequenceInfo = pHcr->nonPcwSideinfo.pEscapeSequenceInfo;
  codewordOffset = pHcr->nonPcwSideinfo.codewordOffset;
  pSta = pHcr->nonPcwSideinfo.pSta;

  escapePrefixUp =
      (pEscapeSequenceInfo[codewordOffset] & MASK_ESCAPE_PREFIX_UP) >>
      LSB_ESCAPE_PREFIX_UP;

  for (; pRemainingBitsInSegment[segmentOffset] > 0;
       pRemainingBitsInSegment[segmentOffset] -= 1) {
    carryBit = HcrGetABitFromBitstream(
        bs, pHcr->decInOut.bitstreamAnchor, &pLeftStartOfSegment[segmentOffset],
        &pRightStartOfSegment[segmentOffset], readDirection);

    if (carryBit == 1) {
      escapePrefixUp += 1;
      if (escapePrefixUp > 8) {
        return BODY_SIGN_ESC__ESC_PREFIX;
      }

      pEscapeSequenceInfo[codewordOffset] &=
          ~MASK_ESCAPE_PREFIX_UP;
      escapePrefixUp <<= LSB_ESCAPE_PREFIX_UP;
      pEscapeSequenceInfo[codewordOffset] |=
          escapePrefixUp;
      escapePrefixUp >>= LSB_ESCAPE_PREFIX_UP;
    } else {
      pRemainingBitsInSegment[segmentOffset] -= 1;
      escapePrefixUp +=
          4;

      pEscapeSequenceInfo[codewordOffset] &=
          ~MASK_ESCAPE_PREFIX_UP;
      escapePrefixUp <<= LSB_ESCAPE_PREFIX_UP;
      pEscapeSequenceInfo[codewordOffset] |=
          escapePrefixUp;
      escapePrefixUp >>= LSB_ESCAPE_PREFIX_UP;

      pEscapeSequenceInfo[codewordOffset] &=
          ~MASK_ESCAPE_PREFIX_DOWN;
      escapePrefixUp <<= LSB_ESCAPE_PREFIX_DOWN;
      pEscapeSequenceInfo[codewordOffset] |=
          escapePrefixUp;

      pSta[codewordOffset] = BODY_SIGN_ESC__ESC_WORD;
      pHcr->nonPcwSideinfo.pState =
          aStateConstant2State[pSta[codewordOffset]];
      break;
    }
  }

  if (pRemainingBitsInSegment[segmentOffset] <= 0) {
    ClearBitFromBitfield(&(pHcr->nonPcwSideinfo.pState), segmentOffset,
                         pSegmentBitfield);

    if (pRemainingBitsInSegment[segmentOffset] < 0) {
      return BODY_SIGN_ESC__ESC_PREFIX;
    }
  }

  return STOP_THIS_STATE;
}

UINT Hcr_State_BODY_SIGN_ESC__ESC_WORD(HANDLE_FDK_BITSTREAM bs, void *ptr) {
  H_HCR_INFO pHcr = (H_HCR_INFO)ptr;
  SCHAR *pRemainingBitsInSegment;
  INT *pLeftStartOfSegment;
  INT *pRightStartOfSegment;
  UCHAR readDirection;
  UINT *pSegmentBitfield;
  UINT *pCodewordBitfield;
  UINT segmentOffset;

  FIXP_DBL *pResultBase;
  USHORT *iResultPointer;
  UINT *pEscapeSequenceInfo;
  UINT codewordOffset;

  UINT escapeWord;
  UINT escapePrefixDown;
  UINT escapePrefixUp;
  UCHAR carryBit;
  UINT iQSC;
  INT sign;
  UINT flagA;
  UINT flagB;
  SCHAR *pSta;

  pRemainingBitsInSegment = pHcr->segmentInfo.pRemainingBitsInSegment;
  pLeftStartOfSegment = pHcr->segmentInfo.pLeftStartOfSegment;
  pRightStartOfSegment = pHcr->segmentInfo.pRightStartOfSegment;
  readDirection = pHcr->segmentInfo.readDirection;
  pSegmentBitfield = pHcr->segmentInfo.pSegmentBitfield;
  pCodewordBitfield = pHcr->segmentInfo.pCodewordBitfield;
  segmentOffset = pHcr->segmentInfo.segmentOffset;

  pResultBase = pHcr->nonPcwSideinfo.pResultBase;
  iResultPointer = pHcr->nonPcwSideinfo.iResultPointer;
  pEscapeSequenceInfo = pHcr->nonPcwSideinfo.pEscapeSequenceInfo;
  codewordOffset = pHcr->nonPcwSideinfo.codewordOffset;
  pSta = pHcr->nonPcwSideinfo.pSta;

  escapeWord = pEscapeSequenceInfo[codewordOffset] & MASK_ESCAPE_WORD;
  escapePrefixDown =
      (pEscapeSequenceInfo[codewordOffset] & MASK_ESCAPE_PREFIX_DOWN) >>
      LSB_ESCAPE_PREFIX_DOWN;

  for (; pRemainingBitsInSegment[segmentOffset] > 0;
       pRemainingBitsInSegment[segmentOffset] -= 1) {
    carryBit = HcrGetABitFromBitstream(
        bs, pHcr->decInOut.bitstreamAnchor, &pLeftStartOfSegment[segmentOffset],
        &pRightStartOfSegment[segmentOffset], readDirection);

    escapeWord <<=
        1;
    escapeWord = escapeWord | carryBit;

    escapePrefixDown -= 1;

    pEscapeSequenceInfo[codewordOffset] &=
        ~MASK_ESCAPE_PREFIX_DOWN;
    escapePrefixDown <<= LSB_ESCAPE_PREFIX_DOWN;
    pEscapeSequenceInfo[codewordOffset] |=
        escapePrefixDown;
    escapePrefixDown >>= LSB_ESCAPE_PREFIX_DOWN;

    pEscapeSequenceInfo[codewordOffset] &=
        ~MASK_ESCAPE_WORD;
    pEscapeSequenceInfo[codewordOffset] |=
        escapeWord;

    if (escapePrefixDown == 0) {
      pRemainingBitsInSegment[segmentOffset] -= 1;

      iQSC = iResultPointer[codewordOffset];
      sign = (pResultBase[iQSC] >= (FIXP_DBL)0)
                 ? 1
                 : -1;

      escapePrefixUp =
          (pEscapeSequenceInfo[codewordOffset] & MASK_ESCAPE_PREFIX_UP) >>
          LSB_ESCAPE_PREFIX_UP;

      pResultBase[iQSC] =
          (FIXP_DBL)(sign * (((INT)1 << escapePrefixUp) + (INT)escapeWord));

      flagA = pEscapeSequenceInfo[codewordOffset] & MASK_FLAG_A;
      flagB = pEscapeSequenceInfo[codewordOffset] & MASK_FLAG_B;

      pEscapeSequenceInfo[codewordOffset] = 0;

      if (flagA != 0) {
        pEscapeSequenceInfo[codewordOffset] &= ~MASK_FLAG_A;

        if (flagB == 0) {
          ClearBitFromBitfield(&(pHcr->nonPcwSideinfo.pState), segmentOffset,
                               pCodewordBitfield);
        } else {
          iQSC++;
          iResultPointer[codewordOffset] = iQSC;

          pSta[codewordOffset] = BODY_SIGN_ESC__ESC_PREFIX;
          pHcr->nonPcwSideinfo.pState =
              aStateConstant2State[pSta[codewordOffset]];
        }
      } else {
        ClearBitFromBitfield(
            &(pHcr->nonPcwSideinfo.pState), segmentOffset,
            pCodewordBitfield);
      }
      break;
    }
  }

  if (pRemainingBitsInSegment[segmentOffset] <= 0) {
    ClearBitFromBitfield(&(pHcr->nonPcwSideinfo.pState), segmentOffset,
                         pSegmentBitfield);

    if (pRemainingBitsInSegment[segmentOffset] < 0) {
      return BODY_SIGN_ESC__ESC_WORD;
    }
  }

  return STOP_THIS_STATE;
}
