/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Noiseless coder module

*******************************************************************************/

#include "dyn_bits.h"
#include "bit_cnt.h"
#include "psy_const.h"
#include "aacenc_pns.h"
#include "aacEnc_ram.h"
#include "aacEnc_rom.h"

typedef INT (*lookUpTable)[CODE_BOOK_ESC_NDX + 1];

static INT FDKaacEnc_getSideInfoBits(const SECTION_INFO* const huffsection,
                                     const SHORT* const sideInfoTab,
                                     const INT useHCR) {
  INT sideInfoBits;

  if (useHCR &&
      ((huffsection->codeBook == 11) || (huffsection->codeBook >= 16))) {
    sideInfoBits = 5;
  } else {
    sideInfoBits = sideInfoTab[huffsection->sfbCnt];
  }

  return (sideInfoBits);
}

static void FDKaacEnc_buildBitLookUp(
    const SHORT* const quantSpectrum, const INT maxSfb,
    const INT* const sfbOffset, const UINT* const sfbMax,
    INT bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1],
    SECTION_INFO* const huffsection) {
  INT i, sfbWidth;

  for (i = 0; i < maxSfb; i++) {
    huffsection[i].sfbCnt = 1;
    huffsection[i].sfbStart = i;
    huffsection[i].sectionBits = INVALID_BITCOUNT;
    huffsection[i].codeBook = -1;
    sfbWidth = sfbOffset[i + 1] - sfbOffset[i];
    FDKaacEnc_bitCount(quantSpectrum + sfbOffset[i], sfbWidth, sfbMax[i],
                       bitLookUp[i]);
  }
}

static inline INT FDKaacEnc_findBestBook(const INT* const bc, INT* const book,
                                         const INT useVCB11) {
  INT minBits = INVALID_BITCOUNT, j;

  int end = CODE_BOOK_ESC_NDX;

  for (j = 0; j <= end; j++) {
    if (bc[j] < minBits) {
      minBits = bc[j];
      *book = j;
    }
  }
  return (minBits);
}

static inline INT FDKaacEnc_findMinMergeBits(const INT* const bc1,
                                             const INT* const bc2,
                                             const INT useVCB11) {
  INT minBits = INVALID_BITCOUNT, j;

  DWORD_ALIGNED(bc1);
  DWORD_ALIGNED(bc2);

  for (j = 0; j <= CODE_BOOK_ESC_NDX; j++) {
    minBits = fixMin(minBits, bc1[j] + bc2[j]);
  }
  return (minBits);
}

static inline void FDKaacEnc_mergeBitLookUp(INT* const RESTRICT bc1,
                                            const INT* const RESTRICT bc2) {
  int j;

  for (j = 0; j <= CODE_BOOK_ESC_NDX; j++) {
    FDK_ASSERT(INVALID_BITCOUNT == 0x1FFFFFFF);
    bc1[j] = fixMin(bc1[j] + bc2[j], INVALID_BITCOUNT);
  }
}

static inline INT FDKaacEnc_findMaxMerge(const INT* const mergeGainLookUp,
                                         const SECTION_INFO* const huffsection,
                                         const INT maxSfb, INT* const maxNdx) {
  INT i, maxMergeGain = 0;
  int lastMaxNdx = 0;

  for (i = 0; i + huffsection[i].sfbCnt < maxSfb; i += huffsection[i].sfbCnt) {
    if (mergeGainLookUp[i] > maxMergeGain) {
      maxMergeGain = mergeGainLookUp[i];
      lastMaxNdx = i;
    }
  }
  *maxNdx = lastMaxNdx;
  return (maxMergeGain);
}

static inline INT FDKaacEnc_CalcMergeGain(
    const SECTION_INFO* const huffsection,
    const INT bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1],
    const SHORT* const sideInfoTab, const INT ndx1, const INT ndx2,
    const INT useVCB11) {
  INT MergeGain, MergeBits, SplitBits;

  MergeBits =
      sideInfoTab[huffsection[ndx1].sfbCnt + huffsection[ndx2].sfbCnt] +
      FDKaacEnc_findMinMergeBits(bitLookUp[ndx1], bitLookUp[ndx2], useVCB11);
  SplitBits =
      huffsection[ndx1].sectionBits +
      huffsection[ndx2].sectionBits;
  MergeGain = SplitBits - MergeBits;

  if ((huffsection[ndx1].codeBook == CODE_BOOK_PNS_NO) ||
      (huffsection[ndx2].codeBook == CODE_BOOK_PNS_NO) ||
      (huffsection[ndx1].codeBook == CODE_BOOK_IS_OUT_OF_PHASE_NO) ||
      (huffsection[ndx2].codeBook == CODE_BOOK_IS_OUT_OF_PHASE_NO) ||
      (huffsection[ndx1].codeBook == CODE_BOOK_IS_IN_PHASE_NO) ||
      (huffsection[ndx2].codeBook == CODE_BOOK_IS_IN_PHASE_NO)) {
    MergeGain = -1;
  }

  return (MergeGain);
}

static void FDKaacEnc_gmStage0(
    SECTION_INFO* const RESTRICT huffsection,
    const INT bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1], const INT maxSfb,
    const INT* const noiseNrg, const INT* const isBook) {
  INT i;

  for (i = 0; i < maxSfb; i++) {
    if (huffsection[i].sectionBits == INVALID_BITCOUNT) {
      if (noiseNrg[i] != NO_NOISE_PNS) {
        huffsection[i].codeBook = CODE_BOOK_PNS_NO;
        huffsection[i].sectionBits = 0;
      } else if (isBook[i]) {
        huffsection[i].codeBook = isBook[i];
        huffsection[i].sectionBits = 0;
      } else {
        huffsection[i].sectionBits =
            FDKaacEnc_findBestBook(bitLookUp[i], &(huffsection[i].codeBook),
                                   0);
      }
    }
  }
}

static void FDKaacEnc_gmStage1(
    SECTION_INFO* const RESTRICT huffsection,
    INT bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1], const INT maxSfb,
    const SHORT* const sideInfoTab, const INT useVCB11) {
  INT mergeStart = 0, mergeEnd;

  do {
    for (mergeEnd = mergeStart + 1; mergeEnd < maxSfb; mergeEnd++) {
      if (huffsection[mergeStart].codeBook != huffsection[mergeEnd].codeBook)
        break;

      huffsection[mergeStart].sfbCnt++;
      huffsection[mergeStart].sectionBits += huffsection[mergeEnd].sectionBits;

      FDKaacEnc_mergeBitLookUp(bitLookUp[mergeStart], bitLookUp[mergeEnd]);
    }

    huffsection[mergeStart].sectionBits += FDKaacEnc_getSideInfoBits(
        &huffsection[mergeStart], sideInfoTab, useVCB11);
    huffsection[mergeEnd - 1].sfbStart =
        huffsection[mergeStart].sfbStart;

    mergeStart = mergeEnd;

  } while (mergeStart < maxSfb);
}

static inline void FDKaacEnc_gmStage2(
    SECTION_INFO* const RESTRICT huffsection,
    INT* const RESTRICT mergeGainLookUp,
    INT bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1], const INT maxSfb,
    const SHORT* const sideInfoTab, const INT useVCB11) {
  INT i;

  for (i = 0; i + huffsection[i].sfbCnt < maxSfb; i += huffsection[i].sfbCnt) {
    mergeGainLookUp[i] =
        FDKaacEnc_CalcMergeGain(huffsection, bitLookUp, sideInfoTab, i,
                                i + huffsection[i].sfbCnt, useVCB11);
  }

  while (TRUE) {
    INT maxMergeGain, maxNdx, maxNdxNext, maxNdxLast;

    maxMergeGain =
        FDKaacEnc_findMaxMerge(mergeGainLookUp, huffsection, maxSfb, &maxNdx);

    if (maxMergeGain <= 0) break;

    maxNdxNext = maxNdx + huffsection[maxNdx].sfbCnt;

    huffsection[maxNdx].sfbCnt += huffsection[maxNdxNext].sfbCnt;
    huffsection[maxNdx].sectionBits +=
        huffsection[maxNdxNext].sectionBits - maxMergeGain;

    FDKaacEnc_mergeBitLookUp(bitLookUp[maxNdx], bitLookUp[maxNdxNext]);

    if (maxNdx != 0) {
      maxNdxLast = huffsection[maxNdx - 1].sfbStart;
      mergeGainLookUp[maxNdxLast] = FDKaacEnc_CalcMergeGain(
          huffsection, bitLookUp, sideInfoTab, maxNdxLast, maxNdx, useVCB11);
    }
    maxNdxNext = maxNdx + huffsection[maxNdx].sfbCnt;

    huffsection[maxNdxNext - 1].sfbStart = huffsection[maxNdx].sfbStart;

    if (maxNdxNext < maxSfb)
      mergeGainLookUp[maxNdx] = FDKaacEnc_CalcMergeGain(
          huffsection, bitLookUp, sideInfoTab, maxNdx, maxNdxNext, useVCB11);
  }
}

static void FDKaacEnc_noiselessCounter(
    SECTION_DATA* const RESTRICT sectionData, INT mergeGainLookUp[MAX_SFB_LONG],
    INT bitLookUp[MAX_SFB_LONG][CODE_BOOK_ESC_NDX + 1],
    const SHORT* const quantSpectrum, const UINT* const maxValueInSfb,
    const INT* const sfbOffset, const INT blockType, const INT* const noiseNrg,
    const INT* const isBook, const INT useVCB11) {
  INT grpNdx, i;
  const SHORT* sideInfoTab = NULL;
  SECTION_INFO* huffsection;

  switch (blockType) {
    case LONG_WINDOW:
    case START_WINDOW:
    case STOP_WINDOW:
    default:
      sideInfoTab = FDKaacEnc_sideInfoTabLong;
      break;
    case SHORT_WINDOW:
      sideInfoTab = FDKaacEnc_sideInfoTabShort;
      break;
  }

  sectionData->noOfSections = 0;
  sectionData->huffmanBits = 0;
  sectionData->sideInfoBits = 0;

  if (sectionData->maxSfbPerGroup == 0) return;

  for (grpNdx = 0; grpNdx < sectionData->sfbCnt;
       grpNdx += sectionData->sfbPerGroup) {
    huffsection = sectionData->huffsection + sectionData->noOfSections;

    FDKaacEnc_buildBitLookUp(quantSpectrum, sectionData->maxSfbPerGroup,
                             sfbOffset + grpNdx, maxValueInSfb + grpNdx,
                             bitLookUp, huffsection);

    FDKaacEnc_gmStage0(huffsection, bitLookUp, sectionData->maxSfbPerGroup,
                       noiseNrg + grpNdx, isBook + grpNdx);

    FDKaacEnc_gmStage1(huffsection, bitLookUp, sectionData->maxSfbPerGroup,
                       sideInfoTab, useVCB11);

    FDKaacEnc_gmStage2(huffsection, mergeGainLookUp, bitLookUp,
                       sectionData->maxSfbPerGroup, sideInfoTab, useVCB11);

    for (i = 0; i < sectionData->maxSfbPerGroup; i += huffsection[i].sfbCnt) {
      if ((huffsection[i].codeBook == CODE_BOOK_PNS_NO) ||
          (huffsection[i].codeBook == CODE_BOOK_IS_OUT_OF_PHASE_NO) ||
          (huffsection[i].codeBook == CODE_BOOK_IS_IN_PHASE_NO)) {
        huffsection[i].sectionBits = 0;
      } else {

        FDKaacEnc_findBestBook(bitLookUp[i], &(huffsection[i].codeBook),
                               useVCB11);

        sectionData->huffmanBits +=
            huffsection[i].sectionBits -
            FDKaacEnc_getSideInfoBits(&huffsection[i], sideInfoTab, useVCB11);
      }

      huffsection[i].sfbStart += grpNdx;

      sectionData->sideInfoBits +=
          FDKaacEnc_getSideInfoBits(&huffsection[i], sideInfoTab, useVCB11);
      sectionData->huffsection[sectionData->noOfSections++] = huffsection[i];
    }
  }
}

static void FDKaacEnc_scfCount(const INT* const scalefacGain,
                               const UINT* const maxValueInSfb,
                               SECTION_DATA* const RESTRICT sectionData,
                               const INT* const isScale) {
  INT i, j, k, m, n;

  INT lastValScf = 0;
  INT deltaScf = 0;
  INT found = 0;
  INT scfSkipCounter = 0;
  INT lastValIs = 0;

  sectionData->scalefacBits = 0;

  if (scalefacGain == NULL) return;

  sectionData->firstScf = 0;

  for (i = 0; i < sectionData->noOfSections; i++) {
    if (sectionData->huffsection[i].codeBook != CODE_BOOK_ZERO_NO) {
      sectionData->firstScf = sectionData->huffsection[i].sfbStart;
      lastValScf = scalefacGain[sectionData->firstScf];
      break;
    }
  }

  for (i = 0; i < sectionData->noOfSections; i++) {
    if ((sectionData->huffsection[i].codeBook ==
         CODE_BOOK_IS_OUT_OF_PHASE_NO) ||
        (sectionData->huffsection[i].codeBook == CODE_BOOK_IS_IN_PHASE_NO)) {
      for (j = sectionData->huffsection[i].sfbStart;
           j < sectionData->huffsection[i].sfbStart +
                   sectionData->huffsection[i].sfbCnt;
           j++) {
        INT deltaIs = isScale[j] - lastValIs;
        lastValIs = isScale[j];
        sectionData->scalefacBits +=
            FDKaacEnc_bitCountScalefactorDelta(deltaIs);
      }
    }
    else if ((sectionData->huffsection[i].codeBook != CODE_BOOK_ZERO_NO) &&
             (sectionData->huffsection[i].codeBook != CODE_BOOK_PNS_NO)) {
      INT tmp = sectionData->huffsection[i].sfbStart +
                sectionData->huffsection[i].sfbCnt;
      for (j = sectionData->huffsection[i].sfbStart; j < tmp; j++) {
        if (maxValueInSfb[j] == 0) {
          found = 0;
          if (scfSkipCounter == 0) {
            if (j == (tmp - 1))
              found = 0;
            else {
              for (k = (j + 1); k < tmp; k++) {
                if (maxValueInSfb[k] != 0) {
                  found = 1;
                  if ((fixp_abs(scalefacGain[k] - lastValScf)) <=
                      CODE_BOOK_SCF_LAV)
                    deltaScf = 0;
                  else {
                    deltaScf = lastValScf - scalefacGain[j];
                    lastValScf = scalefacGain[j];
                    scfSkipCounter = 0;
                  }
                  break;
                }
                scfSkipCounter++;
              }
            }

            for (m = (i + 1); (m < sectionData->noOfSections) && (found == 0);
                 m++) {
              if ((sectionData->huffsection[m].codeBook != CODE_BOOK_ZERO_NO) &&
                  (sectionData->huffsection[m].codeBook != CODE_BOOK_PNS_NO)) {
                INT end = sectionData->huffsection[m].sfbStart +
                          sectionData->huffsection[m].sfbCnt;
                for (n = sectionData->huffsection[m].sfbStart; n < end; n++) {
                  if (maxValueInSfb[n] != 0) {
                    found = 1;
                    if (fixp_abs(scalefacGain[n] - lastValScf) <=
                        CODE_BOOK_SCF_LAV)
                      deltaScf = 0;
                    else {
                      deltaScf = lastValScf - scalefacGain[j];
                      lastValScf = scalefacGain[j];
                      scfSkipCounter = 0;
                    }
                    break;
                  }
                  scfSkipCounter++;
                }
              }
            }
            if (found == 0) {
              deltaScf = 0;
              scfSkipCounter = 0;
            }
          } else {
            deltaScf = 0;
            scfSkipCounter--;
          }
        } else {
          deltaScf = lastValScf - scalefacGain[j];
          lastValScf = scalefacGain[j];
        }
        sectionData->scalefacBits +=
            FDKaacEnc_bitCountScalefactorDelta(deltaScf);
      }
    }
  }
}

static void FDKaacEnc_noiseCount(SECTION_DATA* const RESTRICT sectionData,
                                 const INT* const noiseNrg) {
  INT noisePCMFlag = TRUE;
  INT lastValPns = 0, deltaPns;
  int i, j;

  sectionData->noiseNrgBits = 0;

  for (i = 0; i < sectionData->noOfSections; i++) {
    if (sectionData->huffsection[i].codeBook == CODE_BOOK_PNS_NO) {
      int sfbStart = sectionData->huffsection[i].sfbStart;
      int sfbEnd = sfbStart + sectionData->huffsection[i].sfbCnt;
      for (j = sfbStart; j < sfbEnd; j++) {
        if (noisePCMFlag) {
          sectionData->noiseNrgBits += PNS_PCM_BITS;
          lastValPns = noiseNrg[j];
          noisePCMFlag = FALSE;
        } else {
          deltaPns = noiseNrg[j] - lastValPns;
          lastValPns = noiseNrg[j];
          sectionData->noiseNrgBits +=
              FDKaacEnc_bitCountScalefactorDelta(deltaPns);
        }
      }
    }
  }
}

INT FDKaacEnc_dynBitCount(BITCNTR_STATE* const hBC,
                          const SHORT* const quantSpectrum,
                          const UINT* const maxValueInSfb,
                          const INT* const scalefac, const INT blockType,
                          const INT sfbCnt, const INT maxSfbPerGroup,
                          const INT sfbPerGroup, const INT* const sfbOffset,
                          SECTION_DATA* const RESTRICT sectionData,
                          const INT* const noiseNrg, const INT* const isBook,
                          const INT* const isScale, const UINT syntaxFlags) {
  sectionData->blockType = blockType;
  sectionData->sfbCnt = sfbCnt;
  sectionData->sfbPerGroup = sfbPerGroup;
  sectionData->noOfGroups = sfbCnt / sfbPerGroup;
  sectionData->maxSfbPerGroup = maxSfbPerGroup;

  FDKaacEnc_noiselessCounter(sectionData, hBC->mergeGainLookUp,
                             (lookUpTable)hBC->bitLookUp, quantSpectrum,
                             maxValueInSfb, sfbOffset, blockType, noiseNrg,
                             isBook, (syntaxFlags & AC_ER_VCB11) ? 1 : 0);

  FDKaacEnc_scfCount(scalefac, maxValueInSfb, sectionData, isScale);

  FDKaacEnc_noiseCount(sectionData, noiseNrg);

  return (sectionData->huffmanBits + sectionData->sideInfoBits +
          sectionData->scalefacBits + sectionData->noiseNrgBits);
}

INT FDKaacEnc_BCNew(BITCNTR_STATE** phBC, UCHAR* dynamic_RAM) {
  BITCNTR_STATE* hBC = GetRam_aacEnc_BitCntrState();

  if (hBC) {
    *phBC = hBC;
    hBC->bitLookUp = GetRam_aacEnc_BitLookUp(0, dynamic_RAM);
    hBC->mergeGainLookUp = GetRam_aacEnc_MergeGainLookUp(0, dynamic_RAM);
    if (hBC->bitLookUp == 0 || hBC->mergeGainLookUp == 0) {
      return 1;
    }
  }
  return (hBC == 0) ? 1 : 0;
}

void FDKaacEnc_BCClose(BITCNTR_STATE** phBC) {
  if (*phBC != NULL) {
    FreeRam_aacEnc_BitCntrState(phBC);
  }
}
