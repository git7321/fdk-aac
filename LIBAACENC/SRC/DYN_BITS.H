/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Noiseless coder module

*******************************************************************************/

#ifndef DYN_BITS_H
#define DYN_BITS_H

#include "common_fix.h"

#include "psy_const.h"
#include "aacenc_tns.h"

#define MAX_SECTIONS MAX_GROUPED_SFB
#define SECT_ESC_VAL_LONG 31
#define SECT_ESC_VAL_SHORT 7
#define CODE_BOOK_BITS 4
#define SECT_BITS_LONG 5
#define SECT_BITS_SHORT 3
#define PNS_PCM_BITS 9

typedef struct {
  INT codeBook;
  INT sfbStart;
  INT sfbCnt;
  INT sectionBits;
} SECTION_INFO;

typedef struct {
  INT blockType;
  INT noOfGroups;
  INT sfbCnt;
  INT maxSfbPerGroup;
  INT sfbPerGroup;
  INT noOfSections;
  SECTION_INFO huffsection[MAX_SECTIONS];
  INT sideInfoBits;
  INT huffmanBits;
  INT scalefacBits;
  INT noiseNrgBits;
  INT firstScf;
} SECTION_DATA;

struct BITCNTR_STATE {
  INT* bitLookUp;
  INT* mergeGainLookUp;
};

INT FDKaacEnc_BCNew(BITCNTR_STATE** phBC, UCHAR* dynamic_RAM);

void FDKaacEnc_BCClose(BITCNTR_STATE** phBC);

INT FDKaacEnc_dynBitCount(BITCNTR_STATE* const hBC,
                          const SHORT* const quantSpectrum,
                          const UINT* const maxValueInSfb,
                          const INT* const scalefac, const INT blockType,
                          const INT sfbCnt, const INT maxSfbPerGroup,
                          const INT sfbPerGroup, const INT* const sfbOffset,
                          SECTION_DATA* const RESTRICT sectionData,
                          const INT* const noiseNrg, const INT* const isBook,
                          const INT* const isScale, const UINT syntaxFlags);

#endif
