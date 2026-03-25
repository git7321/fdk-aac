/**************************** AAC decoder library ******************************

   Author(s):   Robert Weidner (DSP Solutions)

   Description: HCR Decoder: Common defines and structures; defines for
                switching error-generator, -detector, and -concealment;

*******************************************************************************/

#ifndef AACDEC_HCR_TYPES_H
#define AACDEC_HCR_TYPES_H

#include "FDK_bitstream.h"
#include "overlapadd.h"

#define LINES_PER_UNIT 4

#define MAX_SFB_HCR \
  (((1024 / 8) / LINES_PER_UNIT) * 8)
#define NUMBER_OF_UNIT_GROUPS (LINES_PER_UNIT * 8)
#define LINES_PER_UNIT_GROUP (1024 / NUMBER_OF_UNIT_GROUPS)

#define FROM_LEFT_TO_RIGHT 0
#define FROM_RIGHT_TO_LEFT 1

#define MAX_CB_PAIRS 23
#define MAX_HCR_SETS 14

#define ESCAPE_VALUE 16
#define POSITION_OF_FLAG_A 21
#define POSITION_OF_FLAG_B 20

#define MAX_CB 32

#define MAX_CB_CHECK \
  32

#define NUMBER_OF_BIT_IN_WORD 32

#define THIRTYTWO_LOG_DIV_TWO_LOG 5
#define EIGHT_LOG_DIV_TWO_LOG 3
#define FOUR_LOG_DIV_TWO_LOG 2

#define CPE_TOP_LENGTH 12288
#define SCE_TOP_LENGTH 6144
#define LEN_OF_LONGEST_CW_TOP_LENGTH 49

#define Q_VALUE_INVALID \
  8192
#define HCR_DIRAC 500

#define MASK_LEFT 0xFFF000
#define MASK_RIGHT 0xFFF
#define CLR_BIT_10 0x3FF
#define TEST_BIT_10 0x400

#define LEFT_OFFSET 12

#define ERROR_LORSD 0
#define ERROR_LOLC 0

#define ERROR_PCW_BODY \
  0
#define ERROR_PCW_BODY_SIGN \
  0
#define ERROR_PCW_BODY_SIGN_ESC \
  0

#define ERROR_PCW_BODY_ONLY_TOO_LONG \
  0
#define ERROR_PCW_BODY_SIGN_TOO_LONG \
  0
#define ERROR_PCW_BODY_SIGN_ESC_TOO_LONG \
  0

#define MODULO_DIVISOR_HCR 30
#define HCR_ERROR_CONCEALMENT \
  1

#define SEGMENT_OVERRIDE_ERR_PCW_BODY \
  0x80000000
#define SEGMENT_OVERRIDE_ERR_PCW_BODY_SIGN \
  0x40000000
#define SEGMENT_OVERRIDE_ERR_PCW_BODY_SIGN_ESC \
  0x20000000
#define EXTENDED_SORTED_COUNTER_OVERFLOW \
  0x10000000
#define TOO_MANY_PCW_BODY_BITS_DECODED \
  0x00080000
#define TOO_MANY_PCW_BODY_SIGN_BITS_DECODED \
  0x00040000
#define TOO_MANY_PCW_BODY_SIGN_ESC_BITS_DECODED \
  0x00020000

#define STATE_ERROR_BODY_ONLY \
  0x00008000
#define STATE_ERROR_BODY_SIGN__BODY \
  0x00004000
#define STATE_ERROR_BODY_SIGN__SIGN \
  0x00002000
#define STATE_ERROR_BODY_SIGN_ESC__BODY \
  0x00001000
#define STATE_ERROR_BODY_SIGN_ESC__SIGN \
  0x00000800
#define STATE_ERROR_BODY_SIGN_ESC__ESC_PREFIX \
  0x00000400
#define STATE_ERROR_BODY_SIGN_ESC__ESC_WORD \
  0x00000200
#define HCR_SI_LENGTHS_FAILURE \
  0x00000100
#define NUM_SECT_OUT_OF_RANGE_SHORT_BLOCK \
  0x00000080
#define NUM_SECT_OUT_OF_RANGE_LONG_BLOCK \
  0x00000040
#define LINE_IN_SECT_OUT_OF_RANGE_SHORT_BLOCK \
  0x00000020
#define CB_OUT_OF_RANGE_SHORT_BLOCK \
  0x00000010
#define LINE_IN_SECT_OUT_OF_RANGE_LONG_BLOCK \
  0x00000008
#define CB_OUT_OF_RANGE_LONG_BLOCK \
  0x00000004
#define LAV_VIOLATION \
  0x00000002
#define BIT_IN_SEGMENTATION_ERROR \
  0x00000001

#define HCR_FATAL_PCW_ERROR_MASK 0x100E01FC

typedef enum { PCW_BODY, PCW_BODY_SIGN, PCW_BODY_SIGN_ESC } PCW_TYPE;

typedef struct {
  UINT errorLog;
  SPECTRAL_PTR pQuantizedSpectralCoefficientsBase;
  int quantizedSpectralCoefficientsIdx;
  SHORT lengthOfReorderedSpectralData;
  SHORT numSection;
  SHORT *pNumLineInSect;
  INT bitstreamAnchor;
  SCHAR lengthOfLongestCodeword;
  UCHAR *pCodebook;
} HCR_INPUT_OUTPUT;

typedef struct {
  const UCHAR *pMinOfCbPair;
  const UCHAR *pMaxOfCbPair;
} HCR_CB_PAIRS;

typedef struct {
  const USHORT *pLargestAbsVal;
  const UCHAR *pMaxCwLength;
  const UCHAR *pCbDimension;
  const UCHAR *pCbDimShift;
  const UCHAR *pCbSign;
  const UCHAR *pCbPriority;
} HCR_TABLE_INFO;

typedef struct {
  UINT numSegment;
  UINT pSegmentBitfield[((1024 >> 1) / NUMBER_OF_BIT_IN_WORD + 1)];
  UINT pCodewordBitfield[((1024 >> 1) / NUMBER_OF_BIT_IN_WORD + 1)];
  UINT segmentOffset;
  INT pLeftStartOfSegment[1024 >> 1];
  INT pRightStartOfSegment[1024 >> 1];
  SCHAR pRemainingBitsInSegment[1024 >> 1];
  UCHAR readDirection;
  UCHAR numWordForBitfield;
  USHORT pNumBitValidInLastWord;
} HCR_SEGMENT_INFO;

typedef struct {
  UINT numCodeword;
  UINT numSortedSection;
  USHORT pNumCodewordInSection[MAX_SFB_HCR];
  USHORT pNumSortedCodewordInSection[MAX_SFB_HCR];
  USHORT pNumExtendedSortedCodewordInSection[MAX_SFB_HCR + MAX_HCR_SETS];
  int numExtendedSortedCodewordInSectionIdx;
  USHORT pNumExtendedSortedSectionsInSets[MAX_HCR_SETS];
  int numExtendedSortedSectionsInSetsIdx;
  USHORT pReorderOffset[MAX_SFB_HCR];
  UCHAR pSortedCodebook[MAX_SFB_HCR];

  UCHAR pExtendedSortedCodebook[MAX_SFB_HCR + MAX_HCR_SETS];
  int extendedSortedCodebookIdx;
  UCHAR pMaxLenOfCbInExtSrtSec[MAX_SFB_HCR + MAX_HCR_SETS];
  int maxLenOfCbInExtSrtSecIdx;
  UCHAR pCodebookSwitch[MAX_SFB_HCR];
} HCR_SECTION_INFO;

typedef UINT (*STATEFUNC)(HANDLE_FDK_BITSTREAM, void *);

typedef struct {
  FIXP_DBL
  *pResultBase;
  UINT iNode[1024 >> 2];
  USHORT
  iResultPointer[1024 >> 2];
  UINT pEscapeSequenceInfo[1024 >> 2];
  UINT codewordOffset;
  STATEFUNC pState;
  UCHAR pCodebook[1024 >> 2];
  UCHAR pCntSign[1024 >> 2];
  SCHAR pSta[1024 >> 2];
} HCR_NON_PCW_SIDEINFO;

typedef struct {
  HCR_INPUT_OUTPUT decInOut;
  HCR_SEGMENT_INFO segmentInfo;
  HCR_SECTION_INFO sectionInfo;
  HCR_NON_PCW_SIDEINFO nonPcwSideinfo;
} CErHcrInfo;

typedef CErHcrInfo *H_HCR_INFO;

#endif /* AACDEC_HCR_TYPES_H */
