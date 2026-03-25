/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Huffman Bitcounter & coder

*******************************************************************************/

#include "bit_cnt.h"
#include "aacEnc_ram.h"

#define HI_LTAB(a) (a >> 16)
#define LO_LTAB(a) (a & 0xffff)

static void FDKaacEnc_count1_2_3_4_5_6_7_8_9_10_11(const SHORT *const values,
                                                   const INT width,
                                                   INT *RESTRICT bitCount) {
  INT i;
  INT bc1_2, bc3_4, bc5_6, bc7_8, bc9_10, bc11, sc;
  INT t0, t1, t2, t3;
  bc1_2 = 0;
  bc3_4 = 0;
  bc5_6 = 0;
  bc7_8 = 0;
  bc9_10 = 0;
  bc11 = 0;
  sc = 0;

  DWORD_ALIGNED(values);

  for (i = 0; i < width; i += 4) {
    t0 = values[i + 0];
    t1 = values[i + 1];
    t2 = values[i + 2];
    t3 = values[i + 3];

    bc1_2 += (INT)FDKaacEnc_huff_ltab1_2[t0 + 1][t1 + 1][t2 + 1][t3 + 1];
    bc5_6 += (INT)FDKaacEnc_huff_ltab5_6[t0 + 4][t1 + 4] +
             (INT)FDKaacEnc_huff_ltab5_6[t2 + 4][t3 + 4];

    t0 = fixp_abs(t0);
    sc += (t0 > 0);
    t1 = fixp_abs(t1);
    sc += (t1 > 0);
    t2 = fixp_abs(t2);
    sc += (t2 > 0);
    t3 = fixp_abs(t3);
    sc += (t3 > 0);

    bc3_4 += (INT)FDKaacEnc_huff_ltab3_4[t0][t1][t2][t3];
    bc7_8 += (INT)FDKaacEnc_huff_ltab7_8[t0][t1] +
             (INT)FDKaacEnc_huff_ltab7_8[t2][t3];
    bc9_10 += (INT)FDKaacEnc_huff_ltab9_10[t0][t1] +
              (INT)FDKaacEnc_huff_ltab9_10[t2][t3];
    bc11 +=
        (INT)FDKaacEnc_huff_ltab11[t0][t1] + (INT)FDKaacEnc_huff_ltab11[t2][t3];
  }
  bitCount[1] = HI_LTAB(bc1_2);
  bitCount[2] = LO_LTAB(bc1_2);
  bitCount[3] = HI_LTAB(bc3_4) + sc;
  bitCount[4] = LO_LTAB(bc3_4) + sc;
  bitCount[5] = HI_LTAB(bc5_6);
  bitCount[6] = LO_LTAB(bc5_6);
  bitCount[7] = HI_LTAB(bc7_8) + sc;
  bitCount[8] = LO_LTAB(bc7_8) + sc;
  bitCount[9] = HI_LTAB(bc9_10) + sc;
  bitCount[10] = LO_LTAB(bc9_10) + sc;
  bitCount[11] = bc11 + sc;
}

static void FDKaacEnc_count3_4_5_6_7_8_9_10_11(const SHORT *const values,
                                               const INT width,
                                               INT *RESTRICT bitCount) {
  INT i;
  INT bc3_4, bc5_6, bc7_8, bc9_10, bc11, sc;
  INT t0, t1, t2, t3;

  bc3_4 = 0;
  bc5_6 = 0;
  bc7_8 = 0;
  bc9_10 = 0;
  bc11 = 0;
  sc = 0;

  DWORD_ALIGNED(values);

  for (i = 0; i < width; i += 4) {
    t0 = values[i + 0];
    t1 = values[i + 1];
    t2 = values[i + 2];
    t3 = values[i + 3];

    bc5_6 += (INT)FDKaacEnc_huff_ltab5_6[t0 + 4][t1 + 4] +
             (INT)FDKaacEnc_huff_ltab5_6[t2 + 4][t3 + 4];

    t0 = fixp_abs(t0);
    sc += (t0 > 0);
    t1 = fixp_abs(t1);
    sc += (t1 > 0);
    t2 = fixp_abs(t2);
    sc += (t2 > 0);
    t3 = fixp_abs(t3);
    sc += (t3 > 0);

    bc3_4 += (INT)FDKaacEnc_huff_ltab3_4[t0][t1][t2][t3];
    bc7_8 += (INT)FDKaacEnc_huff_ltab7_8[t0][t1] +
             (INT)FDKaacEnc_huff_ltab7_8[t2][t3];
    bc9_10 += (INT)FDKaacEnc_huff_ltab9_10[t0][t1] +
              (INT)FDKaacEnc_huff_ltab9_10[t2][t3];
    bc11 +=
        (INT)FDKaacEnc_huff_ltab11[t0][t1] + (INT)FDKaacEnc_huff_ltab11[t2][t3];
  }

  bitCount[1] = INVALID_BITCOUNT;
  bitCount[2] = INVALID_BITCOUNT;
  bitCount[3] = HI_LTAB(bc3_4) + sc;
  bitCount[4] = LO_LTAB(bc3_4) + sc;
  bitCount[5] = HI_LTAB(bc5_6);
  bitCount[6] = LO_LTAB(bc5_6);
  bitCount[7] = HI_LTAB(bc7_8) + sc;
  bitCount[8] = LO_LTAB(bc7_8) + sc;
  bitCount[9] = HI_LTAB(bc9_10) + sc;
  bitCount[10] = LO_LTAB(bc9_10) + sc;
  bitCount[11] = bc11 + sc;
}

static void FDKaacEnc_count5_6_7_8_9_10_11(const SHORT *const values,
                                           const INT width,
                                           INT *RESTRICT bitCount) {
  INT i;
  INT bc5_6, bc7_8, bc9_10, bc11, sc;
  INT t0, t1, t2, t3;
  bc5_6 = 0;
  bc7_8 = 0;
  bc9_10 = 0;
  bc11 = 0;
  sc = 0;

  DWORD_ALIGNED(values);

  for (i = 0; i < width; i += 4) {
    t0 = values[i + 0];
    t1 = values[i + 1];
    t2 = values[i + 2];
    t3 = values[i + 3];

    bc5_6 += (INT)FDKaacEnc_huff_ltab5_6[t0 + 4][t1 + 4] +
             (INT)FDKaacEnc_huff_ltab5_6[t2 + 4][t3 + 4];

    t0 = fixp_abs(t0);
    sc += (t0 > 0);
    t1 = fixp_abs(t1);
    sc += (t1 > 0);
    t2 = fixp_abs(t2);
    sc += (t2 > 0);
    t3 = fixp_abs(t3);
    sc += (t3 > 0);

    bc7_8 += (INT)FDKaacEnc_huff_ltab7_8[t0][t1] +
             (INT)FDKaacEnc_huff_ltab7_8[t2][t3];
    bc9_10 += (INT)FDKaacEnc_huff_ltab9_10[t0][t1] +
              (INT)FDKaacEnc_huff_ltab9_10[t2][t3];
    bc11 +=
        (INT)FDKaacEnc_huff_ltab11[t0][t1] + (INT)FDKaacEnc_huff_ltab11[t2][t3];
  }
  bitCount[1] = INVALID_BITCOUNT;
  bitCount[2] = INVALID_BITCOUNT;
  bitCount[3] = INVALID_BITCOUNT;
  bitCount[4] = INVALID_BITCOUNT;
  bitCount[5] = HI_LTAB(bc5_6);
  bitCount[6] = LO_LTAB(bc5_6);
  bitCount[7] = HI_LTAB(bc7_8) + sc;
  bitCount[8] = LO_LTAB(bc7_8) + sc;
  bitCount[9] = HI_LTAB(bc9_10) + sc;
  bitCount[10] = LO_LTAB(bc9_10) + sc;
  bitCount[11] = bc11 + sc;
}

static void FDKaacEnc_count7_8_9_10_11(const SHORT *const values,
                                       const INT width,
                                       INT *RESTRICT bitCount) {
  INT i;
  INT bc7_8, bc9_10, bc11, sc;
  INT t0, t1, t2, t3;

  bc7_8 = 0;
  bc9_10 = 0;
  bc11 = 0;
  sc = 0;

  DWORD_ALIGNED(values);

  for (i = 0; i < width; i += 4) {
    t0 = values[i + 0];
    t1 = values[i + 1];
    t2 = values[i + 2];
    t3 = values[i + 3];

    t0 = fixp_abs(t0);
    sc += (t0 > 0);
    t1 = fixp_abs(t1);
    sc += (t1 > 0);
    t2 = fixp_abs(t2);
    sc += (t2 > 0);
    t3 = fixp_abs(t3);
    sc += (t3 > 0);

    bc7_8 += (INT)FDKaacEnc_huff_ltab7_8[t0][t1] +
             (INT)FDKaacEnc_huff_ltab7_8[t2][t3];
    bc9_10 += (INT)FDKaacEnc_huff_ltab9_10[t0][t1] +
              (INT)FDKaacEnc_huff_ltab9_10[t2][t3];
    bc11 +=
        (INT)FDKaacEnc_huff_ltab11[t0][t1] + (INT)FDKaacEnc_huff_ltab11[t2][t3];
  }

  bitCount[1] = INVALID_BITCOUNT;
  bitCount[2] = INVALID_BITCOUNT;
  bitCount[3] = INVALID_BITCOUNT;
  bitCount[4] = INVALID_BITCOUNT;
  bitCount[5] = INVALID_BITCOUNT;
  bitCount[6] = INVALID_BITCOUNT;
  bitCount[7] = HI_LTAB(bc7_8) + sc;
  bitCount[8] = LO_LTAB(bc7_8) + sc;
  bitCount[9] = HI_LTAB(bc9_10) + sc;
  bitCount[10] = LO_LTAB(bc9_10) + sc;
  bitCount[11] = bc11 + sc;
}

static void FDKaacEnc_count9_10_11(const SHORT *const values, const INT width,
                                   INT *RESTRICT bitCount) {
  INT i;
  INT bc9_10, bc11, sc;
  INT t0, t1, t2, t3;

  bc9_10 = 0;
  bc11 = 0;
  sc = 0;

  DWORD_ALIGNED(values);

  for (i = 0; i < width; i += 4) {
    t0 = values[i + 0];
    t1 = values[i + 1];
    t2 = values[i + 2];
    t3 = values[i + 3];

    t0 = fixp_abs(t0);
    sc += (t0 > 0);
    t1 = fixp_abs(t1);
    sc += (t1 > 0);
    t2 = fixp_abs(t2);
    sc += (t2 > 0);
    t3 = fixp_abs(t3);
    sc += (t3 > 0);

    bc9_10 += (INT)FDKaacEnc_huff_ltab9_10[t0][t1] +
              (INT)FDKaacEnc_huff_ltab9_10[t2][t3];
    bc11 +=
        (INT)FDKaacEnc_huff_ltab11[t0][t1] + (INT)FDKaacEnc_huff_ltab11[t2][t3];
  }

  bitCount[1] = INVALID_BITCOUNT;
  bitCount[2] = INVALID_BITCOUNT;
  bitCount[3] = INVALID_BITCOUNT;
  bitCount[4] = INVALID_BITCOUNT;
  bitCount[5] = INVALID_BITCOUNT;
  bitCount[6] = INVALID_BITCOUNT;
  bitCount[7] = INVALID_BITCOUNT;
  bitCount[8] = INVALID_BITCOUNT;
  bitCount[9] = HI_LTAB(bc9_10) + sc;
  bitCount[10] = LO_LTAB(bc9_10) + sc;
  bitCount[11] = bc11 + sc;
}

static void FDKaacEnc_count11(const SHORT *const values, const INT width,
                              INT *RESTRICT bitCount) {
  INT i;
  INT bc11, sc;
  INT t0, t1, t2, t3;

  bc11 = 0;
  sc = 0;

  DWORD_ALIGNED(values);

  for (i = 0; i < width; i += 4) {
    t0 = values[i + 0];
    t1 = values[i + 1];
    t2 = values[i + 2];
    t3 = values[i + 3];

    t0 = fixp_abs(t0);
    sc += (t0 > 0);
    t1 = fixp_abs(t1);
    sc += (t1 > 0);
    t2 = fixp_abs(t2);
    sc += (t2 > 0);
    t3 = fixp_abs(t3);
    sc += (t3 > 0);

    bc11 +=
        (INT)FDKaacEnc_huff_ltab11[t0][t1] + (INT)FDKaacEnc_huff_ltab11[t2][t3];
  }

  bitCount[1] = INVALID_BITCOUNT;
  bitCount[2] = INVALID_BITCOUNT;
  bitCount[3] = INVALID_BITCOUNT;
  bitCount[4] = INVALID_BITCOUNT;
  bitCount[5] = INVALID_BITCOUNT;
  bitCount[6] = INVALID_BITCOUNT;
  bitCount[7] = INVALID_BITCOUNT;
  bitCount[8] = INVALID_BITCOUNT;
  bitCount[9] = INVALID_BITCOUNT;
  bitCount[10] = INVALID_BITCOUNT;
  bitCount[11] = bc11 + sc;
}

static void FDKaacEnc_countEsc(const SHORT *const values, const INT width,
                               INT *RESTRICT bitCount) {
  INT i;
  INT bc11, ec, sc;
  INT t0, t1, t00, t01;

  bc11 = 0;
  sc = 0;
  ec = 0;
  for (i = 0; i < width; i += 2) {
    t0 = fixp_abs(values[i + 0]);
    t1 = fixp_abs(values[i + 1]);

    sc += (t0 > 0) + (t1 > 0);

    t00 = fixMin(t0, 16);
    t01 = fixMin(t1, 16);
    bc11 += (INT)FDKaacEnc_huff_ltab11[t00][t01];

    if (t0 >= 16) {
      ec += 5;
      while ((t0 >>= 1) >= 16) ec += 2;
    }

    if (t1 >= 16) {
      ec += 5;
      while ((t1 >>= 1) >= 16) ec += 2;
    }
  }

  for (i = 0; i < 11; i++) bitCount[i] = INVALID_BITCOUNT;

  bitCount[11] = bc11 + sc + ec;
}

typedef void (*COUNT_FUNCTION)(const SHORT *const values, const INT width,
                               INT *RESTRICT bitCount);

static const COUNT_FUNCTION countFuncTable[CODE_BOOK_ESC_LAV + 1] = {

    FDKaacEnc_count1_2_3_4_5_6_7_8_9_10_11, /* 0  */
    FDKaacEnc_count1_2_3_4_5_6_7_8_9_10_11, /* 1  */
    FDKaacEnc_count3_4_5_6_7_8_9_10_11,     /* 2  */
    FDKaacEnc_count5_6_7_8_9_10_11,         /* 3  */
    FDKaacEnc_count5_6_7_8_9_10_11,         /* 4  */
    FDKaacEnc_count7_8_9_10_11,             /* 5  */
    FDKaacEnc_count7_8_9_10_11,             /* 6  */
    FDKaacEnc_count7_8_9_10_11,             /* 7  */
    FDKaacEnc_count9_10_11,                 /* 8  */
    FDKaacEnc_count9_10_11,                 /* 9  */
    FDKaacEnc_count9_10_11,                 /* 10 */
    FDKaacEnc_count9_10_11,                 /* 11 */
    FDKaacEnc_count9_10_11,                 /* 12 */
    FDKaacEnc_count11,                      /* 13 */
    FDKaacEnc_count11,                      /* 14 */
    FDKaacEnc_count11,                      /* 15 */
    FDKaacEnc_countEsc                      /* 16 */
};

INT FDKaacEnc_bitCount(const SHORT *const values, const INT width,
                       const INT maxVal, INT *const RESTRICT bitCount) {
  bitCount[0] = (maxVal == 0) ? 0 : INVALID_BITCOUNT;

  countFuncTable[fixMin(maxVal, (INT)CODE_BOOK_ESC_LAV)](values, width,
                                                         bitCount);

  return (0);
}

INT FDKaacEnc_countValues(SHORT *RESTRICT values, INT width, INT codeBook) {
  INT i, t0, t1, t2, t3;
  INT bitCnt = 0;

  switch (codeBook) {
    case CODE_BOOK_ZERO_NO:
      break;

    case CODE_BOOK_1_NO:
      for (i = 0; i < width; i += 4) {
        t0 = values[i + 0];
        t1 = values[i + 1];
        t2 = values[i + 2];
        t3 = values[i + 3];
        bitCnt +=
            HI_LTAB(FDKaacEnc_huff_ltab1_2[t0 + 1][t1 + 1][t2 + 1][t3 + 1]);
      }
      break;

    case CODE_BOOK_2_NO:
      for (i = 0; i < width; i += 4) {
        t0 = values[i + 0];
        t1 = values[i + 1];
        t2 = values[i + 2];
        t3 = values[i + 3];
        bitCnt +=
            LO_LTAB(FDKaacEnc_huff_ltab1_2[t0 + 1][t1 + 1][t2 + 1][t3 + 1]);
      }
      break;

    case CODE_BOOK_3_NO:
      for (i = 0; i < width; i += 4) {
        t0 = fixp_abs(values[i + 0]);
        bitCnt += (t0 > 0);
        t1 = fixp_abs(values[i + 1]);
        bitCnt += (t1 > 0);
        t2 = fixp_abs(values[i + 2]);
        bitCnt += (t2 > 0);
        t3 = fixp_abs(values[i + 3]);
        bitCnt += (t3 > 0);
        bitCnt += HI_LTAB(FDKaacEnc_huff_ltab3_4[t0][t1][t2][t3]);
      }
      break;

    case CODE_BOOK_4_NO:
      for (i = 0; i < width; i += 4) {
        t0 = fixp_abs(values[i + 0]);
        bitCnt += (t0 > 0);
        t1 = fixp_abs(values[i + 1]);
        bitCnt += (t1 > 0);
        t2 = fixp_abs(values[i + 2]);
        bitCnt += (t2 > 0);
        t3 = fixp_abs(values[i + 3]);
        bitCnt += (t3 > 0);
        bitCnt += LO_LTAB(FDKaacEnc_huff_ltab3_4[t0][t1][t2][t3]);
      }
      break;

    case CODE_BOOK_5_NO:
      for (i = 0; i < width; i += 4) {
        t0 = values[i + 0];
        t1 = values[i + 1];
        t2 = values[i + 2];
        t3 = values[i + 3];
        bitCnt += HI_LTAB(FDKaacEnc_huff_ltab5_6[t0 + 4][t1 + 4]) +
                  HI_LTAB(FDKaacEnc_huff_ltab5_6[t2 + 4][t3 + 4]);
      }
      break;

    case CODE_BOOK_6_NO:
      for (i = 0; i < width; i += 4) {
        t0 = values[i + 0];
        t1 = values[i + 1];
        t2 = values[i + 2];
        t3 = values[i + 3];
        bitCnt += LO_LTAB(FDKaacEnc_huff_ltab5_6[t0 + 4][t1 + 4]) +
                  LO_LTAB(FDKaacEnc_huff_ltab5_6[t2 + 4][t3 + 4]);
      }
      break;

    case CODE_BOOK_7_NO:
      for (i = 0; i < width; i += 4) {
        t0 = fixp_abs(values[i + 0]);
        bitCnt += (t0 > 0);
        t1 = fixp_abs(values[i + 1]);
        bitCnt += (t1 > 0);
        t2 = fixp_abs(values[i + 2]);
        bitCnt += (t2 > 0);
        t3 = fixp_abs(values[i + 3]);
        bitCnt += (t3 > 0);
        bitCnt += HI_LTAB(FDKaacEnc_huff_ltab7_8[t0][t1]) +
                  HI_LTAB(FDKaacEnc_huff_ltab7_8[t2][t3]);
      }
      break;

    case CODE_BOOK_8_NO:
      for (i = 0; i < width; i += 4) {
        t0 = fixp_abs(values[i + 0]);
        bitCnt += (t0 > 0);
        t1 = fixp_abs(values[i + 1]);
        bitCnt += (t1 > 0);
        t2 = fixp_abs(values[i + 2]);
        bitCnt += (t2 > 0);
        t3 = fixp_abs(values[i + 3]);
        bitCnt += (t3 > 0);
        bitCnt += LO_LTAB(FDKaacEnc_huff_ltab7_8[t0][t1]) +
                  LO_LTAB(FDKaacEnc_huff_ltab7_8[t2][t3]);
      }
      break;

    case CODE_BOOK_9_NO:
      for (i = 0; i < width; i += 4) {
        t0 = fixp_abs(values[i + 0]);
        bitCnt += (t0 > 0);
        t1 = fixp_abs(values[i + 1]);
        bitCnt += (t1 > 0);
        t2 = fixp_abs(values[i + 2]);
        bitCnt += (t2 > 0);
        t3 = fixp_abs(values[i + 3]);
        bitCnt += (t3 > 0);
        bitCnt += HI_LTAB(FDKaacEnc_huff_ltab9_10[t0][t1]) +
                  HI_LTAB(FDKaacEnc_huff_ltab9_10[t2][t3]);
      }
      break;

    case CODE_BOOK_10_NO:
      for (i = 0; i < width; i += 4) {
        t0 = fixp_abs(values[i + 0]);
        bitCnt += (t0 > 0);
        t1 = fixp_abs(values[i + 1]);
        bitCnt += (t1 > 0);
        t2 = fixp_abs(values[i + 2]);
        bitCnt += (t2 > 0);
        t3 = fixp_abs(values[i + 3]);
        bitCnt += (t3 > 0);
        bitCnt += LO_LTAB(FDKaacEnc_huff_ltab9_10[t0][t1]) +
                  LO_LTAB(FDKaacEnc_huff_ltab9_10[t2][t3]);
      }
      break;

    case CODE_BOOK_ESC_NO:
      for (i = 0; i < width; i += 2) {
        t0 = fixp_abs(values[i + 0]);
        bitCnt += (t0 > 0);
        t1 = fixp_abs(values[i + 1]);
        bitCnt += (t1 > 0);
        bitCnt += (INT)FDKaacEnc_huff_ltab11[fixMin(t0, 16)][fixMin(t1, 16)];
        if (t0 >= 16) {
          bitCnt += 5;
          while ((t0 >>= 1) >= 16) bitCnt += 2;
        }
        if (t1 >= 16) {
          bitCnt += 5;
          while ((t1 >>= 1) >= 16) bitCnt += 2;
        }
      }
      break;

    default:
      break;
  }

  return (bitCnt);
}

INT FDKaacEnc_codeValues(SHORT *RESTRICT values, INT width, INT codeBook,
                         HANDLE_FDK_BITSTREAM hBitstream) {
  INT i, t0, t1, t2, t3, t00, t01;
  INT codeWord, codeLength;
  INT sign, signLength;

  DWORD_ALIGNED(values);

  switch (codeBook) {
    case CODE_BOOK_ZERO_NO:
      break;

    case CODE_BOOK_1_NO:
      for (i = 0; i < width; i += 4) {
        t0 = values[i + 0] + 1;
        t1 = values[i + 1] + 1;
        t2 = values[i + 2] + 1;
        t3 = values[i + 3] + 1;
        codeWord = FDKaacEnc_huff_ctab1[t0][t1][t2][t3];
        codeLength = HI_LTAB(FDKaacEnc_huff_ltab1_2[t0][t1][t2][t3]);
        FDKwriteBits(hBitstream, codeWord, codeLength);
      }
      break;

    case CODE_BOOK_2_NO:
      for (i = 0; i < width; i += 4) {
        t0 = values[i + 0] + 1;
        t1 = values[i + 1] + 1;
        t2 = values[i + 2] + 1;
        t3 = values[i + 3] + 1;
        codeWord = FDKaacEnc_huff_ctab2[t0][t1][t2][t3];
        codeLength = LO_LTAB(FDKaacEnc_huff_ltab1_2[t0][t1][t2][t3]);
        FDKwriteBits(hBitstream, codeWord, codeLength);
      }
      break;

    case CODE_BOOK_3_NO:
      for (i = 0; i < (width >> 2); i++) {
        sign = 0;
        signLength = 0;
        int index[4];
        for (int j = 0; j < 4; j++) {
          int ti = *values++;
          int zero = (ti == 0) ? 0 : 1;
          signLength += zero;
          sign = (sign << zero) + ((UINT)ti >> 31);
          index[j] = fixp_abs(ti);
        }
        codeWord = FDKaacEnc_huff_ctab3[index[0]][index[1]][index[2]][index[3]];
        codeLength = HI_LTAB(
            FDKaacEnc_huff_ltab3_4[index[0]][index[1]][index[2]][index[3]]);
        FDKwriteBits(hBitstream, (codeWord << signLength) | sign,
                     codeLength + signLength);
      }
      break;

    case CODE_BOOK_4_NO:
      for (i = 0; i < width; i += 4) {
        sign = 0;
        signLength = 0;
        int index[4];
        for (int j = 0; j < 4; j++) {
          int ti = *values++;
          int zero = (ti == 0) ? 0 : 1;
          signLength += zero;
          sign = (sign << zero) + ((UINT)ti >> 31);
          index[j] = fixp_abs(ti);
        }
        codeWord = FDKaacEnc_huff_ctab4[index[0]][index[1]][index[2]][index[3]];
        codeLength = LO_LTAB(
            FDKaacEnc_huff_ltab3_4[index[0]][index[1]][index[2]][index[3]]);
        FDKwriteBits(hBitstream, (codeWord << signLength) | sign,
                     codeLength + signLength);
      }
      break;

    case CODE_BOOK_5_NO:
      for (i = 0; i < (width >> 2); i++) {
        t0 = *values++ + 4;
        t1 = *values++ + 4;
        t2 = *values++ + 4;
        t3 = *values++ + 4;
        codeWord = FDKaacEnc_huff_ctab5[t0][t1];
        codeLength =
            HI_LTAB(FDKaacEnc_huff_ltab5_6[t2][t3]);
        codeWord = (codeWord << codeLength) + FDKaacEnc_huff_ctab5[t2][t3];
        codeLength += HI_LTAB(FDKaacEnc_huff_ltab5_6[t0][t1]);
        FDKwriteBits(hBitstream, codeWord, codeLength);
      }
      break;

    case CODE_BOOK_6_NO:
      for (i = 0; i < (width >> 2); i++) {
        t0 = *values++ + 4;
        t1 = *values++ + 4;
        t2 = *values++ + 4;
        t3 = *values++ + 4;
        codeWord = FDKaacEnc_huff_ctab6[t0][t1];
        codeLength =
            LO_LTAB(FDKaacEnc_huff_ltab5_6[t2][t3]);
        codeWord = (codeWord << codeLength) + FDKaacEnc_huff_ctab6[t2][t3];
        codeLength += LO_LTAB(FDKaacEnc_huff_ltab5_6[t0][t1]);
        FDKwriteBits(hBitstream, codeWord, codeLength);
      }
      break;

    case CODE_BOOK_7_NO:
      for (i = 0; i < (width >> 1); i++) {
        t0 = *values++;
        sign = ((UINT)t0 >> 31);
        t0 = fixp_abs(t0);
        signLength = (t0 == 0) ? 0 : 1;
        t1 = *values++;
        INT zero = (t1 == 0) ? 0 : 1;
        signLength += zero;
        sign = (sign << zero) + ((UINT)t1 >> 31);
        t1 = fixp_abs(t1);
        codeWord = FDKaacEnc_huff_ctab7[t0][t1];
        codeLength = HI_LTAB(FDKaacEnc_huff_ltab7_8[t0][t1]);
        FDKwriteBits(hBitstream, (codeWord << signLength) | sign,
                     codeLength + signLength);
      }
      break;

    case CODE_BOOK_8_NO:
      for (i = 0; i < (width >> 1); i++) {
        t0 = *values++;
        sign = ((UINT)t0 >> 31);
        t0 = fixp_abs(t0);
        signLength = (t0 == 0) ? 0 : 1;
        t1 = *values++;
        INT zero = (t1 == 0) ? 0 : 1;
        signLength += zero;
        sign = (sign << zero) + ((UINT)t1 >> 31);
        t1 = fixp_abs(t1);
        codeWord = FDKaacEnc_huff_ctab8[t0][t1];
        codeLength = LO_LTAB(FDKaacEnc_huff_ltab7_8[t0][t1]);
        FDKwriteBits(hBitstream, (codeWord << signLength) | sign,
                     codeLength + signLength);
      }
      break;

    case CODE_BOOK_9_NO:
      for (i = 0; i < (width >> 1); i++) {
        t0 = *values++;
        sign = ((UINT)t0 >> 31);
        t0 = fixp_abs(t0);
        signLength = (t0 == 0) ? 0 : 1;
        t1 = *values++;
        INT zero = (t1 == 0) ? 0 : 1;
        signLength += zero;
        sign = (sign << zero) + ((UINT)t1 >> 31);
        t1 = fixp_abs(t1);
        codeWord = FDKaacEnc_huff_ctab9[t0][t1];
        codeLength = HI_LTAB(FDKaacEnc_huff_ltab9_10[t0][t1]);
        FDKwriteBits(hBitstream, (codeWord << signLength) | sign,
                     codeLength + signLength);
      }
      break;

    case CODE_BOOK_10_NO:
      for (i = 0; i < (width >> 1); i++) {
        t0 = *values++;
        sign = ((UINT)t0 >> 31);
        t0 = fixp_abs(t0);
        signLength = (t0 == 0) ? 0 : 1;
        t1 = *values++;
        INT zero = (t1 == 0) ? 0 : 1;
        signLength += zero;
        sign = (sign << zero) + ((UINT)t1 >> 31);
        t1 = fixp_abs(t1);
        codeWord = FDKaacEnc_huff_ctab10[t0][t1];
        codeLength = LO_LTAB(FDKaacEnc_huff_ltab9_10[t0][t1]);
        FDKwriteBits(hBitstream, (codeWord << signLength) | sign,
                     codeLength + signLength);
      }
      break;

    case CODE_BOOK_ESC_NO:
      for (i = 0; i < (width >> 1); i++) {
        t0 = *values++;
        sign = ((UINT)t0 >> 31);
        t0 = fixp_abs(t0);
        signLength = (t0 == 0) ? 0 : 1;
        t1 = *values++;
        INT zero = (t1 == 0) ? 0 : 1;
        signLength += zero;
        sign = (sign << zero) + ((UINT)t1 >> 31);
        t1 = fixp_abs(t1);

        t00 = fixMin(t0, 16);
        t01 = fixMin(t1, 16);

        codeWord = FDKaacEnc_huff_ctab11[t00][t01];
        codeLength = (INT)FDKaacEnc_huff_ltab11[t00][t01];
        FDKwriteBits(hBitstream, (codeWord << signLength) | sign,
                     codeLength + signLength);
        for (int j = 0; j < 2; j++) {
          if (t0 >= 16) {
            INT n = 4, p = t0;
            for (; (p >>= 1) >= 16;) n++;
            FDKwriteBits(hBitstream,
                         (((1 << (n - 3)) - 2) << n) | (t0 - (1 << n)),
                         n + n - 3);
          }
          t0 = t1;
        }
      }
      break;

    default:
      break;
  }
  return (0);
}

INT FDKaacEnc_codeScalefactorDelta(INT delta, HANDLE_FDK_BITSTREAM hBitstream) {
  INT codeWord, codeLength;

  if (fixp_abs(delta) > CODE_BOOK_SCF_LAV) return (1);

  codeWord = FDKaacEnc_huff_ctabscf[delta + CODE_BOOK_SCF_LAV];
  codeLength = (INT)FDKaacEnc_huff_ltabscf[delta + CODE_BOOK_SCF_LAV];
  FDKwriteBits(hBitstream, codeWord, codeLength);
  return (0);
}
