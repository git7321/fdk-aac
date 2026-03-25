/*********************** MPEG surround decoder library *************************

   Author(s):   Matthias Hildenbrand

   Description: USAC MPS212 Transient Steering Decorrelator (TSD)

*******************************************************************************/

#include "sac_tsd.h"

#define TSD_START_BAND (7)
#define SIZE_S (4)
#define SIZE_C (5)

RAM_ALIGN
LNK_SECTION_CONSTDATA
static const UCHAR nBitsTsdCW_32slots[32] = {
    5,  9,  13, 16, 18, 20, 22, 24, 25, 26, 27, 28, 29, 29, 30, 30,
    30, 29, 29, 28, 27, 26, 25, 24, 22, 20, 18, 16, 13, 9,  5,  0};

RAM_ALIGN
LNK_SECTION_CONSTDATA
static const UCHAR nBitsTsdCW_64slots[64] = {
    6,  11, 16, 20, 23, 27, 30, 33, 35, 38, 40, 42, 44, 46, 48, 49,
    51, 52, 53, 55, 56, 57, 58, 58, 59, 60, 60, 60, 61, 61, 61, 61,
    61, 61, 61, 60, 60, 60, 59, 58, 58, 57, 56, 55, 53, 52, 51, 49,
    48, 46, 44, 42, 40, 38, 35, 33, 30, 27, 23, 20, 16, 11, 6,  0};

RAM_ALIGN
LNK_SECTION_CONSTDATA
static const FIXP_DPK phiTsd[8] = {
    {{(FIXP_DBL)0x7fffffff, (FIXP_DBL)0x00000000}},
    {{(FIXP_DBL)0x5a82799a, (FIXP_DBL)0x5a82799a}},
    {{(FIXP_DBL)0x00000000, (FIXP_DBL)0x7fffffff}},
    {{(FIXP_DBL)0xa57d8666, (FIXP_DBL)0x5a82799a}},
    {{(FIXP_DBL)0x80000000, (FIXP_DBL)0x00000000}},
    {{(FIXP_DBL)0xa57d8666, (FIXP_DBL)0xa57d8666}},
    {{(FIXP_DBL)0x00000000, (FIXP_DBL)0x80000000}},
    {{(FIXP_DBL)0x5a82799a, (FIXP_DBL)0xa57d8666}}};

static void longmult1(USHORT a[], USHORT b, USHORT d[], int len) {
  int k;
  ULONG tmp;
  ULONG b0 = (ULONG)b;

  tmp = ((ULONG)a[0]) * b0;
  d[0] = (USHORT)tmp;

  for (k = 1; k < len; k++) {
    tmp = (tmp >> 16) + ((ULONG)a[k]) * b0;
    d[k] = (USHORT)tmp;
  }
}

static void longdiv(USHORT b[], USHORT a, USHORT d[], USHORT *pr, int len) {
  ULONG r;
  ULONG tmp;
  int k;

  FDK_ASSERT(a != 0);
  r = 0;

  for (k = len - 1; k >= 0; k--) {
    tmp = ((ULONG)b[k]) + (r << 16);

    if (tmp) {
      d[k] = (USHORT)(tmp / a);
      r = tmp - d[k] * a;
    } else {
      d[k] = 0;
    }
  }
  *pr = (USHORT)r;
}

static void longsub(USHORT a[], USHORT b[], int lena, int lenb) {
  int h;
  LONG carry = 0;

  FDK_ASSERT(lena >= lenb);
  for (h = 0; h < lenb; h++) {
    carry += ((LONG)a[h]) - ((LONG)b[h]);
    a[h] = (USHORT)carry;
    carry = carry >> 16;
  }

  for (; h < lena; h++) {
    carry = ((LONG)a[h]) + carry;
    a[h] = (USHORT)carry;
    carry = carry >> 16;
  }

  FDK_ASSERT(carry ==
             0);
  return;
}

static int longcompare(USHORT a[], USHORT b[], int len) {
  int i;

  for (i = len - 1; i > 0; i--) {
    if (a[i] != b[i]) break;
  }
  return (a[i] >= b[i]) ? 1 : 0;
}

FDK_INLINE int isTrSlot(const TSD_DATA *pTsdData, const int ts) {
  return (pTsdData->bsTsdTrPhaseData[ts] >= 0);
}

int TsdRead(HANDLE_FDK_BITSTREAM hBs, const int numSlots, TSD_DATA *pTsdData) {
  int nBitsTrSlots = 0;
  int bsTsdNumTrSlots;
  const UCHAR *nBitsTsdCW_tab = NULL;

  switch (numSlots) {
    case 32:
      nBitsTrSlots = 4;
      nBitsTsdCW_tab = nBitsTsdCW_32slots;
      break;
    case 64:
      nBitsTrSlots = 5;
      nBitsTsdCW_tab = nBitsTsdCW_64slots;
      break;
    default:
      return 1;
  }

  pTsdData->bsTsdEnable = FDKreadBit(hBs);
  if (!pTsdData->bsTsdEnable) {
    return 0;
  }

  pTsdData->numSlots = numSlots;

  bsTsdNumTrSlots = FDKreadBits(hBs, nBitsTrSlots);

  {
    int nBitsTsdCW = (int)nBitsTsdCW_tab[bsTsdNumTrSlots];
    SCHAR *phaseData = pTsdData->bsTsdTrPhaseData;
    int p = bsTsdNumTrSlots + 1;
    int k, h;
    USHORT s[SIZE_S] = {0};
    USHORT c[SIZE_C] = {0};
    USHORT r[1];

    for (k = 0; k < numSlots; k++) {
      phaseData[k] = -1;
    }

    for (h = (SIZE_S - 1); h >= 0; h--) {
      if (nBitsTsdCW > h * 16) {
        s[h] = (USHORT)FDKreadBits(hBs, nBitsTsdCW - h * 16);
        nBitsTsdCW = h * 16;
      }
    }

    k = numSlots - 1;
    c[0] = k - p + 1;
    for (h = 2; h <= p; h++) {
      longmult1(c, (k - p + h), c, 5);
      longdiv(c, h, c, r, 5);
      FDK_ASSERT(*r == 0);
    }

    for (; k >= 0; k--) {
      if (p > k) {
        for (; k >= 0; k--) {
          phaseData[k] = 1;
        }
        break;
      }
      if (longcompare(s, c, 4)) {
        longsub(s, c, 4, 4);
        phaseData[k] = 1;
        if (p == 1) {
          break;
        }
        longmult1(c, p, c, 5);
        p--;
      } else {
        longmult1(c, (k - p), c, 5);
      }
      longdiv(c, k, c, r, 5);
      FDK_ASSERT(*r == 0);
    }

    for (k = 0; k < numSlots; k++) {
      if (phaseData[k] == 1) {
        phaseData[k] = FDKreadBits(hBs, 3);
      }
    }
  }

  return 0;
}

void TsdGenerateNonTr(const int numHybridBands, const TSD_DATA *pTsdData,
                      const int ts, FIXP_DBL *pVdirectReal,
                      FIXP_DBL *pVdirectImag, FIXP_DBL *pVnonTrReal,
                      FIXP_DBL *pVnonTrImag, FIXP_DBL **ppDecorrInReal,
                      FIXP_DBL **ppDecorrInImag) {
  int k = 0;

  if (!isTrSlot(pTsdData, ts)) {
    *ppDecorrInReal = pVdirectReal;
    *ppDecorrInImag = pVdirectImag;
    return;
  }

  for (; k < TSD_START_BAND; k++) {
    pVnonTrReal[k] = pVdirectReal[k];
    pVnonTrImag[k] = pVdirectImag[k];
  }
  for (; k < numHybridBands; k++) {
    pVnonTrReal[k] = (FIXP_DBL)0;
    pVnonTrImag[k] = (FIXP_DBL)0;
  }
  *ppDecorrInReal = pVnonTrReal;
  *ppDecorrInImag = pVnonTrImag;
}

void TsdApply(const int numHybridBands, const TSD_DATA *pTsdData, int *pTsdTs,
              const FIXP_DBL *pVdirectReal, const FIXP_DBL *pVdirectImag,
              FIXP_DBL *pDnonTrReal, FIXP_DBL *pDnonTrImag) {
  const int ts = *pTsdTs;

  if (isTrSlot(pTsdData, ts)) {
    int k;
    const FIXP_DPK *phi = &phiTsd[pTsdData->bsTsdTrPhaseData[ts]];
    FDK_ASSERT((pTsdData->bsTsdTrPhaseData[ts] >= 0) &&
               (pTsdData->bsTsdTrPhaseData[ts] < 8));

    for (k = TSD_START_BAND; k < numHybridBands; k++) {
      FIXP_DBL tempReal, tempImag;
      cplxMultDiv2(&tempReal, &tempImag, pVdirectReal[k], pVdirectImag[k],
                   *phi);
      pDnonTrReal[k] = SATURATE_LEFT_SHIFT(
          (pDnonTrReal[k] >> 2) + (tempReal >> 1), 2, DFRACT_BITS);
      pDnonTrImag[k] = SATURATE_LEFT_SHIFT(
          (pDnonTrImag[k] >> 2) + (tempImag >> 1), 2, DFRACT_BITS);
    }
  }

  *pTsdTs = (ts + 1) & (MAX_TSD_TIME_SLOTS - 1);
  return;
}
