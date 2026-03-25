/******************* Library for basic calculation routines ********************

   Author(s):

   Description:

*******************************************************************************/

#include "dct.h"
#include "FDK_tools_rom.h"
#include "fft.h"

void dct_getTables(const FIXP_WTP **ptwiddle, const FIXP_STP **sin_twiddle,
                   int *sin_step, int length) {
  const FIXP_WTP *twiddle;
  int ld2_length;

  ld2_length = DFRACT_BITS - 1 - fNormz((FIXP_DBL)length) - 1;

  switch ((length) >> (ld2_length - 1)) {
    case 0x4:
      *sin_twiddle = SineTable1024;
      *sin_step = 1 << (10 - ld2_length);
      twiddle = windowSlopes[0][0][ld2_length - 1];
      break;
    case 0x7:
      *sin_twiddle = SineTable480;
      *sin_step = 1 << (8 - ld2_length);
      twiddle = windowSlopes[0][1][ld2_length];
      break;
    case 0x6:
      *sin_twiddle = SineTable384;
      *sin_step = 1 << (8 - ld2_length);
      twiddle = windowSlopes[0][2][ld2_length];
      break;
    case 0x5:
      *sin_twiddle = SineTable80;
      *sin_step = 1 << (6 - ld2_length);
      twiddle = windowSlopes[0][3][ld2_length];
      break;
    default:
      *sin_twiddle = NULL;
      *sin_step = 0;
      twiddle = NULL;
      break;
  }

  if (ptwiddle != NULL) {
    FDK_ASSERT(twiddle != NULL);
    *ptwiddle = twiddle;
  }

  FDK_ASSERT(*sin_step > 0);
}

#if !defined(FUNCTION_dct_III)
void dct_III(FIXP_DBL *pDat,
             FIXP_DBL *tmp,
             int L,
             int *pDat_e) {
  const FIXP_WTP *sin_twiddle;
  int i;
  FIXP_DBL xr, accu1, accu2;
  int inc, index;
  int M = L >> 1;

  FDK_ASSERT(L % 4 == 0);
  dct_getTables(NULL, &sin_twiddle, &inc, L);
  inc >>= 1;

  FIXP_DBL *pTmp_0 = &tmp[2];
  FIXP_DBL *pTmp_1 = &tmp[(M - 1) * 2];

  index = 4 * inc;

  for (i = 1; i<M>> 1; i++, pTmp_0 += 2, pTmp_1 -= 2) {
    FIXP_DBL accu3, accu4, accu5, accu6;

    cplxMultDiv2(&accu2, &accu1, pDat[L - i], pDat[i], sin_twiddle[i * inc]);
    cplxMultDiv2(&accu4, &accu3, pDat[M + i], pDat[M - i],
                 sin_twiddle[(M - i) * inc]);
    accu3 >>= 1;
    accu4 >>= 1;

    if (2 * i < (M / 2)) {
      cplxMultDiv2(&accu6, &accu5, (accu3 - (accu1 >> 1)),
                   ((accu2 >> 1) + accu4), sin_twiddle[index]);
    } else {
      cplxMultDiv2(&accu6, &accu5, ((accu2 >> 1) + accu4),
                   (accu3 - (accu1 >> 1)), sin_twiddle[index]);
      accu6 = -accu6;
    }
    xr = (accu1 >> 1) + accu3;
    pTmp_0[0] = (xr >> 1) - accu5;
    pTmp_1[0] = (xr >> 1) + accu5;

    xr = (accu2 >> 1) - accu4;
    pTmp_0[1] = (xr >> 1) - accu6;
    pTmp_1[1] = -((xr >> 1) + accu6);

    if (2 * i < ((M / 2) - 1)) {
      index += 4 * inc;
    } else if (2 * i >= ((M / 2))) {
      index -= 4 * inc;
    }
  }

  xr = fMultDiv2(pDat[M], sin_twiddle[M * inc].v.re);
  tmp[0] = ((pDat[0] >> 1) + xr) >> 1;
  tmp[1] = ((pDat[0] >> 1) - xr) >> 1;

  cplxMultDiv2(&accu2, &accu1, pDat[L - (M / 2)], pDat[M / 2],
               sin_twiddle[M * inc / 2]);
  tmp[M] = accu1 >> 1;
  tmp[M + 1] = accu2 >> 1;

  fft(M, tmp, pDat_e);

  pTmp_1 = &tmp[L];
  for (i = M >> 1; i--;) {
    FIXP_DBL tmp1, tmp2, tmp3, tmp4;
    tmp1 = *tmp++;
    tmp2 = *tmp++;
    tmp3 = *--pTmp_1;
    tmp4 = *--pTmp_1;
    *pDat++ = tmp1;
    *pDat++ = tmp3;
    *pDat++ = tmp2;
    *pDat++ = tmp4;
  }

  *pDat_e += 2;
}

void dst_III(FIXP_DBL *pDat,
             FIXP_DBL *tmp,
             int L,
             int *pDat_e) {
  int L2 = L >> 1;
  int i;
  FIXP_DBL t;

  for (i = 0; i < L2; i++) {
    t = pDat[i];
    pDat[i] = pDat[L - 1 - i];
    pDat[L - 1 - i] = t;
  }

  dct_III(pDat, tmp, L, pDat_e);

  for (i = 1; i < L; i += 2) pDat[i] = -pDat[i];
}

#endif

#if !defined(FUNCTION_dct_II)
void dct_II(
    FIXP_DBL *pDat,
    FIXP_DBL *tmp,
    int L,
    int *pDat_e) {
  const FIXP_WTP *sin_twiddle;
  FIXP_DBL accu1, accu2;
  FIXP_DBL *pTmp_0, *pTmp_1;

  int i;
  int inc, index = 0;
  int M = L >> 1;

  FDK_ASSERT(L % 4 == 0);
  dct_getTables(NULL, &sin_twiddle, &inc, L);
  inc >>= 1;

  {
    for (i = 0; i < M; i++) {
      tmp[i] = pDat[2 * i] >> 2;
      tmp[L - 1 - i] = pDat[2 * i + 1] >> 2;
    }
  }

  fft(M, tmp, pDat_e);

  pTmp_0 = &tmp[2];
  pTmp_1 = &tmp[(M - 1) * 2];

  index = inc * 4;

  for (i = 1; i<M>> 1; i++, pTmp_0 += 2, pTmp_1 -= 2) {
    FIXP_DBL a1, a2;
    FIXP_DBL accu3, accu4;

    a1 = ((pTmp_0[1] >> 1) + (pTmp_1[1] >> 1));
    a2 = ((pTmp_1[0] >> 1) - (pTmp_0[0] >> 1));

    if (2 * i < (M / 2)) {
      cplxMultDiv2(&accu1, &accu2, a2, a1, sin_twiddle[index]);
    } else {
      cplxMultDiv2(&accu1, &accu2, a1, a2, sin_twiddle[index]);
      accu1 = -accu1;
    }
    accu1 <<= 1;
    accu2 <<= 1;

    a1 = ((pTmp_0[0] >> 1) + (pTmp_1[0] >> 1));
    a2 = ((pTmp_0[1] >> 1) - (pTmp_1[1] >> 1));

    cplxMult(&accu3, &accu4, (accu1 + a2), (a1 + accu2), sin_twiddle[i * inc]);
    pDat[L - i] = -accu3;
    pDat[i] = accu4;

    cplxMult(&accu3, &accu4, (accu1 - a2), (a1 - accu2),
             sin_twiddle[(M - i) * inc]);
    pDat[M + i] = -accu3;
    pDat[M - i] = accu4;

    if (2 * i < ((M / 2) - 1)) {
      index += 4 * inc;
    } else if (2 * i >= ((M / 2))) {
      index -= 4 * inc;
    }
  }

  cplxMult(&accu1, &accu2, tmp[M], tmp[M + 1], sin_twiddle[(M / 2) * inc]);
  pDat[L - (M / 2)] = accu2;
  pDat[M / 2] = accu1;

  pDat[0] = tmp[0] + tmp[1];
  pDat[M] = fMult(tmp[0] - tmp[1],
                  sin_twiddle[M * inc].v.re);

  *pDat_e += 2;
}
#endif

#if !defined(FUNCTION_dct_IV)

void dct_IV(FIXP_DBL *pDat, int L, int *pDat_e) {
  int sin_step = 0;
  int M = L >> 1;

  const FIXP_WTP *twiddle;
  const FIXP_STP *sin_twiddle;

  FDK_ASSERT(L >= 4);

  FDK_ASSERT(L >= 4);

  dct_getTables(&twiddle, &sin_twiddle, &sin_step, L);

  {
    FIXP_DBL *RESTRICT pDat_0 = &pDat[0];
    FIXP_DBL *RESTRICT pDat_1 = &pDat[L - 2];
    int i;

    for (i = 0; i < M - 1; i += 2, pDat_0 += 2, pDat_1 -= 2) {
      FIXP_DBL accu1, accu2, accu3, accu4;

      accu1 = pDat_1[1];
      accu2 = pDat_0[0];
      accu3 = pDat_0[1];
      accu4 = pDat_1[0];

      cplxMultDiv2(&accu1, &accu2, accu1, accu2, twiddle[i]);
      cplxMultDiv2(&accu3, &accu4, accu4, accu3, twiddle[i + 1]);

      pDat_0[0] = accu2 >> 1;
      pDat_0[1] = accu1 >> 1;
      pDat_1[0] = accu4 >> 1;
      pDat_1[1] = -(accu3 >> 1);
    }
    if (M & 1) {
      FIXP_DBL accu1, accu2;

      accu1 = pDat_1[1];
      accu2 = pDat_0[0];

      cplxMultDiv2(&accu1, &accu2, accu1, accu2, twiddle[i]);

      pDat_0[0] = accu2 >> 1;
      pDat_0[1] = accu1 >> 1;
    }
  }

  fft(M, pDat, pDat_e);

  {
    FIXP_DBL *RESTRICT pDat_0 = &pDat[0];
    FIXP_DBL *RESTRICT pDat_1 = &pDat[L - 2];
    FIXP_DBL accu1, accu2, accu3, accu4;
    int idx, i;

    accu1 = pDat_1[0];
    accu2 = pDat_1[1];

    pDat_1[1] = -pDat_0[1];

    for (idx = sin_step, i = 1; i<(M + 1)>> 1; i++, idx += sin_step) {
      FIXP_STP twd = sin_twiddle[idx];
      cplxMult(&accu3, &accu4, accu1, accu2, twd);
      pDat_0[1] = accu3;
      pDat_1[0] = accu4;

      pDat_0 += 2;
      pDat_1 -= 2;

      cplxMult(&accu3, &accu4, pDat_0[1], pDat_0[0], twd);

      accu1 = pDat_1[0];
      accu2 = pDat_1[1];

      pDat_1[1] = -accu3;
      pDat_0[0] = accu4;
    }

    if ((M & 1) == 0) {
      accu1 = fMult(accu1, WTC(0x5a82799a));
      accu2 = fMult(accu2, WTC(0x5a82799a));

      pDat_1[0] = accu1 + accu2;
      pDat_0[1] = accu1 - accu2;
    }
  }

  *pDat_e += 2;
}
#endif /* defined (FUNCTION_dct_IV) */

#if !defined(FUNCTION_dst_IV)
void dst_IV(FIXP_DBL *pDat, int L, int *pDat_e) {
  int sin_step = 0;
  int M = L >> 1;

  const FIXP_WTP *twiddle;
  const FIXP_STP *sin_twiddle;

  FDK_ASSERT(L >= 4);

  FDK_ASSERT(L >= 4);

  dct_getTables(&twiddle, &sin_twiddle, &sin_step, L);

  {
    FIXP_DBL *RESTRICT pDat_0 = &pDat[0];
    FIXP_DBL *RESTRICT pDat_1 = &pDat[L - 2];
    int i;

    for (i = 0; i < M - 1; i += 2, pDat_0 += 2, pDat_1 -= 2) {
      FIXP_DBL accu1, accu2, accu3, accu4;

      accu1 = pDat_1[1] >> 1;
      accu2 = -(pDat_0[0] >> 1);
      accu3 = pDat_0[1] >> 1;
      accu4 = -(pDat_1[0] >> 1);

      cplxMultDiv2(&accu1, &accu2, accu1, accu2, twiddle[i]);
      cplxMultDiv2(&accu3, &accu4, accu4, accu3, twiddle[i + 1]);

      pDat_0[0] = accu2;
      pDat_0[1] = accu1;
      pDat_1[0] = accu4;
      pDat_1[1] = -accu3;
    }
    if (M & 1) {
      FIXP_DBL accu1, accu2;

      accu1 = pDat_1[1];
      accu2 = -pDat_0[0];

      cplxMultDiv2(&accu1, &accu2, accu1, accu2, twiddle[i]);

      pDat_0[0] = accu2 >> 1;
      pDat_0[1] = accu1 >> 1;
    }
  }

  fft(M, pDat, pDat_e);

  {
    FIXP_DBL *RESTRICT pDat_0;
    FIXP_DBL *RESTRICT pDat_1;
    FIXP_DBL accu1, accu2, accu3, accu4;
    int idx, i;

    pDat_0 = &pDat[0];
    pDat_1 = &pDat[L - 2];

    accu1 = pDat_1[0];
    accu2 = pDat_1[1];

    pDat_1[1] = -pDat_0[0];
    pDat_0[0] = pDat_0[1];

    for (idx = sin_step, i = 1; i<(M + 1)>> 1; i++, idx += sin_step) {
      FIXP_STP twd = sin_twiddle[idx];

      cplxMult(&accu3, &accu4, accu1, accu2, twd);
      pDat_1[0] = -accu3;
      pDat_0[1] = -accu4;

      pDat_0 += 2;
      pDat_1 -= 2;

      cplxMult(&accu3, &accu4, pDat_0[1], pDat_0[0], twd);

      accu1 = pDat_1[0];
      accu2 = pDat_1[1];

      pDat_0[0] = accu3;
      pDat_1[1] = -accu4;
    }

    if ((M & 1) == 0) {
      accu1 = fMult(accu1, WTC(0x5a82799a));
      accu2 = fMult(accu2, WTC(0x5a82799a));

      pDat_0[1] = -accu1 - accu2;
      pDat_1[0] = accu2 - accu1;
    }
  }

  *pDat_e += 2;
}
#endif /* !defined(FUNCTION_dst_IV) */
