#include "fft_rad2.h"
#include "FDK_tools_rom.h"

#define W_PiFOURTH STC(0x5a82799a)
#ifndef SUMDIFF_PIFOURTH
#define SUMDIFF_PIFOURTH(diff, sum, a, b) \
  {                                       \
    FIXP_DBL wa, wb;                      \
    wa = fMultDiv2(a, W_PiFOURTH);        \
    wb = fMultDiv2(b, W_PiFOURTH);        \
    diff = wb - wa;                       \
    sum = wb + wa;                        \
  }
#define SUMDIFF_PIFOURTH16(diff, sum, a, b)       \
  {                                               \
    FIXP_SGL wa, wb;                              \
    wa = FX_DBL2FX_SGL(fMultDiv2(a, W_PiFOURTH)); \
    wb = FX_DBL2FX_SGL(fMultDiv2(b, W_PiFOURTH)); \
    diff = wb - wa;                               \
    sum = wb + wa;                                \
  }
#endif

#define SCALEFACTOR2048 10
#define SCALEFACTOR1024 9
#define SCALEFACTOR512 8
#define SCALEFACTOR256 7
#define SCALEFACTOR128 6
#define SCALEFACTOR64 5
#define SCALEFACTOR32 4
#define SCALEFACTOR16 3
#define SCALEFACTOR8 2
#define SCALEFACTOR4 1
#define SCALEFACTOR2 1

#define SCALEFACTOR3 1
#define SCALEFACTOR5 1
#define SCALEFACTOR6 (SCALEFACTOR2 + SCALEFACTOR3 + 2)
#define SCALEFACTOR7 2
#define SCALEFACTOR9 2
#define SCALEFACTOR10 5
#define SCALEFACTOR12 3
#define SCALEFACTOR15 3
#define SCALEFACTOR18 (SCALEFACTOR2 + SCALEFACTOR9 + 2)
#define SCALEFACTOR20 (SCALEFACTOR4 + SCALEFACTOR5 + 2)
#define SCALEFACTOR21 (SCALEFACTOR3 + SCALEFACTOR7 + 2)
#define SCALEFACTOR24 (SCALEFACTOR2 + SCALEFACTOR12 + 2)
#define SCALEFACTOR30 (SCALEFACTOR2 + SCALEFACTOR15 + 2)
#define SCALEFACTOR40 (SCALEFACTOR5 + SCALEFACTOR8 + 2)
#define SCALEFACTOR48 (SCALEFACTOR4 + SCALEFACTOR12 + 2)
#define SCALEFACTOR60 (SCALEFACTOR4 + SCALEFACTOR15 + 2)
#define SCALEFACTOR80 (SCALEFACTOR5 + SCALEFACTOR16 + 2)
#define SCALEFACTOR96 (SCALEFACTOR3 + SCALEFACTOR32 + 2)
#define SCALEFACTOR120 (SCALEFACTOR8 + SCALEFACTOR15 + 2)
#define SCALEFACTOR160 (SCALEFACTOR10 + SCALEFACTOR16 + 2)
#define SCALEFACTOR168 (SCALEFACTOR21 + SCALEFACTOR8 + 2)
#define SCALEFACTOR192 (SCALEFACTOR12 + SCALEFACTOR16 + 2)
#define SCALEFACTOR240 (SCALEFACTOR16 + SCALEFACTOR15 + 2)
#define SCALEFACTOR320 (SCALEFACTOR10 + SCALEFACTOR32 + 2)
#define SCALEFACTOR336 (SCALEFACTOR21 + SCALEFACTOR16 + 2)
#define SCALEFACTOR384 (SCALEFACTOR12 + SCALEFACTOR32 + 2)
#define SCALEFACTOR480 (SCALEFACTOR32 + SCALEFACTOR15 + 2)

#include "fft.h"

#ifndef FUNCTION_fft2

static inline void fft2(FIXP_DBL *RESTRICT pDat) {
  FIXP_DBL r1, i1;
  FIXP_DBL r2, i2;

  r1 = pDat[2];
  r2 = pDat[0];

  i1 = pDat[3];
  i2 = pDat[1];

  pDat[0] = (r2 + r1) >> 1;
  pDat[2] = (r2 - r1) >> 1;

  pDat[1] = (i2 + i1) >> 1;
  pDat[3] = (i2 - i1) >> 1;
}
#endif

#define C31 (STC(0x91261468))

#ifndef FUNCTION_fft3
static inline void fft3(FIXP_DBL *RESTRICT pDat) {
  FIXP_DBL r1, r2;
  FIXP_DBL s1, s2;
  FIXP_DBL pD;

  r1 = pDat[2] + pDat[4];
  r2 = fMultDiv2((pDat[2] - pDat[4]), C31);
  pD = pDat[0] >> 1;
  pDat[0] = pD + (r1 >> 1);
  r1 = pD - (r1 >> 2);

  s1 = pDat[3] + pDat[5];
  s2 = fMultDiv2((pDat[3] - pDat[5]), C31);
  pD = pDat[1] >> 1;
  pDat[1] = pD + (s1 >> 1);
  s1 = pD - (s1 >> 2);

  pDat[2] = r1 - s2;
  pDat[4] = r1 + s2;
  pDat[3] = s1 + r2;
  pDat[5] = s1 - r2;
}
#endif

#define F5C(x) STC(x)

#define C51 (F5C(0x79bc3854))
#define C52 (F5C(0x9d839db0))
#define C53 (F5C(0xd18053ce))
#define C54 (F5C(0x478dde64))
#define C55 (F5C(0xb0000001))

static inline void fft5(FIXP_DBL *RESTRICT pDat) {
  FIXP_DBL r1, r2, r3, r4;
  FIXP_DBL s1, s2, s3, s4;
  FIXP_DBL t;

  r1 = (pDat[2] + pDat[8]) >> 1;
  r4 = (pDat[2] - pDat[8]) >> 1;
  r3 = (pDat[4] + pDat[6]) >> 1;
  r2 = (pDat[4] - pDat[6]) >> 1;
  t = fMult((r1 - r3), C54);
  r1 = r1 + r3;
  pDat[0] = (pDat[0] >> 1) + r1;
  r1 = pDat[0] + (fMultDiv2(r1, C55) << (2));
  r3 = r1 - t;
  r1 = r1 + t;
  t = fMult((r4 + r2), C51);
  r4 = t + (fMultDiv2(r4, C52) << (2));
  r2 = t + fMult(r2, C53);

  s1 = (pDat[3] + pDat[9]) >> 1;
  s4 = (pDat[3] - pDat[9]) >> 1;
  s3 = (pDat[5] + pDat[7]) >> 1;
  s2 = (pDat[5] - pDat[7]) >> 1;
  t = fMult((s1 - s3), C54);
  s1 = s1 + s3;
  pDat[1] = (pDat[1] >> 1) + s1;
  s1 = pDat[1] + (fMultDiv2(s1, C55) << (2));
  s3 = s1 - t;
  s1 = s1 + t;
  t = fMult((s4 + s2), C51);
  s4 = t + (fMultDiv2(s4, C52) << (2));
  s2 = t + fMult(s2, C53);

  pDat[2] = r1 + s2;
  pDat[8] = r1 - s2;
  pDat[4] = r3 - s4;
  pDat[6] = r3 + s4;

  pDat[3] = s1 - r2;
  pDat[9] = s1 + r2;
  pDat[5] = s3 + r4;
  pDat[7] = s3 - r4;
}

#define F5C(x) STC(x)

#define C51 (F5C(0x79bc3854))
#define C52 (F5C(0x9d839db0))
#define C53 (F5C(0xd18053ce))
#define C54 (F5C(0x478dde64))
#define C55 (F5C(0xb0000001))
static void fft10(FIXP_DBL *x)
{
  FIXP_DBL t;
  FIXP_DBL x0, x1, x2, x3, x4;
  FIXP_DBL r1, r2, r3, r4;
  FIXP_DBL s1, s2, s3, s4;
  FIXP_DBL y00, y01, y02, y03, y04, y05, y06, y07, y08, y09;
  FIXP_DBL y10, y11, y12, y13, y14, y15, y16, y17, y18, y19;

  const int s = 1;

  x0 = (x[s * 0] >> SCALEFACTOR10);
  x1 = (x[s * 4] >> SCALEFACTOR10);
  x2 = (x[s * 8] >> SCALEFACTOR10);
  x3 = (x[s * 12] >> SCALEFACTOR10);
  x4 = (x[s * 16] >> SCALEFACTOR10);

  r1 = (x3 + x2);
  r4 = (x3 - x2);
  r3 = (x1 + x4);
  r2 = (x1 - x4);
  t = fMult((r1 - r3), C54);
  r1 = (r1 + r3);
  y00 = (x0 + r1);
  r1 = (y00 + ((fMult(r1, C55) << 1)));
  r3 = (r1 - t);
  r1 = (r1 + t);
  t = fMult((r4 + r2), C51);
  r4 = (t + (fMult(r4, C52) << 1));
  r2 = (t + fMult(r2, C53));

  x0 = (x[s * 0 + 1] >> SCALEFACTOR10);
  x1 = (x[s * 4 + 1] >> SCALEFACTOR10);
  x2 = (x[s * 8 + 1] >> SCALEFACTOR10);
  x3 = (x[s * 12 + 1] >> SCALEFACTOR10);
  x4 = (x[s * 16 + 1] >> SCALEFACTOR10);

  s1 = (x3 + x2);
  s4 = (x3 - x2);
  s3 = (x1 + x4);
  s2 = (x1 - x4);
  t = fMult((s1 - s3), C54);
  s1 = (s1 + s3);
  y01 = (x0 + s1);
  s1 = (y01 + (fMult(s1, C55) << 1));
  s3 = (s1 - t);
  s1 = (s1 + t);
  t = fMult((s4 + s2), C51);
  s4 = (t + (fMult(s4, C52) << 1));
  s2 = (t + fMult(s2, C53));

  y04 = (r1 + s2);
  y16 = (r1 - s2);
  y08 = (r3 - s4);
  y12 = (r3 + s4);

  y05 = (s1 - r2);
  y17 = (s1 + r2);
  y09 = (s3 + r4);
  y13 = (s3 - r4);

  x0 = (x[s * 10] >> SCALEFACTOR10);
  x1 = (x[s * 2] >> SCALEFACTOR10);
  x2 = (x[s * 6] >> SCALEFACTOR10);
  x3 = (x[s * 14] >> SCALEFACTOR10);
  x4 = (x[s * 18] >> SCALEFACTOR10);

  r1 = (x1 + x4);
  r4 = (x1 - x4);
  r3 = (x3 + x2);
  r2 = (x3 - x2);
  t = fMult((r1 - r3), C54);
  r1 = (r1 + r3);
  y02 = (x0 + r1);
  r1 = (y02 + ((fMult(r1, C55) << 1)));
  r3 = (r1 - t);
  r1 = (r1 + t);
  t = fMult(((r4 + r2)), C51);
  r4 = (t + (fMult(r4, C52) << 1));
  r2 = (t + fMult(r2, C53));

  x0 = (x[s * 10 + 1] >> SCALEFACTOR10);
  x1 = (x[s * 2 + 1] >> SCALEFACTOR10);
  x2 = (x[s * 6 + 1] >> SCALEFACTOR10);
  x3 = (x[s * 14 + 1] >> SCALEFACTOR10);
  x4 = (x[s * 18 + 1] >> SCALEFACTOR10);

  s1 = (x1 + x4);
  s4 = (x1 - x4);
  s3 = (x3 + x2);
  s2 = (x3 - x2);
  t = fMult((s1 - s3), C54);
  s1 = (s1 + s3);
  y03 = (x0 + s1);
  s1 = (y03 + (fMult(s1, C55) << 1));
  s3 = (s1 - t);
  s1 = (s1 + t);
  t = fMult((s4 + s2), C51);
  s4 = (t + (fMult(s4, C52) << 1));
  s2 = (t + fMult(s2, C53));

  y06 = (r1 + s2);
  y18 = (r1 - s2);
  y10 = (r3 - s4);
  y14 = (r3 + s4);

  y07 = (s1 - r2);
  y19 = (s1 + r2);
  y11 = (s3 + r4);
  y15 = (s3 - r4);

  x[s * 0] = (y00 + y02);
  x[s * 0 + 1] = (y01 + y03);
  x[s * 10] = (y00 - y02);
  x[s * 10 + 1] = (y01 - y03);

  x[s * 4] = (y04 + y06);
  x[s * 4 + 1] = (y05 + y07);
  x[s * 14] = (y04 - y06);
  x[s * 14 + 1] = (y05 - y07);

  x[s * 8] = (y08 + y10);
  x[s * 8 + 1] = (y09 + y11);
  x[s * 18] = (y08 - y10);
  x[s * 18 + 1] = (y09 - y11);

  x[s * 12] = (y12 + y14);
  x[s * 12 + 1] = (y13 + y15);
  x[s * 2] = (y12 - y14);
  x[s * 2 + 1] = (y13 - y15);

  x[s * 16] = (y16 + y18);
  x[s * 16 + 1] = (y17 + y19);
  x[s * 6] = (y16 - y18);
  x[s * 6 + 1] = (y17 - y19);
}

#ifndef FUNCTION_fft12
#define FUNCTION_fft12

#undef C31
#define C31 (STC(0x91261468))

static inline void fft12(FIXP_DBL *pInput) {
  FIXP_DBL aDst[24];
  FIXP_DBL *pSrc, *pDst;
  int i;

  pSrc = pInput;
  pDst = aDst;
  FIXP_DBL r1, r2, s1, s2, pD;

  r1 = pSrc[8] + pSrc[16];
  r2 = fMultDiv2((pSrc[8] - pSrc[16]), C31);
  pD = pSrc[0] >> 1;
  pDst[0] = (pD + (r1 >> 1)) >> 1;
  r1 = pD - (r1 >> 2);

  s1 = pSrc[9] + pSrc[17];
  s2 = fMultDiv2((pSrc[9] - pSrc[17]), C31);
  pD = pSrc[1] >> 1;
  pDst[1] = (pD + (s1 >> 1)) >> 1;
  s1 = pD - (s1 >> 2);

  pDst[2] = (r1 - s2) >> 1;
  pDst[3] = (s1 + r2) >> 1;
  pDst[4] = (r1 + s2) >> 1;
  pDst[5] = (s1 - r2) >> 1;
  pSrc += 2;
  pDst += 6;

  const FIXP_STB *pVecRe = RotVectorReal12;
  const FIXP_STB *pVecIm = RotVectorImag12;
  FIXP_DBL re, im;
  FIXP_STB vre, vim;
  for (i = 0; i < 2; i++) {

    r1 = pSrc[8] + pSrc[16];
    r2 = fMultDiv2((pSrc[8] - pSrc[16]), C31);
    pD = pSrc[0] >> 1;
    pDst[0] = (pD + (r1 >> 1)) >> 1;
    r1 = pD - (r1 >> 2);

    s1 = pSrc[9] + pSrc[17];
    s2 = fMultDiv2((pSrc[9] - pSrc[17]), C31);
    pD = pSrc[1] >> 1;
    pDst[1] = (pD + (s1 >> 1)) >> 1;
    s1 = pD - (s1 >> 2);

    re = (r1 - s2) >> 0;
    im = (s1 + r2) >> 0;
    vre = *pVecRe++;
    vim = *pVecIm++;
    cplxMultDiv2(&pDst[3], &pDst[2], im, re, vre, vim);

    re = (r1 + s2) >> 0;
    im = (s1 - r2) >> 0;
    vre = *pVecRe++;
    vim = *pVecIm++;
    cplxMultDiv2(&pDst[5], &pDst[4], im, re, vre, vim);

    pDst += 6;
    pSrc += 2;
  }
  r1 = pSrc[8] + pSrc[16];
  r2 = fMultDiv2((pSrc[8] - pSrc[16]), C31);
  pD = pSrc[0] >> 1;
  pDst[0] = (pD + (r1 >> 1)) >> 1;
  r1 = pD - (r1 >> 2);

  s1 = pSrc[9] + pSrc[17];
  s2 = fMultDiv2((pSrc[9] - pSrc[17]), C31);
  pD = pSrc[1] >> 1;
  pDst[1] = (pD + (s1 >> 1)) >> 1;
  s1 = pD - (s1 >> 2);

  pDst[2] = (s1 + r2) >> 1;
  pDst[3] = (s2 - r1) >> 1;
  pDst[4] = -((r1 + s2) >> 1);
  pDst[5] = (r2 - s1) >> 1;

  pSrc = aDst;
  pDst = pInput;
  for (i = 0; i < 3; i++) {
    FIXP_DBL a00, a10, a20, a30, tmp0, tmp1;

    a00 = (pSrc[0] + pSrc[12]) >> 1;
    a10 = (pSrc[6] + pSrc[18]) >> 1;
    a20 = (pSrc[1] + pSrc[13]) >> 1;
    a30 = (pSrc[7] + pSrc[19]) >> 1;

    pDst[0] = a00 + a10;
    pDst[1] = a20 + a30;

    tmp0 = a00 - pSrc[12];
    tmp1 = a20 - pSrc[13];

    pDst[12] = a00 - a10;
    pDst[13] = a20 - a30;

    a10 = a10 - pSrc[18];
    a30 = a30 - pSrc[19];

    pDst[6] = tmp0 + a30;
    pDst[18] = tmp0 - a30;
    pDst[7] = tmp1 - a10;
    pDst[19] = tmp1 + a10;

    pSrc += 2;
    pDst += 2;
  }
}
#endif

#ifndef FUNCTION_fft15

#define N3 3
#define N5 5
#define N6 6
#define N15 15

static inline void fft15(FIXP_DBL *pInput) {
  FIXP_DBL aDst[2 * N15];
  FIXP_DBL aDst1[2 * N15];
  int i, k, l;

  {
    const FIXP_DBL *pSrc = pInput;
    FIXP_DBL *RESTRICT pDst = aDst;
    for (i = 0, l = 0, k = 0; i < N5; i++, k += 6) {
      pDst[k + 0] = pSrc[l];
      pDst[k + 1] = pSrc[l + 1];
      l += 2 * N5;
      if (l >= (2 * N15)) l -= (2 * N15);

      pDst[k + 2] = pSrc[l];
      pDst[k + 3] = pSrc[l + 1];
      l += 2 * N5;
      if (l >= (2 * N15)) l -= (2 * N15);
      pDst[k + 4] = pSrc[l];
      pDst[k + 5] = pSrc[l + 1];
      l += (2 * N5) + (2 * N3);
      if (l >= (2 * N15)) l -= (2 * N15);

      FIXP_DBL r1, r2, r3;
      FIXP_DBL s1, s2;
      r1 = pDst[k + 2] + pDst[k + 4];
      r2 = fMult((pDst[k + 2] - pDst[k + 4]), C31);
      s1 = pDst[k + 0];
      pDst[k + 0] = (s1 + r1) >> 2;
      r1 = s1 - (r1 >> 1);

      s1 = pDst[k + 3] + pDst[k + 5];
      s2 = fMult((pDst[k + 3] - pDst[k + 5]), C31);
      r3 = pDst[k + 1];
      pDst[k + 1] = (r3 + s1) >> 2;
      s1 = r3 - (s1 >> 1);

      pDst[k + 2] = (r1 - s2) >> 2;
      pDst[k + 4] = (r1 + s2) >> 2;
      pDst[k + 3] = (s1 + r2) >> 2;
      pDst[k + 5] = (s1 - r2) >> 2;
    }
  }
  {
    const FIXP_DBL *pSrc = aDst;
    FIXP_DBL *RESTRICT pDst = aDst1;
    for (i = 0, l = 0, k = 0; i < N3; i++, k += 10) {
      l = 2 * i;
      pDst[k + 0] = pSrc[l + 0];
      pDst[k + 1] = pSrc[l + 1];
      pDst[k + 2] = pSrc[l + 0 + (2 * N3)];
      pDst[k + 3] = pSrc[l + 1 + (2 * N3)];
      pDst[k + 4] = pSrc[l + 0 + (4 * N3)];
      pDst[k + 5] = pSrc[l + 1 + (4 * N3)];
      pDst[k + 6] = pSrc[l + 0 + (6 * N3)];
      pDst[k + 7] = pSrc[l + 1 + (6 * N3)];
      pDst[k + 8] = pSrc[l + 0 + (8 * N3)];
      pDst[k + 9] = pSrc[l + 1 + (8 * N3)];
      fft5(&pDst[k]);
    }
  }
  {
    const FIXP_DBL *pSrc = aDst1;
    FIXP_DBL *RESTRICT pDst = pInput;
    for (i = 0, l = 0, k = 0; i < N3; i++, k += 10) {
      pDst[k + 0] = pSrc[l];
      pDst[k + 1] = pSrc[l + 1];
      l += (2 * N6);
      if (l >= (2 * N15)) l -= (2 * N15);
      pDst[k + 2] = pSrc[l];
      pDst[k + 3] = pSrc[l + 1];
      l += (2 * N6);
      if (l >= (2 * N15)) l -= (2 * N15);
      pDst[k + 4] = pSrc[l];
      pDst[k + 5] = pSrc[l + 1];
      l += (2 * N6);
      if (l >= (2 * N15)) l -= (2 * N15);
      pDst[k + 6] = pSrc[l];
      pDst[k + 7] = pSrc[l + 1];
      l += (2 * N6);
      if (l >= (2 * N15)) l -= (2 * N15);
      pDst[k + 8] = pSrc[l];
      pDst[k + 9] = pSrc[l + 1];
      l += 2;
    }
  }
}
#endif

#if !defined(__arm__)
#define SHIFT_A >> 1
#define SHIFT_B
#else
#define SHIFT_A
#define SHIFT_B >> 1
#endif

#ifndef FUNCTION_fft_16
#define FUNCTION_DATA_fft_16_w16
static const FIXP_STP fft16_w16[2] = {STCP(0x7641af3d, 0x30fbc54d),
                                      STCP(0x30fbc54d, 0x7641af3d)};

LNK_SECTION_CODE_L1
inline void fft_16(FIXP_DBL *RESTRICT x) {
  FIXP_DBL vr, ur;
  FIXP_DBL vr2, ur2;
  FIXP_DBL vr3, ur3;
  FIXP_DBL vr4, ur4;
  FIXP_DBL vi, ui;
  FIXP_DBL vi2, ui2;
  FIXP_DBL vi3, ui3;

  vr = (x[0] >> 1) + (x[16] >> 1);
  ur = (x[1] >> 1) + (x[17] >> 1);
  vi = (x[8] SHIFT_A) + (x[24] SHIFT_A);
  ui = (x[9] SHIFT_A) + (x[25] SHIFT_A);
  x[0] = vr + (vi SHIFT_B);
  x[1] = ur + (ui SHIFT_B);

  vr2 = (x[4] >> 1) + (x[20] >> 1);
  ur2 = (x[5] >> 1) + (x[21] >> 1);

  x[4] = vr - (vi SHIFT_B);
  x[5] = ur - (ui SHIFT_B);
  vr -= x[16];
  vi = (vi SHIFT_B)-x[24];
  ur -= x[17];
  ui = (ui SHIFT_B)-x[25];

  vr3 = (x[2] >> 1) + (x[18] >> 1);
  ur3 = (x[3] >> 1) + (x[19] >> 1);

  x[2] = ui + vr;
  x[3] = ur - vi;

  vr4 = (x[6] >> 1) + (x[22] >> 1);
  ur4 = (x[7] >> 1) + (x[23] >> 1);

  x[6] = vr - ui;
  x[7] = vi + ur;

  vi2 = (x[12] SHIFT_A) + (x[28] SHIFT_A);
  ui2 = (x[13] SHIFT_A) + (x[29] SHIFT_A);
  x[8] = vr2 + (vi2 SHIFT_B);
  x[9] = ur2 + (ui2 SHIFT_B);
  x[12] = vr2 - (vi2 SHIFT_B);
  x[13] = ur2 - (ui2 SHIFT_B);
  vr2 -= x[20];
  ur2 -= x[21];
  vi2 = (vi2 SHIFT_B)-x[28];
  ui2 = (ui2 SHIFT_B)-x[29];

  vi = (x[10] SHIFT_A) + (x[26] SHIFT_A);
  ui = (x[11] SHIFT_A) + (x[27] SHIFT_A);

  x[10] = ui2 + vr2;
  x[11] = ur2 - vi2;

  vi3 = (x[14] SHIFT_A) + (x[30] SHIFT_A);
  ui3 = (x[15] SHIFT_A) + (x[31] SHIFT_A);

  x[14] = vr2 - ui2;
  x[15] = vi2 + ur2;

  x[16] = vr3 + (vi SHIFT_B);
  x[17] = ur3 + (ui SHIFT_B);
  x[20] = vr3 - (vi SHIFT_B);
  x[21] = ur3 - (ui SHIFT_B);
  vr3 -= x[18];
  ur3 -= x[19];
  vi = (vi SHIFT_B)-x[26];
  ui = (ui SHIFT_B)-x[27];
  x[18] = ui + vr3;
  x[19] = ur3 - vi;

  x[24] = vr4 + (vi3 SHIFT_B);
  x[28] = vr4 - (vi3 SHIFT_B);
  x[25] = ur4 + (ui3 SHIFT_B);
  x[29] = ur4 - (ui3 SHIFT_B);
  vr4 -= x[22];
  ur4 -= x[23];

  x[22] = vr3 - ui;
  x[23] = vi + ur3;

  vi3 = (vi3 SHIFT_B)-x[30];
  ui3 = (ui3 SHIFT_B)-x[31];
  x[26] = ui3 + vr4;
  x[30] = vr4 - ui3;
  x[27] = ur4 - vi3;
  x[31] = vi3 + ur4;

  vr = x[8];
  vi = x[9];
  ur = x[0] >> 1;
  ui = x[1] >> 1;
  x[0] = ur + (vr >> 1);
  x[1] = ui + (vi >> 1);
  x[8] = ur - (vr >> 1);
  x[9] = ui - (vi >> 1);

  vr = x[13];
  vi = x[12];
  ur = x[4] >> 1;
  ui = x[5] >> 1;
  x[4] = ur + (vr >> 1);
  x[5] = ui - (vi >> 1);
  x[12] = ur - (vr >> 1);
  x[13] = ui + (vi >> 1);

  vr = x[24];
  vi = x[25];
  ur = x[16] >> 1;
  ui = x[17] >> 1;
  x[16] = ur + (vr >> 1);
  x[17] = ui + (vi >> 1);
  x[24] = ur - (vr >> 1);
  x[25] = ui - (vi >> 1);

  vr = x[29];
  vi = x[28];
  ur = x[20] >> 1;
  ui = x[21] >> 1;
  x[20] = ur + (vr >> 1);
  x[21] = ui - (vi >> 1);
  x[28] = ur - (vr >> 1);
  x[29] = ui + (vi >> 1);

  SUMDIFF_PIFOURTH(vi, vr, x[10], x[11])
  ur = x[2];
  ui = x[3];
  x[2] = (ur >> 1) + vr;
  x[3] = (ui >> 1) + vi;
  x[10] = (ur >> 1) - vr;
  x[11] = (ui >> 1) - vi;

  SUMDIFF_PIFOURTH(vr, vi, x[14], x[15])
  ur = x[6];
  ui = x[7];
  x[6] = (ur >> 1) + vr;
  x[7] = (ui >> 1) - vi;
  x[14] = (ur >> 1) - vr;
  x[15] = (ui >> 1) + vi;

  SUMDIFF_PIFOURTH(vi, vr, x[26], x[27])
  ur = x[18];
  ui = x[19];
  x[18] = (ur >> 1) + vr;
  x[19] = (ui >> 1) + vi;
  x[26] = (ur >> 1) - vr;
  x[27] = (ui >> 1) - vi;

  SUMDIFF_PIFOURTH(vr, vi, x[30], x[31])
  ur = x[22];
  ui = x[23];
  x[22] = (ur >> 1) + vr;
  x[23] = (ui >> 1) - vi;
  x[30] = (ur >> 1) - vr;
  x[31] = (ui >> 1) + vi;

  vr = x[16];
  vi = x[17];
  ur = x[0] >> 1;
  ui = x[1] >> 1;
  x[0] = ur + (vr >> 1);
  x[1] = ui + (vi >> 1);
  x[16] = ur - (vr >> 1);
  x[17] = ui - (vi >> 1);

  vi = x[24];
  vr = x[25];
  ur = x[8] >> 1;
  ui = x[9] >> 1;
  x[8] = ur + (vr >> 1);
  x[9] = ui - (vi >> 1);
  x[24] = ur - (vr >> 1);
  x[25] = ui + (vi >> 1);

  cplxMultDiv2(&vi, &vr, x[19], x[18], fft16_w16[0]);
  ur = x[2];
  ui = x[3];
  x[2] = (ur >> 1) + vr;
  x[3] = (ui >> 1) + vi;
  x[18] = (ur >> 1) - vr;
  x[19] = (ui >> 1) - vi;

  cplxMultDiv2(&vr, &vi, x[27], x[26], fft16_w16[0]);
  ur = x[10];
  ui = x[11];
  x[10] = (ur >> 1) + vr;
  x[11] = (ui >> 1) - vi;
  x[26] = (ur >> 1) - vr;
  x[27] = (ui >> 1) + vi;

  SUMDIFF_PIFOURTH(vi, vr, x[20], x[21])
  ur = x[4];
  ui = x[5];
  x[4] = (ur >> 1) + vr;
  x[5] = (ui >> 1) + vi;
  x[20] = (ur >> 1) - vr;
  x[21] = (ui >> 1) - vi;

  SUMDIFF_PIFOURTH(vr, vi, x[28], x[29])
  ur = x[12];
  ui = x[13];
  x[12] = (ur >> 1) + vr;
  x[13] = (ui >> 1) - vi;
  x[28] = (ur >> 1) - vr;
  x[29] = (ui >> 1) + vi;

  cplxMultDiv2(&vi, &vr, x[23], x[22], fft16_w16[1]);
  ur = x[6];
  ui = x[7];
  x[6] = (ur >> 1) + vr;
  x[7] = (ui >> 1) + vi;
  x[22] = (ur >> 1) - vr;
  x[23] = (ui >> 1) - vi;

  cplxMultDiv2(&vr, &vi, x[31], x[30], fft16_w16[1]);
  ur = x[14];
  ui = x[15];
  x[14] = (ur >> 1) + vr;
  x[15] = (ui >> 1) - vi;
  x[30] = (ur >> 1) - vr;
  x[31] = (ui >> 1) + vi;
}
#endif

#ifndef FUNCTION_fft_32
static const FIXP_STP fft32_w32[6] = {
    STCP(0x7641af3d, 0x30fbc54d), STCP(0x30fbc54d, 0x7641af3d),
    STCP(0x7d8a5f40, 0x18f8b83c), STCP(0x6a6d98a4, 0x471cece7),
    STCP(0x471cece7, 0x6a6d98a4), STCP(0x18f8b83c, 0x7d8a5f40)};
#define W_PiFOURTH STC(0x5a82799a)

LNK_SECTION_CODE_L1
inline void fft_32(FIXP_DBL *const _x) {
  {
    FIXP_DBL *const x = _x;
    FIXP_DBL vi, ui;
    FIXP_DBL vi2, ui2;
    FIXP_DBL vi3, ui3;
    FIXP_DBL vr, ur;
    FIXP_DBL vr2, ur2;
    FIXP_DBL vr3, ur3;
    FIXP_DBL vr4, ur4;

    vr = (x[0] + x[32]) >> 1;     /* Re A + Re B */
    ur = (x[1] + x[33]) >> 1;     /* Im A + Im B */
    vi = (x[16] + x[48]) SHIFT_A; /* Re C + Re D */
    ui = (x[17] + x[49]) SHIFT_A; /* Im C + Im D */

    x[0] = vr + (vi SHIFT_B); /* Re A' = ReA + ReB +ReC + ReD */
    x[1] = ur + (ui SHIFT_B); /* Im A' = sum of imag values */

    vr2 = (x[4] + x[36]) >> 1; /* Re A + Re B */
    ur2 = (x[5] + x[37]) >> 1; /* Im A + Im B */

    x[4] = vr - (vi SHIFT_B); /* Re C' = -(ReC+ReD) + (ReA+ReB) */
    x[5] = ur - (ui SHIFT_B); /* Im C' = -Im C -Im D +Im A +Im B */

    vr -= x[32];             /* Re A - Re B */
    ur -= x[33];             /* Im A - Im B */
    vi = (vi SHIFT_B)-x[48]; /* Re C - Re D */
    ui = (ui SHIFT_B)-x[49]; /* Im C - Im D */

    vr3 = (x[2] + x[34]) >> 1; /* Re A + Re B */
    ur3 = (x[3] + x[35]) >> 1; /* Im A + Im B */

    x[2] = ui + vr; /* Re B' = Im C - Im D  + Re A - Re B */
    x[3] = ur - vi; /* Im B'= -Re C + Re D + Im A - Im B */

    vr4 = (x[6] + x[38]) >> 1; /* Re A + Re B */
    ur4 = (x[7] + x[39]) >> 1; /* Im A + Im B */

    x[6] = vr - ui; /* Re D' = -Im C + Im D + Re A - Re B */
    x[7] = vi + ur; /* Im D'= Re C - Re D + Im A - Im B */

    vi = (x[20] + x[52]) SHIFT_A; /* Re C + Re D */
    ui = (x[21] + x[53]) SHIFT_A; /* Im C + Im D */

    x[16] = vr2 + (vi SHIFT_B); /* Re A' = ReA + ReB +ReC + ReD */
    x[17] = ur2 + (ui SHIFT_B); /* Im A' = sum of imag values */
    x[20] = vr2 - (vi SHIFT_B); /* Re C' = -(ReC+ReD) + (ReA+ReB) */
    x[21] = ur2 - (ui SHIFT_B); /* Im C' = -Im C -Im D +Im A +Im B */

    vr2 -= x[36];            /* Re A - Re B */
    ur2 -= x[37];            /* Im A - Im B */
    vi = (vi SHIFT_B)-x[52]; /* Re C - Re D */
    ui = (ui SHIFT_B)-x[53]; /* Im C - Im D */

    vi2 = (x[18] + x[50]) SHIFT_A; /* Re C + Re D */
    ui2 = (x[19] + x[51]) SHIFT_A; /* Im C + Im D */

    x[18] = ui + vr2; /* Re B' = Im C - Im D  + Re A - Re B */
    x[19] = ur2 - vi; /* Im B'= -Re C + Re D + Im A - Im B */

    vi3 = (x[22] + x[54]) SHIFT_A; /* Re C + Re D */
    ui3 = (x[23] + x[55]) SHIFT_A; /* Im C + Im D */

    x[22] = vr2 - ui; /* Re D' = -Im C + Im D + Re A - Re B */
    x[23] = vi + ur2; /* Im D'= Re C - Re D + Im A - Im B */

    x[32] = vr3 + (vi2 SHIFT_B); /* Re A' = ReA + ReB +ReC + ReD */
    x[33] = ur3 + (ui2 SHIFT_B); /* Im A' = sum of imag values */
    x[36] = vr3 - (vi2 SHIFT_B); /* Re C' = -(ReC+ReD) + (ReA+ReB) */
    x[37] = ur3 - (ui2 SHIFT_B); /* Im C' = -Im C -Im D +Im A +Im B */

    vr3 -= x[34];              /* Re A - Re B */
    ur3 -= x[35];              /* Im A - Im B */
    vi2 = (vi2 SHIFT_B)-x[50]; /* Re C - Re D */
    ui2 = (ui2 SHIFT_B)-x[51]; /* Im C - Im D */

    x[34] = ui2 + vr3; /* Re B' = Im C - Im D  + Re A - Re B */
    x[35] = ur3 - vi2; /* Im B'= -Re C + Re D + Im A - Im B */

    x[48] = vr4 + (vi3 SHIFT_B); /* Re A' = ReA + ReB +ReC + ReD */
    x[52] = vr4 - (vi3 SHIFT_B); /* Re C' = -(ReC+ReD) + (ReA+ReB) */
    x[49] = ur4 + (ui3 SHIFT_B); /* Im A' = sum of imag values */
    x[53] = ur4 - (ui3 SHIFT_B); /* Im C' = -Im C -Im D +Im A +Im B */

    vr4 -= x[38]; /* Re A - Re B */
    ur4 -= x[39]; /* Im A - Im B */

    x[38] = vr3 - ui2; /* Re D' = -Im C + Im D + Re A - Re B */
    x[39] = vi2 + ur3; /* Im D'= Re C - Re D + Im A - Im B */

    vi3 = (vi3 SHIFT_B)-x[54]; /* Re C - Re D */
    ui3 = (ui3 SHIFT_B)-x[55]; /* Im C - Im D */

    x[50] = ui3 + vr4; /* Re B' = Im C - Im D  + Re A - Re B */
    x[54] = vr4 - ui3; /* Re D' = -Im C + Im D + Re A - Re B */
    x[51] = ur4 - vi3; /* Im B'= -Re C + Re D + Im A - Im B */
    x[55] = vi3 + ur4; /* Im D'= Re C - Re D + Im A - Im B */

    vr = (x[8] + x[40]) >> 1;     /* Re A + Re B */
    ur = (x[9] + x[41]) >> 1;     /* Im A + Im B */
    vi = (x[24] + x[56]) SHIFT_A; /* Re C + Re D */
    ui = (x[25] + x[57]) SHIFT_A; /* Im C + Im D */

    x[8] = vr + (vi SHIFT_B); /* Re A' = ReA + ReB +ReC + ReD */
    x[9] = ur + (ui SHIFT_B); /* Im A' = sum of imag values */

    vr2 = (x[12] + x[44]) >> 1; /* Re A + Re B */
    ur2 = (x[13] + x[45]) >> 1; /* Im A + Im B */

    x[12] = vr - (vi SHIFT_B); /* Re C' = -(ReC+ReD) + (ReA+ReB) */
    x[13] = ur - (ui SHIFT_B); /* Im C' = -Im C -Im D +Im A +Im B */

    vr -= x[40];             /* Re A - Re B */
    ur -= x[41];             /* Im A - Im B */
    vi = (vi SHIFT_B)-x[56]; /* Re C - Re D */
    ui = (ui SHIFT_B)-x[57]; /* Im C - Im D */

    vr3 = (x[10] + x[42]) >> 1; /* Re A + Re B */
    ur3 = (x[11] + x[43]) >> 1; /* Im A + Im B */

    x[10] = ui + vr; /* Re B' = Im C - Im D  + Re A - Re B */
    x[11] = ur - vi; /* Im B'= -Re C + Re D + Im A - Im B */

    vr4 = (x[14] + x[46]) >> 1; /* Re A + Re B */
    ur4 = (x[15] + x[47]) >> 1; /* Im A + Im B */

    x[14] = vr - ui; /* Re D' = -Im C + Im D + Re A - Re B */
    x[15] = vi + ur; /* Im D'= Re C - Re D + Im A - Im B */

    vi = (x[28] + x[60]) SHIFT_A; /* Re C + Re D */
    ui = (x[29] + x[61]) SHIFT_A; /* Im C + Im D */

    x[24] = vr2 + (vi SHIFT_B); /* Re A' = ReA + ReB +ReC + ReD */
    x[28] = vr2 - (vi SHIFT_B); /* Re C' = -(ReC+ReD) + (ReA+ReB) */
    x[25] = ur2 + (ui SHIFT_B); /* Im A' = sum of imag values */
    x[29] = ur2 - (ui SHIFT_B); /* Im C' = -Im C -Im D +Im A +Im B */

    vr2 -= x[44];            /* Re A - Re B */
    ur2 -= x[45];            /* Im A - Im B */
    vi = (vi SHIFT_B)-x[60]; /* Re C - Re D */
    ui = (ui SHIFT_B)-x[61]; /* Im C - Im D */

    vi2 = (x[26] + x[58]) SHIFT_A; /* Re C + Re D */
    ui2 = (x[27] + x[59]) SHIFT_A; /* Im C + Im D */

    x[26] = ui + vr2; /* Re B' = Im C - Im D  + Re A - Re B */
    x[27] = ur2 - vi; /* Im B'= -Re C + Re D + Im A - Im B */

    vi3 = (x[30] + x[62]) SHIFT_A; /* Re C + Re D */
    ui3 = (x[31] + x[63]) SHIFT_A; /* Im C + Im D */

    x[30] = vr2 - ui; /* Re D' = -Im C + Im D + Re A - Re B */
    x[31] = vi + ur2; /* Im D'= Re C - Re D + Im A - Im B */

    x[40] = vr3 + (vi2 SHIFT_B); /* Re A' = ReA + ReB +ReC + ReD */
    x[44] = vr3 - (vi2 SHIFT_B); /* Re C' = -(ReC+ReD) + (ReA+ReB) */
    x[41] = ur3 + (ui2 SHIFT_B); /* Im A' = sum of imag values */
    x[45] = ur3 - (ui2 SHIFT_B); /* Im C' = -Im C -Im D +Im A +Im B */

    vr3 -= x[42];              /* Re A - Re B */
    ur3 -= x[43];              /* Im A - Im B */
    vi2 = (vi2 SHIFT_B)-x[58]; /* Re C - Re D */
    ui2 = (ui2 SHIFT_B)-x[59]; /* Im C - Im D */

    x[42] = ui2 + vr3; /* Re B' = Im C - Im D  + Re A - Re B */
    x[43] = ur3 - vi2; /* Im B'= -Re C + Re D + Im A - Im B */

    x[56] = vr4 + (vi3 SHIFT_B); /* Re A' = ReA + ReB +ReC + ReD */
    x[60] = vr4 - (vi3 SHIFT_B); /* Re C' = -(ReC+ReD) + (ReA+ReB) */
    x[57] = ur4 + (ui3 SHIFT_B); /* Im A' = sum of imag values */
    x[61] = ur4 - (ui3 SHIFT_B); /* Im C' = -Im C -Im D +Im A +Im B */

    vr4 -= x[46]; /* Re A - Re B */
    ur4 -= x[47]; /* Im A - Im B */

    x[46] = vr3 - ui2; /* Re D' = -Im C + Im D + Re A - Re B */
    x[47] = vi2 + ur3; /* Im D'= Re C - Re D + Im A - Im B */

    vi3 = (vi3 SHIFT_B)-x[62]; /* Re C - Re D */
    ui3 = (ui3 SHIFT_B)-x[63]; /* Im C - Im D */

    x[58] = ui3 + vr4; /* Re B' = Im C - Im D  + Re A - Re B */
    x[62] = vr4 - ui3; /* Re D' = -Im C + Im D + Re A - Re B */
    x[59] = ur4 - vi3; /* Im B'= -Re C + Re D + Im A - Im B */
    x[63] = vi3 + ur4; /* Im D'= Re C - Re D + Im A - Im B */
  }

  {
    FIXP_DBL *xt = _x;

    int j = 4;
    do {
      FIXP_DBL vi, ui, vr, ur;

      vr = xt[8];
      vi = xt[9];
      ur = xt[0] >> 1;
      ui = xt[1] >> 1;
      xt[0] = ur + (vr >> 1);
      xt[1] = ui + (vi >> 1);
      xt[8] = ur - (vr >> 1);
      xt[9] = ui - (vi >> 1);

      vr = xt[13];
      vi = xt[12];
      ur = xt[4] >> 1;
      ui = xt[5] >> 1;
      xt[4] = ur + (vr >> 1);
      xt[5] = ui - (vi >> 1);
      xt[12] = ur - (vr >> 1);
      xt[13] = ui + (vi >> 1);

      SUMDIFF_PIFOURTH(vi, vr, xt[10], xt[11])
      ur = xt[2];
      ui = xt[3];
      xt[2] = (ur >> 1) + vr;
      xt[3] = (ui >> 1) + vi;
      xt[10] = (ur >> 1) - vr;
      xt[11] = (ui >> 1) - vi;

      SUMDIFF_PIFOURTH(vr, vi, xt[14], xt[15])
      ur = xt[6];
      ui = xt[7];

      xt[6] = (ur >> 1) + vr;
      xt[7] = (ui >> 1) - vi;
      xt[14] = (ur >> 1) - vr;
      xt[15] = (ui >> 1) + vi;
      xt += 16;
    } while (--j != 0);
  }

  {
    FIXP_DBL *const x = _x;
    FIXP_DBL vi, ui, vr, ur;

    vr = x[16];
    vi = x[17];
    ur = x[0] >> 1;
    ui = x[1] >> 1;
    x[0] = ur + (vr >> 1);
    x[1] = ui + (vi >> 1);
    x[16] = ur - (vr >> 1);
    x[17] = ui - (vi >> 1);

    vi = x[24];
    vr = x[25];
    ur = x[8] >> 1;
    ui = x[9] >> 1;
    x[8] = ur + (vr >> 1);
    x[9] = ui - (vi >> 1);
    x[24] = ur - (vr >> 1);
    x[25] = ui + (vi >> 1);

    vr = x[48];
    vi = x[49];
    ur = x[32] >> 1;
    ui = x[33] >> 1;
    x[32] = ur + (vr >> 1);
    x[33] = ui + (vi >> 1);
    x[48] = ur - (vr >> 1);
    x[49] = ui - (vi >> 1);

    vi = x[56];
    vr = x[57];
    ur = x[40] >> 1;
    ui = x[41] >> 1;
    x[40] = ur + (vr >> 1);
    x[41] = ui - (vi >> 1);
    x[56] = ur - (vr >> 1);
    x[57] = ui + (vi >> 1);

    cplxMultDiv2(&vi, &vr, x[19], x[18], fft32_w32[0]);
    ur = x[2];
    ui = x[3];
    x[2] = (ur >> 1) + vr;
    x[3] = (ui >> 1) + vi;
    x[18] = (ur >> 1) - vr;
    x[19] = (ui >> 1) - vi;

    cplxMultDiv2(&vr, &vi, x[27], x[26], fft32_w32[0]);
    ur = x[10];
    ui = x[11];
    x[10] = (ur >> 1) + vr;
    x[11] = (ui >> 1) - vi;
    x[26] = (ur >> 1) - vr;
    x[27] = (ui >> 1) + vi;

    cplxMultDiv2(&vi, &vr, x[51], x[50], fft32_w32[0]);
    ur = x[34];
    ui = x[35];
    x[34] = (ur >> 1) + vr;
    x[35] = (ui >> 1) + vi;
    x[50] = (ur >> 1) - vr;
    x[51] = (ui >> 1) - vi;

    cplxMultDiv2(&vr, &vi, x[59], x[58], fft32_w32[0]);
    ur = x[42];
    ui = x[43];
    x[42] = (ur >> 1) + vr;
    x[43] = (ui >> 1) - vi;
    x[58] = (ur >> 1) - vr;
    x[59] = (ui >> 1) + vi;

    SUMDIFF_PIFOURTH(vi, vr, x[20], x[21])
    ur = x[4];
    ui = x[5];
    x[4] = (ur >> 1) + vr;
    x[5] = (ui >> 1) + vi;
    x[20] = (ur >> 1) - vr;
    x[21] = (ui >> 1) - vi;

    SUMDIFF_PIFOURTH(vr, vi, x[28], x[29])
    ur = x[12];
    ui = x[13];
    x[12] = (ur >> 1) + vr;
    x[13] = (ui >> 1) - vi;
    x[28] = (ur >> 1) - vr;
    x[29] = (ui >> 1) + vi;

    SUMDIFF_PIFOURTH(vi, vr, x[52], x[53])
    ur = x[36];
    ui = x[37];
    x[36] = (ur >> 1) + vr;
    x[37] = (ui >> 1) + vi;
    x[52] = (ur >> 1) - vr;
    x[53] = (ui >> 1) - vi;

    SUMDIFF_PIFOURTH(vr, vi, x[60], x[61])
    ur = x[44];
    ui = x[45];
    x[44] = (ur >> 1) + vr;
    x[45] = (ui >> 1) - vi;
    x[60] = (ur >> 1) - vr;
    x[61] = (ui >> 1) + vi;

    cplxMultDiv2(&vi, &vr, x[23], x[22], fft32_w32[1]);
    ur = x[6];
    ui = x[7];
    x[6] = (ur >> 1) + vr;
    x[7] = (ui >> 1) + vi;
    x[22] = (ur >> 1) - vr;
    x[23] = (ui >> 1) - vi;

    cplxMultDiv2(&vr, &vi, x[31], x[30], fft32_w32[1]);
    ur = x[14];
    ui = x[15];
    x[14] = (ur >> 1) + vr;
    x[15] = (ui >> 1) - vi;
    x[30] = (ur >> 1) - vr;
    x[31] = (ui >> 1) + vi;

    cplxMultDiv2(&vi, &vr, x[55], x[54], fft32_w32[1]);
    ur = x[38];
    ui = x[39];
    x[38] = (ur >> 1) + vr;
    x[39] = (ui >> 1) + vi;
    x[54] = (ur >> 1) - vr;
    x[55] = (ui >> 1) - vi;

    cplxMultDiv2(&vr, &vi, x[63], x[62], fft32_w32[1]);
    ur = x[46];
    ui = x[47];

    x[46] = (ur >> 1) + vr;
    x[47] = (ui >> 1) - vi;
    x[62] = (ur >> 1) - vr;
    x[63] = (ui >> 1) + vi;

    vr = x[32];
    vi = x[33];
    ur = x[0] >> 1;
    ui = x[1] >> 1;
    x[0] = ur + (vr >> 1);
    x[1] = ui + (vi >> 1);
    x[32] = ur - (vr >> 1);
    x[33] = ui - (vi >> 1);

    vi = x[48];
    vr = x[49];
    ur = x[16] >> 1;
    ui = x[17] >> 1;
    x[16] = ur + (vr >> 1);
    x[17] = ui - (vi >> 1);
    x[48] = ur - (vr >> 1);
    x[49] = ui + (vi >> 1);

    cplxMultDiv2(&vi, &vr, x[35], x[34], fft32_w32[2]);
    ur = x[2];
    ui = x[3];
    x[2] = (ur >> 1) + vr;
    x[3] = (ui >> 1) + vi;
    x[34] = (ur >> 1) - vr;
    x[35] = (ui >> 1) - vi;

    cplxMultDiv2(&vr, &vi, x[51], x[50], fft32_w32[2]);
    ur = x[18];
    ui = x[19];
    x[18] = (ur >> 1) + vr;
    x[19] = (ui >> 1) - vi;
    x[50] = (ur >> 1) - vr;
    x[51] = (ui >> 1) + vi;

    cplxMultDiv2(&vi, &vr, x[37], x[36], fft32_w32[0]);
    ur = x[4];
    ui = x[5];
    x[4] = (ur >> 1) + vr;
    x[5] = (ui >> 1) + vi;
    x[36] = (ur >> 1) - vr;
    x[37] = (ui >> 1) - vi;

    cplxMultDiv2(&vr, &vi, x[53], x[52], fft32_w32[0]);
    ur = x[20];
    ui = x[21];
    x[20] = (ur >> 1) + vr;
    x[21] = (ui >> 1) - vi;
    x[52] = (ur >> 1) - vr;
    x[53] = (ui >> 1) + vi;

    cplxMultDiv2(&vi, &vr, x[39], x[38], fft32_w32[3]);
    ur = x[6];
    ui = x[7];
    x[6] = (ur >> 1) + vr;
    x[7] = (ui >> 1) + vi;
    x[38] = (ur >> 1) - vr;
    x[39] = (ui >> 1) - vi;

    cplxMultDiv2(&vr, &vi, x[55], x[54], fft32_w32[3]);
    ur = x[22];
    ui = x[23];
    x[22] = (ur >> 1) + vr;
    x[23] = (ui >> 1) - vi;
    x[54] = (ur >> 1) - vr;
    x[55] = (ui >> 1) + vi;

    SUMDIFF_PIFOURTH(vi, vr, x[40], x[41])
    ur = x[8];
    ui = x[9];
    x[8] = (ur >> 1) + vr;
    x[9] = (ui >> 1) + vi;
    x[40] = (ur >> 1) - vr;
    x[41] = (ui >> 1) - vi;

    SUMDIFF_PIFOURTH(vr, vi, x[56], x[57])
    ur = x[24];
    ui = x[25];
    x[24] = (ur >> 1) + vr;
    x[25] = (ui >> 1) - vi;
    x[56] = (ur >> 1) - vr;
    x[57] = (ui >> 1) + vi;

    cplxMultDiv2(&vi, &vr, x[43], x[42], fft32_w32[4]);
    ur = x[10];
    ui = x[11];

    x[10] = (ur >> 1) + vr;
    x[11] = (ui >> 1) + vi;
    x[42] = (ur >> 1) - vr;
    x[43] = (ui >> 1) - vi;

    cplxMultDiv2(&vr, &vi, x[59], x[58], fft32_w32[4]);
    ur = x[26];
    ui = x[27];
    x[26] = (ur >> 1) + vr;
    x[27] = (ui >> 1) - vi;
    x[58] = (ur >> 1) - vr;
    x[59] = (ui >> 1) + vi;

    cplxMultDiv2(&vi, &vr, x[45], x[44], fft32_w32[1]);
    ur = x[12];
    ui = x[13];
    x[12] = (ur >> 1) + vr;
    x[13] = (ui >> 1) + vi;
    x[44] = (ur >> 1) - vr;
    x[45] = (ui >> 1) - vi;

    cplxMultDiv2(&vr, &vi, x[61], x[60], fft32_w32[1]);
    ur = x[28];
    ui = x[29];
    x[28] = (ur >> 1) + vr;
    x[29] = (ui >> 1) - vi;
    x[60] = (ur >> 1) - vr;
    x[61] = (ui >> 1) + vi;

    cplxMultDiv2(&vi, &vr, x[47], x[46], fft32_w32[5]);
    ur = x[14];
    ui = x[15];
    x[14] = (ur >> 1) + vr;
    x[15] = (ui >> 1) + vi;
    x[46] = (ur >> 1) - vr;
    x[47] = (ui >> 1) - vi;

    cplxMultDiv2(&vr, &vi, x[63], x[62], fft32_w32[5]);
    ur = x[30];
    ui = x[31];
    x[30] = (ur >> 1) + vr;
    x[31] = (ui >> 1) - vi;
    x[62] = (ur >> 1) - vr;
    x[63] = (ui >> 1) + vi;
  }
}
#endif

#define noFFT_APPLY_ROT_VECTOR_HQ

#ifndef FUNCTION_fft_apply_rot_vector__FIXP_DBL
static inline void fft_apply_rot_vector(FIXP_DBL *RESTRICT pData, const int cl,
                                        const int l, const FIXP_STB *pVecRe,
                                        const FIXP_STB *pVecIm) {
  FIXP_DBL re, im;
  FIXP_STB vre, vim;

  int i, c;

  for (i = 0; i < cl; i++) {
    re = pData[2 * i];
    im = pData[2 * i + 1];

    pData[2 * i] = re >> 2;     /* * 0.25 */
    pData[2 * i + 1] = im >> 2; /* * 0.25 */
  }
  for (; i < l; i += cl) {
    re = pData[2 * i];
    im = pData[2 * i + 1];

    pData[2 * i] = re >> 2;     /* * 0.25 */
    pData[2 * i + 1] = im >> 2; /* * 0.25 */

    for (c = i + 1; c < i + cl; c++) {
      re = pData[2 * c] >> 1;
      im = pData[2 * c + 1] >> 1;
      vre = *pVecRe++;
      vim = *pVecIm++;

      cplxMultDiv2(&pData[2 * c + 1], &pData[2 * c], im, re, vre, vim);
    }
  }
}
#endif

#ifndef FUNCTION_fftN2_func
static inline void fftN2_func(FIXP_DBL *pInput, const int length,
                              const int dim1, const int dim2,
                              void (*const fft1)(FIXP_DBL *),
                              void (*const fft2)(FIXP_DBL *),
                              const FIXP_STB *RotVectorReal,
                              const FIXP_STB *RotVectorImag, FIXP_DBL *aDst,
                              FIXP_DBL *aDst2) {
  FIXP_DBL *pSrc, *pDst, *pDstOut;
  int i;

  FDK_ASSERT(length == dim1 * dim2);

  pSrc = pInput;
  pDst = aDst;
  for (i = 0; i < dim2; i++) {
    for (int j = 0; j < dim1; j++) {
      pDst[2 * j] = pSrc[2 * j * dim2];
      pDst[2 * j + 1] = pSrc[2 * j * dim2 + 1];
    }

#ifndef FFT_TWO_STAGE_SWITCH_CASE
    fft1(pDst);
#else
    switch (dim1) {
      case 2:
        fft2(pDst);
        break;
      case 3:
        fft3(pDst);
        break;
      case 4:
        fft_4(pDst);
        break;
      case 12:
        fft12(pDst);
        break;
      case 16:
        fft_16(pDst);
        break;
      case 32:
        fft_32(pDst);
        break;
    }
#endif
    pSrc += 2;
    pDst = pDst + 2 * dim1;
  }

  pSrc = aDst;
  fft_apply_rot_vector(pSrc, dim1, length, RotVectorReal, RotVectorImag);

  pSrc = aDst;
  pDst = aDst2;
  pDstOut = pInput;
  for (i = 0; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      pDst[2 * j] = pSrc[2 * j * dim1];
      pDst[2 * j + 1] = pSrc[2 * j * dim1 + 1];
    }

#ifndef FFT_TWO_STAGE_SWITCH_CASE
    fft2(pDst);
#else
    switch (dim2) {
      case 4:
        fft_4(pDst);
        break;
      case 9:
        fft9(pDst);
        break;
      case 12:
        fft12(pDst);
        break;
      case 15:
        fft15(pDst);
        break;
      case 16:
        fft_16(pDst);
        break;
      case 32:
        fft_32(pDst);
        break;
    }
#endif

    for (int j = 0; j < dim2; j++) {
      pDstOut[2 * j * dim1] = pDst[2 * j];
      pDstOut[2 * j * dim1 + 1] = pDst[2 * j + 1];
    }
    pSrc += 2;
    pDstOut += 2;
  }
}
#endif

#define fftN2(DATA_TYPE, pInput, length, dim1, dim2, fft_func1, fft_func2, \
              RotVectorReal, RotVectorImag)                                \
  {                                                                        \
    C_AALLOC_SCRATCH_START(aDst, DATA_TYPE, 2 * length)                    \
    C_AALLOC_SCRATCH_START(aDst2, DATA_TYPE, 2 * dim2)                     \
    fftN2_func(pInput, length, dim1, dim2, fft_func1, fft_func2,           \
               RotVectorReal, RotVectorImag, aDst, aDst2);                 \
    C_AALLOC_SCRATCH_END(aDst2, DATA_TYPE, 2 * dim2)                       \
    C_AALLOC_SCRATCH_END(aDst, DATA_TYPE, 2 * length)                      \
  }

#ifndef FUNCTION_fft6
static inline void fft6(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 6, 2, 3, fft2, fft3, RotVectorReal6, RotVectorImag6);
}
#endif

#ifndef FUNCTION_fft12
static inline void fft12(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 12, 3, 4, fft3, fft_4, RotVectorReal12,
        RotVectorImag12);
}
#endif

#ifndef FUNCTION_fft20
static inline void fft20(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 20, 4, 5, fft_4, fft5, RotVectorReal20,
        RotVectorImag20);
}
#endif

static inline void fft24(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 24, 2, 12, fft2, fft12, RotVectorReal24,
        RotVectorImag24);
}

static inline void fft48(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 48, 4, 12, fft_4, fft12, RotVectorReal48,
        RotVectorImag48);
}

#ifndef FUNCTION_fft60
static inline void fft60(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 60, 4, 15, fft_4, fft15, RotVectorReal60,
        RotVectorImag60);
}
#endif

#ifndef FUNCTION_fft80
static inline void fft80(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 80, 5, 16, fft5, fft_16, RotVectorReal80,
        RotVectorImag80);
}
#endif

#ifndef FUNCTION_fft96
static inline void fft96(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 96, 3, 32, fft3, fft_32, RotVectorReal96,
        RotVectorImag96);
}
#endif

#ifndef FUNCTION_fft120
static inline void fft120(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 120, 8, 15, fft_8, fft15, RotVectorReal120,
        RotVectorImag120);
}
#endif

#ifndef FUNCTION_fft192
static inline void fft192(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 192, 16, 12, fft_16, fft12, RotVectorReal192,
        RotVectorImag192);
}
#endif

#ifndef FUNCTION_fft240
static inline void fft240(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 240, 16, 15, fft_16, fft15, RotVectorReal240,
        RotVectorImag240);
}
#endif

#ifndef FUNCTION_fft384
static inline void fft384(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 384, 12, 32, fft12, fft_32, RotVectorReal384,
        RotVectorImag384);
}
#endif

#ifndef FUNCTION_fft480
static inline void fft480(FIXP_DBL *pInput) {
  fftN2(FIXP_DBL, pInput, 480, 32, 15, fft_32, fft15, RotVectorReal480,
        RotVectorImag480);
}
#endif

void fft(int length, FIXP_DBL *pInput, INT *pScalefactor) {
  C_ALLOC_ALIGNED_CHECK(pInput);

  if (length == 32) {
    fft_32(pInput);
    *pScalefactor += SCALEFACTOR32;
  } else {
    switch (length) {
      case 16:
        fft_16(pInput);
        *pScalefactor += SCALEFACTOR16;
        break;
      case 8:
        fft_8(pInput);
        *pScalefactor += SCALEFACTOR8;
        break;
      case 2:
        fft2(pInput);
        *pScalefactor += SCALEFACTOR2;
        break;
      case 3:
        fft3(pInput);
        *pScalefactor += SCALEFACTOR3;
        break;
      case 4:
        fft_4(pInput);
        *pScalefactor += SCALEFACTOR4;
        break;
      case 5:
        fft5(pInput);
        *pScalefactor += SCALEFACTOR5;
        break;
      case 6:
        fft6(pInput);
        *pScalefactor += SCALEFACTOR6;
        break;
      case 10:
        fft10(pInput);
        *pScalefactor += SCALEFACTOR10;
        break;
      case 12:
        fft12(pInput);
        *pScalefactor += SCALEFACTOR12;
        break;
      case 15:
        fft15(pInput);
        *pScalefactor += SCALEFACTOR15;
        break;
      case 20:
        fft20(pInput);
        *pScalefactor += SCALEFACTOR20;
        break;
      case 24:
        fft24(pInput);
        *pScalefactor += SCALEFACTOR24;
        break;
      case 48:
        fft48(pInput);
        *pScalefactor += SCALEFACTOR48;
        break;
      case 60:
        fft60(pInput);
        *pScalefactor += SCALEFACTOR60;
        break;
      case 64:
        dit_fft(pInput, 6, SineTable512, 512);
        *pScalefactor += SCALEFACTOR64;
        break;
      case 80:
        fft80(pInput);
        *pScalefactor += SCALEFACTOR80;
        break;
      case 96:
        fft96(pInput);
        *pScalefactor += SCALEFACTOR96;
        break;
      case 120:
        fft120(pInput);
        *pScalefactor += SCALEFACTOR120;
        break;
      case 128:
        dit_fft(pInput, 7, SineTable512, 512);
        *pScalefactor += SCALEFACTOR128;
        break;
      case 192:
        fft192(pInput);
        *pScalefactor += SCALEFACTOR192;
        break;
      case 240:
        fft240(pInput);
        *pScalefactor += SCALEFACTOR240;
        break;
      case 256:
        dit_fft(pInput, 8, SineTable512, 512);
        *pScalefactor += SCALEFACTOR256;
        break;
      case 384:
        fft384(pInput);
        *pScalefactor += SCALEFACTOR384;
        break;
      case 480:
        fft480(pInput);
        *pScalefactor += SCALEFACTOR480;
        break;
      case 512:
        dit_fft(pInput, 9, SineTable512, 512);
        *pScalefactor += SCALEFACTOR512;
        break;
      default:
        FDK_ASSERT(0);
        break;
    }
  }
}

void ifft(int length, FIXP_DBL *pInput, INT *scalefactor) {
  switch (length) {
    default:
      FDK_ASSERT(0);
      break;
  }
}
