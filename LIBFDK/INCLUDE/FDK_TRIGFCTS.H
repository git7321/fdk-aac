#if !defined(FDK_TRIGFCTS_H)
#define FDK_TRIGFCTS_H

#include "common_fix.h"
#include "FDK_tools_rom.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif
#define Q(format) ((FIXP_DBL)(((LONG)1) << (format)))

#define Q_ATANINP (25)
#define Q_ATANOUT (30)
#define ATI_SF ((DFRACT_BITS - 1) - Q_ATANINP)
#define ATI_SCALE ((float)(1 << ATI_SF))
#define ATO_SF ((DFRACT_BITS - 1) - Q_ATANOUT)
#define ATO_SCALE ((float)(1 << ATO_SF))
#define Q_ATAN2OUT (29)
#define AT2O_SF ((DFRACT_BITS - 1) - Q_ATAN2OUT)
#define AT2O_SCALE ((float)(1 << AT2O_SF))

FIXP_DBL fixp_atan(FIXP_DBL x);
FIXP_DBL fixp_atan2(FIXP_DBL y, FIXP_DBL x);

FIXP_DBL fixp_cos(FIXP_DBL x, int scale);
FIXP_DBL fixp_sin(FIXP_DBL x, int scale);

#define FIXP_COS_SIN

#include "FDK_tools_rom.h"

#define SINETAB SineTable512
#define LD 9

#ifndef FUNCTION_inline_fixp_cos_sin
#define FUNCTION_inline_fixp_cos_sin

static inline FIXP_DBL fixp_sin_cos_residual_inline(FIXP_DBL x, int scale,
                                                    FIXP_DBL *sine,
                                                    FIXP_DBL *cosine) {
  FIXP_DBL residual;
  int s;
  int shift = (31 - scale - LD - 1);
  int ssign = 1;
  int csign = 1;

  residual = fMult(x, FL2FXCONST_DBL(1.0 / M_PI));
  s = ((LONG)residual) >> shift;

  residual &= ((1 << shift) - 1);
  residual = fMult(residual, FL2FXCONST_DBL(M_PI / 4.0)) << 2;
  residual <<= scale;

  if (s & ((1 << LD) << 1)) {
    ssign = -ssign;
  }
  if ((s + (1 << LD)) & ((1 << LD) << 1)) {
    csign = -csign;
  }

  s = fAbs(s);

  s &= (((1 << LD) << 1) - 1);

  if (s > (1 << LD)) {
    s = ((1 << LD) << 1) - s;
  }

  {
    LONG sl, cl;
    if (s > (1 << (LD - 1))) {
      FIXP_STP tmp;
      s = (1 << LD) - s;
      tmp = SINETAB[s];
      sl = (LONG)tmp.v.re;
      cl = (LONG)tmp.v.im;
    } else {
      FIXP_STP tmp;
      tmp = SINETAB[s];
      sl = (LONG)tmp.v.im;
      cl = (LONG)tmp.v.re;
    }

#ifdef SINETABLE_16BIT
    *sine = (FIXP_DBL)((sl * ssign) << (DFRACT_BITS - FRACT_BITS));
    *cosine = (FIXP_DBL)((cl * csign) << (DFRACT_BITS - FRACT_BITS));
#else
    *sine = (FIXP_DBL)(sl * ssign) >> 1;
    *cosine = (FIXP_DBL)(cl * csign) >> 1;
#endif
  }

  return residual;
}

static inline void inline_fixp_cos_sin(FIXP_DBL x1, FIXP_DBL x2,
                                       const int scale, FIXP_DBL *out) {
  FIXP_DBL residual, error0, error1, sine, cosine;
  residual = fixp_sin_cos_residual_inline(x1, scale, &sine, &cosine);
  error0 = fMultDiv2(sine, residual);
  error1 = fMultDiv2(cosine, residual);

#ifdef SINETABLE_16BIT
  *out++ = cosine - (error0 << 1);
  *out++ = sine + (error1 << 1);
#else
  *out++ = SATURATE_LEFT_SHIFT(cosine - (error0 << 1), 1, DFRACT_BITS);
  *out++ = SATURATE_LEFT_SHIFT(sine + (error1 << 1), 1, DFRACT_BITS);
#endif

  residual = fixp_sin_cos_residual_inline(x2, scale, &sine, &cosine);
  error0 = fMultDiv2(sine, residual);
  error1 = fMultDiv2(cosine, residual);

#ifdef SINETABLE_16BIT
  *out++ = cosine - (error0 << 1);
  *out++ = sine + (error1 << 1);
#else
  *out++ = SATURATE_LEFT_SHIFT(cosine - (error0 << 1), 1, DFRACT_BITS);
  *out++ = SATURATE_LEFT_SHIFT(sine + (error1 << 1), 1, DFRACT_BITS);
#endif
}
#endif

#endif /* !defined(FDK_TRIGFCTS_H) */
