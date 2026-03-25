#ifndef FIXPOINT_MATH_H
#define FIXPOINT_MATH_H

#include "common_fix.h"
#include "scale.h"

#define LD_DATA_SCALING (64.0f)
#define LD_DATA_SHIFT 6

#define MAX_LD_PRECISION 10
#define LD_PRECISION 10

#ifndef LDCOEFF_16BIT
LNK_SECTION_CONSTDATA_L1
static const FIXP_DBL ldCoeff[MAX_LD_PRECISION] = {
    FL2FXCONST_DBL(-1.0),       FL2FXCONST_DBL(-1.0 / 2.0),
    FL2FXCONST_DBL(-1.0 / 3.0), FL2FXCONST_DBL(-1.0 / 4.0),
    FL2FXCONST_DBL(-1.0 / 5.0), FL2FXCONST_DBL(-1.0 / 6.0),
    FL2FXCONST_DBL(-1.0 / 7.0), FL2FXCONST_DBL(-1.0 / 8.0),
    FL2FXCONST_DBL(-1.0 / 9.0), FL2FXCONST_DBL(-1.0 / 10.0)};
#else
LNK_SECTION_CONSTDATA_L1
static const FIXP_SGL ldCoeff[MAX_LD_PRECISION] = {
    FL2FXCONST_SGL(-1.0),       FL2FXCONST_SGL(-1.0 / 2.0),
    FL2FXCONST_SGL(-1.0 / 3.0), FL2FXCONST_SGL(-1.0 / 4.0),
    FL2FXCONST_SGL(-1.0 / 5.0), FL2FXCONST_SGL(-1.0 / 6.0),
    FL2FXCONST_SGL(-1.0 / 7.0), FL2FXCONST_SGL(-1.0 / 8.0),
    FL2FXCONST_SGL(-1.0 / 9.0), FL2FXCONST_SGL(-1.0 / 10.0)};
#endif

#define SQRT_BITS 7
#define SQRT_VALUES (128 + 2)
#define SQRT_BITS_MASK 0x7f
#define SQRT_FRACT_BITS_MASK 0x007FFFFF

extern const FIXP_DBL invSqrtTab[SQRT_VALUES];

#if defined(__x86__)
#include "x86/fixpoint_math_x86.h"
#endif

#if !defined(FUNCTION_fIsLessThan)
FDK_INLINE INT fIsLessThan(FIXP_DBL a_m, INT a_e, FIXP_DBL b_m, INT b_e) {
  INT n;

  n = fixnorm_D(a_m);
  a_m <<= n;
  a_e -= n;

  n = fixnorm_D(b_m);
  b_m <<= n;
  b_e -= n;

  if (a_m == (FIXP_DBL)0) a_e = b_e;
  if (b_m == (FIXP_DBL)0) b_e = a_e;

  if (a_e > b_e) {
    return ((b_m >> fMin(a_e - b_e, DFRACT_BITS - 1)) > a_m);
  } else {
    return ((a_m >> fMin(b_e - a_e, DFRACT_BITS - 1)) < b_m);
  }
}

FDK_INLINE INT fIsLessThan(FIXP_SGL a_m, INT a_e, FIXP_SGL b_m, INT b_e) {
  INT n;

  n = fixnorm_S(a_m);
  a_m <<= n;
  a_e -= n;

  n = fixnorm_S(b_m);
  b_m <<= n;
  b_e -= n;

  if (a_m == (FIXP_SGL)0) a_e = b_e;
  if (b_m == (FIXP_SGL)0) b_e = a_e;

  if (a_e > b_e) {
    return ((b_m >> fMin(a_e - b_e, FRACT_BITS - 1)) > a_m);
  } else {
    return ((a_m >> fMin(b_e - a_e, FRACT_BITS - 1)) < b_m);
  }
}
#endif

#define CalcLdData(op) fLog2(op, 0)

void LdDataVector(FIXP_DBL *srcVector, FIXP_DBL *destVector, INT number);

extern const UINT exp2_tab_long[32];
extern const UINT exp2w_tab_long[32];
extern const UINT exp2x_tab_long[32];

LNK_SECTION_CODE_L1
FDK_INLINE FIXP_DBL CalcInvLdData(const FIXP_DBL x) {
  int set_zero = (x < FL2FXCONST_DBL(-31.0 / 64.0)) ? 0 : 1;
  int set_max = (x >= FL2FXCONST_DBL(31.0 / 64.0)) | (x == FL2FXCONST_DBL(0.0));

  FIXP_SGL frac = (FIXP_SGL)((LONG)x & 0x3FF);
  UINT index3 = (UINT)(LONG)(x >> 10) & 0x1F;
  UINT index2 = (UINT)(LONG)(x >> 15) & 0x1F;
  UINT index1 = (UINT)(LONG)(x >> 20) & 0x1F;
  int exp = fMin(31, ((x > FL2FXCONST_DBL(0.0f)) ? (31 - (int)(x >> 25))
                                                 : (int)(-(x >> 25))));

  UINT lookup1 = exp2_tab_long[index1] * set_zero;
  UINT lookup2 = exp2w_tab_long[index2];
  UINT lookup3 = exp2x_tab_long[index3];
  UINT lookup3f =
      lookup3 + (UINT)(LONG)fMultDiv2((FIXP_DBL)(0x0016302F), (FIXP_SGL)frac);

  UINT lookup12 = (UINT)(LONG)fMult((FIXP_DBL)lookup1, (FIXP_DBL)lookup2);
  UINT lookup = (UINT)(LONG)fMult((FIXP_DBL)lookup12, (FIXP_DBL)lookup3f);

  FIXP_DBL retVal = (lookup << 3) >> exp;

  if (set_max) {
    retVal = (FIXP_DBL)MAXVAL_DBL;
  }

  return retVal;
}

void InitLdInt();
FIXP_DBL CalcLdInt(INT i);

extern const USHORT sqrt_tab[49];

inline FIXP_DBL sqrtFixp_lookup(FIXP_DBL x) {
  UINT y = (INT)x;
  UCHAR is_zero = (y == 0);
  INT zeros = fixnormz_D(y) & 0x1e;
  y <<= zeros;
  UINT idx = (y >> 26) - 16;
  USHORT frac = (y >> 10) & 0xffff;
  USHORT nfrac = 0xffff ^ frac;
  UINT t = (UINT)nfrac * sqrt_tab[idx] + (UINT)frac * sqrt_tab[idx + 1];
  t = t >> (zeros >> 1);
  return (is_zero ? 0 : t);
}

inline FIXP_DBL sqrtFixp_lookup(FIXP_DBL x, INT *x_e) {
  UINT y = (INT)x;
  INT e;

  if (x == (FIXP_DBL)0) {
    return x;
  }

  e = fixnormz_D(y);
  y <<= e;
  e = *x_e - e + 2;

  if (e & 1) {
    y >>= 1;
    e++;
  }
  UINT idx = (y >> 26) - 16;
  USHORT frac = (y >> 10) & 0xffff;
  USHORT nfrac = 0xffff ^ frac;
  UINT t = (UINT)nfrac * sqrt_tab[idx] + (UINT)frac * sqrt_tab[idx + 1];

  *x_e = e >> 1;
  return (FIXP_DBL)(LONG)(t >> 1);
}

void InitInvSqrtTab();

#ifndef FUNCTION_invSqrtNorm2
static FDK_FORCEINLINE FIXP_DBL invSqrtNorm2(FIXP_DBL op, INT *shift) {
  FIXP_DBL val = op;
  FIXP_DBL reg1, reg2;

  if (val == FL2FXCONST_DBL(0.0)) {
    *shift = 16;
    return ((LONG)MAXVAL_DBL);
  }

#define INVSQRTNORM2_LINEAR_INTERPOLATE
#define INVSQRTNORM2_LINEAR_INTERPOLATE_HQ

  FDK_ASSERT(val > FL2FXCONST_DBL(0.0));
  *shift = fNormz(val) - 1;
  val <<= *shift;
  *shift += 2;

#if defined(INVSQRTNORM2_LINEAR_INTERPOLATE)
  INT index =
      (INT)(val >> (DFRACT_BITS - 1 - (SQRT_BITS + 1))) & SQRT_BITS_MASK;
  FIXP_DBL Fract =
      (FIXP_DBL)(((INT)val & SQRT_FRACT_BITS_MASK) << (SQRT_BITS + 1));
  FIXP_DBL diff = invSqrtTab[index + 1] - invSqrtTab[index];
  reg1 = invSqrtTab[index] + (fMultDiv2(diff, Fract) << 1);
#if defined(INVSQRTNORM2_LINEAR_INTERPOLATE_HQ)
  if (Fract != (FIXP_DBL)0) {
    Fract = fMultDiv2(Fract, (FIXP_DBL)((ULONG)0x80000000 - (ULONG)Fract)) << 1;
    diff = diff - (invSqrtTab[index + 2] - invSqrtTab[index + 1]);
    reg1 = fMultAddDiv2(reg1, Fract, diff);
  }
#endif
#else
#error \
    "Either define INVSQRTNORM2_NEWTON_ITERATE or INVSQRTNORM2_LINEAR_INTERPOLATE"
#endif
  if (*shift & 0x00000001) {
    reg2 = (FIXP_DBL)0x5A827999;
    reg1 = fMultDiv2(reg1, reg2) << 2;
  }

  *shift = *shift >> 1;

  return (reg1);
}
#endif

#ifndef FUNCTION_sqrtFixp
static FDK_FORCEINLINE FIXP_DBL sqrtFixp(FIXP_DBL op) {
  INT tmp_exp = 0;
  FIXP_DBL tmp_inv = invSqrtNorm2(op, &tmp_exp);

  FDK_ASSERT(tmp_exp > 0);
  return ((FIXP_DBL)(fMultDiv2((op << (tmp_exp - 1)), tmp_inv) << 2));
}
#endif

#ifndef FUNCTION_invFixp
static inline FIXP_DBL invFixp(FIXP_DBL op) {
  if ((op == (FIXP_DBL)0x00000000) || (op == (FIXP_DBL)0x00000001)) {
    return ((LONG)MAXVAL_DBL);
  }
  INT tmp_exp;
  FIXP_DBL tmp_inv = invSqrtNorm2(op, &tmp_exp);
  FDK_ASSERT((31 - (2 * tmp_exp + 1)) >= 0);
  int shift = 31 - (2 * tmp_exp + 1);
  tmp_inv = fPow2Div2(tmp_inv);
  if (shift) {
    tmp_inv = ((tmp_inv >> (shift - 1)) + (FIXP_DBL)1) >> 1;
  }
  return tmp_inv;
}

static inline FIXP_DBL invFixp(FIXP_DBL op_m, int *op_e) {
  if ((op_m == (FIXP_DBL)0x00000000) || (op_m == (FIXP_DBL)0x00000001)) {
    *op_e = 31 - *op_e;
    return ((LONG)MAXVAL_DBL);
  }

  INT tmp_exp;
  FIXP_DBL tmp_inv = invSqrtNorm2(op_m, &tmp_exp);

  *op_e = (tmp_exp << 1) - *op_e + 1;
  return fPow2Div2(tmp_inv);
}
#endif

#ifndef FUNCTION_schur_div

FIXP_DBL schur_div(FIXP_DBL num, FIXP_DBL denum, INT count);

#endif

FIXP_DBL mul_dbl_sgl_rnd(const FIXP_DBL op1, const FIXP_SGL op2);

#ifndef FUNCTION_fMultNorm
FIXP_DBL fMultNorm(FIXP_DBL f1, FIXP_DBL f2, INT *result_e);

inline FIXP_DBL fMultNorm(FIXP_DBL f1, FIXP_DBL f2) {
  FIXP_DBL m;
  INT e;

  m = fMultNorm(f1, f2, &e);

  m = scaleValueSaturate(m, e);

  return m;
}

inline FIXP_DBL fMultNorm(FIXP_DBL f1_m, INT f1_e, FIXP_DBL f2_m, INT f2_e,
                          INT result_e) {
  FIXP_DBL m;
  INT e;

  m = fMultNorm(f1_m, f2_m, &e);
  m = scaleValueSaturate(m, e + f1_e + f2_e - result_e);

  return m;
}
#endif

#ifndef FUNCTION_fMultI
inline INT fMultI(FIXP_DBL a, INT b) {
  FIXP_DBL m, mi;
  INT m_e;

  m = fMultNorm(a, (FIXP_DBL)b, &m_e);

  if (m_e < (INT)0) {
    if (m_e > (INT)-DFRACT_BITS) {
      m = m >> ((-m_e) - 1);
      mi = (m + (FIXP_DBL)1) >> 1;
    } else {
      mi = (FIXP_DBL)0;
    }
  } else {
    mi = scaleValueSaturate(m, m_e);
  }

  return ((INT)mi);
}
#endif

#ifndef FUNCTION_fMultIfloor
inline INT fMultIfloor(FIXP_DBL a, INT b) {
  FIXP_DBL m, mi;
  INT m_e;

  m = fMultNorm(a, (FIXP_DBL)b, &m_e);

  if (m_e < (INT)0) {
    if (m_e > (INT)-DFRACT_BITS) {
      mi = m >> (-m_e);
    } else {
      mi = (FIXP_DBL)0;
      if (m < (FIXP_DBL)0) {
        mi = (FIXP_DBL)-1;
      }
    }
  } else {
    mi = scaleValueSaturate(m, m_e);
  }

  return ((INT)mi);
}
#endif

#ifndef FUNCTION_fMultIceil
inline INT fMultIceil(FIXP_DBL a, INT b) {
  FIXP_DBL m, mi;
  INT m_e;

  m = fMultNorm(a, (FIXP_DBL)b, &m_e);

  if (m_e < (INT)0) {
    if (m_e > (INT) - (DFRACT_BITS - 1)) {
      mi = (m >> (-m_e));
      if ((LONG)m & ((1 << (-m_e)) - 1)) {
        mi = mi + (FIXP_DBL)1;
      }
    } else {
      if (m > (FIXP_DBL)0) {
        mi = (FIXP_DBL)1;
      } else {
        if ((m_e == -(DFRACT_BITS - 1)) && (m == (FIXP_DBL)MINVAL_DBL)) {
          mi = (FIXP_DBL)-1;
        } else {
          mi = (FIXP_DBL)0;
        }
      }
    }
  } else {
    mi = scaleValueSaturate(m, m_e);
  }

  return ((INT)mi);
}
#endif

#ifndef FUNCTION_fDivNorm
FIXP_DBL fDivNorm(FIXP_DBL num, FIXP_DBL denom, INT *result_e);
FIXP_DBL fDivNorm(FIXP_DBL num, FIXP_DBL denom);
FIXP_DBL fDivNormSigned(FIXP_DBL L_num, FIXP_DBL L_denum, INT *result_e);
FIXP_DBL fDivNormSigned(FIXP_DBL num, FIXP_DBL denom);
#endif

inline FIXP_DBL fAdjust(FIXP_DBL a_m, INT *pA_e) {
  INT shift;

  shift = fNorm(a_m) - 1;
  *pA_e -= shift;

  return scaleValue(a_m, shift);
}

#ifndef FUNCTION_fAddNorm
inline FIXP_DBL fAddNorm(FIXP_DBL a_m, INT a_e, FIXP_DBL b_m, INT b_e,
                         INT *pResult_e) {
  INT result_e;
  FIXP_DBL result_m;

  if (a_m == (FIXP_DBL)0) {
    *pResult_e = b_e;
    return b_m;
  }
  if (b_m == (FIXP_DBL)0) {
    *pResult_e = a_e;
    return a_m;
  }

  a_m = fAdjust(a_m, &a_e);
  b_m = fAdjust(b_m, &b_e);

  if (a_e > b_e) {
    result_m = a_m + (b_m >> fMin(a_e - b_e, DFRACT_BITS - 1));
    result_e = a_e;
  } else {
    result_m = (a_m >> fMin(b_e - a_e, DFRACT_BITS - 1)) + b_m;
    result_e = b_e;
  }

  *pResult_e = result_e;
  return result_m;
}

inline FIXP_DBL fAddNorm(FIXP_DBL a_m, INT a_e, FIXP_DBL b_m, INT b_e,
                         INT result_e) {
  FIXP_DBL result_m;

  a_m = scaleValue(a_m, a_e - result_e);
  b_m = scaleValue(b_m, b_e - result_e);

  result_m = a_m + b_m;

  return result_m;
}
#endif

FIXP_DBL fDivNormHighPrec(FIXP_DBL L_num, FIXP_DBL L_denum, INT *result_e);

#ifndef FUNCTION_fPow
FIXP_DBL f2Pow(const FIXP_DBL exp_m, const INT exp_e, INT *result_e);
FIXP_DBL f2Pow(const FIXP_DBL exp_m, const INT exp_e);
FIXP_DBL fLdPow(FIXP_DBL baseLd_m, INT baseLd_e, FIXP_DBL exp_m, INT exp_e,
                INT *result_e);
FIXP_DBL fLdPow(FIXP_DBL baseLd_m, INT baseLd_e, FIXP_DBL exp_m, INT exp_e);
FIXP_DBL fPow(FIXP_DBL base_m, INT base_e, FIXP_DBL exp_m, INT exp_e,
              INT *result_e);
FIXP_DBL fPowInt(FIXP_DBL base_m, INT base_e, INT N, INT *result_e);
#endif

#ifndef FUNCTION_fLog2
FIXP_DBL CalcLog2(FIXP_DBL arg, INT arg_e, INT *result_e);
FDK_INLINE FIXP_DBL fLog2(FIXP_DBL x_m, INT x_e, INT *result_e) {
  FIXP_DBL result_m;
  if (x_m <= FL2FXCONST_DBL(0.0f)) {
    *result_e = DFRACT_BITS - 1;
    return FL2FXCONST_DBL(-1.0f);
  }
  {
    FIXP_DBL x2_m;
    {
      INT b_norm;

      b_norm = fNormz(x_m) - 1;
      x2_m = x_m << b_norm;
      x_e = x_e - b_norm;
    }

    x2_m = -(x2_m + FL2FXCONST_DBL(-1.0));
    {
      FIXP_DBL px2_m;
      result_m = FL2FXCONST_DBL(0.0);
      px2_m = x2_m;
      for (int i = 0; i < LD_PRECISION; i++) {
        result_m = fMultAddDiv2(result_m, ldCoeff[i], px2_m);
        px2_m = fMult(px2_m, x2_m);
      }
    }
    result_m =
        fMultAddDiv2(result_m, result_m,
                     FL2FXCONST_DBL(2.0 * 0.4426950408889634073599246810019));

    if (x_e != 0) {
      int enorm;

      enorm = DFRACT_BITS - fNorm((FIXP_DBL)x_e);
      result_m = (result_m >> (enorm - 1)) +
                 ((FIXP_DBL)x_e << (DFRACT_BITS - 1 - enorm));

      *result_e = enorm;
    } else {
      *result_e = 1;
    }
  }

  return result_m;
}

FDK_INLINE FIXP_DBL fLog2(FIXP_DBL x_m, INT x_e) {
  if (x_m <= FL2FXCONST_DBL(0.0f)) {
    x_m = FL2FXCONST_DBL(-1.0f);
  } else {
    INT result_e;
    x_m = fLog2(x_m, x_e, &result_e);
    x_m = scaleValue(x_m, result_e - LD_DATA_SHIFT);
  }
  return x_m;
}

#endif

#ifndef FUNCTION_fAddSaturate
inline FIXP_SGL fAddSaturate(const FIXP_SGL a, const FIXP_SGL b) {
  LONG sum;

  sum = (LONG)(SHORT)a + (LONG)(SHORT)b;
  sum = fMax(fMin((INT)sum, (INT)MAXVAL_SGL), (INT)MINVAL_SGL);
  return (FIXP_SGL)(SHORT)sum;
}

inline FIXP_DBL fAddSaturate(const FIXP_DBL a, const FIXP_DBL b) {
  LONG sum;

  sum = (LONG)(a >> 1) + (LONG)(b >> 1);
  sum = fMax(fMin((INT)sum, (INT)(MAXVAL_DBL >> 1)), (INT)(MINVAL_DBL >> 1));
  return (FIXP_DBL)(LONG)(sum << 1);
}
#endif

INT fixp_floorToInt(FIXP_DBL f_inp, INT sf);
FIXP_DBL fixp_floor(FIXP_DBL f_inp, INT sf);

INT fixp_ceilToInt(FIXP_DBL f_inp, INT sf);
FIXP_DBL fixp_ceil(FIXP_DBL f_inp, INT sf);

INT fixp_truncateToInt(FIXP_DBL f_inp, INT sf);
FIXP_DBL fixp_truncate(FIXP_DBL f_inp, INT sf);

INT fixp_roundToInt(FIXP_DBL f_inp, INT sf);
FIXP_DBL fixp_round(FIXP_DBL f_inp, INT sf);

extern const FIXP_DBL invCount[80];

LNK_SECTION_INITCODE
inline void InitInvInt(void) {}

inline FIXP_DBL GetInvInt(int intValue) {
  return invCount[fMin(fMax(intValue, 0), 80 - 1)];
}

#endif /* FIXPOINT_MATH_H */
