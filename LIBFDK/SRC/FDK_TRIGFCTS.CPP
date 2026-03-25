/******************* Library for basic calculation routines ********************

   Author(s):   Haricharan Lakshman, Manuel Jander

   Description: Trigonometric functions fixed point fractional implementation.

*******************************************************************************/

#include "FDK_trigFcts.h"
#include "fixpoint_math.h"

#define IMPROVE_ATAN2_ACCURACY 1
#define MINSFTAB 7
#define MAXSFTAB 25

#if IMPROVE_ATAN2_ACCURACY
static const FIXP_DBL f_atan_expand_range[MAXSFTAB - (MINSFTAB - 1)] = {
    FL2FXCONST_DBL(7.775862990872099e-001),
    FL2FXCONST_DBL(7.814919928673978e-001),
    FL2FXCONST_DBL(7.834450483314648e-001),
    FL2FXCONST_DBL(7.844216021392089e-001),
    FL2FXCONST_DBL(7.849098823026687e-001),
    FL2FXCONST_DBL(7.851540227918509e-001),
    FL2FXCONST_DBL(7.852760930873737e-001),
    FL2FXCONST_DBL(7.853371282415015e-001),
    FL2FXCONST_DBL(7.853676458193612e-001),
    FL2FXCONST_DBL(7.853829046083906e-001),
    FL2FXCONST_DBL(7.853905340029177e-001),
    FL2FXCONST_DBL(7.853943487001828e-001),
    FL2FXCONST_DBL(7.853962560488155e-001),
    FL2FXCONST_DBL(7.853972097231319e-001),
    FL2FXCONST_DBL(7.853976865602901e-001),
    FL2FXCONST_DBL(7.853979249788692e-001),
    FL2FXCONST_DBL(7.853980441881587e-001),
    FL2FXCONST_DBL(7.853981037928035e-001),
    FL2FXCONST_DBL(7.853981335951259e-001)
};
#endif

FIXP_DBL fixp_atan2(FIXP_DBL y, FIXP_DBL x) {
  FIXP_DBL q;
  FIXP_DBL at;
  FIXP_DBL at2;
  FIXP_DBL ret = FL2FXCONST_DBL(-1.0f);
  INT sf, sfo, stf;

  if (y > FL2FXCONST_DBL(0.0f)) {
    if (x > FL2FXCONST_DBL(0.0f)) {
      q = fDivNormHighPrec(y, x, &sf);
    } else if (x < FL2FXCONST_DBL(0.0f)) {
      q = -fDivNormHighPrec(y, -x, &sf);
    } else {
      q = FL2FXCONST_DBL(+1.0f);
      sf = 0;
    }
  } else if (y < FL2FXCONST_DBL(0.0f)) {
    if (x > FL2FXCONST_DBL(0.0f)) {
      q = -fDivNormHighPrec(-y, x, &sf);
    } else if (x < FL2FXCONST_DBL(0.0f)) {
      q = fDivNormHighPrec(-y, -x, &sf);
    } else {
      q = FL2FXCONST_DBL(-1.0f);
      sf = 0;
    }
  } else {
    q = FL2FXCONST_DBL(0.0f);
    sf = 0;
  }
  sfo = sf;

  if (sfo > ATI_SF) {

#if IMPROVE_ATAN2_ACCURACY
    if (sfo > MAXSFTAB) sfo = MAXSFTAB;
#endif

    if (q > FL2FXCONST_DBL(0.0f)) {
#if IMPROVE_ATAN2_ACCURACY
      at = +f_atan_expand_range[sfo - ATI_SF - 1];
#else
      at = FL2FXCONST_DBL(+M_PI / 2 / ATO_SCALE);
#endif
    } else if (q < FL2FXCONST_DBL(0.0f)) {
#if IMPROVE_ATAN2_ACCURACY
      at = -f_atan_expand_range[sfo - ATI_SF - 1];
#else
      at = FL2FXCONST_DBL(-M_PI / 2 / ATO_SCALE);
#endif
    } else {
      at = FL2FXCONST_DBL(0.0f);
    }
  } else {
    stf = sfo - ATI_SF;
    if (stf > 0)
      q = q << (INT)fMin(stf, DFRACT_BITS - 1);
    else
      q = q >> (INT)fMin(-stf, DFRACT_BITS - 1);
    at = fixp_atan(q);
  }

  at2 = at >> (AT2O_SF - ATO_SF);
  if (x > FL2FXCONST_DBL(0.0f)) {
    ret = at2;
  } else if (x < FL2FXCONST_DBL(0.0f)) {
    if (y >= FL2FXCONST_DBL(0.0f)) {
      ret = at2 + FL2FXCONST_DBL(M_PI / AT2O_SCALE);
    } else {
      ret = at2 - FL2FXCONST_DBL(M_PI / AT2O_SCALE);
    }
  } else {
    if (y > FL2FXCONST_DBL(0.0f)) {
      ret = FL2FXCONST_DBL(+M_PI / 2 / AT2O_SCALE);
    } else if (y < FL2FXCONST_DBL(0.0f)) {
      ret = FL2FXCONST_DBL(-M_PI / 2 / AT2O_SCALE);
    } else if (y == FL2FXCONST_DBL(0.0f)) {
      ret = FL2FXCONST_DBL(0.0f);
    }
  }
  return ret;
}

FIXP_DBL fixp_atan(FIXP_DBL x) {
  INT sign;
  FIXP_DBL result, temp;

  FIXP_DBL P281 = (FIXP_DBL)0x00013000;
  FIXP_DBL ONEP571 = (FIXP_DBL)0x6487ef00;

  if (x < FIXP_DBL(0)) {
    sign = 1;
    x = -x;
  } else {
    sign = 0;
  }
  FDK_ASSERT(FL2FXCONST_DBL(1.0 / 64.0) == Q(Q_ATANINP));
  if (x < FL2FXCONST_DBL(1.0 / 64.0))
  {
    x <<= ATI_SF;
    FIXP_DBL x2 = fPow2(x);
    temp = fMultAddDiv2((FL2FXCONST_DBL(0.1449824901444650f) >> 1), x2,
                        FL2FXCONST_DBL(-0.0382544649702990));
    temp = fMultAddDiv2((FL2FXCONST_DBL(-0.3205332923816640f) >> 2), x2, temp);
    temp = fMultAddDiv2((FL2FXCONST_DBL(0.9991334482227801f) >> 3), x2, temp);
    result = fMult(x, (temp << 2));
  } else if (x < FL2FXCONST_DBL(1.28 / 64.0)) {
    FIXP_DBL delta_fix;
    FIXP_DBL PI_BY_4 = FL2FXCONST_DBL(3.1415926 / 4.0) >> 1;

    delta_fix = (x - FL2FXCONST_DBL(1.0 / 64.0)) << 5;
    result = PI_BY_4 + (delta_fix >> 1) - (fPow2Div2(delta_fix));
  } else {
    INT res_e;

    temp = fPow2Div2(x);
    temp = temp + P281;
    result = fDivNorm(x, temp, &res_e);
    result = scaleValue(result,
                        (Q_ATANOUT - Q_ATANINP + 18 - DFRACT_BITS + 1) + res_e);
    result = ONEP571 - result;
  }
  if (sign) {
    result = -result;
  }

  return (result);
}

#include "FDK_tools_rom.h"

FIXP_DBL fixp_cos(FIXP_DBL x, int scale) {
  FIXP_DBL residual, error, sine, cosine;

  residual = fixp_sin_cos_residual_inline(x, scale, &sine, &cosine);
  error = fMult(sine, residual);

#ifdef SINETABLE_16BIT
  return cosine - error;
#else
  return SATURATE_LEFT_SHIFT(cosine - error, 1, DFRACT_BITS);
#endif
}

FIXP_DBL fixp_sin(FIXP_DBL x, int scale) {
  FIXP_DBL residual, error, sine, cosine;

  residual = fixp_sin_cos_residual_inline(x, scale, &sine, &cosine);
  error = fMult(cosine, residual);

#ifdef SINETABLE_16BIT
  return sine + error;
#else
  return SATURATE_LEFT_SHIFT(sine + error, 1, DFRACT_BITS);
#endif
}

void fixp_cos_sin(FIXP_DBL x, int scale, FIXP_DBL *cos, FIXP_DBL *sin) {
  FIXP_DBL residual, error0, error1, sine, cosine;

  residual = fixp_sin_cos_residual_inline(x, scale, &sine, &cosine);
  error0 = fMult(sine, residual);
  error1 = fMult(cosine, residual);

#ifdef SINETABLE_16BIT
  *cos = cosine - error0;
  *sin = sine + error1;
#else
  *cos = SATURATE_LEFT_SHIFT(cosine - error0, 1, DFRACT_BITS);
  *sin = SATURATE_LEFT_SHIFT(sine + error1, 1, DFRACT_BITS);
#endif
}
