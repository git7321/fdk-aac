/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef TRANSCENDENT_H
#define TRANSCENDENT_H

#include "sbrdecoder.h"
#include "sbr_rom.h"

static inline FIXP_SGL FDK_getNumOctavesDiv8(INT a, /*!< lower band */
                                             INT b) /*!< upper band */
{
  return ((SHORT)((LONG)(CalcLdInt(b) - CalcLdInt(a)) >> (FRACT_BITS - 3)));
}

inline void FDK_add_MantExp(FIXP_SGL a_m, /*!< Mantissa of 1st operand a */
                            SCHAR a_e,    /*!< Exponent of 1st operand a */
                            FIXP_SGL b_m, /*!< Mantissa of 2nd operand b */
                            SCHAR b_e,    /*!< Exponent of 2nd operand b */
                            FIXP_SGL *ptrSum_m, /*!< Mantissa of result */
                            SCHAR *ptrSum_e)    /*!< Exponent of result */
{
  FIXP_DBL accu;
  int shift;
  int shiftAbs;

  FIXP_DBL shiftedMantissa;
  FIXP_DBL otherMantissa;

  shift = (int)(a_e - b_e);

  shiftAbs = (shift > 0) ? shift : -shift;
  shiftAbs = (shiftAbs < DFRACT_BITS - 1) ? shiftAbs : DFRACT_BITS - 1;
  shiftedMantissa = (shift > 0) ? (FX_SGL2FX_DBL(b_m) >> shiftAbs)
                                : (FX_SGL2FX_DBL(a_m) >> shiftAbs);
  otherMantissa = (shift > 0) ? FX_SGL2FX_DBL(a_m) : FX_SGL2FX_DBL(b_m);
  *ptrSum_e = (shift > 0) ? a_e : b_e;

  accu = (shiftedMantissa >> 1) + (otherMantissa >> 1);

  if ((accu >= (FL2FXCONST_DBL(0.5f) - (FIXP_DBL)1)) ||
      (accu <= FL2FXCONST_DBL(-0.5f)))
    *ptrSum_e += 1;
  else
    accu = (shiftedMantissa + otherMantissa);

  *ptrSum_m = FX_DBL2FX_SGL(accu);
}

inline void FDK_add_MantExp(FIXP_DBL a,       /*!< Mantissa of 1st operand a */
                            SCHAR a_e,        /*!< Exponent of 1st operand a */
                            FIXP_DBL b,       /*!< Mantissa of 2nd operand b */
                            SCHAR b_e,        /*!< Exponent of 2nd operand b */
                            FIXP_DBL *ptrSum, /*!< Mantissa of result */
                            SCHAR *ptrSum_e)  /*!< Exponent of result */
{
  FIXP_DBL accu;
  int shift;
  int shiftAbs;

  FIXP_DBL shiftedMantissa;
  FIXP_DBL otherMantissa;

  shift = (int)(a_e - b_e);

  shiftAbs = (shift > 0) ? shift : -shift;
  shiftAbs = (shiftAbs < DFRACT_BITS - 1) ? shiftAbs : DFRACT_BITS - 1;
  shiftedMantissa = (shift > 0) ? (b >> shiftAbs) : (a >> shiftAbs);
  otherMantissa = (shift > 0) ? a : b;
  *ptrSum_e = (shift > 0) ? a_e : b_e;

  accu = (shiftedMantissa >> 1) + (otherMantissa >> 1);

  if ((accu >= (FL2FXCONST_DBL(0.5f) - (FIXP_DBL)1)) ||
      (accu <= FL2FXCONST_DBL(-0.5f)))
    *ptrSum_e += 1;
  else
    accu = (shiftedMantissa + otherMantissa);

  *ptrSum = accu;
}

static inline void FDK_divide_MantExp(
    FIXP_SGL a_m,          /*!< Mantissa of dividend a */
    SCHAR a_e,             /*!< Exponent of dividend a */
    FIXP_SGL b_m,          /*!< Mantissa of divisor b */
    SCHAR b_e,             /*!< Exponent of divisor b */
    FIXP_SGL *ptrResult_m, /*!< Mantissa of quotient a/b */
    SCHAR *ptrResult_e)    /*!< Exponent of quotient a/b */

{
  int preShift, postShift, index, shift;
  FIXP_DBL ratio_m;
  FIXP_SGL bInv_m = FL2FXCONST_SGL(0.0f);

  preShift = CntLeadingZeros(FX_SGL2FX_DBL(b_m));

  shift = (FRACT_BITS - 2 - INV_TABLE_BITS - preShift);

  index = (shift < 0) ? (LONG)b_m << (-shift) : (LONG)b_m >> shift;

  index &= (1 << (INV_TABLE_BITS + 1)) - 1;

  index--;

  index = index >> 1;

  bInv_m = (index < 0) ? bInv_m : FDK_sbrDecoder_invTable[index];

  ratio_m = (index < 0) ? FX_SGL2FX_DBL(a_m >> 1) : fMultDiv2(bInv_m, a_m);

  postShift = CntLeadingZeros(ratio_m) - 1;

  *ptrResult_m = FX_DBL2FX_SGL(ratio_m << postShift);
  *ptrResult_e = a_e - b_e + 1 + preShift - postShift;
}

static inline void FDK_divide_MantExp(
    FIXP_DBL a_m,          /*!< Mantissa of dividend a */
    SCHAR a_e,             /*!< Exponent of dividend a */
    FIXP_DBL b_m,          /*!< Mantissa of divisor b */
    SCHAR b_e,             /*!< Exponent of divisor b */
    FIXP_DBL *ptrResult_m, /*!< Mantissa of quotient a/b */
    SCHAR *ptrResult_e)    /*!< Exponent of quotient a/b */

{
  int preShift, postShift, index, shift;
  FIXP_DBL ratio_m;
  FIXP_SGL bInv_m = FL2FXCONST_SGL(0.0f);

  preShift = CntLeadingZeros(b_m);

  shift = (DFRACT_BITS - 2 - INV_TABLE_BITS - preShift);

  index = (shift < 0) ? (LONG)b_m << (-shift) : (LONG)b_m >> shift;

  index &= (1 << (INV_TABLE_BITS + 1)) - 1;

  index--;

  index = index >> 1;

  bInv_m = (index < 0) ? bInv_m : FDK_sbrDecoder_invTable[index];

  ratio_m = (index < 0) ? (a_m >> 1) : fMultDiv2(bInv_m, a_m);

  postShift = CntLeadingZeros(ratio_m) - 1;

  *ptrResult_m = ratio_m << postShift;
  *ptrResult_e = a_e - b_e + 1 + preShift - postShift;
}

static inline void FDK_sqrt_MantExp(
    FIXP_DBL *mantissa,
    SCHAR *exponent, const SCHAR *destScale) {
  FIXP_DBL input_m = *mantissa;
  int input_e = (int)*exponent;
  FIXP_DBL result = FL2FXCONST_DBL(0.0f);
  int result_e = -FRACT_BITS;

  result = sqrtFixp_lookup(input_m, &input_e);
  result_e = input_e;

  if (exponent == destScale) {
    *mantissa = result;
    *exponent = result_e;
  } else {
    int shift = result_e - *destScale;
    *mantissa = (shift >= 0) ? result << (INT)fixMin(DFRACT_BITS - 1, shift)
                             : result >> (INT)fixMin(DFRACT_BITS - 1, -shift);
    *exponent = *destScale;
  }
}

#endif
