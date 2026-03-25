/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "sbr_misc.h"

void FDKsbrEnc_Shellsort_fract(FIXP_DBL *in, INT n) {
  FIXP_DBL v;
  INT i, j;
  INT inc = 1;

  do
    inc = 3 * inc + 1;
  while (inc <= n);

  do {
    inc = inc / 3;
    for (i = inc + 1; i <= n; i++) {
      v = in[i - 1];
      j = i;
      while (in[j - inc - 1] > v) {
        in[j - 1] = in[j - inc - 1];
        j -= inc;
        if (j <= inc) break;
      }
      in[j - 1] = v;
    }
  } while (inc > 1);
}

void FDKsbrEnc_Shellsort_int(INT *in, INT n) {
  INT i, j, v;
  INT inc = 1;

  do
    inc = 3 * inc + 1;
  while (inc <= n);

  do {
    inc = inc / 3;
    for (i = inc + 1; i <= n; i++) {
      v = in[i - 1];
      j = i;
      while (in[j - inc - 1] > v) {
        in[j - 1] = in[j - inc - 1];
        j -= inc;
        if (j <= inc) break;
      }
      in[j - 1] = v;
    }
  } while (inc > 1);
}

void FDKsbrEnc_AddVecLeft(INT *dst, INT *length_dst, INT *src, INT length_src) {
  INT i;

  for (i = length_src - 1; i >= 0; i--)
    FDKsbrEnc_AddLeft(dst, length_dst, src[i]);
}

void FDKsbrEnc_AddLeft(INT *vector, INT *length_vector, INT value) {
  INT i;

  for (i = *length_vector; i > 0; i--) vector[i] = vector[i - 1];
  vector[0] = value;
  (*length_vector)++;
}

void FDKsbrEnc_AddRight(INT *vector, INT *length_vector, INT value) {
  vector[*length_vector] = value;
  (*length_vector)++;
}

void FDKsbrEnc_AddVecRight(INT *dst, INT *length_dst, INT *src,
                           INT length_src) {
  INT i;
  for (i = 0; i < length_src; i++) FDKsbrEnc_AddRight(dst, length_dst, src[i]);
}

FIXP_DBL FDKsbrEnc_LSI_divide_scale_fract(FIXP_DBL num, FIXP_DBL denom,
                                          FIXP_DBL scale) {
  FIXP_DBL tmp = FL2FXCONST_DBL(0.0f);
  if (num != FL2FXCONST_DBL(0.0f)) {
    INT shiftCommon;
    INT shiftNum = CountLeadingBits(num);
    INT shiftDenom = CountLeadingBits(denom);
    INT shiftScale = CountLeadingBits(scale);

    num = num << shiftNum;
    scale = scale << shiftScale;

    tmp = fMultDiv2(num, scale);

    if (denom > (tmp >> fixMin(shiftNum + shiftScale - 1, (DFRACT_BITS - 1)))) {
      denom = denom << shiftDenom;
      tmp = schur_div(tmp, denom, 15);
      shiftCommon =
          fixMin((shiftNum - shiftDenom + shiftScale - 1), (DFRACT_BITS - 1));
      if (shiftCommon < 0)
        tmp <<= -shiftCommon;
      else
        tmp >>= shiftCommon;
    } else {
      tmp = (FIXP_DBL)MAXVAL_DBL;
    }
  }

  return (tmp);
}
