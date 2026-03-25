/******************* Library for basic calculation routines ********************

   Author(s):   M. Lohwasser

   Description: fixed point abs definitions

*******************************************************************************/

#ifndef AUTOCORR2ND_H
#define AUTOCORR2ND_H

#include "common_fix.h"

typedef struct {
  FIXP_DBL r00r;
  FIXP_DBL r11r;
  FIXP_DBL r22r;
  FIXP_DBL r01r;
  FIXP_DBL r02r;
  FIXP_DBL r12r;
  FIXP_DBL r01i;
  FIXP_DBL r02i;
  FIXP_DBL r12i;
  FIXP_DBL det;
  int det_scale;
} ACORR_COEFS;

#define LPC_ORDER 2

INT autoCorr2nd_real(
    ACORR_COEFS *ac,
    const FIXP_DBL *reBuffer,
    const int len
);

INT autoCorr2nd_cplx(
    ACORR_COEFS *ac,
    const FIXP_DBL *reBuffer,
    const FIXP_DBL *imBuffer,
    const int len
);

#endif /* AUTOCORR2ND_H */
