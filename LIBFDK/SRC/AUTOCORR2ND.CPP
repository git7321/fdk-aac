/******************* Library for basic calculation routines ********************

   Author(s):   M. Lohwasser

   Description: auto-correlation functions

*******************************************************************************/

#include "autocorr2nd.h"

#if !defined(FUNCTION_autoCorr2nd_real)
INT autoCorr2nd_real(
    ACORR_COEFS *ac,
    const FIXP_DBL *reBuffer,
    const int len
) {
  int j, autoCorrScaling, mScale;

  FIXP_DBL accu1, accu2, accu3, accu4, accu5;

  const FIXP_DBL *pReBuf;

  const FIXP_DBL *realBuf = reBuffer;

  const int len_scale = fMax(DFRACT_BITS - fNormz((FIXP_DBL)(len / 2)), 1);
  pReBuf = realBuf - 2;
  accu5 =
      ((fMultDiv2(pReBuf[0], pReBuf[2]) + fMultDiv2(pReBuf[1], pReBuf[3])) >>
       len_scale);
  pReBuf++;

  accu1 = fPow2Div2(pReBuf[0]) >> len_scale;
  accu3 = fMultDiv2(pReBuf[0], pReBuf[1]) >> len_scale;
  pReBuf++;

  for (j = (len - 2) >> 1; j != 0; j--, pReBuf += 2) {
    accu1 += ((fPow2Div2(pReBuf[0]) + fPow2Div2(pReBuf[1])) >> len_scale);

    accu3 +=
        ((fMultDiv2(pReBuf[0], pReBuf[1]) + fMultDiv2(pReBuf[1], pReBuf[2])) >>
         len_scale);

    accu5 +=
        ((fMultDiv2(pReBuf[0], pReBuf[2]) + fMultDiv2(pReBuf[1], pReBuf[3])) >>
         len_scale);
  }

  accu2 = (fPow2Div2(realBuf[-2]) >> len_scale);
  accu2 += accu1;

  accu1 += (fPow2Div2(realBuf[len - 2]) >> len_scale);

  accu4 = (fMultDiv2(realBuf[-1], realBuf[-2]) >> len_scale);
  accu4 += accu3;

  accu3 += (fMultDiv2(realBuf[len - 1], realBuf[len - 2]) >> len_scale);

  mScale = CntLeadingZeros(
               (accu1 | accu2 | fAbs(accu3) | fAbs(accu4) | fAbs(accu5))) -
           1;
  autoCorrScaling = mScale - 1 - len_scale;

  ac->r11r = accu1 << mScale;
  ac->r22r = accu2 << mScale;
  ac->r01r = accu3 << mScale;
  ac->r12r = accu4 << mScale;
  ac->r02r = accu5 << mScale;

  ac->det = (fMultDiv2(ac->r11r, ac->r22r) - fMultDiv2(ac->r12r, ac->r12r));
  mScale = CountLeadingBits(fAbs(ac->det));

  ac->det <<= mScale;
  ac->det_scale = mScale - 1;

  return autoCorrScaling;
}
#endif

#if !defined(FUNCTION_autoCorr2nd_cplx)
INT autoCorr2nd_cplx(
    ACORR_COEFS *ac,
    const FIXP_DBL *reBuffer,
    const FIXP_DBL *imBuffer,
    const int len
) {
  int j, autoCorrScaling, mScale;

  FIXP_DBL accu0, accu1, accu2, accu3, accu4, accu5, accu6, accu7, accu8;

  const FIXP_DBL *pReBuf, *pImBuf;

  const FIXP_DBL *realBuf = reBuffer;
  const FIXP_DBL *imagBuf = imBuffer;

  const int len_scale = fMax(DFRACT_BITS - fNormz((FIXP_DBL)len), 1);
  accu1 = accu3 = accu5 = accu7 = accu8 = FL2FXCONST_DBL(0.0f);

  pReBuf = realBuf - 2, pImBuf = imagBuf - 2;
  accu7 +=
      ((fMultDiv2(pReBuf[2], pReBuf[0]) + fMultDiv2(pImBuf[2], pImBuf[0])) >>
       len_scale);
  accu8 +=
      ((fMultDiv2(pImBuf[2], pReBuf[0]) - fMultDiv2(pReBuf[2], pImBuf[0])) >>
       len_scale);

  pReBuf = realBuf - 1, pImBuf = imagBuf - 1;
  for (j = (len - 1); j != 0; j--, pReBuf++, pImBuf++) {
    accu1 += ((fPow2Div2(pReBuf[0]) + fPow2Div2(pImBuf[0])) >> len_scale);
    accu3 +=
        ((fMultDiv2(pReBuf[0], pReBuf[1]) + fMultDiv2(pImBuf[0], pImBuf[1])) >>
         len_scale);
    accu5 +=
        ((fMultDiv2(pImBuf[1], pReBuf[0]) - fMultDiv2(pReBuf[1], pImBuf[0])) >>
         len_scale);
    accu7 +=
        ((fMultDiv2(pReBuf[2], pReBuf[0]) + fMultDiv2(pImBuf[2], pImBuf[0])) >>
         len_scale);
    accu8 +=
        ((fMultDiv2(pImBuf[2], pReBuf[0]) - fMultDiv2(pReBuf[2], pImBuf[0])) >>
         len_scale);
  }

  accu2 = ((fPow2Div2(realBuf[-2]) + fPow2Div2(imagBuf[-2])) >> len_scale);
  accu2 += accu1;

  accu1 += ((fPow2Div2(realBuf[len - 2]) + fPow2Div2(imagBuf[len - 2])) >>
            len_scale);
  accu0 = ((fPow2Div2(realBuf[len - 1]) + fPow2Div2(imagBuf[len - 1])) >>
           len_scale) -
          ((fPow2Div2(realBuf[-1]) + fPow2Div2(imagBuf[-1])) >> len_scale);
  accu0 += accu1;

  accu4 = ((fMultDiv2(realBuf[-1], realBuf[-2]) +
            fMultDiv2(imagBuf[-1], imagBuf[-2])) >>
           len_scale);
  accu4 += accu3;

  accu3 += ((fMultDiv2(realBuf[len - 1], realBuf[len - 2]) +
             fMultDiv2(imagBuf[len - 1], imagBuf[len - 2])) >>
            len_scale);

  accu6 = ((fMultDiv2(imagBuf[-1], realBuf[-2]) -
            fMultDiv2(realBuf[-1], imagBuf[-2])) >>
           len_scale);
  accu6 += accu5;

  accu5 += ((fMultDiv2(imagBuf[len - 1], realBuf[len - 2]) -
             fMultDiv2(realBuf[len - 1], imagBuf[len - 2])) >>
            len_scale);

  mScale =
      CntLeadingZeros((accu0 | accu1 | accu2 | fAbs(accu3) | fAbs(accu4) |
                       fAbs(accu5) | fAbs(accu6) | fAbs(accu7) | fAbs(accu8))) -
      1;
  autoCorrScaling = mScale - 1 - len_scale;

  ac->r00r = (FIXP_DBL)accu0 << mScale;
  ac->r11r = (FIXP_DBL)accu1 << mScale;
  ac->r22r = (FIXP_DBL)accu2 << mScale;
  ac->r01r = (FIXP_DBL)accu3 << mScale;
  ac->r12r = (FIXP_DBL)accu4 << mScale;
  ac->r01i = (FIXP_DBL)accu5 << mScale;
  ac->r12i = (FIXP_DBL)accu6 << mScale;
  ac->r02r = (FIXP_DBL)accu7 << mScale;
  ac->r02i = (FIXP_DBL)accu8 << mScale;

  ac->det =
      (fMultDiv2(ac->r11r, ac->r22r) >> 1) -
      ((fMultDiv2(ac->r12r, ac->r12r) + fMultDiv2(ac->r12i, ac->r12i)) >> 1);
  mScale = CntLeadingZeros(fAbs(ac->det)) - 1;

  ac->det <<= mScale;
  ac->det_scale = mScale - 2;

  return autoCorrScaling;
}

#endif /* FUNCTION_autoCorr2nd_cplx */
