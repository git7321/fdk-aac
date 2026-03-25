/*********************** MPEG surround encoder library *************************

   Author(s):   M. Luis Valero

   Description: Enhanced Time Domain Downmix

*******************************************************************************/

#include "sacenc_dmx_tdom_enh.h"
#include "FDK_matrixCalloc.h"
#include "FDK_trigFcts.h"
#include "fixpoint_math.h"

#define PI_FLT 3.1415926535897931f
#define ALPHA_FLT 0.0001f

#define PI_E (2)
#define PI_M (FL2FXCONST_DBL(PI_FLT / (1 << PI_E)))

#define ALPHA_E (13)
#define ALPHA_M (FL2FXCONST_DBL(ALPHA_FLT * (1 << ALPHA_E)))

enum { L = 0, R = 1 };

typedef struct T_ENHANCED_TIME_DOMAIN_DMX {
  int maxFramelength;

  int framelength;

  FIXP_DBL prev_gain_m[2];
  INT prev_gain_e;
  FIXP_DBL prev_H1_m[2];
  INT prev_H1_e;

  FIXP_DBL *sinusWindow_m;
  SCHAR sinusWindow_e;

  FIXP_DBL prev_Left_m;
  INT prev_Left_e;
  FIXP_DBL prev_Right_m;
  INT prev_Right_e;
  FIXP_DBL prev_XNrg_m;
  INT prev_XNrg_e;

  FIXP_DBL lin_bbCld_weight_m;
  INT lin_bbCld_weight_e;
  FIXP_DBL gain_weight_m[2];
  INT gain_weight_e;

} ENHANCED_TIME_DOMAIN_DMX;

static void calculateRatio(const FIXP_DBL sqrt_linCld_m,
                           const INT sqrt_linCld_e, const FIXP_DBL lin_Cld_m,
                           const INT lin_Cld_e, const FIXP_DBL Icc_m,
                           const INT Icc_e, FIXP_DBL G_m[2], INT *G_e);

static void calculateDmxGains(const FIXP_DBL lin_Cld_m, const INT lin_Cld_e,
                              const FIXP_DBL lin_Cld2_m, const INT lin_Cld2_e,
                              const FIXP_DBL Icc_m, const INT Icc_e,
                              const FIXP_DBL G_m[2], const INT G_e,
                              FIXP_DBL H1_m[2], INT *pH1_e);

static FIXP_DBL invSqrtNorm2(const FIXP_DBL op_m, const INT op_e,
                             INT *const result_e) {
  FIXP_DBL src_m = op_m;
  int src_e = op_e;

  if (src_e & 1) {
    src_m >>= 1;
    src_e += 1;
  }

  src_m = invSqrtNorm2(src_m, result_e);
  *result_e = (*result_e) - (src_e >> 1);

  return src_m;
}

static FIXP_DBL sqrtFixp(const FIXP_DBL op_m, const INT op_e,
                         INT *const result_e) {
  FIXP_DBL src_m = op_m;
  int src_e = op_e;

  if (src_e & 1) {
    src_m >>= 1;
    src_e += 1;
  }

  *result_e = (src_e >> 1);
  return sqrtFixp(src_m);
}

static FIXP_DBL fixpAdd(const FIXP_DBL src1_m, const INT src1_e,
                        const FIXP_DBL src2_m, const INT src2_e,
                        INT *const dst_e) {
  FIXP_DBL dst_m;

  if (src1_m == FL2FXCONST_DBL(0.f)) {
    *dst_e = src2_e;
    dst_m = src2_m;
  } else if (src2_m == FL2FXCONST_DBL(0.f)) {
    *dst_e = src1_e;
    dst_m = src1_m;
  } else {
    *dst_e = fixMax(src1_e, src2_e) + 1;
    dst_m =
        scaleValue(src1_m, fixMax((src1_e - (*dst_e)), -(DFRACT_BITS - 1))) +
        scaleValue(src2_m, fixMax((src2_e - (*dst_e)), -(DFRACT_BITS - 1)));
  }
  return dst_m;
}

static void fixpAddNorm(const FIXP_DBL value1, const INT q1,
                        FIXP_DBL *const pValue2, INT *const pQ2) {
  const int headroom1 = fNormz(fixp_abs(value1)) - 1;
  const int headroom2 = fNormz(fixp_abs(*pValue2)) - 1;
  int resultScale = fixMax(q1 - headroom1, (*pQ2) - headroom2);

  if ((value1 != FL2FXCONST_DBL(0.f)) && (*pValue2 != FL2FXCONST_DBL(0.f))) {
    resultScale++;
  }

  *pValue2 =
      scaleValue(value1, q1 - resultScale) +
      scaleValue(*pValue2, fixMax(-(DFRACT_BITS - 1), ((*pQ2) - resultScale)));
  *pQ2 = (*pValue2 != (FIXP_DBL)0) ? resultScale : DFRACT_BITS - 1;
}

FDK_SACENC_ERROR fdk_sacenc_open_enhancedTimeDomainDmx(
    HANDLE_ENHANCED_TIME_DOMAIN_DMX *phEnhancedTimeDmx, const INT framelength) {
  FDK_SACENC_ERROR error = SACENC_OK;
  HANDLE_ENHANCED_TIME_DOMAIN_DMX hEnhancedTimeDmx = NULL;

  if (NULL == phEnhancedTimeDmx) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDK_ALLOCATE_MEMORY_1D(hEnhancedTimeDmx, 1, ENHANCED_TIME_DOMAIN_DMX);
    FDK_ALLOCATE_MEMORY_1D(hEnhancedTimeDmx->sinusWindow_m, 1 + framelength,
                           FIXP_DBL);
    hEnhancedTimeDmx->maxFramelength = framelength;
    *phEnhancedTimeDmx = hEnhancedTimeDmx;
  }
  return error;

bail:
  fdk_sacenc_close_enhancedTimeDomainDmx(&hEnhancedTimeDmx);
  return ((SACENC_OK == error) ? SACENC_MEMORY_ERROR : error);
}

FDK_SACENC_ERROR fdk_sacenc_init_enhancedTimeDomainDmx(
    HANDLE_ENHANCED_TIME_DOMAIN_DMX hEnhancedTimeDmx,
    const FIXP_DBL *const pInputGain_m, const INT inputGain_e,
    const FIXP_DBL outputGain_m, const INT outputGain_e,
    const INT framelength) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (hEnhancedTimeDmx == NULL) {
    error = SACENC_INVALID_HANDLE;
  } else {
    int smp;
    if (framelength > hEnhancedTimeDmx->maxFramelength) {
      error = SACENC_INIT_ERROR;
      goto bail;
    }

    hEnhancedTimeDmx->framelength = framelength;

    INT deltax_e;
    FIXP_DBL deltax_m;

    deltax_m = fDivNormHighPrec(
        PI_M, (FIXP_DBL)(2 * hEnhancedTimeDmx->framelength), &deltax_e);
    deltax_m = scaleValue(deltax_m, PI_E + deltax_e - (DFRACT_BITS - 1) - 1);
    deltax_e = 1;

    for (smp = 0; smp < hEnhancedTimeDmx->framelength + 1; smp++) {
      hEnhancedTimeDmx->sinusWindow_m[smp] =
          fMult(ALPHA_M, fPow2(fixp_sin(smp * deltax_m, deltax_e)));
    }
    hEnhancedTimeDmx->sinusWindow_e = -ALPHA_E;

    hEnhancedTimeDmx->prev_Left_m = hEnhancedTimeDmx->prev_Right_m =
        hEnhancedTimeDmx->prev_XNrg_m = FL2FXCONST_DBL(0.f);
    hEnhancedTimeDmx->prev_Left_e = hEnhancedTimeDmx->prev_Right_e =
        hEnhancedTimeDmx->prev_XNrg_e = DFRACT_BITS - 1;

    hEnhancedTimeDmx->lin_bbCld_weight_m =
        fDivNormHighPrec(fPow2(pInputGain_m[L]), fPow2(pInputGain_m[R]),
                         &hEnhancedTimeDmx->lin_bbCld_weight_e);

    hEnhancedTimeDmx->gain_weight_m[L] = fMult(pInputGain_m[L], outputGain_m);
    hEnhancedTimeDmx->gain_weight_m[R] = fMult(pInputGain_m[R], outputGain_m);
    hEnhancedTimeDmx->gain_weight_e =
        -fNorm(fixMax(hEnhancedTimeDmx->gain_weight_m[L],
                      hEnhancedTimeDmx->gain_weight_m[R]));

    hEnhancedTimeDmx->gain_weight_m[L] = scaleValue(
        hEnhancedTimeDmx->gain_weight_m[L], -hEnhancedTimeDmx->gain_weight_e);
    hEnhancedTimeDmx->gain_weight_m[R] = scaleValue(
        hEnhancedTimeDmx->gain_weight_m[R], -hEnhancedTimeDmx->gain_weight_e);
    hEnhancedTimeDmx->gain_weight_e += inputGain_e + outputGain_e;

    hEnhancedTimeDmx->prev_gain_m[L] = hEnhancedTimeDmx->gain_weight_m[L] >> 1;
    hEnhancedTimeDmx->prev_gain_m[R] = hEnhancedTimeDmx->gain_weight_m[R] >> 1;
    hEnhancedTimeDmx->prev_gain_e = hEnhancedTimeDmx->gain_weight_e + 1;

    hEnhancedTimeDmx->prev_H1_m[L] =
        scaleValue(hEnhancedTimeDmx->gain_weight_m[L], -4);
    hEnhancedTimeDmx->prev_H1_m[R] =
        scaleValue(hEnhancedTimeDmx->gain_weight_m[R], -4);
    hEnhancedTimeDmx->prev_H1_e = 2 + 2 + hEnhancedTimeDmx->gain_weight_e;
  }
bail:
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_apply_enhancedTimeDomainDmx(
    HANDLE_ENHANCED_TIME_DOMAIN_DMX hEnhancedTimeDmx,
    const INT_PCM *const *const inputTime, INT_PCM *const outputTimeDmx,
    const INT InputDelay) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((NULL == hEnhancedTimeDmx) || (NULL == inputTime) ||
      (NULL == inputTime[L]) || (NULL == inputTime[R]) ||
      (NULL == outputTimeDmx)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    int smp;
    FIXP_DBL lin_bbCld_m, lin_Cld_m, bbCorr_m, sqrt_linCld_m, G_m[2], H1_m[2],
        gainLeft_m, gainRight_m;
    FIXP_DBL bbNrgLeft_m, bbNrgRight_m, bbXNrg_m, nrgLeft_m, nrgRight_m, nrgX_m;
    INT lin_bbCld_e, lin_Cld_e, bbCorr_e, sqrt_linCld_e, G_e, H1_e;
    INT bbNrgLeft_e, bbNrgRight_e, bbXNrg_e, nrgLeft_e, nrgRight_e, nrgX_e;

    const int granuleLength = fixMin(128, hEnhancedTimeDmx->framelength);
    int granuleShift =
        (granuleLength > 1)
            ? ((DFRACT_BITS - 1) - fNorm((FIXP_DBL)(granuleLength - 1)))
            : 0;
    granuleShift = fixMax(
        3, granuleShift +
               1);

    smp = 0;

    bbNrgLeft_m = bbNrgRight_m = bbXNrg_m = (FIXP_DBL)(1);
    bbNrgLeft_e = bbNrgRight_e = bbXNrg_e = 0;

    do {
      const int offset = smp;
      FIXP_DBL partialL, partialR, partialX;
      partialL = partialR = partialX = FL2FXCONST_DBL(0.f);

      int in_margin = FDKmin(
          getScalefactorPCM(
              &inputTime[L][offset],
              fixMin(offset + granuleLength, hEnhancedTimeDmx->framelength) -
                  offset,
              1),
          getScalefactorPCM(
              &inputTime[R][offset],
              fixMin(offset + granuleLength, hEnhancedTimeDmx->framelength) -
                  offset,
              1));

      for (smp = offset;
           smp < fixMin(offset + granuleLength, hEnhancedTimeDmx->framelength);
           smp++) {
        FIXP_PCM inputL =
            scaleValue((FIXP_PCM)inputTime[L][smp], in_margin - 1);
        FIXP_PCM inputR =
            scaleValue((FIXP_PCM)inputTime[R][smp], in_margin - 1);

        partialL += fPow2Div2(inputL) >> (granuleShift - 3);
        partialR += fPow2Div2(inputR) >> (granuleShift - 3);
        partialX += fMultDiv2(inputL, inputR) >> (granuleShift - 3);
      }

      fixpAddNorm(partialL, granuleShift - 2 * in_margin, &bbNrgLeft_m,
                  &bbNrgLeft_e);
      fixpAddNorm(partialR, granuleShift - 2 * in_margin, &bbNrgRight_m,
                  &bbNrgRight_e);
      fixpAddNorm(partialX, granuleShift - 2 * in_margin, &bbXNrg_m, &bbXNrg_e);
    } while (smp < hEnhancedTimeDmx->framelength);

    nrgLeft_m =
        fixpAdd(hEnhancedTimeDmx->prev_Left_m, hEnhancedTimeDmx->prev_Left_e,
                bbNrgLeft_m, bbNrgLeft_e, &nrgLeft_e);
    nrgRight_m =
        fixpAdd(hEnhancedTimeDmx->prev_Right_m, hEnhancedTimeDmx->prev_Right_e,
                bbNrgRight_m, bbNrgRight_e, &nrgRight_e);
    nrgX_m =
        fixpAdd(hEnhancedTimeDmx->prev_XNrg_m, hEnhancedTimeDmx->prev_XNrg_e,
                bbXNrg_m, bbXNrg_e, &nrgX_e);

    lin_bbCld_m = fMult(hEnhancedTimeDmx->lin_bbCld_weight_m,
                        fDivNorm(nrgLeft_m, nrgRight_m, &lin_bbCld_e));
    lin_bbCld_e +=
        hEnhancedTimeDmx->lin_bbCld_weight_e + nrgLeft_e - nrgRight_e;

    bbCorr_m = fMult(nrgX_m, invSqrtNorm2(fMult(nrgLeft_m, nrgRight_m),
                                          nrgLeft_e + nrgRight_e, &bbCorr_e));
    bbCorr_e += nrgX_e;

    hEnhancedTimeDmx->prev_Left_m = bbNrgLeft_m;
    hEnhancedTimeDmx->prev_Left_e = bbNrgLeft_e;
    hEnhancedTimeDmx->prev_Right_m = bbNrgRight_m;
    hEnhancedTimeDmx->prev_Right_e = bbNrgRight_e;
    hEnhancedTimeDmx->prev_XNrg_m = bbXNrg_m;
    hEnhancedTimeDmx->prev_XNrg_e = bbXNrg_e;

    lin_Cld_m = sqrtFixp(lin_bbCld_m, lin_bbCld_e, &lin_Cld_e);
    sqrt_linCld_m = sqrtFixp(lin_Cld_m, lin_Cld_e, &sqrt_linCld_e);

    calculateRatio(sqrt_linCld_m, sqrt_linCld_e, lin_Cld_m, lin_Cld_e, bbCorr_m,
                   bbCorr_e, G_m, &G_e);

    calculateDmxGains(lin_Cld_m, lin_Cld_e, lin_bbCld_m, lin_bbCld_e, bbCorr_m,
                      bbCorr_e, G_m, G_e, H1_m, &H1_e);

    H1_m[L] = fMult(H1_m[L], hEnhancedTimeDmx->gain_weight_m[L]);
    H1_m[R] = fMult(H1_m[R], hEnhancedTimeDmx->gain_weight_m[R]);
    H1_e += hEnhancedTimeDmx->gain_weight_e;

    gainLeft_m = hEnhancedTimeDmx->prev_gain_m[L];
    gainRight_m = hEnhancedTimeDmx->prev_gain_m[R];

    INT intermediate_gain_e =
        +hEnhancedTimeDmx->sinusWindow_e + H1_e - hEnhancedTimeDmx->prev_gain_e;

    for (smp = 0; smp < hEnhancedTimeDmx->framelength; smp++) {
      const INT N = hEnhancedTimeDmx->framelength;
      FIXP_DBL intermediate_gainLeft_m, intermediate_gainRight_m, tmp;

      intermediate_gainLeft_m =
          scaleValue((fMult(hEnhancedTimeDmx->sinusWindow_m[smp], H1_m[L]) +
                      fMult(hEnhancedTimeDmx->sinusWindow_m[N - smp],
                            hEnhancedTimeDmx->prev_H1_m[L])),
                     intermediate_gain_e);
      intermediate_gainRight_m =
          scaleValue((fMult(hEnhancedTimeDmx->sinusWindow_m[smp], H1_m[R]) +
                      fMult(hEnhancedTimeDmx->sinusWindow_m[N - smp],
                            hEnhancedTimeDmx->prev_H1_m[R])),
                     intermediate_gain_e);

      gainLeft_m = intermediate_gainLeft_m +
                   fMult(FL2FXCONST_DBL(1.f - ALPHA_FLT), gainLeft_m);
      gainRight_m = intermediate_gainRight_m +
                    fMult(FL2FXCONST_DBL(1.f - ALPHA_FLT), gainRight_m);

      tmp = fMultDiv2(gainLeft_m, (FIXP_PCM)inputTime[L][smp + InputDelay]) +
            fMultDiv2(gainRight_m, (FIXP_PCM)inputTime[R][smp + InputDelay]);
      outputTimeDmx[smp] = (INT_PCM)SATURATE_SHIFT(
          tmp,
          -(hEnhancedTimeDmx->prev_gain_e + 1 - (DFRACT_BITS - SAMPLE_BITS)),
          SAMPLE_BITS);
    }

    hEnhancedTimeDmx->prev_gain_m[L] = gainLeft_m;
    hEnhancedTimeDmx->prev_gain_m[R] = gainRight_m;

    hEnhancedTimeDmx->prev_H1_m[L] = H1_m[L];
    hEnhancedTimeDmx->prev_H1_m[R] = H1_m[R];
    hEnhancedTimeDmx->prev_H1_e = H1_e;
  }

  return error;
}

static void calculateRatio(const FIXP_DBL sqrt_linCld_m,
                           const INT sqrt_linCld_e, const FIXP_DBL lin_Cld_m,
                           const INT lin_Cld_e, const FIXP_DBL Icc_m,
                           const INT Icc_e, FIXP_DBL G_m[2], INT *G_e) {
#define G_SCALE_FACTOR (2)

  if (Icc_m >= FL2FXCONST_DBL(0.f)) {
    G_m[0] = G_m[1] = FL2FXCONST_DBL(1.f / (float)(1 << G_SCALE_FACTOR));
    G_e[0] = G_SCALE_FACTOR;
  } else {
    const FIXP_DBL max_gain_factor =
        FL2FXCONST_DBL(2.f / (float)(1 << G_SCALE_FACTOR));
    FIXP_DBL tmp1_m, tmp2_m, numerator_m, denominator_m, r_m, r4_m, q;
    INT tmp1_e, tmp2_e, numerator_e, denominator_e, r_e, r4_e;

    tmp1_m =
        fixpAdd(lin_Cld_m, lin_Cld_e, FL2FXCONST_DBL(1.f / 2.f), 1, &tmp1_e);

    tmp2_m = fMult(Icc_m, sqrt_linCld_m);
    tmp2_e = 1 + Icc_e + sqrt_linCld_e;
    numerator_m = fixpAdd(tmp1_m, tmp1_e, tmp2_m, tmp2_e, &numerator_e);
    denominator_m = fixpAdd(tmp1_m, tmp1_e, -tmp2_m, tmp2_e, &denominator_e);

    if ((numerator_m > FL2FXCONST_DBL(0.f)) &&
        (denominator_m > FL2FXCONST_DBL(0.f))) {
      r_m = fDivNorm(numerator_m, denominator_m, &r_e);
      r_e += numerator_e - denominator_e;

      r4_m = sqrtFixp(r_m, r_e, &r4_e);
      r4_m = sqrtFixp(r4_m, r4_e, &r4_e);

      r4_e -= G_SCALE_FACTOR;

      q = ((r4_e >= 0) && (r4_m >= (max_gain_factor >> r4_e)))
              ? max_gain_factor
              : scaleValue(r4_m, r4_e);
    } else {
      q = FL2FXCONST_DBL(0.f);
    }

    G_m[0] = max_gain_factor - q;
    G_m[1] = q;

    *G_e = G_SCALE_FACTOR;
  }
}

static void calculateDmxGains(const FIXP_DBL lin_Cld_m, const INT lin_Cld_e,
                              const FIXP_DBL lin_Cld2_m, const INT lin_Cld2_e,
                              const FIXP_DBL Icc_m, const INT Icc_e,
                              const FIXP_DBL G_m[2], const INT G_e,
                              FIXP_DBL H1_m[2], INT *pH1_e) {
#define H1_SCALE_FACTOR (2)
  const FIXP_DBL max_gain_factor =
      FL2FXCONST_DBL(2.f / (float)(1 << H1_SCALE_FACTOR));

  FIXP_DBL nrgRight_m, nrgLeft_m, crossNrg_m, inv_weight_num_m,
      inv_weight_denom_m, inverse_weight_m, inverse_weight_limited;
  INT nrgRight_e, nrgLeft_e, crossNrg_e, inv_weight_num_e, inv_weight_denom_e,
      inverse_weight_e;

  nrgRight_m = fixpAdd(lin_Cld2_m, lin_Cld2_e, FL2FXCONST_DBL(1.f / 2.f), 1,
                       &nrgRight_e);
  nrgRight_m = invSqrtNorm2(nrgRight_m, nrgRight_e, &nrgRight_e);

  nrgLeft_m = fMult(lin_Cld_m, nrgRight_m);
  nrgLeft_e = lin_Cld_e + nrgRight_e;

  crossNrg_m = sqrtFixp(fMult(nrgLeft_m, nrgRight_m), nrgLeft_e + nrgRight_e,
                        &crossNrg_e);

  inv_weight_num_m =
      fixpAdd(nrgRight_m, nrgRight_e, nrgLeft_m, nrgLeft_e, &inv_weight_num_e);

  inv_weight_denom_m =
      fixpAdd(fMult(fPow2(G_m[0]), nrgLeft_m), 2 * G_e + nrgLeft_e,
              fMult(fPow2(G_m[1]), nrgRight_m), 2 * G_e + nrgRight_e,
              &inv_weight_denom_e);

  inv_weight_denom_m =
      fixpAdd(fMult(fMult(fMult(G_m[0], G_m[1]), crossNrg_m), Icc_m),
              1 + 2 * G_e + crossNrg_e + Icc_e, inv_weight_denom_m,
              inv_weight_denom_e, &inv_weight_denom_e);

  if (inv_weight_denom_m > FL2FXCONST_DBL(0.f)) {
    inverse_weight_m =
        fDivNorm(inv_weight_num_m, inv_weight_denom_m, &inverse_weight_e);
    inverse_weight_m =
        sqrtFixp(inverse_weight_m,
                 inverse_weight_e + inv_weight_num_e - inv_weight_denom_e,
                 &inverse_weight_e);
    inverse_weight_e -= H1_SCALE_FACTOR;

    inverse_weight_limited =
        ((inverse_weight_e >= 0) &&
         (inverse_weight_m >= (max_gain_factor >> inverse_weight_e)))
            ? max_gain_factor
            : scaleValue(inverse_weight_m, inverse_weight_e);
  } else {
    inverse_weight_limited = max_gain_factor;
  }

  H1_m[0] = fMult(G_m[0], inverse_weight_limited);
  H1_m[1] = fMult(G_m[1], inverse_weight_limited);

  *pH1_e = G_e + H1_SCALE_FACTOR;
}

FDK_SACENC_ERROR fdk_sacenc_close_enhancedTimeDomainDmx(
    HANDLE_ENHANCED_TIME_DOMAIN_DMX *phEnhancedTimeDmx) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (phEnhancedTimeDmx == NULL) {
    error = SACENC_INVALID_HANDLE;
  } else {
    if (*phEnhancedTimeDmx != NULL) {
      if ((*phEnhancedTimeDmx)->sinusWindow_m != NULL) {
        FDK_FREE_MEMORY_1D((*phEnhancedTimeDmx)->sinusWindow_m);
      }
      FDK_FREE_MEMORY_1D(*phEnhancedTimeDmx);
    }
  }
  return error;
}
