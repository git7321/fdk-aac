/*********************** MPEG surround encoder library *************************

   Author(s):   Josef Hoepfl

   Description: Encoder Library Interface
                vector functions

*******************************************************************************/

#include "sacenc_vectorfunctions.h"

FIXP_DBL sumUpCplxPow2(const FIXP_DPK *const x, const INT scaleMode,
                       const INT inScaleFactor, INT *const outScaleFactor,
                       const INT n) {
  int i, cs;

  if (scaleMode == SUM_UP_DYNAMIC_SCALE) {
    FIXP_DBL maxVal = FL2FXCONST_DBL(0.0f);
    for (i = 0; i < n; i++) {
      maxVal |= fAbs(x[i].v.re);
      maxVal |= fAbs(x[i].v.im);
    }
    cs = inScaleFactor - fixMax(0, CntLeadingZeros(maxVal) - 1);
  } else {
    cs = inScaleFactor;
  }

  *outScaleFactor = 2 * cs + 2;
  cs = fixMax(fixMin(cs, DFRACT_BITS - 1), -(DFRACT_BITS - 1));
  FIXP_DBL re, im, sum;

  re = im = sum = FL2FXCONST_DBL(0.0);
  if (cs < 0) {
    cs = -cs;
    for (i = 0; i < n; i++) {
      re += fPow2Div2(x[i].v.re << cs);
      im += fPow2Div2(x[i].v.im << cs);
    }
  } else {
    cs = 2 * cs;
    for (i = 0; i < n; i++) {
      re += fPow2Div2(x[i].v.re) >> cs;
      im += fPow2Div2(x[i].v.im) >> cs;
    }
  }

  sum = (re >> 1) + (im >> 1);

  return (sum);
}

FIXP_DBL sumUpCplxPow2Dim2(const FIXP_DPK *const *const x, const INT scaleMode,
                           const INT inScaleFactor, INT *const outScaleFactor,
                           const INT sDim1, const INT nDim1, const INT sDim2,
                           const INT nDim2) {
  int i, j, cs;

  if (scaleMode == SUM_UP_DYNAMIC_SCALE) {
    FIXP_DBL maxVal = FL2FXCONST_DBL(0.0f);
    for (i = sDim1; i < nDim1; i++) {
      for (j = sDim2; j < nDim2; j++) {
        maxVal |= fAbs(x[i][j].v.re);
        maxVal |= fAbs(x[i][j].v.im);
      }
    }
    cs = inScaleFactor - fixMax(0, CntLeadingZeros(maxVal) - 1);
  } else {
    cs = inScaleFactor;
  }

  *outScaleFactor = 2 * cs + 2;

  cs = fixMax(fixMin(cs, DFRACT_BITS - 1), -(DFRACT_BITS - 1));

  FIXP_DBL re, im, sum;

  re = im = sum = FL2FXCONST_DBL(0.0);
  if (cs < 0) {
    cs = -cs;
    for (i = sDim1; i < nDim1; i++) {
      for (j = sDim2; j < nDim2; j++) {
        re += fPow2Div2(x[i][j].v.re << cs);
        im += fPow2Div2(x[i][j].v.im << cs);
      }
    }
  } else {
    cs = 2 * cs;
    for (i = sDim1; i < nDim1; i++) {
      for (j = sDim2; j < nDim2; j++) {
        re += fPow2Div2(x[i][j].v.re) >> cs;
        im += fPow2Div2(x[i][j].v.im) >> cs;
      }
    }
  }

  sum = (re >> 1) + (im >> 1);

  return (sum);
}

void copyCplxVec(FIXP_DPK *const Z, const FIXP_DPK *const X, const INT n) {
  FDKmemmove(Z, X, sizeof(FIXP_DPK) * n);
}

void setCplxVec(FIXP_DPK *const Z, const FIXP_DBL a, const INT n) {
  int i;

  for (i = 0; i < n; i++) {
    Z[i].v.re = a;
    Z[i].v.im = a;
  }
}

void cplx_cplxScalarProduct(FIXP_DPK *const Z, const FIXP_DPK *const *const X,
                            const FIXP_DPK *const *const Y, const INT scaleX,
                            const INT scaleY, INT *const scaleZ,
                            const INT sDim1, const INT nDim1, const INT sDim2,
                            const INT nDim2) {
  int i, j, sx, sy;
  FIXP_DBL xre, yre, xim, yim, re, im;

  sx = fixMax(fixMin(scaleX, DFRACT_BITS - 1), -(DFRACT_BITS - 1));
  sy = fixMax(fixMin(scaleY, DFRACT_BITS - 1), -(DFRACT_BITS - 1));

  *scaleZ = sx + sy + 2;

  re = (FIXP_DBL)0;
  im = (FIXP_DBL)0;
  if ((sx < 0) && (sy < 0)) {
    sx = -sx;
    sy = -sy;
    for (i = sDim1; i < nDim1; i++) {
      for (j = sDim2; j < nDim2; j++) {
        xre = X[i][j].v.re << sx;
        xim = X[i][j].v.im << sx;
        yre = Y[i][j].v.re << sy;
        yim = Y[i][j].v.im << sy;
        re += fMultDiv2(xre, yre) + fMultDiv2(xim, yim);
        im += fMultDiv2(xim, yre) - fMultDiv2(xre, yim);
      }
    }
  } else if ((sx >= 0) && (sy >= 0)) {
    for (i = sDim1; i < nDim1; i++) {
      for (j = sDim2; j < nDim2; j++) {
        xre = X[i][j].v.re;
        xim = X[i][j].v.im;
        yre = Y[i][j].v.re;
        yim = Y[i][j].v.im;
        re += (fMultDiv2(xre, yre) + fMultDiv2(xim, yim)) >> (sx + sy);
        im += (fMultDiv2(xim, yre) - fMultDiv2(xre, yim)) >> (sx + sy);
      }
    }
  } else if ((sx < 0) && (sy >= 0)) {
    sx = -sx;
    for (i = sDim1; i < nDim1; i++) {
      for (j = sDim2; j < nDim2; j++) {
        xre = X[i][j].v.re << sx;
        xim = X[i][j].v.im << sx;
        yre = Y[i][j].v.re;
        yim = Y[i][j].v.im;
        re += (fMultDiv2(xre, yre) + fMultDiv2(xim, yim)) >> sy;
        im += (fMultDiv2(xim, yre) - fMultDiv2(xre, yim)) >> sy;
      }
    }
  } else {
    sy = -sy;
    for (i = sDim1; i < nDim1; i++) {
      for (j = sDim2; j < nDim2; j++) {
        xre = X[i][j].v.re;
        xim = X[i][j].v.im;
        yre = Y[i][j].v.re << sy;
        yim = Y[i][j].v.im << sy;
        re += (fMultDiv2(xre, yre) + fMultDiv2(xim, yim)) >> sx;
        im += (fMultDiv2(xim, yre) - fMultDiv2(xre, yim)) >> sx;
      }
    }
  }

  Z->v.re = re >> 1;
  Z->v.im = im >> 1;
}

void FDKcalcCorrelationVec(FIXP_DBL *const z, const FIXP_DBL *const pr12,
                           const FIXP_DBL *const p1, const FIXP_DBL *const p2,
                           const INT n) {
  int i, s;
  FIXP_DBL p12, cor;

  for (i = 0; i < n; i++) {
    p12 = fMult(p1[i], p2[i]);
    if (p12 > FL2FXCONST_DBL(0.0f)) {
      p12 = invSqrtNorm2(p12, &s);
      cor = fMult(pr12[i], p12);
      z[i] = SATURATE_LEFT_SHIFT(cor, s, DFRACT_BITS);
    } else {
      z[i] = (FIXP_DBL)MAXVAL_DBL;
    }
  }
}

void calcCoherenceVec(FIXP_DBL *const z, const FIXP_DBL *const p12r,
                      const FIXP_DBL *const p12i, const FIXP_DBL *const p1,
                      const FIXP_DBL *const p2, const INT scaleP12,
                      const INT scaleP, const INT n) {
  int i, s, s1, s2;
  FIXP_DBL coh, p12, p12ri;

  for (i = 0; i < n; i++) {
    s2 = fixMin(fixMax(0, CountLeadingBits(p12r[i]) - 1),
                fixMax(0, CountLeadingBits(p12i[i]) - 1));
    p12ri = sqrtFixp(fPow2Div2(p12r[i] << s2) + fPow2Div2(p12i[i] << s2));
    s1 = fixMin(fixMax(0, CountLeadingBits(p1[i]) - 1),
                fixMax(0, CountLeadingBits(p2[i]) - 1));
    p12 = fMultDiv2(p1[i] << s1, p2[i] << s1);

    if (p12 > FL2FXCONST_DBL(0.0f)) {
      p12 = invSqrtNorm2(p12, &s);
      coh = fMult(p12ri, p12);
      s = fixMax(fixMin((scaleP12 - scaleP + s + s1 - s2), DFRACT_BITS - 1),
                 -(DFRACT_BITS - 1));
      if (s < 0) {
        z[i] = coh >> (-s);
      } else {
        z[i] = SATURATE_LEFT_SHIFT(coh, s, DFRACT_BITS);
      }
    } else {
      z[i] = (FIXP_DBL)MAXVAL_DBL;
    }
  }
}

void addWeightedCplxVec(FIXP_DPK *const *const Z, const FIXP_DBL *const a,
                        const FIXP_DPK *const *const X, const FIXP_DBL *const b,
                        const FIXP_DPK *const *const Y, const INT scale,
                        INT *const scaleCh1, const INT scaleCh2,
                        const UCHAR *const pParameterBand2HybridBandOffset,
                        const INT nParameterBands, const INT nTimeSlots,
                        const INT startTimeSlot) {
  int pb, j, i;
  int cs, s1, s2;

  cs = fixMax(*scaleCh1, scaleCh2);
  s1 = cs - (*scaleCh1);
  s2 = cs - scaleCh2;

  *scaleCh1 = cs;

  for (j = 0, pb = 0; pb < nParameterBands; pb++) {
    FIXP_DBL aPb, bPb;
    aPb = a[pb], bPb = b[pb];
    for (; j < pParameterBand2HybridBandOffset[pb]; j++) {
      for (i = startTimeSlot; i < nTimeSlots; i++) {
        Z[j][i].v.re = ((fMultDiv2(aPb, X[j][i].v.re) >> s1) +
                        (fMultDiv2(bPb, Y[j][i].v.re) >> s2))
                       << (scale + 1);
        Z[j][i].v.im = ((fMultDiv2(aPb, X[j][i].v.im) >> s1) +
                        (fMultDiv2(bPb, Y[j][i].v.im) >> s2))
                       << (scale + 1);
      }
    }
  }
}

void FDKcalcPbScaleFactor(const FIXP_DPK *const *const x,
                          const UCHAR *const pParameterBand2HybridBandOffset,
                          INT *const outScaleFactor, const INT startTimeSlot,
                          const INT nTimeSlots, const INT nParamBands) {
  int i, j, pb;

  for (j = 0, pb = 0; pb < nParamBands; pb++) {
    FIXP_DBL maxVal = FL2FXCONST_DBL(0.0f);
    for (; j < pParameterBand2HybridBandOffset[pb]; j++) {
      for (i = startTimeSlot; i < nTimeSlots; i++) {
        maxVal |= fAbs(x[i][j].v.re);
        maxVal |= fAbs(x[i][j].v.im);
      }
    }
    outScaleFactor[pb] = -fixMax(0, CntLeadingZeros(maxVal) - 1);
  }
}

INT FDKcalcScaleFactor(const FIXP_DBL *const x, const FIXP_DBL *const y,
                       const INT n) {
  int i;

  FIXP_DBL maxVal = FL2FXCONST_DBL(0.0f);
  if (x != NULL) {
    for (i = 0; i < n; i++) {
      maxVal |= fAbs(x[i]);
    }
  }

  if (y != NULL) {
    for (i = 0; i < n; i++) {
      maxVal |= fAbs(y[i]);
    }
  }

  if (maxVal == (FIXP_DBL)0)
    return (-(DFRACT_BITS - 1));
  else
    return (-CountLeadingBits(maxVal));
}

INT FDKcalcScaleFactorDPK(const FIXP_DPK *RESTRICT x, const INT startBand,
                          const INT bands) {
  INT qs, clz;
  FIXP_DBL maxVal = FL2FXCONST_DBL(0.0f);

  for (qs = startBand; qs < bands; qs++) {
    maxVal |= fAbs(x[qs].v.re);
    maxVal |= fAbs(x[qs].v.im);
  }

  clz = -fixMax(0, CntLeadingZeros(maxVal) - 1);

  return (clz);
}
