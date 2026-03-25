/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifdef __ANDROID__
#include "log/log.h"
#endif

#include "lpp_tran.h"
#include "sbr_ram.h"
#include "sbr_rom.h"
#include "genericStds.h"
#include "autocorr2nd.h"
#include "HFgen_preFlat.h"

#define LPC_SCALE_FACTOR 2

static FIXP_DBL mapInvfMode(INVF_MODE mode, INVF_MODE prevMode,
                            WHITENING_FACTORS whFactors) {
  switch (mode) {
    case INVF_LOW_LEVEL:
      if (prevMode == INVF_OFF)
        return whFactors.transitionLevel;
      else
        return whFactors.lowLevel;

    case INVF_MID_LEVEL:
      return whFactors.midLevel;

    case INVF_HIGH_LEVEL:
      return whFactors.highLevel;

    default:
      if (prevMode == INVF_LOW_LEVEL)
        return whFactors.transitionLevel;
      else
        return whFactors.off;
  }
}

static void inverseFilteringLevelEmphasis(
    HANDLE_SBR_LPP_TRANS hLppTrans,
    UCHAR nInvfBands,
    INVF_MODE *sbr_invf_mode,
    INVF_MODE *sbr_invf_mode_prev,
    FIXP_DBL *bwVector
) {
  for (int i = 0; i < nInvfBands; i++) {
    FIXP_DBL accu;
    FIXP_DBL bwTmp = mapInvfMode(sbr_invf_mode[i], sbr_invf_mode_prev[i],
                                 hLppTrans->pSettings->whFactors);

    if (bwTmp < hLppTrans->bwVectorOld[i]) {
      accu = fMultDiv2(FL2FXCONST_DBL(0.75f), bwTmp) +
             fMultDiv2(FL2FXCONST_DBL(0.25f), hLppTrans->bwVectorOld[i]);
    } else {
      accu = fMultDiv2(FL2FXCONST_DBL(0.90625f), bwTmp) +
             fMultDiv2(FL2FXCONST_DBL(0.09375f), hLppTrans->bwVectorOld[i]);
    }

    if (accu<FL2FXCONST_DBL(0.015625f)>> 1) {
      bwVector[i] = FL2FXCONST_DBL(0.0f);
    } else {
      bwVector[i] = fixMin(accu << 1, FL2FXCONST_DBL(0.99609375f));
    }
  }
}

#define ACDET_EXP \
  (2 * (DFRACT_BITS + sbrScaleFactor->lb_scale + 10 - ac.det_scale))
#define AC_EXP (-sbrScaleFactor->lb_scale + LPC_SCALE_FACTOR)
#define ALPHA_EXP (-sbrScaleFactor->lb_scale + LPC_SCALE_FACTOR + 1)
#define QMFOUT_EXP ((SAMPLE_BITS - 15) - sbrScaleFactor->lb_scale)

static inline void calc_qmfBufferReal(FIXP_DBL **qmfBufferReal,
                                      const FIXP_DBL *const lowBandReal,
                                      const int startSample,
                                      const int stopSample, const UCHAR hiBand,
                                      const int dynamicScale,
                                      const FIXP_SGL a0r, const FIXP_SGL a1r) {
  const int dynscale = fixMax(0, dynamicScale - 1) + 1;
  const int rescale = -fixMin(0, dynamicScale - 1) + 1;
  const int descale =
      fixMin(DFRACT_BITS - 1, LPC_SCALE_FACTOR + dynamicScale + rescale);

  for (int i = 0; i < stopSample - startSample; i++) {
    FIXP_DBL accu;

    accu = fMultDiv2(a1r, lowBandReal[i]) + fMultDiv2(a0r, lowBandReal[i + 1]);
    accu = (lowBandReal[i + 2] >> descale) + (accu >> dynscale);

    qmfBufferReal[i + startSample][hiBand] =
        SATURATE_LEFT_SHIFT(accu, rescale, DFRACT_BITS);
  }
}

void lppTransposer(
    HANDLE_SBR_LPP_TRANS hLppTrans,
    QMF_SCALE_FACTOR *sbrScaleFactor,
    FIXP_DBL **qmfBufferReal,

    FIXP_DBL *degreeAlias,
    FIXP_DBL **qmfBufferImag,
    const int useLP, const int fPreWhitening, const int v_k_master0,
    const int timeStep,
    const int firstSlotOffs,
    const int lastSlotOffs,
    const int nInvfBands,
    INVF_MODE *sbr_invf_mode,
    INVF_MODE *sbr_invf_mode_prev
) {
  INT bwIndex[MAX_NUM_PATCHES];
  FIXP_DBL bwVector[MAX_NUM_PATCHES];
  FIXP_DBL preWhiteningGains[(64) / 2];
  int preWhiteningGains_exp[(64) / 2];

  int i;
  int loBand, start, stop;
  TRANSPOSER_SETTINGS *pSettings = hLppTrans->pSettings;
  PATCH_PARAM *patchParam = pSettings->patchParam;
  int patch;

  FIXP_SGL alphar[LPC_ORDER], a0r, a1r;
  FIXP_SGL alphai[LPC_ORDER], a0i = 0, a1i = 0;
  FIXP_SGL bw = FL2FXCONST_SGL(0.0f);

  int autoCorrLength;

  FIXP_DBL k1, k1_below = 0, k1_below2 = 0;

  ACORR_COEFS ac;
  int startSample;
  int stopSample;
  int stopSampleClear;

  int comLowBandScale;
  int ovLowBandShift;
  int lowBandShift;

  alphai[0] = FL2FXCONST_SGL(0.0f);
  alphai[1] = FL2FXCONST_SGL(0.0f);

  startSample = firstSlotOffs * timeStep;
  stopSample = pSettings->nCols + lastSlotOffs * timeStep;
  FDK_ASSERT((lastSlotOffs * timeStep) <= pSettings->overlap);

  inverseFilteringLevelEmphasis(hLppTrans, nInvfBands, sbr_invf_mode,
                                sbr_invf_mode_prev, bwVector);

  stopSampleClear = stopSample;

  autoCorrLength = pSettings->nCols + pSettings->overlap;

  if (pSettings->noOfPatches > 0) {
    int targetStopBand =
        patchParam[pSettings->noOfPatches - 1].targetStartBand +
        patchParam[pSettings->noOfPatches - 1].numBandsInPatch;

    int memSize = ((64) - targetStopBand) * sizeof(FIXP_DBL);

    if (!useLP) {
      for (i = startSample; i < stopSampleClear; i++) {
        FDKmemclear(&qmfBufferReal[i][targetStopBand], memSize);
        FDKmemclear(&qmfBufferImag[i][targetStopBand], memSize);
      }
    } else {
      for (i = startSample; i < stopSampleClear; i++) {
        FDKmemclear(&qmfBufferReal[i][targetStopBand], memSize);
      }
    }
  }
#ifdef __ANDROID__
  else {
    android_errorWriteLog(0x534e4554, "112160868");
  }
#endif

  FDKmemclear(bwIndex, sizeof(bwIndex));

  comLowBandScale =
      fixMin(sbrScaleFactor->ov_lb_scale, sbrScaleFactor->lb_scale);

  ovLowBandShift = sbrScaleFactor->ov_lb_scale - comLowBandScale;
  lowBandShift = sbrScaleFactor->lb_scale - comLowBandScale;

  if (fPreWhitening) {
    sbrDecoder_calculateGainVec(
        qmfBufferReal, qmfBufferImag,
        DFRACT_BITS - 1 - 16 -
            sbrScaleFactor->ov_lb_scale,
        DFRACT_BITS - 1 - 16 -
            sbrScaleFactor->lb_scale,
        pSettings->overlap, preWhiteningGains, preWhiteningGains_exp,
        v_k_master0, startSample, stopSample);
  }

  if (!useLP) {
    start = pSettings->lbStartPatching;
    stop = pSettings->lbStopPatching;
  } else {
    start = fixMax(1, pSettings->lbStartPatching - 2);
    stop = patchParam[0].targetStartBand;
  }

  for (loBand = start; loBand < stop; loBand++) {
    FIXP_DBL lowBandReal[(((1024) / (32) * (4) / 2) + (3 * (4))) + LPC_ORDER];
    FIXP_DBL *plowBandReal = lowBandReal;
    FIXP_DBL **pqmfBufferReal =
        qmfBufferReal + firstSlotOffs * timeStep;
    FIXP_DBL lowBandImag[(((1024) / (32) * (4) / 2) + (3 * (4))) + LPC_ORDER];
    FIXP_DBL *plowBandImag = lowBandImag;
    FIXP_DBL **pqmfBufferImag =
        qmfBufferImag + firstSlotOffs * timeStep;
    int resetLPCCoeffs = 0;
    int dynamicScale = DFRACT_BITS - 1 - LPC_SCALE_FACTOR;
    int acDetScale = 0;

    for (i = 0;
         i < LPC_ORDER + firstSlotOffs * timeStep;
         i++) {
      *plowBandReal++ = hLppTrans->lpcFilterStatesRealLegSBR[i][loBand];
      if (!useLP)
        *plowBandImag++ = hLppTrans->lpcFilterStatesImagLegSBR[i][loBand];
    }

    if (!useLP) {
      for (i = 0;
           i < pSettings->nCols + pSettings->overlap - firstSlotOffs * timeStep;
           i++) {
        *plowBandReal++ = (*pqmfBufferReal++)[loBand];
        *plowBandImag++ = (*pqmfBufferImag++)[loBand];
      }
    } else {
      FDK_ASSERT((pSettings->overlap & 1) == 0);
      for (i = 0; i < ((pSettings->nCols + pSettings->overlap -
                        firstSlotOffs * timeStep) >>
                       1);
           i++) {
        *plowBandReal++ = (*pqmfBufferReal++)[loBand];
        *plowBandReal++ = (*pqmfBufferReal++)[loBand];
      }
      if (pSettings->nCols & 1) {
        *plowBandReal++ = (*pqmfBufferReal++)[loBand];
      }
    }

    dynamicScale =
        fixMin(dynamicScale,
               getScalefactor(lowBandReal, LPC_ORDER + pSettings->overlap) +
                   ovLowBandShift);
    dynamicScale =
        fixMin(dynamicScale,
               getScalefactor(&lowBandReal[LPC_ORDER + pSettings->overlap],
                              pSettings->nCols) +
                   lowBandShift);
    if (!useLP) {
      dynamicScale =
          fixMin(dynamicScale,
                 getScalefactor(lowBandImag, LPC_ORDER + pSettings->overlap) +
                     ovLowBandShift);
      dynamicScale =
          fixMin(dynamicScale,
                 getScalefactor(&lowBandImag[LPC_ORDER + pSettings->overlap],
                                pSettings->nCols) +
                     lowBandShift);
    }

    if (dynamicScale == 0) {
      for (i = 0; i < (LPC_ORDER + pSettings->overlap + pSettings->nCols);
           i++) {
        lowBandReal[i] = fixMax(lowBandReal[i], (FIXP_DBL)0x80000001);
      }
      if (!useLP) {
        for (i = 0; i < (LPC_ORDER + pSettings->overlap + pSettings->nCols);
             i++) {
          lowBandImag[i] = fixMax(lowBandImag[i], (FIXP_DBL)0x80000001);
        }
      }
    } else {
      dynamicScale =
          fixMax(0, dynamicScale -
                        1);
    }

    scaleValues(&lowBandReal[0], LPC_ORDER + pSettings->overlap,
                dynamicScale - ovLowBandShift);
    scaleValues(&lowBandReal[LPC_ORDER + pSettings->overlap], pSettings->nCols,
                dynamicScale - lowBandShift);

    if (!useLP) {
      scaleValues(&lowBandImag[0], LPC_ORDER + pSettings->overlap,
                  dynamicScale - ovLowBandShift);
      scaleValues(&lowBandImag[LPC_ORDER + pSettings->overlap],
                  pSettings->nCols, dynamicScale - lowBandShift);
    }

    if (!useLP) {
      acDetScale += autoCorr2nd_cplx(&ac, lowBandReal + LPC_ORDER,
                                     lowBandImag + LPC_ORDER, autoCorrLength);
    } else {
      acDetScale +=
          autoCorr2nd_real(&ac, lowBandReal + LPC_ORDER, autoCorrLength);
    }

    acDetScale += 2 * (comLowBandScale + dynamicScale);
    acDetScale *= 2;
    acDetScale += ac.det_scale;

    if (acDetScale > 126) {
      resetLPCCoeffs = 1;
    }

    alphar[1] = FL2FXCONST_SGL(0.0f);
    if (!useLP) alphai[1] = FL2FXCONST_SGL(0.0f);

    if (ac.det != FL2FXCONST_DBL(0.0f)) {
      FIXP_DBL tmp, absTmp, absDet;

      absDet = fixp_abs(ac.det);

      if (!useLP) {
        tmp = (fMultDiv2(ac.r01r, ac.r12r) >> (LPC_SCALE_FACTOR - 1)) -
              ((fMultDiv2(ac.r01i, ac.r12i) + fMultDiv2(ac.r02r, ac.r11r)) >>
               (LPC_SCALE_FACTOR - 1));
      } else {
        tmp = (fMultDiv2(ac.r01r, ac.r12r) >> (LPC_SCALE_FACTOR - 1)) -
              (fMultDiv2(ac.r02r, ac.r11r) >> (LPC_SCALE_FACTOR - 1));
      }
      absTmp = fixp_abs(tmp);

      {
        INT scale;
        FIXP_DBL result = fDivNorm(absTmp, absDet, &scale);
        scale = scale + ac.det_scale;

        if ((scale > 0) && (result >= (FIXP_DBL)MAXVAL_DBL >> scale)) {
          resetLPCCoeffs = 1;
        } else {
          alphar[1] = FX_DBL2FX_SGL(scaleValueSaturate(result, scale));
          if ((tmp < FL2FX_DBL(0.0f)) ^ (ac.det < FL2FX_DBL(0.0f))) {
            alphar[1] = -alphar[1];
          }
        }
      }

      if (!useLP) {
        tmp = (fMultDiv2(ac.r01i, ac.r12r) >> (LPC_SCALE_FACTOR - 1)) +
              ((fMultDiv2(ac.r01r, ac.r12i) -
                (FIXP_DBL)fMultDiv2(ac.r02i, ac.r11r)) >>
               (LPC_SCALE_FACTOR - 1));

        absTmp = fixp_abs(tmp);

        {
          INT scale;
          FIXP_DBL result = fDivNorm(absTmp, absDet, &scale);
          scale = scale + ac.det_scale;

          if ((scale > 0) &&
              (result >= (FIXP_DBL)MAXVAL_DBL >>
               scale)) {
            resetLPCCoeffs = 1;
          } else {
            alphai[1] = FX_DBL2FX_SGL(scaleValueSaturate(result, scale));
            if ((tmp < FL2FX_DBL(0.0f)) ^ (ac.det < FL2FX_DBL(0.0f))) {
              alphai[1] = -alphai[1];
            }
          }
        }
      }
    }

    alphar[0] = FL2FXCONST_SGL(0.0f);
    if (!useLP) alphai[0] = FL2FXCONST_SGL(0.0f);

    if (ac.r11r != FL2FXCONST_DBL(0.0f)) {
      FIXP_DBL tmp, absTmp;

      if (!useLP) {
        tmp = (ac.r01r >> (LPC_SCALE_FACTOR + 1)) +
              (fMultDiv2(alphar[1], ac.r12r) + fMultDiv2(alphai[1], ac.r12i));
      } else {
        if (ac.r01r >= FL2FXCONST_DBL(0.0f))
          tmp = (ac.r01r >> (LPC_SCALE_FACTOR + 1)) +
                fMultDiv2(alphar[1], ac.r12r);
        else
          tmp = -((-ac.r01r) >> (LPC_SCALE_FACTOR + 1)) +
                fMultDiv2(alphar[1], ac.r12r);
      }

      absTmp = fixp_abs(tmp);

      if (absTmp >= (ac.r11r >> 1)) {
        resetLPCCoeffs = 1;
      } else {
        INT scale;
        FIXP_DBL result = fDivNorm(absTmp, fixp_abs(ac.r11r), &scale);
        alphar[0] = FX_DBL2FX_SGL(scaleValueSaturate(result, scale + 1));

        if ((tmp > FL2FX_DBL(0.0f)) ^ (ac.r11r < FL2FX_DBL(0.0f)))
          alphar[0] = -alphar[0];
      }

      if (!useLP) {
        tmp = (ac.r01i >> (LPC_SCALE_FACTOR + 1)) +
              (fMultDiv2(alphai[1], ac.r12r) - fMultDiv2(alphar[1], ac.r12i));

        absTmp = fixp_abs(tmp);

        if (absTmp >= (ac.r11r >> 1)) {
          resetLPCCoeffs = 1;
        } else {
          INT scale;
          FIXP_DBL result = fDivNorm(absTmp, fixp_abs(ac.r11r), &scale);
          alphai[0] = FX_DBL2FX_SGL(scaleValueSaturate(result, scale + 1));
          if ((tmp > FL2FX_DBL(0.0f)) ^ (ac.r11r < FL2FX_DBL(0.0f)))
            alphai[0] = -alphai[0];
        }
      }
    }

    if (!useLP) {
      if ((fMultDiv2(alphar[0], alphar[0]) + fMultDiv2(alphai[0], alphai[0])) >=
          FL2FXCONST_DBL(0.5f))
        resetLPCCoeffs = 1;
      if ((fMultDiv2(alphar[1], alphar[1]) + fMultDiv2(alphai[1], alphai[1])) >=
          FL2FXCONST_DBL(0.5f))
        resetLPCCoeffs = 1;
    }

    if (resetLPCCoeffs) {
      alphar[0] = FL2FXCONST_SGL(0.0f);
      alphar[1] = FL2FXCONST_SGL(0.0f);
      if (!useLP) {
        alphai[0] = FL2FXCONST_SGL(0.0f);
        alphai[1] = FL2FXCONST_SGL(0.0f);
      }
    }

    if (useLP) {
      if (ac.r11r == FL2FXCONST_DBL(0.0f)) {
        k1 = FL2FXCONST_DBL(0.0f);
      } else {
        if (fixp_abs(ac.r01r) >= fixp_abs(ac.r11r)) {
          if (fMultDiv2(ac.r01r, ac.r11r) < FL2FX_DBL(0.0f)) {
            k1 = (FIXP_DBL)MAXVAL_DBL;
          } else {
            k1 = (FIXP_DBL)(MINVAL_DBL + 1);
          }
        } else {
          INT scale;
          FIXP_DBL result =
              fDivNorm(fixp_abs(ac.r01r), fixp_abs(ac.r11r), &scale);
          k1 = scaleValueSaturate(result, scale);

          if (!((ac.r01r < FL2FX_DBL(0.0f)) ^ (ac.r11r < FL2FX_DBL(0.0f)))) {
            k1 = -k1;
          }
        }
      }
      if ((loBand > 1) && (loBand < v_k_master0)) {
        FIXP_DBL deg =
            (FIXP_DBL)MAXVAL_DBL - fPow2(k1_below);
        degreeAlias[loBand] = FL2FXCONST_DBL(0.0f);
        if (((loBand & 1) == 0) && (k1 < FL2FXCONST_DBL(0.0f))) {
          if (k1_below < FL2FXCONST_DBL(0.0f)) {
            degreeAlias[loBand] = (FIXP_DBL)MAXVAL_DBL;
            if (k1_below2 >
                FL2FXCONST_DBL(0.0f)) {
              degreeAlias[loBand - 1] = deg;
            }
          } else if (k1_below2 >
                     FL2FXCONST_DBL(0.0f)) {
            degreeAlias[loBand] = deg;
          }
        }
        if (((loBand & 1) == 1) && (k1 > FL2FXCONST_DBL(0.0f))) {
          if (k1_below > FL2FXCONST_DBL(0.0f)) {
            degreeAlias[loBand] = (FIXP_DBL)MAXVAL_DBL;
            if (k1_below2 <
                FL2FXCONST_DBL(0.0f)) {
              degreeAlias[loBand - 1] = deg;
            }
          } else if (k1_below2 <
                     FL2FXCONST_DBL(0.0f)) {
            degreeAlias[loBand] = deg;
          }
        }
      }
      k1_below2 = k1_below;
      k1_below = k1;
    }

    patch = 0;

    while (patch < pSettings->noOfPatches) {

      int hiBand = loBand + patchParam[patch].targetBandOffs;

      if (loBand < patchParam[patch].sourceStartBand ||
          loBand >= patchParam[patch].sourceStopBand
      ) {
        patch++;
        continue;
      }

      FDK_ASSERT(hiBand < (64));

      while (hiBand >= pSettings->bwBorders[bwIndex[patch]] &&
             bwIndex[patch] < MAX_NUM_PATCHES - 1) {
        bwIndex[patch]++;
      }

      bw = FX_DBL2FX_SGL(bwVector[bwIndex[patch]]);

      a0r = FX_DBL2FX_SGL(
          fMult(bw, alphar[0]));

      if (!useLP) a0i = FX_DBL2FX_SGL(fMult(bw, alphai[0]));
      bw = FX_DBL2FX_SGL(fPow2(bw));
      a1r = FX_DBL2FX_SGL(fMult(bw, alphar[1]));
      if (!useLP) a1i = FX_DBL2FX_SGL(fMult(bw, alphai[1]));

      if (bw <= FL2FXCONST_SGL(0.0f)) {
        if (!useLP) {
          int descale =
              fixMin(DFRACT_BITS - 1, (LPC_SCALE_FACTOR + dynamicScale));
          for (i = startSample; i < stopSample; i++) {
            FIXP_DBL accu1, accu2;
            accu1 = lowBandReal[LPC_ORDER + i] >> descale;
            accu2 = lowBandImag[LPC_ORDER + i] >> descale;
            if (fPreWhitening) {
              accu1 = scaleValueSaturate(
                  fMultDiv2(accu1, preWhiteningGains[loBand]),
                  preWhiteningGains_exp[loBand] + 1);
              accu2 = scaleValueSaturate(
                  fMultDiv2(accu2, preWhiteningGains[loBand]),
                  preWhiteningGains_exp[loBand] + 1);
            }
            qmfBufferReal[i][hiBand] = accu1;
            qmfBufferImag[i][hiBand] = accu2;
          }
        } else {
          int descale =
              fixMin(DFRACT_BITS - 1, (LPC_SCALE_FACTOR + dynamicScale));
          for (i = startSample; i < stopSample; i++) {
            qmfBufferReal[i][hiBand] = lowBandReal[LPC_ORDER + i] >> descale;
          }
        }
      } else {

        if (!useLP) {
          const int dynscale = fixMax(0, dynamicScale - 2) + 1;
          const int rescale = -fixMin(0, dynamicScale - 2) + 1;
          const int descale = fixMin(DFRACT_BITS - 1,
                                     LPC_SCALE_FACTOR + dynamicScale + rescale);

          for (i = startSample; i < stopSample; i++) {
            FIXP_DBL accu1, accu2;

            accu1 = ((fMultDiv2(a0r, lowBandReal[LPC_ORDER + i - 1]) -
                      fMultDiv2(a0i, lowBandImag[LPC_ORDER + i - 1])) >>
                     1) +
                    ((fMultDiv2(a1r, lowBandReal[LPC_ORDER + i - 2]) -
                      fMultDiv2(a1i, lowBandImag[LPC_ORDER + i - 2])) >>
                     1);
            accu2 = ((fMultDiv2(a0i, lowBandReal[LPC_ORDER + i - 1]) +
                      fMultDiv2(a0r, lowBandImag[LPC_ORDER + i - 1])) >>
                     1) +
                    ((fMultDiv2(a1i, lowBandReal[LPC_ORDER + i - 2]) +
                      fMultDiv2(a1r, lowBandImag[LPC_ORDER + i - 2])) >>
                     1);

            accu1 =
                (lowBandReal[LPC_ORDER + i] >> descale) + (accu1 >> dynscale);
            accu2 =
                (lowBandImag[LPC_ORDER + i] >> descale) + (accu2 >> dynscale);
            if (fPreWhitening) {
              qmfBufferReal[i][hiBand] = scaleValueSaturate(
                  fMultDiv2(accu1, preWhiteningGains[loBand]),
                  preWhiteningGains_exp[loBand] + 1 + rescale);
              qmfBufferImag[i][hiBand] = scaleValueSaturate(
                  fMultDiv2(accu2, preWhiteningGains[loBand]),
                  preWhiteningGains_exp[loBand] + 1 + rescale);
            } else {
              qmfBufferReal[i][hiBand] =
                  SATURATE_LEFT_SHIFT(accu1, rescale, DFRACT_BITS);
              qmfBufferImag[i][hiBand] =
                  SATURATE_LEFT_SHIFT(accu2, rescale, DFRACT_BITS);
            }
          }
        } else {
          FDK_ASSERT(dynamicScale >= 0);
          calc_qmfBufferReal(
              qmfBufferReal, &(lowBandReal[LPC_ORDER + startSample - 2]),
              startSample, stopSample, hiBand, dynamicScale, a0r, a1r);
        }
      }

      patch++;
    }
  }

  if (useLP) {
    for (loBand = pSettings->lbStartPatching;
         loBand < pSettings->lbStopPatching; loBand++) {
      patch = 0;
      while (patch < pSettings->noOfPatches) {
        UCHAR hiBand = loBand + patchParam[patch].targetBandOffs;

        if (loBand < patchParam[patch].sourceStartBand ||
            loBand >= patchParam[patch].sourceStopBand ||
            hiBand >= (64)
        ) {
          patch++;
          continue;
        }

        if (hiBand != patchParam[patch].targetStartBand)
          degreeAlias[hiBand] = degreeAlias[loBand];

        patch++;
      }
    }
  }

  for (i = 0; i < nInvfBands; i++) {
    hLppTrans->bwVectorOld[i] = bwVector[i];
  }

  sbrScaleFactor->hb_scale = comLowBandScale - (LPC_SCALE_FACTOR);
}

void lppTransposerHBE(
    HANDLE_SBR_LPP_TRANS hLppTrans,
    HANDLE_HBE_TRANSPOSER hQmfTransposer,
    QMF_SCALE_FACTOR *sbrScaleFactor,
    FIXP_DBL **qmfBufferReal,
    FIXP_DBL **qmfBufferImag,
    const int timeStep,
    const int firstSlotOffs,
    const int lastSlotOffs,
    const int nInvfBands,
    INVF_MODE *sbr_invf_mode,
    INVF_MODE *sbr_invf_mode_prev
) {
  INT bwIndex;
  FIXP_DBL bwVector[MAX_NUM_PATCHES_HBE];

  int i;
  int loBand, start, stop;
  TRANSPOSER_SETTINGS *pSettings = hLppTrans->pSettings;
  PATCH_PARAM *patchParam = pSettings->patchParam;

  FIXP_SGL alphar[LPC_ORDER], a0r, a1r;
  FIXP_SGL alphai[LPC_ORDER], a0i = 0, a1i = 0;
  FIXP_SGL bw = FL2FXCONST_SGL(0.0f);

  int autoCorrLength;

  ACORR_COEFS ac;
  int startSample;
  int stopSample;
  int stopSampleClear;

  int comBandScale;
  int ovLowBandShift;
  int lowBandShift;

  alphai[0] = FL2FXCONST_SGL(0.0f);
  alphai[1] = FL2FXCONST_SGL(0.0f);

  startSample = firstSlotOffs * timeStep;
  stopSample = pSettings->nCols + lastSlotOffs * timeStep;

  inverseFilteringLevelEmphasis(hLppTrans, nInvfBands, sbr_invf_mode,
                                sbr_invf_mode_prev, bwVector);

  stopSampleClear = stopSample;

  autoCorrLength = pSettings->nCols + pSettings->overlap;

  if (pSettings->noOfPatches > 0) {
    int targetStopBand =
        patchParam[pSettings->noOfPatches - 1].targetStartBand +
        patchParam[pSettings->noOfPatches - 1].numBandsInPatch;

    int memSize = ((64) - targetStopBand) * sizeof(FIXP_DBL);

    for (i = startSample; i < stopSampleClear; i++) {
      FDKmemclear(&qmfBufferReal[i][targetStopBand], memSize);
      FDKmemclear(&qmfBufferImag[i][targetStopBand], memSize);
    }
  }
#ifdef __ANDROID__
  else {
    android_errorWriteLog(0x534e4554, "112160868");
  }
#endif

  comBandScale = sbrScaleFactor->hb_scale;

  ovLowBandShift = sbrScaleFactor->hb_scale - comBandScale;
  lowBandShift = sbrScaleFactor->hb_scale - comBandScale;

  start = hQmfTransposer->startBand;
  stop = hQmfTransposer->stopBand;

  for (loBand = start; loBand < stop; loBand++) {
    bwIndex = 0;

    FIXP_DBL lowBandReal[(((1024) / (32) * (4) / 2) + (3 * (4))) + LPC_ORDER];
    FIXP_DBL lowBandImag[(((1024) / (32) * (4) / 2) + (3 * (4))) + LPC_ORDER];

    int resetLPCCoeffs = 0;
    int dynamicScale = DFRACT_BITS - 1 - LPC_SCALE_FACTOR;
    int acDetScale = 0;

    for (i = 0; i < LPC_ORDER; i++) {
      lowBandReal[i] = hLppTrans->lpcFilterStatesRealHBE[i][loBand];
      lowBandImag[i] = hLppTrans->lpcFilterStatesImagHBE[i][loBand];
    }

    for (; i < LPC_ORDER + firstSlotOffs * timeStep; i++) {
      lowBandReal[i] = hLppTrans->lpcFilterStatesRealHBE[i][loBand];
      lowBandImag[i] = hLppTrans->lpcFilterStatesImagHBE[i][loBand];
    }

    for (i = firstSlotOffs * timeStep;
         i < pSettings->nCols + pSettings->overlap; i++) {
      lowBandReal[i + LPC_ORDER] = qmfBufferReal[i][loBand];
      lowBandImag[i + LPC_ORDER] = qmfBufferImag[i][loBand];
    }

    for (i = 0; i < LPC_ORDER + pSettings->overlap; i++) {
      hLppTrans->lpcFilterStatesRealHBE[i][loBand] =
          qmfBufferReal[pSettings->nCols - LPC_ORDER + i][loBand];
      hLppTrans->lpcFilterStatesImagHBE[i][loBand] =
          qmfBufferImag[pSettings->nCols - LPC_ORDER + i][loBand];
    }

    dynamicScale =
        fixMin(dynamicScale,
               getScalefactor(lowBandReal, LPC_ORDER + pSettings->overlap) +
                   ovLowBandShift);
    dynamicScale =
        fixMin(dynamicScale,
               getScalefactor(&lowBandReal[LPC_ORDER + pSettings->overlap],
                              pSettings->nCols) +
                   lowBandShift);
    dynamicScale =
        fixMin(dynamicScale,
               getScalefactor(lowBandImag, LPC_ORDER + pSettings->overlap) +
                   ovLowBandShift);
    dynamicScale =
        fixMin(dynamicScale,
               getScalefactor(&lowBandImag[LPC_ORDER + pSettings->overlap],
                              pSettings->nCols) +
                   lowBandShift);

    dynamicScale =
        dynamicScale - 1;

    scaleValues(&lowBandReal[0], LPC_ORDER + pSettings->overlap,
                dynamicScale - ovLowBandShift);
    scaleValues(&lowBandReal[LPC_ORDER + pSettings->overlap], pSettings->nCols,
                dynamicScale - lowBandShift);
    scaleValues(&lowBandImag[0], LPC_ORDER + pSettings->overlap,
                dynamicScale - ovLowBandShift);
    scaleValues(&lowBandImag[LPC_ORDER + pSettings->overlap], pSettings->nCols,
                dynamicScale - lowBandShift);

    acDetScale += autoCorr2nd_cplx(&ac, lowBandReal + LPC_ORDER,
                                   lowBandImag + LPC_ORDER, autoCorrLength);

    acDetScale += 2 * (comBandScale + dynamicScale);
    acDetScale *= 2;
    acDetScale += ac.det_scale;

    if (acDetScale > 126) {
      resetLPCCoeffs = 1;
    }

    alphar[1] = FL2FXCONST_SGL(0.0f);
    alphai[1] = FL2FXCONST_SGL(0.0f);

    if (ac.det != FL2FXCONST_DBL(0.0f)) {
      FIXP_DBL tmp, absTmp, absDet;

      absDet = fixp_abs(ac.det);

      tmp = (fMultDiv2(ac.r01r, ac.r12r) >> (LPC_SCALE_FACTOR - 1)) -
            ((fMultDiv2(ac.r01i, ac.r12i) + fMultDiv2(ac.r02r, ac.r11r)) >>
             (LPC_SCALE_FACTOR - 1));
      absTmp = fixp_abs(tmp);

      {
        INT scale;
        FIXP_DBL result = fDivNorm(absTmp, absDet, &scale);
        scale = scale + ac.det_scale;

        if ((scale > 0) && (result >= (FIXP_DBL)MAXVAL_DBL >> scale)) {
          resetLPCCoeffs = 1;
        } else {
          alphar[1] = FX_DBL2FX_SGL(scaleValueSaturate(result, scale));
          if ((tmp < FL2FX_DBL(0.0f)) ^ (ac.det < FL2FX_DBL(0.0f))) {
            alphar[1] = -alphar[1];
          }
        }
      }

      tmp = (fMultDiv2(ac.r01i, ac.r12r) >> (LPC_SCALE_FACTOR - 1)) +
            ((fMultDiv2(ac.r01r, ac.r12i) -
              (FIXP_DBL)fMultDiv2(ac.r02i, ac.r11r)) >>
             (LPC_SCALE_FACTOR - 1));

      absTmp = fixp_abs(tmp);

      {
        INT scale;
        FIXP_DBL result = fDivNorm(absTmp, absDet, &scale);
        scale = scale + ac.det_scale;

        if ((scale > 0) &&
            (result >= (FIXP_DBL)MAXVAL_DBL >> scale)) {
          resetLPCCoeffs = 1;
        } else {
          alphai[1] = FX_DBL2FX_SGL(scaleValueSaturate(result, scale));
          if ((tmp < FL2FX_DBL(0.0f)) ^ (ac.det < FL2FX_DBL(0.0f))) {
            alphai[1] = -alphai[1];
          }
        }
      }
    }

    alphar[0] = FL2FXCONST_SGL(0.0f);
    alphai[0] = FL2FXCONST_SGL(0.0f);

    if (ac.r11r != FL2FXCONST_DBL(0.0f)) {
      FIXP_DBL tmp, absTmp;

      tmp = (ac.r01r >> (LPC_SCALE_FACTOR + 1)) +
            (fMultDiv2(alphar[1], ac.r12r) + fMultDiv2(alphai[1], ac.r12i));

      absTmp = fixp_abs(tmp);

      if (absTmp >= (ac.r11r >> 1)) {
        resetLPCCoeffs = 1;
      } else {
        INT scale;
        FIXP_DBL result = fDivNorm(absTmp, fixp_abs(ac.r11r), &scale);
        alphar[0] = FX_DBL2FX_SGL(scaleValueSaturate(result, scale + 1));

        if ((tmp > FL2FX_DBL(0.0f)) ^ (ac.r11r < FL2FX_DBL(0.0f)))
          alphar[0] = -alphar[0];
      }

      tmp = (ac.r01i >> (LPC_SCALE_FACTOR + 1)) +
            (fMultDiv2(alphai[1], ac.r12r) - fMultDiv2(alphar[1], ac.r12i));

      absTmp = fixp_abs(tmp);

      if (absTmp >= (ac.r11r >> 1)) {
        resetLPCCoeffs = 1;
      } else {
        INT scale;
        FIXP_DBL result = fDivNorm(absTmp, fixp_abs(ac.r11r), &scale);
        alphai[0] = FX_DBL2FX_SGL(scaleValueSaturate(result, scale + 1));
        if ((tmp > FL2FX_DBL(0.0f)) ^ (ac.r11r < FL2FX_DBL(0.0f))) {
          alphai[0] = -alphai[0];
        }
      }
    }

    if ((fMultDiv2(alphar[0], alphar[0]) + fMultDiv2(alphai[0], alphai[0])) >=
        FL2FXCONST_DBL(0.5f)) {
      resetLPCCoeffs = 1;
    }
    if ((fMultDiv2(alphar[1], alphar[1]) + fMultDiv2(alphai[1], alphai[1])) >=
        FL2FXCONST_DBL(0.5f)) {
      resetLPCCoeffs = 1;
    }

    if (resetLPCCoeffs) {
      alphar[0] = FL2FXCONST_SGL(0.0f);
      alphar[1] = FL2FXCONST_SGL(0.0f);
      alphai[0] = FL2FXCONST_SGL(0.0f);
      alphai[1] = FL2FXCONST_SGL(0.0f);
    }

    while (bwIndex < MAX_NUM_PATCHES - 1 &&
           loBand >= pSettings->bwBorders[bwIndex]) {
      bwIndex++;
    }

    bw = FX_DBL2FX_SGL(bwVector[bwIndex]);

    a0r = FX_DBL2FX_SGL(
        fMult(bw, alphar[0]));
    a0i = FX_DBL2FX_SGL(fMult(bw, alphai[0]));
    bw = FX_DBL2FX_SGL(fPow2(bw));
    a1r = FX_DBL2FX_SGL(fMult(bw, alphar[1]));
    a1i = FX_DBL2FX_SGL(fMult(bw, alphai[1]));

    if (bw <= FL2FXCONST_SGL(0.0f)) {
      int descale = fixMin(DFRACT_BITS - 1, (LPC_SCALE_FACTOR + dynamicScale));
      for (i = startSample; i < stopSample; i++) {
        qmfBufferReal[i][loBand] = lowBandReal[LPC_ORDER + i] >> descale;
        qmfBufferImag[i][loBand] = lowBandImag[LPC_ORDER + i] >> descale;
      }
    } else {

      int descale = fixMin(DFRACT_BITS - 1, (LPC_SCALE_FACTOR + dynamicScale));
      dynamicScale +=
          1;

      for (i = startSample; i < stopSample; i++) {
        FIXP_DBL accu1, accu2;

        accu1 = (fMultDiv2(a0r, lowBandReal[LPC_ORDER + i - 1]) -
                 fMultDiv2(a0i, lowBandImag[LPC_ORDER + i - 1]) +
                 fMultDiv2(a1r, lowBandReal[LPC_ORDER + i - 2]) -
                 fMultDiv2(a1i, lowBandImag[LPC_ORDER + i - 2])) >>
                dynamicScale;
        accu2 = (fMultDiv2(a0i, lowBandReal[LPC_ORDER + i - 1]) +
                 fMultDiv2(a0r, lowBandImag[LPC_ORDER + i - 1]) +
                 fMultDiv2(a1i, lowBandReal[LPC_ORDER + i - 2]) +
                 fMultDiv2(a1r, lowBandImag[LPC_ORDER + i - 2])) >>
                dynamicScale;

        qmfBufferReal[i][loBand] =
            (lowBandReal[LPC_ORDER + i] >> descale) + (accu1 << (1 + 1));
        qmfBufferImag[i][loBand] =
            (lowBandImag[LPC_ORDER + i] >> descale) + (accu2 << (1 + 1));
      }
    }
  }

  for (i = 0; i < nInvfBands; i++) {
    hLppTrans->bwVectorOld[i] = bwVector[i];
  }

  sbrScaleFactor->hb_scale = comBandScale - (LPC_SCALE_FACTOR);
}

SBR_ERROR
createLppTransposer(
    HANDLE_SBR_LPP_TRANS hs,
    TRANSPOSER_SETTINGS *pSettings,
    const int highBandStartSb,
    UCHAR *v_k_master,
    const int numMaster,
    const int usb,
    const int timeSlots,
    const int nCols,
    UCHAR *noiseBandTable,
    const int noNoiseBands,
    UINT fs,
    const int chan,
    const int overlap) {
  hs->pSettings = pSettings;

  pSettings->nCols = nCols;
  pSettings->overlap = overlap;

  switch (timeSlots) {
    case 15:
    case 16:
      break;

    default:
      return SBRDEC_UNSUPPORTED_CONFIG;
  }

  if (chan == 0) {
    hs->pSettings->nCols = nCols;

    return resetLppTransposer(hs, highBandStartSb, v_k_master, numMaster,
                              noiseBandTable, noNoiseBands, usb, fs);
  }
  return SBRDEC_OK;
}

static int findClosestEntry(UCHAR goalSb, UCHAR *v_k_master, UCHAR numMaster,
                            UCHAR direction) {
  int index;

  if (goalSb <= v_k_master[0]) return v_k_master[0];

  if (goalSb >= v_k_master[numMaster]) return v_k_master[numMaster];

  if (direction) {
    index = 0;
    while (v_k_master[index] < goalSb) {
      index++;
    }
  } else {
    index = numMaster;
    while (v_k_master[index] > goalSb) {
      index--;
    }
  }

  return v_k_master[index];
}

SBR_ERROR
resetLppTransposer(
    HANDLE_SBR_LPP_TRANS hLppTrans,
    UCHAR highBandStartSb,
    UCHAR *v_k_master,
    UCHAR numMaster,
    UCHAR *noiseBandTable,
    UCHAR noNoiseBands,
    UCHAR usb,
    UINT fs
) {
  TRANSPOSER_SETTINGS *pSettings = hLppTrans->pSettings;
  PATCH_PARAM *patchParam = pSettings->patchParam;

  int i, patch;
  int targetStopBand;
  int sourceStartBand;
  int patchDistance;
  int numBandsInPatch;

  int lsb = v_k_master[0];
  int xoverOffset = highBandStartSb -
                    lsb;
  int startFreqHz;

  int desiredBorder;

  usb = fixMin(usb, v_k_master[numMaster]);

  if (pSettings->nCols == 64) {
    if (lsb < 4) {
      return SBRDEC_UNSUPPORTED_CONFIG;
    }
  } else if (lsb - SHIFT_START_SB < 4) {
    return SBRDEC_UNSUPPORTED_CONFIG;
  }

  desiredBorder = (((2048000 * 2) / fs) + 1) >> 1;

  desiredBorder = findClosestEntry(desiredBorder, v_k_master, numMaster,
                                   1);

  sourceStartBand = SHIFT_START_SB + xoverOffset;
  targetStopBand = lsb + xoverOffset;

  patch = 0;
  while (targetStopBand < usb) {
    if (patch > MAX_NUM_PATCHES) {
      return SBRDEC_UNSUPPORTED_CONFIG;
    }

    patchParam[patch].guardStartBand = targetStopBand;
    patchParam[patch].targetStartBand = targetStopBand;

    numBandsInPatch =
        desiredBorder - targetStopBand;

    if (numBandsInPatch >= lsb - sourceStartBand) {
      patchDistance =
          targetStopBand - sourceStartBand;
      patchDistance =
          patchDistance & ~1;
      numBandsInPatch =
          lsb - (targetStopBand -
                 patchDistance);
      numBandsInPatch = findClosestEntry(targetStopBand + numBandsInPatch,
                                         v_k_master, numMaster, 0) -
                        targetStopBand;
    }

    if (pSettings->nCols == 64) {
      if (numBandsInPatch == 0 && sourceStartBand == SHIFT_START_SB) {
        return SBRDEC_UNSUPPORTED_CONFIG;
      }
    }

    patchDistance =
        numBandsInPatch + targetStopBand - lsb;
    patchDistance = (patchDistance + 1) &
                    ~1;

    if (numBandsInPatch > 0) {
      patchParam[patch].sourceStartBand = targetStopBand - patchDistance;
      patchParam[patch].targetBandOffs = patchDistance;
      patchParam[patch].numBandsInPatch = numBandsInPatch;
      patchParam[patch].sourceStopBand =
          patchParam[patch].sourceStartBand + numBandsInPatch;

      targetStopBand += patchParam[patch].numBandsInPatch;
      patch++;
    }

    sourceStartBand = SHIFT_START_SB;

    if (desiredBorder - targetStopBand < 3)
    {
      desiredBorder = usb;
    }
  }

  patch--;

  if ((patch > 0) && (patchParam[patch].numBandsInPatch < 3)) {
    patch--;
    targetStopBand =
        patchParam[patch].targetStartBand + patchParam[patch].numBandsInPatch;
  }

  if (patch >= MAX_NUM_PATCHES) {
    return SBRDEC_UNSUPPORTED_CONFIG;
  }

  pSettings->noOfPatches = patch + 1;

  pSettings->lbStartPatching = targetStopBand;
  pSettings->lbStopPatching = 0;
  for (patch = 0; patch < pSettings->noOfPatches; patch++) {
    pSettings->lbStartPatching =
        fixMin(pSettings->lbStartPatching, patchParam[patch].sourceStartBand);
    pSettings->lbStopPatching =
        fixMax(pSettings->lbStopPatching, patchParam[patch].sourceStopBand);
  }

  for (i = 0; i < noNoiseBands; i++) {
    pSettings->bwBorders[i] = noiseBandTable[i + 1];
  }
  for (; i < MAX_NUM_NOISE_VALUES; i++) {
    pSettings->bwBorders[i] = 255;
  }

  startFreqHz =
      ((lsb + xoverOffset) * fs) >> 7;

  for (i = 1; i < NUM_WHFACTOR_TABLE_ENTRIES; i++) {
    if (startFreqHz < FDK_sbrDecoder_sbr_whFactorsIndex[i]) break;
  }
  i--;

  pSettings->whFactors.off = FDK_sbrDecoder_sbr_whFactorsTable[i][0];
  pSettings->whFactors.transitionLevel =
      FDK_sbrDecoder_sbr_whFactorsTable[i][1];
  pSettings->whFactors.lowLevel = FDK_sbrDecoder_sbr_whFactorsTable[i][2];
  pSettings->whFactors.midLevel = FDK_sbrDecoder_sbr_whFactorsTable[i][3];
  pSettings->whFactors.highLevel = FDK_sbrDecoder_sbr_whFactorsTable[i][4];

  return SBRDEC_OK;
}
