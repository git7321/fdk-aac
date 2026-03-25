/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Dec subband processing

*******************************************************************************/

#include "sac_stp.h"
#include "sac_calcM1andM2.h"
#include "sac_bitdec.h"
#include "FDK_matrixCalloc.h"
#include "sac_rom.h"

#define SF_FREQ_DOMAIN_HEADROOM (2 * (1))

#define BP_GF_START 6
#define BP_GF_SIZE 25
#define HP_SIZE 9
#define STP_UPDATE_ENERGY_RATE 32

#define SF_WET 5
#define SF_DRY \
  3
#define SF_DRY_NRG \
  (4 - 1)
#define SF_WET_NRG \
  (4 - 1)
#define SF_PRODUCT_BP_GF 13
#define SF_PRODUCT_BP_GF_GF 26
#define SF_SCALE 2

#define SF_SCALE_LD64 FL2FXCONST_DBL(0.03125)
#define STP_LPF_COEFF1__FDK FL2FXCONST_DBL(0.950f)
#define ONE_MINUS_STP_LPF_COEFF1__FDK FL2FXCONST_DBL(0.05f)
#define STP_LPF_COEFF2__FDK FL2FXCONST_DBL(0.450f)
#define ONE_MINUS_STP_LPF_COEFF2__FDK \
  FL2FXCONST_DBL(1.0f - 0.450f)
#define STP_SCALE_LIMIT__FDK \
  FL2FXCONST_DBL(2.82f / (float)(1 << SF_SCALE))
#define ONE_DIV_STP_SCALE_LIMIT__FDK \
  FL2FXCONST_DBL(1.0f / 2.82f / (float)(1 << SF_SCALE))
#define ABS_THR__FDK \
  FL2FXCONST_DBL(ABS_THR / \
                 ((float)(1 << (22 + 22 - 26))))
#define ABS_THR2__FDK                      \
  FL2FXCONST_DBL(ABS_THR * 32.0f * 32.0f / \
                 ((float)(1 << (22 + 22 - 26))))
#define STP_SCALE_LIMIT_HI \
  FL2FXCONST_DBL(3.02222222222 / (1 << SF_SCALE))
#define STP_SCALE_LIMIT_LO \
  FL2FXCONST_DBL(0.28289992119 / (1 << SF_SCALE))
#define STP_SCALE_LIMIT_HI_LD64 \
  FL2FXCONST_DBL(0.04986280452)
#define STP_SCALE_LIMIT_LO_LD64 \
  FL2FXCONST_DBL(0.05692613500)

#define CALC_WET_SCALE(dryIdx, wetIdx)                                         \
  if ((DryEnerLD64[dryIdx] - STP_SCALE_LIMIT_HI_LD64) > WetEnerLD64[wetIdx]) { \
    scale[wetIdx] = STP_SCALE_LIMIT_HI;                                        \
  } else if (DryEnerLD64[dryIdx] <                                             \
             (WetEnerLD64[wetIdx] - STP_SCALE_LIMIT_LO_LD64)) {                \
    scale[wetIdx] = STP_SCALE_LIMIT_LO;                                        \
  } else {                                                                     \
    tmp = ((DryEnerLD64[dryIdx] - WetEnerLD64[wetIdx]) >> 1) - SF_SCALE_LD64;  \
    scale[wetIdx] = CalcInvLdData(tmp);                                        \
  }

struct STP_DEC {
  FIXP_DBL runDryEner[MAX_INPUT_CHANNELS];
  FIXP_DBL runWetEner[MAX_OUTPUT_CHANNELS];
  FIXP_DBL oldDryEnerLD64[MAX_INPUT_CHANNELS];
  FIXP_DBL oldWetEnerLD64[MAX_OUTPUT_CHANNELS];
  FIXP_DBL prev_tp_scale[MAX_OUTPUT_CHANNELS];
  const FIXP_CFG *BP;
  const FIXP_CFG *BP_GF;
  int update_old_ener;
};

inline void combineSignalCplx(FIXP_DBL *hybOutputRealDry,
                              FIXP_DBL *hybOutputImagDry,
                              FIXP_DBL *hybOutputRealWet,
                              FIXP_DBL *hybOutputImagWet, int bands) {
  int n;

  for (n = bands - 1; n >= 0; n--) {
    *hybOutputRealDry = fAddSaturate(*hybOutputRealDry, *hybOutputRealWet);
    *hybOutputImagDry = fAddSaturate(*hybOutputImagDry, *hybOutputImagWet);
    hybOutputRealDry++, hybOutputRealWet++;
    hybOutputImagDry++, hybOutputImagWet++;
  }
}

inline void combineSignalCplxScale1(FIXP_DBL *hybOutputRealDry,
                                    FIXP_DBL *hybOutputImagDry,
                                    FIXP_DBL *hybOutputRealWet,
                                    FIXP_DBL *hybOutputImagWet,
                                    const FIXP_CFG *pBP, FIXP_DBL scaleX,
                                    int bands) {
  int n;
  FIXP_DBL scaleY;
  for (n = bands - 1; n >= 0; n--) {
    scaleY = fMult(scaleX, *pBP);
    *hybOutputRealDry = SATURATE_LEFT_SHIFT(
        (*hybOutputRealDry >> SF_SCALE) + fMult(*hybOutputRealWet, scaleY),
        SF_SCALE, DFRACT_BITS);
    *hybOutputImagDry = SATURATE_LEFT_SHIFT(
        (*hybOutputImagDry >> SF_SCALE) + fMult(*hybOutputImagWet, scaleY),
        SF_SCALE, DFRACT_BITS);
    hybOutputRealDry++, hybOutputRealWet++;
    hybOutputImagDry++, hybOutputImagWet++;
    pBP++;
  }
}

inline void combineSignalCplxScale2(FIXP_DBL *hybOutputRealDry,
                                    FIXP_DBL *hybOutputImagDry,
                                    FIXP_DBL *hybOutputRealWet,
                                    FIXP_DBL *hybOutputImagWet, FIXP_DBL scaleX,
                                    int bands) {
  int n;

  for (n = bands - 1; n >= 0; n--) {
    *hybOutputRealDry = SATURATE_LEFT_SHIFT(
        (*hybOutputRealDry >> SF_SCALE) + fMult(*hybOutputRealWet, scaleX),
        SF_SCALE, DFRACT_BITS);
    *hybOutputImagDry = SATURATE_LEFT_SHIFT(
        (*hybOutputImagDry >> SF_SCALE) + fMult(*hybOutputImagWet, scaleX),
        SF_SCALE, DFRACT_BITS);
    hybOutputRealDry++, hybOutputRealWet++;
    hybOutputImagDry++, hybOutputImagWet++;
  }
}

SACDEC_ERROR subbandTPCreate(HANDLE_STP_DEC *hStpDec) {
  HANDLE_STP_DEC self = NULL;
  FDK_ALLOCATE_MEMORY_1D(self, 1, struct STP_DEC)
  if (hStpDec != NULL) {
    *hStpDec = self;
  }

  return MPS_OK;
bail:
  return MPS_OUTOFMEMORY;
}

SACDEC_ERROR subbandTPInit(HANDLE_STP_DEC self) {
  SACDEC_ERROR err = MPS_OK;
  int ch;

  for (ch = 0; ch < MAX_OUTPUT_CHANNELS; ch++) {
    self->prev_tp_scale[ch] = FL2FXCONST_DBL(1.0f / (1 << SF_SCALE));
    self->oldWetEnerLD64[ch] = FL2FXCONST_DBL(0.0);
  }
  for (ch = 0; ch < MAX_INPUT_CHANNELS; ch++) {
    self->oldDryEnerLD64[ch] = FL2FXCONST_DBL(0.0);
  }

  self->BP = BP__FDK;
  self->BP_GF = BP_GF__FDK;

  self->update_old_ener = 0;

  return err;
}

void subbandTPDestroy(HANDLE_STP_DEC *hStpDec) {
  if (hStpDec != NULL) {
    FDK_FREE_MEMORY_1D(*hStpDec);
  }
}

SACDEC_ERROR subbandTPApply(spatialDec *self, const SPATIAL_BS_FRAME *frame) {
  FIXP_DBL *qmfOutputRealDry[MAX_OUTPUT_CHANNELS];
  FIXP_DBL *qmfOutputImagDry[MAX_OUTPUT_CHANNELS];
  FIXP_DBL *qmfOutputRealWet[MAX_OUTPUT_CHANNELS];
  FIXP_DBL *qmfOutputImagWet[MAX_OUTPUT_CHANNELS];

  FIXP_DBL DryEner[MAX_INPUT_CHANNELS];
  FIXP_DBL scale[MAX_OUTPUT_CHANNELS];

  FIXP_DBL DryEnerLD64[MAX_INPUT_CHANNELS];
  FIXP_DBL WetEnerLD64[MAX_OUTPUT_CHANNELS];

  FIXP_DBL DryEner0 = FL2FXCONST_DBL(0.0f);
  FIXP_DBL WetEnerX, damp, tmp;
  FIXP_DBL dmxReal0, dmxImag0;
  int skipChannels[MAX_OUTPUT_CHANNELS];
  int n, ch, cplxBands, cplxHybBands;
  int dry_scale_dmx, wet_scale_dmx;
  int i_LF, i_RF;
  HANDLE_STP_DEC hStpDec;
  const FIXP_CFG *pBP;

  int nrgScale = (2 * self->clipProtectGainSF__FDK);

  hStpDec = self->hStpDec;

  FDK_ASSERT(SF_DRY >= 1);
  {
    cplxBands = BP_GF_SIZE;
    cplxHybBands = self->hybridBands;
    if (self->treeConfig == TREE_212) {
      dry_scale_dmx = 2;
    } else {
      dry_scale_dmx = (2 * SF_DRY) - 2;
    }
    wet_scale_dmx = 2;
  }

  for (ch = 0; ch < self->numOutputChannels; ch++) {
    qmfOutputRealDry[ch] = &self->hybOutputRealDry__FDK[ch][7];
    qmfOutputRealWet[ch] = &self->hybOutputRealWet__FDK[ch][7];
    qmfOutputImagDry[ch] = &self->hybOutputImagDry__FDK[ch][7];
    qmfOutputImagWet[ch] = &self->hybOutputImagWet__FDK[ch][7];
  }

  FDKmemset(skipChannels, 0, self->numOutputChannels * sizeof(int));
  FDKmemset(scale, 0, self->numOutputChannels * sizeof(FIXP_DBL));

  if (hStpDec->update_old_ener == STP_UPDATE_ENERGY_RATE) {
    hStpDec->update_old_ener = 1;
    for (ch = 0; ch < self->numInputChannels; ch++) {
      hStpDec->oldDryEnerLD64[ch] =
          CalcLdData(fAddSaturate(hStpDec->runDryEner[ch], ABS_THR__FDK));
    }
    for (ch = 0; ch < self->numOutputChannels; ch++) {
      if (self->treeConfig == TREE_212)
        hStpDec->oldWetEnerLD64[ch] =
            CalcLdData(fAddSaturate(hStpDec->runWetEner[ch], ABS_THR__FDK));
      else
        hStpDec->oldWetEnerLD64[ch] =
            CalcLdData(fAddSaturate(hStpDec->runWetEner[ch], ABS_THR2__FDK));
    }
  } else {
    hStpDec->update_old_ener++;
  }

  switch (self->treeConfig) {
    case TREE_212:
      i_LF = 0;
      i_RF = 1;
      break;
    default:
      return MPS_WRONG_TREECONFIG;
  }

  pBP = hStpDec->BP_GF - BP_GF_START;
  switch (self->treeConfig) {
    case TREE_212:
      INT sMin, sNorm, sReal, sImag;

      sReal = fMin(getScalefactor(&qmfOutputRealDry[i_LF][BP_GF_START],
                                  cplxBands - BP_GF_START),
                   getScalefactor(&qmfOutputRealDry[i_RF][BP_GF_START],
                                  cplxBands - BP_GF_START));
      sImag = fMin(getScalefactor(&qmfOutputImagDry[i_LF][BP_GF_START],
                                  cplxBands - BP_GF_START),
                   getScalefactor(&qmfOutputImagDry[i_RF][BP_GF_START],
                                  cplxBands - BP_GF_START));
      sMin = fMin(sReal, sImag) - 1;

      for (n = BP_GF_START; n < cplxBands; n++) {
        dmxReal0 = scaleValue(qmfOutputRealDry[i_LF][n], sMin) +
                   scaleValue(qmfOutputRealDry[i_RF][n], sMin);
        dmxImag0 = scaleValue(qmfOutputImagDry[i_LF][n], sMin) +
                   scaleValue(qmfOutputImagDry[i_RF][n], sMin);

        DryEner0 += (fMultDiv2(fPow2Div2(dmxReal0), pBP[n]) +
                     fMultDiv2(fPow2Div2(dmxImag0), pBP[n])) >>
                    SF_DRY_NRG;
      }

      sNorm = SF_FREQ_DOMAIN_HEADROOM + SF_DRY_NRG + dry_scale_dmx -
              (2 * sMin) + nrgScale;
      DryEner0 = scaleValueSaturate(
          DryEner0, fMax(fMin(sNorm, DFRACT_BITS - 1), -(DFRACT_BITS - 1)));
      break;
    default:;
  }
  DryEner[0] = DryEner0;

  for (ch = 0; ch < self->numInputChannels; ch++) {
    if (self->treeConfig != TREE_212) DryEner[ch] = DryEner[ch] << nrgScale;
    hStpDec->runDryEner[ch] =
        fMult(STP_LPF_COEFF1__FDK, hStpDec->runDryEner[ch]) +
        fMult(ONE_MINUS_STP_LPF_COEFF1__FDK, DryEner[ch]);
    if (DryEner[ch] != FL2FXCONST_DBL(0.0f)) {
      DryEnerLD64[ch] =
          fixMax((CalcLdData(DryEner[ch]) - hStpDec->oldDryEnerLD64[ch]),
                 FL2FXCONST_DBL(-0.484375f));
    } else {
      DryEnerLD64[ch] = FL2FXCONST_DBL(-0.484375f);
    }
  }
  for (; ch < MAX_INPUT_CHANNELS; ch++) {
    DryEnerLD64[ch] = FL2FXCONST_DBL(-0.484375f);
  }

  pBP = hStpDec->BP_GF - BP_GF_START;
  for (ch = 0; ch < self->numOutputChannels; ch++) {
    if (skipChannels[ch]) {
      continue;
    }

    WetEnerX = FL2FXCONST_DBL(0.0f);

    if (self->treeConfig == TREE_212) {
      INT sMin, sNorm;

      sMin = fMin(getScalefactor(&qmfOutputRealWet[ch][BP_GF_START],
                                 cplxBands - BP_GF_START),
                  getScalefactor(&qmfOutputImagWet[ch][BP_GF_START],
                                 cplxBands - BP_GF_START));

      for (n = BP_GF_START; n < cplxBands; n++) {
        WetEnerX +=
            (fMultDiv2(fPow2Div2(scaleValue(qmfOutputRealWet[ch][n], sMin)),
                       pBP[n]) +
             fMultDiv2(fPow2Div2(scaleValue(qmfOutputImagWet[ch][n], sMin)),
                       pBP[n])) >>
            SF_WET_NRG;
      }
      sNorm = SF_FREQ_DOMAIN_HEADROOM + SF_WET_NRG + wet_scale_dmx -
              (2 * sMin) + nrgScale;
      WetEnerX = scaleValueSaturate(
          WetEnerX, fMax(fMin(sNorm, DFRACT_BITS - 1), -(DFRACT_BITS - 1)));
    } else
      FDK_ASSERT(self->treeConfig == TREE_212);

    hStpDec->runWetEner[ch] =
        fMult(STP_LPF_COEFF1__FDK, hStpDec->runWetEner[ch]) +
        fMult(ONE_MINUS_STP_LPF_COEFF1__FDK, WetEnerX);

    if (WetEnerX == FL2FXCONST_DBL(0.0f)) {
      WetEnerLD64[ch] = FL2FXCONST_DBL(-0.484375f);
    } else {
      WetEnerLD64[ch] =
          fixMax((CalcLdData(WetEnerX) - hStpDec->oldWetEnerLD64[ch]),
                 FL2FXCONST_DBL(-0.484375f));
    }
  }

  switch (self->treeConfig) {
    case TREE_212:
      if (DryEner[0] != FL2FXCONST_DBL(0.0f)) {
        CALC_WET_SCALE(0, i_LF);
        CALC_WET_SCALE(0, i_RF);
      }
      break;
    default:;
  }

  damp = FL2FXCONST_DBL(0.1f / (1 << SF_SCALE));
  for (ch = 0; ch < self->numOutputChannels; ch++) {
    scale[ch] = damp + fMult(FL2FXCONST_DBL(0.9f), scale[ch]);

    if (scale[ch] > STP_SCALE_LIMIT__FDK) {
      scale[ch] = STP_SCALE_LIMIT__FDK;
    }
    if (scale[ch] < ONE_DIV_STP_SCALE_LIMIT__FDK) {
      scale[ch] = ONE_DIV_STP_SCALE_LIMIT__FDK;
    }

    scale[ch] =
        fMult(STP_LPF_COEFF2__FDK, scale[ch]) +
        fMult(ONE_MINUS_STP_LPF_COEFF2__FDK, hStpDec->prev_tp_scale[ch]);
    hStpDec->prev_tp_scale[ch] = scale[ch];
  }

  FDK_ASSERT((HP_SIZE - 3 + 10 - 1) == PC_NUM_HYB_BANDS);
  const SCHAR *channlIndex = row2channelSTP[self->treeConfig];

  for (ch = 0; ch < self->numOutputChannels; ch++) {
    int no_scaling;

    no_scaling = !frame->tempShapeEnableChannelSTP[channlIndex[ch]];
    if (no_scaling) {
      combineSignalCplx(
          &self->hybOutputRealDry__FDK[ch][self->tp_hybBandBorder],
          &self->hybOutputImagDry__FDK[ch][self->tp_hybBandBorder],
          &self->hybOutputRealWet__FDK[ch][self->tp_hybBandBorder],
          &self->hybOutputImagWet__FDK[ch][self->tp_hybBandBorder],
          cplxHybBands - self->tp_hybBandBorder);

    } else {
      FIXP_DBL scaleX;
      scaleX = scale[ch];
      pBP = hStpDec->BP - self->tp_hybBandBorder;
      combineSignalCplxScale1(
          &self->hybOutputRealDry__FDK[ch][self->tp_hybBandBorder],
          &self->hybOutputImagDry__FDK[ch][self->tp_hybBandBorder],
          &self->hybOutputRealWet__FDK[ch][self->tp_hybBandBorder],
          &self->hybOutputImagWet__FDK[ch][self->tp_hybBandBorder],
          &pBP[self->tp_hybBandBorder], scaleX,
          (HP_SIZE - 3 + 10 - 1) - self->tp_hybBandBorder);

      {
        combineSignalCplxScale2(
            &self->hybOutputRealDry__FDK[ch][HP_SIZE - 3 + 10 - 1],
            &self->hybOutputImagDry__FDK[ch][HP_SIZE - 3 + 10 - 1],
            &self->hybOutputRealWet__FDK[ch][HP_SIZE - 3 + 10 - 1],
            &self->hybOutputImagWet__FDK[ch][HP_SIZE - 3 + 10 - 1], scaleX,
            cplxHybBands - (HP_SIZE - 3 + 10 - 1));
      }
    }
  }

  return (SACDEC_ERROR)MPS_OK;
  ;
}
