/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Dec parameter smoothing

*******************************************************************************/

#include "sac_dec.h"
#include "sac_bitdec.h"
#include "sac_smoothing.h"
#include "sac_rom.h"

static FIXP_DBL calcFilterCoeff__FDK(spatialDec *self, int ps,
                                     const SPATIAL_BS_FRAME *frame) {
  int dSlots;
  FIXP_DBL delta;

  dSlots = frame->paramSlot[ps] - self->smoothState->prevParamSlot;

  if (dSlots <= 0) {
    dSlots += self->timeSlots;
  }

  delta = fDivNorm(dSlots, self->smgTime[ps]);

  return delta;
}

static int getSmoothOnOff(spatialDec *self, int ps, int pb) {
  int smoothBand = 0;

  smoothBand = self->smgData[ps][pb];

  return smoothBand;
}

void SpatialDecSmoothM1andM2(spatialDec *self, const SPATIAL_BS_FRAME *frame,
                             int ps) {
  FIXP_DBL delta__FDK;
  FIXP_DBL one_minus_delta__FDK;

  int pb, row, col;
  int residualBands = 0;

  if (self->residualCoding) {
    int i;
    int boxes = self->numOttBoxes;
    for (i = 0; i < boxes; i++) {
      if (self->residualBands[i] > residualBands) {
        residualBands = self->residualBands[i];
      }
    }
  }

  delta__FDK = calcFilterCoeff__FDK(self, ps, frame);
  if (delta__FDK == (FIXP_DBL)MAXVAL_DBL)
    one_minus_delta__FDK = FL2FXCONST_DBL(0.0f);
  else if (delta__FDK == FL2FXCONST_DBL(0.0f))
    one_minus_delta__FDK = (FIXP_DBL)MAXVAL_DBL;
  else
    one_minus_delta__FDK = (FL2FXCONST_DBL(0.5f) - (delta__FDK >> 1)) << 1;

  for (pb = 0; pb < self->numParameterBands; pb++) {
    int smoothBand;

    smoothBand = getSmoothOnOff(self, ps, pb);

    if (smoothBand && (pb >= residualBands)) {
      for (row = 0; row < self->numM2rows; row++) {
        for (col = 0; col < self->numVChannels; col++) {
          self->M2Real__FDK[row][col][pb] =
              ((fMultDiv2(delta__FDK, self->M2Real__FDK[row][col][pb]) +
                fMultDiv2(one_minus_delta__FDK,
                          self->M2RealPrev__FDK[row][col][pb]))
               << 1);
          if (0 || (self->phaseCoding == 3)) {
            self->M2Imag__FDK[row][col][pb] =
                ((fMultDiv2(delta__FDK, self->M2Imag__FDK[row][col][pb]) +
                  fMultDiv2(one_minus_delta__FDK,
                            self->M2ImagPrev__FDK[row][col][pb]))
                 << 1);
          }
        }
      }
    }
  }
  self->smoothState->prevParamSlot = frame->paramSlot[ps];
}

void initParameterSmoothing(spatialDec *self) {
  self->smoothState->prevParamSlot = 0;
}

void SpatialDecSmoothOPD(spatialDec *self, const SPATIAL_BS_FRAME *frame,
                         int ps) {
  int pb;
  int dSlots;
  FIXP_DBL delta__FDK;
  FIXP_DBL one_minus_delta__FDK;
  FIXP_DBL *phaseLeftSmooth__FDK = self->smoothState->opdLeftState__FDK;
  FIXP_DBL *phaseRightSmooth__FDK = self->smoothState->opdRightState__FDK;
  int quantCoarse;

  quantCoarse = frame->IPDLosslessData[0].bsQuantCoarseXXX[ps];

  if (frame->OpdSmoothingMode == 0) {
    FDKmemcpy(phaseLeftSmooth__FDK, self->PhaseLeft__FDK,
              self->numParameterBands * sizeof(FIXP_DBL));
    FDKmemcpy(phaseRightSmooth__FDK, self->PhaseRight__FDK,
              self->numParameterBands * sizeof(FIXP_DBL));
  } else {
    if (ps == 0) {
      dSlots = frame->paramSlot[ps] + 1;
    } else {
      dSlots = frame->paramSlot[ps] - frame->paramSlot[ps - 1];
    }

    delta__FDK = (FIXP_DBL)((INT)(FL2FXCONST_DBL(0.0078125f)) * dSlots);

    if (delta__FDK == (FIXP_DBL)MAXVAL_DBL)
      one_minus_delta__FDK = FL2FXCONST_DBL(0.0f);
    else if (delta__FDK == FL2FXCONST_DBL(0.0f))
      one_minus_delta__FDK = (FIXP_DBL)MAXVAL_DBL;
    else
      one_minus_delta__FDK = (FL2FXCONST_DBL(0.5f) - (delta__FDK >> 1)) << 1;

    for (pb = 0; pb < self->numParameterBands; pb++) {
      FIXP_DBL tmpL, tmpR, tmp;

      tmpL = self->PhaseLeft__FDK[pb];
      tmpR = self->PhaseRight__FDK[pb];

      while (tmpL > phaseLeftSmooth__FDK[pb] + PI__IPD) tmpL -= PI__IPD << 1;
      while (tmpL < phaseLeftSmooth__FDK[pb] - PI__IPD) tmpL += PI__IPD << 1;
      while (tmpR > phaseRightSmooth__FDK[pb] + PI__IPD) tmpR -= PI__IPD << 1;
      while (tmpR < phaseRightSmooth__FDK[pb] - PI__IPD) tmpR += PI__IPD << 1;

      phaseLeftSmooth__FDK[pb] =
          fMult(delta__FDK, tmpL) +
          fMult(one_minus_delta__FDK, phaseLeftSmooth__FDK[pb]);
      phaseRightSmooth__FDK[pb] =
          fMult(delta__FDK, tmpR) +
          fMult(one_minus_delta__FDK, phaseRightSmooth__FDK[pb]);

      tmp = (((tmpL >> 1) - (tmpR >> 1)) - ((phaseLeftSmooth__FDK[pb] >> 1) -
                                            (phaseRightSmooth__FDK[pb] >> 1)))
            << 1;
      while (tmp > PI__IPD) tmp -= PI__IPD << 1;
      while (tmp < -PI__IPD) tmp += PI__IPD << 1;
      if (fixp_abs(tmp) > fMult((quantCoarse ? FL2FXCONST_DBL(50.f / 180.f)
                                             : FL2FXCONST_DBL(25.f / 180.f)),
                                PI__IPD)) {
        phaseLeftSmooth__FDK[pb] = tmpL;
        phaseRightSmooth__FDK[pb] = tmpR;
      }

      while (phaseLeftSmooth__FDK[pb] > PI__IPD << 1)
        phaseLeftSmooth__FDK[pb] -= PI__IPD << 1;
      while (phaseLeftSmooth__FDK[pb] < (FIXP_DBL)0)
        phaseLeftSmooth__FDK[pb] += PI__IPD << 1;
      while (phaseRightSmooth__FDK[pb] > PI__IPD << 1)
        phaseRightSmooth__FDK[pb] -= PI__IPD << 1;
      while (phaseRightSmooth__FDK[pb] < (FIXP_DBL)0)
        phaseRightSmooth__FDK[pb] += PI__IPD << 1;

      self->PhaseLeft__FDK[pb] = phaseLeftSmooth__FDK[pb];
      self->PhaseRight__FDK[pb] = phaseRightSmooth__FDK[pb];
    }
  }
  return;
}
