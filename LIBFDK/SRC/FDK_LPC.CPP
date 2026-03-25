/******************* Library for basic calculation routines ********************

   Author(s):   Manuel Jander

   Description: LPC related functions

*******************************************************************************/

#include "FDK_lpc.h"

static SCHAR order_ld[LPC_MAX_ORDER] = {
    0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3};

#ifndef FUNCTION_CLpc_SynthesisLattice_SGL
void CLpc_SynthesisLattice(FIXP_DBL *signal, const int signal_size,
                           const int signal_e, const int signal_e_out,
                           const int inc, const FIXP_SGL *coeff,
                           const int order, FIXP_DBL *state) {
  int i, j;
  FIXP_DBL *pSignal;
  int shift;

  FDK_ASSERT(order <= LPC_MAX_ORDER);
  FDK_ASSERT(order > 0);

  if (inc == -1)
    pSignal = &signal[signal_size - 1];
  else
    pSignal = &signal[0];

  shift = -order_ld[order - 1];

  for (i = signal_size; i != 0; i--) {
    FIXP_DBL *pState = state + order - 1;
    const FIXP_SGL *pCoeff = coeff + order - 1;
    FIXP_DBL tmp;

    tmp = scaleValue(*pSignal, shift + signal_e) -
          fMultDiv2(*pCoeff--, *pState--);
    for (j = order - 1; j != 0; j--) {
      tmp = fMultSubDiv2(tmp, pCoeff[0], pState[0]);
      pState[1] = pState[0] + (fMultDiv2(*pCoeff--, tmp) << 2);
      pState--;
    }

    *pSignal = scaleValueSaturate(tmp, -shift - signal_e_out);

    pState[1] = tmp << 1;
    pSignal += inc;
  }
}
#endif

#ifndef FUNCTION_CLpc_SynthesisLattice_DBL
void CLpc_SynthesisLattice(FIXP_DBL *signal, const int signal_size,
                           const int signal_e, const int signal_e_out,
                           const int inc, const FIXP_DBL *coeff,
                           const int order, FIXP_DBL *state) {
  int i, j;
  FIXP_DBL *pSignal;

  FDK_ASSERT(order <= LPC_MAX_ORDER);
  FDK_ASSERT(order > 0);

  if (inc == -1)
    pSignal = &signal[signal_size - 1];
  else
    pSignal = &signal[0];

  FDK_ASSERT(signal_size > 0);
  for (i = signal_size; i != 0; i--) {
    FIXP_DBL *pState = state + order - 1;
    const FIXP_DBL *pCoeff = coeff + order - 1;
    FIXP_DBL tmp, accu;

    accu =
        fMultSubDiv2(scaleValue(*pSignal, signal_e - 1), *pCoeff--, *pState--);
    tmp = SATURATE_LEFT_SHIFT_ALT(accu, 1, DFRACT_BITS);

    for (j = order - 1; j != 0; j--) {
      accu = fMultSubDiv2(tmp >> 1, pCoeff[0], pState[0]);
      tmp = SATURATE_LEFT_SHIFT_ALT(accu, 1, DFRACT_BITS);

      accu = fMultAddDiv2(pState[0] >> 1, *pCoeff--, tmp);
      pState[1] = SATURATE_LEFT_SHIFT_ALT(accu, 1, DFRACT_BITS);

      pState--;
    }

    *pSignal = scaleValue(tmp, -signal_e_out);

    pState[1] = tmp;
    pSignal += inc;
  }
}

#endif

void CLpc_Synthesis(FIXP_DBL *signal, const int signal_size, const int signal_e,
                    const int inc, const FIXP_LPC_TNS *lpcCoeff_m,
                    const int lpcCoeff_e, const int order, FIXP_DBL *state,
                    int *pStateIndex) {
  int i, j;
  FIXP_DBL *pSignal;
  int stateIndex = *pStateIndex;

  FIXP_LPC_TNS coeff[2 * LPC_MAX_ORDER];
  FDKmemcpy(&coeff[0], lpcCoeff_m, order * sizeof(FIXP_LPC_TNS));
  FDKmemcpy(&coeff[order], lpcCoeff_m, order * sizeof(FIXP_LPC_TNS));

  FDK_ASSERT(order <= LPC_MAX_ORDER);
  FDK_ASSERT(stateIndex < order);

  if (inc == -1)
    pSignal = &signal[signal_size - 1];
  else
    pSignal = &signal[0];

  for (i = 0; i < signal_size; i++) {
    FIXP_DBL x;
    const FIXP_LPC_TNS *pCoeff = coeff + order - stateIndex;

    x = scaleValue(*pSignal, -(lpcCoeff_e + 1));
    for (j = 0; j < order; j++) {
      x -= fMultDiv2(state[j], pCoeff[j]);
    }
    x = SATURATE_SHIFT(x, -lpcCoeff_e - 1, DFRACT_BITS);

    stateIndex = ((stateIndex - 1) < 0) ? (order - 1) : (stateIndex - 1);
    state[stateIndex] = x;

    *pSignal = scaleValue(x, signal_e);
    pSignal += inc;
  }

  *pStateIndex = stateIndex;
}
void CLpc_Synthesis(FIXP_DBL *signal, const int signal_size, const int signal_e,
                    const int inc, const FIXP_LPC *lpcCoeff_m,
                    const int lpcCoeff_e, const int order, FIXP_DBL *state,
                    int *pStateIndex) {
  int i, j;
  FIXP_DBL *pSignal;
  int stateIndex = *pStateIndex;

  FIXP_LPC coeff[2 * LPC_MAX_ORDER];
  FDKmemcpy(&coeff[0], lpcCoeff_m, order * sizeof(FIXP_LPC));
  FDKmemcpy(&coeff[order], lpcCoeff_m, order * sizeof(FIXP_LPC));

  FDK_ASSERT(order <= LPC_MAX_ORDER);
  FDK_ASSERT(stateIndex < order);

  if (inc == -1)
    pSignal = &signal[signal_size - 1];
  else
    pSignal = &signal[0];

  for (i = 0; i < signal_size; i++) {
    FIXP_DBL x;
    const FIXP_LPC *pCoeff = coeff + order - stateIndex;

    x = scaleValue(*pSignal, -(lpcCoeff_e + 1));
    for (j = 0; j < order; j++) {
      x -= fMultDiv2(state[j], pCoeff[j]);
    }
    x = SATURATE_SHIFT(x, -lpcCoeff_e - 1, DFRACT_BITS);

    stateIndex = ((stateIndex - 1) < 0) ? (order - 1) : (stateIndex - 1);
    state[stateIndex] = x;

    *pSignal = scaleValue(x, signal_e);
    pSignal += inc;
  }

  *pStateIndex = stateIndex;
}

void CLpc_Analysis(FIXP_DBL *RESTRICT signal, const int signal_size,
                   const FIXP_LPC lpcCoeff_m[], const int lpcCoeff_e,
                   const int order, FIXP_DBL *RESTRICT filtState,
                   int *filtStateIndex) {
  int stateIndex;
  INT i, j, shift = lpcCoeff_e + 1;
  FIXP_DBL tmp;

  if (order <= 0) {
    return;
  }
  if (filtStateIndex != NULL) {
    stateIndex = *filtStateIndex;
  } else {
    stateIndex = 0;
  }

  FIXP_LPC coeff[2 * LPC_MAX_ORDER];
  FIXP_LPC *pCoeff;
  FDKmemcpy(&coeff[0], lpcCoeff_m, order * sizeof(FIXP_LPC));
  FDKmemcpy(&coeff[order], lpcCoeff_m, order * sizeof(FIXP_LPC));

  FDK_ASSERT(shift >= 0);

  for (j = 0; j < signal_size; j++) {
    pCoeff = &coeff[(order - stateIndex)];

    tmp = signal[j] >> shift;
    for (i = 0; i < order; i++) {
      tmp = fMultAddDiv2(tmp, pCoeff[i], filtState[i]);
    }

    stateIndex =
        ((stateIndex - 1) < 0) ? (stateIndex - 1 + order) : (stateIndex - 1);
    filtState[stateIndex] = signal[j];

    signal[j] = tmp << shift;
  }

  if (filtStateIndex != NULL) {
    *filtStateIndex = stateIndex;
  }
}

INT CLpc_ParcorToLpc(const FIXP_LPC_TNS reflCoeff[], FIXP_LPC_TNS LpcCoeff[],
                     INT numOfCoeff, FIXP_DBL workBuffer[]) {
  INT i, j;
  INT shiftval,
      par2LpcShiftVal = 6;
  FIXP_DBL maxVal = (FIXP_DBL)0;

  workBuffer[0] = FX_LPC_TNS2FX_DBL(reflCoeff[0]) >> par2LpcShiftVal;
  for (i = 1; i < numOfCoeff; i++) {
    for (j = 0; j < i / 2; j++) {
      FIXP_DBL tmp1, tmp2;

      tmp1 = workBuffer[j];
      tmp2 = workBuffer[i - 1 - j];
      workBuffer[j] += fMult(reflCoeff[i], tmp2);
      workBuffer[i - 1 - j] += fMult(reflCoeff[i], tmp1);
    }
    if (i & 1) {
      workBuffer[j] += fMult(reflCoeff[i], workBuffer[j]);
    }

    workBuffer[i] = FX_LPC_TNS2FX_DBL(reflCoeff[i]) >> par2LpcShiftVal;
  }

  for (i = 0; i < numOfCoeff; i++) {
    maxVal = fMax(maxVal, fAbs(workBuffer[i]));
  }

  shiftval = fMin(fNorm(maxVal), par2LpcShiftVal);

  for (i = 0; i < numOfCoeff; i++) {
    LpcCoeff[i] = FX_DBL2FX_LPC_TNS(workBuffer[i] << shiftval);
  }

  return (par2LpcShiftVal - shiftval);
}
INT CLpc_ParcorToLpc(const FIXP_LPC reflCoeff[], FIXP_LPC LpcCoeff[],
                     INT numOfCoeff, FIXP_DBL workBuffer[]) {
  INT i, j;
  INT shiftval,
      par2LpcShiftVal = 6;
  FIXP_DBL maxVal = (FIXP_DBL)0;

  workBuffer[0] = FX_LPC2FX_DBL(reflCoeff[0]) >> par2LpcShiftVal;
  for (i = 1; i < numOfCoeff; i++) {
    for (j = 0; j < i / 2; j++) {
      FIXP_DBL tmp1, tmp2;

      tmp1 = workBuffer[j];
      tmp2 = workBuffer[i - 1 - j];
      workBuffer[j] += fMult(reflCoeff[i], tmp2);
      workBuffer[i - 1 - j] += fMult(reflCoeff[i], tmp1);
    }
    if (i & 1) {
      workBuffer[j] += fMult(reflCoeff[i], workBuffer[j]);
    }

    workBuffer[i] = FX_LPC2FX_DBL(reflCoeff[i]) >> par2LpcShiftVal;
  }

  for (i = 0; i < numOfCoeff; i++) {
    maxVal = fMax(maxVal, fAbs(workBuffer[i]));
  }

  shiftval = fMin(fNorm(maxVal), par2LpcShiftVal);

  for (i = 0; i < numOfCoeff; i++) {
    LpcCoeff[i] = FX_DBL2FX_LPC(workBuffer[i] << shiftval);
  }

  return (par2LpcShiftVal - shiftval);
}

void CLpc_AutoToParcor(FIXP_DBL acorr[], const int acorr_e,
                       FIXP_LPC reflCoeff[], const int numOfCoeff,
                       FIXP_DBL *pPredictionGain_m, INT *pPredictionGain_e) {
  INT i, j, scale = 0;
  FIXP_DBL parcorWorkBuffer[LPC_MAX_ORDER];

  FIXP_DBL *workBuffer = parcorWorkBuffer;
  FIXP_DBL autoCorr_0 = acorr[0];

  FDKmemclear(reflCoeff, numOfCoeff * sizeof(FIXP_LPC));

  if (autoCorr_0 == FL2FXCONST_DBL(0.0)) {
    if (pPredictionGain_m != NULL) {
      *pPredictionGain_m = FL2FXCONST_DBL(0.5f);
      *pPredictionGain_e = 1;
    }
    return;
  }

  FDKmemcpy(workBuffer, acorr + 1, numOfCoeff * sizeof(FIXP_DBL));
  for (i = 0; i < numOfCoeff; i++) {
    LONG sign = ((LONG)workBuffer[0] >> (DFRACT_BITS - 1));
    FIXP_DBL tmp = (FIXP_DBL)((LONG)workBuffer[0] ^ sign);

    if (acorr[0] < tmp) break;

    tmp = (FIXP_DBL)((LONG)schur_div(tmp, acorr[0], FRACT_BITS) ^ (~sign));

    reflCoeff[i] = FX_DBL2FX_LPC(tmp);

    for (j = numOfCoeff - i - 1; j >= 0; j--) {
      FIXP_DBL accu1 = fMult(tmp, acorr[j]);
      FIXP_DBL accu2 = fMult(tmp, workBuffer[j]);
      workBuffer[j] += accu1;
      acorr[j] += accu2;
    }
    if (acorr[0] == (FIXP_DBL)0) break;

    workBuffer++;
  }

  if (pPredictionGain_m != NULL) {
    if (acorr[0] > (FIXP_DBL)0) {
      *pPredictionGain_m = fDivNormSigned(autoCorr_0, acorr[0], &scale);
      *pPredictionGain_e = scale;
    } else {
      *pPredictionGain_m = (FIXP_DBL)0;
      *pPredictionGain_e = 0;
    }
  }
}
