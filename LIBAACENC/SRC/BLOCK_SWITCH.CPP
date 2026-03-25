/**************************** AAC encoder library ******************************

   Author(s):   M. Werner, Tobias Chalupka

   Description: Block switching

*******************************************************************************/

#include "block_switch.h"
#include "genericStds.h"

#define LOWOV_WINDOW _LOWOV_WINDOW

static FIXP_DBL FDKaacEnc_GetWindowEnergy(const FIXP_DBL in[],
                                          const INT blSwWndIdx);

static void FDKaacEnc_CalcWindowEnergy(
    BLOCK_SWITCHING_CONTROL *RESTRICT blockSwitchingControl, INT windowLen,
    const INT_PCM *pTimeSignal);

static const INT blockType2windowShape[2][5] = {
    {SINE_WINDOW, KBD_WINDOW, WRONG_WINDOW, SINE_WINDOW, KBD_WINDOW},
    {KBD_WINDOW, SINE_WINDOW, SINE_WINDOW, KBD_WINDOW, WRONG_WINDOW}};

#ifndef SINETABLE_16BIT

static const FIXP_DBL hiPassCoeff[BLOCK_SWITCHING_IIR_LEN] = {
    FL2FXCONST_DBL(-0.5095), FL2FXCONST_DBL(0.7548)};

static const FIXP_DBL accWindowNrgFac =
    FL2FXCONST_DBL(0.3f);
static const FIXP_DBL oneMinusAccWindowNrgFac = FL2FXCONST_DBL(0.7f);
static const FIXP_DBL invAttackRatio =
    FL2FXCONST_DBL(0.1f);
static const FIXP_DBL minAttackNrg =
    (FL2FXCONST_DBL(1e+6f * NORM_PCM_ENERGY) >>
     BLOCK_SWITCH_ENERGY_SHIFT);

#else

static const FIXP_SGL hiPassCoeff[BLOCK_SWITCHING_IIR_LEN] = {
    FL2FXCONST_SGL(-0.5095), FL2FXCONST_SGL(0.7548)};
static const FIXP_DBL accWindowNrgFac =
    FL2FXCONST_DBL(0.3f);
static const FIXP_SGL oneMinusAccWindowNrgFac = FL2FXCONST_SGL(0.7f);
static const FIXP_SGL invAttackRatio =
    FL2FXCONST_SGL(0.1f);
static const FIXP_DBL minAttackNrg =
    (FL2FXCONST_DBL(1e+6f * NORM_PCM_ENERGY) >>
     BLOCK_SWITCH_ENERGY_SHIFT);

#endif

void FDKaacEnc_InitBlockSwitching(
    BLOCK_SWITCHING_CONTROL *blockSwitchingControl, INT isLowDelay) {
  FDKmemclear(blockSwitchingControl, sizeof(BLOCK_SWITCHING_CONTROL));

  if (isLowDelay) {
    blockSwitchingControl->nBlockSwitchWindows = 4;
    blockSwitchingControl->allowShortFrames = 0;
    blockSwitchingControl->allowLookAhead = 0;
  } else {
    blockSwitchingControl->nBlockSwitchWindows = 8;
    blockSwitchingControl->allowShortFrames = 1;
    blockSwitchingControl->allowLookAhead = 1;
  }

  blockSwitchingControl->noOfGroups = MAX_NO_OF_GROUPS;
  blockSwitchingControl->lastWindowSequence = LONG_WINDOW;
  blockSwitchingControl->windowShape =
      blockType2windowShape[blockSwitchingControl->allowShortFrames]
                           [blockSwitchingControl->lastWindowSequence];
}

static const INT suggestedGroupingTable[TRANS_FAC][MAX_NO_OF_GROUPS] = {
    {1, 3, 3, 1},
    {1, 1, 3, 3},
    {2, 1, 3, 2},
    {3, 1, 3, 1},
    {3, 1, 1, 3},
    {3, 2, 1, 2},
    {3, 3, 1, 1},
    {3, 3, 1, 1}};

static const INT chgWndSq[2][N_BLOCKTYPES] = {
    {LONG_WINDOW, STOP_WINDOW, WRONG_WINDOW, LONG_WINDOW,
                   STOP_WINDOW, WRONG_WINDOW},
    {START_WINDOW, LOWOV_WINDOW, WRONG_WINDOW, START_WINDOW,
                   LOWOV_WINDOW, WRONG_WINDOW}};

static const INT chgWndSqLkAhd[2][2][N_BLOCKTYPES] = {
    {
        {LONG_WINDOW, SHORT_WINDOW, STOP_WINDOW, LONG_WINDOW, WRONG_WINDOW,
         WRONG_WINDOW},
        {START_WINDOW, SHORT_WINDOW, SHORT_WINDOW, START_WINDOW,
                       WRONG_WINDOW, WRONG_WINDOW}},
    {{LONG_WINDOW, SHORT_WINDOW, SHORT_WINDOW, LONG_WINDOW,
                    WRONG_WINDOW, WRONG_WINDOW},
                   {START_WINDOW, SHORT_WINDOW, SHORT_WINDOW,
                                  START_WINDOW, WRONG_WINDOW,
                                  WRONG_WINDOW}}
};

int FDKaacEnc_BlockSwitching(BLOCK_SWITCHING_CONTROL *blockSwitchingControl,
                             const INT granuleLength, const int isLFE,
                             const INT_PCM *pTimeSignal) {
  UINT i;
  FIXP_DBL enM1, enMax;

  UINT nBlockSwitchWindows = blockSwitchingControl->nBlockSwitchWindows;

  if (isLFE) {
    blockSwitchingControl->lastWindowSequence = LONG_WINDOW;
    blockSwitchingControl->windowShape = SINE_WINDOW;
    blockSwitchingControl->noOfGroups = 1;
    blockSwitchingControl->groupLen[0] = 1;

    return (0);
  };

  blockSwitchingControl->lastattack = blockSwitchingControl->attack;
  blockSwitchingControl->lastAttackIndex = blockSwitchingControl->attackIndex;

  FDKmemcpy(blockSwitchingControl->windowNrg[0],
            blockSwitchingControl->windowNrg[1],
            sizeof(blockSwitchingControl->windowNrg[0]));
  FDKmemcpy(blockSwitchingControl->windowNrgF[0],
            blockSwitchingControl->windowNrgF[1],
            sizeof(blockSwitchingControl->windowNrgF[0]));

  if (blockSwitchingControl->allowShortFrames) {

    FDKmemclear(blockSwitchingControl->groupLen,
                sizeof(blockSwitchingControl->groupLen));

    blockSwitchingControl->noOfGroups = MAX_NO_OF_GROUPS;

    FDKmemcpy(blockSwitchingControl->groupLen,
              suggestedGroupingTable[blockSwitchingControl->lastAttackIndex],
              sizeof(blockSwitchingControl->groupLen));

    if (blockSwitchingControl->attack == TRUE)
      blockSwitchingControl->maxWindowNrg =
          FDKaacEnc_GetWindowEnergy(blockSwitchingControl->windowNrg[0],
                                    blockSwitchingControl->lastAttackIndex);
    else
      blockSwitchingControl->maxWindowNrg = FL2FXCONST_DBL(0.0);
  }

  FDKaacEnc_CalcWindowEnergy(
      blockSwitchingControl,
      granuleLength >> (nBlockSwitchWindows == 4 ? 2 : 3), pTimeSignal);

  blockSwitchingControl->attack = FALSE;

  enMax = FL2FXCONST_DBL(0.0f);
  enM1 = blockSwitchingControl->windowNrgF[0][nBlockSwitchWindows - 1];

  for (i = 0; i < nBlockSwitchWindows; i++) {
    FIXP_DBL tmp =
        fMultDiv2(oneMinusAccWindowNrgFac, blockSwitchingControl->accWindowNrg);
    blockSwitchingControl->accWindowNrg = fMultAdd(tmp, accWindowNrgFac, enM1);

    if (fMult(blockSwitchingControl->windowNrgF[1][i], invAttackRatio) >
        blockSwitchingControl->accWindowNrg) {
      blockSwitchingControl->attack = TRUE;
      blockSwitchingControl->attackIndex = i;
    }
    enM1 = blockSwitchingControl->windowNrgF[1][i];
    enMax = fixMax(enMax, enM1);
  }

  if (enMax < minAttackNrg) blockSwitchingControl->attack = FALSE;

  if ((blockSwitchingControl->attack == FALSE) &&
      (blockSwitchingControl->lastattack == TRUE)) {
    if (((blockSwitchingControl->windowNrgF[0][nBlockSwitchWindows - 1] >> 4) >
         fMult((FIXP_DBL)(10 << (DFRACT_BITS - 1 - 4)),
               blockSwitchingControl->windowNrgF[1][1])) &&
        (blockSwitchingControl->lastAttackIndex ==
         (INT)nBlockSwitchWindows - 1)) {
      blockSwitchingControl->attack = TRUE;
      blockSwitchingControl->attackIndex = 0;
    }
  }

  if (blockSwitchingControl->allowLookAhead) {
    blockSwitchingControl->lastWindowSequence =
        chgWndSqLkAhd[blockSwitchingControl->lastattack]
                     [blockSwitchingControl->attack]
                     [blockSwitchingControl->lastWindowSequence];
  } else {
    blockSwitchingControl->lastWindowSequence =
        chgWndSq[blockSwitchingControl->attack]
                [blockSwitchingControl->lastWindowSequence];
  }

  blockSwitchingControl->windowShape =
      blockType2windowShape[blockSwitchingControl->allowShortFrames]
                           [blockSwitchingControl->lastWindowSequence];

  return (0);
}

static FIXP_DBL FDKaacEnc_GetWindowEnergy(const FIXP_DBL in[],
                                          const INT blSwWndIdx) {
  return in[blSwWndIdx];
}

static void FDKaacEnc_CalcWindowEnergy(
    BLOCK_SWITCHING_CONTROL *RESTRICT blockSwitchingControl, INT windowLen,
    const INT_PCM *pTimeSignal) {
  INT i;
  UINT w;

#ifndef SINETABLE_16BIT
  const FIXP_DBL hiPassCoeff0 = hiPassCoeff[0];
  const FIXP_DBL hiPassCoeff1 = hiPassCoeff[1];
#else
  const FIXP_SGL hiPassCoeff0 = hiPassCoeff[0];
  const FIXP_SGL hiPassCoeff1 = hiPassCoeff[1];
#endif

  FIXP_DBL temp_iirState0 = blockSwitchingControl->iirStates[0];
  FIXP_DBL temp_iirState1 = blockSwitchingControl->iirStates[1];

  for (w = 0; w < blockSwitchingControl->nBlockSwitchWindows; w++) {
    ULONG temp_windowNrg = 0x0;
    ULONG temp_windowNrgF = 0x0;

    for (i = 0; i < windowLen; i++) {
      FIXP_DBL tempUnfiltered, t1, t2;
#if SAMPLE_BITS == DFRACT_BITS
      tempUnfiltered = (FIXP_DBL)*pTimeSignal++ >> 1;
#else
      tempUnfiltered = (FIXP_DBL)*pTimeSignal++
                       << (DFRACT_BITS - SAMPLE_BITS - 1);
#endif
      t1 = fMultDiv2(hiPassCoeff1, tempUnfiltered - temp_iirState0);
      t2 = fMultDiv2(hiPassCoeff0, temp_iirState1);
      temp_iirState0 = tempUnfiltered;
      temp_iirState1 = (t1 - t2) << 1;
      temp_windowNrg += (LONG)fPow2Div2(temp_iirState0) >>
                        (BLOCK_SWITCH_ENERGY_SHIFT - 1 - 2);
      temp_windowNrgF += (LONG)fPow2Div2(temp_iirState1) >>
                         (BLOCK_SWITCH_ENERGY_SHIFT - 1 - 2);
    }
    blockSwitchingControl->windowNrg[1][w] =
        (LONG)fMin(temp_windowNrg, (UINT)MAXVAL_DBL);
    blockSwitchingControl->windowNrgF[1][w] =
        (LONG)fMin(temp_windowNrgF, (UINT)MAXVAL_DBL);
  }
  blockSwitchingControl->iirStates[0] = temp_iirState0;
  blockSwitchingControl->iirStates[1] = temp_iirState1;
}

static const UCHAR synchronizedBlockTypeTable[5][5] = {
    {LONG_WINDOW, START_WINDOW, SHORT_WINDOW, STOP_WINDOW,
                        LOWOV_WINDOW},
    {START_WINDOW, START_WINDOW, SHORT_WINDOW, SHORT_WINDOW, LOWOV_WINDOW},
    {SHORT_WINDOW, SHORT_WINDOW, SHORT_WINDOW, SHORT_WINDOW, WRONG_WINDOW},
    {STOP_WINDOW, SHORT_WINDOW, SHORT_WINDOW, STOP_WINDOW, LOWOV_WINDOW},
    {LOWOV_WINDOW, LOWOV_WINDOW, WRONG_WINDOW, LOWOV_WINDOW, LOWOV_WINDOW},
};

int FDKaacEnc_SyncBlockSwitching(
    BLOCK_SWITCHING_CONTROL *blockSwitchingControlLeft,
    BLOCK_SWITCHING_CONTROL *blockSwitchingControlRight, const INT nChannels,
    const INT commonWindow) {
  UCHAR patchType = LONG_WINDOW;

  if (nChannels == 2 && commonWindow == TRUE) {
    patchType = synchronizedBlockTypeTable[patchType][blockSwitchingControlLeft
                                                          ->lastWindowSequence];
    patchType = synchronizedBlockTypeTable[patchType][blockSwitchingControlRight
                                                          ->lastWindowSequence];

    if (patchType == WRONG_WINDOW) return -1;

    blockSwitchingControlLeft->lastWindowSequence = patchType;
    blockSwitchingControlRight->lastWindowSequence = patchType;

    blockSwitchingControlLeft->windowShape =
        blockType2windowShape[blockSwitchingControlLeft->allowShortFrames]
                             [blockSwitchingControlLeft->lastWindowSequence];
    blockSwitchingControlRight->windowShape =
        blockType2windowShape[blockSwitchingControlLeft->allowShortFrames]
                             [blockSwitchingControlRight->lastWindowSequence];
  }

  if (blockSwitchingControlLeft->allowShortFrames) {
    int i;

    if (nChannels == 2) {
      if (commonWindow == TRUE) {
        int windowSequenceLeftOld =
            blockSwitchingControlLeft->lastWindowSequence;
        int windowSequenceRightOld =
            blockSwitchingControlRight->lastWindowSequence;

        if (patchType != SHORT_WINDOW) {
          blockSwitchingControlLeft->noOfGroups = 1;
          blockSwitchingControlRight->noOfGroups = 1;
          blockSwitchingControlLeft->groupLen[0] = 1;
          blockSwitchingControlRight->groupLen[0] = 1;

          for (i = 1; i < MAX_NO_OF_GROUPS; i++) {
            blockSwitchingControlLeft->groupLen[i] = 0;
            blockSwitchingControlRight->groupLen[i] = 0;
          }
        }

        else {
          if ((windowSequenceLeftOld == SHORT_WINDOW) &&
              (windowSequenceRightOld == SHORT_WINDOW)) {
            if (blockSwitchingControlLeft->maxWindowNrg >
                blockSwitchingControlRight->maxWindowNrg) {
              blockSwitchingControlRight->noOfGroups =
                  blockSwitchingControlLeft->noOfGroups;
              for (i = 0; i < MAX_NO_OF_GROUPS; i++) {
                blockSwitchingControlRight->groupLen[i] =
                    blockSwitchingControlLeft->groupLen[i];
              }
            } else {
              blockSwitchingControlLeft->noOfGroups =
                  blockSwitchingControlRight->noOfGroups;
              for (i = 0; i < MAX_NO_OF_GROUPS; i++) {
                blockSwitchingControlLeft->groupLen[i] =
                    blockSwitchingControlRight->groupLen[i];
              }
            }
          } else if ((windowSequenceLeftOld == SHORT_WINDOW) &&
                     (windowSequenceRightOld != SHORT_WINDOW)) {
            blockSwitchingControlRight->noOfGroups =
                blockSwitchingControlLeft->noOfGroups;
            for (i = 0; i < MAX_NO_OF_GROUPS; i++) {
              blockSwitchingControlRight->groupLen[i] =
                  blockSwitchingControlLeft->groupLen[i];
            }
          } else if ((windowSequenceRightOld == SHORT_WINDOW) &&
                     (windowSequenceLeftOld != SHORT_WINDOW)) {
            blockSwitchingControlLeft->noOfGroups =
                blockSwitchingControlRight->noOfGroups;
            for (i = 0; i < MAX_NO_OF_GROUPS; i++) {
              blockSwitchingControlLeft->groupLen[i] =
                  blockSwitchingControlRight->groupLen[i];
            }
          } else {
            blockSwitchingControlLeft->noOfGroups =
                blockSwitchingControlRight->noOfGroups = 2;
            blockSwitchingControlLeft->groupLen[0] =
                blockSwitchingControlRight->groupLen[0] = 4;
            blockSwitchingControlLeft->groupLen[1] =
                blockSwitchingControlRight->groupLen[1] = 4;
          }
        }
      } else {
        if (blockSwitchingControlLeft->lastWindowSequence != SHORT_WINDOW) {
          blockSwitchingControlLeft->noOfGroups = 1;
          blockSwitchingControlLeft->groupLen[0] = 1;
          for (i = 1; i < MAX_NO_OF_GROUPS; i++) {
            blockSwitchingControlLeft->groupLen[i] = 0;
          }
        }
        if (blockSwitchingControlRight->lastWindowSequence != SHORT_WINDOW) {
          blockSwitchingControlRight->noOfGroups = 1;
          blockSwitchingControlRight->groupLen[0] = 1;
          for (i = 1; i < MAX_NO_OF_GROUPS; i++) {
            blockSwitchingControlRight->groupLen[i] = 0;
          }
        }
      }
    } else {
      if (blockSwitchingControlLeft->lastWindowSequence != SHORT_WINDOW) {
        blockSwitchingControlLeft->noOfGroups = 1;
        blockSwitchingControlLeft->groupLen[0] = 1;

        for (i = 1; i < MAX_NO_OF_GROUPS; i++) {
          blockSwitchingControlLeft->groupLen[i] = 0;
        }
      }
    }
  }

  if (!blockSwitchingControlLeft->allowShortFrames) {
    if (blockSwitchingControlLeft->lastWindowSequence != LONG_WINDOW &&
        blockSwitchingControlLeft->lastWindowSequence != STOP_WINDOW) {
      blockSwitchingControlLeft->lastWindowSequence = LONG_WINDOW;
      blockSwitchingControlLeft->windowShape = LOL_WINDOW;
    }
  }
  if (nChannels == 2) {
    if (!blockSwitchingControlRight->allowShortFrames) {
      if (blockSwitchingControlRight->lastWindowSequence != LONG_WINDOW &&
          blockSwitchingControlRight->lastWindowSequence != STOP_WINDOW) {
        blockSwitchingControlRight->lastWindowSequence = LONG_WINDOW;
        blockSwitchingControlRight->windowShape = LOL_WINDOW;
      }
    }
  }

  return 0;
}
