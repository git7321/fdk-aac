/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Block switching

*******************************************************************************/

#ifndef BLOCK_SWITCH_H
#define BLOCK_SWITCH_H

#include "common_fix.h"

#include "psy_const.h"

#define BLOCK_SWITCH_WINDOWS 8

#define BLOCK_SWITCHING_IIR_LEN \
  2
#define BLOCK_SWITCH_ENERGY_SHIFT \
  7

#define LAST_WINDOW 0
#define THIS_WINDOW 1

typedef struct {
  INT lastWindowSequence;
  INT windowShape;
  INT lastWindowShape;
  UINT nBlockSwitchWindows;
  INT attack;
  INT lastattack;
  INT attackIndex;
  INT lastAttackIndex;
  INT allowShortFrames;
  INT allowLookAhead;
  INT noOfGroups;
  INT groupLen[MAX_NO_OF_GROUPS];
  FIXP_DBL maxWindowNrg;

  FIXP_DBL
  windowNrg[2][BLOCK_SWITCH_WINDOWS];
  FIXP_DBL windowNrgF[2][BLOCK_SWITCH_WINDOWS];
  FIXP_DBL accWindowNrg;

  FIXP_DBL iirStates[BLOCK_SWITCHING_IIR_LEN];

} BLOCK_SWITCHING_CONTROL;

void FDKaacEnc_InitBlockSwitching(
    BLOCK_SWITCHING_CONTROL *blockSwitchingControl, INT isLowDelay);

int FDKaacEnc_BlockSwitching(BLOCK_SWITCHING_CONTROL *blockSwitchingControl,
                             const INT granuleLength, const int isLFE,
                             const INT_PCM *pTimeSignal);

int FDKaacEnc_SyncBlockSwitching(
    BLOCK_SWITCHING_CONTROL *blockSwitchingControlLeft,
    BLOCK_SWITCHING_CONTROL *blockSwitchingControlRight, const INT noOfChannels,
    const INT commonWindow);

#endif /* #ifndef BLOCK_SWITCH_H */
