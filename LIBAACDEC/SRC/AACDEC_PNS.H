/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: perceptual noise substitution tool

*******************************************************************************/

#ifndef AACDEC_PNS_H
#define AACDEC_PNS_H

#include "common_fix.h"

#define NO_OFBANDS ((8 * 16))

typedef struct {
  UCHAR correlated[NO_OFBANDS];
} CPnsInterChannelData;

typedef struct {
  CPnsInterChannelData *pPnsInterChannelData;
  UCHAR pnsUsed[NO_OFBANDS];
  int CurrentEnergy;
  UCHAR PnsActive;
  INT *currentSeed;
  INT *randomSeed;
} CPnsData;

void CPns_UpdateNoiseState(CPnsData *pPnsData, INT *currentSeed,
                           INT *randomSeed);

void CPns_ResetData(CPnsData *pPnsData,
                    CPnsInterChannelData *pPnsInterChannelData);

#endif /* #ifndef AACDEC_PNS_H */
