/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Psychoaccoustic data

*******************************************************************************/

#ifndef PSY_DATA_H
#define PSY_DATA_H

#include "block_switch.h"
#include "mdct.h"

typedef shouldBeUnion {
  FIXP_DBL Long[MAX_GROUPED_SFB];
  FIXP_DBL Short[TRANS_FAC][MAX_SFB_SHORT];
}
SFB_THRESHOLD;

typedef shouldBeUnion {
  FIXP_DBL Long[MAX_GROUPED_SFB];
  FIXP_DBL Short[TRANS_FAC][MAX_SFB_SHORT];
}
SFB_ENERGY;

typedef shouldBeUnion {
  FIXP_DBL Long[MAX_GROUPED_SFB];
  FIXP_DBL Short[TRANS_FAC][MAX_SFB_SHORT];
}
SFB_LD_ENERGY;

typedef shouldBeUnion {
  INT Long[MAX_GROUPED_SFB];
  INT Short[TRANS_FAC][MAX_SFB_SHORT];
}
SFB_MAX_SCALE;

typedef struct {
  INT_PCM* psyInputBuffer;
  FIXP_DBL overlapAddBuffer[3 * 512 / 2];

  mdct_t mdctPers;
  BLOCK_SWITCHING_CONTROL blockSwitchingControl;
  FIXP_DBL sfbThresholdnm1[MAX_SFB];
  INT mdctScalenm1;
  INT calcPreEcho;
  INT isLFE;
} PSY_STATIC;

typedef struct {
  FIXP_DBL* mdctSpectrum;
  SFB_THRESHOLD sfbThreshold;
  SFB_ENERGY sfbEnergy;
  SFB_LD_ENERGY sfbEnergyLdData;
  SFB_MAX_SCALE sfbMaxScaleSpec;
  SFB_ENERGY sfbEnergyMS;
  FIXP_DBL sfbEnergyMSLdData[MAX_GROUPED_SFB];
  SFB_ENERGY sfbSpreadEnergy;
  INT mdctScale;
  INT groupedSfbOffset[MAX_GROUPED_SFB + 1];
  INT sfbActive;
  INT lowpassLine;
} PSY_DATA;

#endif /* PSY_DATA_H */
