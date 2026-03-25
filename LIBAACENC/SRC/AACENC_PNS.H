/**************************** AAC encoder library ******************************

   Author(s):   M. Lohwasser

   Description: pns.h

*******************************************************************************/

#ifndef AACENC_PNS_H
#define AACENC_PNS_H

#include "common_fix.h"
#include "pnsparam.h"

#define NO_NOISE_PNS FDK_INT_MIN

typedef struct {
  NOISEPARAMS np;
  FIXP_DBL minCorrelationEnergy;
  FIXP_DBL noiseCorrelationThresh;
  INT usePns;
} PNS_CONFIG;

typedef struct {
  FIXP_SGL noiseFuzzyMeasure[MAX_GROUPED_SFB];
  FIXP_DBL noiseEnergyCorrelation[MAX_GROUPED_SFB];
  INT pnsFlag[MAX_GROUPED_SFB];
} PNS_DATA;

#endif /* AACENC_PNS_H */
