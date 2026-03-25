/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Perceptual entropie module

*******************************************************************************/

#ifndef LINE_PE_H
#define LINE_PE_H

#include "common_fix.h"
#include "psy_const.h"

#define PE_CONSTPART_SHIFT FRACT_BITS

typedef struct {
  INT sfbNLines[MAX_GROUPED_SFB];
  INT sfbPe[MAX_GROUPED_SFB];
  INT sfbConstPart[MAX_GROUPED_SFB];
  INT sfbNActiveLines[MAX_GROUPED_SFB];
  INT pe;
  INT constPart;
  INT nActiveLines;
} PE_CHANNEL_DATA;

typedef struct {
  PE_CHANNEL_DATA peChannelData[(2)];
  INT pe;
  INT constPart;
  INT nActiveLines;
  INT offset;
} PE_DATA;

void FDKaacEnc_prepareSfbPe(PE_CHANNEL_DATA *RESTRICT const peChanData,
                            const FIXP_DBL *RESTRICT const sfbEnergyLdData,
                            const FIXP_DBL *RESTRICT const sfbThresholdLdData,
                            const FIXP_DBL *RESTRICT const sfbFormFactorLdData,
                            const INT *RESTRICT const sfbOffset,
                            const INT sfbCnt, const INT sfbPerGroup,
                            const INT maxSfbPerGroup);

void FDKaacEnc_calcSfbPe(PE_CHANNEL_DATA *RESTRICT const peChanData,
                         const FIXP_DBL *RESTRICT const sfbEnergyLdData,
                         const FIXP_DBL *RESTRICT const sfbThresholdLdData,
                         const INT sfbCnt, const INT sfbPerGroup,
                         const INT maxSfbPerGroup,
                         const INT *RESTRICT const isBook,
                         const INT *RESTRICT const isScale);

#endif
