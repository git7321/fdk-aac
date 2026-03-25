/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Short block grouping

*******************************************************************************/

#ifndef GRP_DATA_H
#define GRP_DATA_H

#include "common_fix.h"
#include "psy_data.h"

void FDKaacEnc_groupShortData(FIXP_DBL *mdctSpectrum,
                              SFB_THRESHOLD *sfbThreshold,
                              SFB_ENERGY *sfbEnergy,
                              SFB_ENERGY *sfbEnergyMS,
                              SFB_ENERGY *sfbSpreadEnergy, const INT sfbCnt,
                              const INT sfbActive, const INT *sfbOffset,
                              const FIXP_DBL *sfbMinSnrLdData,
                              INT *groupedSfbOffset,
                              INT *maxSfbPerGroup,
                              FIXP_DBL *groupedSfbMinSnrLdData,
                              const INT noOfGroups, const INT *groupLen,
                              const INT granuleLength);

#endif /* _INTERFACE_H */
