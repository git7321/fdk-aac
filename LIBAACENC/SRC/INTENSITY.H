/**************************** AAC encoder library ******************************

   Author(s):   A. Horndasch (code originally from lwr and rtb) / Josef Hoepfl
(FDK)

   Description: intensity stereo prototype

*******************************************************************************/

#ifndef INTENSITY_H
#define INTENSITY_H

#include "aacenc_pns.h"
#include "common_fix.h"

void FDKaacEnc_IntensityStereoProcessing(
    FIXP_DBL *sfbEnergyLeft, FIXP_DBL *sfbEnergyRight,
    FIXP_DBL *mdctSpectrumLeft, FIXP_DBL *mdctSpectrumRight,
    FIXP_DBL *sfbThresholdLeft, FIXP_DBL *sfbThresholdRight,
    FIXP_DBL *sfbThresholdLdDataRight, FIXP_DBL *sfbSpreadEnLeft,
    FIXP_DBL *sfbSpreadEnRight, FIXP_DBL *sfbEnergyLdDataLeft,
    FIXP_DBL *sfbEnergyLdDataRight, INT *msDigest, INT *msMask,
    const INT sfbCnt, const INT sfbPerGroup, const INT maxSfbPerGroup,
    const INT *sfbOffset, const INT allowIS, INT *isBook, INT *isScale,
    PNS_DATA *RESTRICT pnsData[2]);

#endif /* INTENSITY_H */
