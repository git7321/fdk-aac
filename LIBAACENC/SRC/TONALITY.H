/**************************** AAC encoder library ******************************

   Author(s):   M. Lohwasser

   Description: Calculate tonality index

*******************************************************************************/

#ifndef TONALITY_H
#define TONALITY_H

#include "common_fix.h"
#include "chaosmeasure.h"

void FDKaacEnc_CalculateFullTonality(FIXP_DBL *RESTRICT spectrum,
                                     INT *RESTRICT sfbMaxScaleSpec,
                                     FIXP_DBL *RESTRICT sfbEnergyLD64,
                                     FIXP_SGL *RESTRICT sfbTonality, INT sfbCnt,
                                     const INT *sfbOffset, INT usePns);

#endif /* TONALITY_H */
