/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Chaos measure calculation

*******************************************************************************/

#ifndef CHAOSMEASURE_H
#define CHAOSMEASURE_H

#include "common_fix.h"
#include "psy_const.h"

void FDKaacEnc_CalculateChaosMeasure(FIXP_DBL *paMDCTDataNM0, INT numberOfLines,
                                     FIXP_DBL *chaosMeasure);

#endif /* CHAOSMEASURE_H */
