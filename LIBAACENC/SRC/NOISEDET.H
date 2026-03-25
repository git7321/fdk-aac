/**************************** AAC encoder library ******************************

   Author(s):   M. Lohwasser

   Description: noisedet.h

*******************************************************************************/

#ifndef NOISEDET_H
#define NOISEDET_H

#include "common_fix.h"

#include "pnsparam.h"
#include "psy_data.h"

void FDKaacEnc_noiseDetect(FIXP_DBL *mdctSpectrum, INT *sfbMaxScaleSpec,
                           INT sfbActive, const INT *sfbOffset,
                           FIXP_SGL noiseFuzzyMeasure[], NOISEPARAMS *np,
                           FIXP_SGL *sfbtonality);

#endif /* NOISEDET_H */
