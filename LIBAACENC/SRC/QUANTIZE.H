/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Quantization

*******************************************************************************/

#ifndef QUANTIZE_H
#define QUANTIZE_H

#include "common_fix.h"

#define MAX_QUANT 8191

void FDKaacEnc_QuantizeSpectrum(INT sfbCnt, INT maxSfbPerGroup, INT sfbPerGroup,
                                const INT *sfbOffset,
                                const FIXP_DBL *mdctSpectrum, INT globalGain,
                                const INT *scalefactors,
                                SHORT *quantizedSpectrum, INT dZoneQuantEnable);

FIXP_DBL FDKaacEnc_calcSfbDist(const FIXP_DBL *mdctSpectrum,
                               SHORT *quantSpectrum, INT noOfLines, INT gain,
                               INT dZoneQuantEnable);

void FDKaacEnc_calcSfbQuantEnergyAndDist(FIXP_DBL *mdctSpectrum,
                                         SHORT *quantSpectrum, INT noOfLines,
                                         INT gain, FIXP_DBL *en,
                                         FIXP_DBL *dist);

#endif /* QUANTIZE_H */
