/**************************** SBR decoder library ******************************

   Author(s):   Manuel Jander, Matthias Hildenbrand

   Description: QMF frequency pre whitening for SBR

*******************************************************************************/

#include "common_fix.h"

#ifndef HFGEN_PREFLAT_H
#define HFGEN_PREFLAT_H

#define GAIN_VEC_EXP 6

void sbrDecoder_calculateGainVec(FIXP_DBL **sourceBufferReal,
                                 FIXP_DBL **sourceBufferImag,
                                 int sourceBuffer_e_overlap,
                                 int sourceBuffer_e_current, int overlap,
                                 FIXP_DBL GainVec[], int GainVec_exp[],
                                 const int numBands, const int startSample,
                                 const int stopSample);

#endif /* __HFGEN_PREFLAT_H */
