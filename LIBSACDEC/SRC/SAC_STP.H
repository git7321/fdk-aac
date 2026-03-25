/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Dec subband processing

*******************************************************************************/

#ifndef SAC_STP_H
#define SAC_STP_H

#include "sac_dec.h"

SACDEC_ERROR subbandTPCreate(HANDLE_STP_DEC *hStpDec);

SACDEC_ERROR subbandTPInit(HANDLE_STP_DEC self);

SACDEC_ERROR subbandTPApply(spatialDec *self, const SPATIAL_BS_FRAME *frame);
void subbandTPDestroy(HANDLE_STP_DEC *hStpDec);

#endif
