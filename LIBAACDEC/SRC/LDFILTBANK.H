/**************************** AAC decoder library ******************************

   Author(s):

   Description: low delay filterbank interface

*******************************************************************************/

#ifndef LDFILTBANK_H
#define LDFILTBANK_H

#include "common_fix.h"
#include "aac_rom.h"

int InvMdctTransformLowDelay_fdk(FIXP_DBL *mdctdata_m, const int mdctdata_e,
                                 PCM_DEC *mdctOut, FIXP_DBL *fs_buffer,
                                 const int frameLength);

#endif
