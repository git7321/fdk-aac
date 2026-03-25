/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Dec guided envelope shaping

*******************************************************************************/

#ifndef SAC_RESHAPEBBENV_H
#define SAC_RESHAPEBBENV_H

#include "sac_dec_interface.h"

#define BB_ENV_SIZE 9

void initBBEnv(spatialDec *self, int initStatesFlag);
void SpatialDecReshapeBBEnv(spatialDec *self, const SPATIAL_BS_FRAME *frame,
                            int ts);

#endif
