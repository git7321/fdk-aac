/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Dec M1 and M2 calculation

*******************************************************************************/

#ifndef SAC_CALCM1ANDM2_H
#define SAC_CALCM1ANDM2_H

#include "sac_dec.h"

#define SCALE_PARAM_M1 3

#define SCALE_PARAM_CALC_M2 (3)
#define SCALE_PARAM_M2_515X (3)
#define SCALE_PARAM_M2_525 (SCALE_PARAM_M1 + HRG_SF + 1 - SCALE_PARAM_CALC_M2)
#define SCALE_PARAM_M2_212_PRED (3)
#define SCALE_DATA_APPLY_M2 (1)
#define SCALE_DATA_APPLY_M2_PC (2)

SACDEC_ERROR initM1andM2(spatialDec* self, int initStatesFlag,
                         int configChanged);

int SpatialDecGetResidualIndex(spatialDec* self, int row);

SACDEC_ERROR SpatialDecCalculateM1andM2(spatialDec* self, int ps,
                                        const SPATIAL_BS_FRAME* frame);

#endif /* SAC_CALCM1ANDM2_H */
