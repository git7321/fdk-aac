/*********************** MPEG surround encoder library *************************

   Author(s):   Markus Multrus

   Description: Encoder Library Interface
                constants to MPEG-4 spatial encoder lib

*******************************************************************************/

#ifndef SACENC_CONST_H
#define SACENC_CONST_H

/* Includes ******************************************************************/
#include "machine_type.h"

/* Defines *******************************************************************/
#define NUM_QMF_BANDS 64
#define MAX_QMF_BANDS 128

#define SACENC_MAX_NUM_BOXES 1
#define SACENC_MAX_INPUT_CHANNELS 2
#define SACENC_MAX_OUTPUT_CHANNELS 1

#define SACENC_FLOAT_EPSILON (1e-9f)

#endif /* SACENC_CONST_H */
