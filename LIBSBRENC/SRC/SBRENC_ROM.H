/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef SBRENC_ROM_H
#define SBRENC_ROM_H

#include "sbr_def.h"
#include "sbr_encoder.h"

#include "ps_main.h"

extern const INT v_Huff_envelopeLevelC10T[121];
extern const UCHAR v_Huff_envelopeLevelL10T[121];
extern const INT v_Huff_envelopeLevelC10F[121];
extern const UCHAR v_Huff_envelopeLevelL10F[121];
extern const INT bookSbrEnvBalanceC10T[49];
extern const UCHAR bookSbrEnvBalanceL10T[49];
extern const INT bookSbrEnvBalanceC10F[49];
extern const UCHAR bookSbrEnvBalanceL10F[49];
extern const INT v_Huff_envelopeLevelC11T[63];
extern const UCHAR v_Huff_envelopeLevelL11T[63];
extern const INT v_Huff_envelopeLevelC11F[63];
extern const UCHAR v_Huff_envelopeLevelL11F[63];
extern const INT bookSbrEnvBalanceC11T[25];
extern const UCHAR bookSbrEnvBalanceL11T[25];
extern const INT bookSbrEnvBalanceC11F[25];
extern const UCHAR bookSbrEnvBalanceL11F[25];
extern const INT v_Huff_NoiseLevelC11T[63];
extern const UCHAR v_Huff_NoiseLevelL11T[63];
extern const INT bookSbrNoiseBalanceC11T[25];
extern const UCHAR bookSbrNoiseBalanceL11T[25];

extern const sbrTuningTable_t sbrTuningTable[];
extern const int sbrTuningTableSize;

extern const psTuningTable_t psTuningTable[4];

#endif
