/*********************** MPEG surround encoder library *************************

   Author(s):   Christian Goettlinger

   Description: Encoder Library Interfac
                gain management of the encoder

*******************************************************************************/

#ifndef SACENC_STATICGAIN_H
#define SACENC_STATICGAIN_H

#include "common_fix.h"
#include "sacenc_lib.h"
#include "sacenc_const.h"
#include "sacenc_bitstream.h"

#define FIXP_GAIN FIXP_DBL
#define MAXVAL_GAIN ((FIXP_DBL)MAXVAL_DBL)

struct STATIC_GAIN_CONFIG {
  MP4SPACEENC_MODE encMode;
  MP4SPACEENC_DMX_GAIN fixedGainDMX;
  INT preGainFactorDb;
};

typedef struct STATIC_GAIN_CONFIG *HANDLE_STATIC_GAIN_CONFIG;
typedef struct STATIC_GAIN *HANDLE_STATIC_GAIN;

FDK_SACENC_ERROR fdk_sacenc_staticGain_OpenConfig(
    HANDLE_STATIC_GAIN_CONFIG *phStaticGainConfig);

FDK_SACENC_ERROR fdk_sacenc_staticGain_InitDefaultConfig(
    HANDLE_STATIC_GAIN_CONFIG hStaticGainConfig);

FDK_SACENC_ERROR fdk_sacenc_staticGain_CloseConfig(
    HANDLE_STATIC_GAIN_CONFIG *phStaticGainConfig);

FDK_SACENC_ERROR fdk_sacenc_staticGain_Open(HANDLE_STATIC_GAIN *phStaticGain);

FDK_SACENC_ERROR fdk_sacenc_staticGain_Init(
    HANDLE_STATIC_GAIN hStaticGain,
    const HANDLE_STATIC_GAIN_CONFIG hStaticGainConfig, INT *const scale);

FDK_SACENC_ERROR fdk_sacenc_staticGain_Close(HANDLE_STATIC_GAIN *phStaticGain);

FDK_SACENC_ERROR fdk_sacenc_staticPostGain_ApplyFDK(
    const HANDLE_STATIC_GAIN hStaticGain, INT_PCM *const pOutputSamples,
    const INT nOutputSamples, const INT scale);

FIXP_GAIN *fdk_sacenc_getPreGainPtrFDK(HANDLE_STATIC_GAIN hStaticGain);

FIXP_GAIN fdk_sacenc_getPostGainFDK(HANDLE_STATIC_GAIN hStaticGain);

FIXEDGAINDMXCONFIG fdk_sacenc_staticGain_GetDmxGain(
    const HANDLE_STATIC_GAIN hStaticGain);

FDK_SACENC_ERROR fdk_sacenc_staticGain_SetDmxGain(
    HANDLE_STATIC_GAIN_CONFIG hStaticGainCfg,
    const MP4SPACEENC_DMX_GAIN dmxGain);

FDK_SACENC_ERROR fdk_sacenc_staticGain_SetEncMode(
    HANDLE_STATIC_GAIN_CONFIG hStaticGainCfg, const MP4SPACEENC_MODE encMode);

#endif /* SACENC_STATICGAIN_H */
