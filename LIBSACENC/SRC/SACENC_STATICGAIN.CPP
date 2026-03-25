/*********************** MPEG surround encoder library *************************

   Author(s):   Christian Goettlinger

   Description: Encoder Library Interface
                gain management of the encoder

*******************************************************************************/

#include "sacenc_staticgain.h"

#define MP4SPACEENC_DMX_GAIN_DEFAULT SACENC_DMXGAIN_3_dB
#define GAINCF_SF (4)
#define GAINCT1(x) FL2FXCONST_DBL(x)
#define GAINCF(x) FL2FXCONST_DBL(x)

#define GAINCT2(x) FL2FXCONST_DBL(x)
#define FX_DBL2FX_GAIN(x) (x)

struct STATIC_GAIN {
  MP4SPACEENC_MODE encMode;
  MP4SPACEENC_DMX_GAIN fixedGainDMX;
  INT preGainFactorDb;

  FIXP_GAIN PostGain__FDK;
  FIXP_GAIN pPreGain__FDK[SACENC_MAX_INPUT_CHANNELS];
};

static const FIXP_DBL preGainFactorTable__FDK[41] = {
    GAINCF(6.2500000931e-003), GAINCF(7.0126154460e-003),
    GAINCF(7.8682834283e-003), GAINCF(8.8283596560e-003),
    GAINCF(9.9055822939e-003), GAINCF(1.1114246212e-002),
    GAINCF(1.2470389716e-002), GAINCF(1.3992006890e-002),
    GAINCF(1.5699289739e-002), GAINCF(1.7614893615e-002),
    GAINCF(1.9764235243e-002), GAINCF(2.2175837308e-002),
    GAINCF(2.4881698191e-002), GAINCF(2.7917724103e-002),
    GAINCF(3.1324200332e-002), GAINCF(3.5146333277e-002),
    GAINCF(3.9434835315e-002), GAINCF(4.4246610254e-002),
    GAINCF(4.9645513296e-002), GAINCF(5.5703181773e-002),
    GAINCF(6.2500000000e-002), GAINCF(7.0126153529e-002),
    GAINCF(7.8682839870e-002), GAINCF(8.8283598423e-002),
    GAINCF(9.9055826664e-002), GAINCF(1.1114246398e-001),
    GAINCF(1.2470389158e-001), GAINCF(1.3992007077e-001),
    GAINCF(1.5699289739e-001), GAINCF(1.7614893615e-001),
    GAINCF(1.9764235616e-001), GAINCF(2.2175836563e-001),
    GAINCF(2.4881698191e-001), GAINCF(2.7917724848e-001),
    GAINCF(3.1324201822e-001), GAINCF(3.5146331787e-001),
    GAINCF(3.9434835315e-001), GAINCF(4.4246610999e-001),
    GAINCF(4.9645513296e-001), GAINCF(5.5703181028e-001),
    GAINCF(6.2500000000e-001)};

static const FIXP_GAIN dmxGainTable__FDK[] = {
    GAINCT2(0.84089650f),
    GAINCT2(0.70710706f),
    GAINCT2(0.59460385f),
    GAINCT2(0.50000000f),
    GAINCT2(0.42044825f),
    GAINCT2(0.35355341f),
    GAINCT2(0.25000000f)};

FDK_SACENC_ERROR fdk_sacenc_staticGain_OpenConfig(
    HANDLE_STATIC_GAIN_CONFIG *phStaticGainConfig) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == phStaticGainConfig) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDK_ALLOCATE_MEMORY_1D(*phStaticGainConfig, 1, struct STATIC_GAIN_CONFIG);
  }
  return error;

bail:
  fdk_sacenc_staticGain_CloseConfig(phStaticGainConfig);
  return ((SACENC_OK == error) ? SACENC_MEMORY_ERROR : error);
}

FDK_SACENC_ERROR fdk_sacenc_staticGain_InitDefaultConfig(
    HANDLE_STATIC_GAIN_CONFIG hStaticGainConfig) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == hStaticGainConfig) {
    error = SACENC_INVALID_HANDLE;
  } else {
    hStaticGainConfig->encMode = SACENC_INVALID_MODE;

    hStaticGainConfig->fixedGainDMX = MP4SPACEENC_DMX_GAIN_DEFAULT;
    hStaticGainConfig->preGainFactorDb = 0;
  }
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_staticGain_CloseConfig(
    HANDLE_STATIC_GAIN_CONFIG *phStaticGainConfig) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((phStaticGainConfig == NULL) || (*phStaticGainConfig == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDKfree(*phStaticGainConfig);
    *phStaticGainConfig = NULL;
  }
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_staticGain_Open(HANDLE_STATIC_GAIN *phStaticGain) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == phStaticGain) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDK_ALLOCATE_MEMORY_1D(*phStaticGain, 1, struct STATIC_GAIN);
  }
  return error;

bail:
  return ((SACENC_OK == error) ? SACENC_MEMORY_ERROR : error);
}

FDK_SACENC_ERROR fdk_sacenc_staticGain_Init(
    HANDLE_STATIC_GAIN hStaticGain,
    const HANDLE_STATIC_GAIN_CONFIG hStaticGainConfig, INT *const scale) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((hStaticGain == NULL) || (hStaticGainConfig == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    hStaticGain->encMode = hStaticGainConfig->encMode;
    hStaticGain->fixedGainDMX = hStaticGainConfig->fixedGainDMX;
    hStaticGain->preGainFactorDb = hStaticGainConfig->preGainFactorDb;

    if ((hStaticGain->preGainFactorDb < -20) ||
        (hStaticGain->preGainFactorDb > 20)) {
      error = SACENC_INVALID_CONFIG;
      goto bail;
    }

    FIXP_DBL fPreGainFactor__FDK;

    if (hStaticGain->preGainFactorDb == 0) {
      fPreGainFactor__FDK = (FIXP_DBL)MAXVAL_DBL;
      *scale = 0;
    } else {
      int s;
      fPreGainFactor__FDK =
          preGainFactorTable__FDK[hStaticGain->preGainFactorDb + 20];
      s = fixMax(0, CntLeadingZeros(fPreGainFactor__FDK) - 1);
      fPreGainFactor__FDK = fPreGainFactor__FDK << (s);
      *scale = GAINCF_SF - s;
    }

    if (hStaticGain->fixedGainDMX == 0)
      hStaticGain->PostGain__FDK = MAXVAL_GAIN;
    else
      hStaticGain->PostGain__FDK =
          dmxGainTable__FDK[hStaticGain->fixedGainDMX - 1];

    FDKmemclear(
        hStaticGain->pPreGain__FDK,
        sizeof(hStaticGain->pPreGain__FDK));

    if (hStaticGain->encMode == SACENC_212) {
      hStaticGain->pPreGain__FDK[0] =
          FX_DBL2FX_GAIN(fPreGainFactor__FDK);
      hStaticGain->pPreGain__FDK[1] =
          FX_DBL2FX_GAIN(fPreGainFactor__FDK);
    } else {
      error = SACENC_INVALID_CONFIG;
    }
  }

bail:

  return error;
}

FDK_SACENC_ERROR fdk_sacenc_staticGain_Close(HANDLE_STATIC_GAIN *phStaticGain) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((phStaticGain == NULL) || (*phStaticGain == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDKfree(*phStaticGain);
    *phStaticGain = NULL;
  }
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_staticPostGain_ApplyFDK(
    const HANDLE_STATIC_GAIN hStaticGain, INT_PCM *const pOutputSamples,
    const INT nOutputSamples, const INT scale) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == hStaticGain) {
    error = SACENC_INVALID_HANDLE;
  } else {
    int i;
    FIXP_GAIN postGain = hStaticGain->PostGain__FDK;

    if (scale < 0) {
      if (postGain == MAXVAL_GAIN) {
        for (i = 0; i < nOutputSamples; i++) {
          pOutputSamples[i] = pOutputSamples[i] >> (-scale);
        }
      } else {
        for (i = 0; i < nOutputSamples; i++) {
          pOutputSamples[i] = FX_DBL2FX_PCM(
              fMult(postGain, FX_PCM2FX_DBL(pOutputSamples[i])) >> (-scale));
        }
      }
    } else {
      if (postGain == MAXVAL_GAIN) {
        for (i = 0; i < nOutputSamples; i++) {
          pOutputSamples[i] = FX_DBL2FX_PCM(SATURATE_LEFT_SHIFT(
              FX_PCM2FX_DBL(pOutputSamples[i]), scale, DFRACT_BITS));
        }
      } else {
        for (i = 0; i < nOutputSamples; i++) {
          pOutputSamples[i] = FX_DBL2FX_PCM(SATURATE_LEFT_SHIFT(
              fMult(postGain, FX_PCM2FX_DBL(pOutputSamples[i])), scale,
              DFRACT_BITS));
        }
      }
    }
  }
  return error;
}

FIXP_GAIN *fdk_sacenc_getPreGainPtrFDK(HANDLE_STATIC_GAIN hStaticGain) {
  return ((hStaticGain == NULL) ? NULL : hStaticGain->pPreGain__FDK);
}

FIXP_GAIN fdk_sacenc_getPostGainFDK(HANDLE_STATIC_GAIN hStaticGain) {
  return (hStaticGain->PostGain__FDK);
}

FIXEDGAINDMXCONFIG fdk_sacenc_staticGain_GetDmxGain(
    const HANDLE_STATIC_GAIN hStaticGain) {
  FIXEDGAINDMXCONFIG dmxGain = FIXEDGAINDMX_INVALID;

  switch (hStaticGain->fixedGainDMX) {
    case 0:
      dmxGain = FIXEDGAINDMX_0;
      break;
    case 1:
      dmxGain = FIXEDGAINDMX_1;
      break;
    case 2:
      dmxGain = FIXEDGAINDMX_2;
      break;
    case 3:
      dmxGain = FIXEDGAINDMX_3;
      break;
    case 4:
      dmxGain = FIXEDGAINDMX_4;
      break;
    case 5:
      dmxGain = FIXEDGAINDMX_5;
      break;
    case 6:
      dmxGain = FIXEDGAINDMX_6;
      break;
    case 7:
      dmxGain = FIXEDGAINDMX_7;
      break;
    default:
      dmxGain = FIXEDGAINDMX_INVALID;
  }
  return dmxGain;
}

FDK_SACENC_ERROR fdk_sacenc_staticGain_SetDmxGain(
    HANDLE_STATIC_GAIN_CONFIG hStaticGainCfg,
    const MP4SPACEENC_DMX_GAIN dmxGain) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == hStaticGainCfg) {
    error = SACENC_INVALID_HANDLE;
  } else {
    hStaticGainCfg->fixedGainDMX = dmxGain;
  }
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_staticGain_SetEncMode(
    HANDLE_STATIC_GAIN_CONFIG hStaticGainCfg, const MP4SPACEENC_MODE encMode) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == hStaticGainCfg) {
    error = SACENC_INVALID_HANDLE;
  } else {
    hStaticGainCfg->encMode = encMode;
  }
  return error;
}
