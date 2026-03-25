/**************************** PCM utility library ******************************

   Author(s):   Christian Griebel

   Description:

*******************************************************************************/

#ifndef PCMDMX_LIB_H
#define PCMDMX_LIB_H

#include "machine_type.h"
#include "common_fix.h"
#include "FDK_audio.h"
#include "FDK_bitstream.h"

typedef enum {
  PCMDMX_OK = 0x0,
  PCMDMX_UNSUPPORTED =
      0x1,
  pcm_dmx_fatal_error_start,
  PCMDMX_OUT_OF_MEMORY,
  pcm_dmx_fatal_error_end,
  PCMDMX_INVALID_HANDLE,
  PCMDMX_INVALID_ARGUMENT,
  PCMDMX_INVALID_CH_CONFIG,
  PCMDMX_INVALID_MODE,
  PCMDMX_UNKNOWN_PARAM,
  PCMDMX_UNABLE_TO_SET_PARAM,
  PCMDMX_CORRUPT_ANC_DATA,
  PCMDMX_OUTPUT_BUFFER_TOO_SMALL
} PCMDMX_ERROR;

#define PCMDMX_IS_FATAL_ERROR(err)          \
  ((((err) >= pcm_dmx_fatal_error_start) && \
    ((err) <= pcm_dmx_fatal_error_end))     \
       ? 1                                  \
       : 0)

typedef enum {
  DMX_PROFILE_SETTING =
      0x01,
  DMX_BS_DATA_EXPIRY_FRAME =
      0x10,
  DMX_BS_DATA_DELAY =
      0x11,
  MIN_NUMBER_OF_OUTPUT_CHANNELS =
      0x20,
  MAX_NUMBER_OF_OUTPUT_CHANNELS =
      0x21,
  DMX_DUAL_CHANNEL_MODE =
      0x30,
  DMX_PSEUDO_SURROUND_MODE =
      0x31
} PCMDMX_PARAM;

typedef enum {
  DMX_PRFL_STANDARD =
      0x0,
  DMX_PRFL_MATRIX_MIX =
      0x1,
  DMX_PRFL_FORCE_MATRIX_MIX =
      0x2,
  DMX_PRFL_ARIB_JAPAN =
      0x3
} DMX_PROFILE_TYPE;

typedef enum {
  NEVER_DO_PS_DMX =
      -1,
  AUTO_PS_DMX = 0,
  FORCE_PS_DMX =
      1
} PSEUDO_SURROUND_MODE;

typedef enum {
  STEREO_MODE = 0x0,
  CH1_MODE = 0x1,
  CH2_MODE = 0x2,
  MIXED_MODE = 0x3
} DUAL_CHANNEL_MODE;

#define DMX_PCM FIXP_DBL
#define DMX_PCMF FIXP_DBL
#define DMX_PCM_BITS DFRACT_BITS
#define FX_DMX2FX_PCM(x) FX_DBL2FX_PCM((FIXP_DBL)(x))

typedef struct PCM_DMX_INSTANCE *HANDLE_PCM_DOWNMIX;

#define PCMDMX_RESET_PARAMS (1)
#define PCMDMX_RESET_BS_DATA (2)
#define PCMDMX_RESET_FULL (PCMDMX_RESET_PARAMS | PCMDMX_RESET_BS_DATA)

#ifdef __cplusplus
extern "C" {
#endif

PCMDMX_ERROR pcmDmx_Open(HANDLE_PCM_DOWNMIX *pSelf);
PCMDMX_ERROR pcmDmx_SetParam(HANDLE_PCM_DOWNMIX self, const PCMDMX_PARAM param,
                             const INT value);
PCMDMX_ERROR pcmDmx_GetParam(HANDLE_PCM_DOWNMIX self, const PCMDMX_PARAM param,
                             INT *const pValue);
PCMDMX_ERROR pcmDmx_Parse(HANDLE_PCM_DOWNMIX self,
                          HANDLE_FDK_BITSTREAM hBitStream, UINT ancDataBits,
                          int isMpeg2);
PCMDMX_ERROR pcmDmx_ReadDvbAncData(HANDLE_PCM_DOWNMIX self, UCHAR *pAncDataBuf,
                                   UINT ancDataBytes, int isMpeg2);

PCMDMX_ERROR pcmDmx_SetMatrixMixdownFromPce(HANDLE_PCM_DOWNMIX self,
                                            int matrixMixdownPresent,
                                            int matrixMixdownIdx,
                                            int pseudoSurroundEnable);
PCMDMX_ERROR pcmDmx_Reset(HANDLE_PCM_DOWNMIX self, UINT flags);
PCMDMX_ERROR pcmDmx_ApplyFrame(HANDLE_PCM_DOWNMIX self, DMX_PCM *pPcmBuf,
                               const int pcmBufSize, UINT frameSize,
                               INT *nChannels, INT fInterleaved,
                               AUDIO_CHANNEL_TYPE channelType[],
                               UCHAR channelIndices[],
                               const FDK_channelMapDescr *const mapDescr,
                               INT *pDmxOutScale);
PCMDMX_ERROR pcmDmx_Close(HANDLE_PCM_DOWNMIX *pSelf);
PCMDMX_ERROR pcmDmx_GetLibInfo(LIB_INFO *info);

#ifdef __cplusplus
}
#endif

#endif /* PCMDMX_LIB_H */
