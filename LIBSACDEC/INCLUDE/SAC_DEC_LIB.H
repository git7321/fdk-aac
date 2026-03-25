/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: Space Decoder

*******************************************************************************/

#ifndef SAC_DEC_LIB_H
#define SAC_DEC_LIB_H

#include "common_fix.h"
#include "FDK_audio.h"
#include "sac_dec_errorcodes.h"
#include "FDK_bitstream.h"
#include "FDK_qmf_domain.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  SAC_INTERFACE_QMF =
      0,
  SAC_INTERFACE_TIME,
  SAC_INTERFACE_AUTO
} SAC_INPUT_CONFIG;

typedef enum {
  SACDEC_OUT_MODE_NORMAL =
      0,
  SACDEC_OUT_MODE_BINAURAL,
  SACDEC_OUT_MODE_STEREO,
  SACDEC_OUT_MODE_6CHANNEL
} SAC_DEC_OUTPUT_MODE;

typedef enum {
  SAC_BINAURAL_HRTF_KEMAR = 0,
  SAC_BINAURAL_HRTF_VAST,
  SAC_BINAURAL_HRTF_MPSVT,
  SAC_BINAURAL_SINGLE_HRTFS
} SAC_BINAURAL_HRTF_MODEL;

typedef enum {
  SAC_INSTANCE_NOT_FULL_AVAILABLE =
      0,
  SAC_INSTANCE_FULL_AVAILABLE
} SAC_INSTANCE_AVAIL;

typedef enum {
  SACDEC_OUTPUT_MODE = 0x0001,
  SACDEC_BLIND_ENABLE =
      0x0002,
  SACDEC_PARTIALLY_COMPLEX =
      0x0003,
  SACDEC_INTERFACE =
      0x0004,
  SACDEC_BS_DELAY = 0x0005,
  SACDEC_BINAURAL_QUALITY =
      0x0102,
  SACDEC_BINAURAL_DISTANCE = 0x0103,
  SACDEC_BINAURAL_DIALOG_CLARITY =
      0x0104,
  SACDEC_BINAURAL_FRONT_ANGLE = 0x0105,
  SACDEC_BINAURAL_BACK_ANGLE = 0x0106,
  SACDEC_BINAURAL_PRESET = 0x0107,

  SACDEC_BS_INTERRUPTION =
      0x0200,
  SACDEC_CLEAR_HISTORY = 0x0201,
  SACDEC_CONCEAL_NUM_KEEP_FRAMES =
      0x0301,
  SACDEC_CONCEAL_FADE_OUT_SLOPE_LENGTH =
      0x0302,
  SACDEC_CONCEAL_FADE_IN_SLOPE_LENGTH =
      0x0303,
  SACDEC_CONCEAL_NUM_RELEASE_FRAMES =
      0x0304
} SACDEC_PARAM;

#define PCM_MPS LONG

typedef struct MpegSurroundDecoder CMpegSurroundDecoder;

SAC_INSTANCE_AVAIL
mpegSurroundDecoder_IsFullMpegSurroundDecoderInstanceAvailable(
    CMpegSurroundDecoder *pMpegSurroundDecoder);

SACDEC_ERROR mpegSurroundDecoder_Open(
    CMpegSurroundDecoder **pMpegSurroundDecoder, int stereoConfigIndex,
    HANDLE_FDK_QMF_DOMAIN pQmfDomain);
SACDEC_ERROR mpegSurroundDecoder_Init(
    CMpegSurroundDecoder *pMpegSurroundDecoder);
SACDEC_ERROR mpegSurroundDecoder_Config(
    CMpegSurroundDecoder *pMpegSurroundDecoder, HANDLE_FDK_BITSTREAM hBs,
    AUDIO_OBJECT_TYPE coreCodec, INT samplingRate, INT frameSize,
    INT numChannels, INT stereoConfigIndex, INT coreSbrFrameLengthIndex,
    INT configBytes, const UCHAR configMode, UCHAR *configChanged);
SACDEC_ERROR
mpegSurroundDecoder_ConfigureQmfDomain(
    CMpegSurroundDecoder *pMpegSurroundDecoder,
    SAC_INPUT_CONFIG sac_dec_interface, UINT coreSamplingRate,
    AUDIO_OBJECT_TYPE coreCodec);

int mpegSurroundDecoder_ParseNoHeader(
    CMpegSurroundDecoder *pMpegSurroundDecoder, HANDLE_FDK_BITSTREAM hBs,
    int *pMpsDataBits, int fGlobalIndependencyFlag);
int mpegSurroundDecoder_Parse(CMpegSurroundDecoder *pMpegSurroundDecoder,
                              HANDLE_FDK_BITSTREAM hBs, int *pMpsDataBits,
                              AUDIO_OBJECT_TYPE coreCodec, int sampleRate,
                              int frameSize, int fGlobalIndependencyFlag);
int mpegSurroundDecoder_Apply(CMpegSurroundDecoder *pMpegSurroundDecoder,
                              PCM_MPS *input, PCM_MPS *pTimeData,
                              const int timeDataSize, int timeDataFrameSize,
                              int *nChannels, int *frameSize, int sampleRate,
                              AUDIO_OBJECT_TYPE coreCodec,
                              AUDIO_CHANNEL_TYPE channelType[],
                              UCHAR channelIndices[],
                              const FDK_channelMapDescr *const mapDescr,
                              const INT inDataHeadroom, INT *outDataHeadroom);

void mpegSurroundDecoder_Close(CMpegSurroundDecoder *pMpegSurroundDecoder);

SACDEC_ERROR mpegSurroundDecoder_FreeMem(
    CMpegSurroundDecoder *pMpegSurroundDecoder);
SACDEC_ERROR mpegSurroundDecoder_SetParam(
    CMpegSurroundDecoder *pMpegSurroundDecoder, const SACDEC_PARAM param,
    const INT value);

int mpegSurroundDecoder_GetLibInfo(LIB_INFO *libInfo);

UINT mpegSurroundDecoder_GetDelay(const CMpegSurroundDecoder *self);

SACDEC_ERROR mpegSurroundDecoder_IsPseudoLR(
    CMpegSurroundDecoder *pMpegSurroundDecoder, int *bsPseudoLr);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* #ifndef SAC_DEC_LIB_H */
