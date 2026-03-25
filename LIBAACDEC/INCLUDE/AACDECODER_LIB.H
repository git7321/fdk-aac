#ifndef AACDECODER_LIB_H
#define AACDECODER_LIB_H

#include "..\..\libsys\include\machine_type.h"
#include "..\..\libsys\include\FDK_audio.h"

#define AACDECODER_LIB_VL0 3
#define AACDECODER_LIB_VL1 2
#define AACDECODER_LIB_VL2 0

#include "..\..\libsys\include\genericStds.h"
typedef enum {
  AAC_DEC_OK =
      0x0000,
  AAC_DEC_OUT_OF_MEMORY =
      0x0002,
  AAC_DEC_UNKNOWN =
      0x0005,

  aac_dec_sync_error_start = 0x1000,
  AAC_DEC_TRANSPORT_SYNC_ERROR = 0x1001,
  AAC_DEC_NOT_ENOUGH_BITS = 0x1002,
  aac_dec_sync_error_end = 0x1FFF,

  aac_dec_init_error_start = 0x2000,
  AAC_DEC_INVALID_HANDLE =
      0x2001,
  AAC_DEC_UNSUPPORTED_AOT =
      0x2002,
  AAC_DEC_UNSUPPORTED_FORMAT =
      0x2003,
  AAC_DEC_UNSUPPORTED_ER_FORMAT =
      0x2004,
  AAC_DEC_UNSUPPORTED_EPCONFIG =
      0x2005,
  AAC_DEC_UNSUPPORTED_MULTILAYER =
      0x2006,
  AAC_DEC_UNSUPPORTED_CHANNELCONFIG =
      0x2007,
  AAC_DEC_UNSUPPORTED_SAMPLINGRATE = 0x2008,
  AAC_DEC_INVALID_SBR_CONFIG =
      0x2009,
  AAC_DEC_SET_PARAM_FAIL = 0x200A,
  AAC_DEC_NEED_TO_RESTART = 0x200B,
  AAC_DEC_OUTPUT_BUFFER_TOO_SMALL =
      0x200C,
  aac_dec_init_error_end = 0x2FFF,

  aac_dec_decode_error_start = 0x4000,
  AAC_DEC_TRANSPORT_ERROR =
      0x4001,
  AAC_DEC_PARSE_ERROR = 0x4002,
  AAC_DEC_UNSUPPORTED_EXTENSION_PAYLOAD =
      0x4003,
  AAC_DEC_DECODE_FRAME_ERROR = 0x4004,
  AAC_DEC_CRC_ERROR = 0x4005,
  AAC_DEC_INVALID_CODE_BOOK = 0x4006,
  AAC_DEC_UNSUPPORTED_PREDICTION =
      0x4007,
  AAC_DEC_UNSUPPORTED_CCE = 0x4008,
  AAC_DEC_UNSUPPORTED_LFE = 0x4009,
  AAC_DEC_UNSUPPORTED_GAIN_CONTROL_DATA =
      0x400A,
  AAC_DEC_UNSUPPORTED_SBA =
      0x400B,
  AAC_DEC_TNS_READ_ERROR = 0x400C,
  AAC_DEC_RVLC_ERROR =
      0x400D,
  aac_dec_decode_error_end = 0x4FFF,
  aac_dec_anc_data_error_start = 0x8000,
  AAC_DEC_ANC_DATA_ERROR =
      0x8001,
  AAC_DEC_TOO_SMALL_ANC_BUFFER = 0x8002,
  AAC_DEC_TOO_MANY_ANC_ELEMENTS = 0x8003,
  aac_dec_anc_data_error_end = 0x8FFF

} AAC_DECODER_ERROR;

#define IS_INIT_ERROR(err)                                                    \
  ((((err) >= aac_dec_init_error_start) && ((err) <= aac_dec_init_error_end)) \
       ? 1                                                                    \
       : 0)
#define IS_DECODE_ERROR(err)                 \
  ((((err) >= aac_dec_decode_error_start) && \
    ((err) <= aac_dec_decode_error_end))     \
       ? 1                                   \
       : 0)
#define IS_OUTPUT_VALID(err) (((err) == AAC_DEC_OK) || IS_DECODE_ERROR(err))

typedef enum {
  AAC_MD_PROFILE_MPEG_STANDARD =
      0,
  AAC_MD_PROFILE_MPEG_LEGACY =
      1,
  AAC_MD_PROFILE_MPEG_LEGACY_PRIO =
      2,
  AAC_MD_PROFILE_ARIB_JAPAN =
      3
} AAC_MD_PROFILE;

typedef enum {
  AAC_DRC_PARAMETER_HANDLING_DISABLED = -1,
  AAC_DRC_PARAMETER_HANDLING_ENABLED =
      0,
  AAC_DRC_PRESENTATION_MODE_1_DEFAULT =
      1,
  AAC_DRC_PRESENTATION_MODE_2_DEFAULT =
      2
} AAC_DRC_DEFAULT_PRESENTATION_MODE_OPTIONS;

typedef enum {
  AAC_PCM_DUAL_CHANNEL_OUTPUT_MODE =
      0x0002,
  AAC_PCM_OUTPUT_CHANNEL_MAPPING =
      0x0003,
  AAC_PCM_LIMITER_ENABLE =
      0x0004,
  AAC_PCM_LIMITER_ATTACK_TIME = 0x0005,
  AAC_PCM_LIMITER_RELEAS_TIME = 0x0006,
  AAC_PCM_MIN_OUTPUT_CHANNELS =
      0x0011,
  AAC_PCM_MAX_OUTPUT_CHANNELS =
      0x0012,
  AAC_METADATA_PROFILE =
      0x0020,
  AAC_METADATA_EXPIRY_TIME = 0x0021,

  AAC_CONCEAL_METHOD = 0x0100,
  AAC_DRC_BOOST_FACTOR =
      0x0200,
  AAC_DRC_ATTENUATION_FACTOR = 0x0201,
  AAC_DRC_REFERENCE_LEVEL =
      0x0202,
  AAC_DRC_HEAVY_COMPRESSION =
      0x0203,
  AAC_DRC_DEFAULT_PRESENTATION_MODE =
      0x0204,
  AAC_DRC_ENC_TARGET_LEVEL =
      0x0205,
  AAC_UNIDRC_SET_EFFECT = 0x0206,
  AAC_UNIDRC_ALBUM_MODE =
      0x0207,
  AAC_QMF_LOWPOWER =
      0x0300,
  AAC_TPDEC_CLEAR_BUFFER =
      0x0603
} AACDEC_PARAM;

typedef struct {
  INT sampleRate;
  INT frameSize;
  INT numChannels;
  AUDIO_CHANNEL_TYPE
  *pChannelType;
  UCHAR *pChannelIndices;
  INT aacSampleRate;
  INT profile;
  AUDIO_OBJECT_TYPE
  aot;
  INT channelConfig;
  INT bitRate;
  INT aacSamplesPerFrame;
  INT aacNumChannels;
  AUDIO_OBJECT_TYPE extAot;
  INT extSamplingRate;

  UINT outputDelay;
  UINT flags;

  SCHAR epConfig;
  INT numLostAccessUnits;

  INT64 numTotalBytes;
  INT64
  numBadBytes;
  INT64
  numTotalAccessUnits;
  INT64 numBadAccessUnits;

  SCHAR drcProgRefLev;
  SCHAR
  drcPresMode;
  INT outputLoudness;

} CStreamInfo;

typedef struct AAC_DECODER_INSTANCE
    *HANDLE_AACDECODER;

#ifdef __cplusplus
extern "C" {
#endif

LINKSPEC(AAC_DECODER_ERROR) aacDecoder_AncDataInit(HANDLE_AACDECODER self,
                                                    UCHAR *buffer, int size);

LINKSPEC(AAC_DECODER_ERROR) aacDecoder_AncDataGet(HANDLE_AACDECODER self,
                                                   int index, UCHAR **ptr,
                                                   int *size);

LINKSPEC(AAC_DECODER_ERROR) aacDecoder_SetParam(const HANDLE_AACDECODER self,
                                                 const AACDEC_PARAM param,
                                                 const INT value);

LINKSPEC(AAC_DECODER_ERROR)
aacDecoder_GetFreeBytes(const HANDLE_AACDECODER self, UINT *pFreeBytes);

LINKSPEC(HANDLE_AACDECODER) aacDecoder_Open(TRANSPORT_TYPE transportFmt,
                                             UINT nrOfLayers);

LINKSPEC(AAC_DECODER_ERROR) aacDecoder_ConfigRaw(HANDLE_AACDECODER self,
                                                  UCHAR *conf[],
                                                  const UINT length[]);

LINKSPEC(AAC_DECODER_ERROR) aacDecoder_RawISOBMFFData(HANDLE_AACDECODER self,
                                                       UCHAR *buffer,
                                                       UINT length);

LINKSPEC(AAC_DECODER_ERROR) aacDecoder_Fill(HANDLE_AACDECODER self,
                                             UCHAR *pBuffer[],
                                             const UINT bufferSize[],
                                             UINT *bytesValid);

#define AACDEC_CONCEAL 1
#define AACDEC_FLUSH 2
#define AACDEC_INTR 4
#define AACDEC_CLRHIST 8

LINKSPEC(AAC_DECODER_ERROR) aacDecoder_DecodeFrame(HANDLE_AACDECODER self,
                                                    INT_PCM *pTimeData,
                                                    const INT timeDataSize,
                                                    const UINT flags);

LINKSPEC(void) aacDecoder_Close(HANDLE_AACDECODER self);
LINKSPEC(CStreamInfo *) aacDecoder_GetStreamInfo(HANDLE_AACDECODER self);
LINKSPEC(INT) aacDecoder_GetLibInfo(LIB_INFO *info);

#ifdef __cplusplus
}
#endif

#endif /* AACDECODER_LIB_H */
