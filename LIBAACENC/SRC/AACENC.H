/**************************** AAC encoder library ******************************

   Author(s):   M. Schug / A. Groeschel

   Description: fast aac coder interface library functions

*******************************************************************************/

#ifndef AACENC_H
#define AACENC_H

#include "common_fix.h"
#include "FDK_audio.h"

#include "tpenc_lib.h"

#include "sbr_encoder.h"

#define MIN_BUFSIZE_PER_EFF_CHAN 6144

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  AAC_ENC_OK = 0x0000,

  AAC_ENC_UNKNOWN = 0x0002,

  /* initialization errors */
  aac_enc_init_error_start = 0x2000,
  AAC_ENC_INVALID_HANDLE = 0x2020,
  AAC_ENC_INVALID_FRAME_LENGTH =
      0x2080, /*!< Invalid frame length (must be 1024 or 960). */
  AAC_ENC_INVALID_N_CHANNELS =
      0x20e0, /*!< Invalid amount of audio input channels. */
  AAC_ENC_INVALID_SFB_TABLE = 0x2140, /*!< Internal encoder error. */

  AAC_ENC_UNSUPPORTED_AOT =
      0x3000, /*!< The Audio Object Type (AOT) is not supported. */
  AAC_ENC_UNSUPPORTED_FILTERBANK =
      0x3010, /*!< Filterbank type is not supported. */
  AAC_ENC_UNSUPPORTED_BITRATE =
      0x3020, /*!< The chosen bitrate is not supported. */
  AAC_ENC_UNSUPPORTED_BITRATE_MODE =
      0x3028, /*!< Unsupported bit rate mode (CBR or VBR). */
  AAC_ENC_UNSUPPORTED_ANC_BITRATE =
      0x3040, /*!< Unsupported ancillay bitrate. */
  AAC_ENC_UNSUPPORTED_ANC_MODE = 0x3060,
  AAC_ENC_UNSUPPORTED_TRANSPORT_TYPE =
      0x3080, /*!< The bitstream format is not supported. */
  AAC_ENC_UNSUPPORTED_ER_FORMAT =
      0x30a0, /*!< The error resilience tool format is not supported. */
  AAC_ENC_UNSUPPORTED_EPCONFIG =
      0x30c0, /*!< The error protection format is not supported. */
  AAC_ENC_UNSUPPORTED_CHANNELCONFIG =
      0x30e0, /*!< The channel configuration (either number or arrangement) is
                 not supported. */
  AAC_ENC_UNSUPPORTED_SAMPLINGRATE =
      0x3100, /*!< Sample rate of audio input is not supported. */
  AAC_ENC_NO_MEMORY = 0x3120,               /*!< Could not allocate memory. */
  AAC_ENC_PE_INIT_TABLE_NOT_FOUND = 0x3140, /*!< Internal encoder error.    */

  aac_enc_init_error_end,

  /* encode errors */
  aac_enc_error_start = 0x4000,
  AAC_ENC_QUANT_ERROR = 0x4020, /*!< Too many bits used in quantization. */
  AAC_ENC_WRITTEN_BITS_ERROR =
      0x4040, /*!< Unexpected number of written bits, differs to
                   calculated number of bits.                      */
  AAC_ENC_PNS_TABLE_ERROR = 0x4060,      /*!< PNS level out of range.      */
  AAC_ENC_GLOBAL_GAIN_TOO_HIGH = 0x4080, /*!< Internal quantizer error. */
  AAC_ENC_BITRES_TOO_LOW = 0x40a0, /*!< Too few  bits in bit reservoir.       */
  AAC_ENC_BITRES_TOO_HIGH = 0x40a1, /*!< Too many bits in bit reservoir.      */
  AAC_ENC_INVALID_CHANNEL_BITRATE = 0x4100,
  AAC_ENC_INVALID_ELEMENTINFO_TYPE = 0x4120, /*!< Internal encoder error. */

  AAC_ENC_WRITE_SCAL_ERROR = 0x41e0, /*!< Error writing scalefacData. */
  AAC_ENC_WRITE_SEC_ERROR = 0x4200,  /*!< Error writing sectionData.  */
  AAC_ENC_WRITE_SPEC_ERROR = 0x4220, /*!< Error writing spectralData. */
  aac_enc_error_end

} AAC_ENCODER_ERROR;
/*-------------------------- defines --------------------------------------*/

#define ANC_DATA_BUFFERSIZE 1024 /* ancBuffer size */

#define MAX_TOTAL_EXT_PAYLOADS ((((8)) * (1)) + (2 + 2))

typedef enum {
  AACENC_BR_MODE_INVALID = -1, /*!< Invalid bitrate mode. */
  AACENC_BR_MODE_CBR = 0,      /*!< Constant bitrate mode.      */
  AACENC_BR_MODE_VBR_1 = 1, /*!< Variable bitrate mode, very low bitrate.    */
  AACENC_BR_MODE_VBR_2 = 2, /*!< Variable bitrate mode, low bitrate.    */
  AACENC_BR_MODE_VBR_3 = 3, /*!< Variable bitrate mode, medium bitrate.    */
  AACENC_BR_MODE_VBR_4 = 4, /*!< Variable bitrate mode, high bitrate.    */
  AACENC_BR_MODE_VBR_5 = 5, /*!< Variable bitrate mode, very high bitrate.    */
  AACENC_BR_MODE_FF = 6,    /*!< Fixed frame mode.       */
  AACENC_BR_MODE_SFR = 7    /*!< Superframe mode.       */

} AACENC_BITRATE_MODE;

#define AACENC_BR_MODE_IS_VBR(brMode) ((brMode >= 1) && (brMode <= 5))

typedef enum {

  CH_ORDER_MPEG =
      0, /*!< MPEG channel ordering (e. g. 5.1: C, L, R, SL, SR, LFE)       */
  CH_ORDER_WAV, /*!< WAV fileformat channel ordering (e. g. 5.1: L, R, C, LFE,
                   SL, SR) */
  CH_ORDER_WG4  /*!< WG4 fileformat channel ordering (e. g. 5.1: L, R, SL, SR, C, LFE) */

} CHANNEL_ORDER;

/*-------------------- structure definitions ------------------------------*/

struct AACENC_CONFIG {
  INT sampleRate;     /* encoder sample rate */
  INT bitRate;        /* encoder bit rate in bits/sec */
  INT ancDataBitRate; /* additional bits consumed by anc data or sbr have to be
                         consiedered while configuration */

  INT nSubFrames; /* number of frames in super frame (not ADTS/LATM subframes !)
                   */
  AUDIO_OBJECT_TYPE audioObjectType; /* Audio Object Type  */

  INT averageBits;                 /* encoder bit rate in bits/superframe */
  AACENC_BITRATE_MODE bitrateMode; /* encoder bitrate mode (CBR/VBR) */
  INT nChannels;                   /* number of channels to process */
  CHANNEL_ORDER channelOrder;      /* input Channel ordering scheme. */
  INT bandWidth;                   /* targeted audio bandwidth in Hz */
  CHANNEL_MODE channelMode;        /* encoder channel mode configuration */
  INT framelength;                 /* used frame size */

  UINT syntaxFlags; /* bitstreams syntax configuration */
  SCHAR epConfig;   /* error protection configuration */

  INT anc_Rate; /* ancillary rate, 0 (disabled), -1 (default) else desired rate
                 */
  UINT maxAncBytesPerAU;
  INT minBitsPerFrame; /* minimum number of bits in AU */
  INT maxBitsPerFrame; /* maximum number of bits in AU */

  INT audioMuxVersion; /* audio mux version in loas/latm transport format */

  UINT sbrRatio; /* sbr sampling rate ratio: dual- or single-rate */

  UCHAR useTns; /* flag: use temporal noise shaping */
  UCHAR usePns; /* flag: use perceptual noise substitution */
  UCHAR useIS;  /* flag: use intensity coding */
  UCHAR useMS;  /* flag: use ms stereo tool */

  UCHAR useRequant; /* flag: use afterburner */

  UINT downscaleFactor;
};

typedef struct {
  UCHAR *pData;              /* pointer to extension payload data */
  UINT dataSize;             /* extension payload data size in bits */
  EXT_PAYLOAD_TYPE dataType; /* extension payload data type */
  INT associatedChElement; /* number of the channel element the data is assigned
                              to */
} AACENC_EXT_PAYLOAD;

typedef struct AAC_ENC *HANDLE_AAC_ENC;

INT FDKaacEnc_CalcBitsPerFrame(const INT bitRate, const INT frameLength,
                               const INT samplingRate);

INT FDKaacEnc_CalcBitrate(const INT bitsPerFrame, const INT frameLength,
                          const INT samplingRate);

INT FDKaacEnc_LimitBitrate(HANDLE_TRANSPORTENC hTpEnc, AUDIO_OBJECT_TYPE aot,
                           INT coreSamplingRate, INT frameLength, INT nChannels,
                           INT nChannelsEff, INT bitRate, INT averageBits,
                           INT *pAverageBitsPerFrame,
                           AACENC_BITRATE_MODE bitrateMode, INT nSubFrames);

INT FDKaacEnc_GetBitReservoirState(const HANDLE_AAC_ENC hAacEncoder);

INT FDKaacEnc_GetVBRBitrate(AACENC_BITRATE_MODE bitrateMode,
                            CHANNEL_MODE channelMode);

AACENC_BITRATE_MODE FDKaacEnc_AdjustVBRBitrateMode(
    AACENC_BITRATE_MODE bitrateMode, INT bitrate, CHANNEL_MODE channelMode);

void FDKaacEnc_AacInitDefaultConfig(AACENC_CONFIG *config);

AAC_ENCODER_ERROR FDKaacEnc_Open(
    HANDLE_AAC_ENC
        *phAacEnc, /* pointer to an encoder handle, initialized on return */
    const INT nElements, /* number of maximal elements in instance to support */
    const INT nChannels, /* number of maximal channels in instance to support */
    const INT nSubFrames); /* support superframing in instance */

AAC_ENCODER_ERROR FDKaacEnc_Initialize(
    HANDLE_AAC_ENC
        hAacEncoder, /* pointer to an encoder handle, initialized on return */
    AACENC_CONFIG *config, /* pre-initialized config struct */
    HANDLE_TRANSPORTENC hTpEnc, ULONG initFlags);

AAC_ENCODER_ERROR FDKaacEnc_EncodeFrame(
    HANDLE_AAC_ENC hAacEnc, /* encoder handle */
    HANDLE_TRANSPORTENC hTpEnc, INT_PCM *inputBuffer,
    const UINT inputBufferBufSize, INT *numOutBytes,
    AACENC_EXT_PAYLOAD extPayload[MAX_TOTAL_EXT_PAYLOADS]);

void FDKaacEnc_Close(HANDLE_AAC_ENC *phAacEnc); /* encoder handle */

#ifdef __cplusplus
}
#endif

#endif /* AACENC_H */
