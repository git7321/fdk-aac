/*********************** MPEG surround encoder library *************************

   Author(s):   Max Neuendorf

   Description: Encoder Library Interface
                Encoder API

*******************************************************************************/

#ifndef SACENC_LIB_H
#define SACENC_LIB_H

#include "machine_type.h"
#include "FDK_audio.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  SACENC_OK = 0x00000000, /*!< No error happened. All fine. */
  SACENC_INVALID_HANDLE =
      0x00000080, /*!< Handle passed to function call was invalid. */
  SACENC_MEMORY_ERROR = 0x00000800, /*!< Memory allocation failed. */
  SACENC_INIT_ERROR = 0x00008000,   /*!< General initialization error. */
  SACENC_ENCODE_ERROR =
      0x00080000, /*!< The encoding process was interrupted by an unexpected
                     error. */
  SACENC_PARAM_ERROR = 0x00800000,           /*!< Invalid runtime parameter. */
  SACENC_UNSUPPORTED_PARAMETER = 0x00800001, /*!< Parameter not available. */
  SACENC_INVALID_CONFIG = 0x00800002,        /*!< Configuration not provided. */
  SACENC_UNKNOWN_ERROR = 0x08000000          /*!< Unknown error. */

} FDK_SACENC_ERROR;

typedef enum {
  SACENC_INVALID_MODE = 0,
  SACENC_212 = 8,
  SACENC_ESCAPE = 15

} MP4SPACEENC_MODE;

typedef enum {
  SACENC_BANDS_INVALID = 0,
  SACENC_BANDS_4 = 4,
  SACENC_BANDS_5 = 5,
  SACENC_BANDS_7 = 7,
  SACENC_BANDS_9 = 9,
  SACENC_BANDS_12 = 12,
  SACENC_BANDS_15 = 15,
  SACENC_BANDS_23 = 23

} MP4SPACEENC_BANDS_CONFIG;

typedef enum {
  SACENC_QUANTMODE_INVALID = -1,
  SACENC_QUANTMODE_FINE = 0,
  SACENC_QUANTMODE_EBQ1 = 1,
  SACENC_QUANTMODE_EBQ2 = 2,
  SACENC_QUANTMODE_RSVD3 = 3

} MP4SPACEENC_QUANTMODE;

typedef enum {
  SACENC_DMXGAIN_INVALID = -1,
  SACENC_DMXGAIN_0_dB = 0,
  SACENC_DMXGAIN_1_5_dB = 1,
  SACENC_DMXGAIN_3_dB = 2,
  SACENC_DMXGAIN_4_5_dB = 3,
  SACENC_DMXGAIN_6_dB = 4,
  SACENC_DMXGAIN_7_5_dB = 5,
  SACENC_DMXGAIN_9_dB = 6,
  SACENC_DMXGAIN_12_dB = 7

} MP4SPACEENC_DMX_GAIN;

typedef enum {
  SACENC_LOWDELAY, /*!< Configure lowdelay MPEG Surround.
                        - 0: Disable Lowdelay. (default)
                        - 1: Enable Lowdelay.
                        - 2: Enable Lowdelay including keep frame. */

  SACENC_ENC_MODE, /*!< Configure encoder tree mode. See ::MP4SPACEENC_MODE for
                      available values. */

  SACENC_SAMPLERATE, /*!< Configure encoder sampling rate. */

  SACENC_FRAME_TIME_SLOTS, /*!< Configure number of slots per spatial frame. */

  SACENC_PARAM_BANDS, /*!< Configure number of parameter bands. See
                         ::MP4SPACEENC_BANDS_CONFIG for available values. */

  SACENC_TIME_DOM_DMX, /*!< Configure time domain downmix.
                            - 0: No time domain downmix. (default)
                            - 1: Static time domain downmix.
                            - 2: Enhanced time domain downmix, stereo to mono
                          only. */

  SACENC_DMX_GAIN, /*!< Configure downmix gain. See ::MP4SPACEENC_DMX_GAIN for
                      available values. */

  SACENC_COARSE_QUANT, /*!< Use coarse parameter quantization.
                            - 0: No (default)
                            - 1: Yes */

  SACENC_QUANT_MODE, /*!< Configure quanitzation mode. See
                        ::MP4SPACEENC_QUANTMODE for available values. */

  SACENC_TIME_ALIGNMENT, /*!< Configure time alignment in samples. */

  SACENC_INDEPENDENCY_COUNT, /*!< Configure the independency count. (count == 0
                                means independencyFlag == 1) */

  SACENC_INDEPENDENCY_FACTOR, /*!< How often should we set the independency flag
                               */

  SACENC_NONE /*!< ------ */

} SPACEENC_PARAM;

typedef struct {
  INT nSscSizeBits; /*!< Number of valid bits in pSsc buffer. */
  UCHAR *pSsc;      /*!< SpatialSpecificConfig buffer in binary format. */

} MPEG4SPACEENC_SSCBUF;

typedef struct {
  INT nSampleRate;         /*!< Configured sampling rate.*/
  INT nSamplesFrame;       /*!< Frame length in samples. */
  INT nTotalInputChannels; /*!< Number of expected audio input channels. */
  INT nDmxDelay;           /*!< Delay of the downmixed signal. */
  INT nCodecDelay;         /*!< Delay of the whole en-/decoded signal, including
                              core-coder delay. */
  INT nDecoderDelay;       /*!< Delay added by the MP4SPACE decoder. */
  INT nPayloadDelay;       /*!< Delay of the payload. */
  INT nDiscardOutFrames; /*!< Number of dmx frames to discard for alignment with
                            bitstream. */

  MPEG4SPACEENC_SSCBUF
  *pSscBuf; /*!< Pointer to Spatial Specific Config structure. */

} MP4SPACEENC_INFO;

typedef struct MP4SPACE_ENCODER *HANDLE_MP4SPACE_ENCODER;

typedef struct {
  INT nInputSamples; /*!< Number of valid input audio samples (multiple of input
                        channels). */
  UINT inputBufferSizePerChannel; /*!< Size of input buffer (input audio
                                     samples) per channel. */
  UINT isInputInterleaved; /*!< Indicates if input audio samples are represented
                              in blocks or interleaved:
                                - 0 : in blocks.
                                - 1 : interleaved. */

} SACENC_InArgs;

typedef struct {
  INT nOutputBits;    /*!< Number of valid payload bits generated during
                         FDK_sacenc_encode(). */
  INT nOutputSamples; /*!< Number of valid output audio samples generated during
                         FDK_sacenc_encode(). */
  UINT nSamplesConsumed; /*!< Number of input audio samples consumed in
                            FDK_sacenc_encode(). */

} SACENC_OutArgs;

FDK_SACENC_ERROR FDK_sacenc_open(HANDLE_MP4SPACE_ENCODER *phMp4SpaceEnc);

FDK_SACENC_ERROR FDK_sacenc_init(HANDLE_MP4SPACE_ENCODER hMp4SpaceEnc,
                                 const INT dmxDelay);

FDK_SACENC_ERROR FDK_sacenc_close(HANDLE_MP4SPACE_ENCODER *phMp4SpaceEnc);

FDK_SACENC_ERROR FDK_sacenc_encode(const HANDLE_MP4SPACE_ENCODER hMp4SpaceEnc,
                                   const FDK_bufDescr *inBufDesc,
                                   const FDK_bufDescr *outBufDesc,
                                   const SACENC_InArgs *inargs,
                                   SACENC_OutArgs *outargs);

FDK_SACENC_ERROR FDK_sacenc_getInfo(const HANDLE_MP4SPACE_ENCODER hMp4SpaceEnc,
                                    MP4SPACEENC_INFO *const pInfo);

FDK_SACENC_ERROR FDK_sacenc_setParam(HANDLE_MP4SPACE_ENCODER hMp4SpaceEnc,
                                     const SPACEENC_PARAM param,
                                     const UINT value);

FDK_SACENC_ERROR FDK_sacenc_getLibInfo(LIB_INFO *info);

#ifdef __cplusplus
}
#endif

#endif /* SACENC_LIB_H */
