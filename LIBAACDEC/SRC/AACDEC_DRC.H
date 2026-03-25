/**************************** AAC decoder library ******************************

   Author(s):   Christian Griebel

   Description: Dynamic range control (DRC) decoder tool for AAC

*******************************************************************************/

#ifndef AACDEC_DRC_H
#define AACDEC_DRC_H

#include "tp_data.h"

#include "aacdec_drc_types.h"
#include "channel.h"
#include "FDK_bitstream.h"

#define AACDEC_DRC_GAIN_SCALING (11)
#define AACDEC_DRC_GAIN_INIT_VALUE \
  (FL2FXCONST_DBL(                 \
      1.0f / (1 << AACDEC_DRC_GAIN_SCALING)))

#define AACDEC_DRC_DFLT_EXPIRY_FRAMES \
  (0)

typedef enum {
  DRC_CUT_SCALE = 0,
  DRC_BOOST_SCALE,
  TARGET_REF_LEVEL,
  DRC_BS_DELAY,
  DRC_DATA_EXPIRY_FRAME,
  APPLY_HEAVY_COMPRESSION,
  DEFAULT_PRESENTATION_MODE,
  ENCODER_TARGET_LEVEL,
  MAX_OUTPUT_CHANNELS
} AACDEC_DRC_PARAM;

void aacDecoder_drcDisable(HANDLE_AAC_DRC self);

void aacDecoder_drcReset(HANDLE_AAC_DRC self);

void aacDecoder_drcInit(HANDLE_AAC_DRC self);

void aacDecoder_drcInitChannelData(CDrcChannelData *pDrcChannel);

AAC_DECODER_ERROR aacDecoder_drcSetParam(HANDLE_AAC_DRC self,
                                         AACDEC_DRC_PARAM param, INT value);

int aacDecoder_drcMarkPayload(HANDLE_AAC_DRC self, HANDLE_FDK_BITSTREAM hBs,
                              AACDEC_DRC_PAYLOAD_TYPE type);

int aacDecoder_drcProlog(
    HANDLE_AAC_DRC self, HANDLE_FDK_BITSTREAM hBs,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[],
    UCHAR pceInstanceTag, UCHAR channelMapping[], int validChannels);

void aacDecoder_drcApply(HANDLE_AAC_DRC self, void *pSbrDec,
                         CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                         CDrcChannelData *pDrcDat, FIXP_DBL *extGain, int ch,
                         int aacFrameSize, int bSbrPresent);

int aacDecoder_drcEpilog(
    HANDLE_AAC_DRC self, HANDLE_FDK_BITSTREAM hBs,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[],
    UCHAR pceInstanceTag, UCHAR channelMapping[], int validChannels);

void aacDecoder_drcGetInfo(HANDLE_AAC_DRC self, SCHAR *pPresMode,
                           SCHAR *pProgRefLevel);

INT applyDrcLevelNormalization(HANDLE_AAC_DRC hDrcInfo, PCM_DEC *samplesIn,
                               FIXP_DBL *pGain, FIXP_DBL *pGainPerSample,
                               const INT gain_scale, const UINT gain_delay,
                               const UINT nSamples, const UINT channels,
                               const UINT stride, const UINT limiterEnabled);

#endif /* AACDEC_DRC_H */
