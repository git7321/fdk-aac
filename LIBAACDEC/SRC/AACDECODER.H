/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description:

*******************************************************************************/

#ifndef AACDECODER_H
#define AACDECODER_H

#include "common_fix.h"
#include "FDK_bitstream.h"
#include "channel.h"
#include "tpdec_lib.h"
#include "FDK_audio.h"
#include "block.h"
#include "genericStds.h"
#include "FDK_qmf_domain.h"
#include "sbrdecoder.h"
#include "aacdec_drc.h"
#include "pcmdmx_lib.h"
#include "FDK_drcDecLib.h"
#include "limiter.h"
#include "FDK_delay.h"

#define TIME_DATA_FLUSH_SIZE (128)
#define TIME_DATA_FLUSH_SIZE_SF (7)

#define AACDEC_MAX_NUM_PREROLL_AU_USAC (3)
#if (AACDEC_MAX_NUM_PREROLL_AU < 3)
#undef AACDEC_MAX_NUM_PREROLL_AU
#define AACDEC_MAX_NUM_PREROLL_AU AACDEC_MAX_NUM_PREROLL_AU_USAC
#endif

typedef struct AAC_DECODER_INSTANCE *HANDLE_AACDECODER;

enum { L = 0, R = 1 };

typedef struct {
  unsigned char *buffer;
  int bufferSize;
  int offset[8];
  int nrElements;
} CAncData;

typedef enum { NOT_DEFINED = -1, MODE_HQ = 0, MODE_LP = 1 } QMF_MODE;

typedef struct {
  int bsDelay;
} SBR_PARAMS;

enum {
  AACDEC_FLUSH_OFF = 0,
  AACDEC_RSV60_CFG_CHANGE_ATSC_FLUSH_ON = 1,
  AACDEC_RSV60_DASH_IPF_ATSC_FLUSH_ON = 2,
  AACDEC_USAC_DASH_IPF_FLUSH_ON = 3
};

enum {
  AACDEC_BUILD_UP_OFF = 0,
  AACDEC_RSV60_BUILD_UP_ON = 1,
  AACDEC_RSV60_BUILD_UP_ON_IN_BAND = 2,
  AACDEC_USAC_BUILD_UP_ON = 3,
  AACDEC_RSV60_BUILD_UP_IDLE = 4,
  AACDEC_RSV60_BUILD_UP_IDLE_IN_BAND = 5
};

#define AACDEC_CROSSFADE_BITMASK_OFF \
  ((UCHAR)0)
#define AACDEC_CROSSFADE_BITMASK_PREROLL \
  ((UCHAR)1 << 1)

typedef struct {
  USAC_EXT_ELEMENT_TYPE usacExtElementType[(3)];
  UINT usacExtElementDefaultLength[(3)];
  UCHAR usacExtElementPayloadFrag[(3)];
} CUsacCoreExtensions;

struct AAC_DECODER_INSTANCE {
  INT aacChannels;
  INT ascChannels[(1 * 1)];
  INT blockNumber;

  INT nrOfLayers;
  INT outputInterleaved;
  INT aacOutDataHeadroom;

  HANDLE_TRANSPORTDEC hInput;

  SamplingRateInfo
      samplingRateInfo[(1 * 1)];

  UCHAR
  frameOK;

  UINT flags[(1 * 1)];
  UINT elFlags[(3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) +
                1)];

  MP4_ELEMENT_ID elements[(3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) +
                           1)];
  UCHAR elTags[(3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) +
                1)];
  UCHAR chMapping[((8) * 2)];

  AUDIO_CHANNEL_TYPE channelType[(8)];
  UCHAR channelIndices[(8)];

  FDK_channelMapDescr mapDescr;
  UCHAR chMapIndex;
  INT sbrDataLen;

  CProgramConfig pce;
  CStreamInfo
      streamInfo;
  CAacDecoderChannelInfo
      *pAacDecoderChannelInfo[(8)];
  CAacDecoderStaticChannelInfo
      *pAacDecoderStaticChannelInfo[(8)];

  FIXP_DBL *workBufferCore1;
  FIXP_DBL *workBufferCore2;
  PCM_DEC *pTimeData2;
  INT timeData2Size;

  CpePersistentData *cpeStaticData[(
      3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) +
      1)];

  CConcealParams concealCommonData;
  CConcealmentMethod concealMethodUser;

  CUsacCoreExtensions usacCoreExt;
  UINT numUsacElements[(1 * 1)];
  UCHAR usacStereoConfigIndex[(3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1)];
  const CSUsacConfig *pUsacConfig[(1 * 1)];
  INT nbDiv;

  UCHAR useLdQmfTimeAlign;

  INT aacChannelsPrev;
  AUDIO_CHANNEL_TYPE channelTypePrev[(8)];
  UCHAR
  channelIndicesPrev[(8)];

  UCHAR
  downscaleFactor;
  UCHAR downscaleFactorInBS;

  HANDLE_SBRDECODER hSbrDecoder;
  UCHAR sbrEnabled;
  UCHAR sbrEnabledPrev;
  UCHAR psPossible;
  SBR_PARAMS sbrParams;

  UCHAR *pDrmBsBuffer;
  USHORT drmBsBufferSize;
  FDK_QMF_DOMAIN
  qmfDomain;

  QMF_MODE qmfModeCurr;
  QMF_MODE qmfModeUser;

  HANDLE_AAC_DRC hDrcInfo;
  INT metadataExpiry;

  void *pMpegSurroundDecoder;
  UCHAR mpsEnableUser;
  UCHAR mpsEnableCurr;
  UCHAR mpsApplicable;
  SCHAR mpsOutputMode;
  INT mpsOutChannelsLast;
  INT mpsFrameSizeLast;
  CAncData ancData;

  HANDLE_PCM_DOWNMIX hPcmUtils;

  TDLimiterPtr hLimiter;
  UCHAR limiterEnableUser;
  UCHAR limiterEnableCurr;

  FIXP_DBL extGain[1];
  UINT extGainDelay;

  HANDLE_DRC_DECODER hUniDrcDecoder;
  UCHAR multibandDrcPresent;
  UCHAR numTimeSlots;
  UINT loudnessInfoSetPosition[3];
  SCHAR defaultTargetLoudness;

  PCM_DEC
  *pTimeDataFlush[((8) * 2)];
  UCHAR flushStatus;
  SCHAR flushCnt;
  UCHAR buildUpStatus;
  SCHAR buildUpCnt;
  UCHAR hasAudioPreRoll;
  UINT prerollAULength[AACDEC_MAX_NUM_PREROLL_AU + 1];
  INT accessUnit;
  UCHAR applyCrossfade;

  FDK_SignalDelay usacResidualDelay;
};

#define AAC_DEBUG_EXTHLP \
  "\
--- AAC-Core ---\n\
    0x00010000 Header data\n\
    0x00020000 CRC data\n\
    0x00040000 Channel info\n\
    0x00080000 Section data\n\
    0x00100000 Scalefactor data\n\
    0x00200000 Pulse data\n\
    0x00400000 Tns data\n\
    0x00800000 Quantized spectrum\n\
    0x01000000 Requantized spectrum\n\
    0x02000000 Time output\n\
    0x04000000 Fatal errors\n\
    0x08000000 Buffer fullness\n\
    0x10000000 Average bitrate\n\
    0x20000000 Synchronization\n\
    0x40000000 Concealment\n\
    0x7FFF0000 all AAC-Core-Info\n\
"

void CAacDecoder_SyncQmfMode(HANDLE_AACDECODER self);

void CAacDecoder_SignalInterruption(HANDLE_AACDECODER self);

#ifdef __cplusplus
extern "C" {
#endif

AAC_DECODER_ERROR CAacDecoder_AncDataInit(CAncData *ancData,
                                          unsigned char *buffer, int size);

AAC_DECODER_ERROR CAacDecoder_AncDataGet(CAncData *ancData, int index,
                                         unsigned char **ptr, int *size);

HANDLE_AACDECODER CAacDecoder_Open(TRANSPORT_TYPE bsFormat);

AAC_DECODER_ERROR CAacDecoder_Init(HANDLE_AACDECODER self,
                                              const CSAudioSpecificConfig *asc,
                                              UCHAR configMode,
                                              UCHAR *configChanged);
AAC_DECODER_ERROR CAacDecoder_DecodeFrame(
    HANDLE_AACDECODER self, const UINT flags, PCM_DEC *pTimeData,
    const INT timeDataSize, const int timeDataChannelOffset);

AAC_DECODER_ERROR CAacDecoder_FreeMem(HANDLE_AACDECODER self,
                                                 const int subStreamIndex);

AAC_DECODER_ERROR CAacDecoder_PrepareCrossFade(
    const PCM_DEC *pTimeData, PCM_DEC **pTimeDataFlush, const INT numChannels,
    const INT frameSize, const INT interleaved);

AAC_DECODER_ERROR CAacDecoder_ApplyCrossFade(
    PCM_DEC *pTimeData, PCM_DEC **pTimeDataFlush, const INT numChannels,
    const INT frameSize, const INT interleaved);

AAC_DECODER_ERROR CAacDecoder_CtrlCFGChange(HANDLE_AACDECODER self,
                                                       UCHAR flushStatus,
                                                       SCHAR flushCnt,
                                                       UCHAR buildUpStatus,
                                                       SCHAR buildUpCnt);

AAC_DECODER_ERROR CAacDecoder_PreRollExtensionPayloadParse(
    HANDLE_AACDECODER self, UINT *numPrerollAU, UINT *prerollAUOffset,
    UINT *prerollAULength);

void CAacDecoder_Close(HANDLE_AACDECODER self);
CStreamInfo * CAacDecoder_GetStreamInfo(HANDLE_AACDECODER self);

#ifdef __cplusplus
}
#endif

#endif /* #ifndef AACDECODER_H */
