
/**************************** AAC encoder library ******************************

   Author(s):   M. Lohwasser

   Description: FDK HE-AAC Encoder interface library functions

*******************************************************************************/

#include "aacenc_lib.h"
#include "FDK_audio.h"
#include "aacenc.h"
#include "aacEnc_ram.h"
#include "FDK_core.h"

#define AACENCODER_LIB_VL0 4
#define AACENCODER_LIB_VL1 0
#define AACENCODER_LIB_VL2 1
#define AACENCODER_LIB_TITLE "AAC Encoder"
#define AACENCODER_LIB_BUILD_DATE ""
#define AACENCODER_LIB_BUILD_TIME ""

#include "pcm_utils.h"
#include "sbr_encoder.h"
#include "../src/sbrenc_ram.h"
#include "channel_map.h"

#include "psy_const.h"
#include "bitenc.h"

#include "tpenc_lib.h"

#include "metadata_main.h"
#include "mps_main.h"
#include "sacenc_lib.h"

#define SBL(fl) \
  (fl / \
   8)
#define BSLA(fl) \
  (4 * SBL(fl) + SBL(fl) / 2)
#define DELAY_AAC(fl) (fl + BSLA(fl))
#define DELAY_AACLD(fl) (fl)
#define DELAY_AACELD(fl) \
  ((fl) / 2)

#define MAX_DS_DELAY (100)
#define INPUTBUFFER_SIZE \
  (2 * (1024) + MAX_DS_DELAY + 1537)

#define DEFAULT_HEADER_PERIOD_REPETITION_RATE \
  10

enum {
  ENC_MODE_FLAG_AAC = 0x0001,
  ENC_MODE_FLAG_SBR = 0x0002,
  ENC_MODE_FLAG_PS = 0x0004,
  ENC_MODE_FLAG_SAC = 0x0008,
  ENC_MODE_FLAG_META = 0x0010
};

typedef struct {
  AUDIO_OBJECT_TYPE userAOT;
  UINT userSamplerate;
  UINT nChannels;
  CHANNEL_MODE userChannelMode;
  UINT userBitrate;
  UINT userBitrateMode;
  UINT userBandwidth;
  UINT userAfterburner;
  UINT userFramelength;
  UINT userAncDataRate;
  UINT userPeakBitrate;

  UCHAR userTns;
  UCHAR userPns;
  UCHAR userIntensity;

  TRANSPORT_TYPE userTpType;
  UCHAR userTpSignaling;
  UCHAR userTpNsubFrames;
  UCHAR userTpAmxv;
  UCHAR userTpProtection;
  UCHAR userTpHeaderPeriod;

  UCHAR userErTools;
  UINT userPceAdditions;

  UCHAR userMetaDataMode;

  UCHAR userSbrEnabled;
  UINT userSbrRatio;

  UINT userDownscaleFactor;

} USER_PARAM;

typedef struct {
  UCHAR data[(1)][(8)][MAX_PAYLOAD_SIZE];
  UINT dataSize[(1)][(8)];
} SBRENC_EXT_PAYLOAD;

typedef struct AACENC_CONFIG *HANDLE_AACENC_CONFIG;

struct AACENCODER {
  USER_PARAM extParam;
  CODER_CONFIG coderConfig;

  AACENC_CONFIG aacConfig;
  HANDLE_AAC_ENC hAacEnc;
  HANDLE_SBR_ENCODER hEnvEnc;
  SBRENC_EXT_PAYLOAD *pSbrPayload;

  HANDLE_FDK_METADATA_ENCODER hMetadataEnc;
  INT metaDataAllowed;

  HANDLE_MPS_ENCODER hMpsEnc;
  HANDLE_TRANSPORTENC hTpEnc;

  INT_PCM
  *inputBuffer;
  UCHAR *outBuffer;

  INT inputBufferSize;
  INT inputBufferSizePerChannel;
  INT outBufferInBytes;

  INT inputBufferOffset;

  INT nSamplesToRead;
  INT nSamplesRead;
  INT nZerosAppended;
  INT nDelay;
  INT nDelayCore;

  AACENC_EXT_PAYLOAD extPayload[MAX_TOTAL_EXT_PAYLOADS];
  ULONG InitFlags;

  INT nMaxAacElements;
  INT nMaxAacChannels;
  INT nMaxSbrElements;
  INT nMaxSbrChannels;

  UINT encoder_modis;
  UINT CAPF_tpEnc;
};

typedef struct {
  ULONG nChannels;
  ULONG samplingRate;
  ULONG bitrateRange;

  UCHAR sbrMode;
  CHANNEL_MODE chMode;

} ELD_SBR_CONFIGURATOR;

static const ELD_SBR_CONFIGURATOR eldSbrAutoConfigTab[] = {
    {1, 48000, 0, 2, MODE_1},      {1, 48000, 64000, 0, MODE_1},

    {1, 44100, 0, 2, MODE_1},      {1, 44100, 64000, 0, MODE_1},

    {1, 32000, 0, 2, MODE_1},      {1, 32000, 28000, 1, MODE_1},
    {1, 32000, 56000, 0, MODE_1},

    {1, 24000, 0, 1, MODE_1},      {1, 24000, 40000, 0, MODE_1},

    {1, 16000, 0, 1, MODE_1},      {1, 16000, 28000, 0, MODE_1},

    {1, 15999, 0, 0, MODE_1},

    {2, 48000, 0, 2, MODE_2},      {2, 48000, 44000, 2, MODE_2},
    {2, 48000, 128000, 0, MODE_2},

    {2, 44100, 0, 2, MODE_2},      {2, 44100, 44000, 2, MODE_2},
    {2, 44100, 128000, 0, MODE_2},

    {2, 32000, 0, 2, MODE_2},      {2, 32000, 32000, 2, MODE_2},
    {2, 32000, 68000, 1, MODE_2},  {2, 32000, 96000, 0, MODE_2},

    {2, 24000, 0, 1, MODE_2},      {2, 24000, 48000, 1, MODE_2},
    {2, 24000, 80000, 0, MODE_2},

    {2, 16000, 0, 1, MODE_2},      {2, 16000, 32000, 1, MODE_2},
    {2, 16000, 64000, 0, MODE_2},

    {2, 15999, 0, 0, MODE_2}

};

static const ELD_SBR_CONFIGURATOR *eldSbrConfigurator(const ULONG nChannels,
                                                      const ULONG samplingRate,
                                                      const ULONG bitrate) {
  int i;
  const ELD_SBR_CONFIGURATOR *pSetup = NULL;

  for (i = 0;
       i < (int)(sizeof(eldSbrAutoConfigTab) / sizeof(ELD_SBR_CONFIGURATOR));
       i++) {
    if ((nChannels == eldSbrAutoConfigTab[i].nChannels) &&
        (samplingRate <= eldSbrAutoConfigTab[i].samplingRate) &&
        (bitrate >= eldSbrAutoConfigTab[i].bitrateRange)) {
      pSetup = &eldSbrAutoConfigTab[i];
    }
  }

  return pSetup;
}

static inline INT isSbrActive(const HANDLE_AACENC_CONFIG hAacConfig) {
  INT sbrUsed = 0;

  if ((hAacConfig->audioObjectType == AOT_SBR) ||
      (hAacConfig->audioObjectType == AOT_PS) ||
      (hAacConfig->audioObjectType == AOT_MP2_SBR)) {
    sbrUsed = 1;
  }
  if (hAacConfig->audioObjectType == AOT_ER_AAC_ELD &&
      (hAacConfig->syntaxFlags & AC_SBR_PRESENT)) {
    sbrUsed = 1;
  }

  return (sbrUsed);
}

static inline INT isPsActive(const AUDIO_OBJECT_TYPE audioObjectType) {
  INT psUsed = 0;

  if (audioObjectType == AOT_PS) {
    psUsed = 1;
  }

  return (psUsed);
}

static CHANNEL_MODE GetCoreChannelMode(
    const CHANNEL_MODE channelMode, const AUDIO_OBJECT_TYPE audioObjectType) {
  CHANNEL_MODE mappedChannelMode = channelMode;
  if ((isPsActive(audioObjectType) && (channelMode == MODE_2)) ||
      (channelMode == MODE_212)) {
    mappedChannelMode = MODE_1;
  }
  return mappedChannelMode;
}

static SBR_PS_SIGNALING getSbrSignalingMode(
    const AUDIO_OBJECT_TYPE audioObjectType, const TRANSPORT_TYPE transportType,
    const UCHAR transportSignaling, const UINT sbrRatio)

{
  SBR_PS_SIGNALING sbrSignaling;

  if (transportType == TT_UNKNOWN || sbrRatio == 0) {
    sbrSignaling = SIG_UNKNOWN;
    return sbrSignaling;
  } else {
    sbrSignaling =
        SIG_EXPLICIT_HIERARCHICAL;
  }

  if ((audioObjectType == AOT_AAC_LC) || (audioObjectType == AOT_SBR) ||
      (audioObjectType == AOT_PS) || (audioObjectType == AOT_MP2_AAC_LC) ||
      (audioObjectType == AOT_MP2_SBR)) {
    switch (transportType) {
      case TT_MP4_ADIF:
      case TT_MP4_ADTS:
        sbrSignaling = SIG_IMPLICIT;
        break;

      case TT_MP4_RAW:
      case TT_MP4_LATM_MCP1:
      case TT_MP4_LATM_MCP0:
      case TT_MP4_LOAS:
      default:
        if (transportSignaling == 0xFF) {
          sbrSignaling = SIG_EXPLICIT_HIERARCHICAL;
        } else {
          sbrSignaling = (SBR_PS_SIGNALING)transportSignaling;
        }
        break;
    }
  }

  return sbrSignaling;
}

static inline INT getAssociatedChElement(SBR_ELEMENT_INFO *elInfoSbr,
                                         CHANNEL_MAPPING *channelMapping) {
  ELEMENT_INFO *elInfo = channelMapping->elInfo;
  INT nElements = channelMapping->nElements;
  INT associatedChElement = -1;
  int i;

  for (i = 0; i < nElements; i++) {
    if (elInfoSbr->elType == elInfo[i].elType &&
        elInfoSbr->instanceTag == elInfo[i].instanceTag) {
      associatedChElement = i;
      break;
    }
  }

  return associatedChElement;
}

H_ALLOC_MEM(_AacEncoder, AACENCODER)
C_ALLOC_MEM(_AacEncoder, struct AACENCODER, 1)

static void FDKaacEnc_MapConfig(CODER_CONFIG *const cc,
                                const USER_PARAM *const extCfg,
                                const SBR_PS_SIGNALING sbrSignaling,
                                const HANDLE_AACENC_CONFIG hAacConfig) {
  AUDIO_OBJECT_TYPE transport_AOT = AOT_NULL_OBJECT;
  FDKmemclear(cc, sizeof(CODER_CONFIG));

  cc->flags = 0;

  cc->samplesPerFrame = hAacConfig->framelength;
  cc->samplingRate = hAacConfig->sampleRate;
  cc->extSamplingRate = extCfg->userSamplerate;

  switch (hAacConfig->audioObjectType) {
    case AOT_MP2_AAC_LC:
      transport_AOT = AOT_AAC_LC;
      break;
    case AOT_MP2_SBR:
      transport_AOT = AOT_SBR;
      cc->flags |= CC_SBR;
      break;
    default:
      transport_AOT = hAacConfig->audioObjectType;
  }

  if (hAacConfig->audioObjectType == AOT_ER_AAC_ELD) {
    cc->flags |= (hAacConfig->syntaxFlags & AC_SBR_PRESENT) ? CC_SBR : 0;
    cc->flags |= (hAacConfig->syntaxFlags & AC_LD_MPS) ? CC_SAC : 0;
  }

  if ((transport_AOT == AOT_SBR) || (transport_AOT == AOT_PS)) {
    cc->aot = AOT_AAC_LC;
  } else {
    cc->aot = transport_AOT;
  }

  if (sbrSignaling == SIG_IMPLICIT) {
    cc->extAOT = AOT_NULL_OBJECT;
  } else {
    if ((sbrSignaling == SIG_EXPLICIT_BW_COMPATIBLE) &&
        ((transport_AOT == AOT_SBR) || (transport_AOT == AOT_PS))) {
      cc->extAOT = AOT_SBR;
    } else {
      cc->extAOT = transport_AOT;
    }
  }

  if ((transport_AOT == AOT_SBR) || (transport_AOT == AOT_PS)) {
    cc->sbrPresent = 1;
    if (transport_AOT == AOT_PS) {
      cc->psPresent = 1;
    }
  }
  cc->sbrSignaling = sbrSignaling;

  if (hAacConfig->downscaleFactor > 1) {
    cc->downscaleSamplingRate = cc->samplingRate;
    cc->samplingRate *= hAacConfig->downscaleFactor;
    cc->extSamplingRate *= hAacConfig->downscaleFactor;
  }

  cc->bitRate = hAacConfig->bitRate;
  cc->noChannels = hAacConfig->nChannels;
  cc->flags |= CC_IS_BASELAYER;
  cc->channelMode = hAacConfig->channelMode;

  cc->nSubFrames = (hAacConfig->nSubFrames > 1 && extCfg->userTpNsubFrames == 1)
                       ? hAacConfig->nSubFrames
                       : extCfg->userTpNsubFrames;

  cc->flags |= (extCfg->userTpProtection) ? CC_PROTECTION : 0;

  if (extCfg->userTpHeaderPeriod != 0xFF) {
    cc->headerPeriod = extCfg->userTpHeaderPeriod;
  } else {
    switch (extCfg->userTpType) {
      case TT_MP4_ADTS:
      case TT_MP4_LOAS:
      case TT_MP4_LATM_MCP1:
        cc->headerPeriod = DEFAULT_HEADER_PERIOD_REPETITION_RATE;
        break;
      default:
        cc->headerPeriod = 0;
    }
  }

  switch (hAacConfig->audioObjectType) {
    case AOT_MP2_AAC_LC:
    case AOT_MP2_SBR:
      cc->flags &= ~CC_MPEG_ID;
      cc->extAOT = AOT_NULL_OBJECT;
      break;
    default:
      cc->flags |= CC_MPEG_ID;
  }

  cc->flags |= (hAacConfig->syntaxFlags & AC_ER_VCB11) ? CC_VCB11 : 0;
  cc->flags |= (hAacConfig->syntaxFlags & AC_ER_HCR) ? CC_HCR : 0;
  cc->flags |= (hAacConfig->syntaxFlags & AC_ER_RVLC) ? CC_RVLC : 0;

  if ((extCfg->userPceAdditions & 0x1) && (hAacConfig->epConfig == -1) &&
      ((cc->channelMode == MODE_1_2_2) || (cc->channelMode == MODE_1_2_2_1))) {
    cc->matrixMixdownA = ((extCfg->userPceAdditions >> 1) & 0x3) + 1;
    cc->flags |= (extCfg->userPceAdditions >> 3) & 0x1 ? CC_PSEUDO_SURROUND : 0;
  } else {
    cc->matrixMixdownA = 0;
  }

  cc->channelConfigZero = 0;
}

static AACENC_ERROR validateBufDesc(const AACENC_BufDesc *pBufDesc) {
  AACENC_ERROR err = AACENC_OK;

  if (pBufDesc != NULL) {
    int i;
    if ((pBufDesc->bufferIdentifiers == NULL) || (pBufDesc->bufSizes == NULL) ||
        (pBufDesc->bufElSizes == NULL) || (pBufDesc->bufs == NULL)) {
      err = AACENC_UNSUPPORTED_PARAMETER;
      goto bail;
    }
    for (i = 0; i < pBufDesc->numBufs; i++) {
      if (pBufDesc->bufs[i] == NULL) {
        err = AACENC_UNSUPPORTED_PARAMETER;
        goto bail;
      }
    }
  } else {
    err = AACENC_INVALID_HANDLE;
  }
bail:
  return err;
}

static INT getBufDescIdx(const AACENC_BufDesc *pBufDesc,
                         const AACENC_BufferIdentifier identifier) {
  INT i, idx = -1;

  if (pBufDesc != NULL) {
    for (i = 0; i < pBufDesc->numBufs; i++) {
      if ((AACENC_BufferIdentifier)pBufDesc->bufferIdentifiers[i] ==
          identifier) {
        idx = i;
        break;
      }
    }
  }
  return idx;
}

AAC_ENCODER_ERROR aacEncDefaultConfig(HANDLE_AACENC_CONFIG hAacConfig,
                                      USER_PARAM *config) {
  FDKaacEnc_AacInitDefaultConfig(hAacConfig);
  FDKmemclear(config, sizeof(USER_PARAM));

  config->nChannels = hAacConfig->nChannels;
  config->userAOT = hAacConfig->audioObjectType = AOT_AAC_LC;
  config->userSamplerate = hAacConfig->sampleRate;
  config->userChannelMode = hAacConfig->channelMode;
  config->userBitrate = hAacConfig->bitRate;
  config->userBitrateMode = hAacConfig->bitrateMode;
  config->userPeakBitrate = (UINT)-1;
  config->userBandwidth = hAacConfig->bandWidth;
  config->userTns = hAacConfig->useTns;
  config->userPns = hAacConfig->usePns;
  config->userIntensity = hAacConfig->useIS;
  config->userAfterburner = hAacConfig->useRequant;
  config->userFramelength = (UINT)-1;

  config->userDownscaleFactor = 1;

  config->userTpType = TT_UNKNOWN;
  config->userTpAmxv = 0;
  config->userTpSignaling = 0xFF;
  config->userTpNsubFrames = 1;
  config->userTpProtection = 0;
  config->userTpHeaderPeriod = 0xFF;
  config->userPceAdditions = 0;
  config->userMetaDataMode = 0;

  config->userAncDataRate = 0;

  config->userSbrRatio = 0;

  config->userSbrEnabled = (UCHAR)-1;

  return AAC_ENC_OK;
}

static void aacEncDistributeSbrBits(CHANNEL_MAPPING *channelMapping,
                                    SBR_ELEMENT_INFO *sbrElInfo, INT bitRate) {
  INT codebits = bitRate;
  int el;

  for (el = 0; el < channelMapping->nElements; el++) {
    sbrElInfo[el].ChannelIndex[0] = channelMapping->elInfo[el].ChannelIndex[0];
    sbrElInfo[el].ChannelIndex[1] = channelMapping->elInfo[el].ChannelIndex[1];
    sbrElInfo[el].elType = channelMapping->elInfo[el].elType;
    sbrElInfo[el].bitRate =
        fMultIfloor(channelMapping->elInfo[el].relativeBits, bitRate);
    sbrElInfo[el].instanceTag = channelMapping->elInfo[el].instanceTag;
    sbrElInfo[el].nChannelsInEl = channelMapping->elInfo[el].nChannelsInEl;
    sbrElInfo[el].fParametricStereo = 0;
    sbrElInfo[el].fDualMono = 0;

    codebits -= sbrElInfo[el].bitRate;
  }
  sbrElInfo[0].bitRate += codebits;
}

static INT aacEncoder_LimitBitrate(const HANDLE_TRANSPORTENC hTpEnc,
                                   const INT samplingRate,
                                   const INT frameLength, const INT nChannels,
                                   const CHANNEL_MODE channelMode, INT bitRate,
                                   const INT nSubFrames, const INT sbrActive,
                                   const INT sbrDownSampleRate,
                                   const UINT syntaxFlags,
                                   const AUDIO_OBJECT_TYPE aot) {
  INT coreSamplingRate;
  CHANNEL_MAPPING cm;

  FDKaacEnc_InitChannelMapping(channelMode, CH_ORDER_MPEG, &cm);

  if (sbrActive) {
    coreSamplingRate =
        samplingRate >>
        (sbrEncoder_IsSingleRatePossible(aot) ? (sbrDownSampleRate - 1) : 1);
  } else {
    coreSamplingRate = samplingRate;
  }

  bitRate = FDKaacEnc_LimitBitrate(hTpEnc, aot, coreSamplingRate, frameLength,
                                   nChannels, cm.nChannelsEff, bitRate, -1,
                                   NULL, AACENC_BR_MODE_INVALID, nSubFrames);

  if (sbrActive) {
    int numIterations = 0;
    INT initialBitrate, adjustedBitrate;
    adjustedBitrate = bitRate;

    do {
      int e;
      SBR_ELEMENT_INFO sbrElInfo[((8))];
      FDK_ASSERT(cm.nElements <= ((8)));

      initialBitrate = adjustedBitrate;

      aacEncDistributeSbrBits(&cm, sbrElInfo, initialBitrate);

      for (e = 0; e < cm.nElements; e++) {
        INT sbrElementBitRateIn, sbrBitRateOut;

        if (cm.elInfo[e].elType != ID_SCE && cm.elInfo[e].elType != ID_CPE) {
          continue;
        }
        sbrElementBitRateIn = sbrElInfo[e].bitRate;

        sbrBitRateOut = sbrEncoder_LimitBitRate(sbrElementBitRateIn,
                                                cm.elInfo[e].nChannelsInEl,
                                                coreSamplingRate, aot);

        if (sbrBitRateOut == 0) {
          return 0;
        }

        if (sbrElementBitRateIn != sbrBitRateOut) {
          if (sbrElementBitRateIn < sbrBitRateOut) {
            adjustedBitrate = fMax(initialBitrate,
                                   (INT)fDivNorm((FIXP_DBL)(sbrBitRateOut + 8),
                                                 cm.elInfo[e].relativeBits));
            break;
          }

          if (sbrElementBitRateIn > sbrBitRateOut) {
            adjustedBitrate = fMin(initialBitrate,
                                   (INT)fDivNorm((FIXP_DBL)(sbrBitRateOut - 8),
                                                 cm.elInfo[e].relativeBits));
            break;
          }
        }
      }

      numIterations++;

    } while ((initialBitrate != adjustedBitrate) &&
             (numIterations <= cm.nElements));

    bitRate = (initialBitrate == adjustedBitrate) ? adjustedBitrate : 0;
  }

  if ((aot == AOT_ER_AAC_ELD) && (syntaxFlags & AC_LD_MPS) &&
      (channelMode == MODE_1)) {
    bitRate = FDK_MpegsEnc_GetClosestBitRate(
        aot, MODE_212, samplingRate, (sbrActive) ? sbrDownSampleRate : 0,
        bitRate);
  }

  return bitRate;
}

static INT FDKaacEnc_GetCBRBitrate(const HANDLE_AACENC_CONFIG hAacConfig,
                                   const INT userSbrRatio) {
  INT bitrate = FDKaacEnc_GetChannelModeConfiguration(hAacConfig->channelMode)
                    ->nChannelsEff *
                hAacConfig->sampleRate;

  if (isPsActive(hAacConfig->audioObjectType)) {
    bitrate = 1 * bitrate;
  } else if (isSbrActive(hAacConfig)) {
    if ((userSbrRatio == 2) ||
        ((userSbrRatio == 0) &&
         (hAacConfig->audioObjectType != AOT_ER_AAC_ELD))) {
      bitrate = (bitrate + (bitrate >> 2)) >> 1;
    }
    if ((userSbrRatio == 1) ||
        ((userSbrRatio == 0) &&
         (hAacConfig->audioObjectType == AOT_ER_AAC_ELD))) {
      bitrate = (bitrate + (bitrate >> 3));
    }
  } else {
    bitrate = bitrate + (bitrate >> 1);
  }

  return bitrate;
}

static AACENC_ERROR FDKaacEnc_AdjustEncSettings(HANDLE_AACENCODER hAacEncoder,
                                                USER_PARAM *config) {
  AACENC_ERROR err = AACENC_OK;

  HANDLE_AACENC_CONFIG hAacConfig = &hAacEncoder->aacConfig;

  hAacConfig->sampleRate = config->userSamplerate;
  if (config->userDownscaleFactor > 1) {
    hAacConfig->useTns = 0;
    hAacConfig->usePns = 0;
    hAacConfig->useIS = 0;
  } else {
    hAacConfig->useTns = config->userTns;
    hAacConfig->usePns = config->userPns;
    hAacConfig->useIS = config->userIntensity;
  }

  hAacConfig->audioObjectType = config->userAOT;
  hAacConfig->channelMode =
      GetCoreChannelMode(config->userChannelMode, hAacConfig->audioObjectType);
  hAacConfig->nChannels =
      FDKaacEnc_GetChannelModeConfiguration(hAacConfig->channelMode)->nChannels;
  hAacConfig->bitrateMode = (AACENC_BITRATE_MODE)config->userBitrateMode;
  hAacConfig->bandWidth = config->userBandwidth;
  hAacConfig->useRequant = config->userAfterburner;

  hAacConfig->anc_Rate = config->userAncDataRate;
  hAacConfig->syntaxFlags = 0;
  hAacConfig->epConfig = -1;

  if (hAacConfig->audioObjectType != AOT_ER_AAC_ELD &&
      config->userDownscaleFactor > 1) {
    return AACENC_INVALID_CONFIG;
  }
  if (config->userDownscaleFactor > 1 && config->userSbrEnabled == 1) {
    return AACENC_INVALID_CONFIG;
  }
  if (config->userDownscaleFactor > 1 && config->userChannelMode == 128) {
    return AACENC_INVALID_CONFIG;
  }

  if (config->userTpType == TT_MP4_LATM_MCP1 ||
      config->userTpType == TT_MP4_LATM_MCP0 ||
      config->userTpType == TT_MP4_LOAS) {
    hAacConfig->audioMuxVersion = config->userTpAmxv;
  } else {
    hAacConfig->audioMuxVersion = -1;
  }

  switch (config->userAOT) {
    case AOT_MP2_AAC_LC:
    case AOT_MP2_SBR:
      hAacConfig->usePns = 0;
      FDK_FALLTHROUGH;
    case AOT_AAC_LC:
    case AOT_SBR:
    case AOT_PS:
      config->userTpType =
          (config->userTpType != TT_UNKNOWN) ? config->userTpType : TT_MP4_ADTS;
      hAacConfig->framelength = (config->userFramelength != (UINT)-1)
                                    ? config->userFramelength
                                    : 1024;
      if (hAacConfig->framelength != 1024 && hAacConfig->framelength != 960) {
        return AACENC_INVALID_CONFIG;
      }
      break;
    case AOT_ER_AAC_LD:
      hAacConfig->epConfig = 0;
      hAacConfig->syntaxFlags |= AC_ER | AC_LD;
      hAacConfig->syntaxFlags |=
          ((config->userErTools & 0x1) ? AC_ER_VCB11 : 0);
      hAacConfig->syntaxFlags |= ((config->userErTools & 0x2) ? AC_ER_HCR : 0);
      hAacConfig->syntaxFlags |= ((config->userErTools & 0x4) ? AC_ER_RVLC : 0);
      config->userTpType =
          (config->userTpType != TT_UNKNOWN) ? config->userTpType : TT_MP4_LOAS;
      hAacConfig->framelength =
          (config->userFramelength != (UINT)-1) ? config->userFramelength : 512;
      if (hAacConfig->framelength != 512 && hAacConfig->framelength != 480) {
        return AACENC_INVALID_CONFIG;
      }
      break;
    case AOT_ER_AAC_ELD:
      hAacConfig->epConfig = 0;
      hAacConfig->syntaxFlags |= AC_ER | AC_ELD;
      hAacConfig->syntaxFlags |=
          ((config->userErTools & 0x1) ? AC_ER_VCB11 : 0);
      hAacConfig->syntaxFlags |= ((config->userErTools & 0x2) ? AC_ER_HCR : 0);
      hAacConfig->syntaxFlags |= ((config->userErTools & 0x4) ? AC_ER_RVLC : 0);
      hAacConfig->syntaxFlags |=
          ((config->userSbrEnabled == 1) ? AC_SBR_PRESENT : 0);
      hAacConfig->syntaxFlags |=
          ((config->userChannelMode == MODE_212) ? AC_LD_MPS : 0);
      config->userTpType =
          (config->userTpType != TT_UNKNOWN) ? config->userTpType : TT_MP4_LOAS;
      hAacConfig->framelength =
          (config->userFramelength != (UINT)-1) ? config->userFramelength : 512;

      hAacConfig->downscaleFactor = config->userDownscaleFactor;

      switch (config->userDownscaleFactor) {
        case 1:
          break;
        case 2:
        case 4:
          hAacConfig->syntaxFlags |= AC_ELD_DOWNSCALE;
          break;
        default:
          return AACENC_INVALID_CONFIG;
      }

      if (hAacConfig->framelength != 512 && hAacConfig->framelength != 480 &&
          hAacConfig->framelength != 256 && hAacConfig->framelength != 240 &&
          hAacConfig->framelength != 128 && hAacConfig->framelength != 120) {
        return AACENC_INVALID_CONFIG;
      }
      break;
    default:
      break;
  }

  if ((config->userSbrRatio == 0) && (isSbrActive(hAacConfig))) {
    if (hAacConfig->audioObjectType == AOT_ER_AAC_ELD) {
      hAacConfig->sbrRatio = ((hAacConfig->syntaxFlags & AC_LD_MPS) &&
                              (hAacConfig->sampleRate >= 27713))
                                 ? 2
                                 : 1;
    } else {
      hAacConfig->sbrRatio = 2;
    }
  } else {
    hAacConfig->sbrRatio = isSbrActive(hAacConfig) ? config->userSbrRatio : 0;
  }

  hAacConfig->bitRate = config->userBitrate;

  switch (hAacConfig->bitrateMode) {
    case AACENC_BR_MODE_CBR:
      if (config->userBitrate == (UINT)-1) {
        hAacConfig->bitRate =
            FDKaacEnc_GetCBRBitrate(hAacConfig, config->userSbrRatio);
      }
      hAacConfig->averageBits = -1;
      break;
    case AACENC_BR_MODE_VBR_1:
    case AACENC_BR_MODE_VBR_2:
    case AACENC_BR_MODE_VBR_3:
    case AACENC_BR_MODE_VBR_4:
    case AACENC_BR_MODE_VBR_5:
      if ((INT)config->userPeakBitrate != -1) {
        hAacConfig->bitrateMode = FDKaacEnc_AdjustVBRBitrateMode(
            hAacConfig->bitrateMode, config->userPeakBitrate,
            hAacConfig->channelMode);
      }
      hAacConfig->bitRate = FDKaacEnc_GetVBRBitrate(hAacConfig->bitrateMode,
                                                    hAacConfig->channelMode);
      break;
    default:
      return AACENC_INVALID_CONFIG;
  }

  switch (hAacConfig->bitrateMode) {
    case AACENC_BR_MODE_VBR_1:
    case AACENC_BR_MODE_VBR_2:
    case AACENC_BR_MODE_VBR_3:
    case AACENC_BR_MODE_VBR_4:
    case AACENC_BR_MODE_VBR_5:
    case AACENC_BR_MODE_CBR:
      if ((INT)config->userPeakBitrate != -1) {
        hAacConfig->maxBitsPerFrame =
            (FDKaacEnc_CalcBitsPerFrame(
                 fMax(hAacConfig->bitRate, (INT)config->userPeakBitrate),
                 hAacConfig->framelength, hAacConfig->sampleRate) +
             7) &
            ~7;
      } else {
        hAacConfig->maxBitsPerFrame = -1;
      }
      if (hAacConfig->audioMuxVersion == 2) {
        hAacConfig->minBitsPerFrame =
            fMin(32 * 8, FDKaacEnc_CalcBitsPerFrame(hAacConfig->bitRate,
                                                    hAacConfig->framelength,
                                                    hAacConfig->sampleRate)) &
            ~7;
      }
      break;
    default:
      return AACENC_INVALID_CONFIG;
  }

  if ((config->userTpNsubFrames > 1)) {
    int maxFrameLength = 8 * hAacEncoder->outBufferInBytes;
    switch (config->userTpType) {
      case TT_MP4_LOAS:
        maxFrameLength =
            fMin(maxFrameLength, 8 * (1 << 13)) / config->userTpNsubFrames;
        break;
      case TT_MP4_ADTS:
        maxFrameLength = fMin(maxFrameLength, 8 * ((1 << 13) - 1)) /
                         config->userTpNsubFrames;
        break;
      default:
        maxFrameLength = -1;
    }
    if (maxFrameLength != -1) {
      if (hAacConfig->maxBitsPerFrame > maxFrameLength) {
        return AACENC_INVALID_CONFIG;
      } else if (hAacConfig->maxBitsPerFrame == -1) {
        hAacConfig->maxBitsPerFrame = maxFrameLength;
      }
    }
  }

  if ((hAacConfig->audioObjectType == AOT_ER_AAC_ELD) &&
      !(hAacConfig->syntaxFlags & AC_ELD_DOWNSCALE) &&
      (config->userSbrEnabled == (UCHAR)-1) && (config->userSbrRatio == 0) &&
      ((hAacConfig->syntaxFlags & AC_LD_MPS) == 0)) {
    const ELD_SBR_CONFIGURATOR *pConfig = NULL;

    if (NULL !=
        (pConfig = eldSbrConfigurator(
             FDKaacEnc_GetChannelModeConfiguration(hAacConfig->channelMode)
                 ->nChannels,
             hAacConfig->sampleRate, hAacConfig->bitRate))) {
      hAacConfig->syntaxFlags |= (pConfig->sbrMode == 0) ? 0 : AC_SBR_PRESENT;
      hAacConfig->syntaxFlags |= (pConfig->chMode == MODE_212) ? AC_LD_MPS : 0;
      hAacConfig->channelMode =
          GetCoreChannelMode(pConfig->chMode, hAacConfig->audioObjectType);
      hAacConfig->nChannels =
          FDKaacEnc_GetChannelModeConfiguration(hAacConfig->channelMode)
              ->nChannels;
      hAacConfig->sbrRatio =
          (pConfig->sbrMode == 0) ? 0 : (pConfig->sbrMode == 1) ? 1 : 2;
    }
  }

  {
    UCHAR tpSignaling =
        getSbrSignalingMode(hAacConfig->audioObjectType, config->userTpType,
                            config->userTpSignaling, hAacConfig->sbrRatio);

    if ((hAacConfig->audioObjectType == AOT_AAC_LC ||
         hAacConfig->audioObjectType == AOT_SBR ||
         hAacConfig->audioObjectType == AOT_PS) &&
        (config->userTpType == TT_MP4_LATM_MCP1 ||
         config->userTpType == TT_MP4_LATM_MCP0 ||
         config->userTpType == TT_MP4_LOAS) &&
        (tpSignaling == 1) && (config->userTpAmxv == 0)) {
      return AACENC_INVALID_CONFIG;
    }

    if ((hAacConfig->audioObjectType == AOT_AAC_LC ||
         hAacConfig->audioObjectType == AOT_SBR ||
         hAacConfig->audioObjectType == AOT_PS) &&
        (tpSignaling == 0) && (hAacConfig->sbrRatio == 1)) {
      return AACENC_INVALID_CONFIG;
    }
  }

  switch (hAacConfig->bitrateMode) {
    case AACENC_BR_MODE_CBR:
    case AACENC_BR_MODE_VBR_1:
    case AACENC_BR_MODE_VBR_2:
    case AACENC_BR_MODE_VBR_3:
    case AACENC_BR_MODE_VBR_4:
    case AACENC_BR_MODE_VBR_5:
      if (0 >= (hAacConfig->bitRate = aacEncoder_LimitBitrate(
                    NULL, hAacConfig->sampleRate, hAacConfig->framelength,
                    hAacConfig->nChannels, hAacConfig->channelMode,
                    hAacConfig->bitRate, hAacConfig->nSubFrames,
                    isSbrActive(hAacConfig), hAacConfig->sbrRatio,
                    hAacConfig->syntaxFlags, hAacConfig->audioObjectType))) {
        return AACENC_INVALID_CONFIG;
      }
      break;
    default:
      break;
  }

  if (AACENC_BR_MODE_IS_VBR(hAacConfig->bitrateMode)
      || (hAacConfig->useTns == 0))
  {
    hAacConfig->usePns = 0;
  }

  if (hAacConfig->epConfig >= 0) {
    hAacConfig->syntaxFlags |= AC_ER;
    if (((INT)hAacConfig->channelMode < 1) ||
        ((INT)hAacConfig->channelMode > 14)) {
      return AACENC_INVALID_CONFIG;
    }
  }

  if ((hAacConfig->syntaxFlags & AC_LD_MPS) == 0) {
    if (FDKaacEnc_DetermineEncoderMode(&hAacConfig->channelMode,
                                       hAacConfig->nChannels) != AAC_ENC_OK) {
      return AACENC_INVALID_CONFIG;
    }
  }

  if ((hAacConfig->nChannels > hAacEncoder->nMaxAacChannels) ||
      ((FDKaacEnc_GetChannelModeConfiguration(hAacConfig->channelMode)
            ->nChannelsEff > hAacEncoder->nMaxSbrChannels) &&
       isSbrActive(hAacConfig))) {
    return AACENC_INVALID_CONFIG;
  }

  switch (hAacConfig->audioObjectType) {
    case AOT_AAC_LC:
    case AOT_SBR:
    case AOT_PS:
    case AOT_MP2_AAC_LC:
    case AOT_MP2_SBR:
      hAacEncoder->metaDataAllowed = 1;
      if (!((((INT)hAacConfig->channelMode >= 1) &&
             ((INT)hAacConfig->channelMode <= 14)) ||
            (MODE_7_1_REAR_SURROUND == hAacConfig->channelMode) ||
            (MODE_7_1_FRONT_CENTER == hAacConfig->channelMode))) {
        config->userMetaDataMode = 0;
      }
      break;
    default:
      hAacEncoder->metaDataAllowed = 0;
  }

  return err;
}

static INT aacenc_SbrCallback(void *self, HANDLE_FDK_BITSTREAM hBs,
                              const INT sampleRateIn, const INT sampleRateOut,
                              const INT samplesPerFrame,
                              const AUDIO_OBJECT_TYPE coreCodec,
                              const MP4_ELEMENT_ID elementID,
                              const INT elementIndex, const UCHAR harmonicSbr,
                              const UCHAR stereoConfigIndex,
                              const UCHAR configMode, UCHAR *configChanged,
                              const INT downscaleFactor) {
  HANDLE_AACENCODER hAacEncoder = (HANDLE_AACENCODER)self;

  sbrEncoder_GetHeader(hAacEncoder->hEnvEnc, hBs, elementIndex, 0);

  return 0;
}

INT aacenc_SscCallback(void *self, HANDLE_FDK_BITSTREAM hBs,
                       const AUDIO_OBJECT_TYPE coreCodec,
                       const INT samplingRate, const INT frameSize,
                       const INT numChannels, const INT stereoConfigIndex,
                       const INT coreSbrFrameLengthIndex, const INT configBytes,
                       const UCHAR configMode, UCHAR *configChanged) {
  HANDLE_AACENCODER hAacEncoder = (HANDLE_AACENCODER)self;

  return (FDK_MpegsEnc_WriteSpatialSpecificConfig(hAacEncoder->hMpsEnc, hBs));
}

static AACENC_ERROR aacEncInit(HANDLE_AACENCODER hAacEncoder, ULONG InitFlags,
                               USER_PARAM *config) {
  AACENC_ERROR err = AACENC_OK;

  INT aacBufferOffset = 0;
  HANDLE_SBR_ENCODER *hSbrEncoder = &hAacEncoder->hEnvEnc;
  HANDLE_AACENC_CONFIG hAacConfig = &hAacEncoder->aacConfig;

  hAacEncoder->nZerosAppended = 0;

  INT frameLength = hAacConfig->framelength;

  if ((InitFlags & AACENC_INIT_CONFIG)) {
    CHANNEL_MODE prevChMode = hAacConfig->channelMode;

    if ((err = FDKaacEnc_AdjustEncSettings(hAacEncoder, config)) != AACENC_OK) {
      return err;
    }
    frameLength = hAacConfig->framelength;

    if ((prevChMode != hAacConfig->channelMode) && isSbrActive(hAacConfig)) {
      InitFlags |= AACENC_INIT_STATES;
    }
  }

  if (InitFlags == AACENC_INIT_ALL) {
    FDKmemclear(hAacEncoder->inputBuffer,
                sizeof(INT_PCM) * hAacEncoder->inputBufferSize);
  }

  if ((InitFlags & AACENC_INIT_CONFIG)) {
    aacBufferOffset = 0;
    switch (hAacConfig->audioObjectType) {
      case AOT_ER_AAC_LD:
        hAacEncoder->nDelay = DELAY_AACLD(hAacConfig->framelength);
        break;
      case AOT_ER_AAC_ELD:
        hAacEncoder->nDelay = DELAY_AACELD(hAacConfig->framelength);
        break;
      default:
        hAacEncoder->nDelay =
            DELAY_AAC(hAacConfig->framelength);
    }

    hAacConfig->ancDataBitRate = 0;
  }

  if ((NULL != hAacEncoder->hEnvEnc) && isSbrActive(hAacConfig) &&
      ((InitFlags & AACENC_INIT_CONFIG) || (InitFlags & AACENC_INIT_STATES))) {
    INT sbrError;
    UINT initFlag = 0;
    SBR_ELEMENT_INFO sbrElInfo[(8)];
    CHANNEL_MAPPING channelMapping;
    CHANNEL_MODE channelMode = isPsActive(hAacConfig->audioObjectType)
                                   ? config->userChannelMode
                                   : hAacConfig->channelMode;
    INT numChannels = isPsActive(hAacConfig->audioObjectType)
                          ? config->nChannels
                          : hAacConfig->nChannels;

    if (FDKaacEnc_InitChannelMapping(channelMode, hAacConfig->channelOrder,
                                     &channelMapping) != AAC_ENC_OK) {
      return AACENC_INIT_ERROR;
    }

    if (channelMapping.nElements > (8)) {
      return AACENC_INIT_ERROR;
    }

    aacEncDistributeSbrBits(&channelMapping, sbrElInfo, hAacConfig->bitRate);

    initFlag += (InitFlags & AACENC_INIT_STATES) ? 1 : 0;

    sbrError = sbrEncoder_Init(
        *hSbrEncoder, sbrElInfo, channelMapping.nElements,
        hAacEncoder->inputBuffer, hAacEncoder->inputBufferSizePerChannel,
        &hAacConfig->bandWidth, &aacBufferOffset, &numChannels,
        hAacConfig->syntaxFlags, &hAacConfig->sampleRate, &hAacConfig->sbrRatio,
        &frameLength, hAacConfig->audioObjectType, &hAacEncoder->nDelay,
        (hAacConfig->audioObjectType == AOT_ER_AAC_ELD) ? 1 : TRANS_FAC,
        (config->userTpHeaderPeriod != 0xFF)
            ? config->userTpHeaderPeriod
            : DEFAULT_HEADER_PERIOD_REPETITION_RATE,
        initFlag);

    if ((sbrError) || (numChannels != hAacConfig->nChannels)) {
      return AACENC_INIT_SBR_ERROR;
    }

    if (numChannels == 1) {
      hAacConfig->channelMode = MODE_1;
    }

    if (hAacConfig->usePns) {
      hAacConfig->usePns = 0;
    }

    hAacConfig->ancDataBitRate = sbrEncoder_GetEstimateBitrate(*hSbrEncoder);
  }

  if ((hAacEncoder->hMpsEnc != NULL) && (hAacConfig->syntaxFlags & AC_LD_MPS)) {
    int coreCoderDelay = DELAY_AACELD(hAacConfig->framelength);

    if (isSbrActive(hAacConfig)) {
      coreCoderDelay = hAacConfig->sbrRatio * coreCoderDelay +
                       sbrEncoder_GetInputDataDelay(*hSbrEncoder);
    }

    if (MPS_ENCODER_OK !=
        FDK_MpegsEnc_Init(hAacEncoder->hMpsEnc, hAacConfig->audioObjectType,
                          config->userSamplerate, hAacConfig->bitRate,
                          isSbrActive(hAacConfig) ? hAacConfig->sbrRatio : 0,
                          frameLength,
                          hAacEncoder->inputBufferSizePerChannel,
                          coreCoderDelay)) {
      return AACENC_INIT_MPS_ERROR;
    }
  }
  hAacEncoder->nDelay =
      fMax(FDK_MpegsEnc_GetDelay(hAacEncoder->hMpsEnc), hAacEncoder->nDelay);

  if ((InitFlags & AACENC_INIT_TRANSPORT)) {
    UINT flags = 0;

    FDKaacEnc_MapConfig(
        &hAacEncoder->coderConfig, config,
        getSbrSignalingMode(hAacConfig->audioObjectType, config->userTpType,
                            config->userTpSignaling, hAacConfig->sbrRatio),
        hAacConfig);

    if (config->userTpAmxv != 0) {
      flags |= TP_FLAG_LATM_AMV;
    }
    FDKmemclear(hAacEncoder->outBuffer,
                hAacEncoder->outBufferInBytes * sizeof(UCHAR));

    if (transportEnc_Init(hAacEncoder->hTpEnc, hAacEncoder->outBuffer,
                          hAacEncoder->outBufferInBytes, config->userTpType,
                          &hAacEncoder->coderConfig, flags) != 0) {
      return AACENC_INIT_TP_ERROR;
    }
  }

  if ((InitFlags & AACENC_INIT_CONFIG) || (InitFlags & AACENC_INIT_STATES)) {
    if (FDKaacEnc_Initialize(
            hAacEncoder->hAacEnc, hAacConfig, hAacEncoder->hTpEnc,
            (InitFlags & AACENC_INIT_STATES) ? 1 : 0) != AAC_ENC_OK) {
      return AACENC_INIT_AAC_ERROR;
    }
  }

  if (hAacEncoder->hMetadataEnc && (hAacEncoder->metaDataAllowed != 0) &&
      ((InitFlags & AACENC_INIT_CONFIG) || (InitFlags & AACENC_INIT_STATES))) {
    INT inputDataDelay = DELAY_AAC(hAacConfig->framelength);

    if (isSbrActive(hAacConfig) && hSbrEncoder != NULL) {
      inputDataDelay = hAacConfig->sbrRatio * inputDataDelay +
                       sbrEncoder_GetInputDataDelay(*hSbrEncoder);
    }

    if (FDK_MetadataEnc_Init(hAacEncoder->hMetadataEnc,
                             ((InitFlags & AACENC_INIT_STATES) ? 1 : 0),
                             config->userMetaDataMode, inputDataDelay,
                             frameLength, config->userSamplerate,
                             config->nChannels, config->userChannelMode,
                             hAacConfig->channelOrder) != 0) {
      return AACENC_INIT_META_ERROR;
    }

    hAacEncoder->nDelay += FDK_MetadataEnc_GetDelay(hAacEncoder->hMetadataEnc);
  }

  if ((hAacEncoder->hMpsEnc != NULL) && (hAacConfig->syntaxFlags & AC_LD_MPS)) {
    hAacEncoder->nDelayCore =
        hAacEncoder->nDelay -
        fMax(0, FDK_MpegsEnc_GetDecDelay(hAacEncoder->hMpsEnc));
  } else if (isSbrActive(hAacConfig) && hSbrEncoder != NULL) {
    hAacEncoder->nDelayCore =
        hAacEncoder->nDelay -
        fMax(0, sbrEncoder_GetSbrDecDelay(hAacEncoder->hEnvEnc));
  } else {
    hAacEncoder->nDelayCore = hAacEncoder->nDelay;
  }

  if ((InitFlags & AACENC_INIT_CONFIG)) {
    hAacEncoder->inputBufferOffset = aacBufferOffset;

    hAacEncoder->nSamplesToRead = frameLength * config->nChannels;
  }

  return AACENC_OK;
}

LINKSPEC(AACENC_ERROR) aacEncOpen(HANDLE_AACENCODER *phAacEncoder, const UINT encModules,
                        const UINT maxChannels) {
  AACENC_ERROR err = AACENC_OK;
  HANDLE_AACENCODER hAacEncoder = NULL;

  if (phAacEncoder == NULL) {
    err = AACENC_INVALID_HANDLE;
    goto bail;
  }

  hAacEncoder = Get_AacEncoder();

  if (hAacEncoder == NULL) {
    err = AACENC_MEMORY_ERROR;
    goto bail;
  }

  FDKmemclear(hAacEncoder, sizeof(AACENCODER));

  if (encModules == 0) {
    C_ALLOC_SCRATCH_START(_pLibInfo, LIB_INFO, FDK_MODULE_LAST)
    LIB_INFO(*pLibInfo)
    [FDK_MODULE_LAST] = (LIB_INFO(*)[FDK_MODULE_LAST])_pLibInfo;
    FDKinitLibInfo(*pLibInfo);
    aacEncGetLibInfo(*pLibInfo);

    hAacEncoder->encoder_modis = ENC_MODE_FLAG_AAC;
    if (FDKlibInfo_getCapabilities(*pLibInfo, FDK_SBRENC) & CAPF_SBR_HQ) {
      hAacEncoder->encoder_modis |= ENC_MODE_FLAG_SBR;
    }
    if (FDKlibInfo_getCapabilities(*pLibInfo, FDK_SBRENC) & CAPF_SBR_PS_MPEG) {
      hAacEncoder->encoder_modis |= ENC_MODE_FLAG_PS;
    }
    if (FDKlibInfo_getCapabilities(*pLibInfo, FDK_AACENC) & CAPF_AAC_DRC) {
      hAacEncoder->encoder_modis |= ENC_MODE_FLAG_META;
    }
    hAacEncoder->encoder_modis |= ENC_MODE_FLAG_SAC;

    C_ALLOC_SCRATCH_END(_pLibInfo, LIB_INFO, FDK_MODULE_LAST)
  } else {
    hAacEncoder->encoder_modis = encModules;
  }

  if (maxChannels == 0) {
    hAacEncoder->nMaxAacChannels = (8);
    hAacEncoder->nMaxSbrChannels = (8);
  } else {
    hAacEncoder->nMaxAacChannels = (maxChannels & 0x00FF);
    if ((hAacEncoder->encoder_modis & ENC_MODE_FLAG_SBR)) {
      hAacEncoder->nMaxSbrChannels = (maxChannels & 0xFF00)
                                         ? (maxChannels >> 8)
                                         : hAacEncoder->nMaxAacChannels;
    }

    if ((hAacEncoder->nMaxAacChannels > (8)) ||
        (hAacEncoder->nMaxSbrChannels > (8))) {
      err = AACENC_INVALID_CONFIG;
      goto bail;
    }
  }

  hAacEncoder->nMaxAacElements = fixMin(((8)), hAacEncoder->nMaxAacChannels);
  hAacEncoder->nMaxSbrElements = fixMin((8), hAacEncoder->nMaxSbrChannels);

  if (hAacEncoder->encoder_modis & (ENC_MODE_FLAG_SBR | ENC_MODE_FLAG_PS))
    hAacEncoder->inputBufferSizePerChannel = INPUTBUFFER_SIZE;
  else
    hAacEncoder->inputBufferSizePerChannel = (1024);

  hAacEncoder->inputBufferSize =
      hAacEncoder->nMaxAacChannels * hAacEncoder->inputBufferSizePerChannel;

  if (NULL == (hAacEncoder->inputBuffer = (INT_PCM *)FDKcalloc(
                   hAacEncoder->inputBufferSize, sizeof(INT_PCM)))) {
    err = AACENC_MEMORY_ERROR;
    goto bail;
  }

  if (hAacEncoder->encoder_modis & ENC_MODE_FLAG_SBR) {
    if (sbrEncoder_Open(
            &hAacEncoder->hEnvEnc, hAacEncoder->nMaxSbrElements,
            hAacEncoder->nMaxSbrChannels,
            (hAacEncoder->encoder_modis & ENC_MODE_FLAG_PS) ? 1 : 0)) {
      err = AACENC_MEMORY_ERROR;
      goto bail;
    }

    if (NULL == (hAacEncoder->pSbrPayload = (SBRENC_EXT_PAYLOAD *)FDKcalloc(
                     1, sizeof(SBRENC_EXT_PAYLOAD)))) {
      err = AACENC_MEMORY_ERROR;
      goto bail;
    }
  }

  if (FDKaacEnc_Open(&hAacEncoder->hAacEnc, hAacEncoder->nMaxAacElements,
                     hAacEncoder->nMaxAacChannels, (1)) != AAC_ENC_OK) {
    err = AACENC_MEMORY_ERROR;
    goto bail;
  }

  hAacEncoder->outBufferInBytes =
      1 << (DFRACT_BITS - CntLeadingZeros(fixMax(
                              1, ((1) * hAacEncoder->nMaxAacChannels * 6144) >>
                                     3)));
  if (NULL == (hAacEncoder->outBuffer = (UCHAR *)FDKcalloc(
                   hAacEncoder->outBufferInBytes, sizeof(UCHAR)))) {
    err = AACENC_MEMORY_ERROR;
    goto bail;
  }

  if (hAacEncoder->encoder_modis & ENC_MODE_FLAG_META) {
    if (FDK_MetadataEnc_Open(&hAacEncoder->hMetadataEnc,
                             (UINT)hAacEncoder->nMaxAacChannels)) {
      err = AACENC_MEMORY_ERROR;
      goto bail;
    }
  }

  if (hAacEncoder->encoder_modis & ENC_MODE_FLAG_SAC) {
    if (MPS_ENCODER_OK != FDK_MpegsEnc_Open(&hAacEncoder->hMpsEnc)) {
      err = AACENC_MEMORY_ERROR;
      goto bail;
    }
  }

  if (transportEnc_Open(&hAacEncoder->hTpEnc) != 0) {
    err = AACENC_MEMORY_ERROR;
    goto bail;
  } else {
    C_ALLOC_SCRATCH_START(_pLibInfo, LIB_INFO, FDK_MODULE_LAST)

    LIB_INFO(*pLibInfo)
    [FDK_MODULE_LAST] = (LIB_INFO(*)[FDK_MODULE_LAST])_pLibInfo;

    FDKinitLibInfo(*pLibInfo);
    transportEnc_GetLibInfo(*pLibInfo);

    hAacEncoder->CAPF_tpEnc = FDKlibInfo_getCapabilities(*pLibInfo, FDK_TPENC);

    C_ALLOC_SCRATCH_END(_pLibInfo, LIB_INFO, FDK_MODULE_LAST)
  }
  if (transportEnc_RegisterSbrCallback(hAacEncoder->hTpEnc, aacenc_SbrCallback,
                                       hAacEncoder) != 0) {
    err = AACENC_INIT_TP_ERROR;
    goto bail;
  }
  if (transportEnc_RegisterSscCallback(hAacEncoder->hTpEnc, aacenc_SscCallback,
                                       hAacEncoder) != 0) {
    err = AACENC_INIT_TP_ERROR;
    goto bail;
  }

  aacEncDefaultConfig(&hAacEncoder->aacConfig, &hAacEncoder->extParam);

  hAacEncoder->coderConfig.headerPeriod =
      hAacEncoder->extParam.userTpHeaderPeriod;

  hAacEncoder->InitFlags = AACENC_INIT_ALL;

  *phAacEncoder = hAacEncoder;

  return err;

bail:
  aacEncClose(&hAacEncoder);

  return err;
}

LINKSPEC(AACENC_ERROR) aacEncClose(HANDLE_AACENCODER *phAacEncoder) {
  AACENC_ERROR err = AACENC_OK;

  if (phAacEncoder == NULL) {
    err = AACENC_INVALID_HANDLE;
    goto bail;
  }

  if (*phAacEncoder != NULL) {
    HANDLE_AACENCODER hAacEncoder = *phAacEncoder;

    if (hAacEncoder->inputBuffer != NULL) {
      FDKfree(hAacEncoder->inputBuffer);
      hAacEncoder->inputBuffer = NULL;
    }
    if (hAacEncoder->outBuffer != NULL) {
      FDKfree(hAacEncoder->outBuffer);
      hAacEncoder->outBuffer = NULL;
    }

    if (hAacEncoder->hEnvEnc) {
      sbrEncoder_Close(&hAacEncoder->hEnvEnc);
    }
    if (hAacEncoder->pSbrPayload != NULL) {
      FDKfree(hAacEncoder->pSbrPayload);
      hAacEncoder->pSbrPayload = NULL;
    }
    if (hAacEncoder->hAacEnc) {
      FDKaacEnc_Close(&hAacEncoder->hAacEnc);
    }

    transportEnc_Close(&hAacEncoder->hTpEnc);

    if (hAacEncoder->hMetadataEnc) {
      FDK_MetadataEnc_Close(&hAacEncoder->hMetadataEnc);
    }
    if (hAacEncoder->hMpsEnc) {
      FDK_MpegsEnc_Close(&hAacEncoder->hMpsEnc);
    }

    Free_AacEncoder(phAacEncoder);
  }

bail:
  return err;
}

LINKSPEC(AACENC_ERROR) aacEncEncode(const HANDLE_AACENCODER hAacEncoder,
                          const AACENC_BufDesc *inBufDesc,
                          const AACENC_BufDesc *outBufDesc,
                          const AACENC_InArgs *inargs,
                          AACENC_OutArgs *outargs) {
  AACENC_ERROR err = AACENC_OK;
  INT i, nBsBytes = 0;
  INT outBytes[(1)];
  int nExtensions = 0;
  int ancDataExtIdx = -1;

  if (hAacEncoder == NULL) {
    err = AACENC_INVALID_HANDLE;
    goto bail;
  }

  if (hAacEncoder->InitFlags != 0) {
    err =
        aacEncInit(hAacEncoder, hAacEncoder->InitFlags, &hAacEncoder->extParam);

    if (err != AACENC_OK) {
      goto bail;
    }
    hAacEncoder->InitFlags = AACENC_INIT_NONE;
  }

  if (outargs != NULL) {
    FDKmemclear(outargs, sizeof(AACENC_OutArgs));
  }

  if (outBufDesc != NULL) {
    for (i = 0; i < outBufDesc->numBufs; i++) {
      if (outBufDesc->bufs[i] != NULL) {
        FDKmemclear(outBufDesc->bufs[i], outBufDesc->bufSizes[i]);
      }
    }
  }

  if ((inBufDesc == NULL) && (outBufDesc == NULL) && (inargs == NULL) &&
      (outargs == NULL)) {
    goto bail;
  }

  if ((inargs == NULL) || (outargs == NULL) ||
      ((AACENC_OK != validateBufDesc(inBufDesc)) &&
       (inargs->numInSamples > 0)) ||
      (AACENC_OK != validateBufDesc(outBufDesc))) {
    err = AACENC_UNSUPPORTED_PARAMETER;
    goto bail;
  }

  FDKmemclear(outBytes, hAacEncoder->aacConfig.nSubFrames * sizeof(INT));

  if ((inBufDesc != NULL) && (inargs->numInSamples > 0) &&
      (getBufDescIdx(inBufDesc, IN_AUDIO_DATA) != -1)) {
    INT idx = getBufDescIdx(inBufDesc, IN_AUDIO_DATA);
    INT newSamples =
        fixMax(0, fixMin(inargs->numInSamples, hAacEncoder->nSamplesToRead -
                                                   hAacEncoder->nSamplesRead));
    INT_PCM *pIn =
        hAacEncoder->inputBuffer +
        hAacEncoder->inputBufferOffset / hAacEncoder->aacConfig.nChannels +
        hAacEncoder->nSamplesRead / hAacEncoder->extParam.nChannels;
    newSamples -=
        (newSamples %
         hAacEncoder->extParam
             .nChannels);

    if (inBufDesc->bufElSizes[idx] == (INT)sizeof(INT_PCM)) {
      FDK_deinterleave((INT_PCM *)inBufDesc->bufs[idx], pIn,
                       hAacEncoder->extParam.nChannels,
                       newSamples / hAacEncoder->extParam.nChannels,
                       hAacEncoder->inputBufferSizePerChannel);
    } else if (inBufDesc->bufElSizes[idx] > (INT)sizeof(INT_PCM)) {
      FDK_deinterleave((LONG *)inBufDesc->bufs[idx], pIn,
                       hAacEncoder->extParam.nChannels,
                       newSamples / hAacEncoder->extParam.nChannels,
                       hAacEncoder->inputBufferSizePerChannel);
    } else {
      FDK_deinterleave((SHORT *)inBufDesc->bufs[idx], pIn,
                       hAacEncoder->extParam.nChannels,
                       newSamples / hAacEncoder->extParam.nChannels,
                       hAacEncoder->inputBufferSizePerChannel);
    }
    hAacEncoder->nSamplesRead += newSamples;

    outargs->numInSamples = newSamples;
  }

  if (hAacEncoder->nSamplesRead < hAacEncoder->nSamplesToRead) {
    if (inargs->numInSamples == -1) {
      if ((hAacEncoder->nZerosAppended < hAacEncoder->nDelay)) {
        int nZeros = (hAacEncoder->nSamplesToRead - hAacEncoder->nSamplesRead) /
                     hAacEncoder->extParam.nChannels;

        FDK_ASSERT(nZeros >= 0);

        if (nZeros) {
          INT_PCM *pIn =
              hAacEncoder->inputBuffer +
              hAacEncoder->inputBufferOffset /
                  hAacEncoder->aacConfig.nChannels +
              hAacEncoder->nSamplesRead / hAacEncoder->extParam.nChannels;
          for (i = 0; i < (int)hAacEncoder->extParam.nChannels; i++) {
            FDKmemclear(pIn + i * hAacEncoder->inputBufferSizePerChannel,
                        sizeof(INT_PCM) * nZeros);
          }
          hAacEncoder->nZerosAppended += nZeros;
          hAacEncoder->nSamplesRead = hAacEncoder->nSamplesToRead;
        }
      } else {
        err = AACENC_ENCODE_EOF;
        goto bail;
      }
    } else {
      goto bail;
    }
  }

  FDKmemclear(hAacEncoder->extPayload,
              sizeof(AACENC_EXT_PAYLOAD) * MAX_TOTAL_EXT_PAYLOADS);
  for (i = 0; i < MAX_TOTAL_EXT_PAYLOADS; i++) {
    hAacEncoder->extPayload[i].associatedChElement = -1;
  }
  if (hAacEncoder->pSbrPayload != NULL) {
    FDKmemclear(hAacEncoder->pSbrPayload, sizeof(*hAacEncoder->pSbrPayload));
  }

  if ((hAacEncoder->hMetadataEnc != NULL) &&
      (hAacEncoder->metaDataAllowed != 0)) {
    const AACENC_MetaData *pMetaData = NULL;
    AACENC_EXT_PAYLOAD *pMetaDataExtPayload = NULL;
    UINT nMetaDataExtensions = 0;
    INT matrix_mixdown_idx = 0;

    if (getBufDescIdx(inBufDesc, IN_METADATA_SETUP) != -1) {
      pMetaData =
          (AACENC_MetaData *)
              inBufDesc->bufs[getBufDescIdx(inBufDesc, IN_METADATA_SETUP)];
    }

    FDK_MetadataEnc_Process(
        hAacEncoder->hMetadataEnc,
        hAacEncoder->inputBuffer + hAacEncoder->inputBufferOffset /
                                       hAacEncoder->coderConfig.noChannels,
        hAacEncoder->inputBufferSizePerChannel, hAacEncoder->nSamplesRead,
        pMetaData, &pMetaDataExtPayload, &nMetaDataExtensions,
        &matrix_mixdown_idx);

    for (i = 0; i < (INT)nMetaDataExtensions;
         i++) {
      hAacEncoder->extPayload[nExtensions++] = pMetaDataExtPayload[i];
    }

    if ((matrix_mixdown_idx != -1) &&
        ((hAacEncoder->extParam.userChannelMode == MODE_1_2_2) ||
         (hAacEncoder->extParam.userChannelMode == MODE_1_2_2_1))) {
      UINT pceValue = (UINT)((0 << 3) | ((matrix_mixdown_idx & 0x3) << 1) | 1);
      if (hAacEncoder->extParam.userPceAdditions != pceValue) {
        hAacEncoder->extParam.userPceAdditions = pceValue;
        hAacEncoder->InitFlags |= AACENC_INIT_TRANSPORT;
      }
    }
  }

  if ((hAacEncoder->hMpsEnc != NULL) &&
      (hAacEncoder->aacConfig.syntaxFlags & AC_LD_MPS)) {
    AACENC_EXT_PAYLOAD mpsExtensionPayload;
    FDKmemclear(&mpsExtensionPayload, sizeof(AACENC_EXT_PAYLOAD));

    if (MPS_ENCODER_OK !=
        FDK_MpegsEnc_Process(
            hAacEncoder->hMpsEnc,
            hAacEncoder->inputBuffer + hAacEncoder->inputBufferOffset /
                                           hAacEncoder->coderConfig.noChannels,
            hAacEncoder->nSamplesRead, &mpsExtensionPayload)) {
      err = AACENC_ENCODE_ERROR;
      goto bail;
    }

    if ((mpsExtensionPayload.pData != NULL) &&
        ((mpsExtensionPayload.dataSize != 0))) {
      hAacEncoder->extPayload[nExtensions++] = mpsExtensionPayload;
    }
  }

  if ((NULL != hAacEncoder->hEnvEnc) && (NULL != hAacEncoder->pSbrPayload) &&
      isSbrActive(&hAacEncoder->aacConfig)) {
    INT nPayload = 0;

    if (sbrEncoder_EncodeFrame(hAacEncoder->hEnvEnc, hAacEncoder->inputBuffer,
                               hAacEncoder->inputBufferSizePerChannel,
                               hAacEncoder->pSbrPayload->dataSize[nPayload],
                               hAacEncoder->pSbrPayload->data[nPayload])) {
      err = AACENC_ENCODE_ERROR;
      goto bail;
    } else {
      for (i = 0; i < (8); i++) {
        if (hAacEncoder->pSbrPayload->dataSize[nPayload][i] > 0) {
          hAacEncoder->extPayload[nExtensions].pData =
              hAacEncoder->pSbrPayload->data[nPayload][i];
          {
            hAacEncoder->extPayload[nExtensions].dataSize =
                hAacEncoder->pSbrPayload->dataSize[nPayload][i];
            hAacEncoder->extPayload[nExtensions].associatedChElement =
                getAssociatedChElement(
                    &hAacEncoder->hEnvEnc->sbrElement[i]->elInfo,
                    &hAacEncoder->hAacEnc->channelMapping);
            if (hAacEncoder->extPayload[nExtensions].associatedChElement ==
                -1) {
              err = AACENC_ENCODE_ERROR;
              goto bail;
            }
          }
          hAacEncoder->extPayload[nExtensions].dataType =
              EXT_SBR_DATA;
          nExtensions++;
          FDK_ASSERT(nExtensions <= MAX_TOTAL_EXT_PAYLOADS);
        }
      }
      nPayload++;
    }
  }

  if ((inargs->numAncBytes > 0) &&
      (getBufDescIdx(inBufDesc, IN_ANCILLRY_DATA) != -1)) {
    INT idx = getBufDescIdx(inBufDesc, IN_ANCILLRY_DATA);
    hAacEncoder->extPayload[nExtensions].dataSize = inargs->numAncBytes * 8;
    hAacEncoder->extPayload[nExtensions].pData = (UCHAR *)inBufDesc->bufs[idx];
    hAacEncoder->extPayload[nExtensions].dataType = EXT_DATA_ELEMENT;
    hAacEncoder->extPayload[nExtensions].associatedChElement = -1;
    ancDataExtIdx = nExtensions;
    nExtensions++;
  }

  if (FDKaacEnc_EncodeFrame(hAacEncoder->hAacEnc, hAacEncoder->hTpEnc,
                            hAacEncoder->inputBuffer,
                            hAacEncoder->inputBufferSizePerChannel, outBytes,
                            hAacEncoder->extPayload) != AAC_ENC_OK) {
    err = AACENC_ENCODE_ERROR;
    goto bail;
  }

  if (ancDataExtIdx >= 0) {
    outargs->numAncBytes =
        inargs->numAncBytes -
        (hAacEncoder->extPayload[ancDataExtIdx].dataSize >> 3);
  }

  hAacEncoder->nSamplesRead -= hAacEncoder->nSamplesToRead;

  if (isSbrActive(&hAacEncoder->aacConfig)) {
    sbrEncoder_UpdateBuffers(hAacEncoder->hEnvEnc, hAacEncoder->inputBuffer,
                             hAacEncoder->inputBufferSizePerChannel);
  }

  if ((outBufDesc != NULL) && (outBufDesc->numBufs >= 1)) {
    INT bsIdx = getBufDescIdx(outBufDesc, OUT_BITSTREAM_DATA);
    INT auIdx = getBufDescIdx(outBufDesc, OUT_AU_SIZES);

    for (i = 0, nBsBytes = 0; i < hAacEncoder->aacConfig.nSubFrames; i++) {
      nBsBytes += outBytes[i];

      if (auIdx != -1) {
        ((INT *)outBufDesc->bufs[auIdx])[i] = outBytes[i];
      }
    }

    if ((bsIdx != -1) && (outBufDesc->bufSizes[bsIdx] >= nBsBytes)) {
      FDKmemcpy(outBufDesc->bufs[bsIdx], hAacEncoder->outBuffer,
                sizeof(UCHAR) * nBsBytes);
      outargs->numOutBytes = nBsBytes;
      outargs->bitResState =
          FDKaacEnc_GetBitReservoirState(hAacEncoder->hAacEnc);
    } else {
      err = AACENC_ENCODE_ERROR;
      goto bail;
    }
  }

bail:
  if (err == AACENC_ENCODE_ERROR) {
    hAacEncoder->InitFlags = AACENC_INIT_ALL;
  }

  return err;
}

static AAC_ENCODER_ERROR aacEncGetConf(HANDLE_AACENCODER hAacEncoder,
                                       UINT *size, UCHAR *confBuffer) {
  FDK_BITSTREAM tmpConf;
  UINT confType;
  UCHAR buf[64];
  int err;

  FDKinitBitStream(&tmpConf, buf, 64, 0, BS_WRITER);

  err = transportEnc_GetConf(hAacEncoder->hTpEnc, &hAacEncoder->coderConfig,
                             &tmpConf, &confType);

  FDKbyteAlign(&tmpConf, 0);

  if (FDKgetValidBits(&tmpConf) > ((*size) << 3)) return AAC_ENC_UNKNOWN;

  FDKfetchBuffer(&tmpConf, confBuffer, size);

  if (err != 0)
    return AAC_ENC_UNKNOWN;
  else
    return AAC_ENC_OK;
}

LINKSPEC(AACENC_ERROR) aacEncGetLibInfo(LIB_INFO *info) {
  int i = 0;

  if (info == NULL) {
    return AACENC_INVALID_HANDLE;
  }

  FDK_toolsGetLibInfo(info);
  transportEnc_GetLibInfo(info);
  sbrEncoder_GetLibInfo(info);
  FDK_MpegsEnc_GetLibInfo(info);

  for (i = 0; i < FDK_MODULE_LAST; i++) {
    if (info[i].module_id == FDK_NONE) break;
  }
  if (i == FDK_MODULE_LAST) {
    return AACENC_INIT_ERROR;
  }

  info[i].module_id = FDK_AACENC;
  info[i].build_date = AACENCODER_LIB_BUILD_DATE;
  info[i].build_time = AACENCODER_LIB_BUILD_TIME;
  info[i].title = AACENCODER_LIB_TITLE;
  info[i].version =
      LIB_VERSION(AACENCODER_LIB_VL0, AACENCODER_LIB_VL1, AACENCODER_LIB_VL2);
  ;
  LIB_VERSION_STRING(&info[i]);

  info[i].flags = 0 | CAPF_AAC_1024 | CAPF_AAC_LC | CAPF_AAC_512 |
                  CAPF_AAC_480 | CAPF_AAC_DRC | CAPF_AAC_ELD_DOWNSCALE;

  return AACENC_OK;
}

LINKSPEC(AACENC_ERROR) aacEncoder_SetParam(const HANDLE_AACENCODER hAacEncoder,
                                 const AACENC_PARAM param, const UINT value) {
  AACENC_ERROR err = AACENC_OK;
  USER_PARAM *settings = &hAacEncoder->extParam;

  if (hAacEncoder == NULL) {
    err = AACENC_INVALID_HANDLE;
    goto bail;
  }

  switch (param) {
    case AACENC_AOT:
      if (settings->userAOT != (AUDIO_OBJECT_TYPE)value) {
        switch (value) {
          case AOT_PS:
            if (!(hAacEncoder->encoder_modis & (ENC_MODE_FLAG_PS))) {
              err = AACENC_INVALID_CONFIG;
              goto bail;
            }
            FDK_FALLTHROUGH;
          case AOT_SBR:
          case AOT_MP2_SBR:
            if (!(hAacEncoder->encoder_modis & (ENC_MODE_FLAG_SBR))) {
              err = AACENC_INVALID_CONFIG;
              goto bail;
            }
            FDK_FALLTHROUGH;
          case AOT_AAC_LC:
          case AOT_MP2_AAC_LC:
          case AOT_ER_AAC_LD:
          case AOT_ER_AAC_ELD:
            if (!(hAacEncoder->encoder_modis & (ENC_MODE_FLAG_AAC))) {
              err = AACENC_INVALID_CONFIG;
              goto bail;
            }
            break;
          default:
            err = AACENC_INVALID_CONFIG;
            goto bail;
        }
        settings->userAOT = (AUDIO_OBJECT_TYPE)value;
        hAacEncoder->InitFlags |=
            AACENC_INIT_CONFIG | AACENC_INIT_STATES | AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_BITRATE:
      if (settings->userBitrate != value) {
        settings->userBitrate = value;
        hAacEncoder->InitFlags |= AACENC_INIT_CONFIG | AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_BITRATEMODE:
      if (settings->userBitrateMode != value) {
        switch (value) {
          case 0:
          case 1:
          case 2:
          case 3:
          case 4:
          case 5:
            settings->userBitrateMode = value;
            hAacEncoder->InitFlags |=
                AACENC_INIT_CONFIG | AACENC_INIT_TRANSPORT;
            break;
          default:
            err = AACENC_INVALID_CONFIG;
            break;
        }
      }
      break;
    case AACENC_SAMPLERATE:
      if (settings->userSamplerate != value) {
        if (!((value == 8000) || (value == 11025) || (value == 12000) ||
              (value == 16000) || (value == 22050) || (value == 24000) ||
              (value == 32000) || (value == 44100) || (value == 48000) ||
              (value == 64000) || (value == 88200) || (value == 96000))) {
          err = AACENC_INVALID_CONFIG;
          break;
        }
        settings->userSamplerate = value;
        hAacEncoder->nSamplesRead = 0;
        hAacEncoder->InitFlags |=
            AACENC_INIT_CONFIG | AACENC_INIT_STATES | AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_CHANNELMODE:
      if (settings->userChannelMode != (CHANNEL_MODE)value) {
        if (((CHANNEL_MODE)value == MODE_212) &&
            (NULL != hAacEncoder->hMpsEnc)) {
          settings->userChannelMode = (CHANNEL_MODE)value;
          settings->nChannels = 2;
        } else {
          const CHANNEL_MODE_CONFIG_TAB *pConfig =
              FDKaacEnc_GetChannelModeConfiguration((CHANNEL_MODE)value);
          if (pConfig == NULL) {
            err = AACENC_INVALID_CONFIG;
            break;
          }
          if ((pConfig->nElements > hAacEncoder->nMaxAacElements) ||
              (pConfig->nChannelsEff > hAacEncoder->nMaxAacChannels)) {
            err = AACENC_INVALID_CONFIG;
            break;
          }

          settings->userChannelMode = (CHANNEL_MODE)value;
          settings->nChannels = pConfig->nChannels;
        }
        hAacEncoder->nSamplesRead = 0;
        hAacEncoder->InitFlags |= AACENC_INIT_CONFIG | AACENC_INIT_TRANSPORT;
        if (!((value >= 1) && (value <= 6))) {
          hAacEncoder->InitFlags |= AACENC_INIT_STATES;
        }
      }
      break;
    case AACENC_BANDWIDTH:
      if (settings->userBandwidth != value) {
        settings->userBandwidth = value;
        hAacEncoder->InitFlags |= AACENC_INIT_CONFIG;
      }
      break;
    case AACENC_CHANNELORDER:
      if (hAacEncoder->aacConfig.channelOrder != (CHANNEL_ORDER)value) {
        if (!((value == 0) || (value == 1) || (value == 2))) {
          err = AACENC_INVALID_CONFIG;
          break;
        }
        hAacEncoder->aacConfig.channelOrder = (CHANNEL_ORDER)value;
        hAacEncoder->nSamplesRead = 0;
        hAacEncoder->InitFlags |=
            AACENC_INIT_CONFIG | AACENC_INIT_STATES | AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_AFTERBURNER:
      if (settings->userAfterburner != value) {
        if (!((value == 0) || (value == 1))) {
          err = AACENC_INVALID_CONFIG;
          break;
        }
        settings->userAfterburner = value;
        hAacEncoder->InitFlags |= AACENC_INIT_CONFIG;
      }
      break;
    case AACENC_GRANULE_LENGTH:
      if (settings->userFramelength != value) {
        switch (value) {
          case 1024:
          case 512:
          case 480:
          case 256:
          case 240:
          case 128:
          case 120:
            if ((value << 1) == 480 || (value << 1) == 512) {
              settings->userDownscaleFactor = 2;
            } else if ((value << 2) == 480 || (value << 2) == 512) {
              settings->userDownscaleFactor = 4;
            }
            settings->userFramelength = value;
            hAacEncoder->InitFlags |=
                AACENC_INIT_CONFIG | AACENC_INIT_TRANSPORT;
            break;
          default:
            err = AACENC_INVALID_CONFIG;
            break;
        }
      }
      break;
    case AACENC_SBR_RATIO:
      if (settings->userSbrRatio != value) {
        if (!((value == 0) || (value == 1) || (value == 2))) {
          err = AACENC_INVALID_CONFIG;
          break;
        }
        settings->userSbrRatio = value;
        hAacEncoder->InitFlags |=
            AACENC_INIT_CONFIG | AACENC_INIT_STATES | AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_SBR_MODE:
      if ((settings->userSbrEnabled != value) &&
          (NULL != hAacEncoder->hEnvEnc)) {
        settings->userSbrEnabled = value;
        hAacEncoder->InitFlags |=
            AACENC_INIT_CONFIG | AACENC_INIT_STATES | AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_TRANSMUX:
      if (settings->userTpType != (TRANSPORT_TYPE)value) {
        TRANSPORT_TYPE type = (TRANSPORT_TYPE)value;
        UINT flags = hAacEncoder->CAPF_tpEnc;

        if (!(((type == TT_MP4_ADIF) && (flags & CAPF_ADIF)) ||
              ((type == TT_MP4_ADTS) && (flags & CAPF_ADTS)) ||
              ((type == TT_MP4_LATM_MCP0) &&
               ((flags & CAPF_LATM) && (flags & CAPF_RAWPACKETS))) ||
              ((type == TT_MP4_LATM_MCP1) &&
               ((flags & CAPF_LATM) && (flags & CAPF_RAWPACKETS))) ||
              ((type == TT_MP4_LOAS) && (flags & CAPF_LOAS)) ||
              ((type == TT_MP4_RAW) && (flags & CAPF_RAWPACKETS)))) {
          err = AACENC_INVALID_CONFIG;
          break;
        }
        settings->userTpType = (TRANSPORT_TYPE)value;
        hAacEncoder->InitFlags |= AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_SIGNALING_MODE:
      if (settings->userTpSignaling != value) {
        if (!((value == 0) || (value == 1) || (value == 2))) {
          err = AACENC_INVALID_CONFIG;
          break;
        }
        settings->userTpSignaling = value;
        hAacEncoder->InitFlags |= AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_PROTECTION:
      if (settings->userTpProtection != value) {
        if (!((value == 0) || (value == 1))) {
          err = AACENC_INVALID_CONFIG;
          break;
        }
        settings->userTpProtection = value;
        hAacEncoder->InitFlags |= AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_HEADER_PERIOD:
      if (settings->userTpHeaderPeriod != value) {
        if (!(((INT)value >= 0) && (value <= 255))) {
          err = AACENC_INVALID_CONFIG;
          break;
        }
        settings->userTpHeaderPeriod = value;
        hAacEncoder->InitFlags |= AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_AUDIOMUXVER:
      if (settings->userTpAmxv != value) {
        if (!((value == 0) || (value == 1) || (value == 2))) {
          err = AACENC_INVALID_CONFIG;
          break;
        }
        settings->userTpAmxv = value;
        hAacEncoder->InitFlags |= AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_TPSUBFRAMES:
      if (settings->userTpNsubFrames != value) {
        if (!((value >= 1) && (value <= 4))) {
          err = AACENC_INVALID_CONFIG;
          break;
        }
        settings->userTpNsubFrames = value;
        hAacEncoder->InitFlags |= AACENC_INIT_TRANSPORT;
      }
      break;
    case AACENC_ANCILLARY_BITRATE:
      if (settings->userAncDataRate != value) {
        settings->userAncDataRate = value;
      }
      break;
    case AACENC_CONTROL_STATE:
      if (hAacEncoder->InitFlags != value) {
        if (value & AACENC_RESET_INBUFFER) {
          hAacEncoder->nSamplesRead = 0;
        }
        hAacEncoder->InitFlags = value;
      }
      break;
    case AACENC_METADATA_MODE:
      if ((UINT)settings->userMetaDataMode != value) {
        if (!(((INT)value >= 0) && ((INT)value <= 3))) {
          err = AACENC_INVALID_CONFIG;
          break;
        }
        settings->userMetaDataMode = value;
        hAacEncoder->InitFlags |= AACENC_INIT_CONFIG;
      }
      break;
    case AACENC_PEAK_BITRATE:
      if (settings->userPeakBitrate != value) {
        settings->userPeakBitrate = value;
        hAacEncoder->InitFlags |= AACENC_INIT_CONFIG | AACENC_INIT_TRANSPORT;
      }
      break;
    default:
      err = AACENC_UNSUPPORTED_PARAMETER;
      break;
  }

bail:
  return err;
}

LINKSPEC(UINT) aacEncoder_GetParam(const HANDLE_AACENCODER hAacEncoder,
                         const AACENC_PARAM param) {
  UINT value = 0;
  USER_PARAM *settings = &hAacEncoder->extParam;

  if (hAacEncoder == NULL) {
    goto bail;
  }

  switch (param) {
    case AACENC_AOT:
      value = (UINT)hAacEncoder->aacConfig.audioObjectType;
      break;
    case AACENC_BITRATE:
      switch (hAacEncoder->aacConfig.bitrateMode) {
        case AACENC_BR_MODE_CBR:
          value = (UINT)hAacEncoder->aacConfig.bitRate;
          break;
        default:
          value = (UINT)-1;
      }
      break;
    case AACENC_BITRATEMODE:
      value = (UINT)((hAacEncoder->aacConfig.bitrateMode != AACENC_BR_MODE_FF)
                         ? hAacEncoder->aacConfig.bitrateMode
                         : AACENC_BR_MODE_CBR);
      break;
    case AACENC_SAMPLERATE:
      value = (UINT)hAacEncoder->coderConfig.extSamplingRate;
      break;
    case AACENC_CHANNELMODE:
      if ((MODE_1 == hAacEncoder->aacConfig.channelMode) &&
          (hAacEncoder->aacConfig.syntaxFlags & AC_LD_MPS)) {
        value = MODE_212;
      } else {
        value = (UINT)hAacEncoder->aacConfig.channelMode;
      }
      break;
    case AACENC_BANDWIDTH:
      value = (UINT)hAacEncoder->aacConfig.bandWidth;
      break;
    case AACENC_CHANNELORDER:
      value = (UINT)hAacEncoder->aacConfig.channelOrder;
      break;
    case AACENC_AFTERBURNER:
      value = (UINT)hAacEncoder->aacConfig.useRequant;
      break;
    case AACENC_GRANULE_LENGTH:
      value = (UINT)hAacEncoder->aacConfig.framelength;
      break;
    case AACENC_SBR_RATIO:
      value = isSbrActive(&hAacEncoder->aacConfig)
                  ? hAacEncoder->aacConfig.sbrRatio
                  : 0;
      break;
    case AACENC_SBR_MODE:
      value =
          (UINT)(hAacEncoder->aacConfig.syntaxFlags & AC_SBR_PRESENT) ? 1 : 0;
      break;
    case AACENC_TRANSMUX:
      value = (UINT)settings->userTpType;
      break;
    case AACENC_SIGNALING_MODE:
      value = (UINT)getSbrSignalingMode(
          hAacEncoder->aacConfig.audioObjectType, settings->userTpType,
          settings->userTpSignaling, hAacEncoder->aacConfig.sbrRatio);
      break;
    case AACENC_PROTECTION:
      value = (UINT)settings->userTpProtection;
      break;
    case AACENC_HEADER_PERIOD:
      value = (UINT)hAacEncoder->coderConfig.headerPeriod;
      break;
    case AACENC_AUDIOMUXVER:
      value = (UINT)hAacEncoder->aacConfig.audioMuxVersion;
      break;
    case AACENC_TPSUBFRAMES:
      value = (UINT)settings->userTpNsubFrames;
      break;
    case AACENC_ANCILLARY_BITRATE:
      value = (UINT)hAacEncoder->aacConfig.anc_Rate;
      break;
    case AACENC_CONTROL_STATE:
      value = (UINT)hAacEncoder->InitFlags;
      break;
    case AACENC_METADATA_MODE:
      value = (hAacEncoder->metaDataAllowed == 0)
                  ? 0
                  : (UINT)settings->userMetaDataMode;
      break;
    case AACENC_PEAK_BITRATE:
      value = (UINT)-1;
      if (((INT)hAacEncoder->extParam.userPeakBitrate != -1)) {
        value =
            (UINT)(fMax((INT)hAacEncoder->extParam.userPeakBitrate,
                        hAacEncoder->aacConfig
                            .bitRate));
      }
      break;

    default:
      break;
  }

bail:
  return value;
}

LINKSPEC(AACENC_ERROR) aacEncInfo(const HANDLE_AACENCODER hAacEncoder,
                        AACENC_InfoStruct *pInfo) {
  AACENC_ERROR err = AACENC_OK;

  if ((hAacEncoder == NULL) || (pInfo == NULL)) {
    err = AACENC_INVALID_HANDLE;
    goto bail;
  }

  FDKmemclear(pInfo, sizeof(AACENC_InfoStruct));
  pInfo->confSize = 64;

  pInfo->maxOutBufBytes = ((hAacEncoder->nMaxAacChannels * 6144) + 7) >> 3;
  pInfo->maxAncBytes = hAacEncoder->aacConfig.maxAncBytesPerAU;
  pInfo->inBufFillLevel =
      hAacEncoder->nSamplesRead / hAacEncoder->extParam.nChannels;
  pInfo->inputChannels = hAacEncoder->extParam.nChannels;
  pInfo->frameLength =
      hAacEncoder->nSamplesToRead / hAacEncoder->extParam.nChannels;
  pInfo->nDelay = hAacEncoder->nDelay;
  pInfo->nDelayCore = hAacEncoder->nDelayCore;

  if (aacEncGetConf(hAacEncoder, &pInfo->confSize, &pInfo->confBuf[0]) !=
      AAC_ENC_OK) {
    err = AACENC_INIT_ERROR;
    goto bail;
  }
bail:
  return err;
}
