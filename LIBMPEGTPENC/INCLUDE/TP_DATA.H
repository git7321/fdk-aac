/******************* MPEG transport format encoder library *********************

   Author(s):   Manuel Jander

   Description: MPEG Transport data tables

*******************************************************************************/

#ifndef TP_DATA_H
#define TP_DATA_H

#include "machine_type.h"
#include "FDK_audio.h"
#include "FDK_bitstream.h"

#define TP_USAC_MAX_SPEAKERS (24)

#define TP_USAC_MAX_EXT_ELEMENTS ((24))

#define TP_USAC_MAX_ELEMENTS ((24) + TP_USAC_MAX_EXT_ELEMENTS)

#define TP_USAC_MAX_CONFIG_LEN \
  512

#define TPDEC_USAC_NUM_CONFIG_CHANGE_FRAMES \
  (1)

enum {
  TPDEC_FLUSH_OFF = 0,
  TPDEC_RSV60_CFG_CHANGE_ATSC_FLUSH_ON = 1,
  TPDEC_RSV60_DASH_IPF_ATSC_FLUSH_ON = 2,
  TPDEC_USAC_DASH_IPF_FLUSH_ON = 3
};

enum {
  TPDEC_BUILD_UP_OFF = 0,
  TPDEC_RSV60_BUILD_UP_ON = 1,
  TPDEC_RSV60_BUILD_UP_ON_IN_BAND = 2,
  TPDEC_USAC_BUILD_UP_ON = 3,
  TPDEC_RSV60_BUILD_UP_IDLE = 4,
  TPDEC_RSV60_BUILD_UP_IDLE_IN_BAND = 5
};

#define PC_FSB_CHANNELS_MAX 16
#define PC_LFE_CHANNELS_MAX 4
#define PC_ASSOCDATA_MAX 8
#define PC_CCEL_MAX 16
#define PC_COMMENTLENGTH 256
#define PC_NUM_HEIGHT_LAYER 3

typedef struct {
  UCHAR ElementInstanceTag;
  UCHAR Profile;
  UCHAR SamplingFrequencyIndex;
  UCHAR NumFrontChannelElements;
  UCHAR NumSideChannelElements;
  UCHAR NumBackChannelElements;
  UCHAR NumLfeChannelElements;
  UCHAR NumAssocDataElements;
  UCHAR NumValidCcElements;

  UCHAR MonoMixdownPresent;
  UCHAR MonoMixdownElementNumber;

  UCHAR StereoMixdownPresent;
  UCHAR StereoMixdownElementNumber;

  UCHAR MatrixMixdownIndexPresent;
  UCHAR MatrixMixdownIndex;
  UCHAR PseudoSurroundEnable;

  UCHAR FrontElementIsCpe[PC_FSB_CHANNELS_MAX];
  UCHAR FrontElementTagSelect[PC_FSB_CHANNELS_MAX];
  UCHAR FrontElementHeightInfo[PC_FSB_CHANNELS_MAX];

  UCHAR SideElementIsCpe[PC_FSB_CHANNELS_MAX];
  UCHAR SideElementTagSelect[PC_FSB_CHANNELS_MAX];
  UCHAR SideElementHeightInfo[PC_FSB_CHANNELS_MAX];

  UCHAR BackElementIsCpe[PC_FSB_CHANNELS_MAX];
  UCHAR BackElementTagSelect[PC_FSB_CHANNELS_MAX];
  UCHAR BackElementHeightInfo[PC_FSB_CHANNELS_MAX];

  UCHAR LfeElementTagSelect[PC_LFE_CHANNELS_MAX];

  UCHAR AssocDataElementTagSelect[PC_ASSOCDATA_MAX];

  UCHAR CcElementIsIndSw[PC_CCEL_MAX];
  UCHAR ValidCcElementTagSelect[PC_CCEL_MAX];

  UCHAR CommentFieldBytes;
  UCHAR Comment[PC_COMMENTLENGTH];

  UCHAR isValid;
  UCHAR
  NumChannels;
  UCHAR NumEffectiveChannels;
  UCHAR elCounter;
} CProgramConfig;

typedef enum {
  ASCEXT_UNKOWN = -1,
  ASCEXT_SBR = 0x2b7,
  ASCEXT_PS = 0x548,
  ASCEXT_MPS = 0x76a,
  ASCEXT_SAOC = 0x7cb,
  ASCEXT_LDMPS = 0x7cc

} TP_ASC_EXTENSION_ID;

typedef struct {
  UINT m_frameLengthFlag;
  UINT m_dependsOnCoreCoder;
  UINT m_coreCoderDelay;

  UINT m_extensionFlag;
  UINT m_extensionFlag3;

  UINT m_layer;
  UINT m_numOfSubFrame;
  UINT m_layerLength;

} CSGaSpecificConfig;

typedef enum {
  ELDEXT_TERM = 0x0,
  ELDEXT_SAOC = 0x1,
  ELDEXT_LDSAC = 0x2,
  ELDEXT_DOWNSCALEINFO = 0x3
} ASC_ELD_EXT_TYPE;

typedef struct {
  UCHAR m_frameLengthFlag;

  UCHAR m_sbrPresentFlag;
  UCHAR
  m_useLdQmfTimeAlign;
  UCHAR m_sbrSamplingRate;
  UCHAR m_sbrCrcFlag;
  UINT m_downscaledSamplingFrequency;

} CSEldSpecificConfig;

typedef struct {
  USAC_EXT_ELEMENT_TYPE usacExtElementType;
  USHORT usacExtElementConfigLength;
  USHORT usacExtElementDefaultLength;
  UCHAR usacExtElementPayloadFrag;
  UCHAR usacExtElementHasAudioPreRoll;
} CSUsacExtElementConfig;

typedef struct {
  MP4_ELEMENT_ID usacElementType;
  UCHAR m_noiseFilling;
  UCHAR m_harmonicSBR;
  UCHAR m_interTes;
  UCHAR m_pvc;
  UCHAR m_stereoConfigIndex;
  CSUsacExtElementConfig extElement;
} CSUsacElementConfig;

typedef struct {
  UCHAR m_frameLengthFlag;
  UCHAR m_coreSbrFrameLengthIndex;
  UCHAR m_sbrRatioIndex;
  UCHAR m_nUsacChannels;
  UCHAR m_channelConfigurationIndex;
  UINT m_usacNumElements;
  CSUsacElementConfig element[TP_USAC_MAX_ELEMENTS];

  UCHAR numAudioChannels;
  UCHAR m_usacConfigExtensionPresent;
  UCHAR elementLengthPresent;
  UCHAR UsacConfig[TP_USAC_MAX_CONFIG_LEN];
  USHORT UsacConfigBits;
} CSUsacConfig;

typedef struct {
  union {
    CSGaSpecificConfig
        m_gaSpecificConfig;
    CSEldSpecificConfig m_eldSpecificConfig;
    CSUsacConfig m_usacConfig;
  } m_sc;

  CProgramConfig m_progrConfigElement;

  AUDIO_OBJECT_TYPE m_aot;
  UINT m_samplingFrequency;
  UINT m_samplesPerFrame;
  UINT m_directMapping;

  AUDIO_OBJECT_TYPE m_extensionAudioObjectType;
  UINT m_extensionSamplingFrequency;

  SCHAR m_channelConfiguration;

  SCHAR m_epConfig;
  SCHAR m_vcb11Flag;
  SCHAR m_rvlcFlag;
  SCHAR m_hcrFlag;

  SCHAR m_sbrPresentFlag;
  SCHAR
  m_psPresentFlag;
  UCHAR m_samplingFrequencyIndex;
  UCHAR m_extensionSamplingFrequencyIndex;
  SCHAR m_extensionChannelConfiguration;

  UCHAR
  configMode;
  UCHAR AacConfigChanged;
  UCHAR SbrConfigChanged;
  UCHAR SacConfigChanged;

  UCHAR
  config[TP_USAC_MAX_CONFIG_LEN];
  UINT configBits;

} CSAudioSpecificConfig;

typedef struct {
  SCHAR flushCnt;
  UCHAR flushStatus;
  SCHAR buildUpCnt;
  UCHAR buildUpStatus;
  UCHAR cfgChanged;
  UCHAR contentChanged;
  UCHAR forceCfgChange;
} CCtrlCFGChange;

typedef INT (*cbUpdateConfig_t)(void *, const CSAudioSpecificConfig *,
                                const UCHAR configMode, UCHAR *configChanged);
typedef INT (*cbFreeMem_t)(void *, const CSAudioSpecificConfig *);
typedef INT (*cbCtrlCFGChange_t)(void *, const CCtrlCFGChange *);
typedef INT (*cbSsc_t)(void *, HANDLE_FDK_BITSTREAM,
                       const AUDIO_OBJECT_TYPE coreCodec,
                       const INT samplingRate, const INT frameSize,
                       const INT numChannels, const INT stereoConfigIndex,
                       const INT coreSbrFrameLengthIndex, const INT configBytes,
                       const UCHAR configMode, UCHAR *configChanged);

typedef INT (*cbSbr_t)(void *self, HANDLE_FDK_BITSTREAM hBs,
                       const INT sampleRateIn, const INT sampleRateOut,
                       const INT samplesPerFrame,
                       const AUDIO_OBJECT_TYPE coreCodec,
                       const MP4_ELEMENT_ID elementID, const INT elementIndex,
                       const UCHAR harmonicSbr, const UCHAR stereoConfigIndex,
                       const UCHAR configMode, UCHAR *configChanged,
                       const INT downscaleFactor);

typedef INT (*cbUsac_t)(void *self, HANDLE_FDK_BITSTREAM hBs);

typedef INT (*cbUniDrc_t)(void *self, HANDLE_FDK_BITSTREAM hBs,
                          const INT fullPayloadLength, const INT payloadType,
                          const INT subStreamIndex, const INT payloadStart,
                          const AUDIO_OBJECT_TYPE);

typedef struct {
  cbUpdateConfig_t cbUpdateConfig;
  void *cbUpdateConfigData;
  cbFreeMem_t cbFreeMem;
  void *cbFreeMemData;
  cbCtrlCFGChange_t cbCtrlCFGChange;
  void *cbCtrlCFGChangeData;
  cbSsc_t cbSsc;
  void *cbSscData;
  cbSbr_t cbSbr;
  void *cbSbrData;
  cbUsac_t cbUsac;
  void *cbUsacData;
  cbUniDrc_t cbUniDrc;
  void *cbUniDrcData;
} CSTpCallBacks;

static const UINT SamplingRateTable[] = {
    96000, 88200, 64000, 48000, 44100, 32000, 24000, 22050, 16000, 12000, 11025,
    8000,  7350,  0,     0,     57600, 51200, 40000, 38400, 34150, 28800, 25600,
    20000, 19200, 17075, 14400, 12800, 9600,  0,     0,     0,     0};

static inline int getSamplingRateIndex(UINT samplingRate, UINT nBits) {
  UINT sf_index;
  UINT tableSize = (1 << nBits) - 1;

  for (sf_index = 0; sf_index < tableSize; sf_index++) {
    if (SamplingRateTable[sf_index] == samplingRate) break;
  }

  if (sf_index > tableSize) {
    return tableSize - 1;
  }

  return sf_index;
}

static inline int getNumberOfTotalChannels(int channelConfig) {
  switch (channelConfig) {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
      return channelConfig;
    case 7:
    case 12:
    case 14:
      return 8;
    case 11:
      return 7;
    case 13:
      return 24;
    default:
      return 0;
  }
}

static inline int getNumberOfEffectiveChannels(
    const int
        channelConfig) {
  const int n[] = {0, 1, 2, 3, 4, 5, 5, 7, 0, 0, 0, 6, 7, 22, 7, 0};
  return n[channelConfig];
}

#endif /* TP_DATA_H */
