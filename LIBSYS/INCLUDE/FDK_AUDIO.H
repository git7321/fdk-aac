#ifndef FDK_AUDIO_H
#define FDK_AUDIO_H

#include "machine_type.h"
#include "genericStds.h"
#include "syslib_channelMapDescr.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  FF_UNKNOWN = -1,
  FF_RAW = 0,

  FF_MP4_3GPP = 3,
  FF_MP4_MP4F = 4,

  FF_RAWPACKETS = 5
} FILE_FORMAT;

typedef enum {
  TT_UNKNOWN = -1,
  TT_MP4_RAW = 0,
  TT_MP4_ADIF = 1,
  TT_MP4_ADTS = 2,

  TT_MP4_LATM_MCP1 = 6,
  TT_MP4_LATM_MCP0 = 7,

  TT_MP4_LOAS = 10,

  TT_DRM = 12
} TRANSPORT_TYPE;

#define TT_IS_PACKET(x) \
  (((x) == TT_MP4_RAW) || ((x) == TT_DRM) || ((x) == TT_MP4_LATM_MCP0) || \
   ((x) == TT_MP4_LATM_MCP1))

typedef enum {
  AOT_NONE = -1,
  AOT_NULL_OBJECT = 0,
  AOT_AAC_MAIN = 1,
  AOT_AAC_LC = 2,
  AOT_AAC_SSR = 3,
  AOT_AAC_LTP = 4,
  AOT_SBR = 5,
  AOT_AAC_SCAL = 6,
  AOT_TWIN_VQ = 7,
  AOT_CELP = 8,
  AOT_HVXC = 9,
  AOT_RSVD_10 = 10,
  AOT_RSVD_11 = 11,
  AOT_TTSI = 12,
  AOT_MAIN_SYNTH = 13,
  AOT_WAV_TAB_SYNTH = 14,
  AOT_GEN_MIDI = 15,
  AOT_ALG_SYNTH_AUD_FX = 16,
  AOT_ER_AAC_LC = 17,
  AOT_RSVD_18 = 18,
  AOT_ER_AAC_LTP = 19,
  AOT_ER_AAC_SCAL = 20,
  AOT_ER_TWIN_VQ = 21,
  AOT_ER_BSAC = 22,
  AOT_ER_AAC_LD = 23,
  AOT_ER_CELP = 24,
  AOT_ER_HVXC = 25,
  AOT_ER_HILN = 26,
  AOT_ER_PARA = 27,
  AOT_RSVD_28 = 28,
  AOT_PS = 29,
  AOT_MPEGS = 30,

  AOT_ESCAPE = 31,

  AOT_MP3ONMP4_L1 = 32,
  AOT_MP3ONMP4_L2 = 33,
  AOT_MP3ONMP4_L3 = 34,
  AOT_RSVD_35 = 35,
  AOT_RSVD_36 = 36,
  AOT_AAC_SLS = 37,
  AOT_SLS = 38,
  AOT_ER_AAC_ELD = 39,

  AOT_USAC = 42,
  AOT_SAOC = 43,
  AOT_LD_MPEGS = 44,

  AOT_MP2_AAC_LC = 129,
  AOT_MP2_SBR = 132,

  AOT_DRM_AAC = 143,
  AOT_DRM_SBR = 144,
  AOT_DRM_MPEG_PS =
      145,
  AOT_DRM_SURROUND =
      146,
  AOT_DRM_USAC = 147
} AUDIO_OBJECT_TYPE;

#define CAN_DO_PS(aot) \
  ((aot) == AOT_AAC_LC || (aot) == AOT_SBR || (aot) == AOT_PS || \
   (aot) == AOT_ER_BSAC || (aot) == AOT_DRM_AAC)

#define IS_USAC(aot) ((aot) == AOT_USAC)

#define IS_LOWDELAY(aot) ((aot) == AOT_ER_AAC_LD || (aot) == AOT_ER_AAC_ELD)

typedef enum {
  MODE_INVALID = -1,
  MODE_UNKNOWN = 0,
  MODE_1 = 1,
  MODE_2 = 2,
  MODE_1_2 = 3,
  MODE_1_2_1 = 4,
  MODE_1_2_2 = 5,
  MODE_1_2_2_1 = 6,
  MODE_1_2_2_2_1 = 7,

  MODE_6_1 = 11,
  MODE_7_1_BACK = 12,
  MODE_7_1_TOP_FRONT = 14,

  MODE_7_1_REAR_SURROUND = 33,
  MODE_7_1_FRONT_CENTER = 34,

  MODE_212 = 128
} CHANNEL_MODE;

typedef enum {
  ACT_NONE = 0x00,
  ACT_FRONT = 0x01,
  ACT_SIDE = 0x02,
  ACT_BACK = 0x03,
  ACT_LFE = 0x04,

  ACT_TOP =
      0x10,
  ACT_FRONT_TOP = 0x11,
  ACT_SIDE_TOP = 0x12,
  ACT_BACK_TOP = 0x13,

  ACT_BOTTOM =
      0x20,
  ACT_FRONT_BOTTOM = 0x21,
  ACT_SIDE_BOTTOM = 0x22,
  ACT_BACK_BOTTOM = 0x23

} AUDIO_CHANNEL_TYPE;

typedef enum {
  SIG_UNKNOWN = -1,
  SIG_IMPLICIT = 0,
  SIG_EXPLICIT_BW_COMPATIBLE = 1,
  SIG_EXPLICIT_HIERARCHICAL = 2

} SBR_PS_SIGNALING;

#define AC_ER_VCB11 \
  0x000001
#define AC_ER_RVLC \
  0x000002
#define AC_ER_HCR \
  0x000004
#define AC_SCALABLE 0x000008
#define AC_ELD 0x000010
#define AC_LD 0x000020
#define AC_ER 0x000040
#define AC_BSAC 0x000080
#define AC_USAC 0x000100
#define AC_RSV603DA 0x000200
#define AC_HDAAC 0x000400
#define AC_RSVD50 0x004000
#define AC_SBR_PRESENT 0x008000
#define AC_SBRCRC \
  0x010000
#define AC_PS_PRESENT 0x020000
#define AC_MPS_PRESENT \
  0x040000
#define AC_DRM 0x080000
#define AC_INDEP 0x100000
#define AC_MPEGD_RES 0x200000
#define AC_SAOC_PRESENT 0x400000
#define AC_DAB 0x800000
#define AC_ELD_DOWNSCALE 0x1000000
#define AC_LD_MPS 0x2000000
#define AC_DRC_PRESENT \
  0x4000000
#define AC_USAC_SCFGI3 \
  0x8000000
#define AC_CM_DET_CFG_CHANGE \
  0x000001
#define AC_CM_ALLOC_MEM \
  0x000002

#define AC_EL_USAC_TW 0x000001
#define AC_EL_USAC_NOISE 0x000002
#define AC_EL_USAC_ITES 0x000004
#define AC_EL_USAC_PVC \
  0x000008
#define AC_EL_USAC_MPS212 0x000010
#define AC_EL_USAC_LFE 0x000020
#define AC_EL_USAC_CP_POSSIBLE \
  0x000040
#define AC_EL_ENHANCED_NOISE 0x000080
#define AC_EL_IGF_AFTER_TNS 0x000100
#define AC_EL_IGF_INDEP_TILING 0x000200
#define AC_EL_IGF_USE_ENF 0x000400
#define AC_EL_FULLBANDLPD 0x000800
#define AC_EL_LPDSTEREOIDX 0x001000
#define AC_EL_LFE 0x002000

#define CC_MPEG_ID 0x00100000
#define CC_IS_BASELAYER 0x00200000
#define CC_PROTECTION 0x00400000
#define CC_SBR 0x00800000
#define CC_SBRCRC 0x00010000
#define CC_SAC 0x00020000
#define CC_RVLC 0x01000000
#define CC_VCB11 0x02000000
#define CC_HCR 0x04000000
#define CC_PSEUDO_SURROUND 0x08000000
#define CC_USAC_NOISE 0x10000000
#define CC_USAC_TW 0x20000000
#define CC_USAC_HBE 0x40000000

typedef struct {
  AUDIO_OBJECT_TYPE aot;
  AUDIO_OBJECT_TYPE extAOT;
  CHANNEL_MODE channelMode;
  UCHAR channelConfigZero;
  INT samplingRate;
  INT extSamplingRate;
  INT downscaleSamplingRate;
  INT bitRate;
  int samplesPerFrame;
  int noChannels;
  int bitsFrame;
  int nSubFrames;
  int BSACnumOfSubFrame;
  int BSAClayerLength;
  UINT flags;
  UCHAR matrixMixdownA;
  UCHAR headerPeriod;

  UCHAR stereoConfigIndex;
  UCHAR sbrMode;
  SBR_PS_SIGNALING sbrSignaling;

  UCHAR rawConfig[64];
  int rawConfigBits;

  UCHAR sbrPresent;
  UCHAR psPresent;
} CODER_CONFIG;

#define USAC_ID_BIT 16

typedef enum {
  ID_NONE = -1,
  ID_SCE = 0,
  ID_CPE = 1,
  ID_CCE = 2,
  ID_LFE = 3,
  ID_DSE = 4,
  ID_PCE = 5,
  ID_FIL = 6,
  ID_END = 7,
  ID_EXT = 8,
  ID_SCAL = 9,
  ID_USAC_SCE = 0 + USAC_ID_BIT,
  ID_USAC_CPE = 1 + USAC_ID_BIT,
  ID_USAC_LFE = 2 + USAC_ID_BIT,
  ID_USAC_EXT = 3 + USAC_ID_BIT,
  ID_USAC_END = 4 + USAC_ID_BIT,
  ID_LAST
} MP4_ELEMENT_ID;

typedef enum {
  ID_CONFIG_EXT_FILL = 0,
  ID_CONFIG_EXT_DOWNMIX = 1,
  ID_CONFIG_EXT_LOUDNESS_INFO = 2,
  ID_CONFIG_EXT_AUDIOSCENE_INFO = 3,
  ID_CONFIG_EXT_HOA_MATRIX = 4,
  ID_CONFIG_EXT_SIG_GROUP_INFO = 6
} CONFIG_EXT_ID;

#define IS_CHANNEL_ELEMENT(elementId) \
  ((elementId) == ID_SCE || (elementId) == ID_CPE || (elementId) == ID_LFE || \
   (elementId) == ID_USAC_SCE || (elementId) == ID_USAC_CPE || \
   (elementId) == ID_USAC_LFE)

#define IS_MP4_CHANNEL_ELEMENT(elementId) \
  ((elementId) == ID_SCE || (elementId) == ID_CPE || (elementId) == ID_LFE)

#define EXT_ID_BITS 4

typedef enum {
  EXT_FIL = 0x00,
  EXT_FILL_DATA = 0x01,
  EXT_DATA_ELEMENT = 0x02,
  EXT_DATA_LENGTH = 0x03,
  EXT_UNI_DRC = 0x04,
  EXT_LDSAC_DATA = 0x09,
  EXT_SAOC_DATA = 0x0a,
  EXT_DYNAMIC_RANGE = 0x0b,
  EXT_SAC_DATA = 0x0c,
  EXT_SBR_DATA = 0x0d,
  EXT_SBR_DATA_CRC = 0x0e
} EXT_PAYLOAD_TYPE;

#define IS_USAC_CHANNEL_ELEMENT(elementId) \
  ((elementId) == ID_USAC_SCE || (elementId) == ID_USAC_CPE || \
   (elementId) == ID_USAC_LFE)

typedef enum {
  ID_EXT_ELE_FILL = 0x00,
  ID_EXT_ELE_MPEGS = 0x01,
  ID_EXT_ELE_SAOC = 0x02,
  ID_EXT_ELE_AUDIOPREROLL = 0x03,
  ID_EXT_ELE_UNI_DRC = 0x04,
  ID_EXT_ELE_OBJ_METADATA = 0x05,
  ID_EXT_ELE_SAOC_3D = 0x06,
  ID_EXT_ELE_HOA = 0x07,
  ID_EXT_ELE_FMT_CNVRTR = 0x08,
  ID_EXT_ELE_MCT = 0x09,
  ID_EXT_ELE_ENHANCED_OBJ_METADATA = 0x0d,
  ID_EXT_ELE_VR_METADATA = 0x81,
  ID_EXT_ELE_UNKNOWN = 0xFF
} USAC_EXT_ELEMENT_TYPE;

typedef enum {
  TC_NOTHING = 0,
  TC_RAW_ADTS = 2,
  TC_RAW_LATM_MCP1 = 6,
  TC_RAW_SDC = 21

} TP_CONFIG_TYPE;

typedef enum {
  FDK_NONE = 0,
  FDK_TOOLS = 1,
  FDK_SYSLIB = 2,
  FDK_AACDEC = 3,
  FDK_AACENC = 4,
  FDK_SBRDEC = 5,
  FDK_SBRENC = 6,
  FDK_TPDEC = 7,
  FDK_TPENC = 8,
  FDK_MPSDEC = 9,
  FDK_MPEGFILEREAD = 10,
  FDK_MPEGFILEWRITE = 11,
  FDK_PCMDMX = 31,
  FDK_MPSENC = 34,
  FDK_TDLIMIT = 35,
  FDK_UNIDRCDEC = 38,

  FDK_MODULE_LAST

} FDK_MODULE_ID;

#define CAPF_AAC_LC 0x00000001
#define CAPF_ER_AAC_LD \
  0x00000002
#define CAPF_ER_AAC_SCAL 0x00000004
#define CAPF_ER_AAC_LC \
  0x00000008
#define CAPF_AAC_480 \
  0x00000010
#define CAPF_AAC_512 \
  0x00000020
#define CAPF_AAC_960 \
  0x00000040
#define CAPF_AAC_1024 \
  0x00000080
#define CAPF_AAC_HCR \
  0x00000100
#define CAPF_AAC_VCB11 \
  0x00000200
#define CAPF_AAC_RVLC \
  0x00000400
#define CAPF_AAC_MPEG4 0x00000800
#define CAPF_AAC_DRC \
  0x00001000
#define CAPF_AAC_CONCEALMENT \
  0x00002000
#define CAPF_AAC_DRM_BSFORMAT \
  0x00004000
#define CAPF_ER_AAC_ELD \
  0x00008000
#define CAPF_ER_AAC_BSAC \
  0x00010000
#define CAPF_AAC_ELD_DOWNSCALE \
  0x00040000
#define CAPF_AAC_USAC_LP \
  0x00100000
#define CAPF_AAC_USAC \
  0x00200000
#define CAPF_ER_AAC_ELDV2 \
  0x00800000
#define CAPF_AAC_UNIDRC \
  0x01000000

#define CAPF_ADTS \
  0x00000001
#define CAPF_ADIF \
  0x00000002
#define CAPF_LATM \
  0x00000004
#define CAPF_LOAS \
  0x00000008
#define CAPF_RAWPACKETS \
  0x00000010
#define CAPF_DRM \
  0x00000020
#define CAPF_RSVD50 \
  0x00000040
#define CAPF_SBR_LP \
  0x00000001
#define CAPF_SBR_HQ \
  0x00000002
#define CAPF_SBR_DRM_BS \
  0x00000004
#define CAPF_SBR_CONCEALMENT \
  0x00000008
#define CAPF_SBR_DRC \
  0x00000010
#define CAPF_SBR_PS_MPEG \
  0x00000020
#define CAPF_SBR_PS_DRM \
  0x00000040
#define CAPF_SBR_ELD_DOWNSCALE \
  0x00000080
#define CAPF_SBR_HBEHQ \
  0x00000100

#define CAPF_DMX_BLIND \
  0x00000001
#define CAPF_DMX_PCE \
  0x00000002
#define CAPF_DMX_ARIB \
  0x00000004
#define CAPF_DMX_DVB \
  0x00000008
#define CAPF_DMX_CH_EXP \
  0x00000010
#define CAPF_DMX_6_CH \
  0x00000020
#define CAPF_DMX_8_CH \
  0x00000040
#define CAPF_DMX_24_CH \
  0x00000080
#define CAPF_LIMITER \
  0x00002000

#define CAPF_MPS_STD \
  0x00000001
#define CAPF_MPS_LD \
  0x00000002
#define CAPF_MPS_USAC \
  0x00000004
#define CAPF_MPS_HQ \
  0x00000010
#define CAPF_MPS_LP \
  0x00000020
#define CAPF_MPS_BLIND \
  0x00000040
#define CAPF_MPS_BINAURAL \
  0x00000080
#define CAPF_MPS_2CH_OUT \
  0x00000100
#define CAPF_MPS_6CH_OUT \
  0x00000200
#define CAPF_MPS_8CH_OUT \
  0x00000400
#define CAPF_MPS_1CH_IN \
  0x00001000
#define CAPF_MPS_2CH_IN \
  0x00002000
#define CAPF_MPS_6CH_IN \
  0x00004000

#define LIB_VERSION(lev0, lev1, lev2)                      \
  ((lev0 << 24 & 0xff000000) | (lev1 << 16 & 0x00ff0000) | \
   (lev2 << 8 & 0x0000ff00))

#define LIB_VERSION_STRING(info)                                               \
  FDKsprintf((info)->versionStr, "%d.%d.%d", (((info)->version >> 24) & 0xff), \
             (((info)->version >> 16) & 0xff),                                 \
             (((info)->version >> 8) & 0xff))

typedef struct LIB_INFO {
  const char* title;
  const char* build_date;
  const char* build_time;
  FDK_MODULE_ID module_id;
  INT version;
  UINT flags;
  char versionStr[32];
} LIB_INFO;

#ifdef __cplusplus
#define FDK_AUDIO_INLINE inline
#else
#define FDK_AUDIO_INLINE
#endif

static FDK_AUDIO_INLINE void FDKinitLibInfo(LIB_INFO* info) {
  int i;

  for (i = 0; i < FDK_MODULE_LAST; i++) {
    info[i].module_id = FDK_NONE;
  }
}

static FDK_AUDIO_INLINE UINT
FDKlibInfo_getCapabilities(const LIB_INFO* info, FDK_MODULE_ID module_id) {
  int i;

  for (i = 0; i < FDK_MODULE_LAST; i++) {
    if (info[i].module_id == module_id) {
      return info[i].flags;
    }
  }
  return 0;
}

static FDK_AUDIO_INLINE INT FDKlibInfo_lookup(const LIB_INFO* info,
                                              FDK_MODULE_ID module_id) {
  int i = -1;

  for (i = 0; i < FDK_MODULE_LAST; i++) {
    if (info[i].module_id == module_id) return -1;
    if (info[i].module_id == FDK_NONE) break;
  }
  if (i == FDK_MODULE_LAST) return -1;

  return i;
}

typedef struct FDK_bufDescr {
  void** ppBase;
  UINT* pBufSize;
  UINT* pEleSize;
  UINT* pBufType;
  UINT numBufs;

} FDK_bufDescr;

#define FDK_BUF_TYPE_MASK_IO ((UINT)0x03 << 30)
#define FDK_BUF_TYPE_MASK_DESCR ((UINT)0x3F << 16)
#define FDK_BUF_TYPE_MASK_ID ((UINT)0xFF)

#define FDK_BUF_TYPE_INPUT ((UINT)0x1 << 30)
#define FDK_BUF_TYPE_OUTPUT ((UINT)0x2 << 30)

#define FDK_BUF_TYPE_PCM_DATA ((UINT)0x1 << 16)
#define FDK_BUF_TYPE_ANC_DATA ((UINT)0x2 << 16)
#define FDK_BUF_TYPE_BS_DATA ((UINT)0x4 << 16)

#ifdef __cplusplus
}
#endif

#endif /* FDK_AUDIO_H */
