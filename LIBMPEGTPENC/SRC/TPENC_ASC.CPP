/******************* MPEG transport format encoder library *********************

   Author(s):

   Description:

*******************************************************************************/

#include "tp_data.h"
#include "tpenc_lib.h"
#include "tpenc_asc.h"
#include "FDK_bitstream.h"
#include "genericStds.h"
#include "FDK_crc.h"

#define PCE_HEIGHT_EXT_SYNC (0xAC)
#define HEIGHT_NORMAL 0
#define HEIGHT_TOP 1
#define HEIGHT_BOTTOM 2
#define MAX_FRONT_ELEMENTS 8
#define MAX_SIDE_ELEMENTS 3
#define MAX_BACK_ELEMENTS 4

typedef struct {
  UCHAR
  num_front_height_channel_elements[2];
  UCHAR num_side_height_channel_elements[2];
  UCHAR num_back_height_channel_elements[2];
} PCE_HEIGHT_NUM;

typedef struct {
  UCHAR num_front_channel_elements;
  UCHAR num_side_channel_elements;
  UCHAR num_back_channel_elements;
  UCHAR num_lfe_channel_elements;
  const MP4_ELEMENT_ID
      *pEl_type;
  const PCE_HEIGHT_NUM *pHeight_num;
} PCE_CONFIGURATION;

typedef struct {
  CHANNEL_MODE channel_mode;
  PCE_CONFIGURATION
  pce_configuration;

} CHANNEL_CONFIGURATION;

static const MP4_ELEMENT_ID elType_1[] = {ID_SCE};
static const MP4_ELEMENT_ID elType_2[] = {ID_CPE};
static const MP4_ELEMENT_ID elType_1_2[] = {ID_SCE, ID_CPE};
static const MP4_ELEMENT_ID elType_1_2_1[] = {ID_SCE, ID_CPE, ID_SCE};
static const MP4_ELEMENT_ID elType_1_2_2[] = {ID_SCE, ID_CPE, ID_CPE};
static const MP4_ELEMENT_ID elType_1_2_2_1[] = {ID_SCE, ID_CPE, ID_CPE, ID_LFE};
static const MP4_ELEMENT_ID elType_1_2_2_2_1[] = {ID_SCE, ID_CPE, ID_CPE,
                                                  ID_CPE, ID_LFE};
static const MP4_ELEMENT_ID elType_6_1[] = {ID_SCE, ID_CPE, ID_CPE, ID_SCE,
                                            ID_LFE};
static const MP4_ELEMENT_ID elType_7_1_back[] = {ID_SCE, ID_CPE, ID_CPE, ID_CPE,
                                                 ID_LFE};
static const MP4_ELEMENT_ID elType_7_1_top_front[] = {ID_SCE, ID_CPE, ID_CPE,
                                                      ID_LFE, ID_CPE};
static const MP4_ELEMENT_ID elType_7_1_rear_surround[] = {
    ID_SCE, ID_CPE, ID_CPE, ID_CPE, ID_LFE};
static const MP4_ELEMENT_ID elType_7_1_front_center[] = {ID_SCE, ID_CPE, ID_CPE,
                                                         ID_CPE, ID_LFE};

static const PCE_HEIGHT_NUM heightNum_7_1_top_front = {{1, 0}, {0, 0}, {0, 0}};

static const CHANNEL_CONFIGURATION pceConfigTab[] = {
    {MODE_1,
     {1, 0, 0, 0, elType_1,
      NULL}},
    {MODE_2,
     {1, 0, 0, 0, elType_2,
      NULL}},
    {MODE_1_2,
     {2, 0, 0, 0, elType_1_2,
      NULL}},
    {MODE_1_2_1,
     {2, 0, 1, 0, elType_1_2_1,
      NULL}},
    {MODE_1_2_2,
     {2, 0, 1, 0, elType_1_2_2,
      NULL}},
    {MODE_1_2_2_1,
     {2, 0, 1, 1, elType_1_2_2_1,
      NULL}},
    {MODE_1_2_2_2_1,
     {3, 0, 1, 1, elType_1_2_2_2_1,
      NULL}},

    {MODE_6_1,
     {2, 0, 2, 1, elType_6_1,
      NULL}},
    {MODE_7_1_BACK,
     {2, 0, 2, 1, elType_7_1_back,
      NULL}},
    {MODE_7_1_TOP_FRONT,
     {3, 0, 1, 1, elType_7_1_top_front, &heightNum_7_1_top_front}},

    {MODE_7_1_REAR_SURROUND,
     {2, 0, 2, 1, elType_7_1_rear_surround,
      NULL}},
    {MODE_7_1_FRONT_CENTER,
     {3, 0, 1, 1, elType_7_1_front_center,
      NULL}}
};

static const PCE_CONFIGURATION *getPceEntry(const CHANNEL_MODE channel_mode) {
  UINT i;
  const PCE_CONFIGURATION *pce_config = NULL;

  for (i = 0; i < (sizeof(pceConfigTab) / sizeof(CHANNEL_CONFIGURATION)); i++) {
    if (pceConfigTab[i].channel_mode == channel_mode) {
      pce_config = &pceConfigTab[i].pce_configuration;
      break;
    }
  }

  return pce_config;
}

int getChannelConfig(const CHANNEL_MODE channel_mode,
                     const UCHAR channel_config_zero) {
  INT chan_config = 0;

  if (channel_config_zero != 0) {
    chan_config = 0;
  } else {
    switch (channel_mode) {
      case MODE_1:
        chan_config = 1;
        break;
      case MODE_2:
        chan_config = 2;
        break;
      case MODE_1_2:
        chan_config = 3;
        break;
      case MODE_1_2_1:
        chan_config = 4;
        break;
      case MODE_1_2_2:
        chan_config = 5;
        break;
      case MODE_1_2_2_1:
        chan_config = 6;
        break;
      case MODE_1_2_2_2_1:
        chan_config = 7;
        break;
      case MODE_6_1:
        chan_config = 11;
        break;
      case MODE_7_1_BACK:
        chan_config = 12;
        break;
      case MODE_7_1_TOP_FRONT:
        chan_config = 14;
        break;
      default:
        chan_config = 0;
    }
  }

  return chan_config;
}

CHANNEL_MODE transportEnc_GetChannelMode(int noChannels) {
  CHANNEL_MODE chMode;

  if (noChannels <= 8 && noChannels > 0)
    chMode = (CHANNEL_MODE)(
        (noChannels == 8) ? 7
                          : noChannels);
  else
    chMode = MODE_UNKNOWN;

  return chMode;
}

int transportEnc_writePCE(HANDLE_FDK_BITSTREAM hBs, CHANNEL_MODE channelMode,
                          INT sampleRate, int instanceTagPCE, int profile,
                          int matrixMixdownA, int pseudoSurroundEnable,
                          UINT alignAnchor) {
  int sampleRateIndex, i;
  const PCE_CONFIGURATION *config = NULL;
  const MP4_ELEMENT_ID *pEl_list = NULL;
  UCHAR cpeCnt = 0, sceCnt = 0, lfeCnt = 0, frntCnt = 0, sdCnt = 0, bckCnt = 0,
        isCpe = 0, tag = 0, normalFrontEnd = 0, normalSideEnd = 0,
        normalBackEnd = 0, topFrontEnd = 0, topSideEnd = 0, topBackEnd = 0,
        bottomFrontEnd = 0, bottomSideEnd = 0;
#ifdef FDK_ASSERT_ENABLE
  UCHAR bottomBackEnd = 0;
#endif
  enum elementDepth { FRONT, SIDE, BACK } elDepth;

  sampleRateIndex = getSamplingRateIndex(sampleRate, 4);
  if (sampleRateIndex == 15) {
    return -1;
  }

  if ((config = getPceEntry(channelMode)) == NULL) {
    return -1;
  }

  FDK_ASSERT(config->num_front_channel_elements <= MAX_FRONT_ELEMENTS);
  FDK_ASSERT(config->num_side_channel_elements <= MAX_SIDE_ELEMENTS);
  FDK_ASSERT(config->num_back_channel_elements <= MAX_BACK_ELEMENTS);

  UCHAR frontIsCpe[MAX_FRONT_ELEMENTS] = {0},
        frontTag[MAX_FRONT_ELEMENTS] = {0}, sideIsCpe[MAX_SIDE_ELEMENTS] = {0},
        sideTag[MAX_SIDE_ELEMENTS] = {0}, backIsCpe[MAX_BACK_ELEMENTS] = {0},
        backTag[MAX_BACK_ELEMENTS] = {0};

  FDKwriteBits(hBs, instanceTagPCE, 4);
  FDKwriteBits(hBs, profile, 2);
  FDKwriteBits(hBs, sampleRateIndex, 4);

  FDKwriteBits(hBs, config->num_front_channel_elements,
               4);
  FDKwriteBits(hBs, config->num_side_channel_elements,
               4);
  FDKwriteBits(hBs, config->num_back_channel_elements,
               4);
  FDKwriteBits(hBs, config->num_lfe_channel_elements,
               2);

  FDKwriteBits(hBs, 0, 3);
  FDKwriteBits(hBs, 0, 4);
  FDKwriteBits(hBs, 0, 1);
  FDKwriteBits(hBs, 0, 1);

  if (matrixMixdownA != 0 &&
      ((channelMode == MODE_1_2_2) || (channelMode == MODE_1_2_2_1))) {
    FDKwriteBits(hBs, 1, 1);
    FDKwriteBits(hBs, (matrixMixdownA - 1) & 0x3, 2);
    FDKwriteBits(hBs, (pseudoSurroundEnable) ? 1 : 0,
                 1);
  } else {
    FDKwriteBits(hBs, 0, 1);
  }

  if (config->pHeight_num != NULL) {
    normalFrontEnd = config->num_front_channel_elements -
                     config->pHeight_num->num_front_height_channel_elements[0] -
                     config->pHeight_num->num_front_height_channel_elements[1];
    normalSideEnd = normalFrontEnd + config->num_side_channel_elements -
                    config->pHeight_num->num_side_height_channel_elements[0] -
                    config->pHeight_num->num_side_height_channel_elements[1];
    normalBackEnd = normalSideEnd + config->num_back_channel_elements -
                    config->pHeight_num->num_back_height_channel_elements[0] -
                    config->pHeight_num->num_back_height_channel_elements[1];

    topFrontEnd =
        normalBackEnd + config->num_lfe_channel_elements +
        config->pHeight_num->num_front_height_channel_elements[0];
    topSideEnd =
        topFrontEnd + config->pHeight_num->num_side_height_channel_elements[0];
    topBackEnd =
        topSideEnd + config->pHeight_num->num_back_height_channel_elements[0];

    bottomFrontEnd =
        topBackEnd + config->pHeight_num->num_front_height_channel_elements[1];
    bottomSideEnd = bottomFrontEnd +
                    config->pHeight_num->num_side_height_channel_elements[1];
#ifdef FDK_ASSERT_ENABLE
    bottomBackEnd = bottomSideEnd +
                    config->pHeight_num->num_back_height_channel_elements[1];
#endif

  } else {
    normalFrontEnd = config->num_front_channel_elements;
    normalSideEnd = normalFrontEnd + config->num_side_channel_elements;
    normalBackEnd = normalSideEnd + config->num_back_channel_elements;
  }

  pEl_list = config->pEl_type;

  for (i = 0; i < config->num_front_channel_elements +
                      config->num_side_channel_elements +
                      config->num_back_channel_elements +
                      config->num_lfe_channel_elements;
       i++) {
    if (*pEl_list == ID_LFE) {
      pEl_list++;
      continue;
    }
    isCpe = (*pEl_list++ == ID_CPE) ? 1 : 0;
    tag = (isCpe) ? cpeCnt++ : sceCnt++;

    if (i < normalFrontEnd)
      elDepth = FRONT;
    else if (i < normalSideEnd)
      elDepth = SIDE;
    else if (i < normalBackEnd)
      elDepth = BACK;
    else if (i < topFrontEnd)
      elDepth = FRONT;
    else if (i < topSideEnd)
      elDepth = SIDE;
    else if (i < topBackEnd)
      elDepth = BACK;
    else if (i < bottomFrontEnd)
      elDepth = FRONT;
    else if (i < bottomSideEnd)
      elDepth = SIDE;
    else {
      elDepth = BACK;
      FDK_ASSERT(i < bottomBackEnd);
    }

    switch (elDepth) {
      case FRONT:
        FDK_ASSERT(frntCnt < config->num_front_channel_elements);
        frontIsCpe[frntCnt] = isCpe;
        frontTag[frntCnt++] = tag;
        break;
      case SIDE:
        FDK_ASSERT(sdCnt < config->num_side_channel_elements);
        sideIsCpe[sdCnt] = isCpe;
        sideTag[sdCnt++] = tag;
        break;
      case BACK:
        FDK_ASSERT(bckCnt < config->num_back_channel_elements);
        backIsCpe[bckCnt] = isCpe;
        backTag[bckCnt++] = tag;
        break;
    }
  }

  for (i = 0; i < config->num_front_channel_elements; i++) {
    FDKwriteBits(hBs, frontIsCpe[i], 1);
    FDKwriteBits(hBs, frontTag[i], 4);
  }
  for (i = 0; i < config->num_side_channel_elements; i++) {
    FDKwriteBits(hBs, sideIsCpe[i], 1);
    FDKwriteBits(hBs, sideTag[i], 4);
  }
  for (i = 0; i < config->num_back_channel_elements; i++) {
    FDKwriteBits(hBs, backIsCpe[i], 1);
    FDKwriteBits(hBs, backTag[i], 4);
  }
  for (i = 0; i < config->num_lfe_channel_elements; i++) {
    FDKwriteBits(hBs, lfeCnt++, 4);
  }

  FDKbyteAlign(hBs, alignAnchor);

  if (config->pHeight_num != NULL) {

    INT commentBytes =
        1
        + ((((config->num_front_channel_elements +
              config->num_side_channel_elements +
              config->num_back_channel_elements)
             << 1) +
            7) >>
           3)
        + 1;

    FDKwriteBits(hBs, commentBytes, 8);

    FDK_CRCINFO crcInfo;
    INT crcReg;

    FDKcrcInit(&crcInfo, 0x07, 0xFF, 8);
    crcReg = FDKcrcStartReg(&crcInfo, hBs, 0);

    FDKwriteBits(hBs, PCE_HEIGHT_EXT_SYNC, 8);

    for (i = 0;
         i < config->num_front_channel_elements -
                 config->pHeight_num->num_front_height_channel_elements[0] -
                 config->pHeight_num->num_front_height_channel_elements[1];
         i++)
      FDKwriteBits(hBs, HEIGHT_NORMAL, 2);
    for (i = 0; i < config->pHeight_num->num_front_height_channel_elements[0];
         i++)
      FDKwriteBits(hBs, HEIGHT_TOP, 2);
    for (i = 0; i < config->pHeight_num->num_front_height_channel_elements[1];
         i++)
      FDKwriteBits(hBs, HEIGHT_BOTTOM, 2);

    for (i = 0;
         i < config->num_side_channel_elements -
                 config->pHeight_num->num_side_height_channel_elements[0] -
                 config->pHeight_num->num_side_height_channel_elements[1];
         i++)
      FDKwriteBits(hBs, HEIGHT_NORMAL, 2);
    for (i = 0; i < config->pHeight_num->num_side_height_channel_elements[0];
         i++)
      FDKwriteBits(hBs, HEIGHT_TOP, 2);
    for (i = 0; i < config->pHeight_num->num_side_height_channel_elements[1];
         i++)
      FDKwriteBits(hBs, HEIGHT_BOTTOM, 2);

    for (i = 0;
         i < config->num_back_channel_elements -
                 config->pHeight_num->num_back_height_channel_elements[0] -
                 config->pHeight_num->num_back_height_channel_elements[1];
         i++)
      FDKwriteBits(hBs, HEIGHT_NORMAL, 2);
    for (i = 0; i < config->pHeight_num->num_back_height_channel_elements[0];
         i++)
      FDKwriteBits(hBs, HEIGHT_TOP, 2);
    for (i = 0; i < config->pHeight_num->num_back_height_channel_elements[1];
         i++)
      FDKwriteBits(hBs, HEIGHT_BOTTOM, 2);

    FDKbyteAlign(hBs, alignAnchor);

    FDKcrcEndReg(&crcInfo, hBs, crcReg);
    FDKwriteBits(hBs, FDKcrcGetCRC(&crcInfo), 8);

  } else {
    FDKwriteBits(hBs, 0,
                 8);
  }

  return 0;
}

int transportEnc_GetPCEBits(CHANNEL_MODE channelMode, int matrixMixdownA,
                            int bits) {
  const PCE_CONFIGURATION *config = NULL;

  if ((config = getPceEntry(channelMode)) == NULL) {
    return -1;
  }

  bits +=
      4 + 2 + 4;
  bits += 4 + 4 + 4 + 2;
  bits += 3 + 4;
  bits += 1 + 1 + 1;

  if (matrixMixdownA != 0 &&
      ((channelMode == MODE_1_2_2) || (channelMode == MODE_1_2_2_1))) {
    bits += 3;
  }

  bits += (1 + 4) * (INT)config->num_front_channel_elements;
  bits += (1 + 4) * (INT)config->num_side_channel_elements;
  bits += (1 + 4) * (INT)config->num_back_channel_elements;
  bits += (4) * (INT)config->num_lfe_channel_elements;

  if ((bits % 8) != 0) {
    bits += (8 - (bits % 8));
  }

  bits += 8;

  if (config->pHeight_num != NULL) {

    bits +=
        8
        +
        ((config->num_front_channel_elements +
          config->num_side_channel_elements + config->num_back_channel_elements)
         << 1)
        + 8;

    if ((bits % 8) != 0) {
      bits += (8 - (bits % 8));
    }
  }

  return bits;
}

static void writeAot(HANDLE_FDK_BITSTREAM hBitstreamBuffer,
                     AUDIO_OBJECT_TYPE aot) {
  int tmp = (int)aot;

  if (tmp > 31) {
    FDKwriteBits(hBitstreamBuffer, AOT_ESCAPE, 5);
    FDKwriteBits(hBitstreamBuffer, tmp - 32, 6);
  } else {
    FDKwriteBits(hBitstreamBuffer, tmp, 5);
  }
}

static void writeSampleRate(HANDLE_FDK_BITSTREAM hBs, int sampleRate,
                            int nBits) {
  int srIdx = getSamplingRateIndex(sampleRate, nBits);

  FDKwriteBits(hBs, srIdx, nBits);
  if (srIdx == (1 << nBits) - 1) {
    FDKwriteBits(hBs, sampleRate, 24);
  }
}

static int transportEnc_writeGASpecificConfig(HANDLE_FDK_BITSTREAM asc,
                                              CODER_CONFIG *config, int extFlg,
                                              UINT alignAnchor) {
  int aot = config->aot;
  int samplesPerFrame = config->samplesPerFrame;

  FDKwriteBits(asc,
               ((samplesPerFrame == 960 || samplesPerFrame == 480) ? 1 : 0),
               1);
  FDKwriteBits(asc, 0,
               1);
  FDKwriteBits(asc, extFlg,
               1);

  if (getChannelConfig(config->channelMode, config->channelConfigZero) == 0) {
    transportEnc_writePCE(asc, config->channelMode, config->samplingRate, 0, 1,
                          config->matrixMixdownA,
                          (config->flags & CC_PSEUDO_SURROUND) ? 1 : 0,
                          alignAnchor);
  }
  if ((aot == AOT_AAC_SCAL) || (aot == AOT_ER_AAC_SCAL)) {
    FDKwriteBits(asc, 0, 3);
  }
  if (extFlg) {
    if (aot == AOT_ER_BSAC) {
      FDKwriteBits(asc, config->BSACnumOfSubFrame, 5);
      FDKwriteBits(asc, config->BSAClayerLength, 11);
    }
    if ((aot == AOT_ER_AAC_LC) || (aot == AOT_ER_AAC_LTP) ||
        (aot == AOT_ER_AAC_SCAL) || (aot == AOT_ER_AAC_LD)) {
      FDKwriteBits(asc, (config->flags & CC_VCB11) ? 1 : 0,
                   1);
      FDKwriteBits(asc, (config->flags & CC_RVLC) ? 1 : 0,
                   1);
      FDKwriteBits(asc, (config->flags & CC_HCR) ? 1 : 0,
                   1);
    }
    FDKwriteBits(asc, 0, 1);
  }
  return 0;
}

static int transportEnc_writeELDSpecificConfig(HANDLE_FDK_BITSTREAM hBs,
                                               CODER_CONFIG *config,
                                               int epConfig,
                                               CSTpCallBacks *cb) {
  UINT frameLengthFlag = 0;
  switch (config->samplesPerFrame) {
    case 512:
    case 256:
    case 128:
    case 64:
      frameLengthFlag = 0;
      break;
    case 480:
    case 240:
    case 160:
    case 120:
    case 60:
      frameLengthFlag = 1;
      break;
  }

  FDKwriteBits(hBs, frameLengthFlag, 1);

  FDKwriteBits(hBs, (config->flags & CC_VCB11) ? 1 : 0, 1);
  FDKwriteBits(hBs, (config->flags & CC_RVLC) ? 1 : 0, 1);
  FDKwriteBits(hBs, (config->flags & CC_HCR) ? 1 : 0, 1);

  FDKwriteBits(hBs, (config->flags & CC_SBR) ? 1 : 0, 1);
  if ((config->flags & CC_SBR)) {
    FDKwriteBits(hBs, (config->samplingRate == config->extSamplingRate) ? 0 : 1,
                 1);
    FDKwriteBits(hBs, (config->flags & CC_SBRCRC) ? 1 : 0, 1);

    if (cb->cbSbr != NULL) {
      const PCE_CONFIGURATION *pPce;
      int e, sbrElementIndex = 0;

      pPce = getPceEntry(config->channelMode);

      for (e = 0; e < pPce->num_front_channel_elements +
                          pPce->num_side_channel_elements +
                          pPce->num_back_channel_elements +
                          pPce->num_lfe_channel_elements;
           e++) {
        if ((pPce->pEl_type[e] == ID_SCE) || (pPce->pEl_type[e] == ID_CPE)) {
          cb->cbSbr(cb->cbSbrData, hBs, 0, 0, 0, config->aot, pPce->pEl_type[e],
                    sbrElementIndex, 0, 0, 0, NULL, 1);
          sbrElementIndex++;
        }
      }
    }
  }

  if ((config->flags & CC_SAC) && (cb->cbSsc != NULL)) {
    FDKwriteBits(hBs, ELDEXT_LDSAC, 4);

    const INT eldExtLen =
        (cb->cbSsc(cb->cbSscData, NULL, config->aot, config->extSamplingRate, 0,
                   0, 0, 0, 0, 0, NULL) +
         7) >>
        3;
    INT cnt = eldExtLen;

    if (cnt < 0xF) {
      FDKwriteBits(hBs, cnt, 4);
    } else {
      FDKwriteBits(hBs, 0xF, 4);
      cnt -= 0xF;

      if (cnt < 0xFF) {
        FDKwriteBits(hBs, cnt, 8);
      } else {
        FDKwriteBits(hBs, 0xFF, 8);
        cnt -= 0xFF;

        FDK_ASSERT(cnt <= 0xFFFF);
        FDKwriteBits(hBs, cnt, 16);
      }
    }

    cb->cbSsc(cb->cbSscData, hBs, config->aot, config->extSamplingRate, 0, 0, 0,
              0, 0, 0, NULL);
  }

  if (config->downscaleSamplingRate != 0 &&
      config->downscaleSamplingRate != config->extSamplingRate) {

    int eldExtLenDsc = 1;
    int downscaleSamplingRate = config->downscaleSamplingRate;
    FDKwriteBits(hBs, ELDEXT_DOWNSCALEINFO, 4);

    if ((downscaleSamplingRate != 96000) && (downscaleSamplingRate != 88200) &&
        (downscaleSamplingRate != 64000) && (downscaleSamplingRate != 48000) &&
        (downscaleSamplingRate != 44100) && (downscaleSamplingRate != 32000) &&
        (downscaleSamplingRate != 24000) && (downscaleSamplingRate != 22050) &&
        (downscaleSamplingRate != 16000) && (downscaleSamplingRate != 12000) &&
        (downscaleSamplingRate != 11025) && (downscaleSamplingRate != 8000) &&
        (downscaleSamplingRate != 7350)) {
      eldExtLenDsc = 4;
    }

    FDKwriteBits(hBs, eldExtLenDsc, 4);
    writeSampleRate(hBs, downscaleSamplingRate, 4);
    FDKwriteBits(hBs, 0x0, 4);
  }

  FDKwriteBits(hBs, ELDEXT_TERM, 4);

  return 0;
}

static int transportEnc_writeUsacSpecificConfig(HANDLE_FDK_BITSTREAM hBs,
                                                int extFlag, CODER_CONFIG *cc,
                                                CSTpCallBacks *cb) {
  FDK_BITSTREAM usacConf;
  int usacConfigBits = cc->rawConfigBits;

  if ((usacConfigBits <= 0) ||
      ((usacConfigBits + 7) / 8 > (int)sizeof(cc->rawConfig))) {
    return TRANSPORTENC_UNSUPPORTED_FORMAT;
  }
  FDKinitBitStream(&usacConf, cc->rawConfig, BUFSIZE_DUMMY_VALUE,
                   usacConfigBits, BS_READER);

  for (; usacConfigBits > 0; usacConfigBits--) {
    UINT tmp = FDKreadBit(&usacConf);
    FDKwriteBits(hBs, tmp, 1);
  }
  FDKsyncCache(hBs);

  return TRANSPORTENC_OK;
}

int transportEnc_writeASC(HANDLE_FDK_BITSTREAM asc, CODER_CONFIG *config,
                          CSTpCallBacks *cb) {
  UINT extFlag = 0;
  int err;
  int epConfig = 0;

  UINT alignAnchor = FDKgetValidBits(asc);

  switch (config->aot) {
    case AOT_ER_AAC_LC:
    case AOT_ER_AAC_LTP:
    case AOT_ER_AAC_SCAL:
    case AOT_ER_TWIN_VQ:
    case AOT_ER_BSAC:
    case AOT_ER_AAC_LD:
    case AOT_ER_AAC_ELD:
    case AOT_USAC:
      extFlag = 1;
      break;
    default:
      break;
  }

  if (config->sbrSignaling == SIG_EXPLICIT_HIERARCHICAL && config->sbrPresent)
    writeAot(asc, config->extAOT);
  else
    writeAot(asc, config->aot);

  writeSampleRate(asc, config->samplingRate, 4);

  if (config->channelMode == MODE_INVALID) {
    config->channelMode = transportEnc_GetChannelMode(config->noChannels);
    if (config->channelMode == MODE_INVALID) return -1;
  }

  FDKwriteBits(
      asc, getChannelConfig(config->channelMode, config->channelConfigZero), 4);

  if (config->sbrSignaling == SIG_EXPLICIT_HIERARCHICAL && config->sbrPresent) {
    writeSampleRate(asc, config->extSamplingRate, 4);
    writeAot(asc, config->aot);
  }

  switch (config->aot) {
    case AOT_AAC_MAIN:
    case AOT_AAC_LC:
    case AOT_AAC_SSR:
    case AOT_AAC_LTP:
    case AOT_AAC_SCAL:
    case AOT_TWIN_VQ:
    case AOT_ER_AAC_LC:
    case AOT_ER_AAC_LTP:
    case AOT_ER_AAC_SCAL:
    case AOT_ER_TWIN_VQ:
    case AOT_ER_BSAC:
    case AOT_ER_AAC_LD:
      err =
          transportEnc_writeGASpecificConfig(asc, config, extFlag, alignAnchor);
      if (err) return err;
      break;
    case AOT_ER_AAC_ELD:
      err = transportEnc_writeELDSpecificConfig(asc, config, epConfig, cb);
      if (err) return err;
      break;
    case AOT_USAC:
      err = transportEnc_writeUsacSpecificConfig(asc, extFlag, config, cb);
      if (err) {
        return err;
      }
      break;
    default:
      return -1;
  }

  switch (config->aot) {
    case AOT_ER_AAC_LC:
    case AOT_ER_AAC_LTP:
    case AOT_ER_AAC_SCAL:
    case AOT_ER_TWIN_VQ:
    case AOT_ER_BSAC:
    case AOT_ER_AAC_LD:
    case AOT_ER_CELP:
    case AOT_ER_HVXC:
    case AOT_ER_HILN:
    case AOT_ER_PARA:
    case AOT_ER_AAC_ELD:
      FDKwriteBits(asc, 0, 2);
      break;
    default:
      break;
  }

  if (config->sbrSignaling == SIG_EXPLICIT_BW_COMPATIBLE) {
    TP_ASC_EXTENSION_ID ascExtId = ASCEXT_UNKOWN;

    if (config->sbrPresent) {
      ascExtId = ASCEXT_SBR;
      FDKwriteBits(asc, ascExtId, 11);
      writeAot(asc, config->extAOT);
      FDKwriteBits(asc, 1, 1);
      writeSampleRate(asc, config->extSamplingRate, 4);
      if (config->psPresent) {
        ascExtId = ASCEXT_PS;
        FDKwriteBits(asc, ascExtId, 11);
        FDKwriteBits(asc, 1, 1);
      }
    }
  }

  FDKsyncCache(asc);

  return 0;
}
