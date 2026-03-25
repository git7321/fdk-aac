/**************************** AAC encoder library ******************************

   Author(s):   M. Schug / A. Groeschel

   Description: bandwidth expert

*******************************************************************************/

#include "channel_map.h"
#include "bandwidth.h"
#include "aacEnc_ram.h"

typedef struct {
  INT chanBitRate;
  INT bandWidthMono;
  INT bandWidth2AndMoreChan;
} BANDWIDTH_TAB;

static const BANDWIDTH_TAB bandWidthTable[] = {
    {0, 6000, 6000},       {12000, 9000, 9000},   {20000, 12000, 12000},
    {28000, 14000, 14000}, {40000, 16000, 16000}, {56000, 18000, 18000},
    {72000, 20000, 20000}, {96000, 20000, 20000}, {576001, 20000, 20000}};

static const BANDWIDTH_TAB bandWidthTable_LD_22050[] = {
    {8000, 6000, 6000},    {12000, 8000, 8000},   {16000, 10000, 10000},
    {24000, 12000, 12000}, {32000, 16000, 16000}, {40000, 18000, 18000},
    {48000, 20000, 20000}, {282241, 20000, 20000}};

static const BANDWIDTH_TAB bandWidthTable_LD_24000[] = {
    {8000, 6000, 6000},    {12000, 8000, 8000},   {16000, 10000, 10000},
    {24000, 12000, 12000}, {32000, 16000, 16000}, {40000, 18000, 18000},
    {48000, 20000, 20000}, {64000, 20000, 20000}, {307201, 20000, 20000}};

static const BANDWIDTH_TAB bandWidthTable_LD_32000[] = {
    {8000, 6000, 6000},    {12000, 8000, 8000},   {24000, 10000, 10000},
    {32000, 12000, 12000}, {40000, 16000, 16000}, {48000, 18000, 18000},
    {64000, 20000, 20000}, {76000, 20000, 20000}, {409601, 20000, 20000}};

static const BANDWIDTH_TAB bandWidthTable_LD_44100[] = {
    {8000, 6000, 6000},     {24000, 8000, 8000},   {32000, 10000, 10000},
    {40000, 12000, 12000},  {48000, 14000, 14000}, {56000, 16000, 16000},
    {64000, 18000, 18000},  {80000, 18000, 18000}, {96000, 20000, 20000},
    {128000, 20000, 20000}, {564481, 20000, 20000}};

static const BANDWIDTH_TAB bandWidthTable_LD_48000[] = {
    {8000, 6000, 6000},     {24000, 8000, 8000},   {32000, 10000, 10000},
    {40000, 12000, 12000},  {48000, 14000, 14000}, {56000, 16000, 16000},
    {64000, 18000, 18000},  {80000, 18000, 18000}, {96000, 20000, 20000},
    {128000, 20000, 20000}, {614401, 20000, 20000}};

typedef struct {
  AACENC_BITRATE_MODE bitrateMode;
  int bandWidthMono;
  int bandWidth2AndMoreChan;
} BANDWIDTH_TAB_VBR;

static const BANDWIDTH_TAB_VBR bandWidthTableVBR[] = {
    {AACENC_BR_MODE_CBR, 0, 0},
    {AACENC_BR_MODE_VBR_1, 18000, 18000},
    {AACENC_BR_MODE_VBR_2, 18000, 18000},
    {AACENC_BR_MODE_VBR_3, 18000, 18000},
    {AACENC_BR_MODE_VBR_4, 19000, 19000},
    {AACENC_BR_MODE_VBR_5, 20000, 20000},
    {AACENC_BR_MODE_SFR, 0, 0},
    {AACENC_BR_MODE_FF, 0, 0}
};

static INT GetBandwidthEntry(const INT frameLength, const INT sampleRate,
                             const INT chanBitRate, const INT entryNo) {
  INT bandwidth = -1;
  const BANDWIDTH_TAB *pBwTab = NULL;
  INT bwTabSize = 0;

  switch (frameLength) {
    case 960:
    case 1024:
      pBwTab = bandWidthTable;
      bwTabSize = sizeof(bandWidthTable) / sizeof(BANDWIDTH_TAB);
      break;
    case 120:
    case 128:
    case 240:
    case 256:
    case 480:
    case 512:
      switch (sampleRate) {
        case 8000:
        case 11025:
        case 12000:
        case 16000:
        case 22050:
          pBwTab = bandWidthTable_LD_22050;
          bwTabSize = sizeof(bandWidthTable_LD_22050) / sizeof(BANDWIDTH_TAB);
          break;
        case 24000:
          pBwTab = bandWidthTable_LD_24000;
          bwTabSize = sizeof(bandWidthTable_LD_24000) / sizeof(BANDWIDTH_TAB);
          break;
        case 32000:
          pBwTab = bandWidthTable_LD_32000;
          bwTabSize = sizeof(bandWidthTable_LD_32000) / sizeof(BANDWIDTH_TAB);
          break;
        case 44100:
          pBwTab = bandWidthTable_LD_44100;
          bwTabSize = sizeof(bandWidthTable_LD_44100) / sizeof(BANDWIDTH_TAB);
          break;
        case 48000:
        case 64000:
        case 88200:
        case 96000:
          pBwTab = bandWidthTable_LD_48000;
          bwTabSize = sizeof(bandWidthTable_LD_48000) / sizeof(BANDWIDTH_TAB);
          break;
      }
      break;
    default:
      pBwTab = NULL;
      bwTabSize = 0;
  }

  if (pBwTab != NULL) {
    int i;
    for (i = 0; i < bwTabSize - 1; i++) {
      if (chanBitRate >= pBwTab[i].chanBitRate &&
          chanBitRate < pBwTab[i + 1].chanBitRate) {
        switch (frameLength) {
          case 960:
          case 1024:
            bandwidth = (entryNo == 0) ? pBwTab[i].bandWidthMono
                                       : pBwTab[i].bandWidth2AndMoreChan;
            break;
          case 120:
          case 128:
          case 240:
          case 256:
          case 480:
          case 512: {
            INT q_res = 0;
            INT startBw = (entryNo == 0) ? pBwTab[i].bandWidthMono
                                         : pBwTab[i].bandWidth2AndMoreChan;
            INT endBw = (entryNo == 0) ? pBwTab[i + 1].bandWidthMono
                                       : pBwTab[i + 1].bandWidth2AndMoreChan;
            INT startBr = pBwTab[i].chanBitRate;
            INT endBr = pBwTab[i + 1].chanBitRate;

            FIXP_DBL bwFac_fix =
                fDivNorm(chanBitRate - startBr, endBr - startBr, &q_res);
            bandwidth =
                (INT)scaleValue(fMult(bwFac_fix, (FIXP_DBL)(endBw - startBw)),
                 q_res) +
                startBw;
          } break;
          default:
            bandwidth = -1;
        }
        break;
      }
    }
  }

  return bandwidth;
}

AAC_ENCODER_ERROR FDKaacEnc_DetermineBandWidth(
    const INT proposedBandWidth, const INT bitrate,
    const AACENC_BITRATE_MODE bitrateMode, const INT sampleRate,
    const INT frameLength, const CHANNEL_MAPPING *const cm,
    const CHANNEL_MODE encoderMode, INT *const bandWidth) {
  AAC_ENCODER_ERROR ErrorStatus = AAC_ENC_OK;
  INT chanBitRate = bitrate / cm->nChannelsEff;

  switch (bitrateMode) {
    case AACENC_BR_MODE_VBR_1:
    case AACENC_BR_MODE_VBR_2:
    case AACENC_BR_MODE_VBR_3:
    case AACENC_BR_MODE_VBR_4:
    case AACENC_BR_MODE_VBR_5:
      if (proposedBandWidth != 0) {
        *bandWidth = proposedBandWidth;
      } else {
        switch (encoderMode) {
          case MODE_1:
            *bandWidth = bandWidthTableVBR[bitrateMode].bandWidthMono;
            break;
          case MODE_2:
          case MODE_1_2:
          case MODE_1_2_1:
          case MODE_1_2_2:
          case MODE_1_2_2_1:
          case MODE_6_1:
          case MODE_1_2_2_2_1:
          case MODE_7_1_REAR_SURROUND:
          case MODE_7_1_FRONT_CENTER:
          case MODE_7_1_BACK:
          case MODE_7_1_TOP_FRONT:
            *bandWidth = bandWidthTableVBR[bitrateMode].bandWidth2AndMoreChan;
            break;
          default:
            return AAC_ENC_UNSUPPORTED_CHANNELCONFIG;
        }
      }
      break;
    case AACENC_BR_MODE_CBR:
    case AACENC_BR_MODE_SFR:
    case AACENC_BR_MODE_FF:

      if (proposedBandWidth != 0) {
        *bandWidth = fMin(proposedBandWidth, fMin(20000, sampleRate >> 1));
      } else {
        int entryNo = 0;

        switch (encoderMode) {
          case MODE_1:
            entryNo = 0;
            break;

          case MODE_2:
          case MODE_1_2:
          case MODE_1_2_1:
          case MODE_1_2_2:
          case MODE_1_2_2_1:
          case MODE_6_1:
          case MODE_1_2_2_2_1:
          case MODE_7_1_REAR_SURROUND:
          case MODE_7_1_FRONT_CENTER:
          case MODE_7_1_BACK:
          case MODE_7_1_TOP_FRONT:
            entryNo = 1;
            break;

          default:
            return AAC_ENC_UNSUPPORTED_CHANNELCONFIG;
        }
        *bandWidth =
            GetBandwidthEntry(frameLength, sampleRate, chanBitRate, entryNo);

        if (*bandWidth == -1) {
          switch (frameLength) {
            case 120:
            case 128:
            case 240:
            case 256:
              *bandWidth = 20000;
              break;
            default:
              ErrorStatus = AAC_ENC_INVALID_CHANNEL_BITRATE;
          }
        }
      }
      break;
    default:
      *bandWidth = 0;
      return AAC_ENC_UNSUPPORTED_BITRATE_MODE;
  }
  *bandWidth = fMin(*bandWidth, sampleRate / 2);

  return ErrorStatus;
}
