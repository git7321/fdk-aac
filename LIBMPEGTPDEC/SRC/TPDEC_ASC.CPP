/******************* MPEG transport format decoder library *********************

   Author(s):   Daniel Homm

   Description:

*******************************************************************************/

#include "tpdec_lib.h"
#include "tp_data.h"
#include "FDK_crc.h"
#include "common_fix.h"

static const MP4_ELEMENT_ID channel_configuration_0[] = {ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_1[] = {ID_SCE, ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_2[] = {ID_CPE, ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_3[] = {ID_SCE, ID_CPE,
                                                         ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_4[] = {ID_SCE, ID_CPE, ID_SCE,
                                                         ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_5[] = {ID_SCE, ID_CPE, ID_CPE,
                                                         ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_6[] = {ID_SCE, ID_CPE, ID_CPE,
                                                         ID_LFE, ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_7[] = {
    ID_SCE, ID_CPE, ID_CPE, ID_CPE, ID_LFE, ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_8[] = {
    ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_9[] = {
    ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_10[] = {
    ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_11[] = {
    ID_SCE, ID_CPE, ID_CPE, ID_SCE, ID_LFE, ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_12[] = {
    ID_SCE, ID_CPE, ID_CPE, ID_CPE, ID_LFE, ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_13[] = {
    ID_SCE, ID_CPE, ID_CPE, ID_CPE, ID_CPE, ID_SCE, ID_LFE, ID_LFE, ID_SCE,
    ID_CPE, ID_CPE, ID_SCE, ID_CPE, ID_SCE, ID_SCE, ID_CPE, ID_NONE};
static const MP4_ELEMENT_ID channel_configuration_14[] = {
    ID_SCE, ID_CPE, ID_CPE, ID_LFE, ID_CPE, ID_NONE};

static const MP4_ELEMENT_ID *channel_configuration_array[] = {
    channel_configuration_0,  channel_configuration_1,
    channel_configuration_2,  channel_configuration_3,
    channel_configuration_4,  channel_configuration_5,
    channel_configuration_6,  channel_configuration_7,
    channel_configuration_8,  channel_configuration_9,
    channel_configuration_10, channel_configuration_11,
    channel_configuration_12, channel_configuration_13,
    channel_configuration_14};

#define TP_USAC_MAX_CHANNEL_CONFIGURATION_INDEX (13)
#define SC_CHANNEL_CONFIG_TAB_SIZE (TP_USAC_MAX_CHANNEL_CONFIGURATION_INDEX + 1)

typedef struct {
  SCHAR nCh;
  SCHAR nSCE;
  SCHAR nCPE;
  SCHAR nLFE;
} SC_CHANNEL_CONFIG;

static const SC_CHANNEL_CONFIG sc_chan_config_tab[SC_CHANNEL_CONFIG_TAB_SIZE] =
    {
        {0, 0, 0, 0},
        {1, 1, 0, 0},
        {2, 0, 1, 0},
        {3, 1, 1, 0},
        {4, 2, 1, 0},
        {5, 1, 2, 0},
        {6, 1, 2, 1},
        {8, 1, 3, 1},
        {2, 2, 0, 0},
        {3, 1, 1, 0},
        {4, 0, 2, 0},
        {7, 2, 2, 1},
        {8, 1, 3, 1},
        {24, 6, 8, 2}
};

void CProgramConfig_Reset(CProgramConfig *pPce) { pPce->elCounter = 0; }

void CProgramConfig_Init(CProgramConfig *pPce) {
  FDKmemclear(pPce, sizeof(CProgramConfig));
  pPce->SamplingFrequencyIndex = 0xf;
}

int CProgramConfig_IsValid(const CProgramConfig *pPce) {
  return ((pPce->isValid) ? 1 : 0);
}

#define PCE_HEIGHT_EXT_SYNC (0xAC)

static int CProgramConfig_ReadHeightExt(CProgramConfig *pPce,
                                        HANDLE_FDK_BITSTREAM bs,
                                        int *const bytesAvailable,
                                        const UINT alignmentAnchor) {
  int err = 0;
  FDK_CRCINFO crcInfo;
  INT crcReg;
  FDKcrcInit(&crcInfo, 0x07, 0xFF, 8);
  crcReg = FDKcrcStartReg(&crcInfo, bs, 0);
  UINT startAnchor = FDKgetValidBits(bs);

  FDK_ASSERT(pPce != NULL);
  FDK_ASSERT(bs != NULL);
  FDK_ASSERT(bytesAvailable != NULL);

  if ((startAnchor >= 24) && (*bytesAvailable >= 3) &&
      (FDKreadBits(bs, 8) == PCE_HEIGHT_EXT_SYNC)) {
    int i;

    for (i = 0; i < pPce->NumFrontChannelElements; i++) {
      if ((pPce->FrontElementHeightInfo[i] = (UCHAR)FDKreadBits(bs, 2)) >=
          PC_NUM_HEIGHT_LAYER) {
        err = -2;
      }
    }
    for (i = 0; i < pPce->NumSideChannelElements; i++) {
      if ((pPce->SideElementHeightInfo[i] = (UCHAR)FDKreadBits(bs, 2)) >=
          PC_NUM_HEIGHT_LAYER) {
        err = -2;
      }
    }
    for (i = 0; i < pPce->NumBackChannelElements; i++) {
      if ((pPce->BackElementHeightInfo[i] = (UCHAR)FDKreadBits(bs, 2)) >=
          PC_NUM_HEIGHT_LAYER) {
        err = -2;
      }
    }
    FDKbyteAlign(bs, alignmentAnchor);

    FDKcrcEndReg(&crcInfo, bs, crcReg);
    if ((USHORT)FDKreadBits(bs, 8) != FDKcrcGetCRC(&crcInfo)) {
      err = -1;
    }
    if (err != 0) {
      FDKmemclear(pPce->FrontElementHeightInfo,
                  sizeof(pPce->FrontElementHeightInfo));
      FDKmemclear(pPce->SideElementHeightInfo,
                  sizeof(pPce->SideElementHeightInfo));
      FDKmemclear(pPce->BackElementHeightInfo,
                  sizeof(pPce->BackElementHeightInfo));
    }
  } else {
    FDKpushBack(bs, (INT)startAnchor - (INT)FDKgetValidBits(bs));
  }

  *bytesAvailable -= ((INT)startAnchor - (INT)FDKgetValidBits(bs)) >> 3;

  return (err);
}

static int CProgramConfig_Check(CProgramConfig *pPce) {
  INT i;
  INT err = 0;
  INT numBackChannels[3] = {0};
  INT numSideChannels[3] = {0};
  INT numFrontChannels[3] = {0};
  UCHAR *pCpeFront = pPce->FrontElementIsCpe;
  UCHAR *pCpeSide = pPce->SideElementIsCpe;
  UCHAR *pCpeBack = pPce->BackElementIsCpe;
  UCHAR *pHeight;

  pHeight = pPce->BackElementHeightInfo;
  for (i = 0; i < pPce->NumBackChannelElements; i++) {
    numBackChannels[*pHeight] += pPce->BackElementIsCpe[i] ? 2 : 1;
    pHeight++;
  }
  pHeight = pPce->SideElementHeightInfo;
  for (i = 0; i < pPce->NumSideChannelElements; i++) {
    numSideChannels[*pHeight] += pPce->SideElementIsCpe[i] ? 2 : 1;
    pHeight++;
  }
  pHeight = pPce->FrontElementHeightInfo;
  for (i = 0; i < pPce->NumFrontChannelElements; i++) {
    numFrontChannels[*pHeight] += pPce->FrontElementIsCpe[i] ? 2 : 1;
    pHeight++;
  }

  for (i = 0; i < 3; i++) {
    if (numFrontChannels[i] & 1) {
      if (*pCpeFront++ == ID_CPE) {
        err = 1;
        goto bail;
      }
      numFrontChannels[i]--;
    }
    while (numFrontChannels[i] > 0) {
      if (*pCpeFront++ == ID_SCE) {
        if (*pCpeFront++ == ID_CPE) {
          err = 1;
          goto bail;
        }
      }
      numFrontChannels[i] -= 2;
    };

    if (i != 1) {
      if (numSideChannels[i] & 1) {
        err = 1;
        goto bail;
      }
      while (numSideChannels[i] > 0) {
        if (*pCpeSide++ == ID_SCE) {
          if (*pCpeSide++ == ID_CPE) {
            err = 1;
            goto bail;
          }
        }
        numSideChannels[i] -= 2;
      };
    }

    while (numBackChannels[i] > 1) {
      if (*pCpeBack++ == ID_SCE) {
        if (*pCpeBack++ == ID_CPE) {
          err = 1;
          goto bail;
        }
      }
      numBackChannels[i] -= 2;
    };
    if (numBackChannels[i]) {
      if (*pCpeBack++ == ID_CPE) {
        err = 1;
        goto bail;
      }
    }
  }

bail:

  return err;
}

void CProgramConfig_Read(CProgramConfig *pPce, HANDLE_FDK_BITSTREAM bs,
                         UINT alignmentAnchor) {
  int i;
  int commentBytes;
  UCHAR tag, isCpe;
  UCHAR checkElementTagSelect[3][PC_FSB_CHANNELS_MAX] = {{0}};

  pPce->isValid = 1;
  pPce->NumEffectiveChannels = 0;
  pPce->NumChannels = 0;
  pPce->ElementInstanceTag = (UCHAR)FDKreadBits(bs, 4);
  pPce->Profile = (UCHAR)FDKreadBits(bs, 2);
  pPce->SamplingFrequencyIndex = (UCHAR)FDKreadBits(bs, 4);
  pPce->NumFrontChannelElements = (UCHAR)FDKreadBits(bs, 4);
  pPce->NumSideChannelElements = (UCHAR)FDKreadBits(bs, 4);
  pPce->NumBackChannelElements = (UCHAR)FDKreadBits(bs, 4);
  pPce->NumLfeChannelElements = (UCHAR)FDKreadBits(bs, 2);
  pPce->NumAssocDataElements = (UCHAR)FDKreadBits(bs, 3);
  pPce->NumValidCcElements = (UCHAR)FDKreadBits(bs, 4);

  if ((pPce->MonoMixdownPresent = (UCHAR)FDKreadBits(bs, 1)) != 0) {
    pPce->MonoMixdownElementNumber = (UCHAR)FDKreadBits(bs, 4);
  }

  if ((pPce->StereoMixdownPresent = (UCHAR)FDKreadBits(bs, 1)) != 0) {
    pPce->StereoMixdownElementNumber = (UCHAR)FDKreadBits(bs, 4);
  }

  if ((pPce->MatrixMixdownIndexPresent = (UCHAR)FDKreadBits(bs, 1)) != 0) {
    pPce->MatrixMixdownIndex = (UCHAR)FDKreadBits(bs, 2);
    pPce->PseudoSurroundEnable = (UCHAR)FDKreadBits(bs, 1);
  }

  for (i = 0; i < pPce->NumFrontChannelElements; i++) {
    pPce->FrontElementIsCpe[i] = isCpe = (UCHAR)FDKreadBits(bs, 1);
    pPce->FrontElementTagSelect[i] = tag = (UCHAR)FDKreadBits(bs, 4);
    pPce->NumChannels += pPce->FrontElementIsCpe[i] ? 2 : 1;

    if (checkElementTagSelect[isCpe][tag] == 0) {
      checkElementTagSelect[isCpe][tag] = 1;
    } else {
      pPce->isValid = 0;
    }
  }

  for (i = 0; i < pPce->NumSideChannelElements; i++) {
    pPce->SideElementIsCpe[i] = isCpe = (UCHAR)FDKreadBits(bs, 1);
    pPce->SideElementTagSelect[i] = tag = (UCHAR)FDKreadBits(bs, 4);
    pPce->NumChannels += pPce->SideElementIsCpe[i] ? 2 : 1;

    if (checkElementTagSelect[isCpe][tag] == 0) {
      checkElementTagSelect[isCpe][tag] = 1;
    } else {
      pPce->isValid = 0;
    }
  }

  for (i = 0; i < pPce->NumBackChannelElements; i++) {
    pPce->BackElementIsCpe[i] = isCpe = (UCHAR)FDKreadBits(bs, 1);
    pPce->BackElementTagSelect[i] = tag = (UCHAR)FDKreadBits(bs, 4);
    pPce->NumChannels += pPce->BackElementIsCpe[i] ? 2 : 1;

    if (checkElementTagSelect[isCpe][tag] == 0) {
      checkElementTagSelect[isCpe][tag] = 1;
    } else {
      pPce->isValid = 0;
    }
  }

  pPce->NumEffectiveChannels = pPce->NumChannels;

  for (i = 0; i < pPce->NumLfeChannelElements; i++) {
    pPce->LfeElementTagSelect[i] = tag = (UCHAR)FDKreadBits(bs, 4);
    pPce->NumChannels += 1;

    if (checkElementTagSelect[2][tag] == 0) {
      checkElementTagSelect[2][tag] = 1;
    } else {
      pPce->isValid = 0;
    }
  }

  for (i = 0; i < pPce->NumAssocDataElements; i++) {
    pPce->AssocDataElementTagSelect[i] = (UCHAR)FDKreadBits(bs, 4);
  }

  for (i = 0; i < pPce->NumValidCcElements; i++) {
    pPce->CcElementIsIndSw[i] = (UCHAR)FDKreadBits(bs, 1);
    pPce->ValidCcElementTagSelect[i] = (UCHAR)FDKreadBits(bs, 4);
  }

  FDKbyteAlign(bs, alignmentAnchor);

  pPce->CommentFieldBytes = (UCHAR)FDKreadBits(bs, 8);
  commentBytes = pPce->CommentFieldBytes;

  if (CProgramConfig_ReadHeightExt(pPce, bs, &commentBytes, alignmentAnchor)) {
    pPce->isValid = 0;
  }

  if (CProgramConfig_Check(pPce)) {
    pPce->isValid = 0;
  }

  for (i = 0; i < commentBytes; i++) {
    UCHAR text;

    text = (UCHAR)FDKreadBits(bs, 8);

    if (i < PC_COMMENTLENGTH) {
      pPce->Comment[i] = text;
    }
  }
}

int CProgramConfig_Compare(const CProgramConfig *const pPce1,
                           const CProgramConfig *const pPce2) {
  int result = 0;

  if (FDKmemcmp(pPce1, pPce2, sizeof(CProgramConfig)) !=
      0) {
    result = -1;

    if (pPce1->NumChannels ==
        pPce2->NumChannels) {
      result = 1;

      if (pPce1->NumFrontChannelElements != pPce2->NumFrontChannelElements) {
        result = 2;
      } else {
        int el, numCh1 = 0, numCh2 = 0;
        for (el = 0; el < pPce1->NumFrontChannelElements; el += 1) {
          if (pPce1->FrontElementHeightInfo[el] !=
              pPce2->FrontElementHeightInfo[el]) {
            result = 2;
            break;
          }
          numCh1 += pPce1->FrontElementIsCpe[el] ? 2 : 1;
          numCh2 += pPce2->FrontElementIsCpe[el] ? 2 : 1;
        }
        if (numCh1 != numCh2) {
          result = 2;
        }
      }
      if (pPce1->NumSideChannelElements != pPce2->NumSideChannelElements) {
        result = 2;
      } else {
        int el, numCh1 = 0, numCh2 = 0;
        for (el = 0; el < pPce1->NumSideChannelElements; el += 1) {
          if (pPce1->SideElementHeightInfo[el] !=
              pPce2->SideElementHeightInfo[el]) {
            result = 2;
            break;
          }
          numCh1 += pPce1->SideElementIsCpe[el] ? 2 : 1;
          numCh2 += pPce2->SideElementIsCpe[el] ? 2 : 1;
        }
        if (numCh1 != numCh2) {
          result = 2;
        }
      }
      if (pPce1->NumBackChannelElements != pPce2->NumBackChannelElements) {
        result = 2;
      } else {
        int el, numCh1 = 0, numCh2 = 0;
        for (el = 0; el < pPce1->NumBackChannelElements; el += 1) {
          if (pPce1->BackElementHeightInfo[el] !=
              pPce2->BackElementHeightInfo[el]) {
            result = 2;
            break;
          }
          numCh1 += pPce1->BackElementIsCpe[el] ? 2 : 1;
          numCh2 += pPce2->BackElementIsCpe[el] ? 2 : 1;
        }
        if (numCh1 != numCh2) {
          result = 2;
        }
      }
      if (pPce1->NumLfeChannelElements != pPce2->NumLfeChannelElements) {
        result = 2;
      }
    }
  }

  return result;
}

void CProgramConfig_GetDefault(CProgramConfig *pPce, const UINT channelConfig) {
  FDK_ASSERT(pPce != NULL);

  CProgramConfig_Init(pPce);
  pPce->Profile =
      1;

  switch (channelConfig) {
    case 32:
      pPce->NumFrontChannelElements = 2;
      pPce->FrontElementIsCpe[0] = 0;
      pPce->FrontElementIsCpe[1] = 1;
      pPce->NumSideChannelElements = 1;
      pPce->SideElementIsCpe[0] = 1;
      pPce->NumBackChannelElements = 1;
      pPce->BackElementIsCpe[0] = 1;
      pPce->NumLfeChannelElements = 1;
      pPce->NumChannels = 8;
      pPce->NumEffectiveChannels = 7;
      pPce->isValid = 1;
      break;
    case 12:
      pPce->BackElementIsCpe[1] = 1;
      pPce->NumChannels += 1;
      pPce->NumEffectiveChannels += 1;
      FDK_FALLTHROUGH;
    case 11:
      pPce->NumFrontChannelElements += 2;
      pPce->FrontElementIsCpe[0] = 0;
      pPce->FrontElementIsCpe[1] = 1;
      pPce->NumBackChannelElements += 2;
      pPce->BackElementIsCpe[0] = 1;
      pPce->BackElementIsCpe[1] += 0;
      pPce->NumLfeChannelElements += 1;
      pPce->NumChannels += 7;
      pPce->NumEffectiveChannels += 6;
      pPce->isValid = 1;
      break;
    case 14:
      pPce->FrontElementHeightInfo[2] = 1;
      FDK_FALLTHROUGH;
    case 7:
      pPce->NumFrontChannelElements += 1;
      pPce->FrontElementIsCpe[2] = 1;
      pPce->NumChannels += 2;
      pPce->NumEffectiveChannels += 2;
      FDK_FALLTHROUGH;
    case 6:
      pPce->NumLfeChannelElements += 1;
      pPce->NumChannels += 1;
      FDK_FALLTHROUGH;
    case 5:
    case 4:
      pPce->NumBackChannelElements += 1;
      pPce->BackElementIsCpe[0] = (channelConfig > 4) ? 1 : 0;
      pPce->NumChannels += (channelConfig > 4) ? 2 : 1;
      pPce->NumEffectiveChannels += (channelConfig > 4) ? 2 : 1;
      FDK_FALLTHROUGH;
    case 3:
      pPce->NumFrontChannelElements += 1;
      pPce->FrontElementIsCpe[1] = 1;
      pPce->NumChannels += 2;
      pPce->NumEffectiveChannels += 2;
      FDK_FALLTHROUGH;
    case 1:
      pPce->NumFrontChannelElements += 1;
      pPce->FrontElementIsCpe[0] = 0;
      pPce->NumChannels += 1;
      pPce->NumEffectiveChannels += 1;
      pPce->isValid = 1;
      break;
    case 2:
      pPce->NumFrontChannelElements = 1;
      pPce->FrontElementIsCpe[0] = 1;
      pPce->NumChannels += 2;
      pPce->NumEffectiveChannels += 2;
      pPce->isValid = 1;
      break;
    default:
      pPce->isValid = 0;
      break;
  }

  if (pPce->isValid) {
    int el, elTagSce = 0, elTagCpe = 0;

    for (el = 0; el < pPce->NumFrontChannelElements; el += 1) {
      pPce->FrontElementTagSelect[el] =
          (pPce->FrontElementIsCpe[el]) ? elTagCpe++ : elTagSce++;
    }
    for (el = 0; el < pPce->NumSideChannelElements; el += 1) {
      pPce->SideElementTagSelect[el] =
          (pPce->SideElementIsCpe[el]) ? elTagCpe++ : elTagSce++;
    }
    for (el = 0; el < pPce->NumBackChannelElements; el += 1) {
      pPce->BackElementTagSelect[el] =
          (pPce->BackElementIsCpe[el]) ? elTagCpe++ : elTagSce++;
    }
    elTagSce = 0;
    for (el = 0; el < pPce->NumLfeChannelElements; el += 1) {
      pPce->LfeElementTagSelect[el] = elTagSce++;
    }
  }
}

static void getImplicitAudioChannelTypeAndIndex(AUDIO_CHANNEL_TYPE *chType,
                                                UCHAR *chIndex,
                                                UINT channelConfig,
                                                UINT index) {
  if (index < 3) {
    *chType = ACT_FRONT;
    *chIndex = index;
  } else {
    switch (channelConfig) {
      case 4:
      case 5:
      case 6:
        switch (index) {
          case 3:
          case 4:
            *chType = ACT_BACK;
            *chIndex = index - 3;
            break;
          case 5:
            *chType = ACT_LFE;
            *chIndex = 0;
            break;
        }
        break;
      case 7:
        switch (index) {
          case 3:
          case 4:
            *chType = ACT_FRONT;
            *chIndex = index;
            break;
          case 5:
          case 6:
            *chType = ACT_BACK;
            *chIndex = index - 5;
            break;
          case 7:
            *chType = ACT_LFE;
            *chIndex = 0;
            break;
        }
        break;
      case 11:
        if (index < 6) {
          *chType = ACT_BACK;
          *chIndex = index - 3;
        } else {
          *chType = ACT_LFE;
          *chIndex = 0;
        }
        break;
      case 12:
        if (index < 7) {
          *chType = ACT_BACK;
          *chIndex = index - 3;
        } else {
          *chType = ACT_LFE;
          *chIndex = 0;
        }
        break;
      case 14:
        switch (index) {
          case 3:
          case 4:
            *chType = ACT_BACK;
            *chIndex = index - 3;
            break;
          case 5:
            *chType = ACT_LFE;
            *chIndex = 0;
            break;
          case 6:
          case 7:
            *chType = ACT_FRONT_TOP;
            *chIndex = index - 6;
            break;
        }
        break;
      default:
        *chType = ACT_NONE;
        break;
    }
  }
}

int CProgramConfig_LookupElement(CProgramConfig *pPce, UINT channelConfig,
                                 const UINT tag, const UINT channelIdx,
                                 UCHAR chMapping[], AUDIO_CHANNEL_TYPE chType[],
                                 UCHAR chIndex[], const UINT chDescrLen,
                                 UCHAR *elMapping, MP4_ELEMENT_ID elList[],
                                 MP4_ELEMENT_ID elType) {
  if (channelConfig > 0) {
    if (IS_CHANNEL_ELEMENT(elType)) {
      *elMapping = pPce->elCounter;
      if (elList[pPce->elCounter] != elType &&
          !IS_USAC_CHANNEL_ELEMENT(elType)) {
        if ((channelConfig == 2) &&
            (elType == ID_SCE)) {
          channelConfig = 1;
        } else if ((elList[pPce->elCounter] == ID_LFE) &&
                   (elType ==
                    ID_SCE)) {
          ;
        } else {
          return 0;
        }
      }
      getImplicitAudioChannelTypeAndIndex(
          &chType[channelIdx], &chIndex[channelIdx], channelConfig, channelIdx);
      if (elType == ID_CPE || elType == ID_USAC_CPE) {
        chType[channelIdx + 1] = chType[channelIdx];
        chIndex[channelIdx + 1] = chIndex[channelIdx] + 1;
      }
      pPce->elCounter++;
    }
    return 1;
  } else {
    if ((!pPce->isValid) || (pPce->NumChannels > chDescrLen)) {
      if (IS_USAC_CHANNEL_ELEMENT(elType)) {
        *elMapping = pPce->elCounter++;
      } else if (IS_MP4_CHANNEL_ELEMENT(elType)) {
        elList[pPce->elCounter] = elType;
        *elMapping = pPce->elCounter++;
      }
    } else {
      int isCpe = 0, i;
      int ec[PC_NUM_HEIGHT_LAYER] = {0};
      int cc[PC_NUM_HEIGHT_LAYER] = {0};
      int fc[PC_NUM_HEIGHT_LAYER] = {0};
      int sc[PC_NUM_HEIGHT_LAYER] = {0};
      int bc[PC_NUM_HEIGHT_LAYER] = {0};
      int lc = 0;

      switch (elType) {
        case ID_CPE:
          isCpe = 1;
          FDK_FALLTHROUGH;
        case ID_SCE:
          for (i = 0; i < pPce->NumFrontChannelElements; i++) {
            int heightLayer = pPce->FrontElementHeightInfo[i];
            if (isCpe == pPce->FrontElementIsCpe[i] &&
                pPce->FrontElementTagSelect[i] == tag) {
              int h, elIdx = ec[heightLayer], chIdx = cc[heightLayer];
              AUDIO_CHANNEL_TYPE aChType =
                  (AUDIO_CHANNEL_TYPE)((heightLayer << 4) | ACT_FRONT);
              for (h = heightLayer - 1; h >= 0; h -= 1) {
                int el;
                for (el = 0; el < pPce->NumFrontChannelElements; el += 1) {
                  if (pPce->FrontElementHeightInfo[el] == h) {
                    elIdx += 1;
                    chIdx += (pPce->FrontElementIsCpe[el]) ? 2 : 1;
                  }
                }
                for (el = 0; el < pPce->NumSideChannelElements; el += 1) {
                  if (pPce->SideElementHeightInfo[el] == h) {
                    elIdx += 1;
                    chIdx += (pPce->SideElementIsCpe[el]) ? 2 : 1;
                  }
                }
                for (el = 0; el < pPce->NumBackChannelElements; el += 1) {
                  if (pPce->BackElementHeightInfo[el] == h) {
                    elIdx += 1;
                    chIdx += (pPce->BackElementIsCpe[el]) ? 2 : 1;
                  }
                }
                if (h == 0) {
                  elIdx += pPce->NumLfeChannelElements;
                  chIdx += pPce->NumLfeChannelElements;
                }
              }
              chMapping[chIdx] = channelIdx;
              chType[chIdx] = aChType;
              chIndex[chIdx] = fc[heightLayer];
              if (isCpe) {
                chMapping[chIdx + 1] = channelIdx + 1;
                chType[chIdx + 1] = aChType;
                chIndex[chIdx + 1] = fc[heightLayer] + 1;
              }
              *elMapping = elIdx;
              return 1;
            }
            ec[heightLayer] += 1;
            if (pPce->FrontElementIsCpe[i]) {
              cc[heightLayer] += 2;
              fc[heightLayer] += 2;
            } else {
              cc[heightLayer] += 1;
              fc[heightLayer] += 1;
            }
          }
          for (i = 0; i < pPce->NumSideChannelElements; i++) {
            int heightLayer = pPce->SideElementHeightInfo[i];
            if (isCpe == pPce->SideElementIsCpe[i] &&
                pPce->SideElementTagSelect[i] == tag) {
              int h, elIdx = ec[heightLayer], chIdx = cc[heightLayer];
              AUDIO_CHANNEL_TYPE aChType =
                  (AUDIO_CHANNEL_TYPE)((heightLayer << 4) | ACT_SIDE);
              for (h = heightLayer - 1; h >= 0; h -= 1) {
                int el;
                for (el = 0; el < pPce->NumFrontChannelElements; el += 1) {
                  if (pPce->FrontElementHeightInfo[el] == h) {
                    elIdx += 1;
                    chIdx += (pPce->FrontElementIsCpe[el]) ? 2 : 1;
                  }
                }
                for (el = 0; el < pPce->NumSideChannelElements; el += 1) {
                  if (pPce->SideElementHeightInfo[el] == h) {
                    elIdx += 1;
                    chIdx += (pPce->SideElementIsCpe[el]) ? 2 : 1;
                  }
                }
                for (el = 0; el < pPce->NumBackChannelElements; el += 1) {
                  if (pPce->BackElementHeightInfo[el] == h) {
                    elIdx += 1;
                    chIdx += (pPce->BackElementIsCpe[el]) ? 2 : 1;
                  }
                }
                if (h ==
                    0) {
                  elIdx += pPce->NumLfeChannelElements;
                  chIdx += pPce->NumLfeChannelElements;
                }
              }
              chMapping[chIdx] = channelIdx;
              chType[chIdx] = aChType;
              chIndex[chIdx] = sc[heightLayer];
              if (isCpe) {
                chMapping[chIdx + 1] = channelIdx + 1;
                chType[chIdx + 1] = aChType;
                chIndex[chIdx + 1] = sc[heightLayer] + 1;
              }
              *elMapping = elIdx;
              return 1;
            }
            ec[heightLayer] += 1;
            if (pPce->SideElementIsCpe[i]) {
              cc[heightLayer] += 2;
              sc[heightLayer] += 2;
            } else {
              cc[heightLayer] += 1;
              sc[heightLayer] += 1;
            }
          }
          for (i = 0; i < pPce->NumBackChannelElements; i++) {
            int heightLayer = pPce->BackElementHeightInfo[i];
            if (isCpe == pPce->BackElementIsCpe[i] &&
                pPce->BackElementTagSelect[i] == tag) {
              int h, elIdx = ec[heightLayer], chIdx = cc[heightLayer];
              AUDIO_CHANNEL_TYPE aChType =
                  (AUDIO_CHANNEL_TYPE)((heightLayer << 4) | ACT_BACK);
              for (h = heightLayer - 1; h >= 0; h -= 1) {
                int el;
                for (el = 0; el < pPce->NumFrontChannelElements; el += 1) {
                  if (pPce->FrontElementHeightInfo[el] == h) {
                    elIdx += 1;
                    chIdx += (pPce->FrontElementIsCpe[el]) ? 2 : 1;
                  }
                }
                for (el = 0; el < pPce->NumSideChannelElements; el += 1) {
                  if (pPce->SideElementHeightInfo[el] == h) {
                    elIdx += 1;
                    chIdx += (pPce->SideElementIsCpe[el]) ? 2 : 1;
                  }
                }
                for (el = 0; el < pPce->NumBackChannelElements; el += 1) {
                  if (pPce->BackElementHeightInfo[el] == h) {
                    elIdx += 1;
                    chIdx += (pPce->BackElementIsCpe[el]) ? 2 : 1;
                  }
                }
                if (h == 0) {
                  elIdx += pPce->NumLfeChannelElements;
                  chIdx += pPce->NumLfeChannelElements;
                }
              }
              chMapping[chIdx] = channelIdx;
              chType[chIdx] = aChType;
              chIndex[chIdx] = bc[heightLayer];
              if (isCpe) {
                chMapping[chIdx + 1] = channelIdx + 1;
                chType[chIdx + 1] = aChType;
                chIndex[chIdx + 1] = bc[heightLayer] + 1;
              }
              *elMapping = elIdx;
              return 1;
            }
            ec[heightLayer] += 1;
            if (pPce->BackElementIsCpe[i]) {
              cc[heightLayer] += 2;
              bc[heightLayer] += 2;
            } else {
              cc[heightLayer] += 1;
              bc[heightLayer] += 1;
            }
          }
          break;

        case ID_LFE: {
          for (i = 0; i < pPce->NumFrontChannelElements; i += 1) {
            int heightLayer = pPce->FrontElementHeightInfo[i];
            ec[heightLayer] += 1;
            cc[heightLayer] += (pPce->FrontElementIsCpe[i]) ? 2 : 1;
          }
          for (i = 0; i < pPce->NumSideChannelElements; i += 1) {
            int heightLayer = pPce->SideElementHeightInfo[i];
            ec[heightLayer] += 1;
            cc[heightLayer] += (pPce->SideElementIsCpe[i]) ? 2 : 1;
          }
          for (i = 0; i < pPce->NumBackChannelElements; i += 1) {
            int heightLayer = pPce->BackElementHeightInfo[i];
            ec[heightLayer] += 1;
            cc[heightLayer] += (pPce->BackElementIsCpe[i]) ? 2 : 1;
          }

          for (i = 0; i < pPce->NumLfeChannelElements; i++) {
            int elIdx =
                ec[0];
            int chIdx = cc[0];
            if (pPce->LfeElementTagSelect[i] == tag) {
              chMapping[chIdx] = channelIdx;
              *elMapping = elIdx;
              chType[chIdx] = ACT_LFE;
              chIndex[chIdx] = lc;
              return 1;
            }
            ec[0] += 1;
            cc[0] += 1;
            lc += 1;
          }
        } break;

        case ID_CCE:
          for (i = 0; i < pPce->NumValidCcElements; i++) {
            if (pPce->ValidCcElementTagSelect[i] == tag) {
              return 1;
            }
          }
          break;
        case ID_DSE:
          for (i = 0; i < pPce->NumAssocDataElements; i++) {
            if (pPce->AssocDataElementTagSelect[i] == tag) {
              return 1;
            }
          }
          break;
        default:
          return 0;
      }
      return 0;
    }
  }

  return 1;
}

#define SPEAKER_PLANE_NORMAL 0
#define SPEAKER_PLANE_TOP 1
#define SPEAKER_PLANE_BOTTOM 2

void CProgramConfig_GetChannelDescription(const UINT chConfig,
                                          const CProgramConfig *pPce,
                                          AUDIO_CHANNEL_TYPE chType[],
                                          UCHAR chIndex[]) {
  FDK_ASSERT(chType != NULL);
  FDK_ASSERT(chIndex != NULL);

  if ((chConfig == 0) && (pPce != NULL)) {
    if (pPce->isValid) {
      int spkPlane, chIdx = 0;
      for (spkPlane = SPEAKER_PLANE_NORMAL; spkPlane <= SPEAKER_PLANE_BOTTOM;
           spkPlane += 1) {
        int elIdx, grpChIdx = 0;
        for (elIdx = 0; elIdx < pPce->NumFrontChannelElements; elIdx += 1) {
          if (pPce->FrontElementHeightInfo[elIdx] == spkPlane) {
            chType[chIdx] = (AUDIO_CHANNEL_TYPE)((spkPlane << 4) | ACT_FRONT);
            chIndex[chIdx++] = grpChIdx++;
            if (pPce->FrontElementIsCpe[elIdx]) {
              chType[chIdx] = (AUDIO_CHANNEL_TYPE)((spkPlane << 4) | ACT_FRONT);
              chIndex[chIdx++] = grpChIdx++;
            }
          }
        }
        grpChIdx = 0;
        for (elIdx = 0; elIdx < pPce->NumSideChannelElements; elIdx += 1) {
          if (pPce->SideElementHeightInfo[elIdx] == spkPlane) {
            chType[chIdx] = (AUDIO_CHANNEL_TYPE)((spkPlane << 4) | ACT_SIDE);
            chIndex[chIdx++] = grpChIdx++;
            if (pPce->SideElementIsCpe[elIdx]) {
              chType[chIdx] = (AUDIO_CHANNEL_TYPE)((spkPlane << 4) | ACT_SIDE);
              chIndex[chIdx++] = grpChIdx++;
            }
          }
        }
        grpChIdx = 0;
        for (elIdx = 0; elIdx < pPce->NumBackChannelElements; elIdx += 1) {
          if (pPce->BackElementHeightInfo[elIdx] == spkPlane) {
            chType[chIdx] = (AUDIO_CHANNEL_TYPE)((spkPlane << 4) | ACT_BACK);
            chIndex[chIdx++] = grpChIdx++;
            if (pPce->BackElementIsCpe[elIdx]) {
              chType[chIdx] = (AUDIO_CHANNEL_TYPE)((spkPlane << 4) | ACT_BACK);
              chIndex[chIdx++] = grpChIdx++;
            }
          }
        }
        grpChIdx = 0;
        if (spkPlane == SPEAKER_PLANE_NORMAL) {
          for (elIdx = 0; elIdx < pPce->NumLfeChannelElements; elIdx += 1) {
            chType[chIdx] = ACT_LFE;
            chIndex[chIdx++] = grpChIdx++;
          }
        }
      }
    }
  } else {
    int chIdx;
    for (chIdx = 0; chIdx < getNumberOfTotalChannels(chConfig); chIdx += 1) {
      getImplicitAudioChannelTypeAndIndex(&chType[chIdx], &chIndex[chIdx],
                                          chConfig, chIdx);
    }
  }
}

int CProgramConfig_GetPceChMap(const CProgramConfig *pPce, UCHAR pceChMap[],
                               const UINT pceChMapLen) {
  const UCHAR *nElements = &pPce->NumFrontChannelElements;
  const UCHAR *elHeight[3], *elIsCpe[3];
  unsigned chIdx, plane, grp, offset, totCh[3], numCh[3][4];

  FDK_ASSERT(pPce != NULL);
  FDK_ASSERT(pceChMap != NULL);

  FDKmemclear(totCh, 3 * sizeof(unsigned));
  FDKmemclear(numCh, 3 * 4 * sizeof(unsigned));

  elHeight[0] = pPce->FrontElementHeightInfo;
  elIsCpe[0] = pPce->FrontElementIsCpe;
  elHeight[1] = pPce->SideElementHeightInfo;
  elIsCpe[1] = pPce->SideElementIsCpe;
  elHeight[2] = pPce->BackElementHeightInfo;
  elIsCpe[2] = pPce->BackElementIsCpe;

  for (plane = 0; plane <= SPEAKER_PLANE_BOTTOM; plane += 1) {
    for (grp = 0; grp < 3; grp += 1) {
      unsigned el;
      for (el = 0; el < nElements[grp]; el += 1) {
        if (elHeight[grp][el] == plane) {
          unsigned elCh = elIsCpe[grp][el] ? 2 : 1;
          numCh[plane][grp] += elCh;
          totCh[plane] += elCh;
        }
      }
    }
    if (plane == SPEAKER_PLANE_NORMAL) {
      unsigned elCh = pPce->NumLfeChannelElements;
      numCh[plane][grp] += elCh;
      totCh[plane] += elCh;
    }
  }
  chIdx = totCh[SPEAKER_PLANE_NORMAL] + totCh[SPEAKER_PLANE_TOP] +
          totCh[SPEAKER_PLANE_BOTTOM];
  if (chIdx > pceChMapLen) {
    return -1;
  }

  offset = grp = 0;
  unsigned grpThresh = numCh[SPEAKER_PLANE_NORMAL][grp];
  for (chIdx = 0; chIdx < totCh[SPEAKER_PLANE_NORMAL]; chIdx += 1) {
    while ((chIdx >= grpThresh) && (grp < 3)) {
      offset += numCh[1][grp] + numCh[2][grp];
      grp += 1;
      grpThresh += numCh[SPEAKER_PLANE_NORMAL][grp];
    }
    pceChMap[chIdx] = chIdx + offset;
  }
  offset = 0;
  for (grp = 0; grp < 4; grp += 1) {
    offset += numCh[SPEAKER_PLANE_NORMAL][grp];
    for (plane = SPEAKER_PLANE_TOP; plane <= SPEAKER_PLANE_BOTTOM; plane += 1) {
      unsigned mapCh;
      for (mapCh = 0; mapCh < numCh[plane][grp]; mapCh += 1) {
        pceChMap[chIdx++] = offset;
        offset += 1;
      }
    }
  }
  return 0;
}

int CProgramConfig_GetElementTable(const CProgramConfig *pPce,
                                   MP4_ELEMENT_ID elList[],
                                   const INT elListSize, UCHAR *pChMapIdx) {
  int i, el = 0;

  FDK_ASSERT(elList != NULL);
  FDK_ASSERT(pChMapIdx != NULL);
  FDK_ASSERT(pPce != NULL);

  *pChMapIdx = 0;

  if ((elListSize <
       pPce->NumFrontChannelElements + pPce->NumSideChannelElements +
           pPce->NumBackChannelElements + pPce->NumLfeChannelElements) ||
      (pPce->NumChannels == 0)) {
    return 0;
  }

  for (i = 0; i < pPce->NumFrontChannelElements; i += 1) {
    elList[el++] = (pPce->FrontElementIsCpe[i]) ? ID_CPE : ID_SCE;
  }

  for (i = 0; i < pPce->NumSideChannelElements; i += 1) {
    elList[el++] = (pPce->SideElementIsCpe[i]) ? ID_CPE : ID_SCE;
  }

  for (i = 0; i < pPce->NumBackChannelElements; i += 1) {
    elList[el++] = (pPce->BackElementIsCpe[i]) ? ID_CPE : ID_SCE;
  }

  for (i = 0; i < pPce->NumLfeChannelElements; i += 1) {
    elList[el++] = ID_LFE;
  }

  switch (pPce->NumChannels) {
    case 1:
    case 2:
      *pChMapIdx = pPce->NumChannels;
      break;
    case 3:
    case 4:
    case 5:
    case 6: {
      C_ALLOC_SCRATCH_START(tmpPce, CProgramConfig, 1);
      CProgramConfig_GetDefault(tmpPce, pPce->NumChannels);
      *pChMapIdx = (!(CProgramConfig_Compare(pPce, tmpPce) & 0xE))
                       ? pPce->NumChannels
                       : 0;
      C_ALLOC_SCRATCH_END(tmpPce, CProgramConfig, 1);
    } break;
    case 7: {
      C_ALLOC_SCRATCH_START(tmpPce, CProgramConfig, 1);
      CProgramConfig_GetDefault(tmpPce, 11);
      *pChMapIdx = (!(CProgramConfig_Compare(pPce, tmpPce) & 0xE)) ? 11 : 0;
      C_ALLOC_SCRATCH_END(tmpPce, CProgramConfig, 1);
    } break;
    case 8: {
      UCHAR testCfg[4] = {32, 14, 12, 7};
      C_ALLOC_SCRATCH_START(tmpPce, CProgramConfig, 1);
      for (i = 0; i < 4; i += 1) {
        CProgramConfig_GetDefault(tmpPce, testCfg[i]);
        if (!(CProgramConfig_Compare(pPce, tmpPce) & 0xE)) {
          *pChMapIdx = (testCfg[i] == 32) ? 12 : testCfg[i];
        }
      }
      C_ALLOC_SCRATCH_END(tmpPce, CProgramConfig, 1);
    } break;
    default:
      *pChMapIdx = 0;
      break;
  }

  return el;
}

static AUDIO_OBJECT_TYPE getAOT(HANDLE_FDK_BITSTREAM bs) {
  int tmp = 0;

  tmp = FDKreadBits(bs, 5);
  if (tmp == AOT_ESCAPE) {
    int tmp2 = FDKreadBits(bs, 6);
    tmp = 32 + tmp2;
  }

  return (AUDIO_OBJECT_TYPE)tmp;
}

static INT getSampleRate(HANDLE_FDK_BITSTREAM bs, UCHAR *index, int nBits) {
  INT sampleRate;
  int idx;

  idx = FDKreadBits(bs, nBits);
  if (idx == (1 << nBits) - 1) {
    if (FDKgetValidBits(bs) < 24) {
      return 0;
    }
    sampleRate = FDKreadBits(bs, 24);
  } else {
    sampleRate = SamplingRateTable[idx];
  }

  *index = idx;

  return sampleRate;
}

static TRANSPORTDEC_ERROR GaSpecificConfig_Parse(CSGaSpecificConfig *self,
                                                 CSAudioSpecificConfig *asc,
                                                 HANDLE_FDK_BITSTREAM bs,
                                                 UINT ascStartAnchor) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;

  self->m_frameLengthFlag = FDKreadBits(bs, 1);

  self->m_dependsOnCoreCoder = FDKreadBits(bs, 1);

  if (self->m_dependsOnCoreCoder) self->m_coreCoderDelay = FDKreadBits(bs, 14);

  self->m_extensionFlag = FDKreadBits(bs, 1);

  if (asc->m_channelConfiguration == 0) {
    CProgramConfig_Read(&asc->m_progrConfigElement, bs, ascStartAnchor);
  }

  if ((asc->m_aot == AOT_AAC_SCAL) || (asc->m_aot == AOT_ER_AAC_SCAL)) {
    self->m_layer = FDKreadBits(bs, 3);
  }

  if (self->m_extensionFlag) {
    if (asc->m_aot == AOT_ER_BSAC) {
      self->m_numOfSubFrame = FDKreadBits(bs, 5);
      self->m_layerLength = FDKreadBits(bs, 11);
    }

    if ((asc->m_aot == AOT_ER_AAC_LC) || (asc->m_aot == AOT_ER_AAC_LTP) ||
        (asc->m_aot == AOT_ER_AAC_SCAL) || (asc->m_aot == AOT_ER_AAC_LD)) {
      asc->m_vcb11Flag = FDKreadBits(bs, 1);
      asc->m_rvlcFlag =
          FDKreadBits(bs, 1);
      asc->m_hcrFlag = FDKreadBits(bs, 1);
    }

    self->m_extensionFlag3 = FDKreadBits(bs, 1);
  }
  return (ErrorStatus);
}

static INT skipSbrHeader(HANDLE_FDK_BITSTREAM hBs, int isUsac) {
  INT dflt_header_extra1, dflt_header_extra2, bitsToSkip = 0;

  if (!isUsac) {
    bitsToSkip = 6;
    FDKpushFor(hBs, 6);
  }
  bitsToSkip += 8;
  FDKpushFor(hBs, 8);
  bitsToSkip += 2;
  dflt_header_extra1 = FDKreadBit(hBs);
  dflt_header_extra2 = FDKreadBit(hBs);
  bitsToSkip += 5 * dflt_header_extra1 + 6 * dflt_header_extra2;
  FDKpushFor(hBs, 5 * dflt_header_extra1 + 6 * dflt_header_extra2);

  return bitsToSkip;
}

static INT ld_sbr_header(CSAudioSpecificConfig *asc, const INT dsFactor,
                         HANDLE_FDK_BITSTREAM hBs, CSTpCallBacks *cb) {
  const int channelConfiguration = asc->m_channelConfiguration;
  int i = 0, j = 0;
  INT error = 0;
  MP4_ELEMENT_ID element = ID_NONE;

  if (channelConfiguration < 0 ||
      channelConfiguration > (INT)(sizeof(channel_configuration_array) /
                                       sizeof(MP4_ELEMENT_ID **) -
                                   1)) {
    return TRANSPORTDEC_PARSE_ERROR;
  }

  while ((element = channel_configuration_array[channelConfiguration][j]) !=
         ID_NONE) {
    if (element == ID_SCE || element == ID_CPE || element == ID_LFE) {
      error |= cb->cbSbr(
          cb->cbSbrData, hBs, asc->m_samplingFrequency / dsFactor,
          asc->m_extensionSamplingFrequency / dsFactor,
          asc->m_samplesPerFrame / dsFactor, AOT_ER_AAC_ELD, element, i++, 0, 0,
          asc->configMode, &asc->SbrConfigChanged, dsFactor);
      if (error != TRANSPORTDEC_OK) {
        goto bail;
      }
    }
    j++;
  }
bail:
  return error;
}

static TRANSPORTDEC_ERROR EldSpecificConfig_Parse(CSAudioSpecificConfig *asc,
                                                  HANDLE_FDK_BITSTREAM hBs,
                                                  CSTpCallBacks *cb) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;
  CSEldSpecificConfig *esc = &asc->m_sc.m_eldSpecificConfig;
  UINT eldExtType;
  int eldExtLen, len, cnt, ldSbrLen = 0, eldExtLenSum, numSbrHeader = 0,
                           sbrIndex, eldExtCnt = 0;

  unsigned char downscale_fill_nibble;

  FDKmemclear(esc, sizeof(CSEldSpecificConfig));

  esc->m_frameLengthFlag = FDKreadBits(hBs, 1);
  if (esc->m_frameLengthFlag) {
    asc->m_samplesPerFrame = 480;
  } else {
    asc->m_samplesPerFrame = 512;
  }

  asc->m_vcb11Flag = FDKreadBits(hBs, 1);
  asc->m_rvlcFlag = FDKreadBits(hBs, 1);
  asc->m_hcrFlag = FDKreadBits(hBs, 1);

  esc->m_sbrPresentFlag = FDKreadBits(hBs, 1);

  if (esc->m_sbrPresentFlag == 1) {
    esc->m_sbrSamplingRate =
        FDKreadBits(hBs, 1);
    esc->m_sbrCrcFlag = FDKreadBits(hBs, 1);

    asc->m_extensionSamplingFrequency = asc->m_samplingFrequency
                                        << esc->m_sbrSamplingRate;

    if (cb->cbSbr != NULL) {
      switch (asc->m_channelConfiguration) {
        case 1:
        case 2:
          numSbrHeader = 1;
          break;
        case 3:
          numSbrHeader = 2;
          break;
        case 4:
        case 5:
        case 6:
          numSbrHeader = 3;
          break;
        case 7:
        case 11:
        case 12:
        case 14:
          numSbrHeader = 4;
          break;
        default:
          numSbrHeader = 0;
          break;
      }
      for (sbrIndex = 0; sbrIndex < numSbrHeader; sbrIndex++) {
        ldSbrLen += skipSbrHeader(hBs, 0);
      }
    } else {
      return TRANSPORTDEC_UNSUPPORTED_FORMAT;
    }
  }
  esc->m_useLdQmfTimeAlign = 0;

  eldExtLenSum = FDKgetValidBits(hBs);
  esc->m_downscaledSamplingFrequency = asc->m_samplingFrequency;
  while (((eldExtType = FDKreadBits(hBs, 4)) != ELDEXT_TERM) &&
         ((INT)FDKgetValidBits(hBs) >= 0) && (eldExtCnt++ < 15)) {
    eldExtLen = len = FDKreadBits(hBs, 4);
    if (len == 0xf) {
      len = FDKreadBits(hBs, 8);
      eldExtLen += len;

      if (len == 0xff) {
        len = FDKreadBits(hBs, 16);
        eldExtLen += len;
      }
    }

    switch (eldExtType) {
      case ELDEXT_LDSAC:
        esc->m_useLdQmfTimeAlign = 1;
        if (cb->cbSsc != NULL) {
          ErrorStatus = (TRANSPORTDEC_ERROR)cb->cbSsc(
              cb->cbSscData, hBs, asc->m_aot,
              asc->m_samplingFrequency << esc->m_sbrSamplingRate,
              asc->m_samplesPerFrame << esc->m_sbrSamplingRate,
              asc->m_channelConfiguration, 1,
              -1,
              eldExtLen, asc->configMode, &asc->SacConfigChanged);
          if (ErrorStatus != TRANSPORTDEC_OK) {
            return TRANSPORTDEC_PARSE_ERROR;
          }
          if (esc->m_downscaledSamplingFrequency != asc->m_samplingFrequency) {
            return TRANSPORTDEC_UNSUPPORTED_FORMAT;
          }
          break;
        }

        FDK_FALLTHROUGH;
      default:
        for (cnt = 0; cnt < eldExtLen; cnt++) {
          FDKreadBits(hBs, 8);
        }
        break;

      case ELDEXT_DOWNSCALEINFO:
        UCHAR tmpDownscaleFreqIdx;
        esc->m_downscaledSamplingFrequency =
            getSampleRate(hBs, &tmpDownscaleFreqIdx, 4);
        if (esc->m_downscaledSamplingFrequency == 0 ||
            esc->m_downscaledSamplingFrequency > 96000) {
          return TRANSPORTDEC_PARSE_ERROR;
        }
        downscale_fill_nibble = FDKreadBits(hBs, 4);
        if (downscale_fill_nibble != 0x0) {
          return TRANSPORTDEC_PARSE_ERROR;
        }
        if (esc->m_useLdQmfTimeAlign == 1) {
          return TRANSPORTDEC_UNSUPPORTED_FORMAT;
        }
        break;
    }
  }
  if (eldExtType != ELDEXT_TERM) {
    return TRANSPORTDEC_PARSE_ERROR;
  }

  if ((INT)FDKgetValidBits(hBs) < 0) {
    return TRANSPORTDEC_PARSE_ERROR;
  }

  if (esc->m_sbrPresentFlag == 1 && numSbrHeader != 0) {
    INT dsFactor = 1;
    if (esc->m_downscaledSamplingFrequency != 0) {
      if (asc->m_samplingFrequency % esc->m_downscaledSamplingFrequency != 0) {
        return TRANSPORTDEC_UNSUPPORTED_FORMAT;
      }
      dsFactor = asc->m_samplingFrequency / esc->m_downscaledSamplingFrequency;
      if (dsFactor != 1 && (dsFactor)&1) {
        return TRANSPORTDEC_UNSUPPORTED_FORMAT;
      }
      if (dsFactor != 1 && dsFactor != 2 && dsFactor != 4) {
        dsFactor = 1;
      }
      if ((INT)asc->m_samplesPerFrame % dsFactor != 0) {
        return TRANSPORTDEC_UNSUPPORTED_FORMAT;
      }
    }
    eldExtLenSum = eldExtLenSum - FDKgetValidBits(hBs);
    FDKpushBack(hBs, eldExtLenSum + ldSbrLen);
    if (0 != ld_sbr_header(asc, dsFactor, hBs, cb)) {
      return TRANSPORTDEC_PARSE_ERROR;
    }
    FDKpushFor(hBs, eldExtLenSum);
  }
  return (ErrorStatus);
}

static UINT StoreConfigAsBitstream(
    HANDLE_FDK_BITSTREAM hBs, const INT configSize_bits,
    UCHAR *configTargetBuffer, const USHORT configTargetBufferSize_bytes) {
  FDK_BITSTREAM usacConf;
  UINT const nBits = fAbs(configSize_bits);
  UINT j, tmp;

  if (nBits > 8 * (UINT)configTargetBufferSize_bytes) {
    return 1;
  }
  FDKmemclear(configTargetBuffer, configTargetBufferSize_bytes);

  FDKinitBitStream(&usacConf, configTargetBuffer, configTargetBufferSize_bytes,
                   nBits, BS_WRITER);
  if (configSize_bits < 0) {
    FDKpushBack(hBs, nBits);
  }
  for (j = nBits; j > 31; j -= 32) {
    tmp = FDKreadBits(hBs, 32);
    FDKwriteBits(&usacConf, tmp, 32);
  }
  if (j > 0) {
    tmp = FDKreadBits(hBs, j);
    FDKwriteBits(&usacConf, tmp, j);
  }
  FDKsyncCache(&usacConf);
  if (configSize_bits > 0) {
    FDKpushBack(hBs, nBits);
  }

  return 0;
}

static const USHORT usacFrameLength[8] = {768, 1024, 2048, 2048, 4096, 0, 0, 0};
static const UCHAR sbrRatioIndex[8] = {0, 0, 2, 3, 1, 0, 0, 0};

static TRANSPORTDEC_ERROR extElementConfig(CSUsacExtElementConfig *extElement,
                                           HANDLE_FDK_BITSTREAM hBs,
                                           const CSTpCallBacks *cb,
                                           const UCHAR numSignalsInGroup,
                                           const UINT coreFrameLength,
                                           const int subStreamIndex,
                                           const AUDIO_OBJECT_TYPE aot) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;
  UINT usacExtElementType = escapedValue(hBs, 4, 8, 16);

  if (aot == AOT_USAC) {
    switch (usacExtElementType) {
      case ID_EXT_ELE_FILL:
      case ID_EXT_ELE_MPEGS:
      case ID_EXT_ELE_SAOC:
      case ID_EXT_ELE_AUDIOPREROLL:
      case ID_EXT_ELE_UNI_DRC:
        break;
      default:
        usacExtElementType = ID_EXT_ELE_UNKNOWN;
        break;
    }
  }

  int usacExtElementConfigLength = escapedValue(hBs, 4, 8, 16);
  extElement->usacExtElementConfigLength = (USHORT)usacExtElementConfigLength;
  INT bsAnchor;

  if (FDKreadBit(hBs))
    extElement->usacExtElementDefaultLength = escapedValue(hBs, 8, 16, 0) + 1;
  else
    extElement->usacExtElementDefaultLength = 0;

  extElement->usacExtElementPayloadFrag = FDKreadBit(hBs);

  bsAnchor = (INT)FDKgetValidBits(hBs);

  switch (usacExtElementType) {
    case ID_EXT_ELE_UNKNOWN:
    case ID_EXT_ELE_FILL:
      break;
    case ID_EXT_ELE_AUDIOPREROLL:
      extElement->usacExtElementHasAudioPreRoll = 1;
      break;
    case ID_EXT_ELE_UNI_DRC: {
      if (cb->cbUniDrc != NULL) {
        ErrorStatus = (TRANSPORTDEC_ERROR)cb->cbUniDrc(
            cb->cbUniDrcData, hBs, usacExtElementConfigLength,
            0,
            subStreamIndex, 0, aot);
        if (ErrorStatus != TRANSPORTDEC_OK) {
          return ErrorStatus;
        }
      }
    } break;
    default:
      usacExtElementType = ID_EXT_ELE_UNKNOWN;
      break;
  }
  extElement->usacExtElementType = (USAC_EXT_ELEMENT_TYPE)usacExtElementType;

  {
    INT left_bits = (usacExtElementConfigLength << 3) -
                    (bsAnchor - (INT)FDKgetValidBits(hBs));
    if (left_bits >= 0) {
      FDKpushFor(hBs, left_bits);
    } else {
      ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
    }
  }

  return ErrorStatus;
}

static TRANSPORTDEC_ERROR configExtension(CSUsacConfig *usc,
                                          HANDLE_FDK_BITSTREAM hBs,
                                          const CSTpCallBacks *cb) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;

  int numConfigExtensions;
  UINT usacConfigExtType;
  int usacConfigExtLength;
  int loudnessInfoSetIndex =
      -1;
  int tmp_subStreamIndex = 0;
  AUDIO_OBJECT_TYPE tmp_aot = AOT_USAC;

  numConfigExtensions = (int)escapedValue(hBs, 2, 4, 8) + 1;
  for (int confExtIdx = 0; confExtIdx < numConfigExtensions; confExtIdx++) {
    INT nbits;
    int loudnessInfoSetConfigExtensionPosition = FDKgetValidBits(hBs);
    usacConfigExtType = escapedValue(hBs, 4, 8, 16);
    usacConfigExtLength = (int)escapedValue(hBs, 4, 8, 16);

    nbits = (INT)FDKgetValidBits(hBs);

    if (nbits < usacConfigExtLength * 8) {
      return TRANSPORTDEC_PARSE_ERROR;
    }

    switch (usacConfigExtType) {
      case ID_CONFIG_EXT_FILL:
        for (int i = 0; i < usacConfigExtLength; i++) {
          if (FDKreadBits(hBs, 8) != 0xa5) {
            return TRANSPORTDEC_PARSE_ERROR;
          }
        }
        break;
      case ID_CONFIG_EXT_LOUDNESS_INFO: {
        if (cb->cbUniDrc != NULL) {
          ErrorStatus = (TRANSPORTDEC_ERROR)cb->cbUniDrc(
              cb->cbUniDrcData, hBs, usacConfigExtLength,
              1,
              tmp_subStreamIndex, loudnessInfoSetConfigExtensionPosition,
              tmp_aot);
          if (ErrorStatus != TRANSPORTDEC_OK) {
            return ErrorStatus;
          }
          loudnessInfoSetIndex = confExtIdx;
        }
      } break;
      default:
        break;
    }

    usacConfigExtLength =
        8 * usacConfigExtLength - (nbits - (INT)FDKgetValidBits(hBs));
    if (usacConfigExtLength < 0) {
      return TRANSPORTDEC_PARSE_ERROR;
    }
    FDKpushFor(hBs, usacConfigExtLength);
  }

  if (loudnessInfoSetIndex == -1 && cb->cbUniDrc != NULL) {
    ErrorStatus = (TRANSPORTDEC_ERROR)cb->cbUniDrc(
        cb->cbUniDrcData, NULL, 0, 1, tmp_subStreamIndex,
        0, tmp_aot);
    if (ErrorStatus != TRANSPORTDEC_OK) {
      return ErrorStatus;
    }
  }

  return ErrorStatus;
}

static TRANSPORTDEC_ERROR UsacRsv60DecoderConfig_Parse(
    CSAudioSpecificConfig *asc, HANDLE_FDK_BITSTREAM hBs,
    const CSTpCallBacks *cb) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;
  CSUsacConfig *usc = &asc->m_sc.m_usacConfig;
  int i, numberOfElements;
  int channelElementIdx =
      0;
  SC_CHANNEL_CONFIG sc_chan_config = {0, 0, 0, 0};
  int uniDrcElement =
      -1;

  numberOfElements = (int)escapedValue(hBs, 4, 8, 16) + 1;
  usc->m_usacNumElements = numberOfElements;
  if (numberOfElements > TP_USAC_MAX_ELEMENTS) {
    return TRANSPORTDEC_UNSUPPORTED_FORMAT;
  }
  usc->m_nUsacChannels = 0;
  usc->m_channelConfigurationIndex = asc->m_channelConfiguration;

  if (asc->m_aot == AOT_USAC) {
    sc_chan_config = sc_chan_config_tab[usc->m_channelConfigurationIndex];

    if (sc_chan_config.nCh > (SCHAR)TP_USAC_MAX_SPEAKERS) {
      return TRANSPORTDEC_PARSE_ERROR;
    }
  }

  for (i = 0; i < numberOfElements; i++) {
    MP4_ELEMENT_ID usacElementType = (MP4_ELEMENT_ID)(
        FDKreadBits(hBs, 2) | USAC_ID_BIT);
    usc->element[i].usacElementType = usacElementType;

    if (asc->m_aot == AOT_USAC) {
      switch (usacElementType) {
        case ID_USAC_SCE:
          sc_chan_config.nSCE--;
          break;
        case ID_USAC_CPE:
          sc_chan_config.nCPE--;
          break;
        case ID_USAC_LFE:
          sc_chan_config.nLFE--;
          break;
        default:
          break;
      }
      if (usc->m_channelConfigurationIndex) {
        if (sc_chan_config.nCPE < 0 || sc_chan_config.nSCE < 0 ||
            sc_chan_config.nLFE < 0) {
          return TRANSPORTDEC_PARSE_ERROR;
        }
      }
    }

    switch (usacElementType) {
      case ID_USAC_SCE:
        if (FDKreadBit(hBs)) {
          return TRANSPORTDEC_UNSUPPORTED_FORMAT;
        }
        usc->element[i].m_noiseFilling = FDKreadBits(hBs, 1);
        if (usc->m_sbrRatioIndex > 0) {
          if (cb->cbSbr == NULL) {
            return TRANSPORTDEC_UNKOWN_ERROR;
          }
          usc->element[i].m_harmonicSBR = FDKreadBit(hBs);
          usc->element[i].m_interTes = FDKreadBit(hBs);
          usc->element[i].m_pvc = FDKreadBit(hBs);
          if (cb->cbSbr(cb->cbSbrData, hBs, asc->m_samplingFrequency,
                        asc->m_extensionSamplingFrequency,
                        asc->m_samplesPerFrame, asc->m_aot, ID_SCE,
                        channelElementIdx, usc->element[i].m_harmonicSBR,
                        usc->element[i].m_stereoConfigIndex, asc->configMode,
                        &asc->SbrConfigChanged, 1)) {
            return TRANSPORTDEC_PARSE_ERROR;
          }
        }
        usc->m_nUsacChannels += 1;
        channelElementIdx++;
        break;

      case ID_USAC_CPE:
        if (FDKreadBit(hBs)) {
          return TRANSPORTDEC_UNSUPPORTED_FORMAT;
        }
        usc->element[i].m_noiseFilling = FDKreadBits(hBs, 1);
        if (usc->m_sbrRatioIndex > 0) {
          if (cb->cbSbr == NULL) return TRANSPORTDEC_UNKOWN_ERROR;
          usc->element[i].m_harmonicSBR = FDKreadBit(hBs);
          usc->element[i].m_interTes = FDKreadBit(hBs);
          usc->element[i].m_pvc = FDKreadBit(hBs);
          {
            INT bitsToSkip = skipSbrHeader(hBs, 1);
            usc->element[i].m_stereoConfigIndex = FDKreadBits(hBs, 2);
            FDKpushBack(hBs, bitsToSkip + 2);
          }
          {
            MP4_ELEMENT_ID el_type =
                (usc->element[i].m_stereoConfigIndex == 1 ||
                 usc->element[i].m_stereoConfigIndex == 2)
                    ? ID_SCE
                    : ID_CPE;
            if (cb->cbSbr(cb->cbSbrData, hBs, asc->m_samplingFrequency,
                          asc->m_extensionSamplingFrequency,
                          asc->m_samplesPerFrame, asc->m_aot, el_type,
                          channelElementIdx, usc->element[i].m_harmonicSBR,
                          usc->element[i].m_stereoConfigIndex, asc->configMode,
                          &asc->SbrConfigChanged, 1)) {
              return TRANSPORTDEC_PARSE_ERROR;
            }
          }

          usc->element[i].m_stereoConfigIndex =
              FDKreadBits(hBs, 2);

          if (usc->element[i].m_stereoConfigIndex > 0) {
            if (cb->cbSsc != NULL) {
              int samplesPerFrame = asc->m_samplesPerFrame;

              if (usc->m_sbrRatioIndex == 1) samplesPerFrame <<= 2;
              if (usc->m_sbrRatioIndex == 2)
                samplesPerFrame = (samplesPerFrame * 8) / 3;
              if (usc->m_sbrRatioIndex == 3) samplesPerFrame <<= 1;

              if (cb->cbSsc(cb->cbSscData, hBs, asc->m_aot,
                            asc->m_extensionSamplingFrequency, samplesPerFrame,
                            1,
                            usc->element[i].m_stereoConfigIndex,
                            usc->m_coreSbrFrameLengthIndex,
                            0,
                            asc->configMode, &asc->SacConfigChanged)) {
                return TRANSPORTDEC_PARSE_ERROR;
              }
            } else {
              return TRANSPORTDEC_UNKOWN_ERROR;
            }
          }
        } else {
          usc->element[i].m_stereoConfigIndex = 0;
        }
        usc->m_nUsacChannels += 2;

        channelElementIdx++;
        break;

      case ID_USAC_LFE:
        usc->element[i].m_noiseFilling = 0;
        usc->m_nUsacChannels += 1;
        if (usc->m_sbrRatioIndex > 0) {
          if (cb->cbSbr == NULL) return ErrorStatus = TRANSPORTDEC_UNKOWN_ERROR;
          usc->element[i].m_harmonicSBR = (UCHAR)0;
          usc->element[i].m_interTes = (UCHAR)0;
          usc->element[i].m_pvc = (UCHAR)0;
          if (cb->cbSbr(cb->cbSbrData, hBs, asc->m_samplingFrequency,
                        asc->m_extensionSamplingFrequency,
                        asc->m_samplesPerFrame, asc->m_aot, ID_LFE,
                        channelElementIdx, usc->element[i].m_harmonicSBR,
                        usc->element[i].m_stereoConfigIndex, asc->configMode,
                        &asc->SbrConfigChanged, 1)) {
            return ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
          }
        }
        channelElementIdx++;
        break;

      case ID_USAC_EXT:
        ErrorStatus = extElementConfig(&usc->element[i].extElement, hBs, cb, 0,
                                       asc->m_samplesPerFrame, 0, asc->m_aot);
        if (usc->element[i].extElement.usacExtElementType ==
            ID_EXT_ELE_UNI_DRC) {
          uniDrcElement = i;
        }

        if (ErrorStatus) {
          return ErrorStatus;
        }
        break;

      default:
        return TRANSPORTDEC_PARSE_ERROR;
    }
  }

  if (asc->m_aot == AOT_USAC) {
    if (usc->m_channelConfigurationIndex) {
      if (sc_chan_config.nCPE | sc_chan_config.nSCE | sc_chan_config.nLFE) {
        return TRANSPORTDEC_PARSE_ERROR;
      }
    } else {
      if ((INT)(-2 * sc_chan_config.nCPE - sc_chan_config.nSCE -
                sc_chan_config.nLFE) < (INT)usc->numAudioChannels) {
        return TRANSPORTDEC_PARSE_ERROR;
      }
    }
  }

  if (uniDrcElement == -1 && cb->cbUniDrc != NULL) {
    int subStreamIndex = 0;
    ErrorStatus = (TRANSPORTDEC_ERROR)cb->cbUniDrc(
        cb->cbUniDrcData, NULL, 0, 0, subStreamIndex, 0,
        asc->m_aot);
    if (ErrorStatus != TRANSPORTDEC_OK) {
      return ErrorStatus;
    }
  }

  return ErrorStatus;
}

static TRANSPORTDEC_ERROR UsacConfig_SetCoreSbrFrameLengthIndex(
    CSAudioSpecificConfig *asc, int coreSbrFrameLengthIndex) {
  int sbrRatioIndex_val;

  if (coreSbrFrameLengthIndex > 4) {
    return TRANSPORTDEC_PARSE_ERROR;
  }
  asc->m_sc.m_usacConfig.m_coreSbrFrameLengthIndex = coreSbrFrameLengthIndex;
  asc->m_samplesPerFrame = usacFrameLength[coreSbrFrameLengthIndex];
  sbrRatioIndex_val = sbrRatioIndex[coreSbrFrameLengthIndex];
  asc->m_sc.m_usacConfig.m_sbrRatioIndex = sbrRatioIndex_val;

  if (sbrRatioIndex_val > 0) {
    asc->m_sbrPresentFlag = 1;
    asc->m_extensionSamplingFrequency = asc->m_samplingFrequency;
    asc->m_extensionSamplingFrequencyIndex = asc->m_samplingFrequencyIndex;
    switch (sbrRatioIndex_val) {
      case 1:
        asc->m_samplingFrequency >>= 2;
        asc->m_samplesPerFrame >>= 2;
        break;
      case 2:
        asc->m_samplingFrequency = (asc->m_samplingFrequency * 3) / 8;
        asc->m_samplesPerFrame = (asc->m_samplesPerFrame * 3) / 8;
        break;
      case 3:
        asc->m_samplingFrequency >>= 1;
        asc->m_samplesPerFrame >>= 1;
        break;
      default:
        return TRANSPORTDEC_PARSE_ERROR;
    }
    asc->m_samplingFrequencyIndex =
        getSamplingRateIndex(asc->m_samplingFrequency, 4);
  }

  return TRANSPORTDEC_OK;
}

static TRANSPORTDEC_ERROR UsacConfig_Parse(CSAudioSpecificConfig *asc,
                                           HANDLE_FDK_BITSTREAM hBs,
                                           CSTpCallBacks *cb) {
  int usacSamplingFrequency, channelConfigurationIndex, coreSbrFrameLengthIndex;
  TRANSPORTDEC_ERROR err = TRANSPORTDEC_OK;

  INT nbits = (INT)FDKgetValidBits(hBs);

  usacSamplingFrequency = getSampleRate(hBs, &asc->m_samplingFrequencyIndex, 5);
  if (usacSamplingFrequency == 0 || usacSamplingFrequency > 96000) {
    return TRANSPORTDEC_PARSE_ERROR;
  }
  asc->m_samplingFrequency = (UINT)usacSamplingFrequency;

  coreSbrFrameLengthIndex = FDKreadBits(hBs, 3);
  if (UsacConfig_SetCoreSbrFrameLengthIndex(asc, coreSbrFrameLengthIndex) !=
      TRANSPORTDEC_OK) {
    return TRANSPORTDEC_PARSE_ERROR;
  }

  channelConfigurationIndex = FDKreadBits(hBs, 5);
  if (channelConfigurationIndex > 2) {
    return TRANSPORTDEC_PARSE_ERROR;
  }

  if (channelConfigurationIndex == 0) {
    return TRANSPORTDEC_PARSE_ERROR;
  }
  asc->m_channelConfiguration = channelConfigurationIndex;

  err = UsacRsv60DecoderConfig_Parse(asc, hBs, cb);
  if (err != TRANSPORTDEC_OK) {
    return err;
  }

  if (FDKreadBits(hBs, 1)) {
    err = configExtension(&asc->m_sc.m_usacConfig, hBs, cb);
    if (err != TRANSPORTDEC_OK) {
      return err;
    }
  } else if (cb->cbUniDrc != NULL) {
    err = (TRANSPORTDEC_ERROR)cb->cbUniDrc(
        cb->cbUniDrcData, NULL, 0, 1, 0, 0, asc->m_aot);
    if (err != TRANSPORTDEC_OK) {
      return err;
    }
  }

  if ((channelConfigurationIndex > 0) &&
      (sc_chan_config_tab[channelConfigurationIndex].nCh !=
       asc->m_sc.m_usacConfig.m_nUsacChannels)) {
    return TRANSPORTDEC_PARSE_ERROR;
  }

  INT configSize_bits = (INT)FDKgetValidBits(hBs) - nbits;
  if (StoreConfigAsBitstream(hBs, configSize_bits,
                             asc->m_sc.m_usacConfig.UsacConfig,
                             TP_USAC_MAX_CONFIG_LEN)) {
    return TRANSPORTDEC_PARSE_ERROR;
  }
  asc->m_sc.m_usacConfig.UsacConfigBits = fAbs(configSize_bits);

  return err;
}

static TRANSPORTDEC_ERROR AudioSpecificConfig_ExtensionParse(
    CSAudioSpecificConfig *self, HANDLE_FDK_BITSTREAM bs, CSTpCallBacks *cb) {
  TP_ASC_EXTENSION_ID lastAscExt, ascExtId = ASCEXT_UNKOWN;
  INT bitsAvailable = (INT)FDKgetValidBits(bs);

  while (bitsAvailable >= 11) {
    lastAscExt = ascExtId;
    ascExtId = (TP_ASC_EXTENSION_ID)FDKreadBits(bs, 11);
    bitsAvailable -= 11;

    switch (ascExtId) {
      case ASCEXT_SBR:
        if ((self->m_extensionAudioObjectType != AOT_SBR) &&
            (bitsAvailable >= 5)) {
          self->m_extensionAudioObjectType = getAOT(bs);

          if ((self->m_extensionAudioObjectType == AOT_SBR) ||
              (self->m_extensionAudioObjectType ==
               AOT_ER_BSAC)) {
            self->m_sbrPresentFlag = FDKreadBits(bs, 1);
            if (self->m_aot == AOT_USAC && self->m_sbrPresentFlag > 0 &&
                self->m_sc.m_usacConfig.m_sbrRatioIndex == 0) {
              return TRANSPORTDEC_PARSE_ERROR;
            }

            if (self->m_sbrPresentFlag == 1) {
              self->m_extensionSamplingFrequency = getSampleRate(
                  bs, &self->m_extensionSamplingFrequencyIndex, 4);

              if (self->m_extensionSamplingFrequency == 0 ||
                  self->m_extensionSamplingFrequency > 96000) {
                return TRANSPORTDEC_PARSE_ERROR;
              }
            }
            if (self->m_extensionAudioObjectType == AOT_ER_BSAC) {
              self->m_extensionChannelConfiguration = FDKreadBits(bs, 4);
            }
          }
          bitsAvailable = (INT)FDKgetValidBits(bs);
        }
        break;
      case ASCEXT_PS:
        if ((lastAscExt == ASCEXT_SBR) &&
            (self->m_extensionAudioObjectType == AOT_SBR) &&
            (bitsAvailable > 0)) {
          self->m_psPresentFlag = FDKreadBits(bs, 1);
          bitsAvailable -= 1;
        }
        break;
      case ASCEXT_MPS:
        if (self->m_extensionAudioObjectType == AOT_MPEGS) break;
        FDK_FALLTHROUGH;
      case ASCEXT_LDMPS:
        if ((ascExtId == ASCEXT_LDMPS) &&
            (self->m_extensionAudioObjectType == AOT_LD_MPEGS))
          break;
        if (bitsAvailable >= 1) {
          bitsAvailable -= 1;
          if (FDKreadBits(bs, 1)) {
            int sscLen = FDKreadBits(bs, 8);
            bitsAvailable -= 8;
            if (sscLen == 0xFF) {
              sscLen += FDKreadBits(bs, 16);
              bitsAvailable -= 16;
            }
            FDKpushFor(bs, sscLen);
            bitsAvailable -= sscLen * 8;
          }
        }
        break;
      case ASCEXT_SAOC:
        if ((ascExtId == ASCEXT_SAOC) &&
            (self->m_extensionAudioObjectType == AOT_SAOC))
          break;
        if (FDKreadBits(bs, 1)) {
          int saocscLen = FDKreadBits(bs, 8);
          bitsAvailable -= 8;
          if (saocscLen == 0xFF) {
            saocscLen += FDKreadBits(bs, 16);
            bitsAvailable -= 16;
          }
          FDKpushFor(bs, saocscLen);
          bitsAvailable -= saocscLen * 8;
        }
        break;
      default:
        return TRANSPORTDEC_OK;
    }
  }

  return TRANSPORTDEC_OK;
}

void AudioSpecificConfig_Init(CSAudioSpecificConfig *asc) {
  FDKmemclear(asc, sizeof(CSAudioSpecificConfig));

  asc->m_aot = AOT_NONE;
  asc->m_samplingFrequencyIndex = 0xf;
  asc->m_epConfig = -1;
  asc->m_extensionAudioObjectType = AOT_NULL_OBJECT;
  CProgramConfig_Init(&asc->m_progrConfigElement);
}

TRANSPORTDEC_ERROR AudioSpecificConfig_Parse(
    CSAudioSpecificConfig *self, HANDLE_FDK_BITSTREAM bs,
    int fExplicitBackwardCompatible, CSTpCallBacks *cb, UCHAR configMode,
    UCHAR configChanged, AUDIO_OBJECT_TYPE m_aot) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;
  UINT ascStartAnchor = FDKgetValidBits(bs);
  int frameLengthFlag = -1;

  AudioSpecificConfig_Init(self);

  self->configMode = configMode;
  self->AacConfigChanged = configChanged;
  self->SbrConfigChanged = configChanged;
  self->SacConfigChanged = configChanged;

  if (m_aot != AOT_NULL_OBJECT) {
    self->m_aot = m_aot;
  } else {
    self->m_aot = getAOT(bs);
    self->m_samplingFrequency =
        getSampleRate(bs, &self->m_samplingFrequencyIndex, 4);
    if (self->m_samplingFrequency <= 0 ||
        (self->m_samplingFrequency > 96000 && self->m_aot != 39) ||
        self->m_samplingFrequency > 4 * 96000) {
      return TRANSPORTDEC_PARSE_ERROR;
    }

    self->m_channelConfiguration = FDKreadBits(bs, 4);

    if ((self->m_channelConfiguration == 0) &&
        ((self->m_aot == AOT_ER_AAC_LC) || (self->m_aot == AOT_ER_AAC_LTP) ||
         (self->m_aot == AOT_ER_AAC_LD) || (self->m_aot == AOT_ER_AAC_SCAL) ||
         (self->m_aot == AOT_ER_AAC_ELD))) {
      return TRANSPORTDEC_UNSUPPORTED_FORMAT;
    }
    if ((self->m_channelConfiguration > 2) &&
        ((self->m_aot == AOT_AAC_SCAL) || (self->m_aot == AOT_ER_AAC_SCAL))) {
      return TRANSPORTDEC_UNSUPPORTED_FORMAT;
    }

    self->m_sbrPresentFlag = 0;
    self->m_psPresentFlag = 0;

    if (self->m_aot == AOT_SBR || self->m_aot == AOT_PS) {
      self->m_extensionAudioObjectType = AOT_SBR;

      self->m_sbrPresentFlag = 1;
      if (self->m_aot == AOT_PS) {
        self->m_psPresentFlag = 1;
      }

      self->m_extensionSamplingFrequency =
          getSampleRate(bs, &self->m_extensionSamplingFrequencyIndex, 4);
      if (self->m_extensionSamplingFrequency == 0 ||
          self->m_extensionSamplingFrequency > 96000) {
        return TRANSPORTDEC_PARSE_ERROR;
      }
      self->m_aot = getAOT(bs);

      switch (self->m_aot) {
        case AOT_AAC_LC:
          break;
        case AOT_ER_BSAC:
          break;
        default:
          return TRANSPORTDEC_UNSUPPORTED_FORMAT;
      }

      if (self->m_aot == AOT_ER_BSAC) {
        self->m_extensionChannelConfiguration = FDKreadBits(bs, 4);
      }
    } else {
      self->m_extensionAudioObjectType = AOT_NULL_OBJECT;
    }
  }

  switch (self->m_aot) {
    case AOT_AAC_LC:
    case AOT_AAC_SCAL:
    case AOT_ER_AAC_LC:
    case AOT_ER_AAC_LD:
    case AOT_ER_AAC_SCAL:
    case AOT_ER_BSAC:
      if ((ErrorStatus = GaSpecificConfig_Parse(&self->m_sc.m_gaSpecificConfig,
                                                self, bs, ascStartAnchor)) !=
          TRANSPORTDEC_OK) {
        return (ErrorStatus);
      }
      frameLengthFlag = self->m_sc.m_gaSpecificConfig.m_frameLengthFlag;
      break;
    case AOT_MPEGS:
      if (cb->cbSsc != NULL) {
        if (cb->cbSsc(cb->cbSscData, bs, self->m_aot, self->m_samplingFrequency,
                      self->m_samplesPerFrame, self->m_channelConfiguration, 1,
                      -1,
                      0,
                      self->configMode, &self->SacConfigChanged)) {
          return TRANSPORTDEC_UNSUPPORTED_FORMAT;
        }
      } else {
        return TRANSPORTDEC_UNSUPPORTED_FORMAT;
      }
      break;
    case AOT_ER_AAC_ELD:
      if ((ErrorStatus = EldSpecificConfig_Parse(self, bs, cb)) !=
          TRANSPORTDEC_OK) {
        return (ErrorStatus);
      }
      frameLengthFlag = self->m_sc.m_eldSpecificConfig.m_frameLengthFlag;
      self->m_sbrPresentFlag = self->m_sc.m_eldSpecificConfig.m_sbrPresentFlag;
      self->m_extensionSamplingFrequency =
          (self->m_sc.m_eldSpecificConfig.m_sbrSamplingRate + 1) *
          self->m_samplingFrequency;
      break;
    case AOT_USAC:
      if ((ErrorStatus = UsacConfig_Parse(self, bs, cb)) != TRANSPORTDEC_OK) {
        return (ErrorStatus);
      }
      break;

    default:
      return TRANSPORTDEC_UNSUPPORTED_FORMAT;
  }

  switch (self->m_aot) {
    case AOT_AAC_LC:
    case AOT_AAC_SCAL:
    case AOT_ER_AAC_LC:
    case AOT_ER_AAC_SCAL:
    case AOT_ER_BSAC:
      if (!frameLengthFlag)
        self->m_samplesPerFrame = 1024;
      else
        self->m_samplesPerFrame = 960;
      break;
    case AOT_ER_AAC_LD:
      if (!frameLengthFlag)
        self->m_samplesPerFrame = 512;
      else
        self->m_samplesPerFrame = 480;
      break;
    default:
      break;
  }

  switch (self->m_aot) {
    case AOT_ER_AAC_LC:
    case AOT_ER_AAC_LD:
    case AOT_ER_AAC_ELD:
    case AOT_ER_AAC_SCAL:
    case AOT_ER_CELP:
    case AOT_ER_HVXC:
    case AOT_ER_BSAC:
      self->m_epConfig = FDKreadBits(bs, 2);

      if (self->m_epConfig > 1) {
        return TRANSPORTDEC_UNSUPPORTED_FORMAT;
      }
      break;
    default:
      break;
  }

  if (fExplicitBackwardCompatible &&
      (self->m_aot == AOT_AAC_LC || self->m_aot == AOT_ER_AAC_LD ||
       self->m_aot == AOT_ER_BSAC)) {
    ErrorStatus = AudioSpecificConfig_ExtensionParse(self, bs, cb);
  }

  if ((ErrorStatus == TRANSPORTDEC_OK) && (self->m_aot == AOT_USAC)) {
    INT configSize_bits = (INT)FDKgetValidBits(bs) - (INT)ascStartAnchor;
    if (StoreConfigAsBitstream(bs, configSize_bits, self->config,
                               TP_USAC_MAX_CONFIG_LEN)) {
      return TRANSPORTDEC_PARSE_ERROR;
    }
    self->configBits = fAbs(configSize_bits);
  }

  return (ErrorStatus);
}

static TRANSPORTDEC_ERROR Drm_xHEAACDecoderConfig(
    CSAudioSpecificConfig *asc, HANDLE_FDK_BITSTREAM hBs, int audioMode,
    CSTpCallBacks *cb
) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;
  CSUsacConfig *usc = &asc->m_sc.m_usacConfig;
  int elemIdx = 0;

  usc->element[elemIdx].m_stereoConfigIndex = 0;

  usc->m_usacNumElements = 1;

  switch (audioMode) {
    case 0:
      usc->element[elemIdx].usacElementType = ID_USAC_SCE;
      usc->m_nUsacChannels = 1;
      usc->element[elemIdx].m_noiseFilling = FDKreadBits(hBs, 1);
      if (usc->m_sbrRatioIndex > 0) {
        if (cb == NULL) {
          return ErrorStatus;
        }
        if (cb->cbSbr != NULL) {
          usc->element[elemIdx].m_harmonicSBR = FDKreadBit(hBs);
          usc->element[elemIdx].m_interTes = FDKreadBit(hBs);
          usc->element[elemIdx].m_pvc = FDKreadBit(hBs);
          if (cb->cbSbr(cb->cbSbrData, hBs, asc->m_samplingFrequency,
                        asc->m_extensionSamplingFrequency,
                        asc->m_samplesPerFrame, asc->m_aot, ID_SCE, elemIdx,
                        usc->element[elemIdx].m_harmonicSBR,
                        usc->element[elemIdx].m_stereoConfigIndex,
                        asc->configMode, &asc->SbrConfigChanged, 1)) {
            return ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
          }
        }
      }
      break;
    case 2:
      usc->element[elemIdx].usacElementType = ID_USAC_CPE;
      usc->m_nUsacChannels = 2;
      usc->element[elemIdx].m_noiseFilling = FDKreadBits(hBs, 1);
      if (usc->m_sbrRatioIndex > 0) {
        usc->element[elemIdx].m_harmonicSBR = FDKreadBit(hBs);
        usc->element[elemIdx].m_interTes = FDKreadBit(hBs);
        usc->element[elemIdx].m_pvc = FDKreadBit(hBs);
        {
          INT bitsToSkip = skipSbrHeader(hBs, 1);
          usc->element[elemIdx].m_stereoConfigIndex = FDKreadBits(hBs, 2);
          FDKpushBack(hBs, bitsToSkip + 2);
        }
        if ((usc->element[elemIdx].m_stereoConfigIndex > 1) &&
            usc->element[elemIdx].m_harmonicSBR) {
          return ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
        }
        if ((usc->m_sbrRatioIndex == 1) &&
            (usc->element[elemIdx].m_stereoConfigIndex != 1)) {
          return ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
        }
        if (cb == NULL) {
          return ErrorStatus;
        }
        {
          MP4_ELEMENT_ID el_type =
              (usc->element[elemIdx].m_stereoConfigIndex == 1 ||
               usc->element[elemIdx].m_stereoConfigIndex == 2)
                  ? ID_SCE
                  : ID_CPE;
          if (cb->cbSbr == NULL) return ErrorStatus = TRANSPORTDEC_UNKOWN_ERROR;
          if (cb->cbSbr(cb->cbSbrData, hBs, asc->m_samplingFrequency,
                        asc->m_extensionSamplingFrequency,
                        asc->m_samplesPerFrame, asc->m_aot, el_type, elemIdx,
                        usc->element[elemIdx].m_harmonicSBR,
                        usc->element[elemIdx].m_stereoConfigIndex,
                        asc->configMode, &asc->SbrConfigChanged, 1)) {
            return ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
          }
        }
        FDKreadBits(hBs, 2);
        if (usc->element[elemIdx].m_stereoConfigIndex > 0) {
          if (cb->cbSsc != NULL) {
            int samplesPerFrame = asc->m_samplesPerFrame;

            if (usc->m_sbrRatioIndex == 1) samplesPerFrame <<= 2;
            if (usc->m_sbrRatioIndex == 2)
              samplesPerFrame = (samplesPerFrame * 8) / 3;
            if (usc->m_sbrRatioIndex == 3) samplesPerFrame <<= 1;

            ErrorStatus = (TRANSPORTDEC_ERROR)cb->cbSsc(
                cb->cbSscData, hBs,
                AOT_DRM_USAC,
                asc->m_extensionSamplingFrequency, samplesPerFrame,
                1,
                usc->element[elemIdx].m_stereoConfigIndex,
                usc->m_coreSbrFrameLengthIndex, 0,
                asc->configMode, &asc->SacConfigChanged);
          } else {
          }
        }
      }
      break;
    default:
      return TRANSPORTDEC_PARSE_ERROR;
  }

  return ErrorStatus;
}

TRANSPORTDEC_ERROR Drm_xHEAACStaticConfig(
    CSAudioSpecificConfig *asc, HANDLE_FDK_BITSTREAM bs, int audioMode,
    CSTpCallBacks *cb
) {
  int coreSbrFrameLengthIndexDrm = FDKreadBits(bs, 2);
  if (UsacConfig_SetCoreSbrFrameLengthIndex(
          asc, coreSbrFrameLengthIndexDrm + 1) != TRANSPORTDEC_OK) {
    return TRANSPORTDEC_PARSE_ERROR;
  }

  asc->m_channelConfiguration = (audioMode) ? 2 : 1;

  if (Drm_xHEAACDecoderConfig(asc, bs, audioMode, cb) != TRANSPORTDEC_OK) {
    return TRANSPORTDEC_PARSE_ERROR;
  }

  return TRANSPORTDEC_OK;
}

const UCHAR mapSr2MPEGIdx[8] = {
    0x1b,
    0x09,
    0x08,
    0x17,
    0x06,
    0x05,
    0x12,
    0x03
};

TRANSPORTDEC_ERROR DrmRawSdcAudioConfig_Parse(
    CSAudioSpecificConfig *self, HANDLE_FDK_BITSTREAM bs,
    CSTpCallBacks *cb,
    UCHAR configMode, UCHAR configChanged) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;

  AudioSpecificConfig_Init(self);

  if ((INT)FDKgetValidBits(bs) < 16) {
    ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
    goto bail;
  } else {
    int audioCoding, audioMode, cSamplingFreq, coderField, sfIdx, sbrFlag;

    self->configMode = configMode;
    self->AacConfigChanged = configChanged;
    self->SbrConfigChanged = configChanged;
    self->SacConfigChanged = configChanged;

    audioCoding = FDKreadBits(bs, 2);
    sbrFlag = FDKreadBits(bs, 1);
    audioMode = FDKreadBits(bs, 2);
    cSamplingFreq = FDKreadBits(bs, 3);

    FDKreadBits(bs, 2);
    coderField = FDKreadBits(bs, 5);
    FDKreadBits(bs, 1);

    if (audioCoding == 3) {
      sfIdx = (int)mapSr2MPEGIdx[cSamplingFreq];
      sbrFlag = 0;
    } else {
      switch (cSamplingFreq) {
        case 0:
          sfIdx = 11;
          break;
        case 1:
          sfIdx = 9;
          break;
        case 2:
          sfIdx = 8;
          break;
        case 3:
          sfIdx = 6;
          break;
        case 5:
          sfIdx = 3;
          break;
        case 4:
        case 6:
        case 7:
        default:
          ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
          goto bail;
      }
    }

    self->m_samplingFrequencyIndex = sfIdx;
    self->m_samplingFrequency = SamplingRateTable[sfIdx];

    if (sbrFlag) {
      UINT i;
      int tmp = -1;
      self->m_sbrPresentFlag = 1;
      self->m_extensionAudioObjectType = AOT_SBR;
      self->m_extensionSamplingFrequency = self->m_samplingFrequency << 1;
      for (i = 0;
           i < (sizeof(SamplingRateTable) / sizeof(SamplingRateTable[0]));
           i++) {
        if (SamplingRateTable[i] == self->m_extensionSamplingFrequency) {
          tmp = i;
          break;
        }
      }
      self->m_extensionSamplingFrequencyIndex = tmp;
    }

    switch (audioCoding) {
      case 0:
        if ((coderField >> 2) && (audioMode != 1)) {
          self->m_aot = AOT_DRM_SURROUND;
        } else {
          self->m_aot = AOT_DRM_AAC;
        }
        switch (audioMode) {
          case 1:
            self->m_psPresentFlag = 1;
            FDK_FALLTHROUGH;
          case 0:
            self->m_channelConfiguration = 1;
            break;
          case 2:
            self->m_channelConfiguration = 2;
            break;
          default:
            ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
            goto bail;
        }
        self->m_vcb11Flag = 1;
        self->m_hcrFlag = 1;
        self->m_samplesPerFrame = 960;
        self->m_epConfig = 1;
        break;
      case 1:
        self->m_aot = AOT_ER_CELP;
        self->m_channelConfiguration = 1;
        break;
      case 2:
        self->m_aot = AOT_ER_HVXC;
        self->m_channelConfiguration = 1;
        break;
      case 3:
      {
        self->m_aot = AOT_USAC;
      }
        switch (audioMode) {
          case 0:
          case 2:
            ErrorStatus = Drm_xHEAACStaticConfig(self, bs, audioMode, cb);
            break;
          default:
            ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
            goto bail;
        }
        break;
      default:
        ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
        self->m_aot = AOT_NONE;
        break;
    }

    if (self->m_psPresentFlag && !self->m_sbrPresentFlag) {
      ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
      goto bail;
    }
  }

bail:
  return (ErrorStatus);
}
