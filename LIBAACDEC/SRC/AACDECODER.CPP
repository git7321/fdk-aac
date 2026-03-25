/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description:

*******************************************************************************/

#include "aacdecoder.h"
#include "aac_rom.h"
#include "aac_ram.h"
#include "channel.h"
#include "FDK_audio.h"

#include "aacdec_pns.h"
#include "sbrdecoder.h"
#include "sac_dec_lib.h"
#include "aacdec_hcr.h"
#include "rvlc.h"
#include "usacdec_lpd.h"
#include "ac_arith_coder.h"
#include "tpdec_lib.h"
#include "conceal.h"
#include "FDK_crc.h"

#define PS_IS_EXPLICITLY_DISABLED(aot, flags) \
  (((aot) == AOT_DRM_AAC) && !(flags & AC_PS_PRESENT))

#define IS_STEREO_SBR(el_id, stereoConfigIndex)            \
  (((el_id) == ID_USAC_CPE && (stereoConfigIndex) == 0) || \
   ((el_id) == ID_USAC_CPE && (stereoConfigIndex) == 3))

void CAacDecoder_SyncQmfMode(HANDLE_AACDECODER self) {
  FDK_ASSERT(
      !((self->flags[0] & AC_MPS_PRESENT) && (self->flags[0] & AC_PS_PRESENT)));

  self->qmfModeCurr = self->qmfModeUser;

  if (IS_USAC(self->streamInfo.aot)) {
    self->qmfModeCurr = MODE_HQ;
  }

  if (self->qmfModeCurr == NOT_DEFINED) {
    if ((IS_LOWDELAY(self->streamInfo.aot) &&
         (self->flags[0] & AC_MPS_PRESENT)) ||
        ((self->streamInfo.aacNumChannels == 1) &&
         ((CAN_DO_PS(self->streamInfo.aot) &&
           !(self->flags[0] & AC_MPS_PRESENT)) ||
          (IS_USAC(self->streamInfo.aot))))) {
      self->qmfModeCurr = MODE_HQ;
    } else {
      self->qmfModeCurr = MODE_LP;
    }
  }

  if (self->mpsEnableCurr) {
    if (IS_LOWDELAY(self->streamInfo.aot) &&
        (self->qmfModeCurr == MODE_LP)) {
      self->qmfModeCurr = MODE_HQ;
    }
    switch (mpegSurroundDecoder_SetParam(
        (CMpegSurroundDecoder *)self->pMpegSurroundDecoder,
        SACDEC_PARTIALLY_COMPLEX, self->qmfModeCurr == MODE_LP)) {
      case MPS_OK:
        break;
      case MPS_INVALID_PARAMETER: {
        LIB_INFO libInfo[FDK_MODULE_LAST];
        UINT mpsCaps;

        FDKinitLibInfo(libInfo);
        mpegSurroundDecoder_GetLibInfo(libInfo);
        mpsCaps = FDKlibInfo_getCapabilities(libInfo, FDK_MPSDEC);

        if (((mpsCaps & CAPF_MPS_LP) && (self->qmfModeCurr == MODE_LP)) ||
            ((mpsCaps & CAPF_MPS_HQ) &&
             (self->qmfModeCurr ==
              MODE_HQ))) {
          break;
        }
      }
        FDK_FALLTHROUGH;
      default:
        if (self->qmfModeUser == NOT_DEFINED) {
          self->qmfModeCurr =
              (self->qmfModeCurr == MODE_LP) ? MODE_HQ : MODE_LP;
        } else {
          self->mpsEnableCurr = 0;
        }
    }
  }

  sbrDecoder_SetParam(self->hSbrDecoder, SBR_QMF_MODE,
                      (self->qmfModeCurr == MODE_LP));
  self->psPossible =
      ((CAN_DO_PS(self->streamInfo.aot) &&
        !PS_IS_EXPLICITLY_DISABLED(self->streamInfo.aot, self->flags[0]) &&
        self->streamInfo.aacNumChannels == 1 &&
        !(self->flags[0] & AC_MPS_PRESENT))) &&
      self->qmfModeCurr == MODE_HQ;
  FDK_ASSERT(!((self->flags[0] & AC_MPS_PRESENT) && self->psPossible));
}

void CAacDecoder_SignalInterruption(HANDLE_AACDECODER self) {
  if (self->flags[0] & (AC_USAC | AC_RSVD50 | AC_RSV603DA)) {
    int i;

    for (i = 0; i < fMin(self->aacChannels, (8)); i++) {
      if (self->pAacDecoderStaticChannelInfo
              [i]) {
        self->pAacDecoderStaticChannelInfo[i]->hArCo->m_numberLinesPrev = 0;
      }
    }
  }
}

static int CAacDecoder_GetELChannels(MP4_ELEMENT_ID type,
                                     UCHAR usacStereoConfigIndex) {
  int el_channels = 0;

  switch (type) {
    case ID_USAC_CPE:
      if (usacStereoConfigIndex == 1) {
        el_channels = 1;
      } else {
        el_channels = 2;
      }
      break;
    case ID_CPE:
      el_channels = 2;
      break;
    case ID_USAC_SCE:
    case ID_USAC_LFE:
    case ID_SCE:
    case ID_LFE:
      el_channels = 1;
      break;
    default:
      el_channels = 0;
      break;
  }

  return el_channels;
}

static AAC_DECODER_ERROR CAacDecoder_AncDataReset(CAncData *ancData) {
  int i;
  for (i = 0; i < 8; i++) {
    ancData->offset[i] = 0;
  }
  ancData->nrElements = 0;

  return AAC_DEC_OK;
}

AAC_DECODER_ERROR CAacDecoder_AncDataInit(CAncData *ancData,
                                          unsigned char *buffer, int size) {
  if (size >= 0) {
    ancData->buffer = buffer;
    ancData->bufferSize = size;

    CAacDecoder_AncDataReset(ancData);

    return AAC_DEC_OK;
  }

  return AAC_DEC_ANC_DATA_ERROR;
}

AAC_DECODER_ERROR CAacDecoder_AncDataGet(CAncData *ancData, int index,
                                         unsigned char **ptr, int *size) {
  AAC_DECODER_ERROR error = AAC_DEC_OK;

  *ptr = NULL;
  *size = 0;

  if (index >= 0 && index < 8 - 1 && index < ancData->nrElements) {
    *ptr = &ancData->buffer[ancData->offset[index]];
    *size = ancData->offset[index + 1] - ancData->offset[index];
  }

  return error;
}

static AAC_DECODER_ERROR CAacDecoder_AncDataParse(CAncData *ancData,
                                                  HANDLE_FDK_BITSTREAM hBs,
                                                  const int ancBytes) {
  AAC_DECODER_ERROR error = AAC_DEC_OK;
  int readBytes = 0;

  if (ancData->buffer != NULL) {
    if (ancBytes > 0) {
      int offset = ancData->offset[ancData->nrElements];

      if ((offset + ancBytes) > ancData->bufferSize) {
        error = AAC_DEC_TOO_SMALL_ANC_BUFFER;
      } else if (ancData->nrElements >= 8 - 1) {
        error = AAC_DEC_TOO_MANY_ANC_ELEMENTS;
      } else {
        int i;

        for (i = 0; i < ancBytes; i++) {
          ancData->buffer[i + offset] = FDKreadBits(hBs, 8);
          readBytes++;
        }

        ancData->nrElements++;
        ancData->offset[ancData->nrElements] =
            ancBytes + ancData->offset[ancData->nrElements - 1];
      }
    }
  }

  readBytes = ancBytes - readBytes;

  if (readBytes > 0) {
    FDKpushFor(hBs, readBytes << 3);
  }

  return error;
}

static AAC_DECODER_ERROR CDataStreamElement_Read(HANDLE_AACDECODER self,
                                                 HANDLE_FDK_BITSTREAM bs,
                                                 UCHAR *elementInstanceTag,
                                                 UINT alignmentAnchor) {
  AAC_DECODER_ERROR error = AAC_DEC_OK;
  UINT dseBits;
  INT dataStart;
  int dataByteAlignFlag, count;

  FDK_ASSERT(self != NULL);

  int crcReg = transportDec_CrcStartReg(self->hInput, 0);

  *elementInstanceTag = FDKreadBits(bs, 4);
  dataByteAlignFlag = FDKreadBits(bs, 1);

  count = FDKreadBits(bs, 8);

  if (count == 255) {
    count += FDKreadBits(bs, 8);
  }
  dseBits = count * 8;

  if (dataByteAlignFlag) {
    FDKbyteAlign(bs, alignmentAnchor);
  }

  dataStart = (INT)FDKgetValidBits(bs);

  error = CAacDecoder_AncDataParse(&self->ancData, bs, count);
  transportDec_CrcEndReg(self->hInput, crcReg);

  {
    FDKpushBack(bs, dataStart - (INT)FDKgetValidBits(bs));

    aacDecoder_drcMarkPayload(self->hDrcInfo, bs, DVB_DRC_ANC_DATA);
  }

  {
    PCMDMX_ERROR dmxErr = PCMDMX_OK;
    FDKpushBack(bs, dataStart - (INT)FDKgetValidBits(bs));

    dmxErr = pcmDmx_Parse(self->hPcmUtils, bs, dseBits, 0);
    if (error == AAC_DEC_OK && dmxErr != PCMDMX_OK) {
      error = AAC_DEC_UNKNOWN;
    }
  }

  FDKpushBiDirectional(bs, (INT)FDKgetValidBits(bs) - dataStart + (INT)dseBits);

  return error;
}

static INT findElementInstanceTag(
    INT elementTag, MP4_ELEMENT_ID elementId,
    CAacDecoderChannelInfo **pAacDecoderChannelInfo, INT nChannels,
    MP4_ELEMENT_ID *pElementIdTab, INT nElements) {
  int el, chCnt = 0;

  for (el = 0; el < nElements; el++) {
    switch (pElementIdTab[el]) {
      case ID_CPE:
      case ID_SCE:
      case ID_LFE:
        if ((elementTag == pAacDecoderChannelInfo[chCnt]->ElementInstanceTag) &&
            (elementId == pElementIdTab[el])) {
          return 1;
        }
        chCnt += (pElementIdTab[el] == ID_CPE) ? 2 : 1;
        break;
      default:
        break;
    }
    if (chCnt >= nChannels) break;
    if (pElementIdTab[el] == ID_END) break;
  }

  return 0;
}

static INT validateElementInstanceTags(
    CProgramConfig *pce, CAacDecoderChannelInfo **pAacDecoderChannelInfo,
    INT nChannels, MP4_ELEMENT_ID *pElementIdTab, INT nElements) {
  if (nChannels >= pce->NumChannels) {
    for (int el = 0; el < pce->NumFrontChannelElements; el++) {
      if (!findElementInstanceTag(pce->FrontElementTagSelect[el],
                                  pce->FrontElementIsCpe[el] ? ID_CPE : ID_SCE,
                                  pAacDecoderChannelInfo, nChannels,
                                  pElementIdTab, nElements)) {
        return 0;
      }
    }
    for (int el = 0; el < pce->NumSideChannelElements; el++) {
      if (!findElementInstanceTag(pce->SideElementTagSelect[el],
                                  pce->SideElementIsCpe[el] ? ID_CPE : ID_SCE,
                                  pAacDecoderChannelInfo, nChannels,
                                  pElementIdTab, nElements)) {
        return 0;
      }
    }
    for (int el = 0; el < pce->NumBackChannelElements; el++) {
      if (!findElementInstanceTag(pce->BackElementTagSelect[el],
                                  pce->BackElementIsCpe[el] ? ID_CPE : ID_SCE,
                                  pAacDecoderChannelInfo, nChannels,
                                  pElementIdTab, nElements)) {
        return 0;
      }
    }
    for (int el = 0; el < pce->NumLfeChannelElements; el++) {
      if (!findElementInstanceTag(pce->LfeElementTagSelect[el], ID_LFE,
                                  pAacDecoderChannelInfo, nChannels,
                                  pElementIdTab, nElements)) {
        return 0;
      }
    }
  } else {
    return 0;
  }

  return 1;
}

static int CProgramConfigElement_Read(HANDLE_FDK_BITSTREAM bs,
                                      HANDLE_TRANSPORTDEC pTp,
                                      CProgramConfig *pce,
                                      const UINT channelConfig,
                                      const UINT alignAnchor) {
  int pceStatus = 0;
  int crcReg;

  C_ALLOC_SCRATCH_START(tmpPce, CProgramConfig, 1);

  CProgramConfig_Init(tmpPce);

  crcReg = transportDec_CrcStartReg(pTp, 0);

  CProgramConfig_Read(tmpPce, bs, alignAnchor);

  transportDec_CrcEndReg(pTp, crcReg);

  if (CProgramConfig_IsValid(tmpPce) && (tmpPce->Profile == 1)) {
    if (!CProgramConfig_IsValid(pce) && (channelConfig > 0)) {
      CProgramConfig_GetDefault(pce, channelConfig);
    }

    if (CProgramConfig_IsValid(pce)) {
      switch (CProgramConfig_Compare(pce, tmpPce)) {
        case 1:
          FDKmemcpy(pce, tmpPce,
                    sizeof(CProgramConfig));
          pceStatus = 1;
          break;
        case 2:
        case -1:
          pceStatus = -1;
          break;
        case 0:
        default:
          break;
      }
    }
  }

  C_ALLOC_SCRATCH_END(tmpPce, CProgramConfig, 1);

  return pceStatus;
}

AAC_DECODER_ERROR CAacDecoder_PrepareCrossFade(
    const PCM_DEC *pTimeData, PCM_DEC **pTimeDataFlush, const INT numChannels,
    const INT frameSize, const INT interleaved) {
  int i, ch, s1, s2;
  AAC_DECODER_ERROR ErrorStatus;

  ErrorStatus = AAC_DEC_OK;

  if (interleaved) {
    s1 = 1;
    s2 = numChannels;
  } else {
    s1 = frameSize;
    s2 = 1;
  }

  for (ch = 0; ch < numChannels; ch++) {
    const PCM_DEC *pIn = &pTimeData[ch * s1];
    for (i = 0; i < TIME_DATA_FLUSH_SIZE; i++) {
      pTimeDataFlush[ch][i] = *pIn;
      pIn += s2;
    }
  }

  return ErrorStatus;
}

AAC_DECODER_ERROR CAacDecoder_ApplyCrossFade(
    PCM_DEC *pTimeData, PCM_DEC **pTimeDataFlush, const INT numChannels,
    const INT frameSize, const INT interleaved) {
  int i, ch, s1, s2;
  AAC_DECODER_ERROR ErrorStatus;

  ErrorStatus = AAC_DEC_OK;

  if (interleaved) {
    s1 = 1;
    s2 = numChannels;
  } else {
    s1 = frameSize;
    s2 = 1;
  }

  for (ch = 0; ch < numChannels; ch++) {
    PCM_DEC *pIn = &pTimeData[ch * s1];
    for (i = 0; i < TIME_DATA_FLUSH_SIZE; i++) {
      FIXP_SGL alpha = (FIXP_SGL)i
                       << (FRACT_BITS - 1 - TIME_DATA_FLUSH_SIZE_SF);
      FIXP_DBL time = PCM_DEC2FIXP_DBL(*pIn);
      FIXP_DBL timeFlush = PCM_DEC2FIXP_DBL(pTimeDataFlush[ch][i]);

      *pIn = FIXP_DBL2PCM_DEC(timeFlush - fMult(timeFlush, alpha) +
                              fMult(time, alpha));
      pIn += s2;
    }
  }

  return ErrorStatus;
}

AAC_DECODER_ERROR CAacDecoder_PreRollExtensionPayloadParse(
    HANDLE_AACDECODER self, UINT *numPrerollAU, UINT *prerollAUOffset,
    UINT *prerollAULength) {
  FDK_BITSTREAM bs;
  HANDLE_FDK_BITSTREAM hBs;
  AAC_DECODER_ERROR ErrorStatus;

  INT auStartAnchor;
  UINT independencyFlag;
  UINT extPayloadPresentFlag;
  UINT useDefaultLengthFlag;
  UINT configLength = 0;
  UINT preRollPossible = 1;
  UINT i;
  UCHAR configChanged = 0;
  UCHAR config[TP_USAC_MAX_CONFIG_LEN] = {0};
  UCHAR
  implicitExplicitCfgDiff = 0;

  ErrorStatus = AAC_DEC_OK;

  hBs = transportDec_GetBitstream(self->hInput, 0);
  bs = *hBs;

  auStartAnchor = (INT)FDKgetValidBits(hBs);
  if (auStartAnchor <= 0) {
    ErrorStatus = AAC_DEC_NOT_ENOUGH_BITS;
    goto bail;
  }

  FDKreadBit(hBs);

  extPayloadPresentFlag = FDKreadBits(hBs, 1);
  if (!extPayloadPresentFlag) {
    preRollPossible = 0;
  }

  useDefaultLengthFlag = FDKreadBits(hBs, 1);
  if (useDefaultLengthFlag) {
    preRollPossible = 0;
  }

  if (preRollPossible) {
    escapedValue(hBs, 8, 16, 0);

    configLength = escapedValue(hBs, 4, 4, 8);
  }

  if (preRollPossible == 0) {
    if (self->flushStatus != AACDEC_FLUSH_OFF) {
      self->flushStatus = AACDEC_FLUSH_OFF;
      ErrorStatus = AAC_DEC_PARSE_ERROR;
    }
    goto bail;
  }

  if (self->flags[0] & AC_USAC) {
    if (configLength > 0) {
      for (i = 0; i < configLength; i++) {
        config[i] = FDKreadBits(hBs, 8);
      }
      TRANSPORTDEC_ERROR terr;
      terr = transportDec_InBandConfig(self->hInput, config, configLength,
                                       self->buildUpStatus, &configChanged, 0,
                                       &implicitExplicitCfgDiff);
      if (terr != TRANSPORTDEC_OK) {
        ErrorStatus = AAC_DEC_PARSE_ERROR;
        goto bail;
      }
    }
  }

  if ((self->streamInfo.numChannels == 0) && !implicitExplicitCfgDiff &&
      (self->flags[0] & AC_USAC)) {
    self->buildUpStatus = AACDEC_USAC_BUILD_UP_ON;
    if (self->flushStatus != AACDEC_FLUSH_OFF) {
      self->flushStatus = AACDEC_FLUSH_OFF;
      ErrorStatus = AAC_DEC_PARSE_ERROR;
      goto bail;
    }
  }

  if (self->flags[0] & AC_USAC) {
    if (self->buildUpStatus == AACDEC_USAC_BUILD_UP_ON) {
      UCHAR applyCrossfade = FDKreadBit(hBs);
      if (applyCrossfade) {
        self->applyCrossfade |= AACDEC_CROSSFADE_BITMASK_PREROLL;
      } else {
        self->applyCrossfade &= ~AACDEC_CROSSFADE_BITMASK_PREROLL;
      }
      FDKreadBit(hBs);
      *numPrerollAU = escapedValue(hBs, 2, 4, 0);
      if (*numPrerollAU > AACDEC_MAX_NUM_PREROLL_AU_USAC) {
        *numPrerollAU = 0;
        ErrorStatus = AAC_DEC_PARSE_ERROR;
        goto bail;
      }
    }
  }

  for (i = 0; i < *numPrerollAU; i++) {
    prerollAULength[i] = escapedValue(hBs, 16, 16, 0);
    if (prerollAULength[i] > 0) {
      prerollAUOffset[i] = auStartAnchor - (INT)FDKgetValidBits(hBs);
      independencyFlag = FDKreadBit(hBs);
      if (i == 0 && !independencyFlag) {
        *numPrerollAU = 0;
        ErrorStatus = AAC_DEC_PARSE_ERROR;
        goto bail;
      }
      FDKpushFor(hBs, prerollAULength[i] * 8 - 1);
      self->prerollAULength[i] = (prerollAULength[i] * 8) + prerollAUOffset[i];
    } else {
      *numPrerollAU = 0;
      ErrorStatus = AAC_DEC_PARSE_ERROR;
      goto bail;
    }
  }

bail:

  *hBs = bs;

  return ErrorStatus;
}

static AAC_DECODER_ERROR CAacDecoder_ExtPayloadParse(
    HANDLE_AACDECODER self, HANDLE_FDK_BITSTREAM hBs, int *count,
    MP4_ELEMENT_ID previous_element, int elIndex, int fIsFillElement) {
  AAC_DECODER_ERROR error = AAC_DEC_OK;
  EXT_PAYLOAD_TYPE extension_type;
  int bytes = (*count) >> 3;
  int crcFlag = 0;

  if (*count < 4) {
    return AAC_DEC_PARSE_ERROR;
  } else if ((INT)FDKgetValidBits(hBs) < *count) {
    return AAC_DEC_DECODE_FRAME_ERROR;
  }

  extension_type =
      (EXT_PAYLOAD_TYPE)FDKreadBits(hBs, 4);
  *count -= 4;

  if ((self->flags[0] & AC_ELD) && ((extension_type == EXT_SBR_DATA_CRC) ||
                                    (extension_type == EXT_SBR_DATA))) {
    extension_type = EXT_FIL;
  }

  switch (extension_type) {
    case EXT_DYNAMIC_RANGE: {
      INT readBits =
          aacDecoder_drcMarkPayload(self->hDrcInfo, hBs, MPEG_DRC_EXT_DATA);

      if (readBits > *count) {
        error = AAC_DEC_PARSE_ERROR;
      }
      *count -= readBits;
    } break;
    case EXT_LDSAC_DATA:
    case EXT_SAC_DATA:
      {
        int err, mpsSampleRate, mpsFrameSize;

        if (self->flags[0] & AC_PS_PRESENT) {
          error = AAC_DEC_PARSE_ERROR;
          goto bail;
        }

        if (self->streamInfo.extSamplingRate != 0) {
          mpsSampleRate = self->streamInfo.extSamplingRate;
          mpsFrameSize = self->streamInfo.aacSamplesPerFrame *
                         (self->streamInfo.extSamplingRate /
                          self->streamInfo.aacSampleRate);
        } else {
          mpsSampleRate = self->streamInfo.aacSampleRate;
          mpsFrameSize = self->streamInfo.aacSamplesPerFrame;
        }
        if (!(self->flags[0] & (AC_USAC | AC_RSVD50 | AC_ELD))) {
          self->mpsEnableCurr = self->mpsEnableUser;
        }
        if (self->mpsEnableCurr) {
          if (!self->qmfDomain.globalConf.qmfDomainExplicitConfig) {
            if (mpegSurroundDecoder_IsFullMpegSurroundDecoderInstanceAvailable(
                    (CMpegSurroundDecoder *)self->pMpegSurroundDecoder) ==
                SAC_INSTANCE_NOT_FULL_AVAILABLE) {
              if (mpegSurroundDecoder_Open(
                      (CMpegSurroundDecoder **)&self->pMpegSurroundDecoder, -1,
                      &self->qmfDomain)) {
                return AAC_DEC_OUT_OF_MEMORY;
              }
            }
          }
          err = mpegSurroundDecoder_Parse(
              (CMpegSurroundDecoder *)self->pMpegSurroundDecoder, hBs, count,
              self->streamInfo.aot, mpsSampleRate, mpsFrameSize,
              self->flags[0] & AC_INDEP);
          if (err == MPS_OK) {
            self->flags[0] |= AC_MPS_PRESENT;
          } else {
            error = AAC_DEC_PARSE_ERROR;
          }
        }
        FDKpushFor(hBs, *count);
        *count = 0;
      }
      break;

    case EXT_SBR_DATA_CRC:
      crcFlag = 1;
      FDK_FALLTHROUGH;
    case EXT_SBR_DATA:
      if (IS_CHANNEL_ELEMENT(previous_element)) {
        SBR_ERROR sbrError;
        UCHAR configMode = 0;
        UCHAR configChanged = 0;

        CAacDecoder_SyncQmfMode(self);

        configMode |= AC_CM_ALLOC_MEM;

        sbrError = sbrDecoder_InitElement(
            self->hSbrDecoder, self->streamInfo.aacSampleRate,
            self->streamInfo.extSamplingRate,
            self->streamInfo.aacSamplesPerFrame, self->streamInfo.aot,
            previous_element, elIndex,
            2,
            0, configMode, &configChanged, self->downscaleFactor);

        if (sbrError == SBRDEC_OK) {
          sbrError = sbrDecoder_Parse(self->hSbrDecoder, hBs,
                                      self->pDrmBsBuffer, self->drmBsBufferSize,
                                      count, *count, crcFlag, previous_element,
                                      elIndex, self->flags[0], self->elFlags);
          if ((sbrError == SBRDEC_OK) || (sbrError == SBRDEC_PARSE_ERROR)) {
            self->sbrEnabled = 1;
          }
        } else {
          self->sbrEnabled = 0;
        }
        if (fIsFillElement) {
          FDKpushBiDirectional(hBs, *count);
          *count = 0;
        } else {
          if (sbrError != SBRDEC_OK) {
            self->frameOK = 0;
          }
        }
      } else {
        error = AAC_DEC_PARSE_ERROR;
      }
      break;

    case EXT_FILL_DATA: {
      int temp;

      temp = FDKreadBits(hBs, 4);
      bytes--;
      if (temp != 0) {
        error = AAC_DEC_PARSE_ERROR;
        break;
      }
      while (bytes > 0) {
        temp = FDKreadBits(hBs, 8);
        bytes--;
        if (temp != 0xa5) {
          error = AAC_DEC_PARSE_ERROR;
          break;
        }
      }
      *count = bytes << 3;
    } break;

    case EXT_DATA_ELEMENT: {
      int dataElementVersion;

      dataElementVersion = FDKreadBits(hBs, 4);
      *count -= 4;
      if (dataElementVersion == 0)
      {
        int temp, dataElementLength = 0;
        do {
          temp = FDKreadBits(hBs, 8);
          *count -= 8;
          dataElementLength += temp;
        } while (temp == 255);

        CAacDecoder_AncDataParse(&self->ancData, hBs, dataElementLength);
        *count -= (dataElementLength << 3);
      } else {
        error = AAC_DEC_PARSE_ERROR;
        goto bail;
      }
    } break;

    case EXT_DATA_LENGTH:
      if (!fIsFillElement
          &&
          (self->flags[0] &
           AC_ER))
      {
        int bitCnt, len = FDKreadBits(hBs, 4);
        *count -= 4;

        if (len == 15) {
          int add_len = FDKreadBits(hBs, 8);
          *count -= 8;
          len += add_len;

          if (add_len == 255) {
            len += FDKreadBits(hBs, 16);
            *count -= 16;
          }
        }
        len <<= 3;
        bitCnt = len;

        if ((EXT_PAYLOAD_TYPE)FDKreadBits(hBs, 4) == EXT_DATA_LENGTH) {
          error = AAC_DEC_PARSE_ERROR;
          goto bail;
        } else {
          FDKpushBack(hBs, 4);

          error = CAacDecoder_ExtPayloadParse(
              self, hBs, &bitCnt, previous_element, elIndex,
              1);

          *count -= len - bitCnt;
        }
        break;
      }
      FDK_FALLTHROUGH;

    case EXT_FIL:

    default:
      FDKpushFor(hBs, *count);
      *count = 0;
      break;
  }

bail:
  if ((error != AAC_DEC_OK) &&
      fIsFillElement) {
    FDKpushBiDirectional(hBs, *count);
    *count = 0;
    error = AAC_DEC_OK;
  }
  return error;
}

static AAC_DECODER_ERROR aacDecoder_ParseExplicitMpsAndSbr(
    HANDLE_AACDECODER self, HANDLE_FDK_BITSTREAM bs,
    const MP4_ELEMENT_ID previous_element, const int previous_element_index,
    const int element_index, const int el_cnt[]) {
  AAC_DECODER_ERROR ErrorStatus = AAC_DEC_OK;
  INT bitCnt = 0;

  bitCnt = transportDec_GetAuBitsRemaining(self->hInput, 0);

  if ((self->flags[0] & AC_SBR_PRESENT) &&
      (self->flags[0] & (AC_USAC | AC_RSVD50 | AC_ELD | AC_DRM))) {
    SBR_ERROR err = SBRDEC_OK;
    int chElIdx, numChElements = el_cnt[ID_SCE] + el_cnt[ID_CPE] +
                                 el_cnt[ID_LFE] + el_cnt[ID_USAC_SCE] +
                                 el_cnt[ID_USAC_CPE] + el_cnt[ID_USAC_LFE];
    INT bitCntTmp = bitCnt;

    if (self->flags[0] & AC_USAC) {
      chElIdx = numChElements - 1;
    } else {
      chElIdx = 0;
    }

    for (; chElIdx < numChElements; chElIdx += 1) {
      MP4_ELEMENT_ID sbrType;
      SBR_ERROR errTmp;
      if (self->flags[0] & (AC_USAC)) {
        FDK_ASSERT((self->elements[element_index] == ID_USAC_SCE) ||
                   (self->elements[element_index] == ID_USAC_CPE));
        sbrType = IS_STEREO_SBR(self->elements[element_index],
                                self->usacStereoConfigIndex[element_index])
                      ? ID_CPE
                      : ID_SCE;
      } else
        sbrType = self->elements[chElIdx];
      errTmp = sbrDecoder_Parse(self->hSbrDecoder, bs, self->pDrmBsBuffer,
                                self->drmBsBufferSize, &bitCnt, -1,
                                self->flags[0] & AC_SBRCRC, sbrType, chElIdx,
                                self->flags[0], self->elFlags);
      if (errTmp != SBRDEC_OK) {
        err = errTmp;
        bitCntTmp = bitCnt;
        bitCnt = 0;
      }
    }
    switch (err) {
      case SBRDEC_PARSE_ERROR:
        FDKpushFor(bs, bitCntTmp);
        bitCnt = 0;
        break;
      case SBRDEC_OK:
        self->sbrEnabled = 1;
        break;
      default:
        self->frameOK = 0;
        break;
    }
  }

  if ((bitCnt > 0) && (self->flags[0] & (AC_USAC | AC_RSVD50))) {
    if ((self->flags[0] & AC_MPS_PRESENT) ||
        (self->elFlags[element_index] & AC_EL_USAC_MPS212)) {
      int err;

      err = mpegSurroundDecoder_ParseNoHeader(
          (CMpegSurroundDecoder *)self->pMpegSurroundDecoder, bs, &bitCnt,
          self->flags[0] & AC_INDEP);
      if (err != MPS_OK) {
        self->frameOK = 0;
        ErrorStatus = AAC_DEC_PARSE_ERROR;
      }
    }
  }

  if (self->flags[0] & AC_DRM) {
    if ((bitCnt = (INT)FDKgetValidBits(bs)) != 0) {
      FDKpushBiDirectional(bs, bitCnt);
    }
  }

  if (!(self->flags[0] & (AC_USAC | AC_RSVD50 | AC_DRM))) {
    while (bitCnt > 7) {
      ErrorStatus = CAacDecoder_ExtPayloadParse(
          self, bs, &bitCnt, previous_element, previous_element_index, 0);
      if (ErrorStatus != AAC_DEC_OK) {
        self->frameOK = 0;
        ErrorStatus = AAC_DEC_PARSE_ERROR;
        break;
      }
    }
  }
  return ErrorStatus;
}

static void CStreamInfoInit(CStreamInfo *pStreamInfo) {
  pStreamInfo->aacSampleRate = 0;
  pStreamInfo->profile = -1;
  pStreamInfo->aot = AOT_NONE;

  pStreamInfo->channelConfig = -1;
  pStreamInfo->bitRate = 0;
  pStreamInfo->aacSamplesPerFrame = 0;

  pStreamInfo->extAot = AOT_NONE;
  pStreamInfo->extSamplingRate = 0;

  pStreamInfo->flags = 0;

  pStreamInfo->epConfig = -1;

  pStreamInfo->numChannels = 0;
  pStreamInfo->sampleRate = 0;
  pStreamInfo->frameSize = 0;

  pStreamInfo->outputDelay = 0;

  pStreamInfo->drcProgRefLev =
      -1;
  pStreamInfo->drcPresMode = -1;

  pStreamInfo->outputLoudness = -1;
}

HANDLE_AACDECODER CAacDecoder_Open(
    TRANSPORT_TYPE bsFormat)
{
  HANDLE_AACDECODER self;

  self = GetAacDecoder();
  if (self == NULL) {
    goto bail;
  }

  FDK_QmfDomain_ClearRequested(&self->qmfDomain.globalConf);

  self->streamInfo.pChannelIndices = self->channelIndices;
  self->streamInfo.pChannelType = self->channelType;
  self->downscaleFactor = 1;
  self->downscaleFactorInBS = 1;

  CAacDecoder_AncDataInit(&self->ancData, NULL, 0);
  CStreamInfoInit(&self->streamInfo);
  CProgramConfig_Init(&self->pce);

  CConcealment_InitCommonData(&self->concealCommonData);
  self->concealMethodUser = ConcealMethodNone;

  self->hDrcInfo = GetDrcInfo();
  if (self->hDrcInfo == NULL) {
    goto bail;
  }
  aacDecoder_drcInit(self->hDrcInfo);
  aacDecoder_drcSetParam(self->hDrcInfo, DRC_BS_DELAY,
                         CConcealment_GetDelay(&self->concealCommonData));
  self->workBufferCore1 = (FIXP_DBL *)GetWorkBufferCore1();

  self->workBufferCore2 = GetWorkBufferCore2();
  if (self->workBufferCore2 == NULL) goto bail;

  self->pTimeData2 = GetWorkBufferCore5();
  self->timeData2Size = GetRequiredMemWorkBufferCore5();
  if (self->pTimeData2 == NULL) {
    goto bail;
  }

  return self;

bail:
  CAacDecoder_Close(self);

  return NULL;
}

static void CAacDecoder_DeInit(HANDLE_AACDECODER self,
                               const int subStreamIndex) {
  int ch;
  int aacChannelOffset = 0, aacChannels = (8);
  int numElements = (3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1),
      elementOffset = 0;

  if (self == NULL) return;

  {
    self->ascChannels[0] = 0;
    self->elements[0] = ID_END;
  }

  for (ch = aacChannelOffset; ch < aacChannelOffset + aacChannels; ch++) {
    if (self->pAacDecoderChannelInfo[ch] != NULL) {
      if (self->pAacDecoderChannelInfo[ch]->pComStaticData != NULL) {
        if (self->pAacDecoderChannelInfo[ch]
                ->pComStaticData->pWorkBufferCore1 != NULL) {
          if (ch == aacChannelOffset) {
            FreeWorkBufferCore1(&self->pAacDecoderChannelInfo[ch]
                                     ->pComStaticData->pWorkBufferCore1);
          }
        }
        if (self->pAacDecoderChannelInfo[ch]
                ->pComStaticData->cplxPredictionData != NULL) {
          FreeCplxPredictionData(&self->pAacDecoderChannelInfo[ch]
                                      ->pComStaticData->cplxPredictionData);
        }
        if (ch < (8) - 1) {
          if ((self->pAacDecoderChannelInfo[ch + 1] != NULL) &&
              (self->pAacDecoderChannelInfo[ch + 1]->pComStaticData ==
               self->pAacDecoderChannelInfo[ch]->pComStaticData)) {
            self->pAacDecoderChannelInfo[ch + 1]->pComStaticData = NULL;
          }
        }
        FDKfree(self->pAacDecoderChannelInfo[ch]->pComStaticData);
        self->pAacDecoderChannelInfo[ch]->pComStaticData = NULL;
      }
      if (self->pAacDecoderChannelInfo[ch]->pComData != NULL) {
        if (ch < (8) - 1) {
          if ((self->pAacDecoderChannelInfo[ch + 1] != NULL) &&
              (self->pAacDecoderChannelInfo[ch + 1]->pComData ==
               self->pAacDecoderChannelInfo[ch]->pComData)) {
            self->pAacDecoderChannelInfo[ch + 1]->pComData = NULL;
          }
        }
        if (ch == aacChannelOffset) {
          FreeWorkBufferCore6(
              (SCHAR **)&self->pAacDecoderChannelInfo[ch]->pComData);
        } else {
          FDKafree(self->pAacDecoderChannelInfo[ch]->pComData);
        }
        self->pAacDecoderChannelInfo[ch]->pComData = NULL;
      }
    }
    if (self->pAacDecoderStaticChannelInfo[ch] != NULL) {
      if (self->pAacDecoderStaticChannelInfo[ch]->pOverlapBuffer != NULL) {
        FreeOverlapBuffer(
            &self->pAacDecoderStaticChannelInfo[ch]->pOverlapBuffer);
      }
      if (self->pAacDecoderStaticChannelInfo[ch]->hArCo != NULL) {
        CArco_Destroy(self->pAacDecoderStaticChannelInfo[ch]->hArCo);
      }
      FreeAacDecoderStaticChannelInfo(&self->pAacDecoderStaticChannelInfo[ch]);
    }
    if (self->pAacDecoderChannelInfo[ch] != NULL) {
      FreeAacDecoderChannelInfo(&self->pAacDecoderChannelInfo[ch]);
    }
  }

  {
    int el;
    for (el = elementOffset; el < elementOffset + numElements; el++) {
      if (self->cpeStaticData[el] != NULL) {
        FreeCpePersistentData(&self->cpeStaticData[el]);
      }
    }
  }

  FDK_Delay_Destroy(&self->usacResidualDelay);

  self->aacChannels = 0;
  self->streamInfo.aacSampleRate = 0;
  self->streamInfo.sampleRate = 0;
  self->samplingRateInfo[subStreamIndex].samplingRate = 0;
}

static void CAacDecoder_AcceptFlags(HANDLE_AACDECODER self,
                                    const CSAudioSpecificConfig *asc,
                                    UINT flags, UINT *elFlags, int streamIndex,
                                    int elementOffset) {
  FDKmemcpy(self->elFlags, elFlags, sizeof(self->elFlags));

  self->flags[streamIndex] = flags;
}

AAC_DECODER_ERROR CAacDecoder_CtrlCFGChange(HANDLE_AACDECODER self,
                                                         UCHAR flushStatus,
                                                         SCHAR flushCnt,
                                                         UCHAR buildUpStatus,
                                                         SCHAR buildUpCnt) {
  AAC_DECODER_ERROR err = AAC_DEC_OK;

  self->flushStatus = flushStatus;
  self->flushCnt = flushCnt;
  self->buildUpStatus = buildUpStatus;
  self->buildUpCnt = buildUpCnt;

  return (err);
}

AAC_DECODER_ERROR CAacDecoder_FreeMem(HANDLE_AACDECODER self,
                                                   const int subStreamIndex) {
  AAC_DECODER_ERROR err = AAC_DEC_OK;

  CAacDecoder_DeInit(self, subStreamIndex);

  return (err);
}

void CAacDecoder_Close(HANDLE_AACDECODER self) {
  if (self == NULL) return;

  CAacDecoder_DeInit(self, 0);

  {
    int ch;
    for (ch = 0; ch < (8); ch++) {
      if (self->pTimeDataFlush[ch] != NULL) {
        FreeTimeDataFlush(&self->pTimeDataFlush[ch]);
      }
    }
  }

  if (self->hDrcInfo) {
    FreeDrcInfo(&self->hDrcInfo);
  }

  if (self->workBufferCore1 != NULL) {
    FreeWorkBufferCore1((CWorkBufferCore1 **)&self->workBufferCore1);
  }

  if (self->workBufferCore2 != NULL) {
    FreeWorkBufferCore2(&self->workBufferCore2);
  }
  if (self->pTimeData2 != NULL) {
    FreeWorkBufferCore5(&self->pTimeData2);
  }

  FDK_QmfDomain_Close(&self->qmfDomain);

  FreeAacDecoder(&self);
}

AAC_DECODER_ERROR CAacDecoder_Init(HANDLE_AACDECODER self, const CSAudioSpecificConfig *asc,
                 UCHAR configMode, UCHAR *configChanged) {
  AAC_DECODER_ERROR err = AAC_DEC_OK;
  INT ascChannels, ascChanged = 0;
  AACDEC_RENDER_MODE initRenderMode = AACDEC_RENDER_INVALID;
  SCHAR usacStereoConfigIndex = -1;
  int usacResidualDelayCompSamples = 0;
  int elementOffset, aacChannelsOffset, aacChannelsOffsetIdx;
  const int streamIndex = 0;
  INT flushChannels = 0;

  UINT flags;
  UINT elFlags[(3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1)];

  UCHAR sbrEnabled = self->sbrEnabled;
  UCHAR sbrEnabledPrev = self->sbrEnabledPrev;
  UCHAR mpsEnableCurr = self->mpsEnableCurr;

  if (!self) return AAC_DEC_INVALID_HANDLE;

  UCHAR downscaleFactor = self->downscaleFactor;
  UCHAR downscaleFactorInBS = self->downscaleFactorInBS;

  self->aacOutDataHeadroom = (3);

  switch (asc->m_aot) {
    case AOT_AAC_LC:
      self->streamInfo.profile = 1;
      FDK_FALLTHROUGH;
    case AOT_ER_AAC_SCAL:
      if (asc->m_sc.m_gaSpecificConfig.m_layer > 0) {
        return AAC_DEC_UNSUPPORTED_FORMAT;
      }
      FDK_FALLTHROUGH;
    case AOT_SBR:
    case AOT_PS:
    case AOT_ER_AAC_LC:
    case AOT_ER_AAC_LD:
    case AOT_DRM_AAC:
    case AOT_DRM_SURROUND:
      initRenderMode = AACDEC_RENDER_IMDCT;
      break;
    case AOT_ER_AAC_ELD:
      initRenderMode = AACDEC_RENDER_ELDFB;
      break;
    case AOT_USAC:
      initRenderMode = AACDEC_RENDER_IMDCT;
      break;
    default:
      return AAC_DEC_UNSUPPORTED_AOT;
  }

  if (CProgramConfig_IsValid(&self->pce) && (asc->m_channelConfiguration > 0)) {
    int pceCmpResult;
    C_ALLOC_SCRATCH_START(tmpPce, CProgramConfig, 1);

    CProgramConfig_GetDefault(tmpPce, asc->m_channelConfiguration);
    pceCmpResult = CProgramConfig_Compare(&self->pce, tmpPce);
    if ((pceCmpResult < 0)
        ||
        (pceCmpResult > 1)) {
      CProgramConfig_Init(&self->pce);
    }
    C_ALLOC_SCRATCH_END(tmpPce, CProgramConfig, 1);
  } else {
    CProgramConfig_Init(&self->pce);
  }

  switch (asc->m_channelConfiguration) {
    case 0:
      switch (asc->m_aot) {
        case AOT_USAC:
          self->chMapIndex = 0;
          ascChannels = asc->m_sc.m_usacConfig.m_nUsacChannels;
          break;
        default:
          if (CProgramConfig_IsValid(&asc->m_progrConfigElement)) {
            ascChannels = asc->m_progrConfigElement.NumChannels;
            if (ascChannels > 0) {
              int el_tmp;
              FDKmemcpy(&self->pce, &asc->m_progrConfigElement,
                        sizeof(CProgramConfig));
              el_tmp = CProgramConfig_GetElementTable(
                  &asc->m_progrConfigElement, self->elements, (((8)) + (8)),
                  &self->chMapIndex);
              for (; el_tmp < (3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1);
                   el_tmp++) {
                self->elements[el_tmp] = ID_NONE;
              }
            } else {
              return AAC_DEC_UNSUPPORTED_CHANNELCONFIG;
            }
          } else {
            self->chMapIndex = 0;
            return AAC_DEC_UNSUPPORTED_CHANNELCONFIG;
          }
          break;
      }
      break;
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
      ascChannels = asc->m_channelConfiguration;
      break;
    case 11:
      ascChannels = 7;
      break;
    case 7:
    case 12:
    case 14:
      ascChannels = 8;
      break;
    default:
      return AAC_DEC_UNSUPPORTED_CHANNELCONFIG;
  }

  if (asc->m_aot == AOT_USAC) {
    flushChannels = fMin(ascChannels, (8));
    INT numChannel;
    pcmDmx_GetParam(self->hPcmUtils, MIN_NUMBER_OF_OUTPUT_CHANNELS,
                    &numChannel);
    flushChannels = fMin(fMax(numChannel, flushChannels), (8));
  }

  if (IS_USAC(asc->m_aot)) {
    for (int el = 0; el < (INT)asc->m_sc.m_usacConfig.m_usacNumElements; el++) {
      if (asc->m_sc.m_usacConfig.element[el].m_stereoConfigIndex == 1) {
        ascChannels--;
      } else if (asc->m_sc.m_usacConfig.element[el].m_stereoConfigIndex == 2) {
        if (asc->m_sc.m_usacConfig.element[el].m_harmonicSBR) {
          usacResidualDelayCompSamples += asc->m_samplesPerFrame;
        }
        if (asc->m_sc.m_usacConfig.m_coreSbrFrameLengthIndex == 4) {
          usacResidualDelayCompSamples +=
              6 * 16;
        }
      }
    }
  }

  aacChannelsOffset = 0;
  aacChannelsOffsetIdx = 0;
  elementOffset = 0;
  if ((ascChannels <= 0) || (ascChannels > (8)) ||
      (asc->m_channelConfiguration > AACDEC_MAX_CH_CONF)) {
    return AAC_DEC_UNSUPPORTED_CHANNELCONFIG;
  }

  flags = 0;
  { FDKmemclear(elFlags, sizeof(elFlags)); }

  if ((asc->m_channelConfiguration > 0) || IS_USAC(asc->m_aot)) {
    if (IS_USAC(asc->m_aot)) {
      self->pUsacConfig[streamIndex] = &asc->m_sc.m_usacConfig;

      if (self->pUsacConfig[streamIndex]->m_usacNumElements >
          (3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1)) {
        goto bail;
      }

      if (self->numUsacElements[streamIndex] !=
          asc->m_sc.m_usacConfig.m_usacNumElements) {
        ascChanged = 1;
      }

      if (configMode & AC_CM_ALLOC_MEM) {
        self->numUsacElements[streamIndex] =
            asc->m_sc.m_usacConfig.m_usacNumElements;
      }

      mpsEnableCurr = 0;
      for (int _el = 0;
           _el < (int)self->pUsacConfig[streamIndex]->m_usacNumElements;
           _el++) {
        int el = _el + elementOffset;
        if (self->elements[el] !=
            self->pUsacConfig[streamIndex]->element[_el].usacElementType) {
          ascChanged = 1;
        }
        if (self->usacStereoConfigIndex[el] !=
            asc->m_sc.m_usacConfig.element[_el].m_stereoConfigIndex) {
          ascChanged = 1;
        }
        if (configMode & AC_CM_ALLOC_MEM) {
          self->elements[el] =
              self->pUsacConfig[streamIndex]->element[_el].usacElementType;
          self->usacStereoConfigIndex[el] =
              asc->m_sc.m_usacConfig.element[_el].m_stereoConfigIndex;
          if (self->elements[el] == ID_USAC_CPE) {
            mpsEnableCurr |= self->usacStereoConfigIndex[el] ? 1 : 0;
          }
        }

        elFlags[el] |= (asc->m_sc.m_usacConfig.element[_el].m_noiseFilling)
                           ? AC_EL_USAC_NOISE
                           : 0;
        elFlags[el] |=
            (asc->m_sc.m_usacConfig.element[_el].m_stereoConfigIndex > 0)
                ? AC_EL_USAC_MPS212
                : 0;
        elFlags[el] |= (asc->m_sc.m_usacConfig.element[_el].m_interTes)
                           ? AC_EL_USAC_ITES
                           : 0;
        elFlags[el] |=
            (asc->m_sc.m_usacConfig.element[_el].m_pvc) ? AC_EL_USAC_PVC : 0;
        elFlags[el] |=
            (asc->m_sc.m_usacConfig.element[_el].usacElementType == ID_USAC_LFE)
                ? AC_EL_USAC_LFE
                : 0;
        elFlags[el] |=
            (asc->m_sc.m_usacConfig.element[_el].usacElementType == ID_USAC_LFE)
                ? AC_EL_LFE
                : 0;
        if ((asc->m_sc.m_usacConfig.element[_el].usacElementType ==
             ID_USAC_CPE) &&
            ((self->usacStereoConfigIndex[el] == 0))) {
          elFlags[el] |= AC_EL_USAC_CP_POSSIBLE;
        }
      }

      self->hasAudioPreRoll = 0;
      if (self->pUsacConfig[streamIndex]->m_usacNumElements) {
        self->hasAudioPreRoll = asc->m_sc.m_usacConfig.element[0]
                                    .extElement.usacExtElementHasAudioPreRoll;
      }
      if (configMode & AC_CM_ALLOC_MEM) {
        self->elements[elementOffset +
                       self->pUsacConfig[streamIndex]->m_usacNumElements] =
            ID_END;
      }
    } else {
      int i;
      for (i = 0; i < AACDEC_CH_ELEMENTS_TAB_SIZE; i++) {
        self->elements[i] = elementsTab[asc->m_channelConfiguration - 1][i];
      }
      for (; i < (3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1); i++) {
        self->elements[i] = ID_NONE;
      }
    }

    {
      int ch;

      for (ch = 0; ch < ascChannels; ch++) {
        self->chMapping[ch] = ch;
      }
      for (; ch < (8); ch++) {
        self->chMapping[ch] = 255;
      }
    }

    self->chMapIndex = asc->m_channelConfiguration;
  } else {
    if (CProgramConfig_IsValid(&asc->m_progrConfigElement)) {
      pcmDmx_SetMatrixMixdownFromPce(
          self->hPcmUtils, asc->m_progrConfigElement.MatrixMixdownIndexPresent,
          asc->m_progrConfigElement.MatrixMixdownIndex,
          asc->m_progrConfigElement.PseudoSurroundEnable);
    }
  }

  self->streamInfo.channelConfig = asc->m_channelConfiguration;

  if (self->streamInfo.aot != asc->m_aot) {
    if (configMode & AC_CM_ALLOC_MEM) {
      self->streamInfo.aot = asc->m_aot;
    }
    ascChanged = 1;
  }

  if (asc->m_aot == AOT_ER_AAC_ELD &&
      asc->m_sc.m_eldSpecificConfig.m_downscaledSamplingFrequency != 0) {
    if (self->samplingRateInfo[0].samplingRate !=
            asc->m_sc.m_eldSpecificConfig.m_downscaledSamplingFrequency ||
        self->samplingRateInfo[0].samplingRate * self->downscaleFactor !=
            asc->m_samplingFrequency) {
      downscaleFactorInBS =
          asc->m_samplingFrequency /
          asc->m_sc.m_eldSpecificConfig.m_downscaledSamplingFrequency;
      if ((downscaleFactorInBS == 1 || downscaleFactorInBS == 2 ||
           (downscaleFactorInBS == 3 &&
            asc->m_sc.m_eldSpecificConfig.m_frameLengthFlag) ||
           downscaleFactorInBS == 4) &&
          ((asc->m_samplingFrequency %
            asc->m_sc.m_eldSpecificConfig.m_downscaledSamplingFrequency) ==
           0)) {
        downscaleFactor = downscaleFactorInBS;
      } else {
        downscaleFactorInBS = 1;
        downscaleFactor = 1;
      }
    }
  } else {
    downscaleFactorInBS = 1;
    downscaleFactor = 1;
  }

  if (self->downscaleFactorInBS != downscaleFactorInBS) {
    if (configMode & AC_CM_ALLOC_MEM) {
      self->downscaleFactorInBS = downscaleFactorInBS;
      self->downscaleFactor = downscaleFactor;
    }
    ascChanged = 1;
  }

  if ((INT)asc->m_samplesPerFrame % downscaleFactor != 0) {
    return AAC_DEC_UNSUPPORTED_SAMPLINGRATE;
  }

  self->streamInfo.bitRate = 0;

  if (asc->m_aot == AOT_ER_AAC_ELD) {
    if (self->useLdQmfTimeAlign !=
        asc->m_sc.m_eldSpecificConfig.m_useLdQmfTimeAlign) {
      ascChanged = 1;
    }
    if (configMode & AC_CM_ALLOC_MEM) {
      self->useLdQmfTimeAlign =
          asc->m_sc.m_eldSpecificConfig.m_useLdQmfTimeAlign;
    }
    if (sbrEnabled != asc->m_sbrPresentFlag) {
      ascChanged = 1;
    }
  }

  self->streamInfo.extAot = asc->m_extensionAudioObjectType;
  if (self->streamInfo.extSamplingRate !=
      (INT)asc->m_extensionSamplingFrequency) {
    ascChanged = 1;
  }
  if (configMode & AC_CM_ALLOC_MEM) {
    self->streamInfo.extSamplingRate = asc->m_extensionSamplingFrequency;
  }
  flags |= (asc->m_sbrPresentFlag) ? AC_SBR_PRESENT : 0;
  flags |= (asc->m_psPresentFlag) ? AC_PS_PRESENT : 0;
  if (asc->m_sbrPresentFlag) {
    sbrEnabled = 1;
    sbrEnabledPrev = 1;
  } else {
    sbrEnabled = 0;
    sbrEnabledPrev = 0;
  }
  if (sbrEnabled && asc->m_extensionSamplingFrequency) {
    if (downscaleFactor != 1 && (downscaleFactor)&1) {
      return AAC_DEC_UNSUPPORTED_SAMPLINGRATE;
    }
    if (configMode & AC_CM_ALLOC_MEM) {
      self->streamInfo.extSamplingRate =
          self->streamInfo.extSamplingRate / self->downscaleFactor;
    }
  }
  if ((asc->m_aot == AOT_AAC_LC) && (asc->m_sbrPresentFlag == 1) &&
      (asc->m_extensionSamplingFrequency > (2 * asc->m_samplingFrequency))) {
    return AAC_DEC_UNSUPPORTED_SAMPLINGRATE;
  }

  flags |= (asc->m_vcb11Flag) ? AC_ER_VCB11 : 0;
  flags |= (asc->m_rvlcFlag) ? AC_ER_RVLC : 0;
  flags |= (asc->m_hcrFlag) ? AC_ER_HCR : 0;

  if (asc->m_aot == AOT_ER_AAC_ELD) {
    mpsEnableCurr = 0;
    flags |= AC_ELD;
    flags |= (asc->m_sbrPresentFlag)
                 ? AC_SBR_PRESENT
                 : 0;
    flags |= (asc->m_sc.m_eldSpecificConfig.m_sbrCrcFlag) ? AC_SBRCRC : 0;
    flags |= (asc->m_sc.m_eldSpecificConfig.m_useLdQmfTimeAlign)
                 ? AC_MPS_PRESENT
                 : 0;
    if (self->mpsApplicable) {
      mpsEnableCurr = asc->m_sc.m_eldSpecificConfig.m_useLdQmfTimeAlign;
    }
  }
  flags |= (asc->m_aot == AOT_ER_AAC_LD) ? AC_LD : 0;
  flags |= (asc->m_epConfig >= 0) ? AC_ER : 0;

  if (asc->m_aot == AOT_USAC) {
    flags |= AC_USAC;
    flags |= (asc->m_sc.m_usacConfig.element[0].m_stereoConfigIndex > 0)
                 ? AC_MPS_PRESENT
                 : 0;
  }
  if (asc->m_aot == AOT_DRM_AAC) {
    flags |= AC_DRM | AC_SBRCRC | AC_SCALABLE;
  }
  if (asc->m_aot == AOT_DRM_SURROUND) {
    flags |= AC_DRM | AC_SBRCRC | AC_SCALABLE | AC_MPS_PRESENT;
    FDK_ASSERT(!asc->m_psPresentFlag);
  }
  if ((asc->m_aot == AOT_AAC_SCAL) || (asc->m_aot == AOT_ER_AAC_SCAL)) {
    flags |= AC_SCALABLE;
  }

  if ((asc->m_epConfig >= 0) && (asc->m_channelConfiguration <= 0)) {
    return AAC_DEC_UNSUPPORTED_ER_FORMAT;
  }

  self->streamInfo.epConfig = asc->m_epConfig;

  if (asc->m_epConfig > 1) return AAC_DEC_UNSUPPORTED_ER_FORMAT;

  if ((self->samplingRateInfo[streamIndex].samplingRate !=
       asc->m_samplingFrequency) ||
      (self->streamInfo.aacSamplesPerFrame !=
       (INT)asc->m_samplesPerFrame / downscaleFactor)) {
    AAC_DECODER_ERROR error;

    ascChanged = 1;

    if (configMode & AC_CM_ALLOC_MEM) {
      error = getSamplingRateInfo(
          &self->samplingRateInfo[streamIndex], asc->m_samplesPerFrame,
          asc->m_samplingFrequencyIndex, asc->m_samplingFrequency);
      if (error != AAC_DEC_OK) {
        return error;
      }
      self->streamInfo.aacSampleRate =
          self->samplingRateInfo[0].samplingRate / self->downscaleFactor;
      self->streamInfo.aacSamplesPerFrame =
          asc->m_samplesPerFrame / self->downscaleFactor;
      if (self->streamInfo.aacSampleRate <= 0) {
        return AAC_DEC_UNSUPPORTED_SAMPLINGRATE;
      }
    }
  }

  if (self->ascChannels[streamIndex] != ascChannels) {
    ascChanged = 1;
  }

  if (configMode & AC_CM_DET_CFG_CHANGE) {
    if (ascChanged != 0) {
      *configChanged = 1;
    }

    CAacDecoder_AcceptFlags(self, asc, flags, elFlags, streamIndex,
                            elementOffset);

    return err;
  }

  switch (asc->m_aot) {
    case AOT_USAC:
      if (sbrEnabled) {
        for (int _el = 0;
             _el < (int)self->pUsacConfig[streamIndex]->m_usacNumElements;
             _el++) {
          int el = elementOffset + _el;
          if (IS_USAC_CHANNEL_ELEMENT(self->elements[el])) {
            if (usacStereoConfigIndex < 0) {
              usacStereoConfigIndex = self->usacStereoConfigIndex[el];
            } else {
              if ((usacStereoConfigIndex != self->usacStereoConfigIndex[el]) ||
                  (self->usacStereoConfigIndex[el] > 0)) {
                goto bail;
              }
            }
          }
        }

        if (usacStereoConfigIndex < 0) {
          goto bail;
        }

        if (usacStereoConfigIndex == 3) {
          flags |= AC_USAC_SCFGI3;
        }
      }
      break;
    default:
      break;
  }

  if (*configChanged) {
    switch (asc->m_aot) {
      case AOT_USAC:
        if (sbrEnabled) {
          const UCHAR map_sbrRatio_2_nAnaBands[] = {16, 24, 32};

          FDK_ASSERT(asc->m_sc.m_usacConfig.m_sbrRatioIndex > 0);
          FDK_ASSERT(streamIndex == 0);

          self->qmfDomain.globalConf.nInputChannels_requested = ascChannels;
          self->qmfDomain.globalConf.nOutputChannels_requested =
              (usacStereoConfigIndex == 1) ? 2 : ascChannels;
          self->qmfDomain.globalConf.flags_requested = 0;
          self->qmfDomain.globalConf.nBandsAnalysis_requested =
              map_sbrRatio_2_nAnaBands[asc->m_sc.m_usacConfig.m_sbrRatioIndex -
                                       1];
          self->qmfDomain.globalConf.nBandsSynthesis_requested = 64;
          self->qmfDomain.globalConf.nQmfTimeSlots_requested =
              (asc->m_sc.m_usacConfig.m_sbrRatioIndex == 1) ? 64 : 32;
          self->qmfDomain.globalConf.nQmfOvTimeSlots_requested =
              (asc->m_sc.m_usacConfig.m_sbrRatioIndex == 1) ? 12 : 6;
          self->qmfDomain.globalConf.nQmfProcBands_requested = 64;
          self->qmfDomain.globalConf.nQmfProcChannels_requested = 1;
          self->qmfDomain.globalConf.parkChannel =
              (usacStereoConfigIndex == 3) ? 1 : 0;
          self->qmfDomain.globalConf.parkChannel_requested =
              (usacStereoConfigIndex == 3) ? 1 : 0;
          self->qmfDomain.globalConf.qmfDomainExplicitConfig = 1;
        }
        break;
      case AOT_ER_AAC_ELD:
        if (mpsEnableCurr &&
            asc->m_sc.m_eldSpecificConfig.m_useLdQmfTimeAlign) {
          SAC_INPUT_CONFIG sac_interface = (sbrEnabled && self->hSbrDecoder)
                                               ? SAC_INTERFACE_QMF
                                               : SAC_INTERFACE_TIME;
          mpegSurroundDecoder_ConfigureQmfDomain(
              (CMpegSurroundDecoder *)self->pMpegSurroundDecoder, sac_interface,
              (UINT)self->streamInfo.aacSampleRate, asc->m_aot);
          self->qmfDomain.globalConf.qmfDomainExplicitConfig = 1;
        }
        break;
      default:
        self->qmfDomain.globalConf.qmfDomainExplicitConfig =
            0;
        break;
    }

    {
      int ch = aacChannelsOffset;
      for (int _ch = 0; _ch < ascChannels; _ch++) {
        if (ch >= (8)) {
          goto bail;
        }
        self->pAacDecoderChannelInfo[ch] = GetAacDecoderChannelInfo(ch);
        if (self->pAacDecoderChannelInfo[ch] == NULL) {
          goto bail;
        }
        ch++;
      }

      int chIdx = aacChannelsOffsetIdx;
      ch = aacChannelsOffset;
      int _numElements;
      _numElements = (((8)) + (8));
      if (flags & (AC_RSV603DA | AC_USAC)) {
        _numElements = (int)asc->m_sc.m_usacConfig.m_usacNumElements;
      }
      for (int _el = 0; _el < _numElements; _el++) {
        int el_channels = 0;
        int el = elementOffset + _el;

        if (flags &
            (AC_ER | AC_LD | AC_ELD | AC_RSV603DA | AC_USAC | AC_RSVD50)) {
          if (ch >= ascChannels) {
            break;
          }
        }

        switch (self->elements[el]) {
          case ID_SCE:
          case ID_CPE:
          case ID_LFE:
          case ID_USAC_SCE:
          case ID_USAC_CPE:
          case ID_USAC_LFE:

            el_channels = CAacDecoder_GetELChannels(
                self->elements[el], self->usacStereoConfigIndex[el]);

            {
              self->pAacDecoderChannelInfo[ch]->pComStaticData =
                  (CAacDecoderCommonStaticData *)FDKcalloc(
                      1, sizeof(CAacDecoderCommonStaticData));
              if (self->pAacDecoderChannelInfo[ch]->pComStaticData == NULL) {
                goto bail;
              }
              if (ch == aacChannelsOffset) {
                self->pAacDecoderChannelInfo[ch]->pComData =
                    (CAacDecoderCommonData *)GetWorkBufferCore6();
                self->pAacDecoderChannelInfo[ch]
                    ->pComStaticData->pWorkBufferCore1 = GetWorkBufferCore1();
              } else {
                self->pAacDecoderChannelInfo[ch]->pComData =
                    (CAacDecoderCommonData *)FDKaalloc(
                        sizeof(CAacDecoderCommonData), ALIGNMENT_DEFAULT);
                self->pAacDecoderChannelInfo[ch]
                    ->pComStaticData->pWorkBufferCore1 =
                    self->pAacDecoderChannelInfo[aacChannelsOffset]
                        ->pComStaticData->pWorkBufferCore1;
              }
              if ((self->pAacDecoderChannelInfo[ch]->pComData == NULL) ||
                  (self->pAacDecoderChannelInfo[ch]
                       ->pComStaticData->pWorkBufferCore1 == NULL)) {
                goto bail;
              }
              self->pAacDecoderChannelInfo[ch]->pDynData =
                  &(self->pAacDecoderChannelInfo[ch]
                        ->pComData->pAacDecoderDynamicData[0]);
              self->pAacDecoderChannelInfo[ch]->pSpectralCoefficient =
                  (SPECTRAL_PTR)&self->workBufferCore2[ch * 1024];

              if (el_channels == 2) {
                if (ch >= (8) - 1) {
                  return AAC_DEC_UNSUPPORTED_CHANNELCONFIG;
                }
                self->pAacDecoderChannelInfo[ch + 1]->pComData =
                    self->pAacDecoderChannelInfo[ch]->pComData;
                self->pAacDecoderChannelInfo[ch + 1]->pComStaticData =
                    self->pAacDecoderChannelInfo[ch]->pComStaticData;
                self->pAacDecoderChannelInfo[ch + 1]
                    ->pComStaticData->pWorkBufferCore1 =
                    self->pAacDecoderChannelInfo[ch]
                        ->pComStaticData->pWorkBufferCore1;
                self->pAacDecoderChannelInfo[ch + 1]->pDynData =
                    &(self->pAacDecoderChannelInfo[ch]
                          ->pComData->pAacDecoderDynamicData[1]);
                self->pAacDecoderChannelInfo[ch + 1]->pSpectralCoefficient =
                    (SPECTRAL_PTR)&self->workBufferCore2[(ch + 1) * 1024];
              }

              ch += el_channels;
            }
            chIdx += el_channels;
            break;

          default:
            break;
        }

        if (self->elements[el] == ID_END) {
          break;
        }

        el++;
      }

      chIdx = aacChannelsOffsetIdx;
      ch = aacChannelsOffset;
      for (int _ch = 0; _ch < ascChannels; _ch++) {
        {
          self->pAacDecoderStaticChannelInfo[ch] =
              GetAacDecoderStaticChannelInfo(ch);
          if (self->pAacDecoderStaticChannelInfo[ch] == NULL) {
            goto bail;
          }
          self->pAacDecoderStaticChannelInfo[ch]->pOverlapBuffer =
              GetOverlapBuffer(ch);
          if (self->pAacDecoderStaticChannelInfo[ch]->pOverlapBuffer == NULL) {
            goto bail;
          }
          if (flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA)) {
            self->pAacDecoderStaticChannelInfo[ch]->hArCo = CArco_Create();
            if (self->pAacDecoderStaticChannelInfo[ch]->hArCo == NULL) {
              goto bail;
            }
          }

          if (!(flags & (AC_USAC | AC_RSV603DA))) {
            CPns_UpdateNoiseState(
                &self->pAacDecoderChannelInfo[ch]->data.aac.PnsData,
                &self->pAacDecoderStaticChannelInfo[ch]->pnsCurrentSeed,
                self->pAacDecoderChannelInfo[ch]->pComData->pnsRandomSeed);
          }
          ch++;
        }
        chIdx++;
      }

      if (flags & AC_USAC) {
        for (int _ch = 0; _ch < flushChannels; _ch++) {
          ch = aacChannelsOffset + _ch;
          if (self->pTimeDataFlush[ch] == NULL) {
            self->pTimeDataFlush[ch] = GetTimeDataFlush(ch);
            if (self->pTimeDataFlush[ch] == NULL) {
              goto bail;
            }
          }
        }
      }

      if (flags & (AC_USAC | AC_RSV603DA)) {
        int complexStereoPredPossible = 0;
        ch = aacChannelsOffset;
        chIdx = aacChannelsOffsetIdx;
        for (int _el2 = 0; _el2 < (int)asc->m_sc.m_usacConfig.m_usacNumElements;
             _el2++) {
          int el2 = elementOffset + _el2;
          int elCh = 0, ch2;

          if ((self->elements[el2] == ID_USAC_CPE) &&
              !(self->usacStereoConfigIndex[el2] == 1)) {
            elCh = 2;
          } else if (IS_CHANNEL_ELEMENT(self->elements[el2])) {
            elCh = 1;
          }

          if (elFlags[el2] & AC_EL_USAC_CP_POSSIBLE) {
            complexStereoPredPossible = 1;
            if (self->cpeStaticData[el2] == NULL) {
              self->cpeStaticData[el2] = GetCpePersistentData();
              if (self->cpeStaticData[el2] == NULL) {
                goto bail;
              }
            }
          }

          for (ch2 = 0; ch2 < elCh; ch2++) {
            self->pAacDecoderStaticChannelInfo[ch]->pCpeStaticData =
                self->cpeStaticData[el2];
            if (self->pAacDecoderStaticChannelInfo[ch]->pCpeStaticData !=
                NULL) {
              self->pAacDecoderStaticChannelInfo[ch]
                  ->pCpeStaticData->jointStereoPersistentData
                  .spectralCoeffs[ch2] =
                  self->pAacDecoderStaticChannelInfo[ch]
                      ->concealmentInfo.spectralCoefficient;
              self->pAacDecoderStaticChannelInfo[ch]
                  ->pCpeStaticData->jointStereoPersistentData.specScale[ch2] =
                  self->pAacDecoderStaticChannelInfo[ch]
                      ->concealmentInfo.specScale;
              self->pAacDecoderStaticChannelInfo[ch]
                  ->pCpeStaticData->jointStereoPersistentData.scratchBuffer =
                  (FIXP_DBL *)self->pTimeData2;
            }
            chIdx++;
            ch++;
          }
          if (complexStereoPredPossible && (elCh == 2)) {
            if (self->pAacDecoderChannelInfo[ch - 1]
                    ->pComStaticData->cplxPredictionData == NULL) {
              self->pAacDecoderChannelInfo[ch - 1]
                  ->pComStaticData->cplxPredictionData =
                  GetCplxPredictionData();
            }
            if (self->pAacDecoderChannelInfo[ch - 1]
                    ->pComStaticData->cplxPredictionData == NULL) {
              goto bail;
            }
          }
          if (elCh > 0) {
            self->pAacDecoderStaticChannelInfo[ch - elCh]->nfRandomSeed =
                (ULONG)0x3039;
            if (self->elements[el2] == ID_USAC_CPE) {
              if (asc->m_sc.m_usacConfig.element[el2].m_stereoConfigIndex !=
                  1) {
                self->pAacDecoderStaticChannelInfo[ch - elCh + 1]
                    ->nfRandomSeed = (ULONG)0x10932;
              }
            }
          }
        }
      }

      if (ascChannels != self->aacChannels) {
        self->aacChannels = aacChannelsOffset + ch;
      }
    }

    if (usacResidualDelayCompSamples) {
      INT delayErr = FDK_Delay_Create(&self->usacResidualDelay,
                                      (USHORT)usacResidualDelayCompSamples, 1);
      if (delayErr) {
        goto bail;
      }
    }

    self->ascChannels[streamIndex] = ascChannels;
    self->aacChannelsPrev = 0;
  }

  if (self->pAacDecoderChannelInfo[0] != NULL) {
    self->pDrmBsBuffer = self->pAacDecoderChannelInfo[0]
                             ->pComStaticData->pWorkBufferCore1->DrmBsBuffer;
    self->drmBsBufferSize = DRM_BS_BUFFER_SIZE;
  }

  if (*configChanged) {
    int ch = 0;
    int chIdx = 0;
    for (int _ch = 0; _ch < self->ascChannels[streamIndex]; _ch++) {
      switch (self->streamInfo.aot) {
        case AOT_ER_AAC_ELD:
        case AOT_ER_AAC_LD:
          self->pAacDecoderChannelInfo[ch]->granuleLength =
              self->streamInfo.aacSamplesPerFrame;
          break;
        default:
          self->pAacDecoderChannelInfo[ch]->granuleLength =
              self->streamInfo.aacSamplesPerFrame / 8;
          break;
      }
      self->pAacDecoderChannelInfo[ch]->renderMode = initRenderMode;

      mdct_init(&self->pAacDecoderStaticChannelInfo[ch]->IMdct,
                self->pAacDecoderStaticChannelInfo[ch]->pOverlapBuffer,
                OverlapBufferSize);

      self->pAacDecoderStaticChannelInfo[ch]->last_core_mode = FD_LONG;
      self->pAacDecoderStaticChannelInfo[ch]->last_lpd_mode = 255;

      self->pAacDecoderStaticChannelInfo[ch]->last_tcx_pitch = L_DIV;

      aacDecoder_drcInitChannelData(
          &self->pAacDecoderStaticChannelInfo[ch]->drcData);

      pcmDmx_Reset(self->hPcmUtils, PCMDMX_RESET_BS_DATA);

      CConcealment_InitChannelData(
          &self->pAacDecoderStaticChannelInfo[ch]->concealmentInfo,
          &self->concealCommonData, initRenderMode,
          self->streamInfo.aacSamplesPerFrame);
      ch++;
      chIdx++;
    }
  }

  if (*configChanged) {
    int drcDecSampleRate, drcDecFrameSize;

    if (self->streamInfo.extSamplingRate != 0) {
      drcDecSampleRate = self->streamInfo.extSamplingRate;
      drcDecFrameSize = (self->streamInfo.aacSamplesPerFrame *
                         self->streamInfo.extSamplingRate) /
                        self->streamInfo.aacSampleRate;
    } else {
      drcDecSampleRate = self->streamInfo.aacSampleRate;
      drcDecFrameSize = self->streamInfo.aacSamplesPerFrame;
    }

    if (FDK_drcDec_Init(self->hUniDrcDecoder, drcDecFrameSize, drcDecSampleRate,
                        self->aacChannels) != 0)
      goto bail;
  }

  if (*configChanged) {
    if (asc->m_aot == AOT_USAC) {
      aacDecoder_drcDisable(self->hDrcInfo);
    }
  }

  if (asc->m_aot == AOT_USAC) {
    pcmLimiter_SetAttack(self->hLimiter, (5));
    pcmLimiter_SetThreshold(self->hLimiter, FL2FXCONST_DBL(0.89125094f));
  }

  CAacDecoder_AcceptFlags(self, asc, flags, elFlags, streamIndex,
                          elementOffset);
  self->sbrEnabled = sbrEnabled;
  self->sbrEnabledPrev = sbrEnabledPrev;
  self->mpsEnableCurr = mpsEnableCurr;

  self->streamInfo.flags = self->flags[0];

  return err;

bail:
  CAacDecoder_DeInit(self, 0);
  return AAC_DEC_OUT_OF_MEMORY;
}

AAC_DECODER_ERROR CAacDecoder_DecodeFrame(
    HANDLE_AACDECODER self, const UINT flags, PCM_DEC *pTimeData,
    const INT timeDataSize, const int timeDataChannelOffset) {
  AAC_DECODER_ERROR ErrorStatus = AAC_DEC_OK;

  CProgramConfig *pce;
  HANDLE_FDK_BITSTREAM bs = transportDec_GetBitstream(self->hInput, 0);

  MP4_ELEMENT_ID type = ID_NONE;
  INT aacChannels = 0;
  const int streamIndex = 0;

  INT auStartAnchor = (INT)FDKgetValidBits(
      bs);

  INT checkSampleRate = self->streamInfo.aacSampleRate;

  INT CConceal_TDFading_Applied[(8)] = {
      0};

  if (self->aacChannels <= 0) {
    return AAC_DEC_UNSUPPORTED_CHANNELCONFIG;
  }

  if ((transportDec_GetAuBitsRemaining(self->hInput, 0) < 15) &&
      (flags & (AACDEC_CONCEAL | AACDEC_FLUSH)) == 0) {
    self->frameOK = 0;
    ErrorStatus = AAC_DEC_DECODE_FRAME_ERROR;
  }

  pce = &self->pce;
  CProgramConfig_Reset(pce);

  CAacDecoder_AncDataReset(&self->ancData);
  if (!(flags & (AACDEC_CONCEAL | AACDEC_FLUSH)) &&
      !(self->flags[0] & (AC_USAC | AC_RSV603DA))) {
    int ch;
    if (self->streamInfo.channelConfig == 0) {
      for (ch = 0; ch < (8); ch++) {
        self->chMapping[ch] = 255;
      }
      if (!CProgramConfig_IsValid(pce)) {
        int el;
        for (el = 0; el < (3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1);
             el++) {
          self->elements[el] = ID_NONE;
        }
      }
    }
  }

  if (self->downscaleFactor > 1 && (self->flags[0] & AC_ELD)) {
    self->flags[0] |= AC_ELD_DOWNSCALE;
  } else {
    self->flags[0] &= ~AC_ELD_DOWNSCALE;
  }
  if (self->downscaleFactorInBS > 1 &&
      (self->flags[0] & AC_ELD_DOWNSCALE) == 0) {
    checkSampleRate =
        self->streamInfo.aacSampleRate / self->downscaleFactorInBS;
  }

  if (self->streamInfo.aacSampleRate <= 0) {
    return AAC_DEC_UNSUPPORTED_SAMPLINGRATE;
  }
  switch (checkSampleRate) {
    case 96000:
    case 88200:
    case 64000:
    case 16000:
    case 12000:
    case 11025:
    case 8000:
    case 7350:
    case 48000:
    case 44100:
    case 32000:
    case 24000:
    case 22050:
      break;
    default:
      if (!(self->flags[0] & (AC_USAC | AC_RSVD50 | AC_RSV603DA))) {
        return AAC_DEC_UNSUPPORTED_SAMPLINGRATE;
      }
      break;
  }

  if (flags & AACDEC_CLRHIST) {
    if (!(self->flags[0] & AC_USAC)) {
      int ch;
      for (ch = 0; ch < self->aacChannels; ch++) {
        CConcealment_InitChannelData(
            &self->pAacDecoderStaticChannelInfo[ch]->concealmentInfo,
            &self->concealCommonData,
            self->pAacDecoderChannelInfo[0]->renderMode,
            self->streamInfo.aacSamplesPerFrame);
        FDKmemclear(self->pAacDecoderStaticChannelInfo[ch]->pOverlapBuffer,
                    OverlapBufferSize * sizeof(FIXP_DBL));
      }
      if (self->streamInfo.channelConfig > 0) {
        CProgramConfig_Init(pce);
      }
    }
  }

  int pceRead = 0;

  INT hdaacDecoded = 0;
  MP4_ELEMENT_ID previous_element =
      ID_END;
  UCHAR previous_element_index = 0;
  int element_count =
      0;
  int channel_element_count = 0;
  MP4_ELEMENT_ID
  channel_elements[(3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) +
                    1)];
  int el_cnt[ID_LAST] = {0};
  int element_count_prev_streams =
      0;

  while ((type != ID_END) && (!(flags & (AACDEC_CONCEAL | AACDEC_FLUSH))) &&
         self->frameOK) {
    int el_channels;

    if (!(self->flags[0] &
          (AC_USAC | AC_RSVD50 | AC_RSV603DA | AC_ELD | AC_SCALABLE | AC_ER)))
      type = (MP4_ELEMENT_ID)FDKreadBits(bs, 3);
    else {
      if (element_count >= (3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1)) {
        self->frameOK = 0;
        ErrorStatus = AAC_DEC_PARSE_ERROR;
        break;
      }
      type = self->elements[element_count];
    }

    if ((self->flags[streamIndex] & (AC_USAC | AC_RSVD50) &&
         element_count == 0) ||
        (self->flags[streamIndex] & AC_RSV603DA)) {
      self->flags[streamIndex] &= ~AC_INDEP;

      if (FDKreadBit(bs)) {
        self->flags[streamIndex] |= AC_INDEP;
      }

      int ch = aacChannels;
      for (int chIdx = aacChannels; chIdx < self->ascChannels[streamIndex];
           chIdx++) {
        {
          if (ch >= self->aacChannels) {
            return AAC_DEC_UNKNOWN;
          }

          if ((CConcealment_GetLastFrameOk(
                   &self->pAacDecoderStaticChannelInfo[ch]->concealmentInfo,
                   1) == 0) &&
              !(self->flags[streamIndex] & AC_INDEP)) {
            self->frameOK = 0;
          }
          ch++;
        }
      }
    }

    if ((INT)FDKgetValidBits(bs) < 0) {
      self->frameOK = 0;
    }

    switch (type) {
      case ID_SCE:
      case ID_CPE:
      case ID_LFE:
      case ID_USAC_SCE:
      case ID_USAC_CPE:
      case ID_USAC_LFE:
        if (element_count >= (3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1)) {
          self->frameOK = 0;
          ErrorStatus = AAC_DEC_PARSE_ERROR;
          break;
        }

        el_channels = CAacDecoder_GetELChannels(
            type, self->usacStereoConfigIndex[element_count]);

        {
          int totalAscChannels = 0;

          for (int i = 0; i < (1 * 1); i++) {
            totalAscChannels += self->ascChannels[i];
          }
          if ((el_cnt[type] >= (totalAscChannels >> (el_channels - 1))) ||
              (aacChannels > (totalAscChannels - el_channels))) {
            ErrorStatus = AAC_DEC_DECODE_FRAME_ERROR;
            self->frameOK = 0;
            break;
          }
        }

        if (!(self->flags[streamIndex] & (AC_USAC | AC_RSVD50 | AC_RSV603DA))) {
          int ch;
          for (ch = 0; ch < el_channels; ch += 1) {
            CPns_ResetData(&self->pAacDecoderChannelInfo[aacChannels + ch]
                                ->data.aac.PnsData,
                           &self->pAacDecoderChannelInfo[aacChannels + ch]
                                ->pComData->pnsInterChannelData);
          }
        }

        if (self->frameOK) {
          ErrorStatus = CChannelElement_Read(
              bs, &self->pAacDecoderChannelInfo[aacChannels],
              &self->pAacDecoderStaticChannelInfo[aacChannels],
              self->streamInfo.aot, &self->samplingRateInfo[streamIndex],
              self->flags[streamIndex], self->elFlags[element_count],
              self->streamInfo.aacSamplesPerFrame, el_channels,
              self->streamInfo.epConfig, self->hInput);
          if (ErrorStatus != AAC_DEC_OK) {
            self->frameOK = 0;
          }
        }

        if (self->frameOK) {
          if (CProgramConfig_LookupElement(
                  pce, self->streamInfo.channelConfig,
                  self->pAacDecoderChannelInfo[aacChannels]->ElementInstanceTag,
                  aacChannels, self->chMapping, self->channelType,
                  self->channelIndices, (8), &previous_element_index,
                  self->elements, type)) {
            channel_elements[channel_element_count++] = type;
            aacChannels += el_channels;
          } else {
            self->frameOK = 0;
          }
          if (self->frameOK &&
              ((self->flags[streamIndex] & AC_SBR_PRESENT) ||
               (self->sbrEnabled == 1)) &&
              !(self->flags[streamIndex] &
                AC_USAC)
          ) {
            SBR_ERROR sbrError;
            UCHAR configMode = 0;
            UCHAR configChanged = 0;
            configMode |= AC_CM_ALLOC_MEM;

            sbrError = sbrDecoder_InitElement(
                self->hSbrDecoder, self->streamInfo.aacSampleRate,
                self->streamInfo.extSamplingRate,
                self->streamInfo.aacSamplesPerFrame, self->streamInfo.aot, type,
                previous_element_index, 2,
                0, configMode, &configChanged, self->downscaleFactor);
            if (sbrError != SBRDEC_OK) {
              self->sbrEnabled = 0;
            }
          }
        }

        el_cnt[type]++;
        if (self->frameOK && (self->flags[streamIndex] & AC_USAC) &&
            (type == ID_USAC_CPE || type == ID_USAC_SCE)) {
          ErrorStatus = aacDecoder_ParseExplicitMpsAndSbr(
              self, bs, previous_element, previous_element_index, element_count,
              el_cnt);
          if (ErrorStatus != AAC_DEC_OK) {
            self->frameOK = 0;
          }
        }
        break;

      case ID_CCE:
        if (el_cnt[type] > self->ascChannels[streamIndex]) {
          ErrorStatus = AAC_DEC_DECODE_FRAME_ERROR;
          self->frameOK = 0;
          break;
        }

        if (self->frameOK) {
          CAacDecoderCommonData commonData;
          CAacDecoderCommonStaticData commonStaticData;
          CWorkBufferCore1 workBufferCore1;
          commonStaticData.pWorkBufferCore1 = &workBufferCore1;
          C_AALLOC_SCRATCH_START(mdctSpec, FIXP_DBL, 1024);
          CAacDecoderChannelInfo tmpAacDecoderChannelInfo,
              *pTmpAacDecoderChannelInfo;

          FDKmemclear(mdctSpec, 1024 * sizeof(FIXP_DBL));

          tmpAacDecoderChannelInfo.pDynData = commonData.pAacDecoderDynamicData;
          tmpAacDecoderChannelInfo.pComData = &commonData;
          tmpAacDecoderChannelInfo.pComStaticData = &commonStaticData;
          tmpAacDecoderChannelInfo.pSpectralCoefficient =
              (SPECTRAL_PTR)mdctSpec;
          tmpAacDecoderChannelInfo.granuleLength =
              self->streamInfo.aacSamplesPerFrame / 8;
          CPns_ResetData(
              &tmpAacDecoderChannelInfo.data.aac.PnsData,
              &tmpAacDecoderChannelInfo.pComData->pnsInterChannelData);
          pTmpAacDecoderChannelInfo = &tmpAacDecoderChannelInfo;
          ErrorStatus = CChannelElement_Read(
              bs, &pTmpAacDecoderChannelInfo, NULL, self->streamInfo.aot,
              &self->samplingRateInfo[streamIndex], self->flags[streamIndex],
              AC_EL_GA_CCE, self->streamInfo.aacSamplesPerFrame, 1,
              self->streamInfo.epConfig, self->hInput);

          C_AALLOC_SCRATCH_END(mdctSpec, FIXP_DBL, 1024);

          if (ErrorStatus) {
            self->frameOK = 0;
          }

          if (self->frameOK) {
            if (CProgramConfig_LookupElement(
                    pce, self->streamInfo.channelConfig,
                    pTmpAacDecoderChannelInfo->ElementInstanceTag, 0,
                    self->chMapping, self->channelType, self->channelIndices,
                    (8), &previous_element_index, self->elements, type)) {
            } else {
              self->frameOK = 0;
            }
          }
        }
        el_cnt[type]++;
        break;

      case ID_DSE: {
        UCHAR element_instance_tag;

        CDataStreamElement_Read(self, bs, &element_instance_tag, auStartAnchor);

        if (!CProgramConfig_LookupElement(
                pce, self->streamInfo.channelConfig, element_instance_tag, 0,
                self->chMapping, self->channelType, self->channelIndices, (8),
                &previous_element_index, self->elements, type)) {
        }
      } break;

      case ID_PCE: {
        int result = CProgramConfigElement_Read(bs, self->hInput, pce,
                                                self->streamInfo.channelConfig,
                                                auStartAnchor);
        if (result < 0) {
          ErrorStatus = AAC_DEC_PARSE_ERROR;
          self->frameOK = 0;
        } else if (result > 1) {
          int elIdx = CProgramConfig_GetElementTable(
              pce, self->elements, (((8)) + (8)), &self->chMapIndex);
          for (; elIdx < (3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1);
               elIdx++) {
            self->elements[elIdx] = ID_NONE;
          }
          self->ascChannels[streamIndex] = pce->NumChannels;
          if (element_count != 0) {
            self->frameOK = 0;
          }
        }
        pceRead = (result >= 0) ? 1 : 0;
      } break;

      case ID_FIL: {
        int bitCnt = FDKreadBits(bs, 4);

        if (bitCnt == 15) {
          int esc_count = FDKreadBits(bs, 8);
          bitCnt = esc_count + 14;
        }

        bitCnt <<= 3;

        while (bitCnt > 0) {
          ErrorStatus = CAacDecoder_ExtPayloadParse(
              self, bs, &bitCnt, previous_element, previous_element_index, 1);
          if (ErrorStatus != AAC_DEC_OK) {
            self->frameOK = 0;
            break;
          }
        }
      } break;

      case ID_EXT:
        if (element_count >= (3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1)) {
          self->frameOK = 0;
          ErrorStatus = AAC_DEC_PARSE_ERROR;
          break;
        }

        ErrorStatus = aacDecoder_ParseExplicitMpsAndSbr(
            self, bs, previous_element, previous_element_index, element_count,
            el_cnt);
        break;

      case ID_USAC_EXT: {
        if ((element_count - element_count_prev_streams) >=
            TP_USAC_MAX_ELEMENTS) {
          self->frameOK = 0;
          ErrorStatus = AAC_DEC_PARSE_ERROR;
          break;
        }
        int usacExtElementPayloadLength;

        if (FDKreadBit(bs)) {
          if (FDKreadBit(bs)) {
            usacExtElementPayloadLength =
                self->pUsacConfig[streamIndex]
                    ->element[element_count - element_count_prev_streams]
                    .extElement.usacExtElementDefaultLength;
          } else {
            usacExtElementPayloadLength = FDKreadBits(bs, 8);
            if (usacExtElementPayloadLength == (UINT)(1 << 8) - 1) {
              UINT valueAdd = FDKreadBits(bs, 16);
              usacExtElementPayloadLength += (INT)valueAdd - 2;
            }
          }
          if (usacExtElementPayloadLength > 0) {
            int usacExtBitPos;

            if (self->pUsacConfig[streamIndex]
                    ->element[element_count - element_count_prev_streams]
                    .extElement.usacExtElementPayloadFrag) {
              FDKreadBit(bs);
              FDKreadBit(bs);
            }

            usacExtBitPos = (INT)FDKgetValidBits(bs);

            USAC_EXT_ELEMENT_TYPE usacExtElementType =
                self->pUsacConfig[streamIndex]
                    ->element[element_count - element_count_prev_streams]
                    .extElement.usacExtElementType;

            switch (usacExtElementType) {
              case ID_EXT_ELE_UNI_DRC:
                if (streamIndex == 0) {
                  int drcErr;

                  drcErr = FDK_drcDec_ReadUniDrcGain(self->hUniDrcDecoder, bs);
                  if (drcErr != 0) {
                    ErrorStatus = AAC_DEC_PARSE_ERROR;
                  }
                }
                break;

              default:
                break;
            }

            usacExtBitPos = (usacExtElementPayloadLength * 8) -
                            (usacExtBitPos - (INT)FDKgetValidBits(bs));
            if (usacExtBitPos < 0) {
              self->frameOK = 0;
              ErrorStatus = AAC_DEC_PARSE_ERROR;
            }
            FDKpushBiDirectional(bs, usacExtBitPos);
          }
        }
      } break;
      case ID_END:
      case ID_USAC_END:
        break;

      default:
        ErrorStatus = AAC_DEC_DECODE_FRAME_ERROR;
        self->frameOK = 0;
        break;
    }

    previous_element = type;
    element_count++;
  }

  if (!(self->flags[streamIndex] &
        (AC_USAC | AC_RSVD50 | AC_RSV603DA | AC_BSAC | AC_LD | AC_ELD | AC_ER |
         AC_SCALABLE)) &&
      (self->streamInfo.channelConfig == 0) && pce->isValid &&
      (ErrorStatus == AAC_DEC_OK) && self->frameOK &&
      !(flags & (AACDEC_CONCEAL | AACDEC_FLUSH))) {
    if (!validateElementInstanceTags(
            &self->pce, self->pAacDecoderChannelInfo, aacChannels,
            channel_elements,
            fMin(channel_element_count, (int)(sizeof(channel_elements) /
                                              sizeof(*channel_elements))))) {
      ErrorStatus = AAC_DEC_DECODE_FRAME_ERROR;
      self->frameOK = 0;
    }
  }

  if (!(flags & (AACDEC_CONCEAL | AACDEC_FLUSH))) {
    if (!(self->flags[streamIndex] & (AC_RSVD50 | AC_USAC)) ||
        (self->prerollAULength[self->accessUnit])
    ) {
      FDKbyteAlign(bs, auStartAnchor);
    }

    if (transportDec_GetAuBitsTotal(self->hInput, 0) > 0) {
      INT unreadBits = transportDec_GetAuBitsRemaining(self->hInput, 0);
      if (self->prerollAULength[self->accessUnit]) {
        unreadBits = unreadBits - transportDec_GetAuBitsTotal(self->hInput, 0) +
                     (INT)self->prerollAULength[self->accessUnit];
      }
      if (((self->flags[streamIndex] & (AC_RSVD50 | AC_USAC)) &&
           ((unreadBits < 0) || (unreadBits > 7)) &&
           !(self->prerollAULength[self->accessUnit])) ||
          ((!(self->flags[streamIndex] & (AC_RSVD50 | AC_USAC)) ||
            (self->prerollAULength[self->accessUnit])) &&
           (unreadBits != 0))) {
        if ((((unreadBits < 0) || (unreadBits > 7)) && self->frameOK) &&
            ((transportDec_GetFormat(self->hInput) == TT_DRM) &&
             (self->flags[streamIndex] & AC_USAC))) {
          self->frameOK = 1;
        } else {
          self->frameOK = 0;
        }

        if (ErrorStatus == AAC_DEC_OK && self->frameOK == 0) {
          ErrorStatus = AAC_DEC_PARSE_ERROR;
        }
        FDKpushBiDirectional(bs, unreadBits);
      }
    }

    if (self->frameOK && type != ID_END) {
      if (ErrorStatus == AAC_DEC_OK) {
        ErrorStatus = AAC_DEC_PARSE_ERROR;
      }
      self->frameOK = 0;
    }
  }

  if (!(flags & (AACDEC_CONCEAL | AACDEC_FLUSH)) && self->frameOK) {
    channel_elements[channel_element_count++] = ID_END;
  }
  element_count = 0;
  aacChannels = 0;
  type = ID_NONE;
  previous_element_index = 0;

  while (type != ID_END &&
         element_count < (3 * ((8) * 2) + (((8) * 2)) / 2 + 4 * (1) + 1)) {
    int el_channels;

    if ((flags & (AACDEC_CONCEAL | AACDEC_FLUSH)) || !self->frameOK) {
      channel_elements[element_count] = self->elements[element_count];
      if (channel_elements[element_count] == ID_NONE) {
        channel_elements[element_count] = ID_END;
      }
    }

    if (self->flags[streamIndex] & (AC_USAC | AC_RSV603DA | AC_BSAC)) {
      type = self->elements[element_count];
    } else {
      type = channel_elements[element_count];
    }

    if (!(flags & (AACDEC_CONCEAL | AACDEC_FLUSH)) && self->frameOK) {
      switch (type) {
        case ID_SCE:
        case ID_CPE:
        case ID_LFE:
        case ID_USAC_SCE:
        case ID_USAC_CPE:
        case ID_USAC_LFE:

          el_channels = CAacDecoder_GetELChannels(
              type, self->usacStereoConfigIndex[element_count]);

          if (!hdaacDecoded) {
            if (self->pAacDecoderStaticChannelInfo[aacChannels]
                    ->pCpeStaticData != NULL) {
              self->pAacDecoderStaticChannelInfo[aacChannels]
                  ->pCpeStaticData->jointStereoPersistentData.scratchBuffer =
                  (FIXP_DBL *)pTimeData;
            }
            CChannelElement_Decode(
                &self->pAacDecoderChannelInfo[aacChannels],
                &self->pAacDecoderStaticChannelInfo[aacChannels],
                &self->samplingRateInfo[streamIndex], self->flags[streamIndex],
                self->elFlags[element_count], el_channels);
          }
          aacChannels += el_channels;
          break;
        case ID_NONE:
          type = ID_END;
          break;
        default:
          break;
      }
    }
    element_count++;
  }

  if ((aacChannels == 0 || aacChannels > self->aacChannels) &&
      !(flags & (AACDEC_CONCEAL | AACDEC_FLUSH))) {
    if (ErrorStatus == AAC_DEC_OK) {
      ErrorStatus = AAC_DEC_DECODE_FRAME_ERROR;
    }
    self->frameOK = 0;
    aacChannels = 0;
  }

  if (!(flags & (AACDEC_CONCEAL | AACDEC_FLUSH))) {
    if (TRANSPORTDEC_OK != transportDec_CrcCheck(self->hInput)) {
      ErrorStatus = AAC_DEC_CRC_ERROR;
      self->frameOK = 0;
    }
  }

  if ((self->frameOK == 0) && (ErrorStatus == AAC_DEC_OK)) {
    ErrorStatus = AAC_DEC_DECODE_FRAME_ERROR;
  }

  if (self->frameOK && (flags & AACDEC_FLUSH)) {
    aacChannels = self->aacChannelsPrev;
    FDKmemcpy(self->channelType, self->channelTypePrev,
              (8) * sizeof(AUDIO_CHANNEL_TYPE));
    FDKmemcpy(self->channelIndices, self->channelIndicesPrev,
              (8) * sizeof(UCHAR));
    self->sbrEnabled = self->sbrEnabledPrev;
  } else {
    if (self->frameOK && !(flags & AACDEC_CONCEAL)) {
      self->aacChannelsPrev = aacChannels;
      FDKmemcpy(self->channelTypePrev, self->channelType,
                (8) * sizeof(AUDIO_CHANNEL_TYPE));
      FDKmemcpy(self->channelIndicesPrev, self->channelIndices,
                (8) * sizeof(UCHAR));
      self->sbrEnabledPrev = self->sbrEnabled;
    } else {
      if (self->aacChannels > 0) {
        if ((self->buildUpStatus == AACDEC_RSV60_BUILD_UP_ON) ||
            (self->buildUpStatus == AACDEC_RSV60_BUILD_UP_ON_IN_BAND) ||
            (self->buildUpStatus == AACDEC_USAC_BUILD_UP_ON)) {
          aacChannels = self->aacChannels;
          self->aacChannelsPrev = aacChannels;
        } else {
          aacChannels = self->aacChannelsPrev;
        }
        FDKmemcpy(self->channelType, self->channelTypePrev,
                  (8) * sizeof(AUDIO_CHANNEL_TYPE));
        FDKmemcpy(self->channelIndices, self->channelIndicesPrev,
                  (8) * sizeof(UCHAR));
        self->sbrEnabled = self->sbrEnabledPrev;
      }
    }
  }

  self->streamInfo.aacNumChannels = aacChannels;

  if (aacChannels == 0) {
    ErrorStatus = AAC_DEC_UNKNOWN;
  }

  if (pceRead == 1 && CProgramConfig_IsValid(pce)) {
    pcmDmx_SetMatrixMixdownFromPce(
        self->hPcmUtils, pce->MatrixMixdownIndexPresent,
        pce->MatrixMixdownIndex, pce->PseudoSurroundEnable);
  }

  if (!IS_OUTPUT_VALID(ErrorStatus)) {
    return ErrorStatus;
  }

  if ((self->streamInfo.channelConfig == 0) && !pce->isValid) {
    self->chMapIndex = (aacChannels < 7) ? aacChannels : 0;
  }

  {
    int c, cIdx;
    int mapped, fCopyChMap = 1;
    UCHAR drcChMap[(8)];

    if ((self->streamInfo.channelConfig == 0) && CProgramConfig_IsValid(pce)) {
      UCHAR tmpChMap[(8)];
      if (CProgramConfig_GetPceChMap(pce, tmpChMap, (8)) == 0) {
        for (c = 0; c < aacChannels; c += 1) {
          drcChMap[c] =
              (self->chMapping[c] == 255) ? 255 : tmpChMap[self->chMapping[c]];
        }
        fCopyChMap = 0;
      }
    }
    if (fCopyChMap != 0) {
      FDKmemcpy(drcChMap, self->chMapping, (8) * sizeof(UCHAR));
    }

    mapped = aacDecoder_drcProlog(
        self->hDrcInfo, bs, self->pAacDecoderStaticChannelInfo,
        pce->ElementInstanceTag, drcChMap, aacChannels);
    if (mapped > 0) {
      if (!(self->flags[streamIndex] & (AC_USAC | AC_RSV603DA))) {
        self->flags[streamIndex] |= AC_DRC_PRESENT;
      } else {
        ErrorStatus = AAC_DEC_UNSUPPORTED_FORMAT;
      }
    }
    if (self->flags[streamIndex] & (AC_USAC | AC_RSV603DA)) {
      aacDecoder_drcDisable(self->hDrcInfo);
    }

    UCHAR Reverse_chMapping[((8) * 2)];
    for (c = 0; c < aacChannels; c++) {
      int d;
      for (d = 0; d < aacChannels - 1; d++) {
        if (self->chMapping[d] == c) {
          break;
        }
      }
      Reverse_chMapping[c] = d;
    }

    int el;
    int el_channels;
    c = 0;
    cIdx = 0;
    el_channels = 0;
    for (el = 0; el < element_count; el++) {
      int frameOk_butConceal =
          0;
      int concealApplyReturnCode;

      if (self->flags[streamIndex] & (AC_USAC | AC_RSV603DA | AC_BSAC)) {
        type = self->elements[el];
      } else {
        type = channel_elements[el];
      }

      {
        int nElementChannels;

        nElementChannels =
            CAacDecoder_GetELChannels(type, self->usacStereoConfigIndex[el]);

        el_channels += nElementChannels;

        if (nElementChannels == 0) {
          continue;
        }
      }

      int offset;
      int elCh = 0;
      for (; cIdx < el_channels; c++, cIdx++, elCh++) {
        if (c >= aacChannels) {
          return AAC_DEC_UNKNOWN;
        }

        CAacDecoderChannelInfo *pAacDecoderChannelInfo =
            self->pAacDecoderChannelInfo[c];
        CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo =
            self->pAacDecoderStaticChannelInfo[c];

        {
          pAacDecoderStaticChannelInfo =
              self->pAacDecoderStaticChannelInfo[Reverse_chMapping[c]];
          offset =
              FDK_chMapDescr_getMapValue(
                  &self->mapDescr, Reverse_chMapping[cIdx], self->chMapIndex) *
              timeDataChannelOffset;
        }

        if (flags & AACDEC_FLUSH) {
          FDKmemclear(pAacDecoderChannelInfo->pSpectralCoefficient,
                      sizeof(FIXP_DBL) * self->streamInfo.aacSamplesPerFrame);
        }

        {
          CIcsInfo *pIcsInfo = &pAacDecoderChannelInfo->icsInfo;
          CConcealmentInfo *hConcealmentInfo =
              &pAacDecoderStaticChannelInfo->concealmentInfo;
          const int mute_release_active =
              (self->frameOK && !(flags & AACDEC_CONCEAL)) &&
              ((hConcealmentInfo->concealState >= ConcealState_Mute) &&
               (hConcealmentInfo->cntValidFrames + 1 <=
                hConcealmentInfo->pConcealParams->numMuteReleaseFrames));
          const int icsIsInvalid = (GetScaleFactorBandsTransmitted(pIcsInfo) >
                                    GetScaleFactorBandsTotal(pIcsInfo));
          const int icsInfoUsedinFadeOut =
              !(pAacDecoderChannelInfo->renderMode == AACDEC_RENDER_LPD &&
                pAacDecoderStaticChannelInfo->last_lpd_mode == 0);
          if (icsInfoUsedinFadeOut && icsIsInvalid && !mute_release_active) {
            self->frameOK = 0;
          }
        }

        {
          CAacDecoderChannelInfo **ppAacDecoderChannelInfo =
              &pAacDecoderChannelInfo;
          CAacDecoderStaticChannelInfo **ppAacDecoderStaticChannelInfo =
              &pAacDecoderStaticChannelInfo;
          {
            concealApplyReturnCode = CConcealment_Apply(
                &(*ppAacDecoderStaticChannelInfo)->concealmentInfo,
                *ppAacDecoderChannelInfo, *ppAacDecoderStaticChannelInfo,
                &self->samplingRateInfo[streamIndex],
                self->streamInfo.aacSamplesPerFrame,
                pAacDecoderStaticChannelInfo->last_lpd_mode,
                (self->frameOK && !(flags & AACDEC_CONCEAL)),
                self->flags[streamIndex]);
          }
        }
        if (concealApplyReturnCode == -1) {
          frameOk_butConceal = 1;
        }

        if (flags & (AACDEC_INTR)) {
          aacDecoder_drcInitChannelData(&pAacDecoderStaticChannelInfo->drcData);
        }
        if (flags & (AACDEC_CLRHIST)) {
          if (!(self->flags[0] & AC_USAC)) {
            aacDecoder_drcInitChannelData(
                &pAacDecoderStaticChannelInfo->drcData);
          }
        }

        self->extGain[0] = (FIXP_DBL)AACDEC_DRC_GAIN_SCALING;

        aacDecoder_drcApply(
            self->hDrcInfo, self->hSbrDecoder, pAacDecoderChannelInfo,
            &pAacDecoderStaticChannelInfo->drcData, self->extGain, c,
            self->streamInfo.aacSamplesPerFrame, self->sbrEnabled

        );

        if (timeDataSize < timeDataChannelOffset * self->aacChannels) {
          ErrorStatus = AAC_DEC_OUTPUT_BUFFER_TOO_SMALL;
          break;
        }
        if (self->flushStatus && (self->flushCnt > 0) &&
            !(flags & AACDEC_CONCEAL)) {
          FDKmemclear(pTimeData + offset,
                      sizeof(PCM_DEC) * self->streamInfo.aacSamplesPerFrame);
        } else
          switch (pAacDecoderChannelInfo->renderMode) {
            case AACDEC_RENDER_IMDCT:

              CBlock_FrequencyToTime(
                  pAacDecoderStaticChannelInfo, pAacDecoderChannelInfo,
                  pTimeData + offset, self->streamInfo.aacSamplesPerFrame,
                  (self->frameOK && !(flags & AACDEC_CONCEAL) &&
                   !frameOk_butConceal),
                  pAacDecoderChannelInfo->pComStaticData->pWorkBufferCore1
                      ->mdctOutTemp,
                  self->aacOutDataHeadroom, self->elFlags[el], elCh);

              self->extGainDelay = self->streamInfo.aacSamplesPerFrame;
              break;
            case AACDEC_RENDER_ELDFB: {
              CBlock_FrequencyToTimeLowDelay(
                  pAacDecoderStaticChannelInfo, pAacDecoderChannelInfo,
                  pTimeData + offset, self->streamInfo.aacSamplesPerFrame);
              self->extGainDelay =
                  (self->streamInfo.aacSamplesPerFrame * 2 -
                   self->streamInfo.aacSamplesPerFrame / 2 - 1) /
                  2;
            } break;
            case AACDEC_RENDER_LPD:

              CLpd_RenderTimeSignal(
                  pAacDecoderStaticChannelInfo, pAacDecoderChannelInfo,
                  pTimeData + offset, self->streamInfo.aacSamplesPerFrame,
                  &self->samplingRateInfo[streamIndex],
                  (self->frameOK && !(flags & AACDEC_CONCEAL) &&
                   !frameOk_butConceal),
                  self->aacOutDataHeadroom, flags, self->flags[streamIndex]);

              self->extGainDelay = self->streamInfo.aacSamplesPerFrame;
              break;
            default:
              ErrorStatus = AAC_DEC_UNKNOWN;
              break;
          }
        if (!CConceal_TDFading_Applied[c]) {
          CConceal_TDFading_Applied[c] = CConcealment_TDFading(
              self->streamInfo.aacSamplesPerFrame,
              &self->pAacDecoderStaticChannelInfo[c], self->aacOutDataHeadroom,
              pTimeData + offset, 0);
          if (c + 1 < (8) && c < aacChannels - 1) {
            self->pAacDecoderStaticChannelInfo[c + 1]
                ->concealmentInfo.TDNoiseSeed =
                self->pAacDecoderStaticChannelInfo[c]
                    ->concealmentInfo.TDNoiseSeed;
          }
        }
      }
    }

    if (self->flags[streamIndex] & AC_USAC) {
      int bsPseudoLr = 0;
      mpegSurroundDecoder_IsPseudoLR(
          (CMpegSurroundDecoder *)self->pMpegSurroundDecoder, &bsPseudoLr);
      if ((aacChannels == 2) && bsPseudoLr) {
        int i, offset2;
        const FIXP_SGL invSqrt2 = FL2FXCONST_SGL(0.707106781186547f);
        PCM_DEC *pTD = pTimeData;

        offset2 = timeDataChannelOffset;

        for (i = 0; i < self->streamInfo.aacSamplesPerFrame; i++) {
          FIXP_DBL L = PCM_DEC2FIXP_DBL(pTD[0]);
          FIXP_DBL R = PCM_DEC2FIXP_DBL(pTD[offset2]);
          L = fMult(L, invSqrt2);
          R = fMult(R, invSqrt2);
          pTD[0] = L + R;
          pTD[offset2] = L - R;
          pTD++;
        }
      }
    }

    mapped = aacDecoder_drcEpilog(
        self->hDrcInfo, bs, self->pAacDecoderStaticChannelInfo,
        pce->ElementInstanceTag, drcChMap, aacChannels);
    if (mapped > 0) {
      if (!(self->flags[streamIndex] & (AC_USAC | AC_RSV603DA))) {
        self->flags[streamIndex] |= AC_DRC_PRESENT;
      } else {
        ErrorStatus = AAC_DEC_UNSUPPORTED_FORMAT;
      }
    }
    if (self->flags[streamIndex] & (AC_USAC | AC_RSV603DA)) {
      aacDecoder_drcDisable(self->hDrcInfo);
    }
  }

  self->streamInfo.outputDelay +=
      CConcealment_GetDelay(&self->concealCommonData) *
      self->streamInfo.aacSamplesPerFrame;

  aacDecoder_drcGetInfo(self->hDrcInfo, &self->streamInfo.drcPresMode,
                        &self->streamInfo.drcProgRefLev);

  if (!(self->flags[0] & AC_RSV603DA)) {
    AUDIO_CHANNEL_TYPE types[(8)];
    UCHAR idx[(8)];
    int c;
    int mapValue;

    FDK_ASSERT(sizeof(self->channelType) == sizeof(types));
    FDK_ASSERT(sizeof(self->channelIndices) == sizeof(idx));

    FDKmemcpy(types, self->channelType, sizeof(types));
    FDKmemcpy(idx, self->channelIndices, sizeof(idx));

    for (c = 0; c < aacChannels; c++) {
      mapValue =
          FDK_chMapDescr_getMapValue(&self->mapDescr, c, self->chMapIndex);
      self->channelType[mapValue] = types[c];
      self->channelIndices[mapValue] = idx[c];
    }
  }

  self->blockNumber++;

  return ErrorStatus;
}

CStreamInfo * CAacDecoder_GetStreamInfo(HANDLE_AACDECODER self) {
  if (!self) {
    return NULL;
  }
  return &self->streamInfo;
}
