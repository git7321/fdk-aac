/**************************** AAC decoder library ******************************

   Author(s):   Christian Griebel

   Description: Dynamic range control (DRC) decoder tool for AAC

*******************************************************************************/

#include "aacdec_drc.h"
#include "channelinfo.h"
#include "aac_rom.h"
#include "sbrdecoder.h"

#define DRC_PARAMETER_BITS (7)
#define DRC_MAX_QUANT_STEPS (1 << DRC_PARAMETER_BITS)
#define DRC_MAX_QUANT_FACTOR (DRC_MAX_QUANT_STEPS - 1)
#define DRC_PARAM_QUANT_STEP \
  (FL2FXCONST_DBL(1.0f / (float)DRC_MAX_QUANT_FACTOR))
#define DRC_PARAM_SCALE (1)
#define DRC_SCALING_MAX \
  ((FIXP_DBL)((INT)(DRC_PARAM_QUANT_STEP >> DRC_PARAM_SCALE) * (INT)127))

#define DRC_BLOCK_LEN (1024)
#define DRC_BAND_MULT (4)
#define DRC_BLOCK_LEN_DIV_BAND_MULT (DRC_BLOCK_LEN / DRC_BAND_MULT)

#define MAX_REFERENCE_LEVEL (127)

#define DRC_HEAVY_THRESHOLD_DB (10)

#define DVB_ANC_DATA_SYNC_BYTE (0xBC)

#define OFF 0
#define ON 1

static INT convert_drcParam(FIXP_DBL param_dbl) {
  LONG param_long;

  param_long = (LONG)param_dbl >> 7;
  param_long = param_long * (INT)DRC_MAX_QUANT_FACTOR;
  param_long >>= 31 - 7 - DRC_PARAM_SCALE - 1;
  param_long += 1;
  param_long >>= 1;

  return (INT)param_long;
}

void aacDecoder_drcDisable(HANDLE_AAC_DRC self) {
  self->enable = 0;
  self->applyExtGain = 0;
  self->progRefLevelPresent = 0;
}

void aacDecoder_drcReset(HANDLE_AAC_DRC self) {
  self->applyExtGain = 0;
  self->additionalGainPrev = AACDEC_DRC_GAIN_INIT_VALUE;
  self->additionalGainFilterState = AACDEC_DRC_GAIN_INIT_VALUE;
  self->additionalGainFilterState1 = AACDEC_DRC_GAIN_INIT_VALUE;
}

void aacDecoder_drcInit(HANDLE_AAC_DRC self) {
  CDrcParams *pParams;

  if (self == NULL) {
    return;
  }

  self->enable = OFF;
  self->numThreads = 0;

  pParams = &self->params;
  pParams->bsDelayEnable = 0;
  pParams->cut = FL2FXCONST_DBL(0.0f);
  pParams->usrCut = FL2FXCONST_DBL(0.0f);
  pParams->boost = FL2FXCONST_DBL(0.0f);
  pParams->usrBoost = FL2FXCONST_DBL(0.0f);
  pParams->targetRefLevel = 96;
  pParams->expiryFrame = AACDEC_DRC_DFLT_EXPIRY_FRAMES;
  pParams->applyHeavyCompression = OFF;
  pParams->usrApplyHeavyCompression = OFF;

  pParams->defaultPresentationMode = DISABLED_PARAMETER_HANDLING;
  pParams->encoderTargetLevel = MAX_REFERENCE_LEVEL;

  self->update = 1;
  self->numOutChannels = 0;
  self->prevAacNumChannels = 0;

  self->progRefLevel = pParams->targetRefLevel;
  self->progRefLevelPresent = 0;
  self->presMode = -1;

  aacDecoder_drcReset(self);
}

void aacDecoder_drcInitChannelData(CDrcChannelData *pDrcChData) {
  if (pDrcChData != NULL) {
    pDrcChData->expiryCount = 0;
    pDrcChData->numBands = 1;
    pDrcChData->bandTop[0] = DRC_BLOCK_LEN_DIV_BAND_MULT - 1;
    pDrcChData->drcValue[0] = 0;
    pDrcChData->drcInterpolationScheme = 0;
    pDrcChData->drcDataType = UNKNOWN_PAYLOAD;
  }
}

AAC_DECODER_ERROR aacDecoder_drcSetParam(HANDLE_AAC_DRC self,
                                         AACDEC_DRC_PARAM param, INT value) {
  AAC_DECODER_ERROR ErrorStatus = AAC_DEC_OK;

  switch (param) {
    case DRC_CUT_SCALE:
      if ((value < 0) || (value > DRC_MAX_QUANT_FACTOR)) {
        return AAC_DEC_SET_PARAM_FAIL;
      }
      if (self == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      }
      self->params.usrCut = (FIXP_DBL)(
          (INT)(DRC_PARAM_QUANT_STEP >> DRC_PARAM_SCALE) * (INT)value);
      self->update = 1;
      break;
    case DRC_BOOST_SCALE:
      if ((value < 0) || (value > DRC_MAX_QUANT_FACTOR)) {
        return AAC_DEC_SET_PARAM_FAIL;
      }
      if (self == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      }
      self->params.usrBoost = (FIXP_DBL)(
          (INT)(DRC_PARAM_QUANT_STEP >> DRC_PARAM_SCALE) * (INT)value);
      self->update = 1;
      break;
    case TARGET_REF_LEVEL:
      if (value > MAX_REFERENCE_LEVEL || value < -MAX_REFERENCE_LEVEL) {
        return AAC_DEC_SET_PARAM_FAIL;
      }
      if (self == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      }
      if (value < 0) {
        self->params.targetRefLevel = -1;
      } else {
        if (self->params.targetRefLevel != (SCHAR)value) {
          self->params.targetRefLevel = (SCHAR)value;
          self->progRefLevel = (SCHAR)value;
        }
        self->update = 1;
      }
      break;
    case APPLY_HEAVY_COMPRESSION:
      if ((value != OFF) && (value != ON)) {
        return AAC_DEC_SET_PARAM_FAIL;
      }
      if (self == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      }
      self->params.usrApplyHeavyCompression = (UCHAR)value;
      self->update = 1;
      break;
    case DEFAULT_PRESENTATION_MODE:
      if (value < AAC_DRC_PARAMETER_HANDLING_DISABLED ||
          value > AAC_DRC_PRESENTATION_MODE_2_DEFAULT) {
        return AAC_DEC_SET_PARAM_FAIL;
      }
      if (self == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      }
      self->params.defaultPresentationMode =
          (AACDEC_DRC_PARAMETER_HANDLING)value;
      self->update = 1;
      break;
    case ENCODER_TARGET_LEVEL:
      if (value > MAX_REFERENCE_LEVEL || value < 0) {
        return AAC_DEC_SET_PARAM_FAIL;
      }
      if (self == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      }
      self->params.encoderTargetLevel = (UCHAR)value;
      self->update = 1;
      break;
    case DRC_BS_DELAY:
      if (value < 0 || value > 1) {
        return AAC_DEC_SET_PARAM_FAIL;
      }
      if (self == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      }
      self->params.bsDelayEnable = value;
      break;
    case DRC_DATA_EXPIRY_FRAME:
      if (self == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      }
      self->params.expiryFrame = (value > 0) ? (UINT)value : 0;
      break;
    case MAX_OUTPUT_CHANNELS:
      if (self == NULL) {
        return AAC_DEC_INVALID_HANDLE;
      }
      self->numOutChannels = (INT)value;
      self->update = 1;
      break;
    default:
      return AAC_DEC_SET_PARAM_FAIL;
  }

  return ErrorStatus;
}

static int parseExcludedChannels(UINT *excludedChnsMask,
                                 HANDLE_FDK_BITSTREAM bs) {
  UINT excludeMask = 0;
  UINT i, j;
  int bitCnt = 9;

  for (i = 0, j = 1; i < 7; i++, j <<= 1) {
    if (FDKreadBits(bs, 1)) {
      excludeMask |= j;
    }
  }

  while (FDKreadBits(bs, 1)) {
    for (i = 0; i < 7; i++, j <<= 1) {
      if (FDKreadBits(bs, 1)) {
        excludeMask |= j;
      }
    }
    bitCnt += 9;
    FDK_ASSERT(j < (UINT)-1);
  }

  *excludedChnsMask = excludeMask;

  return (bitCnt);
}

int aacDecoder_drcMarkPayload(HANDLE_AAC_DRC self, HANDLE_FDK_BITSTREAM bs,
                              AACDEC_DRC_PAYLOAD_TYPE type) {
  UINT bsStartPos;
  int i, numBands = 1, bitCnt = 0;

  if (self == NULL) {
    return 0;
  }

  bsStartPos = FDKgetValidBits(bs);

  switch (type) {
    case MPEG_DRC_EXT_DATA: {
      bitCnt = 4;

      if (FDKreadBits(bs, 1)) {
        FDKreadBits(bs, 8);
        bitCnt += 8;
      }

      if (FDKreadBits(bs, 1)) {
        FDKreadBits(bs, 7);
        bitCnt += 8;
        while (FDKreadBits(bs, 1)) {
          FDKreadBits(bs, 7);
          bitCnt += 8;
        }
      }

      if (FDKreadBits(bs, 1)) {
        numBands += FDKreadBits(bs, 4);
        FDKreadBits(bs, 4);
        bitCnt += 8;
        for (i = 0; i < numBands; i++) {
          FDKreadBits(bs, 8);
          bitCnt += 8;
        }
      }

      if (FDKreadBits(bs, 1)) {
        FDKreadBits(bs, 8);
        bitCnt += 8;
      }

      for (i = 0; i < numBands; i++) {
        FDKreadBits(bs, 8);
        bitCnt += 8;
      }

      if ((self->numPayloads < MAX_DRC_THREADS) &&
          ((INT)FDKgetValidBits(bs) >= 0)) {
        self->drcPayloadPosition[self->numPayloads++] = bsStartPos;
      }
    } break;

    case DVB_DRC_ANC_DATA:
      bitCnt += 8;
      if (FDKreadBits(bs, 8) == DVB_ANC_DATA_SYNC_BYTE) {
        int dmxLevelsPresent, compressionPresent;
        int coarseGrainTcPresent, fineGrainTcPresent;

        FDKreadBits(
            bs,
            8);
        bitCnt += 8;

        FDKreadBits(bs, 3);
        dmxLevelsPresent =
            FDKreadBits(bs, 1);
        FDKreadBits(bs, 1);
        compressionPresent =
            FDKreadBits(bs, 1);
        coarseGrainTcPresent =
            FDKreadBits(bs, 1);
        fineGrainTcPresent =
            FDKreadBits(bs, 1);
        bitCnt += 8;

        if (dmxLevelsPresent) {
          FDKreadBits(bs, 8);
          bitCnt += 8;
        }
        if (compressionPresent) {
          FDKreadBits(bs, 16);
          bitCnt += 16;
        }
        if (coarseGrainTcPresent) {
          FDKreadBits(bs, 16);
          bitCnt += 16;
        }
        if (fineGrainTcPresent) {
          FDKreadBits(bs, 16);
          bitCnt += 16;
        }
        if (!self->dvbAncDataAvailable && ((INT)FDKgetValidBits(bs) >= 0)) {
          self->dvbAncDataPosition = bsStartPos;
          self->dvbAncDataAvailable = 1;
        }
      }
      break;

    default:
      break;
  }

  return (bitCnt);
}

static int aacDecoder_drcParse(HANDLE_FDK_BITSTREAM bs, CDrcPayload *pDrcBs,
                               UINT payloadPosition) {
  int i, numBands;

  FDKpushBiDirectional(bs, (INT)FDKgetValidBits(bs) - (INT)payloadPosition);

  if (FDKreadBits(bs, 1)) {
    pDrcBs->pceInstanceTag = FDKreadBits(bs, 4);
    FDKreadBits(bs, 4);
  } else {
    pDrcBs->pceInstanceTag = -1;
  }

  if (FDKreadBits(bs, 1)) {
    parseExcludedChannels(&pDrcBs->excludedChnsMask, bs);
  } else {
    pDrcBs->excludedChnsMask = 0;
  }

  numBands = 1;
  if (FDKreadBits(bs, 1))
  {
    numBands += FDKreadBits(bs, 4);
    pDrcBs->channelData.drcInterpolationScheme =
        FDKreadBits(bs, 4);
    for (i = 0; i < numBands; i++) {
      pDrcBs->channelData.bandTop[i] = FDKreadBits(bs, 8);
    }
  } else {
    pDrcBs->channelData.bandTop[0] = DRC_BLOCK_LEN_DIV_BAND_MULT -
                                     1;
    ;
  }

  pDrcBs->channelData.numBands = numBands;

  if (FDKreadBits(bs, 1))
  {
    pDrcBs->progRefLevel = FDKreadBits(bs, 7);
    FDKreadBits(bs, 1);
  } else {
    pDrcBs->progRefLevel = -1;
  }

  for (i = 0; i < numBands; i++) {
    pDrcBs->channelData.drcValue[i] = FDKreadBits(bs, 1)
                                      << 7;
    pDrcBs->channelData.drcValue[i] |=
        FDKreadBits(bs, 7) & 0x7F;
  }

  pDrcBs->channelData.drcDataType = MPEG_DRC_EXT_DATA;

  return (1);
}

#define DVB_COMPRESSION_SCALE (8)

static int aacDecoder_drcReadCompression(HANDLE_FDK_BITSTREAM bs,
                                         CDrcPayload *pDrcBs,
                                         UINT payloadPosition) {
  int foundDrcData = 0;
  int dmxLevelsPresent, compressionPresent;

  FDKpushBiDirectional(bs, (INT)FDKgetValidBits(bs) - (INT)payloadPosition);

  if (FDKgetValidBits(bs) < 24) {
    return 0;
  }

  if (FDKreadBits(bs, 8) != DVB_ANC_DATA_SYNC_BYTE) {
    return 0;
  }

  if (FDKreadBits(bs, 2) != 3) {
    return 0;
  }
  FDKreadBits(bs, 2);
  pDrcBs->presMode = FDKreadBits(bs, 2);
  FDKreadBits(bs, 1);
  if (FDKreadBits(bs, 1) != 0) {
    return 0;
  }

  if (FDKreadBits(bs, 3) != 0) {
    return 0;
  }
  dmxLevelsPresent = FDKreadBits(bs, 1);
  FDKreadBits(bs, 1);
  compressionPresent =
      FDKreadBits(bs, 1);
  FDKreadBits(bs, 1);
  FDKreadBits(bs, 1);

  if (dmxLevelsPresent) {
    FDKreadBits(bs, 8);
  }

  if (compressionPresent) {
    UCHAR compressionOn, compressionValue;

    if (FDKreadBits(bs, 7) != 0) {
      return 0;
    }
    compressionOn = (UCHAR)FDKreadBits(bs, 1);
    compressionValue = (UCHAR)FDKreadBits(bs, 8);

    if (compressionOn) {
      pDrcBs->channelData.numBands = 1;
      pDrcBs->channelData.drcValue[0] =
          compressionValue;
      pDrcBs->channelData.bandTop[0] =
          DRC_BLOCK_LEN_DIV_BAND_MULT -
          1;
      pDrcBs->pceInstanceTag = -1;
      pDrcBs->progRefLevel = -1;
      pDrcBs->channelData.drcDataType =
          DVB_DRC_ANC_DATA;
      foundDrcData = 1;
    }
  }

  return (foundDrcData);
}

static int aacDecoder_drcExtractAndMap(
    HANDLE_AAC_DRC self, HANDLE_FDK_BITSTREAM hBs,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[],
    UCHAR pceInstanceTag,
    UCHAR channelMapping[],
    int validChannels) {
  CDrcPayload threadBs[MAX_DRC_THREADS];
  CDrcPayload *validThreadBs[MAX_DRC_THREADS];
  CDrcParams *pParams;
  UINT backupBsPosition;
  int result = 0;
  int i, thread, validThreads = 0;

  FDK_ASSERT(self != NULL);
  FDK_ASSERT(hBs != NULL);
  FDK_ASSERT(pAacDecoderStaticChannelInfo != NULL);

  pParams = &self->params;

  self->numThreads = 0;
  backupBsPosition = FDKgetValidBits(hBs);

  for (i = 0; i < self->numPayloads && self->numThreads < MAX_DRC_THREADS;
       i++) {
    FDKmemclear(&threadBs[self->numThreads], sizeof(CDrcPayload));
    threadBs[self->numThreads].channelData.bandTop[0] =
        DRC_BLOCK_LEN_DIV_BAND_MULT - 1;

    self->numThreads += aacDecoder_drcParse(hBs, &threadBs[self->numThreads],
                                            self->drcPayloadPosition[i]);
  }
  self->numPayloads = 0;

  if (self->dvbAncDataAvailable &&
      self->numThreads < MAX_DRC_THREADS) {

    FDKmemclear(&threadBs[self->numThreads], sizeof(CDrcPayload));
    threadBs[self->numThreads].channelData.bandTop[0] =
        DRC_BLOCK_LEN_DIV_BAND_MULT - 1;

    self->numThreads += aacDecoder_drcReadCompression(
        hBs, &threadBs[self->numThreads], self->dvbAncDataPosition);
  }
  self->dvbAncDataAvailable = 0;

  FDKpushBiDirectional(hBs, (INT)FDKgetValidBits(hBs) - (INT)backupBsPosition);

  for (thread = 0; thread < self->numThreads; thread++) {
    CDrcPayload *pThreadBs = &threadBs[thread];
    int numExclChns = 0;

    switch ((AACDEC_DRC_PAYLOAD_TYPE)pThreadBs->channelData.drcDataType) {
      default:
        continue;
      case MPEG_DRC_EXT_DATA:
      case DVB_DRC_ANC_DATA:
        break;
    }

    if (pThreadBs->pceInstanceTag >= 0) {
      if (pThreadBs->pceInstanceTag != pceInstanceTag) {
        continue;
      }
    }

    if (pThreadBs->excludedChnsMask > 0) {
      INT exclMask = pThreadBs->excludedChnsMask;
      int ch;
      for (ch = 0; ch < validChannels; ch++) {
        numExclChns += exclMask & 0x1;
        exclMask >>= 1;
      }
    }
    if (numExclChns < validChannels) {
      validThreadBs[validThreads] = pThreadBs;
      validThreads++;
    }
  }

  for (thread = 0; thread < validThreads; thread++) {
    CDrcPayload *pThreadBs = validThreadBs[thread];
    INT exclMask = pThreadBs->excludedChnsMask;
    AACDEC_DRC_PAYLOAD_TYPE drcPayloadType =
        (AACDEC_DRC_PAYLOAD_TYPE)pThreadBs->channelData.drcDataType;
    int ch;

    if (pThreadBs->progRefLevel >= 0) {
      self->progRefLevel = pThreadBs->progRefLevel;
      self->progRefLevelPresent = 1;
      self->prlExpiryCount = 0;
    }

    if (drcPayloadType == DVB_DRC_ANC_DATA) {
      self->presMode = pThreadBs->presMode;
    }

    for (ch = 0; ch < validChannels; ch++) {
      AACDEC_DRC_PAYLOAD_TYPE prvPayloadType = UNKNOWN_PAYLOAD;
      int mapedChannel = channelMapping[ch];

      if ((mapedChannel >= validChannels) ||
          ((exclMask & (1 << mapedChannel)) != 0))
        continue;

      if ((pParams->expiryFrame <= 0) ||
          (pAacDecoderStaticChannelInfo[ch]->drcData.expiryCount <
           pParams->expiryFrame)) {
        prvPayloadType =
            (AACDEC_DRC_PAYLOAD_TYPE)pAacDecoderStaticChannelInfo[ch]
                ->drcData.drcDataType;
      }
      if (((drcPayloadType == MPEG_DRC_EXT_DATA) &&
           (prvPayloadType != DVB_DRC_ANC_DATA)) ||
          ((drcPayloadType == DVB_DRC_ANC_DATA) &&
           (pParams->applyHeavyCompression ==
            ON))) {
        pAacDecoderStaticChannelInfo[ch]->drcData = pThreadBs->channelData;
        result = 1;
      }
    }
  }

  if ((pParams->expiryFrame > 0) &&
      (self->prlExpiryCount++ >
       pParams->expiryFrame)) {
    self->progRefLevelPresent = 0;
    self->progRefLevel = pParams->targetRefLevel;
    self->prlExpiryCount = 0;
  }

  return result;
}

void aacDecoder_drcApply(HANDLE_AAC_DRC self, void *pSbrDec,
                         CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                         CDrcChannelData *pDrcChData, FIXP_DBL *extGain,
                         int ch,
                         int aacFrameSize, int bSbrPresent) {
  int band, bin, numBands;
  int bottom = 0;
  int modifyBins = 0;

  FIXP_DBL max_mantissa;
  INT max_exponent;

  FIXP_DBL norm_mantissa = FL2FXCONST_DBL(0.5f);
  INT norm_exponent = 1;

  FIXP_DBL fact_mantissa[MAX_DRC_BANDS];
  INT fact_exponent[MAX_DRC_BANDS];

  CDrcParams *pParams = &self->params;

  FIXP_DBL *pSpectralCoefficient =
      SPEC_LONG(pAacDecoderChannelInfo->pSpectralCoefficient);
  CIcsInfo *pIcsInfo = &pAacDecoderChannelInfo->icsInfo;
  SHORT *pSpecScale = pAacDecoderChannelInfo->specScale;

  int winSeq = pIcsInfo->WindowSequence;

  if ((pParams->expiryFrame > 0) &&
      (++pDrcChData->expiryCount >
       pParams->expiryFrame)) {
    aacDecoder_drcInitChannelData(pDrcChData);
  }

  if (self->enable != ON) {
    sbrDecoder_drcDisable((HANDLE_SBRDECODER)pSbrDec, ch);
    if (extGain != NULL) {
      INT gainScale = (INT)*extGain;
      if (gainScale >= 0 && gainScale <= DFRACT_BITS) {
        *extGain = scaleValue(norm_mantissa, norm_exponent - gainScale);
      } else {
        FDK_ASSERT(0);
      }
    }
    return;
  }

  numBands = pDrcChData->numBands;

  if (pParams->targetRefLevel >= 0) {
    norm_mantissa =
        fLdPow(FL2FXCONST_DBL(-1.0),
               0,
               (FIXP_DBL)((INT)(FL2FXCONST_DBL(1.0f / 24.0) >> 3) *
                          (INT)(pParams->targetRefLevel - self->progRefLevel)),
               3, &norm_exponent);
  }
  if (extGain != NULL) {
    INT gainScale = (INT)*extGain;
    if (gainScale >= 0 && gainScale <= DFRACT_BITS) {
      *extGain = scaleValue(norm_mantissa, norm_exponent - gainScale);
    } else {
      FDK_ASSERT(0);
    }
  }
  norm_mantissa = FL2FXCONST_DBL(0.5f);
  norm_exponent = 1;

  for (band = 0; band < numBands; band++) {
    UCHAR drcVal = pDrcChData->drcValue[band];

    fact_mantissa[band] = FL2FXCONST_DBL(0.5f);
    fact_exponent[band] = 1;

    if ((pParams->applyHeavyCompression == ON) &&
        ((AACDEC_DRC_PAYLOAD_TYPE)pDrcChData->drcDataType ==
         DVB_DRC_ANC_DATA)) {
      INT compressionFactorVal_e;
      int valX, valY;

      valX = drcVal >> 4;
      valY = drcVal & 0x0F;

      if (drcVal != 0x7F) {
        fact_mantissa[band] = fPowInt(
            FL2FXCONST_DBL(0.95483867181),
            0, valY, &compressionFactorVal_e);

        fact_mantissa[band] =
            fMult(FL2FXCONST_DBL(0.99990790084), fact_mantissa[band]);

        fact_exponent[band] =
            DVB_COMPRESSION_SCALE - valX + compressionFactorVal_e;
      }
    } else if ((AACDEC_DRC_PAYLOAD_TYPE)pDrcChData->drcDataType ==
               MPEG_DRC_EXT_DATA) {

      if ((drcVal & 0x7F) > 0) {
        FIXP_DBL tParamVal = (drcVal & 0x80) ? -pParams->cut : pParams->boost;

        fact_mantissa[band] = f2Pow(
            (FIXP_DBL)((INT)fMult(FL2FXCONST_DBL(1.0f / 192.0f), tParamVal) *
                       (drcVal & 0x7F)),
            3 + DRC_PARAM_SCALE, &fact_exponent[band]);
      }
    }

    fact_mantissa[band] = fMult(fact_mantissa[band], norm_mantissa);
    fact_exponent[band] += norm_exponent;
  }

  {
    int res;

    max_mantissa = FL2FXCONST_DBL(0.0f);
    max_exponent = 0;
    for (band = 0; band < numBands; band++) {
      max_mantissa = fixMax(max_mantissa, fact_mantissa[band]);
      max_exponent = fixMax(max_exponent, fact_exponent[band]);
    }

    res = CntLeadingZeros(max_mantissa) - 1;

    if (((pDrcChData->bandTop[fMax(0, numBands - 1)] + 1) << 2) < aacFrameSize)
      res = 0;

    if (res > 0) {
      res = fixMin(res, max_exponent);
      max_exponent -= res;

      for (band = 0; band < numBands; band++) {
        fact_mantissa[band] <<= res;
        fact_exponent[band] -= res;
      }
    }

    for (band = 0; band < numBands; band++) {
      if (fact_exponent[band] < max_exponent) {
        fact_mantissa[band] >>= max_exponent - fact_exponent[band];
      }
      if (fact_mantissa[band] != FL2FXCONST_DBL(0.5f)) {
        modifyBins = 1;
      }
    }
    if (max_exponent != 1) {
      modifyBins = 1;
    }
  }

  if (!bSbrPresent) {
    bottom = 0;

    if (!modifyBins) {
      max_exponent -= 1;
    } else {
      for (band = 0; band < numBands; band++) {
        int top = fixMin((int)((pDrcChData->bandTop[band] + 1) << 2),
                         aacFrameSize);

        for (bin = bottom; bin < top; bin++) {
          pSpectralCoefficient[bin] =
              fMult(pSpectralCoefficient[bin], fact_mantissa[band]);
        }

        bottom = top;
      }
    }

    if (max_exponent > 0) {
      for (bin = bottom; bin < aacFrameSize; bin += 1) {
        pSpectralCoefficient[bin] >>= max_exponent;
      }
    }

    pSpecScale[0] += max_exponent;

    if (winSeq == BLOCK_SHORT) {
      int win;
      for (win = 1; win < 8; win++) {
        pSpecScale[win] += max_exponent;
      }
    }
  } else {
    HANDLE_SBRDECODER hSbrDecoder = (HANDLE_SBRDECODER)pSbrDec;
    numBands = pDrcChData->numBands;

    sbrDecoder_drcFeedChannel(hSbrDecoder, ch, numBands, fact_mantissa,
                              max_exponent, pDrcChData->drcInterpolationScheme,
                              winSeq, pDrcChData->bandTop);
  }

  return;
}

static void aacDecoder_drcParameterHandling(HANDLE_AAC_DRC self,
                                            INT aacNumChannels,
                                            SCHAR prevDrcProgRefLevel,
                                            SCHAR prevDrcPresMode) {
  int isDownmix, isMonoDownmix, isStereoDownmix;
  int dDmx, dHr;
  AACDEC_DRC_PARAMETER_HANDLING drcParameterHandling;
  CDrcParams *p;

  FDK_ASSERT(self != NULL);

  p = &self->params;

  if (self->progRefLevel != prevDrcProgRefLevel) self->update = 1;

  if (self->presMode != prevDrcPresMode) self->update = 1;

  if (self->prevAacNumChannels != aacNumChannels) self->update = 1;

  if (!self->update) {
    return;
  }

  isDownmix = (aacNumChannels > self->numOutChannels);
  isDownmix = (isDownmix && (self->numOutChannels > 0));
  isMonoDownmix = (isDownmix && (self->numOutChannels == 1));
  isStereoDownmix = (isDownmix && (self->numOutChannels == 2));

  if ((self->presMode == 1) || (self->presMode == 2)) {
    drcParameterHandling = (AACDEC_DRC_PARAMETER_HANDLING)self->presMode;
  } else {
    drcParameterHandling = p->defaultPresentationMode;
  }

  p->cut = p->usrCut;
  p->boost = p->usrBoost;
  p->applyHeavyCompression = p->usrApplyHeavyCompression;

  switch (drcParameterHandling) {
    case DISABLED_PARAMETER_HANDLING:
    default:
      break;

    case ENABLED_PARAMETER_HANDLING:
      if (isDownmix) {
        FIXP_DBL dmxTmp;
        int e_log, e_mult;
        dmxTmp = fDivNorm(self->numOutChannels,
                          aacNumChannels);
        dmxTmp = fLog2(dmxTmp, 0, &e_log);
        dmxTmp = fMultNorm(
            dmxTmp, FL2FXCONST_DBL(4.0f * 20.0f * 0.30103f / (float)(1 << 5)),
            &e_mult);
        dDmx = (int)scaleValue(dmxTmp, e_log + e_mult + 5 - (DFRACT_BITS - 1));
      } else {
        dDmx = 0;
      }

      if (p->targetRefLevel >= 0) {
        dHr = p->targetRefLevel + dDmx - self->progRefLevel;
      } else {
        dHr = dDmx;
      }

      if (dHr < 0) {
        int eHr = fixMin(p->encoderTargetLevel - self->progRefLevel, 0);
        if (eHr <
            dHr) {
          FIXP_DBL calcFactor_norm;
          INT calcFactor;
          calcFactor_norm =
              fDivNorm(-dHr, -eHr);
          calcFactor_norm = calcFactor_norm >> DRC_PARAM_SCALE;
          calcFactor = convert_drcParam(
              calcFactor_norm);
          calcFactor_norm = (FIXP_DBL)(
              (INT)(DRC_PARAM_QUANT_STEP >> DRC_PARAM_SCALE) * calcFactor);
          p->cut = (calcFactor_norm > p->cut)
                       ? calcFactor_norm
                       : p->cut;
        } else {
          p->cut = DRC_SCALING_MAX;
          if ((dHr - eHr) <=
              -4 * DRC_HEAVY_THRESHOLD_DB) {
            p->applyHeavyCompression = ON;
          }
        }
      } else {
      }
      break;

    case DRC_PRESENTATION_MODE_1:
      if ((p->targetRefLevel >= 0) &&
          (p->targetRefLevel <
           124)) {
        p->applyHeavyCompression = ON;
      } else {
        if (isMonoDownmix || isStereoDownmix) {
          p->cut = DRC_SCALING_MAX;
        }
      }
      break;

    case DRC_PRESENTATION_MODE_2:
      if ((p->targetRefLevel >= 0) &&
          (p->targetRefLevel <
           124)) {
        if (isMonoDownmix) {
          p->applyHeavyCompression = ON;
        } else {
          p->applyHeavyCompression = OFF;
          p->cut = DRC_SCALING_MAX;
        }
      } else {
        p->applyHeavyCompression = OFF;
        if (isMonoDownmix || isStereoDownmix) {
          p->cut = DRC_SCALING_MAX;
        }
      }
      break;
  }

  if (p->applyHeavyCompression == ON) {
    p->boost = DRC_SCALING_MAX;
    p->cut = DRC_SCALING_MAX;
  }

  self->enable = ((p->boost > (FIXP_DBL)0) || (p->cut > (FIXP_DBL)0) ||
                  (p->applyHeavyCompression == ON) || (p->targetRefLevel >= 0));

  self->prevAacNumChannels = aacNumChannels;
  self->update = 0;
}

int aacDecoder_drcProlog(
    HANDLE_AAC_DRC self, HANDLE_FDK_BITSTREAM hBs,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[],
    UCHAR pceInstanceTag,
    UCHAR channelMapping[],
    int validChannels) {
  int result = 0;

  if (self == NULL) {
    return -1;
  }

  if (!self->params.bsDelayEnable) {
    INT prevPRL, prevPM = 0;
    prevPRL = self->progRefLevel;
    prevPM = self->presMode;

    result = aacDecoder_drcExtractAndMap(
        self, hBs, pAacDecoderStaticChannelInfo, pceInstanceTag, channelMapping,
        validChannels);

    if (result < 0) {
      return result;
    }

    aacDecoder_drcParameterHandling(self, validChannels, prevPRL, prevPM);
  }

  return result;
}

int aacDecoder_drcEpilog(
    HANDLE_AAC_DRC self, HANDLE_FDK_BITSTREAM hBs,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[],
    UCHAR pceInstanceTag,
    UCHAR channelMapping[],
    int validChannels) {
  int result = 0;

  if (self == NULL) {
    return -1;
  }

  if (self->params.bsDelayEnable) {
    INT prevPRL, prevPM = 0;
    prevPRL = self->progRefLevel;
    prevPM = self->presMode;

    result = aacDecoder_drcExtractAndMap(
        self, hBs, pAacDecoderStaticChannelInfo, pceInstanceTag, channelMapping,
        validChannels);

    if (result < 0) {
      return result;
    }

    aacDecoder_drcParameterHandling(self, validChannels, prevPRL, prevPM);
  }

  return result;
}

void aacDecoder_drcGetInfo(HANDLE_AAC_DRC self, SCHAR *pPresMode,
                           SCHAR *pProgRefLevel) {
  if (self != NULL) {
    if (pPresMode != NULL) {
      *pPresMode = self->presMode;
    }
    if (pProgRefLevel != NULL) {
      if (self->progRefLevelPresent) {
        *pProgRefLevel = self->progRefLevel;
      } else {
        *pProgRefLevel = -1;
      }
    }
  }
}

INT applyDrcLevelNormalization(HANDLE_AAC_DRC hDrcInfo, PCM_DEC *samplesIn,
                               FIXP_DBL *pGain, FIXP_DBL *pGainPerSample,
                               const INT gain_scale, const UINT gain_delay,
                               const UINT nSamples, const UINT channels,
                               const UINT stride, const UINT limiterEnabled) {
  UINT i;
  INT additionalGain_scaling;
  FIXP_DBL additionalGain;

  FDK_ASSERT(gain_delay <= nSamples);

  FIXP_DBL additionalGainSmoothState = hDrcInfo->additionalGainFilterState;
  FIXP_DBL additionalGainSmoothState1 = hDrcInfo->additionalGainFilterState1;

  if (!gain_delay) {
    additionalGain = pGain[0];

    additionalGain_scaling =
        fMin(gain_scale, CntLeadingZeros(additionalGain) - 1);
    additionalGain = scaleValue(additionalGain, additionalGain_scaling);

    additionalGain_scaling -= gain_scale;

    if (additionalGain_scaling) {
      scaleValuesSaturate(samplesIn, channels * nSamples,
                          -additionalGain_scaling);
    }

    if (limiterEnabled) {
      FDK_ASSERT(pGainPerSample != NULL);

      for (i = 0; i < nSamples; i++) {
        pGainPerSample[i] = additionalGain;
      }
    } else {
      for (i = 0; i < channels * nSamples; i++) {
        samplesIn[i] = FIXP_DBL2PCM_DEC(fMult(samplesIn[i], additionalGain));
      }
    }
  } else {
    UINT inc;
    FIXP_DBL additionalGainUnfiltered;

    inc = (stride == 1) ? channels : 1;

    for (i = 0; i < nSamples; i++) {
      if (i < gain_delay) {
        additionalGainUnfiltered = hDrcInfo->additionalGainPrev;
      } else {
        additionalGainUnfiltered = pGain[0];
      }

      static const FIXP_SGL b[] = {FL2FXCONST_SGL(0.015466 * 2.0),
                                   FL2FXCONST_SGL(0.015466 * 2.0)};
      static const FIXP_SGL a[] = {(FIXP_SGL)MAXVAL_SGL,
                                   FL2FXCONST_SGL(-0.96907)};

      additionalGain = -fMult(additionalGainSmoothState, a[1]) +
                       fMultDiv2(additionalGainUnfiltered, b[0]) +
                       fMultDiv2(additionalGainSmoothState1, b[1]);
      additionalGainSmoothState1 = additionalGainUnfiltered;
      additionalGainSmoothState = additionalGain;

      additionalGain_scaling =
          fMin(gain_scale, CntLeadingZeros(additionalGain) - 1);
      additionalGain = scaleValue(additionalGain, additionalGain_scaling);

      additionalGain_scaling -= gain_scale;

      if (limiterEnabled) {
        FDK_ASSERT(stride == 1);
        FDK_ASSERT(pGainPerSample != NULL);

        if (additionalGain_scaling) {
          scaleValuesSaturate(samplesIn, channels, -additionalGain_scaling);
        }

        pGainPerSample[i] = additionalGain;
      } else {
        if (additionalGain_scaling) {
          for (UINT k = 0; k < channels; k++) {
            scaleValuesSaturate(&samplesIn[k * stride], 1,
                                -additionalGain_scaling);
          }
        }

        for (UINT k = 0; k < channels; k++) {
          samplesIn[k * stride] =
              FIXP_DBL2PCM_DEC(fMult(samplesIn[k * stride], additionalGain));
        }
      }

      samplesIn += inc;
    }
  }

  hDrcInfo->additionalGainPrev = pGain[0];
  hDrcInfo->additionalGainFilterState = additionalGainSmoothState;
  hDrcInfo->additionalGainFilterState1 = additionalGainSmoothState1;

  return (AACDEC_DRC_GAIN_SCALING);
}
