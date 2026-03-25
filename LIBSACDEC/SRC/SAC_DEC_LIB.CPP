/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Decoder Library Interface

*******************************************************************************/

#include "sac_dec_lib.h"
#include "sac_dec_interface.h"
#include "sac_dec.h"
#include "sac_bitdec.h"
#include "FDK_matrixCalloc.h"

#define MPS_DATA_BUFFER_SIZE (2048)

typedef enum {
  MPEGS_ANCTYPE_FRAME = 0,
  MPEGS_ANCTYPE_HEADER_AND_FRAME = 1,
  MPEGS_ANCTYPE_RESERVED_1 = 2,
  MPEGS_ANCTYPE_RESERVED_2 = 3
} MPEGS_ANCTYPE;

typedef enum {
  MPEGS_CONTINUE = 0,
  MPEGS_STOP = 1,
  MPEGS_START = 2,
  MPEGS_START_STOP =
      3
} MPEGS_ANCSTARTSTOP;

typedef enum {
  MPEGS_SYNC_LOST =
      0,
  MPEGS_SYNC_FOUND = 1,
  MPEGS_SYNC_COMPLETE = 2
} MPEGS_SYNCSTATE;

typedef enum {
  MPEGS_OPMODE_EMM = 0,
  MPEGS_OPMODE_MPS_PAYLOAD = 1,
  MPEGS_OPMODE_NO_MPS_PAYLOAD = 2
} MPEGS_OPMODE;

typedef enum {
  MPEGS_INIT_OK = 0x00000000,
  MPEGS_INIT_ENFORCE_REINIT =
      0x00000001,

  MPEGS_INIT_CHANGE_OUTPUT_MODE =
      0x00000010,
  MPEGS_INIT_CHANGE_PARTIALLY_COMPLEX =
      0x00000020,
  MPEGS_INIT_CHANGE_TIME_FREQ_INTERFACE =
      0x00000040,
  MPEGS_INIT_CHANGE_HEADER = 0x00000080,

  MPEGS_INIT_ERROR_PAYLOAD =
      0x00000100,

  MPEGS_INIT_BS_INTERRUPTION =
      0x00001000,
  MPEGS_INIT_CLEAR_HISTORY =
      0x00002000,

  MPEGS_INIT_CHANGE_CONCEAL_PARAMS = 0x00100000,

  MPEGS_INIT_CHANGE_BYPASS_MODE =
      0x01000000,

  MPEGS_INIT_ERROR_ANC_TYPE = 0x10000000,
  MPEGS_INIT_ERROR_ANC_STARTSTOP =
      0x20000000
} MPEGS_INIT_FLAGS;

struct MpegSurroundDecoder {
  HANDLE_FDK_QMF_DOMAIN pQmfDomain;
  UCHAR mpsData[MPS_DATA_BUFFER_SIZE];
  INT mpsDataBits;
  SPATIAL_SPECIFIC_CONFIG spatialSpecificConfig[1];
  spatialDec *pSpatialDec;
  SPATIAL_SPECIFIC_CONFIG
  spatialSpecificConfigBackup;

  UCHAR mpegSurroundDecoderLevel;
  UCHAR mpegSurroundSscIsGlobalCfg;
  UCHAR mpegSurroundUseTimeInterface;

  SPATIAL_BS_FRAME
  bsFrames[1];
  BS_LL_STATE llState;
  UCHAR bsFrameParse;
  UCHAR bsFrameDecode;
  UCHAR bsFrameDelay;

  SPATIALDEC_PARAM mpegSurroundUserParams;

  SPATIAL_DEC_UPMIX_TYPE upmixType;
  int initFlags[1];
  MPEGS_ANCSTARTSTOP ancStartStopPrev;
  MPEGS_SYNCSTATE fOnSync[1];

  SPATIAL_DEC_CONFIG decConfig;
};

SACDEC_ERROR
static sscCheckOutOfBand(const SPATIAL_SPECIFIC_CONFIG *pSsc,
                         const INT coreCodec, const INT sampleRate,
                         const INT frameSize);

static SACDEC_ERROR sscParseCheck(const SPATIAL_SPECIFIC_CONFIG *pSsc);

static int mpegSurroundDecoder_GetNrOfQmfBands(
    const SPATIAL_SPECIFIC_CONFIG *pSsc, UINT sampleRate) {
  UINT samplingFrequency = sampleRate;
  int qmfBands = 64;

  if (pSsc != NULL) {
    switch (pSsc->coreCodec) {
      case AOT_USAC:
        if ((pSsc->stereoConfigIndex == 3)) {
          static const UCHAR mapIdx2QmfBands[3] = {24, 32, 16};
          FDK_ASSERT((pSsc->coreSbrFrameLengthIndex >= 2) &&
                     (pSsc->coreSbrFrameLengthIndex <= 4));
          qmfBands = mapIdx2QmfBands[pSsc->coreSbrFrameLengthIndex - 2];
        }
        return qmfBands;
      default:
        samplingFrequency = pSsc->samplingFreq;
        break;
    }
  }

  if (samplingFrequency < 27713) {
    qmfBands = 32;
  }
  if (samplingFrequency > 55426) {
    qmfBands = 128;
  }

  return qmfBands;
}

static int mpegSurroundDecoder_CalcInitFlags(SPATIAL_SPECIFIC_CONFIG *pSsc1,
                                             SPATIAL_SPECIFIC_CONFIG *pSsc2,
                                             int upmixTypeFlag,
                                             int binauralQualityFlag,
                                             int partiallyComplexFlag,
                                             int *ctrlFlags) {
  if (pSsc1->coreCodec != pSsc2->coreCodec) {
    *ctrlFlags |= MASK_MPEGS_INIT_ALL_STATES;
    *ctrlFlags |= MASK_MPEGS_INIT_ALL_PARAMS;
  } else {
    if ((partiallyComplexFlag) || (pSsc1->treeConfig != pSsc2->treeConfig) ||
        (pSsc1->samplingFreq != pSsc2->samplingFreq)) {
      *ctrlFlags |= MPEGS_INIT_STATES_ANA_QMF_FILTER;
      *ctrlFlags |= MPEGS_INIT_STATES_ANA_HYB_FILTER;
    }

    if ((upmixTypeFlag) || (partiallyComplexFlag) ||
        (pSsc1->treeConfig != pSsc2->treeConfig) ||
        (pSsc1->samplingFreq != pSsc2->samplingFreq) ||
        (pSsc1->bsFixedGainDMX != pSsc2->bsFixedGainDMX)) {
      *ctrlFlags |= MPEGS_INIT_STATES_SYN_QMF_FILTER;
    }

    if ((upmixTypeFlag) || (partiallyComplexFlag) ||
        (pSsc1->treeConfig != pSsc2->treeConfig) ||
        (pSsc1->samplingFreq != pSsc2->samplingFreq) ||
        (pSsc1->decorrConfig != pSsc2->decorrConfig)) {
      *ctrlFlags |= MPEGS_INIT_STATES_DECORRELATOR;
    }

    if ((upmixTypeFlag) || (binauralQualityFlag) ||
        (pSsc1->treeConfig != pSsc2->treeConfig) ||
        (pSsc1->samplingFreq != pSsc2->samplingFreq))

    {
      *ctrlFlags |= MPEGS_INIT_STATES_M1M2;
    }

    if ((upmixTypeFlag) || (pSsc1->treeConfig != pSsc2->treeConfig) ||
        (pSsc1->tempShapeConfig != pSsc2->tempShapeConfig)) {
      *ctrlFlags |= MPEGS_INIT_STATES_GES;
    }

    if ((upmixTypeFlag) || (binauralQualityFlag) || (partiallyComplexFlag) ||
        (pSsc1->samplingFreq != pSsc2->samplingFreq) ||
        (pSsc1->nTimeSlots != pSsc2->nTimeSlots)) {
      *ctrlFlags |= MPEGS_INIT_STATES_REVERB;
    }

    if (*ctrlFlags & MPEGS_INIT_CONFIG) {
      *ctrlFlags |= MPEGS_INIT_STATES_PARAM;
    }
  }

  return MPS_OK;
}

static void updateMpegSurroundDecoderStatus(
    CMpegSurroundDecoder *pMpegSurroundDecoder, int initFlags,
    MPEGS_SYNCSTATE fOnSync, MPEGS_ANCSTARTSTOP ancStartStopPrev) {
  pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode] |=
      initFlags;
  if ((pMpegSurroundDecoder->mpegSurroundSscIsGlobalCfg != 0) &&
      (pMpegSurroundDecoder->fOnSync[pMpegSurroundDecoder->bsFrameDecode] >=
       MPEGS_SYNC_FOUND) &&
      (fOnSync < MPEGS_SYNC_FOUND)) {
    pMpegSurroundDecoder->fOnSync[pMpegSurroundDecoder->bsFrameDecode] =
        MPEGS_SYNC_FOUND;
  } else {
    pMpegSurroundDecoder->fOnSync[pMpegSurroundDecoder->bsFrameDecode] =
        fOnSync;
  }
  pMpegSurroundDecoder->ancStartStopPrev = ancStartStopPrev;
}

static SACDEC_ERROR mpegSurroundDecoder_Create(
    CMpegSurroundDecoder **pMpegSurroundDecoder, int stereoConfigIndex,
    HANDLE_FDK_QMF_DOMAIN pQmfDomain);

SAC_INSTANCE_AVAIL
mpegSurroundDecoder_IsFullMpegSurroundDecoderInstanceAvailable(
    CMpegSurroundDecoder *pMpegSurroundDecoder) {
  SAC_INSTANCE_AVAIL instanceAvailable = SAC_INSTANCE_NOT_FULL_AVAILABLE;

  if (pMpegSurroundDecoder->pSpatialDec != NULL) {
    instanceAvailable = SAC_INSTANCE_FULL_AVAILABLE;
  }

  return instanceAvailable;
}

SACDEC_ERROR mpegSurroundDecoder_Open(
    CMpegSurroundDecoder **pMpegSurroundDecoder, int stereoConfigIndex,
    HANDLE_FDK_QMF_DOMAIN pQmfDomain) {
  SACDEC_ERROR error;

  error = mpegSurroundDecoder_Create(pMpegSurroundDecoder, stereoConfigIndex,
                                     pQmfDomain);

  return error;
}

static SACDEC_ERROR check_UParam_Build_DecConfig(
    SPATIALDEC_PARAM const *pUserParams, SPATIAL_DEC_CONFIG const *pDecConfig,
    const SPATIAL_SPECIFIC_CONFIG *pSsc, SPATIAL_DEC_UPMIX_TYPE *pUpmixType) {
  int dmxChannels, outChannels, maxNumOutChannels;

  FDK_ASSERT(pUserParams != NULL);
  FDK_ASSERT(pUpmixType != NULL);

  switch (pSsc->treeConfig) {
    case SPATIALDEC_MODE_RSVD7:
      dmxChannels = 1;
      outChannels = 2;
      break;
    default:
      return MPS_UNSUPPORTED_CONFIG;
  }

  switch (pDecConfig->binauralMode) {
    case BINAURAL_NONE:
      break;
    default:
      return MPS_UNSUPPORTED_CONFIG;
  }

  switch (pDecConfig->decoderMode) {
    case EXT_HQ_ONLY:
      break;
    default:
      return MPS_UNSUPPORTED_CONFIG;
  }

  switch (pDecConfig->maxNumOutputChannels) {
    case OUTPUT_CHANNELS_DEFAULT:
      switch (pDecConfig->decoderLevel) {
        case DECODER_LEVEL_0:
          maxNumOutChannels = 2;
          break;
        default:
          return MPS_UNSUPPORTED_CONFIG;
      }
      break;
    case OUTPUT_CHANNELS_2_0:
      maxNumOutChannels = 2;
      break;
    default:
      return MPS_UNSUPPORTED_CONFIG;
  }

  if (pUserParams->blindEnable == 1) {
    return MPS_UNSUPPORTED_CONFIG;
  }
  {
    switch ((SAC_DEC_OUTPUT_MODE)pUserParams->outputMode) {
      case SACDEC_OUT_MODE_NORMAL:
        if (maxNumOutChannels >= outChannels) {
          *pUpmixType = UPMIX_TYPE_NORMAL;
        } else {
          { *pUpmixType = UPMIX_TYPE_BYPASS; }
        }
        break;
      case SACDEC_OUT_MODE_STEREO:
        if (dmxChannels == 1) {
          if (outChannels == 2) {
            *pUpmixType = UPMIX_TYPE_NORMAL;
          }
        } else {
          *pUpmixType = UPMIX_TYPE_BYPASS;
        }
        break;
      case SACDEC_OUT_MODE_6CHANNEL:
        if (outChannels > 6) {
          { *pUpmixType = UPMIX_TYPE_BYPASS; }
        } else {
          *pUpmixType = UPMIX_TYPE_NORMAL;
        }
        break;
      default:
        return MPS_UNSUPPORTED_CONFIG;
    }
  }

  return MPS_OK;
}

static SACDEC_ERROR initMpegSurroundDecoder(
    CMpegSurroundDecoder *pMpegSurroundDecoder) {
  SACDEC_ERROR err;
  int initFlags = MPEGS_INIT_NONE, initFlagsDec;
  int upmixTypeCurr = pMpegSurroundDecoder->upmixType;

  FDK_ASSERT(pMpegSurroundDecoder != NULL);

  SPATIAL_SPECIFIC_CONFIG *const pSSCinput =
      &pMpegSurroundDecoder->spatialSpecificConfigBackup;
  SPATIAL_SPECIFIC_CONFIG *const pSSCtarget =
      &pMpegSurroundDecoder
           ->spatialSpecificConfig[pMpegSurroundDecoder->bsFrameDecode];
  initFlagsDec =
      pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode];

  if (pSSCinput->coreCodec != AOT_USAC) {
    err = sscParseCheck(pSSCinput);
    if (err != MPS_OK) goto bail;
  }

  err = check_UParam_Build_DecConfig(
      &pMpegSurroundDecoder->mpegSurroundUserParams,
      &pMpegSurroundDecoder->decConfig, pSSCinput,
      &pMpegSurroundDecoder->upmixType);
  if (err != MPS_OK) goto bail;

  if (initFlagsDec & MPEGS_INIT_CHANGE_HEADER) {
    initFlags |= MPEGS_INIT_CONFIG;
  }
  if (initFlagsDec & MPEGS_INIT_CLEAR_HISTORY) {
    initFlags |= MASK_MPEGS_INIT_ALL_STATES;
  }
  if (initFlagsDec & MPEGS_INIT_CHANGE_CONCEAL_PARAMS) {
    initFlags |= MPEGS_INIT_PARAMS_ERROR_CONCEALMENT;
  }

  if (initFlagsDec & MPEGS_INIT_ENFORCE_REINIT) {
    initFlags |= MASK_MPEGS_INIT_ALL_STATES;
    initFlags |= MASK_MPEGS_INIT_ALL_PARAMS;
  } else {
    mpegSurroundDecoder_CalcInitFlags(
        pSSCtarget, pSSCinput,
        (upmixTypeCurr !=
         pMpegSurroundDecoder->upmixType),
        0, (initFlagsDec & MPEGS_INIT_CHANGE_PARTIALLY_COMPLEX) ? 1 : 0,
        &initFlags);
  }

  {
    int nrOfQmfBands;
    FDKmemcpy(pSSCtarget, pSSCinput, sizeof(SPATIAL_SPECIFIC_CONFIG));

    nrOfQmfBands = mpegSurroundDecoder_GetNrOfQmfBands(
        pSSCtarget, pSSCtarget->samplingFreq);
    err = FDK_SpatialDecInit(
        pMpegSurroundDecoder->pSpatialDec,
        &pMpegSurroundDecoder->bsFrames[pMpegSurroundDecoder->bsFrameDecode],
        pSSCtarget, nrOfQmfBands, pMpegSurroundDecoder->upmixType,
        &pMpegSurroundDecoder->mpegSurroundUserParams, initFlags);

    if (err != MPS_OK) goto bail;

    if (err == MPS_OK) {
      initFlagsDec = MPEGS_INIT_OK;
      {
        pMpegSurroundDecoder->fOnSync[pMpegSurroundDecoder->bsFrameDecode] =
            MPEGS_SYNC_FOUND;
      }
    }
  }

bail:
  pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode] =
      initFlagsDec;
  return err;
}

SACDEC_ERROR mpegSurroundDecoder_Init(
    CMpegSurroundDecoder *pMpegSurroundDecoder) {
  SACDEC_ERROR err = MPS_OK;

  if (pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode]) {
    err = initMpegSurroundDecoder(pMpegSurroundDecoder);
  }
  return err;
}

static SACDEC_ERROR mpegSurroundDecoder_Create(
    CMpegSurroundDecoder **pMpegSurroundDecoder, int stereoConfigIndex,
    HANDLE_FDK_QMF_DOMAIN pQmfDomain) {
  SACDEC_ERROR err = MPS_OK;
  CMpegSurroundDecoder *sacDec = NULL;
  spatialDec *self = NULL;

  static const SPATIAL_DEC_CONFIG decConfig = {
      (CFG_LEVEL)(0), EXT_HQ_ONLY, OUTPUT_CHANNELS_DEFAULT, BINAURAL_NONE};

  if (*pMpegSurroundDecoder == NULL) {
    FDK_ALLOCATE_MEMORY_1D(*pMpegSurroundDecoder, 1, CMpegSurroundDecoder)

    for (int i = 0; i < 1; i++) {
      err = SpatialDecCreateBsFrame(&(*pMpegSurroundDecoder)->bsFrames[i],
                                    &(*pMpegSurroundDecoder)->llState);
      if (err != MPS_OK) {
        sacDec = *pMpegSurroundDecoder;
        goto bail;
      }
    }
    (*pMpegSurroundDecoder)->pQmfDomain = pQmfDomain;

    (*pMpegSurroundDecoder)->bsFrameDelay = 1;
    (*pMpegSurroundDecoder)->bsFrameParse = 0;
    (*pMpegSurroundDecoder)->bsFrameDecode = 0;

    return err;
  } else {
    sacDec = *pMpegSurroundDecoder;
  }

  if (sacDec->pSpatialDec == NULL) {
    if ((self = FDK_SpatialDecOpen(&decConfig, stereoConfigIndex)) == NULL) {
      err = MPS_OUTOFMEMORY;
      goto bail;
    }
  } else {
    self = sacDec->pSpatialDec;
  }

  self->pQmfDomain = sacDec->pQmfDomain;

  sacDec->pSpatialDec = self;

  sacDec->mpegSurroundUserParams.outputMode = SACDEC_OUT_MODE_NORMAL;
  sacDec->mpegSurroundUserParams.blindEnable = 0;
  sacDec->mpegSurroundUserParams.bypassMode = 0;
  sacDec->mpegSurroundUserParams.concealMethod = 1;
  sacDec->mpegSurroundUserParams.concealNumKeepFrames = 10;
  sacDec->mpegSurroundUserParams.concealFadeOutSlopeLength = 5;
  sacDec->mpegSurroundUserParams.concealFadeInSlopeLength = 5;
  sacDec->mpegSurroundUserParams.concealNumReleaseFrames = 3;
  sacDec->mpegSurroundSscIsGlobalCfg = 0;
  sacDec->mpegSurroundUseTimeInterface = 1;
  sacDec->mpegSurroundDecoderLevel = decConfig.decoderLevel;

  sacDec->upmixType = UPMIX_TYPE_NORMAL;

  updateMpegSurroundDecoderStatus(sacDec, MPEGS_INIT_ENFORCE_REINIT,
                                  MPEGS_SYNC_LOST, MPEGS_STOP);

  *pMpegSurroundDecoder = sacDec;
  sacDec->decConfig = decConfig;

  SpatialDecInitParserContext(sacDec->pSpatialDec);

  return err;

bail:
  if (sacDec != NULL) {
    mpegSurroundDecoder_Close(sacDec);
  }
  *pMpegSurroundDecoder = NULL;
  if (err == MPS_OK) {
    return MPS_OUTOFMEMORY;
  } else {
    return err;
  }
}

SACDEC_ERROR mpegSurroundDecoder_Config(
    CMpegSurroundDecoder *pMpegSurroundDecoder, HANDLE_FDK_BITSTREAM hBs,
    AUDIO_OBJECT_TYPE coreCodec, INT samplingRate, INT frameSize,
    INT numChannels, INT stereoConfigIndex, INT coreSbrFrameLengthIndex,
    INT configBytes, const UCHAR configMode, UCHAR *configChanged) {
  SACDEC_ERROR err = MPS_OK;
  INT nInputChannels;
  SPATIAL_SPECIFIC_CONFIG spatialSpecificConfig;
  SPATIAL_SPECIFIC_CONFIG *pSsc =
      &pMpegSurroundDecoder->spatialSpecificConfigBackup;

  switch (coreCodec) {
    case AOT_DRM_USAC:
    case AOT_USAC:
      if (configMode == AC_CM_DET_CFG_CHANGE) {
        err = SpatialDecParseMps212Config(
            hBs, &spatialSpecificConfig, samplingRate, coreCodec,
            stereoConfigIndex, coreSbrFrameLengthIndex);
        nInputChannels = spatialSpecificConfig.nInputChannels;
        pSsc = &spatialSpecificConfig;
      } else {
        err = SpatialDecParseMps212Config(
            hBs, &pMpegSurroundDecoder->spatialSpecificConfigBackup,
            samplingRate, coreCodec, stereoConfigIndex,
            coreSbrFrameLengthIndex);
        nInputChannels =
            pMpegSurroundDecoder->spatialSpecificConfigBackup.nInputChannels;
      }
      if ((err == MPS_OK) && (numChannels != nInputChannels)) {
        err = MPS_PARSE_ERROR;
      }
      break;
    case AOT_ER_AAC_ELD:
    case AOT_ER_AAC_LD:
      if (configMode == AC_CM_DET_CFG_CHANGE) {
        err = SpatialDecParseSpecificConfig(hBs, &spatialSpecificConfig,
                                            configBytes, coreCodec);
        nInputChannels = spatialSpecificConfig.nInputChannels;
        pSsc = &spatialSpecificConfig;
      } else {
        err = SpatialDecParseSpecificConfig(
            hBs, &pMpegSurroundDecoder->spatialSpecificConfigBackup,
            configBytes, coreCodec);
        nInputChannels =
            pMpegSurroundDecoder->spatialSpecificConfigBackup.nInputChannels;
      }
      if ((err == MPS_OK) && (numChannels > 0) &&
          (numChannels != nInputChannels)) {
        err = MPS_PARSE_ERROR;
      }
      break;
    default:
      err = MPS_UNSUPPORTED_FORMAT;
      break;
  }

  if (err != MPS_OK) {
    goto bail;
  }

  err = sscCheckOutOfBand(pSsc, coreCodec, samplingRate, frameSize);

  if (err != MPS_OK) {
    goto bail;
  }

  if (configMode & AC_CM_DET_CFG_CHANGE) {
    return err;
  }

  if (configMode & AC_CM_ALLOC_MEM) {
    if (*configChanged) {
      err = mpegSurroundDecoder_Open(&pMpegSurroundDecoder, stereoConfigIndex,
                                     NULL);
      if (err) {
        return err;
      }
    }
  }

  {
    SPATIAL_SPECIFIC_CONFIG *sscParse =
        &pMpegSurroundDecoder
             ->spatialSpecificConfig[pMpegSurroundDecoder->bsFrameParse];

    if (FDK_SpatialDecCompareSpatialSpecificConfigHeader(
            &pMpegSurroundDecoder->spatialSpecificConfigBackup, sscParse)) {
      pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameParse] |=
          MPEGS_INIT_CHANGE_HEADER;
      if (pMpegSurroundDecoder->pSpatialDec == NULL) {
        err = MPS_NOTOK;
        goto bail;
      }
      SpatialDecInitParserContext(pMpegSurroundDecoder->pSpatialDec);
      pMpegSurroundDecoder->pSpatialDec->pConfigCurrent =
          &pMpegSurroundDecoder
             ->spatialSpecificConfig[pMpegSurroundDecoder->bsFrameDecode];
    }
  }

  if (err == MPS_OK) {
    pMpegSurroundDecoder->mpegSurroundSscIsGlobalCfg = 1;
  }

bail:
  return err;
}

static MPEGS_OPMODE mpegSurroundOperationMode(
    CMpegSurroundDecoder *pMpegSurroundDecoder, int mpsDataBits) {
  MPEGS_OPMODE mode;

  {
    if ((mpsDataBits > 0) &&
        (pMpegSurroundDecoder->mpegSurroundUserParams.blindEnable == 0)) {
      mode = MPEGS_OPMODE_MPS_PAYLOAD;
    } else {
      mode = MPEGS_OPMODE_NO_MPS_PAYLOAD;
      updateMpegSurroundDecoderStatus(pMpegSurroundDecoder,
                                      MPEGS_INIT_ERROR_PAYLOAD, MPEGS_SYNC_LOST,
                                      MPEGS_STOP);
    }
  }

  return (mode);
}

static SACDEC_ERROR sscParseCheck(const SPATIAL_SPECIFIC_CONFIG *pSsc) {
  if (pSsc->samplingFreq > 96000) return MPS_PARSE_ERROR;
  if (pSsc->samplingFreq < 8000) return MPS_PARSE_ERROR;

  if ((pSsc->treeConfig < 0) || (pSsc->treeConfig > 7)) {
    return MPS_PARSE_ERROR;
  }

  if ((pSsc->quantMode < 0) || (pSsc->quantMode > 2)) {
    return MPS_PARSE_ERROR;
  }

  return MPS_OK;
}

SACDEC_ERROR checkTimeSlots(int frameLength, int qmfBands, int timeSlots) {
  int len;
  int maxFrameLength;

  if (qmfBands == 64) {
    switch (frameLength) {
      case 960:
      case 1920:
        maxFrameLength = 3840;
        break;
      case 1024:
      case 2048:
        maxFrameLength = 4096;
        break;
      case 512:
      case 1152:
        maxFrameLength = 4608;
        break;
      default:
        return MPS_PARSE_ERROR;
    }
  } else if (qmfBands == 32) {
    switch (frameLength) {
      case 960:
      case 1920:
        maxFrameLength = 1920;
        break;
      case 512:
      case 1024:
      case 2048:
        maxFrameLength = 2048;
        break;
      case 1152:
        maxFrameLength = 2304;
        break;
      default:
        return MPS_PARSE_ERROR;
    }
  } else if (qmfBands == 128) {
    switch (frameLength) {
      case 1920:
        maxFrameLength = 7680;
        break;
      case 1024:
        maxFrameLength = 9216;
        break;
      case 2048:
        maxFrameLength = 8192;
        break;
      case 512:
      case 960:
      case 1152:
      default:
        return MPS_PARSE_ERROR;
    }
  } else {
    return MPS_PARSE_ERROR;
  }

  len = frameLength;

  while (len <= maxFrameLength) {
    if (len == timeSlots * qmfBands) {
      return MPS_OK;
    }
    len += frameLength;
  }
  return MPS_PARSE_ERROR;
}

static SACDEC_ERROR sscCheckInBand(SPATIAL_SPECIFIC_CONFIG *pSsc,
                                   int frameLength, int sampleRate) {
  SACDEC_ERROR err = MPS_OK;
  int qmfBands;

  FDK_ASSERT(pSsc != NULL);

  if (sscParseCheck(pSsc) != MPS_OK) {
    err = MPS_PARSE_ERROR;
  }

  if (pSsc->samplingFreq != sampleRate) {
    err = MPS_PARSE_ERROR;
  }

  qmfBands = mpegSurroundDecoder_GetNrOfQmfBands(pSsc, pSsc->samplingFreq);

  if (checkTimeSlots(frameLength, qmfBands, pSsc->nTimeSlots) != MPS_OK) {
    err = MPS_PARSE_ERROR;
  }

  return err;
}

SACDEC_ERROR
mpegSurroundDecoder_ConfigureQmfDomain(
    CMpegSurroundDecoder *pMpegSurroundDecoder,
    SAC_INPUT_CONFIG sac_dec_interface, UINT coreSamplingRate,
    AUDIO_OBJECT_TYPE coreCodec) {
  SACDEC_ERROR err = MPS_OK;
  FDK_QMF_DOMAIN_GC *pGC = NULL;

  if (pMpegSurroundDecoder == NULL) {
    return MPS_INVALID_HANDLE;
  }

  FDK_ASSERT(pMpegSurroundDecoder->pSpatialDec);

  pGC = &pMpegSurroundDecoder->pQmfDomain->globalConf;
  if (pMpegSurroundDecoder->mpegSurroundSscIsGlobalCfg) {
    SPATIAL_SPECIFIC_CONFIG *pSSC =
        &pMpegSurroundDecoder->spatialSpecificConfigBackup;
    if (sac_dec_interface == SAC_INTERFACE_TIME) {
      pGC->nBandsAnalysis_requested = mpegSurroundDecoder_GetNrOfQmfBands(
          pSSC, coreSamplingRate);
      pGC->nBandsSynthesis_requested = pGC->nBandsAnalysis_requested;
      pGC->nInputChannels_requested =
          fMax((UINT)pSSC->nInputChannels, (UINT)pGC->nInputChannels_requested);
    }
    pGC->nOutputChannels_requested =
        fMax((UINT)pSSC->nOutputChannels, (UINT)pGC->nOutputChannels_requested);
  } else {
    if (sac_dec_interface == SAC_INTERFACE_TIME) {
      pGC->nBandsAnalysis_requested = mpegSurroundDecoder_GetNrOfQmfBands(
          NULL, coreSamplingRate);
      pGC->nBandsSynthesis_requested = pGC->nBandsAnalysis_requested;
      pGC->nInputChannels_requested =
          pMpegSurroundDecoder->pSpatialDec->createParams.maxNumInputChannels;
    }
    pGC->nOutputChannels_requested =
        pMpegSurroundDecoder->pSpatialDec->createParams.maxNumOutputChannels;
  }
  pGC->nQmfProcBands_requested = 64;
  pGC->nQmfProcChannels_requested =
      fMin((INT)pGC->nInputChannels_requested,
           pMpegSurroundDecoder->pSpatialDec->createParams.maxNumInputChannels);

  if (coreCodec == AOT_ER_AAC_ELD) {
    pGC->flags_requested |= QMF_FLAG_MPSLDFB;
    pGC->flags_requested &= ~QMF_FLAG_CLDFB;
  }

  return err;
}

SACDEC_ERROR
sscCheckOutOfBand(const SPATIAL_SPECIFIC_CONFIG *pSsc, const INT coreCodec,
                  const INT sampleRate, const INT frameSize) {
  FDK_ASSERT(pSsc != NULL);
  int qmfBands = 0;

  if (sscParseCheck(pSsc) != MPS_OK) {
    return MPS_PARSE_ERROR;
  }

  switch (coreCodec) {
    case AOT_USAC:
    case AOT_DRM_USAC:
      if (pSsc->samplingFreq >= 55426) {
        return MPS_PARSE_ERROR;
      }
      break;
    case AOT_ER_AAC_LD:
    case AOT_ER_AAC_ELD:
      if (pSsc->samplingFreq != sampleRate) {
        return MPS_PARSE_ERROR;
      }

      if (pSsc->samplingFreq > 48000) {
        return MPS_PARSE_ERROR;
      }

      qmfBands = mpegSurroundDecoder_GetNrOfQmfBands(pSsc, pSsc->samplingFreq);
      switch (frameSize) {
        case 480:
          if (!((qmfBands == 32) && (pSsc->nTimeSlots == 15))) {
            return MPS_PARSE_ERROR;
          }
          break;
        case 960:
          if (!((qmfBands == 64) && (pSsc->nTimeSlots == 15))) {
            return MPS_PARSE_ERROR;
          }
          break;
        case 512:
          if (!(((qmfBands == 32) && (pSsc->nTimeSlots == 16)) ||
                ((qmfBands == 64) && (pSsc->nTimeSlots == 8)))) {
            return MPS_PARSE_ERROR;
          }
          break;
        case 1024:
          if (!((qmfBands == 64) && (pSsc->nTimeSlots == 16))) {
            return MPS_PARSE_ERROR;
          }
          break;
        default:
          return MPS_PARSE_ERROR;
      }
      break;
    default:
      return MPS_PARSE_ERROR;
      break;
  }

  return MPS_OK;
}

int mpegSurroundDecoder_ParseNoHeader(
    CMpegSurroundDecoder *pMpegSurroundDecoder, HANDLE_FDK_BITSTREAM hBs,
    int *pMpsDataBits, int fGlobalIndependencyFlag) {
  SACDEC_ERROR err = MPS_OK;
  SPATIAL_SPECIFIC_CONFIG *sscParse;
  int bitsAvail, numSacBits;

  if (pMpegSurroundDecoder == NULL || hBs == NULL) {
    return MPS_INVALID_HANDLE;
  }

  sscParse = &pMpegSurroundDecoder
                  ->spatialSpecificConfig[pMpegSurroundDecoder->bsFrameParse];

  bitsAvail = FDKgetValidBits(hBs);

  if (pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameParse]) {
    FDKmemcpy(sscParse, &pMpegSurroundDecoder->spatialSpecificConfigBackup,
              sizeof(SPATIAL_SPECIFIC_CONFIG));
    pMpegSurroundDecoder->fOnSync[pMpegSurroundDecoder->bsFrameParse] =
        MPEGS_SYNC_FOUND;
  }

  if (bitsAvail <= 0) {
    err = MPS_PARSE_ERROR;
  } else {
    err = SpatialDecParseFrameData(
        pMpegSurroundDecoder->pSpatialDec,
        &pMpegSurroundDecoder->bsFrames[pMpegSurroundDecoder->bsFrameParse],
        hBs, sscParse, (UPMIXTYPE)pMpegSurroundDecoder->upmixType,
        fGlobalIndependencyFlag);
    if (err == MPS_OK) {
      pMpegSurroundDecoder->bsFrames[pMpegSurroundDecoder->bsFrameParse]
          .newBsData = 1;
    }
  }

  numSacBits = bitsAvail - (INT)FDKgetValidBits(hBs);

  if (numSacBits > bitsAvail) {
    pMpegSurroundDecoder->bsFrames[pMpegSurroundDecoder->bsFrameParse]
        .newBsData = 0;
    err = MPS_PARSE_ERROR;
  }

  *pMpsDataBits -= numSacBits;

  return err;
}

static int isValidAncType(CMpegSurroundDecoder *pMpegSurroundDecoder,
                          int ancType) {
  int ret = 1;

  if ((ancType != MPEGS_ANCTYPE_HEADER_AND_FRAME) &&
      (ancType != MPEGS_ANCTYPE_FRAME)) {
    ret = 0;
  }

  if (ret == 0) {
    updateMpegSurroundDecoderStatus(pMpegSurroundDecoder,
                                    MPEGS_INIT_ERROR_PAYLOAD, MPEGS_SYNC_LOST,
                                    MPEGS_STOP);
  }

  return (ret);
}

static int isValidAncStartStop(CMpegSurroundDecoder *pMpegSurroundDecoder,
                               int ancStartStop) {
  int ret = 1;

  switch (ancStartStop) {
    case MPEGS_START:
      if ((pMpegSurroundDecoder->ancStartStopPrev == MPEGS_START) ||
          (pMpegSurroundDecoder->ancStartStopPrev == MPEGS_CONTINUE)) {
        ret = 0;
      }
      break;

    case MPEGS_STOP:
      if ((pMpegSurroundDecoder->ancStartStopPrev == MPEGS_STOP) ||
          (pMpegSurroundDecoder->ancStartStopPrev == MPEGS_START_STOP)) {
        ret = 0;
      }
      break;

    case MPEGS_CONTINUE:
    case MPEGS_START_STOP:
      break;
  }

  if (ret == 0) {
    updateMpegSurroundDecoderStatus(pMpegSurroundDecoder,
                                    MPEGS_INIT_ERROR_PAYLOAD, MPEGS_SYNC_LOST,
                                    MPEGS_STOP);
  } else {
    pMpegSurroundDecoder->ancStartStopPrev = (MPEGS_ANCSTARTSTOP)ancStartStop;
  }

  return (ret);
}

int mpegSurroundDecoder_Parse(CMpegSurroundDecoder *pMpegSurroundDecoder,
                              HANDLE_FDK_BITSTREAM hBs, int *pMpsDataBits,
                              AUDIO_OBJECT_TYPE coreCodec, int sampleRate,
                              int frameSize, int fGlobalIndependencyFlag) {
  SACDEC_ERROR err = MPS_OK;
  SPATIAL_SPECIFIC_CONFIG *sscParse;
  SPATIAL_BS_FRAME *bsFrame;
  HANDLE_FDK_BITSTREAM hMpsBsData = NULL;
  FDK_BITSTREAM mpsBsData;
  int mpsDataBits = *pMpsDataBits;
  int mpsBsBits;
  MPEGS_ANCTYPE ancType;
  MPEGS_ANCSTARTSTOP ancStartStop;

  if (pMpegSurroundDecoder == NULL) {
    return MPS_INVALID_HANDLE;
  }

  FDK_ASSERT(pMpegSurroundDecoder->pSpatialDec);

  mpsBsBits = (INT)FDKgetValidBits(hBs);

  sscParse = &pMpegSurroundDecoder
                  ->spatialSpecificConfig[pMpegSurroundDecoder->bsFrameParse];
  bsFrame = &pMpegSurroundDecoder->bsFrames[pMpegSurroundDecoder->bsFrameParse];

  {
    ancType = (MPEGS_ANCTYPE)FDKreadBits(hBs, 2);
    ancStartStop = (MPEGS_ANCSTARTSTOP)FDKreadBits(hBs, 2);
    mpsDataBits -= 4;

    if (isValidAncType(pMpegSurroundDecoder, ancType)) {
      if (isValidAncStartStop(pMpegSurroundDecoder, ancStartStop)) {
        switch (ancStartStop) {
          case MPEGS_START:
            if (mpsDataBits > MPS_DATA_BUFFER_SIZE * 8) {
              err = MPS_NOTOK;
              goto bail;
            }
            for (int i = 0; i < mpsDataBits / 8; i++) {
              pMpegSurroundDecoder->mpsData[i] = FDKreadBits(hBs, 8);
            }
            pMpegSurroundDecoder->mpsDataBits = mpsDataBits;
            break;

          case MPEGS_CONTINUE:
          case MPEGS_STOP:
            if ((pMpegSurroundDecoder->mpsDataBits + mpsDataBits) >
                MPS_DATA_BUFFER_SIZE * 8) {
              err = MPS_NOTOK;
              goto bail;
            }
            for (int i = 0; i < mpsDataBits / 8; i++) {
              pMpegSurroundDecoder
                  ->mpsData[(pMpegSurroundDecoder->mpsDataBits / 8) + i] =
                  FDKreadBits(hBs, 8);
            }
            pMpegSurroundDecoder->mpsDataBits += mpsDataBits;
            FDKinitBitStream(&mpsBsData, pMpegSurroundDecoder->mpsData,
                             MAX_BUFSIZE_BYTES,
                             pMpegSurroundDecoder->mpsDataBits, BS_READER);
            hMpsBsData = &mpsBsData;
            break;

          case MPEGS_START_STOP:
            pMpegSurroundDecoder->mpsDataBits = mpsDataBits;
            hMpsBsData = hBs;
            break;

          default:
            FDK_ASSERT(0);
        }

        if ((ancStartStop == MPEGS_STOP) ||
            (ancStartStop == MPEGS_START_STOP)) {
          switch (ancType) {
            case MPEGS_ANCTYPE_HEADER_AND_FRAME: {
              int parseResult, bitsRead;
              SPATIAL_SPECIFIC_CONFIG spatialSpecificConfigTmp =
                  pMpegSurroundDecoder->spatialSpecificConfigBackup;

              bitsRead = (INT)FDKgetValidBits(hMpsBsData);

              err = SpatialDecParseSpecificConfigHeader(
                  hMpsBsData,
                  &pMpegSurroundDecoder->spatialSpecificConfigBackup, coreCodec,
                  pMpegSurroundDecoder->upmixType);

              bitsRead = (bitsRead - (INT)FDKgetValidBits(hMpsBsData));
              parseResult = ((err == MPS_OK) ? bitsRead : -bitsRead);

              if (parseResult < 0) {
                parseResult = -parseResult;
                err = MPS_PARSE_ERROR;
              } else if (err == MPS_OK) {
                err = sscCheckInBand(
                    &pMpegSurroundDecoder->spatialSpecificConfigBackup,
                    frameSize, sampleRate);
              }
              if (err != MPS_OK) {
                pMpegSurroundDecoder->spatialSpecificConfigBackup =
                    spatialSpecificConfigTmp;
                break;
              }

              pMpegSurroundDecoder->mpsDataBits -= parseResult;

              if (FDK_SpatialDecCompareSpatialSpecificConfigHeader(
                      &pMpegSurroundDecoder->spatialSpecificConfigBackup,
                      sscParse) == MPS_UNEQUAL_SSC) {
                pMpegSurroundDecoder
                    ->initFlags[pMpegSurroundDecoder->bsFrameParse] |=
                    MPEGS_INIT_CHANGE_HEADER;
                SpatialDecInitParserContext(pMpegSurroundDecoder->pSpatialDec);
                pMpegSurroundDecoder->mpegSurroundSscIsGlobalCfg = 0;
              }
            }
              FDK_FALLTHROUGH;
            case MPEGS_ANCTYPE_FRAME:

              if (pMpegSurroundDecoder
                      ->initFlags[pMpegSurroundDecoder->bsFrameParse] &
                  MPEGS_INIT_ERROR_PAYLOAD) {
                err = MPS_PARSE_ERROR;
                break;
              }

              if (pMpegSurroundDecoder
                      ->initFlags[pMpegSurroundDecoder->bsFrameParse]) {
                FDKmemcpy(sscParse,
                          &pMpegSurroundDecoder->spatialSpecificConfigBackup,
                          sizeof(SPATIAL_SPECIFIC_CONFIG));
                pMpegSurroundDecoder
                    ->fOnSync[pMpegSurroundDecoder->bsFrameParse] =
                    MPEGS_SYNC_FOUND;
              }

              if (pMpegSurroundDecoder
                      ->fOnSync[pMpegSurroundDecoder->bsFrameParse] >=
                  MPEGS_SYNC_FOUND) {
                int nbits = 0, bitsAvail;

                if (err != MPS_OK) {
                  break;
                }

                bitsAvail = FDKgetValidBits(hMpsBsData);

                if (bitsAvail <= 0) {
                  err = MPS_PARSE_ERROR;
                } else {
                  err = SpatialDecParseFrameData(
                      pMpegSurroundDecoder->pSpatialDec, bsFrame, hMpsBsData,
                      sscParse, (UPMIXTYPE)pMpegSurroundDecoder->upmixType,
                      fGlobalIndependencyFlag);
                  if (err == MPS_OK) {
                    bsFrame->newBsData = 1;
                  }
                }

                nbits = bitsAvail - (INT)FDKgetValidBits(hMpsBsData);

                if ((nbits > bitsAvail) ||
                    (nbits > pMpegSurroundDecoder->mpsDataBits) ||
                    (pMpegSurroundDecoder->mpsDataBits > nbits + 7 &&
                     !IS_LOWDELAY(coreCodec))) {
                  bsFrame->newBsData = 0;
                  err = MPS_PARSE_ERROR;
                  break;
                }
                pMpegSurroundDecoder->mpsDataBits -= nbits;
              }
              break;

            default:
              err = MPS_NOTOK;
              break;
          }

          if (err == MPS_OK) {
            pMpegSurroundDecoder->ancStartStopPrev = ancStartStop;
          } else {
            updateMpegSurroundDecoderStatus(pMpegSurroundDecoder,
                                            MPEGS_INIT_ERROR_PAYLOAD,
                                            MPEGS_SYNC_LOST, MPEGS_STOP);
            pMpegSurroundDecoder->mpsDataBits = 0;
          }
        }
      }
    }
  }

bail:

  *pMpsDataBits -= (mpsBsBits - (INT)FDKgetValidBits(hBs));

  return err;
}

int mpegSurroundDecoder_Apply(CMpegSurroundDecoder *pMpegSurroundDecoder,
                              PCM_MPS *input, PCM_MPS *pTimeData,
                              const int timeDataSize, int timeDataFrameSize,
                              int *nChannels, int *frameSize, int sampleRate,
                              AUDIO_OBJECT_TYPE coreCodec,
                              AUDIO_CHANNEL_TYPE channelType[],
                              UCHAR channelIndices[],
                              const FDK_channelMapDescr *const mapDescr,
                              const INT inDataHeadroom, INT *outDataHeadroom) {
  SACDEC_ERROR err = MPS_OK;
  PCM_MPS *pTimeOut = pTimeData;
  PCM_MPS *TDinput = NULL;
  UINT initControlFlags = 0, controlFlags = 0;
  int timeDataRequiredSize = 0;
  int newData;

  if (pMpegSurroundDecoder == NULL) {
    return MPS_INVALID_HANDLE;
  }

  FDK_ASSERT(pMpegSurroundDecoder->pSpatialDec);

  if (!FDK_chMapDescr_isValid(mapDescr)) {
    return MPS_INVALID_HANDLE;
  }

  if ((*nChannels <= 0) || (*nChannels > 2)) {
    return MPS_NOTOK;
  }

  pMpegSurroundDecoder->pSpatialDec->sacInDataHeadroom = inDataHeadroom;
  *outDataHeadroom = (INT)(8);

  pMpegSurroundDecoder->pSpatialDec->pConfigCurrent =
      &pMpegSurroundDecoder
           ->spatialSpecificConfig[pMpegSurroundDecoder->bsFrameDecode];
  newData = pMpegSurroundDecoder->bsFrames[pMpegSurroundDecoder->bsFrameParse]
                .newBsData;

  switch (mpegSurroundOperationMode(pMpegSurroundDecoder, 1000)) {
    case MPEGS_OPMODE_MPS_PAYLOAD:
      if (pMpegSurroundDecoder
              ->initFlags[pMpegSurroundDecoder->bsFrameDecode]) {
        err = initMpegSurroundDecoder(pMpegSurroundDecoder);
      }

      if (err == MPS_OK) {
        if ((pMpegSurroundDecoder
                 ->fOnSync[pMpegSurroundDecoder->bsFrameDecode] !=
             MPEGS_SYNC_COMPLETE) &&
            (pMpegSurroundDecoder->bsFrames[pMpegSurroundDecoder->bsFrameDecode]
                 .bsIndependencyFlag == 1)) {
          pMpegSurroundDecoder->fOnSync[pMpegSurroundDecoder->bsFrameDecode] =
              MPEGS_SYNC_COMPLETE;
        }
      } else {
        pMpegSurroundDecoder->fOnSync[pMpegSurroundDecoder->bsFrameDecode] =
            MPEGS_SYNC_FOUND;
        controlFlags |= MPEGS_CONCEAL;
        err = MPS_OK;
      }
      if (pMpegSurroundDecoder->fOnSync[pMpegSurroundDecoder->bsFrameDecode] !=
          MPEGS_SYNC_COMPLETE) {
        controlFlags |= MPEGS_CONCEAL;
      }
      break;

    case MPEGS_OPMODE_NO_MPS_PAYLOAD:
      controlFlags |= MPEGS_CONCEAL;
      break;

    default:
      err = MPS_NOTOK;
  }

  if (err != MPS_OK) {
    goto bail;
  }

  if (pMpegSurroundDecoder->mpegSurroundUserParams.bypassMode) {
    controlFlags |= MPEGS_BYPASSMODE;
  }

  if (pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode]) {
    int startWithDfltCfg = 0;
    if (pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode] &
        MPEGS_INIT_ENFORCE_REINIT) {
      if (FDK_SpatialDecInitDefaultSpatialSpecificConfig(
              &pMpegSurroundDecoder->spatialSpecificConfigBackup, coreCodec,
              *nChannels, sampleRate,
              *frameSize /
                  mpegSurroundDecoder_GetNrOfQmfBands(NULL, sampleRate),
              pMpegSurroundDecoder->mpegSurroundDecoderLevel,
              pMpegSurroundDecoder->mpegSurroundUserParams.blindEnable)) {
        err = MPS_NOTOK;
        goto bail;
      }

      if (FDK_SpatialDecCompareSpatialSpecificConfigHeader(
              &pMpegSurroundDecoder->spatialSpecificConfigBackup,
              &pMpegSurroundDecoder->spatialSpecificConfig
                   [pMpegSurroundDecoder->bsFrameDecode]) == MPS_UNEQUAL_SSC) {
        pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode] |=
            MPEGS_INIT_CHANGE_HEADER;
        SpatialDecInitParserContext(pMpegSurroundDecoder->pSpatialDec);
      }

      startWithDfltCfg = 1;
    }

    err = initMpegSurroundDecoder(pMpegSurroundDecoder);

    if (startWithDfltCfg) {
      pMpegSurroundDecoder->fOnSync[pMpegSurroundDecoder->bsFrameDecode] =
          MPEGS_SYNC_LOST;
    }

    controlFlags |= MPEGS_CONCEAL;

    if (err != MPS_OK) {
      goto bail;
    }
  }

  initControlFlags = controlFlags;

  if (pMpegSurroundDecoder->pQmfDomain->globalConf.nBandsAnalysis == 0) {
    err = MPS_UNSUPPORTED_FORMAT;
    goto bail;
  }
  timeDataRequiredSize =
      (timeDataFrameSize *
       pMpegSurroundDecoder->pSpatialDec->numOutputChannelsAT *
       pMpegSurroundDecoder->pQmfDomain->globalConf.nBandsSynthesis) /
      pMpegSurroundDecoder->pQmfDomain->globalConf.nBandsAnalysis;
  if (timeDataSize < timeDataRequiredSize) {
    err = MPS_OUTPUT_BUFFER_TOO_SMALL;
    goto bail;
  }

  if ((pMpegSurroundDecoder->pSpatialDec->pConfigCurrent->syntaxFlags &
       SACDEC_SYNTAX_USAC) &&
      (pMpegSurroundDecoder->pSpatialDec->stereoConfigIndex > 1)) {
    FDK_ASSERT(timeDataRequiredSize >= timeDataFrameSize * *nChannels);
    int timeDataFrameSizeOut =
        (timeDataFrameSize *
         pMpegSurroundDecoder->pQmfDomain->globalConf.nBandsSynthesis) /
        pMpegSurroundDecoder->pQmfDomain->globalConf.nBandsAnalysis;
    TDinput = pTimeData + timeDataFrameSizeOut - timeDataFrameSize;
    for (int i = *nChannels - 1; i >= 0; i--) {
      FDKmemmove(pTimeData + (i + 1) * timeDataFrameSizeOut - timeDataFrameSize,
                 pTimeData + timeDataFrameSize * i,
                 sizeof(PCM_MPS) * timeDataFrameSize);
      FDKmemclear(pTimeData + i * timeDataFrameSizeOut,
                  sizeof(PCM_MPS) * (timeDataFrameSizeOut - timeDataFrameSize));
    }
  } else {
    if (pMpegSurroundDecoder->mpegSurroundUseTimeInterface) {
      FDKmemcpy(input, pTimeData,
                sizeof(PCM_MPS) * (*nChannels) * (*frameSize));
      TDinput = input;
    }
  }

  err = SpatialDecApplyFrame(
      pMpegSurroundDecoder->pSpatialDec,
      &pMpegSurroundDecoder->bsFrames[pMpegSurroundDecoder->bsFrameDecode],
      pMpegSurroundDecoder->mpegSurroundUseTimeInterface ? INPUTMODE_TIME
                                                         : INPUTMODE_QMF_SBR,
      TDinput, NULL, NULL, pTimeOut, *frameSize, &controlFlags, *nChannels,
      mapDescr);
  *nChannels = pMpegSurroundDecoder->pSpatialDec->numOutputChannelsAT;

  if (err !=
      MPS_OK) {
    updateMpegSurroundDecoderStatus(pMpegSurroundDecoder,
                                    MPEGS_INIT_ENFORCE_REINIT, MPEGS_SYNC_LOST,
                                    MPEGS_STOP);
    *frameSize =
        0;
  } else {
    if (((controlFlags & MPEGS_CONCEAL) &&
         !(initControlFlags & MPEGS_CONCEAL)) ||
        (pMpegSurroundDecoder->pSpatialDec->errInt !=
         MPS_OK)) {
      updateMpegSurroundDecoderStatus(pMpegSurroundDecoder,
                                      MPEGS_INIT_ERROR_PAYLOAD, MPEGS_SYNC_LOST,
                                      MPEGS_STOP);
    }
  }

  if ((err == MPS_OK) && !(controlFlags & MPEGS_BYPASSMODE) &&
      !(pMpegSurroundDecoder->upmixType == UPMIX_TYPE_BYPASS)) {
    SpatialDecChannelProperties(pMpegSurroundDecoder->pSpatialDec, channelType,
                                channelIndices, mapDescr);
  }

bail:

  if (newData) {
    pMpegSurroundDecoder->pSpatialDec->numParameterSetsPrev =
        pMpegSurroundDecoder->bsFrames[pMpegSurroundDecoder->bsFrameDecode]
            .numParameterSets;
  }

  return (err);
}

SACDEC_ERROR mpegSurroundDecoder_FreeMem(
    CMpegSurroundDecoder *pMpegSurroundDecoder) {
  SACDEC_ERROR err = MPS_OK;

  if (pMpegSurroundDecoder != NULL) {
    FDK_SpatialDecClose(pMpegSurroundDecoder->pSpatialDec);
    pMpegSurroundDecoder->pSpatialDec = NULL;
  }

  return err;
}

void mpegSurroundDecoder_Close(CMpegSurroundDecoder *pMpegSurroundDecoder) {
  if (pMpegSurroundDecoder != NULL) {
    FDK_SpatialDecClose(pMpegSurroundDecoder->pSpatialDec);
    pMpegSurroundDecoder->pSpatialDec = NULL;

    for (int i = 0; i < 1; i++) {
      SpatialDecCloseBsFrame(&pMpegSurroundDecoder->bsFrames[i]);
    }

    FDK_FREE_MEMORY_1D(pMpegSurroundDecoder);
  }
}

#define SACDEC_VL0 2
#define SACDEC_VL1 1
#define SACDEC_VL2 0

int mpegSurroundDecoder_GetLibInfo(LIB_INFO *info) {
  int i;

  if (info == NULL) {
    return -1;
  }

  for (i = 0; i < FDK_MODULE_LAST; i++) {
    if (info[i].module_id == FDK_NONE) break;
  }
  if (i == FDK_MODULE_LAST) return -1;

  info += i;

  info->module_id = FDK_MPSDEC;
  info->build_date = "";
  info->build_time = "";
  info->title = "MPEG Surround Decoder";
  info->version = LIB_VERSION(SACDEC_VL0, SACDEC_VL1, SACDEC_VL2);
  LIB_VERSION_STRING(info);
  info->flags = 0 | CAPF_MPS_LD | CAPF_MPS_USAC | CAPF_MPS_HQ |
                CAPF_MPS_1CH_IN | CAPF_MPS_2CH_OUT;

  return 0;
}

SACDEC_ERROR mpegSurroundDecoder_SetParam(
    CMpegSurroundDecoder *pMpegSurroundDecoder, const SACDEC_PARAM param,
    const INT value) {
  SACDEC_ERROR err = MPS_OK;
  SPATIALDEC_PARAM *pUserParams = NULL;

  if (pMpegSurroundDecoder != NULL) {
    pUserParams = &pMpegSurroundDecoder->mpegSurroundUserParams;
  } else {
    err = MPS_INVALID_HANDLE;
  }

  switch (param) {
    case SACDEC_OUTPUT_MODE:
      switch ((SAC_DEC_OUTPUT_MODE)value) {
        case SACDEC_OUT_MODE_NORMAL:
        case SACDEC_OUT_MODE_STEREO:
          break;
        default:
          err = MPS_INVALID_PARAMETER;
      }
      if (err == MPS_OK) {
        if (0) {
          err = MPS_INVALID_PARAMETER;
        } else if (pUserParams->outputMode != (UCHAR)value) {
          pUserParams->outputMode = (UCHAR)value;
          pMpegSurroundDecoder
              ->initFlags[pMpegSurroundDecoder->bsFrameDecode] |=
              MPEGS_INIT_CHANGE_OUTPUT_MODE;
        }
      }
      break;

    case SACDEC_INTERFACE:
      if (value < 0 || value > 1) {
        err = MPS_INVALID_PARAMETER;
      }
      if (err != MPS_OK) {
        goto bail;
      }
      if (pMpegSurroundDecoder->mpegSurroundUseTimeInterface != (UCHAR)value) {
        pMpegSurroundDecoder->mpegSurroundUseTimeInterface = (UCHAR)value;
        pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode] |=
            MPEGS_INIT_CHANGE_TIME_FREQ_INTERFACE;
      }
      break;

    case SACDEC_BS_INTERRUPTION:
      if ((err == MPS_OK) && (value != 0)) {
        updateMpegSurroundDecoderStatus(pMpegSurroundDecoder,
                                        MPEGS_INIT_BS_INTERRUPTION,
                                        MPEGS_SYNC_LOST, MPEGS_STOP);
      }
      break;

    case SACDEC_CLEAR_HISTORY:
      if ((err == MPS_OK) && (value != 0)) {
        updateMpegSurroundDecoderStatus(pMpegSurroundDecoder,
                                        MPEGS_INIT_CLEAR_HISTORY,
                                        MPEGS_SYNC_LOST, MPEGS_STOP);
      }
      break;

    case SACDEC_CONCEAL_NUM_KEEP_FRAMES:
      if (value < 0) {
        err = MPS_INVALID_PARAMETER;
      }
      if (err != MPS_OK) {
        goto bail;
      }
      if (pUserParams->concealNumKeepFrames != (UINT)value) {
        pUserParams->concealNumKeepFrames = (UINT)value;
        pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode] |=
            MPEGS_INIT_CHANGE_CONCEAL_PARAMS;
      }
      break;

    case SACDEC_CONCEAL_FADE_OUT_SLOPE_LENGTH:
      if (value < 0) {
        err = MPS_INVALID_PARAMETER;
      }
      if (err != MPS_OK) {
        goto bail;
      }
      if (pUserParams->concealFadeOutSlopeLength != (UINT)value) {
        pUserParams->concealFadeOutSlopeLength = (UINT)value;
        pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode] |=
            MPEGS_INIT_CHANGE_CONCEAL_PARAMS;
      }
      break;

    case SACDEC_CONCEAL_FADE_IN_SLOPE_LENGTH:
      if (value < 0) {
        err = MPS_INVALID_PARAMETER;
      }
      if (err != MPS_OK) {
        goto bail;
      }
      if (pUserParams->concealFadeInSlopeLength != (UINT)value) {
        pUserParams->concealFadeInSlopeLength = (UINT)value;
        pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode] |=
            MPEGS_INIT_CHANGE_CONCEAL_PARAMS;
      }
      break;

    case SACDEC_CONCEAL_NUM_RELEASE_FRAMES:
      if (value < 0) {
        err = MPS_INVALID_PARAMETER;
      }
      if (err != MPS_OK) {
        goto bail;
      }
      if (pUserParams->concealNumReleaseFrames != (UINT)value) {
        pUserParams->concealNumReleaseFrames = (UINT)value;
        pMpegSurroundDecoder->initFlags[pMpegSurroundDecoder->bsFrameDecode] |=
            MPEGS_INIT_CHANGE_CONCEAL_PARAMS;
      }
      break;

    default:
      err = MPS_INVALID_PARAMETER;
      break;
  }

bail:
  return err;
}

SACDEC_ERROR mpegSurroundDecoder_IsPseudoLR(
    CMpegSurroundDecoder *pMpegSurroundDecoder, int *bsPseudoLr) {
  if (pMpegSurroundDecoder != NULL) {
    const SPATIAL_SPECIFIC_CONFIG *sscDecode =
        &pMpegSurroundDecoder
             ->spatialSpecificConfig[pMpegSurroundDecoder->bsFrameDecode];
    *bsPseudoLr = (int)sscDecode->bsPseudoLr;
    return MPS_OK;
  } else
    return MPS_INVALID_HANDLE;
}

UINT mpegSurroundDecoder_GetDelay(const CMpegSurroundDecoder *self) {
  INT outputDelay = 0;

  if (self != NULL) {
    const SPATIAL_SPECIFIC_CONFIG *sscDecode =
        &self->spatialSpecificConfig[self->bsFrameDecode];
    AUDIO_OBJECT_TYPE coreCodec = sscDecode->coreCodec;

    if (coreCodec > AOT_NULL_OBJECT) {
      if (IS_LOWDELAY(coreCodec)) {
        outputDelay += 256;
      } else if (!IS_USAC(coreCodec)) {
        outputDelay += 320 + 257;
        if (self->mpegSurroundUseTimeInterface) {
          outputDelay += 320 + 384;
        }
      }
    }
  }

  return (outputDelay);
}
