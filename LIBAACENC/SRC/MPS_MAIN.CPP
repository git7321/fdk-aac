/**************************** AAC encoder library ******************************

   Author(s):   Markus Lohwasser

   Description: Mpeg Surround library interface functions

*******************************************************************************/

#include "mps_main.h"
#include "sacenc_lib.h"

struct MPS_ENCODER {
  HANDLE_MP4SPACE_ENCODER hSacEncoder;

  AUDIO_OBJECT_TYPE audioObjectType;

  FDK_bufDescr inBufDesc;
  FDK_bufDescr outBufDesc;
  SACENC_InArgs inargs;
  SACENC_OutArgs outargs;

  void *pInBuffer[1];
  UINT pInBufferSize[1];
  UINT pInBufferElSize[1];
  UINT pInBufferType[1];

  void *pOutBuffer[2];
  UINT pOutBufferSize[2];
  UINT pOutBufferElSize[2];
  UINT pOutBufferType[2];

  UCHAR sacOutBuffer[1024];
};

struct MPS_CONFIG_TAB {
  AUDIO_OBJECT_TYPE audio_object_type;
  CHANNEL_MODE channel_mode;
  ULONG sbr_ratio;
  ULONG sampling_rate;
  ULONG bitrate_min;
  ULONG bitrate_max;
};

static const MPS_CONFIG_TAB mpsConfigTab[] = {
    {AOT_ER_AAC_ELD, MODE_212, 0, 16000, 16000, 39999},
    {AOT_ER_AAC_ELD, MODE_212, 0, 22050, 16000, 49999},
    {AOT_ER_AAC_ELD, MODE_212, 0, 24000, 16000, 61999},
    {AOT_ER_AAC_ELD, MODE_212, 0, 32000, 20000, 84999},
    {AOT_ER_AAC_ELD, MODE_212, 0, 44100, 50000, 192000},
    {AOT_ER_AAC_ELD, MODE_212, 0, 48000, 62000, 192000},

    {AOT_ER_AAC_ELD, MODE_212, 1, 16000, 18000, 31999},
    {AOT_ER_AAC_ELD, MODE_212, 1, 22050, 18000, 31999},
    {AOT_ER_AAC_ELD, MODE_212, 1, 24000, 20000, 64000},

    {AOT_ER_AAC_ELD, MODE_212, 2, 32000, 18000, 64000},
    {AOT_ER_AAC_ELD, MODE_212, 2, 44100, 21000, 64000},
    {AOT_ER_AAC_ELD, MODE_212, 2, 48000, 26000, 64000}

};

static INT FDK_MpegsEnc_WriteFrameHeader(HANDLE_MPS_ENCODER hMpsEnc,
                                         UCHAR *const pOutputBuffer,
                                         const int outputBufferSize);

MPS_ENCODER_ERROR FDK_MpegsEnc_Open(HANDLE_MPS_ENCODER *phMpsEnc) {
  MPS_ENCODER_ERROR error = MPS_ENCODER_OK;
  HANDLE_MPS_ENCODER hMpsEnc = NULL;

  if (phMpsEnc == NULL) {
    error = MPS_ENCODER_INVALID_HANDLE;
    goto bail;
  }

  if (NULL ==
      (hMpsEnc = (HANDLE_MPS_ENCODER)FDKcalloc(1, sizeof(MPS_ENCODER)))) {
    error = MPS_ENCODER_MEMORY_ERROR;
    goto bail;
  }
  FDKmemclear(hMpsEnc, sizeof(MPS_ENCODER));

  if (SACENC_OK != FDK_sacenc_open(&hMpsEnc->hSacEncoder)) {
    error = MPS_ENCODER_MEMORY_ERROR;
    goto bail;
  }

  *phMpsEnc = hMpsEnc;

bail:
  if (error != MPS_ENCODER_OK) {
    FDK_MpegsEnc_Close(&hMpsEnc);
  }
  return error;
}

MPS_ENCODER_ERROR FDK_MpegsEnc_Close(HANDLE_MPS_ENCODER *phMpsEnc) {
  MPS_ENCODER_ERROR error = MPS_ENCODER_OK;

  if (phMpsEnc == NULL) {
    error = MPS_ENCODER_INVALID_HANDLE;
    goto bail;
  }

  if (*phMpsEnc != NULL) {
    FDK_sacenc_close(&(*phMpsEnc)->hSacEncoder);
    FDKfree(*phMpsEnc);
    *phMpsEnc = NULL;
  }
bail:
  return error;
}

MPS_ENCODER_ERROR FDK_MpegsEnc_Init(HANDLE_MPS_ENCODER hMpsEnc,
                                    const AUDIO_OBJECT_TYPE audioObjectType,
                                    const UINT samplingrate, const UINT bitrate,
                                    const UINT sbrRatio, const UINT framelength,
                                    const UINT inputBufferSizePerChannel,
                                    const UINT coreCoderDelay) {
  MPS_ENCODER_ERROR error = MPS_ENCODER_OK;
  const UINT fs_low = 27713;
  const UINT fs_high = 55426;
  UINT nTimeSlots = 0, nQmfBandsLd = 0;

  if (hMpsEnc == NULL) {
    error = MPS_ENCODER_INVALID_HANDLE;
    goto bail;
  }

  switch (sbrRatio) {
    case 1:
      if (!(samplingrate < fs_low)) {
        error = MPS_ENCODER_INIT_ERROR;
        goto bail;
      }
      break;
    case 2:
      if (!((samplingrate >= fs_low) && (samplingrate < fs_high))) {
        error = MPS_ENCODER_INIT_ERROR;
        goto bail;
      }
      break;
    case 0:
    default:;
  }

  nQmfBandsLd =
      (samplingrate < fs_low) ? 5 : ((samplingrate > fs_high) ? 7 : 6);
  nTimeSlots = framelength >> nQmfBandsLd;

  if (framelength != (nTimeSlots << nQmfBandsLd)) {
    error = MPS_ENCODER_INIT_ERROR;
    goto bail;
  }

  if ((INT)bitrate != FDK_MpegsEnc_GetClosestBitRate(audioObjectType, MODE_212,
                                                     samplingrate, sbrRatio,
                                                     bitrate)) {
    error = MPS_ENCODER_INIT_ERROR;
    goto bail;
  }

  switch (audioObjectType) {
    case AOT_ER_AAC_ELD: {
      const UINT noInterFrameCoding = 0;

      if ((SACENC_OK !=
           FDK_sacenc_setParam(hMpsEnc->hSacEncoder, SACENC_LOWDELAY,
                               (noInterFrameCoding == 1) ? 1 : 2)) ||
          (SACENC_OK != FDK_sacenc_setParam(hMpsEnc->hSacEncoder,
                                            SACENC_ENC_MODE, SACENC_212)) ||
          (SACENC_OK != FDK_sacenc_setParam(hMpsEnc->hSacEncoder,
                                            SACENC_SAMPLERATE, samplingrate)) ||
          (SACENC_OK != FDK_sacenc_setParam(hMpsEnc->hSacEncoder,
                                            SACENC_FRAME_TIME_SLOTS,
                                            nTimeSlots)) ||
          (SACENC_OK != FDK_sacenc_setParam(hMpsEnc->hSacEncoder,
                                            SACENC_PARAM_BANDS,
                                            SACENC_BANDS_15)) ||
          (SACENC_OK !=
           FDK_sacenc_setParam(hMpsEnc->hSacEncoder, SACENC_TIME_DOM_DMX, 2)) ||
          (SACENC_OK !=
           FDK_sacenc_setParam(hMpsEnc->hSacEncoder, SACENC_COARSE_QUANT, 0)) ||
          (SACENC_OK != FDK_sacenc_setParam(hMpsEnc->hSacEncoder,
                                            SACENC_QUANT_MODE,
                                            SACENC_QUANTMODE_FINE)) ||
          (SACENC_OK != FDK_sacenc_setParam(hMpsEnc->hSacEncoder,
                                            SACENC_TIME_ALIGNMENT, 0)) ||
          (SACENC_OK != FDK_sacenc_setParam(hMpsEnc->hSacEncoder,
                                            SACENC_INDEPENDENCY_FACTOR, 20))) {
        error = MPS_ENCODER_INIT_ERROR;
        goto bail;
      }
      break;
    }
    default:
      error = MPS_ENCODER_INIT_ERROR;
      goto bail;
  }

  if (SACENC_OK != FDK_sacenc_init(hMpsEnc->hSacEncoder, coreCoderDelay)) {
    error = MPS_ENCODER_INIT_ERROR;
  }

  hMpsEnc->audioObjectType = audioObjectType;

  hMpsEnc->inBufDesc.ppBase = (void **)&hMpsEnc->pInBuffer;
  hMpsEnc->inBufDesc.pBufSize = hMpsEnc->pInBufferSize;
  hMpsEnc->inBufDesc.pEleSize = hMpsEnc->pInBufferElSize;
  hMpsEnc->inBufDesc.pBufType = hMpsEnc->pInBufferType;
  hMpsEnc->inBufDesc.numBufs = 1;

  hMpsEnc->outBufDesc.ppBase = (void **)&hMpsEnc->pOutBuffer;
  hMpsEnc->outBufDesc.pBufSize = hMpsEnc->pOutBufferSize;
  hMpsEnc->outBufDesc.pEleSize = hMpsEnc->pOutBufferElSize;
  hMpsEnc->outBufDesc.pBufType = hMpsEnc->pOutBufferType;
  hMpsEnc->outBufDesc.numBufs = 2;

  hMpsEnc->pInBuffer[0] = NULL;
  hMpsEnc->pInBufferSize[0] = 0;
  hMpsEnc->pInBufferElSize[0] = sizeof(INT_PCM);
  hMpsEnc->pInBufferType[0] = (FDK_BUF_TYPE_INPUT | FDK_BUF_TYPE_PCM_DATA);

  hMpsEnc->pOutBuffer[0] = NULL;
  hMpsEnc->pOutBufferSize[0] = 0;
  hMpsEnc->pOutBufferElSize[0] = sizeof(INT_PCM);
  hMpsEnc->pOutBufferType[0] = (FDK_BUF_TYPE_OUTPUT | FDK_BUF_TYPE_PCM_DATA);

  hMpsEnc->pOutBuffer[1] = NULL;
  hMpsEnc->pOutBufferSize[1] = 0;
  hMpsEnc->pOutBufferElSize[1] = sizeof(UCHAR);
  hMpsEnc->pOutBufferType[1] = (FDK_BUF_TYPE_OUTPUT | FDK_BUF_TYPE_BS_DATA);

  hMpsEnc->inargs.isInputInterleaved = 0;
  hMpsEnc->inargs.inputBufferSizePerChannel = inputBufferSizePerChannel;

bail:
  return error;
}

MPS_ENCODER_ERROR FDK_MpegsEnc_Process(HANDLE_MPS_ENCODER hMpsEnc,
                                       INT_PCM *const pAudioSamples,
                                       const INT nAudioSamples,
                                       AACENC_EXT_PAYLOAD *pMpsExtPayload) {
  MPS_ENCODER_ERROR error = MPS_ENCODER_OK;

  if (hMpsEnc == NULL) {
    error = MPS_ENCODER_INVALID_HANDLE;
  } else {
    int sacHeaderFlag = 1;
    int sacOutBufferOffset = 0;

    if (hMpsEnc->audioObjectType == AOT_ER_AAC_ELD) {
      sacHeaderFlag = 0;
    }

    hMpsEnc->sacOutBuffer[0] = (sacHeaderFlag == 0) ? 0x3 : 0x7;
    sacOutBufferOffset += 1;

    if (sacHeaderFlag) {
      sacOutBufferOffset += FDK_MpegsEnc_WriteFrameHeader(
          hMpsEnc, &hMpsEnc->sacOutBuffer[sacOutBufferOffset],
          sizeof(hMpsEnc->sacOutBuffer) - sacOutBufferOffset);
    }

    hMpsEnc->pInBuffer[0] = (void *)pAudioSamples;
    hMpsEnc->inargs.nInputSamples = nAudioSamples;

    hMpsEnc->pOutBuffer[0] = (void *)pAudioSamples;
    hMpsEnc->pOutBufferSize[0] = sizeof(INT_PCM) * nAudioSamples / 2;

    hMpsEnc->pOutBuffer[1] = (void *)&hMpsEnc->sacOutBuffer[sacOutBufferOffset];
    hMpsEnc->pOutBufferSize[1] =
        sizeof(hMpsEnc->sacOutBuffer) - sacOutBufferOffset;

    if (SACENC_OK != FDK_sacenc_encode(hMpsEnc->hSacEncoder,
                                       &hMpsEnc->inBufDesc,
                                       &hMpsEnc->outBufDesc, &hMpsEnc->inargs,
                                       &hMpsEnc->outargs)) {
      error = MPS_ENCODER_ENCODE_ERROR;
      goto bail;
    }

    pMpsExtPayload->pData = (UCHAR *)hMpsEnc->sacOutBuffer;
    pMpsExtPayload->dataSize =
        hMpsEnc->outargs.nOutputBits + 8 * (sacOutBufferOffset - 1);
    pMpsExtPayload->dataType = EXT_LDSAC_DATA;
    pMpsExtPayload->associatedChElement = -1;
  }

bail:
  return error;
}

INT FDK_MpegsEnc_WriteSpatialSpecificConfig(HANDLE_MPS_ENCODER hMpsEnc,
                                            HANDLE_FDK_BITSTREAM hBs) {
  INT sscBits = 0;

  if (NULL != hMpsEnc) {
    MP4SPACEENC_INFO mp4SpaceEncoderInfo;
    FDK_sacenc_getInfo(hMpsEnc->hSacEncoder, &mp4SpaceEncoderInfo);

    if (hBs != NULL) {
      int i;
      int writtenBits = 0;
      for (i = 0; i<mp4SpaceEncoderInfo.pSscBuf->nSscSizeBits>> 3; i++) {
        FDKwriteBits(hBs, mp4SpaceEncoderInfo.pSscBuf->pSsc[i], 8);
        writtenBits += 8;
      }
      FDKwriteBits(hBs, mp4SpaceEncoderInfo.pSscBuf->pSsc[i],
                   mp4SpaceEncoderInfo.pSscBuf->nSscSizeBits - writtenBits);
    }

    sscBits = mp4SpaceEncoderInfo.pSscBuf->nSscSizeBits;
  }

  return sscBits;
}

static INT FDK_MpegsEnc_WriteFrameHeader(HANDLE_MPS_ENCODER hMpsEnc,
                                         UCHAR *const pOutputBuffer,
                                         const int outputBufferSize) {
  const int sacTimeAlignFlag = 0;

  int numBits = 0;

  if ((NULL != hMpsEnc) && (NULL != pOutputBuffer)) {
    UINT alignAnchor, cnt;
    FDK_BITSTREAM Bs;
    FDKinitBitStream(&Bs, pOutputBuffer, outputBufferSize, 0, BS_WRITER);

    cnt = (FDK_MpegsEnc_WriteSpatialSpecificConfig(hMpsEnc, NULL) + 7) >> 3;

    FDKwriteBits(&Bs, sacTimeAlignFlag, 1);

    if (cnt < 127) {
      FDKwriteBits(&Bs, cnt, 7);
    } else {
      FDKwriteBits(&Bs, 127, 7);
      FDKwriteBits(&Bs, cnt - 127, 16);
    }

    alignAnchor = FDKgetValidBits(&Bs);
    FDK_MpegsEnc_WriteSpatialSpecificConfig(hMpsEnc, &Bs);
    FDKbyteAlign(&Bs, alignAnchor);

    if (sacTimeAlignFlag) {
      FDK_ASSERT(1);
    }

    numBits = FDKgetValidBits(&Bs);
  }

  return ((numBits + 7) >> 3);
}

INT FDK_MpegsEnc_GetClosestBitRate(const AUDIO_OBJECT_TYPE audioObjectType,
                                   const CHANNEL_MODE channelMode,
                                   const UINT samplingrate, const UINT sbrRatio,
                                   const UINT bitrate) {
  unsigned int i;
  int targetBitrate = -1;

  for (i = 0; i < sizeof(mpsConfigTab) / sizeof(MPS_CONFIG_TAB); i++) {
    if ((mpsConfigTab[i].audio_object_type == audioObjectType) &&
        (mpsConfigTab[i].channel_mode == channelMode) &&
        (mpsConfigTab[i].sbr_ratio == sbrRatio) &&
        (mpsConfigTab[i].sampling_rate == samplingrate)) {
      targetBitrate = fMin(fMax(bitrate, mpsConfigTab[i].bitrate_min),
                           mpsConfigTab[i].bitrate_max);
    }
  }

  return targetBitrate;
}

INT FDK_MpegsEnc_GetDelay(HANDLE_MPS_ENCODER hMpsEnc) {
  INT delay = 0;

  if (NULL != hMpsEnc) {
    MP4SPACEENC_INFO mp4SpaceEncoderInfo;
    FDK_sacenc_getInfo(hMpsEnc->hSacEncoder, &mp4SpaceEncoderInfo);
    delay = mp4SpaceEncoderInfo.nCodecDelay;
  }

  return delay;
}

INT FDK_MpegsEnc_GetDecDelay(HANDLE_MPS_ENCODER hMpsEnc) {
  INT delay = 0;

  if (NULL != hMpsEnc) {
    MP4SPACEENC_INFO mp4SpaceEncoderInfo;
    FDK_sacenc_getInfo(hMpsEnc->hSacEncoder, &mp4SpaceEncoderInfo);
    delay = mp4SpaceEncoderInfo.nDecoderDelay;
  }

  return delay;
}

MPS_ENCODER_ERROR FDK_MpegsEnc_GetLibInfo(LIB_INFO *info) {
  MPS_ENCODER_ERROR error = MPS_ENCODER_OK;

  if (NULL == info) {
    error = MPS_ENCODER_INVALID_HANDLE;
  } else if (SACENC_OK != FDK_sacenc_getLibInfo(info)) {
    error = MPS_ENCODER_INIT_ERROR;
  }

  return error;
}
