/**************************** AAC encoder library ******************************

   Author(s):   M. Schug / A. Groeschel

   Description: fast aac coder functions

*******************************************************************************/

#include "aacenc.h"

#include "bitenc.h"
#include "interface.h"
#include "psy_configuration.h"
#include "psy_main.h"
#include "qc_main.h"
#include "bandwidth.h"
#include "channel_map.h"
#include "tns_func.h"
#include "aacEnc_ram.h"

#include "genericStds.h"

#define BITRES_MIN \
  300
#define BITRES_MAX_LD 4000
#define BITRES_MIN_LD 500
#define BITRATE_MAX_LD 70000
#define BITRATE_MIN_LD 12000

INT FDKaacEnc_CalcBitsPerFrame(const INT bitRate, const INT frameLength,
                               const INT samplingRate) {
  int shift = 0;
  while ((frameLength & ~((1 << (shift + 1)) - 1)) == frameLength &&
         (samplingRate & ~((1 << (shift + 1)) - 1)) == samplingRate) {
    shift++;
  }

  return (bitRate * (frameLength >> shift)) / (samplingRate >> shift);
}

INT FDKaacEnc_CalcBitrate(const INT bitsPerFrame, const INT frameLength,
                          const INT samplingRate) {
  int shift = 0;
  while ((frameLength & ~((1 << (shift + 1)) - 1)) == frameLength &&
         (samplingRate & ~((1 << (shift + 1)) - 1)) == samplingRate) {
    shift++;
  }

  return (bitsPerFrame * (samplingRate >> shift)) / (frameLength >> shift);
}

static AAC_ENCODER_ERROR FDKaacEnc_InitCheckAncillary(
    INT bitRate, INT framelength, INT ancillaryRate, INT *ancillaryBitsPerFrame,
    INT sampleRate);

INT FDKaacEnc_LimitBitrate(HANDLE_TRANSPORTENC hTpEnc, AUDIO_OBJECT_TYPE aot,
                           INT coreSamplingRate, INT frameLength, INT nChannels,
                           INT nChannelsEff, INT bitRate, INT averageBits,
                           INT *pAverageBitsPerFrame,
                           AACENC_BITRATE_MODE bitrateMode, INT nSubFrames) {
  INT transportBits, prevBitRate, averageBitsPerFrame, minBitrate = 0, iter = 0;
  INT minBitsPerFrame = 40 * nChannels;
  if (isLowDelay(aot)) {
    minBitrate = 8000 * nChannelsEff;
  }

  do {
    prevBitRate = bitRate;
    averageBitsPerFrame =
        FDKaacEnc_CalcBitsPerFrame(bitRate, frameLength, coreSamplingRate) /
        nSubFrames;

    if (pAverageBitsPerFrame != NULL) {
      *pAverageBitsPerFrame = averageBitsPerFrame;
    }

    if (hTpEnc != NULL) {
      transportBits = transportEnc_GetStaticBits(hTpEnc, averageBitsPerFrame);
    } else {
      transportBits = 208;
    }

    bitRate = fMax(bitRate,
                   fMax(minBitrate,
                        FDKaacEnc_CalcBitrate((minBitsPerFrame + transportBits),
                                              frameLength, coreSamplingRate)));
    FDK_ASSERT(bitRate >= 0);

    bitRate = fMin(bitRate, FDKaacEnc_CalcBitrate(
                                (nChannelsEff * MIN_BUFSIZE_PER_EFF_CHAN),
                                frameLength, coreSamplingRate));
    FDK_ASSERT(bitRate >= 0);

  } while (prevBitRate != bitRate && iter++ < 3);

  return bitRate;
}

typedef struct {
  AACENC_BITRATE_MODE bitrateMode;
  int chanBitrate[2];
} CONFIG_TAB_ENTRY_VBR;

static const CONFIG_TAB_ENTRY_VBR configTabVBR[] = {
    {AACENC_BR_MODE_CBR, {0, 0}},
    {AACENC_BR_MODE_VBR_1, {32000, 20000}},
    {AACENC_BR_MODE_VBR_2, {40000, 32000}},
    {AACENC_BR_MODE_VBR_3, {56000, 48000}},
    {AACENC_BR_MODE_VBR_4, {72000, 64000}},
    {AACENC_BR_MODE_VBR_5, {112000, 96000}}};

INT FDKaacEnc_GetVBRBitrate(AACENC_BITRATE_MODE bitrateMode,
                            CHANNEL_MODE channelMode) {
  INT bitrate = 0;
  INT monoStereoMode = 0;

  if (FDKaacEnc_GetMonoStereoMode(channelMode) == EL_MODE_STEREO) {
    monoStereoMode = 1;
  }

  switch (bitrateMode) {
    case AACENC_BR_MODE_VBR_1:
    case AACENC_BR_MODE_VBR_2:
    case AACENC_BR_MODE_VBR_3:
    case AACENC_BR_MODE_VBR_4:
    case AACENC_BR_MODE_VBR_5:
      bitrate = configTabVBR[bitrateMode].chanBitrate[monoStereoMode];
      break;
    case AACENC_BR_MODE_INVALID:
    case AACENC_BR_MODE_CBR:
    case AACENC_BR_MODE_SFR:
    case AACENC_BR_MODE_FF:
    default:
      bitrate = 0;
      break;
  }

  bitrate *= FDKaacEnc_GetChannelModeConfiguration(channelMode)->nChannelsEff;

  return bitrate;
}

AACENC_BITRATE_MODE FDKaacEnc_AdjustVBRBitrateMode(
    AACENC_BITRATE_MODE bitrateMode, INT bitrate, CHANNEL_MODE channelMode) {
  AACENC_BITRATE_MODE newBitrateMode = bitrateMode;

  if (bitrate != -1) {
    const INT monoStereoMode =
        (FDKaacEnc_GetMonoStereoMode(channelMode) == EL_MODE_STEREO) ? 1 : 0;
    const INT nChannelsEff =
        FDKaacEnc_GetChannelModeConfiguration(channelMode)->nChannelsEff;
    newBitrateMode = AACENC_BR_MODE_INVALID;

    for (int idx = (int)(sizeof(configTabVBR) / sizeof(*configTabVBR)) - 1;
         idx >= 0; idx--) {
      if (bitrate >=
          configTabVBR[idx].chanBitrate[monoStereoMode] * nChannelsEff) {
        if (configTabVBR[idx].chanBitrate[monoStereoMode] * nChannelsEff <
            FDKaacEnc_GetVBRBitrate(bitrateMode, channelMode)) {
          newBitrateMode = configTabVBR[idx].bitrateMode;
        } else {
          newBitrateMode = bitrateMode;
        }
        break;
      }
    }
  }

  return AACENC_BR_MODE_IS_VBR(newBitrateMode) ? newBitrateMode
                                               : AACENC_BR_MODE_INVALID;
}

static INT FDKaacEnc_EncBitresToTpBitres(const HANDLE_AAC_ENC hAacEnc) {
  INT transportBitreservoir = 0;

  switch (hAacEnc->bitrateMode) {
    case AACENC_BR_MODE_CBR:
      transportBitreservoir =
          hAacEnc->qcKernel->bitResTot;
      break;
    case AACENC_BR_MODE_VBR_1:
    case AACENC_BR_MODE_VBR_2:
    case AACENC_BR_MODE_VBR_3:
    case AACENC_BR_MODE_VBR_4:
    case AACENC_BR_MODE_VBR_5:
      transportBitreservoir = FDK_INT_MAX;
      break;
    case AACENC_BR_MODE_SFR:
      transportBitreservoir = 0;
      break;
    default:
    case AACENC_BR_MODE_INVALID:
      transportBitreservoir = 0;
  }

  if (hAacEnc->config->audioMuxVersion == 2) {
    transportBitreservoir =
        MIN_BUFSIZE_PER_EFF_CHAN * hAacEnc->channelMapping.nChannelsEff;
  }

  return transportBitreservoir;
}

INT FDKaacEnc_GetBitReservoirState(const HANDLE_AAC_ENC hAacEncoder) {
  return FDKaacEnc_EncBitresToTpBitres(hAacEncoder);
}

void FDKaacEnc_AacInitDefaultConfig(AACENC_CONFIG *config) {
  FDKmemclear(config, sizeof(AACENC_CONFIG));

  config->anc_Rate = 0;
  config->ancDataBitRate = 0;

  config->bitRate = -1;
  config->averageBits =
      -1;
  config->bitrateMode =
      AACENC_BR_MODE_CBR;
  config->bandWidth = 0;
  config->useTns = TNS_ENABLE_MASK;
  config->usePns =
      1;
  config->useIS = 1;
  config->useMS = 1;
  config->framelength = -1;
  config->syntaxFlags = 0;
  config->epConfig = -1;
  config->nSubFrames = 1;
  config->channelOrder = CH_ORDER_MPEG;
  config->channelMode = MODE_UNKNOWN;
  config->minBitsPerFrame = -1;
  config->maxBitsPerFrame = -1;
  config->audioMuxVersion = -1;
  config->downscaleFactor =
      1;
}

AAC_ENCODER_ERROR FDKaacEnc_Open(HANDLE_AAC_ENC *phAacEnc, const INT nElements,
                                 const INT nChannels, const INT nSubFrames) {
  AAC_ENCODER_ERROR ErrorStatus;
  AAC_ENC *hAacEnc = NULL;
  UCHAR *dynamicRAM = NULL;

  if (phAacEnc == NULL) {
    return AAC_ENC_INVALID_HANDLE;
  }

  hAacEnc = GetRam_aacEnc_AacEncoder();
  if (hAacEnc == NULL) {
    ErrorStatus = AAC_ENC_NO_MEMORY;
    goto bail;
  }
  FDKmemclear(hAacEnc, sizeof(AAC_ENC));

  if (NULL == (hAacEnc->dynamic_RAM = GetAACdynamic_RAM())) {
    ErrorStatus = AAC_ENC_NO_MEMORY;
    goto bail;
  }
  dynamicRAM = (UCHAR *)hAacEnc->dynamic_RAM;

  ErrorStatus =
      FDKaacEnc_PsyNew(&hAacEnc->psyKernel, nElements, nChannels, dynamicRAM);
  if (ErrorStatus != AAC_ENC_OK) goto bail;

  ErrorStatus = FDKaacEnc_PsyOutNew(hAacEnc->psyOut, nElements, nChannels,
                                    nSubFrames, dynamicRAM);
  if (ErrorStatus != AAC_ENC_OK) goto bail;

  ErrorStatus = FDKaacEnc_QCOutNew(hAacEnc->qcOut, nElements, nChannels,
                                   nSubFrames, dynamicRAM);
  if (ErrorStatus != AAC_ENC_OK) goto bail;

  ErrorStatus = FDKaacEnc_QCNew(&hAacEnc->qcKernel, nElements, dynamicRAM);
  if (ErrorStatus != AAC_ENC_OK) goto bail;

  hAacEnc->maxChannels = nChannels;
  hAacEnc->maxElements = nElements;
  hAacEnc->maxFrames = nSubFrames;

bail:
  *phAacEnc = hAacEnc;
  return ErrorStatus;
}

AAC_ENCODER_ERROR FDKaacEnc_Initialize(
    HANDLE_AAC_ENC hAacEnc,
    AACENC_CONFIG *config,
    HANDLE_TRANSPORTENC hTpEnc, ULONG initFlags) {
  AAC_ENCODER_ERROR ErrorStatus;
  INT psyBitrate, tnsMask;
  CHANNEL_MAPPING *cm = NULL;

  INT mbfac_e, qbw;
  FIXP_DBL mbfac, bw_ratio;
  QC_INIT qcInit;
  INT averageBitsPerFrame = 0;
  const CHANNEL_MODE prevChannelMode = hAacEnc->encoderMode;

  if (config == NULL) return AAC_ENC_INVALID_HANDLE;

  if (config->nChannels < 1 || config->nChannels > (8)) {
    return AAC_ENC_UNSUPPORTED_CHANNELCONFIG;
  }

  switch (config->sampleRate) {
    case 8000:
    case 11025:
    case 12000:
    case 16000:
    case 22050:
    case 24000:
    case 32000:
    case 44100:
    case 48000:
    case 64000:
    case 88200:
    case 96000:
      break;
    default:
      return AAC_ENC_UNSUPPORTED_SAMPLINGRATE;
  }

  if (config->bitRate == -1) {
    return AAC_ENC_UNSUPPORTED_BITRATE;
  }

  if (FDKaacEnc_LimitBitrate(
          hTpEnc, config->audioObjectType, config->sampleRate,
          config->framelength, config->nChannels,
          FDKaacEnc_GetChannelModeConfiguration(config->channelMode)
              ->nChannelsEff,
          config->bitRate, config->averageBits, &averageBitsPerFrame,
          config->bitrateMode, config->nSubFrames) != config->bitRate &&
      !(AACENC_BR_MODE_IS_VBR(config->bitrateMode))) {
    return AAC_ENC_UNSUPPORTED_BITRATE;
  }

  if (config->syntaxFlags & AC_ER_VCB11) {
    return AAC_ENC_UNSUPPORTED_ER_FORMAT;
  }
  if (config->syntaxFlags & AC_ER_HCR) {
    return AAC_ENC_UNSUPPORTED_ER_FORMAT;
  }

  switch (config->framelength) {
    case 1024:
      if (isLowDelay(config->audioObjectType)) {
        return AAC_ENC_INVALID_FRAME_LENGTH;
      }
      break;
    case 128:
    case 256:
    case 512:
    case 120:
    case 240:
    case 480:
      if (!isLowDelay(config->audioObjectType)) {
        return AAC_ENC_INVALID_FRAME_LENGTH;
      }
      break;
    default:
      return AAC_ENC_INVALID_FRAME_LENGTH;
  }

  if (config->anc_Rate != 0) {
    ErrorStatus = FDKaacEnc_InitCheckAncillary(
        config->bitRate, config->framelength, config->anc_Rate,
        &hAacEnc->ancillaryBitsPerFrame, config->sampleRate);
    if (ErrorStatus != AAC_ENC_OK) goto bail;

    config->ancDataBitRate +=
        FDKaacEnc_CalcBitrate(hAacEnc->ancillaryBitsPerFrame,
                              config->framelength, config->sampleRate);
  }

  config->maxAncBytesPerAU =
      fMin((256), fMax(0, FDKaacEnc_CalcBitsPerFrame(
                              (config->bitRate - (config->nChannels * 8000)),
                              config->framelength, config->sampleRate) >>
                              3));

  hAacEnc->config = config;

  hAacEnc->bitrateMode = config->bitrateMode;

  hAacEnc->encoderMode = config->channelMode;

  ErrorStatus = FDKaacEnc_InitChannelMapping(
      hAacEnc->encoderMode, config->channelOrder, &hAacEnc->channelMapping);
  if (ErrorStatus != AAC_ENC_OK) goto bail;

  cm = &hAacEnc->channelMapping;

  ErrorStatus = FDKaacEnc_DetermineBandWidth(
      config->bandWidth, config->bitRate - config->ancDataBitRate,
      hAacEnc->bitrateMode, config->sampleRate, config->framelength, cm,
      hAacEnc->encoderMode, &hAacEnc->config->bandWidth);
  if (ErrorStatus != AAC_ENC_OK) goto bail;

  hAacEnc->bandwidth90dB = (INT)hAacEnc->config->bandWidth;

  tnsMask = config->useTns ? TNS_ENABLE_MASK : 0x0;
  psyBitrate = config->bitRate - config->ancDataBitRate;

  if ((hAacEnc->encoderMode != prevChannelMode) || (initFlags != 0)) {
    ErrorStatus = FDKaacEnc_psyInit(hAacEnc->psyKernel, hAacEnc->psyOut,
                                    hAacEnc->maxFrames, hAacEnc->maxChannels,
                                    config->audioObjectType, cm);
    if (ErrorStatus != AAC_ENC_OK) goto bail;
  }

  ErrorStatus = FDKaacEnc_psyMainInit(
      hAacEnc->psyKernel, config->audioObjectType, cm, config->sampleRate,
      config->framelength, psyBitrate, tnsMask, hAacEnc->bandwidth90dB,
      config->usePns, config->useIS, config->useMS, config->syntaxFlags,
      initFlags);
  if (ErrorStatus != AAC_ENC_OK) goto bail;

  ErrorStatus = FDKaacEnc_QCOutInit(hAacEnc->qcOut, hAacEnc->maxFrames, cm);
  if (ErrorStatus != AAC_ENC_OK) goto bail;

  qcInit.channelMapping = &hAacEnc->channelMapping;
  qcInit.sceCpe = 0;

  if (AACENC_BR_MODE_IS_VBR(config->bitrateMode)) {
    qcInit.averageBits = (averageBitsPerFrame + 7) & ~7;
    qcInit.bitRes = MIN_BUFSIZE_PER_EFF_CHAN * cm->nChannelsEff;
    qcInit.maxBits = MIN_BUFSIZE_PER_EFF_CHAN * cm->nChannelsEff;
    qcInit.maxBits = (config->maxBitsPerFrame != -1)
                         ? fixMin(qcInit.maxBits, config->maxBitsPerFrame)
                         : qcInit.maxBits;
    qcInit.maxBits = fixMax(qcInit.maxBits, (averageBitsPerFrame + 7) & ~7);
    qcInit.minBits =
        (config->minBitsPerFrame != -1) ? config->minBitsPerFrame : 0;
    qcInit.minBits = fixMin(qcInit.minBits, averageBitsPerFrame & ~7);
  } else {
    INT bitreservoir = -1;
    if (isLowDelay(config->audioObjectType)) {
      INT brPerChannel = config->bitRate / config->nChannels;
      brPerChannel = fMin(BITRATE_MAX_LD, fMax(BITRATE_MIN_LD, brPerChannel));

      FIXP_DBL slope = fDivNorm(
          (brPerChannel - BITRATE_MIN_LD),
          BITRATE_MAX_LD - BITRATE_MIN_LD);
      bitreservoir = fMultI(slope, (INT)(BITRES_MAX_LD - BITRES_MIN_LD)) +
                     BITRES_MIN_LD;
      bitreservoir = bitreservoir & ~7;
    }

    int maxBitres;
    qcInit.averageBits = (averageBitsPerFrame + 7) & ~7;
    maxBitres =
        (MIN_BUFSIZE_PER_EFF_CHAN * cm->nChannelsEff) - qcInit.averageBits;
    qcInit.bitRes =
        (bitreservoir != -1) ? fMin(bitreservoir, maxBitres) : maxBitres;

    qcInit.maxBits = fixMin(MIN_BUFSIZE_PER_EFF_CHAN * cm->nChannelsEff,
                            ((averageBitsPerFrame + 7) & ~7) + qcInit.bitRes);
    qcInit.maxBits = (config->maxBitsPerFrame != -1)
                         ? fixMin(qcInit.maxBits, config->maxBitsPerFrame)
                         : qcInit.maxBits;
    qcInit.maxBits =
        fixMin(MIN_BUFSIZE_PER_EFF_CHAN * cm->nChannelsEff,
               fixMax(qcInit.maxBits, (averageBitsPerFrame + 7 + 8) & ~7));

    qcInit.minBits = fixMax(
        0, ((averageBitsPerFrame - 1) & ~7) - qcInit.bitRes -
               transportEnc_GetStaticBits(
                   hTpEnc, ((averageBitsPerFrame + 7) & ~7) + qcInit.bitRes));
    qcInit.minBits = (config->minBitsPerFrame != -1)
                         ? fixMax(qcInit.minBits, config->minBitsPerFrame)
                         : qcInit.minBits;
    qcInit.minBits = fixMin(
        qcInit.minBits, (averageBitsPerFrame -
                         transportEnc_GetStaticBits(hTpEnc, qcInit.maxBits)) &
                            ~7);
  }

  qcInit.sampleRate = config->sampleRate;
  qcInit.isLowDelay = isLowDelay(config->audioObjectType) ? 1 : 0;
  qcInit.nSubFrames = config->nSubFrames;
  qcInit.padding.paddingRest = config->sampleRate;

  if (qcInit.maxBits - qcInit.averageBits >=
      ((qcInit.isLowDelay) ? BITRES_MIN_LD : BITRES_MIN) * config->nChannels) {
    qcInit.bitResMode = AACENC_BR_MODE_FULL;
  } else if (qcInit.maxBits > qcInit.averageBits) {
    qcInit.bitResMode = AACENC_BR_MODE_REDUCED;
  } else {
    qcInit.bitResMode = AACENC_BR_MODE_DISABLED;
  }

  switch (config->channelMode) {
    case MODE_1_2:
    case MODE_1_2_1:
    case MODE_1_2_2:
    case MODE_1_2_2_1:
    case MODE_6_1:
    case MODE_1_2_2_2_1:
    case MODE_7_1_BACK:
    case MODE_7_1_TOP_FRONT:
    case MODE_7_1_REAR_SURROUND:
    case MODE_7_1_FRONT_CENTER:
      qcInit.bitDistributionMode = 0;
      break;
    case MODE_1:
    case MODE_2:
    default:
      qcInit.bitDistributionMode = 1;
  }

  bw_ratio =
      fDivNorm((FIXP_DBL)(10 * config->framelength * hAacEnc->bandwidth90dB),
               (FIXP_DBL)(config->sampleRate), &qbw);
  qcInit.meanPe =
      fMax((INT)scaleValue(bw_ratio, qbw + 1 - (DFRACT_BITS - 1)), 1);

  mbfac = fDivNorm(qcInit.maxBits, qcInit.averageBits / qcInit.nSubFrames,
                   &mbfac_e);
  qcInit.maxBitFac = scaleValue(mbfac, -(DFRACT_BITS - 1 - 24 - mbfac_e));

  switch (config->bitrateMode) {
    case AACENC_BR_MODE_CBR:
      qcInit.bitrateMode = QCDATA_BR_MODE_CBR;
      break;
    case AACENC_BR_MODE_VBR_1:
      qcInit.bitrateMode = QCDATA_BR_MODE_VBR_1;
      break;
    case AACENC_BR_MODE_VBR_2:
      qcInit.bitrateMode = QCDATA_BR_MODE_VBR_2;
      break;
    case AACENC_BR_MODE_VBR_3:
      qcInit.bitrateMode = QCDATA_BR_MODE_VBR_3;
      break;
    case AACENC_BR_MODE_VBR_4:
      qcInit.bitrateMode = QCDATA_BR_MODE_VBR_4;
      break;
    case AACENC_BR_MODE_VBR_5:
      qcInit.bitrateMode = QCDATA_BR_MODE_VBR_5;
      break;
    case AACENC_BR_MODE_SFR:
      qcInit.bitrateMode = QCDATA_BR_MODE_SFR;
      break;
    case AACENC_BR_MODE_FF:
      qcInit.bitrateMode = QCDATA_BR_MODE_FF;
      break;
    default:
      ErrorStatus = AAC_ENC_UNSUPPORTED_BITRATE_MODE;
      goto bail;
  }

  qcInit.invQuant = (config->useRequant) ? 2 : 0;

  if (isLowDelay(config->audioObjectType)) {
    qcInit.maxIterations = 2;
  } else {
    qcInit.maxIterations = 5;
  }

  qcInit.bitrate = config->bitRate - config->ancDataBitRate;

  qcInit.staticBits = transportEnc_GetStaticBits(
      hTpEnc, qcInit.averageBits / qcInit.nSubFrames);

  ErrorStatus = FDKaacEnc_QCInit(hAacEnc->qcKernel, &qcInit, initFlags);
  if (ErrorStatus != AAC_ENC_OK) goto bail;

  switch (hAacEnc->config->audioObjectType) {
    case AOT_MP2_AAC_LC:
      hAacEnc->aot = AOT_AAC_LC;
      break;
    case AOT_MP2_SBR:
      hAacEnc->aot = AOT_SBR;
      break;
    default:
      hAacEnc->aot = hAacEnc->config->audioObjectType;
  }

  return AAC_ENC_OK;

bail:

  return ErrorStatus;
}

AAC_ENCODER_ERROR FDKaacEnc_EncodeFrame(
    HANDLE_AAC_ENC hAacEnc,
    HANDLE_TRANSPORTENC hTpEnc, INT_PCM *RESTRICT inputBuffer,
    const UINT inputBufferBufSize, INT *nOutBytes,
    AACENC_EXT_PAYLOAD extPayload[MAX_TOTAL_EXT_PAYLOADS]) {
  AAC_ENCODER_ERROR ErrorStatus;
  int el, n, c = 0;
  UCHAR extPayloadUsed[MAX_TOTAL_EXT_PAYLOADS];

  CHANNEL_MAPPING *cm = &hAacEnc->channelMapping;

  PSY_OUT *psyOut = hAacEnc->psyOut[c];
  QC_OUT *qcOut = hAacEnc->qcOut[c];

  FDKmemclear(extPayloadUsed, MAX_TOTAL_EXT_PAYLOADS * sizeof(UCHAR));

  qcOut->elementExtBits = 0;
  qcOut->staticBits = 0;
  qcOut->totalNoRedPe = 0;

  for (el = 0; el < cm->nElements; el++) {
    ELEMENT_INFO elInfo = cm->elInfo[el];

    if ((elInfo.elType == ID_SCE) || (elInfo.elType == ID_CPE) ||
        (elInfo.elType == ID_LFE)) {
      int ch;

      for (ch = 0; ch < elInfo.nChannelsInEl; ch++) {
        PSY_OUT_CHANNEL *psyOutChan =
            psyOut->psyOutElement[el]->psyOutChannel[ch];
        QC_OUT_CHANNEL *qcOutChan = qcOut->qcElement[el]->qcOutChannel[ch];

        psyOutChan->mdctSpectrum = qcOutChan->mdctSpectrum;
        psyOutChan->sfbSpreadEnergy = qcOutChan->sfbSpreadEnergy;
        psyOutChan->sfbEnergy = qcOutChan->sfbEnergy;
        psyOutChan->sfbEnergyLdData = qcOutChan->sfbEnergyLdData;
        psyOutChan->sfbMinSnrLdData = qcOutChan->sfbMinSnrLdData;
        psyOutChan->sfbThresholdLdData = qcOutChan->sfbThresholdLdData;
      }

      ErrorStatus = FDKaacEnc_psyMain(
          elInfo.nChannelsInEl, hAacEnc->psyKernel->psyElement[el],
          hAacEnc->psyKernel->psyDynamic, hAacEnc->psyKernel->psyConf,
          psyOut->psyOutElement[el], inputBuffer, inputBufferBufSize,
          cm->elInfo[el].ChannelIndex, cm->nChannels);

      if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;

      ErrorStatus = FDKaacEnc_QCMainPrepare(
          &elInfo, hAacEnc->qcKernel->hAdjThr->adjThrStateElem[el],
          psyOut->psyOutElement[el], qcOut->qcElement[el], hAacEnc->aot,
          hAacEnc->config->syntaxFlags, hAacEnc->config->epConfig);

      if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;

      qcOut->qcElement[el]->extBitsUsed = 0;
      qcOut->qcElement[el]->nExtensions = 0;
      FDKmemclear(&qcOut->qcElement[el]->extension,
                  (1) * sizeof(QC_OUT_EXTENSION));

      for (n = 0; n < MAX_TOTAL_EXT_PAYLOADS; n++) {
        if (!extPayloadUsed[n] && (extPayload[n].associatedChElement == el) &&
            (extPayload[n].dataSize > 0) && (extPayload[n].pData != NULL)) {
          int idx = qcOut->qcElement[el]->nExtensions++;

          qcOut->qcElement[el]->extension[idx].type =
              extPayload[n].dataType;
          qcOut->qcElement[el]->extension[idx].nPayloadBits =
              extPayload[n].dataSize;
          qcOut->qcElement[el]->extension[idx].pPayload = extPayload[n].pData;
          qcOut->qcElement[el]->extBitsUsed += FDKaacEnc_writeExtensionData(
              NULL, &qcOut->qcElement[el]->extension[idx], 0, 0,
              hAacEnc->config->syntaxFlags, hAacEnc->aot,
              hAacEnc->config->epConfig);
          extPayloadUsed[n] = 1;
        }
      }

      qcOut->elementExtBits += qcOut->qcElement[el]->extBitsUsed;
      qcOut->staticBits += qcOut->qcElement[el]->staticBitsUsed;

      qcOut->totalNoRedPe += qcOut->qcElement[el]->peData.pe;
    }
  }

  qcOut->nExtensions = 0;
  qcOut->globalExtBits = 0;

  FDKmemclear(&qcOut->extension, (2 + 2) * sizeof(QC_OUT_EXTENSION));

  for (n = 0; n < MAX_TOTAL_EXT_PAYLOADS; n++) {
    if (!extPayloadUsed[n] && (extPayload[n].associatedChElement == -1) &&
        (extPayload[n].pData != NULL)) {
      UINT payloadBits = 0;

      if (extPayload[n].dataType == EXT_DATA_ELEMENT) {
        if (hAacEnc->ancillaryBitsPerFrame) {
          payloadBits = hAacEnc->ancillaryBitsPerFrame;
        } else {
          if ((extPayload[n].dataSize >> 3) <=
              hAacEnc->config->maxAncBytesPerAU) {
            payloadBits = extPayload[n].dataSize;
          }
        }
        payloadBits = fixMin(extPayload[n].dataSize, payloadBits);
      } else {
        payloadBits = extPayload[n].dataSize;
      }

      if (payloadBits > 0) {
        int idx = qcOut->nExtensions++;

        qcOut->extension[idx].type =
            extPayload[n].dataType;
        qcOut->extension[idx].nPayloadBits = payloadBits;
        qcOut->extension[idx].pPayload = extPayload[n].pData;
        qcOut->globalExtBits += FDKaacEnc_writeExtensionData(
            NULL, &qcOut->extension[idx], 0, 0, hAacEnc->config->syntaxFlags,
            hAacEnc->aot, hAacEnc->config->epConfig);
        if (extPayload[n].dataType == EXT_DATA_ELEMENT) {
          extPayload[n].dataSize -= payloadBits;
        }
        extPayloadUsed[n] = 1;
      }
    }
  }

  if (!(hAacEnc->config->syntaxFlags & (AC_SCALABLE | AC_ER))) {
    qcOut->globalExtBits += EL_ID_BITS;
  }

  {
    INT totalBits = 0;
    INT avgTotalBits = 0;

    {
      FDKaacEnc_AdjustBitrate(
          hAacEnc->qcKernel, cm, &avgTotalBits, hAacEnc->config->bitRate,
          hAacEnc->config->sampleRate, hAacEnc->config->framelength);

      avgTotalBits *= hAacEnc->config->nSubFrames;
    }

    hAacEnc->qcKernel->globHdrBits = transportEnc_GetStaticBits(
        hTpEnc, avgTotalBits + hAacEnc->qcKernel->bitResTot);

    ErrorStatus = FDKaacEnc_QCMain(
        hAacEnc->qcKernel, hAacEnc->psyOut, hAacEnc->qcOut, avgTotalBits, cm,
        hAacEnc->aot, hAacEnc->config->syntaxFlags, hAacEnc->config->epConfig);

    if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;

    ErrorStatus = FDKaacEnc_updateFillBits(
        cm, hAacEnc->qcKernel, hAacEnc->qcKernel->elementBits, hAacEnc->qcOut);
    if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;

    ErrorStatus = FDKaacEnc_FinalizeBitConsumption(
        cm, hAacEnc->qcKernel, qcOut, qcOut->qcElement, hTpEnc, hAacEnc->aot,
        hAacEnc->config->syntaxFlags, hAacEnc->config->epConfig);
    if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;
    totalBits += qcOut->totalBits;

    FDKaacEnc_updateBitres(cm, hAacEnc->qcKernel, hAacEnc->qcOut);

    if (TRANSPORTENC_OK !=
        transportEnc_WriteAccessUnit(hTpEnc, totalBits,
                                     FDKaacEnc_EncBitresToTpBitres(hAacEnc),
                                     cm->nChannelsEff)) {
      return AAC_ENC_UNKNOWN;
    }

    ErrorStatus = FDKaacEnc_WriteBitstream(
        hTpEnc, cm, qcOut, psyOut, hAacEnc->qcKernel, hAacEnc->aot,
        hAacEnc->config->syntaxFlags, hAacEnc->config->epConfig);

    if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;

    if (TRANSPORTENC_OK != transportEnc_GetFrame(hTpEnc, nOutBytes)) {
      return AAC_ENC_UNKNOWN;
    }
  }

  return AAC_ENC_OK;
}

void FDKaacEnc_Close(HANDLE_AAC_ENC *phAacEnc)
{
  if (*phAacEnc == NULL) {
    return;
  }
  AAC_ENC *hAacEnc = (AAC_ENC *)*phAacEnc;

  if (hAacEnc->dynamic_RAM != NULL) FreeAACdynamic_RAM(&hAacEnc->dynamic_RAM);

  FDKaacEnc_PsyClose(&hAacEnc->psyKernel, hAacEnc->psyOut);

  FDKaacEnc_QCClose(&hAacEnc->qcKernel, hAacEnc->qcOut);

  FreeRam_aacEnc_AacEncoder(phAacEnc);
}

#define MAX_ANCRATE 19200

static AAC_ENCODER_ERROR FDKaacEnc_InitCheckAncillary(
    INT bitRate, INT framelength, INT ancillaryRate, INT *ancillaryBitsPerFrame,
    INT sampleRate) {
  if (ancillaryRate < -1) return AAC_ENC_UNSUPPORTED_ANC_BITRATE;

  if ((ancillaryRate != (-1)) && (ancillaryRate != 0)) {
    if ((ancillaryRate >= MAX_ANCRATE) ||
        ((ancillaryRate * 20) > (bitRate * 3))) {
      return AAC_ENC_UNSUPPORTED_ANC_BITRATE;
    }
  } else if (ancillaryRate == -1) {
    if (bitRate >= (MAX_ANCRATE * 10)) {
      ancillaryRate = (MAX_ANCRATE - 1);
    } else {
      ancillaryRate = bitRate / 10;
    }
  }

  *ancillaryBitsPerFrame =
      FDKaacEnc_CalcBitsPerFrame(ancillaryRate, framelength, sampleRate) & ~0x7;

  return AAC_ENC_OK;
}
