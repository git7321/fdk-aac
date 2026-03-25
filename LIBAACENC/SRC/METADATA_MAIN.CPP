/**************************** AAC encoder library ******************************

   Author(s):   V. Bacigalupo

   Description: Metadata Encoder library interface functions

*******************************************************************************/

#include "metadata_main.h"
#include "metadata_compressor.h"
#include "FDK_bitstream.h"
#include "FDK_audio.h"
#include "genericStds.h"

#define MAX_DRC_BANDS (1 << 4)
#define MAX_DRC_FRAMELEN (2 * 1024)
#define MAX_DELAY_FRAMES (3)

typedef struct AAC_METADATA {
  struct {
    UCHAR prog_ref_level_present;
    SCHAR prog_ref_level;

    UCHAR dyn_rng_sgn[MAX_DRC_BANDS];
    UCHAR dyn_rng_ctl[MAX_DRC_BANDS];

    UCHAR drc_bands_present;
    UCHAR drc_band_incr;
    UCHAR drc_band_top[MAX_DRC_BANDS];
    UCHAR drc_interpolation_scheme;
    AACENC_METADATA_DRC_PROFILE drc_profile;
    INT drc_TargetRefLevel;

    UCHAR excluded_chns_present;
    UCHAR exclude_mask[2];
  } mpegDrc;

  struct {
    UCHAR compression_on;
    UCHAR compression_value;
    AACENC_METADATA_DRC_PROFILE comp_profile;
    INT comp_TargetRefLevel;
    INT timecode_coarse_status;
    INT timecode_fine_status;

    UCHAR extAncDataStatus;

    struct {
      UCHAR ext_downmix_lvl_status;
      UCHAR ext_downmix_gain_status;
      UCHAR ext_lfe_downmix_status;
      UCHAR
      ext_dmix_a_idx;
      UCHAR
      ext_dmix_b_idx;
      UCHAR dmx_gain_5_sgn;
      UCHAR dmx_gain_5_idx;
      UCHAR dmx_gain_2_sgn;
      UCHAR dmx_gain_2_idx;
      UCHAR ext_dmix_lfe_idx;

    } extAncData;

  } etsiAncData;

  SCHAR centerMixLevel;
  SCHAR
  surroundMixLevel;
  UCHAR WritePCEMixDwnIdx;
  UCHAR DmxLvl_On;

  UCHAR dolbySurroundMode;
  UCHAR drcPresentationMode;

  UCHAR
  metadataMode;

} AAC_METADATA;

typedef struct FDK_METADATA_ENCODER {
  INT metadataMode;
  HDRC_COMP hDrcComp;
  AACENC_MetaData submittedMetaData;

  INT nAudioDataDelay;
  INT nMetaDataDelay;
  INT nChannels;
  CHANNEL_MODE channelMode;

  INT_PCM* pAudioDelayBuffer;

  AAC_METADATA metaDataBuffer[MAX_DELAY_FRAMES];
  INT metaDataDelayIdx;

  UCHAR drcInfoPayload[12];
  UCHAR drcDsePayload[8];

  INT matrix_mixdown_idx;

  AACENC_EXT_PAYLOAD exPayload[2];
  INT nExtensions;

  UINT maxChannels;

  INT finalizeMetaData;
  INT initializeMetaData;
} FDK_METADATA_ENCODER;

static const AACENC_MetaData defaultMetaDataSetup = {
    AACENC_METADATA_DRC_NONE,
    AACENC_METADATA_DRC_NOT_PRESENT,
    -(31 << 16),
    -(23 << 16),
    0,
    -(23 << 16),
    0,
    0,
    0,
    0,
    0,
    0,
    {0, 0, 0, 0, 0, 0, 0, 0, 0}
};

static const FIXP_DBL dmxTable[8] = {
    ((FIXP_DBL)MAXVAL_DBL), FL2FXCONST_DBL(0.841f), FL2FXCONST_DBL(0.707f),
    FL2FXCONST_DBL(0.596f), FL2FXCONST_DBL(0.500f), FL2FXCONST_DBL(0.422f),
    FL2FXCONST_DBL(0.355f), FL2FXCONST_DBL(0.000f)};

#define FL2DMXLFE(a) FL2FXCONST_DBL((a) / (1 << LFE_LEV_SCALE))
static const FIXP_DBL dmxLfeTable[16] = {
    FL2DMXLFE(3.162f), FL2DMXLFE(2.000f), FL2DMXLFE(1.679f), FL2DMXLFE(1.413f),
    FL2DMXLFE(1.189f), FL2DMXLFE(1.000f), FL2DMXLFE(0.841f), FL2DMXLFE(0.707f),
    FL2DMXLFE(0.596f), FL2DMXLFE(0.500f), FL2DMXLFE(0.316f), FL2DMXLFE(0.178f),
    FL2DMXLFE(0.100f), FL2DMXLFE(0.032f), FL2DMXLFE(0.010f), FL2DMXLFE(0.000f)};

static const UCHAR surmix2matrix_mixdown_idx[8] = {0, 0, 0, 1, 1, 2, 2, 3};

static FDK_METADATA_ERROR WriteMetadataPayload(
    const HANDLE_FDK_METADATA_ENCODER hMetaData,
    const AAC_METADATA* const pMetadata);

static INT WriteDynamicRangeInfoPayload(const AAC_METADATA* const pMetadata,
                                        UCHAR* const pExtensionPayload);

static INT WriteEtsiAncillaryDataPayload(const AAC_METADATA* const pMetadata,
                                         UCHAR* const pExtensionPayload);

static FDK_METADATA_ERROR CompensateAudioDelay(
    HANDLE_FDK_METADATA_ENCODER hMetaDataEnc, INT_PCM* const pAudioSamples,
    const UINT audioSamplesBufSize, const INT nAudioSamples);

static FDK_METADATA_ERROR LoadSubmittedMetadata(
    const AACENC_MetaData* const hMetadata, const INT nChannels,
    const INT metadataMode, AAC_METADATA* const pAacMetaData);

static FDK_METADATA_ERROR ProcessCompressor(AAC_METADATA* pMetadata,
                                            HDRC_COMP hDrcComp,
                                            const INT_PCM* const pSamples,
                                            const UINT samplesBufSize,
                                            const INT nSamples);

static DRC_PROFILE convertProfile(AACENC_METADATA_DRC_PROFILE aacProfile) {
  DRC_PROFILE drcProfile = DRC_NONE;

  switch (aacProfile) {
    case AACENC_METADATA_DRC_NONE:
      drcProfile = DRC_NONE;
      break;
    case AACENC_METADATA_DRC_FILMSTANDARD:
      drcProfile = DRC_FILMSTANDARD;
      break;
    case AACENC_METADATA_DRC_FILMLIGHT:
      drcProfile = DRC_FILMLIGHT;
      break;
    case AACENC_METADATA_DRC_MUSICSTANDARD:
      drcProfile = DRC_MUSICSTANDARD;
      break;
    case AACENC_METADATA_DRC_MUSICLIGHT:
      drcProfile = DRC_MUSICLIGHT;
      break;
    case AACENC_METADATA_DRC_SPEECH:
      drcProfile = DRC_SPEECH;
      break;
    case AACENC_METADATA_DRC_NOT_PRESENT:
      drcProfile = DRC_NOT_PRESENT;
      break;
    default:
      drcProfile = DRC_NONE;
      break;
  }
  return drcProfile;
}

static UCHAR dialnorm2progreflvl(const INT d) {
  return ((UCHAR)fMax(0, fMin((-d + (1 << 13)) >> 14, 127)));
}

static INT progreflvl2dialnorm(const UCHAR p) {
  return -((INT)(p << (16 - 2)));
}

static SCHAR encodeDmxLvls(const SCHAR cmixlev, const SCHAR surmixlev) {
  SCHAR dmxLvls = 0;
  dmxLvls |= 0x80 | (cmixlev << 4);
  dmxLvls |= 0x08 | surmixlev;

  return dmxLvls;
}

static void encodeDynrng(INT gain, UCHAR* const dyn_rng_ctl,
                         UCHAR* const dyn_rng_sgn) {
  if (gain < 0) {
    *dyn_rng_sgn = 1;
    gain = -gain;
  } else {
    *dyn_rng_sgn = 0;
  }
  gain = fMin(gain, (127 << 14));

  *dyn_rng_ctl = (UCHAR)((gain + (1 << 13)) >> 14);
}

static INT decodeDynrng(const UCHAR dyn_rng_ctl, const UCHAR dyn_rng_sgn) {
  INT tmp = ((INT)dyn_rng_ctl << (16 - 2));
  if (dyn_rng_sgn) tmp = -tmp;

  return tmp;
}

static UCHAR encodeCompr(INT gain) {
  UCHAR x, y;
  INT tmp;

  tmp = ((3156476 - gain) * 15 + 197283) / 394566;

  if (tmp >= 240) {
    return 0xFF;
  } else if (tmp < 0) {
    return 0;
  } else {
    x = tmp / 15;
    y = tmp % 15;
  }

  return (x << 4) | y;
}

static INT decodeCompr(const UCHAR compr) {
  INT gain;
  SCHAR x = compr >> 4;
  UCHAR y = (compr & 0x0F);

  gain = (INT)(
      scaleValue((FIXP_DBL)(((LONG)FL2FXCONST_DBL(6.0206f / 128.f) * (8 - x) -
                             (LONG)FL2FXCONST_DBL(0.4014f / 128.f) * y)),
                 -(DFRACT_BITS - 1 - 7 - 16)));

  return gain;
}

FDK_METADATA_ERROR FDK_MetadataEnc_Open(HANDLE_FDK_METADATA_ENCODER* phMetaData,
                                        const UINT maxChannels) {
  FDK_METADATA_ERROR err = METADATA_OK;
  HANDLE_FDK_METADATA_ENCODER hMetaData = NULL;

  if (phMetaData == NULL) {
    err = METADATA_INVALID_HANDLE;
    goto bail;
  }

  if (NULL == (hMetaData = (HANDLE_FDK_METADATA_ENCODER)FDKcalloc(
                   1, sizeof(FDK_METADATA_ENCODER)))) {
    err = METADATA_MEMORY_ERROR;
    goto bail;
  }
  FDKmemclear(hMetaData, sizeof(FDK_METADATA_ENCODER));

  if (NULL == (hMetaData->pAudioDelayBuffer = (INT_PCM*)FDKcalloc(
                   maxChannels * MAX_DRC_FRAMELEN, sizeof(INT_PCM)))) {
    err = METADATA_MEMORY_ERROR;
    goto bail;
  }
  FDKmemclear(hMetaData->pAudioDelayBuffer,
              maxChannels * MAX_DRC_FRAMELEN * sizeof(INT_PCM));
  hMetaData->maxChannels = maxChannels;

  if (FDK_DRC_Generator_Open(&hMetaData->hDrcComp) != 0) {
    err = METADATA_MEMORY_ERROR;
    goto bail;
  }
  hMetaData->channelMode = MODE_UNKNOWN;

  *phMetaData = hMetaData;

  return err;

bail:
  FDK_MetadataEnc_Close(&hMetaData);
  return err;
}

FDK_METADATA_ERROR FDK_MetadataEnc_Close(
    HANDLE_FDK_METADATA_ENCODER* phMetaData) {
  FDK_METADATA_ERROR err = METADATA_OK;

  if (phMetaData == NULL) {
    err = METADATA_INVALID_HANDLE;
    goto bail;
  }

  if (*phMetaData != NULL) {
    FDK_DRC_Generator_Close(&(*phMetaData)->hDrcComp);
    FDKfree((*phMetaData)->pAudioDelayBuffer);
    FDKfree(*phMetaData);
    *phMetaData = NULL;
  }
bail:
  return err;
}

FDK_METADATA_ERROR FDK_MetadataEnc_Init(
    HANDLE_FDK_METADATA_ENCODER hMetaData, const INT resetStates,
    const INT metadataMode, const INT audioDelay, const UINT frameLength,
    const UINT sampleRate, const UINT nChannels, const CHANNEL_MODE channelMode,
    const CHANNEL_ORDER channelOrder) {
  FDK_METADATA_ERROR err = METADATA_OK;
  int i, nFrames, delay;

  if (hMetaData == NULL) {
    err = METADATA_INVALID_HANDLE;
    goto bail;
  }

  for (nFrames = 0, delay = audioDelay - (INT)frameLength; delay > 0;
       delay -= (INT)frameLength, nFrames++)
    ;

  if ((nChannels > (8)) || (nChannels > hMetaData->maxChannels) ||
      ((-delay) > MAX_DRC_FRAMELEN) || nFrames >= MAX_DELAY_FRAMES) {
    err = METADATA_INIT_ERROR;
    goto bail;
  }

  FDKmemcpy(&hMetaData->submittedMetaData, &defaultMetaDataSetup,
            sizeof(AACENC_MetaData));

  hMetaData->finalizeMetaData =
      0;
  hMetaData->initializeMetaData =
      0;

  if (resetStates || (hMetaData->nAudioDataDelay != -delay) ||
      (hMetaData->channelMode != channelMode)) {
    if (resetStates || (hMetaData->channelMode == MODE_UNKNOWN)) {
      FDKmemclear(hMetaData->pAudioDelayBuffer,
                  hMetaData->maxChannels * MAX_DRC_FRAMELEN * sizeof(INT_PCM));
    } else {
      FDK_channelMapDescr mapDescrPrev, mapDescr;
      int c, src[2] = {-1, -1}, dst[2] = {-1, -1};

      if (channelOrder == CH_ORDER_WG4) {
        FDK_chMapDescr_init(&mapDescrPrev, FDK_mapInfoTabWg4,
                            FDK_mapInfoTabLenWg4, 0);
        FDK_chMapDescr_init(&mapDescr, FDK_mapInfoTabWg4,
                            FDK_mapInfoTabLenWg4, 0);
      } else {
        FDK_chMapDescr_init(&mapDescrPrev, NULL, 0,
                            (channelOrder == CH_ORDER_MPEG) ? 1 : 0);
        FDK_chMapDescr_init(&mapDescr, NULL, 0,
                            (channelOrder == CH_ORDER_MPEG) ? 1 : 0);
      }

      switch (channelMode) {
        case MODE_1:
          if ((INT)nChannels != 2) {
            src[0] = FDK_chMapDescr_getMapValue(&mapDescrPrev, 0,
                                                hMetaData->channelMode);
            dst[0] = FDK_chMapDescr_getMapValue(&mapDescr, 0, channelMode);
          }
          break;
        case MODE_2:
        case MODE_1_2:
        case MODE_1_2_1:
        case MODE_1_2_2:
        case MODE_1_2_2_1:
          if (hMetaData->nChannels >= 2) {
            src[0] = FDK_chMapDescr_getMapValue(
                &mapDescrPrev, ((hMetaData->channelMode == 2) ? 0 : 1),
                hMetaData->channelMode);
            src[1] = FDK_chMapDescr_getMapValue(
                &mapDescrPrev, ((hMetaData->channelMode == 2) ? 1 : 2),
                hMetaData->channelMode);
            dst[0] = FDK_chMapDescr_getMapValue(
                &mapDescr, ((channelMode == 2) ? 0 : 1), channelMode);
            dst[1] = FDK_chMapDescr_getMapValue(
                &mapDescr, ((channelMode == 2) ? 1 : 2), channelMode);
          }
          break;
        default:;
      }
      C_ALLOC_SCRATCH_START(scratch_audioDelayBuffer, INT_PCM, (8));
      FDKmemclear(scratch_audioDelayBuffer, (8) * sizeof(INT_PCM));

      i = (hMetaData->nChannels > (INT)nChannels)
              ? 0
              : hMetaData->nAudioDataDelay - 1;
      do {
        for (c = 0; c < 2; c++) {
          if (src[c] != -1 && dst[c] != -1) {
            scratch_audioDelayBuffer[dst[c]] =
                hMetaData->pAudioDelayBuffer[i * hMetaData->nChannels + src[c]];
          }
        }
        FDKmemcpy(&hMetaData->pAudioDelayBuffer[i * nChannels],
                  scratch_audioDelayBuffer, nChannels * sizeof(INT_PCM));
        i += (hMetaData->nChannels > (INT)nChannels) ? 1 : -1;
      } while ((i < hMetaData->nAudioDataDelay) && (i >= 0));

      C_ALLOC_SCRATCH_END(scratch_audioDelayBuffer, INT_PCM, (8));
    }
    FDKmemclear(hMetaData->metaDataBuffer, sizeof(hMetaData->metaDataBuffer));
    hMetaData->metaDataDelayIdx = 0;
    hMetaData->initializeMetaData =
        1;
  } else {
    if ((hMetaData->metadataMode == 0) && (metadataMode != 0)) {
      for (i = 0;
           i < (int)(sizeof(hMetaData->metaDataBuffer) / sizeof(AAC_METADATA));
           i++) {
        LoadSubmittedMetadata(&hMetaData->submittedMetaData, nChannels, 0,
                              &hMetaData->metaDataBuffer[i]);
      }
    }

    if ((hMetaData->metadataMode != 0) && (metadataMode == 0)) {
      hMetaData->finalizeMetaData = hMetaData->metadataMode;
    }
  }

  hMetaData->nAudioDataDelay = -delay;
  hMetaData->nMetaDataDelay = nFrames;
  hMetaData->nChannels = nChannels;
  hMetaData->channelMode = channelMode;
  hMetaData->metadataMode = metadataMode;

  if ((metadataMode == 1) || (metadataMode == 2)) {
    if (FDK_DRC_Generator_Initialize(hMetaData->hDrcComp, DRC_NONE, DRC_NONE,
                                     frameLength, sampleRate, channelMode,
                                     channelOrder, 1) != 0) {
      err = METADATA_INIT_ERROR;
    }
  }
bail:
  return err;
}

static FDK_METADATA_ERROR ProcessCompressor(AAC_METADATA* pMetadata,
                                            HDRC_COMP hDrcComp,
                                            const INT_PCM* const pSamples,
                                            const UINT samplesBufSize,
                                            const INT nSamples) {
  FDK_METADATA_ERROR err = METADATA_OK;

  INT dynrng, compr;
  INT dmxGain5, dmxGain2;
  DRC_PROFILE profileDrc;
  DRC_PROFILE profileComp;

  if ((pMetadata == NULL) || (hDrcComp == NULL)) {
    err = METADATA_INVALID_HANDLE;
    return err;
  }

  profileDrc = convertProfile(pMetadata->mpegDrc.drc_profile);
  profileComp = convertProfile(pMetadata->etsiAncData.comp_profile);

  if ((profileDrc != FDK_DRC_Generator_getDrcProfile(hDrcComp)) ||
      (profileComp != FDK_DRC_Generator_getCompProfile(hDrcComp))) {
    FDK_DRC_Generator_setDrcProfile(hDrcComp, profileDrc, profileComp);
  }

  if (profileComp == DRC_NONE) {
    pMetadata->etsiAncData.compression_value = 0x80;
  }

  dynrng = decodeDynrng(pMetadata->mpegDrc.dyn_rng_ctl[0],
                        pMetadata->mpegDrc.dyn_rng_sgn[0]);
  compr = decodeCompr(pMetadata->etsiAncData.compression_value);

  dmxGain5 = decodeDynrng(pMetadata->etsiAncData.extAncData.dmx_gain_5_idx,
                          pMetadata->etsiAncData.extAncData.dmx_gain_5_sgn);
  dmxGain2 = decodeDynrng(pMetadata->etsiAncData.extAncData.dmx_gain_2_idx,
                          pMetadata->etsiAncData.extAncData.dmx_gain_2_sgn);

  if (FDK_DRC_Generator_Calc(
          hDrcComp, pSamples, samplesBufSize,
          progreflvl2dialnorm(pMetadata->mpegDrc.prog_ref_level),
          pMetadata->mpegDrc.drc_TargetRefLevel,
          pMetadata->etsiAncData.comp_TargetRefLevel,
          dmxTable[pMetadata->centerMixLevel],
          dmxTable[pMetadata->surroundMixLevel],
          dmxTable[pMetadata->etsiAncData.extAncData.ext_dmix_a_idx],
          dmxTable[pMetadata->etsiAncData.extAncData.ext_dmix_b_idx],
          pMetadata->etsiAncData.extAncData.ext_lfe_downmix_status
              ? dmxLfeTable[pMetadata->etsiAncData.extAncData.ext_dmix_lfe_idx]
              : (FIXP_DBL)0,
          dmxGain5, dmxGain2, &dynrng, &compr) != 0) {
    err = METADATA_ENCODE_ERROR;
    goto bail;
  }

  pMetadata->mpegDrc.drc_band_incr = 0;
  encodeDynrng(dynrng, pMetadata->mpegDrc.dyn_rng_ctl,
               pMetadata->mpegDrc.dyn_rng_sgn);
  pMetadata->etsiAncData.compression_value = encodeCompr(compr);

bail:
  return err;
}

FDK_METADATA_ERROR FDK_MetadataEnc_Process(
    HANDLE_FDK_METADATA_ENCODER hMetaDataEnc, INT_PCM* const pAudioSamples,
    const UINT audioSamplesBufSize, const INT nAudioSamples,
    const AACENC_MetaData* const pMetadata,
    AACENC_EXT_PAYLOAD** ppMetaDataExtPayload, UINT* nMetaDataExtensions,
    INT* matrix_mixdown_idx) {
  FDK_METADATA_ERROR err = METADATA_OK;
  int metaDataDelayWriteIdx, metaDataDelayReadIdx, metadataMode;

  metaDataDelayWriteIdx = hMetaDataEnc->metaDataDelayIdx;
  metadataMode = hMetaDataEnc->metadataMode;

  hMetaDataEnc->metaDataDelayIdx++;
  if (hMetaDataEnc->metaDataDelayIdx > hMetaDataEnc->nMetaDataDelay)
    hMetaDataEnc->metaDataDelayIdx = 0;

  metaDataDelayReadIdx = hMetaDataEnc->metaDataDelayIdx;

  if (pMetadata != NULL) {
    FDKmemcpy(&hMetaDataEnc->submittedMetaData, pMetadata,
              sizeof(AACENC_MetaData));
  }

  if ((hMetaDataEnc->finalizeMetaData != 0) &&
      (hMetaDataEnc->metadataMode == 0)) {
    FDKmemcpy(&hMetaDataEnc->submittedMetaData, &defaultMetaDataSetup,
              sizeof(AACENC_MetaData));
    metadataMode = hMetaDataEnc->finalizeMetaData;
    hMetaDataEnc->finalizeMetaData = 0;
  }

  if ((err = LoadSubmittedMetadata(
           &hMetaDataEnc->submittedMetaData, hMetaDataEnc->nChannels,
           metadataMode,
           &hMetaDataEnc->metaDataBuffer[metaDataDelayWriteIdx])) !=
      METADATA_OK) {
    goto bail;
  }

  if ((hMetaDataEnc->metaDataBuffer[metaDataDelayWriteIdx].metadataMode == 1) ||
      (hMetaDataEnc->metaDataBuffer[metaDataDelayWriteIdx].metadataMode == 2)) {
    if ((err = ProcessCompressor(
             &hMetaDataEnc->metaDataBuffer[metaDataDelayWriteIdx],
             hMetaDataEnc->hDrcComp, pAudioSamples, audioSamplesBufSize,
             nAudioSamples)) != METADATA_OK) {
      LoadSubmittedMetadata(
          &hMetaDataEnc->submittedMetaData, hMetaDataEnc->nChannels,
          metadataMode, &hMetaDataEnc->metaDataBuffer[metaDataDelayWriteIdx]);
    }
  }

  if ((hMetaDataEnc->initializeMetaData != 0) &&
      (hMetaDataEnc->metadataMode != 0)) {
    int i;
    for (i = 0;
         i < (int)(sizeof(hMetaDataEnc->metaDataBuffer) / sizeof(AAC_METADATA));
         i++) {
      if (i != metaDataDelayWriteIdx) {
        FDKmemcpy(&hMetaDataEnc->metaDataBuffer[i],
                  &hMetaDataEnc->metaDataBuffer[metaDataDelayWriteIdx],
                  sizeof(hMetaDataEnc->metaDataBuffer[metaDataDelayWriteIdx]));
      }
    }
    hMetaDataEnc->initializeMetaData = 0;
  }

  FDK_ASSERT(metaDataDelayReadIdx < MAX_DELAY_FRAMES);
  if ((err = WriteMetadataPayload(
           hMetaDataEnc,
           &hMetaDataEnc->metaDataBuffer[metaDataDelayReadIdx])) !=
      METADATA_OK) {
    goto bail;
  }

  *ppMetaDataExtPayload = hMetaDataEnc->exPayload;
  *nMetaDataExtensions = hMetaDataEnc->nExtensions;
  *matrix_mixdown_idx = hMetaDataEnc->matrix_mixdown_idx;

bail:
  err = CompensateAudioDelay(hMetaDataEnc, pAudioSamples, audioSamplesBufSize,
                             nAudioSamples / hMetaDataEnc->nChannels);

  return err;
}

static FDK_METADATA_ERROR CompensateAudioDelay(
    HANDLE_FDK_METADATA_ENCODER hMetaDataEnc, INT_PCM* const pAudioSamples,
    const UINT audioSamplesBufSize, const INT nAudioSamples) {
  FDK_METADATA_ERROR err = METADATA_OK;

  if (hMetaDataEnc->nAudioDataDelay) {
    C_ALLOC_SCRATCH_START(scratch_audioDelayBuffer, INT_PCM, 1024);

    for (int c = 0; c < hMetaDataEnc->nChannels; c++) {
      int M = 1024;
      INT_PCM* pAudioSamples2 = pAudioSamples + c * audioSamplesBufSize;
      int delayIdx = hMetaDataEnc->nAudioDataDelay;

      do {
        M = fMin(M, delayIdx);
        delayIdx -= M;

        FDKmemcpy(&scratch_audioDelayBuffer[0],
                  &pAudioSamples2[(nAudioSamples - M)], sizeof(INT_PCM) * M);
        FDKmemmove(&pAudioSamples2[M], &pAudioSamples2[0],
                   sizeof(INT_PCM) * (nAudioSamples - M));
        FDKmemcpy(
            &pAudioSamples2[0],
            &hMetaDataEnc->pAudioDelayBuffer[delayIdx +
                                             c * hMetaDataEnc->nAudioDataDelay],
            sizeof(INT_PCM) * M);
        FDKmemcpy(
            &hMetaDataEnc->pAudioDelayBuffer[delayIdx +
                                             c * hMetaDataEnc->nAudioDataDelay],
            &scratch_audioDelayBuffer[0], sizeof(INT_PCM) * M);

      } while (delayIdx > 0);
    }

    C_ALLOC_SCRATCH_END(scratch_audioDelayBuffer, INT_PCM, 1024);
  }

  return err;
}

static FDK_METADATA_ERROR WriteMetadataPayload(
    const HANDLE_FDK_METADATA_ENCODER hMetaData,
    const AAC_METADATA* const pMetadata) {
  FDK_METADATA_ERROR err = METADATA_OK;

  if ((hMetaData == NULL) || (pMetadata == NULL)) {
    err = METADATA_INVALID_HANDLE;
    goto bail;
  }

  hMetaData->nExtensions = 0;
  hMetaData->matrix_mixdown_idx = -1;

  if (pMetadata->metadataMode != 0) {
    if ((pMetadata->metadataMode == 1) || (pMetadata->metadataMode == 2)) {
      hMetaData->exPayload[hMetaData->nExtensions].pData =
          hMetaData->drcInfoPayload;
      hMetaData->exPayload[hMetaData->nExtensions].dataType = EXT_DYNAMIC_RANGE;
      hMetaData->exPayload[hMetaData->nExtensions].associatedChElement = -1;

      hMetaData->exPayload[hMetaData->nExtensions].dataSize =
          WriteDynamicRangeInfoPayload(
              pMetadata, hMetaData->exPayload[hMetaData->nExtensions].pData);

      hMetaData->nExtensions++;
    }

    if (pMetadata->WritePCEMixDwnIdx) {
      hMetaData->matrix_mixdown_idx =
          surmix2matrix_mixdown_idx[pMetadata->surroundMixLevel];
    }

    if ((pMetadata->metadataMode == 2) ||
        (pMetadata->metadataMode == 3))
    {
      hMetaData->exPayload[hMetaData->nExtensions].pData =
          hMetaData->drcDsePayload;
      hMetaData->exPayload[hMetaData->nExtensions].dataType = EXT_DATA_ELEMENT;
      hMetaData->exPayload[hMetaData->nExtensions].associatedChElement = -1;

      hMetaData->exPayload[hMetaData->nExtensions].dataSize =
          WriteEtsiAncillaryDataPayload(
              pMetadata, hMetaData->exPayload[hMetaData->nExtensions].pData);

      hMetaData->nExtensions++;
    }

  }

bail:
  return err;
}

static INT WriteDynamicRangeInfoPayload(const AAC_METADATA* const pMetadata,
                                        UCHAR* const pExtensionPayload) {
  const INT pce_tag_present = 0;
  const INT prog_ref_lev_res_bits = 0;
  INT i, drc_num_bands = 1;

  FDK_BITSTREAM bsWriter;
  FDKinitBitStream(&bsWriter, pExtensionPayload, 16, 0, BS_WRITER);

  FDKwriteBits(&bsWriter, pce_tag_present, 1);
  if (pce_tag_present) {
    FDKwriteBits(&bsWriter, 0x0, 4);
    FDKwriteBits(&bsWriter, 0x0, 4);
  }

  FDKwriteBits(&bsWriter, (pMetadata->mpegDrc.excluded_chns_present) ? 1 : 0,
               1);

  FDKwriteBits(&bsWriter, (pMetadata->mpegDrc.drc_bands_present) ? 1 : 0,
               1);
  if (pMetadata->mpegDrc.drc_bands_present) {
    FDKwriteBits(&bsWriter, pMetadata->mpegDrc.drc_band_incr,
                 4);
    FDKwriteBits(&bsWriter, pMetadata->mpegDrc.drc_interpolation_scheme,
                 4);
    drc_num_bands += pMetadata->mpegDrc.drc_band_incr;
    for (i = 0; i < drc_num_bands; i++) {
      FDKwriteBits(&bsWriter, pMetadata->mpegDrc.drc_band_top[i],
                   8);
    }
  }

  FDKwriteBits(&bsWriter, pMetadata->mpegDrc.prog_ref_level_present,
               1);
  if (pMetadata->mpegDrc.prog_ref_level_present) {
    FDKwriteBits(&bsWriter, pMetadata->mpegDrc.prog_ref_level,
                 7);
    FDKwriteBits(&bsWriter, prog_ref_lev_res_bits,
                 1);
  }

  for (i = 0; i < drc_num_bands; i++) {
    FDKwriteBits(&bsWriter, (pMetadata->mpegDrc.dyn_rng_sgn[i]) ? 1 : 0,
                 1);
    FDKwriteBits(&bsWriter, pMetadata->mpegDrc.dyn_rng_ctl[i],
                 7);
  }

  return FDKgetValidBits(&bsWriter);
}

static INT WriteEtsiAncillaryDataPayload(const AAC_METADATA* const pMetadata,
                                         UCHAR* const pExtensionPayload) {
  FDK_BITSTREAM bsWriter;
  FDKinitBitStream(&bsWriter, pExtensionPayload, 16, 0, BS_WRITER);

  FDKwriteBits(&bsWriter, 0xBC, 8);

  FDKwriteBits(&bsWriter, 0x3, 2);
  FDKwriteBits(&bsWriter, pMetadata->dolbySurroundMode,
               2);
  FDKwriteBits(&bsWriter, pMetadata->drcPresentationMode,
               2);
  FDKwriteBits(&bsWriter, 0x0, 1);
  FDKwriteBits(&bsWriter, 0x0, 1);

  FDKwriteBits(&bsWriter, 0, 3);
  FDKwriteBits(&bsWriter, (pMetadata->DmxLvl_On) ? 1 : 0,
               1);
  FDKwriteBits(&bsWriter, pMetadata->etsiAncData.extAncDataStatus,
               1);
  FDKwriteBits(&bsWriter, (pMetadata->etsiAncData.compression_on) ? 1 : 0,
               1);
  FDKwriteBits(&bsWriter,
               (pMetadata->etsiAncData.timecode_coarse_status) ? 1 : 0,
               1);
  FDKwriteBits(&bsWriter, (pMetadata->etsiAncData.timecode_fine_status) ? 1 : 0,
               1);

  if (pMetadata->DmxLvl_On) {
    FDKwriteBits(
        &bsWriter,
        encodeDmxLvls(pMetadata->centerMixLevel, pMetadata->surroundMixLevel),
        8);
  }

  if (pMetadata->etsiAncData.compression_on) {
    FDKwriteBits(&bsWriter, 0x01, 8);
    FDKwriteBits(&bsWriter, pMetadata->etsiAncData.compression_value,
                 8);
  }

  if (pMetadata->etsiAncData.timecode_coarse_status) {
    FDKwriteBits(&bsWriter, 0x0, 16);
  }

  if (pMetadata->etsiAncData.timecode_fine_status) {
    FDKwriteBits(&bsWriter, 0x0, 16);
  }

  if (pMetadata->etsiAncData.extAncDataStatus) {
    FDKwriteBits(&bsWriter, 0, 1);
    FDKwriteBits(&bsWriter,
                 pMetadata->etsiAncData.extAncData.ext_downmix_lvl_status,
                 1);
    FDKwriteBits(&bsWriter,
                 pMetadata->etsiAncData.extAncData.ext_downmix_gain_status,
                 1);
    FDKwriteBits(&bsWriter,
                 pMetadata->etsiAncData.extAncData.ext_lfe_downmix_status,
                 1);
    FDKwriteBits(&bsWriter, 0, 4);

    if (pMetadata->etsiAncData.extAncData.ext_downmix_lvl_status) {
      FDKwriteBits(&bsWriter, pMetadata->etsiAncData.extAncData.ext_dmix_a_idx,
                   3);
      FDKwriteBits(&bsWriter, pMetadata->etsiAncData.extAncData.ext_dmix_b_idx,
                   3);
      FDKwriteBits(&bsWriter, 0, 2);
    }

    if (pMetadata->etsiAncData.extAncData.ext_downmix_gain_status) {
      FDKwriteBits(&bsWriter, pMetadata->etsiAncData.extAncData.dmx_gain_5_sgn,
                   1);
      FDKwriteBits(&bsWriter, pMetadata->etsiAncData.extAncData.dmx_gain_5_idx,
                   6);
      FDKwriteBits(&bsWriter, 0, 1);
      FDKwriteBits(&bsWriter, pMetadata->etsiAncData.extAncData.dmx_gain_2_sgn,
                   1);
      FDKwriteBits(&bsWriter, pMetadata->etsiAncData.extAncData.dmx_gain_2_idx,
                   6);
      FDKwriteBits(&bsWriter, 0, 1);
    }

    if (pMetadata->etsiAncData.extAncData.ext_lfe_downmix_status) {
      FDKwriteBits(&bsWriter,
                   pMetadata->etsiAncData.extAncData.ext_dmix_lfe_idx,
                   4);
      FDKwriteBits(&bsWriter, 0, 4);
    }
  }

  return FDKgetValidBits(&bsWriter);
}

static FDK_METADATA_ERROR LoadSubmittedMetadata(
    const AACENC_MetaData* const hMetadata, const INT nChannels,
    const INT metadataMode, AAC_METADATA* const pAacMetaData) {
  FDK_METADATA_ERROR err = METADATA_OK;

  if (pAacMetaData == NULL) {
    err = METADATA_INVALID_HANDLE;
  } else {
    FDKmemclear(pAacMetaData, sizeof(AAC_METADATA));

    if (hMetadata != NULL) {
      pAacMetaData->mpegDrc.drc_profile = hMetadata->drc_profile;
      pAacMetaData->etsiAncData.comp_profile = hMetadata->comp_profile;
      pAacMetaData->mpegDrc.drc_TargetRefLevel = hMetadata->drc_TargetRefLevel;
      pAacMetaData->etsiAncData.comp_TargetRefLevel =
          hMetadata->comp_TargetRefLevel;
      pAacMetaData->mpegDrc.prog_ref_level_present =
          hMetadata->prog_ref_level_present;
      pAacMetaData->mpegDrc.prog_ref_level =
          dialnorm2progreflvl(hMetadata->prog_ref_level);

      pAacMetaData->centerMixLevel = hMetadata->centerMixLevel;
      pAacMetaData->surroundMixLevel = hMetadata->surroundMixLevel;
      pAacMetaData->WritePCEMixDwnIdx = hMetadata->PCE_mixdown_idx_present;
      pAacMetaData->DmxLvl_On = hMetadata->ETSI_DmxLvl_present;

      pAacMetaData->etsiAncData.compression_on =
          (hMetadata->comp_profile == AACENC_METADATA_DRC_NOT_PRESENT ? 0 : 1);

      if (pAacMetaData->mpegDrc.drc_profile ==
          AACENC_METADATA_DRC_NOT_PRESENT) {
        pAacMetaData->mpegDrc.drc_profile =
            AACENC_METADATA_DRC_NONE;
      }

      if (nChannels == 2) {
        pAacMetaData->dolbySurroundMode =
            hMetadata->dolbySurroundMode;
      } else {
        pAacMetaData->dolbySurroundMode = 0;
      }

      pAacMetaData->drcPresentationMode = hMetadata->drcPresentationMode;
      if (pAacMetaData->drcPresentationMode == 1) {
        pAacMetaData->mpegDrc.drc_TargetRefLevel =
            fMax(-(31 << 16), pAacMetaData->mpegDrc.drc_TargetRefLevel);
        pAacMetaData->etsiAncData.comp_TargetRefLevel = fMax(
            -(20 << 16),
            pAacMetaData->etsiAncData.comp_TargetRefLevel);
      }
      if (pAacMetaData->drcPresentationMode == 2) {
        pAacMetaData->mpegDrc.drc_TargetRefLevel =
            fMax(-(23 << 16), pAacMetaData->mpegDrc.drc_TargetRefLevel);
        pAacMetaData->etsiAncData.comp_TargetRefLevel =
            fMax(-(23 << 16), pAacMetaData->etsiAncData.comp_TargetRefLevel);
      }
      if (pAacMetaData->etsiAncData.comp_profile ==
          AACENC_METADATA_DRC_NOT_PRESENT) {
        if (pAacMetaData->drcPresentationMode != 0) {
          pAacMetaData->mpegDrc.drc_TargetRefLevel =
              fMax(pAacMetaData->etsiAncData.comp_TargetRefLevel,
                   pAacMetaData->mpegDrc.drc_TargetRefLevel);
        }
      }

      pAacMetaData->etsiAncData.timecode_coarse_status =
          0;
      pAacMetaData->etsiAncData.timecode_fine_status =
          0;

      pAacMetaData->etsiAncData.extAncDataStatus =
          ((hMetadata->ExtMetaData.extAncDataEnable == 1) ? 1 : 0);

      if (pAacMetaData->etsiAncData.extAncDataStatus) {
        pAacMetaData->etsiAncData.extAncData.ext_downmix_lvl_status =
            (hMetadata->ExtMetaData.extDownmixLevelEnable ? 1 : 0);
        pAacMetaData->etsiAncData.extAncData.ext_downmix_gain_status =
            (hMetadata->ExtMetaData.dmxGainEnable ? 1 : 0);
        pAacMetaData->etsiAncData.extAncData.ext_lfe_downmix_status =
            (hMetadata->ExtMetaData.lfeDmxEnable ? 1 : 0);

        pAacMetaData->etsiAncData.extAncData.ext_dmix_a_idx =
            hMetadata->ExtMetaData.extDownmixLevel_A;
        pAacMetaData->etsiAncData.extAncData.ext_dmix_b_idx =
            hMetadata->ExtMetaData.extDownmixLevel_B;

        if (pAacMetaData->etsiAncData.extAncData.ext_downmix_gain_status) {
          encodeDynrng(hMetadata->ExtMetaData.dmxGain5,
                       &(pAacMetaData->etsiAncData.extAncData.dmx_gain_5_idx),
                       &(pAacMetaData->etsiAncData.extAncData.dmx_gain_5_sgn));
          encodeDynrng(hMetadata->ExtMetaData.dmxGain2,
                       &(pAacMetaData->etsiAncData.extAncData.dmx_gain_2_idx),
                       &(pAacMetaData->etsiAncData.extAncData.dmx_gain_2_sgn));
        } else {
          encodeDynrng(1 << 16,
                       &(pAacMetaData->etsiAncData.extAncData.dmx_gain_5_idx),
                       &(pAacMetaData->etsiAncData.extAncData.dmx_gain_5_sgn));
          encodeDynrng(1 << 16,
                       &(pAacMetaData->etsiAncData.extAncData.dmx_gain_2_idx),
                       &(pAacMetaData->etsiAncData.extAncData.dmx_gain_2_sgn));
        }

        if (pAacMetaData->etsiAncData.extAncData.ext_lfe_downmix_status) {
          pAacMetaData->etsiAncData.extAncData.ext_dmix_lfe_idx =
              hMetadata->ExtMetaData.lfeDmxLevel;
        } else {
          pAacMetaData->etsiAncData.extAncData.ext_dmix_lfe_idx =
              15;
        }
      } else {
        pAacMetaData->etsiAncData.extAncData.ext_downmix_lvl_status = 0;
        pAacMetaData->etsiAncData.extAncData.ext_downmix_gain_status = 0;
        pAacMetaData->etsiAncData.extAncData.ext_lfe_downmix_status = 0;

        pAacMetaData->etsiAncData.extAncData.ext_dmix_a_idx = 7;
        pAacMetaData->etsiAncData.extAncData.ext_dmix_b_idx = 7;

        encodeDynrng(1 << 16,
                     &(pAacMetaData->etsiAncData.extAncData.dmx_gain_5_idx),
                     &(pAacMetaData->etsiAncData.extAncData.dmx_gain_5_sgn));
        encodeDynrng(1 << 16,
                     &(pAacMetaData->etsiAncData.extAncData.dmx_gain_2_idx),
                     &(pAacMetaData->etsiAncData.extAncData.dmx_gain_2_sgn));

        pAacMetaData->etsiAncData.extAncData.ext_dmix_lfe_idx =
            15;
      }

      pAacMetaData->metadataMode = metadataMode;
    } else {
      pAacMetaData->metadataMode = 0;
    }
  }

  return err;
}

INT FDK_MetadataEnc_GetDelay(HANDLE_FDK_METADATA_ENCODER hMetadataEnc) {
  INT delay = 0;

  if (hMetadataEnc != NULL) {
    delay = hMetadataEnc->nAudioDataDelay;
  }

  return delay;
}
