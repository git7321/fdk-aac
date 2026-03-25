/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "env_extr.h"
#include "sbr_ram.h"
#include "sbr_rom.h"
#include "huff_dec.h"

#include "psbitdec.h"

#define DRM_PARAMETRIC_STEREO 0
#define EXTENSION_ID_PS_CODING 2

static int extractPvcFrameInfo(
    HANDLE_FDK_BITSTREAM hBs,
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_frame_data,
    HANDLE_SBR_PREV_FRAME_DATA h_prev_frame_data,
    UCHAR pvc_mode_last,
    const UINT flags);
static int extractFrameInfo(HANDLE_FDK_BITSTREAM hBs,
                            HANDLE_SBR_HEADER_DATA hHeaderData,
                            HANDLE_SBR_FRAME_DATA h_frame_data,
                            const UINT nrOfChannels, const UINT flags);

static int sbrGetPvcEnvelope(HANDLE_SBR_HEADER_DATA hHeaderData,
                             HANDLE_SBR_FRAME_DATA h_frame_data,
                             HANDLE_FDK_BITSTREAM hBs, const UINT flags,
                             const UINT pvcMode);
static int sbrGetEnvelope(HANDLE_SBR_HEADER_DATA hHeaderData,
                          HANDLE_SBR_FRAME_DATA h_frame_data,
                          HANDLE_FDK_BITSTREAM hBs, const UINT flags);

static void sbrGetDirectionControlData(HANDLE_SBR_FRAME_DATA hFrameData,
                                       HANDLE_FDK_BITSTREAM hBs,
                                       const UINT flags, const int bs_pvc_mode);

static void sbrGetNoiseFloorData(HANDLE_SBR_HEADER_DATA hHeaderData,
                                 HANDLE_SBR_FRAME_DATA h_frame_data,
                                 HANDLE_FDK_BITSTREAM hBs);

static int checkFrameInfo(FRAME_INFO *pFrameInfo, int numberOfTimeSlots,
                          int overlap, int timeStep);

typedef struct SR_MAPPING {
  UINT fsRangeLo;
  UINT fsMapped;
} SR_MAPPING;

static const SR_MAPPING stdSampleRatesMapping[] = {
    {0, 8000},      {9391, 11025},  {11502, 12000}, {13856, 16000},
    {18783, 22050}, {23004, 24000}, {27713, 32000}, {37566, 44100},
    {46009, 48000}, {55426, 64000}, {75132, 88200}, {92017, 96000}};
static const SR_MAPPING stdSampleRatesMappingUsac[] = {
    {0, 16000},     {18783, 22050}, {23004, 24000}, {27713, 32000},
    {35777, 40000}, {42000, 44100}, {46009, 48000}, {55426, 64000},
    {75132, 88200}, {92017, 96000}};

UINT sbrdec_mapToStdSampleRate(UINT fs,
                               UINT isUsac)
{
  UINT fsMapped = fs, tableSize = 0;
  const SR_MAPPING *mappingTable;
  int i;

  if (!isUsac) {
    mappingTable = stdSampleRatesMapping;
    tableSize = sizeof(stdSampleRatesMapping) / sizeof(SR_MAPPING);
  } else {
    mappingTable = stdSampleRatesMappingUsac;
    tableSize = sizeof(stdSampleRatesMappingUsac) / sizeof(SR_MAPPING);
  }

  for (i = tableSize - 1; i >= 0; i--) {
    if (fs >= mappingTable[i].fsRangeLo) {
      fsMapped = mappingTable[i].fsMapped;
      break;
    }
  }

  return (fsMapped);
}

SBR_ERROR
initHeaderData(HANDLE_SBR_HEADER_DATA hHeaderData, const int sampleRateIn,
               const int sampleRateOut, const INT downscaleFactor,
               const int samplesPerFrame, const UINT flags,
               const int setDefaultHdr) {
  HANDLE_FREQ_BAND_DATA hFreq = &hHeaderData->freqBandData;
  SBR_ERROR sbrError = SBRDEC_OK;
  int numAnalysisBands;
  int sampleRateProc;

  if (!(flags & (SBRDEC_SYNTAX_USAC | SBRDEC_SYNTAX_RSVD50))) {
    sampleRateProc =
        sbrdec_mapToStdSampleRate(sampleRateOut * downscaleFactor, 0);
  } else {
    sampleRateProc = sampleRateOut * downscaleFactor;
  }

  if (sampleRateIn == sampleRateOut) {
    hHeaderData->sbrProcSmplRate = sampleRateProc << 1;
    numAnalysisBands = 32;
  } else {
    hHeaderData->sbrProcSmplRate = sampleRateProc;
    if ((sampleRateOut >> 1) == sampleRateIn) {
      numAnalysisBands = 32;
    } else if ((sampleRateOut >> 2) == sampleRateIn) {
      numAnalysisBands = 16;
    } else if ((sampleRateOut * 3) >> 3 == (sampleRateIn * 8) >> 3) {
      numAnalysisBands = 24;
    } else {
      sbrError = SBRDEC_UNSUPPORTED_CONFIG;
      goto bail;
    }
  }
  numAnalysisBands /= downscaleFactor;

  if (setDefaultHdr) {
    hHeaderData->syncState = SBR_NOT_INITIALIZED;
    hHeaderData->status = 0;
    hHeaderData->frameErrorFlag = 0;

    hHeaderData->bs_info.ampResolution = 1;
    hHeaderData->bs_info.xover_band = 0;
    hHeaderData->bs_info.sbr_preprocessing = 0;
    hHeaderData->bs_info.pvc_mode = 0;

    hHeaderData->bs_data.startFreq = 5;
    hHeaderData->bs_data.stopFreq = 0;
    hHeaderData->bs_data.freqScale =
        0;
    hHeaderData->bs_data.alterScale = 1;
    hHeaderData->bs_data.noise_bands = 2;
    hHeaderData->bs_data.limiterBands = 2;
    hHeaderData->bs_data.limiterGains = 2;
    hHeaderData->bs_data.interpolFreq = 1;
    hHeaderData->bs_data.smoothingLength = 1;

    if (sampleRateOut * downscaleFactor >= 96000) {
      hHeaderData->bs_data.startFreq =
          4;
      hHeaderData->bs_data.stopFreq = 3;
    } else if (sampleRateOut * downscaleFactor >
               24000) {
      hHeaderData->bs_data.startFreq =
          7;
      hHeaderData->bs_data.stopFreq = 3;
    }
  }

  if ((sampleRateOut >> 2) == sampleRateIn) {
    hHeaderData->timeStep = 4;
  } else {
    hHeaderData->timeStep = (flags & SBRDEC_ELD_GRID) ? 1 : 2;
  }

  hFreq->freqBandTable[0] = hFreq->freqBandTableLo;
  hFreq->freqBandTable[1] = hFreq->freqBandTableHi;

  hHeaderData->numberTimeSlots =
      (samplesPerFrame / numAnalysisBands) >> (hHeaderData->timeStep - 1);
  if (hHeaderData->numberTimeSlots > (16)) {
    sbrError = SBRDEC_UNSUPPORTED_CONFIG;
  }

  hHeaderData->numberOfAnalysisBands = numAnalysisBands;
  if ((sampleRateOut >> 2) == sampleRateIn) {
    hHeaderData->numberTimeSlots <<= 1;
  }

bail:
  return sbrError;
}

void initSbrPrevFrameData(
    HANDLE_SBR_PREV_FRAME_DATA
        h_prev_data,
    int timeSlots)
{
  int i;

  for (i = 0; i < MAX_FREQ_COEFFS; i++)
    h_prev_data->sfb_nrg_prev[i] = (FIXP_DBL)0;
  for (i = 0; i < MAX_NOISE_COEFFS; i++)
    h_prev_data->prevNoiseLevel[i] = (FIXP_DBL)0;
  for (i = 0; i < MAX_INVF_BANDS; i++) h_prev_data->sbr_invf_mode[i] = INVF_OFF;

  h_prev_data->stopPos = timeSlots;
  h_prev_data->coupling = COUPLING_OFF;
  h_prev_data->ampRes = 0;

  FDKmemclear(&h_prev_data->prevFrameInfo, sizeof(h_prev_data->prevFrameInfo));
}

SBR_HEADER_STATUS
sbrGetHeaderData(HANDLE_SBR_HEADER_DATA hHeaderData, HANDLE_FDK_BITSTREAM hBs,
                 const UINT flags, const int fIsSbrData,
                 const UCHAR configMode) {
  SBR_HEADER_DATA_BS *pBsData;
  SBR_HEADER_DATA_BS lastHeader;
  SBR_HEADER_DATA_BS_INFO lastInfo;
  int headerExtra1 = 0, headerExtra2 = 0;

  if (configMode & AC_CM_DET_CFG_CHANGE) {
    if (!(flags & (SBRDEC_SYNTAX_RSVD50 | SBRDEC_SYNTAX_USAC))) {
      FDKreadBits(hBs, 1);
    }
    FDKpushFor(hBs, 8);
    if (!(flags & (SBRDEC_SYNTAX_RSVD50 | SBRDEC_SYNTAX_USAC))) {
      FDKreadBits(hBs, 3);
      FDKreadBits(hBs, 2);
    }
    headerExtra1 = FDKreadBit(hBs);
    headerExtra2 = FDKreadBit(hBs);
    FDKpushFor(hBs, 5 * headerExtra1 + 6 * headerExtra2);

    return HEADER_OK;
  }

  lastHeader = hHeaderData->bs_data;
  lastInfo = hHeaderData->bs_info;

  if ((flags & (SBRDEC_SYNTAX_RSVD50 | SBRDEC_SYNTAX_USAC)) && !fIsSbrData) {
    pBsData = &hHeaderData->bs_dflt;
  } else {
    pBsData = &hHeaderData->bs_data;
  }

  if (!(flags & (SBRDEC_SYNTAX_RSVD50 | SBRDEC_SYNTAX_USAC))) {
    hHeaderData->bs_info.ampResolution = FDKreadBits(hBs, 1);
  }

  pBsData->startFreq = FDKreadBits(hBs, 4);
  pBsData->stopFreq = FDKreadBits(hBs, 4);

  if (!(flags & (SBRDEC_SYNTAX_RSVD50 | SBRDEC_SYNTAX_USAC))) {
    hHeaderData->bs_info.xover_band = FDKreadBits(hBs, 3);
    FDKreadBits(hBs, 2);
  }

  headerExtra1 = FDKreadBits(hBs, 1);
  headerExtra2 = FDKreadBits(hBs, 1);

  if (headerExtra1) {
    pBsData->freqScale = FDKreadBits(hBs, 2);
    pBsData->alterScale = FDKreadBits(hBs, 1);
    pBsData->noise_bands = FDKreadBits(hBs, 2);
  } else {
    pBsData->freqScale = 2;
    pBsData->alterScale = 1;
    pBsData->noise_bands = 2;
  }

  if (headerExtra2) {
    pBsData->limiterBands = FDKreadBits(hBs, 2);
    pBsData->limiterGains = FDKreadBits(hBs, 2);
    pBsData->interpolFreq = FDKreadBits(hBs, 1);
    pBsData->smoothingLength = FDKreadBits(hBs, 1);
  } else {
    pBsData->limiterBands = 2;
    pBsData->limiterGains = 2;
    pBsData->interpolFreq = 1;
    pBsData->smoothingLength = 1;
  }

  if (hHeaderData->syncState < SBR_HEADER ||
      lastHeader.startFreq != pBsData->startFreq ||
      lastHeader.stopFreq != pBsData->stopFreq ||
      lastHeader.freqScale != pBsData->freqScale ||
      lastHeader.alterScale != pBsData->alterScale ||
      lastHeader.noise_bands != pBsData->noise_bands ||
      lastInfo.xover_band != hHeaderData->bs_info.xover_band) {
    return HEADER_RESET;
  }

  return HEADER_OK;
}

int sbrGetSyntheticCodedData(HANDLE_SBR_HEADER_DATA hHeaderData,
                             HANDLE_SBR_FRAME_DATA hFrameData,
                             HANDLE_FDK_BITSTREAM hBs, const UINT flags) {
  int i, bitsRead = 0;

  int add_harmonic_flag = FDKreadBits(hBs, 1);
  bitsRead++;

  if (add_harmonic_flag) {
    int nSfb = hHeaderData->freqBandData.nSfb[1];
    for (i = 0; i < ADD_HARMONICS_FLAGS_SIZE; i++) {
      int readBits = fMin(32, nSfb);
      nSfb -= readBits;
      if (readBits > 0) {
        hFrameData->addHarmonics[i] = FDKreadBits(hBs, readBits)
                                      << (32 - readBits);
      } else {
        hFrameData->addHarmonics[i] = 0;
      }

      bitsRead += readBits;
    }
    if (flags & SBRDEC_SYNTAX_USAC) {
      if (hHeaderData->bs_info.pvc_mode) {
        int bs_sinusoidal_position = 31;
        if (FDKreadBit(hBs)) {
          bs_sinusoidal_position = FDKreadBits(hBs, 5);
        }
        hFrameData->sinusoidal_position = bs_sinusoidal_position;
      }
    }
  } else {
    for (i = 0; i < ADD_HARMONICS_FLAGS_SIZE; i++)
      hFrameData->addHarmonics[i] = 0;
  }

  return (bitsRead);
}

static int extractExtendedData(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_FDK_BITSTREAM hBs
    ,
    HANDLE_PS_DEC hParametricStereoDec
) {
  INT nBitsLeft;
  int extended_data;
  int i, frameOk = 1;

  extended_data = FDKreadBits(hBs, 1);

  if (extended_data) {
    int cnt;
    int bPsRead = 0;

    cnt = FDKreadBits(hBs, 4);
    if (cnt == (1 << 4) - 1) cnt += FDKreadBits(hBs, 8);

    nBitsLeft = 8 * cnt;

    if (nBitsLeft > (INT)FDKgetValidBits(hBs)) {
      nBitsLeft = (INT)FDKgetValidBits(hBs);
      frameOk = 0;
    }

    while (nBitsLeft > 7) {
      int extension_id = FDKreadBits(hBs, 2);
      nBitsLeft -= 2;

      switch (extension_id) {
        case EXTENSION_ID_PS_CODING:

          if (hParametricStereoDec != NULL) {
            if (bPsRead &&
                !hParametricStereoDec->bsData[hParametricStereoDec->bsReadSlot]
                     .mpeg.bPsHeaderValid) {
              cnt = nBitsLeft >> 3;
              for (i = 0; i < cnt; i++) FDKreadBits(hBs, 8);
              nBitsLeft -= cnt * 8;
            } else {
              nBitsLeft -=
                  (INT)ReadPsData(hParametricStereoDec, hBs, nBitsLeft);
              bPsRead = 1;
            }
          }

          break;

        default:
          cnt = nBitsLeft >> 3;
          for (i = 0; i < cnt; i++) FDKreadBits(hBs, 8);
          nBitsLeft -= cnt * 8;
          break;
      }
    }

    if (nBitsLeft < 0) {
      frameOk = 0;
      goto bail;
    } else {
      FDKreadBits(hBs, nBitsLeft);
    }
  }

bail:
  return (frameOk);
}

int sbrGetChannelElement(HANDLE_SBR_HEADER_DATA hHeaderData,
                         HANDLE_SBR_FRAME_DATA hFrameDataLeft,
                         HANDLE_SBR_FRAME_DATA hFrameDataRight,
                         HANDLE_SBR_PREV_FRAME_DATA hFrameDataLeftPrev,
                         UCHAR pvc_mode_last, HANDLE_FDK_BITSTREAM hBs,
                         HANDLE_PS_DEC hParametricStereoDec, const UINT flags,
                         const int overlap) {
  int i, bs_coupling = COUPLING_OFF;
  const int nCh = (hFrameDataRight == NULL) ? 1 : 2;

  if (!(flags & (SBRDEC_SYNTAX_USAC | SBRDEC_SYNTAX_RSVD50))) {
    if (FDKreadBits(hBs, 1)) {
      FDKreadBits(hBs, 4);
      if ((flags & SBRDEC_SYNTAX_SCAL) || (nCh == 2)) {
        FDKreadBits(hBs, 4);
      }
    }
  }

  if (nCh == 2) {
    bs_coupling = FDKreadBits(hBs, 1);
    if (bs_coupling) {
      hFrameDataLeft->coupling = COUPLING_LEVEL;
      hFrameDataRight->coupling = COUPLING_BAL;
    } else {
      hFrameDataLeft->coupling = COUPLING_OFF;
      hFrameDataRight->coupling = COUPLING_OFF;
    }
  } else {
    if (flags & SBRDEC_SYNTAX_SCAL) {
      FDKreadBits(hBs, 1);
    }
    hFrameDataLeft->coupling = COUPLING_OFF;
  }

  if (flags & (SBRDEC_SYNTAX_USAC | SBRDEC_SYNTAX_RSVD50)) {
    if (flags & SBRDEC_USAC_HARMONICSBR) {
      hFrameDataLeft->sbrPatchingMode = FDKreadBit(hBs);
      if (hFrameDataLeft->sbrPatchingMode == 0) {
        hFrameDataLeft->sbrOversamplingFlag = FDKreadBit(hBs);
        if (FDKreadBit(hBs)) { /* sbrPitchInBinsFlag */
          hFrameDataLeft->sbrPitchInBins = FDKreadBits(hBs, 7);
        } else {
          hFrameDataLeft->sbrPitchInBins = 0;
        }
      } else {
        hFrameDataLeft->sbrOversamplingFlag = 0;
        hFrameDataLeft->sbrPitchInBins = 0;
      }

      if (nCh == 2) {
        if (bs_coupling) {
          hFrameDataRight->sbrPatchingMode = hFrameDataLeft->sbrPatchingMode;
          hFrameDataRight->sbrOversamplingFlag =
              hFrameDataLeft->sbrOversamplingFlag;
          hFrameDataRight->sbrPitchInBins = hFrameDataLeft->sbrPitchInBins;
        } else {
          hFrameDataRight->sbrPatchingMode = FDKreadBit(hBs);
          if (hFrameDataRight->sbrPatchingMode == 0) {
            hFrameDataRight->sbrOversamplingFlag = FDKreadBit(hBs);
            if (FDKreadBit(hBs)) { /* sbrPitchInBinsFlag */
              hFrameDataRight->sbrPitchInBins = FDKreadBits(hBs, 7);
            } else {
              hFrameDataRight->sbrPitchInBins = 0;
            }
          } else {
            hFrameDataRight->sbrOversamplingFlag = 0;
            hFrameDataRight->sbrPitchInBins = 0;
          }
        }
      }
    } else {
      if (nCh == 2) {
        hFrameDataRight->sbrPatchingMode = 1;
        hFrameDataRight->sbrOversamplingFlag = 0;
        hFrameDataRight->sbrPitchInBins = 0;
      }

      hFrameDataLeft->sbrPatchingMode = 1;
      hFrameDataLeft->sbrOversamplingFlag = 0;
      hFrameDataLeft->sbrPitchInBins = 0;
    }
  } else {
    if (nCh == 2) {
      hFrameDataRight->sbrPatchingMode = 1;
      hFrameDataRight->sbrOversamplingFlag = 0;
      hFrameDataRight->sbrPitchInBins = 0;
    }

    hFrameDataLeft->sbrPatchingMode = 1;
    hFrameDataLeft->sbrOversamplingFlag = 0;
    hFrameDataLeft->sbrPitchInBins = 0;
  }

  if (hHeaderData->bs_info.pvc_mode) {
    FDK_ASSERT(nCh == 1);
    if (!extractPvcFrameInfo(hBs, hHeaderData, hFrameDataLeft,
                             hFrameDataLeftPrev, pvc_mode_last, flags))
      return 0;

    if (!checkFrameInfo(&hFrameDataLeft->frameInfo,
                        hHeaderData->numberTimeSlots, overlap,
                        hHeaderData->timeStep))
      return 0;
  } else {
    if (!extractFrameInfo(hBs, hHeaderData, hFrameDataLeft, 1, flags)) return 0;

    if (!checkFrameInfo(&hFrameDataLeft->frameInfo,
                        hHeaderData->numberTimeSlots, overlap,
                        hHeaderData->timeStep))
      return 0;
  }
  if (nCh == 2) {
    if (hFrameDataLeft->coupling) {
      FDKmemcpy(&hFrameDataRight->frameInfo, &hFrameDataLeft->frameInfo,
                sizeof(FRAME_INFO));
      hFrameDataRight->ampResolutionCurrentFrame =
          hFrameDataLeft->ampResolutionCurrentFrame;
    } else {
      if (!extractFrameInfo(hBs, hHeaderData, hFrameDataRight, 2, flags))
        return 0;

      if (!checkFrameInfo(&hFrameDataRight->frameInfo,
                          hHeaderData->numberTimeSlots, overlap,
                          hHeaderData->timeStep))
        return 0;
    }
  }

  sbrGetDirectionControlData(hFrameDataLeft, hBs, flags,
                             hHeaderData->bs_info.pvc_mode);
  if (nCh == 2) {
    sbrGetDirectionControlData(hFrameDataRight, hBs, flags, 0);
  }

  for (i = 0; i < hHeaderData->freqBandData.nInvfBands; i++) {
    hFrameDataLeft->sbr_invf_mode[i] = (INVF_MODE)FDKreadBits(hBs, 2);
  }
  if (nCh == 2) {
    if (hFrameDataLeft->coupling) {
      for (i = 0; i < hHeaderData->freqBandData.nInvfBands; i++) {
        hFrameDataRight->sbr_invf_mode[i] = hFrameDataLeft->sbr_invf_mode[i];
      }
    } else {
      for (i = 0; i < hHeaderData->freqBandData.nInvfBands; i++) {
        hFrameDataRight->sbr_invf_mode[i] = (INVF_MODE)FDKreadBits(hBs, 2);
      }
    }
  }

  if (nCh == 1) {
    if (hHeaderData->bs_info.pvc_mode) {
      if (!sbrGetPvcEnvelope(hHeaderData, hFrameDataLeft, hBs, flags,
                             hHeaderData->bs_info.pvc_mode))
        return 0;
    } else if (!sbrGetEnvelope(hHeaderData, hFrameDataLeft, hBs, flags))
      return 0;

    sbrGetNoiseFloorData(hHeaderData, hFrameDataLeft, hBs);
  } else if (hFrameDataLeft->coupling) {
    if (!sbrGetEnvelope(hHeaderData, hFrameDataLeft, hBs, flags)) {
      return 0;
    }

    sbrGetNoiseFloorData(hHeaderData, hFrameDataLeft, hBs);

    if (!sbrGetEnvelope(hHeaderData, hFrameDataRight, hBs, flags)) {
      return 0;
    }
    sbrGetNoiseFloorData(hHeaderData, hFrameDataRight, hBs);
  } else {

    if (!sbrGetEnvelope(hHeaderData, hFrameDataLeft, hBs, flags)) return 0;

    if (!sbrGetEnvelope(hHeaderData, hFrameDataRight, hBs, flags)) return 0;

    sbrGetNoiseFloorData(hHeaderData, hFrameDataLeft, hBs);

    sbrGetNoiseFloorData(hHeaderData, hFrameDataRight, hBs);
  }

  sbrGetSyntheticCodedData(hHeaderData, hFrameDataLeft, hBs, flags);
  if (nCh == 2) {
    sbrGetSyntheticCodedData(hHeaderData, hFrameDataRight, hBs, flags);
  }

  if (!(flags & (SBRDEC_SYNTAX_USAC | SBRDEC_SYNTAX_RSVD50))) {
    if (!extractExtendedData(hHeaderData, hBs, hParametricStereoDec)) {
      return 0;
    }
  }

  return 1;
}

void sbrGetDirectionControlData(
    HANDLE_SBR_FRAME_DATA h_frame_data,
    HANDLE_FDK_BITSTREAM hBs,
    const UINT flags, const int bs_pvc_mode)

{
  int i;
  int indepFlag = 0;

  if (flags & (SBRDEC_SYNTAX_USAC | SBRDEC_SYNTAX_RSVD50)) {
    indepFlag = flags & SBRDEC_USAC_INDEP;
  }

  if (bs_pvc_mode == 0) {
    i = 0;
    if (indepFlag) {
      h_frame_data->domain_vec[i++] = 0;
    }
    for (; i < h_frame_data->frameInfo.nEnvelopes; i++) {
      h_frame_data->domain_vec[i] = FDKreadBits(hBs, 1);
    }
  }

  i = 0;
  if (indepFlag) {
    h_frame_data->domain_vec_noise[i++] = 0;
  }
  for (; i < h_frame_data->frameInfo.nNoiseEnvelopes; i++) {
    h_frame_data->domain_vec_noise[i] = FDKreadBits(hBs, 1);
  }
}

void sbrGetNoiseFloorData(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_frame_data,
    HANDLE_FDK_BITSTREAM hBs)
{
  int i, j;
  int delta;
  COUPLING_MODE coupling;
  int noNoiseBands = hHeaderData->freqBandData.nNfb;

  Huffman hcb_noiseF;
  Huffman hcb_noise;
  int envDataTableCompFactor;

  coupling = h_frame_data->coupling;

  if (coupling == COUPLING_BAL) {
    hcb_noise = (Huffman)&FDK_sbrDecoder_sbr_huffBook_NoiseBalance11T;
    hcb_noiseF =
        (Huffman)&FDK_sbrDecoder_sbr_huffBook_EnvBalance11F;
    envDataTableCompFactor = 1;
  } else {
    hcb_noise = (Huffman)&FDK_sbrDecoder_sbr_huffBook_NoiseLevel11T;
    hcb_noiseF =
        (Huffman)&FDK_sbrDecoder_sbr_huffBook_EnvLevel11F;
    envDataTableCompFactor = 0;
  }

  for (i = 0; i < h_frame_data->frameInfo.nNoiseEnvelopes; i++) {
    if (h_frame_data->domain_vec_noise[i] == 0) {
      if (coupling == COUPLING_BAL) {
        h_frame_data->sbrNoiseFloorLevel[i * noNoiseBands] =
            (FIXP_SGL)(((int)FDKreadBits(hBs, 5)) << envDataTableCompFactor);
      } else {
        h_frame_data->sbrNoiseFloorLevel[i * noNoiseBands] =
            (FIXP_SGL)(int)FDKreadBits(hBs, 5);
      }

      for (j = 1; j < noNoiseBands; j++) {
        delta = DecodeHuffmanCW(hcb_noiseF, hBs);
        h_frame_data->sbrNoiseFloorLevel[i * noNoiseBands + j] =
            (FIXP_SGL)(delta << envDataTableCompFactor);
      }
    } else {
      for (j = 0; j < noNoiseBands; j++) {
        delta = DecodeHuffmanCW(hcb_noise, hBs);
        h_frame_data->sbrNoiseFloorLevel[i * noNoiseBands + j] =
            (FIXP_SGL)(delta << envDataTableCompFactor);
      }
    }
  }
}

static const UCHAR mapNsMode2ns[2][2] = {
    {16, 4},
    {12, 3}
};

static int sbrGetPvcEnvelope(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_frame_data,
    HANDLE_FDK_BITSTREAM hBs,
    const UINT flags, const UINT pvcMode) {
  int divMode, nsMode;
  int indepFlag = flags & SBRDEC_USAC_INDEP;
  UCHAR *pvcID = h_frame_data->pvcID;

  divMode = FDKreadBits(hBs, PVC_DIVMODE_BITS);
  nsMode = FDKreadBit(hBs);
  FDK_ASSERT((pvcMode == 1) || (pvcMode == 2));
  h_frame_data->ns = mapNsMode2ns[pvcMode - 1][nsMode];

  if (divMode <= 3) {
    int i, k = 1, sum_length = 0, reuse_pcvID;

    indepFlag ? (reuse_pcvID = 0) : (reuse_pcvID = FDKreadBit(hBs));
    if (reuse_pcvID) {
      pvcID[0] = hHeaderData->pvcIDprev;
    } else {
      pvcID[0] = FDKreadBits(hBs, PVC_PVCID_BITS);
    }

    for (i = 0; i < divMode; i++) {
      int length, numBits = 4;

      if (sum_length >= 13) {
        numBits = 1;
      } else if (sum_length >= 11) {
        numBits = 2;
      } else if (sum_length >= 7) {
        numBits = 3;
      }

      length = FDKreadBits(hBs, numBits);
      sum_length += length + 1;
      if (sum_length >= PVC_NTIMESLOT) {
        return 0;
      }
      for (; length--; k++) {
        pvcID[k] = pvcID[k - 1];
      }
      pvcID[k++] = FDKreadBits(hBs, PVC_PVCID_BITS);
    }
    for (; k < 16; k++) {
      pvcID[k] = pvcID[k - 1];
    }
  } else {
    int num_grid_info, fixed_length, grid_info, j, k = 0;

    divMode -= 4;
    num_grid_info = 2 << divMode;
    fixed_length = 8 >> divMode;
    FDK_ASSERT(num_grid_info * fixed_length == PVC_NTIMESLOT);

    indepFlag ? (grid_info = 1) : (grid_info = FDKreadBit(hBs));
    if (grid_info) {
      pvcID[k++] = FDKreadBits(hBs, PVC_PVCID_BITS);
    } else {
      pvcID[k++] = hHeaderData->pvcIDprev;
    }
    j = fixed_length - 1;
    for (; j--; k++) {
      pvcID[k] = pvcID[k - 1];
    }
    num_grid_info--;

    for (; num_grid_info--;) {
      j = fixed_length;
      grid_info = FDKreadBit(hBs);
      if (grid_info) {
        pvcID[k++] = FDKreadBits(hBs, PVC_PVCID_BITS);
        j--;
      }
      for (; j--; k++) {
        pvcID[k] = pvcID[k - 1];
      }
    }
  }

  hHeaderData->pvcIDprev = pvcID[PVC_NTIMESLOT - 1];

  h_frame_data->iTESactive = (UCHAR)0;

  return 1;
}
static int sbrGetEnvelope(
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_frame_data,
    HANDLE_FDK_BITSTREAM hBs,
    const UINT flags) {
  int i, j;
  UCHAR no_band[MAX_ENVELOPES];
  int delta = 0;
  int offset = 0;
  COUPLING_MODE coupling = h_frame_data->coupling;
  int ampRes = hHeaderData->bs_info.ampResolution;
  int nEnvelopes = h_frame_data->frameInfo.nEnvelopes;
  int envDataTableCompFactor;
  int start_bits, start_bits_balance;
  Huffman hcb_t, hcb_f;

  h_frame_data->nScaleFactors = 0;

  if ((h_frame_data->frameInfo.frameClass == 0) && (nEnvelopes == 1)) {
    if (flags & SBRDEC_ELD_GRID)
      ampRes = h_frame_data->ampResolutionCurrentFrame;
    else
      ampRes = 0;
  }
  h_frame_data->ampResolutionCurrentFrame = ampRes;

  if (ampRes == 1) {
    start_bits = 6;
    start_bits_balance = 5;
  } else {
    start_bits = 7;
    start_bits_balance = 6;
  }

  for (i = 0; i < nEnvelopes; i++) {
    no_band[i] =
        hHeaderData->freqBandData.nSfb[h_frame_data->frameInfo.freqRes[i]];
    h_frame_data->nScaleFactors += no_band[i];
  }
  if (h_frame_data->nScaleFactors > MAX_NUM_ENVELOPE_VALUES) return 0;

  if (coupling == COUPLING_BAL) {
    envDataTableCompFactor = 1;
    if (ampRes == 0) {
      hcb_t = (Huffman)&FDK_sbrDecoder_sbr_huffBook_EnvBalance10T;
      hcb_f = (Huffman)&FDK_sbrDecoder_sbr_huffBook_EnvBalance10F;
    } else {
      hcb_t = (Huffman)&FDK_sbrDecoder_sbr_huffBook_EnvBalance11T;
      hcb_f = (Huffman)&FDK_sbrDecoder_sbr_huffBook_EnvBalance11F;
    }
  } else {
    envDataTableCompFactor = 0;
    if (ampRes == 0) {
      hcb_t = (Huffman)&FDK_sbrDecoder_sbr_huffBook_EnvLevel10T;
      hcb_f = (Huffman)&FDK_sbrDecoder_sbr_huffBook_EnvLevel10F;
    } else {
      hcb_t = (Huffman)&FDK_sbrDecoder_sbr_huffBook_EnvLevel11T;
      hcb_f = (Huffman)&FDK_sbrDecoder_sbr_huffBook_EnvLevel11F;
    }
  }

  h_frame_data->iTESactive = (UCHAR)0;
  for (j = 0, offset = 0; j < nEnvelopes; j++) {
    if (h_frame_data->domain_vec[j] == 0) {
      if (coupling == COUPLING_BAL) {
        h_frame_data->iEnvelope[offset] =
            (FIXP_SGL)(((int)FDKreadBits(hBs, start_bits_balance))
                       << envDataTableCompFactor);
      } else {
        h_frame_data->iEnvelope[offset] =
            (FIXP_SGL)(int)FDKreadBits(hBs, start_bits);
      }
    }

    for (i = (1 - h_frame_data->domain_vec[j]); i < no_band[j]; i++) {
      if (h_frame_data->domain_vec[j] == 0) {
        delta = DecodeHuffmanCW(hcb_f, hBs);
      } else {
        delta = DecodeHuffmanCW(hcb_t, hBs);
      }

      h_frame_data->iEnvelope[offset + i] =
          (FIXP_SGL)(delta << envDataTableCompFactor);
    }
    if ((flags & SBRDEC_SYNTAX_USAC) && (flags & SBRDEC_USAC_ITES)) {
      int bs_temp_shape = FDKreadBit(hBs);
      FDK_ASSERT(j < 8);
      h_frame_data->iTESactive |= (UCHAR)(bs_temp_shape << j);
      if (bs_temp_shape) {
        h_frame_data->interTempShapeMode[j] =
            FDKread2Bits(hBs);
      } else {
        h_frame_data->interTempShapeMode[j] = 0;
      }
    }
    offset += no_band[j];
  }

#if ENV_EXP_FRACT
  for (i = 0; i < h_frame_data->nScaleFactors; i++) {
    h_frame_data->iEnvelope[i] <<= ENV_EXP_FRACT;
  }
#endif

  return 1;
}

static int generateFixFixOnly(FRAME_INFO *hSbrFrameInfo, int tranPosInternal,
                              int numberTimeSlots, const UINT flags) {
  int nEnv, i, tranIdx;
  const int *pTable;

  if (tranPosInternal >= numberTimeSlots) {
    return 0;
  }

  switch (numberTimeSlots) {
    case 8:
      pTable = FDK_sbrDecoder_envelopeTable_8[tranPosInternal];
      break;
    case 15:
      pTable = FDK_sbrDecoder_envelopeTable_15[tranPosInternal];
      break;
    case 16:
      pTable = FDK_sbrDecoder_envelopeTable_16[tranPosInternal];
      break;
    default:
      return 0;
  }

  nEnv = pTable[0];
  for (i = 1; i < nEnv; i++) hSbrFrameInfo->borders[i] = pTable[i + 2];
  hSbrFrameInfo->borders[0] = 0;
  hSbrFrameInfo->borders[nEnv] = numberTimeSlots;
  hSbrFrameInfo->nEnvelopes = nEnv;

  tranIdx = hSbrFrameInfo->tranEnv = pTable[1];

  hSbrFrameInfo->bordersNoise[0] = 0;
  hSbrFrameInfo->bordersNoise[1] =
      hSbrFrameInfo->borders[tranIdx ? tranIdx : 1];
  hSbrFrameInfo->bordersNoise[2] = numberTimeSlots;
  hSbrFrameInfo->nNoiseEnvelopes = 2;

  return 1;
}

static int extractLowDelayGrid(
    HANDLE_FDK_BITSTREAM hBitBuf,
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA
        h_frame_data,
    int timeSlots, const UINT flags) {
  FRAME_INFO *pFrameInfo = &h_frame_data->frameInfo;
  INT numberTimeSlots = hHeaderData->numberTimeSlots;
  INT temp = 0, k;

  h_frame_data->frameInfo.frameClass = 0;

  switch (timeSlots) {
    case 8:
      temp = FDKreadBits(hBitBuf, 3);
      break;

    case 16:
    case 15:
      temp = FDKreadBits(hBitBuf, 4);
      break;

    default:
      return 0;
  }

  if (temp >= timeSlots) {
    return 0;
  }

  if (!generateFixFixOnly(pFrameInfo, temp, numberTimeSlots, flags)) {
    return 0;
  }

  for (k = 0; k < pFrameInfo->nEnvelopes; k++) {
    pFrameInfo->freqRes[k] =
        (UCHAR)FDKreadBits(hBitBuf, 1);
  }

  return 1;
}

int extractPvcFrameInfo(
    HANDLE_FDK_BITSTREAM hBs,
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_frame_data,
    HANDLE_SBR_PREV_FRAME_DATA h_prev_frame_data,
    UCHAR pvc_mode_last,
    const UINT flags) {
  FRAME_INFO *pFrameInfo = &h_frame_data->frameInfo;
  FRAME_INFO *pPrevFrameInfo = &h_prev_frame_data->prevFrameInfo;
  int bs_var_len_hf, bs_noise_position;
  bs_noise_position = FDKreadBits(hBs, 4);
  bs_var_len_hf = FDKreadBit(hBs);
  pFrameInfo->noisePosition = bs_noise_position;
  pFrameInfo->tranEnv = -1;

  pFrameInfo->nEnvelopes = 1;
  pFrameInfo->nNoiseEnvelopes = 1;
  pFrameInfo->freqRes[0] = 0;

  if (bs_var_len_hf) {
    pFrameInfo->varLength = FDKreadBits(hBs, 2) + 1;
    if (pFrameInfo->varLength > 3) {
      pFrameInfo->varLength =
          0;
      return 0;
    }
  } else {
    pFrameInfo->varLength = 0;
  }

  if (bs_noise_position) {
    pFrameInfo->nEnvelopes = 2;
    pFrameInfo->nNoiseEnvelopes = 2;
    FDKmemclear(pFrameInfo->freqRes, sizeof(pFrameInfo->freqRes));
  }

  if (hHeaderData->bs_info.pvc_mode > 0) {

    FDK_ASSERT((pFrameInfo->nEnvelopes == 1) || (pFrameInfo->nEnvelopes == 2));

    if (pPrevFrameInfo->nEnvelopes > 0) {
      pFrameInfo->borders[0] =
          pPrevFrameInfo->borders[pPrevFrameInfo->nEnvelopes] - PVC_NTIMESLOT;
      FDK_ASSERT(pFrameInfo->borders[0] <= 3);
    } else {
      pFrameInfo->borders[0] = 0;
    }

    pFrameInfo->borders[pFrameInfo->nEnvelopes] = 16 + pFrameInfo->varLength;

    if (pFrameInfo->nEnvelopes == 2) {
      pFrameInfo->borders[1] = pFrameInfo->noisePosition;
    }

    if (pvc_mode_last == 0) {
      pFrameInfo->pvcBorders[0] = pFrameInfo->borders[0];
    } else {
      pFrameInfo->pvcBorders[0] = 0;
    }
    if (pFrameInfo->nEnvelopes == 2) {
      pFrameInfo->pvcBorders[1] = pFrameInfo->borders[1];
    }
    pFrameInfo->pvcBorders[pFrameInfo->nEnvelopes] = 16;

    for (INT i = 0; i <= pFrameInfo->nNoiseEnvelopes; i++) {
      pFrameInfo->bordersNoise[i] = pFrameInfo->borders[i];
    }

    pFrameInfo->tranEnv = -1;
  }
  return 1;
}

int extractFrameInfo(
    HANDLE_FDK_BITSTREAM hBs,
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA h_frame_data,
    const UINT nrOfChannels, const UINT flags) {
  FRAME_INFO *pFrameInfo = &h_frame_data->frameInfo;
  int numberTimeSlots = hHeaderData->numberTimeSlots;
  int pointer_bits = 0, nEnv = 0, b = 0, border, i, n = 0, k, p, aL, aR, nL, nR,
      temp = 0, staticFreqRes;
  UCHAR frameClass;

  if (flags & SBRDEC_ELD_GRID) {
    frameClass = FDKreadBits(hBs, 1);
    if (frameClass == 1) {
      pFrameInfo->frameClass = frameClass;
      int err = 1;
      err = extractLowDelayGrid(hBs, hHeaderData, h_frame_data, numberTimeSlots,
                                flags);
      return err;
    }
  } else {
    frameClass = FDKreadBits(hBs, 2);
  }

  switch (frameClass) {
    case 0:
      temp = FDKreadBits(hBs, 2);
      nEnv = (int)(1 << temp);

      if ((flags & SBRDEC_ELD_GRID) && (nEnv == 1))
        h_frame_data->ampResolutionCurrentFrame =
            FDKreadBits(hBs, 1);

      staticFreqRes = FDKreadBits(hBs, 1);

      if (flags & (SBRDEC_SYNTAX_USAC | SBRDEC_SYNTAX_RSVD50)) {
        if (nEnv > MAX_ENVELOPES_USAC) return 0;
      } else

        b = nEnv + 1;
      switch (nEnv) {
        case 1:
          switch (numberTimeSlots) {
            case 15:
              FDKmemcpy(pFrameInfo, &FDK_sbrDecoder_sbr_frame_info1_15,
                        sizeof(FRAME_INFO));
              break;
            case 16:
              FDKmemcpy(pFrameInfo, &FDK_sbrDecoder_sbr_frame_info1_16,
                        sizeof(FRAME_INFO));
              break;
            default:
              FDK_ASSERT(0);
          }
          break;
        case 2:
          switch (numberTimeSlots) {
            case 15:
              FDKmemcpy(pFrameInfo, &FDK_sbrDecoder_sbr_frame_info2_15,
                        sizeof(FRAME_INFO));
              break;
            case 16:
              FDKmemcpy(pFrameInfo, &FDK_sbrDecoder_sbr_frame_info2_16,
                        sizeof(FRAME_INFO));
              break;
            default:
              FDK_ASSERT(0);
          }
          break;
        case 4:
          switch (numberTimeSlots) {
            case 15:
              FDKmemcpy(pFrameInfo, &FDK_sbrDecoder_sbr_frame_info4_15,
                        sizeof(FRAME_INFO));
              break;
            case 16:
              FDKmemcpy(pFrameInfo, &FDK_sbrDecoder_sbr_frame_info4_16,
                        sizeof(FRAME_INFO));
              break;
            default:
              FDK_ASSERT(0);
          }
          break;
        case 8:
#if (MAX_ENVELOPES >= 8)
          switch (numberTimeSlots) {
            case 15:
              FDKmemcpy(pFrameInfo, &FDK_sbrDecoder_sbr_frame_info8_15,
                        sizeof(FRAME_INFO));
              break;
            case 16:
              FDKmemcpy(pFrameInfo, &FDK_sbrDecoder_sbr_frame_info8_16,
                        sizeof(FRAME_INFO));
              break;
            default:
              FDK_ASSERT(0);
          }
          break;
#else
          return 0;
#endif
      }
      if (!staticFreqRes) {
        for (i = 0; i < nEnv; i++) pFrameInfo->freqRes[i] = 0;
      }

      break;
    case 1:
    case 2:
      temp = FDKreadBits(hBs, 2);

      n = FDKreadBits(hBs, 2);

      nEnv = n + 1;
      b = nEnv + 1;

      break;
  }

  switch (frameClass) {
    case 1:
      pFrameInfo->borders[0] = 0;
      border = temp + numberTimeSlots;
      i = b - 1;
      pFrameInfo->borders[i] = border;

      for (k = 0; k < n; k++) {
        temp = FDKreadBits(hBs, 2);
        border -= (2 * temp + 2);
        pFrameInfo->borders[--i] = border;
      }

      pointer_bits = DFRACT_BITS - 1 - CountLeadingBits((FIXP_DBL)(n + 1));
      p = FDKreadBits(hBs, pointer_bits);

      if (p > n + 1) return 0;

      pFrameInfo->tranEnv = p ? n + 2 - p : -1;

      for (k = n; k >= 0; k--) {
        pFrameInfo->freqRes[k] = FDKreadBits(hBs, 1);
      }

      if (p == 0 || p == 1)
        pFrameInfo->bordersNoise[1] = pFrameInfo->borders[n];
      else
        pFrameInfo->bordersNoise[1] = pFrameInfo->borders[pFrameInfo->tranEnv];

      break;

    case 2:
      border = temp;
      pFrameInfo->borders[0] = border;

      for (k = 1; k <= n; k++) {
        temp = FDKreadBits(hBs, 2);
        border += (2 * temp + 2);
        pFrameInfo->borders[k] = border;
      }
      pFrameInfo->borders[k] = numberTimeSlots;

      pointer_bits = DFRACT_BITS - 1 - CountLeadingBits((FIXP_DBL)(n + 1));
      p = FDKreadBits(hBs, pointer_bits);
      if (p > n + 1) return 0;

      if (p == 0 || p == 1)
        pFrameInfo->tranEnv = -1;
      else
        pFrameInfo->tranEnv = p - 1;

      for (k = 0; k <= n; k++) {
        pFrameInfo->freqRes[k] = FDKreadBits(hBs, 1);
      }

      switch (p) {
        case 0:
          pFrameInfo->bordersNoise[1] = pFrameInfo->borders[1];
          break;
        case 1:
          pFrameInfo->bordersNoise[1] = pFrameInfo->borders[n];
          break;
        default:
          pFrameInfo->bordersNoise[1] =
              pFrameInfo->borders[pFrameInfo->tranEnv];
          break;
      }

      break;

    case 3:

      aL = FDKreadBits(hBs, 2);

      aR = FDKreadBits(hBs, 2) + numberTimeSlots;

      nL = FDKreadBits(hBs, 2);

      nR = FDKreadBits(hBs, 2);

      nEnv = nL + nR + 1;
      if (nEnv > MAX_ENVELOPES) return 0;
      b = nEnv + 1;

      border = aL;
      pFrameInfo->borders[0] = border;

      for (k = 1; k <= nL; k++) {
        temp = FDKreadBits(hBs, 2);
        border += (2 * temp + 2);
        pFrameInfo->borders[k] = border;
      }

      border = aR;
      i = nEnv;

      pFrameInfo->borders[i] = border;

      for (k = 0; k < nR; k++) {
        temp = FDKreadBits(hBs, 2);
        border -= (2 * temp + 2);
        pFrameInfo->borders[--i] = border;
      }

      pointer_bits =
          DFRACT_BITS - 1 - CountLeadingBits((FIXP_DBL)(nL + nR + 1));
      p = FDKreadBits(hBs, pointer_bits);

      if (p > nL + nR + 1) return 0;

      pFrameInfo->tranEnv = p ? b - p : -1;

      for (k = 0; k < nEnv; k++) {
        pFrameInfo->freqRes[k] = FDKreadBits(hBs, 1);
      }

      pFrameInfo->bordersNoise[0] = aL;

      if (nEnv == 1) {
        pFrameInfo->bordersNoise[1] = aR;
      } else {
        if (p == 0 || p == 1)
          pFrameInfo->bordersNoise[1] = pFrameInfo->borders[nEnv - 1];
        else
          pFrameInfo->bordersNoise[1] =
              pFrameInfo->borders[pFrameInfo->tranEnv];
        pFrameInfo->bordersNoise[2] = aR;
      }
      break;
  }

  pFrameInfo->nEnvelopes = nEnv;

  if (nEnv == 1)
    pFrameInfo->nNoiseEnvelopes = 1;
  else
    pFrameInfo->nNoiseEnvelopes = 2;

  pFrameInfo->frameClass = frameClass;

  if (pFrameInfo->frameClass == 2 || pFrameInfo->frameClass == 1) {
    pFrameInfo->bordersNoise[0] = pFrameInfo->borders[0];
    pFrameInfo->bordersNoise[pFrameInfo->nNoiseEnvelopes] =
        pFrameInfo->borders[nEnv];
  }

  return 1;
}

static int checkFrameInfo(
    FRAME_INFO *pFrameInfo,
    int numberOfTimeSlots,
    int overlap,
    int timeStep)
{
  int maxPos, i, j;
  int startPos;
  int stopPos;
  int tranEnv;
  int startPosNoise;
  int stopPosNoise;
  int nEnvelopes = pFrameInfo->nEnvelopes;
  int nNoiseEnvelopes = pFrameInfo->nNoiseEnvelopes;

  if (nEnvelopes < 1 || nEnvelopes > MAX_ENVELOPES) return 0;

  if (nNoiseEnvelopes > MAX_NOISE_ENVELOPES) return 0;

  startPos = pFrameInfo->borders[0];
  stopPos = pFrameInfo->borders[nEnvelopes];
  tranEnv = pFrameInfo->tranEnv;
  startPosNoise = pFrameInfo->bordersNoise[0];
  stopPosNoise = pFrameInfo->bordersNoise[nNoiseEnvelopes];

  if (overlap < 0 || overlap > (3 * (4))) {
    return 0;
  }
  if (timeStep < 1 || timeStep > (4)) {
    return 0;
  }
  maxPos = numberOfTimeSlots + (overlap / timeStep);

  if ((startPos < 0) || (startPos >= stopPos)) return 0;
  if (startPos > maxPos - numberOfTimeSlots)
    return 0;
  if (stopPos < numberOfTimeSlots)
    return 0;
  if (stopPos > maxPos) return 0;

  for (i = 0; i < nEnvelopes; i++) {
    if (pFrameInfo->borders[i] >= pFrameInfo->borders[i + 1]) return 0;
  }

  if (tranEnv > nEnvelopes) return 0;

  if (nEnvelopes == 1 && nNoiseEnvelopes > 1) return 0;

  if (startPos != startPosNoise || stopPos != stopPosNoise) return 0;

  for (i = 0; i < nNoiseEnvelopes; i++) {
    if (pFrameInfo->bordersNoise[i] >= pFrameInfo->bordersNoise[i + 1])
      return 0;
  }

  for (i = 0; i < nNoiseEnvelopes; i++) {
    startPosNoise = pFrameInfo->bordersNoise[i];

    for (j = 0; j < nEnvelopes; j++) {
      if (pFrameInfo->borders[j] == startPosNoise) break;
    }
    if (j == nEnvelopes) return 0;
  }

  return 1;
}
