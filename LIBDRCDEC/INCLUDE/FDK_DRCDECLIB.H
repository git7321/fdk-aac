/************************* MPEG-D DRC decoder library **************************

   Author(s):   Bernhard Neugebauer

   Description: MPEG-D DRC Decoder

*******************************************************************************/

#ifndef FDK_DRCDECLIB_H
#define FDK_DRCDECLIB_H

#include "FDK_bitstream.h"
#include "FDK_audio.h"
#include "common_fix.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DRC_DEC_LOUDNESS_NOT_PRESENT (LONG)0x7FFFFFFE

typedef struct s_drc_decoder* HANDLE_DRC_DECODER;
typedef struct s_uni_drc_interface* HANDLE_UNI_DRC_INTERFACE;
typedef struct s_selection_process_output* HANDLE_SEL_PROC_OUTPUT;

typedef enum {
  DRC_DEC_SELECTION = 0x1,
  DRC_DEC_GAIN = 0x2,
  DRC_DEC_ALL = 0x3
} DRC_DEC_FUNCTIONAL_RANGE;

typedef enum {
  DRC_DEC_BOOST,
  DRC_DEC_COMPRESS,
  DRC_DEC_LOUDNESS_NORMALIZATION_ON,
  DRC_DEC_TARGET_LOUDNESS,
  DRC_DEC_EFFECT_TYPE,
  DRC_DEC_EFFECT_TYPE_FALLBACK_CODE,
  DRC_DEC_LOUDNESS_MEASUREMENT_METHOD,
  DRC_DEC_ALBUM_MODE,
  DRC_DEC_DOWNMIX_ID,
  DRC_DEC_TARGET_CHANNEL_COUNT_REQUESTED,
  DRC_DEC_BASE_CHANNEL_COUNT,
  DRC_DEC_FRAME_SIZE,
  DRC_DEC_SAMPLE_RATE,
  DRC_DEC_IS_MULTIBAND_DRC_1,
  DRC_DEC_IS_MULTIBAND_DRC_2,
  DRC_DEC_IS_ACTIVE,
  DRC_DEC_TARGET_CHANNEL_COUNT_SELECTED,
  DRC_DEC_OUTPUT_LOUDNESS
} DRC_DEC_USERPARAM;

typedef enum {
  DRC_DEC_OK = 0,

  DRC_DEC_NOT_OK = -10000,
  DRC_DEC_OUT_OF_MEMORY,
  DRC_DEC_NOT_OPENED,
  DRC_DEC_NOT_READY,
  DRC_DEC_PARAM_OUT_OF_RANGE,
  DRC_DEC_INVALID_PARAM,
  DRC_DEC_UNSUPPORTED_FUNCTION
} DRC_DEC_ERROR;

typedef enum {
  DRC_DEC_TEST_TIME_DOMAIN = -100,
  DRC_DEC_TEST_QMF_DOMAIN,
  DRC_DEC_TEST_STFT_DOMAIN,
  DRC_DEC_CODEC_MODE_UNDEFINED = -1,
  DRC_DEC_MPEG_4_AAC,
  DRC_DEC_MPEG_D_USAC,
  DRC_DEC_MPEG_H_3DA
} DRC_DEC_CODEC_MODE;

typedef enum {
  DRC_DEC_DRC1,
  DRC_DEC_DRC1_DRC2,
  DRC_DEC_DRC2,
  DRC_DEC_DRC3,
  DRC_DEC_DRC2_DRC3
} DRC_DEC_LOCATION;

DRC_DEC_ERROR
FDK_drcDec_Open(HANDLE_DRC_DECODER* phDrcDec,
                const DRC_DEC_FUNCTIONAL_RANGE functionalRange);

DRC_DEC_ERROR
FDK_drcDec_SetCodecMode(HANDLE_DRC_DECODER hDrcDec,
                        const DRC_DEC_CODEC_MODE codecMode);

DRC_DEC_ERROR
FDK_drcDec_Init(HANDLE_DRC_DECODER hDrcDec, const int frameSize,
                const int sampleRate, const int baseChannelCount);

DRC_DEC_ERROR
FDK_drcDec_Close(HANDLE_DRC_DECODER* phDrcDec);

DRC_DEC_ERROR
FDK_drcDec_SetParam(HANDLE_DRC_DECODER hDrcDec,
                    const DRC_DEC_USERPARAM requestType,
                    const FIXP_DBL requestValue);

LONG FDK_drcDec_GetParam(HANDLE_DRC_DECODER hDrcDec,
                         const DRC_DEC_USERPARAM requestType);

DRC_DEC_ERROR
FDK_drcDec_SetInterfaceParameters(HANDLE_DRC_DECODER hDrcDec,
                                  HANDLE_UNI_DRC_INTERFACE uniDrcInterface);

DRC_DEC_ERROR
FDK_drcDec_SetSelectionProcessMpeghParameters_simple(
    HANDLE_DRC_DECODER hDrcDec, const int groupPresetIdRequested,
    const int numGroupIdsRequested, const int* groupIdsRequested);

DRC_DEC_ERROR
FDK_drcDec_SetDownmixInstructions(HANDLE_DRC_DECODER hDrcDec,
                                  const int numDowmixId, const int* downmixId,
                                  const int* targetLayout,
                                  const int* targetChannelCount);

void FDK_drcDec_SetSelectionProcessOutput(
    HANDLE_DRC_DECODER hDrcDec, HANDLE_SEL_PROC_OUTPUT hSelProcOutput);

HANDLE_SEL_PROC_OUTPUT
FDK_drcDec_GetSelectionProcessOutput(HANDLE_DRC_DECODER hDrcDec);

LONG
FDK_drcDec_GetGroupLoudness(HANDLE_SEL_PROC_OUTPUT hSelProcOutput,
                            const int groupID, int* groupLoudnessAvailable);

void FDK_drcDec_SetChannelGains(HANDLE_DRC_DECODER hDrcDec,
                                const int numChannels, const int frameSize,
                                FIXP_DBL* channelGainDb, FIXP_DBL* audioBuffer,
                                const int audioBufferChannelOffset);

DRC_DEC_ERROR
FDK_drcDec_ReadUniDrcConfig(HANDLE_DRC_DECODER hDrcDec,
                            HANDLE_FDK_BITSTREAM hBitstream);

DRC_DEC_ERROR
FDK_drcDec_ReadLoudnessInfoSet(HANDLE_DRC_DECODER hDrcDec,
                               HANDLE_FDK_BITSTREAM hBitstream);

DRC_DEC_ERROR
FDK_drcDec_ReadLoudnessBox(HANDLE_DRC_DECODER hDrcDec,
                           HANDLE_FDK_BITSTREAM hBitstream);

DRC_DEC_ERROR
FDK_drcDec_ReadDownmixInstructions_Box(HANDLE_DRC_DECODER hDrcDec,
                                       HANDLE_FDK_BITSTREAM hBitstream);

DRC_DEC_ERROR
FDK_drcDec_ReadUniDrcInstructions_Box(HANDLE_DRC_DECODER hDrcDec,
                                      HANDLE_FDK_BITSTREAM hBitstream);

DRC_DEC_ERROR
FDK_drcDec_ReadUniDrcCoefficients_Box(HANDLE_DRC_DECODER hDrcDec,
                                      HANDLE_FDK_BITSTREAM hBitstream);

DRC_DEC_ERROR
FDK_drcDec_ReadUniDrcGain(HANDLE_DRC_DECODER hDrcDec,
                          HANDLE_FDK_BITSTREAM hBitstream);

DRC_DEC_ERROR
FDK_drcDec_ReadUniDrc(HANDLE_DRC_DECODER hDrcDec,
                      HANDLE_FDK_BITSTREAM hBitstream);

DRC_DEC_ERROR
FDK_drcDec_Preprocess(HANDLE_DRC_DECODER hDrcDec);

DRC_DEC_ERROR
FDK_drcDec_ProcessTime(HANDLE_DRC_DECODER hDrcDec, const int delaySamples,
                       const DRC_DEC_LOCATION drcLocation,
                       const int channelOffset, const int drcChannelOffset,
                       const int numChannelsProcessed, FIXP_DBL* realBuffer,
                       const int timeDataChannelOffset);

DRC_DEC_ERROR
FDK_drcDec_ProcessFreq(HANDLE_DRC_DECODER hDrcDec, const int delaySamples,
                       const DRC_DEC_LOCATION drcLocation,
                       const int channelOffset, const int drcChannelOffset,
                       const int numChannelsProcessed,
                       const int processSingleTimeslot, FIXP_DBL** realBuffer,
                       FIXP_DBL** imagBuffer);

DRC_DEC_ERROR
FDK_drcDec_ApplyDownmix(HANDLE_DRC_DECODER hDrcDec, int* reverseInChannelMap,
                        int* reverseOutChannelMap, FIXP_DBL* realBuffer,
                        int* pNChannels);

DRC_DEC_ERROR
FDK_drcDec_GetLibInfo(LIB_INFO* info);

#ifdef __cplusplus
}
#endif
#endif
