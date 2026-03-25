/************************* MPEG-D DRC decoder library **************************

   Author(s):   Andreas Hoelzer

   Description: DRC Set Selection

*******************************************************************************/

#ifndef DRCDEC_SELECTIONPROCESS_H
#define DRCDEC_SELECTIONPROCESS_H

#include "drcDec_types.h"
#include "drcDecoder.h"

typedef struct s_drcdec_selection_process* HANDLE_DRC_SELECTION_PROCESS;

#define UNDEFINED_LOUDNESS_VALUE (FIXP_DBL)(MAXVAL_DBL - 1)

typedef enum {
  DRCDEC_SELECTION_PROCESS_NO_ERROR = 0,

  DRCDEC_SELECTION_PROCESS_WARNING = -1000,

  DRCDEC_SELECTION_PROCESS_NOT_OK = -2000,
  DRCDEC_SELECTION_PROCESS_OUTOFMEMORY,
  DRCDEC_SELECTION_PROCESS_INVALID_HANDLE,
  DRCDEC_SELECTION_PROCESS_NOT_SUPPORTED,
  DRCDEC_SELECTION_PROCESS_INVALID_PARAM,
  DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE

} DRCDEC_SELECTION_PROCESS_RETURN;

typedef enum {
  SEL_PROC_TEST_TIME_DOMAIN = -100,
  SEL_PROC_TEST_QMF_DOMAIN,
  SEL_PROC_TEST_STFT_DOMAIN,

  SEL_PROC_CODEC_MODE_UNDEFINED = -1,
  SEL_PROC_MPEG_4_AAC,
  SEL_PROC_MPEG_D_USAC,
  SEL_PROC_MPEG_H_3DA
} SEL_PROC_CODEC_MODE;

typedef enum {
  SEL_PROC_LOUDNESS_NORMALIZATION_ON,
  SEL_PROC_DYNAMIC_RANGE_CONTROL_ON,
  SEL_PROC_TARGET_LOUDNESS,
  SEL_PROC_EFFECT_TYPE,
  SEL_PROC_EFFECT_TYPE_FALLBACK_CODE,
  SEL_PROC_LOUDNESS_MEASUREMENT_METHOD,
  SEL_PROC_ALBUM_MODE,
  SEL_PROC_DOWNMIX_ID,
  SEL_PROC_TARGET_LAYOUT,
  SEL_PROC_TARGET_CHANNEL_COUNT,
  SEL_PROC_BASE_CHANNEL_COUNT,
  SEL_PROC_SAMPLE_RATE,
  SEL_PROC_BOOST,
  SEL_PROC_COMPRESS
} SEL_PROC_USER_PARAM;

typedef struct s_selection_process_output {
  FIXP_DBL outputPeakLevelDb;
  FIXP_DBL loudnessNormalizationGainDb;
  FIXP_DBL outputLoudness;

  UCHAR numSelectedDrcSets;
  SCHAR selectedDrcSetIds[MAX_ACTIVE_DRCS];
  UCHAR selectedDownmixIds[MAX_ACTIVE_DRCS];

  UCHAR activeDownmixId;
  UCHAR baseChannelCount;
  UCHAR targetChannelCount;
  SCHAR targetLayout;
  UCHAR downmixMatrixPresent;
  FIXP_DBL downmixMatrix[8][8];

  FIXP_SGL boost;
  FIXP_SGL compress;
  FIXP_DBL mixingLevel;

} SEL_PROC_OUTPUT, *HANDLE_SEL_PROC_OUTPUT;

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_Create(HANDLE_DRC_SELECTION_PROCESS* phInstance);

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_Delete(HANDLE_DRC_SELECTION_PROCESS* phInstance);

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_Init(HANDLE_DRC_SELECTION_PROCESS hInstance);

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_SetCodecMode(HANDLE_DRC_SELECTION_PROCESS hInstance,
                                     const SEL_PROC_CODEC_MODE codecMode);

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_SetParam(HANDLE_DRC_SELECTION_PROCESS hInstance,
                                 const SEL_PROC_USER_PARAM requestType,
                                 FIXP_DBL requestValue, int* pDiff);

FIXP_DBL
drcDec_SelectionProcess_GetParam(HANDLE_DRC_SELECTION_PROCESS hInstance,
                                 const SEL_PROC_USER_PARAM requestType);

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_SetMpeghParams(
    HANDLE_DRC_SELECTION_PROCESS hInstance, const int numGroupIdsRequested,
    const int* groupIdRequested, const int numGroupPresetIdsRequested,
    const int* groupPresetIdRequested,
    const int* numMembersGroupPresetIdsRequested,
    const int groupPresetIdRequestedPreference, int* pDiff);

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_Process(HANDLE_DRC_SELECTION_PROCESS hInstance,
                                HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
                                HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
                                HANDLE_SEL_PROC_OUTPUT hSelProcOutput);

#endif
