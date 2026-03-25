/************************* MPEG-D DRC decoder library **************************

   Author(s):   Andreas Hoelzer

   Description: DRC Set Selection

*******************************************************************************/

#include "drcDec_selectionProcess.h"
#include "drcDec_tools.h"

typedef enum {
  DETR_NONE = 0,
  DETR_NIGHT = 1,
  DETR_NOISY = 2,
  DETR_LIMITED = 3,
  DETR_LOWLEVEL = 4,
  DETR_DIALOG = 5,
  DETR_GENERAL_COMPR = 6,
  DETR_EXPAND = 7,
  DETR_ARTISTIC = 8,
  DETR_COUNT
} DRC_EFFECT_TYPE_REQUEST;

typedef enum {
  DFRT_EFFECT_TYPE,
  DFRT_DYNAMIC_RANGE,
  DFRT_DRC_CHARACTERISTIC
} DRC_FEATURE_REQUEST_TYPE;

typedef enum {
  MDR_DEFAULT = 0,
  MDR_PROGRAM_LOUDNESS = 1,
  MDR_ANCHOR_LOUDNESS = 2
} METHOD_DEFINITION_REQUEST;

typedef enum {
  MSR_DEFAULT = 0,
  MSR_BS_1770_4 = 1,
  MSR_USER = 2,
  MSR_EXPERT_PANEL = 3,
  MSR_RESERVED_A = 4,
  MSR_RESERVED_B = 5,
  MSR_RESERVED_C = 6,
  MSR_RESERVED_D = 7,
  MSR_RESERVED_E = 8
} MEASUREMENT_SYSTEM_REQUEST;

typedef enum {
  LPR_DEFAULT = 0,
  LPR_OFF = 1,
  LPR_HIGHPASS = 2
} LOUDNESS_PREPROCESSING_REQUEST;

typedef enum {
  DRMRT_SHORT_TERM_LOUDNESS_TO_AVG = 0,
  DRMRT_MOMENTARY_LOUDNESS_TO_AVG = 1,
  DRMRT_TOP_OF_LOUDNESS_RANGE_TO_AVG = 2
} DYN_RANGE_MEASUREMENT_REQUEST_TYPE;

typedef enum {
  TCRT_DOWNMIX_ID = 0,
  TCRT_TARGET_LAYOUT = 1,
  TCRT_TARGET_CHANNEL_COUNT = 2
} TARGET_CONFIG_REQUEST_TYPE;

typedef shouldBeUnion {
  struct {
    UCHAR numRequests;
    UCHAR numRequestsDesired;
    DRC_EFFECT_TYPE_REQUEST request[MAX_REQUESTS_DRC_EFFECT_TYPE];
  } drcEffectType;
  struct {
    DYN_RANGE_MEASUREMENT_REQUEST_TYPE measurementRequestType;
    UCHAR requestedIsRange;
    FIXP_DBL requestValue;
    FIXP_DBL requestValueMin;
    FIXP_DBL requestValueMax;
  } dynamicRange;
  UCHAR drcCharacteristic;
}
DRC_FEATURE_REQUEST;

typedef struct {
  SCHAR baseChannelCount;
  SCHAR baseLayout;
  TARGET_CONFIG_REQUEST_TYPE targetConfigRequestType;
  UCHAR numDownmixIdRequests;
  UCHAR downmixIdRequested[MAX_REQUESTS_DOWNMIX_ID];
  UCHAR targetLayoutRequested;
  UCHAR targetChannelCountRequested;
  LONG audioSampleRate;

  UCHAR loudnessNormalizationOn;
  FIXP_DBL targetLoudness;
  UCHAR albumMode;
  UCHAR peakLimiterPresent;
  UCHAR loudnessDeviationMax;
  METHOD_DEFINITION_REQUEST loudnessMeasurementMethod;
  MEASUREMENT_SYSTEM_REQUEST loudnessMeasurementSystem;
  LOUDNESS_PREPROCESSING_REQUEST loudnessMeasurementPreProc;
  LONG deviceCutOffFrequency;
  FIXP_DBL loudnessNormalizationGainDbMax;
  FIXP_DBL loudnessNormalizationGainModificationDb;
  FIXP_DBL outputPeakLevelMax;

  UCHAR dynamicRangeControlOn;
  UCHAR numDrcFeatureRequests;
  DRC_FEATURE_REQUEST_TYPE drcFeatureRequestType[MAX_REQUESTS_DRC_FEATURE];
  DRC_FEATURE_REQUEST drcFeatureRequest[MAX_REQUESTS_DRC_FEATURE];

  FIXP_SGL boost;
  FIXP_SGL compress;
  UCHAR drcCharacteristicTarget;
} SEL_PROC_INPUT, *HANDLE_SEL_PROC_INPUT;

static DRC_EFFECT_TYPE_REQUEST fallbackEffectTypeRequests[6][5] = {
    {DETR_GENERAL_COMPR, DETR_NOISY, DETR_LIMITED, DETR_LOWLEVEL,
     DETR_DIALOG},
    {DETR_GENERAL_COMPR, DETR_NIGHT, DETR_LIMITED, DETR_LOWLEVEL, DETR_DIALOG},
    {DETR_GENERAL_COMPR, DETR_NIGHT, DETR_NOISY, DETR_LOWLEVEL, DETR_DIALOG},
    {DETR_GENERAL_COMPR, DETR_NOISY, DETR_NIGHT, DETR_LIMITED, DETR_DIALOG},
    {DETR_GENERAL_COMPR, DETR_NIGHT, DETR_NOISY, DETR_LIMITED, DETR_LOWLEVEL},
    {DETR_NIGHT, DETR_NOISY, DETR_LIMITED, DETR_LOWLEVEL, DETR_DIALOG}};

typedef struct {
  UCHAR selectionFlag;
  UCHAR downmixIdRequestIndex;
  FIXP_DBL outputPeakLevel;
  FIXP_DBL loudnessNormalizationGainDbAdjusted;
  FIXP_DBL outputLoudness;
  DRC_INSTRUCTIONS_UNI_DRC* pInst;

} DRCDEC_SELECTION_DATA;

typedef struct {
  UCHAR numData;
  DRCDEC_SELECTION_DATA data[(12 + 1 + 6)];

} DRCDEC_SELECTION;

static int _isError(int x) {
  if (x < DRCDEC_SELECTION_PROCESS_WARNING) {
    return 1;
  }

  return 0;
}

static inline int _compAssign(UCHAR* dest, const UCHAR src) {
  int diff = 0;
  if (*dest != src) diff = 1;
  *dest = src;
  return diff;
}

static inline int _compAssign(SCHAR* dest, const SCHAR src) {
  int diff = 0;
  if (*dest != src) diff = 1;
  *dest = src;
  return diff;
}

static inline int _compAssign(FIXP_DBL* dest, const FIXP_DBL src) {
  int diff = 0;
  if (*dest != src) diff = 1;
  *dest = src;
  return diff;
}

static inline int _compAssign(FIXP_SGL* dest, const FIXP_SGL src) {
  int diff = 0;
  if (*dest != src) diff = 1;
  *dest = src;
  return diff;
}

static inline int _compAssign(TARGET_CONFIG_REQUEST_TYPE* dest, const int src) {
  int diff = 0;
  if (*dest != src) diff = 1;
  *dest = (TARGET_CONFIG_REQUEST_TYPE)src;
  return diff;
}

static inline int _compAssign(METHOD_DEFINITION_REQUEST* dest, const int src) {
  int diff = 0;
  if (*dest != src) diff = 1;
  *dest = (METHOD_DEFINITION_REQUEST)src;
  return diff;
}

static inline int _compAssign(DRC_FEATURE_REQUEST_TYPE* dest, const int src) {
  int diff = 0;
  if (*dest != src) diff = 1;
  *dest = (DRC_FEATURE_REQUEST_TYPE)src;
  return diff;
}

static inline int _compAssign(DRC_EFFECT_TYPE_REQUEST* dest, const int src) {
  int diff = 0;
  if (*dest != src) diff = 1;
  *dest = (DRC_EFFECT_TYPE_REQUEST)src;
  return diff;
}

static DRCDEC_SELECTION_DATA* _drcdec_selection_addNew(
    DRCDEC_SELECTION* pSelection);

static DRCDEC_SELECTION_DATA* _drcdec_selection_add(
    DRCDEC_SELECTION* pSelection, DRCDEC_SELECTION_DATA* pDataIn);

static int _drcdec_selection_clear(DRCDEC_SELECTION* pSelection);

static int _drcdec_selection_getNumber(DRCDEC_SELECTION* pSelection);

static int _drcdec_selection_setNumber(DRCDEC_SELECTION* pSelection, int num);

static DRCDEC_SELECTION_DATA* _drcdec_selection_getAt(
    DRCDEC_SELECTION* pSelection, int at);

static int _swapSelectionAndClear(DRCDEC_SELECTION** ppCandidatesPotential,
                                  DRCDEC_SELECTION** ppCandidatesSelected);

static int _swapSelection(DRCDEC_SELECTION** ppCandidatesPotential,
                          DRCDEC_SELECTION** ppCandidatesSelected);

static DRCDEC_SELECTION_PROCESS_RETURN _initDefaultParams(
    HANDLE_SEL_PROC_INPUT hSelProcInput);

static DRCDEC_SELECTION_PROCESS_RETURN _initCodecModeParams(
    HANDLE_SEL_PROC_INPUT hSelProcInput, const SEL_PROC_CODEC_MODE codecMode);

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetPreSelection(
    SEL_PROC_INPUT* hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
    DRCDEC_SELECTION** ppCandidatesPotential,
    DRCDEC_SELECTION** ppCandidatesSelected, SEL_PROC_CODEC_MODE codecMode);

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetFinalSelection_peakValue0(
    DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected);

static DRCDEC_SELECTION_PROCESS_RETURN _dynamicRangeMeasurement(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, DRC_INSTRUCTIONS_UNI_DRC* pInst,
    UCHAR downmixIdRequested,
    DYN_RANGE_MEASUREMENT_REQUEST_TYPE dynamicRangeMeasurementType,
    int albumMode, int* peakToAveragePresent, FIXP_DBL* peakToAverage);

static DRCDEC_SELECTION_PROCESS_RETURN _channelLayoutToDownmixIdMapping(
    HANDLE_SEL_PROC_INPUT hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig);

static DRCDEC_SELECTION_PROCESS_RETURN _generateVirtualDrcSets(
    HANDLE_SEL_PROC_INPUT hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    SEL_PROC_CODEC_MODE codecMode);

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetRequestSelection(
    SEL_PROC_INPUT* hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
    DRCDEC_SELECTION** ppCandidatesPotential,
    DRCDEC_SELECTION** ppCandidatesSelected);

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetFinalSelection(
    HANDLE_SEL_PROC_INPUT hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    DRCDEC_SELECTION** ppCandidatesPotential,
    DRCDEC_SELECTION** ppCandidatesSelected, SEL_PROC_CODEC_MODE codecMode);

static DRCDEC_SELECTION_PROCESS_RETURN _generateOutputInfo(
    HANDLE_SEL_PROC_INPUT hSelProcInput, HANDLE_SEL_PROC_OUTPUT hSelProcOutput,
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
    DRCDEC_SELECTION_DATA* pSelectionData, SEL_PROC_CODEC_MODE codecMode);

static DRCDEC_SELECTION_PROCESS_RETURN _selectDownmixMatrix(
    HANDLE_SEL_PROC_OUTPUT hSelProcOutput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig);

static DRCDEC_SELECTION_PROCESS_RETURN _getLoudness(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, int albumMode,
    METHOD_DEFINITION_REQUEST measurementMethodRequested,
    MEASUREMENT_SYSTEM_REQUEST measurementSystemRequested,
    FIXP_DBL targetLoudness, int drcSetId, int downmixIdRequested,
    FIXP_DBL* pLoudnessNormalizationGain, FIXP_DBL* pLoudness);

static DRCDEC_SELECTION_PROCESS_RETURN _getMixingLevel(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, int downmixIdRequested,
    int drcSetIdRequested, int albumMode, FIXP_DBL* pMixingLevel);

static DRCDEC_SELECTION_PROCESS_RETURN _getSignalPeakLevel(
    HANDLE_SEL_PROC_INPUT hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, DRC_INSTRUCTIONS_UNI_DRC* pInst,
    int downmixIdRequested, int* explicitPeakInformationPresent,
    FIXP_DBL* signalPeakLevelOut,
    SEL_PROC_CODEC_MODE codecMode);

static DRCDEC_SELECTION_PROCESS_RETURN _extractLoudnessPeakToAverageValue(
    LOUDNESS_INFO* loudnessInfo,
    DYN_RANGE_MEASUREMENT_REQUEST_TYPE dynamicRangeMeasurementType,
    int* pLoudnessPeakToAverageValuePresent,
    FIXP_DBL* pLoudnessPeakToAverageValue);

static DRCDEC_SELECTION_PROCESS_RETURN _selectAlbumLoudness(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
    DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected);

static int _findMethodDefinition(LOUDNESS_INFO* pLoudnessInfo,
                                 int methodDefinition, int startIndex);

struct s_drcdec_selection_process {
  SEL_PROC_CODEC_MODE codecMode;
  SEL_PROC_INPUT selProcInput;
  DRCDEC_SELECTION
  selectionData[2];
};

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_Create(HANDLE_DRC_SELECTION_PROCESS* phInstance) {
  HANDLE_DRC_SELECTION_PROCESS hInstance;
  hInstance = (HANDLE_DRC_SELECTION_PROCESS)FDKcalloc(
      1, sizeof(struct s_drcdec_selection_process));

  if (!hInstance) return DRCDEC_SELECTION_PROCESS_OUTOFMEMORY;

  hInstance->codecMode = SEL_PROC_CODEC_MODE_UNDEFINED;

  *phInstance = hInstance;
  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_Init(HANDLE_DRC_SELECTION_PROCESS hInstance) {
  if (!hInstance) return DRCDEC_SELECTION_PROCESS_NOT_OK;

  _initDefaultParams(&hInstance->selProcInput);
  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_SetCodecMode(HANDLE_DRC_SELECTION_PROCESS hInstance,
                                     const SEL_PROC_CODEC_MODE codecMode) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;

  if (!hInstance) return DRCDEC_SELECTION_PROCESS_NOT_OK;

  switch (codecMode) {
    case SEL_PROC_MPEG_4_AAC:
    case SEL_PROC_MPEG_D_USAC:
    case SEL_PROC_TEST_TIME_DOMAIN:
    case SEL_PROC_TEST_QMF_DOMAIN:
    case SEL_PROC_TEST_STFT_DOMAIN:
      hInstance->codecMode = codecMode;
      break;

    case SEL_PROC_CODEC_MODE_UNDEFINED:
    default:
      return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  retVal = _initCodecModeParams(&(hInstance->selProcInput),
                                hInstance->codecMode = codecMode);

  return retVal;
}

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_SetParam(HANDLE_DRC_SELECTION_PROCESS hInstance,
                                 const SEL_PROC_USER_PARAM requestType,
                                 FIXP_DBL requestValue, int* pDiff) {
  INT requestValueInt = (INT)requestValue;
  int i, diff = 0;
  SEL_PROC_INPUT* pSelProcInput = &(hInstance->selProcInput);

  switch (requestType) {
    case SEL_PROC_LOUDNESS_NORMALIZATION_ON:
      if ((requestValueInt != 0) && (requestValueInt != 1))
        return DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE;
      diff |=
          _compAssign(&pSelProcInput->loudnessNormalizationOn, requestValueInt);
      break;
    case SEL_PROC_TARGET_LOUDNESS:
      if ((requestValue < FL2FXCONST_DBL(-63.0f / (float)(1 << 7))) ||
          (requestValue > (FIXP_DBL)0))
        return DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE;
      if (requestValue >
          FL2FXCONST_DBL(-10.0f /
                         (float)(1 << 7)))
        requestValue = FL2FXCONST_DBL(-10.0f / (float)(1 << 7));
      diff |= _compAssign(&pSelProcInput->targetLoudness, requestValue);
      break;
    case SEL_PROC_EFFECT_TYPE:
      if ((requestValueInt < -1) || (requestValueInt >= DETR_COUNT))
        return DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE;
      if (requestValueInt == -1) {
        diff |= _compAssign(&pSelProcInput->dynamicRangeControlOn, 0);
      } else if (requestValueInt == DETR_NONE) {
        diff |= _compAssign(&pSelProcInput->dynamicRangeControlOn, 1);
        diff |= _compAssign(&pSelProcInput->numDrcFeatureRequests, 0);
      } else {
        diff |= _compAssign(&pSelProcInput->dynamicRangeControlOn, 1);
        diff |= _compAssign(&pSelProcInput->numDrcFeatureRequests, 1);
        diff |= _compAssign(&pSelProcInput->drcFeatureRequestType[0],
                            DFRT_EFFECT_TYPE);
        diff |= _compAssign(&pSelProcInput->drcFeatureRequest[0]
                                 .drcEffectType.numRequestsDesired,
                            1);
        diff |= _compAssign(
            &pSelProcInput->drcFeatureRequest[0].drcEffectType.request[0],
            requestValueInt);
        if ((requestValueInt > DETR_NONE) &&
            (requestValueInt <= DETR_GENERAL_COMPR)) {
          for (i = 0; i < 5; i++) {
            diff |=
                _compAssign(&pSelProcInput->drcFeatureRequest[0]
                                 .drcEffectType.request[i + 1],
                            fallbackEffectTypeRequests[requestValueInt - 1][i]);
          }
          diff |= _compAssign(
              &pSelProcInput->drcFeatureRequest[0].drcEffectType.numRequests,
              6);
        } else {
          diff |= _compAssign(
              &pSelProcInput->drcFeatureRequest[0].drcEffectType.numRequests,
              1);
        }
      }
      break;
    case SEL_PROC_LOUDNESS_MEASUREMENT_METHOD:
      if ((requestValueInt < 0) || (requestValueInt > 2))
        return DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE;
      diff |= _compAssign(&pSelProcInput->loudnessMeasurementMethod,
                          requestValueInt);
      break;
    case SEL_PROC_ALBUM_MODE:
      if ((requestValueInt < 0) || (requestValueInt > 1))
        return DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE;
      diff |= _compAssign(&pSelProcInput->albumMode, requestValueInt);
      break;
    case SEL_PROC_DOWNMIX_ID:
      diff |=
          _compAssign(&pSelProcInput->targetConfigRequestType, TCRT_DOWNMIX_ID);
      if (requestValueInt < 0) {
        diff |= _compAssign(&pSelProcInput->numDownmixIdRequests, 0);
      } else {
        diff |= _compAssign(&pSelProcInput->numDownmixIdRequests, 1);
        diff |=
            _compAssign(&pSelProcInput->downmixIdRequested[0], requestValueInt);
      }
      break;
    case SEL_PROC_TARGET_LAYOUT:
      if ((requestValueInt < 1) || (requestValueInt > 63))
        return DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE;
      diff |= _compAssign(&pSelProcInput->targetConfigRequestType,
                          TCRT_TARGET_LAYOUT);
      diff |=
          _compAssign(&pSelProcInput->targetLayoutRequested, requestValueInt);
      break;
    case SEL_PROC_TARGET_CHANNEL_COUNT:
      if ((requestValueInt < 1) || (requestValueInt > 8))
        return DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE;
      diff |= _compAssign(&pSelProcInput->targetConfigRequestType,
                          TCRT_TARGET_CHANNEL_COUNT);
      diff |= _compAssign(&pSelProcInput->targetChannelCountRequested,
                          requestValueInt);
      break;
    case SEL_PROC_BASE_CHANNEL_COUNT:
      if (requestValueInt < 0)
        return DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE;
      diff |= _compAssign(&pSelProcInput->baseChannelCount, requestValueInt);
      break;
    case SEL_PROC_SAMPLE_RATE:
      if (requestValueInt < 0)
        return DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE;
      diff |= _compAssign(&pSelProcInput->audioSampleRate, requestValueInt);
      break;
    case SEL_PROC_BOOST:
      if ((requestValue < (FIXP_DBL)0) ||
          (requestValue > FL2FXCONST_DBL(1.0f / (float)(1 << 1))))
        return DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE;
      diff |= _compAssign(
          &pSelProcInput->boost,
          FX_DBL2FX_SGL(
              requestValue +
              (FIXP_DBL)(1 << 15)));
      break;
    case SEL_PROC_COMPRESS:
      if ((requestValue < (FIXP_DBL)0) ||
          (requestValue > FL2FXCONST_DBL(1.0f / (float)(1 << 1))))
        return DRCDEC_SELECTION_PROCESS_PARAM_OUT_OF_RANGE;
      diff |= _compAssign(
          &pSelProcInput->compress,
          FX_DBL2FX_SGL(
              requestValue +
              (FIXP_DBL)(1 << 15)));
      break;
    default:
      return DRCDEC_SELECTION_PROCESS_INVALID_PARAM;
  }

  if (pDiff != NULL) {
    *pDiff |= diff;
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

FIXP_DBL
drcDec_SelectionProcess_GetParam(HANDLE_DRC_SELECTION_PROCESS hInstance,
                                 const SEL_PROC_USER_PARAM requestType) {
  SEL_PROC_INPUT* pSelProcInput = &(hInstance->selProcInput);

  switch (requestType) {
    case SEL_PROC_LOUDNESS_NORMALIZATION_ON:
      return (FIXP_DBL)pSelProcInput->loudnessNormalizationOn;
    case SEL_PROC_DYNAMIC_RANGE_CONTROL_ON:
      return (FIXP_DBL)pSelProcInput->dynamicRangeControlOn;
    default:
      return (FIXP_DBL)0;
  }
}

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_Delete(HANDLE_DRC_SELECTION_PROCESS* phInstance) {
  if (phInstance == NULL || *phInstance == NULL)
    return DRCDEC_SELECTION_PROCESS_INVALID_HANDLE;

  FDKfree(*phInstance);
  *phInstance = NULL;
  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

DRCDEC_SELECTION_PROCESS_RETURN
drcDec_SelectionProcess_Process(HANDLE_DRC_SELECTION_PROCESS hInstance,
                                HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
                                HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
                                HANDLE_SEL_PROC_OUTPUT hSelProcOutput) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;
  DRCDEC_SELECTION* pCandidatesSelected;
  DRCDEC_SELECTION* pCandidatesPotential;

  if (hInstance == NULL) return DRCDEC_SELECTION_PROCESS_INVALID_HANDLE;

  pCandidatesSelected = &(hInstance->selectionData[0]);
  pCandidatesPotential = &(hInstance->selectionData[1]);
  _drcdec_selection_setNumber(pCandidatesSelected, 0);
  _drcdec_selection_setNumber(pCandidatesPotential, 0);

  retVal = _generateVirtualDrcSets(&(hInstance->selProcInput), hUniDrcConfig,
                                   hInstance->codecMode);
  if (retVal) return (retVal);

  if (hInstance->selProcInput.baseChannelCount !=
      hUniDrcConfig->channelLayout.baseChannelCount) {
    hInstance->selProcInput.baseChannelCount =
        hUniDrcConfig->channelLayout.baseChannelCount;
  }

  if ((hInstance->selProcInput.targetConfigRequestType != 0) ||
      (hInstance->selProcInput.targetConfigRequestType == 0 &&
       hInstance->selProcInput.numDownmixIdRequests == 0)) {
    retVal = _channelLayoutToDownmixIdMapping(&(hInstance->selProcInput),
                                              hUniDrcConfig);

    if (_isError(retVal)) return (retVal);
  }

  retVal = _drcSetPreSelection(&(hInstance->selProcInput), hUniDrcConfig,
                               hLoudnessInfoSet, &pCandidatesPotential,
                               &pCandidatesSelected, hInstance->codecMode);
  if (retVal) return (retVal);

  if (hInstance->selProcInput.albumMode) {
    _swapSelectionAndClear(&pCandidatesPotential, &pCandidatesSelected);

    retVal = _selectAlbumLoudness(hLoudnessInfoSet, pCandidatesPotential,
                                  pCandidatesSelected);
    if (retVal) return (retVal);

    if (_drcdec_selection_getNumber(pCandidatesSelected) == 0) {
      _swapSelection(&pCandidatesPotential, &pCandidatesSelected);
    }
  }

  _swapSelectionAndClear(&pCandidatesPotential, &pCandidatesSelected);

  retVal = _drcSetRequestSelection(&(hInstance->selProcInput), hUniDrcConfig,
                                   hLoudnessInfoSet, &pCandidatesPotential,
                                   &pCandidatesSelected);
  if (retVal) return (retVal);

  retVal = _drcSetFinalSelection(&(hInstance->selProcInput), hUniDrcConfig,
                                 &pCandidatesPotential, &pCandidatesSelected,
                                 hInstance->codecMode);
  if (retVal) return (retVal);

  retVal = _generateOutputInfo(
      &(hInstance->selProcInput), hSelProcOutput, hUniDrcConfig,
      hLoudnessInfoSet, &(pCandidatesSelected->data[0]), hInstance->codecMode);

  if (_isError(retVal)) return (retVal);

  retVal = _selectDownmixMatrix(hSelProcOutput, hUniDrcConfig);
  if (retVal) return (retVal);

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _initDefaultParams(
    HANDLE_SEL_PROC_INPUT hSelProcInput) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;

  if (hSelProcInput == NULL) return DRCDEC_SELECTION_PROCESS_INVALID_HANDLE;

  hSelProcInput->baseChannelCount = -1;
  hSelProcInput->baseLayout = -1;
  hSelProcInput->targetConfigRequestType = TCRT_DOWNMIX_ID;
  hSelProcInput->numDownmixIdRequests = 0;

  hSelProcInput->albumMode = 0;
  hSelProcInput->peakLimiterPresent = 0;
  hSelProcInput->loudnessNormalizationOn = 1;
  hSelProcInput->targetLoudness = FL2FXCONST_DBL(-24.0f / (float)(1 << 7));
  hSelProcInput->loudnessDeviationMax = DEFAULT_LOUDNESS_DEVIATION_MAX;
  hSelProcInput->loudnessMeasurementMethod = MDR_ANCHOR_LOUDNESS;
  hSelProcInput->loudnessMeasurementSystem = MSR_EXPERT_PANEL;
  hSelProcInput->loudnessMeasurementPreProc = LPR_DEFAULT;
  hSelProcInput->deviceCutOffFrequency = 500;
  hSelProcInput->loudnessNormalizationGainDbMax =
      (FIXP_DBL)MAXVAL_DBL;
  hSelProcInput->loudnessNormalizationGainModificationDb = (FIXP_DBL)0;
  hSelProcInput->outputPeakLevelMax = (FIXP_DBL)0;
  if (hSelProcInput->peakLimiterPresent == 1) {
    hSelProcInput->outputPeakLevelMax = FL2FXCONST_DBL(6.0f / (float)(1 << 7));
  }

  hSelProcInput->dynamicRangeControlOn = 1;
  hSelProcInput->numDrcFeatureRequests = 0;

  hSelProcInput->boost = FL2FXCONST_SGL(1.f / (float)(1 << 1));
  hSelProcInput->compress = FL2FXCONST_SGL(1.f / (float)(1 << 1));
  hSelProcInput->drcCharacteristicTarget = 0;

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _initCodecModeParams(
    HANDLE_SEL_PROC_INPUT hSelProcInput, const SEL_PROC_CODEC_MODE codecMode) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;

  if (hSelProcInput == NULL) return DRCDEC_SELECTION_PROCESS_INVALID_HANDLE;

  switch (codecMode) {
    case SEL_PROC_MPEG_H_3DA:
      hSelProcInput->loudnessDeviationMax = 0;
      hSelProcInput->peakLimiterPresent = 1;
      hSelProcInput->outputPeakLevelMax = (FIXP_DBL)0;
      break;
    case SEL_PROC_MPEG_4_AAC:
    case SEL_PROC_MPEG_D_USAC:
      hSelProcInput->loudnessDeviationMax = DEFAULT_LOUDNESS_DEVIATION_MAX;
      hSelProcInput->peakLimiterPresent = 1;
      hSelProcInput->outputPeakLevelMax =
          FL2FXCONST_DBL(6.0f / (float)(1 << 7));
      break;
    case SEL_PROC_TEST_TIME_DOMAIN:
    case SEL_PROC_TEST_QMF_DOMAIN:
    case SEL_PROC_TEST_STFT_DOMAIN:
      hSelProcInput->loudnessNormalizationOn = 0;
      hSelProcInput->dynamicRangeControlOn = 0;
      break;
    case SEL_PROC_CODEC_MODE_UNDEFINED:
    default:
      hSelProcInput->loudnessDeviationMax = DEFAULT_LOUDNESS_DEVIATION_MAX;
      hSelProcInput->peakLimiterPresent = 0;
  }

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _channelLayoutToDownmixIdMapping(
    HANDLE_SEL_PROC_INPUT hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;

  DOWNMIX_INSTRUCTIONS* pDown = NULL;

  int i;

  hSelProcInput->numDownmixIdRequests = 0;

  switch (hSelProcInput->targetConfigRequestType) {
    case TCRT_DOWNMIX_ID:
      if (hSelProcInput->numDownmixIdRequests == 0) {
        hSelProcInput->downmixIdRequested[0] = 0;
        hSelProcInput->numDownmixIdRequests = 1;
      }

      break;

    case TCRT_TARGET_LAYOUT:
      if (hSelProcInput->targetLayoutRequested == hSelProcInput->baseLayout) {
        hSelProcInput->downmixIdRequested[0] = 0;
        hSelProcInput->numDownmixIdRequests = 1;
      }

      if (hSelProcInput->numDownmixIdRequests == 0) {
        for (i = 0; i < hUniDrcConfig->downmixInstructionsCount; i++) {
          pDown = &(hUniDrcConfig->downmixInstructions[i]);

          if (hSelProcInput->targetLayoutRequested == pDown->targetLayout) {
            hSelProcInput
                ->downmixIdRequested[hSelProcInput->numDownmixIdRequests] =
                pDown->downmixId;
            hSelProcInput->numDownmixIdRequests++;
          }
        }
      }

      if (hSelProcInput->baseLayout == -1) {
        retVal = DRCDEC_SELECTION_PROCESS_WARNING;
      }

      if (hSelProcInput->numDownmixIdRequests == 0) {
        hSelProcInput->downmixIdRequested[0] = 0;
        hSelProcInput->numDownmixIdRequests = 1;
        retVal = DRCDEC_SELECTION_PROCESS_WARNING;
      }

      break;

    case TCRT_TARGET_CHANNEL_COUNT:
      if (hSelProcInput->targetChannelCountRequested ==
          hSelProcInput->baseChannelCount) {
        hSelProcInput->downmixIdRequested[0] = 0;
        hSelProcInput->numDownmixIdRequests = 1;
      }

      if (hSelProcInput->numDownmixIdRequests == 0) {
        for (i = 0; i < hUniDrcConfig->downmixInstructionsCount; i++) {
          pDown = &(hUniDrcConfig->downmixInstructions[i]);

          if (hSelProcInput->targetChannelCountRequested ==
              pDown->targetChannelCount) {
            hSelProcInput
                ->downmixIdRequested[hSelProcInput->numDownmixIdRequests] =
                pDown->downmixId;
            hSelProcInput->numDownmixIdRequests++;
          }
        }
      }

      if (hSelProcInput->baseChannelCount == -1) {
        retVal = DRCDEC_SELECTION_PROCESS_WARNING;
      }

      if (hSelProcInput->numDownmixIdRequests == 0) {
        retVal = DRCDEC_SELECTION_PROCESS_WARNING;
        hSelProcInput->downmixIdRequested[0] = 0;
        hSelProcInput->numDownmixIdRequests = 1;
      }

      break;

    default:
      return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _preSelectionRequirement123(
    int nRequestedDownmixId, DRC_INSTRUCTIONS_UNI_DRC* pDrcInstructionUniDrc,
    int* pMatchFound) {
  int i;
  *pMatchFound = 0;

  for (i = 0; i < pDrcInstructionUniDrc->downmixIdCount; i++) {
    if ((pDrcInstructionUniDrc->downmixId[i] == nRequestedDownmixId) ||
        (pDrcInstructionUniDrc->downmixId[i] == DOWNMIX_ID_ANY_DOWNMIX) ||
        ((pDrcInstructionUniDrc->downmixId[i] == DOWNMIX_ID_BASE_LAYOUT) &&
         (pDrcInstructionUniDrc->drcSetId > 0))) {
      *pMatchFound = 1;
      break;
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _preSelectionRequirement4(
    DRC_INSTRUCTIONS_UNI_DRC* pDrcInstruction, int nDynamicRangeControlOn,
    int* pMatchFound) {
  *pMatchFound = 0;

  if (nDynamicRangeControlOn == 1) {
    if ((pDrcInstruction->drcSetEffect != EB_FADE) &&
        (pDrcInstruction->drcSetEffect != EB_DUCK_OTHER) &&
        (pDrcInstruction->drcSetEffect != EB_DUCK_SELF) &&
        (pDrcInstruction->drcSetEffect != 0 || pDrcInstruction->drcSetId < 0)) {
      *pMatchFound = 1;
    }
  } else {
    *pMatchFound = 1;
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _preSelectionRequirement5(
    DRC_INSTRUCTIONS_UNI_DRC* pDrcInstructionUniDrc,
    DRC_COEFFICIENTS_UNI_DRC* pCoef, int* pMatchFound) {
  int b, i;

  *pMatchFound = 1;

  if (pDrcInstructionUniDrc->drcSetId < 0)
  {
    return DRCDEC_SELECTION_PROCESS_NO_ERROR;
  }

  if (pCoef == NULL)
  {
    *pMatchFound = 0;
    return DRCDEC_SELECTION_PROCESS_NO_ERROR;
  }

  if (pCoef->drcLocation !=
      pDrcInstructionUniDrc
          ->drcLocation)
  {
    *pMatchFound = 0;
    return DRCDEC_SELECTION_PROCESS_NO_ERROR;
  }

  for (i = 0; i < pDrcInstructionUniDrc->nDrcChannelGroups; i++) {
    int indexDrcCoeff = pDrcInstructionUniDrc->gainSetIndexForChannelGroup[i];
    int bandCount = 0;

    if (indexDrcCoeff >= 12) {
      *pMatchFound = 0;
      return DRCDEC_SELECTION_PROCESS_NO_ERROR;
    }

    if (indexDrcCoeff > pCoef->gainSetCount - 1)
    {
      continue;
    }

    GAIN_SET* gainSet = &(pCoef->gainSet[indexDrcCoeff]);
    bandCount = gainSet->bandCount;

    if (bandCount > 4) {
      *pMatchFound = 0;
    }

    for (b = 0; b < bandCount; b++) {
      if ((gainSet->gainSequenceIndex[b] >= 12) ||
          (gainSet->gainSequenceIndex[b] >= pCoef->gainSequenceCount)) {
        *pMatchFound = 0;
        return DRCDEC_SELECTION_PROCESS_NO_ERROR;
      }
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _preSelectionRequirement6(
    DRC_INSTRUCTIONS_UNI_DRC* pDrcInstructionUniDrc, int* pMatchFound) {
  *pMatchFound = 0;

  if (((pDrcInstructionUniDrc->dependsOnDrcSetPresent == 0) &&
       (pDrcInstructionUniDrc->noIndependentUse == 0)) ||
      (pDrcInstructionUniDrc->dependsOnDrcSetPresent == 1)) {
    *pMatchFound = 1;
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _preSelectionRequirement7(
    DRC_INSTRUCTIONS_UNI_DRC* pDrcInstructionUniDrc, int* pMatchFound) {
  *pMatchFound = 1;

  if (pDrcInstructionUniDrc->requiresEq) {
    *pMatchFound = 0;
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static void _setSelectionDataInfo(
    DRCDEC_SELECTION_DATA* pData, FIXP_DBL loudness,
    FIXP_DBL loudnessNormalizationGainDb,
    FIXP_DBL loudnessNormalizationGainDbMax,
    FIXP_DBL loudnessDeviationMax,
    FIXP_DBL signalPeakLevel,
    FIXP_DBL outputPeakLevelMax,
    int applyAdjustment) {
  FIXP_DBL adjustment = 0;

  loudness >>= 1;
  loudnessNormalizationGainDb >>= 1;
  loudnessNormalizationGainDbMax >>= 1;
  loudnessDeviationMax >>= 1;
  signalPeakLevel >>= 1;
  outputPeakLevelMax >>= 1;

  if (applyAdjustment) {
    adjustment =
        fMax((FIXP_DBL)0, signalPeakLevel + loudnessNormalizationGainDb -
                              outputPeakLevelMax);
    adjustment = fMin(adjustment, fMax((FIXP_DBL)0, loudnessDeviationMax));
  }

  pData->loudnessNormalizationGainDbAdjusted = fMin(
      loudnessNormalizationGainDb - adjustment, loudnessNormalizationGainDbMax);
  pData->outputLoudness = loudness + pData->loudnessNormalizationGainDbAdjusted;
  pData->outputPeakLevel =
      signalPeakLevel + pData->loudnessNormalizationGainDbAdjusted;

  pData->loudnessNormalizationGainDbAdjusted = SATURATE_LEFT_SHIFT(
      pData->loudnessNormalizationGainDbAdjusted, 1, DFRACT_BITS);
  pData->outputLoudness =
      SATURATE_LEFT_SHIFT(pData->outputLoudness, 1, DFRACT_BITS);
  pData->outputPeakLevel =
      SATURATE_LEFT_SHIFT(pData->outputPeakLevel, 1, DFRACT_BITS);
}

static int _targetLoudnessInRange(
    DRC_INSTRUCTIONS_UNI_DRC* pDrcInstructionUniDrc, FIXP_DBL targetLoudness) {
  int retVal = 0;

  FIXP_DBL drcSetTargetLoudnessValueUpper =
      ((FIXP_DBL)pDrcInstructionUniDrc->drcSetTargetLoudnessValueUpper)
      << (DFRACT_BITS - 1 - 7);
  FIXP_DBL drcSetTargetLoudnessValueLower =
      ((FIXP_DBL)pDrcInstructionUniDrc->drcSetTargetLoudnessValueLower)
      << (DFRACT_BITS - 1 - 7);

  if (pDrcInstructionUniDrc->drcSetTargetLoudnessPresent &&
      drcSetTargetLoudnessValueUpper >= targetLoudness &&
      drcSetTargetLoudnessValueLower < targetLoudness) {
    retVal = 1;
  }

  return retVal;
}

static int _drcSetIsUsable(HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
                           DRC_INSTRUCTIONS_UNI_DRC* pInst) {
  int usable = 0;
  DRC_COEFFICIENTS_UNI_DRC* pCoef =
      selectDrcCoefficients(hUniDrcConfig, LOCATION_SELECTED);

  if (selectDrcInstructions(hUniDrcConfig, pInst->drcSetId) != pInst) return 0;
  _preSelectionRequirement5(pInst, pCoef, &usable);
  return usable;
}

static DRCDEC_SELECTION_PROCESS_RETURN _preSelectionRequirement8(
    SEL_PROC_INPUT* hSelProcInput, int downmixIdIndex,
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
    DRC_INSTRUCTIONS_UNI_DRC* pDrcInstructionUniDrc,
    DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected, SEL_PROC_CODEC_MODE codecMode) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;
  int explicitPeakInformationPresent;
  FIXP_DBL signalPeakLevel;
  int addToCandidate = 0;

  FIXP_DBL loudnessNormalizationGainDb;
  FIXP_DBL loudness;

  FIXP_DBL loudnessDeviationMax =
      ((FIXP_DBL)hSelProcInput->loudnessDeviationMax) << (DFRACT_BITS - 1 - 7);

  {
    retVal = _getLoudness(hLoudnessInfoSet, hSelProcInput->albumMode,
                          hSelProcInput->loudnessMeasurementMethod,
                          hSelProcInput->loudnessMeasurementSystem,
                          hSelProcInput->targetLoudness,
                          pDrcInstructionUniDrc->drcSetId,
                          hSelProcInput->downmixIdRequested[downmixIdIndex],
                          &loudnessNormalizationGainDb, &loudness);
    if (retVal) return (retVal);
  }

  if (!hSelProcInput->loudnessNormalizationOn) {
    loudnessNormalizationGainDb = (FIXP_DBL)0;
  }

  retVal = _getSignalPeakLevel(
      hSelProcInput, hUniDrcConfig, hLoudnessInfoSet, pDrcInstructionUniDrc,
      hSelProcInput->downmixIdRequested[downmixIdIndex],
      &explicitPeakInformationPresent, &signalPeakLevel, codecMode

  );
  if (retVal) return (retVal);

  if (hSelProcInput->dynamicRangeControlOn) {
    if (explicitPeakInformationPresent == 0) {
      if (pDrcInstructionUniDrc->drcSetTargetLoudnessPresent &&
          ((hSelProcInput->loudnessNormalizationOn &&
            _targetLoudnessInRange(pDrcInstructionUniDrc,
                                   hSelProcInput->targetLoudness)) ||
           !hSelProcInput->loudnessNormalizationOn)) {
        DRCDEC_SELECTION_DATA* pData =
            _drcdec_selection_addNew(pCandidatesSelected);
        if (pData == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

        _setSelectionDataInfo(pData, loudness, loudnessNormalizationGainDb,
                              hSelProcInput->loudnessNormalizationGainDbMax,
                              loudnessDeviationMax, signalPeakLevel,
                              hSelProcInput->outputPeakLevelMax, 0);
        pData->downmixIdRequestIndex = downmixIdIndex;
        pData->pInst = pDrcInstructionUniDrc;
        pData->selectionFlag =
            1;

        if (hSelProcInput->loudnessNormalizationOn) {
          pData->outputPeakLevel =
              hSelProcInput->targetLoudness -
              (((FIXP_DBL)pData->pInst->drcSetTargetLoudnessValueUpper)
               << (DFRACT_BITS - 1 - 7));
        } else {
          pData->outputPeakLevel = (FIXP_DBL)0;
        }
      } else {
        if ((!hSelProcInput->loudnessNormalizationOn) ||
            (!pDrcInstructionUniDrc->drcSetTargetLoudnessPresent) ||
            (hSelProcInput->loudnessNormalizationOn &&
             _targetLoudnessInRange(pDrcInstructionUniDrc,
                                    hSelProcInput->targetLoudness))) {
          addToCandidate = 1;
        }
      }
    } else {
      addToCandidate = 1;
    }

    if (addToCandidate) {
      DRCDEC_SELECTION_DATA* pData =
          _drcdec_selection_addNew(pCandidatesPotential);
      if (pData == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

      _setSelectionDataInfo(pData, loudness, loudnessNormalizationGainDb,
                            hSelProcInput->loudnessNormalizationGainDbMax,
                            loudnessDeviationMax, signalPeakLevel,
                            hSelProcInput->outputPeakLevelMax, 0);
      pData->downmixIdRequestIndex = downmixIdIndex;
      pData->pInst = pDrcInstructionUniDrc;
      pData->selectionFlag = 0;
    }
  } else {
    if (pDrcInstructionUniDrc->drcSetId < 0) {
      DRCDEC_SELECTION_DATA* pData =
          _drcdec_selection_addNew(pCandidatesSelected);
      if (pData == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

      _setSelectionDataInfo(pData, loudness, loudnessNormalizationGainDb,
                            hSelProcInput->loudnessNormalizationGainDbMax,
                            loudnessDeviationMax, signalPeakLevel,
                            hSelProcInput->outputPeakLevelMax, 1);

      pData->downmixIdRequestIndex = downmixIdIndex;
      pData->pInst = pDrcInstructionUniDrc;
      pData->selectionFlag = 0;
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _preSelectionRequirement9(
    SEL_PROC_INPUT* hSelProcInput, DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected) {
  int i;

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    DRCDEC_SELECTION_DATA* pCandidate =
        _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    if (pCandidate->outputPeakLevel <= hSelProcInput->outputPeakLevelMax) {
      if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
        return DRCDEC_SELECTION_PROCESS_NOT_OK;
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetPreSelectionSingleInstruction(
    SEL_PROC_INPUT* hSelProcInput, int downmixIdIndex,
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
    DRC_INSTRUCTIONS_UNI_DRC* pDrcInstructionUniDrc,
    DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected, SEL_PROC_CODEC_MODE codecMode) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;
  int matchFound = 0;
  DRC_COEFFICIENTS_UNI_DRC* pCoef =
      selectDrcCoefficients(hUniDrcConfig, LOCATION_SELECTED);

  retVal = _preSelectionRequirement123(
      hSelProcInput->downmixIdRequested[downmixIdIndex], pDrcInstructionUniDrc,
      &matchFound);

  if (!retVal && matchFound)
    retVal = _preSelectionRequirement4(pDrcInstructionUniDrc,
                                       hSelProcInput->dynamicRangeControlOn,
                                       &matchFound);

  if (!retVal && matchFound)
    retVal =
        _preSelectionRequirement5(pDrcInstructionUniDrc, pCoef, &matchFound);

  if (!retVal && matchFound)
    retVal = _preSelectionRequirement6(pDrcInstructionUniDrc, &matchFound);

  if (!retVal && matchFound)
    retVal = _preSelectionRequirement7(pDrcInstructionUniDrc, &matchFound);

  if (!retVal && matchFound)
    retVal = _preSelectionRequirement8(
        hSelProcInput, downmixIdIndex, hUniDrcConfig, hLoudnessInfoSet,
        pDrcInstructionUniDrc, pCandidatesPotential, pCandidatesSelected,
        codecMode);

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetSelectionAddCandidates(
    SEL_PROC_INPUT* hSelProcInput, DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;
  int nHitCount = 0;
  int i;

  DRCDEC_SELECTION_DATA* pCandidate = NULL;
  DRC_INSTRUCTIONS_UNI_DRC* pDrcInstructionUniDrc = NULL;

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    pDrcInstructionUniDrc = pCandidate->pInst;

    if (_targetLoudnessInRange(pDrcInstructionUniDrc,
                               hSelProcInput->targetLoudness)) {
      nHitCount++;
    }
  }

  if (nHitCount != 0) {
    for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
      pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
      if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

      pDrcInstructionUniDrc = pCandidate->pInst;

      if (_targetLoudnessInRange(pDrcInstructionUniDrc,
                                 hSelProcInput->targetLoudness)) {
        if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
          return DRCDEC_SELECTION_PROCESS_NOT_OK;
      }
    }
  } else {
    FIXP_DBL lowestPeakLevel = MAXVAL_DBL;
    FIXP_DBL peakLevel = 0;

    for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
      pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
      if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

      peakLevel = pCandidate->outputPeakLevel;

      if (peakLevel < lowestPeakLevel) {
        lowestPeakLevel = peakLevel;
      }
    }

    for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
      FIXP_DBL loudnessDeviationMax =
          ((FIXP_DBL)hSelProcInput->loudnessDeviationMax)
          << (DFRACT_BITS - 1 - 7);

      pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
      if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

      peakLevel = pCandidate->outputPeakLevel;

      if (peakLevel == lowestPeakLevel ||
          peakLevel <=
              lowestPeakLevel + FL2FXCONST_DBL(1.0f / (float)(1 << 7))) {
        FIXP_DBL adjustment =
            fMax((FIXP_DBL)0, peakLevel - hSelProcInput->outputPeakLevelMax);
        adjustment = fMin(adjustment, fMax((FIXP_DBL)0, loudnessDeviationMax));

        pCandidate->loudnessNormalizationGainDbAdjusted -= adjustment;
        pCandidate->outputPeakLevel -= adjustment;
        pCandidate->outputLoudness -= adjustment;
        if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
          return DRCDEC_SELECTION_PROCESS_NOT_OK;
      }
    }
  }

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _dependentDrcInstruction(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, DRC_INSTRUCTIONS_UNI_DRC* pInst,
    DRC_INSTRUCTIONS_UNI_DRC** ppDrcInstructionsDependent) {
  int i;
  DRC_INSTRUCTIONS_UNI_DRC* pDependentDrc = NULL;

  for (i = 0; i < hUniDrcConfig->drcInstructionsUniDrcCount; i++) {
    pDependentDrc =
        (DRC_INSTRUCTIONS_UNI_DRC*)&(hUniDrcConfig->drcInstructionsUniDrc[i]);

    if (pDependentDrc->drcSetId == pInst->dependsOnDrcSet) {
      break;
    }
  }

  if (i == hUniDrcConfig->drcInstructionsUniDrcCount) {
    return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  if (pDependentDrc->dependsOnDrcSetPresent == 1) {
    return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  *ppDrcInstructionsDependent = pDependentDrc;

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _selectDrcSetEffectNone(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected) {
  int i;

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    DRCDEC_SELECTION_DATA* pCandidate =
        _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    if ((pCandidate->pInst->drcSetEffect & 0xff) == 0) {
      if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
        return DRCDEC_SELECTION_PROCESS_NOT_OK;
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _selectSingleEffectType(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, DRC_EFFECT_TYPE_REQUEST effectType,
    DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected) {
  int i;
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;
  DRC_INSTRUCTIONS_UNI_DRC* pInst;
  DRC_INSTRUCTIONS_UNI_DRC* pDrcInstructionsDependent;

  if (effectType == DETR_NONE) {
    retVal = _selectDrcSetEffectNone(hUniDrcConfig, pCandidatesPotential,
                                     pCandidatesSelected);
    if (retVal) return (retVal);
  } else {
    int effectBitPosition = 1 << (effectType - 1);

    for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
      DRCDEC_SELECTION_DATA* pCandidate =
          _drcdec_selection_getAt(pCandidatesPotential, i);
      if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

      pInst = pCandidate->pInst;

      if (!pInst->dependsOnDrcSetPresent) {
        if ((pInst->drcSetEffect & effectBitPosition)) {
          if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
            return DRCDEC_SELECTION_PROCESS_NOT_OK;
        }
      } else {
        retVal = _dependentDrcInstruction(hUniDrcConfig, pInst,
                                          &pDrcInstructionsDependent);
        if (retVal) return (retVal);

        if (((pInst->drcSetEffect & effectBitPosition)) ||
            ((pDrcInstructionsDependent->drcSetEffect & effectBitPosition))) {
          if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
            return DRCDEC_SELECTION_PROCESS_NOT_OK;
        }
      }
    }
  }

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _selectEffectTypeFeature(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, DRC_FEATURE_REQUEST drcFeatureRequest,
    DRCDEC_SELECTION** ppCandidatesPotential,
    DRCDEC_SELECTION** ppCandidatesSelected) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;
  int i;
  int desiredEffectTypeFound = 0;

  for (i = 0; i < drcFeatureRequest.drcEffectType.numRequestsDesired; i++) {
    retVal = _selectSingleEffectType(
        hUniDrcConfig, drcFeatureRequest.drcEffectType.request[i],
        *ppCandidatesPotential, *ppCandidatesSelected);
    if (retVal) return (retVal);

    if (_drcdec_selection_getNumber(*ppCandidatesSelected)) {
      desiredEffectTypeFound = 1;
      _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
    }
  }

  if (!desiredEffectTypeFound) {
    for (i = drcFeatureRequest.drcEffectType.numRequestsDesired;
         i < drcFeatureRequest.drcEffectType.numRequests; i++) {
      retVal = _selectSingleEffectType(
          hUniDrcConfig, drcFeatureRequest.drcEffectType.request[i],
          *ppCandidatesPotential, *ppCandidatesSelected);
      if (retVal) return (retVal);

      if (_drcdec_selection_getNumber(*ppCandidatesSelected)) {
        _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
        break;
      }
    }
  }

  _swapSelection(ppCandidatesPotential, ppCandidatesSelected);

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _selectDynamicRange(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
    DRC_FEATURE_REQUEST drcFeatureRequest, UCHAR* pDownmixIdRequested,
    int albumMode, DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* ppCandidatesSelected) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;
  int i;
  int peakToAveragePresent;
  FIXP_DBL peakToAverage;

  FIXP_DBL minVal = MAXVAL_DBL;
  FIXP_DBL val = 0;

  int numSelectedCandidates = _drcdec_selection_getNumber(ppCandidatesSelected);

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    DRCDEC_SELECTION_DATA* pCandidate =
        _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    retVal = _dynamicRangeMeasurement(
        hLoudnessInfoSet, pCandidate->pInst,
        pDownmixIdRequested[pCandidate->downmixIdRequestIndex],
        drcFeatureRequest.dynamicRange.measurementRequestType, albumMode,
        &peakToAveragePresent, &peakToAverage);
    if (retVal) return (retVal);

    if (peakToAveragePresent) {
      if (!drcFeatureRequest.dynamicRange.requestedIsRange) {
        val = fAbs(drcFeatureRequest.dynamicRange.requestValue - peakToAverage);

        if (minVal > val) {
          minVal = val;

          _drcdec_selection_setNumber(ppCandidatesSelected,
                                      numSelectedCandidates);
        }
        if (_drcdec_selection_add(ppCandidatesSelected, pCandidate) == NULL)
          return DRCDEC_SELECTION_PROCESS_NOT_OK;
      } else {
        if ((peakToAverage >= drcFeatureRequest.dynamicRange.requestValueMin) &&
            (peakToAverage <= drcFeatureRequest.dynamicRange.requestValueMax)) {
          if (_drcdec_selection_add(ppCandidatesSelected, pCandidate) == NULL)
            return DRCDEC_SELECTION_PROCESS_NOT_OK;
        }
      }
    }
  }

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _selectSingleDrcCharacteristic(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, int requestedDrcCharacteristic,
    DRCDEC_SELECTION** ppCandidatesPotential,
    DRCDEC_SELECTION** ppCandidatesSelected) {
  int i, j, b;
  int hit = 0;

  DRC_INSTRUCTIONS_UNI_DRC* pInst = NULL;
  DRC_COEFFICIENTS_UNI_DRC* pCoef = NULL;
  GAIN_SET* pGainSet = NULL;

  if (requestedDrcCharacteristic < 1) {
    return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  pCoef = selectDrcCoefficients(hUniDrcConfig, LOCATION_SELECTED);

  if (pCoef == NULL)
  {
    return DRCDEC_SELECTION_PROCESS_NO_ERROR;
  }

  for (i = 0; i < _drcdec_selection_getNumber(*ppCandidatesPotential); i++) {
    DRCDEC_SELECTION_DATA* pCandidate =
        _drcdec_selection_getAt(*ppCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    pInst = pCandidate->pInst;

    hit = 0;

    for (j = 0; j < pInst->nDrcChannelGroups; j++) {
      int bandCount = 0;
      int indexDrcCoeff = pInst->gainSetIndexForChannelGroup[j];

      if (indexDrcCoeff > pCoef->gainSetCount - 1)
      {
        return DRCDEC_SELECTION_PROCESS_NO_ERROR;
      }

      pGainSet = &(pCoef->gainSet[indexDrcCoeff]);
      bandCount = pGainSet->bandCount;

      for (b = 0; b < bandCount; b++) {
        if ((pGainSet->drcCharacteristic[b].isCICP) &&
            (pGainSet->drcCharacteristic[b].cicpIndex ==
             requestedDrcCharacteristic)) {
          hit = 1;
          break;
        }
      }

      if (hit) break;
    }

    if (hit) {
      if (_drcdec_selection_add(*ppCandidatesSelected, pCandidate) == NULL)
        return DRCDEC_SELECTION_PROCESS_NOT_OK;
    }
  }

  if (_drcdec_selection_getNumber(*ppCandidatesSelected)) {
    _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _selectDrcCharacteristic(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, int drcCharacteristicRequested,
    DRCDEC_SELECTION** ppCandidatesPotential,
    DRCDEC_SELECTION** ppCandidatesSelected) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;

  const int secondTry[12] = {0, 2, 3, 4, 5, 6, 5, 9, 10, 7, 8, 10};

  retVal = _selectSingleDrcCharacteristic(
      hUniDrcConfig, drcCharacteristicRequested, ppCandidatesPotential,
      ppCandidatesSelected);
  if (retVal) return (retVal);

  if ((drcCharacteristicRequested <= 11) &&
      (_drcdec_selection_getNumber(*ppCandidatesSelected) == 0)) {
    retVal = _selectSingleDrcCharacteristic(
        hUniDrcConfig, secondTry[drcCharacteristicRequested],
        ppCandidatesPotential, ppCandidatesSelected);
    if (retVal) return (retVal);
  }

  if (_drcdec_selection_getNumber(*ppCandidatesSelected) == 0) {
    if ((drcCharacteristicRequested >= 2) &&
        (drcCharacteristicRequested <= 5)) {
      retVal = _selectSingleDrcCharacteristic(
          hUniDrcConfig, drcCharacteristicRequested - 1, ppCandidatesPotential,
          ppCandidatesSelected);
      if (retVal) return (retVal);
    } else if (drcCharacteristicRequested == 11) {
      retVal = _selectSingleDrcCharacteristic(
          hUniDrcConfig, 9, ppCandidatesPotential, ppCandidatesSelected);
      if (retVal) return (retVal);
    }
  }

  _swapSelection(ppCandidatesPotential, ppCandidatesSelected);

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetFinalSelection_peakValue0(
    DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected) {
  int i;

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    DRCDEC_SELECTION_DATA* pCandidate =
        _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    if (pCandidate->outputPeakLevel <= FIXP_DBL(0)) {
      if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
        return DRCDEC_SELECTION_PROCESS_NOT_OK;
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetFinalSelection_downmixId(
    HANDLE_SEL_PROC_INPUT hSelProcInput,
    DRCDEC_SELECTION** ppCandidatesPotential,
    DRCDEC_SELECTION** ppCandidatesSelected) {
  int i, j;
  DRCDEC_SELECTION_DATA* pCandidate = NULL;
  DRC_INSTRUCTIONS_UNI_DRC* pInst = NULL;

  for (i = 0; i < _drcdec_selection_getNumber(*ppCandidatesPotential); i++) {
    pCandidate = _drcdec_selection_getAt(*ppCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    pInst = pCandidate->pInst;

    for (j = 0; j < pInst->downmixIdCount; j++) {
      if (DOWNMIX_ID_BASE_LAYOUT != pInst->downmixId[j] &&
          DOWNMIX_ID_ANY_DOWNMIX != pInst->downmixId[j] &&
          hSelProcInput
                  ->downmixIdRequested[pCandidate->downmixIdRequestIndex] ==
              pInst->downmixId[j]) {
        if (_drcdec_selection_add(*ppCandidatesSelected, pCandidate) == NULL)
          return DRCDEC_SELECTION_PROCESS_NOT_OK;
      }
    }
  }

  if (_drcdec_selection_getNumber(*ppCandidatesSelected) == 0) {
    _swapSelection(ppCandidatesPotential, ppCandidatesSelected);
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static int _crossSum(int value) {
  int sum = 0;

  while (value != 0) {
    if ((value & 1) == 1) {
      sum++;
    }

    value >>= 1;
  }

  return sum;
}

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetFinalSelection_effectTypes(
    DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected) {
  int i;
  int minNumEffects = 1000;
  int numEffects = 0;
  int effects = 0;
  DRCDEC_SELECTION_DATA* pCandidate = NULL;
  DRC_INSTRUCTIONS_UNI_DRC* pInst = NULL;

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    pInst = pCandidate->pInst;

    effects = pInst->drcSetEffect;
    effects &= 0xffff ^ (EB_GENERAL_COMPR);
    numEffects = _crossSum(effects);

    if (numEffects < minNumEffects) {
      minNumEffects = numEffects;
    }
  }

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    pInst = pCandidate->pInst;

    effects = pInst->drcSetEffect;
    effects &= 0xffff ^ (EB_GENERAL_COMPR);
    numEffects = _crossSum(effects);

    if (numEffects == minNumEffects) {
      if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
        return DRCDEC_SELECTION_PROCESS_NOT_OK;
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _selectSmallestTargetLoudnessValueUpper(
    DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected) {
  int i;
  SCHAR minVal = 0x7F;
  SCHAR val = 0;
  DRCDEC_SELECTION_DATA* pCandidate = NULL;

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    val = pCandidate->pInst->drcSetTargetLoudnessValueUpper;

    if (val < minVal) {
      minVal = val;
    }
  }

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    val = pCandidate->pInst->drcSetTargetLoudnessValueUpper;

    if (val == minVal) {
      if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
        return DRCDEC_SELECTION_PROCESS_NOT_OK;
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetFinalSelection_targetLoudness(
    FIXP_DBL targetLoudness, DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;
  int i;
  DRCDEC_SELECTION_DATA* pCandidate = NULL;

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    if (pCandidate->selectionFlag == 0) {
      if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
        return DRCDEC_SELECTION_PROCESS_NOT_OK;
    }
  }

  if (_drcdec_selection_getNumber(pCandidatesSelected) == 0) {
    retVal = _selectSmallestTargetLoudnessValueUpper(pCandidatesPotential,
                                                     pCandidatesSelected);
    if (retVal) return (retVal);
  }

  if (_drcdec_selection_getNumber(pCandidatesSelected) > 1) {
    DRC_INSTRUCTIONS_UNI_DRC* pDrcInstructionUniDrc = NULL;

    _swapSelectionAndClear(&pCandidatesPotential, &pCandidatesSelected);

    for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
      pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
      if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

      pDrcInstructionUniDrc = pCandidate->pInst;

      if (_targetLoudnessInRange(pDrcInstructionUniDrc, targetLoudness)) {
        if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
          return DRCDEC_SELECTION_PROCESS_NOT_OK;
      }
    }

    if (_drcdec_selection_getNumber(pCandidatesSelected) > 1) {
      _swapSelectionAndClear(&pCandidatesPotential, &pCandidatesSelected);

      retVal = _selectSmallestTargetLoudnessValueUpper(pCandidatesPotential,
                                                       pCandidatesSelected);
      if (retVal) return (retVal);
    }
  }

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetFinalSelection_peakValueLargest(
    DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected) {
  int i;
  FIXP_DBL largestPeakLevel = MINVAL_DBL;
  FIXP_DBL peakLevel = 0;
  DRCDEC_SELECTION_DATA* pCandidate = NULL;

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    peakLevel = pCandidate->outputPeakLevel;

    if (peakLevel > largestPeakLevel) {
      largestPeakLevel = peakLevel;
    }
  }

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    peakLevel = pCandidate->outputPeakLevel;

    if (peakLevel == largestPeakLevel) {
      if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
        return DRCDEC_SELECTION_PROCESS_NOT_OK;
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetFinalSelection_drcSetId(
    DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected) {
  int i;
  int largestId = -1000;
  int id = 0;
  DRCDEC_SELECTION_DATA* pCandidate = NULL;
  DRCDEC_SELECTION_DATA* pCandidateSelected = NULL;

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    pCandidate = _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    id = pCandidate->pInst->drcSetId;

    if (id > largestId) {
      largestId = id;
      pCandidateSelected = pCandidate;
    }
  }

  if (pCandidateSelected != NULL) {
    if (_drcdec_selection_add(pCandidatesSelected, pCandidateSelected) == NULL)
      return DRCDEC_SELECTION_PROCESS_NOT_OK;
  } else {
    return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetFinalSelection(
    HANDLE_SEL_PROC_INPUT hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    DRCDEC_SELECTION** ppCandidatesPotential,
    DRCDEC_SELECTION** ppCandidatesSelected, SEL_PROC_CODEC_MODE codecMode) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;

  if (_drcdec_selection_getNumber(*ppCandidatesPotential) == 0) {
    return DRCDEC_SELECTION_PROCESS_NOT_OK;
  } else if (_drcdec_selection_getNumber(*ppCandidatesPotential) == 1) {
    _swapSelection(ppCandidatesPotential, ppCandidatesSelected);
  } else
  {
    retVal = _drcSetFinalSelection_peakValue0(*ppCandidatesPotential,
                                              *ppCandidatesSelected);
    if (retVal) return (retVal);

    if (_drcdec_selection_getNumber(*ppCandidatesSelected) > 1) {
      _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
      retVal = _drcSetFinalSelection_downmixId(
          hSelProcInput, ppCandidatesPotential, ppCandidatesSelected);
      if (retVal) return (retVal);
    }

    if (_drcdec_selection_getNumber(*ppCandidatesSelected) > 1) {
      _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
      retVal = _drcSetFinalSelection_effectTypes(*ppCandidatesPotential,
                                                 *ppCandidatesSelected);
      if (retVal) return (retVal);
    }

    if (_drcdec_selection_getNumber(*ppCandidatesSelected) > 1) {
      _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
      retVal = _drcSetFinalSelection_targetLoudness(
          hSelProcInput->targetLoudness, *ppCandidatesPotential,
          *ppCandidatesSelected);
      if (retVal) return (retVal);
    }

    if (_drcdec_selection_getNumber(*ppCandidatesSelected) > 1) {
      _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
      retVal = _drcSetFinalSelection_peakValueLargest(*ppCandidatesPotential,
                                                      *ppCandidatesSelected);
      if (retVal) return (retVal);
    }

    if (_drcdec_selection_getNumber(*ppCandidatesSelected) > 1) {
      _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
      retVal = _drcSetFinalSelection_drcSetId(*ppCandidatesPotential,
                                              *ppCandidatesSelected);
      if (retVal) return (retVal);
    }
  }

  if (_drcdec_selection_getNumber(*ppCandidatesSelected) == 0) {
    return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _generateVirtualDrcSets(
    HANDLE_SEL_PROC_INPUT hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    SEL_PROC_CODEC_MODE codecMode) {
  int i;
  int nMixes = hUniDrcConfig->downmixInstructionsCount + 1;
  int index = hUniDrcConfig->drcInstructionsUniDrcCount;
  int indexVirtual = -1;
  DRC_INSTRUCTIONS_UNI_DRC* pDrcInstruction =
      &(hUniDrcConfig->drcInstructionsUniDrc[index]);

  if (codecMode == SEL_PROC_MPEG_H_3DA) {
    nMixes = 1;
  }

  if ((index + nMixes) > (12 + 1 + 6)) {
    return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  FDKmemset(pDrcInstruction, 0, sizeof(DRC_INSTRUCTIONS_UNI_DRC));

  pDrcInstruction->drcSetId = indexVirtual;
  index++;
  indexVirtual--;
  pDrcInstruction->downmixIdCount = 1;

  if ((codecMode == SEL_PROC_MPEG_H_3DA) &&
      (hSelProcInput->numDownmixIdRequests)) {
    pDrcInstruction->downmixId[0] = hSelProcInput->downmixIdRequested[0];
  } else {
    pDrcInstruction->downmixId[0] = DOWNMIX_ID_BASE_LAYOUT;
  }

  for (i = 1; i < nMixes; i++) {
    pDrcInstruction = &(hUniDrcConfig->drcInstructionsUniDrc[index]);
    FDKmemset(pDrcInstruction, 0, sizeof(DRC_INSTRUCTIONS_UNI_DRC));
    pDrcInstruction->drcSetId = indexVirtual;
    pDrcInstruction->downmixId[0] =
        hUniDrcConfig->downmixInstructions[i - 1].downmixId;
    pDrcInstruction->downmixIdCount = 1;
    index++;
    indexVirtual--;
  }

  hUniDrcConfig->drcInstructionsCountInclVirtual =
      hUniDrcConfig->drcInstructionsUniDrcCount + nMixes;

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _generateOutputInfo(
    HANDLE_SEL_PROC_INPUT hSelProcInput, HANDLE_SEL_PROC_OUTPUT hSelProcOutput,
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
    DRCDEC_SELECTION_DATA* pSelectionData, SEL_PROC_CODEC_MODE codecMode) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;

  int i, j;
  int hasDependend = 0;
  int hasFading = 0;
  int hasDucking = 0;
  int selectedDrcSetIds;
  int selectedDownmixIds;
  FIXP_DBL mixingLevel = 0;
  int albumMode = hSelProcInput->albumMode;
  UCHAR* pDownmixIdRequested = hSelProcInput->downmixIdRequested;
  FIXP_SGL boost = hSelProcInput->boost;
  FIXP_SGL compress = hSelProcInput->compress;

  hSelProcOutput->numSelectedDrcSets = 1;
  hSelProcOutput->selectedDrcSetIds[0] = pSelectionData->pInst->drcSetId;
  hSelProcOutput->selectedDownmixIds[0] =
      pSelectionData->pInst->drcApplyToDownmix == 1
          ? pSelectionData->pInst->downmixId[0]
          : 0;
  hSelProcOutput->loudnessNormalizationGainDb =
      pSelectionData->loudnessNormalizationGainDbAdjusted +
      hSelProcInput->loudnessNormalizationGainModificationDb;
  hSelProcOutput->outputPeakLevelDb = pSelectionData->outputPeakLevel;
  hSelProcOutput->outputLoudness = pSelectionData->outputLoudness;

  hSelProcOutput->boost = boost;
  hSelProcOutput->compress = compress;
  hSelProcOutput->baseChannelCount =
      hUniDrcConfig->channelLayout.baseChannelCount;
  hSelProcOutput->targetChannelCount =
      hUniDrcConfig->channelLayout.baseChannelCount;
  hSelProcOutput->activeDownmixId =
      pDownmixIdRequested[pSelectionData->downmixIdRequestIndex];

  _getMixingLevel(hLoudnessInfoSet, *pDownmixIdRequested,
                  hSelProcOutput->selectedDrcSetIds[0], albumMode,
                  &mixingLevel);
  hSelProcOutput->mixingLevel = mixingLevel;

  if (pSelectionData->pInst->dependsOnDrcSetPresent) {
    int dependsOnDrcSetID = pSelectionData->pInst->dependsOnDrcSet;

    for (i = 0; i < hUniDrcConfig->drcInstructionsCountInclVirtual; i++) {
      DRC_INSTRUCTIONS_UNI_DRC* pInst =
          &(hUniDrcConfig->drcInstructionsUniDrc[i]);
      if (!_drcSetIsUsable(hUniDrcConfig, pInst)) continue;

      if (pInst->drcSetId == dependsOnDrcSetID) {
        hSelProcOutput->selectedDrcSetIds[hSelProcOutput->numSelectedDrcSets] =
            hUniDrcConfig->drcInstructionsUniDrc[i].drcSetId;
        hSelProcOutput->selectedDownmixIds[hSelProcOutput->numSelectedDrcSets] =
            hUniDrcConfig->drcInstructionsUniDrc[i].drcApplyToDownmix == 1
                ? hUniDrcConfig->drcInstructionsUniDrc[i].downmixId[0]
                : 0;
        hSelProcOutput->numSelectedDrcSets++;
        hasDependend = 1;
        break;
      }
    }
  }

  if (hSelProcInput->albumMode == 0) {
    for (i = 0; i < hUniDrcConfig->drcInstructionsUniDrcCount; i++) {
      DRC_INSTRUCTIONS_UNI_DRC* pInst =
          &(hUniDrcConfig->drcInstructionsUniDrc[i]);
      if (!_drcSetIsUsable(hUniDrcConfig, pInst)) continue;

      if (pInst->drcSetEffect & EB_FADE) {
        if (pInst->downmixId[0] == DOWNMIX_ID_ANY_DOWNMIX) {
          hSelProcOutput->numSelectedDrcSets = hasDependend + 1;
          hSelProcOutput
              ->selectedDrcSetIds[hSelProcOutput->numSelectedDrcSets] =
              hUniDrcConfig->drcInstructionsUniDrc[i].drcSetId;
          hSelProcOutput
              ->selectedDownmixIds[hSelProcOutput->numSelectedDrcSets] =
              hUniDrcConfig->drcInstructionsUniDrc[i].drcApplyToDownmix == 1
                  ? hUniDrcConfig->drcInstructionsUniDrc[i].downmixId[0]
                  : 0;
          hSelProcOutput->numSelectedDrcSets++;
          hasFading = 1;

        } else {
          retVal = DRCDEC_SELECTION_PROCESS_NOT_OK;
          if (retVal) return DRCDEC_SELECTION_PROCESS_NOT_OK;
        }
      }
    }
  }

  for (i = 0; i < hUniDrcConfig->drcInstructionsUniDrcCount; i++) {
    DRC_INSTRUCTIONS_UNI_DRC* pInst =
        &(hUniDrcConfig->drcInstructionsUniDrc[i]);
    if (!_drcSetIsUsable(hUniDrcConfig, pInst)) continue;

    if (pInst->drcSetEffect & (EB_DUCK_OTHER | EB_DUCK_SELF)) {
      for (j = 0; j < pInst->downmixIdCount; j++) {
        if (pInst->downmixId[j] == hSelProcOutput->activeDownmixId) {
          hSelProcOutput->numSelectedDrcSets =
              hasDependend + 1;

          hSelProcOutput
              ->selectedDrcSetIds[hSelProcOutput->numSelectedDrcSets] =
              hUniDrcConfig->drcInstructionsUniDrc[i].drcSetId;
          hSelProcOutput
              ->selectedDownmixIds[hSelProcOutput->numSelectedDrcSets] = 0;
          hSelProcOutput->numSelectedDrcSets++;
          hasDucking = 1;
        }
      }
    }
  }

  if (!hasDucking) {
    for (i = 0; i < hUniDrcConfig->drcInstructionsUniDrcCount; i++) {
      DRC_INSTRUCTIONS_UNI_DRC* pInst =
          &(hUniDrcConfig->drcInstructionsUniDrc[i]);
      if (!_drcSetIsUsable(hUniDrcConfig, pInst)) continue;

      if (pInst->drcSetEffect & (EB_DUCK_OTHER | EB_DUCK_SELF)) {
        for (j = 0; j < pInst->downmixIdCount; j++) {
          if (pInst->downmixId[j] == DOWNMIX_ID_BASE_LAYOUT) {
            hSelProcOutput->numSelectedDrcSets = hasDependend + hasFading + 1;
            hSelProcOutput
                ->selectedDrcSetIds[hSelProcOutput->numSelectedDrcSets] =
                hUniDrcConfig->drcInstructionsUniDrc[i].drcSetId;
            hSelProcOutput
                ->selectedDownmixIds[hSelProcOutput->numSelectedDrcSets] = 0;
            hSelProcOutput->numSelectedDrcSets++;
          }
        }
      }
    }
  }

  if (hSelProcOutput->numSelectedDrcSets > 3) {
    hSelProcOutput->numSelectedDrcSets = 0;
    return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  if (hSelProcOutput->numSelectedDrcSets == 3) {
    selectedDrcSetIds = hSelProcOutput->selectedDrcSetIds[0];
    selectedDownmixIds = hSelProcOutput->selectedDownmixIds[0];
    hSelProcOutput->selectedDrcSetIds[0] = hSelProcOutput->selectedDrcSetIds[2];
    hSelProcOutput->selectedDownmixIds[0] =
        hSelProcOutput->selectedDownmixIds[2];
    hSelProcOutput->selectedDrcSetIds[2] = selectedDrcSetIds;
    hSelProcOutput->selectedDownmixIds[2] = selectedDownmixIds;
  } else if (hSelProcOutput->numSelectedDrcSets == 2) {
    selectedDrcSetIds = hSelProcOutput->selectedDrcSetIds[0];
    selectedDownmixIds = hSelProcOutput->selectedDownmixIds[0];
    hSelProcOutput->selectedDrcSetIds[0] = hSelProcOutput->selectedDrcSetIds[1];
    hSelProcOutput->selectedDownmixIds[0] =
        hSelProcOutput->selectedDownmixIds[1];
    hSelProcOutput->selectedDrcSetIds[1] = selectedDrcSetIds;
    hSelProcOutput->selectedDownmixIds[1] = selectedDownmixIds;
  }

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _selectDownmixMatrix(
    HANDLE_SEL_PROC_OUTPUT hSelProcOutput,
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig) {
  int i;
  hSelProcOutput->baseChannelCount =
      hUniDrcConfig->channelLayout.baseChannelCount;
  hSelProcOutput->targetChannelCount =
      hUniDrcConfig->channelLayout.baseChannelCount;
  hSelProcOutput->targetLayout = -1;
  hSelProcOutput->downmixMatrixPresent = 0;

  if (hSelProcOutput->activeDownmixId != 0) {
    for (i = 0; i < hUniDrcConfig->downmixInstructionsCount; i++) {
      DOWNMIX_INSTRUCTIONS* pDown = &(hUniDrcConfig->downmixInstructions[i]);
      if (pDown->targetChannelCount > 8) {
        continue;
      }

      if (hSelProcOutput->activeDownmixId == pDown->downmixId) {
        hSelProcOutput->targetChannelCount = pDown->targetChannelCount;
        hSelProcOutput->targetLayout = pDown->targetLayout;

        if (pDown->downmixCoefficientsPresent) {
          int j, k;
          FIXP_DBL downmixOffset = getDownmixOffset(
              pDown, hSelProcOutput->baseChannelCount);

          for (j = 0; j < hSelProcOutput->baseChannelCount; j++) {
            for (k = 0; k < hSelProcOutput->targetChannelCount; k++) {
              hSelProcOutput->downmixMatrix[j][k] =
                  fMultDiv2(
                      downmixOffset,
                      pDown->downmixCoefficient[j + k * hSelProcOutput
                                                            ->baseChannelCount])
                  << 2;
            }
          }

          hSelProcOutput->downmixMatrixPresent = 1;
        }
        break;
      }
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetPreSelection(
    SEL_PROC_INPUT* hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
    DRCDEC_SELECTION** ppCandidatesPotential,
    DRCDEC_SELECTION** ppCandidatesSelected, SEL_PROC_CODEC_MODE codecMode) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;
  int i, j;

  for (i = 0; i < hSelProcInput->numDownmixIdRequests; i++) {
    for (j = 0; j < hUniDrcConfig->drcInstructionsCountInclVirtual; j++) {
      DRC_INSTRUCTIONS_UNI_DRC* pDrcInstruction =
          &(hUniDrcConfig->drcInstructionsUniDrc[j]);
      if (selectDrcInstructions(hUniDrcConfig, pDrcInstruction->drcSetId) !=
          pDrcInstruction)
        continue;

      retVal = _drcSetPreSelectionSingleInstruction(
          hSelProcInput, i, hUniDrcConfig, hLoudnessInfoSet, pDrcInstruction,
          *ppCandidatesPotential, *ppCandidatesSelected, codecMode);
      if (retVal) return DRCDEC_SELECTION_PROCESS_NOT_OK;
    }
  }

  retVal = _preSelectionRequirement9(hSelProcInput, *ppCandidatesPotential,
                                     *ppCandidatesSelected);
  if (retVal) return DRCDEC_SELECTION_PROCESS_NOT_OK;

  if (_drcdec_selection_getNumber(*ppCandidatesSelected) == 0) {
    retVal = _drcSetSelectionAddCandidates(
        hSelProcInput, *ppCandidatesPotential, *ppCandidatesSelected);
    if (retVal) return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _drcSetRequestSelection(
    SEL_PROC_INPUT* hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
    DRCDEC_SELECTION** ppCandidatesPotential,
    DRCDEC_SELECTION** ppCandidatesSelected) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal;
  int i;

  if (_drcdec_selection_getNumber(*ppCandidatesPotential) == 0) {
    retVal = DRCDEC_SELECTION_PROCESS_NOT_OK;
    if (retVal) return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  if (hSelProcInput->dynamicRangeControlOn) {
    if (hSelProcInput->numDrcFeatureRequests == 0) {
      retVal = _selectDrcSetEffectNone(hUniDrcConfig, *ppCandidatesPotential,
                                       *ppCandidatesSelected);
      if (retVal) return (retVal);

      if (_drcdec_selection_getNumber(*ppCandidatesSelected) == 0) {
        DRC_FEATURE_REQUEST fallbackRequest;
        fallbackRequest.drcEffectType.numRequests = 5;
        fallbackRequest.drcEffectType.numRequestsDesired = 5;
        fallbackRequest.drcEffectType.request[0] = DETR_GENERAL_COMPR;
        fallbackRequest.drcEffectType.request[1] = DETR_NIGHT;
        fallbackRequest.drcEffectType.request[2] = DETR_NOISY;
        fallbackRequest.drcEffectType.request[3] = DETR_LIMITED;
        fallbackRequest.drcEffectType.request[4] = DETR_LOWLEVEL;

        retVal = _selectEffectTypeFeature(hUniDrcConfig, fallbackRequest,
                                          ppCandidatesPotential,
                                          ppCandidatesSelected);
        if (retVal) return DRCDEC_SELECTION_PROCESS_NOT_OK;
      }

      _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
    } else {
      for (i = 0; i < hSelProcInput->numDrcFeatureRequests; i++) {
        if (hSelProcInput->drcFeatureRequestType[i] == DFRT_EFFECT_TYPE) {
          retVal = _selectEffectTypeFeature(
              hUniDrcConfig, hSelProcInput->drcFeatureRequest[i],
              ppCandidatesPotential, ppCandidatesSelected);

          _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
          if (retVal) return DRCDEC_SELECTION_PROCESS_NOT_OK;
        }

        else if (hSelProcInput->drcFeatureRequestType[i] ==
                 DFRT_DYNAMIC_RANGE) {
          retVal = _selectDynamicRange(
              hUniDrcConfig, hLoudnessInfoSet,
              hSelProcInput->drcFeatureRequest[i],
              hSelProcInput->downmixIdRequested, hSelProcInput->albumMode,
              *ppCandidatesPotential, *ppCandidatesSelected);

          if (_drcdec_selection_getNumber(*ppCandidatesSelected) > 0) {
            _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
          }
          if (retVal) return DRCDEC_SELECTION_PROCESS_NOT_OK;
        } else if (hSelProcInput->drcFeatureRequestType[i] ==
                   DFRT_DRC_CHARACTERISTIC) {
          retVal = _selectDrcCharacteristic(
              hUniDrcConfig,
              hSelProcInput->drcFeatureRequest[i].drcCharacteristic,
              ppCandidatesPotential, ppCandidatesSelected);

          if (_drcdec_selection_getNumber(*ppCandidatesSelected) > 0) {
            _swapSelectionAndClear(ppCandidatesPotential, ppCandidatesSelected);
          }
          if (retVal) return DRCDEC_SELECTION_PROCESS_NOT_OK;
        }
      }
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _dynamicRangeMeasurement(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, DRC_INSTRUCTIONS_UNI_DRC* pInst,
    UCHAR downmixIdRequested,
    DYN_RANGE_MEASUREMENT_REQUEST_TYPE dynamicRangeMeasurementType,
    int albumMode, int* pPeakToAveragePresent, FIXP_DBL* pPeakToAverage) {
  int i;
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;
  int drcSetId = fMax(0, pInst->drcSetId);

  *pPeakToAveragePresent = 0;

  if (albumMode) {
    for (i = 0; i < hLoudnessInfoSet->loudnessInfoAlbumCount; i++) {
      LOUDNESS_INFO* pLoudnessInfo = &(hLoudnessInfoSet->loudnessInfoAlbum[i]);

      if (drcSetId == pLoudnessInfo->drcSetId) {
        if (downmixIdRequested == pLoudnessInfo->downmixId) {
          retVal = _extractLoudnessPeakToAverageValue(
              pLoudnessInfo, dynamicRangeMeasurementType, pPeakToAveragePresent,
              pPeakToAverage);
          if (retVal) return (retVal);
        }
      }
    }
  }

  if (*pPeakToAveragePresent == 0) {
    for (i = 0; i < hLoudnessInfoSet->loudnessInfoCount; i++) {
      LOUDNESS_INFO* pLoudnessInfo = &(hLoudnessInfoSet->loudnessInfo[i]);

      if (drcSetId == pLoudnessInfo->drcSetId) {
        if (downmixIdRequested == pLoudnessInfo->downmixId) {
          retVal = _extractLoudnessPeakToAverageValue(
              pLoudnessInfo, dynamicRangeMeasurementType, pPeakToAveragePresent,
              pPeakToAverage);
          if (retVal) return (retVal);
        }
      }
    }
  }

  return retVal;
}

static DRCDEC_SELECTION_DATA* _drcdec_selection_addNew(
    DRCDEC_SELECTION* pSelection) {
  if (pSelection->numData < (12 + 1 + 6)) {
    DRCDEC_SELECTION_DATA* pData = &(pSelection->data[pSelection->numData]);
    FDKmemset(pData, 0, sizeof(DRCDEC_SELECTION_DATA));
    pSelection->numData++;

    return pData;
  } else {
    return NULL;
  }
}

static DRCDEC_SELECTION_DATA* _drcdec_selection_add(
    DRCDEC_SELECTION* pSelection, DRCDEC_SELECTION_DATA* pDataIn) {
  if (pSelection->numData < (12 + 1 + 6)) {
    DRCDEC_SELECTION_DATA* pData = &(pSelection->data[pSelection->numData]);
    FDKmemcpy(pData, pDataIn, sizeof(DRCDEC_SELECTION_DATA));
    pSelection->numData++;
    return pData;
  } else {
    return NULL;
  }
}

static int _drcdec_selection_clear(DRCDEC_SELECTION* pSelection) {
  return pSelection->numData = 0;
}

static int _drcdec_selection_getNumber(DRCDEC_SELECTION* pSelection) {
  return pSelection->numData;
}

static int _drcdec_selection_setNumber(DRCDEC_SELECTION* pSelection, int num) {
  if (num >= 0 && num < pSelection->numData) {
    return pSelection->numData = num;
  } else {
    return pSelection->numData;
  }
}

static DRCDEC_SELECTION_DATA* _drcdec_selection_getAt(
    DRCDEC_SELECTION* pSelection, int at) {
  if (at >= 0 && at < (12 + 1 + 6)) {
    return &(pSelection->data[at]);
  } else {
    return NULL;
  }
}

static int _swapSelectionAndClear(DRCDEC_SELECTION** ppCandidatesPotential,
                                  DRCDEC_SELECTION** ppCandidatesSelected) {
  DRCDEC_SELECTION* pTmp = *ppCandidatesPotential;
  *ppCandidatesPotential = *ppCandidatesSelected;
  *ppCandidatesSelected = pTmp;
  _drcdec_selection_clear(*ppCandidatesSelected);
  return 0;
}

static int _swapSelection(DRCDEC_SELECTION** ppCandidatesPotential,
                          DRCDEC_SELECTION** ppCandidatesSelected) {
  DRCDEC_SELECTION* pTmp = *ppCandidatesPotential;
  *ppCandidatesPotential = *ppCandidatesSelected;
  *ppCandidatesSelected = pTmp;
  return 0;
}

static LOUDNESS_INFO* _getLoudnessInfoStructure(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, int drcSetId, int downmixId,
    int albumMode) {
  int i, j;
  int count;

  LOUDNESS_INFO* pLoudnessInfo = NULL;

  if (albumMode) {
    count = hLoudnessInfoSet->loudnessInfoAlbumCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfoAlbum;
  } else {
    count = hLoudnessInfoSet->loudnessInfoCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfo;
  }

  for (i = 0; i < count; i++) {
    if ((pLoudnessInfo[i].drcSetId == drcSetId) &&
        (pLoudnessInfo[i].downmixId == downmixId)) {
      for (j = 0; j < pLoudnessInfo[i].measurementCount; j++) {
        if ((pLoudnessInfo[i].loudnessMeasurement[j].methodDefinition == 1) ||
            (pLoudnessInfo[i].loudnessMeasurement[j].methodDefinition == 2)) {
          return &pLoudnessInfo[i];
        }
      }
    }
  }

  return NULL;
}

static LOUDNESS_INFO* _getApplicableLoudnessInfoStructure(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, int drcSetId,
    int downmixIdRequested, int albumMode) {
  LOUDNESS_INFO* pLoudnessInfo = NULL;

  pLoudnessInfo = _getLoudnessInfoStructure(hLoudnessInfoSet, drcSetId,
                                            downmixIdRequested, albumMode);

  if (pLoudnessInfo == NULL) {
    pLoudnessInfo =
        _getLoudnessInfoStructure(hLoudnessInfoSet, drcSetId, 0x7F, albumMode);
  }

  if (pLoudnessInfo == NULL) {
    pLoudnessInfo = _getLoudnessInfoStructure(hLoudnessInfoSet, 0x3F,
                                              downmixIdRequested, albumMode);
  }

  if (pLoudnessInfo == NULL) {
    pLoudnessInfo = _getLoudnessInfoStructure(hLoudnessInfoSet, 0,
                                              downmixIdRequested, albumMode);
  }

  if (pLoudnessInfo == NULL) {
    pLoudnessInfo =
        _getLoudnessInfoStructure(hLoudnessInfoSet, 0x3F, 0x7F, albumMode);
  }

  if (pLoudnessInfo == NULL) {
    pLoudnessInfo =
        _getLoudnessInfoStructure(hLoudnessInfoSet, 0, 0x7F, albumMode);
  }

  if (pLoudnessInfo == NULL) {
    pLoudnessInfo =
        _getLoudnessInfoStructure(hLoudnessInfoSet, drcSetId, 0, albumMode);
  }

  if (pLoudnessInfo == NULL) {
    pLoudnessInfo =
        _getLoudnessInfoStructure(hLoudnessInfoSet, 0x3F, 0, albumMode);
  }

  if (pLoudnessInfo == NULL) {
    pLoudnessInfo =
        _getLoudnessInfoStructure(hLoudnessInfoSet, 0, 0, albumMode);
  }

  return pLoudnessInfo;
}

typedef struct {
  FIXP_DBL value;
  int order;
} VALUE_ORDER;

void _initValueOrder(VALUE_ORDER* pValue) {
  pValue->value = (FIXP_DBL)0;
  pValue->order = -1;
}

enum {
  MS_BONUS0 = 0,
  MS_BONUS1770,
  MS_BONUSUSER,
  MS_BONUSEXPERT,
  MS_RESA,
  MS_RESB,
  MS_RESC,
  MS_RESD,
  MS_RESE,
  MS_PROGRAMLOUDNESS,
  MS_PEAKLOUDNESS
};

static DRCDEC_SELECTION_PROCESS_RETURN _getMethodValue(
    VALUE_ORDER* pValueOrder, FIXP_DBL value, int measurementSystem,
    int measurementSystemRequested) {
  const int rows = 11;
  const int columns = 12;
  const int pOrdering[rows][columns] = {
      {0, 0, 8, 0, 1, 3, 0, 5, 6, 7, 4, 2},
      {0, 0, 8, 0, 1, 3, 0, 5, 6, 7, 4, 2},
      {0, 0, 1, 0, 8, 5, 0, 2, 3, 4, 6, 7},
      {0, 0, 3, 0, 1, 8, 0, 4, 5, 6, 7, 2},
      {0, 0, 5, 0, 1, 3, 0, 8, 6, 7, 4, 2},
      {0, 0, 5, 0, 1, 3, 0, 6, 8, 7, 4, 2},
      {0, 0, 5, 0, 1, 3, 0, 6, 7, 8, 4, 2},
      {0, 0, 3, 0, 1, 7, 0, 4, 5, 6, 8, 2},
      {0, 0, 1, 0, 7, 5, 0, 2, 3, 4, 6, 8},
      {0, 0, 1, 0, 0, 0, 0, 2, 3, 4, 0, 0},
      {0, 7, 0, 0, 0, 0, 6, 5, 4, 3, 2, 1}
  };

  if (measurementSystemRequested < 0 || measurementSystemRequested >= rows ||
      measurementSystem < 0 || measurementSystem >= columns) {
    return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  if (pOrdering[measurementSystemRequested][measurementSystem] >
      pValueOrder->order) {
    pValueOrder->order =
        pOrdering[measurementSystemRequested][measurementSystem];
    pValueOrder->value = value;
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _getLoudness(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, int albumMode,
    METHOD_DEFINITION_REQUEST measurementMethodRequested,
    MEASUREMENT_SYSTEM_REQUEST measurementSystemRequested,
    FIXP_DBL targetLoudness,
    int drcSetId, int downmixIdRequested,
    FIXP_DBL* pLoudnessNormalizationGain,
    FIXP_DBL* pLoudness)
{
  int index;

  LOUDNESS_INFO* pLoudnessInfo = NULL;
  VALUE_ORDER valueOrder;

  METHOD_DEFINITION_REQUEST requestedMethodDefinition =
      measurementMethodRequested < MDR_ANCHOR_LOUDNESS ? MDR_PROGRAM_LOUDNESS
                                                       : MDR_ANCHOR_LOUDNESS;

  if (measurementMethodRequested > MDR_ANCHOR_LOUDNESS) {
    return DRCDEC_SELECTION_PROCESS_NOT_OK;
  }

  _initValueOrder(&valueOrder);

  *pLoudness = UNDEFINED_LOUDNESS_VALUE;
  *pLoudnessNormalizationGain = (FIXP_DBL)0;

  if (drcSetId < 0) {
    drcSetId = 0;
  }

  pLoudnessInfo = _getApplicableLoudnessInfoStructure(
      hLoudnessInfoSet, drcSetId, downmixIdRequested, albumMode);

  if (albumMode && (pLoudnessInfo == NULL)) {
    pLoudnessInfo = _getApplicableLoudnessInfoStructure(
        hLoudnessInfoSet, drcSetId, downmixIdRequested, 0);
  }

  if (pLoudnessInfo == NULL) {
    return DRCDEC_SELECTION_PROCESS_NO_ERROR;
  }

  index = -1;

  do {
    index = _findMethodDefinition(pLoudnessInfo, requestedMethodDefinition,
                                  index + 1);

    if (index >= 0) {
      _getMethodValue(
          &valueOrder, pLoudnessInfo->loudnessMeasurement[index].methodValue,
          pLoudnessInfo->loudnessMeasurement[index].measurementSystem,
          measurementSystemRequested);
    }
  } while (index >= 0);

  if (valueOrder.order == -1) {
    index = -1;

    do {
      index = _findMethodDefinition(
          pLoudnessInfo,
          requestedMethodDefinition == MDR_PROGRAM_LOUDNESS
              ? MDR_ANCHOR_LOUDNESS
              : MDR_PROGRAM_LOUDNESS,
          index + 1);

      if (index >= 0) {
        _getMethodValue(
            &valueOrder, pLoudnessInfo->loudnessMeasurement[index].methodValue,
            pLoudnessInfo->loudnessMeasurement[index].measurementSystem,
            measurementSystemRequested);
      }
    } while (index >= 0);
  }

  if (valueOrder.order == -1) {
    return DRCDEC_SELECTION_PROCESS_NOT_OK;
  } else {
    *pLoudnessNormalizationGain = targetLoudness - valueOrder.value;
    *pLoudness = valueOrder.value;
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static int _truePeakLevelIsPresent(HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
                                   int drcSetId, int downmixId, int albumMode) {
  int i;
  int count;
  LOUDNESS_INFO* pLoudnessInfo = NULL;

  if (albumMode) {
    count = hLoudnessInfoSet->loudnessInfoAlbumCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfoAlbum;
  } else {
    count = hLoudnessInfoSet->loudnessInfoCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfo;
  }

  for (i = 0; i < count; i++) {
    if ((pLoudnessInfo[i].drcSetId == drcSetId) &&
        (pLoudnessInfo[i].downmixId == downmixId)) {
      if (pLoudnessInfo[i].truePeakLevelPresent) return 1;
    }
  }

  return 0;
}

static DRCDEC_SELECTION_PROCESS_RETURN _getTruePeakLevel(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, int drcSetId, int downmixId,
    int albumMode, FIXP_DBL* pTruePeakLevel) {
  int i;
  int count;
  LOUDNESS_INFO* pLoudnessInfo = NULL;

  if (albumMode) {
    count = hLoudnessInfoSet->loudnessInfoAlbumCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfoAlbum;
  } else {
    count = hLoudnessInfoSet->loudnessInfoCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfo;
  }

  for (i = 0; i < count; i++) {
    if ((pLoudnessInfo[i].drcSetId == drcSetId) &&
        (pLoudnessInfo[i].downmixId == downmixId)) {
      if (pLoudnessInfo[i].truePeakLevelPresent) {
        *pTruePeakLevel = pLoudnessInfo[i].truePeakLevel;
        return DRCDEC_SELECTION_PROCESS_NO_ERROR;
      }
    }
  }

  return DRCDEC_SELECTION_PROCESS_NOT_OK;
}

static int _samplePeakLevelIsPresent(HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
                                     int drcSetId, int downmixId,
                                     int albumMode) {
  int i;
  int count;
  LOUDNESS_INFO* pLoudnessInfo = NULL;

  if (albumMode) {
    count = hLoudnessInfoSet->loudnessInfoAlbumCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfoAlbum;
  } else {
    count = hLoudnessInfoSet->loudnessInfoCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfo;
  }

  for (i = 0; i < count; i++) {
    if ((pLoudnessInfo[i].drcSetId == drcSetId) &&
        (pLoudnessInfo[i].downmixId == downmixId)) {
      if (pLoudnessInfo[i].samplePeakLevelPresent) return 1;
    }
  }

  return 0;
}

static DRCDEC_SELECTION_PROCESS_RETURN _getSamplePeakLevel(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, int drcSetId, int downmixId,
    int albumMode, FIXP_DBL* pSamplePeakLevel
) {
  int i;
  int count;
  LOUDNESS_INFO* pLoudnessInfo = NULL;

  if (albumMode) {
    count = hLoudnessInfoSet->loudnessInfoAlbumCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfoAlbum;
  } else {
    count = hLoudnessInfoSet->loudnessInfoCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfo;
  }

  for (i = 0; i < count; i++) {
    if ((pLoudnessInfo[i].drcSetId == drcSetId) &&
        (pLoudnessInfo[i].downmixId == downmixId)) {
      if (pLoudnessInfo[i].samplePeakLevelPresent) {
        *pSamplePeakLevel = pLoudnessInfo[i].samplePeakLevel;
        return DRCDEC_SELECTION_PROCESS_NO_ERROR;
      }
    }
  }

  return DRCDEC_SELECTION_PROCESS_NOT_OK;
}

static int _limiterPeakTargetIsPresent(
    DRC_INSTRUCTIONS_UNI_DRC* pDrcInstruction, int drcSetId, int downmixId) {
  int i;

  if (pDrcInstruction->limiterPeakTargetPresent) {
    if ((pDrcInstruction->downmixId[0] == downmixId) ||
        (pDrcInstruction->downmixId[0] == 0x7F)) {
      return 1;
    }

    for (i = 0; i < pDrcInstruction->downmixIdCount; i++) {
      if (pDrcInstruction->downmixId[i] == downmixId) {
        return 1;
      }
    }
  }

  return 0;
}

static DRCDEC_SELECTION_PROCESS_RETURN _getLimiterPeakTarget(
    DRC_INSTRUCTIONS_UNI_DRC* pDrcInstruction, int drcSetId, int downmixId,
    FIXP_DBL* pLimiterPeakTarget) {
  int i;

  if (pDrcInstruction->limiterPeakTargetPresent) {
    if ((pDrcInstruction->downmixId[0] == downmixId) ||
        (pDrcInstruction->downmixId[0] == 0x7F)) {
      *pLimiterPeakTarget =
          ((FX_SGL2FX_DBL(pDrcInstruction->limiterPeakTarget) >> 2));
      return DRCDEC_SELECTION_PROCESS_NO_ERROR;
    }

    for (i = 0; i < pDrcInstruction->downmixIdCount; i++) {
      if (pDrcInstruction->downmixId[i] == downmixId) {
        *pLimiterPeakTarget =
            ((FX_SGL2FX_DBL(pDrcInstruction->limiterPeakTarget) >> 2));
        return DRCDEC_SELECTION_PROCESS_NO_ERROR;
      }
    }
  }

  return DRCDEC_SELECTION_PROCESS_NOT_OK;
}

static int _downmixCoefficientsArePresent(HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
                                          int downmixId, int* pIndex) {
  int i;
  *pIndex = -1;

  for (i = 0; i < hUniDrcConfig->downmixInstructionsCount; i++) {
    if (hUniDrcConfig->downmixInstructions[i].downmixId == downmixId) {
      if (hUniDrcConfig->downmixInstructions[i].downmixCoefficientsPresent) {
        if (hUniDrcConfig->downmixInstructions[i].targetChannelCount > 8)
          return 0;
        *pIndex = i;
        return 1;
      }
    }
  }

  return 0;
}

static DRCDEC_SELECTION_PROCESS_RETURN _getSignalPeakLevel(
    HANDLE_SEL_PROC_INPUT hSelProcInput, HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, DRC_INSTRUCTIONS_UNI_DRC* pInst,
    int downmixIdRequested, int* explicitPeakInformationPresent,
    FIXP_DBL* signalPeakLevelOut,
    SEL_PROC_CODEC_MODE codecMode

) {
  DRCDEC_SELECTION_PROCESS_RETURN retVal = DRCDEC_SELECTION_PROCESS_NO_ERROR;

  int albumMode = hSelProcInput->albumMode;

  FIXP_DBL signalPeakLevelTmp = (FIXP_DBL)0;
  FIXP_DBL signalPeakLevel = FIXP_DBL(0);

  int dmxId = downmixIdRequested;

  int drcSetId = pInst->drcSetId;

  if (drcSetId < 0) {
    drcSetId = 0;
  }

  *explicitPeakInformationPresent = 1;

  if (_truePeakLevelIsPresent(hLoudnessInfoSet, drcSetId, dmxId, albumMode)) {
    retVal = _getTruePeakLevel(hLoudnessInfoSet, drcSetId, dmxId, albumMode,
                               &signalPeakLevel);
    if (retVal) return (retVal);
  } else if (_samplePeakLevelIsPresent(hLoudnessInfoSet, drcSetId, dmxId,
                                       albumMode)) {
    retVal = _getSamplePeakLevel(hLoudnessInfoSet, drcSetId, dmxId, albumMode,
                                 &signalPeakLevel);
    if (retVal) return (retVal);
  } else if (_truePeakLevelIsPresent(hLoudnessInfoSet, 0x3F, dmxId,
                                     albumMode)) {
    retVal = _getTruePeakLevel(hLoudnessInfoSet, 0x3F, dmxId, albumMode,
                               &signalPeakLevel);
    if (retVal) return (retVal);
  } else if (_samplePeakLevelIsPresent(hLoudnessInfoSet, 0x3F, dmxId,
                                       albumMode)) {
    retVal = _getSamplePeakLevel(hLoudnessInfoSet, 0x3F, dmxId, albumMode,
                                 &signalPeakLevel);
    if (retVal) return (retVal);
  } else if (_limiterPeakTargetIsPresent(pInst, drcSetId, dmxId)) {
    retVal = _getLimiterPeakTarget(pInst, drcSetId, dmxId, &signalPeakLevel);
    if (retVal) return (retVal);
  } else if (dmxId != 0) {
    int downmixInstructionIndex = 0;
    FIXP_DBL downmixPeakLevelDB = 0;

    *explicitPeakInformationPresent = 0;

    signalPeakLevelTmp = FIXP_DBL(0);

    if (_downmixCoefficientsArePresent(hUniDrcConfig, dmxId,
                                       &downmixInstructionIndex)) {
      FIXP_DBL dB_m;
      int dB_e;
      FIXP_DBL coeff;
      FIXP_DBL sum, maxSum;
      int i, j;
      DOWNMIX_INSTRUCTIONS* pDown =
          &(hUniDrcConfig->downmixInstructions[downmixInstructionIndex]);
      FIXP_DBL downmixOffset = getDownmixOffset(
          pDown, hUniDrcConfig->channelLayout.baseChannelCount);
      maxSum = (FIXP_DBL)0;

      for (i = 0; i < pDown->targetChannelCount; i++) {
        sum = (FIXP_DBL)0;
        for (j = 0; j < hUniDrcConfig->channelLayout.baseChannelCount; j++) {
          coeff = pDown->downmixCoefficient[j + i * hUniDrcConfig->channelLayout
                                                        .baseChannelCount];
          sum += coeff >> 5;
        }
        if (maxSum < sum) maxSum = sum;
      }

      maxSum = fMultDiv2(maxSum, downmixOffset) << 2;

      if (maxSum == FL2FXCONST_DBL(1.0f / (float)(1 << 7))) {
        downmixPeakLevelDB = (FIXP_DBL)0;
      } else {
        dB_m = lin2dB(maxSum, 7, &dB_e);
        downmixPeakLevelDB =
            scaleValue(dB_m, dB_e - 7);
      }
    }

    if (_truePeakLevelIsPresent(hLoudnessInfoSet, drcSetId, 0, albumMode)) {
      retVal = _getTruePeakLevel(hLoudnessInfoSet, drcSetId, 0, albumMode,
                                 &signalPeakLevelTmp);
      if (retVal) return (retVal);
    } else if (_samplePeakLevelIsPresent(hLoudnessInfoSet, drcSetId, 0,
                                         albumMode)) {
      retVal = _getSamplePeakLevel(hLoudnessInfoSet, drcSetId, 0, albumMode,
                                   &signalPeakLevelTmp);
      if (retVal) return (retVal);
    } else if (_truePeakLevelIsPresent(hLoudnessInfoSet, 0x3F, 0, albumMode)) {
      retVal = _getTruePeakLevel(hLoudnessInfoSet, 0x3F, 0, albumMode,
                                 &signalPeakLevelTmp);
      if (retVal) return (retVal);
    } else if (_samplePeakLevelIsPresent(hLoudnessInfoSet, 0x3F, 0,
                                         albumMode)) {
      retVal = _getSamplePeakLevel(hLoudnessInfoSet, 0x3F, 0, albumMode,
                                   &signalPeakLevelTmp);
      if (retVal) return (retVal);
    } else if (_limiterPeakTargetIsPresent(pInst, drcSetId, 0)) {
      retVal = _getLimiterPeakTarget(pInst, drcSetId, 0, &signalPeakLevelTmp);
      if (retVal) return (retVal);
    }

    signalPeakLevel = signalPeakLevelTmp + downmixPeakLevelDB;
  } else {
    signalPeakLevel = FIXP_DBL(0);
    *explicitPeakInformationPresent = FIXP_DBL(0);
  }

  *signalPeakLevelOut = signalPeakLevel;

  return retVal;
}

static DRCDEC_SELECTION_PROCESS_RETURN _extractLoudnessPeakToAverageValue(
    LOUDNESS_INFO* loudnessInfo,
    DYN_RANGE_MEASUREMENT_REQUEST_TYPE dynamicRangeMeasurementType,
    int* pLoudnessPeakToAverageValuePresent,
    FIXP_DBL* pLoudnessPeakToAverageValue) {
  int i;

  VALUE_ORDER valueOrderLoudness;
  VALUE_ORDER valueOrderPeakLoudness;

  _initValueOrder(&valueOrderLoudness);
  _initValueOrder(&valueOrderPeakLoudness);

  LOUDNESS_MEASUREMENT* pLoudnessMeasure = NULL;

  *pLoudnessPeakToAverageValuePresent = 0;

  for (i = 0; i < loudnessInfo->measurementCount; i++) {
    pLoudnessMeasure = &(loudnessInfo->loudnessMeasurement[i]);

    if (pLoudnessMeasure->methodDefinition == MD_PROGRAM_LOUDNESS) {
      _getMethodValue(&valueOrderLoudness, pLoudnessMeasure->methodValue,
                      pLoudnessMeasure->measurementSystem, MS_PROGRAMLOUDNESS);
    }

    if ((dynamicRangeMeasurementType == DRMRT_SHORT_TERM_LOUDNESS_TO_AVG) &&
        (pLoudnessMeasure->methodDefinition == MD_SHORT_TERM_LOUDNESS_MAX)) {
      _getMethodValue(&valueOrderPeakLoudness, pLoudnessMeasure->methodValue,
                      pLoudnessMeasure->measurementSystem, MS_PEAKLOUDNESS);
    }

    if ((dynamicRangeMeasurementType == DRMRT_MOMENTARY_LOUDNESS_TO_AVG) &&
        (pLoudnessMeasure->methodDefinition == MD_MOMENTARY_LOUDNESS_MAX)) {
      _getMethodValue(&valueOrderPeakLoudness, pLoudnessMeasure->methodValue,
                      pLoudnessMeasure->measurementSystem, MS_PEAKLOUDNESS);
    }

    if ((dynamicRangeMeasurementType == DRMRT_TOP_OF_LOUDNESS_RANGE_TO_AVG) &&
        (pLoudnessMeasure->methodDefinition == MD_MAX_OF_LOUDNESS_RANGE)) {
      _getMethodValue(&valueOrderPeakLoudness, pLoudnessMeasure->methodValue,
                      pLoudnessMeasure->measurementSystem, MS_PEAKLOUDNESS);
    }
  }

  if ((valueOrderLoudness.order > -1) && (valueOrderPeakLoudness.order > -1)) {
    *pLoudnessPeakToAverageValue =
        valueOrderPeakLoudness.value - valueOrderLoudness.value;
    *pLoudnessPeakToAverageValuePresent = 1;
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static DRCDEC_SELECTION_PROCESS_RETURN _selectAlbumLoudness(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
    DRCDEC_SELECTION* pCandidatesPotential,
    DRCDEC_SELECTION* pCandidatesSelected) {
  int i, j;

  for (i = 0; i < _drcdec_selection_getNumber(pCandidatesPotential); i++) {
    DRCDEC_SELECTION_DATA* pCandidate =
        _drcdec_selection_getAt(pCandidatesPotential, i);
    if (pCandidate == NULL) return DRCDEC_SELECTION_PROCESS_NOT_OK;

    for (j = 0; j < hLoudnessInfoSet->loudnessInfoAlbumCount; j++) {
      if (pCandidate->pInst->drcSetId ==
          hLoudnessInfoSet->loudnessInfoAlbum[j].drcSetId) {
        if (_drcdec_selection_add(pCandidatesSelected, pCandidate) == NULL)
          return DRCDEC_SELECTION_PROCESS_NOT_OK;
      }
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}

static int _findMethodDefinition(LOUDNESS_INFO* pLoudnessInfo,
                                 int methodDefinition, int startIndex) {
  int i;
  int index = -1;

  for (i = startIndex; i < pLoudnessInfo->measurementCount; i++) {
    if (pLoudnessInfo->loudnessMeasurement[i].methodDefinition ==
        methodDefinition) {
      index = i;
      break;
    }
  }

  return index;
}

static DRCDEC_SELECTION_PROCESS_RETURN _getMixingLevel(
    HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet, int downmixIdRequested,
    int drcSetIdRequested, int albumMode, FIXP_DBL* pMixingLevel) {
  const FIXP_DBL mixingLevelDefault = FL2FXCONST_DBL(85.0f / (float)(1 << 7));

  int i;
  int count;

  LOUDNESS_INFO* pLoudnessInfo = NULL;

  *pMixingLevel = mixingLevelDefault;

  if (drcSetIdRequested < 0) {
    drcSetIdRequested = 0;
  }

  if (albumMode) {
    count = hLoudnessInfoSet->loudnessInfoAlbumCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfoAlbum;
  } else {
    count = hLoudnessInfoSet->loudnessInfoCount;
    pLoudnessInfo = hLoudnessInfoSet->loudnessInfo;
  }

  for (i = 0; i < count; i++) {
    if ((drcSetIdRequested == pLoudnessInfo[i].drcSetId) &&
        ((downmixIdRequested == pLoudnessInfo[i].downmixId) ||
         (DOWNMIX_ID_ANY_DOWNMIX == pLoudnessInfo[i].downmixId))) {
      int index = _findMethodDefinition(&pLoudnessInfo[i], MD_MIXING_LEVEL, 0);

      if (index >= 0) {
        *pMixingLevel = pLoudnessInfo[i].loudnessMeasurement[index].methodValue;
        break;
      }
    }
  }

  return DRCDEC_SELECTION_PROCESS_NO_ERROR;
}
