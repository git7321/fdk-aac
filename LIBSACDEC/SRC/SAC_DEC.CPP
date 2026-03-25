/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Decoder Library

*******************************************************************************/

#include "sac_dec_errorcodes.h"
#include "sac_dec.h"
#include "sac_process.h"
#include "sac_bitdec.h"
#include "sac_smoothing.h"
#include "sac_calcM1andM2.h"
#include "sac_reshapeBBEnv.h"
#include "sac_stp.h"
#include "sac_rom.h"

#include "FDK_decorrelate.h"

#include "FDK_trigFcts.h"
#include "FDK_matrixCalloc.h"

enum {
  APPLY_M2_NONE = 0,
  APPLY_M2 = 1,
  APPLY_M2_MODE212 = 2,
  APPLY_M2_MODE212_Res_PhaseCoding =
      3
};

int FDK_SpatialDecInitDefaultSpatialSpecificConfig(
    SPATIAL_SPECIFIC_CONFIG *pSpatialSpecificConfig,
    AUDIO_OBJECT_TYPE coreCodec, int coreChannels, int samplingFreq,
    int nTimeSlots, int decoderLevel, int isBlind) {
  return SpatialDecDefaultSpecificConfig(pSpatialSpecificConfig, coreCodec,
                                         samplingFreq, nTimeSlots, decoderLevel,
                                         isBlind, coreChannels);
}

int FDK_SpatialDecCompareSpatialSpecificConfigHeader(
    SPATIAL_SPECIFIC_CONFIG *pSsc1, SPATIAL_SPECIFIC_CONFIG *pSsc2) {
  int result = MPS_OK;

  if (FDKmemcmp(pSsc1, pSsc2, sizeof(SPATIAL_SPECIFIC_CONFIG)) != 0) {
    result = MPS_UNEQUAL_SSC;
  }
  return result;
}

static void SpatialDecClearFrameData(
    spatialDec *self,
    SPATIAL_BS_FRAME *bsFrame, const SACDEC_CREATION_PARAMS *const setup) {
  int i;

  FDK_ASSERT(self != NULL);
  FDK_ASSERT(bsFrame != NULL);
  FDK_ASSERT(setup != NULL);

  for (i = 0; i < setup->maxNumOutputChannels;
       i += 1) {
    bsFrame->tempShapeEnableChannelSTP[i] = 0;
    bsFrame->tempShapeEnableChannelGES[i] = 0;
  }

  bsFrame->TsdData->bsTsdEnable = 0;

  bsFrame->numParameterSets = 1;
  bsFrame->paramSlot[0] = self->timeSlots - 1;

  bsFrame->bsSmoothMode[0] = 0;
  initParameterSmoothing(self);

  {
    int resQmfBands, resTimeSlots = (1);

    resQmfBands = setup->maxNumQmfBands;

    for (i = 0; i < setup->bProcResidual
                    ? fMin(setup->maxNumResChannels,
                           setup->maxNumOttBoxes + setup->maxNumInputChannels)
                    : 0;
         i += 1) {
      for (int j = 0; j < resTimeSlots; j += 1) {
        for (int k = 0; k < resQmfBands; k += 1) {
          self->qmfResidualReal__FDK[i][j][k] = FL2FXCONST_DBL(0.0f);
          self->qmfResidualImag__FDK[i][j][k] = FL2FXCONST_DBL(0.0f);
        }
      }
    }
  }

  return;
}

spatialDec *FDK_SpatialDecOpen(const SPATIAL_DEC_CONFIG *config,
                               int stereoConfigIndex) {
  int i;
  int lfSize, hfSize;
  spatialDec *self = NULL;
  SACDEC_CREATION_PARAMS setup;

  switch (config->decoderLevel) {
    case DECODER_LEVEL_0:
      setup.maxNumInputChannels = 1;
      setup.maxNumOutputChannels = 2;
      setup.maxNumQmfBands = 64;
      setup.maxNumXChannels = 2;
      setup.maxNumVChannels = 2;
      setup.maxNumDecorChannels = 1;
      setup.bProcResidual = 1;
      setup.maxNumResidualChannels = 0;
      setup.maxNumOttBoxes = 1;
      setup.maxNumParams = setup.maxNumInputChannels + setup.maxNumOttBoxes;
      break;
    default:
      return NULL;
  }

  setup.maxNumResChannels = 1;

  {
    switch (config->maxNumOutputChannels) {
      case OUTPUT_CHANNELS_2_0:
        setup.maxNumOutputChannels = fMin(setup.maxNumOutputChannels, 2);
        break;
      case OUTPUT_CHANNELS_DEFAULT:
      default:
        break;
    }
  }

  setup.maxNumHybridBands = SacGetHybridSubbands(setup.maxNumQmfBands);

  switch (config->decoderMode) {
    case EXT_HQ_ONLY:
      setup.maxNumCmplxQmfBands = setup.maxNumQmfBands;
      setup.maxNumCmplxHybBands = setup.maxNumHybridBands;
      break;
    default:
      setup.maxNumCmplxQmfBands = fixMax(PC_NUM_BANDS, setup.maxNumQmfBands);
      setup.maxNumCmplxHybBands =
          fixMax(PC_NUM_HYB_BANDS, setup.maxNumHybridBands);
      break;
  }

  FDK_ALLOCATE_MEMORY_1D_INT(self, 1, spatialDec, SECT_DATA_L2)

  self->createParams = setup;

  FDK_ALLOCATE_MEMORY_1D(self->param2hyb, MAX_PARAMETER_BANDS + 1, int)

  FDK_ALLOCATE_MEMORY_1D(self->numOttBands, setup.maxNumOttBoxes, int)

  FDK_ALLOCATE_MEMORY_1D(self->smgTime, MAX_PARAMETER_SETS, int)
  FDK_ALLOCATE_MEMORY_2D(self->smgData, MAX_PARAMETER_SETS, MAX_PARAMETER_BANDS,
                         UCHAR)

  FDK_ALLOCATE_MEMORY_3D(self->ottCLD__FDK, setup.maxNumOttBoxes,
                         MAX_PARAMETER_SETS, MAX_PARAMETER_BANDS, SCHAR)
  FDK_ALLOCATE_MEMORY_3D(self->ottICC__FDK, setup.maxNumOttBoxes,
                         MAX_PARAMETER_SETS, MAX_PARAMETER_BANDS, SCHAR)
  FDK_ALLOCATE_MEMORY_3D(self->ottIPD__FDK, setup.maxNumOttBoxes,
                         MAX_PARAMETER_SETS, MAX_PARAMETER_BANDS, SCHAR)

  FDK_ALLOCATE_MEMORY_2D(self->ottCLDidxPrev, setup.maxNumOttBoxes,
                         MAX_PARAMETER_BANDS, SCHAR)
  FDK_ALLOCATE_MEMORY_2D(self->ottICCidxPrev, setup.maxNumOttBoxes,
                         MAX_PARAMETER_BANDS, SCHAR)
  FDK_ALLOCATE_MEMORY_3D(self->ottICCdiffidx, setup.maxNumOttBoxes,
                         MAX_PARAMETER_SETS, MAX_PARAMETER_BANDS, SCHAR)
  FDK_ALLOCATE_MEMORY_2D(self->ottIPDidxPrev, setup.maxNumOttBoxes,
                         MAX_PARAMETER_BANDS, SCHAR)
  FDK_ALLOCATE_MEMORY_2D(self->arbdmxGainIdxPrev, setup.maxNumInputChannels,
                         MAX_PARAMETER_BANDS, SCHAR)
  FDK_ALLOCATE_MEMORY_2D(self->cmpOttCLDidxPrev, setup.maxNumOttBoxes,
                         MAX_PARAMETER_BANDS, SCHAR)
  FDK_ALLOCATE_MEMORY_2D(self->cmpOttICCidxPrev, setup.maxNumOttBoxes,
                         MAX_PARAMETER_BANDS, SCHAR)
  FDK_ALLOCATE_MEMORY_3D(self->outIdxData, setup.maxNumOttBoxes,
                         MAX_PARAMETER_SETS, MAX_PARAMETER_BANDS, SCHAR)

  FDK_ALLOCATE_MEMORY_3D(self->arbdmxGain__FDK, setup.maxNumInputChannels,
                         MAX_PARAMETER_SETS, MAX_PARAMETER_BANDS, SCHAR)
  FDK_ALLOCATE_MEMORY_1D(self->arbdmxAlpha__FDK, setup.maxNumInputChannels,
                         FIXP_DBL)
  FDK_ALLOCATE_MEMORY_1D(self->arbdmxAlphaPrev__FDK, setup.maxNumInputChannels,
                         FIXP_DBL)
  FDK_ALLOCATE_MEMORY_2D(self->cmpArbdmxGainIdxPrev, setup.maxNumInputChannels,
                         MAX_PARAMETER_BANDS, SCHAR)

  FDK_ALLOCATE_MEMORY_2D(self->cmpOttIPDidxPrev, setup.maxNumOttBoxes,
                         MAX_PARAMETER_BANDS, SCHAR)

  FDK_ALLOCATE_MEMORY_3D_INT(self->M2Real__FDK, setup.maxNumOutputChannels,
                             setup.maxNumVChannels, MAX_PARAMETER_BANDS,
                             FIXP_DBL, SECT_DATA_L2)
  FDK_ALLOCATE_MEMORY_3D(self->M2Imag__FDK, setup.maxNumOutputChannels,
                         setup.maxNumVChannels, MAX_PARAMETER_BANDS, FIXP_DBL)

  FDK_ALLOCATE_MEMORY_3D_INT(self->M2RealPrev__FDK, setup.maxNumOutputChannels,
                             setup.maxNumVChannels, MAX_PARAMETER_BANDS,
                             FIXP_DBL, SECT_DATA_L2)
  FDK_ALLOCATE_MEMORY_3D(self->M2ImagPrev__FDK, setup.maxNumOutputChannels,
                         setup.maxNumVChannels, MAX_PARAMETER_BANDS, FIXP_DBL)

  FDK_ALLOCATE_MEMORY_2D_INT_ALIGNED(
      self->qmfInputReal__FDK, setup.maxNumInputChannels, setup.maxNumQmfBands,
      FIXP_DBL, SECT_DATA_L2)
  FDK_ALLOCATE_MEMORY_2D_INT_ALIGNED(
      self->qmfInputImag__FDK, setup.maxNumInputChannels,
      setup.maxNumCmplxQmfBands, FIXP_DBL, SECT_DATA_L2)

  FDK_ALLOCATE_MEMORY_2D_INT(self->hybInputReal__FDK, setup.maxNumInputChannels,
                             setup.maxNumHybridBands, FIXP_DBL, SECT_DATA_L2)
  FDK_ALLOCATE_MEMORY_2D_INT(self->hybInputImag__FDK, setup.maxNumInputChannels,
                             setup.maxNumCmplxHybBands, FIXP_DBL, SECT_DATA_L2)

  if (setup.bProcResidual) {
    FDK_ALLOCATE_MEMORY_1D(self->qmfResidualReal__FDK, setup.maxNumResChannels,
                           FIXP_DBL **)
    FDK_ALLOCATE_MEMORY_1D(self->qmfResidualImag__FDK, setup.maxNumResChannels,
                           FIXP_DBL **)

    FDK_ALLOCATE_MEMORY_1D(self->hybResidualReal__FDK, setup.maxNumResChannels,
                           FIXP_DBL *)
    FDK_ALLOCATE_MEMORY_1D(self->hybResidualImag__FDK, setup.maxNumResChannels,
                           FIXP_DBL *)

    for (i = 0; i < setup.maxNumResChannels; i++) {
      int resQmfBands = (config->decoderMode == EXT_LP_ONLY)
                            ? PC_NUM_BANDS
                            : setup.maxNumQmfBands;
      int resHybBands = (config->decoderMode == EXT_LP_ONLY)
                            ? PC_NUM_HYB_BANDS
                            : setup.maxNumHybridBands;
      FDK_ALLOCATE_MEMORY_2D_INT_ALIGNED(self->qmfResidualReal__FDK[i], (1),
                                         resQmfBands, FIXP_DBL, SECT_DATA_L1)
      FDK_ALLOCATE_MEMORY_2D_INT_ALIGNED(self->qmfResidualImag__FDK[i], (1),
                                         resQmfBands, FIXP_DBL, SECT_DATA_L1)

      FDK_ALLOCATE_MEMORY_1D(self->hybResidualReal__FDK[i],
                             setup.maxNumHybridBands, FIXP_DBL)
      FDK_ALLOCATE_MEMORY_1D(self->hybResidualImag__FDK[i], resHybBands,
                             FIXP_DBL)
    }
  }

  FDK_ALLOCATE_MEMORY_2D_INT(self->wReal__FDK, setup.maxNumVChannels,
                             setup.maxNumHybridBands, FIXP_DBL, SECT_DATA_L2)
  FDK_ALLOCATE_MEMORY_2D_INT(self->wImag__FDK, setup.maxNumVChannels,
                             setup.maxNumCmplxHybBands, FIXP_DBL, SECT_DATA_L2)

  FDK_ALLOCATE_MEMORY_2D_INT(self->hybOutputRealDry__FDK,
                             setup.maxNumOutputChannels,
                             setup.maxNumHybridBands, FIXP_DBL, SECT_DATA_L2)
  FDK_ALLOCATE_MEMORY_2D_INT(self->hybOutputImagDry__FDK,
                             setup.maxNumOutputChannels,
                             setup.maxNumCmplxHybBands, FIXP_DBL, SECT_DATA_L2)

  FDK_ALLOCATE_MEMORY_2D_INT(self->hybOutputRealWet__FDK,
                             setup.maxNumOutputChannels,
                             setup.maxNumHybridBands, FIXP_DBL, SECT_DATA_L2)
  FDK_ALLOCATE_MEMORY_2D_INT(self->hybOutputImagWet__FDK,
                             setup.maxNumOutputChannels,
                             setup.maxNumCmplxHybBands, FIXP_DBL, SECT_DATA_L2)

  FDK_ALLOCATE_MEMORY_1D(self->hybridSynthesis, setup.maxNumOutputChannels,
                         FDK_SYN_HYB_FILTER)

  FDK_ALLOCATE_MEMORY_1D(
      self->hybridAnalysis,
      setup.bProcResidual ? setup.maxNumInputChannels + setup.maxNumResChannels
                          : setup.maxNumInputChannels,
      FDK_ANA_HYB_FILTER)

  lfSize = 2 * BUFFER_LEN_LF * MAX_QMF_BANDS_TO_HYBRID;
  {
    hfSize =
        BUFFER_LEN_HF * ((setup.maxNumQmfBands - MAX_QMF_BANDS_TO_HYBRID) +
                         (setup.maxNumCmplxQmfBands - MAX_QMF_BANDS_TO_HYBRID));
  }

  FDK_ALLOCATE_MEMORY_2D_INT(self->pHybridAnaStatesLFdmx,
                             setup.maxNumInputChannels, lfSize, FIXP_DBL,
                             SECT_DATA_L2) {
    FDK_ALLOCATE_MEMORY_2D(self->pHybridAnaStatesHFdmx,
                           setup.maxNumInputChannels, hfSize, FIXP_DBL)
  }

  for (i = 0; i < setup.maxNumInputChannels; i++) {
    FIXP_DBL *pHybridAnaStatesHFdmx;

    pHybridAnaStatesHFdmx = self->pHybridAnaStatesHFdmx[i];

    FDKhybridAnalysisOpen(&self->hybridAnalysis[i],
                          self->pHybridAnaStatesLFdmx[i],
                          lfSize * sizeof(FIXP_DBL), pHybridAnaStatesHFdmx,
                          hfSize * sizeof(FIXP_DBL));
  }
  if (setup.bProcResidual) {
    lfSize = 2 * BUFFER_LEN_LF * MAX_QMF_BANDS_TO_HYBRID;
    hfSize = BUFFER_LEN_HF *
             ((((config->decoderMode == EXT_LP_ONLY) ? PC_NUM_BANDS
                                                     : setup.maxNumQmfBands) -
               MAX_QMF_BANDS_TO_HYBRID) +
              (setup.maxNumCmplxQmfBands - MAX_QMF_BANDS_TO_HYBRID));

    FDK_ALLOCATE_MEMORY_2D_INT(self->pHybridAnaStatesLFres,
                               setup.maxNumResChannels, lfSize, FIXP_DBL,
                               SECT_DATA_L2)
    FDK_ALLOCATE_MEMORY_2D(self->pHybridAnaStatesHFres, setup.maxNumResChannels,
                           hfSize, FIXP_DBL)

    for (i = setup.maxNumInputChannels;
         i < (setup.maxNumInputChannels + setup.maxNumResChannels); i++) {
      FDKhybridAnalysisOpen(
          &self->hybridAnalysis[i],
          self->pHybridAnaStatesLFres[i - setup.maxNumInputChannels],
          lfSize * sizeof(FIXP_DBL),
          self->pHybridAnaStatesHFres[i - setup.maxNumInputChannels],
          hfSize * sizeof(FIXP_DBL));
    }
  }

  FDK_ALLOCATE_MEMORY_1D(self->smoothState, 1, SMOOTHING_STATE)
  FDK_ALLOCATE_MEMORY_1D(self->reshapeBBEnvState, 1, RESHAPE_BBENV_STATE)

  FDK_ALLOCATE_MEMORY_1D(self->apDecor, setup.maxNumDecorChannels, DECORR_DEC)
  FDK_ALLOCATE_MEMORY_2D_INT(self->pDecorBufferCplx, setup.maxNumDecorChannels,
                             (2 * ((825) + (373))), FIXP_DBL, SECT_DATA_L2)

  for (i = 0; i < setup.maxNumDecorChannels; i++) {
    if (FDKdecorrelateOpen(&self->apDecor[i], self->pDecorBufferCplx[i],
                           (2 * ((825) + (373))))) {
      goto bail;
    }
  }

  if (subbandTPCreate(&self->hStpDec) != MPS_OK) {
    goto bail;
  }

  self->decoderLevel = config->decoderLevel;
  self->decoderMode = config->decoderMode;
  self->binauralMode = config->binauralMode;

  self->partiallyComplex = (config->decoderMode != EXT_HQ_ONLY) ? 1 : 0;

  SpatialDecConcealment_Init(&self->concealInfo, MPEGS_CONCEAL_RESET_ALL);

  return self;

bail:
  FDK_SpatialDecClose(self);

  return NULL;
}

static int isValidConfig(spatialDec const *const self,
                         const SPATIAL_DEC_UPMIX_TYPE upmixType,
                         SPATIALDEC_PARAM const *const pUserParams,
                         const AUDIO_OBJECT_TYPE coreAot) {
  UPMIXTYPE nUpmixType;

  FDK_ASSERT(self != NULL);
  FDK_ASSERT(pUserParams != NULL);

  nUpmixType = (UPMIXTYPE)upmixType;

  switch (nUpmixType) {
    case UPMIXTYPE_BYPASS:
      break;
    case UPMIXTYPE_NORMAL:
      break;
    default:
      return 0;
  }

  return 1;
}

static SACDEC_ERROR CheckLevelTreeUpmixType(
    const SACDEC_CREATION_PARAMS *const pCreateParams,
    const SPATIAL_SPECIFIC_CONFIG *const pSsc, const int decoderLevel,
    const UPMIXTYPE upmixType) {
  SACDEC_ERROR err = MPS_OK;
  int nOutputChannels, treeConfig;

  FDK_ASSERT(pCreateParams != NULL);
  FDK_ASSERT(pSsc != NULL);

  treeConfig = pSsc->treeConfig;

  switch (decoderLevel) {
    case 0: {
      if (treeConfig != SPATIALDEC_MODE_RSVD7) {
        err = MPS_INVALID_TREECONFIG;
        goto bail;
      }
      break;
    }
    default:
      err = MPS_INVALID_PARAMETER;
      goto bail;
  }

  switch (upmixType) {
    case UPMIXTYPE_BYPASS:
      nOutputChannels = pSsc->nInputChannels;
      break;
    default:
      nOutputChannels = pSsc->nOutputChannels;
      break;
  }

  if ((pSsc->nInputChannels > pCreateParams->maxNumInputChannels) ||
      (nOutputChannels > pCreateParams->maxNumOutputChannels) ||
      (pSsc->nOttBoxes > pCreateParams->maxNumOttBoxes)) {
    err = MPS_INVALID_PARAMETER;
  }

bail:
  return err;
}

void SpatialDecInitParserContext(spatialDec *self) {
  int i, j;

  for (i = 0; i < self->createParams.maxNumOttBoxes; i += 1) {
    for (j = 0; j < MAX_PARAMETER_BANDS; j++) {
      self->ottCLDidxPrev[i][j] = 0;
      self->ottICCidxPrev[i][j] = 0;
      self->cmpOttCLDidxPrev[i][j] = 0;
      self->cmpOttICCidxPrev[i][j] = 0;
    }
  }
  for (i = 0; i < self->createParams.maxNumInputChannels; i++) {
    for (j = 0; j < MAX_PARAMETER_BANDS; j++) {
      self->arbdmxGainIdxPrev[i][j] = 0;
      self->cmpArbdmxGainIdxPrev[i][j] = 0;
    }
  }
}

SACDEC_ERROR FDK_SpatialDecInit(spatialDec *self, SPATIAL_BS_FRAME *frame,
                                SPATIAL_SPECIFIC_CONFIG *pSpatialSpecificConfig,
                                int nQmfBands,
                                SPATIAL_DEC_UPMIX_TYPE const upmixType,
                                SPATIALDEC_PARAM *pUserParams, UINT initFlags) {
  SACDEC_ERROR err = MPS_OK;
  int nCh, i, j, k;
  int maxQmfBands;
  int bypassMode = 0;

  self->useFDreverb = 0;

  if (!isValidConfig(self, upmixType, pUserParams,
                     pSpatialSpecificConfig->coreCodec)) {
    return MPS_INVALID_PARAMETER;
  }

  err = CheckLevelTreeUpmixType(&self->createParams, pSpatialSpecificConfig,
                                self->decoderLevel, (UPMIXTYPE)upmixType);
  if (err != MPS_OK) {
    goto bail;
  }

  self->upmixType = (UPMIXTYPE)upmixType;

  if (initFlags & MPEGS_INIT_PARAMS_ERROR_CONCEALMENT) {
    err = SpatialDecConcealment_SetParam(
        &self->concealInfo, SAC_DEC_CONCEAL_METHOD, pUserParams->concealMethod);
    if (err != MPS_OK) {
      goto bail;
    }
    err = SpatialDecConcealment_SetParam(&self->concealInfo,
                                         SAC_DEC_CONCEAL_NUM_KEEP_FRAMES,
                                         pUserParams->concealNumKeepFrames);
    if (err != MPS_OK) {
      goto bail;
    }
    err = SpatialDecConcealment_SetParam(
        &self->concealInfo, SAC_DEC_CONCEAL_FADE_OUT_SLOPE_LENGTH,
        pUserParams->concealFadeOutSlopeLength);
    if (err != MPS_OK) {
      goto bail;
    }
    err = SpatialDecConcealment_SetParam(&self->concealInfo,
                                         SAC_DEC_CONCEAL_FADE_IN_SLOPE_LENGTH,
                                         pUserParams->concealFadeInSlopeLength);
    if (err != MPS_OK) {
      goto bail;
    }
    err = SpatialDecConcealment_SetParam(&self->concealInfo,
                                         SAC_DEC_CONCEAL_NUM_RELEASE_FRAMES,
                                         pUserParams->concealNumReleaseFrames);
    if (err != MPS_OK) {
      goto bail;
    }
  }

  if (initFlags &
      MPEGS_INIT_STATES_ERROR_CONCEALMENT) {
    SpatialDecConcealment_Init(&self->concealInfo, MPEGS_CONCEAL_RESET_STATE);
  }

  bypassMode |= pUserParams->bypassMode;
  bypassMode |= ((self->upmixType == UPMIXTYPE_BYPASS) ? 1 : 0);

  switch (nQmfBands) {
    case 16:
    case 24:
    case 32:
      self->staticDecScale = 21;
      break;
    case 64:
      self->staticDecScale = 22;
      break;
    default:
      return MPS_INVALID_PARAMETER;
  }

  self->numParameterSetsPrev = 1;

  self->qmfBands = nQmfBands;

  self->bShareDelayWithSBR = 0;

  err = SpatialDecDecodeHeader(self, pSpatialSpecificConfig);
  if (err != MPS_OK) {
    goto bail;
  }

  self->stereoConfigIndex = pSpatialSpecificConfig->stereoConfigIndex;

  if (initFlags & MPEGS_INIT_STATES_ANA_QMF_FILTER) {
    self->qmfInputDelayBufPos = 0;
    self->pc_filterdelay = 1;
  }

  maxQmfBands = self->qmfBands;

  if (initFlags & MPEGS_INIT_STATES_PARAM) {
    initParameterSmoothing(self);
  }

  initBBEnv(self, (initFlags & MPEGS_INIT_STATES_GES) ? 1 : 0);

  if (!isTwoChMode(self->upmixType) && !bypassMode) {
    self->staticDecScale += self->clipProtectGainSF__FDK;
  }

  {
    UINT flags = 0;
    INT initStatesFlag = (initFlags & MPEGS_INIT_STATES_ANA_QMF_FILTER) ? 1 : 0;
    INT useLdFilter =
        (self->pConfigCurrent->syntaxFlags & SACDEC_SYNTAX_LD) ? 1 : 0;

    flags = self->pQmfDomain->globalConf.flags_requested;
    flags &= (~(UINT)QMF_FLAG_LP);

    if (initStatesFlag)
      flags &= ~QMF_FLAG_KEEP_STATES;
    else
      flags |= QMF_FLAG_KEEP_STATES;

    if (useLdFilter)
      flags |= QMF_FLAG_MPSLDFB;
    else
      flags &= ~QMF_FLAG_MPSLDFB;

    self->pQmfDomain->globalConf.flags_requested = flags;
    FDK_QmfDomain_Configure(self->pQmfDomain);

    for (nCh = 0; nCh < self->numOutputChannelsAT; nCh++) {
      int outputScale = 0, outputGain_e = 0, scale = -(8) + (1);
      FIXP_DBL outputGain_m = getChGain(self, nCh, &outputGain_e);

      if (!isTwoChMode(self->upmixType) && !bypassMode) {
        outputScale +=
            self->clipProtectGainSF__FDK;
      }

      scale += outputScale;

      qmfChangeOutScalefactor(&self->pQmfDomain->QmfDomainOut[nCh].fb, scale);
      qmfChangeOutGain(&self->pQmfDomain->QmfDomainOut[nCh].fb, outputGain_m,
                       outputGain_e);
    }
  }

  for (nCh = 0; nCh < self->numOutputChannelsAT; nCh++) {
    FDKhybridSynthesisInit(&self->hybridSynthesis[nCh], THREE_TO_TEN,
                           self->qmfBands, maxQmfBands);
  }

  for (nCh = 0; nCh < self->createParams.maxNumInputChannels; nCh++) {
    FDKhybridAnalysisInit(
        &self->hybridAnalysis[nCh], THREE_TO_TEN, self->qmfBands, maxQmfBands,
        (initFlags & MPEGS_INIT_STATES_ANA_HYB_FILTER) ? 1 : 0);
  }
  for (; nCh < (self->createParams.bProcResidual
                    ? (self->createParams.maxNumInputChannels +
                       self->createParams.maxNumResChannels)
                    : self->createParams.maxNumInputChannels);
       nCh++) {
    FDKhybridAnalysisInit(&self->hybridAnalysis[nCh], THREE_TO_TEN, maxQmfBands,
                          maxQmfBands, 0);
  }

  {
    for (k = 0; k < self->numDecorSignals; k++) {
      int errCode, idec;
      FDK_DECORR_TYPE decorrType = DECORR_PS;
      decorrType = DECORR_LD;
      if (self->pConfigCurrent->syntaxFlags &
          (SACDEC_SYNTAX_USAC | SACDEC_SYNTAX_RSVD50)) {
        decorrType =
            ((self->treeConfig == TREE_212) && (self->decorrType == DECORR_PS))
                ? DECORR_PS
                : DECORR_USAC;
      }
      {
        idec = k;
        if (self->pConfigCurrent->syntaxFlags & SACDEC_SYNTAX_LD) {
          if (self->treeConfig == TREE_212 && k == 0) {
            idec = 2;
          }
        }
      }
      errCode = FDKdecorrelateInit(
          &self->apDecor[k], self->hybridBands, decorrType, DUCKER_AUTOMATIC,
          self->decorrConfig, idec, 0,
          0, 0,
          (initFlags & MPEGS_INIT_STATES_DECORRELATOR) ? 1 : 0);
      if (errCode) return MPS_NOTOK;
    }
  }

  err = initM1andM2(self, (initFlags & MPEGS_INIT_STATES_M1M2) ? 1 : 0,
                    (initFlags & MPEGS_INIT_CONFIG) ? 1 : 0);
  if (err != MPS_OK) return err;

  if (initFlags & MPEGS_INIT_STATES_PARAM) {
    for (i = 0; i < self->createParams.maxNumOttBoxes; i += 1) {
      for (k = 0; k < MAX_PARAMETER_SETS; k += 1) {
        for (j = 0; j < MAX_PARAMETER_BANDS; j += 1) {
          self->ottICCdiffidx[i][k][j] = 0;
        }
      }
    }
    self->smoothState->prevSmgTime = smgTimeTable[2];
    FDKmemclear(self->smoothState->prevSmgData,
                MAX_PARAMETER_BANDS * sizeof(UCHAR));
    FDKmemclear(self->smoothState->opdLeftState__FDK,
                MAX_PARAMETER_BANDS * sizeof(FIXP_DBL));
    FDKmemclear(self->smoothState->opdRightState__FDK,
                MAX_PARAMETER_BANDS * sizeof(FIXP_DBL));
  }

  self->prevTimeSlot = -1;
  self->curTimeSlot =
      MAX_TIME_SLOTS + 1;
  self->curPs = 0;

  subbandTPInit(self->hStpDec);

bail:
  return err;
}

void SpatialDecChannelProperties(spatialDec *self,
                                 AUDIO_CHANNEL_TYPE channelType[],
                                 UCHAR channelIndices[],
                                 const FDK_channelMapDescr *const mapDescr) {
  if ((self == NULL) || (channelType == NULL) || (channelIndices == NULL) ||
      (mapDescr == NULL)) {
    return;
  }

  if (self->numOutputChannelsAT !=
      treePropertyTable[self->treeConfig].numOutputChannels) {
    int ch;
    for (ch = 0; ch < self->numOutputChannelsAT; ch += 1) {
      channelType[ch] = ACT_FRONT;
      channelIndices[ch] = ch;
    }
  } else {
    switch (self->treeConfig) {
      case TREE_212:
        channelType[0] = ACT_FRONT;
        channelIndices[0] = 0;
        channelType[1] = ACT_FRONT;
        channelIndices[1] = 1;
        break;
      default:;
    }
  }
}

void FDK_SpatialDecClose(spatialDec *self) {
  if (self) {
    int k;

    if (self->apDecor != NULL) {
      for (k = 0; k < self->createParams.maxNumDecorChannels; k++) {
        FDKdecorrelateClose(&(self->apDecor[k]));
      }
      FDK_FREE_MEMORY_1D(self->apDecor);
    }
    if (self->pDecorBufferCplx != NULL) {
      FDK_FREE_MEMORY_2D(self->pDecorBufferCplx);
    }

    subbandTPDestroy(&self->hStpDec);

    FDK_FREE_MEMORY_1D(self->reshapeBBEnvState);
    FDK_FREE_MEMORY_1D(self->smoothState);

    FDK_FREE_MEMORY_2D(self->pHybridAnaStatesLFdmx);
    FDK_FREE_MEMORY_2D(self->pHybridAnaStatesHFdmx);
    FDK_FREE_MEMORY_2D(self->pHybridAnaStatesLFres);
    FDK_FREE_MEMORY_2D(self->pHybridAnaStatesHFres);
    FDK_FREE_MEMORY_1D(self->hybridAnalysis);

    FDK_FREE_MEMORY_1D(self->hybridSynthesis);

    FDK_FREE_MEMORY_2D(self->hybOutputImagWet__FDK);
    FDK_FREE_MEMORY_2D(self->hybOutputRealWet__FDK);

    FDK_FREE_MEMORY_2D(self->hybOutputImagDry__FDK);
    FDK_FREE_MEMORY_2D(self->hybOutputRealDry__FDK);

    FDK_FREE_MEMORY_2D(self->wImag__FDK);
    FDK_FREE_MEMORY_2D(self->wReal__FDK);

    if (self->createParams.bProcResidual) {
      int i;

      for (i = 0; i < self->createParams.maxNumResChannels; i++) {
        if (self->hybResidualImag__FDK != NULL)
          FDK_FREE_MEMORY_1D(self->hybResidualImag__FDK[i]);
        if (self->hybResidualReal__FDK != NULL)
          FDK_FREE_MEMORY_1D(self->hybResidualReal__FDK[i]);
        if (self->qmfResidualImag__FDK != NULL)
          FDK_FREE_MEMORY_2D_ALIGNED(self->qmfResidualImag__FDK[i]);
        if (self->qmfResidualReal__FDK != NULL)
          FDK_FREE_MEMORY_2D_ALIGNED(self->qmfResidualReal__FDK[i]);
      }

      FDK_FREE_MEMORY_1D(self->hybResidualImag__FDK);
      FDK_FREE_MEMORY_1D(self->hybResidualReal__FDK);

      FDK_FREE_MEMORY_1D(self->qmfResidualImag__FDK);
      FDK_FREE_MEMORY_1D(self->qmfResidualReal__FDK);

    }

    FDK_FREE_MEMORY_2D(self->hybInputImag__FDK);
    FDK_FREE_MEMORY_2D(self->hybInputReal__FDK);

    FDK_FREE_MEMORY_2D_ALIGNED(self->qmfInputImag__FDK);
    FDK_FREE_MEMORY_2D_ALIGNED(self->qmfInputReal__FDK);

    FDK_FREE_MEMORY_3D(self->M2ImagPrev__FDK);

    FDK_FREE_MEMORY_3D(self->M2RealPrev__FDK);

    FDK_FREE_MEMORY_3D(self->M2Imag__FDK);

    FDK_FREE_MEMORY_3D(self->M2Real__FDK);

    FDK_FREE_MEMORY_1D(self->arbdmxAlphaPrev__FDK);
    FDK_FREE_MEMORY_1D(self->arbdmxAlpha__FDK);

    FDK_FREE_MEMORY_3D(self->arbdmxGain__FDK);

    FDK_FREE_MEMORY_3D(self->ottIPD__FDK);
    FDK_FREE_MEMORY_3D(self->ottICC__FDK);
    FDK_FREE_MEMORY_3D(self->ottCLD__FDK);

    FDK_FREE_MEMORY_2D(self->ottCLDidxPrev);
    FDK_FREE_MEMORY_2D(self->ottICCidxPrev);
    FDK_FREE_MEMORY_3D(self->ottICCdiffidx);
    FDK_FREE_MEMORY_2D(self->ottIPDidxPrev);
    FDK_FREE_MEMORY_2D(self->arbdmxGainIdxPrev);

    FDK_FREE_MEMORY_2D(self->cmpOttCLDidxPrev);
    FDK_FREE_MEMORY_2D(self->cmpOttICCidxPrev);
    FDK_FREE_MEMORY_3D(self->outIdxData);
    FDK_FREE_MEMORY_2D(self->cmpOttIPDidxPrev);
    FDK_FREE_MEMORY_2D(self->cmpArbdmxGainIdxPrev);

    FDK_FREE_MEMORY_2D(self->smgData);
    FDK_FREE_MEMORY_1D(self->smgTime);

    FDK_FREE_MEMORY_1D(self->numOttBands);

    FDK_FREE_MEMORY_1D(self->param2hyb);

    FDK_FREE_MEMORY_1D(self);
  }

  return;
}

static void SpatialDecApplyBypass(spatialDec *self, FIXP_DBL **hybInputReal,
                                  FIXP_DBL **hybInputImag,
                                  FIXP_DBL **hybOutputReal,
                                  FIXP_DBL **hybOutputImag,
                                  const int numInputChannels) {
  int complexHybBands;

  complexHybBands = self->hybridBands;

  {
    int ch;
    int rf = -1, lf = -1, cf = -1;

    switch (self->treeConfig) {
      case TREE_212:
        lf = 0;
        rf = 1;
        break;
      default:;
    }

    switch (numInputChannels) {
      case 1:
        if (cf > 0) {
          FDKmemcpy(hybOutputReal[cf], hybInputReal[0],
                    self->hybridBands * sizeof(FIXP_DBL));
          FDKmemcpy(hybOutputImag[cf], hybInputImag[0],
                    complexHybBands * sizeof(FIXP_DBL));
        } else {
          FDKmemcpy(hybOutputReal[lf], hybInputReal[0],
                    self->hybridBands * sizeof(FIXP_DBL));
          FDKmemcpy(hybOutputReal[rf], hybInputReal[0],
                    self->hybridBands * sizeof(FIXP_DBL));
          FDKmemcpy(hybOutputImag[lf], hybInputImag[0],
                    complexHybBands * sizeof(FIXP_DBL));
          FDKmemcpy(hybOutputImag[rf], hybInputImag[0],
                    complexHybBands * sizeof(FIXP_DBL));
        }
        break;
      case 2:
        FDK_ASSERT(lf != -1);
        FDK_ASSERT(rf != -1);
        FDKmemcpy(hybOutputReal[lf], hybInputReal[0],
                  self->hybridBands * sizeof(FIXP_DBL));
        FDKmemcpy(hybOutputReal[rf], hybInputReal[1],
                  self->hybridBands * sizeof(FIXP_DBL));
        FDKmemcpy(hybOutputImag[lf], hybInputImag[0],
                  complexHybBands * sizeof(FIXP_DBL));
        FDKmemcpy(hybOutputImag[rf], hybInputImag[1],
                  complexHybBands * sizeof(FIXP_DBL));
        break;
    }
    for (ch = 0; ch < self->numOutputChannelsAT; ch++) {
      if (ch == lf || ch == rf || ch == cf) {
        continue;
      }
      FDKmemclear(hybOutputReal[ch], self->hybridBands * sizeof(FIXP_DBL));
      FDKmemclear(hybOutputImag[ch], complexHybBands * sizeof(FIXP_DBL));
    }
  }
}

static SACDEC_ERROR SpatialDecSetInternalError(spatialDec *self,
                                               int *bypassMode,
                                               SACDEC_ERROR err) {
  *bypassMode = 1;

  if (self->errInt == MPS_OK) {
    self->errInt = err;
  }

  return MPS_OK;
}

static SACDEC_ERROR SpatialDecApplyParameterSets(
    spatialDec *self, const SPATIAL_BS_FRAME *frame, SPATIALDEC_INPUT_MODE mode,
    PCM_MPS *inData,
    FIXP_DBL **qmfInDataReal,
    FIXP_DBL **qmfInDataImag,
    UINT nSamples, UINT controlFlags, int numInputChannels,
    const FDK_channelMapDescr *const mapDescr) {
  SACDEC_ERROR err = MPS_OK;

  FIXP_SGL alpha = FL2FXCONST_SGL(0.0);

  int ts;
  int ch;
  int hyb;

  int prevSlot = self->prevTimeSlot;
  int ps = self->curPs;
  int ts_io = 0;
  int bypassMode = (controlFlags & MPEGS_BYPASSMODE) ? 1 : 0;

  bypassMode |= ((self->upmixType == UPMIXTYPE_BYPASS) ? 1 : 0);

  for (ts = self->curTimeSlot;
       ts <= fixMin(self->curTimeSlot + (int)nSamples / self->qmfBands - 1,
                    self->timeSlots - 1);
       ts++, ts_io++) {
    int currSlot = frame->paramSlot[ps];

    err = (currSlot < ts) ? MPS_WRONG_PARAMETERSETS : MPS_OK;
    if (err != MPS_OK) {
      err = SpatialDecSetInternalError(self, &bypassMode, err);
    }

    if (ts == prevSlot + 1) {
      if (bypassMode == 0) {
        err = SpatialDecCalculateM1andM2(
            self, ps, frame);
        if (err != MPS_OK) {
          err = SpatialDecSetInternalError(self, &bypassMode, err);
        }
      }

      if ((ps == 0) && (self->bOverwriteM1M2prev != 0)) {
        SpatialDecBufferMatrices(self);
        self->bOverwriteM1M2prev = 0;
      }

      if (bypassMode == 0) {
        SpatialDecSmoothM1andM2(
            self, frame,
            ps);
      }
    }

    if (bypassMode == 0) {
      alpha = FX_DBL2FX_SGL(fDivNorm(ts - prevSlot, currSlot - prevSlot));
    }

    switch (mode) {
      case INPUTMODE_QMF_SBR:
        if (self->pConfigCurrent->syntaxFlags & SACDEC_SYNTAX_LD)
          self->bShareDelayWithSBR = 0;
        else
          self->bShareDelayWithSBR = 1;
        SpatialDecFeedQMF(
            self, qmfInDataReal, qmfInDataImag, ts_io, bypassMode,
            self->qmfInputReal__FDK, self->qmfInputImag__FDK,
            (bypassMode) ? numInputChannels : self->numInputChannels);
        break;
      case INPUTMODE_TIME:
        self->bShareDelayWithSBR = 0;
        SpatialDecQMFAnalysis(
            self, inData, ts_io, bypassMode, self->qmfInputReal__FDK,
            self->qmfInputImag__FDK,
            (bypassMode) ? numInputChannels : self->numInputChannels);
        break;
      default:
        break;
    }

    if ((self->pConfigCurrent->syntaxFlags & SACDEC_SYNTAX_USAC) &&
        self->residualCoding) {
      int offset;
      ch = 1;

      offset = self->pQmfDomain->globalConf.nBandsSynthesis *
               self->pQmfDomain->globalConf.nQmfTimeSlots;

      {
        const PCM_MPS *inSamples =
            &inData[ts * self->pQmfDomain->globalConf.nBandsAnalysis];

        CalculateSpaceAnalysisQmf(
            &self->pQmfDomain->QmfDomainIn[ch].fb, inSamples + (ch * offset),
            self->qmfResidualReal__FDK[0][0], self->qmfResidualImag__FDK[0][0]);

        if (!isTwoChMode(self->upmixType) && !bypassMode) {
          int i;
          FIXP_DBL *RESTRICT self_qmfResidualReal__FDK_0_0 =
              &self->qmfResidualReal__FDK[0][0][0];
          FIXP_DBL *RESTRICT self_qmfResidualImag__FDK_0_0 =
              &self->qmfResidualImag__FDK[0][0][0];

          if ((self->pQmfDomain->globalConf.nBandsAnalysis == 24) &&
              !(self->stereoConfigIndex == 3)) {
            for (i = 0; i < self->qmfBands; i++) {
              self_qmfResidualReal__FDK_0_0[i] =
                  fMult(scaleValueSaturate(self_qmfResidualReal__FDK_0_0[i],
                                           1 + self->sacInDataHeadroom - (1)),
                        self->clipProtectGain__FDK);
              self_qmfResidualImag__FDK_0_0[i] =
                  fMult(scaleValueSaturate(self_qmfResidualImag__FDK_0_0[i],
                                           1 + self->sacInDataHeadroom - (1)),
                        self->clipProtectGain__FDK);
            }
          } else {
            for (i = 0; i < self->qmfBands; i++) {
              self_qmfResidualReal__FDK_0_0[i] =
                  fMult(scaleValueSaturate(self_qmfResidualReal__FDK_0_0[i],
                                           self->sacInDataHeadroom - (1)),
                        self->clipProtectGain__FDK);
              self_qmfResidualImag__FDK_0_0[i] =
                  fMult(scaleValueSaturate(self_qmfResidualImag__FDK_0_0[i],
                                           self->sacInDataHeadroom - (1)),
                        self->clipProtectGain__FDK);
            }
          }
        }
      }
    }

    SpatialDecHybridAnalysis(
        self,
        self->qmfInputReal__FDK, self->qmfInputImag__FDK,
        self->hybInputReal__FDK, self->hybInputImag__FDK, ts, numInputChannels);

    if (bypassMode) {
      SpatialDecApplyBypass(
          self, self->hybInputReal__FDK,
          self->hybInputImag__FDK,
          self->hybOutputRealDry__FDK,
          self->hybOutputImagDry__FDK, numInputChannels);
    } else
    {
      FIXP_DBL *pxReal[MAX_NUM_XCHANNELS] = {NULL};
      FIXP_DBL *pxImag[MAX_NUM_XCHANNELS] = {NULL};

      SpatialDecCreateX(self,
                        self->hybInputReal__FDK,
                        self->hybInputImag__FDK, pxReal, pxImag);

      {
        SpatialDecApplyM1_CreateW_Mode212(
            self, frame, pxReal, pxImag,
            self->wReal__FDK,
            self->wImag__FDK);
      }
      if (err != MPS_OK) goto bail;

      int applyM2Config = APPLY_M2_NONE;

      applyM2Config = APPLY_M2;
      if ((self->pConfigCurrent->syntaxFlags &
           (SACDEC_SYNTAX_USAC | SACDEC_SYNTAX_RSVD50)) &&
          (self->tempShapeConfig != 1) && (self->tempShapeConfig != 2)) {
        if (self->phaseCoding == 3)
          applyM2Config = APPLY_M2_MODE212_Res_PhaseCoding;
        else
          applyM2Config = APPLY_M2_MODE212;
      }

      switch (applyM2Config) {
        case APPLY_M2_MODE212: {
          err = SpatialDecApplyM2_Mode212(
              self, ps, alpha, self->wReal__FDK, self->wImag__FDK,
              self->hybOutputRealDry__FDK, self->hybOutputImagDry__FDK);
        } break;
        case APPLY_M2_MODE212_Res_PhaseCoding:
          err = SpatialDecApplyM2_Mode212_ResidualsPlusPhaseCoding(
              self, ps, alpha, self->wReal__FDK, self->wImag__FDK,
              self->hybOutputRealDry__FDK, self->hybOutputImagDry__FDK);
          break;
        case APPLY_M2:
          err = SpatialDecApplyM2(
              self, ps, alpha, self->wReal__FDK, self->wImag__FDK,
              self->hybOutputRealDry__FDK, self->hybOutputImagDry__FDK,
              self->hybOutputRealWet__FDK, self->hybOutputImagWet__FDK);
          break;
        default:
          err = MPS_APPLY_M2_ERROR;
          goto bail;
      }

      if (err != MPS_OK) goto bail;

      if ((self->tempShapeConfig == 2) && (!isTwoChMode(self->upmixType))) {
        SpatialDecReshapeBBEnv(self, frame,
                               ts);
      }

      if ((self->tempShapeConfig == 1) && (!isTwoChMode(self->upmixType))) {
        for (ch = 0; ch < self->numOutputChannels; ch++) {
          for (hyb = 0; hyb < self->tp_hybBandBorder; hyb++) {
            self->hybOutputRealDry__FDK[ch][hyb] =
                fAddSaturate(self->hybOutputRealDry__FDK[ch][hyb],
                             self->hybOutputRealWet__FDK[ch][hyb]);
            self->hybOutputImagDry__FDK[ch][hyb] =
                fAddSaturate(self->hybOutputImagDry__FDK[ch][hyb],
                             self->hybOutputImagWet__FDK[ch][hyb]);
          }
        }
        err = subbandTPApply(
            self, frame);
        if (err != MPS_OK) goto bail;
      }
      else {
        if ((self->tempShapeConfig == 1) || (self->tempShapeConfig == 2)) {
          int nHybBands;
          nHybBands = self->hybridBands;

          for (ch = 0; ch < self->numOutputChannels; ch++) {
            FIXP_DBL *RESTRICT pRealDry = self->hybOutputRealDry__FDK[ch];
            FIXP_DBL *RESTRICT pImagDry = self->hybOutputImagDry__FDK[ch];
            FIXP_DBL *RESTRICT pRealWet = self->hybOutputRealWet__FDK[ch];
            FIXP_DBL *RESTRICT pImagWet = self->hybOutputImagWet__FDK[ch];
            for (hyb = 0; hyb < nHybBands; hyb++) {
              pRealDry[hyb] = fAddSaturate(pRealDry[hyb], pRealWet[hyb]);
              pImagDry[hyb] = fAddSaturate(pImagDry[hyb], pImagWet[hyb]);
            }
            for (; hyb < self->hybridBands; hyb++) {
              pRealDry[hyb] = fAddSaturate(pRealDry[hyb], pRealWet[hyb]);
            }
          }
        }
      }
    }

    if ((self->phaseCoding == 1) && (bypassMode == 0)) {
      SpatialDecApplyPhase(
          self, alpha, (ts == currSlot)
      );
    }

    err = SpatialDecSynthesis(
        self, ts_io,
        self->hybOutputRealDry__FDK,
        self->hybOutputImagDry__FDK, self->timeOut__FDK,
        numInputChannels, mapDescr);

    if (err != MPS_OK) goto bail;

    if (ts == currSlot) {
      SpatialDecBufferMatrices(self);

      prevSlot = currSlot;
      ps++;
    }
  }

  self->prevTimeSlot = prevSlot;
  self->curTimeSlot = ts;
  self->curPs = ps;

bail:

  return err;
}

SACDEC_ERROR SpatialDecApplyFrame(
    spatialDec *self,
    SPATIAL_BS_FRAME *frame,
    SPATIALDEC_INPUT_MODE inputMode, PCM_MPS *inData,
    FIXP_DBL **qmfInDataReal,
    FIXP_DBL **qmfInDataImag,
    PCM_MPS *pcmOutBuf,
    UINT nSamples, UINT *pControlFlags, int numInputChannels,
    const FDK_channelMapDescr *const mapDescr) {
  SACDEC_ERROR err = MPS_OK;

  int fDecAndMapFrameData;
  int controlFlags;

  FDK_ASSERT(self != NULL);
  FDK_ASSERT(pControlFlags != NULL);
  FDK_ASSERT(pcmOutBuf != NULL);
  FDK_ASSERT(self->sacInDataHeadroom >= (1));

  self->errInt = err;

  controlFlags = *pControlFlags;

  if ((self->pConfigCurrent->syntaxFlags & SACDEC_SYNTAX_USAC) &&
      (self->stereoConfigIndex > 1)) {
    numInputChannels =
        1;
  }

  if (numInputChannels != self->numInputChannels) {
    controlFlags |= MPEGS_CONCEAL;
    if (numInputChannels > self->createParams.maxNumInputChannels) {
      return MPS_INVALID_PARAMETER;
    }
  }

  self->timeOut__FDK = pcmOutBuf;

  fDecAndMapFrameData = frame->newBsData;

  if (((fDecAndMapFrameData ==
        0)
       && (self->curTimeSlot + (int)nSamples / self->qmfBands >
           self->timeSlots)) ||
      (frame->numParameterSets ==
       0)) {
    fDecAndMapFrameData = 1;
    controlFlags |= MPEGS_CONCEAL;
  }

  if ((fDecAndMapFrameData == 0) &&
      (frame->paramSlot[fMax(0, frame->numParameterSets - 1)] !=
           (self->timeSlots - 1) ||
       self->curTimeSlot >
           frame->paramSlot[self->curPs])) {
    fDecAndMapFrameData = 1;
    controlFlags |= MPEGS_CONCEAL;
  }

  SpatialDecConcealment_UpdateState(
      &self->concealInfo,
      (controlFlags & MPEGS_CONCEAL)
          ? 0
          : 1);

  if (fDecAndMapFrameData) {
    frame->newBsData = 0;
    self->prevTimeSlot = -1;
    self->curTimeSlot = 0;
    self->curPs = 0;

    if (controlFlags & MPEGS_CONCEAL) {
      SpatialDecClearFrameData(self, frame, &self->createParams);
    }

    {
      err = SpatialDecDecodeFrame(self, frame);
    }

    if (err != MPS_OK) {
      controlFlags |= MPEGS_CONCEAL;
      if (self->errInt == MPS_OK) {
        self->errInt = err;
      }
    }
  }

  err = SpatialDecApplyParameterSets(
      self, frame, inputMode, inData, qmfInDataReal, qmfInDataImag, nSamples,
      controlFlags | ((err == MPS_OK) ? 0 : MPEGS_BYPASSMODE), numInputChannels,
      mapDescr);
  if (err != MPS_OK) {
    goto bail;
  }

bail:

  *pControlFlags = controlFlags;

  return err;
}
