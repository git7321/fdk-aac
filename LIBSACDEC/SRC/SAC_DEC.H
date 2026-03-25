/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Decoder Library structures

*******************************************************************************/

#ifndef SAC_DEC_H
#define SAC_DEC_H

#include "common_fix.h"
#include "sac_dec_interface.h"
#include "FDK_qmf_domain.h"
#include "sac_qmf.h"
#include "FDK_bitstream.h"
#include "sac_calcM1andM2.h"
#include "FDK_hybrid.h"
#include "FDK_decorrelate.h"
#include "sac_reshapeBBEnv.h"
#include "sac_dec_conceal.h"
#include "sac_tsd.h"

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

#define ICCdefault 0
#define IPDdefault 0
#define arbdmxGainDefault 0
#define CPCdefault 10
#define tttCLD1default 15
#define tttCLD2default 0

#define IS_HQ_ONLY(aot)                                                      \
  ((aot) == AOT_ER_AAC_LD || (aot) == AOT_ER_AAC_ELD || (aot) == AOT_USAC || \
   (aot) == AOT_RSVD50)

#define SCONST(x) FL2FXCONST_DBL(x)

#define PC_NUM_BANDS (8)
#define PC_NUM_HYB_BANDS (PC_NUM_BANDS - 3 + 10)
#define ABS_THR (1e-9f * 32768 * 32768)

#define MAX_HYBRID_BANDS (MAX_NUM_QMF_BANDS - 3 + 10)
#define HYBRID_FILTER_DELAY (6)

#define MAX_RESIDUAL_FRAMES (4)
#define MAX_RESIDUAL_BISTREAM \
  (836)
#define MAX_MDCT_COEFFS (1024)
#define SACDEC_RESIDUAL_BS_BUF_SIZE \
  (1024)

#define MAX_NUM_PARAMS (MAX_NUM_OTT + 4 * MAX_NUM_TTT + MAX_INPUT_CHANNELS)
#define MAX_NUM_PARAMETERS (MAX(MAX_NUM_PARAMS, MAX_NUM_OTT))

#define MAX_PARAMETER_SETS (9)
#define MAX_M2_INPUT (MAX_OUTPUT_CHANNELS)
#define MAX_QMF_BANDS_TO_HYBRID \
  (3)
#define PROTO_LEN (13)
#define BUFFER_LEN_LF (PROTO_LEN)
#define BUFFER_LEN_HF ((PROTO_LEN - 1) / 2)

#define MAX_NO_DECORR_CHANNELS (MAX_OUTPUT_CHANNELS)
#define HRTF_AZIMUTHS (5)

#define MAX_NUM_OTT_AT 0

typedef enum {
  UPMIXTYPE_BYPASS = -1,
  UPMIXTYPE_NORMAL = 0
} UPMIXTYPE;

static inline int isTwoChMode(UPMIXTYPE upmixType) {
  int retval = 0;
  return retval;
}

#define MPEGS_BYPASSMODE (0x00000001)
#define MPEGS_CONCEAL (0x00000002)

typedef struct STP_DEC *HANDLE_STP_DEC;

typedef struct {
  SCHAR bsQuantCoarseXXXprev;
  SCHAR bsQuantCoarseXXXprevParse;
} LOSSLESSSTATE;

typedef struct {
  SCHAR bsXXXDataMode[MAX_PARAMETER_SETS];
  SCHAR bsQuantCoarseXXX[MAX_PARAMETER_SETS];
  SCHAR bsFreqResStrideXXX[MAX_PARAMETER_SETS];
  SCHAR nocmpQuantCoarseXXX[MAX_PARAMETER_SETS];
  LOSSLESSSTATE *state;
} LOSSLESSDATA;

struct SPATIAL_BS_FRAME_struct {
  UCHAR bsIndependencyFlag;
  UCHAR newBsData;
  UCHAR numParameterSets;

  INT paramSlot[MAX_PARAMETER_SETS];

  SCHAR cmpOttCLDidx[MAX_NUM_OTT + MAX_NUM_OTT_AT][MAX_PARAMETER_SETS]
                    [MAX_PARAMETER_BANDS];
  SCHAR cmpOttICCidx[MAX_NUM_OTT][MAX_PARAMETER_SETS][MAX_PARAMETER_BANDS];

  UCHAR bsSmoothMode[MAX_PARAMETER_SETS];
  UCHAR bsSmoothTime[MAX_PARAMETER_SETS];
  UCHAR bsFreqResStrideSmg[MAX_PARAMETER_SETS];
  UCHAR bsSmgData[MAX_PARAMETER_SETS]
                 [MAX_PARAMETER_BANDS];

  SCHAR (*cmpArbdmxGainIdx)[MAX_PARAMETER_SETS][MAX_PARAMETER_BANDS];

  LOSSLESSDATA *CLDLosslessData;
  LOSSLESSDATA *ICCLosslessData;

  LOSSLESSDATA *IPDLosslessData;
  SCHAR (*cmpOttIPDidx)[MAX_PARAMETER_SETS][MAX_PARAMETER_BANDS];
  int phaseMode;
  int OpdSmoothingMode;

  UCHAR tempShapeEnableChannelGES[MAX_OUTPUT_CHANNELS];
  UCHAR bsEnvShapeData[MAX_OUTPUT_CHANNELS]
                      [MAX_TIME_SLOTS];

  UCHAR tempShapeEnableChannelSTP[MAX_OUTPUT_CHANNELS];

  TSD_DATA TsdData[1];
};

typedef struct {
  LOSSLESSSTATE CLDLosslessState[MAX_NUM_PARAMETERS];
  LOSSLESSSTATE ICCLosslessState[MAX_NUM_PARAMETERS];
  LOSSLESSSTATE IPDLosslessState[MAX_NUM_PARAMETERS];
} BS_LL_STATE;

typedef struct {
  int prevParamSlot;
  int prevSmgTime;
  UCHAR prevSmgData[MAX_PARAMETER_BANDS];

  FIXP_DBL opdLeftState__FDK[MAX_PARAMETER_BANDS];
  FIXP_DBL opdRightState__FDK[MAX_PARAMETER_BANDS];

} SMOOTHING_STATE;

typedef struct {
  FIXP_DBL alpha__FDK;
  FIXP_DBL beta__FDK;
  FIXP_DBL partNrgPrev__FDK[2 * MAX_OUTPUT_CHANNELS + MAX_INPUT_CHANNELS]
                           [BB_ENV_SIZE];
  FIXP_DBL normNrgPrev__FDK[2 * MAX_OUTPUT_CHANNELS + MAX_INPUT_CHANNELS];
  FIXP_DBL frameNrgPrev__FDK[2 * MAX_OUTPUT_CHANNELS + MAX_INPUT_CHANNELS];
  INT partNrgPrevSF[2 * MAX_OUTPUT_CHANNELS + MAX_INPUT_CHANNELS];
  INT partNrgPrev2SF[2 * MAX_OUTPUT_CHANNELS + MAX_INPUT_CHANNELS];
  INT normNrgPrevSF[2 * MAX_OUTPUT_CHANNELS + MAX_INPUT_CHANNELS];
  INT frameNrgPrevSF[2 * MAX_OUTPUT_CHANNELS + MAX_INPUT_CHANNELS];
} RESHAPE_BBENV_STATE;

typedef struct {
  int maxNumInputChannels;
  int maxNumOutputChannels;
  int maxNumQmfBands;
  int maxNumHybridBands;
  int maxNumXChannels;
  int maxNumVChannels;
  int maxNumDecorChannels;
  int maxNumCmplxQmfBands;
  int maxNumCmplxHybBands;
  int maxNumResChannels;
  int bProcResidual;
  int maxNumResidualChannels;
  int maxNumOttBoxes;
  int maxNumParams;

} SACDEC_CREATION_PARAMS;

struct spatialDec_struct {
  SACDEC_ERROR
  errInt;
  int staticDecScale;

  int samplingFreq;
  CFG_LEVEL decoderLevel;
  CFG_EXTENT decoderMode;
  CFG_BINAURAL binauralMode;

  SACDEC_CREATION_PARAMS createParams;

  int numComplexProcessingBands;

  int treeConfig;

  int numInputChannels;
  int numOutputChannels;
  int numOttBoxes;
  int numM2rows;

  int numOutputChannelsAT;

  int quantMode;
  int arbitraryDownmix;
  int residualCoding;
  UCHAR nrResidualFrame;
  UCHAR nrArbDownmixResidualFrame;
  FDK_BITSTREAM **hResidualBitstreams;
  int tempShapeConfig;
  int decorrType;
  int decorrConfig;
  int envQuantMode;

  FIXP_DBL clipProtectGain__FDK;
  char clipProtectGainSF__FDK;

  int numDirektSignals;
  int wStartResidualIdx;
  int numDecorSignals;
  int numVChannels;
  int numXChannels;

  int timeSlots;
  int curTimeSlot;
  int prevTimeSlot;
  int curPs;
  int frameLength;
  UPMIXTYPE upmixType;
  int partiallyComplex;
  int useFDreverb;

  int bShareDelayWithSBR;

  int tp_hybBandBorder;

  int qmfBands;
  int hybridBands;
  const SCHAR *kernels;

  int TsdTs;

  int *param2hyb;
  int kernels_width[MAX_PARAMETER_BANDS];

  int residualSamplingFreq;
  UCHAR residualPresent[MAX_NUM_OTT + MAX_NUM_TTT];
  UCHAR residualBands[MAX_NUM_OTT + MAX_NUM_TTT];
  UCHAR residualQMFBands[MAX_NUM_OTT + MAX_NUM_TTT];
  SPATIAL_SPECIFIC_CONFIG *pConfigCurrent;

  int arbdmxFramesPerSpatialFrame;
  int arbdmxUpdQMF;

  int numParameterBands;
  int bitstreamParameterBands;
  int *numOttBands;

  UCHAR extendFrame;
  UCHAR numParameterSetsPrev;

  int *smgTime;
  UCHAR **smgData;

  SCHAR **ottCLDidxPrev;
  SCHAR **ottICCidxPrev;
  SCHAR **arbdmxGainIdxPrev;
  SCHAR **ottIPDidxPrev;
  SCHAR ***outIdxData;

  SCHAR **cmpOttCLDidxPrev;
  SCHAR **cmpOttICCidxPrev;
  SCHAR ***ottICCdiffidx;
  SCHAR **cmpOttIPDidxPrev;

  SCHAR **cmpArbdmxGainIdxPrev;

  SCHAR ***ottCLD__FDK;
  SCHAR ***ottICC__FDK;

  SCHAR ***arbdmxGain__FDK;

  FIXP_DBL *arbdmxAlpha__FDK;
  FIXP_DBL *arbdmxAlphaPrev__FDK;

  UCHAR stereoConfigIndex;
  int highRateMode;

  int phaseCoding;

  SCHAR ***ottIPD__FDK;

  FIXP_DBL PhaseLeft__FDK[MAX_PARAMETER_BANDS];
  FIXP_DBL PhaseRight__FDK[MAX_PARAMETER_BANDS];
  FIXP_DBL PhasePrevLeft__FDK[MAX_PARAMETER_BANDS];
  FIXP_DBL PhasePrevRight__FDK[MAX_PARAMETER_BANDS];
  int numOttBandsIPD;

  FIXP_DBL ***M2Real__FDK;
  FIXP_DBL ***M2Imag__FDK;
  FIXP_DBL ***M2RealPrev__FDK;
  FIXP_DBL ***M2ImagPrev__FDK;

  FIXP_DBL ***qmfInputRealDelayBuffer__FDK;
  FIXP_DBL ***qmfInputImagDelayBuffer__FDK;

  int pc_filterdelay;
  int qmfInputDelayBufPos;
  FIXP_DBL **qmfInputReal__FDK;
  FIXP_DBL **qmfInputImag__FDK;

  FIXP_DBL **hybInputReal__FDK;
  FIXP_DBL **hybInputImag__FDK;

  FIXP_DBL **binInputReverb;

  FIXP_DBL binGain, reverbGain;
  FIXP_DBL binCenterGain, reverbCenterGain;

  FIXP_DBL ***qmfResidualReal__FDK;
  FIXP_DBL ***qmfResidualImag__FDK;

  FIXP_DBL **hybResidualReal__FDK;
  FIXP_DBL **hybResidualImag__FDK;

  int qmfOutputRealDryDelayBufPos;
  FIXP_DBL ***qmfOutputRealDryDelayBuffer__FDK;
  FIXP_DBL ***qmfOutputImagDryFilterBuffer__FDK;
  FIXP_DBL *qmfOutputImagDryFilterBufferBase__FDK;

  FIXP_DBL **wReal__FDK;
  FIXP_DBL **wImag__FDK;

  FIXP_DBL **hybOutputRealDry__FDK;
  FIXP_DBL **hybOutputImagDry__FDK;
  FIXP_DBL **hybOutputRealWet__FDK;
  FIXP_DBL **hybOutputImagWet__FDK;
  PCM_MPS *timeOut__FDK;

  HANDLE_FDK_QMF_DOMAIN pQmfDomain;

  FDK_ANA_HYB_FILTER
  *hybridAnalysis;
  FDK_SYN_HYB_FILTER
  *hybridSynthesis;
  FIXP_DBL **
      pHybridAnaStatesLFdmx;
  FIXP_DBL **
      pHybridAnaStatesHFdmx;
  FIXP_DBL **
      pHybridAnaStatesLFres;
  FIXP_DBL **
      pHybridAnaStatesHFres;

  DECORR_DEC *apDecor;
  FIXP_DBL **pDecorBufferCplx;

  SMOOTHING_STATE *smoothState;

  RESHAPE_BBENV_STATE *reshapeBBEnvState;
  SCHAR row2channelDmxGES[MAX_OUTPUT_CHANNELS];

  HANDLE_STP_DEC hStpDec;

  const UCHAR *pActivM2ParamBands;

  int bOverwriteM1M2prev;
  SpatialDecConcealmentInfo concealInfo;

  INT sacInDataHeadroom;
};

#define SACDEC_SYNTAX_MPS 1
#define SACDEC_SYNTAX_USAC 2
#define SACDEC_SYNTAX_RSVD50 4
#define SACDEC_SYNTAX_L2 8
#define SACDEC_SYNTAX_L3 16
#define SACDEC_SYNTAX_LD 32

static inline int GetProcBand(spatialDec_struct *self, int qs) {
  return self->kernels[qs];
}

#endif /* SAC_DEC_H */
