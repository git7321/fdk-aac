/**************************** SBR encoder library ******************************

   Author(s):   Markus Multrus

   Description: PS Wrapper, Downmix header file

*******************************************************************************/

#ifndef PS_MAIN_H
#define PS_MAIN_H

/* Includes ******************************************************************/

#include "sbr_def.h"
#include "qmf.h"
#include "ps_encode.h"
#include "FDK_bitstream.h"
#include "FDK_hybrid.h"

/* Data Types ****************************************************************/
typedef enum {
  PSENC_STEREO_BANDS_INVALID = 0,
  PSENC_STEREO_BANDS_10 = 10,
  PSENC_STEREO_BANDS_20 = 20

} PSENC_STEREO_BANDS_CONFIG;

typedef enum {
  PSENC_NENV_1 = 1,
  PSENC_NENV_2 = 2,
  PSENC_NENV_4 = 4,
  PSENC_NENV_DEFAULT = PSENC_NENV_2,
  PSENC_NENV_MAX = PSENC_NENV_4

} PSENC_NENV_CONFIG;

typedef struct {
  UINT bitrateFrom; /* inclusive */
  UINT bitrateTo;   /* exclusive */
  PSENC_STEREO_BANDS_CONFIG nStereoBands;
  PSENC_NENV_CONFIG nEnvelopes;
  LONG iidQuantErrorThreshold;
} psTuningTable_t;

/* Function / Class Declarations *********************************************/

typedef struct T_PARAMETRIC_STEREO {
  HANDLE_PS_ENCODE hPsEncode;
  PS_OUT psOut[2];

  FIXP_DBL __staticHybridData[HYBRID_READ_OFFSET][MAX_PS_CHANNELS][2]
                             [MAX_HYBRID_BANDS];
  FIXP_DBL
  *pHybridData[HYBRID_READ_OFFSET + HYBRID_FRAMESIZE][MAX_PS_CHANNELS][2];

  FIXP_DBL qmfDelayLines[2][32 >> 1][64];
  int qmfDelayScale;

  INT psDelay;
  UINT maxEnvelopes;
  UCHAR dynBandScale[PS_MAX_BANDS];
  FIXP_DBL maxBandValue[PS_MAX_BANDS];
  SCHAR dmxScale;
  INT initPS;
  INT noQmfSlots;
  INT noQmfBands;

  FIXP_DBL __staticHybAnaStatesLF[MAX_PS_CHANNELS][2 * HYBRID_FILTER_LENGTH *
                                                   HYBRID_MAX_QMF_BANDS];
  FIXP_DBL __staticHybAnaStatesHF[MAX_PS_CHANNELS][2 * HYBRID_FILTER_DELAY *
                                                   (64 - HYBRID_MAX_QMF_BANDS)];
  FDK_ANA_HYB_FILTER fdkHybAnaFilter[MAX_PS_CHANNELS];
  FDK_SYN_HYB_FILTER fdkHybSynFilter;

} PARAMETRIC_STEREO;

typedef struct T_PSENC_CONFIG {
  INT frameSize;
  INT qmfFilterMode;
  INT sbrPsDelay;
  PSENC_STEREO_BANDS_CONFIG nStereoBands;
  PSENC_NENV_CONFIG maxEnvelopes;
  FIXP_DBL iidQuantErrorThreshold;

} PSENC_CONFIG, *HANDLE_PSENC_CONFIG;

typedef struct T_PARAMETRIC_STEREO *HANDLE_PARAMETRIC_STEREO;

FDK_PSENC_ERROR PSEnc_Create(HANDLE_PARAMETRIC_STEREO *phParametricStereo);

FDK_PSENC_ERROR PSEnc_Init(HANDLE_PARAMETRIC_STEREO hParametricStereo,
                           const HANDLE_PSENC_CONFIG hPsEncConfig,
                           INT noQmfSlots, INT noQmfBands, UCHAR *dynamic_RAM);

FDK_PSENC_ERROR PSEnc_Destroy(HANDLE_PARAMETRIC_STEREO *phParametricStereo);

FDK_PSENC_ERROR FDKsbrEnc_PSEnc_ParametricStereoProcessing(
    HANDLE_PARAMETRIC_STEREO hParametricStereo, INT_PCM *samples[2],
    UINT timeInStride, QMF_FILTER_BANK **hQmfAnalysis,
    FIXP_DBL **RESTRICT downmixedRealQmfData,
    FIXP_DBL **RESTRICT downmixedImagQmfData, INT_PCM *downsampledOutSignal,
    HANDLE_QMF_FILTER_BANK sbrSynthQmf, SCHAR *qmfScale, const int sendHeader);

INT FDKsbrEnc_PSEnc_WritePSData(HANDLE_PARAMETRIC_STEREO hParametricStereo,
                                HANDLE_FDK_BITSTREAM hBitstream);

#endif /* PS_MAIN_H */
