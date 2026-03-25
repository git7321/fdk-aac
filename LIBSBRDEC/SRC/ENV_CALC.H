/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef ENV_CALC_H
#define ENV_CALC_H

#include "sbrdecoder.h"
#include "env_extr.h"

typedef struct {
  FIXP_DBL filtBuffer[MAX_FREQ_COEFFS];
  FIXP_DBL filtBufferNoise[MAX_FREQ_COEFFS];
  SCHAR filtBuffer_e[MAX_FREQ_COEFFS];
  SCHAR filtBufferNoise_e;

  int startUp;
  int phaseIndex;
  int prevTranEnv;

  ULONG harmFlagsPrev[ADD_HARMONICS_FLAGS_SIZE];
  UCHAR harmIndex;
  int sbrPatchingMode;

  FIXP_SGL prevSbrNoiseFloorLevel[MAX_NOISE_COEFFS];
  UCHAR prevNNfb;
  UCHAR prevNSfb[2];
  UCHAR prevLoSubband;
  UCHAR prevHiSubband;
  UCHAR prev_ov_highSubband;
  UCHAR *prevFreqBandTable[2];
  UCHAR prevFreqBandTableLo[MAX_FREQ_COEFFS / 2 + 1];
  UCHAR prevFreqBandTableHi[MAX_FREQ_COEFFS + 1];
  UCHAR prevFreqBandTableNoise[MAX_NOISE_COEFFS + 1];
  SCHAR sinusoidal_positionPrev;
  ULONG harmFlagsPrevActive[ADD_HARMONICS_FLAGS_SIZE];
} SBR_CALCULATE_ENVELOPE;

typedef SBR_CALCULATE_ENVELOPE *HANDLE_SBR_CALCULATE_ENVELOPE;

void calculateSbrEnvelope(
    QMF_SCALE_FACTOR *sbrScaleFactor,
    HANDLE_SBR_CALCULATE_ENVELOPE h_sbr_cal_env,
    HANDLE_SBR_HEADER_DATA hHeaderData, HANDLE_SBR_FRAME_DATA hFrameData,
    PVC_DYNAMIC_DATA *pPvcDynamicData, FIXP_DBL **analysBufferReal,
    FIXP_DBL *
        *analysBufferImag,
    const int useLP,
    FIXP_DBL *degreeAlias,
    const UINT flags, const int frameErrorFlag);

SBR_ERROR
createSbrEnvelopeCalc(HANDLE_SBR_CALCULATE_ENVELOPE hSbrCalculateEnvelope,
                      HANDLE_SBR_HEADER_DATA hHeaderData, const int chan,
                      const UINT flags);

int deleteSbrEnvelopeCalc(HANDLE_SBR_CALCULATE_ENVELOPE hSbrCalculateEnvelope);

void resetSbrEnvelopeCalc(HANDLE_SBR_CALCULATE_ENVELOPE hCalEnv);

SBR_ERROR
ResetLimiterBands(UCHAR *limiterBandTable, UCHAR *noLimiterBands,
                  UCHAR *freqBandTable, int noFreqBands,
                  const PATCH_PARAM *patchParam, int noPatches,
                  int limiterBands, UCHAR sbrPatchingMode,
                  int xOverQmf[MAX_NUM_PATCHES], int sbrRatio);

void rescaleSubbandSamples(FIXP_DBL **re, FIXP_DBL **im, int lowSubband,
                           int noSubbands, int start_pos, int next_pos,
                           int shift);

FIXP_DBL maxSubbandSample(FIXP_DBL **analysBufferReal_m,
                          FIXP_DBL **analysBufferImag_m, int lowSubband,
                          int highSubband, int start_pos, int stop_pos);

#endif  // ENV_CALC_H
