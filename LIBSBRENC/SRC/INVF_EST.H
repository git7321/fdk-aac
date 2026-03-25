/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef INVF_EST_H
#define INVF_EST_H

#include "sbr_encoder.h"
#include "sbr_def.h"

#define INVF_SMOOTHING_LENGTH 2

typedef struct {
  const FIXP_DBL *quantStepsSbr;
  const FIXP_DBL *quantStepsOrig;
  const FIXP_DBL *nrgBorders;
  INT numRegionsSbr;
  INT numRegionsOrig;
  INT numRegionsNrg;
  INVF_MODE regionSpace[5][5];
  INVF_MODE regionSpaceTransient[5][5];
  INT EnergyCompFactor[5];

} DETECTOR_PARAMETERS;

typedef struct {
  FIXP_DBL origQuotaMean[INVF_SMOOTHING_LENGTH + 1];
  FIXP_DBL sbrQuotaMean[INVF_SMOOTHING_LENGTH + 1];
  FIXP_DBL origQuotaMeanStrongest[INVF_SMOOTHING_LENGTH + 1];
  FIXP_DBL sbrQuotaMeanStrongest[INVF_SMOOTHING_LENGTH + 1];

  FIXP_DBL origQuotaMeanFilt;
  FIXP_DBL sbrQuotaMeanFilt;
  FIXP_DBL origQuotaMeanStrongestFilt;
  FIXP_DBL sbrQuotaMeanStrongestFilt;

  FIXP_DBL origQuotaMax;
  FIXP_DBL sbrQuotaMax;

  FIXP_DBL avgNrg;
} DETECTOR_VALUES;

typedef struct {
  INT numberOfStrongest;

  INT prevRegionSbr[MAX_NUM_NOISE_VALUES];
  INT prevRegionOrig[MAX_NUM_NOISE_VALUES];

  INT freqBandTableInvFilt[MAX_NUM_NOISE_VALUES];
  INT noDetectorBands;
  INT noDetectorBandsMax;

  const DETECTOR_PARAMETERS *detectorParams;

  INVF_MODE prevInvfMode[MAX_NUM_NOISE_VALUES];
  DETECTOR_VALUES detectorValues[MAX_NUM_NOISE_VALUES];

  FIXP_DBL nrgAvg;
  FIXP_DBL wmQmf[MAX_NUM_NOISE_VALUES];
} SBR_INV_FILT_EST;

typedef SBR_INV_FILT_EST *HANDLE_SBR_INV_FILT_EST;

void FDKsbrEnc_qmfInverseFilteringDetector(HANDLE_SBR_INV_FILT_EST hInvFilt,
                                           FIXP_DBL **quotaMatrix,
                                           FIXP_DBL *nrgVector,
                                           SCHAR *indexVector, INT startIndex,
                                           INT stopIndex, INT transientFlag,
                                           INVF_MODE *infVec);

INT FDKsbrEnc_initInvFiltDetector(HANDLE_SBR_INV_FILT_EST hInvFilt,
                                  INT *freqBandTableDetector,
                                  INT numDetectorBands, UINT useSpeechConfig);

INT FDKsbrEnc_resetInvFiltDetector(HANDLE_SBR_INV_FILT_EST hInvFilt,
                                   INT *freqBandTableDetector,
                                   INT numDetectorBands);

#endif /* _QMF_INV_FILT_H */
