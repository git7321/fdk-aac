/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "invf_est.h"
#include "sbr_misc.h"

#include "genericStds.h"

#define MAX_NUM_REGIONS 10
#define SCALE_FAC_QUO 512.0f
#define SCALE_FAC_NRG 256.0f

#ifndef min
#define min(a, b) (a < b ? a : b)
#endif

#ifndef max
#define max(a, b) (a > b ? a : b)
#endif

static const FIXP_DBL quantStepsSbr[4] = {
    0x00400000, 0x02800000, 0x03800000,
    0x04c00000};
static const FIXP_DBL quantStepsOrig[4] = {
    0x00000000, 0x00c00000, 0x01c00000,
    0x02800000};
static const FIXP_DBL nrgBorders[4] = {
    0x0c800000, 0x0f000000, 0x11800000,
    0x14000000};

static const DETECTOR_PARAMETERS detectorParamsAAC = {
    quantStepsSbr,
    quantStepsOrig,
    nrgBorders,
    4,
    4,
    4,
    {
        {INVF_MID_LEVEL, INVF_LOW_LEVEL, INVF_OFF, INVF_OFF,
         INVF_OFF},
        {INVF_MID_LEVEL, INVF_LOW_LEVEL, INVF_OFF, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_LOW_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_OFF,
         INVF_OFF}
    },
    {
        {INVF_LOW_LEVEL, INVF_LOW_LEVEL, INVF_LOW_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_LOW_LEVEL, INVF_LOW_LEVEL, INVF_LOW_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_MID_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_OFF,
         INVF_OFF}
    },
    {-4, -3, -2, -1,
     0}
};

static const FIXP_DBL hysteresis =
    0x00400000;

static const DETECTOR_PARAMETERS detectorParamsAACSpeech = {
    quantStepsSbr,
    quantStepsOrig,
    nrgBorders,
    4,
    4,
    4,
    {
        {INVF_MID_LEVEL, INVF_MID_LEVEL, INVF_LOW_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_MID_LEVEL, INVF_MID_LEVEL, INVF_LOW_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_MID_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_OFF,
         INVF_OFF}
    },
    {
        {INVF_MID_LEVEL, INVF_MID_LEVEL, INVF_LOW_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_MID_LEVEL, INVF_MID_LEVEL, INVF_LOW_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_MID_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_OFF,
         INVF_OFF},
        {INVF_HIGH_LEVEL, INVF_HIGH_LEVEL, INVF_MID_LEVEL, INVF_OFF,
         INVF_OFF}
    },
    {-4, -3, -2, -1,
     0}
};

typedef const FIXP_DBL FIR_FILTER[5];

static const FIR_FILTER fir_0 = {0x7fffffff, 0x00000000, 0x00000000, 0x00000000,
                                 0x00000000};
static const FIR_FILTER fir_1 = {0x2aaaaa80, 0x555554ff, 0x00000000, 0x00000000,
                                 0x00000000};
static const FIR_FILTER fir_2 = {0x10000000, 0x30000000, 0x40000000, 0x00000000,
                                 0x00000000};
static const FIR_FILTER fir_3 = {0x077f80e8, 0x199999a0, 0x2bb3b240, 0x33333340,
                                 0x00000000};
static const FIR_FILTER fir_4 = {0x04130598, 0x0ebdb000, 0x1becfa60, 0x2697a4c0,
                                 0x2aaaaa80};

static const FIR_FILTER *const fir_table[5] = {&fir_0, &fir_1, &fir_2, &fir_3,
                                               &fir_4};

static void calculateDetectorValues(
    FIXP_DBL **quotaMatrixOrig,
    SCHAR *indexVector,
    FIXP_DBL *nrgVector,
    DETECTOR_VALUES *detectorValues,
    INT startChannel,
    INT stopChannel,
    INT startIndex,
    INT stopIndex,
    INT numberOfStrongest
) {
  INT i, temp, j;

  const FIXP_DBL *filter = *fir_table[INVF_SMOOTHING_LENGTH];
  FIXP_DBL origQuotaMeanStrongest, sbrQuotaMeanStrongest;
  FIXP_DBL origQuota, sbrQuota;
  FIXP_DBL invIndex, invChannel, invTemp;
  FIXP_DBL quotaVecOrig[64], quotaVecSbr[64];

  FDKmemclear(quotaVecOrig, 64 * sizeof(FIXP_DBL));
  FDKmemclear(quotaVecSbr, 64 * sizeof(FIXP_DBL));

  invIndex = GetInvInt(stopIndex - startIndex);
  invChannel = GetInvInt(stopChannel - startChannel);

  detectorValues->avgNrg = FL2FXCONST_DBL(0.0f);
  for (j = startIndex; j < stopIndex; j++) {
    for (i = startChannel; i < stopChannel; i++) {
      quotaVecOrig[i] += fMult(quotaMatrixOrig[j][i], invIndex);

      if (indexVector[i] != -1)
        quotaVecSbr[i] += fMult(quotaMatrixOrig[j][indexVector[i]], invIndex);
    }
    detectorValues->avgNrg += fMult(nrgVector[j], invIndex);
  }

  origQuota = FL2FXCONST_DBL(0.0f);
  sbrQuota = FL2FXCONST_DBL(0.0f);
  for (i = startChannel; i < stopChannel; i++) {
    origQuota += fMultDiv2(quotaVecOrig[i], invChannel);
    sbrQuota += fMultDiv2(quotaVecSbr[i], invChannel);
  }

  FDKsbrEnc_Shellsort_fract(quotaVecOrig + startChannel,
                            stopChannel - startChannel);
  FDKsbrEnc_Shellsort_fract(quotaVecSbr + startChannel,
                            stopChannel - startChannel);

  origQuotaMeanStrongest = FL2FXCONST_DBL(0.0f);
  sbrQuotaMeanStrongest = FL2FXCONST_DBL(0.0f);

  temp = min(stopChannel - startChannel, numberOfStrongest);
  invTemp = GetInvInt(temp);

  for (i = 0; i < temp; i++) {
    origQuotaMeanStrongest +=
        fMultDiv2(quotaVecOrig[i + stopChannel - temp], invTemp);
    sbrQuotaMeanStrongest +=
        fMultDiv2(quotaVecSbr[i + stopChannel - temp], invTemp);
  }

  detectorValues->origQuotaMax = quotaVecOrig[stopChannel - 1];
  detectorValues->sbrQuotaMax = quotaVecSbr[stopChannel - 1];

  FDKmemmove(detectorValues->origQuotaMean, detectorValues->origQuotaMean + 1,
             INVF_SMOOTHING_LENGTH * sizeof(FIXP_DBL));
  FDKmemmove(detectorValues->sbrQuotaMean, detectorValues->sbrQuotaMean + 1,
             INVF_SMOOTHING_LENGTH * sizeof(FIXP_DBL));
  FDKmemmove(detectorValues->origQuotaMeanStrongest,
             detectorValues->origQuotaMeanStrongest + 1,
             INVF_SMOOTHING_LENGTH * sizeof(FIXP_DBL));
  FDKmemmove(detectorValues->sbrQuotaMeanStrongest,
             detectorValues->sbrQuotaMeanStrongest + 1,
             INVF_SMOOTHING_LENGTH * sizeof(FIXP_DBL));

  detectorValues->origQuotaMean[INVF_SMOOTHING_LENGTH] = origQuota << 1;
  detectorValues->sbrQuotaMean[INVF_SMOOTHING_LENGTH] = sbrQuota << 1;
  detectorValues->origQuotaMeanStrongest[INVF_SMOOTHING_LENGTH] =
      origQuotaMeanStrongest << 1;
  detectorValues->sbrQuotaMeanStrongest[INVF_SMOOTHING_LENGTH] =
      sbrQuotaMeanStrongest << 1;

  detectorValues->origQuotaMeanFilt = FL2FXCONST_DBL(0.0f);
  detectorValues->sbrQuotaMeanFilt = FL2FXCONST_DBL(0.0f);
  detectorValues->origQuotaMeanStrongestFilt = FL2FXCONST_DBL(0.0f);
  detectorValues->sbrQuotaMeanStrongestFilt = FL2FXCONST_DBL(0.0f);

  for (i = 0; i < INVF_SMOOTHING_LENGTH + 1; i++) {
    detectorValues->origQuotaMeanFilt +=
        fMult(detectorValues->origQuotaMean[i], filter[i]);
    detectorValues->sbrQuotaMeanFilt +=
        fMult(detectorValues->sbrQuotaMean[i], filter[i]);
    detectorValues->origQuotaMeanStrongestFilt +=
        fMult(detectorValues->origQuotaMeanStrongest[i], filter[i]);
    detectorValues->sbrQuotaMeanStrongestFilt +=
        fMult(detectorValues->sbrQuotaMeanStrongest[i], filter[i]);
  }
}

static INT findRegion(
    FIXP_DBL currVal,
    const FIXP_DBL *borders,
    const INT numBorders
) {
  INT i;

  if (currVal < borders[0]) {
    return 0;
  }

  for (i = 1; i < numBorders; i++) {
    if (currVal >= borders[i - 1] && currVal < borders[i]) {
      return i;
    }
  }

  if (currVal >= borders[numBorders - 1]) {
    return numBorders;
  }

  return 0;
}

static INVF_MODE decisionAlgorithm(
    const DETECTOR_PARAMETERS
        *detectorParams,
    DETECTOR_VALUES *detectorValues,
    INT transientFlag,
    INT *prevRegionSbr,
    INT *prevRegionOrig
) {
  INT invFiltLevel, regionSbr, regionOrig, regionNrg;

  const INT numRegionsSbr = detectorParams->numRegionsSbr;
  const INT numRegionsOrig = detectorParams->numRegionsOrig;
  const INT numRegionsNrg = detectorParams->numRegionsNrg;

  FIXP_DBL quantStepsSbrTmp[MAX_NUM_REGIONS];
  FIXP_DBL quantStepsOrigTmp[MAX_NUM_REGIONS];
  FIXP_DBL origQuotaMeanFilt;
  FIXP_DBL sbrQuotaMeanFilt;
  FIXP_DBL nrg;

  origQuotaMeanFilt =
      (fMultDiv2(FL2FXCONST_DBL(2.f * 0.375f),
                 (FIXP_DBL)(CalcLdData(max(detectorValues->origQuotaMeanFilt,
                                           (FIXP_DBL)1)) +
                            FL2FXCONST_DBL(0.31143075889f))))
      << 0;
  sbrQuotaMeanFilt =
      (fMultDiv2(FL2FXCONST_DBL(2.f * 0.375f),
                 (FIXP_DBL)(CalcLdData(max(detectorValues->sbrQuotaMeanFilt,
                                           (FIXP_DBL)1)) +
                            FL2FXCONST_DBL(0.31143075889f))))
      << 0;
  nrg =
      (fMultDiv2(FL2FXCONST_DBL(2.f * 0.375f),
                 (FIXP_DBL)(CalcLdData(detectorValues->avgNrg + (FIXP_DBL)1) +
                            FL2FXCONST_DBL(0.0625f) + FL2FXCONST_DBL(0.6875f))))
      << 0;

  FDKmemcpy(quantStepsSbrTmp, detectorParams->quantStepsSbr,
            numRegionsSbr * sizeof(FIXP_DBL));
  FDKmemcpy(quantStepsOrigTmp, detectorParams->quantStepsOrig,
            numRegionsOrig * sizeof(FIXP_DBL));

  if (*prevRegionSbr < numRegionsSbr)
    quantStepsSbrTmp[*prevRegionSbr] =
        detectorParams->quantStepsSbr[*prevRegionSbr] + hysteresis;
  if (*prevRegionSbr > 0)
    quantStepsSbrTmp[*prevRegionSbr - 1] =
        detectorParams->quantStepsSbr[*prevRegionSbr - 1] - hysteresis;

  if (*prevRegionOrig < numRegionsOrig)
    quantStepsOrigTmp[*prevRegionOrig] =
        detectorParams->quantStepsOrig[*prevRegionOrig] + hysteresis;
  if (*prevRegionOrig > 0)
    quantStepsOrigTmp[*prevRegionOrig - 1] =
        detectorParams->quantStepsOrig[*prevRegionOrig - 1] - hysteresis;

  regionSbr = findRegion(sbrQuotaMeanFilt, quantStepsSbrTmp, numRegionsSbr);
  regionOrig = findRegion(origQuotaMeanFilt, quantStepsOrigTmp, numRegionsOrig);
  regionNrg = findRegion(nrg, detectorParams->nrgBorders, numRegionsNrg);

  *prevRegionSbr = regionSbr;
  *prevRegionOrig = regionOrig;

  invFiltLevel =
      (transientFlag == 1)
          ? detectorParams->regionSpaceTransient[regionSbr][regionOrig]
          : detectorParams->regionSpace[regionSbr][regionOrig];

  invFiltLevel =
      max(invFiltLevel + detectorParams->EnergyCompFactor[regionNrg], 0);

  return (INVF_MODE)(invFiltLevel);
}

void FDKsbrEnc_qmfInverseFilteringDetector(
    HANDLE_SBR_INV_FILT_EST
        hInvFilt,
    FIXP_DBL **quotaMatrix,
    FIXP_DBL *nrgVector,
    SCHAR *indexVector,
    INT startIndex,
    INT stopIndex,
    INT transientFlag,
    INVF_MODE *infVec
) {
  INT band;

  for (band = 0; band < hInvFilt->noDetectorBands; band++) {
    INT startChannel = hInvFilt->freqBandTableInvFilt[band];
    INT stopChannel = hInvFilt->freqBandTableInvFilt[band + 1];

    calculateDetectorValues(quotaMatrix, indexVector, nrgVector,
                            &hInvFilt->detectorValues[band], startChannel,
                            stopChannel, startIndex, stopIndex,
                            hInvFilt->numberOfStrongest);

    infVec[band] = decisionAlgorithm(
        hInvFilt->detectorParams, &hInvFilt->detectorValues[band],
        transientFlag, &hInvFilt->prevRegionSbr[band],
        &hInvFilt->prevRegionOrig[band]);
  }
}

INT FDKsbrEnc_initInvFiltDetector(
    HANDLE_SBR_INV_FILT_EST
        hInvFilt,
    INT *freqBandTableDetector,
    INT numDetectorBands,
    UINT useSpeechConfig
) {
  INT i;

  FDKmemclear(hInvFilt, sizeof(SBR_INV_FILT_EST));

  hInvFilt->detectorParams =
      (useSpeechConfig) ? &detectorParamsAACSpeech : &detectorParamsAAC;

  hInvFilt->noDetectorBandsMax = numDetectorBands;

  for (i = 0; i < hInvFilt->noDetectorBandsMax; i++) {
    FDKmemclear(&hInvFilt->detectorValues[i], sizeof(DETECTOR_VALUES));
    hInvFilt->prevInvfMode[i] = INVF_OFF;
    hInvFilt->prevRegionOrig[i] = 0;
    hInvFilt->prevRegionSbr[i] = 0;
  }

  FDKsbrEnc_resetInvFiltDetector(hInvFilt, freqBandTableDetector,
                                 hInvFilt->noDetectorBandsMax);

  return (0);
}

INT FDKsbrEnc_resetInvFiltDetector(
    HANDLE_SBR_INV_FILT_EST
        hInvFilt,
    INT *freqBandTableDetector,
    INT numDetectorBands)
{
  hInvFilt->numberOfStrongest = 1;
  FDKmemcpy(hInvFilt->freqBandTableInvFilt, freqBandTableDetector,
            (numDetectorBands + 1) * sizeof(INT));
  hInvFilt->noDetectorBands = numDetectorBands;

  return (0);
}
