/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "psdec.h"
#include "FDK_bitbuffer.h"
#include "sbr_rom.h"
#include "sbr_ram.h"
#include "FDK_tools_rom.h"
#include "genericStds.h"
#include "FDK_trigFcts.h"

#define FRACT_ZERO FRACT_BITS - 1

SBR_ERROR ResetPsDec(HANDLE_PS_DEC h_ps_d);

int CreatePsDec(HANDLE_PS_DEC *h_PS_DEC,
                int aacSamplesPerFrame) {
  SBR_ERROR errorInfo = SBRDEC_OK;
  HANDLE_PS_DEC h_ps_d;
  int i;

  if (*h_PS_DEC == NULL) {
    h_ps_d = GetRam_ps_dec();
    if (h_ps_d == NULL) {
      goto bail;
    }
  } else {
    h_ps_d = *h_PS_DEC;
  }

  FDKhybridAnalysisOpen(&h_ps_d->specificTo.mpeg.hybridAnalysis,
                        h_ps_d->specificTo.mpeg.pHybridAnaStatesLFdmx,
                        sizeof(h_ps_d->specificTo.mpeg.pHybridAnaStatesLFdmx),
                        NULL, 0);

  switch (aacSamplesPerFrame) {
    case 960:
      h_ps_d->noSubSamples = 30;
      break;
    case 1024:
      h_ps_d->noSubSamples = 32;
      break;
    default:
      h_ps_d->noSubSamples = -1;
      break;
  }

  if (h_ps_d->noSubSamples > MAX_NUM_COL || h_ps_d->noSubSamples <= 0) {
    goto bail;
  }
  h_ps_d->noChannels = NO_QMF_CHANNELS;

  h_ps_d->psDecodedPrv = 0;
  h_ps_d->procFrameBased = -1;
  for (i = 0; i < (1) + 1; i++) {
    h_ps_d->bPsDataAvail[i] = ppt_none;
  }
  {
    int error;
    error = FDKdecorrelateOpen(&(h_ps_d->specificTo.mpeg.apDecor),
                               h_ps_d->specificTo.mpeg.decorrBufferCplx,
                               (2 * ((825) + (373))));
    if (error) goto bail;
  }

  for (i = 0; i < (1) + 1; i++) {
    FDKmemclear(&h_ps_d->bsData[i].mpeg, sizeof(MPEG_PS_BS_DATA));
  }

  errorInfo = ResetPsDec(h_ps_d);

  if (errorInfo != SBRDEC_OK) goto bail;

  *h_PS_DEC = h_ps_d;

  return 0;

bail:
  if (h_ps_d != NULL) {
    DeletePsDec(&h_ps_d);
  }

  return -1;
}

int DeletePsDec(HANDLE_PS_DEC *h_PS_DEC)
{
  if (*h_PS_DEC == NULL) {
    return -1;
  }

  {
    HANDLE_PS_DEC h_ps_d = *h_PS_DEC;
    FDKdecorrelateClose(&(h_ps_d->specificTo.mpeg.apDecor));
  }

  FreeRam_ps_dec(h_PS_DEC);

  return 0;
}

SBR_ERROR ResetPsDec(HANDLE_PS_DEC h_ps_d)
{
  SBR_ERROR errorInfo = SBRDEC_OK;
  INT i;

  h_ps_d->specificTo.mpeg.lastUsb = 0;

  FDKhybridAnalysisInit(&h_ps_d->specificTo.mpeg.hybridAnalysis, THREE_TO_TEN,
                        NO_QMF_BANDS_HYBRID20, NO_QMF_BANDS_HYBRID20, 1);

  for (i = 0; i < 2; i++) {
    FDKhybridSynthesisInit(&h_ps_d->specificTo.mpeg.hybridSynthesis[i],
                           THREE_TO_TEN, NO_QMF_CHANNELS, NO_QMF_CHANNELS);
  }
  {
    INT error;
    error = FDKdecorrelateInit(&h_ps_d->specificTo.mpeg.apDecor, 71, DECORR_PS,
                               DUCKER_AUTOMATIC, 0, 0, 0, 0, 1,
                               1);
    if (error) return SBRDEC_NOT_INITIALIZED;
  }

  for (i = 0; i < NO_IID_GROUPS; i++) {
    h_ps_d->specificTo.mpeg.h11rPrev[i] = FL2FXCONST_DBL(0.5f);
    h_ps_d->specificTo.mpeg.h12rPrev[i] = FL2FXCONST_DBL(0.5f);
  }

  FDKmemclear(h_ps_d->specificTo.mpeg.h21rPrev,
              sizeof(h_ps_d->specificTo.mpeg.h21rPrev));
  FDKmemclear(h_ps_d->specificTo.mpeg.h22rPrev,
              sizeof(h_ps_d->specificTo.mpeg.h22rPrev));

  return errorInfo;
}

void PreparePsProcessing(HANDLE_PS_DEC h_ps_d,
                         const FIXP_DBL *const *const rIntBufferLeft,
                         const FIXP_DBL *const *const iIntBufferLeft,
                         const int scaleFactorLowBand) {
  if (h_ps_d->procFrameBased ==
      1)
  {
    int i, j;

    for (i = 0; i < HYBRID_FILTER_DELAY; i++) {
      FIXP_DBL qmfInputData[2][NO_QMF_BANDS_HYBRID20];
      FIXP_DBL hybridOutputData[2][NO_SUB_QMF_CHANNELS];

      for (j = 0; j < NO_QMF_BANDS_HYBRID20; j++) {
        qmfInputData[0][j] =
            scaleValue(rIntBufferLeft[i][j], scaleFactorLowBand);
        qmfInputData[1][j] =
            scaleValue(iIntBufferLeft[i][j], scaleFactorLowBand);
      }

      FDKhybridAnalysisApply(&h_ps_d->specificTo.mpeg.hybridAnalysis,
                             qmfInputData[0], qmfInputData[1],
                             hybridOutputData[0], hybridOutputData[1]);
    }
    h_ps_d->procFrameBased = 0;
  }
}

void initSlotBasedRotation(
    HANDLE_PS_DEC h_ps_d,
    int env, int usb) {
  INT group = 0;
  INT bin = 0;
  INT noIidSteps;

  FIXP_SGL invL;
  FIXP_DBL ScaleL, ScaleR;
  FIXP_DBL Alpha, Beta;
  FIXP_DBL h11r, h12r, h21r, h22r;

  const FIXP_DBL *PScaleFactors;

  if (h_ps_d->bsData[h_ps_d->processSlot].mpeg.bFineIidQ) {
    PScaleFactors = ScaleFactorsFine;
    noIidSteps = NO_IID_STEPS_FINE;
  } else {
    PScaleFactors = ScaleFactors;
    noIidSteps = NO_IID_STEPS;
  }

  for (group = 0; group < NO_IID_GROUPS; group++) {
    bin = bins2groupMap20[group];

    ScaleR = PScaleFactors[noIidSteps + h_ps_d->specificTo.mpeg.pCoef
                                            ->aaIidIndexMapped[env][bin]];
    ScaleL = PScaleFactors[noIidSteps - h_ps_d->specificTo.mpeg.pCoef
                                            ->aaIidIndexMapped[env][bin]];

    Beta = fMult(
        fMult(Alphas[h_ps_d->specificTo.mpeg.pCoef->aaIccIndexMapped[env][bin]],
              (ScaleR - ScaleL)),
        FIXP_SQRT05);
    Alpha =
        Alphas[h_ps_d->specificTo.mpeg.pCoef->aaIccIndexMapped[env][bin]] >> 1;

    {
      FIXP_DBL trigData[4];

      inline_fixp_cos_sin(Beta + Alpha, Beta - Alpha, 2, trigData);
      h11r = fMult(ScaleL, trigData[0]);
      h12r = fMult(ScaleR, trigData[2]);
      h21r = fMult(ScaleL, trigData[1]);
      h22r = fMult(ScaleR, trigData[3]);
    }

    invL = FX_DBL2FX_SGL(GetInvInt(
        h_ps_d->bsData[h_ps_d->processSlot].mpeg.aEnvStartStop[env + 1] -
        h_ps_d->bsData[h_ps_d->processSlot].mpeg.aEnvStartStop[env]));

    h_ps_d->specificTo.mpeg.pCoef->H11r[group] =
        h_ps_d->specificTo.mpeg.h11rPrev[group];
    h_ps_d->specificTo.mpeg.pCoef->H12r[group] =
        h_ps_d->specificTo.mpeg.h12rPrev[group];
    h_ps_d->specificTo.mpeg.pCoef->H21r[group] =
        h_ps_d->specificTo.mpeg.h21rPrev[group];
    h_ps_d->specificTo.mpeg.pCoef->H22r[group] =
        h_ps_d->specificTo.mpeg.h22rPrev[group];

    h_ps_d->specificTo.mpeg.pCoef->DeltaH11r[group] =
        fMult(h11r - h_ps_d->specificTo.mpeg.pCoef->H11r[group], invL);
    h_ps_d->specificTo.mpeg.pCoef->DeltaH12r[group] =
        fMult(h12r - h_ps_d->specificTo.mpeg.pCoef->H12r[group], invL);
    h_ps_d->specificTo.mpeg.pCoef->DeltaH21r[group] =
        fMult(h21r - h_ps_d->specificTo.mpeg.pCoef->H21r[group], invL);
    h_ps_d->specificTo.mpeg.pCoef->DeltaH22r[group] =
        fMult(h22r - h_ps_d->specificTo.mpeg.pCoef->H22r[group], invL);

    h_ps_d->specificTo.mpeg.h11rPrev[group] = h11r;
    h_ps_d->specificTo.mpeg.h12rPrev[group] = h12r;
    h_ps_d->specificTo.mpeg.h21rPrev[group] = h21r;
    h_ps_d->specificTo.mpeg.h22rPrev[group] = h22r;
  }
}

static const UCHAR groupTable[NO_IID_GROUPS + 1] = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
    12, 13, 14, 15, 16, 18, 21, 25, 30, 42, 71};

static void applySlotBasedRotation(
    HANDLE_PS_DEC h_ps_d,

    FIXP_DBL *mHybridRealLeft,
    FIXP_DBL *mHybridImagLeft,

    FIXP_DBL *mHybridRealRight,
    FIXP_DBL *mHybridImagRight
) {
  INT group;
  INT subband;

  PS_DEC_COEFFICIENTS *pCoef = h_ps_d->specificTo.mpeg.pCoef;

  for (group = 0; group < NO_IID_GROUPS; group++) {
    pCoef->H11r[group] += pCoef->DeltaH11r[group];
    pCoef->H12r[group] += pCoef->DeltaH12r[group];
    pCoef->H21r[group] += pCoef->DeltaH21r[group];
    pCoef->H22r[group] += pCoef->DeltaH22r[group];

    const int start = groupTable[group];
    const int stop = groupTable[group + 1];
    for (subband = start; subband < stop; subband++) {
      FIXP_DBL tmpLeft =
          fMultAdd(fMultDiv2(pCoef->H11r[group], mHybridRealLeft[subband]),
                   pCoef->H21r[group], mHybridRealRight[subband]);
      FIXP_DBL tmpRight =
          fMultAdd(fMultDiv2(pCoef->H12r[group], mHybridRealLeft[subband]),
                   pCoef->H22r[group], mHybridRealRight[subband]);
      mHybridRealLeft[subband] = tmpLeft;
      mHybridRealRight[subband] = tmpRight;

      tmpLeft =
          fMultAdd(fMultDiv2(pCoef->H11r[group], mHybridImagLeft[subband]),
                   pCoef->H21r[group], mHybridImagRight[subband]);
      tmpRight =
          fMultAdd(fMultDiv2(pCoef->H12r[group], mHybridImagLeft[subband]),
                   pCoef->H22r[group], mHybridImagRight[subband]);
      mHybridImagLeft[subband] = tmpLeft;
      mHybridImagRight[subband] = tmpRight;
    }
  }
}

void ApplyPsSlot(
    HANDLE_PS_DEC h_ps_d,
    FIXP_DBL **rIntBufferLeft,
    FIXP_DBL **iIntBufferLeft,
    FIXP_DBL *rIntBufferRight,
    FIXP_DBL *iIntBufferRight,
    const int scaleFactorLowBand_no_ov, const int scaleFactorLowBand,
    const int scaleFactorHighBand, const int lsb, const int usb) {
#define NO_HYBRID_DATA_BANDS (71)

  int i;
  FIXP_DBL qmfInputData[2][NO_QMF_BANDS_HYBRID20];
  FIXP_DBL *hybridData[2][2];
  C_ALLOC_SCRATCH_START(pHybridData, FIXP_DBL, 4 * NO_HYBRID_DATA_BANDS);

  hybridData[0][0] =
      pHybridData + 0 * NO_HYBRID_DATA_BANDS;
  hybridData[0][1] =
      pHybridData + 1 * NO_HYBRID_DATA_BANDS;
  hybridData[1][0] =
      pHybridData + 2 * NO_HYBRID_DATA_BANDS;
  hybridData[1][1] =
      pHybridData + 3 * NO_HYBRID_DATA_BANDS;

  for (i = 0; i < NO_QMF_BANDS_HYBRID20; i++) {
    qmfInputData[0][i] = scaleValue(rIntBufferLeft[HYBRID_FILTER_DELAY][i],
                                    scaleFactorLowBand_no_ov);
    qmfInputData[1][i] = scaleValue(iIntBufferLeft[HYBRID_FILTER_DELAY][i],
                                    scaleFactorLowBand_no_ov);
  }

  FDKhybridAnalysisApply(&h_ps_d->specificTo.mpeg.hybridAnalysis,
                         qmfInputData[0], qmfInputData[1], hybridData[0][0],
                         hybridData[0][1]);

  scaleValues(&hybridData[0][0][NO_SUB_QMF_CHANNELS - 2],
              &rIntBufferLeft[0][NO_QMF_BANDS_HYBRID20],
              lsb - NO_QMF_BANDS_HYBRID20, scaleFactorLowBand);
  scaleValues(&hybridData[0][1][NO_SUB_QMF_CHANNELS - 2],
              &iIntBufferLeft[0][NO_QMF_BANDS_HYBRID20],
              lsb - NO_QMF_BANDS_HYBRID20, scaleFactorLowBand);

  scaleValues(&hybridData[0][0][lsb + (NO_SUB_QMF_CHANNELS - 2 -
                                       NO_QMF_BANDS_HYBRID20)],
              &rIntBufferLeft[0][lsb], usb - lsb, scaleFactorHighBand);
  scaleValues(&hybridData[0][1][lsb + (NO_SUB_QMF_CHANNELS - 2 -
                                       NO_QMF_BANDS_HYBRID20)],
              &iIntBufferLeft[0][lsb], usb - lsb, scaleFactorHighBand);

  FDKmemcpy(
      &hybridData[0][0]
                 [usb + (NO_SUB_QMF_CHANNELS - 2 - NO_QMF_BANDS_HYBRID20)],
      &rIntBufferLeft[0][usb], sizeof(FIXP_DBL) * (NO_QMF_CHANNELS - usb));
  FDKmemcpy(
      &hybridData[0][1]
                 [usb + (NO_SUB_QMF_CHANNELS - 2 - NO_QMF_BANDS_HYBRID20)],
      &iIntBufferLeft[0][usb], sizeof(FIXP_DBL) * (NO_QMF_CHANNELS - usb));

  FDKdecorrelateApply(&h_ps_d->specificTo.mpeg.apDecor,
                      &hybridData[0][0][0],
                      &hybridData[0][1][0],
                      &hybridData[1][0][0],
                      &hybridData[1][1][0],
                      0
  );

  applySlotBasedRotation(h_ps_d,
                         &hybridData[0][0][0],
                         &hybridData[0][1][0],
                         &hybridData[1][0][0],
                         &hybridData[1][1][0]
  );

  for (i = 0; i < 2; i++) {
    FDKhybridSynthesisApply(
        &h_ps_d->specificTo.mpeg.hybridSynthesis[i],
        hybridData[i][0],
        hybridData[i][1],
        (i == 0) ? rIntBufferLeft[0]
                 : rIntBufferRight,
        (i == 0) ? iIntBufferLeft[0]
                 : iIntBufferRight
    );
  }

  C_ALLOC_SCRATCH_END(pHybridData, FIXP_DBL, 4 * NO_HYBRID_DATA_BANDS);
}
