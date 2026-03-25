/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "ton_corr.h"
#include "sbrenc_ram.h"
#include "sbr_misc.h"
#include "genericStds.h"
#include "autocorr2nd.h"

#define BAND_V_SIZE 32
#define NUM_V_COMBINE \
  8

void FDKsbrEnc_CalculateTonalityQuotas(
    HANDLE_SBR_TON_CORR_EST hTonCorr,
    FIXP_DBL **RESTRICT
        sourceBufferReal,
    FIXP_DBL **RESTRICT
        sourceBufferImag,
    INT usb,
    INT qmfScale
) {
  INT i, k, r, r2, timeIndex, autoCorrScaling;

  INT startIndexMatrix = hTonCorr->startIndexMatrix;
  INT totNoEst = hTonCorr->numberOfEstimates;
  INT noEstPerFrame = hTonCorr->numberOfEstimatesPerFrame;
  INT move = hTonCorr->move;
  INT noQmfChannels = hTonCorr->noQmfChannels;
  INT buffLen = hTonCorr->bufferLength;
  INT stepSize = hTonCorr->stepSize;
  INT *pBlockLength = hTonCorr->lpcLength;
  INT **RESTRICT signMatrix = hTonCorr->signMatrix;
  FIXP_DBL *RESTRICT nrgVector = hTonCorr->nrgVector;
  FIXP_DBL **RESTRICT quotaMatrix = hTonCorr->quotaMatrix;
  FIXP_DBL *RESTRICT pNrgVectorFreq = hTonCorr->nrgVectorFreq;

  FIXP_DBL *realBuf;
  FIXP_DBL *imagBuf;

  FIXP_DBL alphar[2], alphai[2], fac;

  C_ALLOC_SCRATCH_START(ac, ACORR_COEFS, 1)
  C_ALLOC_SCRATCH_START(realBufRef, FIXP_DBL, 2 * BAND_V_SIZE * NUM_V_COMBINE)
  realBuf = realBufRef;
  imagBuf = realBuf + BAND_V_SIZE * NUM_V_COMBINE;

  FDK_ASSERT(buffLen <= BAND_V_SIZE);
  FDK_ASSERT(sizeof(FIXP_DBL) * NUM_V_COMBINE * BAND_V_SIZE * 2 <
             (1024 * sizeof(FIXP_DBL) - sizeof(ACORR_COEFS)));

  for (i = 0; i < move; i++) {
    FDKmemcpy(quotaMatrix[i], quotaMatrix[i + noEstPerFrame],
              noQmfChannels * sizeof(FIXP_DBL));
    FDKmemcpy(signMatrix[i], signMatrix[i + noEstPerFrame],
              noQmfChannels * sizeof(INT));
  }

  FDKmemmove(nrgVector, nrgVector + noEstPerFrame, move * sizeof(FIXP_DBL));
  FDKmemclear(nrgVector + startIndexMatrix,
              (totNoEst - startIndexMatrix) * sizeof(FIXP_DBL));
  FDKmemclear(pNrgVectorFreq, noQmfChannels * sizeof(FIXP_DBL));

  for (r = 0; r < usb; r++) {
    int blockLength;

    k = hTonCorr->nextSample;
    timeIndex = startIndexMatrix;
    if (realBuf != realBufRef) {
      realBuf -= BAND_V_SIZE;
      imagBuf -= BAND_V_SIZE;
    } else {
      realBuf += BAND_V_SIZE * (NUM_V_COMBINE - 1);
      imagBuf += BAND_V_SIZE * (NUM_V_COMBINE - 1);

      for (i = 0; i < buffLen; i++) {
        int v;
        FIXP_DBL *ptr;
        ptr = realBuf + i;
        for (v = 0; v < NUM_V_COMBINE; v++) {
          ptr[0] = sourceBufferReal[i][r + v];
          ptr[0 + BAND_V_SIZE * NUM_V_COMBINE] = sourceBufferImag[i][r + v];
          ptr -= BAND_V_SIZE;
        }
      }
    }

    blockLength = pBlockLength[0];

    while (k <= buffLen - blockLength) {
      autoCorrScaling = fixMin(
          getScalefactor(&realBuf[k - LPC_ORDER], LPC_ORDER + blockLength),
          getScalefactor(&imagBuf[k - LPC_ORDER], LPC_ORDER + blockLength));
      autoCorrScaling = fixMax(0, autoCorrScaling - 1);

      scaleValues(&realBuf[k - LPC_ORDER], LPC_ORDER + blockLength,
                  autoCorrScaling);
      scaleValues(&imagBuf[k - LPC_ORDER], LPC_ORDER + blockLength,
                  autoCorrScaling);

      autoCorrScaling <<= 1;
      autoCorrScaling +=
          autoCorr2nd_cplx(ac, realBuf + k, imagBuf + k, blockLength);

      if (ac->det == FL2FXCONST_DBL(0.0f)) {
        alphar[1] = alphai[1] = FL2FXCONST_DBL(0.0f);

        alphar[0] = (ac->r01r) >> 2;
        alphai[0] = (ac->r01i) >> 2;

        fac = fMultDiv2(ac->r00r, ac->r11r) >> 1;
      } else {
        alphar[1] = (fMultDiv2(ac->r01r, ac->r12r) >> 1) -
                    (fMultDiv2(ac->r01i, ac->r12i) >> 1) -
                    (fMultDiv2(ac->r02r, ac->r11r) >> 1);
        alphai[1] = (fMultDiv2(ac->r01i, ac->r12r) >> 1) +
                    (fMultDiv2(ac->r01r, ac->r12i) >> 1) -
                    (fMultDiv2(ac->r02i, ac->r11r) >> 1);

        alphar[0] = (fMultDiv2(ac->r01r, ac->det) >> (ac->det_scale + 1)) +
                    fMult(alphar[1], ac->r12r) + fMult(alphai[1], ac->r12i);
        alphai[0] = (fMultDiv2(ac->r01i, ac->det) >> (ac->det_scale + 1)) +
                    fMult(alphai[1], ac->r12r) - fMult(alphar[1], ac->r12i);

        fac = fMultDiv2(ac->r00r, fMult(ac->det, ac->r11r)) >>
              (ac->det_scale + 1);
      }

      if (fac == FL2FXCONST_DBL(0.0f)) {
        quotaMatrix[timeIndex][r] = FL2FXCONST_DBL(0.0f);
        signMatrix[timeIndex][r] = 0;
      } else {
        FIXP_DBL tmp, num, denom;
        INT numShift, denomShift, commonShift;
        INT sign;

        num = fMultDiv2(alphar[0], ac->r01r) + fMultDiv2(alphai[0], ac->r01i) -
              fMultDiv2(alphar[1], fMult(ac->r02r, ac->r11r)) -
              fMultDiv2(alphai[1], fMult(ac->r02i, ac->r11r));
        num = fixp_abs(num);

        denom = (fac >> 1) +
                (fMultDiv2(fac, RELAXATION_FRACT) >> RELAXATION_SHIFT) - num;
        denom = fixp_abs(denom);

        num = fMult(num, RELAXATION_FRACT);

        numShift = CountLeadingBits(num) - 2;
        num = scaleValue(num, numShift);

        denomShift = CountLeadingBits(denom);
        denom = (FIXP_DBL)denom << denomShift;

        if ((num > FL2FXCONST_DBL(0.0f)) && (denom != FL2FXCONST_DBL(0.0f))) {
          commonShift =
              fixMin(numShift - denomShift + RELAXATION_SHIFT, DFRACT_BITS - 1);
          if (commonShift < 0) {
            commonShift = -commonShift;
            tmp = schur_div(num, denom, 16);
            commonShift = fixMin(commonShift, CountLeadingBits(tmp));
            quotaMatrix[timeIndex][r] = tmp << commonShift;
          } else {
            quotaMatrix[timeIndex][r] =
                schur_div(num, denom, 16) >> commonShift;
          }
        } else {
          quotaMatrix[timeIndex][r] = FL2FXCONST_DBL(0.0f);
        }

        if (ac->r11r != FL2FXCONST_DBL(0.0f)) {
          if (((ac->r01r >= FL2FXCONST_DBL(0.0f)) &&
               (ac->r11r >= FL2FXCONST_DBL(0.0f))) ||
              ((ac->r01r < FL2FXCONST_DBL(0.0f)) &&
               (ac->r11r < FL2FXCONST_DBL(0.0f)))) {
            sign = 1;
          } else {
            sign = -1;
          }
        } else {
          sign = 1;
        }

        if (sign < 0) {
          r2 = r;
        } else {
          r2 = r + 1;
        }
        signMatrix[timeIndex][r] = 1 - 2 * (r2 & 0x1);
      }

      nrgVector[timeIndex] +=
          ((ac->r00r) >>
           fixMin(DFRACT_BITS - 1,
                  (2 * qmfScale + autoCorrScaling + SCALE_NRGVEC)));
      pNrgVectorFreq[r] =
          pNrgVectorFreq[r] +
          ((ac->r00r) >>
           fixMin(DFRACT_BITS - 1,
                  (2 * qmfScale + autoCorrScaling + SCALE_NRGVEC)));

      blockLength = pBlockLength[1];
      k += stepSize;
      timeIndex++;
    }
  }

  C_ALLOC_SCRATCH_END(realBufRef, FIXP_DBL, 2 * BAND_V_SIZE * NUM_V_COMBINE)
  C_ALLOC_SCRATCH_END(ac, ACORR_COEFS, 1)
}

void FDKsbrEnc_TonCorrParamExtr(
    HANDLE_SBR_TON_CORR_EST hTonCorr,
    INVF_MODE *infVec,
    FIXP_DBL *noiseLevels,
    INT *missingHarmonicFlag,
    UCHAR *missingHarmonicsIndex,
    UCHAR *envelopeCompensation,
    const SBR_FRAME_INFO *frameInfo,
    UCHAR *transientInfo,
    UCHAR *freqBandTable,
    INT nSfb,
    XPOS_MODE xposType,
    UINT sbrSyntaxFlags) {
  INT band;
  INT transientFlag = transientInfo[1];
  INT transientPos = transientInfo[0];
  INT transientFrame, transientFrameInvfEst;
  INVF_MODE *infVecPtr;

  transientFrame = 0;
  if (hTonCorr->transientNextFrame) {
    transientFrame = 1;
    hTonCorr->transientNextFrame = 0;

    if (transientFlag) {
      if (transientPos + hTonCorr->transientPosOffset >=
          frameInfo->borders[frameInfo->nEnvelopes]) {
        hTonCorr->transientNextFrame = 1;
      }
    }
  } else {
    if (transientFlag) {
      if (transientPos + hTonCorr->transientPosOffset <
          frameInfo->borders[frameInfo->nEnvelopes]) {
        transientFrame = 1;
        hTonCorr->transientNextFrame = 0;
      } else {
        hTonCorr->transientNextFrame = 1;
      }
    }
  }
  transientFrameInvfEst = transientFrame;

  if (hTonCorr->switchInverseFilt)
    FDKsbrEnc_qmfInverseFilteringDetector(
        &hTonCorr->sbrInvFilt, hTonCorr->quotaMatrix, hTonCorr->nrgVector,
        hTonCorr->indexVector, hTonCorr->frameStartIndexInvfEst,
        hTonCorr->numberOfEstimatesPerFrame + hTonCorr->frameStartIndexInvfEst,
        transientFrameInvfEst, infVec);

  if (xposType == XPOS_LC) {
    FDKsbrEnc_SbrMissingHarmonicsDetectorQmf(
        &hTonCorr->sbrMissingHarmonicsDetector, hTonCorr->quotaMatrix,
        hTonCorr->signMatrix, hTonCorr->indexVector, frameInfo, transientInfo,
        missingHarmonicFlag, missingHarmonicsIndex, freqBandTable, nSfb,
        envelopeCompensation, hTonCorr->nrgVectorFreq);
  } else {
    *missingHarmonicFlag = 0;
    FDKmemclear(missingHarmonicsIndex, nSfb * sizeof(UCHAR));
  }

  infVecPtr = hTonCorr->sbrInvFilt.prevInvfMode;

  FDKsbrEnc_sbrNoiseFloorEstimateQmf(
      &hTonCorr->sbrNoiseFloorEstimate, frameInfo, noiseLevels,
      hTonCorr->quotaMatrix, hTonCorr->indexVector, *missingHarmonicFlag,
      hTonCorr->frameStartIndex, hTonCorr->numberOfEstimatesPerFrame,
      transientFrame, infVecPtr, sbrSyntaxFlags);

  for (band = 0; band < hTonCorr->sbrInvFilt.noDetectorBands; band++) {
    hTonCorr->sbrInvFilt.prevInvfMode[band] = infVec[band];
  }
}

static INT findClosestEntry(INT goalSb, UCHAR *v_k_master, INT numMaster,
                            INT direction) {
  INT index;

  if (goalSb <= v_k_master[0]) return v_k_master[0];

  if (goalSb >= v_k_master[numMaster]) return v_k_master[numMaster];

  if (direction) {
    index = 0;
    while (v_k_master[index] < goalSb) {
      index++;
    }
  } else {
    index = numMaster;
    while (v_k_master[index] > goalSb) {
      index--;
    }
  }

  return v_k_master[index];
}

static INT resetPatch(
    HANDLE_SBR_TON_CORR_EST hTonCorr,
    INT xposctrl,
    INT highBandStartSb,
    UCHAR *v_k_master,
    INT numMaster,
    INT fs,
    INT noChannels)
{
  INT patch, k, i;
  INT targetStopBand;

  PATCH_PARAM *patchParam = hTonCorr->patchParam;

  INT sbGuard = hTonCorr->guard;
  INT sourceStartBand;
  INT patchDistance;
  INT numBandsInPatch;

  INT lsb =
      v_k_master[0];
  INT usb = v_k_master[numMaster];
  INT xoverOffset =
      highBandStartSb -
      v_k_master[0];

  INT goalSb;

  if (xposctrl == 1) {
    lsb += xoverOffset;
    xoverOffset = 0;
  }

  goalSb = (INT)((2 * noChannels * 16000 + (fs >> 1)) / fs);
  goalSb = findClosestEntry(goalSb, v_k_master, numMaster,
                            1);

  sourceStartBand = hTonCorr->shiftStartSb + xoverOffset;
  targetStopBand = lsb + xoverOffset;

  patch = 0;
  while (targetStopBand < usb) {
    if (patch >= MAX_NUM_PATCHES) return (1);

    patchParam[patch].guardStartBand = targetStopBand;
    targetStopBand += sbGuard;
    patchParam[patch].targetStartBand = targetStopBand;

    numBandsInPatch =
        goalSb - targetStopBand;

    if (numBandsInPatch >= lsb - sourceStartBand) {
      patchDistance =
          targetStopBand - sourceStartBand;
      patchDistance =
          patchDistance & ~1;
      numBandsInPatch = lsb - (targetStopBand - patchDistance);
      numBandsInPatch = findClosestEntry(targetStopBand + numBandsInPatch,
                                         v_k_master, numMaster, 0) -
                        targetStopBand;
    }

    patchDistance =
        numBandsInPatch + targetStopBand - lsb;
    patchDistance = (patchDistance + 1) &
                    ~1;

    if (numBandsInPatch <= 0) {
      patch--;
    } else {
      patchParam[patch].sourceStartBand = targetStopBand - patchDistance;
      patchParam[patch].targetBandOffs = patchDistance;
      patchParam[patch].numBandsInPatch = numBandsInPatch;
      patchParam[patch].sourceStopBand =
          patchParam[patch].sourceStartBand + numBandsInPatch;

      targetStopBand += patchParam[patch].numBandsInPatch;
    }

    sourceStartBand = hTonCorr->shiftStartSb;

    if (fixp_abs(targetStopBand - goalSb) < 3) {
      goalSb = usb;
    }

    patch++;
  }

  patch--;

  if (patchParam[patch].numBandsInPatch < 3 && patch > 0) {
    patch--;
  }

  hTonCorr->noOfPatches = patch + 1;

  for (k = 0; k < hTonCorr->patchParam[0].guardStartBand; k++)
    hTonCorr->indexVector[k] = k;

  for (i = 0; i < hTonCorr->noOfPatches; i++) {
    INT sourceStart = hTonCorr->patchParam[i].sourceStartBand;
    INT targetStart = hTonCorr->patchParam[i].targetStartBand;
    INT numberOfBands = hTonCorr->patchParam[i].numBandsInPatch;
    INT startGuardBand = hTonCorr->patchParam[i].guardStartBand;

    for (k = 0; k < (targetStart - startGuardBand); k++)
      hTonCorr->indexVector[startGuardBand + k] = -1;

    for (k = 0; k < numberOfBands; k++)
      hTonCorr->indexVector[targetStart + k] = sourceStart + k;
  }

  return (0);
}

INT FDKsbrEnc_CreateTonCorrParamExtr(
    HANDLE_SBR_TON_CORR_EST
        hTonCorr,
    INT chan)
{
  INT i;
  FIXP_DBL *quotaMatrix = GetRam_Sbr_quotaMatrix(chan);
  INT *signMatrix = GetRam_Sbr_signMatrix(chan);

  if ((NULL == quotaMatrix) || (NULL == signMatrix)) {
    goto bail;
  }

  FDKmemclear(hTonCorr, sizeof(SBR_TON_CORR_EST));

  for (i = 0; i < MAX_NO_OF_ESTIMATES; i++) {
    hTonCorr->quotaMatrix[i] = quotaMatrix + (i * 64);
    hTonCorr->signMatrix[i] = signMatrix + (i * 64);
  }

  if (0 != FDKsbrEnc_CreateSbrMissingHarmonicsDetector(
               &hTonCorr->sbrMissingHarmonicsDetector, chan)) {
    goto bail;
  }

  return 0;

bail:
  hTonCorr->quotaMatrix[0] = quotaMatrix;
  hTonCorr->signMatrix[0] = signMatrix;

  FDKsbrEnc_DeleteTonCorrParamExtr(hTonCorr);

  return -1;
}

INT FDKsbrEnc_InitTonCorrParamExtr(
    INT frameSize,
    HANDLE_SBR_TON_CORR_EST
        hTonCorr,
    HANDLE_SBR_CONFIG_DATA
        sbrCfg,
    INT timeSlots,
    INT xposCtrl,
    INT ana_max_level,
    INT noiseBands,
    INT noiseFloorOffset,
    UINT useSpeechConfig)
{
  INT nCols = sbrCfg->noQmfSlots;
  INT fs = sbrCfg->sampleFreq;
  INT noQmfChannels = sbrCfg->noQmfBands;

  INT highBandStartSb = sbrCfg->freqBandTable[LOW_RES][0];
  UCHAR *v_k_master = sbrCfg->v_k_master;
  INT numMaster = sbrCfg->num_Master;

  UCHAR **freqBandTable = sbrCfg->freqBandTable;
  INT *nSfb = sbrCfg->nSfb;

  INT i;

  if (sbrCfg->sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY) {
    switch (timeSlots) {
      case NUMBER_TIME_SLOTS_1920:
        hTonCorr->lpcLength[0] = 8 - LPC_ORDER;
        hTonCorr->lpcLength[1] = 7 - LPC_ORDER;
        hTonCorr->numberOfEstimates = NO_OF_ESTIMATES_LD;
        hTonCorr->numberOfEstimatesPerFrame = 2;
        hTonCorr->frameStartIndexInvfEst = 0;
        hTonCorr->transientPosOffset = FRAME_MIDDLE_SLOT_512LD;
        break;
      case NUMBER_TIME_SLOTS_2048:
        hTonCorr->lpcLength[0] = 8 - LPC_ORDER;
        hTonCorr->lpcLength[1] = 8 - LPC_ORDER;
        hTonCorr->numberOfEstimates = NO_OF_ESTIMATES_LD;
        hTonCorr->numberOfEstimatesPerFrame = 2;
        hTonCorr->frameStartIndexInvfEst = 0;
        hTonCorr->transientPosOffset = FRAME_MIDDLE_SLOT_512LD;
        break;
    }
  } else
    switch (timeSlots) {
      case NUMBER_TIME_SLOTS_2048:
        hTonCorr->lpcLength[0] = 16 - LPC_ORDER;
        hTonCorr->lpcLength[1] = 16 - LPC_ORDER;
        hTonCorr->numberOfEstimates = NO_OF_ESTIMATES_LC;
        hTonCorr->numberOfEstimatesPerFrame = sbrCfg->noQmfSlots / 16;
        hTonCorr->frameStartIndexInvfEst = 0;
        hTonCorr->transientPosOffset = FRAME_MIDDLE_SLOT_2048;
        break;
      case NUMBER_TIME_SLOTS_1920:
        hTonCorr->lpcLength[0] = 15 - LPC_ORDER;
        hTonCorr->lpcLength[1] = 15 - LPC_ORDER;
        hTonCorr->numberOfEstimates = NO_OF_ESTIMATES_LC;
        hTonCorr->numberOfEstimatesPerFrame = sbrCfg->noQmfSlots / 15;
        hTonCorr->frameStartIndexInvfEst = 0;
        hTonCorr->transientPosOffset = FRAME_MIDDLE_SLOT_1920;
        break;
      default:
        return -1;
    }

  hTonCorr->bufferLength = nCols;
  hTonCorr->stepSize =
      hTonCorr->lpcLength[0] + LPC_ORDER;

  hTonCorr->nextSample = LPC_ORDER;
  hTonCorr->move = hTonCorr->numberOfEstimates -
                   hTonCorr->numberOfEstimatesPerFrame;
  if (hTonCorr->move < 0) {
    return -1;
  }
  hTonCorr->startIndexMatrix =
      hTonCorr->numberOfEstimates -
      hTonCorr->numberOfEstimatesPerFrame;
  hTonCorr->frameStartIndex = 0;
  hTonCorr->prevTransientFlag = 0;
  hTonCorr->transientNextFrame = 0;

  hTonCorr->noQmfChannels = noQmfChannels;

  for (i = 0; i < hTonCorr->numberOfEstimates; i++) {
    FDKmemclear(hTonCorr->quotaMatrix[i], sizeof(FIXP_DBL) * noQmfChannels);
    FDKmemclear(hTonCorr->signMatrix[i], sizeof(INT) * noQmfChannels);
  }

  hTonCorr->guard = 0;
  hTonCorr->shiftStartSb = 1;

  if (resetPatch(hTonCorr, xposCtrl, highBandStartSb, v_k_master, numMaster, fs,
                 noQmfChannels))
    return (1);

  if (FDKsbrEnc_InitSbrNoiseFloorEstimate(
          &hTonCorr->sbrNoiseFloorEstimate, ana_max_level, freqBandTable[LO],
          nSfb[LO], noiseBands, noiseFloorOffset, timeSlots, useSpeechConfig))
    return (1);

  if (FDKsbrEnc_initInvFiltDetector(
          &hTonCorr->sbrInvFilt,
          hTonCorr->sbrNoiseFloorEstimate.freqBandTableQmf,
          hTonCorr->sbrNoiseFloorEstimate.noNoiseBands, useSpeechConfig))
    return (1);

  if (FDKsbrEnc_InitSbrMissingHarmonicsDetector(
          &hTonCorr->sbrMissingHarmonicsDetector, fs, frameSize, nSfb[HI],
          noQmfChannels, hTonCorr->numberOfEstimates, hTonCorr->move,
          hTonCorr->numberOfEstimatesPerFrame, sbrCfg->sbrSyntaxFlags))
    return (1);

  return (0);
}

INT FDKsbrEnc_ResetTonCorrParamExtr(
    HANDLE_SBR_TON_CORR_EST hTonCorr,
    INT xposctrl,
    INT highBandStartSb,
    UCHAR *v_k_master,
    INT numMaster,
    INT fs,
    UCHAR *
        *freqBandTable,
    INT *nSfb,
    INT noQmfChannels
) {
  hTonCorr->guard = 0;
  hTonCorr->shiftStartSb = 1;

  if (resetPatch(hTonCorr, xposctrl, highBandStartSb, v_k_master, numMaster, fs,
                 noQmfChannels))
    return (1);

  if (FDKsbrEnc_resetSbrNoiseFloorEstimate(&hTonCorr->sbrNoiseFloorEstimate,
                                           freqBandTable[LO], nSfb[LO]))
    return (1);

  if (FDKsbrEnc_resetInvFiltDetector(
          &hTonCorr->sbrInvFilt,
          hTonCorr->sbrNoiseFloorEstimate.freqBandTableQmf,
          hTonCorr->sbrNoiseFloorEstimate.noNoiseBands))
    return (1);
  if (FDKsbrEnc_ResetSbrMissingHarmonicsDetector(
          &hTonCorr->sbrMissingHarmonicsDetector, nSfb[HI]))
    return (1);

  return (0);
}

void FDKsbrEnc_DeleteTonCorrParamExtr(
    HANDLE_SBR_TON_CORR_EST hTonCorr)
{
  if (hTonCorr) {
    FreeRam_Sbr_quotaMatrix(hTonCorr->quotaMatrix);

    FreeRam_Sbr_signMatrix(hTonCorr->signMatrix);

    FDKsbrEnc_DeleteSbrMissingHarmonicsDetector(
        &hTonCorr->sbrMissingHarmonicsDetector);
  }
}
