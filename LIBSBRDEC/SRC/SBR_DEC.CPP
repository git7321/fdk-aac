/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "sbr_dec.h"
#include "sbr_ram.h"
#include "env_extr.h"
#include "env_calc.h"
#include "scale.h"
#include "FDK_matrixCalloc.h"
#include "hbe.h"
#include "genericStds.h"
#include "sbrdec_drc.h"

static void copyHarmonicSpectrum(int *xOverQmf, FIXP_DBL **qmfReal,
                                 FIXP_DBL **qmfImag, int noCols, int overlap,
                                 KEEP_STATES_SYNCED_MODE keepStatesSynced) {
  int patchBands;
  int patch, band, col, target, sourceBands, i;
  int numPatches = 0;
  int slotOffset = 0;

  FIXP_DBL **ppqmfReal = qmfReal + overlap;
  FIXP_DBL **ppqmfImag = qmfImag + overlap;

  if (keepStatesSynced == KEEP_STATES_SYNCED_NORMAL) {
    slotOffset = noCols - overlap - LPC_ORDER;
  }

  if (keepStatesSynced == KEEP_STATES_SYNCED_OUTDIFF) {
    ppqmfReal = qmfReal;
    ppqmfImag = qmfImag;
  }

  for (i = 1; i < MAX_NUM_PATCHES; i++) {
    if (xOverQmf[i] != 0) {
      numPatches++;
    }
  }

  for (patch = (MAX_STRETCH_HBE - 1); patch < numPatches; patch++) {
    patchBands = xOverQmf[patch + 1] - xOverQmf[patch];
    target = xOverQmf[patch];
    sourceBands = xOverQmf[MAX_STRETCH_HBE - 1] - xOverQmf[MAX_STRETCH_HBE - 2];

    while (patchBands > 0) {
      int numBands = sourceBands;
      int startBand = xOverQmf[MAX_STRETCH_HBE - 1] - 1;
      if (target + numBands >= xOverQmf[patch + 1]) {
        numBands = xOverQmf[patch + 1] - target;
      }
      if ((((target + numBands - 1) % 2) +
           ((xOverQmf[MAX_STRETCH_HBE - 1] - 1) % 2)) %
          2) {
        if (numBands == sourceBands) {
          numBands--;
        } else {
          startBand--;
        }
      }
      if (keepStatesSynced == KEEP_STATES_SYNCED_OUTDIFF) {
        for (col = slotOffset; col < overlap + LPC_ORDER; col++) {
          i = 0;
          for (band = numBands; band > 0; band--) {
            if ((target + band - 1 < 64) &&
                (target + band - 1 < xOverQmf[patch + 1])) {
              ppqmfReal[col][target + band - 1] = ppqmfReal[col][startBand - i];
              ppqmfImag[col][target + band - 1] = ppqmfImag[col][startBand - i];
              i++;
            }
          }
        }
      } else {
        for (col = slotOffset; col < noCols; col++) {
          i = 0;
          for (band = numBands; band > 0; band--) {
            if ((target + band - 1 < 64) &&
                (target + band - 1 < xOverQmf[patch + 1])) {
              ppqmfReal[col][target + band - 1] = ppqmfReal[col][startBand - i];
              ppqmfImag[col][target + band - 1] = ppqmfImag[col][startBand - i];
              i++;
            }
          }
        }
      }
      target += numBands;
      patchBands -= numBands;
    }
  }
}

void sbr_dec(
    HANDLE_SBR_DEC hSbrDec,
    LONG *timeIn,
    LONG *timeOut,
    HANDLE_SBR_DEC hSbrDecRight,
    LONG *timeOutRight,
    const int strideOut,
    HANDLE_SBR_HEADER_DATA hHeaderData,
    HANDLE_SBR_FRAME_DATA hFrameData,
    HANDLE_SBR_PREV_FRAME_DATA
        hPrevFrameData,
    const int applyProcessing,
    HANDLE_PS_DEC h_ps_d, const UINT flags, const int codecFrameSize,
    const INT sbrInDataHeadroom) {
  int i, slot, reserve;
  int saveLbScale;
  int lastSlotOffs;
  FIXP_DBL maxVal;

  LONG *pTimeInQmf = timeIn;

  int ov_len = hSbrDec->LppTrans.pSettings->overlap;
  int noCols = hHeaderData->numberTimeSlots * hHeaderData->timeStep;

  FIXP_DBL *pLowBandReal[(3 * 4) + 2 * ((1024) / (32) * (4) / 2)];
  FIXP_DBL *pLowBandImag[(3 * 4) + 2 * ((1024) / (32) * (4) / 2)];

  FIXP_DBL **pReal = pLowBandReal + ov_len;
  FIXP_DBL **pImag = pLowBandImag + ov_len;

  for (i = 0; i < noCols + ov_len; i++) {
    pLowBandReal[i] = hSbrDec->qmfDomainInCh->hQmfSlotsReal[i];
    pLowBandImag[i] = hSbrDec->qmfDomainInCh->hQmfSlotsImag[i];
  }

  if ((flags & SBRDEC_USAC_HARMONICSBR)) {
    for (i = 0; i < noCols; i++) {
      pLowBandReal[i + noCols + ov_len] = hSbrDec->hQmfHBESlotsReal[i];
      pLowBandImag[i + noCols + ov_len] = hSbrDec->hQmfHBESlotsImag[i];
    }

    hSbrDec->scale_ov = hSbrDec->scale_lb;
    hSbrDec->scale_lb = hSbrDec->scale_hbe;

    pReal += noCols;
    pImag += noCols;
    if (flags & SBRDEC_SKIP_QMF_ANA) {
      FDK_QmfDomain_QmfData2HBE(hSbrDec->qmfDomainInCh,
                                hSbrDec->hQmfHBESlotsReal,
                                hSbrDec->hQmfHBESlotsImag);
    } else {
      for (i = 0; i < noCols; i++) {
        FDKmemcpy(pLowBandReal[ov_len + i], hSbrDec->hQmfHBESlotsReal[i],
                  hHeaderData->numberOfAnalysisBands * sizeof(FIXP_DBL));
        FDKmemcpy(pLowBandImag[ov_len + i], hSbrDec->hQmfHBESlotsImag[i],
                  hHeaderData->numberOfAnalysisBands * sizeof(FIXP_DBL));
      }
    }
  }

  if (flags & SBRDEC_SKIP_QMF_ANA) {
    if (!(flags & SBRDEC_USAC_HARMONICSBR))
      FDK_QmfDomain_WorkBuffer2ProcChannel(hSbrDec->qmfDomainInCh);
  } else {
    C_AALLOC_SCRATCH_START(qmfTemp, FIXP_DBL, 2 * (64));
    qmfAnalysisFiltering(&hSbrDec->qmfDomainInCh->fb, pReal, pImag,
                         &hSbrDec->qmfDomainInCh->scaling, pTimeInQmf,
                         0 + sbrInDataHeadroom, 1, qmfTemp);

    C_AALLOC_SCRATCH_END(qmfTemp, FIXP_DBL, 2 * (64));
  }

  if (!((flags & SBRDEC_USAC_HARMONICSBR) &&
        (hFrameData->sbrPatchingMode == 0))) {
    int nAnalysisBands = hHeaderData->numberOfAnalysisBands;

    if (!(flags & SBRDEC_LOW_POWER)) {
      for (slot = ov_len; slot < noCols + ov_len; slot++) {
        FDKmemclear(&pLowBandReal[slot][nAnalysisBands],
                    ((64) - nAnalysisBands) * sizeof(FIXP_DBL));
        FDKmemclear(&pLowBandImag[slot][nAnalysisBands],
                    ((64) - nAnalysisBands) * sizeof(FIXP_DBL));
      }
    } else {
      for (slot = ov_len; slot < noCols + ov_len; slot++) {
        FDKmemclear(&pLowBandReal[slot][nAnalysisBands],
                    ((64) - nAnalysisBands) * sizeof(FIXP_DBL));
      }
    }
  }

  maxVal = maxSubbandSample(pReal, (flags & SBRDEC_LOW_POWER) ? NULL : pImag, 0,
                            hSbrDec->qmfDomainInCh->fb.no_channels, 0, noCols);

  reserve = fixMax(0, CntLeadingZeros(maxVal) - 1);
  reserve = fixMin(reserve,
                   DFRACT_BITS - 1 - hSbrDec->qmfDomainInCh->scaling.lb_scale);

  rescaleSubbandSamples(pReal, (flags & SBRDEC_LOW_POWER) ? NULL : pImag, 0,
                        hSbrDec->qmfDomainInCh->fb.no_channels, 0, noCols,
                        reserve);

  hSbrDec->qmfDomainInCh->scaling.lb_scale += reserve;

  if ((flags & SBRDEC_USAC_HARMONICSBR)) {
    hSbrDec->scale_hbe = hSbrDec->qmfDomainInCh->scaling.lb_scale;
    hSbrDec->qmfDomainInCh->scaling.lb_scale = hSbrDec->scale_lb;
  }
  saveLbScale = hSbrDec->qmfDomainInCh->scaling.lb_scale;

  if (applyProcessing) {
    UCHAR *borders = hFrameData->frameInfo.borders;
    lastSlotOffs = borders[hFrameData->frameInfo.nEnvelopes] -
                   hHeaderData->numberTimeSlots;

    FIXP_DBL degreeAlias[(64)];
    PVC_DYNAMIC_DATA pvcDynamicData;
    pvcInitFrame(
        &hSbrDec->PvcStaticData, &pvcDynamicData,
        (hHeaderData->frameErrorFlag ? 0 : hHeaderData->bs_info.pvc_mode),
        hFrameData->ns, hHeaderData->timeStep,
        hHeaderData->freqBandData.lowSubband,
        hFrameData->frameInfo.pvcBorders[0], hFrameData->pvcID);

    if (!hHeaderData->frameErrorFlag && (hHeaderData->bs_info.pvc_mode > 0)) {
      pvcDecodeFrame(&hSbrDec->PvcStaticData, &pvcDynamicData, pLowBandReal,
                     pLowBandImag, ov_len,
                     SCALE2EXP(hSbrDec->qmfDomainInCh->scaling.ov_lb_scale),
                     SCALE2EXP(hSbrDec->qmfDomainInCh->scaling.lb_scale));
    }
    pvcEndFrame(&hSbrDec->PvcStaticData, &pvcDynamicData);

    if (flags & SBRDEC_LOW_POWER)
      FDKmemclear(&degreeAlias[hHeaderData->freqBandData.lowSubband],
                  (hHeaderData->freqBandData.highSubband -
                   hHeaderData->freqBandData.lowSubband) *
                      sizeof(FIXP_DBL));

    {
      KEEP_STATES_SYNCED_MODE keepStatesSyncedMode =
          ((flags & SBRDEC_USAC_HARMONICSBR) &&
           (hFrameData->sbrPatchingMode != 0))
              ? KEEP_STATES_SYNCED_NORMAL
              : KEEP_STATES_SYNCED_OFF;

      if (flags & SBRDEC_USAC_HARMONICSBR) {
        if (flags & SBRDEC_QUAD_RATE) {
          pReal -= 32;
          pImag -= 32;
        }

        if ((hSbrDec->savedStates == 0) && (hFrameData->sbrPatchingMode == 1)) {
          for (i = 0; i < LPC_ORDER + ov_len; i++) {
            FDKmemcpy(
                hSbrDec->LppTrans.lpcFilterStatesRealLegSBR[i],
                hSbrDec->codecQMFBufferReal[noCols - LPC_ORDER - ov_len + i],
                hSbrDec->hHBE->noChannels * sizeof(FIXP_DBL));
            FDKmemcpy(
                hSbrDec->LppTrans.lpcFilterStatesImagLegSBR[i],
                hSbrDec->codecQMFBufferImag[noCols - LPC_ORDER - ov_len + i],
                hSbrDec->hHBE->noChannels * sizeof(FIXP_DBL));
          }
        }

        for (i = 0; i < hSbrDec->hHBE->noCols; i++) {
          FDKmemcpy(hSbrDec->codecQMFBufferReal[i], pLowBandReal[ov_len + i],
                    hSbrDec->hHBE->noChannels * sizeof(FIXP_DBL));
          FDKmemcpy(hSbrDec->codecQMFBufferImag[i], pLowBandImag[ov_len + i],
                    hSbrDec->hHBE->noChannels * sizeof(FIXP_DBL));
        }

        QmfTransposerApply(
            hSbrDec->hHBE, pReal, pImag, noCols, pLowBandReal, pLowBandImag,
            hSbrDec->LppTrans.lpcFilterStatesRealHBE,
            hSbrDec->LppTrans.lpcFilterStatesImagHBE,
            hFrameData->sbrPitchInBins, hSbrDec->scale_lb, hSbrDec->scale_hbe,
            &hSbrDec->qmfDomainInCh->scaling.hb_scale, hHeaderData->timeStep,
            borders[0], ov_len, keepStatesSyncedMode);

        if (flags & SBRDEC_QUAD_RATE) {
          int *xOverQmf = GetxOverBandQmfTransposer(hSbrDec->hHBE);

          copyHarmonicSpectrum(xOverQmf, pLowBandReal, pLowBandImag, noCols,
                               ov_len, keepStatesSyncedMode);
        }
      }
    }

    if ((flags & SBRDEC_USAC_HARMONICSBR) &&
        (hFrameData->sbrPatchingMode == 0)) {
      hSbrDec->prev_frame_lSbr = 0;
      hSbrDec->prev_frame_hbeSbr = 1;

      lppTransposerHBE(
          &hSbrDec->LppTrans, hSbrDec->hHBE, &hSbrDec->qmfDomainInCh->scaling,
          pLowBandReal, pLowBandImag, hHeaderData->timeStep, borders[0],
          lastSlotOffs, hHeaderData->freqBandData.nInvfBands,
          hFrameData->sbr_invf_mode, hPrevFrameData->sbr_invf_mode);

    } else {
      if (flags & SBRDEC_USAC_HARMONICSBR) {
        for (i = 0; i < LPC_ORDER + hSbrDec->LppTrans.pSettings->overlap; i++) {
          FDKmemcpy(hSbrDec->LppTrans.lpcFilterStatesRealHBE[i],
                    hSbrDec->qmfDomainInCh
                        ->hQmfSlotsReal[hSbrDec->hHBE->noCols - LPC_ORDER + i],
                    (64) * sizeof(FIXP_DBL));
          FDKmemcpy(hSbrDec->LppTrans.lpcFilterStatesImagHBE[i],
                    hSbrDec->qmfDomainInCh
                        ->hQmfSlotsImag[hSbrDec->hHBE->noCols - LPC_ORDER + i],
                    (64) * sizeof(FIXP_DBL));
        }
      }
      {
        hSbrDec->prev_frame_lSbr = 1;
        hSbrDec->prev_frame_hbeSbr = 0;
      }

      lppTransposer(
          &hSbrDec->LppTrans, &hSbrDec->qmfDomainInCh->scaling, pLowBandReal,
          degreeAlias,
          pLowBandImag, flags & SBRDEC_LOW_POWER,
          hHeaderData->bs_info.sbr_preprocessing,
          hHeaderData->freqBandData.v_k_master[0], hHeaderData->timeStep,
          borders[0], lastSlotOffs, hHeaderData->freqBandData.nInvfBands,
          hFrameData->sbr_invf_mode, hPrevFrameData->sbr_invf_mode);
    }

    if ((hFrameData->sbrPatchingMode !=
         hSbrDec->SbrCalculateEnvelope.sbrPatchingMode)) {
      ResetLimiterBands(hHeaderData->freqBandData.limiterBandTable,
                        &hHeaderData->freqBandData.noLimiterBands,
                        hHeaderData->freqBandData.freqBandTable[0],
                        hHeaderData->freqBandData.nSfb[0],
                        hSbrDec->LppTrans.pSettings->patchParam,
                        hSbrDec->LppTrans.pSettings->noOfPatches,
                        hHeaderData->bs_data.limiterBands,
                        hFrameData->sbrPatchingMode,
                        (flags & SBRDEC_USAC_HARMONICSBR) &&
                                (hFrameData->sbrPatchingMode == 0)
                            ? GetxOverBandQmfTransposer(hSbrDec->hHBE)
                            : NULL,
                        Get41SbrQmfTransposer(hSbrDec->hHBE));

      hSbrDec->SbrCalculateEnvelope.sbrPatchingMode =
          hFrameData->sbrPatchingMode;
    }

    calculateSbrEnvelope(
        &hSbrDec->qmfDomainInCh->scaling, &hSbrDec->SbrCalculateEnvelope,
        hHeaderData, hFrameData, &pvcDynamicData, pLowBandReal, pLowBandImag,
        flags & SBRDEC_LOW_POWER,

        degreeAlias, flags,
        (hHeaderData->frameErrorFlag || hPrevFrameData->frameErrorFlag));

#if (SBRDEC_MAX_HB_FADE_FRAMES > 0)
    if (hHeaderData->frameErrorFlag) {
      if (hSbrDec->highBandFadeCnt < SBRDEC_MAX_HB_FADE_FRAMES) {
        hSbrDec->highBandFadeCnt += 1;
      }
    } else {
      if (hSbrDec->highBandFadeCnt >
          0) {
        hSbrDec->qmfDomainInCh->scaling.hb_scale += hSbrDec->highBandFadeCnt;
        hSbrDec->qmfDomainInCh->scaling.hb_scale =
            fMin(hSbrDec->qmfDomainInCh->scaling.hb_scale, DFRACT_BITS - 1);
        hSbrDec->highBandFadeCnt -= 1;
      }
    }

#endif
    for (i = 0; i < hHeaderData->freqBandData.nInvfBands; i++) {
      hPrevFrameData->sbr_invf_mode[i] = hFrameData->sbr_invf_mode[i];
    }
    hPrevFrameData->coupling = hFrameData->coupling;
    hPrevFrameData->stopPos = borders[hFrameData->frameInfo.nEnvelopes];
    hPrevFrameData->ampRes = hFrameData->ampResolutionCurrentFrame;
    hPrevFrameData->prevSbrPitchInBins = hFrameData->sbrPitchInBins;
    FDKmemcpy(&hPrevFrameData->prevFrameInfo, &hFrameData->frameInfo,
              sizeof(FRAME_INFO));
  } else {
    int rescale;
    int lsb;
    int length;

    hSbrDec->qmfDomainInCh->scaling.hb_scale = saveLbScale;

    rescale = hSbrDec->qmfDomainInCh->scaling.hb_scale -
              hSbrDec->qmfDomainInCh->scaling.ov_lb_scale;
    lsb = hSbrDec->qmfDomainOutCh->fb.lsb;
    length = (hSbrDec->qmfDomainInCh->fb.no_channels - lsb);

    if ((rescale < 0) && (length > 0)) {
      if (!(flags & SBRDEC_LOW_POWER)) {
        for (i = 0; i < ov_len; i++) {
          scaleValues(&pLowBandReal[i][lsb], length, rescale);
          scaleValues(&pLowBandImag[i][lsb], length, rescale);
        }
      } else {
        for (i = 0; i < ov_len; i++) {
          scaleValues(&pLowBandReal[i][lsb], length, rescale);
        }
      }
    }
  }

  if (!(flags & SBRDEC_USAC_HARMONICSBR)) {
    int length = hSbrDec->qmfDomainInCh->fb.lsb;
    if (flags & SBRDEC_SYNTAX_USAC) {
      length = hSbrDec->qmfDomainInCh->fb.no_channels;
    }

    for (i = 0; i < LPC_ORDER + ov_len; i++) {
      if (!(flags & SBRDEC_LOW_POWER)) {
        FDKmemcpy(hSbrDec->LppTrans.lpcFilterStatesRealLegSBR[i],
                  pLowBandReal[noCols - LPC_ORDER + i],
                  length * sizeof(FIXP_DBL));
        FDKmemcpy(hSbrDec->LppTrans.lpcFilterStatesImagLegSBR[i],
                  pLowBandImag[noCols - LPC_ORDER + i],
                  length * sizeof(FIXP_DBL));
      } else
        FDKmemcpy(hSbrDec->LppTrans.lpcFilterStatesRealLegSBR[i],
                  pLowBandReal[noCols - LPC_ORDER + i],
                  length * sizeof(FIXP_DBL));
    }
  }

  if (!(flags & SBRDEC_PS_DECODED)) {
    if (!(flags & SBRDEC_SKIP_QMF_SYN)) {
      int outScalefactor = -(8);

      if (h_ps_d != NULL) {
        h_ps_d->procFrameBased = 1;
      }

      sbrDecoder_drcApply(&hSbrDec->sbrDrcChannel, pLowBandReal,
                          (flags & SBRDEC_LOW_POWER) ? NULL : pLowBandImag,
                          hSbrDec->qmfDomainOutCh->fb.no_col, &outScalefactor);

      qmfChangeOutScalefactor(&hSbrDec->qmfDomainOutCh->fb, outScalefactor);

      {
        HANDLE_FREQ_BAND_DATA hFreq = &hHeaderData->freqBandData;
        int save_usb = hSbrDec->qmfDomainOutCh->fb.usb;

#if (QMF_MAX_SYNTHESIS_BANDS <= 64)
        C_AALLOC_SCRATCH_START(qmfTemp, FIXP_DBL, 2 * QMF_MAX_SYNTHESIS_BANDS);
#else
        C_AALLOC_STACK_START(qmfTemp, FIXP_DBL, 2 * QMF_MAX_SYNTHESIS_BANDS);
#endif
        if (hSbrDec->qmfDomainOutCh->fb.usb < hFreq->ov_highSubband) {
          hSbrDec->qmfDomainOutCh->fb.usb =
              fMin((UINT)hFreq->ov_highSubband,
                   (UINT)hSbrDec->qmfDomainOutCh->fb.no_channels);
        }
        {
          qmfSynthesisFiltering(
              &hSbrDec->qmfDomainOutCh->fb, pLowBandReal,
              (flags & SBRDEC_LOW_POWER) ? NULL : pLowBandImag,
              &hSbrDec->qmfDomainInCh->scaling,
              hSbrDec->LppTrans.pSettings->overlap, timeOut, strideOut,
              qmfTemp);
        }
        hSbrDec->qmfDomainOutCh->fb.usb = save_usb;
        hFreq->ov_highSubband = save_usb;
#if (QMF_MAX_SYNTHESIS_BANDS <= 64)
        C_AALLOC_SCRATCH_END(qmfTemp, FIXP_DBL, 2 * QMF_MAX_SYNTHESIS_BANDS);
#else
        C_AALLOC_STACK_END(qmfTemp, FIXP_DBL, 2 * QMF_MAX_SYNTHESIS_BANDS);
#endif
      }
    }

  } else {
    INT sdiff;
    INT scaleFactorHighBand, scaleFactorLowBand_ov, scaleFactorLowBand_no_ov,
        outScalefactor, outScalefactorR, outScalefactorL;

    HANDLE_QMF_FILTER_BANK synQmf = &hSbrDec->qmfDomainOutCh->fb;
    HANDLE_QMF_FILTER_BANK synQmfRight = &hSbrDecRight->qmfDomainOutCh->fb;

    sdiff = hSbrDec->qmfDomainInCh->scaling.lb_scale -
            reserve;
    scaleFactorHighBand = sdiff - hSbrDec->qmfDomainInCh->scaling.hb_scale;
    scaleFactorLowBand_ov = sdiff - hSbrDec->qmfDomainInCh->scaling.ov_lb_scale;
    scaleFactorLowBand_no_ov = sdiff - hSbrDec->qmfDomainInCh->scaling.lb_scale;

    scaleFactorLowBand_ov =
        fMin(DFRACT_BITS - 1, fMax(-(DFRACT_BITS - 1), scaleFactorLowBand_ov));
    scaleFactorLowBand_no_ov = fMin(
        DFRACT_BITS - 1, fMax(-(DFRACT_BITS - 1), scaleFactorLowBand_no_ov));
    scaleFactorHighBand =
        fMin(DFRACT_BITS - 1, fMax(-(DFRACT_BITS - 1), scaleFactorHighBand));

    if (h_ps_d->procFrameBased == 1)
    {
      FDK_ASSERT(hSbrDec->qmfDomainInCh->pGlobalConf->nBandsSynthesis <=
                 QMF_MAX_SYNTHESIS_BANDS);
      synQmfRight->outScalefactor = synQmf->outScalefactor;
      FDKmemcpy(synQmfRight->FilterStates, synQmf->FilterStates,
                9 * hSbrDec->qmfDomainInCh->pGlobalConf->nBandsSynthesis *
                    sizeof(FIXP_QSS));
    }

    PreparePsProcessing(h_ps_d, pLowBandReal, pLowBandImag,
                        scaleFactorLowBand_ov);

    synQmfRight->no_col = synQmf->no_col;
    synQmfRight->lsb = synQmf->lsb;
    synQmfRight->usb = synQmf->usb;

    int env = 0;

    {
#if (QMF_MAX_SYNTHESIS_BANDS <= 64)
      C_AALLOC_SCRATCH_START(pWorkBuffer, FIXP_DBL,
                             2 * QMF_MAX_SYNTHESIS_BANDS);
#else
      C_AALLOC_STACK_START(pWorkBuffer, FIXP_DBL, 2 * QMF_MAX_SYNTHESIS_BANDS);
#endif

      int maxShift = 0;

      if (hSbrDec->sbrDrcChannel.enable != 0) {
        if (hSbrDec->sbrDrcChannel.prevFact_exp > maxShift) {
          maxShift = hSbrDec->sbrDrcChannel.prevFact_exp;
        }
        if (hSbrDec->sbrDrcChannel.currFact_exp > maxShift) {
          maxShift = hSbrDec->sbrDrcChannel.currFact_exp;
        }
        if (hSbrDec->sbrDrcChannel.nextFact_exp > maxShift) {
          maxShift = hSbrDec->sbrDrcChannel.nextFact_exp;
        }
      }

      FDKmemcpy(&hSbrDecRight->sbrDrcChannel, &hSbrDec->sbrDrcChannel,
                sizeof(SBRDEC_DRC_CHANNEL));

      outScalefactor = maxShift - (8);
      outScalefactorL = outScalefactorR =
          sbrInDataHeadroom + 1;

      for (i = 0; i < synQmf->no_col; i++) {
        FIXP_DBL *rQmfReal = pWorkBuffer;
        FIXP_DBL *rQmfImag = pWorkBuffer + synQmf->no_channels;

        {
          if (i ==
              h_ps_d->bsData[h_ps_d->processSlot].mpeg.aEnvStartStop[env]) {
            initSlotBasedRotation(h_ps_d, env,
                                  hHeaderData->freqBandData.highSubband);
            env++;
          }

          ApplyPsSlot(
              h_ps_d,
              (pLowBandReal + i),
              (pLowBandImag + i),
              rQmfReal,
              rQmfImag,
              scaleFactorLowBand_no_ov,
              (i < hSbrDec->LppTrans.pSettings->overlap)
                  ? scaleFactorLowBand_ov
                  : scaleFactorLowBand_no_ov,
              scaleFactorHighBand, synQmf->lsb, synQmf->usb);
        }

        sbrDecoder_drcApplySlot(
                                &hSbrDecRight->sbrDrcChannel, rQmfReal,
                                rQmfImag, i, synQmfRight->no_col, maxShift);

        sbrDecoder_drcApplySlot(
                                &hSbrDec->sbrDrcChannel, *(pLowBandReal + i),
                                *(pLowBandImag + i), i, synQmf->no_col,
                                maxShift);

        if (!(flags & SBRDEC_SKIP_QMF_SYN)) {
          qmfChangeOutScalefactor(synQmf, outScalefactor);
          qmfChangeOutScalefactor(synQmfRight, outScalefactor);

          qmfSynthesisFilteringSlot(
              synQmfRight, rQmfReal,
              rQmfImag,
              outScalefactorL, outScalefactorL,
              timeOutRight + (i * synQmf->no_channels * strideOut), strideOut,
              pWorkBuffer);

          qmfSynthesisFilteringSlot(
              synQmf, *(pLowBandReal + i),
              *(pLowBandImag + i),
              outScalefactorR, outScalefactorR,
              timeOut + (i * synQmf->no_channels * strideOut), strideOut,
              pWorkBuffer);
        }
      }
#if (QMF_MAX_SYNTHESIS_BANDS <= 64)
      C_AALLOC_SCRATCH_END(pWorkBuffer, FIXP_DBL, 2 * QMF_MAX_SYNTHESIS_BANDS);
#else
      C_AALLOC_STACK_END(pWorkBuffer, FIXP_DBL, 2 * QMF_MAX_SYNTHESIS_BANDS);
#endif
    }
  }

  sbrDecoder_drcUpdateChannel(&hSbrDec->sbrDrcChannel);

  if (!(flags & SBRDEC_SKIP_QMF_SYN)) {
    {
      FDK_QmfDomain_SaveOverlap(hSbrDec->qmfDomainInCh, 0);
      FDK_ASSERT(hSbrDec->qmfDomainInCh->scaling.ov_lb_scale == saveLbScale);
    }
  }

  hSbrDec->savedStates = 0;

  hPrevFrameData->frameErrorFlag = hHeaderData->frameErrorFlag;
  hSbrDec->applySbrProc_old = applyProcessing;
}

SBR_ERROR
createSbrDec(SBR_CHANNEL *hSbrChannel,
             HANDLE_SBR_HEADER_DATA hHeaderData,
             TRANSPOSER_SETTINGS *pSettings,
             const int downsampleFac,
             const UINT qmfFlags,
             const UINT flags, const int overlap,
             int chan,
             int codecFrameSize)

{
  SBR_ERROR err = SBRDEC_OK;
  int timeSlots =
      hHeaderData->numberTimeSlots;
  int noCols =
      timeSlots * hHeaderData->timeStep;
  HANDLE_SBR_DEC hs = &(hSbrChannel->SbrDec);

#if (SBRDEC_MAX_HB_FADE_FRAMES > 0)
  hs->highBandFadeCnt = SBRDEC_MAX_HB_FADE_FRAMES;

#endif
  hs->scale_hbe = 15;
  hs->scale_lb = 15;
  hs->scale_ov = 15;

  hs->prev_frame_lSbr = 0;
  hs->prev_frame_hbeSbr = 0;

  hs->codecFrameSize = codecFrameSize;

  err = createSbrEnvelopeCalc(&hs->SbrCalculateEnvelope, hHeaderData, chan,
                              flags);
  if (err != SBRDEC_OK) {
    return err;
  }

  initSbrPrevFrameData(&hSbrChannel->prevFrameData, timeSlots);

  err = createLppTransposer(
      &hs->LppTrans, pSettings, hHeaderData->freqBandData.lowSubband,
      hHeaderData->freqBandData.v_k_master, hHeaderData->freqBandData.numMaster,
      hHeaderData->freqBandData.highSubband, timeSlots, noCols,
      hHeaderData->freqBandData.freqBandTableNoise,
      hHeaderData->freqBandData.nNfb, hHeaderData->sbrProcSmplRate, chan,
      overlap);
  if (err != SBRDEC_OK) {
    return err;
  }

  if (flags & SBRDEC_USAC_HARMONICSBR) {
    int noChannels, bSbr41 = flags & SBRDEC_QUAD_RATE ? 1 : 0;

    noChannels =
        QMF_SYNTH_CHANNELS /
        ((bSbr41 + 1) * 2);

    hSbrChannel->SbrDec.tmp_memory = (FIXP_DBL **)fdkCallocMatrix2D_aligned(
        noCols, noChannels, sizeof(FIXP_DBL));
    if (hSbrChannel->SbrDec.tmp_memory == NULL) {
      return SBRDEC_MEM_ALLOC_FAILED;
    }

    hSbrChannel->SbrDec.hQmfHBESlotsReal = hSbrChannel->SbrDec.tmp_memory;
    hSbrChannel->SbrDec.hQmfHBESlotsImag =
        (FIXP_DBL **)fdkCallocMatrix2D_aligned(noCols, noChannels,
                                               sizeof(FIXP_DBL));
    if (hSbrChannel->SbrDec.hQmfHBESlotsImag == NULL) {
      return SBRDEC_MEM_ALLOC_FAILED;
    }

    hSbrChannel->SbrDec.codecQMFBufferReal =
        (FIXP_DBL **)fdkCallocMatrix2D_aligned(noCols, noChannels,
                                               sizeof(FIXP_DBL));
    if (hSbrChannel->SbrDec.codecQMFBufferReal == NULL) {
      return SBRDEC_MEM_ALLOC_FAILED;
    }

    hSbrChannel->SbrDec.codecQMFBufferImag =
        (FIXP_DBL **)fdkCallocMatrix2D_aligned(noCols, noChannels,
                                               sizeof(FIXP_DBL));
    if (hSbrChannel->SbrDec.codecQMFBufferImag == NULL) {
      return SBRDEC_MEM_ALLOC_FAILED;
    }

    err = QmfTransposerCreate(&hs->hHBE, codecFrameSize, 0, bSbr41);
    if (err != SBRDEC_OK) {
      return err;
    }
  }

  return err;
}

int deleteSbrDec(SBR_CHANNEL *hSbrChannel) {
  HANDLE_SBR_DEC hs = &hSbrChannel->SbrDec;

  deleteSbrEnvelopeCalc(&hs->SbrCalculateEnvelope);

  if (hs->tmp_memory != NULL) {
    FDK_FREE_MEMORY_2D_ALIGNED(hs->tmp_memory);
  }

  FDK_FREE_MEMORY_2D_ALIGNED(hs->hQmfHBESlotsImag);

  if (hs->hHBE != NULL) QmfTransposerClose(hs->hHBE);

  if (hs->codecQMFBufferReal != NULL) {
    FDK_FREE_MEMORY_2D_ALIGNED(hs->codecQMFBufferReal);
  }

  if (hs->codecQMFBufferImag != NULL) {
    FDK_FREE_MEMORY_2D_ALIGNED(hs->codecQMFBufferImag);
  }

  return 0;
}

SBR_ERROR
resetSbrDec(HANDLE_SBR_DEC hSbrDec, HANDLE_SBR_HEADER_DATA hHeaderData,
            HANDLE_SBR_PREV_FRAME_DATA hPrevFrameData, const int downsampleFac,
            const UINT flags, HANDLE_SBR_FRAME_DATA hFrameData) {
  SBR_ERROR sbrError = SBRDEC_OK;
  int i;
  FIXP_DBL *pLowBandReal[128];
  FIXP_DBL *pLowBandImag[128];
  int useLP = flags & SBRDEC_LOW_POWER;

  int old_lsb = hSbrDec->qmfDomainInCh->fb.lsb;
  int old_usb = hSbrDec->qmfDomainInCh->fb.usb;
  int new_lsb = hHeaderData->freqBandData.lowSubband;
  int l, startBand, stopBand, startSlot, size;

  FIXP_DBL **OverlapBufferReal = hSbrDec->qmfDomainInCh->hQmfSlotsReal;
  FIXP_DBL **OverlapBufferImag = hSbrDec->qmfDomainInCh->hQmfSlotsImag;

  int applySbrProc = (hHeaderData->syncState == SBR_ACTIVE ||
                      (hHeaderData->frameErrorFlag == 0 &&
                       hHeaderData->syncState == SBR_HEADER));
  int applySbrProc_old = hSbrDec->applySbrProc_old;

  if (!applySbrProc) {
    new_lsb = (hSbrDec->qmfDomainInCh->fb).no_channels;
  }
  if (!applySbrProc_old) {
    old_lsb = (hSbrDec->qmfDomainInCh->fb).no_channels;
    old_usb = old_lsb;
  }

  resetSbrEnvelopeCalc(&hSbrDec->SbrCalculateEnvelope);

  FDK_ASSERT(hSbrDec->qmfDomainOutCh != NULL);
  hSbrDec->qmfDomainOutCh->fb.lsb =
      fixMin((INT)hSbrDec->qmfDomainOutCh->fb.no_channels,
             (INT)hHeaderData->freqBandData.lowSubband);
  hSbrDec->qmfDomainOutCh->fb.usb =
      fixMin((INT)hSbrDec->qmfDomainOutCh->fb.no_channels,
             (INT)hHeaderData->freqBandData.highSubband);
  /* Analysis */
  FDK_ASSERT(hSbrDec->qmfDomainInCh != NULL);
  hSbrDec->qmfDomainInCh->fb.lsb = hSbrDec->qmfDomainOutCh->fb.lsb;
  hSbrDec->qmfDomainInCh->fb.usb = hSbrDec->qmfDomainOutCh->fb.usb;

  startBand = old_lsb;
  stopBand = new_lsb;
  startSlot = fMax(0, hHeaderData->timeStep * (hPrevFrameData->stopPos -
                                               hHeaderData->numberTimeSlots));
  size = fMax(0, stopBand - startBand);

  if (!(flags & SBRDEC_SYNTAX_USAC)) {
    if (!useLP) {
      for (l = startSlot; l < hSbrDec->LppTrans.pSettings->overlap; l++) {
        FDKmemclear(&OverlapBufferReal[l][startBand], size * sizeof(FIXP_DBL));
        FDKmemclear(&OverlapBufferImag[l][startBand], size * sizeof(FIXP_DBL));
      }
    } else {
      for (l = startSlot; l < hSbrDec->LppTrans.pSettings->overlap; l++) {
        FDKmemclear(&OverlapBufferReal[l][startBand], size * sizeof(FIXP_DBL));
      }
    }

    startBand = fixMin(old_lsb, new_lsb);
    stopBand = fixMax(old_lsb, new_lsb);
    size = fixMax(0, stopBand - startBand);

    FDKmemclear(&hSbrDec->LppTrans.lpcFilterStatesRealLegSBR[0][startBand],
                size * sizeof(FIXP_DBL));
    FDKmemclear(&hSbrDec->LppTrans.lpcFilterStatesRealLegSBR[1][startBand],
                size * sizeof(FIXP_DBL));
    if (!useLP) {
      FDKmemclear(&hSbrDec->LppTrans.lpcFilterStatesImagLegSBR[0][startBand],
                  size * sizeof(FIXP_DBL));
      FDKmemclear(&hSbrDec->LppTrans.lpcFilterStatesImagLegSBR[1][startBand],
                  size * sizeof(FIXP_DBL));
    }
  }

  if (startSlot != 0) {
    int source_exp, target_exp, delta_exp, target_lsb, target_usb, reserve;
    FIXP_DBL maxVal;

    if (new_lsb > old_lsb) {
      source_exp = SCALE2EXP(hSbrDec->qmfDomainInCh->scaling.ov_hb_scale);
      target_exp = SCALE2EXP(hSbrDec->qmfDomainInCh->scaling.ov_lb_scale);

      startBand = old_lsb;

      if (new_lsb >= old_usb) {
        stopBand = old_usb;
      } else {
        stopBand = new_lsb;
      }

      target_lsb = 0;
      target_usb = old_lsb;
    } else {
      source_exp = SCALE2EXP(hSbrDec->qmfDomainInCh->scaling.ov_lb_scale);
      target_exp = SCALE2EXP(hSbrDec->qmfDomainInCh->scaling.ov_hb_scale);

      startBand = new_lsb;
      stopBand = old_lsb;

      target_lsb = old_lsb;
      target_usb = old_usb;
    }

    maxVal =
        maxSubbandSample(OverlapBufferReal, (useLP) ? NULL : OverlapBufferImag,
                         startBand, stopBand, 0, startSlot);

    reserve = ((LONG)maxVal != 0 ? CntLeadingZeros(maxVal) - 1 : 0);
    reserve = fixMin(
        reserve,
        DFRACT_BITS - 1 -
            EXP2SCALE(
                source_exp));

    if (target_exp - (source_exp - reserve) >= 0) {
      rescaleSubbandSamples(OverlapBufferReal,
                            (useLP) ? NULL : OverlapBufferImag, startBand,
                            stopBand, 0, startSlot, reserve);
      source_exp -= reserve;
    }

    delta_exp = target_exp - source_exp;

    if (delta_exp < 0) {
      startBand = target_lsb;
      stopBand = target_usb;
      delta_exp = -delta_exp;

      if (new_lsb > old_lsb) {
        hSbrDec->qmfDomainInCh->scaling.ov_lb_scale = EXP2SCALE(source_exp);
      } else {
        hSbrDec->qmfDomainInCh->scaling.ov_hb_scale = EXP2SCALE(source_exp);
      }
    }

    FDK_ASSERT(startBand <= stopBand);

    if (!useLP) {
      for (l = 0; l < startSlot; l++) {
        scaleValues(OverlapBufferReal[l] + startBand, stopBand - startBand,
                    -delta_exp);
        scaleValues(OverlapBufferImag[l] + startBand, stopBand - startBand,
                    -delta_exp);
      }
    } else
      for (l = 0; l < startSlot; l++) {
        scaleValues(OverlapBufferReal[l] + startBand, stopBand - startBand,
                    -delta_exp);
      }
  }

  sbrError = resetLppTransposer(
      &hSbrDec->LppTrans, hHeaderData->freqBandData.lowSubband,
      hHeaderData->freqBandData.v_k_master, hHeaderData->freqBandData.numMaster,
      hHeaderData->freqBandData.freqBandTableNoise,
      hHeaderData->freqBandData.nNfb, hHeaderData->freqBandData.highSubband,
      hHeaderData->sbrProcSmplRate);
  if (sbrError != SBRDEC_OK) return sbrError;

  hSbrDec->savedStates = 0;

  if ((flags & SBRDEC_USAC_HARMONICSBR) && applySbrProc) {
    sbrError = QmfTransposerReInit(hSbrDec->hHBE,
                                   hHeaderData->freqBandData.freqBandTable,
                                   hHeaderData->freqBandData.nSfb);
    if (sbrError != SBRDEC_OK) return sbrError;

    for (i = 0; i < LPC_ORDER + hSbrDec->LppTrans.pSettings->overlap; i++) {
      FDKmemcpy(
          hSbrDec->LppTrans.lpcFilterStatesRealLegSBR[i],
          hSbrDec->codecQMFBufferReal[hSbrDec->hHBE->noCols - LPC_ORDER -
                                      hSbrDec->LppTrans.pSettings->overlap + i],
          hSbrDec->hHBE->noChannels * sizeof(FIXP_DBL));
      FDKmemcpy(
          hSbrDec->LppTrans.lpcFilterStatesImagLegSBR[i],
          hSbrDec->codecQMFBufferImag[hSbrDec->hHBE->noCols - LPC_ORDER -
                                      hSbrDec->LppTrans.pSettings->overlap + i],
          hSbrDec->hHBE->noChannels * sizeof(FIXP_DBL));
    }
    hSbrDec->savedStates = 1;

    {
      for (i = 0; i < hSbrDec->LppTrans.pSettings->overlap + LPC_ORDER; i++) {
        pLowBandReal[i] = hSbrDec->LppTrans.lpcFilterStatesRealHBE[i];
        pLowBandImag[i] = hSbrDec->LppTrans.lpcFilterStatesImagHBE[i];
      }

      for (i = 0; i < hSbrDec->hHBE->noCols; i++) {
        pLowBandReal[i + hSbrDec->LppTrans.pSettings->overlap + LPC_ORDER] =
            hSbrDec->codecQMFBufferReal[i];
        pLowBandImag[i + hSbrDec->LppTrans.pSettings->overlap + LPC_ORDER] =
            hSbrDec->codecQMFBufferImag[i];
      }

      if (flags & SBRDEC_QUAD_RATE) {
        if (hFrameData->sbrPatchingMode == 0) {
          int *xOverQmf = GetxOverBandQmfTransposer(hSbrDec->hHBE);

          for (i = 0; i < hSbrDec->hHBE->noCols / 2; i++) {
            pLowBandReal[i + hSbrDec->hHBE->noCols +
                         hSbrDec->LppTrans.pSettings->overlap + LPC_ORDER] =
                hSbrDec->hQmfHBESlotsReal[i];
            pLowBandImag[i + hSbrDec->hHBE->noCols +
                         hSbrDec->LppTrans.pSettings->overlap + LPC_ORDER] =
                hSbrDec->hQmfHBESlotsImag[i];
          }

          QmfTransposerApply(
              hSbrDec->hHBE,
              pLowBandReal + hSbrDec->LppTrans.pSettings->overlap +
                  hSbrDec->hHBE->noCols / 2 + LPC_ORDER,
              pLowBandImag + hSbrDec->LppTrans.pSettings->overlap +
                  hSbrDec->hHBE->noCols / 2 + LPC_ORDER,
              hSbrDec->hHBE->noCols, pLowBandReal, pLowBandImag,
              hSbrDec->LppTrans.lpcFilterStatesRealHBE,
              hSbrDec->LppTrans.lpcFilterStatesImagHBE,
              hPrevFrameData->prevSbrPitchInBins, hSbrDec->scale_lb,
              hSbrDec->scale_hbe, &hSbrDec->qmfDomainInCh->scaling.hb_scale,
              hHeaderData->timeStep, hFrameData->frameInfo.borders[0],
              hSbrDec->LppTrans.pSettings->overlap, KEEP_STATES_SYNCED_OUTDIFF);

          copyHarmonicSpectrum(
              xOverQmf, pLowBandReal, pLowBandImag, hSbrDec->hHBE->noCols,
              hSbrDec->LppTrans.pSettings->overlap, KEEP_STATES_SYNCED_OUTDIFF);
        }
      } else {
        for (i = 0; i < hSbrDec->hHBE->noCols; i++) {
          pLowBandReal[i + hSbrDec->hHBE->noCols +
                       hSbrDec->LppTrans.pSettings->overlap + LPC_ORDER] =
              hSbrDec->hQmfHBESlotsReal[i];
          pLowBandImag[i + hSbrDec->hHBE->noCols +
                       hSbrDec->LppTrans.pSettings->overlap + LPC_ORDER] =
              hSbrDec->hQmfHBESlotsImag[i];
        }

        if (hFrameData->sbrPatchingMode == 0) {
          QmfTransposerApply(
              hSbrDec->hHBE,
              pLowBandReal + hSbrDec->LppTrans.pSettings->overlap + LPC_ORDER,
              pLowBandImag + hSbrDec->LppTrans.pSettings->overlap + LPC_ORDER,
              hSbrDec->hHBE->noCols, pLowBandReal, pLowBandImag,
              hSbrDec->LppTrans.lpcFilterStatesRealHBE,
              hSbrDec->LppTrans.lpcFilterStatesImagHBE,
              0,
              hSbrDec->scale_lb, hSbrDec->scale_lb,
              &hSbrDec->qmfDomainInCh->scaling.hb_scale, hHeaderData->timeStep,
              hFrameData->frameInfo.borders[0],
              hSbrDec->LppTrans.pSettings->overlap, KEEP_STATES_SYNCED_NOOUT);
        }

        QmfTransposerApply(
            hSbrDec->hHBE,
            pLowBandReal + hSbrDec->LppTrans.pSettings->overlap +
                hSbrDec->hHBE->noCols + LPC_ORDER,
            pLowBandImag + hSbrDec->LppTrans.pSettings->overlap +
                hSbrDec->hHBE->noCols + LPC_ORDER,
            hSbrDec->hHBE->noCols, pLowBandReal, pLowBandImag,
            hSbrDec->LppTrans.lpcFilterStatesRealHBE,
            hSbrDec->LppTrans.lpcFilterStatesImagHBE,
            hPrevFrameData->prevSbrPitchInBins, hSbrDec->scale_lb,
            hSbrDec->scale_hbe, &hSbrDec->qmfDomainInCh->scaling.hb_scale,
            hHeaderData->timeStep, hFrameData->frameInfo.borders[0],
            hSbrDec->LppTrans.pSettings->overlap, KEEP_STATES_SYNCED_OUTDIFF);
      }

      if (hFrameData->sbrPatchingMode == 0) {
        for (i = startSlot; i < hSbrDec->LppTrans.pSettings->overlap; i++) {
          FDKmemcpy(hSbrDec->qmfDomainInCh->hQmfSlotsReal[i],
                    hSbrDec->LppTrans.lpcFilterStatesRealHBE[i + LPC_ORDER],
                    (64) * sizeof(FIXP_DBL));
          FDKmemcpy(hSbrDec->qmfDomainInCh->hQmfSlotsImag[i],
                    hSbrDec->LppTrans.lpcFilterStatesImagHBE[i + LPC_ORDER],
                    (64) * sizeof(FIXP_DBL));
        }

        for (i = startSlot; i < hSbrDec->LppTrans.pSettings->overlap; i++) {
          FDKmemcpy(
              hSbrDec->qmfDomainInCh->hQmfSlotsReal[i],
              hSbrDec->codecQMFBufferReal[hSbrDec->hHBE->noCols -
                                          hSbrDec->LppTrans.pSettings->overlap +
                                          i],
              new_lsb * sizeof(FIXP_DBL));
          FDKmemcpy(
              hSbrDec->qmfDomainInCh->hQmfSlotsImag[i],
              hSbrDec->codecQMFBufferImag[hSbrDec->hHBE->noCols -
                                          hSbrDec->LppTrans.pSettings->overlap +
                                          i],
              new_lsb * sizeof(FIXP_DBL));
        }
      }
    }
  }

  {
    int adapt_lb = 0, diff = 0,
        new_scale = hSbrDec->qmfDomainInCh->scaling.ov_lb_scale;

    if ((hSbrDec->qmfDomainInCh->scaling.ov_lb_scale !=
         hSbrDec->qmfDomainInCh->scaling.lb_scale) &&
        startSlot != 0) {
      diff = SCALE2EXP(hSbrDec->qmfDomainInCh->scaling.ov_lb_scale) -
             SCALE2EXP(hSbrDec->qmfDomainInCh->scaling.lb_scale);

      if (diff > 0) {
        adapt_lb = 1;
        diff = -diff;
        new_scale = hSbrDec->qmfDomainInCh->scaling.ov_lb_scale;
      }

      stopBand = new_lsb;
    }

    if (hFrameData->sbrPatchingMode == 1) {
      for (i = 0; i < hSbrDec->LppTrans.pSettings->overlap + LPC_ORDER; i++) {
        scaleValues(hSbrDec->LppTrans.lpcFilterStatesRealLegSBR[i], new_lsb,
                    diff);
        if (!useLP) {
          scaleValues(hSbrDec->LppTrans.lpcFilterStatesImagLegSBR[i], new_lsb,
                      diff);
        }
      }

      if (flags & SBRDEC_SYNTAX_USAC) {
        for (i = startSlot; i < hSbrDec->LppTrans.pSettings->overlap; i++) {
          FDKmemcpy(&OverlapBufferReal[i][old_lsb],
                    &hSbrDec->LppTrans
                         .lpcFilterStatesRealLegSBR[LPC_ORDER + i][old_lsb],
                    fMax(new_lsb - old_lsb, 0) * sizeof(FIXP_DBL));
          if (!useLP) {
            FDKmemcpy(&OverlapBufferImag[i][old_lsb],
                      &hSbrDec->LppTrans
                           .lpcFilterStatesImagLegSBR[LPC_ORDER + i][old_lsb],
                      fMax(new_lsb - old_lsb, 0) * sizeof(FIXP_DBL));
          }
        }
      }

      if (new_lsb > old_lsb) {
        stopBand = old_lsb;
      }
    }
    if ((adapt_lb == 1) && (stopBand > startBand)) {
      for (l = startSlot; l < hSbrDec->LppTrans.pSettings->overlap; l++) {
        scaleValues(OverlapBufferReal[l] + startBand, stopBand - startBand,
                    diff);
        if (!useLP) {
          scaleValues(OverlapBufferImag[l] + startBand, stopBand - startBand,
                      diff);
        }
      }
    }
    hSbrDec->qmfDomainInCh->scaling.ov_lb_scale = new_scale;
  }

  sbrError = ResetLimiterBands(hHeaderData->freqBandData.limiterBandTable,
                               &hHeaderData->freqBandData.noLimiterBands,
                               hHeaderData->freqBandData.freqBandTable[0],
                               hHeaderData->freqBandData.nSfb[0],
                               hSbrDec->LppTrans.pSettings->patchParam,
                               hSbrDec->LppTrans.pSettings->noOfPatches,
                               hHeaderData->bs_data.limiterBands,
                               hFrameData->sbrPatchingMode,
                               GetxOverBandQmfTransposer(hSbrDec->hHBE),
                               Get41SbrQmfTransposer(hSbrDec->hHBE));

  hSbrDec->SbrCalculateEnvelope.sbrPatchingMode = hFrameData->sbrPatchingMode;

  return sbrError;
}
