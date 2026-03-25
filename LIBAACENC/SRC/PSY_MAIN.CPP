/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Psychoaccoustic major function block

*******************************************************************************/

#include "psy_const.h"

#include "block_switch.h"
#include "transform.h"
#include "spreading.h"
#include "pre_echo_control.h"
#include "band_nrg.h"
#include "psy_configuration.h"
#include "psy_data.h"
#include "ms_stereo.h"
#include "interface.h"
#include "psy_main.h"
#include "grp_data.h"
#include "tns_func.h"
#include "pns_func.h"
#include "tonality.h"
#include "aacEnc_ram.h"
#include "intensity.h"

#define FADE_OUT_LEN 6
static const FIXP_DBL fadeOutFactor[FADE_OUT_LEN] = {
    1840644096, 1533870080, 1227096064, 920322048, 613548032, 306774016};

AAC_ENCODER_ERROR FDKaacEnc_PsyNew(PSY_INTERNAL **phpsy, const INT nElements,
                                   const INT nChannels, UCHAR *dynamic_RAM) {
  AAC_ENCODER_ERROR ErrorStatus;
  PSY_INTERNAL *hPsy;
  INT i;

  hPsy = GetRam_aacEnc_PsyInternal();
  *phpsy = hPsy;
  if (hPsy == NULL) {
    ErrorStatus = AAC_ENC_NO_MEMORY;
    goto bail;
  }

  for (i = 0; i < nElements; i++) {
    hPsy->psyElement[i] = GetRam_aacEnc_PsyElement(i);
    if (hPsy->psyElement[i] == NULL) {
      ErrorStatus = AAC_ENC_NO_MEMORY;
      goto bail;
    }
  }

  for (i = 0; i < nChannels; i++) {
    hPsy->pStaticChannels[i] = GetRam_aacEnc_PsyStatic(i);
    if (hPsy->pStaticChannels[i] == NULL) {
      ErrorStatus = AAC_ENC_NO_MEMORY;
      goto bail;
    }
    hPsy->pStaticChannels[i]->psyInputBuffer = GetRam_aacEnc_PsyInputBuffer(i);
    if (hPsy->pStaticChannels[i]->psyInputBuffer == NULL) {
      ErrorStatus = AAC_ENC_NO_MEMORY;
      goto bail;
    }
  }

  hPsy->psyDynamic = GetRam_aacEnc_PsyDynamic(0, dynamic_RAM);

  return AAC_ENC_OK;

bail:
  FDKaacEnc_PsyClose(phpsy, NULL);

  return ErrorStatus;
}

AAC_ENCODER_ERROR FDKaacEnc_PsyOutNew(PSY_OUT **phpsyOut, const INT nElements,
                                      const INT nChannels, const INT nSubFrames,
                                      UCHAR *dynamic_RAM) {
  AAC_ENCODER_ERROR ErrorStatus;
  int n, i;
  int elInc = 0, chInc = 0;

  for (n = 0; n < nSubFrames; n++) {
    phpsyOut[n] = GetRam_aacEnc_PsyOut(n);

    if (phpsyOut[n] == NULL) {
      ErrorStatus = AAC_ENC_NO_MEMORY;
      goto bail;
    }

    for (i = 0; i < nChannels; i++) {
      phpsyOut[n]->pPsyOutChannels[i] = GetRam_aacEnc_PsyOutChannel(chInc++);
      if (NULL == phpsyOut[n]->pPsyOutChannels[i]) {
        ErrorStatus = AAC_ENC_NO_MEMORY;
        goto bail;
      }
    }

    for (i = 0; i < nElements; i++) {
      phpsyOut[n]->psyOutElement[i] = GetRam_aacEnc_PsyOutElements(elInc++);
      if (phpsyOut[n]->psyOutElement[i] == NULL) {
        ErrorStatus = AAC_ENC_NO_MEMORY;
        goto bail;
      }
    }
  }

  return AAC_ENC_OK;

bail:
  FDKaacEnc_PsyClose(NULL, phpsyOut);
  return ErrorStatus;
}

AAC_ENCODER_ERROR FDKaacEnc_psyInitStates(PSY_INTERNAL *hPsy,
                                          PSY_STATIC *psyStatic,
                                          AUDIO_OBJECT_TYPE audioObjectType) {
  FDKmemclear(psyStatic->psyInputBuffer,
              MAX_INPUT_BUFFER_SIZE * sizeof(INT_PCM));

  FDKaacEnc_InitBlockSwitching(&psyStatic->blockSwitchingControl,
                               isLowDelay(audioObjectType));

  return AAC_ENC_OK;
}

AAC_ENCODER_ERROR FDKaacEnc_psyInit(PSY_INTERNAL *hPsy, PSY_OUT **phpsyOut,
                                    const INT nSubFrames,
                                    const INT nMaxChannels,
                                    const AUDIO_OBJECT_TYPE audioObjectType,
                                    CHANNEL_MAPPING *cm) {
  AAC_ENCODER_ERROR ErrorStatus = AAC_ENC_OK;
  int i, ch, n, chInc = 0, resetChannels = 3;

  if ((nMaxChannels > 2) && (cm->nChannels == 2)) {
    chInc = 1;
    FDKaacEnc_psyInitStates(hPsy, hPsy->pStaticChannels[0], audioObjectType);
  }

  if ((nMaxChannels == 2)) {
    resetChannels = 0;
  }

  for (i = 0; i < cm->nElements; i++) {
    for (ch = 0; ch < cm->elInfo[i].nChannelsInEl; ch++) {
      hPsy->psyElement[i]->psyStatic[ch] = hPsy->pStaticChannels[chInc];
      if (cm->elInfo[i].elType != ID_LFE) {
        if (chInc >= resetChannels) {
          FDKaacEnc_psyInitStates(hPsy, hPsy->psyElement[i]->psyStatic[ch],
                                  audioObjectType);
        }
        mdct_init(&(hPsy->psyElement[i]->psyStatic[ch]->mdctPers), NULL, 0);
        hPsy->psyElement[i]->psyStatic[ch]->isLFE = 0;
      } else {
        hPsy->psyElement[i]->psyStatic[ch]->isLFE = 1;
      }
      chInc++;
    }
  }

  for (n = 0; n < nSubFrames; n++) {
    chInc = 0;
    for (i = 0; i < cm->nElements; i++) {
      for (ch = 0; ch < cm->elInfo[i].nChannelsInEl; ch++) {
        phpsyOut[n]->psyOutElement[i]->psyOutChannel[ch] =
            phpsyOut[n]->pPsyOutChannels[chInc++];
      }
    }
  }

  return ErrorStatus;
}

AAC_ENCODER_ERROR FDKaacEnc_psyMainInit(
    PSY_INTERNAL *hPsy, AUDIO_OBJECT_TYPE audioObjectType, CHANNEL_MAPPING *cm,
    INT sampleRate, INT granuleLength, INT bitRate, INT tnsMask, INT bandwidth,
    INT usePns, INT useIS, INT useMS, UINT syntaxFlags, ULONG initFlags) {
  AAC_ENCODER_ERROR ErrorStatus;
  int i, ch;
  int channelsEff = cm->nChannelsEff;
  int tnsChannels = 0;
  FB_TYPE filterBank;

  switch (FDKaacEnc_GetMonoStereoMode(cm->encMode)) {
    case EL_MODE_MONO:
      tnsChannels = 1;
      break;
    case EL_MODE_STEREO:
      tnsChannels = 2;
      break;
    default:
      tnsChannels = 0;
  }

  switch (audioObjectType) {
    default:
      filterBank = FB_LC;
      break;
    case AOT_ER_AAC_LD:
      filterBank = FB_LD;
      break;
    case AOT_ER_AAC_ELD:
      filterBank = FB_ELD;
      break;
  }

  hPsy->granuleLength = granuleLength;

  ErrorStatus = FDKaacEnc_InitPsyConfiguration(
      bitRate / channelsEff, sampleRate, bandwidth, LONG_WINDOW,
      hPsy->granuleLength, useIS, useMS, &(hPsy->psyConf[0]), filterBank);
  if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;

  ErrorStatus = FDKaacEnc_InitTnsConfiguration(
      (bitRate * tnsChannels) / channelsEff, sampleRate, tnsChannels,
      LONG_WINDOW, hPsy->granuleLength, isLowDelay(audioObjectType),
      (syntaxFlags & AC_SBR_PRESENT) ? 1 : 0, &(hPsy->psyConf[0].tnsConf),
      &hPsy->psyConf[0], (INT)(tnsMask & 2), (INT)(tnsMask & 8));

  if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;

  if (granuleLength > 512) {
    ErrorStatus = FDKaacEnc_InitPsyConfiguration(
        bitRate / channelsEff, sampleRate, bandwidth, SHORT_WINDOW,
        hPsy->granuleLength, useIS, useMS, &hPsy->psyConf[1], filterBank);

    if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;

    ErrorStatus = FDKaacEnc_InitTnsConfiguration(
        (bitRate * tnsChannels) / channelsEff, sampleRate, tnsChannels,
        SHORT_WINDOW, hPsy->granuleLength, isLowDelay(audioObjectType),
        (syntaxFlags & AC_SBR_PRESENT) ? 1 : 0, &hPsy->psyConf[1].tnsConf,
        &hPsy->psyConf[1], (INT)(tnsMask & 1), (INT)(tnsMask & 4));

    if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;
  }

  for (i = 0; i < cm->nElements; i++) {
    for (ch = 0; ch < cm->elInfo[i].nChannelsInEl; ch++) {
      if (initFlags) {
        FDKaacEnc_psyInitStates(hPsy, hPsy->psyElement[i]->psyStatic[ch],
                                audioObjectType);
      }

      FDKaacEnc_InitPreEchoControl(
          hPsy->psyElement[i]->psyStatic[ch]->sfbThresholdnm1,
          &hPsy->psyElement[i]->psyStatic[ch]->calcPreEcho,
          hPsy->psyConf[0].sfbCnt, hPsy->psyConf[0].sfbPcmQuantThreshold,
          &hPsy->psyElement[i]->psyStatic[ch]->mdctScalenm1);
    }
  }

  ErrorStatus = FDKaacEnc_InitPnsConfiguration(
      &hPsy->psyConf[0].pnsConf, bitRate / channelsEff, sampleRate, usePns,
      hPsy->psyConf[0].sfbCnt, hPsy->psyConf[0].sfbOffset,
      cm->elInfo[0].nChannelsInEl, (hPsy->psyConf[0].filterbank == FB_LC));
  if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;

  if (granuleLength > 512) {
    ErrorStatus = FDKaacEnc_InitPnsConfiguration(
        &hPsy->psyConf[1].pnsConf, bitRate / channelsEff, sampleRate, usePns,
        hPsy->psyConf[1].sfbCnt, hPsy->psyConf[1].sfbOffset,
        cm->elInfo[1].nChannelsInEl, (hPsy->psyConf[1].filterbank == FB_LC));
    if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;
  }

  return ErrorStatus;
}

AAC_ENCODER_ERROR FDKaacEnc_psyMain(INT channels, PSY_ELEMENT *psyElement,
                                    PSY_DYNAMIC *psyDynamic,
                                    PSY_CONFIGURATION *psyConf,
                                    PSY_OUT_ELEMENT *RESTRICT psyOutElement,
                                    INT_PCM *pInput, const UINT inputBufSize,
                                    INT *chIdx, INT totalChannels) {
  const INT commonWindow = 1;
  INT maxSfbPerGroup[(2)];
  INT mdctSpectrum_e;
  INT ch;
  INT w;
  INT sfb;
  INT line;

  PSY_CONFIGURATION *RESTRICT hPsyConfLong = &psyConf[0];
  PSY_CONFIGURATION *RESTRICT hPsyConfShort = &psyConf[1];
  PSY_OUT_CHANNEL **RESTRICT psyOutChannel = psyOutElement->psyOutChannel;
  FIXP_SGL sfbTonality[(2)][MAX_SFB_LONG];

  PSY_STATIC **RESTRICT psyStatic = psyElement->psyStatic;

  PSY_DATA *RESTRICT psyData[(2)];
  TNS_DATA *RESTRICT tnsData[(2)];
  PNS_DATA *RESTRICT pnsData[(2)];

  INT zeroSpec = TRUE;

  INT blockSwitchingOffset;

  PSY_CONFIGURATION *RESTRICT hThisPsyConf[(2)];
  INT windowLength[(2)];
  INT nWindows[(2)];
  INT wOffset;

  INT maxSfb[(2)];
  INT *pSfbMaxScaleSpec[(2)];
  FIXP_DBL *pSfbEnergy[(2)];
  FIXP_DBL *pSfbSpreadEnergy[(2)];
  FIXP_DBL *pSfbEnergyLdData[(2)];
  FIXP_DBL *pSfbEnergyMS[(2)];
  FIXP_DBL *pSfbThreshold[(2)];

  INT isShortWindow[(2)];

  const INT nTimeSamples = psyConf->granuleLength;

  switch (hPsyConfLong->filterbank) {
    case FB_LC:
      blockSwitchingOffset =
          nTimeSamples + (9 * nTimeSamples / (2 * TRANS_FAC));
      break;
    case FB_LD:
    case FB_ELD:
      blockSwitchingOffset = nTimeSamples;
      break;
    default:
      return AAC_ENC_UNSUPPORTED_FILTERBANK;
  }

  for (ch = 0; ch < channels; ch++) {
    psyData[ch] = &psyDynamic->psyData[ch];
    tnsData[ch] = &psyDynamic->tnsData[ch];
    pnsData[ch] = &psyDynamic->pnsData[ch];

    psyData[ch]->mdctSpectrum = psyOutChannel[ch]->mdctSpectrum;
  }

  if (hPsyConfLong->filterbank != FB_ELD) {
    int err;

    for (ch = 0; ch < channels; ch++) {
      C_ALLOC_SCRATCH_START(pTimeSignal, INT_PCM, (1024))

      FDKmemcpy(pTimeSignal, pInput + chIdx[ch] * inputBufSize,
                nTimeSamples * sizeof(INT_PCM));

      FDKaacEnc_BlockSwitching(&psyStatic[ch]->blockSwitchingControl,
                               nTimeSamples, psyStatic[ch]->isLFE, pTimeSignal);

      FDKmemcpy(psyStatic[ch]->psyInputBuffer + blockSwitchingOffset,
                pTimeSignal,
                (2 * nTimeSamples - blockSwitchingOffset) * sizeof(INT_PCM));

      C_ALLOC_SCRATCH_END(pTimeSignal, INT_PCM, (1024))
    }

    err = FDKaacEnc_SyncBlockSwitching(
        &psyStatic[0]->blockSwitchingControl,
        (channels > 1) ? &psyStatic[1]->blockSwitchingControl : NULL, channels,
        commonWindow);

    if (err) {
      return AAC_ENC_UNSUPPORTED_AOT;
    }

  } else {
    for (ch = 0; ch < channels; ch++) {
      FDKmemcpy(psyStatic[ch]->psyInputBuffer + blockSwitchingOffset,
                pInput + chIdx[ch] * inputBufSize,
                nTimeSamples * sizeof(INT_PCM));
    }
  }

  for (ch = 0; ch < channels; ch++)
    isShortWindow[ch] =
        (psyStatic[ch]->blockSwitchingControl.lastWindowSequence ==
         SHORT_WINDOW);

  for (ch = 0; ch < channels; ch++) {
    if (isShortWindow[ch]) {
      hThisPsyConf[ch] = hPsyConfShort;
      windowLength[ch] = psyConf->granuleLength / TRANS_FAC;
      nWindows[ch] = TRANS_FAC;
      maxSfb[ch] = MAX_SFB_SHORT;

      pSfbMaxScaleSpec[ch] = psyData[ch]->sfbMaxScaleSpec.Short[0];
      pSfbEnergy[ch] = psyData[ch]->sfbEnergy.Short[0];
      pSfbSpreadEnergy[ch] = psyData[ch]->sfbSpreadEnergy.Short[0];
      pSfbEnergyLdData[ch] = psyData[ch]->sfbEnergyLdData.Short[0];
      pSfbEnergyMS[ch] = psyData[ch]->sfbEnergyMS.Short[0];
      pSfbThreshold[ch] = psyData[ch]->sfbThreshold.Short[0];

    } else {
      hThisPsyConf[ch] = hPsyConfLong;
      windowLength[ch] = psyConf->granuleLength;
      nWindows[ch] = 1;
      maxSfb[ch] = MAX_GROUPED_SFB;

      pSfbMaxScaleSpec[ch] = psyData[ch]->sfbMaxScaleSpec.Long;
      pSfbEnergy[ch] = psyData[ch]->sfbEnergy.Long;
      pSfbSpreadEnergy[ch] = psyData[ch]->sfbSpreadEnergy.Long;
      pSfbEnergyLdData[ch] = psyData[ch]->sfbEnergyLdData.Long;
      pSfbEnergyMS[ch] = psyData[ch]->sfbEnergyMS.Long;
      pSfbThreshold[ch] = psyData[ch]->sfbThreshold.Long;
    }
  }

  for (ch = 0; ch < channels; ch++) {
    if (psyStatic[ch]->isLFE) {
      psyData[ch]->sfbActive = hThisPsyConf[ch]->sfbActiveLFE;
      psyData[ch]->lowpassLine = hThisPsyConf[ch]->lowpassLineLFE;
    } else {
      psyData[ch]->sfbActive = hThisPsyConf[ch]->sfbActive;
      psyData[ch]->lowpassLine = hThisPsyConf[ch]->lowpassLine;
    }

    if (hThisPsyConf[ch]->filterbank == FB_ELD) {
      if (FDKaacEnc_Transform_Real_Eld(
              psyStatic[ch]->psyInputBuffer, psyData[ch]->mdctSpectrum,
              psyStatic[ch]->blockSwitchingControl.lastWindowSequence,
              psyStatic[ch]->blockSwitchingControl.windowShape,
              &psyStatic[ch]->blockSwitchingControl.lastWindowShape,
              nTimeSamples, &mdctSpectrum_e, hThisPsyConf[ch]->filterbank,
              psyStatic[ch]->overlapAddBuffer) != 0) {
        return AAC_ENC_UNSUPPORTED_FILTERBANK;
      }
    } else {
      if (FDKaacEnc_Transform_Real(
              psyStatic[ch]->psyInputBuffer, psyData[ch]->mdctSpectrum,
              psyStatic[ch]->blockSwitchingControl.lastWindowSequence,
              psyStatic[ch]->blockSwitchingControl.windowShape,
              &psyStatic[ch]->blockSwitchingControl.lastWindowShape,
              &psyStatic[ch]->mdctPers, nTimeSamples, &mdctSpectrum_e,
              hThisPsyConf[ch]->filterbank) != 0) {
        return AAC_ENC_UNSUPPORTED_FILTERBANK;
      }
    }

    for (w = 0; w < nWindows[ch]; w++) {
      wOffset = w * windowLength[ch];

      FDKmemclear(
          &psyData[ch]->mdctSpectrum[psyData[ch]->lowpassLine + wOffset],
          (windowLength[ch] - psyData[ch]->lowpassLine) * sizeof(FIXP_DBL));

      if ((hPsyConfLong->filterbank != FB_LC) &&
          (psyData[ch]->lowpassLine >= FADE_OUT_LEN)) {
        for (int i = 0; i < FADE_OUT_LEN; i++) {
          psyData[ch]->mdctSpectrum[psyData[ch]->lowpassLine + wOffset -
                                    FADE_OUT_LEN + i] =
              fMult(psyData[ch]->mdctSpectrum[psyData[ch]->lowpassLine +
                                              wOffset - FADE_OUT_LEN + i],
                    fadeOutFactor[i]);
        }
      }

      for (line = 0; (line < psyData[ch]->lowpassLine) && (zeroSpec == TRUE);
           line++) {
        if (psyData[ch]->mdctSpectrum[line + wOffset] != (FIXP_DBL)0) {
          zeroSpec = FALSE;
          break;
        }
      }
    }

    psyData[ch]->mdctScale = mdctSpectrum_e;

    FDKmemmove(psyStatic[ch]->psyInputBuffer,
               psyStatic[ch]->psyInputBuffer + nTimeSamples,
               nTimeSamples * sizeof(INT_PCM));

    FDKmemcpy(psyStatic[ch]->psyInputBuffer + nTimeSamples,
              pInput + (2 * nTimeSamples - blockSwitchingOffset) +
                  chIdx[ch] * inputBufSize,
              (blockSwitchingOffset - nTimeSamples) * sizeof(INT_PCM));
  }

  if (zeroSpec == FALSE) {
    INT minSpecShift = MAX_SHIFT_DBL;
    INT nrgShift = MAX_SHIFT_DBL;
    INT finalShift = MAX_SHIFT_DBL;
    FIXP_DBL currNrg = 0;
    FIXP_DBL maxNrg = 0;

    for (ch = 0; ch < channels; ch++) {
      for (w = 0; w < nWindows[ch]; w++) {
        wOffset = w * windowLength[ch];
        FDKaacEnc_CalcSfbMaxScaleSpec(
            psyData[ch]->mdctSpectrum + wOffset, hThisPsyConf[ch]->sfbOffset,
            pSfbMaxScaleSpec[ch] + w * maxSfb[ch], psyData[ch]->sfbActive);

        for (sfb = 0; sfb < psyData[ch]->sfbActive; sfb++)
          minSpecShift = fixMin(minSpecShift,
                                (pSfbMaxScaleSpec[ch] + w * maxSfb[ch])[sfb]);
      }
    }

    for (ch = 0; ch < channels; ch++) {
      for (w = 0; w < nWindows[ch]; w++) {
        wOffset = w * windowLength[ch];
        currNrg = FDKaacEnc_CheckBandEnergyOptim(
            psyData[ch]->mdctSpectrum + wOffset,
            pSfbMaxScaleSpec[ch] + w * maxSfb[ch], hThisPsyConf[ch]->sfbOffset,
            psyData[ch]->sfbActive, pSfbEnergy[ch] + w * maxSfb[ch],
            pSfbEnergyLdData[ch] + w * maxSfb[ch], minSpecShift - 4);

        maxNrg = fixMax(maxNrg, currNrg);
      }
    }

    if (maxNrg != (FIXP_DBL)0) {
      nrgShift = (CountLeadingBits(maxNrg) >> 1) + (minSpecShift - 4);
    }

    if (isShortWindow[0]) nrgShift--;

    finalShift = fixMin(minSpecShift, nrgShift);

    if (finalShift > psyData[0]->mdctScale + 3)
      finalShift = psyData[0]->mdctScale + 3;

    FDK_ASSERT(finalShift >= 0);

    FIXP_DBL ldShift = finalShift * FL2FXCONST_DBL(2.0 / 64);
    for (ch = 0; ch < channels; ch++) {
      INT maxSfb_ch = maxSfb[ch];
      INT w_maxSfb_ch = 0;
      for (w = 0; w < nWindows[ch]; w++) {
        for (sfb = 0; sfb < psyData[ch]->sfbActive; sfb++) {
          INT scale = fixMax(0, (pSfbMaxScaleSpec[ch] + w_maxSfb_ch)[sfb] - 4);
          scale = fixMin((scale - finalShift) << 1, DFRACT_BITS - 1);
          if (scale >= 0)
            (pSfbEnergy[ch] + w_maxSfb_ch)[sfb] >>= (scale);
          else
            (pSfbEnergy[ch] + w_maxSfb_ch)[sfb] <<= (-scale);
          (pSfbThreshold[ch] + w_maxSfb_ch)[sfb] =
              fMult((pSfbEnergy[ch] + w_maxSfb_ch)[sfb], C_RATIO);
          (pSfbEnergyLdData[ch] + w_maxSfb_ch)[sfb] += ldShift;
        }
        w_maxSfb_ch += maxSfb_ch;
      }
    }

    if (finalShift != 0) {
      for (ch = 0; ch < channels; ch++) {
        INT wLen = windowLength[ch];
        INT lowpassLine = psyData[ch]->lowpassLine;
        wOffset = 0;
        FIXP_DBL *mdctSpectrum = &psyData[ch]->mdctSpectrum[0];
        for (w = 0; w < nWindows[ch]; w++) {
          FIXP_DBL *spectrum = &mdctSpectrum[wOffset];
          for (line = 0; line < lowpassLine; line++) {
            spectrum[line] <<= finalShift;
          }
          wOffset += wLen;

          for (sfb = 0; sfb < psyData[ch]->sfbActive; sfb++)
            (pSfbMaxScaleSpec[ch] + w * maxSfb[ch])[sfb] -= finalShift;
        }
        psyData[ch]->mdctScale -= finalShift;
      }
    }

  } else {
    for (ch = 0; ch < channels; ch++) {
      psyData[ch]->mdctScale =
          0;
      for (w = 0; w < nWindows[ch]; w++) {
        for (sfb = 0; sfb < psyData[ch]->sfbActive; sfb++) {
          (pSfbMaxScaleSpec[ch] + w * maxSfb[ch])[sfb] = 0;
          (pSfbEnergy[ch] + w * maxSfb[ch])[sfb] = (FIXP_DBL)0;
          (pSfbEnergyLdData[ch] + w * maxSfb[ch])[sfb] = FL2FXCONST_DBL(-1.0f);
          (pSfbThreshold[ch] + w * maxSfb[ch])[sfb] = (FIXP_DBL)0;
        }
      }
    }
  }

  if ((channels >= 1) && (psyStatic[0]->isLFE)) {
    tnsData[0]->dataRaw.Long.subBlockInfo.tnsActive[HIFILT] = 0;
    tnsData[0]->dataRaw.Long.subBlockInfo.tnsActive[LOFILT] = 0;
  } else {
    for (ch = 0; ch < channels; ch++) {
      if (!isShortWindow[ch]) {
        FDKaacEnc_CalculateFullTonality(
            psyData[ch]->mdctSpectrum, pSfbMaxScaleSpec[ch],
            pSfbEnergyLdData[ch], sfbTonality[ch], psyData[ch]->sfbActive,
            hThisPsyConf[ch]->sfbOffset, hThisPsyConf[ch]->pnsConf.usePns);
      }
    }

    if (hPsyConfLong->tnsConf.tnsActive || hPsyConfShort->tnsConf.tnsActive) {
      INT tnsActive[TRANS_FAC] = {0};
      INT nrgScaling[2] = {0, 0};
      INT tnsSpecShift = 0;

      for (ch = 0; ch < channels; ch++) {
        for (w = 0; w < nWindows[ch]; w++) {
          wOffset = w * windowLength[ch];
          FDKaacEnc_TnsDetect(
              tnsData[ch], &hThisPsyConf[ch]->tnsConf,
              &psyOutChannel[ch]->tnsInfo, hThisPsyConf[ch]->sfbCnt,
              psyData[ch]->mdctSpectrum + wOffset, w,
              psyStatic[ch]->blockSwitchingControl.lastWindowSequence);
        }
      }

      if (channels == 2) {
        FDKaacEnc_TnsSync(
            tnsData[1], tnsData[0], &psyOutChannel[1]->tnsInfo,
            &psyOutChannel[0]->tnsInfo,

            psyStatic[1]->blockSwitchingControl.lastWindowSequence,
            psyStatic[0]->blockSwitchingControl.lastWindowSequence,
            &hThisPsyConf[1]->tnsConf);
      }

      if (channels >= 1) {
        FDK_ASSERT(1 == commonWindow);
        for (w = 0; w < nWindows[0]; w++) {
          if (isShortWindow[0])
            tnsActive[w] =
                tnsData[0]->dataRaw.Short.subBlockInfo[w].tnsActive[HIFILT] ||
                tnsData[0]->dataRaw.Short.subBlockInfo[w].tnsActive[LOFILT] ||
                tnsData[channels - 1]
                    ->dataRaw.Short.subBlockInfo[w]
                    .tnsActive[HIFILT] ||
                tnsData[channels - 1]
                    ->dataRaw.Short.subBlockInfo[w]
                    .tnsActive[LOFILT];
          else
            tnsActive[w] =
                tnsData[0]->dataRaw.Long.subBlockInfo.tnsActive[HIFILT] ||
                tnsData[0]->dataRaw.Long.subBlockInfo.tnsActive[LOFILT] ||
                tnsData[channels - 1]
                    ->dataRaw.Long.subBlockInfo.tnsActive[HIFILT] ||
                tnsData[channels - 1]
                    ->dataRaw.Long.subBlockInfo.tnsActive[LOFILT];
        }
      }

      for (ch = 0; ch < channels; ch++) {
        if (tnsActive[0] && !isShortWindow[ch]) {
          INT shift = 1;
          for (sfb = 0; sfb < hThisPsyConf[ch]->lowpassLine; sfb++) {
            psyData[ch]->mdctSpectrum[sfb] =
                psyData[ch]->mdctSpectrum[sfb] >> shift;
          }
          for (sfb = 0; sfb < psyData[ch]->sfbActive; sfb++) {
            pSfbThreshold[ch][sfb] >>= (2 * shift);
          }

          psyData[ch]->mdctScale += shift;
        }
      }

      for (ch = 0; ch < channels; ch++) {
        for (w = 0; w < nWindows[ch]; w++) {
          wOffset = w * windowLength[ch];
          FDKaacEnc_TnsEncode(
              &psyOutChannel[ch]->tnsInfo, tnsData[ch],
              hThisPsyConf[ch]->sfbCnt, &hThisPsyConf[ch]->tnsConf,
              hThisPsyConf[ch]->sfbOffset[psyData[ch]->sfbActive],
                  psyData[ch]->mdctSpectrum +
                  wOffset,
              w, psyStatic[ch]->blockSwitchingControl.lastWindowSequence);

          if (tnsActive[w]) {
            FDKaacEnc_CalcSfbMaxScaleSpec(psyData[ch]->mdctSpectrum + wOffset,
                                          hThisPsyConf[ch]->sfbOffset,
                                          pSfbMaxScaleSpec[ch] + w * maxSfb[ch],
                                          psyData[ch]->sfbActive);
          }
        }
      }

      for (ch = 0; ch < channels; ch++) {
        for (w = 0; w < nWindows[ch]; w++) {
          if (tnsActive[w]) {
            if (isShortWindow[ch]) {
              FDKaacEnc_CalcBandEnergyOptimShort(
                  psyData[ch]->mdctSpectrum + w * windowLength[ch],
                  pSfbMaxScaleSpec[ch] + w * maxSfb[ch],
                  hThisPsyConf[ch]->sfbOffset, psyData[ch]->sfbActive,
                  pSfbEnergy[ch] + w * maxSfb[ch]);
            } else {
              nrgScaling[ch] =
                  FDKaacEnc_CalcBandEnergyOptimLong(
                      psyData[ch]->mdctSpectrum, pSfbMaxScaleSpec[ch],
                      hThisPsyConf[ch]->sfbOffset, psyData[ch]->sfbActive,
                      pSfbEnergy[ch], pSfbEnergyLdData[ch]);
              tnsSpecShift =
                  fixMax(tnsSpecShift, nrgScaling[ch]);
            }
          }
        }
      }

      for (ch = 0; ch < channels; ch++) {
        if ((tnsSpecShift != 0) && !isShortWindow[ch]) {
          for (line = 0; line < hThisPsyConf[ch]->lowpassLine; line++) {
            psyData[ch]->mdctSpectrum[line] >>= tnsSpecShift;
          }
          INT scale = (tnsSpecShift - nrgScaling[ch]) << 1;
          for (sfb = 0; sfb < psyData[ch]->sfbActive; sfb++) {
            pSfbEnergyLdData[ch][sfb] -=
                scale * FL2FXCONST_DBL(1.0 / LD_DATA_SCALING);
            pSfbEnergy[ch][sfb] >>= scale;
            pSfbThreshold[ch][sfb] >>= (tnsSpecShift << 1);
          }
          psyData[ch]->mdctScale += tnsSpecShift;
        }
      }

    }
    else {
      FDKmemclear(psyDynamic->tnsData, sizeof(psyDynamic->tnsData));
    }
  }

  for (ch = 0; ch < channels; ch++) {
    INT headroom;

    FIXP_DBL clipEnergy;
    INT energyShift = psyData[ch]->mdctScale * 2;
    INT clipNrgShift = energyShift - THR_SHIFTBITS;
    if (isShortWindow[ch])
      headroom = 6;
    else
      headroom = 0;

    if (clipNrgShift >= 0)
      clipEnergy = hThisPsyConf[ch]->clipEnergy >> clipNrgShift;
    else if (clipNrgShift >= -headroom)
      clipEnergy = hThisPsyConf[ch]->clipEnergy << -clipNrgShift;
    else
      clipEnergy = (FIXP_DBL)MAXVAL_DBL;

    for (w = 0; w < nWindows[ch]; w++) {
      INT i;
      for (i = 0; i < psyData[ch]->sfbActive; i++) {
        *(pSfbThreshold[ch] + w * maxSfb[ch] + i) =
            fixMin(*(pSfbThreshold[ch] + w * maxSfb[ch] + i), clipEnergy);
      }

      FDKaacEnc_SpreadingMax(psyData[ch]->sfbActive,
                             hThisPsyConf[ch]->sfbMaskLowFactor,
                             hThisPsyConf[ch]->sfbMaskHighFactor,
                             pSfbThreshold[ch] + w * maxSfb[ch]);

      energyShift += PCM_QUANT_THR_SCALE;
      if (energyShift >= 0) {
        energyShift = fixMin(DFRACT_BITS - 1, energyShift);
        for (i = 0; i < psyData[ch]->sfbActive; i++) {
          *(pSfbThreshold[ch] + w * maxSfb[ch] + i) = fixMax(
              *(pSfbThreshold[ch] + w * maxSfb[ch] + i) >> THR_SHIFTBITS,
              (hThisPsyConf[ch]->sfbPcmQuantThreshold[i] >> energyShift));
        }
      } else {
        energyShift = fixMin(DFRACT_BITS - 1, -energyShift);
        for (i = 0; i < psyData[ch]->sfbActive; i++) {
          *(pSfbThreshold[ch] + w * maxSfb[ch] + i) = fixMax(
              *(pSfbThreshold[ch] + w * maxSfb[ch] + i) >> THR_SHIFTBITS,
              (hThisPsyConf[ch]->sfbPcmQuantThreshold[i] << energyShift));
        }
      }

      if (!psyStatic[ch]->isLFE) {
        if (psyStatic[ch]->blockSwitchingControl.lastWindowSequence ==
            STOP_WINDOW) {
          for (i = 0; i < psyData[ch]->sfbActive; i++) {
            psyStatic[ch]->sfbThresholdnm1[i] = (FIXP_DBL)MAXVAL_DBL;
          }

          psyStatic[ch]->mdctScalenm1 = 0;
          psyStatic[ch]->calcPreEcho = 0;
        }

        FDKaacEnc_PreEchoControl(
            psyStatic[ch]->sfbThresholdnm1, psyStatic[ch]->calcPreEcho,
            psyData[ch]->sfbActive, hThisPsyConf[ch]->maxAllowedIncreaseFactor,
            hThisPsyConf[ch]->minRemainingThresholdFactor,
            pSfbThreshold[ch] + w * maxSfb[ch], psyData[ch]->mdctScale,
            &psyStatic[ch]->mdctScalenm1);

        psyStatic[ch]->calcPreEcho = 1;

        if (psyStatic[ch]->blockSwitchingControl.lastWindowSequence ==
            START_WINDOW) {
          for (i = 0; i < psyData[ch]->sfbActive; i++) {
            psyStatic[ch]->sfbThresholdnm1[i] = (FIXP_DBL)MAXVAL_DBL;
          }

          psyStatic[ch]->mdctScalenm1 = 0;
          psyStatic[ch]->calcPreEcho = 0;
        }
      }

      FDKmemcpy(pSfbSpreadEnergy[ch] + w * maxSfb[ch],
                pSfbEnergy[ch] + w * maxSfb[ch],
                psyData[ch]->sfbActive * sizeof(FIXP_DBL));

      FDKaacEnc_SpreadingMax(psyData[ch]->sfbActive,
                             hThisPsyConf[ch]->sfbMaskLowFactorSprEn,
                             hThisPsyConf[ch]->sfbMaskHighFactorSprEn,
                             pSfbSpreadEnergy[ch] + w * maxSfb[ch]);
    }
  }

  if (channels == 2) {
    for (w = 0; w < nWindows[1]; w++) {
      wOffset = w * windowLength[1];
      FDKaacEnc_CalcBandNrgMSOpt(
          psyData[0]->mdctSpectrum + wOffset,
          psyData[1]->mdctSpectrum + wOffset,
          pSfbMaxScaleSpec[0] + w * maxSfb[0],
          pSfbMaxScaleSpec[1] + w * maxSfb[1], hThisPsyConf[1]->sfbOffset,
          psyData[0]->sfbActive, pSfbEnergyMS[0] + w * maxSfb[0],
          pSfbEnergyMS[1] + w * maxSfb[1],
          (psyStatic[1]->blockSwitchingControl.lastWindowSequence !=
           SHORT_WINDOW),
          psyData[0]->sfbEnergyMSLdData, psyData[1]->sfbEnergyMSLdData);
    }
  }

  for (ch = 0; ch < channels; ch++) {
    if (isShortWindow[ch]) {
      int sfbGrp;
      int noSfb = psyStatic[ch]->blockSwitchingControl.noOfGroups *
                  hPsyConfShort->sfbCnt;
      FDKaacEnc_groupShortData(
          psyData[ch]->mdctSpectrum, &psyData[ch]->sfbThreshold,
          &psyData[ch]->sfbEnergy, &psyData[ch]->sfbEnergyMS,
          &psyData[ch]->sfbSpreadEnergy, hPsyConfShort->sfbCnt,
          psyData[ch]->sfbActive, hPsyConfShort->sfbOffset,
          hPsyConfShort->sfbMinSnrLdData, psyData[ch]->groupedSfbOffset,
          &maxSfbPerGroup[ch], psyOutChannel[ch]->sfbMinSnrLdData,
          psyStatic[ch]->blockSwitchingControl.noOfGroups,
          psyStatic[ch]->blockSwitchingControl.groupLen,
          psyConf[1].granuleLength);

      for (sfbGrp = 0; sfbGrp < noSfb; sfbGrp += hPsyConfShort->sfbCnt) {
        LdDataVector(&psyData[ch]->sfbEnergy.Long[sfbGrp],
                     &psyOutChannel[ch]->sfbEnergyLdData[sfbGrp],
                     psyData[ch]->sfbActive);
      }

      for (sfbGrp = 0; sfbGrp < noSfb; sfbGrp += hPsyConfShort->sfbCnt) {
        LdDataVector(&psyData[ch]->sfbThreshold.Long[sfbGrp],
                     &psyOutChannel[ch]->sfbThresholdLdData[sfbGrp],
                     psyData[ch]->sfbActive);
        for (sfb = 0; sfb < psyData[ch]->sfbActive; sfb++) {
          psyOutChannel[ch]->sfbThresholdLdData[sfbGrp + sfb] =
              fixMax(psyOutChannel[ch]->sfbThresholdLdData[sfbGrp + sfb],
                     FL2FXCONST_DBL(-0.515625f));
        }
      }

      if (channels == 2) {
        for (sfbGrp = 0; sfbGrp < noSfb; sfbGrp += hPsyConfShort->sfbCnt) {
          LdDataVector(&psyData[ch]->sfbEnergyMS.Long[sfbGrp],
                       &psyData[ch]->sfbEnergyMSLdData[sfbGrp],
                       psyData[ch]->sfbActive);
        }
      }

      FDKmemcpy(psyOutChannel[ch]->sfbOffsets, psyData[ch]->groupedSfbOffset,
                (MAX_GROUPED_SFB + 1) * sizeof(INT));

    } else {
      int i;
      for (sfb = psyData[ch]->sfbActive - 1; sfb >= 0; sfb--) {
        for (line = hPsyConfLong->sfbOffset[sfb + 1] - 1;
             line >= hPsyConfLong->sfbOffset[sfb]; line--) {
          if (psyData[ch]->mdctSpectrum[line] != FL2FXCONST_SGL(0.0f)) break;
        }
        if (line > hPsyConfLong->sfbOffset[sfb]) break;
      }
      maxSfbPerGroup[ch] = sfb + 1;
      maxSfbPerGroup[ch] =
          fixMax(fixMin(5, psyData[ch]->sfbActive), maxSfbPerGroup[ch]);

      FDKmemcpy(psyOutChannel[ch]->sfbEnergyLdData,
                psyData[ch]->sfbEnergyLdData.Long,
                psyData[ch]->sfbActive * sizeof(FIXP_DBL));

      FDKmemcpy(psyOutChannel[ch]->sfbOffsets, hPsyConfLong->sfbOffset,
                (MAX_GROUPED_SFB + 1) * sizeof(INT));

      FDKmemcpy(psyOutChannel[ch]->sfbMinSnrLdData,
                hPsyConfLong->sfbMinSnrLdData,
                psyData[ch]->sfbActive * sizeof(FIXP_DBL));

      LdDataVector(psyData[ch]->sfbThreshold.Long,
                   psyOutChannel[ch]->sfbThresholdLdData,
                   psyData[ch]->sfbActive);
      for (i = 0; i < psyData[ch]->sfbActive; i++) {
        psyOutChannel[ch]->sfbThresholdLdData[i] =
            fixMax(psyOutChannel[ch]->sfbThresholdLdData[i],
                   FL2FXCONST_DBL(-0.515625f));
      }
    }
  }

  for (ch = 0; ch < channels; ch++) {
    FDKmemclear(psyOutChannel[ch]->isBook, MAX_GROUPED_SFB * sizeof(INT));
    FDKmemclear(psyOutChannel[ch]->isScale, MAX_GROUPED_SFB * sizeof(INT));
  }

  for (ch = 0; ch < channels; ch++) {
    INT win = (isShortWindow[ch] ? 1 : 0);
    if (!psyStatic[ch]->isLFE) {
      FDKaacEnc_PnsDetect(
          &(psyConf[0].pnsConf), pnsData[ch],
          psyStatic[ch]->blockSwitchingControl.lastWindowSequence,
          psyData[ch]->sfbActive,
          maxSfbPerGroup[ch],
          psyOutChannel[ch]->sfbThresholdLdData, psyConf[win].sfbOffset,
          psyData[ch]->mdctSpectrum, psyData[ch]->sfbMaxScaleSpec.Long,
          sfbTonality[ch], psyOutChannel[ch]->tnsInfo.order[0][0],
          tnsData[ch]->dataRaw.Long.subBlockInfo.predictionGain[HIFILT],
          tnsData[ch]->dataRaw.Long.subBlockInfo.tnsActive[HIFILT],
          psyOutChannel[ch]->sfbEnergyLdData, psyOutChannel[ch]->noiseNrg);
    }
  }

  if (channels == 2) {
    psyOutElement->toolsInfo.msDigest = MS_NONE;
    psyOutElement->commonWindow = commonWindow;
    if (psyOutElement->commonWindow)
      maxSfbPerGroup[0] = maxSfbPerGroup[1] =
          fixMax(maxSfbPerGroup[0], maxSfbPerGroup[1]);
    if (psyStatic[0]->blockSwitchingControl.lastWindowSequence !=
        SHORT_WINDOW) {
      FDKaacEnc_PreProcessPnsChannelPair(
          psyData[0]->sfbActive, (&psyData[0]->sfbEnergy)->Long,
          (&psyData[1]->sfbEnergy)->Long, psyOutChannel[0]->sfbEnergyLdData,
          psyOutChannel[1]->sfbEnergyLdData, psyData[0]->sfbEnergyMS.Long,
          &(psyConf[0].pnsConf), pnsData[0], pnsData[1]);

      FDKaacEnc_IntensityStereoProcessing(
          psyData[0]->sfbEnergy.Long, psyData[1]->sfbEnergy.Long,
          psyData[0]->mdctSpectrum, psyData[1]->mdctSpectrum,
          psyData[0]->sfbThreshold.Long, psyData[1]->sfbThreshold.Long,
          psyOutChannel[1]->sfbThresholdLdData,
          psyData[0]->sfbSpreadEnergy.Long, psyData[1]->sfbSpreadEnergy.Long,
          psyOutChannel[0]->sfbEnergyLdData, psyOutChannel[1]->sfbEnergyLdData,
          &psyOutElement->toolsInfo.msDigest, psyOutElement->toolsInfo.msMask,
          psyConf[0].sfbCnt, psyConf[0].sfbCnt, maxSfbPerGroup[0],
          psyConf[0].sfbOffset,
          psyConf[0].allowIS && psyOutElement->commonWindow,
          psyOutChannel[1]->isBook, psyOutChannel[1]->isScale, pnsData);

      FDKaacEnc_MsStereoProcessing(
          psyData, psyOutChannel, psyOutChannel[1]->isBook,
          &psyOutElement->toolsInfo.msDigest, psyOutElement->toolsInfo.msMask,
          psyConf[0].allowMS, psyData[0]->sfbActive, psyData[0]->sfbActive,
          maxSfbPerGroup[0], psyOutChannel[0]->sfbOffsets);

      FDKaacEnc_PostProcessPnsChannelPair(
          psyData[0]->sfbActive, &(psyConf[0].pnsConf), pnsData[0], pnsData[1],
          psyOutElement->toolsInfo.msMask, &psyOutElement->toolsInfo.msDigest);

    } else {
      FDKaacEnc_IntensityStereoProcessing(
          psyData[0]->sfbEnergy.Long, psyData[1]->sfbEnergy.Long,
          psyData[0]->mdctSpectrum, psyData[1]->mdctSpectrum,
          psyData[0]->sfbThreshold.Long, psyData[1]->sfbThreshold.Long,
          psyOutChannel[1]->sfbThresholdLdData,
          psyData[0]->sfbSpreadEnergy.Long, psyData[1]->sfbSpreadEnergy.Long,
          psyOutChannel[0]->sfbEnergyLdData, psyOutChannel[1]->sfbEnergyLdData,
          &psyOutElement->toolsInfo.msDigest, psyOutElement->toolsInfo.msMask,
          psyStatic[0]->blockSwitchingControl.noOfGroups *
              hPsyConfShort->sfbCnt,
          psyConf[1].sfbCnt, maxSfbPerGroup[0], psyData[0]->groupedSfbOffset,
          psyConf[0].allowIS && psyOutElement->commonWindow,
          psyOutChannel[1]->isBook, psyOutChannel[1]->isScale, pnsData);

      FDKaacEnc_MsStereoProcessing(
          psyData, psyOutChannel, psyOutChannel[1]->isBook,
          &psyOutElement->toolsInfo.msDigest, psyOutElement->toolsInfo.msMask,
          psyConf[1].allowMS,
          psyStatic[0]->blockSwitchingControl.noOfGroups *
              hPsyConfShort->sfbCnt,
          hPsyConfShort->sfbCnt, maxSfbPerGroup[0],
          psyOutChannel[0]->sfbOffsets);
    }
  }

  for (ch = 0; ch < channels; ch++) {
    if (psyStatic[ch]->isLFE) {
      for (sfb = 0; sfb < psyData[ch]->sfbActive; sfb++) {
        psyOutChannel[ch]->noiseNrg[sfb] = NO_NOISE_PNS;
      }
    } else {
      FDKaacEnc_CodePnsChannel(
          psyData[ch]->sfbActive, &(hThisPsyConf[ch]->pnsConf),
          pnsData[ch]->pnsFlag, psyData[ch]->sfbEnergyLdData.Long,
          psyOutChannel[ch]->noiseNrg,
          psyOutChannel[ch]->sfbThresholdLdData);
    }
  }

  for (ch = 0; ch < channels; ch++) {
    INT mask;
    int grp;
    psyOutChannel[ch]->maxSfbPerGroup = maxSfbPerGroup[ch];
    psyOutChannel[ch]->mdctScale = psyData[ch]->mdctScale;
    if (isShortWindow[ch] == 0) {
      psyOutChannel[ch]->sfbCnt = hPsyConfLong->sfbActive;
      psyOutChannel[ch]->sfbPerGroup = hPsyConfLong->sfbActive;
      psyOutChannel[ch]->lastWindowSequence =
          psyStatic[ch]->blockSwitchingControl.lastWindowSequence;
      psyOutChannel[ch]->windowShape =
          psyStatic[ch]->blockSwitchingControl.windowShape;
    } else {
      INT sfbCnt = psyStatic[ch]->blockSwitchingControl.noOfGroups *
                   hPsyConfShort->sfbCnt;

      psyOutChannel[ch]->sfbCnt = sfbCnt;
      psyOutChannel[ch]->sfbPerGroup = hPsyConfShort->sfbCnt;
      psyOutChannel[ch]->lastWindowSequence = SHORT_WINDOW;
      psyOutChannel[ch]->windowShape = SINE_WINDOW;
    }
    mask = 0;
    for (grp = 0; grp < psyStatic[ch]->blockSwitchingControl.noOfGroups;
         grp++) {
      int j;
      mask <<= 1;
      for (j = 1; j < psyStatic[ch]->blockSwitchingControl.groupLen[grp]; j++) {
        mask = (mask << 1) | 1;
      }
    }
    psyOutChannel[ch]->groupingMask = mask;

    FDKmemcpy(psyOutChannel[ch]->groupLen,
              psyStatic[ch]->blockSwitchingControl.groupLen,
              MAX_NO_OF_GROUPS * sizeof(INT));
    FDKmemcpy(psyOutChannel[ch]->sfbEnergy, (&psyData[ch]->sfbEnergy)->Long,
              MAX_GROUPED_SFB * sizeof(FIXP_DBL));
    FDKmemcpy(psyOutChannel[ch]->sfbSpreadEnergy,
              (&psyData[ch]->sfbSpreadEnergy)->Long,
              MAX_GROUPED_SFB * sizeof(FIXP_DBL));
  }

  return AAC_ENC_OK;
}

void FDKaacEnc_PsyClose(PSY_INTERNAL **phPsyInternal, PSY_OUT **phPsyOut) {
  int n, i;

  if (phPsyInternal != NULL) {
    PSY_INTERNAL *hPsyInternal = *phPsyInternal;

    if (hPsyInternal) {
      for (i = 0; i < (8); i++) {
        if (hPsyInternal->pStaticChannels[i]) {
          if (hPsyInternal->pStaticChannels[i]->psyInputBuffer)
            FreeRam_aacEnc_PsyInputBuffer(
                &hPsyInternal->pStaticChannels[i]
                     ->psyInputBuffer);

          FreeRam_aacEnc_PsyStatic(
              &hPsyInternal->pStaticChannels[i]);
        }
      }

      for (i = 0; i < ((8)); i++) {
        if (hPsyInternal->psyElement[i])
          FreeRam_aacEnc_PsyElement(
              &hPsyInternal->psyElement[i]);
      }

      FreeRam_aacEnc_PsyInternal(phPsyInternal);
    }
  }

  if (phPsyOut != NULL) {
    for (n = 0; n < (1); n++) {
      if (phPsyOut[n]) {
        for (i = 0; i < (8); i++) {
          if (phPsyOut[n]->pPsyOutChannels[i])
            FreeRam_aacEnc_PsyOutChannel(
                &phPsyOut[n]->pPsyOutChannels[i]);
        }

        for (i = 0; i < ((8)); i++) {
          if (phPsyOut[n]->psyOutElement[i])
            FreeRam_aacEnc_PsyOutElements(
                &phPsyOut[n]->psyOutElement[i]);
        }

        FreeRam_aacEnc_PsyOut(&phPsyOut[n]);
      }
    }
  }
}
