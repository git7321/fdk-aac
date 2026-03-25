/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description:

*******************************************************************************/

#include "channel.h"
#include "aacdecoder.h"
#include "block.h"
#include "aacdec_tns.h"
#include "FDK_bitstream.h"
#include "conceal.h"
#include "rvlc.h"
#include "aacdec_hcr.h"
#include "usacdec_lpd.h"
#include "usacdec_fac.h"

static void MapMidSideMaskToPnsCorrelation(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo[2]) {
  int group;

  for (group = 0; group < pAacDecoderChannelInfo[L]->icsInfo.WindowGroups;
       group++) {
    UCHAR groupMask = 1 << group;

    for (UCHAR band = 0; band < pAacDecoderChannelInfo[L]->icsInfo.MaxSfBands;
         band++) {
      if (pAacDecoderChannelInfo[L]->pComData->jointStereoData.MsUsed[band] &
          groupMask) {
        CPns_SetCorrelation(&pAacDecoderChannelInfo[L]->data.aac.PnsData, group,
                            band, 0);

        if (CPns_IsPnsUsed(&pAacDecoderChannelInfo[L]->data.aac.PnsData, group,
                           band) &&
            CPns_IsPnsUsed(&pAacDecoderChannelInfo[R]->data.aac.PnsData, group,
                           band))
          pAacDecoderChannelInfo[L]->pComData->jointStereoData.MsUsed[band] ^=
              groupMask;
      }
    }
  }
}

static void Clean_Complex_Prediction_coefficients(
    CJointStereoPersistentData *pJointStereoPersistentData, int windowGroups,
    const int low_limit, const int high_limit) {
  for (int group = 0; group < windowGroups; group++) {
    for (int sfb = low_limit; sfb < high_limit; sfb++) {
      pJointStereoPersistentData->alpha_q_re_prev[group][sfb] = 0;
      pJointStereoPersistentData->alpha_q_im_prev[group][sfb] = 0;
    }
  }
}

void CChannelElement_Decode(
    CAacDecoderChannelInfo
        *pAacDecoderChannelInfo[2],
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[2],
    SamplingRateInfo *pSamplingRateInfo, UINT flags, UINT elFlags,
    int el_channels) {
  int ch = 0;

  int maxSfBandsL = 0, maxSfBandsR = 0;
  int maybe_jstereo = (el_channels > 1);

  if (flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA) && el_channels == 2) {
    if (pAacDecoderChannelInfo[L]->data.usac.core_mode ||
        pAacDecoderChannelInfo[R]->data.usac.core_mode) {
      maybe_jstereo = 0;
    }
  }

  if (maybe_jstereo) {
    maxSfBandsL =
        GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo[L]->icsInfo);
    maxSfBandsR =
        GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo[R]->icsInfo);

    if (pAacDecoderChannelInfo[L]->pDynData->RawDataInfo.CommonWindow) {
      if (!(flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA))) {
        if (pAacDecoderChannelInfo[L]->data.aac.PnsData.PnsActive ||
            pAacDecoderChannelInfo[R]->data.aac.PnsData.PnsActive) {
          MapMidSideMaskToPnsCorrelation(pAacDecoderChannelInfo);
        }
      }
      if (((flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA)) &&
           (pAacDecoderChannelInfo[L]->pDynData->specificTo.usac.tns_on_lr ==
            1)) ||
          ((flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA)) == 0)) {
        int max_sfb_ste = (INT)(pAacDecoderChannelInfo[L]->icsInfo.max_sfb_ste);

        CJointStereo_ApplyMS(
            pAacDecoderChannelInfo, pAacDecoderStaticChannelInfo,
            pAacDecoderChannelInfo[L]->pSpectralCoefficient,
            pAacDecoderChannelInfo[R]->pSpectralCoefficient,
            pAacDecoderChannelInfo[L]->pDynData->aSfbScale,
            pAacDecoderChannelInfo[R]->pDynData->aSfbScale,
            pAacDecoderChannelInfo[L]->specScale,
            pAacDecoderChannelInfo[R]->specScale,
            GetScaleFactorBandOffsets(&pAacDecoderChannelInfo[L]->icsInfo,
                                      pSamplingRateInfo),
            GetWindowGroupLengthTable(&pAacDecoderChannelInfo[L]->icsInfo),
            GetWindowGroups(&pAacDecoderChannelInfo[L]->icsInfo), max_sfb_ste,
            maxSfBandsL, maxSfBandsR,
            pAacDecoderChannelInfo[L]
                ->pComData->jointStereoData.store_dmx_re_prev,
            &(pAacDecoderChannelInfo[L]
                  ->pComData->jointStereoData.store_dmx_re_prev_e),
            1);

      }
    }

    if (!(flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA))) {
      if ((pAacDecoderChannelInfo[L]->pDynData->RawDataInfo.CommonWindow ==
           1) &&
          (el_channels == 2)) {
        CJointStereo_ApplyIS(
            pAacDecoderChannelInfo,
            GetScaleFactorBandOffsets(&pAacDecoderChannelInfo[L]->icsInfo,
                                      pSamplingRateInfo),
            GetWindowGroupLengthTable(&pAacDecoderChannelInfo[L]->icsInfo),
            GetWindowGroups(&pAacDecoderChannelInfo[L]->icsInfo),
            GetScaleFactorBandsTransmitted(
                &pAacDecoderChannelInfo[L]->icsInfo));
      }
    }
  }

  for (ch = 0; ch < el_channels; ch++) {
    if (pAacDecoderChannelInfo[ch]->renderMode == AACDEC_RENDER_LPD) {
      CLpdChannelStream_Decode(pAacDecoderChannelInfo[ch],
                               pAacDecoderStaticChannelInfo[ch], flags);
    } else {
      UCHAR noSfbs =
          GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo[ch]->icsInfo);
      if ((maybe_jstereo == 1) &&
          (pAacDecoderChannelInfo[L]->pDynData->RawDataInfo.CommonWindow ==
           1)) {
        noSfbs = fMax(maxSfBandsL, maxSfBandsR);
      }
      int CP_active = 0;
      if (elFlags & AC_EL_USAC_CP_POSSIBLE) {
        CP_active = pAacDecoderChannelInfo[ch]
                        ->pComData->jointStereoData.cplx_pred_flag;
      }

      int max_sfb_ste = (INT)(pAacDecoderChannelInfo[L]->icsInfo.max_sfb_ste);

      if (!(CP_active && (max_sfb_ste == noSfbs)) ||
          !(CP_active &&
            !(pAacDecoderChannelInfo[ch]->pDynData->TnsData.Active)) ||
          ((flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA)) &&
           (pAacDecoderChannelInfo[L]->pDynData->specificTo.usac.tns_on_lr ==
            0))) {
        CBlock_ScaleSpectralData(pAacDecoderChannelInfo[ch], noSfbs,
                                 pSamplingRateInfo);

        if ((flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA)) &&
            (pAacDecoderChannelInfo[L]->pDynData->specificTo.usac.tns_on_lr ==
             0)) {
          if (IsLongBlock(&pAacDecoderChannelInfo[ch]->icsInfo)) {
            for (int i = 0; i < noSfbs; i++) {
              pAacDecoderChannelInfo[ch]->pDynData->aSfbScale[i] =
                  pAacDecoderChannelInfo[ch]->specScale[0];
            }
          } else {
            for (int i = 0; i < 8; i++) {
              for (int j = 0; j < noSfbs; j++) {
                pAacDecoderChannelInfo[ch]->pDynData->aSfbScale[i * 16 + j] =
                    pAacDecoderChannelInfo[ch]->specScale[i];
              }
            }
          }
        }
      }
    }
  }

  if (maybe_jstereo) {
    if (pAacDecoderChannelInfo[L]->pDynData->RawDataInfo.CommonWindow) {
    }
    else {
      if (elFlags & AC_EL_USAC_CP_POSSIBLE) {
        FDKmemclear(
            pAacDecoderStaticChannelInfo[L]
                ->pCpeStaticData->jointStereoPersistentData.alpha_q_re_prev,
            JointStereoMaximumGroups * JointStereoMaximumBands * sizeof(SHORT));
        FDKmemclear(
            pAacDecoderStaticChannelInfo[L]
                ->pCpeStaticData->jointStereoPersistentData.alpha_q_im_prev,
            JointStereoMaximumGroups * JointStereoMaximumBands * sizeof(SHORT));
      }
    }
  }

  for (ch = 0; ch < el_channels; ch++) {
    if (pAacDecoderChannelInfo[ch]->renderMode == AACDEC_RENDER_LPD) {
    } else {
      if (!(flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA))) {
        int pnsCh = (ch > 0) ? L : ch;
        CPns_UpdateNoiseState(
            &pAacDecoderChannelInfo[ch]->data.aac.PnsData,
            pAacDecoderChannelInfo[pnsCh]->data.aac.PnsData.currentSeed,
            pAacDecoderChannelInfo[ch]->pComData->pnsRandomSeed);
      }

      if ((!(flags & (AC_USAC))) ||
          ((flags & (AC_USAC)) &&
           (pAacDecoderChannelInfo[L]->pDynData->specificTo.usac.tns_active ==
            1)) ||
          (maybe_jstereo == 0)) {
        ApplyTools(
            pAacDecoderChannelInfo, pSamplingRateInfo, flags, elFlags, ch,
            pAacDecoderChannelInfo[L]->pDynData->RawDataInfo.CommonWindow);
      }
    }
  }

  if (maybe_jstereo) {
    if (pAacDecoderChannelInfo[L]->pDynData->RawDataInfo.CommonWindow) {
      if ((flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA)) &&
          (pAacDecoderChannelInfo[L]->pDynData->specificTo.usac.tns_on_lr ==
           0)) {
        int max_sfb_ste = (INT)(pAacDecoderChannelInfo[L]->icsInfo.max_sfb_ste);

        CJointStereo_ApplyMS(
            pAacDecoderChannelInfo, pAacDecoderStaticChannelInfo,
            pAacDecoderChannelInfo[L]->pSpectralCoefficient,
            pAacDecoderChannelInfo[R]->pSpectralCoefficient,
            pAacDecoderChannelInfo[L]->pDynData->aSfbScale,
            pAacDecoderChannelInfo[R]->pDynData->aSfbScale,
            pAacDecoderChannelInfo[L]->specScale,
            pAacDecoderChannelInfo[R]->specScale,
            GetScaleFactorBandOffsets(&pAacDecoderChannelInfo[L]->icsInfo,
                                      pSamplingRateInfo),
            GetWindowGroupLengthTable(&pAacDecoderChannelInfo[L]->icsInfo),
            GetWindowGroups(&pAacDecoderChannelInfo[L]->icsInfo), max_sfb_ste,
            maxSfBandsL, maxSfBandsR,
            pAacDecoderChannelInfo[L]
                ->pComData->jointStereoData.store_dmx_re_prev,
            &(pAacDecoderChannelInfo[L]
                  ->pComData->jointStereoData.store_dmx_re_prev_e),
            1);
      }
    }
  }

  for (ch = 0; ch < el_channels; ch++) {
    if (elFlags & AC_EL_USAC_CP_POSSIBLE) {
      pAacDecoderStaticChannelInfo[L]
          ->pCpeStaticData->jointStereoPersistentData.clearSpectralCoeffs = 0;
    }
  }

  CRvlc_ElementCheck(pAacDecoderChannelInfo, pAacDecoderStaticChannelInfo,
                     flags, el_channels);
}

void CChannel_CodebookTableInit(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo) {
  int b, w, maxBands, maxWindows;
  int maxSfb = GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo->icsInfo);
  UCHAR *pCodeBook = pAacDecoderChannelInfo->pDynData->aCodeBook;

  if (IsLongBlock(&pAacDecoderChannelInfo->icsInfo)) {
    maxBands = 64;
    maxWindows = 1;
  } else {
    maxBands = 16;
    maxWindows = 8;
  }

  for (w = 0; w < maxWindows; w++) {
    for (b = 0; b < maxSfb; b++) {
      pCodeBook[b] = ESCBOOK;
    }
    for (; b < maxBands; b++) {
      pCodeBook[b] = ZERO_HCB;
    }
    pCodeBook += maxBands;
  }
}

AAC_DECODER_ERROR CChannelElement_Read(
    HANDLE_FDK_BITSTREAM hBs, CAacDecoderChannelInfo *pAacDecoderChannelInfo[],
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[],
    const AUDIO_OBJECT_TYPE aot, SamplingRateInfo *pSamplingRateInfo,
    const UINT flags, const UINT elFlags, const UINT frame_length,
    const UCHAR numberOfChannels, const SCHAR epConfig,
    HANDLE_TRANSPORTDEC pTpDec) {
  AAC_DECODER_ERROR error = AAC_DEC_OK;
  const element_list_t *list;
  int i, ch, decision_bit;
  int crcReg1 = -1, crcReg2 = -1;
  int cplxPred;
  int ind_sw_cce_flag = 0, num_gain_element_lists = 0;

  FDK_ASSERT((numberOfChannels == 1) || (numberOfChannels == 2));

  list = getBitstreamElementList(aot, epConfig, numberOfChannels, 0, elFlags);
  if (list == NULL) {
    error = AAC_DEC_UNSUPPORTED_FORMAT;
    goto bail;
  }

  CTns_Reset(&pAacDecoderChannelInfo[0]->pDynData->TnsData);
  pAacDecoderChannelInfo[0]->pDynData->RawDataInfo.CommonWindow = 0;
  if (flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA)) {
    pAacDecoderChannelInfo[0]->pDynData->specificTo.usac.tns_active = 0;
    pAacDecoderChannelInfo[0]->pDynData->specificTo.usac.tns_on_lr = 0;
  }
  if (numberOfChannels == 2) {
    CTns_Reset(&pAacDecoderChannelInfo[1]->pDynData->TnsData);
    pAacDecoderChannelInfo[1]->pDynData->RawDataInfo.CommonWindow = 0;
  }

  cplxPred = 0;
  if (pAacDecoderStaticChannelInfo != NULL) {
    if (elFlags & AC_EL_USAC_CP_POSSIBLE) {
      pAacDecoderChannelInfo[0]->pComData->jointStereoData.cplx_pred_flag = 0;
      cplxPred = 1;
    }
  }

  if (0 || (flags & (AC_ELD | AC_SCALABLE))) {
    pAacDecoderChannelInfo[0]->pDynData->RawDataInfo.CommonWindow = 1;
    if (numberOfChannels == 2) {
      pAacDecoderChannelInfo[1]->pDynData->RawDataInfo.CommonWindow =
          pAacDecoderChannelInfo[0]->pDynData->RawDataInfo.CommonWindow;
    }
  }

  i = 0;
  ch = 0;
  decision_bit = 0;
  do {
    switch (list->id[i]) {
      case element_instance_tag:
        pAacDecoderChannelInfo[0]->ElementInstanceTag = FDKreadBits(hBs, 4);
        if (numberOfChannels == 2) {
          pAacDecoderChannelInfo[1]->ElementInstanceTag =
              pAacDecoderChannelInfo[0]->ElementInstanceTag;
        }
        break;
      case common_window:
        decision_bit =
            pAacDecoderChannelInfo[ch]->pDynData->RawDataInfo.CommonWindow =
                FDKreadBits(hBs, 1);
        if (numberOfChannels == 2) {
          pAacDecoderChannelInfo[1]->pDynData->RawDataInfo.CommonWindow =
              pAacDecoderChannelInfo[0]->pDynData->RawDataInfo.CommonWindow;
        }
        break;
      case ics_info:
        if (cplxPred) {
          if (pAacDecoderChannelInfo[0]->pDynData->RawDataInfo.CommonWindow) {
            pAacDecoderStaticChannelInfo[0]
                ->pCpeStaticData->jointStereoPersistentData.winSeqPrev =
                pAacDecoderChannelInfo[0]->icsInfo.WindowSequence;
            pAacDecoderStaticChannelInfo[0]
                ->pCpeStaticData->jointStereoPersistentData.winShapePrev =
                pAacDecoderChannelInfo[0]->icsInfo.WindowShape;
          }
        }
        error = IcsRead(hBs, &pAacDecoderChannelInfo[ch]->icsInfo,
                        pSamplingRateInfo, flags);

        if (elFlags & AC_EL_LFE &&
            GetWindowSequence(&pAacDecoderChannelInfo[ch]->icsInfo) !=
                BLOCK_LONG) {
          error = AAC_DEC_PARSE_ERROR;
          break;
        }

        if (numberOfChannels == 2 &&
            pAacDecoderChannelInfo[0]->pDynData->RawDataInfo.CommonWindow) {
          pAacDecoderChannelInfo[1]->icsInfo =
              pAacDecoderChannelInfo[0]->icsInfo;
        }
        break;

      case common_max_sfb:
        if (FDKreadBit(hBs) == 0) {
          error = IcsReadMaxSfb(hBs, &pAacDecoderChannelInfo[1]->icsInfo,
                                pSamplingRateInfo);
        }
        break;

      case ltp_data_present:
        if (FDKreadBits(hBs, 1) != 0) {
          error = AAC_DEC_UNSUPPORTED_PREDICTION;
        }
        break;

      case ms:

        INT max_sfb_ste;
        INT max_sfb_ste_clear;

        max_sfb_ste = GetScaleMaxFactorBandsTransmitted(
            &pAacDecoderChannelInfo[0]->icsInfo,
            &pAacDecoderChannelInfo[1]->icsInfo);

        max_sfb_ste_clear = 64;

        pAacDecoderChannelInfo[0]->icsInfo.max_sfb_ste = (UCHAR)max_sfb_ste;
        pAacDecoderChannelInfo[1]->icsInfo.max_sfb_ste = (UCHAR)max_sfb_ste;

        if (flags & (AC_USAC | AC_RSV603DA) &&
            pAacDecoderChannelInfo[ch]->pDynData->RawDataInfo.CommonWindow ==
                0) {
          Clean_Complex_Prediction_coefficients(
              &pAacDecoderStaticChannelInfo[0]
                   ->pCpeStaticData->jointStereoPersistentData,
              GetWindowGroups(&pAacDecoderChannelInfo[0]->icsInfo), 0, 64);
        }

        if (CJointStereo_Read(
                hBs, &pAacDecoderChannelInfo[0]->pComData->jointStereoData,
                GetWindowGroups(&pAacDecoderChannelInfo[0]->icsInfo),
                max_sfb_ste, max_sfb_ste_clear,
                ((cplxPred == 0) || (pAacDecoderStaticChannelInfo == NULL))
                    ? NULL
                    : &pAacDecoderStaticChannelInfo[0]
                           ->pCpeStaticData->jointStereoPersistentData,
                ((cplxPred == 0) || (pAacDecoderChannelInfo[0] == NULL))
                    ? NULL
                    : pAacDecoderChannelInfo[0]
                          ->pComStaticData->cplxPredictionData,
                cplxPred,
                GetScaleFactorBandsTotal(&pAacDecoderChannelInfo[0]->icsInfo),
                GetWindowSequence(&pAacDecoderChannelInfo[0]->icsInfo),
                flags)) {
          error = AAC_DEC_PARSE_ERROR;
        }

        break;

      case global_gain:
        pAacDecoderChannelInfo[ch]->pDynData->RawDataInfo.GlobalGain =
            (UCHAR)FDKreadBits(hBs, 8);
        break;

      case section_data:
        error = CBlock_ReadSectionData(hBs, pAacDecoderChannelInfo[ch],
                                       pSamplingRateInfo, flags);
        break;

      case scale_factor_data_usac:
        pAacDecoderChannelInfo[ch]->currAliasingSymmetry = 0;
        CChannel_CodebookTableInit(
            pAacDecoderChannelInfo[ch]);
        FDK_FALLTHROUGH;

      case scale_factor_data:
        if (flags & AC_ER_RVLC) {
          CRvlc_Read(pAacDecoderChannelInfo[ch], hBs);
        } else {
          error = CBlock_ReadScaleFactorData(pAacDecoderChannelInfo[ch], hBs,
                                             flags);
        }
        break;

      case pulse:
        if (CPulseData_Read(
                hBs,
                &pAacDecoderChannelInfo[ch]->pDynData->specificTo.aac.PulseData,
                pSamplingRateInfo->ScaleFactorBands_Long,
                (void *)&pAacDecoderChannelInfo[ch]->icsInfo,
                frame_length) != 0) {
          error = AAC_DEC_DECODE_FRAME_ERROR;
        }
        break;
      case tns_data_present:
        CTns_ReadDataPresentFlag(
            hBs, &pAacDecoderChannelInfo[ch]->pDynData->TnsData);
        if (elFlags & AC_EL_LFE &&
            pAacDecoderChannelInfo[ch]->pDynData->TnsData.DataPresent) {
          error = AAC_DEC_PARSE_ERROR;
        }
        break;
      case tns_data:
        error = CTns_Read(hBs, &pAacDecoderChannelInfo[ch]->pDynData->TnsData,
                          &pAacDecoderChannelInfo[ch]->icsInfo, flags);

        break;

      case gain_control_data:
        break;

      case gain_control_data_present:
        if (FDKreadBits(hBs, 1)) {
          error = AAC_DEC_UNSUPPORTED_GAIN_CONTROL_DATA;
        }
        break;

      case tw_data:
        break;
      case common_tw:
        break;
      case tns_data_present_usac:
        if (pAacDecoderChannelInfo[0]->pDynData->specificTo.usac.tns_active) {
          CTns_ReadDataPresentUsac(
              hBs, &pAacDecoderChannelInfo[0]->pDynData->TnsData,
              &pAacDecoderChannelInfo[1]->pDynData->TnsData,
              &pAacDecoderChannelInfo[0]->pDynData->specificTo.usac.tns_on_lr,
              &pAacDecoderChannelInfo[0]->icsInfo, flags, elFlags,
              pAacDecoderChannelInfo[0]->pDynData->RawDataInfo.CommonWindow);
        } else {
          pAacDecoderChannelInfo[0]->pDynData->specificTo.usac.tns_on_lr =
              (UCHAR)1;
        }
        break;
      case core_mode:
        decision_bit = FDKreadBits(hBs, 1);
        pAacDecoderChannelInfo[ch]->data.usac.core_mode = decision_bit;
        if ((ch == 1) && (pAacDecoderChannelInfo[0]->data.usac.core_mode !=
                          pAacDecoderChannelInfo[1]->data.usac.core_mode)) {
          pAacDecoderChannelInfo[0]->pDynData->RawDataInfo.CommonWindow = 0;
          pAacDecoderChannelInfo[1]->pDynData->RawDataInfo.CommonWindow = 0;
        }
        break;
      case tns_active:
        pAacDecoderChannelInfo[0]->pDynData->specificTo.usac.tns_active =
            FDKreadBit(hBs);
        break;
      case noise:
        if (elFlags & AC_EL_USAC_NOISE) {
          pAacDecoderChannelInfo[ch]
              ->pDynData->specificTo.usac.fd_noise_level_and_offset =
              FDKreadBits(hBs, 3 + 5);
        }
        break;
      case lpd_channel_stream:

      {
        error = CLpdChannelStream_Read(
                                       hBs, pAacDecoderChannelInfo[ch],
                                       pAacDecoderStaticChannelInfo[ch],
                                       pSamplingRateInfo, flags);
      }

        pAacDecoderChannelInfo[ch]->renderMode = AACDEC_RENDER_LPD;
        break;
      case fac_data: {
        int fFacDatPresent = FDKreadBit(hBs);

        pAacDecoderChannelInfo[ch]->data.usac.fac_data[0] =
            pAacDecoderChannelInfo[ch]->data.usac.fac_data0;

        if (fFacDatPresent) {
          if (elFlags & AC_EL_LFE) {
            error = AAC_DEC_PARSE_ERROR;
            break;
          }
          if (pAacDecoderStaticChannelInfo[ch]->last_core_mode != LPD ||
              pAacDecoderStaticChannelInfo[ch]->last_lpd_mode != 0) {
            pAacDecoderChannelInfo[ch]->data.usac.core_mode_last = LPD;
            pAacDecoderChannelInfo[ch]->data.usac.lpd_mode_last = 0;
          }
          CLpd_FAC_Read(hBs, pAacDecoderChannelInfo[ch]->data.usac.fac_data[0],
                        pAacDecoderChannelInfo[ch]->data.usac.fac_data_e,
                        CLpd_FAC_getLength(
                            IsLongBlock(&pAacDecoderChannelInfo[ch]->icsInfo),
                            pAacDecoderChannelInfo[ch]->granuleLength),
                        1, 0);
        } else {
          if (pAacDecoderStaticChannelInfo[ch]->last_core_mode == LPD &&
              pAacDecoderStaticChannelInfo[ch]->last_lpd_mode == 0) {
            FDKmemclear(pAacDecoderChannelInfo[ch]->data.usac.fac_data0,
                        LFAC * sizeof(FIXP_DBL));
          }
        }
      } break;
      case esc2_rvlc:
        if (flags & AC_ER_RVLC) {
          CRvlc_Decode(pAacDecoderChannelInfo[ch],
                       pAacDecoderStaticChannelInfo[ch], hBs);
        }
        break;

      case esc1_hcr:
        if (flags & AC_ER_HCR) {
          CHcr_Read(hBs, pAacDecoderChannelInfo[ch],
                    numberOfChannels == 2 ? ID_CPE : ID_SCE);
        }
        break;

      case spectral_data:
        error = CBlock_ReadSpectralData(hBs, pAacDecoderChannelInfo[ch],
                                        pSamplingRateInfo, flags);
        if (flags & AC_ELD) {
          pAacDecoderChannelInfo[ch]->renderMode = AACDEC_RENDER_ELDFB;
        } else {
          if (flags & AC_HDAAC) {
            pAacDecoderChannelInfo[ch]->renderMode = AACDEC_RENDER_INTIMDCT;
          } else {
            pAacDecoderChannelInfo[ch]->renderMode = AACDEC_RENDER_IMDCT;
          }
        }
        break;

      case ac_spectral_data:
        error = CBlock_ReadAcSpectralData(
            hBs, pAacDecoderChannelInfo[ch], pAacDecoderStaticChannelInfo[ch],
            pSamplingRateInfo, frame_length, flags);
        pAacDecoderChannelInfo[ch]->renderMode = AACDEC_RENDER_IMDCT;
        break;

      case coupled_elements: {
        int num_coupled_elements, c;

        ind_sw_cce_flag = FDKreadBit(hBs);
        num_coupled_elements = FDKreadBits(hBs, 3);

        for (c = 0; c < (num_coupled_elements + 1); c++) {
          int cc_target_is_cpe;

          num_gain_element_lists++;
          cc_target_is_cpe = FDKreadBit(hBs);
          FDKreadBits(hBs, 4);

          if (cc_target_is_cpe) {
            int cc_l, cc_r;

            cc_l = FDKreadBit(hBs);
            cc_r = FDKreadBit(hBs);

            if (cc_l && cc_r) {
              num_gain_element_lists++;
            }
          }
        }
        FDKreadBit(hBs);
        FDKreadBit(hBs);
        FDKreadBits(hBs, 2);
      } break;

      case gain_element_lists: {
        const CodeBookDescription *hcb;
        UCHAR *pCodeBook;
        int c;

        hcb = &AACcodeBookDescriptionTable[BOOKSCL];
        pCodeBook = pAacDecoderChannelInfo[ch]->pDynData->aCodeBook;

        for (c = 1; c < num_gain_element_lists; c++) {
          int cge;
          if (ind_sw_cce_flag) {
            cge = 1;
          } else {
            cge = FDKreadBits(hBs, 1);
          }
          if (cge) {
            CBlock_DecodeHuffmanWord(
                hBs, hcb);
          } else {
            int g, sfb;
            for (g = 0;
                 g < GetWindowGroups(&pAacDecoderChannelInfo[ch]->icsInfo);
                 g++) {
              for (sfb = 0; sfb < GetScaleFactorBandsTransmitted(
                                      &pAacDecoderChannelInfo[ch]->icsInfo);
                   sfb++) {
                if (pCodeBook[sfb] != ZERO_HCB) {
                  CBlock_DecodeHuffmanWord(
                      hBs,
                      hcb);
                }
              }
            }
          }
        }
      } break;

      case adtscrc_start_reg1:
        if (pTpDec != NULL) {
          crcReg1 = transportDec_CrcStartReg(pTpDec, 192);
        }
        break;
      case adtscrc_start_reg2:
        if (pTpDec != NULL) {
          crcReg2 = transportDec_CrcStartReg(pTpDec, 128);
        }
        break;
      case adtscrc_end_reg1:
      case drmcrc_end_reg:
        if (pTpDec != NULL) {
          transportDec_CrcEndReg(pTpDec, crcReg1);
          crcReg1 = -1;
        }
        break;
      case adtscrc_end_reg2:
        if (crcReg1 != -1) {
          error = AAC_DEC_DECODE_FRAME_ERROR;
        } else if (pTpDec != NULL) {
          transportDec_CrcEndReg(pTpDec, crcReg2);
          crcReg2 = -1;
        }
        break;
      case drmcrc_start_reg:
        if (pTpDec != NULL) {
          crcReg1 = transportDec_CrcStartReg(pTpDec, 0);
        }
        break;

      case next_channel:
        ch = (ch + 1) % numberOfChannels;
        break;
      case link_sequence:
        list = list->next[decision_bit];
        i = -1;
        break;

      default:
        error = AAC_DEC_UNSUPPORTED_FORMAT;
        break;
    }

    if (error != AAC_DEC_OK) {
      goto bail;
    }

    i++;

  } while (list->id[i] != end_of_sequence);

  for (ch = 0; ch < numberOfChannels; ch++) {
    if (pAacDecoderChannelInfo[ch]->renderMode == AACDEC_RENDER_IMDCT ||
        pAacDecoderChannelInfo[ch]->renderMode == AACDEC_RENDER_ELDFB) {
      UCHAR *band_is_noise =
          pAacDecoderChannelInfo[ch]->pDynData->band_is_noise;
      FDKmemset(band_is_noise, (UCHAR)1, sizeof(UCHAR) * (8 * 16));

      error = CBlock_InverseQuantizeSpectralData(
          pAacDecoderChannelInfo[ch], pSamplingRateInfo, band_is_noise, 1);
      if (error != AAC_DEC_OK) {
        return error;
      }

      if (elFlags & AC_EL_USAC_NOISE) {
        CBlock_ApplyNoise(pAacDecoderChannelInfo[ch], pSamplingRateInfo,
                          &pAacDecoderStaticChannelInfo[ch]->nfRandomSeed,
                          band_is_noise);
      }
    }
  }

bail:
  if (crcReg1 != -1 || crcReg2 != -1) {
    if (error == AAC_DEC_OK) {
      error = AAC_DEC_DECODE_FRAME_ERROR;
    }
    if (crcReg1 != -1) {
      transportDec_CrcEndReg(pTpDec, crcReg1);
    }
    if (crcReg2 != -1) {
      transportDec_CrcEndReg(pTpDec, crcReg2);
    }
  }
  return error;
}
