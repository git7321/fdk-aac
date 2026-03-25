/**************************** AAC decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "rvlcconceal.h"
#include "block.h"
#include "rvlc.h"

static void calcRefValFwd(CErRvlcInfo *pRvlc,
                          CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                          int *refIsFwd, int *refNrgFwd, int *refScfFwd) {
  int band, bnds, group, startBand;
  int idIs, idNrg, idScf;
  int conceal_min, conceal_group_min;
  int MaximumScaleFactorBands;

  if (GetWindowSequence(&pAacDecoderChannelInfo->icsInfo) == BLOCK_SHORT)
    MaximumScaleFactorBands = 16;
  else
    MaximumScaleFactorBands = 64;

  conceal_min = pRvlc->conceal_min % MaximumScaleFactorBands;
  conceal_group_min = pRvlc->conceal_min / MaximumScaleFactorBands;

  idIs = idNrg = idScf = 1;

  *refIsFwd = -SF_OFFSET;
  *refNrgFwd = pAacDecoderChannelInfo->pDynData->RawDataInfo.GlobalGain -
               SF_OFFSET - 90 - 256;
  *refScfFwd =
      pAacDecoderChannelInfo->pDynData->RawDataInfo.GlobalGain - SF_OFFSET;

  startBand = conceal_min - 1;
  for (group = conceal_group_min; group >= 0; group--) {
    for (band = startBand; band >= 0; band--) {
      bnds = 16 * group + band;
      switch (pAacDecoderChannelInfo->pDynData->aCodeBook[bnds]) {
        case ZERO_HCB:
          break;
        case INTENSITY_HCB:
        case INTENSITY_HCB2:
          if (idIs) {
            *refIsFwd =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
            idIs = 0;
          }
          break;
        case NOISE_HCB:
          if (idNrg) {
            *refNrgFwd =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
            idNrg = 0;
          }
          break;
        default:
          if (idScf) {
            *refScfFwd =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
            idScf = 0;
          }
          break;
      }
    }
    startBand = pRvlc->maxSfbTransmitted - 1;
  }
}

static void calcRefValBwd(CErRvlcInfo *pRvlc,
                          CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                          int *refIsBwd, int *refNrgBwd, int *refScfBwd) {
  int band, bnds, group, startBand;
  int idIs, idNrg, idScf;
  int conceal_max, conceal_group_max;
  int MaximumScaleFactorBands;

  if (GetWindowSequence(&pAacDecoderChannelInfo->icsInfo) == BLOCK_SHORT)
    MaximumScaleFactorBands = 16;
  else
    MaximumScaleFactorBands = 64;

  conceal_max = pRvlc->conceal_max % MaximumScaleFactorBands;
  conceal_group_max = pRvlc->conceal_max / MaximumScaleFactorBands;

  idIs = idNrg = idScf = 1;

  *refIsBwd = pRvlc->dpcm_is_last_position - SF_OFFSET;
  *refNrgBwd = pRvlc->rev_global_gain + pRvlc->dpcm_noise_last_position -
               SF_OFFSET - 90 - 256 + pRvlc->dpcm_noise_nrg;
  *refScfBwd = pRvlc->rev_global_gain - SF_OFFSET;

  startBand = conceal_max + 1;

  for (group = conceal_group_max; group < pRvlc->numWindowGroups; group++) {
    for (band = startBand; band < pRvlc->maxSfbTransmitted; band++) {
      bnds = 16 * group + band;
      switch (pAacDecoderChannelInfo->pDynData->aCodeBook[bnds]) {
        case ZERO_HCB:
          break;
        case INTENSITY_HCB:
        case INTENSITY_HCB2:
          if (idIs) {
            *refIsBwd =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
            idIs = 0;
          }
          break;
        case NOISE_HCB:
          if (idNrg) {
            *refNrgBwd =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
            idNrg = 0;
          }
          break;
        default:
          if (idScf) {
            *refScfBwd =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
            idScf = 0;
          }
          break;
      }
    }
    startBand = 0;
  }
}

void BidirectionalEstimation_UseLowerScfOfCurrentFrame(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo) {
  CErRvlcInfo *pRvlc =
      &pAacDecoderChannelInfo->pComData->overlay.aac.erRvlcInfo;
  int band, bnds, startBand, endBand, group;
  int conceal_min, conceal_max;
  int conceal_group_min, conceal_group_max;
  int MaximumScaleFactorBands;

  if (GetWindowSequence(&pAacDecoderChannelInfo->icsInfo) == BLOCK_SHORT) {
    MaximumScaleFactorBands = 16;
  } else {
    MaximumScaleFactorBands = 64;
  }

  if (pRvlc->conceal_min == CONCEAL_MIN_INIT) pRvlc->conceal_min = 0;

  if (pRvlc->conceal_max == CONCEAL_MAX_INIT)
    pRvlc->conceal_max =
        (pRvlc->numWindowGroups - 1) * 16 + pRvlc->maxSfbTransmitted - 1;

  conceal_min = pRvlc->conceal_min % MaximumScaleFactorBands;
  conceal_group_min = pRvlc->conceal_min / MaximumScaleFactorBands;
  conceal_max = pRvlc->conceal_max % MaximumScaleFactorBands;
  conceal_group_max = pRvlc->conceal_max / MaximumScaleFactorBands;

  if (pRvlc->conceal_min == pRvlc->conceal_max) {
    int refIsFwd, refNrgFwd, refScfFwd;
    int refIsBwd, refNrgBwd, refScfBwd;

    bnds = pRvlc->conceal_min;
    calcRefValFwd(pRvlc, pAacDecoderChannelInfo, &refIsFwd, &refNrgFwd,
                  &refScfFwd);
    calcRefValBwd(pRvlc, pAacDecoderChannelInfo, &refIsBwd, &refNrgBwd,
                  &refScfBwd);

    switch (pAacDecoderChannelInfo->pDynData->aCodeBook[bnds]) {
      case ZERO_HCB:
        break;
      case INTENSITY_HCB:
      case INTENSITY_HCB2:
        if (refIsFwd < refIsBwd)
          pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = refIsFwd;
        else
          pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = refIsBwd;
        break;
      case NOISE_HCB:
        if (refNrgFwd < refNrgBwd)
          pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = refNrgFwd;
        else
          pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = refNrgBwd;
        break;
      default:
        if (refScfFwd < refScfBwd)
          pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = refScfFwd;
        else
          pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = refScfBwd;
        break;
    }
  } else {
    pAacDecoderChannelInfo->pComData->overlay.aac
        .aRvlcScfFwd[pRvlc->conceal_max] =
        pAacDecoderChannelInfo->pComData->overlay.aac
            .aRvlcScfBwd[pRvlc->conceal_max];
    pAacDecoderChannelInfo->pComData->overlay.aac
        .aRvlcScfBwd[pRvlc->conceal_min] =
        pAacDecoderChannelInfo->pComData->overlay.aac
            .aRvlcScfFwd[pRvlc->conceal_min];

    startBand = conceal_min;
    if (conceal_group_min == conceal_group_max)
      endBand = conceal_max;
    else
      endBand = pRvlc->maxSfbTransmitted - 1;

    for (group = conceal_group_min; group <= conceal_group_max; group++) {
      for (band = startBand; band <= endBand; band++) {
        bnds = 16 * group + band;
        if (pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds] <
            pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds])
          pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
              pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
        else
          pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
              pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
      }
      startBand = 0;
      if ((group + 1) == conceal_group_max) endBand = conceal_max;
    }
  }

  if (conceal_group_min == 0)
    endBand = conceal_min;
  else
    endBand = pRvlc->maxSfbTransmitted;
  for (group = 0; group <= conceal_group_min; group++) {
    for (band = 0; band < endBand; band++) {
      bnds = 16 * group + band;
      pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
          pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
    }
    if ((group + 1) == conceal_group_min) endBand = conceal_min;
  }

  startBand = conceal_max + 1;
  for (group = conceal_group_max; group < pRvlc->numWindowGroups; group++) {
    for (band = startBand; band < pRvlc->maxSfbTransmitted; band++) {
      bnds = 16 * group + band;
      pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
          pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
    }
    startBand = 0;
  }
}

void BidirectionalEstimation_UseScfOfPrevFrameAsReference(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo) {
  CErRvlcInfo *pRvlc =
      &pAacDecoderChannelInfo->pComData->overlay.aac.erRvlcInfo;
  int band, bnds, startBand, endBand, group;
  int conceal_min, conceal_max;
  int conceal_group_min, conceal_group_max;
  int MaximumScaleFactorBands;
  SHORT commonMin;

  if (GetWindowSequence(&pAacDecoderChannelInfo->icsInfo) == BLOCK_SHORT) {
    MaximumScaleFactorBands = 16;
  } else {
    MaximumScaleFactorBands = 64;
  }

  if (pRvlc->conceal_min == CONCEAL_MIN_INIT) pRvlc->conceal_min = 0;

  if (pRvlc->conceal_max == CONCEAL_MAX_INIT)
    pRvlc->conceal_max =
        (pRvlc->numWindowGroups - 1) * 16 + pRvlc->maxSfbTransmitted - 1;

  conceal_min = pRvlc->conceal_min % MaximumScaleFactorBands;
  conceal_group_min = pRvlc->conceal_min / MaximumScaleFactorBands;
  conceal_max = pRvlc->conceal_max % MaximumScaleFactorBands;
  conceal_group_max = pRvlc->conceal_max / MaximumScaleFactorBands;

  pAacDecoderChannelInfo->pComData->overlay.aac
      .aRvlcScfFwd[pRvlc->conceal_max] =
      pAacDecoderChannelInfo->pComData->overlay.aac
          .aRvlcScfBwd[pRvlc->conceal_max];
  pAacDecoderChannelInfo->pComData->overlay.aac
      .aRvlcScfBwd[pRvlc->conceal_min] =
      pAacDecoderChannelInfo->pComData->overlay.aac
          .aRvlcScfFwd[pRvlc->conceal_min];

  startBand = conceal_min;
  if (conceal_group_min == conceal_group_max)
    endBand = conceal_max;
  else
    endBand = pRvlc->maxSfbTransmitted - 1;

  for (group = conceal_group_min; group <= conceal_group_max; group++) {
    for (band = startBand; band <= endBand; band++) {
      bnds = 16 * group + band;
      switch (pAacDecoderChannelInfo->pDynData->aCodeBook[bnds]) {
        case ZERO_HCB:
          pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = 0;
          break;

        case INTENSITY_HCB:
        case INTENSITY_HCB2:
          if ((pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] == INTENSITY_HCB) ||
              (pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] == INTENSITY_HCB2)) {
            commonMin = fMin(
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds],
                pAacDecoderChannelInfo->pComData->overlay.aac
                    .aRvlcScfBwd[bnds]);
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                fMin(commonMin, pAacDecoderStaticChannelInfo->concealmentInfo
                                    .aRvlcPreviousScaleFactor[bnds]);
          } else {
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = fMin(
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds],
                pAacDecoderChannelInfo->pComData->overlay.aac
                    .aRvlcScfBwd[bnds]);
          }
          break;

        case NOISE_HCB:
          if (pAacDecoderStaticChannelInfo->concealmentInfo
                  .aRvlcPreviousCodebook[bnds] == NOISE_HCB) {
            commonMin = fMin(
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds],
                pAacDecoderChannelInfo->pComData->overlay.aac
                    .aRvlcScfBwd[bnds]);
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                fMin(commonMin, pAacDecoderStaticChannelInfo->concealmentInfo
                                    .aRvlcPreviousScaleFactor[bnds]);
          } else {
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = fMin(
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds],
                pAacDecoderChannelInfo->pComData->overlay.aac
                    .aRvlcScfBwd[bnds]);
          }
          break;

        default:
          if ((pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] != ZERO_HCB) &&
              (pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] != NOISE_HCB) &&
              (pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] != INTENSITY_HCB) &&
              (pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] != INTENSITY_HCB2)) {
            commonMin = fMin(
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds],
                pAacDecoderChannelInfo->pComData->overlay.aac
                    .aRvlcScfBwd[bnds]);
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                fMin(commonMin, pAacDecoderStaticChannelInfo->concealmentInfo
                                    .aRvlcPreviousScaleFactor[bnds]);
          } else {
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = fMin(
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds],
                pAacDecoderChannelInfo->pComData->overlay.aac
                    .aRvlcScfBwd[bnds]);
          }
          break;
      }
    }
    startBand = 0;
    if ((group + 1) == conceal_group_max) endBand = conceal_max;
  }

  if (conceal_group_min == 0)
    endBand = conceal_min;
  else
    endBand = pRvlc->maxSfbTransmitted;
  for (group = 0; group <= conceal_group_min; group++) {
    for (band = 0; band < endBand; band++) {
      bnds = 16 * group + band;
      pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
          pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
    }
    if ((group + 1) == conceal_group_min) endBand = conceal_min;
  }

  startBand = conceal_max + 1;
  for (group = conceal_group_max; group < pRvlc->numWindowGroups; group++) {
    for (band = startBand; band < pRvlc->maxSfbTransmitted; band++) {
      bnds = 16 * group + band;
      pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
          pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
    }
    startBand = 0;
  }
}

void StatisticalEstimation(CAacDecoderChannelInfo *pAacDecoderChannelInfo) {
  CErRvlcInfo *pRvlc =
      &pAacDecoderChannelInfo->pComData->overlay.aac.erRvlcInfo;
  int band, bnds, group;
  int sumIsFwd, sumIsBwd;
  int sumNrgFwd, sumNrgBwd;
  int sumScfFwd, sumScfBwd;
  int useIsFwd, useNrgFwd, useScfFwd;

  sumIsFwd = sumIsBwd = sumNrgFwd = sumNrgBwd = sumScfFwd = sumScfBwd = 0;
  useIsFwd = useNrgFwd = useScfFwd = 0;

  for (group = 0; group < pRvlc->numWindowGroups; group++) {
    for (band = 0; band < pRvlc->maxSfbTransmitted; band++) {
      bnds = 16 * group + band;
      switch (pAacDecoderChannelInfo->pDynData->aCodeBook[bnds]) {
        case ZERO_HCB:
          break;

        case INTENSITY_HCB:
        case INTENSITY_HCB2:
          sumIsFwd +=
              pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
          sumIsBwd +=
              pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
          break;

        case NOISE_HCB:
          sumNrgFwd +=
              pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
          sumNrgBwd +=
              pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
          break;

        default:
          sumScfFwd +=
              pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
          sumScfBwd +=
              pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
          break;
      }
    }
  }

  if (sumIsFwd < sumIsBwd) useIsFwd = 1;
  if (sumNrgFwd < sumNrgBwd) useNrgFwd = 1;
  if (sumScfFwd < sumScfBwd) useScfFwd = 1;

  for (group = 0; group < pRvlc->numWindowGroups; group++) {
    for (band = 0; band < pRvlc->maxSfbTransmitted; band++) {
      bnds = 16 * group + band;
      switch (pAacDecoderChannelInfo->pDynData->aCodeBook[bnds]) {
        case ZERO_HCB:
          break;

        case INTENSITY_HCB:
        case INTENSITY_HCB2:
          if (useIsFwd)
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
          else
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
          break;

        case NOISE_HCB:
          if (useNrgFwd)
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
          else
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
          break;

        default:
          if (useScfFwd)
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
          else
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd[bnds];
          break;
      }
    }
  }
}

void PredictiveInterpolation(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo) {
  CErRvlcInfo *pRvlc =
      &pAacDecoderChannelInfo->pComData->overlay.aac.erRvlcInfo;
  int band, bnds, group;
  SHORT commonMin;

  for (group = 0; group < pRvlc->numWindowGroups; group++) {
    for (band = 0; band < pRvlc->maxSfbTransmitted; band++) {
      bnds = 16 * group + band;
      switch (pAacDecoderChannelInfo->pDynData->aCodeBook[bnds]) {
        case ZERO_HCB:
          pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = 0;
          break;

        case INTENSITY_HCB:
        case INTENSITY_HCB2:
          if ((pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] == INTENSITY_HCB) ||
              (pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] == INTENSITY_HCB2)) {
            commonMin = fMin(
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds],
                pAacDecoderChannelInfo->pComData->overlay.aac
                    .aRvlcScfBwd[bnds]);
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                fMin(commonMin, pAacDecoderStaticChannelInfo->concealmentInfo
                                    .aRvlcPreviousScaleFactor[bnds]);
          } else {
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = -110;
          }
          break;

        case NOISE_HCB:
          if (pAacDecoderStaticChannelInfo->concealmentInfo
                  .aRvlcPreviousCodebook[bnds] == NOISE_HCB) {
            commonMin = fMin(
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds],
                pAacDecoderChannelInfo->pComData->overlay.aac
                    .aRvlcScfBwd[bnds]);
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                fMin(commonMin, pAacDecoderStaticChannelInfo->concealmentInfo
                                    .aRvlcPreviousScaleFactor[bnds]);
          } else {
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = -110;
          }
          break;

        default:
          if ((pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] != ZERO_HCB) &&
              (pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] != NOISE_HCB) &&
              (pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] != INTENSITY_HCB) &&
              (pAacDecoderStaticChannelInfo->concealmentInfo
                   .aRvlcPreviousCodebook[bnds] != INTENSITY_HCB2)) {
            commonMin = fMin(
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds],
                pAacDecoderChannelInfo->pComData->overlay.aac
                    .aRvlcScfBwd[bnds]);
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
                fMin(commonMin, pAacDecoderStaticChannelInfo->concealmentInfo
                                    .aRvlcPreviousScaleFactor[bnds]);
          } else {
            pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] = 0;
          }
          break;
      }
    }
  }
}
