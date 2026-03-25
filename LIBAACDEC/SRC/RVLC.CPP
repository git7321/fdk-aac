/**************************** AAC decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "rvlc.h"
#include "block.h"
#include "aac_rom.h"
#include "rvlcbit.h"
#include "rvlcconceal.h"
#include "aacdec_hcr.h"

static void rvlcInit(CErRvlcInfo *pRvlc,
                     CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                     HANDLE_FDK_BITSTREAM bs) {
  SHORT *pScfEsc = pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfEsc;
  SHORT *pScfFwd = pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd;
  SHORT *pScfBwd = pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd;
  SHORT *pScaleFactor = pAacDecoderChannelInfo->pDynData->aScaleFactor;
  int bnds;

  pAacDecoderChannelInfo->pDynData->specificTo.aac.rvlcIntensityUsed = 0;

  pRvlc->numDecodedEscapeWordsEsc = 0;
  pRvlc->numDecodedEscapeWordsFwd = 0;
  pRvlc->numDecodedEscapeWordsBwd = 0;

  pRvlc->intensity_used = 0;
  pRvlc->errorLogRvlc = 0;

  pRvlc->conceal_max = CONCEAL_MAX_INIT;
  pRvlc->conceal_min = CONCEAL_MIN_INIT;

  pRvlc->conceal_max_esc = CONCEAL_MAX_INIT;
  pRvlc->conceal_min_esc = CONCEAL_MIN_INIT;

  pRvlc->pHuffTreeRvlcEscape = aHuffTreeRvlcEscape;
  pRvlc->pHuffTreeRvlCodewds = aHuffTreeRvlCodewds;

  for (bnds = 0; bnds < RVLC_MAX_SFB; bnds++) {
    pScfFwd[bnds] = 0;
    pScfBwd[bnds] = 0;
    pScfEsc[bnds] = 0;
    pScaleFactor[bnds] = 0;
  }

  FDKsyncCache(bs);
  pRvlc->bsAnchor = (INT)FDKgetValidBits(bs);

  pRvlc->bitstreamIndexRvlFwd =
      0;
  pRvlc->bitstreamIndexRvlBwd =
      pRvlc->length_of_rvlc_sf - 1;

  FDKpushFor(bs, pRvlc->length_of_rvlc_sf);

  if (pRvlc->sf_escapes_present != 0) {
    FDKsyncCache(bs);
    pRvlc->bitstreamIndexEsc = pRvlc->bsAnchor - (INT)FDKgetValidBits(bs);

    FDKpushFor(bs, pRvlc->length_of_rvlc_escapes);
  }
}

static void rvlcCheckIntensityCb(
    CErRvlcInfo *pRvlc, CAacDecoderChannelInfo *pAacDecoderChannelInfo) {
  int group, band, bnds;

  pRvlc->intensity_used = 0;

  for (group = 0; group < pRvlc->numWindowGroups; group++) {
    for (band = 0; band < pRvlc->maxSfbTransmitted; band++) {
      bnds = 16 * group + band;
      if ((pAacDecoderChannelInfo->pDynData->aCodeBook[bnds] ==
           INTENSITY_HCB) ||
          (pAacDecoderChannelInfo->pDynData->aCodeBook[bnds] ==
           INTENSITY_HCB2)) {
        pRvlc->intensity_used = 1;
        break;
      }
    }
  }
}

static SCHAR rvlcDecodeEscapeWord(CErRvlcInfo *pRvlc, HANDLE_FDK_BITSTREAM bs) {
  int i;
  SCHAR value;
  UCHAR carryBit;
  UINT treeNode;
  UINT branchValue;
  UINT branchNode;

  INT *pBitstreamIndexEsc;
  const UINT *pEscTree;

  pEscTree = pRvlc->pHuffTreeRvlcEscape;
  pBitstreamIndexEsc = &(pRvlc->bitstreamIndexEsc);
  treeNode = *pEscTree;

  for (i = MAX_LEN_RVLC_ESCAPE_WORD - 1; i >= 0; i--) {
    carryBit =
        rvlcReadBitFromBitstream(bs,
                                 pRvlc->bsAnchor, pBitstreamIndexEsc, FWD);

    CarryBitToBranchValue(carryBit,
                          treeNode, &branchValue, &branchNode);

    if ((branchNode & TEST_BIT_10) ==
        TEST_BIT_10) {
      value = (SCHAR)branchNode & CLR_BIT_10;
      pRvlc->length_of_rvlc_escapes -= (MAX_LEN_RVLC_ESCAPE_WORD - i);

      if (pRvlc->length_of_rvlc_escapes < 0) {
        pRvlc->errorLogRvlc |= RVLC_ERROR_ALL_ESCAPE_WORDS_INVALID;
        value = -1;
      }

      return value;
    } else {
      treeNode = *(
          pEscTree +
          branchValue);
    }
  }

  pRvlc->errorLogRvlc |= RVLC_ERROR_ALL_ESCAPE_WORDS_INVALID;

  return -1;
}

static void rvlcDecodeEscapes(CErRvlcInfo *pRvlc, SHORT *pEsc,
                              HANDLE_FDK_BITSTREAM bs) {
  SCHAR escWord;
  SCHAR escCnt = 0;
  SHORT *pEscBitCntSum;

  pEscBitCntSum = &(pRvlc->length_of_rvlc_escapes);

  while (*pEscBitCntSum > 0) {
    escWord = rvlcDecodeEscapeWord(pRvlc, bs);

    if (escWord >= 0) {
      pEsc[escCnt] = escWord;
      escCnt++;
    } else {
      pRvlc->errorLogRvlc |= RVLC_ERROR_ALL_ESCAPE_WORDS_INVALID;
      pRvlc->numDecodedEscapeWordsEsc = escCnt;

      return;
    }
  }

  pRvlc->numDecodedEscapeWordsEsc = escCnt;
}

SCHAR decodeRVLCodeword(HANDLE_FDK_BITSTREAM bs, CErRvlcInfo *pRvlc) {
  int i;
  SCHAR value;
  UCHAR carryBit;
  UINT branchValue;
  UINT branchNode;

  const UINT *pRvlCodeTree = pRvlc->pHuffTreeRvlCodewds;
  UCHAR direction = pRvlc->direction;
  INT *pBitstrIndxRvl = pRvlc->pBitstrIndxRvl_RVL;
  UINT treeNode = *pRvlCodeTree;

  for (i = MAX_LEN_RVLC_CODE_WORD - 1; i >= 0; i--) {
    carryBit =
        rvlcReadBitFromBitstream(bs,
                                 pRvlc->bsAnchor, pBitstrIndxRvl, direction);

    CarryBitToBranchValue(carryBit,
                          treeNode, &branchValue, &branchNode);

    if ((branchNode & TEST_BIT_10) ==
        TEST_BIT_10) {
      value = (SCHAR)(branchNode & CLR_BIT_10);
      *pRvlc->pRvlBitCnt_RVL -= (MAX_LEN_RVLC_CODE_WORD - i);

      if (*pRvlc->pRvlBitCnt_RVL < 0) {
        if (direction == FWD) {
          pRvlc->errorLogRvlc |= RVLC_ERROR_RVL_SUM_BIT_COUNTER_BELOW_ZERO_FWD;
        } else {
          pRvlc->errorLogRvlc |= RVLC_ERROR_RVL_SUM_BIT_COUNTER_BELOW_ZERO_BWD;
        }
        value = -1;
      }

      if (value > MAX_ALLOWED_DPCM_INDEX) {
        if (direction == FWD) {
          pRvlc->errorLogRvlc |= RVLC_ERROR_FORBIDDEN_CW_DETECTED_FWD;
        } else {
          pRvlc->errorLogRvlc |= RVLC_ERROR_FORBIDDEN_CW_DETECTED_BWD;
        }
        value = -1;
      }

      return value;
    } else {
      treeNode = *(
          pRvlCodeTree +
          branchValue);
    }
  }

  return -1;
}

static void rvlcDecodeForward(CErRvlcInfo *pRvlc,
                              CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                              HANDLE_FDK_BITSTREAM bs) {
  int band = 0;
  int group = 0;
  int bnds = 0;

  SHORT dpcm;

  SHORT factor =
      pAacDecoderChannelInfo->pDynData->RawDataInfo.GlobalGain - SF_OFFSET;
  SHORT position = -SF_OFFSET;
  SHORT noisenrg = pAacDecoderChannelInfo->pDynData->RawDataInfo.GlobalGain -
                   SF_OFFSET - 90 - 256;

  SHORT *pScfFwd = pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd;
  SHORT *pScfEsc = pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfEsc;
  UCHAR *pEscFwdCnt = &(pRvlc->numDecodedEscapeWordsFwd);

  pRvlc->pRvlBitCnt_RVL = &(pRvlc->length_of_rvlc_sf_fwd);
  pRvlc->pBitstrIndxRvl_RVL = &(pRvlc->bitstreamIndexRvlFwd);

  *pEscFwdCnt = 0;
  pRvlc->direction = FWD;
  pRvlc->noise_used = 0;
  pRvlc->sf_used = 0;
  pRvlc->lastScf = 0;
  pRvlc->lastNrg = 0;
  pRvlc->lastIs = 0;

  rvlcCheckIntensityCb(pRvlc, pAacDecoderChannelInfo);

  for (group = 0; group < pRvlc->numWindowGroups; group++) {
    for (band = 0; band < pRvlc->maxSfbTransmitted; band++) {
      bnds = 16 * group + band;

      switch (pAacDecoderChannelInfo->pDynData->aCodeBook[bnds]) {
        case ZERO_HCB:
          pScfFwd[bnds] = 0;
          break;

        case INTENSITY_HCB2:
        case INTENSITY_HCB:
          dpcm = decodeRVLCodeword(bs, pRvlc);
          if (dpcm < 0) {
            pRvlc->conceal_max = bnds;
            return;
          }
          dpcm -= TABLE_OFFSET;
          if ((dpcm == MIN_RVL) || (dpcm == MAX_RVL)) {
            if (pRvlc->length_of_rvlc_escapes) {
              pRvlc->conceal_max = bnds;
              return;
            } else {
              if (dpcm == MIN_RVL) {
                dpcm -= *pScfEsc++;
              } else {
                dpcm += *pScfEsc++;
              }
              (*pEscFwdCnt)++;
              if (pRvlc->conceal_max_esc == CONCEAL_MAX_INIT) {
                pRvlc->conceal_max_esc = bnds;
              }
            }
          }
          position += dpcm;
          pScfFwd[bnds] = position;
          pRvlc->lastIs = position;
          break;

        case NOISE_HCB:
          if (pRvlc->noise_used == 0) {
            pRvlc->noise_used = 1;
            pRvlc->first_noise_band = bnds;
            noisenrg += pRvlc->dpcm_noise_nrg;
            pScfFwd[bnds] = 100 + noisenrg;
            pRvlc->lastNrg = noisenrg;
          } else {
            dpcm = decodeRVLCodeword(bs, pRvlc);
            if (dpcm < 0) {
              pRvlc->conceal_max = bnds;
              return;
            }
            dpcm -= TABLE_OFFSET;
            if ((dpcm == MIN_RVL) || (dpcm == MAX_RVL)) {
              if (pRvlc->length_of_rvlc_escapes) {
                pRvlc->conceal_max = bnds;
                return;
              } else {
                if (dpcm == MIN_RVL) {
                  dpcm -= *pScfEsc++;
                } else {
                  dpcm += *pScfEsc++;
                }
                (*pEscFwdCnt)++;
                if (pRvlc->conceal_max_esc == CONCEAL_MAX_INIT) {
                  pRvlc->conceal_max_esc = bnds;
                }
              }
            }
            noisenrg += dpcm;
            pScfFwd[bnds] = 100 + noisenrg;
            pRvlc->lastNrg = noisenrg;
          }
          pAacDecoderChannelInfo->data.aac.PnsData.pnsUsed[bnds] = 1;
          break;

        default:
          pRvlc->sf_used = 1;
          dpcm = decodeRVLCodeword(bs, pRvlc);
          if (dpcm < 0) {
            pRvlc->conceal_max = bnds;
            return;
          }
          dpcm -= TABLE_OFFSET;
          if ((dpcm == MIN_RVL) || (dpcm == MAX_RVL)) {
            if (pRvlc->length_of_rvlc_escapes) {
              pRvlc->conceal_max = bnds;
              return;
            } else {
              if (dpcm == MIN_RVL) {
                dpcm -= *pScfEsc++;
              } else {
                dpcm += *pScfEsc++;
              }
              (*pEscFwdCnt)++;
              if (pRvlc->conceal_max_esc == CONCEAL_MAX_INIT) {
                pRvlc->conceal_max_esc = bnds;
              }
            }
          }
          factor += dpcm;
          pScfFwd[bnds] = factor;
          pRvlc->lastScf = factor;
          break;
      }
    }
  }

  if (pRvlc->intensity_used) {
    dpcm = decodeRVLCodeword(bs, pRvlc);
    if (dpcm < 0) {
      pRvlc->conceal_max = bnds;
      return;
    }
    dpcm -= TABLE_OFFSET;
    if ((dpcm == MIN_RVL) || (dpcm == MAX_RVL)) {
      if (pRvlc->length_of_rvlc_escapes) {
        pRvlc->conceal_max = bnds;
        return;
      } else {
        if (dpcm == MIN_RVL) {
          dpcm -= *pScfEsc++;
        } else {
          dpcm += *pScfEsc++;
        }
        (*pEscFwdCnt)++;
        if (pRvlc->conceal_max_esc == CONCEAL_MAX_INIT) {
          pRvlc->conceal_max_esc = bnds;
        }
      }
    }
    pRvlc->dpcm_is_last_position = dpcm;
  }
}

static void rvlcDecodeBackward(CErRvlcInfo *pRvlc,
                               CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                               HANDLE_FDK_BITSTREAM bs) {
  SHORT band, group, dpcm, offset;
  SHORT bnds = pRvlc->maxSfbTransmitted - 1;

  SHORT factor = pRvlc->rev_global_gain - SF_OFFSET;
  SHORT position = pRvlc->dpcm_is_last_position - SF_OFFSET;
  SHORT noisenrg = pRvlc->rev_global_gain + pRvlc->dpcm_noise_last_position -
                   SF_OFFSET - 90 - 256;

  SHORT *pScfBwd = pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfBwd;
  SHORT *pScfEsc = pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfEsc;
  UCHAR escEscCnt = pRvlc->numDecodedEscapeWordsEsc;
  UCHAR *pEscBwdCnt = &(pRvlc->numDecodedEscapeWordsBwd);

  pRvlc->pRvlBitCnt_RVL = &(pRvlc->length_of_rvlc_sf_bwd);
  pRvlc->pBitstrIndxRvl_RVL = &(pRvlc->bitstreamIndexRvlBwd);

  *pEscBwdCnt = 0;
  pRvlc->direction = BWD;
  pScfEsc += escEscCnt - 1;
  pRvlc->firstScf = 0;
  pRvlc->firstNrg = 0;
  pRvlc->firstIs = 0;

  if (pRvlc->intensity_used) {
    dpcm = decodeRVLCodeword(bs, pRvlc);
    if (dpcm < 0) {
      pRvlc->dpcm_is_last_position = 0;
      pRvlc->conceal_min = bnds;
      return;
    }
    dpcm -= TABLE_OFFSET;
    if ((dpcm == MIN_RVL) || (dpcm == MAX_RVL)) {
      if ((pRvlc->length_of_rvlc_escapes) || (*pEscBwdCnt >= escEscCnt)) {
        pRvlc->conceal_min = bnds;
        return;
      } else {
        if (dpcm == MIN_RVL) {
          dpcm -= *pScfEsc--;
        } else {
          dpcm += *pScfEsc--;
        }
        (*pEscBwdCnt)++;
        if (pRvlc->conceal_min_esc == CONCEAL_MIN_INIT) {
          pRvlc->conceal_min_esc = bnds;
        }
      }
    }
    pRvlc->dpcm_is_last_position = dpcm;
  }

  for (group = pRvlc->numWindowGroups - 1; group >= 0; group--) {
    for (band = pRvlc->maxSfbTransmitted - 1; band >= 0; band--) {
      bnds = 16 * group + band;
      if ((band == 0) && (pRvlc->numWindowGroups != 1))
        offset = 16 - pRvlc->maxSfbTransmitted + 1;
      else
        offset = 1;

      switch (pAacDecoderChannelInfo->pDynData->aCodeBook[bnds]) {
        case ZERO_HCB:
          pScfBwd[bnds] = 0;
          break;

        case INTENSITY_HCB2:
        case INTENSITY_HCB:
          dpcm = decodeRVLCodeword(bs, pRvlc);
          if (dpcm < 0) {
            pScfBwd[bnds] = position;
            pRvlc->conceal_min = fMax(0, bnds - offset);
            return;
          }
          dpcm -= TABLE_OFFSET;
          if ((dpcm == MIN_RVL) || (dpcm == MAX_RVL)) {
            if ((pRvlc->length_of_rvlc_escapes) || (*pEscBwdCnt >= escEscCnt)) {
              pScfBwd[bnds] = position;
              pRvlc->conceal_min = fMax(0, bnds - offset);
              return;
            } else {
              if (dpcm == MIN_RVL) {
                dpcm -= *pScfEsc--;
              } else {
                dpcm += *pScfEsc--;
              }
              (*pEscBwdCnt)++;
              if (pRvlc->conceal_min_esc == CONCEAL_MIN_INIT) {
                pRvlc->conceal_min_esc = fMax(0, bnds - offset);
              }
            }
          }
          pScfBwd[bnds] = position;
          position -= dpcm;
          pRvlc->firstIs = position;
          break;

        case NOISE_HCB:
          if (bnds == pRvlc->first_noise_band) {
            pScfBwd[bnds] =
                pRvlc->dpcm_noise_nrg +
                pAacDecoderChannelInfo->pDynData->RawDataInfo.GlobalGain -
                SF_OFFSET - 90 - 256;
            pRvlc->firstNrg = pScfBwd[bnds];
          } else {
            dpcm = decodeRVLCodeword(bs, pRvlc);
            if (dpcm < 0) {
              pScfBwd[bnds] = noisenrg;
              pRvlc->conceal_min = fMax(0, bnds - offset);
              return;
            }
            dpcm -= TABLE_OFFSET;
            if ((dpcm == MIN_RVL) || (dpcm == MAX_RVL)) {
              if ((pRvlc->length_of_rvlc_escapes) ||
                  (*pEscBwdCnt >= escEscCnt)) {
                pScfBwd[bnds] = noisenrg;
                pRvlc->conceal_min = fMax(0, bnds - offset);
                return;
              } else {
                if (dpcm == MIN_RVL) {
                  dpcm -= *pScfEsc--;
                } else {
                  dpcm += *pScfEsc--;
                }
                (*pEscBwdCnt)++;
                if (pRvlc->conceal_min_esc == CONCEAL_MIN_INIT) {
                  pRvlc->conceal_min_esc = fMax(0, bnds - offset);
                }
              }
            }
            pScfBwd[bnds] = noisenrg;
            noisenrg -= dpcm;
            pRvlc->firstNrg = noisenrg;
          }
          break;

        default:
          dpcm = decodeRVLCodeword(bs, pRvlc);
          if (dpcm < 0) {
            pScfBwd[bnds] = factor;
            pRvlc->conceal_min = fMax(0, bnds - offset);
            return;
          }
          dpcm -= TABLE_OFFSET;
          if ((dpcm == MIN_RVL) || (dpcm == MAX_RVL)) {
            if ((pRvlc->length_of_rvlc_escapes) || (*pEscBwdCnt >= escEscCnt)) {
              pScfBwd[bnds] = factor;
              pRvlc->conceal_min = fMax(0, bnds - offset);
              return;
            } else {
              if (dpcm == MIN_RVL) {
                dpcm -= *pScfEsc--;
              } else {
                dpcm += *pScfEsc--;
              }
              (*pEscBwdCnt)++;
              if (pRvlc->conceal_min_esc == CONCEAL_MIN_INIT) {
                pRvlc->conceal_min_esc = fMax(0, bnds - offset);
              }
            }
          }
          pScfBwd[bnds] = factor;
          factor -= dpcm;
          pRvlc->firstScf = factor;
          break;
      }
    }
  }
}

static void rvlcFinalErrorDetection(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo) {
  CErRvlcInfo *pRvlc =
      &pAacDecoderChannelInfo->pComData->overlay.aac.erRvlcInfo;
  UCHAR ErrorStatusComplete = 0;
  UCHAR ErrorStatusLengthFwd = 0;
  UCHAR ErrorStatusLengthBwd = 0;
  UCHAR ErrorStatusLengthEscapes = 0;
  UCHAR ErrorStatusFirstScf = 0;
  UCHAR ErrorStatusLastScf = 0;
  UCHAR ErrorStatusFirstNrg = 0;
  UCHAR ErrorStatusLastNrg = 0;
  UCHAR ErrorStatusFirstIs = 0;
  UCHAR ErrorStatusLastIs = 0;
  UCHAR ErrorStatusForbiddenCwFwd = 0;
  UCHAR ErrorStatusForbiddenCwBwd = 0;
  UCHAR ErrorStatusNumEscapesFwd = 0;
  UCHAR ErrorStatusNumEscapesBwd = 0;
  UCHAR ConcealStatus = 1;
  UCHAR currentBlockType;

  pAacDecoderChannelInfo->pDynData->specificTo.aac.rvlcCurrentScaleFactorOK = 1;

  if (pRvlc->errorLogRvlc & RVLC_ERROR_FORBIDDEN_CW_DETECTED_FWD)
    ErrorStatusForbiddenCwFwd = 1;

  if (pRvlc->errorLogRvlc & RVLC_ERROR_FORBIDDEN_CW_DETECTED_BWD)
    ErrorStatusForbiddenCwBwd = 1;

  if (pRvlc->length_of_rvlc_sf_fwd) ErrorStatusLengthFwd = 1;

  if (pRvlc->length_of_rvlc_sf_bwd) ErrorStatusLengthBwd = 1;

  if (pRvlc->sf_escapes_present)
    if (pRvlc->length_of_rvlc_escapes) ErrorStatusLengthEscapes = 1;

  if (pRvlc->sf_used) {
    if (pRvlc->firstScf !=
        (pAacDecoderChannelInfo->pDynData->RawDataInfo.GlobalGain - SF_OFFSET))
      ErrorStatusFirstScf = 1;

    if (pRvlc->lastScf != (pRvlc->rev_global_gain - SF_OFFSET))
      ErrorStatusLastScf = 1;
  }

  if (pRvlc->noise_used) {
    if (pRvlc->firstNrg !=
        (pAacDecoderChannelInfo->pDynData->RawDataInfo.GlobalGain +
         pRvlc->dpcm_noise_nrg - SF_OFFSET - 90 - 256))
      ErrorStatusFirstNrg = 1;

    if (pRvlc->lastNrg !=
        (pRvlc->rev_global_gain + pRvlc->dpcm_noise_last_position - SF_OFFSET -
         90 - 256))
      ErrorStatusLastNrg = 1;
  }

  if (pRvlc->intensity_used) {
    if (pRvlc->firstIs != (-SF_OFFSET)) ErrorStatusFirstIs = 1;

    if (pRvlc->lastIs != (pRvlc->dpcm_is_last_position - SF_OFFSET))
      ErrorStatusLastIs = 1;
  }

  if ((pRvlc->numDecodedEscapeWordsFwd != pRvlc->numDecodedEscapeWordsEsc) &&
      (pRvlc->conceal_max == CONCEAL_MAX_INIT)) {
    ErrorStatusNumEscapesFwd = 1;
  }

  if ((pRvlc->numDecodedEscapeWordsBwd != pRvlc->numDecodedEscapeWordsEsc) &&
      (pRvlc->conceal_min == CONCEAL_MIN_INIT)) {
    ErrorStatusNumEscapesBwd = 1;
  }

  if (ErrorStatusLengthEscapes ||
      (((pRvlc->conceal_max == CONCEAL_MAX_INIT) &&
        (pRvlc->numDecodedEscapeWordsFwd != pRvlc->numDecodedEscapeWordsEsc) &&
        (ErrorStatusLastScf || ErrorStatusLastNrg || ErrorStatusLastIs))

       &&

       ((pRvlc->conceal_min == CONCEAL_MIN_INIT) &&
        (pRvlc->numDecodedEscapeWordsBwd != pRvlc->numDecodedEscapeWordsEsc) &&
        (ErrorStatusFirstScf || ErrorStatusFirstNrg || ErrorStatusFirstIs))) ||
      ((pRvlc->conceal_max == CONCEAL_MAX_INIT) &&
       ((pRvlc->rev_global_gain - SF_OFFSET - pRvlc->lastScf) < -15)) ||
      ((pRvlc->conceal_min == CONCEAL_MIN_INIT) &&
       ((pAacDecoderChannelInfo->pDynData->RawDataInfo.GlobalGain - SF_OFFSET -
         pRvlc->firstScf) < -15))) {
    if ((pRvlc->conceal_max == CONCEAL_MAX_INIT) ||
        (pRvlc->conceal_min == CONCEAL_MIN_INIT)) {
      pRvlc->conceal_max = 0;
      pRvlc->conceal_min = fMax(
          0, (pRvlc->numWindowGroups - 1) * 16 + pRvlc->maxSfbTransmitted - 1);
    } else {
      pRvlc->conceal_max = fMin(pRvlc->conceal_max, pRvlc->conceal_max_esc);
      pRvlc->conceal_min = fMax(pRvlc->conceal_min, pRvlc->conceal_min_esc);
    }
  }

  ErrorStatusComplete = ErrorStatusLastScf || ErrorStatusFirstScf ||
                        ErrorStatusLastNrg || ErrorStatusFirstNrg ||
                        ErrorStatusLastIs || ErrorStatusFirstIs ||
                        ErrorStatusForbiddenCwFwd ||
                        ErrorStatusForbiddenCwBwd || ErrorStatusLengthFwd ||
                        ErrorStatusLengthBwd || ErrorStatusLengthEscapes ||
                        ErrorStatusNumEscapesFwd || ErrorStatusNumEscapesBwd;

  currentBlockType =
      (GetWindowSequence(&pAacDecoderChannelInfo->icsInfo) == BLOCK_SHORT) ? 0
                                                                           : 1;

  if (!ErrorStatusComplete) {
    int band;
    int group;
    int bnds;
    int lastSfbIndex;

    lastSfbIndex = (pRvlc->numWindowGroups > 1) ? 16 : 64;

    for (group = 0; group < pRvlc->numWindowGroups; group++) {
      for (band = 0; band < pRvlc->maxSfbTransmitted; band++) {
        bnds = 16 * group + band;
        pAacDecoderChannelInfo->pDynData->aScaleFactor[bnds] =
            pAacDecoderStaticChannelInfo->concealmentInfo
                .aRvlcPreviousScaleFactor[bnds] =
                pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfFwd[bnds];
      }
    }

    for (group = 0; group < pRvlc->numWindowGroups; group++) {
      for (band = 0; band < pRvlc->maxSfbTransmitted; band++) {
        bnds = 16 * group + band;
        pAacDecoderStaticChannelInfo->concealmentInfo
            .aRvlcPreviousCodebook[bnds] =
            pAacDecoderChannelInfo->pDynData->aCodeBook[bnds];
      }
      for (; band < lastSfbIndex; band++) {
        bnds = 16 * group + band;
        FDK_ASSERT(bnds >= 0 && bnds < RVLC_MAX_SFB);
        pAacDecoderStaticChannelInfo->concealmentInfo
            .aRvlcPreviousCodebook[bnds] = ZERO_HCB;
      }
    }
  } else {
    int band;
    int group;

    if (((pRvlc->conceal_min != CONCEAL_MIN_INIT) ||
         (pRvlc->conceal_max != CONCEAL_MAX_INIT)) &&
        (pRvlc->conceal_min <= pRvlc->conceal_max) &&
        (pAacDecoderStaticChannelInfo->concealmentInfo.rvlcPreviousBlockType ==
         currentBlockType) &&
        pAacDecoderStaticChannelInfo->concealmentInfo
            .rvlcPreviousScaleFactorOK &&
        pRvlc->sf_concealment && ConcealStatus) {
      BidirectionalEstimation_UseScfOfPrevFrameAsReference(
          pAacDecoderChannelInfo, pAacDecoderStaticChannelInfo);
      ConcealStatus = 0;
    }

    if ((pRvlc->conceal_min <= pRvlc->conceal_max) &&
        ((pRvlc->conceal_min != CONCEAL_MIN_INIT) ||
         (pRvlc->conceal_max != CONCEAL_MAX_INIT)) &&
        !(pAacDecoderStaticChannelInfo->concealmentInfo
              .rvlcPreviousScaleFactorOK &&
          pRvlc->sf_concealment &&
          (pAacDecoderStaticChannelInfo->concealmentInfo
               .rvlcPreviousBlockType == currentBlockType)) &&
        ConcealStatus) {
      BidirectionalEstimation_UseLowerScfOfCurrentFrame(pAacDecoderChannelInfo);
      ConcealStatus = 0;
    }

    if ((pRvlc->conceal_min <= pRvlc->conceal_max) &&
        ((ErrorStatusLastScf && ErrorStatusFirstScf) ||
         (ErrorStatusLastNrg && ErrorStatusFirstNrg) ||
         (ErrorStatusLastIs && ErrorStatusFirstIs)) &&
        !(ErrorStatusForbiddenCwFwd || ErrorStatusForbiddenCwBwd ||
          ErrorStatusLengthEscapes) &&
        ConcealStatus) {
      StatisticalEstimation(pAacDecoderChannelInfo);
      ConcealStatus = 0;
    }

    if ((pRvlc->conceal_min <= pRvlc->conceal_max) &&
        pAacDecoderStaticChannelInfo->concealmentInfo
            .rvlcPreviousScaleFactorOK &&
        pRvlc->sf_concealment &&
        (pAacDecoderStaticChannelInfo->concealmentInfo.rvlcPreviousBlockType ==
         currentBlockType) &&
        ConcealStatus) {
      PredictiveInterpolation(pAacDecoderChannelInfo,
                              pAacDecoderStaticChannelInfo);
      ConcealStatus = 0;
    }

    if (ConcealStatus) {
      for (group = 0; group < pRvlc->numWindowGroups; group++) {
        for (band = 0; band < pRvlc->maxSfbTransmitted; band++) {
          pAacDecoderChannelInfo->pDynData->aScaleFactor[16 * group + band] = 0;
        }
      }
      pAacDecoderChannelInfo->pDynData->specificTo.aac
          .rvlcCurrentScaleFactorOK = 0;
    }
  }
}

void CRvlc_Read(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                HANDLE_FDK_BITSTREAM bs) {
  CErRvlcInfo *pRvlc =
      &pAacDecoderChannelInfo->pComData->overlay.aac.erRvlcInfo;

  int group, band;

  pRvlc->numWindowGroups = GetWindowGroups(&pAacDecoderChannelInfo->icsInfo);
  pRvlc->maxSfbTransmitted =
      GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo->icsInfo);
  pRvlc->noise_used = 0;
  pRvlc->dpcm_noise_nrg = 0;
  pRvlc->dpcm_noise_last_position = 0;
  pRvlc->length_of_rvlc_escapes =
      -1;

  pRvlc->sf_concealment = FDKreadBits(bs, 1);
  pRvlc->rev_global_gain = FDKreadBits(bs, 8);

  if (GetWindowSequence(&pAacDecoderChannelInfo->icsInfo) == BLOCK_SHORT) {
    pRvlc->length_of_rvlc_sf = FDKreadBits(bs, 11);
  } else {
    pRvlc->length_of_rvlc_sf = FDKreadBits(bs, 9);
  }

  for (group = 0; group < pRvlc->numWindowGroups; group++) {
    for (band = 0; band < pRvlc->maxSfbTransmitted; band++) {
      if (pAacDecoderChannelInfo->pDynData->aCodeBook[16 * group + band] ==
          NOISE_HCB) {
        pRvlc->noise_used = 1;
        break;
      }
    }
  }

  if (pRvlc->noise_used)
    pRvlc->dpcm_noise_nrg = FDKreadBits(bs, 9);

  pRvlc->sf_escapes_present = FDKreadBits(bs, 1);

  if (pRvlc->sf_escapes_present) {
    pRvlc->length_of_rvlc_escapes = FDKreadBits(bs, 8);
  }

  if (pRvlc->noise_used) {
    pRvlc->dpcm_noise_last_position = FDKreadBits(bs, 9);
    pRvlc->length_of_rvlc_sf -= 9;
  }

  pRvlc->length_of_rvlc_sf_fwd = pRvlc->length_of_rvlc_sf;
  pRvlc->length_of_rvlc_sf_bwd = pRvlc->length_of_rvlc_sf;
}

void CRvlc_Decode(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                  CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
                  HANDLE_FDK_BITSTREAM bs) {
  CErRvlcInfo *pRvlc =
      &pAacDecoderChannelInfo->pComData->overlay.aac.erRvlcInfo;
  INT bitCntOffst;
  INT saveBitCnt;

  rvlcInit(pRvlc, pAacDecoderChannelInfo, bs);

  saveBitCnt = (INT)FDKgetValidBits(bs);

  if (pRvlc->sf_escapes_present)
    rvlcDecodeEscapes(
        pRvlc, pAacDecoderChannelInfo->pComData->overlay.aac.aRvlcScfEsc, bs);

  rvlcDecodeForward(pRvlc, pAacDecoderChannelInfo, bs);
  rvlcDecodeBackward(pRvlc, pAacDecoderChannelInfo, bs);
  rvlcFinalErrorDetection(pAacDecoderChannelInfo, pAacDecoderStaticChannelInfo);

  pAacDecoderChannelInfo->pDynData->specificTo.aac.rvlcIntensityUsed =
      pRvlc->intensity_used;
  pAacDecoderChannelInfo->data.aac.PnsData.PnsActive = pRvlc->noise_used;

  bitCntOffst = (INT)FDKgetValidBits(bs) - saveBitCnt;
  if (bitCntOffst) {
    FDKpushBiDirectional(bs, bitCntOffst);
  }
}

void CRvlc_ElementCheck(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo[],
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[],
    const UINT flags, const INT elChannels) {
  int ch;

  if (pAacDecoderStaticChannelInfo == NULL) {
    return;
  }

  if ((flags & AC_ER_RVLC) && (elChannels == 2)) {
    if (((pAacDecoderChannelInfo[0]
              ->pDynData->specificTo.aac.rvlcCurrentScaleFactorOK == 0) ||
         (pAacDecoderChannelInfo[1]
              ->pDynData->specificTo.aac.rvlcCurrentScaleFactorOK == 0)) &&
        pAacDecoderChannelInfo[0]->pComData->jointStereoData.MsMaskPresent) {
      pAacDecoderChannelInfo[0]
          ->pDynData->specificTo.aac.rvlcCurrentScaleFactorOK = 0;
      pAacDecoderChannelInfo[1]
          ->pDynData->specificTo.aac.rvlcCurrentScaleFactorOK = 0;
    }

    if ((pAacDecoderChannelInfo[0]
             ->pDynData->specificTo.aac.rvlcCurrentScaleFactorOK == 0) &&
        (pAacDecoderChannelInfo[1]
             ->pDynData->specificTo.aac.rvlcCurrentScaleFactorOK == 1) &&
        (pAacDecoderChannelInfo[1]
             ->pDynData->specificTo.aac.rvlcIntensityUsed == 1)) {
      pAacDecoderChannelInfo[1]
          ->pDynData->specificTo.aac.rvlcCurrentScaleFactorOK = 0;
    }
  }

  for (ch = 0; ch < elChannels; ch++) {
    pAacDecoderStaticChannelInfo[ch]->concealmentInfo.rvlcPreviousBlockType =
        (GetWindowSequence(&pAacDecoderChannelInfo[ch]->icsInfo) == BLOCK_SHORT)
            ? 0
            : 1;
    if (flags & AC_ER_RVLC) {
      pAacDecoderStaticChannelInfo[ch]
          ->concealmentInfo.rvlcPreviousScaleFactorOK =
          pAacDecoderChannelInfo[ch]
              ->pDynData->specificTo.aac.rvlcCurrentScaleFactorOK;
    } else {
      pAacDecoderStaticChannelInfo[ch]
          ->concealmentInfo.rvlcPreviousScaleFactorOK = 0;
    }
  }
}
