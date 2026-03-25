/**************************** AAC encoder library ******************************

   Author(s):   A. Horndasch (code originally from lwr) / Josef Hoepfl (FDK)

   Description: intensity stereo processing

*******************************************************************************/

#include "intensity.h"
#include "interface.h"
#include "psy_configuration.h"
#include "psy_const.h"
#include "qc_main.h"
#include "bit_cnt.h"

#define IS_CORR_THRESH FL2FXCONST_DBL(0.95f)

#define IS_TOTAL_ERROR_THRESH FL2FXCONST_DBL(0.04f)
#define IS_LOCAL_ERROR_THRESH FL2FXCONST_DBL(0.01f)

#define IS_DIRECTION_DEVIATION_THRESH_SF 2
#define IS_DIRECTION_DEVIATION_THRESH \
  FL2FXCONST_DBL(2.0f / (1 << IS_DIRECTION_DEVIATION_THRESH_SF))

#define IS_REGION_MIN_LOUDNESS FL2FXCONST_DBL(0.1f)
#define IS_MIN_SFBS 6
#define IS_LEFT_RIGHT_RATIO_THRESH FL2FXCONST_DBL(0.7f)

#define REAL_SCALE_SF 1
#define OVERALL_LOUDNESS_SF 6
#define MAX_SFB_PER_GROUP_SF 6
#define MDCT_SPEC_SF 6

typedef struct {
  FIXP_DBL corr_thresh;

  FIXP_DBL total_error_thresh;

  FIXP_DBL local_error_thresh;

  FIXP_DBL direction_deviation_thresh;

  FIXP_DBL is_region_min_loudness;

  INT min_is_sfbs;

  FIXP_DBL left_right_ratio_threshold;

} INTENSITY_PARAMETERS;

static INT calcSfbMaxScale(const FIXP_DBL *mdctSpectrum, const INT l1,
                           const INT l2) {
  INT i;
  INT sfbMaxScale;
  FIXP_DBL maxSpc;

  maxSpc = FL2FXCONST_DBL(0.0);
  for (i = l1; i < l2; i++) {
    FIXP_DBL tmp = fixp_abs((FIXP_DBL)mdctSpectrum[i]);
    maxSpc = fixMax(maxSpc, tmp);
  }
  sfbMaxScale = (maxSpc == FL2FXCONST_DBL(0.0)) ? (DFRACT_BITS - 2)
                                                : CntLeadingZeros(maxSpc) - 1;

  return sfbMaxScale;
}

static void FDKaacEnc_initIsParams(INTENSITY_PARAMETERS *isParams) {
  isParams->corr_thresh = IS_CORR_THRESH;
  isParams->total_error_thresh = IS_TOTAL_ERROR_THRESH;
  isParams->local_error_thresh = IS_LOCAL_ERROR_THRESH;
  isParams->direction_deviation_thresh = IS_DIRECTION_DEVIATION_THRESH;
  isParams->is_region_min_loudness = IS_REGION_MIN_LOUDNESS;
  isParams->min_is_sfbs = IS_MIN_SFBS;
  isParams->left_right_ratio_threshold = IS_LEFT_RIGHT_RATIO_THRESH;
}

static void FDKaacEnc_prepareIntensityDecision(
    const FIXP_DBL *sfbEnergyLeft, const FIXP_DBL *sfbEnergyRight,
    const FIXP_DBL *sfbEnergyLdDataLeft, const FIXP_DBL *sfbEnergyLdDataRight,
    const FIXP_DBL *mdctSpectrumLeft, const FIXP_DBL *mdctSpectrumRight,
    const INTENSITY_PARAMETERS *isParams, FIXP_DBL *hrrErr, INT *isMask,
    FIXP_DBL *realScale, FIXP_DBL *normSfbLoudness, const INT sfbCnt,
    const INT sfbPerGroup, const INT maxSfbPerGroup, const INT *sfbOffset) {
  INT j, sfb, sfboffs;
  INT grpCounter;

  FIXP_DBL overallLoudness[MAX_NO_OF_GROUPS];
  FIXP_DBL channelCorr[MAX_GROUPED_SFB];
  FIXP_DBL ml, mr;
  FIXP_DBL prod_lr;
  FIXP_DBL square_l, square_r;
  FIXP_DBL tmp_l, tmp_r;
  FIXP_DBL inv_n;

  FDKmemclear(channelCorr, MAX_GROUPED_SFB * sizeof(FIXP_DBL));
  FDKmemclear(normSfbLoudness, MAX_GROUPED_SFB * sizeof(FIXP_DBL));
  FDKmemclear(overallLoudness, MAX_NO_OF_GROUPS * sizeof(FIXP_DBL));
  FDKmemclear(realScale, MAX_GROUPED_SFB * sizeof(FIXP_DBL));

  for (grpCounter = 0, sfboffs = 0; sfboffs < sfbCnt;
       sfboffs += sfbPerGroup, grpCounter++) {
    overallLoudness[grpCounter] = FL2FXCONST_DBL(0.0f);
    for (sfb = 0; sfb < maxSfbPerGroup; sfb++) {
      INT sL, sR, s;
      FIXP_DBL isValue = sfbEnergyLdDataLeft[sfb + sfboffs] -
                         sfbEnergyLdDataRight[sfb + sfboffs];

      realScale[sfb + sfboffs] = fixMin(
          FL2FXCONST_DBL(60.f / (1 << (REAL_SCALE_SF + LD_DATA_SHIFT))),
          fixMax(FL2FXCONST_DBL(-60.f / (1 << (REAL_SCALE_SF + LD_DATA_SHIFT))),
                 isValue));

      sL = fixMax(0, (CntLeadingZeros(sfbEnergyLeft[sfb + sfboffs]) - 1));
      sR = fixMax(0, (CntLeadingZeros(sfbEnergyRight[sfb + sfboffs]) - 1));
      s = (fixMin(sL, sR) >> 2) << 2;
      normSfbLoudness[sfb + sfboffs] =
          sqrtFixp(sqrtFixp(((sfbEnergyLeft[sfb + sfboffs] << s) >> 1) +
                            ((sfbEnergyRight[sfb + sfboffs] << s) >> 1))) >>
          (s >> 2);

      overallLoudness[grpCounter] +=
          normSfbLoudness[sfb + sfboffs] >> OVERALL_LOUDNESS_SF;
      if ((sfbEnergyLeft[sfb + sfboffs] >=
           fMult(isParams->left_right_ratio_threshold,
                 sfbEnergyRight[sfb + sfboffs])) &&
          (fMult(isParams->left_right_ratio_threshold,
                 sfbEnergyLeft[sfb + sfboffs]) <=
           sfbEnergyRight[sfb + sfboffs])) {
        hrrErr[sfb + sfboffs] = FL2FXCONST_DBL(1.0 / 8.0);
      }
    }
  }

  for (grpCounter = 0, sfboffs = 0; sfboffs < sfbCnt;
       sfboffs += sfbPerGroup, grpCounter++) {
    INT invOverallLoudnessSF;
    FIXP_DBL invOverallLoudness;

    if (overallLoudness[grpCounter] == FL2FXCONST_DBL(0.0)) {
      invOverallLoudness = FL2FXCONST_DBL(0.0);
      invOverallLoudnessSF = 0;
    } else {
      invOverallLoudness =
          fDivNorm((FIXP_DBL)MAXVAL_DBL, overallLoudness[grpCounter],
                   &invOverallLoudnessSF);
      invOverallLoudnessSF =
          invOverallLoudnessSF - OVERALL_LOUDNESS_SF +
          1;
    }
    invOverallLoudnessSF = fixMin(
        fixMax(invOverallLoudnessSF, -(DFRACT_BITS - 1)), DFRACT_BITS - 1);

    for (sfb = 0; sfb < maxSfbPerGroup; sfb++) {
      FIXP_DBL tmp;

      tmp = fMultDiv2((normSfbLoudness[sfb + sfboffs] >> OVERALL_LOUDNESS_SF)
                          << OVERALL_LOUDNESS_SF,
                      invOverallLoudness);

      normSfbLoudness[sfb + sfboffs] = scaleValue(tmp, invOverallLoudnessSF);

      channelCorr[sfb + sfboffs] = FL2FXCONST_DBL(0.0f);

      inv_n = GetInvInt(
          (sfbOffset[sfb + sfboffs + 1] - sfbOffset[sfb + sfboffs]) >> 1);

      if (inv_n > FL2FXCONST_DBL(0.0f)) {
        INT s, sL, sR;

        ml = FL2FXCONST_DBL(0.0f);
        mr = FL2FXCONST_DBL(0.0f);
        prod_lr = FL2FXCONST_DBL(0.0f);
        square_l = FL2FXCONST_DBL(0.0f);
        square_r = FL2FXCONST_DBL(0.0f);

        sL = calcSfbMaxScale(mdctSpectrumLeft, sfbOffset[sfb + sfboffs],
                             sfbOffset[sfb + sfboffs + 1]);
        sR = calcSfbMaxScale(mdctSpectrumRight, sfbOffset[sfb + sfboffs],
                             sfbOffset[sfb + sfboffs + 1]);
        s = fixMin(sL, sR);

        for (j = sfbOffset[sfb + sfboffs]; j < sfbOffset[sfb + sfboffs + 1];
             j++) {
          ml += fMultDiv2((mdctSpectrumLeft[j] << s),
                          inv_n);
          mr += fMultDiv2((mdctSpectrumRight[j] << s),
                          inv_n);
        }
        ml = fMultDiv2(ml, inv_n);
        mr = fMultDiv2(mr, inv_n);

        for (j = sfbOffset[sfb + sfboffs]; j < sfbOffset[sfb + sfboffs + 1];
             j++) {
          tmp_l = fMultDiv2((mdctSpectrumLeft[j] << s), inv_n) -
                  ml;
          tmp_r = fMultDiv2((mdctSpectrumRight[j] << s), inv_n) -
                  mr;

          prod_lr += fMultDiv2(
              tmp_l, tmp_r);
          square_l +=
              fPow2Div2(tmp_l);
          square_r +=
              fPow2Div2(tmp_r);
        }
        prod_lr = prod_lr << 1;
        square_l = square_l << 1;
        square_r = square_r << 1;

        if (square_l > FL2FXCONST_DBL(0.0f) &&
            square_r > FL2FXCONST_DBL(0.0f)) {
          INT channelCorrSF = 0;

          sL = fixMax(0, (CntLeadingZeros(square_l) - 1));
          sR = fixMax(0, (CntLeadingZeros(square_r) - 1));
          s = ((sL + sR) >> 1) << 1;
          sL = fixMin(sL, s);
          sR = s - sL;
          tmp = fMult(square_l << sL, square_r << sR);
          tmp = sqrtFixp(tmp);

          FDK_ASSERT(tmp > FL2FXCONST_DBL(0.0f));

          if (prod_lr < FL2FXCONST_DBL(0.0f)) {
            channelCorr[sfb + sfboffs] =
                -(fDivNorm(-prod_lr, tmp, &channelCorrSF));

          } else {
            channelCorr[sfb + sfboffs] =
                (fDivNorm(prod_lr, tmp, &channelCorrSF));
          }
          channelCorrSF = fixMin(
              fixMax((channelCorrSF + ((sL + sR) >> 1)), -(DFRACT_BITS - 1)),
              DFRACT_BITS - 1);

          if (channelCorrSF < 0) {
            channelCorr[sfb + sfboffs] =
                channelCorr[sfb + sfboffs] >> (-channelCorrSF);
          } else {
            if (fAbs(channelCorr[sfb + sfboffs]) >
                (((FIXP_DBL)MAXVAL_DBL) >> channelCorrSF)) {
              if (channelCorr[sfb + sfboffs] < FL2FXCONST_DBL(0.0f))
                channelCorr[sfb + sfboffs] = -(FIXP_DBL)MAXVAL_DBL;
              else
                channelCorr[sfb + sfboffs] = (FIXP_DBL)MAXVAL_DBL;
            } else {
              channelCorr[sfb + sfboffs] = channelCorr[sfb + sfboffs]
                                           << channelCorrSF;
            }
          }
        }
      }

      if (hrrErr[sfb + sfboffs] == FL2FXCONST_DBL(1.0 / 8.0)) {
        continue;
      }

      hrrErr[sfb + sfboffs] =
          fMultDiv2((FL2FXCONST_DBL(0.25f) - (channelCorr[sfb + sfboffs] >> 2)),
                    normSfbLoudness[sfb + sfboffs]);

      if (fAbs(channelCorr[sfb + sfboffs]) >= isParams->corr_thresh) {
        isMask[sfb + sfboffs] = 1;
      }
    }
  }
}

static void FDKaacEnc_finalizeIntensityDecision(
    const FIXP_DBL *hrrErr, INT *isMask, const FIXP_DBL *realIsScale,
    const FIXP_DBL *normSfbLoudness, const INTENSITY_PARAMETERS *isParams,
    const INT sfbCnt, const INT sfbPerGroup, const INT maxSfbPerGroup) {
  INT sfb, sfboffs, j;
  FIXP_DBL isScaleLast = FL2FXCONST_DBL(0.0f);
  INT isStartValueFound = 0;

  for (sfboffs = 0; sfboffs < sfbCnt; sfboffs += sfbPerGroup) {
    INT startIsSfb = 0;
    INT inIsBlock = 0;
    INT currentIsSfbCount = 0;
    FIXP_DBL overallHrrError = FL2FXCONST_DBL(0.0f);
    FIXP_DBL isRegionLoudness = FL2FXCONST_DBL(0.0f);

    for (sfb = 0; sfb < maxSfbPerGroup; sfb++) {
      if (isMask[sfboffs + sfb] == 1) {
        if (currentIsSfbCount == 0) {
          startIsSfb = sfboffs + sfb;
        }
        if (isStartValueFound == 0) {
          isScaleLast = realIsScale[sfboffs + sfb];
          isStartValueFound = 1;
        }
        inIsBlock = 1;
        currentIsSfbCount++;
        overallHrrError += hrrErr[sfboffs + sfb] >> (MAX_SFB_PER_GROUP_SF - 3);
        isRegionLoudness +=
            normSfbLoudness[sfboffs + sfb] >> MAX_SFB_PER_GROUP_SF;
      } else {
        if (inIsBlock) {
          overallHrrError +=
              hrrErr[sfboffs + sfb] >> (MAX_SFB_PER_GROUP_SF - 3);
          isRegionLoudness +=
              normSfbLoudness[sfboffs + sfb] >> MAX_SFB_PER_GROUP_SF;

          if ((hrrErr[sfboffs + sfb] < (isParams->local_error_thresh >> 3)) &&
              (overallHrrError <
               (isParams->total_error_thresh >> MAX_SFB_PER_GROUP_SF))) {
            currentIsSfbCount++;
            isMask[sfboffs + sfb] = 1;
          } else {
            inIsBlock = 0;
          }
        }
      }
      if (inIsBlock) {
        if (fAbs(isScaleLast - realIsScale[sfboffs + sfb]) <
            (isParams->direction_deviation_thresh >>
             (REAL_SCALE_SF + LD_DATA_SHIFT -
              IS_DIRECTION_DEVIATION_THRESH_SF))) {
          isScaleLast = realIsScale[sfboffs + sfb];
        } else {
          isMask[sfboffs + sfb] = 0;
          inIsBlock = 0;
          currentIsSfbCount--;
        }
      }

      if (currentIsSfbCount > 0 && (!inIsBlock || sfb == maxSfbPerGroup - 1)) {
        if (currentIsSfbCount < isParams->min_is_sfbs ||
            (isRegionLoudness<isParams->is_region_min_loudness>>
             MAX_SFB_PER_GROUP_SF)) {
          for (j = startIsSfb; j <= sfboffs + sfb; j++) {
            isMask[j] = 0;
          }
          isScaleLast = FL2FXCONST_DBL(0.0f);
          isStartValueFound = 0;
          for (j = 0; j < startIsSfb; j++) {
            if (isMask[j] != 0) {
              isScaleLast = realIsScale[j];
              isStartValueFound = 1;
            }
          }
        }
        currentIsSfbCount = 0;
        overallHrrError = FL2FXCONST_DBL(0.0f);
        isRegionLoudness = FL2FXCONST_DBL(0.0f);
      }
    }
  }
}

void FDKaacEnc_IntensityStereoProcessing(
    FIXP_DBL *sfbEnergyLeft, FIXP_DBL *sfbEnergyRight,
    FIXP_DBL *mdctSpectrumLeft, FIXP_DBL *mdctSpectrumRight,
    FIXP_DBL *sfbThresholdLeft, FIXP_DBL *sfbThresholdRight,
    FIXP_DBL *sfbThresholdLdDataRight, FIXP_DBL *sfbSpreadEnLeft,
    FIXP_DBL *sfbSpreadEnRight, FIXP_DBL *sfbEnergyLdDataLeft,
    FIXP_DBL *sfbEnergyLdDataRight, INT *msDigest, INT *msMask,
    const INT sfbCnt, const INT sfbPerGroup, const INT maxSfbPerGroup,
    const INT *sfbOffset, const INT allowIS, INT *isBook, INT *isScale,
    PNS_DATA *RESTRICT pnsData[2]) {
  INT sfb, sfboffs, j;
  FIXP_DBL scale;
  FIXP_DBL lr;
  FIXP_DBL hrrErr[MAX_GROUPED_SFB];
  FIXP_DBL normSfbLoudness[MAX_GROUPED_SFB];
  FIXP_DBL realIsScale[MAX_GROUPED_SFB];
  INTENSITY_PARAMETERS isParams;
  INT isMask[MAX_GROUPED_SFB];

  FDKmemclear((void *)isBook, sfbCnt * sizeof(INT));
  FDKmemclear((void *)isMask, sfbCnt * sizeof(INT));
  FDKmemclear((void *)realIsScale, sfbCnt * sizeof(FIXP_DBL));
  FDKmemclear((void *)isScale, sfbCnt * sizeof(INT));
  FDKmemclear((void *)hrrErr, sfbCnt * sizeof(FIXP_DBL));

  if (!allowIS) return;

  FDKaacEnc_initIsParams(&isParams);

  FDKaacEnc_prepareIntensityDecision(
      sfbEnergyLeft, sfbEnergyRight, sfbEnergyLdDataLeft, sfbEnergyLdDataRight,
      mdctSpectrumLeft, mdctSpectrumRight, &isParams, hrrErr, isMask,
      realIsScale, normSfbLoudness, sfbCnt, sfbPerGroup, maxSfbPerGroup,
      sfbOffset);

  FDKaacEnc_finalizeIntensityDecision(hrrErr, isMask, realIsScale,
                                      normSfbLoudness, &isParams, sfbCnt,
                                      sfbPerGroup, maxSfbPerGroup);

  for (sfb = 0; sfb < sfbCnt; sfb += sfbPerGroup) {
    for (sfboffs = 0; sfboffs < maxSfbPerGroup; sfboffs++) {
      INT sL, sR;
      FIXP_DBL inv_n;
      INT mdct_spec_sf = MDCT_SPEC_SF;

      msMask[sfb + sfboffs] = 0;
      if (isMask[sfb + sfboffs] == 0) {
        continue;
      }

      if ((sfbEnergyLeft[sfb + sfboffs] < sfbThresholdLeft[sfb + sfboffs]) &&
          (fMult(FL2FXCONST_DBL(1.0f / 1.5f), sfbEnergyRight[sfb + sfboffs]) >
           sfbThresholdRight[sfb + sfboffs])) {
        continue;
      }
      if (pnsData[0]) {
        if (pnsData[0]->pnsFlag[sfb + sfboffs]) {
          pnsData[0]->pnsFlag[sfb + sfboffs] = 0;
        }
        if (pnsData[1]->pnsFlag[sfb + sfboffs]) {
          pnsData[1]->pnsFlag[sfb + sfboffs] = 0;
        }
      }

      if (sfbOffset[sfb + sfboffs + 1] - sfbOffset[sfb + sfboffs] >
          1 << mdct_spec_sf) {
        mdct_spec_sf++;
      }

      inv_n = GetInvInt(
          (sfbOffset[sfb + sfboffs + 1] - sfbOffset[sfb + sfboffs]) >>
          1);
      sL = calcSfbMaxScale(mdctSpectrumLeft, sfbOffset[sfb + sfboffs],
                           sfbOffset[sfb + sfboffs + 1]);
      sR = calcSfbMaxScale(mdctSpectrumRight, sfbOffset[sfb + sfboffs],
                           sfbOffset[sfb + sfboffs + 1]);

      lr = FL2FXCONST_DBL(0.0f);
      for (j = sfbOffset[sfb + sfboffs]; j < sfbOffset[sfb + sfboffs + 1]; j++)
        lr += fMultDiv2(
            fMultDiv2(mdctSpectrumLeft[j] << sL, mdctSpectrumRight[j] << sR),
            inv_n);
      lr = lr << 1;

      if (lr < FL2FXCONST_DBL(0.0f)) {
        INT s0, s1, s2;
        FIXP_DBL tmp, d, ed = FL2FXCONST_DBL(0.0f);

        s0 = fixMin(sL, sR);
        for (j = sfbOffset[sfb + sfboffs]; j < sfbOffset[sfb + sfboffs + 1];
             j++) {
          d = ((mdctSpectrumLeft[j] << s0) >> 1) -
              ((mdctSpectrumRight[j] << s0) >> 1);
          ed += fMultDiv2(d, d) >> (mdct_spec_sf - 1);
        }
        msMask[sfb + sfboffs] = 1;
        tmp = fDivNorm(sfbEnergyLeft[sfb + sfboffs], ed, &s1);
        s2 = (s1) + (2 * s0) - 2 - mdct_spec_sf;
        if (s2 & 1) {
          tmp = tmp >> 1;
          s2 = s2 + 1;
        }
        s2 = (s2 >> 1) + 1;
        s2 = fixMin(fixMax(s2, -(DFRACT_BITS - 1)), (DFRACT_BITS - 1));
        scale = sqrtFixp(tmp);
        if (s2 < 0) {
          s2 = -s2;
          for (j = sfbOffset[sfb + sfboffs]; j < sfbOffset[sfb + sfboffs + 1];
               j++) {
            mdctSpectrumLeft[j] = (fMultDiv2(mdctSpectrumLeft[j], scale) -
                                   fMultDiv2(mdctSpectrumRight[j], scale)) >>
                                  s2;
            mdctSpectrumRight[j] = FL2FXCONST_DBL(0.0f);
          }
        } else {
          for (j = sfbOffset[sfb + sfboffs]; j < sfbOffset[sfb + sfboffs + 1];
               j++) {
            mdctSpectrumLeft[j] = (fMultDiv2(mdctSpectrumLeft[j], scale) -
                                   fMultDiv2(mdctSpectrumRight[j], scale))
                                  << s2;
            mdctSpectrumRight[j] = FL2FXCONST_DBL(0.0f);
          }
        }
      } else {
        INT s0, s1, s2;
        FIXP_DBL tmp, s, es = FL2FXCONST_DBL(0.0f);

        s0 = fixMin(sL, sR);
        for (j = sfbOffset[sfb + sfboffs]; j < sfbOffset[sfb + sfboffs + 1];
             j++) {
          s = ((mdctSpectrumLeft[j] << s0) >> 1) +
              ((mdctSpectrumRight[j] << s0) >> 1);
          es += fMultDiv2(s, s) >>
                (mdct_spec_sf -
                 1);
        }
        msMask[sfb + sfboffs] = 0;
        tmp = fDivNorm(sfbEnergyLeft[sfb + sfboffs], es, &s1);
        s2 = (s1) + (2 * s0) - 2 - mdct_spec_sf;
        if (s2 & 1) {
          tmp = tmp >> 1;
          s2 = s2 + 1;
        }
        s2 = (s2 >> 1) + 1;
        s2 = fixMin(fixMax(s2, -(DFRACT_BITS - 1)), (DFRACT_BITS - 1));
        scale = sqrtFixp(tmp);
        if (s2 < 0) {
          s2 = -s2;
          for (j = sfbOffset[sfb + sfboffs]; j < sfbOffset[sfb + sfboffs + 1];
               j++) {
            mdctSpectrumLeft[j] = (fMultDiv2(mdctSpectrumLeft[j], scale) +
                                   fMultDiv2(mdctSpectrumRight[j], scale)) >>
                                  s2;
            mdctSpectrumRight[j] = FL2FXCONST_DBL(0.0f);
          }
        } else {
          for (j = sfbOffset[sfb + sfboffs]; j < sfbOffset[sfb + sfboffs + 1];
               j++) {
            mdctSpectrumLeft[j] = (fMultDiv2(mdctSpectrumLeft[j], scale) +
                                   fMultDiv2(mdctSpectrumRight[j], scale))
                                  << s2;
            mdctSpectrumRight[j] = FL2FXCONST_DBL(0.0f);
          }
        }
      }

      isBook[sfb + sfboffs] = CODE_BOOK_IS_IN_PHASE_NO;

      if (realIsScale[sfb + sfboffs] < FL2FXCONST_DBL(0.0f)) {
        isScale[sfb + sfboffs] =
            (INT)(((realIsScale[sfb + sfboffs] >> 1) -
                   FL2FXCONST_DBL(
                       0.5f / (1 << (REAL_SCALE_SF + LD_DATA_SHIFT + 1)))) >>
                  (DFRACT_BITS - 1 - REAL_SCALE_SF - LD_DATA_SHIFT - 1)) +
            1;
      } else {
        isScale[sfb + sfboffs] =
            (INT)(((realIsScale[sfb + sfboffs] >> 1) +
                   FL2FXCONST_DBL(
                       0.5f / (1 << (REAL_SCALE_SF + LD_DATA_SHIFT + 1)))) >>
                  (DFRACT_BITS - 1 - REAL_SCALE_SF - LD_DATA_SHIFT - 1));
      }

      sfbEnergyRight[sfb + sfboffs] = FL2FXCONST_DBL(0.0f);
      sfbEnergyLdDataRight[sfb + sfboffs] = FL2FXCONST_DBL(-1.0f);
      sfbThresholdRight[sfb + sfboffs] = FL2FXCONST_DBL(0.0f);
      sfbThresholdLdDataRight[sfb + sfboffs] = FL2FXCONST_DBL(-0.515625f);
      sfbSpreadEnRight[sfb + sfboffs] = FL2FXCONST_DBL(0.0f);

      *msDigest = MS_SOME;
    }
  }
}
