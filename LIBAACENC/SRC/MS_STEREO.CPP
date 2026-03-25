/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: MS stereo processing

*******************************************************************************/

#include "ms_stereo.h"
#include "psy_const.h"

void FDKaacEnc_MsStereoProcessing(PSY_DATA *RESTRICT psyData[(2)],
                                  PSY_OUT_CHANNEL *psyOutChannel[2],
                                  const INT *isBook, INT *msDigest,
                                  INT *msMask,
                                  const INT allowMS, const INT sfbCnt,
                                  const INT sfbPerGroup,
                                  const INT maxSfbPerGroup,
                                  const INT *sfbOffset) {
  FIXP_DBL *sfbEnergyLeft =
      psyData[0]->sfbEnergy.Long;
  FIXP_DBL *sfbEnergyRight =
      psyData[1]->sfbEnergy.Long;
  const FIXP_DBL *sfbEnergyMid = psyData[0]->sfbEnergyMS.Long;
  const FIXP_DBL *sfbEnergySide = psyData[1]->sfbEnergyMS.Long;
  FIXP_DBL *sfbThresholdLeft =
      psyData[0]->sfbThreshold.Long;
  FIXP_DBL *sfbThresholdRight =
      psyData[1]->sfbThreshold.Long;

  FIXP_DBL *sfbSpreadEnLeft = psyData[0]->sfbSpreadEnergy.Long;
  FIXP_DBL *sfbSpreadEnRight = psyData[1]->sfbSpreadEnergy.Long;

  FIXP_DBL *sfbEnergyLeftLdData =
      psyOutChannel[0]->sfbEnergyLdData;
  FIXP_DBL *sfbEnergyRightLdData =
      psyOutChannel[1]->sfbEnergyLdData;
  FIXP_DBL *sfbEnergyMidLdData = psyData[0]->sfbEnergyMSLdData;
  FIXP_DBL *sfbEnergySideLdData = psyData[1]->sfbEnergyMSLdData;
  FIXP_DBL *sfbThresholdLeftLdData =
      psyOutChannel[0]->sfbThresholdLdData;
  FIXP_DBL *sfbThresholdRightLdData =
      psyOutChannel[1]->sfbThresholdLdData;

  FIXP_DBL *mdctSpectrumLeft =
      psyData[0]->mdctSpectrum;
  FIXP_DBL *mdctSpectrumRight =
      psyData[1]->mdctSpectrum;

  INT sfb, sfboffs, j;
  FIXP_DBL pnlrLdData, pnmsLdData;
  FIXP_DBL minThresholdLdData;
  FIXP_DBL minThreshold;
  INT useMS;

  INT msMaskTrueSomewhere = 0;
  INT numMsMaskFalse =
      0;

  for (sfb = 0; sfb < sfbCnt; sfb += sfbPerGroup) {
    for (sfboffs = 0; sfboffs < maxSfbPerGroup; sfboffs++) {
      if ((isBook == NULL) ? 1 : (isBook[sfb + sfboffs] == 0)) {
        FIXP_DBL tmp;

        minThresholdLdData = fixMin(sfbThresholdLeftLdData[sfb + sfboffs],
                                    sfbThresholdRightLdData[sfb + sfboffs]);

        tmp = fixMax(sfbEnergyLeftLdData[sfb + sfboffs],
                     sfbThresholdLeftLdData[sfb + sfboffs]);
        pnlrLdData = (sfbThresholdLeftLdData[sfb + sfboffs] >> 1) - (tmp >> 1);
        pnlrLdData = pnlrLdData + (sfbThresholdRightLdData[sfb + sfboffs] >> 1);
        tmp = fixMax(sfbEnergyRightLdData[sfb + sfboffs],
                     sfbThresholdRightLdData[sfb + sfboffs]);
        pnlrLdData = pnlrLdData - (tmp >> 1);

        tmp = fixMax(sfbEnergyMidLdData[sfb + sfboffs], minThresholdLdData);
        pnmsLdData = minThresholdLdData - (tmp >> 1);
        tmp = fixMax(sfbEnergySideLdData[sfb + sfboffs], minThresholdLdData);
        pnmsLdData = pnmsLdData - (tmp >> 1);
        useMS = ((allowMS != 0) && (pnmsLdData > pnlrLdData)) ? 1 : 0;

        if (useMS) {
          msMask[sfb + sfboffs] = 1;
          msMaskTrueSomewhere = 1;
          for (j = sfbOffset[sfb + sfboffs]; j < sfbOffset[sfb + sfboffs + 1];
               j++) {
            FIXP_DBL specL, specR;
            specL = mdctSpectrumLeft[j] >> 1;
            specR = mdctSpectrumRight[j] >> 1;
            mdctSpectrumLeft[j] = specL + specR;
            mdctSpectrumRight[j] = specL - specR;
          }
          minThreshold = fixMin(sfbThresholdLeft[sfb + sfboffs],
                                sfbThresholdRight[sfb + sfboffs]);
          sfbThresholdLeft[sfb + sfboffs] = sfbThresholdRight[sfb + sfboffs] =
              minThreshold;
          sfbThresholdLeftLdData[sfb + sfboffs] =
              sfbThresholdRightLdData[sfb + sfboffs] = minThresholdLdData;
          sfbEnergyLeft[sfb + sfboffs] = sfbEnergyMid[sfb + sfboffs];
          sfbEnergyRight[sfb + sfboffs] = sfbEnergySide[sfb + sfboffs];
          sfbEnergyLeftLdData[sfb + sfboffs] =
              sfbEnergyMidLdData[sfb + sfboffs];
          sfbEnergyRightLdData[sfb + sfboffs] =
              sfbEnergySideLdData[sfb + sfboffs];

          sfbSpreadEnLeft[sfb + sfboffs] = sfbSpreadEnRight[sfb + sfboffs] =
              fixMin(sfbSpreadEnLeft[sfb + sfboffs],
                     sfbSpreadEnRight[sfb + sfboffs]) >>
              1;

        } else {
          msMask[sfb + sfboffs] = 0;
          numMsMaskFalse++;
        }
      }
      else {
        if (msMask[sfb + sfboffs]) {
          msMaskTrueSomewhere = 1;
        }
        numMsMaskFalse = 9;
      }
    }
  }

  if (msMaskTrueSomewhere == 1) {
    if ((numMsMaskFalse == 0) ||
        ((numMsMaskFalse < maxSfbPerGroup) && (numMsMaskFalse < 9))) {
      *msDigest = SI_MS_MASK_ALL;
      for (sfb = 0; sfb < sfbCnt; sfb += sfbPerGroup) {
        for (sfboffs = 0; sfboffs < maxSfbPerGroup; sfboffs++) {
          if (((isBook == NULL) ? 1 : (isBook[sfb + sfboffs] == 0)) &&
              (msMask[sfb + sfboffs] == 0)) {
            msMask[sfb + sfboffs] = 1;
            for (j = sfbOffset[sfb + sfboffs]; j < sfbOffset[sfb + sfboffs + 1];
                 j++) {
              FIXP_DBL specL, specR;
              specL = mdctSpectrumLeft[j] >> 1;
              specR = mdctSpectrumRight[j] >> 1;
              mdctSpectrumLeft[j] = specL + specR;
              mdctSpectrumRight[j] = specL - specR;
            }
            minThreshold = fixMin(sfbThresholdLeft[sfb + sfboffs],
                                  sfbThresholdRight[sfb + sfboffs]);
            sfbThresholdLeft[sfb + sfboffs] = sfbThresholdRight[sfb + sfboffs] =
                minThreshold;
            minThresholdLdData = fixMin(sfbThresholdLeftLdData[sfb + sfboffs],
                                        sfbThresholdRightLdData[sfb + sfboffs]);
            sfbThresholdLeftLdData[sfb + sfboffs] =
                sfbThresholdRightLdData[sfb + sfboffs] = minThresholdLdData;
            sfbEnergyLeft[sfb + sfboffs] = sfbEnergyMid[sfb + sfboffs];
            sfbEnergyRight[sfb + sfboffs] = sfbEnergySide[sfb + sfboffs];
            sfbEnergyLeftLdData[sfb + sfboffs] =
                sfbEnergyMidLdData[sfb + sfboffs];
            sfbEnergyRightLdData[sfb + sfboffs] =
                sfbEnergySideLdData[sfb + sfboffs];

            sfbSpreadEnLeft[sfb + sfboffs] = sfbSpreadEnRight[sfb + sfboffs] =
                fixMin(sfbSpreadEnLeft[sfb + sfboffs],
                       sfbSpreadEnRight[sfb + sfboffs]) >>
                1;
          }
        }
      }
    } else {
      *msDigest = SI_MS_MASK_SOME;
    }
  } else {
    *msDigest = SI_MS_MASK_NONE;
  }
}
