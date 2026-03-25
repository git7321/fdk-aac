/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Short block grouping

*******************************************************************************/

#include "psy_const.h"
#include "interface.h"

static inline FIXP_DBL nrgAddSaturate(const FIXP_DBL a, const FIXP_DBL b) {
  return ((a >= (FIXP_DBL)MAXVAL_DBL - b) ? (FIXP_DBL)MAXVAL_DBL : (a + b));
}

void FDKaacEnc_groupShortData(FIXP_DBL *mdctSpectrum,
                              SFB_THRESHOLD *sfbThreshold,
                              SFB_ENERGY *sfbEnergy,
                              SFB_ENERGY *sfbEnergyMS,
                              SFB_ENERGY *sfbSpreadEnergy, const INT sfbCnt,
                              const INT sfbActive, const INT *sfbOffset,
                              const FIXP_DBL *sfbMinSnrLdData,
                              INT *groupedSfbOffset,
                              INT *maxSfbPerGroup,
                              FIXP_DBL *groupedSfbMinSnrLdData,
                              const INT noOfGroups, const INT *groupLen,
                              const INT granuleLength) {
  INT i, j;
  INT line;
  INT sfb;
  INT grp;
  INT wnd;
  INT offset;
  INT highestSfb;
  INT granuleLength_short = granuleLength / TRANS_FAC;

  C_ALLOC_SCRATCH_START(tmpSpectrum, FIXP_DBL, (1024))

  highestSfb = 0;
  for (wnd = 0; wnd < TRANS_FAC; wnd++) {
    for (sfb = sfbActive - 1; sfb >= highestSfb; sfb--) {
      for (line = sfbOffset[sfb + 1] - 1; line >= sfbOffset[sfb]; line--) {
        if (mdctSpectrum[wnd * granuleLength_short + line] !=
            FL2FXCONST_SPC(0.0))
          break;
      }
      if (line >= sfbOffset[sfb]) break;
    }
    highestSfb = fixMax(highestSfb, sfb);
  }
  highestSfb = highestSfb > 0 ? highestSfb : 0;
  *maxSfbPerGroup = highestSfb + 1;

  i = 0;
  offset = 0;
  for (grp = 0; grp < noOfGroups; grp++) {
    for (sfb = 0; sfb < sfbActive + 1; sfb++) {
      groupedSfbOffset[i++] = offset + sfbOffset[sfb] * groupLen[grp];
    }
    i += sfbCnt - sfb;
    offset += groupLen[grp] * granuleLength_short;
  }
  groupedSfbOffset[i++] = granuleLength;

  i = 0;
  for (grp = 0; grp < noOfGroups; grp++) {
    for (sfb = 0; sfb < sfbActive; sfb++) {
      groupedSfbMinSnrLdData[i++] = sfbMinSnrLdData[sfb];
    }
    i += sfbCnt - sfb;
  }

  wnd = 0;
  i = 0;
  for (grp = 0; grp < noOfGroups; grp++) {
    for (sfb = 0; sfb < sfbActive; sfb++) {
      FIXP_DBL thresh = sfbThreshold->Short[wnd][sfb];
      for (j = 1; j < groupLen[grp]; j++) {
        thresh = nrgAddSaturate(thresh, sfbThreshold->Short[wnd + j][sfb]);
      }
      sfbThreshold->Long[i++] = thresh;
    }
    i += sfbCnt - sfb;
    wnd += groupLen[grp];
  }

  wnd = 0;
  i = 0;
  for (grp = 0; grp < noOfGroups; grp++) {
    for (sfb = 0; sfb < sfbActive; sfb++) {
      FIXP_DBL energy = sfbEnergy->Short[wnd][sfb];
      for (j = 1; j < groupLen[grp]; j++) {
        energy = nrgAddSaturate(energy, sfbEnergy->Short[wnd + j][sfb]);
      }
      sfbEnergy->Long[i++] = energy;
    }
    i += sfbCnt - sfb;
    wnd += groupLen[grp];
  }

  wnd = 0;
  i = 0;
  for (grp = 0; grp < noOfGroups; grp++) {
    for (sfb = 0; sfb < sfbActive; sfb++) {
      FIXP_DBL energy = sfbEnergyMS->Short[wnd][sfb];
      for (j = 1; j < groupLen[grp]; j++) {
        energy = nrgAddSaturate(energy, sfbEnergyMS->Short[wnd + j][sfb]);
      }
      sfbEnergyMS->Long[i++] = energy;
    }
    i += sfbCnt - sfb;
    wnd += groupLen[grp];
  }

  wnd = 0;
  i = 0;
  for (grp = 0; grp < noOfGroups; grp++) {
    for (sfb = 0; sfb < sfbActive; sfb++) {
      FIXP_DBL energy = sfbSpreadEnergy->Short[wnd][sfb];
      for (j = 1; j < groupLen[grp]; j++) {
        energy = nrgAddSaturate(energy, sfbSpreadEnergy->Short[wnd + j][sfb]);
      }
      sfbSpreadEnergy->Long[i++] = energy;
    }
    i += sfbCnt - sfb;
    wnd += groupLen[grp];
  }

  wnd = 0;
  i = 0;
  for (grp = 0; grp < noOfGroups; grp++) {
    for (sfb = 0; sfb < sfbActive; sfb++) {
      int width = sfbOffset[sfb + 1] - sfbOffset[sfb];
      FIXP_DBL *pMdctSpectrum =
          &mdctSpectrum[sfbOffset[sfb]] + wnd * granuleLength_short;
      for (j = 0; j < groupLen[grp]; j++) {
        FIXP_DBL *pTmp = pMdctSpectrum;
        for (line = width; line > 0; line--) {
          tmpSpectrum[i++] = *pTmp++;
        }
        pMdctSpectrum += granuleLength_short;
      }
    }
    i += (groupLen[grp] * (sfbOffset[sfbCnt] - sfbOffset[sfb]));
    wnd += groupLen[grp];
  }

  FDKmemcpy(mdctSpectrum, tmpSpectrum, granuleLength * sizeof(FIXP_DBL));
  C_ALLOC_SCRATCH_END(tmpSpectrum, FIXP_DBL, (1024))
}
