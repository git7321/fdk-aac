/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Perceptual entropie module

*******************************************************************************/

#include "line_pe.h"
#include "sf_estim.h"
#include "bit_cnt.h"

#include "genericStds.h"

static const FIXP_DBL C1LdData =
    FL2FXCONST_DBL(3.0 / LD_DATA_SCALING);
static const FIXP_DBL C2LdData = FL2FXCONST_DBL(
    1.3219281 / LD_DATA_SCALING);
static const FIXP_DBL C3LdData = FL2FXCONST_DBL(0.5593573);

void FDKaacEnc_prepareSfbPe(PE_CHANNEL_DATA *RESTRICT const peChanData,
                            const FIXP_DBL *RESTRICT const sfbEnergyLdData,
                            const FIXP_DBL *RESTRICT const sfbThresholdLdData,
                            const FIXP_DBL *RESTRICT const sfbFormFactorLdData,
                            const INT *RESTRICT const sfbOffset,
                            const INT sfbCnt, const INT sfbPerGroup,
                            const INT maxSfbPerGroup) {
  INT sfbGrp, sfb;
  INT sfbWidth;
  FIXP_DBL avgFormFactorLdData;
  const FIXP_DBL formFacScaling =
      FL2FXCONST_DBL((float)FORM_FAC_SHIFT / LD_DATA_SCALING);

  for (sfbGrp = 0; sfbGrp < sfbCnt; sfbGrp += sfbPerGroup) {
    for (sfb = 0; sfb < maxSfbPerGroup; sfb++) {
      if ((FIXP_DBL)sfbEnergyLdData[sfbGrp + sfb] >
          (FIXP_DBL)sfbThresholdLdData[sfbGrp + sfb]) {
        sfbWidth = sfbOffset[sfbGrp + sfb + 1] - sfbOffset[sfbGrp + sfb];
        avgFormFactorLdData = ((-sfbEnergyLdData[sfbGrp + sfb] >> 1) +
                               (CalcLdInt(sfbWidth) >> 1)) >>
                              1;
        peChanData->sfbNLines[sfbGrp + sfb] = (INT)CalcInvLdData(
            (sfbFormFactorLdData[sfbGrp + sfb] + formFacScaling) +
            avgFormFactorLdData);
        peChanData->sfbNLines[sfbGrp + sfb] =
            fMin(sfbWidth, peChanData->sfbNLines[sfbGrp + sfb]);
      } else {
        peChanData->sfbNLines[sfbGrp + sfb] = 0;
      }
    }
  }
}

void FDKaacEnc_calcSfbPe(PE_CHANNEL_DATA *RESTRICT const peChanData,
                         const FIXP_DBL *RESTRICT const sfbEnergyLdData,
                         const FIXP_DBL *RESTRICT const sfbThresholdLdData,
                         const INT sfbCnt, const INT sfbPerGroup,
                         const INT maxSfbPerGroup,
                         const INT *RESTRICT const isBook,
                         const INT *RESTRICT const isScale) {
  INT sfbGrp, sfb, thisSfb;
  INT nLines;
  FIXP_DBL logDataRatio;
  FIXP_DBL scaleLd = (FIXP_DBL)0;
  INT lastValIs = 0;

  FIXP_DBL pe = 0;
  FIXP_DBL constPart = 0;
  FIXP_DBL nActiveLines = 0;

  FIXP_DBL tmpPe, tmpConstPart, tmpNActiveLines;

  for (sfbGrp = 0; sfbGrp < sfbCnt; sfbGrp += sfbPerGroup) {
    for (sfb = 0; sfb < maxSfbPerGroup; sfb++) {
      tmpPe = (FIXP_DBL)0;
      tmpConstPart = (FIXP_DBL)0;
      tmpNActiveLines = (FIXP_DBL)0;

      thisSfb = sfbGrp + sfb;

      if (sfbEnergyLdData[thisSfb] > sfbThresholdLdData[thisSfb]) {
        logDataRatio = sfbEnergyLdData[thisSfb] - sfbThresholdLdData[thisSfb];
        nLines = peChanData->sfbNLines[thisSfb];

        FIXP_DBL factor = nLines << (LD_DATA_SHIFT + PE_CONSTPART_SHIFT + 1);
        if (logDataRatio >= C1LdData) {
          tmpPe = fMultDiv2(logDataRatio, factor);
          tmpConstPart = fMultDiv2(sfbEnergyLdData[thisSfb] + scaleLd, factor);
        } else {
          tmpPe = fMultDiv2(
              ((FIXP_DBL)C2LdData + fMult(C3LdData, logDataRatio)), factor);
          tmpConstPart =
              fMultDiv2(((FIXP_DBL)C2LdData +
                         fMult(C3LdData, sfbEnergyLdData[thisSfb] + scaleLd)),
                        factor);

          nLines = fMultI(C3LdData, nLines);
        }
        tmpNActiveLines = (FIXP_DBL)nLines;
      } else if (isBook[thisSfb]) {
        INT delta = isScale[thisSfb] - lastValIs;
        lastValIs = isScale[thisSfb];
        peChanData->sfbPe[thisSfb] = FDKaacEnc_bitCountScalefactorDelta(delta)
                                     << PE_CONSTPART_SHIFT;
        peChanData->sfbConstPart[thisSfb] = 0;
        peChanData->sfbNActiveLines[thisSfb] = 0;
      }
      peChanData->sfbPe[thisSfb] = tmpPe;
      peChanData->sfbConstPart[thisSfb] = tmpConstPart;
      peChanData->sfbNActiveLines[thisSfb] = tmpNActiveLines;

      pe += tmpPe;
      constPart += tmpConstPart;
      nActiveLines += tmpNActiveLines;
    }
  }

  peChanData->pe = pe >> PE_CONSTPART_SHIFT;
  peChanData->constPart = constPart >> PE_CONSTPART_SHIFT;

  peChanData->nActiveLines = nActiveLines;
}
