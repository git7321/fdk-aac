/**************************** AAC encoder library ******************************

   Author(s):   M. Lohwasser

   Description: noisedet.c
                Routines for Noise Detection

*******************************************************************************/

#include "noisedet.h"
#include "aacenc_pns.h"
#include "pnsparam.h"

static FIXP_SGL FDKaacEnc_fuzzyIsSmaller(FIXP_DBL testVal, FIXP_DBL refVal,
                                         FIXP_DBL loLim, FIXP_DBL hiLim) {
  if (refVal <= FL2FXCONST_DBL(0.0))
    return (FL2FXCONST_SGL(0.0f));
  else if (testVal >= fMult((hiLim >> 1) + (loLim >> 1), refVal))
    return (FL2FXCONST_SGL(0.0f));
  else
    return ((FIXP_SGL)MAXVAL_SGL);
}

void FDKaacEnc_noiseDetect(FIXP_DBL *RESTRICT mdctSpectrum,
                           INT *RESTRICT sfbMaxScaleSpec, INT sfbActive,
                           const INT *RESTRICT sfbOffset,
                           FIXP_SGL *RESTRICT noiseFuzzyMeasure,
                           NOISEPARAMS *np, FIXP_SGL *RESTRICT sfbtonality)

{
  int i, k, sfb, sfbWidth;
  FIXP_SGL fuzzy, fuzzyTotal;
  FIXP_DBL refVal, testVal;

  for (sfb = 0; sfb < sfbActive; sfb++) {
    fuzzyTotal = (FIXP_SGL)MAXVAL_SGL;
    sfbWidth = sfbOffset[sfb + 1] - sfbOffset[sfb];

    if (sfb < np->startSfb || sfbWidth < np->minSfbWidth) {
      noiseFuzzyMeasure[sfb] = FL2FXCONST_SGL(0.0f);
      continue;
    }

    if ((np->detectionAlgorithmFlags & USE_POWER_DISTRIBUTION) &&
        (fuzzyTotal > FL2FXCONST_SGL(0.5f))) {
      FIXP_DBL fhelp1, fhelp2, fhelp3, fhelp4, maxVal, minVal;
      INT leadingBits = fixMax(
          0, (sfbMaxScaleSpec[sfb] -
              3));

      fhelp1 = fhelp2 = fhelp3 = fhelp4 = FL2FXCONST_DBL(0.0f);
      k = sfbWidth >> 2;

      for (i = sfbOffset[sfb]; i < sfbOffset[sfb] + k; i++) {
        fhelp1 = fPow2AddDiv2(fhelp1, mdctSpectrum[i] << leadingBits);
        fhelp2 = fPow2AddDiv2(fhelp2, mdctSpectrum[i + k] << leadingBits);
        fhelp3 = fPow2AddDiv2(fhelp3, mdctSpectrum[i + 2 * k] << leadingBits);
        fhelp4 = fPow2AddDiv2(fhelp4, mdctSpectrum[i + 3 * k] << leadingBits);
      }

      maxVal = fixMax(fhelp1, fhelp2);
      maxVal = fixMax(maxVal, fhelp3);
      maxVal = fixMax(maxVal, fhelp4);

      minVal = fixMin(fhelp1, fhelp2);
      minVal = fixMin(minVal, fhelp3);
      minVal = fixMin(minVal, fhelp4);

      leadingBits = CountLeadingBits(maxVal);
      testVal = maxVal << leadingBits;
      refVal = minVal << leadingBits;

      testVal = fMultDiv2(testVal, np->powDistPSDcurve[sfb]);

      fuzzy = FDKaacEnc_fuzzyIsSmaller(
          testVal,
          refVal,
          FL2FXCONST_DBL(0.495),
          FL2FXCONST_DBL(0.505));

      fuzzyTotal = fixMin(fuzzyTotal, fuzzy);
    }

    if ((np->detectionAlgorithmFlags & USE_PSYCH_TONALITY) &&
        (fuzzyTotal > FL2FXCONST_SGL(0.5f))) {

      testVal = FX_SGL2FX_DBL(sfbtonality[sfb]) >> 1;
      refVal = np->refTonality;

      fuzzy = FDKaacEnc_fuzzyIsSmaller(
          testVal, refVal, FL2FXCONST_DBL(0.45f),
          FL2FXCONST_DBL(0.55f));

      fuzzyTotal = fixMin(fuzzyTotal, fuzzy);
    }

    noiseFuzzyMeasure[sfb] = fuzzyTotal;
  }
}
