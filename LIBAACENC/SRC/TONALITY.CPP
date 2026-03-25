/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Convert chaos measure to the tonality index

*******************************************************************************/

#include "tonality.h"

#include "chaosmeasure.h"

#if defined(__arm__)
#endif

static const FIXP_DBL normlog =
    (FIXP_DBL)0xd977d949;

static void FDKaacEnc_CalcSfbTonality(FIXP_DBL *RESTRICT spectrum,
                                      INT *RESTRICT sfbMaxScaleSpec,
                                      FIXP_DBL *RESTRICT chaosMeasure,
                                      FIXP_SGL *RESTRICT sfbTonality,
                                      INT sfbCnt, const INT *RESTRICT sfbOffset,
                                      FIXP_DBL *RESTRICT sfbEnergyLD64);

void FDKaacEnc_CalculateFullTonality(FIXP_DBL *RESTRICT spectrum,
                                     INT *RESTRICT sfbMaxScaleSpec,
                                     FIXP_DBL *RESTRICT sfbEnergyLD64,
                                     FIXP_SGL *RESTRICT sfbTonality, INT sfbCnt,
                                     const INT *sfbOffset, INT usePns) {
  INT j;
  INT numberOfLines = sfbOffset[sfbCnt];

  if (usePns) {
    C_ALLOC_SCRATCH_START(chaosMeasurePerLine, FIXP_DBL, (1024))

    FDKaacEnc_CalculateChaosMeasure(spectrum, numberOfLines,
                                    chaosMeasurePerLine);

    FIXP_DBL left = chaosMeasurePerLine[0];
    FIXP_DBL right;
    for (j = 1; j < (numberOfLines - 1); j += 2) {
      right = chaosMeasurePerLine[j];
      right = right - (right >> 2);
      left = right + (left >> 2);
      chaosMeasurePerLine[j] = left;

      right = chaosMeasurePerLine[j + 1];
      right = right - (right >> 2);
      left = right + (left >> 2);
      chaosMeasurePerLine[j + 1] = left;
    }
    if (j == (numberOfLines - 1)) {
      right = chaosMeasurePerLine[j];
      right = right - (right >> 2);
      left = right + (left >> 2);
      chaosMeasurePerLine[j] = left;
    }

    FDKaacEnc_CalcSfbTonality(spectrum, sfbMaxScaleSpec, chaosMeasurePerLine,
                              sfbTonality, sfbCnt, sfbOffset, sfbEnergyLD64);

    C_ALLOC_SCRATCH_END(chaosMeasurePerLine, FIXP_DBL, (1024))
  }
}

static void FDKaacEnc_CalcSfbTonality(FIXP_DBL *RESTRICT spectrum,
                                      INT *RESTRICT sfbMaxScaleSpec,
                                      FIXP_DBL *RESTRICT chaosMeasure,
                                      FIXP_SGL *RESTRICT sfbTonality,
                                      INT sfbCnt, const INT *RESTRICT sfbOffset,
                                      FIXP_DBL *RESTRICT sfbEnergyLD64) {
  INT i;

  for (i = 0; i < sfbCnt; i++) {
    FIXP_DBL chaosMeasureSfbLD64;
    INT shiftBits =
        fixMax(0, sfbMaxScaleSpec[i] -
                      4);

    INT j;
    FIXP_DBL chaosMeasureSfb = FL2FXCONST_DBL(0.0);

    for (j = (sfbOffset[i + 1] - sfbOffset[i]) - 1; j >= 0; j--) {
      FIXP_DBL tmp = (*spectrum++) << shiftBits;
      FIXP_DBL lineNrg = fMultDiv2(tmp, tmp);
      chaosMeasureSfb = fMultAddDiv2(chaosMeasureSfb, lineNrg, *chaosMeasure++);
    }

    if (chaosMeasureSfb != FL2FXCONST_DBL(0.0)) {
      chaosMeasureSfbLD64 = CalcLdData((chaosMeasureSfb)) - sfbEnergyLD64[i];
      chaosMeasureSfbLD64 += FL2FXCONST_DBL(3.0f / 64) -
                             ((FIXP_DBL)(shiftBits) << (DFRACT_BITS - 6));

      if (chaosMeasureSfbLD64 >
          FL2FXCONST_DBL(-0.0519051))
      {
        if (chaosMeasureSfbLD64 <= FL2FXCONST_DBL(0.0))
          sfbTonality[i] =
              FX_DBL2FX_SGL(fMultDiv2(chaosMeasureSfbLD64, normlog) << 7);
        else
          sfbTonality[i] = FL2FXCONST_SGL(0.0);
      } else
        sfbTonality[i] = (FIXP_SGL)MAXVAL_SGL;
    } else
      sfbTonality[i] = (FIXP_SGL)MAXVAL_SGL;
  }
}
