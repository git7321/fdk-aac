/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Quantization

*******************************************************************************/

#include "quantize.h"

#include "aacEnc_rom.h"

static void FDKaacEnc_quantizeLines(INT gain, INT noOfLines,
                                    const FIXP_DBL *mdctSpectrum,
                                    SHORT *quaSpectrum, INT dZoneQuantEnable) {
  int line;
  FIXP_DBL k = FL2FXCONST_DBL(0.0f);
  FIXP_QTD quantizer = FDKaacEnc_quantTableQ[(-gain) & 3];
  INT quantizershift = ((-gain) >> 2) + 1;
  const INT kShift = 16;

  if (dZoneQuantEnable)
    k = FL2FXCONST_DBL(0.23f) >> kShift;
  else
    k = FL2FXCONST_DBL(-0.0946f + 0.5f) >> kShift;

  for (line = 0; line < noOfLines; line++) {
    FIXP_DBL accu = fMultDiv2(mdctSpectrum[line], quantizer);

    if (accu < FL2FXCONST_DBL(0.0f)) {
      accu = -accu;
      INT accuShift = CntLeadingZeros(accu) - 1;
      accu <<= accuShift;
      INT tabIndex =
          (INT)(accu >> (DFRACT_BITS - 2 - MANT_DIGITS)) & (~MANT_SIZE);
      INT totalShift = quantizershift - accuShift + 1;
      accu = fMultDiv2(FDKaacEnc_mTab_3_4[tabIndex],
                       FDKaacEnc_quantTableE[totalShift & 3]);
      totalShift = (16 - 4) - (3 * (totalShift >> 2));
      FDK_ASSERT(totalShift >= 0);
      accu >>= fixMin(totalShift, DFRACT_BITS - 1);
      quaSpectrum[line] =
          (SHORT)(-((LONG)(k + accu) >> (DFRACT_BITS - 1 - 16)));
    } else if (accu > FL2FXCONST_DBL(0.0f)) {
      INT accuShift = CntLeadingZeros(accu) - 1;
      accu <<= accuShift;
      INT tabIndex =
          (INT)(accu >> (DFRACT_BITS - 2 - MANT_DIGITS)) & (~MANT_SIZE);
      INT totalShift = quantizershift - accuShift + 1;
      accu = fMultDiv2(FDKaacEnc_mTab_3_4[tabIndex],
                       FDKaacEnc_quantTableE[totalShift & 3]);
      totalShift = (16 - 4) - (3 * (totalShift >> 2));
      FDK_ASSERT(totalShift >= 0);
      accu >>= fixMin(totalShift, DFRACT_BITS - 1);
      quaSpectrum[line] = (SHORT)((LONG)(k + accu) >> (DFRACT_BITS - 1 - 16));
    } else {
      quaSpectrum[line] = 0;
    }
  }
}

static void FDKaacEnc_invQuantizeLines(INT gain, INT noOfLines,
                                       SHORT *quantSpectrum,
                                       FIXP_DBL *mdctSpectrum)

{
  INT iquantizermod;
  INT iquantizershift;
  INT line;

  iquantizermod = gain & 3;
  iquantizershift = gain >> 2;

  for (line = 0; line < noOfLines; line++) {
    if (quantSpectrum[line] < 0) {
      FIXP_DBL accu;
      INT ex, specExp, tabIndex;
      FIXP_DBL s, t;

      accu = (FIXP_DBL)-quantSpectrum[line];

      ex = CountLeadingBits(accu);
      accu <<= ex;
      specExp = (DFRACT_BITS - 1) - ex;

      FDK_ASSERT(specExp < 14);

      tabIndex = (INT)(accu >> (DFRACT_BITS - 2 - MANT_DIGITS)) & (~MANT_SIZE);

      s = FDKaacEnc_mTab_4_3Elc[tabIndex];

      t = FDKaacEnc_specExpMantTableCombElc[iquantizermod][specExp];

      accu = fMult(s, t);

      specExp = FDKaacEnc_specExpTableComb[iquantizermod][specExp] -
                1;

      if ((-iquantizershift - specExp) < 0)
        accu <<= -(-iquantizershift - specExp);
      else
        accu >>= -iquantizershift - specExp;

      mdctSpectrum[line] = -accu;
    } else if (quantSpectrum[line] > 0) {
      FIXP_DBL accu;
      INT ex, specExp, tabIndex;
      FIXP_DBL s, t;

      accu = (FIXP_DBL)(INT)quantSpectrum[line];

      ex = CountLeadingBits(accu);
      accu <<= ex;
      specExp = (DFRACT_BITS - 1) - ex;

      FDK_ASSERT(specExp < 14);

      tabIndex = (INT)(accu >> (DFRACT_BITS - 2 - MANT_DIGITS)) & (~MANT_SIZE);

      s = FDKaacEnc_mTab_4_3Elc[tabIndex];

      t = FDKaacEnc_specExpMantTableCombElc[iquantizermod][specExp];

      accu = fMult(s, t);

      specExp = FDKaacEnc_specExpTableComb[iquantizermod][specExp] -
                1;

      if ((-iquantizershift - specExp) < 0)
        accu <<= -(-iquantizershift - specExp);
      else
        accu >>= -iquantizershift - specExp;

      mdctSpectrum[line] = accu;
    } else {
      mdctSpectrum[line] = FL2FXCONST_DBL(0.0f);
    }
  }
}

void FDKaacEnc_QuantizeSpectrum(INT sfbCnt, INT maxSfbPerGroup, INT sfbPerGroup,
                                const INT *sfbOffset,
                                const FIXP_DBL *mdctSpectrum, INT globalGain,
                                const INT *scalefactors,
                                SHORT *quantizedSpectrum,
                                INT dZoneQuantEnable) {
  INT sfbOffs, sfb;

  for (sfbOffs = 0; sfbOffs < sfbCnt; sfbOffs += sfbPerGroup)
    for (sfb = 0; sfb < maxSfbPerGroup; sfb++) {
      INT scalefactor = scalefactors[sfbOffs + sfb];

      FDKaacEnc_quantizeLines(
          globalGain - scalefactor,
          sfbOffset[sfbOffs + sfb + 1] - sfbOffset[sfbOffs + sfb],
          mdctSpectrum + sfbOffset[sfbOffs + sfb],
          quantizedSpectrum + sfbOffset[sfbOffs + sfb], dZoneQuantEnable);
    }
}

FIXP_DBL FDKaacEnc_calcSfbDist(const FIXP_DBL *mdctSpectrum,
                               SHORT *quantSpectrum, INT noOfLines, INT gain,
                               INT dZoneQuantEnable) {
  INT i, scale;
  FIXP_DBL xfsf;
  FIXP_DBL diff;
  FIXP_DBL invQuantSpec;

  xfsf = FL2FXCONST_DBL(0.0f);

  for (i = 0; i < noOfLines; i++) {
    FDKaacEnc_quantizeLines(gain, 1, &mdctSpectrum[i], &quantSpectrum[i],
                            dZoneQuantEnable);

    if (fAbs(quantSpectrum[i]) > MAX_QUANT) {
      return FL2FXCONST_DBL(0.0f);
    }
    FDKaacEnc_invQuantizeLines(gain, 1, &quantSpectrum[i], &invQuantSpec);

    diff = fixp_abs(fixp_abs(invQuantSpec) - fixp_abs(mdctSpectrum[i] >> 1));

    scale = CountLeadingBits(diff);
    diff = scaleValue(diff, scale);
    diff = fPow2(diff);
    scale = fixMin(2 * (scale - 1), DFRACT_BITS - 1);

    diff = scaleValue(diff, -scale);

    xfsf = xfsf + diff;
  }

  xfsf = CalcLdData(xfsf);

  return xfsf;
}

void FDKaacEnc_calcSfbQuantEnergyAndDist(FIXP_DBL *mdctSpectrum,
                                         SHORT *quantSpectrum, INT noOfLines,
                                         INT gain, FIXP_DBL *en,
                                         FIXP_DBL *dist) {
  INT i, scale;
  FIXP_DBL invQuantSpec;
  FIXP_DBL diff;

  FIXP_DBL energy = FL2FXCONST_DBL(0.0f);
  FIXP_DBL distortion = FL2FXCONST_DBL(0.0f);

  for (i = 0; i < noOfLines; i++) {
    if (fAbs(quantSpectrum[i]) > MAX_QUANT) {
      *en = FL2FXCONST_DBL(0.0f);
      *dist = FL2FXCONST_DBL(0.0f);
      return;
    }

    FDKaacEnc_invQuantizeLines(gain, 1, &quantSpectrum[i], &invQuantSpec);

    energy += fPow2(invQuantSpec);

    diff = fixp_abs(fixp_abs(invQuantSpec) - fixp_abs(mdctSpectrum[i] >> 1));

    scale = CountLeadingBits(diff);
    diff = scaleValue(diff, scale);
    diff = fPow2(diff);

    scale = fixMin(2 * (scale - 1), DFRACT_BITS - 1);

    diff = scaleValue(diff, -scale);

    distortion += diff;
  }

  *en = CalcLdData(energy) + FL2FXCONST_DBL(0.03125f);
  *dist = CalcLdData(distortion);
}
