/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Chaos measure calculation

*******************************************************************************/

#include "chaosmeasure.h"

static void FDKaacEnc_FDKaacEnc_CalculateChaosMeasurePeakFast(
    FIXP_DBL *RESTRICT paMDCTDataNM0, INT numberOfLines,
    FIXP_DBL *RESTRICT chaosMeasure) {
  INT i, j;

  FIXP_DBL left_0_div2,
      center_0;
  FIXP_DBL left_1_div2, center_1;

  left_0_div2 = (FIXP_DBL)(((LONG)paMDCTDataNM0[0] ^
                            ((LONG)paMDCTDataNM0[0] >> (DFRACT_BITS - 1))) >>
                           1);
  left_1_div2 = (FIXP_DBL)(((LONG)paMDCTDataNM0[1] ^
                            ((LONG)paMDCTDataNM0[1] >> (DFRACT_BITS - 1))) >>
                           1);
  center_0 = (FIXP_DBL)((LONG)paMDCTDataNM0[2] ^
                        ((LONG)paMDCTDataNM0[2] >> (DFRACT_BITS - 1)));
  center_1 = (FIXP_DBL)((LONG)paMDCTDataNM0[3] ^
                        ((LONG)paMDCTDataNM0[3] >> (DFRACT_BITS - 1)));

  for (j = 2; j < numberOfLines - 2; j += 2) {
    FIXP_DBL right_0 =
        (FIXP_DBL)((LONG)paMDCTDataNM0[j + 2] ^
                   ((LONG)paMDCTDataNM0[j + 2] >> (DFRACT_BITS - 1)));
    FIXP_DBL tmp_0 = left_0_div2 + (right_0 >> 1);
    FIXP_DBL right_1 =
        (FIXP_DBL)((LONG)paMDCTDataNM0[j + 3] ^
                   ((LONG)paMDCTDataNM0[j + 3] >> (DFRACT_BITS - 1)));
    FIXP_DBL tmp_1 = left_1_div2 + (right_1 >> 1);

    if (tmp_0 < center_0) {
      INT leadingBits = CntLeadingZeros(center_0) - 1;
      tmp_0 = schur_div(tmp_0 << leadingBits, center_0 << leadingBits, 8);
      tmp_0 = fMult(tmp_0, tmp_0);
    } else {
      tmp_0 = (FIXP_DBL)MAXVAL_DBL;
    }
    chaosMeasure[j + 0] = tmp_0;
    left_0_div2 = center_0 >> 1;
    center_0 = right_0;

    if (tmp_1 < center_1) {
      INT leadingBits = CntLeadingZeros(center_1) - 1;
      tmp_1 = schur_div(tmp_1 << leadingBits, center_1 << leadingBits, 8);
      tmp_1 = fMult(tmp_1, tmp_1);
    } else {
      tmp_1 = (FIXP_DBL)MAXVAL_DBL;
    }

    left_1_div2 = center_1 >> 1;
    center_1 = right_1;
    chaosMeasure[j + 1] = tmp_1;
  }

  chaosMeasure[0] = chaosMeasure[2];
  chaosMeasure[1] = chaosMeasure[2];

  for (i = (numberOfLines - 3); i < numberOfLines; i++)
    chaosMeasure[i] = FL2FXCONST_DBL(0.5);
}

void FDKaacEnc_CalculateChaosMeasure(FIXP_DBL *paMDCTDataNM0, INT numberOfLines,
                                     FIXP_DBL *chaosMeasure)

{
  FDKaacEnc_FDKaacEnc_CalculateChaosMeasurePeakFast(
      paMDCTDataNM0, numberOfLines, chaosMeasure);
}
