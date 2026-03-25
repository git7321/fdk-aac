/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Pre echo control

*******************************************************************************/

#include "pre_echo_control.h"
#include "psy_configuration.h"

void FDKaacEnc_InitPreEchoControl(FIXP_DBL *RESTRICT pbThresholdNm1,
                                  INT *calcPreEcho, INT numPb,
                                  FIXP_DBL *RESTRICT sfbPcmQuantThreshold,
                                  INT *mdctScalenm1) {
  *mdctScalenm1 = PCM_QUANT_THR_SCALE >> 1;

  FDKmemcpy(pbThresholdNm1, sfbPcmQuantThreshold, numPb * sizeof(FIXP_DBL));

  *calcPreEcho = 1;
}

void FDKaacEnc_PreEchoControl(FIXP_DBL *RESTRICT pbThresholdNm1,
                              INT calcPreEcho, INT numPb,
                              INT maxAllowedIncreaseFactor,
                              FIXP_SGL minRemainingThresholdFactor,
                              FIXP_DBL *RESTRICT pbThreshold, INT mdctScale,
                              INT *mdctScalenm1) {
  int i;
  FIXP_DBL tmpThreshold1, tmpThreshold2;
  int scaling;

  if (calcPreEcho == 0) {
    FDKmemcpy(pbThresholdNm1, pbThreshold, numPb * sizeof(FIXP_DBL));
    *mdctScalenm1 = mdctScale;
    return;
  }

  if (mdctScale > *mdctScalenm1) {
    scaling = 2 * (mdctScale - *mdctScalenm1);
    for (i = 0; i < numPb; i++) {
      FDK_ASSERT(scaling >= 0);
      tmpThreshold1 = maxAllowedIncreaseFactor * (pbThresholdNm1[i] >> scaling);
      tmpThreshold2 = fMult(minRemainingThresholdFactor, pbThreshold[i]);

      FIXP_DBL tmp = pbThreshold[i];

      pbThresholdNm1[i] = tmp;

      tmp = fixMin(tmp, tmpThreshold1);
      pbThreshold[i] = fixMax(tmp, tmpThreshold2);
    }
  } else {
    scaling = 2 * (*mdctScalenm1 - mdctScale);
    for (i = 0; i < numPb; i++) {
      tmpThreshold1 = (maxAllowedIncreaseFactor >> 1) * pbThresholdNm1[i];
      tmpThreshold2 = fMult(minRemainingThresholdFactor, pbThreshold[i]);

      pbThresholdNm1[i] = pbThreshold[i];

      FDK_ASSERT(scaling >= 0);
      if ((pbThreshold[i] >> (scaling + 1)) > tmpThreshold1) {
        pbThreshold[i] = tmpThreshold1 << (scaling + 1);
      }
      pbThreshold[i] = fixMax(pbThreshold[i], tmpThreshold2);
    }
  }

  *mdctScalenm1 = mdctScale;
}
