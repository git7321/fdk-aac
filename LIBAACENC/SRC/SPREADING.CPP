/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Spreading of energy

*******************************************************************************/

#include "spreading.h"

void FDKaacEnc_SpreadingMax(const INT pbCnt,
                            const FIXP_DBL *RESTRICT maskLowFactor,
                            const FIXP_DBL *RESTRICT maskHighFactor,
                            FIXP_DBL *RESTRICT pbSpreadEnergy) {
  int i;
  FIXP_DBL delay;

  delay = pbSpreadEnergy[0];
  for (i = 1; i < pbCnt; i++) {
    delay = fixMax(pbSpreadEnergy[i], fMult(maskHighFactor[i], delay));
    pbSpreadEnergy[i] = delay;
  }

  delay = pbSpreadEnergy[pbCnt - 1];
  for (i = pbCnt - 2; i >= 0; i--) {
    delay = fixMax(pbSpreadEnergy[i], fMult(maskLowFactor[i], delay));
    pbSpreadEnergy[i] = delay;
  }
}
