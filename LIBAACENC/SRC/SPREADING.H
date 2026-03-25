/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Spreading of energy and weighted tonality

*******************************************************************************/

#ifndef SPREADING_H
#define SPREADING_H

#include "common_fix.h"

void FDKaacEnc_SpreadingMax(const INT pbCnt,
                            const FIXP_DBL *RESTRICT maskLowFactor,
                            const FIXP_DBL *RESTRICT maskHighFactor,
                            FIXP_DBL *RESTRICT pbSpreadEnergy);

#endif /* #ifndef SPREADING_H */
