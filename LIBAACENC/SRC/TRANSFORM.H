/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: MDCT Transform

*******************************************************************************/

#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "mdct.h"
#include "common_fix.h"

#define WTS0 1
#define WTS1 0
#define WTS2 -2

INT FDKaacEnc_Transform_Real(const INT_PCM* pTimeData,
                             FIXP_DBL* RESTRICT mdctData, const INT blockType,
                             const INT windowShape, INT* prevWindowShape,
                             H_MDCT mdctPers, const INT frameLength,
                             INT* pMdctData_e, INT filterType);

INT FDKaacEnc_Transform_Real_Eld(const INT_PCM* pTimeData,
                                 FIXP_DBL* RESTRICT mdctData,
                                 const INT blockType, const INT windowShape,
                                 INT* prevWindowShape, const INT frameLength,
                                 INT* mdctData_e, INT filterType,
                                 FIXP_DBL* RESTRICT overlapAddBuffer);

#endif /* #!defined (TRANSFORM_H) */
