/**************************** AAC encoder library ******************************

   Author(s):   Alex Goeschel

   Description: Temporal noise shaping

*******************************************************************************/

#ifndef TNS_FUNC_H
#define TNS_FUNC_H

#include "common_fix.h"
#include "psy_configuration.h"

AAC_ENCODER_ERROR FDKaacEnc_InitTnsConfiguration(
    INT bitrate, INT samplerate, INT channels, INT blocktype, INT granuleLength,
    INT isLowDelay, INT ldSbrPresent, TNS_CONFIG *tnsConfig,
    PSY_CONFIGURATION *psyConfig, INT active, INT useTnsPeak);

INT FDKaacEnc_TnsDetect(TNS_DATA *tnsData, const TNS_CONFIG *tC,
                        TNS_INFO *tnsInfo, INT sfbCnt, const FIXP_DBL *spectrum,
                        INT subBlockNumber, INT blockType);

void FDKaacEnc_TnsSync(TNS_DATA *tnsDataDest, const TNS_DATA *tnsDataSrc,
                       TNS_INFO *tnsInfoDest, TNS_INFO *tnsInfoSrc,
                       const INT blockTypeDest, const INT blockTypeSrc,
                       const TNS_CONFIG *tC);

INT FDKaacEnc_TnsEncode(TNS_INFO *tnsInfo, TNS_DATA *tnsData,
                        const INT numOfSfb, const TNS_CONFIG *tC,
                        const INT lowPassLine, FIXP_DBL *spectrum,
                        const INT subBlockNumber, const INT blockType);

#endif /* TNS_FUNC_H */
