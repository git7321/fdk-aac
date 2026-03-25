/**************************** AAC encoder library ******************************

   Author(s):   M. Lohwasser

   Description: pns_func.h

*******************************************************************************/

#ifndef PNS_FUNC_H
#define PNS_FUNC_H

#include "common_fix.h"
#include "aacenc_pns.h"
#include "psy_data.h"

AAC_ENCODER_ERROR FDKaacEnc_InitPnsConfiguration(
    PNS_CONFIG *pnsConf, INT bitRate, INT sampleRate, INT usePns, INT sfbCnt,
    const INT *sfbOffset, const INT numChan, const INT isLC);

void FDKaacEnc_PnsDetect(PNS_CONFIG *pnsConf, PNS_DATA *pnsData,
                         const INT lastWindowSequence, const INT sfbActive,
                         const INT maxSfbPerGroup, FIXP_DBL *sfbThresholdLdData,
                         const INT *sfbOffset, FIXP_DBL *mdctSpectrum,
                         INT *sfbMaxScaleSpec, FIXP_SGL *sfbtonality,
                         int tnsOrder, INT tnsPredictionGain, INT tnsActive,
                         FIXP_DBL *sfbEnergyLdData, INT *noiseNrg);

void FDKaacEnc_CodePnsChannel(const INT sfbActive, PNS_CONFIG *pnsConf,
                              INT *pnsFlag, FIXP_DBL *sfbEnergy, INT *noiseNrg,
                              FIXP_DBL *sfbThreshold);

void FDKaacEnc_PreProcessPnsChannelPair(
    const INT sfbActive, FIXP_DBL *sfbEnergyLeft, FIXP_DBL *sfbEnergyRight,
    FIXP_DBL *sfbEnergyLeftLD, FIXP_DBL *sfbEnergyRightLD,
    FIXP_DBL *sfbEnergyMid, PNS_CONFIG *pnsConfLeft, PNS_DATA *pnsDataLeft,
    PNS_DATA *pnsDataRight);

void FDKaacEnc_PostProcessPnsChannelPair(const INT sfbActive,
                                         PNS_CONFIG *pnsConf,
                                         PNS_DATA *pnsDataLeft,
                                         PNS_DATA *pnsDataRight, INT *msMask,
                                         INT *msDigest);

#endif /* PNS_FUNC_H */
