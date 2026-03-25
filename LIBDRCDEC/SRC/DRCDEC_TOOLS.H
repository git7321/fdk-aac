/************************* MPEG-D DRC decoder library **************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef DRCDEC_TOOLS_H
#define DRCDEC_TOOLS_H

#include "drcDec_types.h"
#include "drcDec_selectionProcess.h"

int getDeltaTmin(const int sampleRate);

DRC_COEFFICIENTS_UNI_DRC* selectDrcCoefficients(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, const int location);

DRC_INSTRUCTIONS_UNI_DRC* selectDrcInstructions(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, const int drcSetId);

DOWNMIX_INSTRUCTIONS* selectDownmixInstructions(
    HANDLE_UNI_DRC_CONFIG hUniDrcConfig, const int downmixId);

DRC_ERROR
deriveDrcChannelGroups(
    const int drcSetEffect,
    const int channelCount,
    const SCHAR* gainSetIndex,
    const DUCKING_MODIFICATION* duckingModificationForChannel,
    UCHAR* nDrcChannelGroups,
    SCHAR* uniqueIndex,
    SCHAR* groupForChannel,
    DUCKING_MODIFICATION* duckingModificationForChannelGroup);

FIXP_DBL
dB2lin(const FIXP_DBL dB_m, const int dB_e, int* pLin_e);

FIXP_DBL
lin2dB(const FIXP_DBL lin_m, const int lin_e, int* pDb_e);

FIXP_DBL
approxDb2lin(const FIXP_DBL dB_m, const int dB_e, int* pLin_e);

int bitstreamContainsMultibandDrc(HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
                                  const int downmixId);

FIXP_DBL
getDownmixOffset(DOWNMIX_INSTRUCTIONS* pDown, int baseChannelCount);

#endif
