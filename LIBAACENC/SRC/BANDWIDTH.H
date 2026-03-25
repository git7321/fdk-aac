/**************************** AAC encoder library ******************************

   Author(s):   M. Schug / A. Groeschel

   Description: bandwidth expert

*******************************************************************************/

#ifndef BANDWIDTH_H
#define BANDWIDTH_H

#include "qc_data.h"

AAC_ENCODER_ERROR FDKaacEnc_DetermineBandWidth(
    const INT proposedBandWidth, const INT bitrate,
    const AACENC_BITRATE_MODE bitrateMode, const INT sampleRate,
    const INT frameLength, const CHANNEL_MAPPING *const cm,
    const CHANNEL_MODE encoderMode, INT *const bandWidth);

#endif /* BANDWIDTH_H */
