/**************************** AAC encoder library ******************************

   Author(s):   A. Groeschel

   Description: channel mapping functionality

*******************************************************************************/

#ifndef CHANNEL_MAP_H
#define CHANNEL_MAP_H

#include "aacenc.h"
#include "psy_const.h"
#include "qc_data.h"

typedef struct {
  CHANNEL_MODE encMode;
  INT nChannels;
  INT nChannelsEff;
  INT nElements;
} CHANNEL_MODE_CONFIG_TAB;

typedef enum { EL_MODE_INVALID = 0, EL_MODE_MONO, EL_MODE_STEREO } ELEMENT_MODE;

AAC_ENCODER_ERROR FDKaacEnc_DetermineEncoderMode(CHANNEL_MODE* mode,
                                                 INT nChannels);

AAC_ENCODER_ERROR FDKaacEnc_InitChannelMapping(CHANNEL_MODE mode,
                                               CHANNEL_ORDER co,
                                               CHANNEL_MAPPING* chMap);

AAC_ENCODER_ERROR FDKaacEnc_InitElementBits(QC_STATE* hQC, CHANNEL_MAPPING* cm,
                                            INT bitrateTot, INT averageBitsTot,
                                            INT maxChannelBits);

ELEMENT_MODE FDKaacEnc_GetMonoStereoMode(const CHANNEL_MODE mode);

const CHANNEL_MODE_CONFIG_TAB* FDKaacEnc_GetChannelModeConfiguration(
    const CHANNEL_MODE mode);

#endif /* CHANNEL_MAP_H */
