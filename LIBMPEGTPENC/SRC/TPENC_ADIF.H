/******************* MPEG transport format encoder library *********************

   Author(s):   Alex Goeschel

   Description: Transport Headers support

*******************************************************************************/

#ifndef TPENC_ADIF_H
#define TPENC_ADIF_H

#include "machine_type.h"
#include "FDK_bitstream.h"

#include "tp_data.h"

typedef struct {
  CHANNEL_MODE cm;
  INT samplingRate;
  INT bitRate;
  int profile;
  int bVariableRate;
  int instanceTag;
  int headerWritten;
  int matrixMixdownA;
  int pseudoSurroundEnable;

} ADIF_INFO;

int adifWrite_EncodeHeader(ADIF_INFO *adif, HANDLE_FDK_BITSTREAM hBitStream,
                           INT adif_buffer_fullness);

int adifWrite_GetHeaderBits(ADIF_INFO *adif);

#endif /* TPENC_ADIF_H */
