/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef CMONDATA_H
#define CMONDATA_H

#include "FDK_bitstream.h"

struct COMMON_DATA {
  INT sbrHdrBits;
  INT sbrDataBits;
  INT sbrFillBits;
  FDK_BITSTREAM sbrBitbuf;
  FDK_BITSTREAM tmpWriteBitbuf;
  INT xOverFreq;
  INT dynBwEnabled;
  INT sbrNumChannels;
  INT dynXOverFreqEnc;
};

typedef struct COMMON_DATA *HANDLE_COMMON_DATA;

#endif
