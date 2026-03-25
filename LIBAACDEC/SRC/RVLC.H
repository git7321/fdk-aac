/**************************** AAC decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef RVLC_H
#define RVLC_H

#include "aacdecoder.h"
#include "channel.h"
#include "rvlc_info.h"

#define RVLC_ERROR_ALL_ESCAPE_WORDS_INVALID \
  0x80000000
#define RVLC_ERROR_RVL_SUM_BIT_COUNTER_BELOW_ZERO_FWD \
  0x40000000
#define RVLC_ERROR_RVL_SUM_BIT_COUNTER_BELOW_ZERO_BWD \
  0x20000000
#define RVLC_ERROR_FORBIDDEN_CW_DETECTED_FWD \
  0x08000000
#define RVLC_ERROR_FORBIDDEN_CW_DETECTED_BWD \
  0x04000000

void CRvlc_Read(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                HANDLE_FDK_BITSTREAM bs);

void CRvlc_Decode(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                  CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
                  HANDLE_FDK_BITSTREAM bs);

void CRvlc_ElementCheck(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo[],
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[],
    const UINT flags, const INT elChannels);

#endif /* RVLC_H */
