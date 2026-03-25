/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description:

*******************************************************************************/

#ifndef CHANNEL_H
#define CHANNEL_H

#include "common_fix.h"

#include "FDK_tools_rom.h"
#include "channelinfo.h"
#include "tpdec_lib.h"

void CChannel_CodebookTableInit(CAacDecoderChannelInfo *pAacDecoderChannelInfo);

void CChannelElement_Decode(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo[2],
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[2],
    SamplingRateInfo *pSamplingRateInfo, UINT flags, UINT elFlags,
    int el_channels);

AAC_DECODER_ERROR CChannelElement_Read(
    HANDLE_FDK_BITSTREAM hBs, CAacDecoderChannelInfo *pAacDecoderChannelInfo[],
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[],
    const AUDIO_OBJECT_TYPE aot, SamplingRateInfo *pSamplingRateInfo,
    const UINT flags, const UINT elFlags, const UINT frame_length,
    const UCHAR numberOfChannels, const SCHAR epConfig,
    HANDLE_TRANSPORTDEC pTpDec);

#endif /* #ifndef CHANNEL_H */
