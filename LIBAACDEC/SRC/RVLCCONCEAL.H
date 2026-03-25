/**************************** AAC decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef RVLCCONCEAL_H
#define RVLCCONCEAL_H

#include "rvlc.h"

void BidirectionalEstimation_UseLowerScfOfCurrentFrame(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo);

void BidirectionalEstimation_UseScfOfPrevFrameAsReference(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo);

void StatisticalEstimation(CAacDecoderChannelInfo *pAacDecoderChannelInfo);

void PredictiveInterpolation(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo);

#endif /* RVLCCONCEAL_H */
