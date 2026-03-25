/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: independent channel concealment

*******************************************************************************/

#ifndef CONCEAL_H
#define CONCEAL_H

#include "channelinfo.h"

#define AACDEC_CONCEAL_PARAM_NOT_SPECIFIED (0xFFFE)

void CConcealment_InitCommonData(CConcealParams *pConcealCommonData);

void CConcealment_InitChannelData(CConcealmentInfo *hConcealmentInfo,
                                  CConcealParams *pConcealCommonData,
                                  AACDEC_RENDER_MODE initRenderMode,
                                  int samplesPerFrame);

CConcealmentMethod CConcealment_GetMethod(CConcealParams *pConcealCommonData);

UINT CConcealment_GetDelay(CConcealParams *pConcealCommonData);

AAC_DECODER_ERROR
CConcealment_SetParams(CConcealParams *concealParams, int method,
                       int fadeOutSlope, int fadeInSlope, int muteRelease,
                       FIXP_DBL comfNoiseLevel);

CConcealmentState CConcealment_GetState(CConcealmentInfo *hConcealmentInfo);

AAC_DECODER_ERROR
CConcealment_SetAttenuation(CConcealParams *concealParams,
                            const SHORT *fadeOutAttenuationVector,
                            const SHORT *fadeInAttenuationVector);

void CConcealment_Store(
    CConcealmentInfo *hConcealmentInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo);

int CConcealment_Apply(
    CConcealmentInfo *hConcealmentInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const int samplesPerFrame,
    const UCHAR lastLpdMode, const int FrameOk, const UINT flags);

int CConcealment_GetLastFrameOk(CConcealmentInfo *hConcealmentInfo,
                                const int fBeforeApply);

INT CConcealment_TDFading(
    int len, CAacDecoderStaticChannelInfo **ppAacDecoderStaticChannelInfo,
    const INT aacOutDataHeadroom, PCM_DEC *pcmdata, PCM_DEC *pcmdata_1);

#endif /* #ifndef CONCEAL_H */
