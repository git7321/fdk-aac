/**************************** AAC decoder library ******************************

   Author(s):   Manuel Jander

   Description: USAC Linear Prediction Domain coding

*******************************************************************************/

#ifndef USACDEC_LPD_H
#define USACDEC_LPD_H

#include "channelinfo.h"

#define OPTIMIZE_AVG_PERFORMANCE

AAC_DECODER_ERROR CLpdChannelStream_Read(
    HANDLE_FDK_BITSTREAM hBs, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, UINT flags);

void CLpdChannelStream_Decode(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo, UINT flags);

AAC_DECODER_ERROR CLpd_RenderTimeSignal(
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo, PCM_DEC *pTimeData,
    INT samplesPerFrame, SamplingRateInfo *pSamplingRateInfo, UINT frameOk,
    const INT aacOutDataHeadroom, UINT flags, UINT strmFlags);

static inline INT CLpd_FAC_getLength(int fNotShortBlock, int fac_length_long) {
  if (fNotShortBlock) {
    return (fac_length_long);
  } else {
    return fac_length_long / 2;
  }
}

void filtLP(const FIXP_DBL *syn, PCM_DEC *syn_out, FIXP_DBL *noise,
            const FIXP_SGL *filt, const INT aacOutDataHeadroom, INT stop,
            int len);

void bass_pf_1sf_delay(FIXP_DBL syn[], const INT T_sf[], FIXP_DBL *pit_gain,
                       const int frame_length, const INT l_frame,
                       const INT l_next, PCM_DEC *synth_out,
                       const INT aacOutDataHeadroom, FIXP_DBL mem_bpf[]);

FDK_INLINE
int UsacRandomSign(ULONG *seed) {
  *seed = (ULONG)((UINT64)(*seed) * 69069 + 5);

  return (int)((*seed) & 0x10000);
}

void CFdp_Reset(CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo);

#endif /* USACDEC_LPD_H */
