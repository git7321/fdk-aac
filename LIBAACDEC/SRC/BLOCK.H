/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: long/short-block decoding

*******************************************************************************/

#ifndef BLOCK_H
#define BLOCK_H

#include "common_fix.h"

#include "channelinfo.h"
#include "FDK_bitstream.h"

void CPns_Read(CPnsData *pPnsData, HANDLE_FDK_BITSTREAM bs,
               const CodeBookDescription *hcb, SHORT *pScaleFactor,
               UCHAR global_gain, int band, int group);

void CPns_Apply(const CPnsData *pPnsData, const CIcsInfo *pIcsInfo,
                SPECTRAL_PTR pSpectrum, const SHORT *pSpecScale,
                const SHORT *pScaleFactor,
                const SamplingRateInfo *pSamplingRateInfo,
                const INT granuleLength, const int channel);

void CBlock_ApplyNoise(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                       SamplingRateInfo *pSamplingRateInfo, ULONG *nfRandomSeed,
                       UCHAR *band_is_noise);

void CTns_ReadDataPresentFlag(HANDLE_FDK_BITSTREAM bs, CTnsData *pTnsData);

void CTns_ReadDataPresentUsac(HANDLE_FDK_BITSTREAM hBs, CTnsData *pTnsData0,
                              CTnsData *pTnsData1, UCHAR *ptns_on_lr,
                              const CIcsInfo *pIcsInfo, const UINT flags,
                              const UINT elFlags, const int fCommonWindow);

AAC_DECODER_ERROR CTns_Read(HANDLE_FDK_BITSTREAM bs, CTnsData *pTnsData,
                            const CIcsInfo *pIcsInfo, const UINT flags);

void CTns_Apply(CTnsData *RESTRICT pTnsData,
                const CIcsInfo *pIcsInfo, SPECTRAL_PTR pSpectralCoefficient,
                const SamplingRateInfo *pSamplingRateInfo,
                const INT granuleLength, const UCHAR nbands,
                const UCHAR igf_active, const UINT flags);

LONG CBlock_GetEscape(HANDLE_FDK_BITSTREAM bs, const LONG q);

AAC_DECODER_ERROR CBlock_ReadScaleFactorData(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo, HANDLE_FDK_BITSTREAM bs,
    const UINT flags);

AAC_DECODER_ERROR CBlock_ReadSpectralData(
    HANDLE_FDK_BITSTREAM bs, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const UINT flags);

AAC_DECODER_ERROR CBlock_ReadAcSpectralData(
    HANDLE_FDK_BITSTREAM hBs, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const UINT frame_length,
    const UINT flags);

AAC_DECODER_ERROR CBlock_ReadSectionData(
    HANDLE_FDK_BITSTREAM bs, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const UINT flags);

void CBlock_ScaleSpectralData(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                              UCHAR maxSfbs,
                              SamplingRateInfo *pSamplingRateInfo);

void ApplyTools(CAacDecoderChannelInfo *pAacDecoderChannelInfo[],
                const SamplingRateInfo *pSamplingRateInfo, const UINT flags,
                const UINT elFlags, const int channel, const int maybe_jstereo);

void CBlock_FrequencyToTime(
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo, PCM_DEC outSamples[],
    const SHORT frameLen, const int frameOk, FIXP_DBL *pWorkBuffer1,
    const INT aacOutDataHeadroom, UINT elFlags, INT elCh);

void CBlock_FrequencyToTimeLowDelay(
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo, PCM_DEC outSamples[],
    const short frameLen);

AAC_DECODER_ERROR CBlock_InverseQuantizeSpectralData(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    SamplingRateInfo *pSamplingRateInfo, UCHAR *band_is_noise,
    UCHAR active_band_search);

FDK_INLINE
int EvaluatePower43(FIXP_DBL *pValue, UINT lsb) {
  FIXP_DBL value;
  UINT freeBits;
  UINT exponent;

  value = *pValue;
  freeBits = fNormz(value);
  exponent = DFRACT_BITS - freeBits;
  FDK_ASSERT(exponent < 14);

  UINT x = (((int)value << freeBits) >> 19);
  UINT tableIndex = (x & 0x0FFF) >> 4;
  FIXP_DBL invQVal;

  x = x & 0x0F;

  UINT r0 = (LONG)InverseQuantTable[tableIndex + 0];
  UINT r1 = (LONG)InverseQuantTable[tableIndex + 1];
  USHORT nx = 16 - x;
  UINT temp = (r0)*nx + (r1)*x;
  invQVal = (FIXP_DBL)temp;

  FDK_ASSERT(lsb < 4);
  *pValue = fMultDiv2(invQVal, MantissaTable[lsb][exponent]);

  return ExponentTable[lsb][exponent] + 1;
}

FIXP_DBL get_gain(const FIXP_DBL *x, const FIXP_DBL *y, int n);

FDK_INLINE int GetScaleFromValue(FIXP_DBL value, unsigned int lsb) {
  if (value != (FIXP_DBL)0) {
    int scale = EvaluatePower43(&value, lsb);
    return CntLeadingZeros(value) - scale - 2;
  } else
    return 0;
}

inline int CBlock_DecodeHuffmanWord(
    HANDLE_FDK_BITSTREAM bs,
    const CodeBookDescription *hcb)
{
  UINT val;
  UINT index = 0;
  const USHORT(*CodeBook)[HuffmanEntries] = hcb->CodeBook;

  while (1) {
    val = CodeBook[index]
                  [FDKreadBits(bs, HuffmanBits)];

    if ((val & 1) == 0) {
      index = val >> 2;
      continue;
    } else {
      if (val & 2) {
        FDKpushBackCache(bs, 1);
      }

      val >>= 2;
      break;
    }
  }

  return val;
}
inline int CBlock_DecodeHuffmanWordCB(
    HANDLE_FDK_BITSTREAM bs,
    const USHORT (
        *CodeBook)[HuffmanEntries])
{
  UINT index = 0;

  while (1) {
    index = CodeBook[index][FDKread2Bits(bs)];
    if (index & 1) break;
    index >>= 2;
  }
  if (index & 2) {
    FDKpushBackCache(bs, 1);
  }
  return index >> 2;
}

#endif /* #ifndef BLOCK_H */
