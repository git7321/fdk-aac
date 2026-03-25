/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: long/short-block decoding

*******************************************************************************/

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "block.h"
#include "aac_rom.h"
#include "FDK_bitstream.h"
#include "scale.h"
#include "FDK_tools_rom.h"
#include "usacdec_fac.h"
#include "usacdec_lpd.h"
#include "usacdec_lpc.h"
#include "FDK_trigFcts.h"
#include "ac_arith_coder.h"
#include "aacdec_hcr.h"
#include "rvlc.h"

#if defined(__arm__)
#include "arm/block_arm.cpp"
#endif

LONG CBlock_GetEscape(HANDLE_FDK_BITSTREAM bs,
                      const LONG q)
{
  if (fAbs(q) != 16) return (q);

  LONG i, off;
  for (i = 4; i < 13; i++) {
    if (FDKreadBit(bs) == 0) break;
  }

  if (i == 13) return (MAX_QUANTIZED_VALUE + 1);

  off = FDKreadBits(bs, i);
  i = off + (1 << i);

  if (q < 0) i = -i;

  return i;
}

AAC_DECODER_ERROR CBlock_ReadScaleFactorData(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo, HANDLE_FDK_BITSTREAM bs,
    UINT flags) {
  int temp;
  int band;
  int group;
  int position = 0;
  int factor = pAacDecoderChannelInfo->pDynData->RawDataInfo
                   .GlobalGain;
  UCHAR *pCodeBook = pAacDecoderChannelInfo->pDynData->aCodeBook;
  SHORT *pScaleFactor = pAacDecoderChannelInfo->pDynData->aScaleFactor;
  const CodeBookDescription *hcb = &AACcodeBookDescriptionTable[BOOKSCL];

  const USHORT(*CodeBook)[HuffmanEntries] = hcb->CodeBook;

  int ScaleFactorBandsTransmitted =
      GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo->icsInfo);
  for (group = 0; group < GetWindowGroups(&pAacDecoderChannelInfo->icsInfo);
       group++) {
    for (band = 0; band < ScaleFactorBandsTransmitted; band++) {
      switch (pCodeBook[band]) {
        case ZERO_HCB:
          pScaleFactor[band] = 0;
          break;

        default:
          if (!((flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA)) && band == 0 &&
                group == 0)) {
            temp = CBlock_DecodeHuffmanWordCB(bs, CodeBook);
            factor += temp - 60;
          }
          pScaleFactor[band] = factor - 100;
          break;

        case INTENSITY_HCB:
        case INTENSITY_HCB2:
          temp = CBlock_DecodeHuffmanWordCB(bs, CodeBook);
          position += temp - 60;
          pScaleFactor[band] = position - 100;
          break;

        case NOISE_HCB:
          if (flags & (AC_MPEGD_RES | AC_USAC | AC_RSVD50 | AC_RSV603DA)) {
            return AAC_DEC_PARSE_ERROR;
          }
          CPns_Read(&pAacDecoderChannelInfo->data.aac.PnsData, bs, hcb,
                    pAacDecoderChannelInfo->pDynData->aScaleFactor,
                    pAacDecoderChannelInfo->pDynData->RawDataInfo.GlobalGain,
                    band, group);
          break;
      }
    }
    pCodeBook += 16;
    pScaleFactor += 16;
  }

  return AAC_DEC_OK;
}

void CBlock_ScaleSpectralData(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                              UCHAR maxSfbs,
                              SamplingRateInfo *pSamplingRateInfo) {
  int band;
  int window;
  const SHORT *RESTRICT pSfbScale = pAacDecoderChannelInfo->pDynData->aSfbScale;
  SHORT *RESTRICT pSpecScale = pAacDecoderChannelInfo->specScale;
  int groupwin, group;
  const SHORT *RESTRICT BandOffsets = GetScaleFactorBandOffsets(
      &pAacDecoderChannelInfo->icsInfo, pSamplingRateInfo);
  SPECTRAL_PTR RESTRICT pSpectralCoefficient =
      pAacDecoderChannelInfo->pSpectralCoefficient;

  FDKmemclear(pSpecScale, 8 * sizeof(SHORT));

  for (window = 0, group = 0;
       group < GetWindowGroups(&pAacDecoderChannelInfo->icsInfo); group++) {
    for (groupwin = 0; groupwin < GetWindowGroupLength(
                                      &pAacDecoderChannelInfo->icsInfo, group);
         groupwin++, window++) {
      int SpecScale_window = pSpecScale[window];
      FIXP_DBL *pSpectrum = SPEC(pSpectralCoefficient, window,
                                 pAacDecoderChannelInfo->granuleLength);

      for (band = 0; band < maxSfbs; band++) {
        SpecScale_window =
            fMax(SpecScale_window, (int)pSfbScale[window * 16 + band]);
      }

      if (pAacDecoderChannelInfo->pDynData->TnsData.Active &&
          pAacDecoderChannelInfo->pDynData->TnsData.NumberOfFilters[window] >
              0) {
        int filter_index, SpecScale_window_tns;
        int tns_start, tns_stop;

        SpecScale_window_tns = 0;
        tns_start = GetMaximumTnsBands(&pAacDecoderChannelInfo->icsInfo,
                                       pSamplingRateInfo->samplingRateIndex);
        tns_stop = 0;
        for (filter_index = 0;
             filter_index < (int)pAacDecoderChannelInfo->pDynData->TnsData
                                .NumberOfFilters[window];
             filter_index++) {
          for (band = pAacDecoderChannelInfo->pDynData->TnsData
                          .Filter[window][filter_index]
                          .StartBand;
               band < pAacDecoderChannelInfo->pDynData->TnsData
                          .Filter[window][filter_index]
                          .StopBand;
               band++) {
            SpecScale_window_tns =
                fMax(SpecScale_window_tns, (int)pSfbScale[window * 16 + band]);
          }
          tns_start =
              fMin(tns_start, (int)pAacDecoderChannelInfo->pDynData->TnsData
                                  .Filter[window][filter_index]
                                  .StartBand);
          tns_stop =
              fMax(tns_stop, (int)pAacDecoderChannelInfo->pDynData->TnsData
                                 .Filter[window][filter_index]
                                 .StopBand);
        }
        SpecScale_window_tns = SpecScale_window_tns +
                               pAacDecoderChannelInfo->pDynData->TnsData.GainLd;
        FDK_ASSERT(tns_stop >= tns_start);
        SpecScale_window_tns -=
            getScalefactor(pSpectrum + BandOffsets[tns_start],
                           BandOffsets[tns_stop] - BandOffsets[tns_start]);
        if (SpecScale_window <= 17) {
          SpecScale_window_tns++;
        }
        SpecScale_window = fMax(SpecScale_window, SpecScale_window_tns);
      }

      pSpecScale[window] = SpecScale_window;

#ifdef FUNCTION_CBlock_ScaleSpectralData_func1

      CBlock_ScaleSpectralData_func1(pSpectrum, maxSfbs, BandOffsets,
                                     SpecScale_window, pSfbScale, window);

#else
      for (band = 0; band < maxSfbs; band++) {
        int scale = fMin(DFRACT_BITS - 1,
                         SpecScale_window - pSfbScale[window * 16 + band]);
        if (scale) {
          FDK_ASSERT(scale > 0);

          int max_index = BandOffsets[band + 1];
          DWORD_ALIGNED(pSpectrum);
          for (int index = BandOffsets[band]; index < max_index; index++) {
            pSpectrum[index] >>= scale;
          }
        }
      }
#endif
    }
  }
}

AAC_DECODER_ERROR CBlock_ReadSectionData(
    HANDLE_FDK_BITSTREAM bs, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const UINT flags) {
  int top, band;
  int sect_len, sect_len_incr;
  int group;
  UCHAR sect_cb;
  UCHAR *pCodeBook = pAacDecoderChannelInfo->pDynData->aCodeBook;
  SHORT *pNumLinesInSec =
      pAacDecoderChannelInfo->pDynData->specificTo.aac.aNumLineInSec4Hcr;
  int numLinesInSecIdx = 0;
  UCHAR *pHcrCodeBook =
      pAacDecoderChannelInfo->pDynData->specificTo.aac.aCodeBooks4Hcr;
  const SHORT *BandOffsets = GetScaleFactorBandOffsets(
      &pAacDecoderChannelInfo->icsInfo, pSamplingRateInfo);
  pAacDecoderChannelInfo->pDynData->specificTo.aac.numberSection = 0;
  AAC_DECODER_ERROR ErrorStatus = AAC_DEC_OK;

  FDKmemclear(pCodeBook, sizeof(UCHAR) * (8 * 16));

  const int nbits =
      (IsLongBlock(&pAacDecoderChannelInfo->icsInfo) == 1) ? 5 : 3;

  int sect_esc_val = (1 << nbits) - 1;

  UCHAR ScaleFactorBandsTransmitted =
      GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo->icsInfo);
  for (group = 0; group < GetWindowGroups(&pAacDecoderChannelInfo->icsInfo);
       group++) {
    for (band = 0; band < ScaleFactorBandsTransmitted;) {
      sect_len = 0;
      if (flags & AC_ER_VCB11) {
        sect_cb = (UCHAR)FDKreadBits(bs, 5);
      } else
        sect_cb = (UCHAR)FDKreadBits(bs, 4);

      if (((flags & AC_ER_VCB11) == 0) || (sect_cb < 11) ||
          ((sect_cb > 11) && (sect_cb < 16))) {
        sect_len_incr = FDKreadBits(bs, nbits);
        while (sect_len_incr == sect_esc_val) {
          sect_len += sect_esc_val;
          sect_len_incr = FDKreadBits(bs, nbits);
        }
      } else {
        sect_len_incr = 1;
      }

      sect_len += sect_len_incr;

      top = band + sect_len;

      if (flags & AC_ER_HCR) {
        if (numLinesInSecIdx >= MAX_SFB_HCR) {
          return AAC_DEC_PARSE_ERROR;
        }
        if (top > (int)GetNumberOfScaleFactorBands(
                      &pAacDecoderChannelInfo->icsInfo, pSamplingRateInfo)) {
          return AAC_DEC_PARSE_ERROR;
        }
        pNumLinesInSec[numLinesInSecIdx] = BandOffsets[top] - BandOffsets[band];
        numLinesInSecIdx++;
        if (sect_cb == BOOKSCL) {
          return AAC_DEC_INVALID_CODE_BOOK;
        } else {
          *pHcrCodeBook++ = sect_cb;
        }
        pAacDecoderChannelInfo->pDynData->specificTo.aac.numberSection++;
      }

      if (IsLongBlock(&(pAacDecoderChannelInfo->icsInfo))) {
        if (top > 64) {
          return AAC_DEC_DECODE_FRAME_ERROR;
        }
      } else {
        if (top + group * 16 > (8 * 16)) {
          return AAC_DEC_DECODE_FRAME_ERROR;
        }
      }

      if ((sect_cb == BOOKSCL) ||
          ((sect_cb == INTENSITY_HCB || sect_cb == INTENSITY_HCB2) &&
           pAacDecoderChannelInfo->pDynData->RawDataInfo.CommonWindow == 0)) {
        return AAC_DEC_INVALID_CODE_BOOK;
      }

      for (; band < top; band++) {
        pCodeBook[group * 16 + band] = sect_cb;
      }
    }
  }

  return ErrorStatus;
}

static inline void InverseQuantizeBand(
    FIXP_DBL *RESTRICT spectrum, const FIXP_DBL *RESTRICT InverseQuantTabler,
    const FIXP_DBL *RESTRICT MantissaTabler,
    const SCHAR *RESTRICT ExponentTabler, INT noLines, INT scale) {
  scale = scale + 1;

  FIXP_DBL *RESTRICT ptr = spectrum;
  FIXP_DBL signedValue;

  for (INT i = noLines; i--;) {
    if ((signedValue = *ptr++) != FL2FXCONST_DBL(0)) {
      FIXP_DBL value = fAbs(signedValue);
      UINT freeBits = CntLeadingZeros(value);
      UINT exponent = 32 - freeBits;

      UINT x = (UINT)(LONG)value << (INT)freeBits;
      x <<= 1;
      UINT tableIndex = x >> 24;
      x = (x >> 20) & 0x0F;

      UINT r0 = (UINT)(LONG)InverseQuantTabler[tableIndex + 0];
      UINT r1 = (UINT)(LONG)InverseQuantTabler[tableIndex + 1];
      UINT temp = (r1 - r0) * x + (r0 << 4);

      value = fMultDiv2((FIXP_DBL)temp, MantissaTabler[exponent]);

      scaleValueInPlace(&value, scale + ExponentTabler[exponent]);

      signedValue = (signedValue < (FIXP_DBL)0) ? -value : value;
      ptr[-1] = signedValue;
    }
  }
}

static inline FIXP_DBL maxabs_D(const FIXP_DBL *pSpectralCoefficient,
                                const int noLines) {
  FIXP_DBL locMax = (FIXP_DBL)0;
  int i;

  DWORD_ALIGNED(pSpectralCoefficient);

  for (i = noLines; i-- > 0;) {
    locMax = fMax(fixp_abs(pSpectralCoefficient[i]), locMax);
  }

  return locMax;
}

AAC_DECODER_ERROR CBlock_InverseQuantizeSpectralData(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    SamplingRateInfo *pSamplingRateInfo, UCHAR *band_is_noise,
    UCHAR active_band_search) {
  int window, group, groupwin, band;
  int ScaleFactorBandsTransmitted =
      GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo->icsInfo);
  UCHAR *RESTRICT pCodeBook = pAacDecoderChannelInfo->pDynData->aCodeBook;
  SHORT *RESTRICT pSfbScale = pAacDecoderChannelInfo->pDynData->aSfbScale;
  SHORT *RESTRICT pScaleFactor = pAacDecoderChannelInfo->pDynData->aScaleFactor;
  const SHORT *RESTRICT BandOffsets = GetScaleFactorBandOffsets(
      &pAacDecoderChannelInfo->icsInfo, pSamplingRateInfo);
  const SHORT total_bands =
      GetScaleFactorBandsTotal(&pAacDecoderChannelInfo->icsInfo);

  FDKmemclear(pAacDecoderChannelInfo->pDynData->aSfbScale,
              (8 * 16) * sizeof(SHORT));

  for (window = 0, group = 0;
       group < GetWindowGroups(&pAacDecoderChannelInfo->icsInfo); group++) {
    for (groupwin = 0; groupwin < GetWindowGroupLength(
                                      &pAacDecoderChannelInfo->icsInfo, group);
         groupwin++, window++) {
      for (band = 0; band < ScaleFactorBandsTransmitted; band++) {
        FIXP_DBL *pSpectralCoefficient =
            SPEC(pAacDecoderChannelInfo->pSpectralCoefficient, window,
                 pAacDecoderChannelInfo->granuleLength) +
            BandOffsets[band];
        FIXP_DBL locMax;

        const int noLines = BandOffsets[band + 1] - BandOffsets[band];
        const int bnds = group * 16 + band;

        if ((pCodeBook[bnds] == ZERO_HCB) ||
            (pCodeBook[bnds] == INTENSITY_HCB) ||
            (pCodeBook[bnds] == INTENSITY_HCB2))
          continue;

        if (pCodeBook[bnds] == NOISE_HCB) {
          pSfbScale[window * 16 + band] = (pScaleFactor[bnds] >> 2) + 1;
          continue;
        }

        locMax = maxabs_D(pSpectralCoefficient, noLines);

        if (active_band_search) {
          if (locMax != FIXP_DBL(0)) {
            band_is_noise[group * 16 + band] = 0;
          }
        }

        if (fixp_abs(locMax) > (FIXP_DBL)MAX_QUANTIZED_VALUE) {
          return AAC_DEC_PARSE_ERROR;
        }

        int msb = pScaleFactor[bnds] >> 2;

        if (locMax != FIXP_DBL(0)) {
          int lsb = pScaleFactor[bnds] & 0x03;

          int scale = EvaluatePower43(&locMax, lsb);

          scale = CntLeadingZeros(locMax) - scale - 2;

          pSfbScale[window * 16 + band] = msb - scale;
          InverseQuantizeBand(pSpectralCoefficient, InverseQuantTable,
                              MantissaTable[lsb], ExponentTable[lsb], noLines,
                              scale);
        } else {
          pSfbScale[window * 16 + band] = msb;
        }
      }

      SHORT start_clear = BandOffsets[ScaleFactorBandsTransmitted];
      SHORT end_clear = BandOffsets[total_bands];
      int diff_clear = (int)(end_clear - start_clear);
      FIXP_DBL *pSpectralCoefficient =
          SPEC(pAacDecoderChannelInfo->pSpectralCoefficient, window,
               pAacDecoderChannelInfo->granuleLength) +
          start_clear;
      FDKmemclear(pSpectralCoefficient, diff_clear * sizeof(FIXP_DBL));

    }
  }

  return AAC_DEC_OK;
}

AAC_DECODER_ERROR CBlock_ReadSpectralData(
    HANDLE_FDK_BITSTREAM bs, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const UINT flags) {
  int index, i;
  const SHORT *RESTRICT BandOffsets = GetScaleFactorBandOffsets(
      &pAacDecoderChannelInfo->icsInfo, pSamplingRateInfo);

  SPECTRAL_PTR pSpectralCoefficient =
      pAacDecoderChannelInfo->pSpectralCoefficient;

  FDK_ASSERT(BandOffsets != NULL);

  FDKmemclear(pSpectralCoefficient, sizeof(SPECTRUM));

  if ((flags & AC_ER_HCR) == 0) {
    int group;
    int groupoffset;
    UCHAR *pCodeBook = pAacDecoderChannelInfo->pDynData->aCodeBook;
    int ScaleFactorBandsTransmitted =
        GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo->icsInfo);
    int granuleLength = pAacDecoderChannelInfo->granuleLength;

    groupoffset = 0;

    int max_group = GetWindowGroups(&pAacDecoderChannelInfo->icsInfo);

    for (group = 0; group < max_group; group++) {
      int max_groupwin =
          GetWindowGroupLength(&pAacDecoderChannelInfo->icsInfo, group);
      int band;

      int bnds = group * 16;

      int bandOffset1 = BandOffsets[0];
      for (band = 0; band < ScaleFactorBandsTransmitted; band++, bnds++) {
        UCHAR currentCB = pCodeBook[bnds];
        int bandOffset0 = bandOffset1;
        bandOffset1 = BandOffsets[band + 1];

        if ((currentCB >= 16) && (currentCB <= 31)) {
          pCodeBook[bnds] = currentCB = 11;
        }
        if (((currentCB != ZERO_HCB) && (currentCB != NOISE_HCB) &&
             (currentCB != INTENSITY_HCB) && (currentCB != INTENSITY_HCB2))) {
          const CodeBookDescription *hcb =
              &AACcodeBookDescriptionTable[currentCB];
          int step = hcb->Dimension;
          int offset = hcb->Offset;
          int bits = hcb->numBits;
          int mask = (1 << bits) - 1;
          const USHORT(*CodeBook)[HuffmanEntries] = hcb->CodeBook;
          int groupwin;

          FIXP_DBL *mdctSpectrum =
              &pSpectralCoefficient[groupoffset * granuleLength];

          if (offset == 0) {
            for (groupwin = 0; groupwin < max_groupwin; groupwin++) {
              for (index = bandOffset0; index < bandOffset1; index += step) {
                int idx = CBlock_DecodeHuffmanWordCB(bs, CodeBook);
                for (i = 0; i < step; i++, idx >>= bits) {
                  FIXP_DBL tmp = (FIXP_DBL)((idx & mask) - offset);
                  if (tmp != FIXP_DBL(0)) tmp = (FDKreadBit(bs)) ? -tmp : tmp;
                  mdctSpectrum[index + i] = tmp;
                }

                if (currentCB == ESCBOOK) {
                  for (int j = 0; j < 2; j++)
                    mdctSpectrum[index + j] = (FIXP_DBL)CBlock_GetEscape(
                        bs, (LONG)mdctSpectrum[index + j]);
                }
              }
              mdctSpectrum += granuleLength;
            }
          } else {
            for (groupwin = 0; groupwin < max_groupwin; groupwin++) {
              for (index = bandOffset0; index < bandOffset1; index += step) {
                int idx = CBlock_DecodeHuffmanWordCB(bs, CodeBook);
                for (i = 0; i < step; i++, idx >>= bits) {
                  mdctSpectrum[index + i] = (FIXP_DBL)((idx & mask) - offset);
                }
                if (currentCB == ESCBOOK) {
                  for (int j = 0; j < 2; j++)
                    mdctSpectrum[index + j] = (FIXP_DBL)CBlock_GetEscape(
                        bs, (LONG)mdctSpectrum[index + j]);
                }
              }
              mdctSpectrum += granuleLength;
            }
          }
        }
      }
      groupoffset += max_groupwin;
    }
  }

  else

  {
    H_HCR_INFO hHcr = &pAacDecoderChannelInfo->pComData->overlay.aac.erHcrInfo;

    int hcrStatus = 0;

    if (pAacDecoderChannelInfo->pDynData->specificTo.aac
            .lenOfReorderedSpectralData != 0) {
      hcrStatus = HcrInit(hHcr, pAacDecoderChannelInfo, pSamplingRateInfo, bs);

      if (hcrStatus != 0) {
        return AAC_DEC_DECODE_FRAME_ERROR;
      }

      hcrStatus =
          HcrDecoder(hHcr, pAacDecoderChannelInfo, pSamplingRateInfo, bs);

      if (hcrStatus != 0) {
#if HCR_ERROR_CONCEALMENT
        HcrMuteErroneousLines(hHcr);
#else
        return AAC_DEC_DECODE_FRAME_ERROR;
#endif
      }

      FDKpushFor(bs, pAacDecoderChannelInfo->pDynData->specificTo.aac
                         .lenOfReorderedSpectralData);
    }
  }

  if (IsLongBlock(&pAacDecoderChannelInfo->icsInfo) &&
      !(flags & (AC_ELD | AC_SCALABLE))) {
    CPulseData_Apply(
        &pAacDecoderChannelInfo->pDynData->specificTo.aac.PulseData,
        GetScaleFactorBandOffsets(&pAacDecoderChannelInfo->icsInfo,
                                  pSamplingRateInfo),
        SPEC_LONG(pSpectralCoefficient));
  }

  return AAC_DEC_OK;
}

static const FIXP_SGL noise_level_tab[8] = {
    FX_DBL2FXCONST_SGL(0x00000000),
    FX_DBL2FXCONST_SGL(0x0cb2ff5e),
    FX_DBL2FXCONST_SGL(0x10000000),
    FX_DBL2FXCONST_SGL(0x1428a2e7),
    FX_DBL2FXCONST_SGL(0x1965febd),
    FX_DBL2FXCONST_SGL(0x20000000),
    FX_DBL2FXCONST_SGL(0x28514606),
    FX_DBL2FXCONST_SGL(0x32cbfd33)};

void CBlock_ApplyNoise(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                       SamplingRateInfo *pSamplingRateInfo, ULONG *nfRandomSeed,
                       UCHAR *band_is_noise) {
  const SHORT *swb_offset = GetScaleFactorBandOffsets(
      &pAacDecoderChannelInfo->icsInfo, pSamplingRateInfo);
  int g, win, gwin, sfb, noiseFillingStartOffset, nfStartOffset_sfb;
  int noise_level = pAacDecoderChannelInfo->pDynData->specificTo.usac
                        .fd_noise_level_and_offset >>
                    5;
  const FIXP_SGL noiseVal_pos = noise_level_tab[noise_level];

  const int noise_offset = (pAacDecoderChannelInfo->pDynData->specificTo.usac
                                .fd_noise_level_and_offset &
                            0x1f) -
                           16;

  int max_sfb =
      GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo->icsInfo);

  noiseFillingStartOffset =
      (GetWindowSequence(&pAacDecoderChannelInfo->icsInfo) == BLOCK_SHORT)
          ? 20
          : 160;
  if (pAacDecoderChannelInfo->granuleLength == 96) {
    noiseFillingStartOffset =
        (3 * noiseFillingStartOffset) /
        4;
  }

  for (sfb = 0; swb_offset[sfb] < noiseFillingStartOffset; sfb++)
    ;
  nfStartOffset_sfb = sfb;

  {
    for (g = 0, win = 0; g < GetWindowGroups(&pAacDecoderChannelInfo->icsInfo);
         g++) {
      int windowGroupLength =
          GetWindowGroupLength(&pAacDecoderChannelInfo->icsInfo, g);
      for (sfb = nfStartOffset_sfb; sfb < max_sfb; sfb++) {
        int bin_start = swb_offset[sfb];
        int bin_stop = swb_offset[sfb + 1];

        int flagN = band_is_noise[g * 16 + sfb];

        if (flagN) {
          pAacDecoderChannelInfo->pDynData->aScaleFactor[g * 16 + sfb] +=
              noise_offset;
          for (gwin = 0; gwin < windowGroupLength; gwin++) {
            pAacDecoderChannelInfo->pDynData
                ->aSfbScale[(win + gwin) * 16 + sfb] += (noise_offset >> 2);
          }
        }

        ULONG seed = *nfRandomSeed;
        int scale =
            (pAacDecoderChannelInfo->pDynData->aScaleFactor[g * 16 + sfb] >>
             2) +
            1;
        int lsb =
            pAacDecoderChannelInfo->pDynData->aScaleFactor[g * 16 + sfb] & 3;
        FIXP_DBL mantissa = MantissaTable[lsb][0];

        for (gwin = 0; gwin < windowGroupLength; gwin++) {
          FIXP_DBL *pSpec =
              SPEC(pAacDecoderChannelInfo->pSpectralCoefficient, win + gwin,
                   pAacDecoderChannelInfo->granuleLength);

          int scale1 = scale - pAacDecoderChannelInfo->pDynData
                                   ->aSfbScale[(win + gwin) * 16 + sfb];
          FIXP_DBL scaled_noiseVal_pos =
              scaleValue(fMultDiv2(noiseVal_pos, mantissa), scale1);
          FIXP_DBL scaled_noiseVal_neg = -scaled_noiseVal_pos;

          if (flagN) {
            for (int bin = bin_start; bin < bin_stop; bin++) {
              seed = (ULONG)(
                  (UINT64)seed * 69069 +
                  5);
              pSpec[bin] =
                  (seed & 0x10000) ? scaled_noiseVal_neg : scaled_noiseVal_pos;
            }
          }
          else {
            for (int bin = bin_start; bin < bin_stop; bin++) {
              if (pSpec[bin] == (FIXP_DBL)0) {
                seed = (ULONG)(
                    (UINT64)seed * 69069 +
                    5);
                pSpec[bin] = (seed & 0x10000) ? scaled_noiseVal_neg
                                              : scaled_noiseVal_pos;
              }
            }
          }
        }
        *nfRandomSeed = seed;
      }
      win += windowGroupLength;
    }
  }
}

AAC_DECODER_ERROR CBlock_ReadAcSpectralData(
    HANDLE_FDK_BITSTREAM hBs, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, const UINT frame_length,
    const UINT flags) {
  AAC_DECODER_ERROR errorAAC = AAC_DEC_OK;
  ARITH_CODING_ERROR error = ARITH_CODER_OK;
  int arith_reset_flag, lg, numWin, win, winLen;
  const SHORT *RESTRICT BandOffsets;

  BandOffsets = GetScaleFactorBandOffsets(&pAacDecoderChannelInfo->icsInfo,
                                          pSamplingRateInfo);
  lg = BandOffsets[GetScaleFactorBandsTransmitted(
      &pAacDecoderChannelInfo->icsInfo)];

  numWin = GetWindowsPerFrame(&pAacDecoderChannelInfo->icsInfo);
  winLen = (IsLongBlock(&pAacDecoderChannelInfo->icsInfo))
               ? (int)frame_length
               : (int)frame_length / numWin;

  if (flags & AC_INDEP) {
    arith_reset_flag = 1;
  } else {
    arith_reset_flag = (USHORT)FDKreadBits(hBs, 1);
  }

  for (win = 0; win < numWin; win++) {
    error =
        CArco_DecodeArithData(pAacDecoderStaticChannelInfo->hArCo, hBs,
                              SPEC(pAacDecoderChannelInfo->pSpectralCoefficient,
                                   win, pAacDecoderChannelInfo->granuleLength),
                              lg, winLen, arith_reset_flag && (win == 0));
    if (error != ARITH_CODER_OK) {
      goto bail;
    }
  }

bail:
  if (error == ARITH_CODER_ERROR) {
    errorAAC = AAC_DEC_PARSE_ERROR;
  }

  return errorAAC;
}

void ApplyTools(CAacDecoderChannelInfo *pAacDecoderChannelInfo[],
                const SamplingRateInfo *pSamplingRateInfo, const UINT flags,
                const UINT elFlags, const int channel,
                const int common_window) {
  if (!(flags & (AC_USAC | AC_RSVD50 | AC_MPEGD_RES | AC_RSV603DA))) {
    CPns_Apply(&pAacDecoderChannelInfo[channel]->data.aac.PnsData,
               &pAacDecoderChannelInfo[channel]->icsInfo,
               pAacDecoderChannelInfo[channel]->pSpectralCoefficient,
               pAacDecoderChannelInfo[channel]->specScale,
               pAacDecoderChannelInfo[channel]->pDynData->aScaleFactor,
               pSamplingRateInfo,
               pAacDecoderChannelInfo[channel]->granuleLength, channel);
  }

  UCHAR nbands =
      GetScaleFactorBandsTransmitted(&pAacDecoderChannelInfo[channel]->icsInfo);

  CTns_Apply(&pAacDecoderChannelInfo[channel]->pDynData->TnsData,
             &pAacDecoderChannelInfo[channel]->icsInfo,
             pAacDecoderChannelInfo[channel]->pSpectralCoefficient,
             pSamplingRateInfo, pAacDecoderChannelInfo[channel]->granuleLength,
             nbands, (elFlags & AC_EL_ENHANCED_NOISE) ? 1 : 0, flags);
}

static int getWindow2Nr(int length, int shape) {
  int nr = 0;

  if (shape == 2) {
    nr = (length * 3) >> 2;
  }

  return nr;
}

FIXP_DBL get_gain(const FIXP_DBL *x, const FIXP_DBL *y, int n) {
  FIXP_DBL corr = (FIXP_DBL)0;
  FIXP_DBL ener = (FIXP_DBL)1;

  int headroom_x = getScalefactor(x, n);
  int headroom_y = getScalefactor(y, n);

  INT width_shift = (INT)(fNormz((FIXP_DBL)n));
  width_shift = 31 - width_shift;

  for (int i = 0; i < n; i++) {
    corr +=
        fMultDiv2((x[i] << headroom_x), (y[i] << headroom_y)) >> width_shift;
    ener += fPow2Div2((y[i] << headroom_y)) >> width_shift;
  }

  int exp_corr = (17 - headroom_x) + (17 - headroom_y) + width_shift + 1;
  int exp_ener = ((17 - headroom_y) << 1) + width_shift + 1;

  int temp_exp = 0;
  FIXP_DBL output = fDivNormSigned(corr, ener, &temp_exp);

  int output_exp = (exp_corr - exp_ener) + temp_exp;

  INT output_shift = 17 - output_exp;
  output_shift = fMin(output_shift, 31);

  output = scaleValue(output, -output_shift);

  return output;
}

void CBlock_FrequencyToTime(
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo, PCM_DEC outSamples[],
    const SHORT frameLen, const int frameOk, FIXP_DBL *pWorkBuffer1,
    const INT aacOutDataHeadroom, UINT elFlags, INT elCh) {
  int fr, fl, tl, nSpec;

#if defined(FDK_ASSERT_ENABLE)
  LONG nSamples;
#endif

  tl = frameLen;
  nSpec = 1;

  switch (pAacDecoderChannelInfo->icsInfo.WindowSequence) {
    default:
    case BLOCK_LONG:
      fl = frameLen;
      fr = frameLen -
           getWindow2Nr(frameLen,
                        GetWindowShape(&pAacDecoderChannelInfo->icsInfo));
      if (pAacDecoderStaticChannelInfo->IMdct.prev_tl == 0) {
        fl = fr;
      }
      break;
    case BLOCK_STOP:
      fl = frameLen >> 3;
      fr = frameLen;
      break;
    case BLOCK_START:
      fl = frameLen;
      fr = frameLen >> 3;
      break;
    case BLOCK_SHORT:
      fl = fr = frameLen >> 3;
      tl >>= 3;
      nSpec = 8;
      break;
  }

  {
    int last_frame_lost = pAacDecoderStaticChannelInfo->last_lpc_lost;

    if (pAacDecoderStaticChannelInfo->last_core_mode == LPD) {
      INT fac_FB = 1;
      if (elFlags & AC_EL_FULLBANDLPD) {
        fac_FB = 2;
      }

      FIXP_DBL *synth;

      if (!(elFlags & AC_EL_LPDSTEREOIDX)) {
        synth = pWorkBuffer1 + ((PIT_MAX_MAX - (1 * L_SUBFR)) * fac_FB);
      } else {
        synth = pWorkBuffer1 + PIT_MAX_MAX * fac_FB;
      }

      int fac_length =
          (pAacDecoderChannelInfo->icsInfo.WindowSequence == BLOCK_SHORT)
              ? (frameLen >> 4)
              : (frameLen >> 3);

      INT pitch[NB_SUBFR_SUPERFR + SYN_SFD];
      FIXP_DBL pit_gain[NB_SUBFR_SUPERFR + SYN_SFD];

      int nbDiv = (elFlags & AC_EL_FULLBANDLPD) ? 2 : 4;
      int lFrame = (elFlags & AC_EL_FULLBANDLPD) ? frameLen / 2 : frameLen;
      int nbSubfr =
          lFrame / (nbDiv * L_SUBFR);
      int LpdSfd = (nbDiv * nbSubfr) >> 1;
      int SynSfd = LpdSfd - BPF_SFD;

      FDKmemclear(
          pitch,
          sizeof(
              pitch));
      FDKmemclear(pit_gain, sizeof(pit_gain));

      if (pAacDecoderStaticChannelInfo->last_lpd_mode == 0 ||
          pAacDecoderStaticChannelInfo->last_lpd_mode == 4) {
        FIXP_DBL fac_buf[LFAC];
        FIXP_LPC *A = pAacDecoderChannelInfo->data.usac.lp_coeff[0];

        if (!frameOk || last_frame_lost ||
            (pAacDecoderChannelInfo->data.usac.fac_data[0] == NULL)) {
          FDKmemclear(fac_buf,
                      pAacDecoderChannelInfo->granuleLength * sizeof(FIXP_DBL));
          pAacDecoderChannelInfo->data.usac.fac_data[0] = fac_buf;
          pAacDecoderChannelInfo->data.usac.fac_data_e[0] = 0;
        }

        INT A_exp;
        {
          for (int i = 0; i < M_LP_FILTER_ORDER; i++) {
            A[i] = FX_DBL2FX_LPC(fixp_cos(
                fMult(pAacDecoderStaticChannelInfo->lpc4_lsf[i],
                      FL2FXCONST_SGL((1 << LSPARG_SCALE) * M_PI / 6400.0)),
                LSF_SCALE - LSPARG_SCALE));
          }

          E_LPC_f_lsp_a_conversion(A, A, &A_exp);
        }

#if defined(FDK_ASSERT_ENABLE)
        nSamples =
#endif
            CLpd_FAC_Acelp2Mdct(
                &pAacDecoderStaticChannelInfo->IMdct, synth,
                SPEC_LONG(pAacDecoderChannelInfo->pSpectralCoefficient),
                pAacDecoderChannelInfo->specScale, nSpec,
                pAacDecoderChannelInfo->data.usac.fac_data[0],
                pAacDecoderChannelInfo->data.usac.fac_data_e[0], fac_length,
                frameLen, tl,
                FDKgetWindowSlope(
                    fr, GetWindowShape(&pAacDecoderChannelInfo->icsInfo)),
                fr, A, A_exp, &pAacDecoderStaticChannelInfo->acelp,
                (FIXP_DBL)0,
                (last_frame_lost || !frameOk), 1,
                pAacDecoderStaticChannelInfo->last_lpd_mode, 0,
                pAacDecoderChannelInfo->currAliasingSymmetry);

      } else {
#if defined(FDK_ASSERT_ENABLE)
        nSamples =
#endif
            imlt_block(
                &pAacDecoderStaticChannelInfo->IMdct, synth,
                SPEC_LONG(pAacDecoderChannelInfo->pSpectralCoefficient),
                pAacDecoderChannelInfo->specScale, nSpec, frameLen, tl,
                FDKgetWindowSlope(
                    fl, GetWindowShape(&pAacDecoderChannelInfo->icsInfo)),
                fl,
                FDKgetWindowSlope(
                    fr, GetWindowShape(&pAacDecoderChannelInfo->icsInfo)),
                fr, (FIXP_DBL)0,
                pAacDecoderChannelInfo->currAliasingSymmetry
                    ? MLT_FLAG_CURR_ALIAS_SYMMETRY
                    : 0);
      }
      FDK_ASSERT(nSamples == frameLen);

      if (!(elFlags & AC_EL_LPDSTEREOIDX)) {
        FDKmemcpy(pitch, pAacDecoderStaticChannelInfo->old_T_pf,
                  SynSfd * sizeof(INT));
        FDKmemcpy(pit_gain, pAacDecoderStaticChannelInfo->old_gain_pf,
                  SynSfd * sizeof(FIXP_DBL));

        for (int i = SynSfd; i < LpdSfd + 3; i++) {
          pitch[i] = L_SUBFR;
          pit_gain[i] = (FIXP_DBL)0;
        }

        if (pAacDecoderStaticChannelInfo->last_lpd_mode == 0) {
          pitch[SynSfd] = pitch[SynSfd - 1];
          pit_gain[SynSfd] = pit_gain[SynSfd - 1];
          if (IsLongBlock(&pAacDecoderChannelInfo->icsInfo)) {
            pitch[SynSfd + 1] = pitch[SynSfd];
            pit_gain[SynSfd + 1] = pit_gain[SynSfd];
          }
        }

        {
          FDKmemcpy(
              pWorkBuffer1, pAacDecoderStaticChannelInfo->old_synth,
              ((PIT_MAX_MAX - (1 * L_SUBFR)) * fac_FB) * sizeof(FIXP_DBL));
        }

        FIXP_DBL *p2_synth = pWorkBuffer1 + (PIT_MAX_MAX * fac_FB);

        for (int i = 0; i < SynSfd + 2; i++) {
          int T = pitch[i];
          FIXP_DBL gain = pit_gain[i];

          if (gain > (FIXP_DBL)0) {
            gain = get_gain(&p2_synth[i * L_SUBFR * fac_FB],
                            &p2_synth[(i * L_SUBFR * fac_FB) - fac_FB * T],
                            L_SUBFR * fac_FB);
            pit_gain[i] = gain;
          }
        }

        bass_pf_1sf_delay(p2_synth, pitch, pit_gain, frameLen,
                          (LpdSfd + 2) * L_SUBFR + BPF_SFD * L_SUBFR,
                          frameLen - (LpdSfd + 4) * L_SUBFR, outSamples,
                          aacOutDataHeadroom,
                          pAacDecoderStaticChannelInfo->mem_bpf);
      }

    } else
    {
      FIXP_DBL *tmp =
          pAacDecoderChannelInfo->pComStaticData->pWorkBufferCore1->mdctOutTemp;
#if defined(FDK_ASSERT_ENABLE)
      nSamples =
#endif
          imlt_block(&pAacDecoderStaticChannelInfo->IMdct, tmp,
                     SPEC_LONG(pAacDecoderChannelInfo->pSpectralCoefficient),
                     pAacDecoderChannelInfo->specScale, nSpec, frameLen, tl,
                     FDKgetWindowSlope(
                         fl, GetWindowShape(&pAacDecoderChannelInfo->icsInfo)),
                     fl,
                     FDKgetWindowSlope(
                         fr, GetWindowShape(&pAacDecoderChannelInfo->icsInfo)),
                     fr, (FIXP_DBL)0,
                     pAacDecoderChannelInfo->currAliasingSymmetry
                         ? MLT_FLAG_CURR_ALIAS_SYMMETRY
                         : 0);

      scaleValuesSaturate(outSamples, tmp, frameLen,
                          MDCT_OUT_HEADROOM - aacOutDataHeadroom);
    }
  }

  FDK_ASSERT(nSamples == frameLen);

  pAacDecoderStaticChannelInfo->last_core_mode =
      (pAacDecoderChannelInfo->icsInfo.WindowSequence == BLOCK_SHORT) ? FD_SHORT
                                                                      : FD_LONG;
  pAacDecoderStaticChannelInfo->last_lpd_mode = 255;
}

#include "ldfiltbank.h"
void CBlock_FrequencyToTimeLowDelay(
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo, PCM_DEC outSamples[],
    const short frameLen) {
  InvMdctTransformLowDelay_fdk(
      SPEC_LONG(pAacDecoderChannelInfo->pSpectralCoefficient),
      pAacDecoderChannelInfo->specScale[0], outSamples,
      pAacDecoderStaticChannelInfo->pOverlapBuffer, frameLen);
}
