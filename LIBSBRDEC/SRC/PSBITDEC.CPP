/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "psbitdec.h"
#include "sbr_rom.h"
#include "huff_dec.h"

SBR_ERROR ResetPsDec(HANDLE_PS_DEC h_ps_d);

static SCHAR decode_huff_cw(
    Huffman h,
    HANDLE_FDK_BITSTREAM hBitBuf,
    int *length)
{
  UCHAR bit = 0;
  SCHAR index = 0;
  UCHAR bitCount = 0;

  while (index >= 0) {
    bit = FDKreadBits(hBitBuf, 1);
    bitCount++;
    index = h[index][bit];
  }
  if (length) {
    *length = bitCount;
  }
  return (index + 64);
}

static SCHAR limitMinMax(SCHAR i, SCHAR min, SCHAR max) {
  if (i < min)
    return min;
  else if (i > max)
    return max;
  else
    return i;
}

static void deltaDecodeArray(
    SCHAR enable, SCHAR *aIndex,
    SCHAR *aPrevFrameIndex,
    SCHAR DtDf, UCHAR nrElements,
    UCHAR stride,
    SCHAR minIdx, SCHAR maxIdx) {
  int i;

  if (enable == 1) {
    if (DtDf == 0) {
      aIndex[0] = 0 + aIndex[0];
      aIndex[0] = limitMinMax(aIndex[0], minIdx, maxIdx);
      for (i = 1; i < nrElements; i++) {
        aIndex[i] = aIndex[i - 1] + aIndex[i];
        aIndex[i] = limitMinMax(aIndex[i], minIdx, maxIdx);
      }
    } else {
      for (i = 0; i < nrElements; i++) {
        aIndex[i] = aPrevFrameIndex[i * stride] + aIndex[i];
        aIndex[i] = limitMinMax(aIndex[i], minIdx, maxIdx);
      }
    }
  } else {
    for (i = 0; i < nrElements; i++) {
      aIndex[i] = 0;
    }
  }
  if (stride == 2) {
    for (i = nrElements * stride - 1; i > 0; i--) {
      aIndex[i] = aIndex[i >> 1];
    }
  }
}

static void map34IndexTo20(SCHAR *aIndex,
                           UCHAR noBins)
{
  aIndex[0] = (2 * aIndex[0] + aIndex[1]) / 3;
  aIndex[1] = (aIndex[1] + 2 * aIndex[2]) / 3;
  aIndex[2] = (2 * aIndex[3] + aIndex[4]) / 3;
  aIndex[3] = (aIndex[4] + 2 * aIndex[5]) / 3;
  aIndex[4] = (aIndex[6] + aIndex[7]) / 2;
  aIndex[5] = (aIndex[8] + aIndex[9]) / 2;
  aIndex[6] = aIndex[10];
  aIndex[7] = aIndex[11];
  aIndex[8] = (aIndex[12] + aIndex[13]) / 2;
  aIndex[9] = (aIndex[14] + aIndex[15]) / 2;
  aIndex[10] = aIndex[16];

  if (noBins == NO_HI_RES_BINS) {
    aIndex[11] = aIndex[17];
    aIndex[12] = aIndex[18];
    aIndex[13] = aIndex[19];
    aIndex[14] = (aIndex[20] + aIndex[21]) / 2;
    aIndex[15] = (aIndex[22] + aIndex[23]) / 2;
    aIndex[16] = (aIndex[24] + aIndex[25]) / 2;
    aIndex[17] = (aIndex[26] + aIndex[27]) / 2;
    aIndex[18] = (aIndex[28] + aIndex[29] + aIndex[30] + aIndex[31]) / 4;
    aIndex[19] = (aIndex[32] + aIndex[33]) / 2;
  }
}

int DecodePs(struct PS_DEC *h_ps_d,
             const UCHAR frameError,
             PS_DEC_COEFFICIENTS *pScratch) {
  MPEG_PS_BS_DATA *pBsData;
  UCHAR gr, env;
  int bPsHeaderValid, bPsDataAvail;

  h_ps_d->specificTo.mpeg.pCoef = pScratch;

  pBsData = &h_ps_d->bsData[h_ps_d->processSlot].mpeg;
  bPsHeaderValid = pBsData->bPsHeaderValid;
  bPsDataAvail =
      (h_ps_d->bPsDataAvail[h_ps_d->processSlot] == ppt_mpeg) ? 1 : 0;

  if ((h_ps_d->psDecodedPrv && !frameError && !bPsDataAvail) ||
      (!h_ps_d->psDecodedPrv &&
       (frameError || !bPsDataAvail || !bPsHeaderValid))) {
    pBsData->bPsHeaderValid = 0;
    h_ps_d->bPsDataAvail[h_ps_d->processSlot] = ppt_none;
    return (0);
  }

  if (frameError ||
      !bPsHeaderValid) {
    pBsData->noEnv = 0;
  }

  for (env = 0; env < pBsData->noEnv; env++) {
    SCHAR *aPrevIidIndex;
    SCHAR *aPrevIccIndex;

    UCHAR noIidSteps = pBsData->bFineIidQ ? NO_IID_STEPS_FINE : NO_IID_STEPS;

    if (env == 0) {
      aPrevIidIndex = h_ps_d->specificTo.mpeg.aIidPrevFrameIndex;
      aPrevIccIndex = h_ps_d->specificTo.mpeg.aIccPrevFrameIndex;
    } else {
      aPrevIidIndex = pBsData->aaIidIndex[env - 1];
      aPrevIccIndex = pBsData->aaIccIndex[env - 1];
    }

    deltaDecodeArray(pBsData->bEnableIid, pBsData->aaIidIndex[env],
                     aPrevIidIndex, pBsData->abIidDtFlag[env],
                     FDK_sbrDecoder_aNoIidBins[pBsData->freqResIid],
                     (pBsData->freqResIid) ? 1 : 2, -noIidSteps, noIidSteps);

    deltaDecodeArray(pBsData->bEnableIcc, pBsData->aaIccIndex[env],
                     aPrevIccIndex, pBsData->abIccDtFlag[env],
                     FDK_sbrDecoder_aNoIccBins[pBsData->freqResIcc],
                     (pBsData->freqResIcc) ? 1 : 2, 0, NO_ICC_STEPS - 1);
  }

  if (pBsData->noEnv == 0) {
    pBsData->noEnv = 1;

    if (pBsData->bEnableIid) {
      pBsData->bFineIidQ = h_ps_d->specificTo.mpeg.bPrevFrameFineIidQ;
      pBsData->freqResIid = h_ps_d->specificTo.mpeg.prevFreqResIid;
      for (gr = 0; gr < NO_HI_RES_IID_BINS; gr++) {
        pBsData->aaIidIndex[pBsData->noEnv - 1][gr] =
            h_ps_d->specificTo.mpeg.aIidPrevFrameIndex[gr];
      }
    } else {
      for (gr = 0; gr < NO_HI_RES_IID_BINS; gr++) {
        pBsData->aaIidIndex[pBsData->noEnv - 1][gr] = 0;
      }
    }

    if (pBsData->bEnableIcc) {
      pBsData->freqResIcc = h_ps_d->specificTo.mpeg.prevFreqResIcc;
      for (gr = 0; gr < NO_HI_RES_ICC_BINS; gr++) {
        pBsData->aaIccIndex[pBsData->noEnv - 1][gr] =
            h_ps_d->specificTo.mpeg.aIccPrevFrameIndex[gr];
      }
    } else {
      for (gr = 0; gr < NO_HI_RES_ICC_BINS; gr++) {
        pBsData->aaIccIndex[pBsData->noEnv - 1][gr] = 0;
      }
    }
  }

  h_ps_d->specificTo.mpeg.bPrevFrameFineIidQ = pBsData->bFineIidQ;

  h_ps_d->specificTo.mpeg.prevFreqResIid = pBsData->freqResIid;

  h_ps_d->specificTo.mpeg.prevFreqResIcc = pBsData->freqResIcc;

  for (gr = 0; gr < NO_HI_RES_IID_BINS; gr++) {
    h_ps_d->specificTo.mpeg.aIidPrevFrameIndex[gr] =
        pBsData->aaIidIndex[pBsData->noEnv - 1][gr];
  }
  for (gr = 0; gr < NO_HI_RES_ICC_BINS; gr++) {
    h_ps_d->specificTo.mpeg.aIccPrevFrameIndex[gr] =
        pBsData->aaIccIndex[pBsData->noEnv - 1][gr];
  }

  h_ps_d->bPsDataAvail[h_ps_d->processSlot] = ppt_none;

  if (pBsData->bFrameClass == 0) {
    pBsData->aEnvStartStop[0] = 0;
    for (env = 1; env < pBsData->noEnv; env++) {
      pBsData->aEnvStartStop[env] =
          (env * h_ps_d->noSubSamples) / pBsData->noEnv;
    }
    pBsData->aEnvStartStop[pBsData->noEnv] = h_ps_d->noSubSamples;
  } else {
    pBsData->aEnvStartStop[0] = 0;

    if (pBsData->aEnvStartStop[pBsData->noEnv] < h_ps_d->noSubSamples) {
      for (gr = 0; gr < NO_HI_RES_IID_BINS; gr++) {
        pBsData->aaIidIndex[pBsData->noEnv][gr] =
            pBsData->aaIidIndex[pBsData->noEnv - 1][gr];
      }
      for (gr = 0; gr < NO_HI_RES_ICC_BINS; gr++) {
        pBsData->aaIccIndex[pBsData->noEnv][gr] =
            pBsData->aaIccIndex[pBsData->noEnv - 1][gr];
      }
      pBsData->noEnv++;
      pBsData->aEnvStartStop[pBsData->noEnv] = h_ps_d->noSubSamples;
    }

    for (env = 1; env < pBsData->noEnv; env++) {
      UCHAR thr;
      thr = (UCHAR)h_ps_d->noSubSamples - (pBsData->noEnv - env);
      if (pBsData->aEnvStartStop[env] > thr) {
        pBsData->aEnvStartStop[env] = thr;
      } else {
        thr = pBsData->aEnvStartStop[env - 1] + 1;
        if (pBsData->aEnvStartStop[env] < thr) {
          pBsData->aEnvStartStop[env] = thr;
        }
      }
    }
  }

  for (env = 0; env < pBsData->noEnv; env++) {
    UCHAR i;
    for (i = 0; i < NO_HI_RES_IID_BINS; i++) {
      h_ps_d->specificTo.mpeg.pCoef->aaIidIndexMapped[env][i] =
          pBsData->aaIidIndex[env][i];
    }
    for (i = 0; i < NO_HI_RES_ICC_BINS; i++) {
      h_ps_d->specificTo.mpeg.pCoef->aaIccIndexMapped[env][i] =
          pBsData->aaIccIndex[env][i];
    }
  }

  for (env = 0; env < pBsData->noEnv; env++) {
    if (pBsData->freqResIid == 2)
      map34IndexTo20(h_ps_d->specificTo.mpeg.pCoef->aaIidIndexMapped[env],
                     NO_HI_RES_IID_BINS);
    if (pBsData->freqResIcc == 2)
      map34IndexTo20(h_ps_d->specificTo.mpeg.pCoef->aaIccIndexMapped[env],
                     NO_HI_RES_ICC_BINS);
  }

  return (1);
}

unsigned int ReadPsData(
    HANDLE_PS_DEC h_ps_d,
    HANDLE_FDK_BITSTREAM hBitBuf,
    int nBitsLeft
) {
  MPEG_PS_BS_DATA *pBsData;

  UCHAR gr, env;
  SCHAR dtFlag;
  INT startbits;
  Huffman CurrentTable;
  SCHAR bEnableHeader;

  if (!h_ps_d) return 0;

  pBsData = &h_ps_d->bsData[h_ps_d->bsReadSlot].mpeg;

  if (h_ps_d->bsReadSlot != h_ps_d->bsLastSlot) {
    FDKmemcpy(pBsData, &h_ps_d->bsData[h_ps_d->bsLastSlot].mpeg,
              sizeof(MPEG_PS_BS_DATA));
  }

  startbits = (INT)FDKgetValidBits(hBitBuf);

  bEnableHeader = (SCHAR)FDKreadBits(hBitBuf, 1);

  if (bEnableHeader) {
    pBsData->bPsHeaderValid = 1;
    pBsData->bEnableIid = (UCHAR)FDKreadBits(hBitBuf, 1);
    if (pBsData->bEnableIid) {
      pBsData->modeIid = (UCHAR)FDKreadBits(hBitBuf, 3);
    }

    pBsData->bEnableIcc = (UCHAR)FDKreadBits(hBitBuf, 1);
    if (pBsData->bEnableIcc) {
      pBsData->modeIcc = (UCHAR)FDKreadBits(hBitBuf, 3);
    }

    pBsData->bEnableExt = (UCHAR)FDKreadBits(hBitBuf, 1);
  }

  pBsData->bFrameClass = (UCHAR)FDKreadBits(hBitBuf, 1);
  if (pBsData->bFrameClass == 0) {
    pBsData->noEnv =
        FDK_sbrDecoder_aFixNoEnvDecode[(UCHAR)FDKreadBits(hBitBuf, 2)];
  } else {
    pBsData->noEnv = 1 + (UCHAR)FDKreadBits(hBitBuf, 2);
    for (env = 1; env < pBsData->noEnv + 1; env++)
      pBsData->aEnvStartStop[env] = ((UCHAR)FDKreadBits(hBitBuf, 5)) + 1;
  }

  if ((pBsData->modeIid > 5) || (pBsData->modeIcc > 5)) {
    h_ps_d->bPsDataAvail[h_ps_d->bsReadSlot] = ppt_none;
    nBitsLeft -= startbits - (INT)FDKgetValidBits(hBitBuf);
    while (nBitsLeft > 0) {
      int i = nBitsLeft;
      if (i > 8) {
        i = 8;
      }
      FDKreadBits(hBitBuf, i);
      nBitsLeft -= i;
    }
    return (UINT)(startbits - (INT)FDKgetValidBits(hBitBuf));
  }

  if (pBsData->modeIid > 2) {
    pBsData->freqResIid = pBsData->modeIid - 3;
    pBsData->bFineIidQ = 1;
  } else {
    pBsData->freqResIid = pBsData->modeIid;
    pBsData->bFineIidQ = 0;
  }

  if (pBsData->modeIcc > 2) {
    pBsData->freqResIcc = pBsData->modeIcc - 3;
  } else {
    pBsData->freqResIcc = pBsData->modeIcc;
  }

  if (pBsData->bEnableIid) {
    for (env = 0; env < pBsData->noEnv; env++) {
      dtFlag = (SCHAR)FDKreadBits(hBitBuf, 1);
      if (!dtFlag) {
        if (pBsData->bFineIidQ)
          CurrentTable = (Huffman)&aBookPsIidFineFreqDecode;
        else
          CurrentTable = (Huffman)&aBookPsIidFreqDecode;
      } else {
        if (pBsData->bFineIidQ)
          CurrentTable = (Huffman)&aBookPsIidFineTimeDecode;
        else
          CurrentTable = (Huffman)&aBookPsIidTimeDecode;
      }

      for (gr = 0; gr < FDK_sbrDecoder_aNoIidBins[pBsData->freqResIid]; gr++)
        pBsData->aaIidIndex[env][gr] =
            decode_huff_cw(CurrentTable, hBitBuf, NULL);
      pBsData->abIidDtFlag[env] = dtFlag;
    }
  }

  if (pBsData->bEnableIcc) {
    for (env = 0; env < pBsData->noEnv; env++) {
      dtFlag = (SCHAR)FDKreadBits(hBitBuf, 1);
      if (!dtFlag)
        CurrentTable = (Huffman)&aBookPsIccFreqDecode;
      else
        CurrentTable = (Huffman)&aBookPsIccTimeDecode;

      for (gr = 0; gr < FDK_sbrDecoder_aNoIccBins[pBsData->freqResIcc]; gr++)
        pBsData->aaIccIndex[env][gr] =
            decode_huff_cw(CurrentTable, hBitBuf, NULL);
      pBsData->abIccDtFlag[env] = dtFlag;
    }
  }

  if (pBsData->bEnableExt) {
    int cnt = FDKreadBits(hBitBuf, PS_EXTENSION_SIZE_BITS);
    if (cnt == (1 << PS_EXTENSION_SIZE_BITS) - 1) {
      cnt += FDKreadBits(hBitBuf, PS_EXTENSION_ESC_COUNT_BITS);
    }
    while (cnt--) FDKreadBits(hBitBuf, 8);
  }

  h_ps_d->bPsDataAvail[h_ps_d->bsReadSlot] = ppt_mpeg;

  return (startbits - (INT)FDKgetValidBits(hBitBuf));
}
