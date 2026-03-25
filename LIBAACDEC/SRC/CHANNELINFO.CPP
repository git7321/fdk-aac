/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: individual channel stream info

*******************************************************************************/

#include "channelinfo.h"
#include "aac_rom.h"
#include "aac_ram.h"
#include "FDK_bitstream.h"

AAC_DECODER_ERROR IcsReadMaxSfb(HANDLE_FDK_BITSTREAM bs, CIcsInfo *pIcsInfo,
                                const SamplingRateInfo *pSamplingRateInfo) {
  AAC_DECODER_ERROR ErrorStatus = AAC_DEC_OK;
  int nbits;

  if (IsLongBlock(pIcsInfo)) {
    nbits = 6;
    pIcsInfo->TotalSfBands = pSamplingRateInfo->NumberOfScaleFactorBands_Long;
  } else {
    nbits = 4;
    pIcsInfo->TotalSfBands = pSamplingRateInfo->NumberOfScaleFactorBands_Short;
  }
  pIcsInfo->MaxSfBands = (UCHAR)FDKreadBits(bs, nbits);

  if (pIcsInfo->MaxSfBands > pIcsInfo->TotalSfBands) {
    ErrorStatus = AAC_DEC_PARSE_ERROR;
  }

  return ErrorStatus;
}

AAC_DECODER_ERROR IcsRead(HANDLE_FDK_BITSTREAM bs, CIcsInfo *pIcsInfo,
                          const SamplingRateInfo *pSamplingRateInfo,
                          const UINT flags) {
  AAC_DECODER_ERROR ErrorStatus = AAC_DEC_OK;

  pIcsInfo->Valid = 0;

  if (flags & AC_ELD) {
    pIcsInfo->WindowSequence = BLOCK_LONG;
    pIcsInfo->WindowShape = 0;
  } else {
    if (!(flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA))) {
      FDKreadBits(bs, 1);
    }
    pIcsInfo->WindowSequence = (BLOCK_TYPE)FDKreadBits(bs, 2);
    pIcsInfo->WindowShape = (UCHAR)FDKreadBits(bs, 1);
    if (flags & AC_LD) {
      if (pIcsInfo->WindowShape) {
        pIcsInfo->WindowShape = 2;
      }
    }
  }

  if ((flags & (AC_ELD | AC_LD)) && pIcsInfo->WindowSequence != BLOCK_LONG) {
    pIcsInfo->WindowSequence = BLOCK_LONG;
    ErrorStatus = AAC_DEC_PARSE_ERROR;
    goto bail;
  }

  ErrorStatus = IcsReadMaxSfb(bs, pIcsInfo, pSamplingRateInfo);
  if (ErrorStatus != AAC_DEC_OK) {
    goto bail;
  }

  if (IsLongBlock(pIcsInfo)) {
    if (!(flags & (AC_ELD | AC_SCALABLE | AC_BSAC | AC_USAC | AC_RSVD50 |
                   AC_RSV603DA)))
    {
      if ((UCHAR)FDKreadBits(bs, 1) != 0)
      {
        ErrorStatus = AAC_DEC_UNSUPPORTED_PREDICTION;
        goto bail;
      }
    }

    pIcsInfo->WindowGroups = 1;
    pIcsInfo->WindowGroupLength[0] = 1;
  } else {
    INT i;
    UINT mask;

    pIcsInfo->ScaleFactorGrouping = (UCHAR)FDKreadBits(bs, 7);

    pIcsInfo->WindowGroups = 0;

    for (i = 0; i < (8 - 1); i++) {
      mask = 1 << (6 - i);
      pIcsInfo->WindowGroupLength[i] = 1;

      if (pIcsInfo->ScaleFactorGrouping & mask) {
        pIcsInfo->WindowGroupLength[pIcsInfo->WindowGroups]++;
      } else {
        pIcsInfo->WindowGroups++;
      }
    }

    pIcsInfo->WindowGroupLength[8 - 1] = 1;
    pIcsInfo->WindowGroups++;
  }

bail:
  if (ErrorStatus == AAC_DEC_OK) pIcsInfo->Valid = 1;

  return ErrorStatus;
}

AAC_DECODER_ERROR getSamplingRateInfo(SamplingRateInfo *t, UINT samplesPerFrame,
                                      UINT samplingRateIndex,
                                      UINT samplingRate) {
  int index = 0;

  if ((samplingRateIndex >= 15) || (samplesPerFrame == 768)) {
    const UINT borders[] = {(UINT)-1, 92017, 75132, 55426, 46009, 37566,
                            27713,    23004, 18783, 13856, 11502, 9391};
    UINT i, samplingRateSearch = samplingRate;

    if (samplesPerFrame == 768) {
      samplingRateSearch = (samplingRate * 4) / 3;
    }

    for (i = 0; i < 11; i++) {
      if (borders[i] > samplingRateSearch &&
          samplingRateSearch >= borders[i + 1]) {
        break;
      }
    }
    samplingRateIndex = i;
  }

  t->samplingRateIndex = samplingRateIndex;
  t->samplingRate = samplingRate;

  switch (samplesPerFrame) {
    case 1024:
      index = 0;
      break;
    case 960:
      index = 1;
      break;
    case 768:
      index = 2;
      break;
    case 512:
      index = 3;
      break;
    case 480:
      index = 4;
      break;

    default:
      return AAC_DEC_UNSUPPORTED_FORMAT;
  }

  t->ScaleFactorBands_Long =
      sfbOffsetTables[index][samplingRateIndex].sfbOffsetLong;
  t->ScaleFactorBands_Short =
      sfbOffsetTables[index][samplingRateIndex].sfbOffsetShort;
  t->NumberOfScaleFactorBands_Long =
      sfbOffsetTables[index][samplingRateIndex].numberOfSfbLong;
  t->NumberOfScaleFactorBands_Short =
      sfbOffsetTables[index][samplingRateIndex].numberOfSfbShort;

  if (t->ScaleFactorBands_Long == NULL ||
      t->NumberOfScaleFactorBands_Long == 0) {
    t->samplingRate = 0;
    return AAC_DEC_UNSUPPORTED_FORMAT;
  }

  FDK_ASSERT((UINT)t->ScaleFactorBands_Long[t->NumberOfScaleFactorBands_Long] ==
             samplesPerFrame);
  FDK_ASSERT(
      t->ScaleFactorBands_Short == NULL ||
      (UINT)t->ScaleFactorBands_Short[t->NumberOfScaleFactorBands_Short] * 8 ==
          samplesPerFrame);

  return AAC_DEC_OK;
}
