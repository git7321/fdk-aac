/******************* MPEG transport format decoder library *********************

   Author(s):   Josef Hoepfl

   Description: ADTS interface

*******************************************************************************/

#include "tpdec_adts.h"
#include "FDK_bitstream.h"

void adtsRead_CrcInit(
    HANDLE_ADTS pAdts)
{
  FDKcrcInit(&pAdts->crcInfo, 0x8005, 0xFFFF, 16);
}

int adtsRead_CrcStartReg(
    HANDLE_ADTS pAdts,
    HANDLE_FDK_BITSTREAM hBs,
    int mBits
) {
  if (pAdts->bs.protection_absent) {
    return 0;
  }

  return (FDKcrcStartReg(&pAdts->crcInfo, hBs, mBits));
}

void adtsRead_CrcEndReg(
    HANDLE_ADTS pAdts,
    HANDLE_FDK_BITSTREAM hBs,
    int reg
) {
  if (pAdts->bs.protection_absent == 0) {
    FDKcrcEndReg(&pAdts->crcInfo, hBs, reg);
  }
}

TRANSPORTDEC_ERROR adtsRead_CrcCheck(HANDLE_ADTS pAdts) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;
  USHORT crc;

  if (pAdts->bs.protection_absent) return TRANSPORTDEC_OK;

  crc = FDKcrcGetCRC(&pAdts->crcInfo);
  if (crc != pAdts->crcReadValue) {
    return (TRANSPORTDEC_CRC_ERROR);
  }

  return (ErrorStatus);
}

#define Adts_Length_SyncWord 12
#define Adts_Length_Id 1
#define Adts_Length_Layer 2
#define Adts_Length_ProtectionAbsent 1
#define Adts_Length_Profile 2
#define Adts_Length_SamplingFrequencyIndex 4
#define Adts_Length_PrivateBit 1
#define Adts_Length_ChannelConfiguration 3
#define Adts_Length_OriginalCopy 1
#define Adts_Length_Home 1
#define Adts_Length_CopyrightIdentificationBit 1
#define Adts_Length_CopyrightIdentificationStart 1
#define Adts_Length_FrameLength 13
#define Adts_Length_BufferFullness 11
#define Adts_Length_NumberOfRawDataBlocksInFrame 2
#define Adts_Length_CrcCheck 16

TRANSPORTDEC_ERROR adtsRead_DecodeHeader(HANDLE_ADTS pAdts,
                                         CSAudioSpecificConfig *pAsc,
                                         HANDLE_FDK_BITSTREAM hBs,
                                         const INT ignoreBufferFullness) {
  INT crcReg;

  INT valBits;
  INT cmp_buffer_fullness;
  int i, adtsHeaderLength;

  STRUCT_ADTS_BS bs;

  CProgramConfig oldPce;
  FDKmemcpy(&oldPce, &pAsc->m_progrConfigElement, sizeof(CProgramConfig));

  valBits = FDKgetValidBits(hBs) + ADTS_SYNCLENGTH;

  if (valBits < ADTS_HEADERLENGTH) {
    return TRANSPORTDEC_NOT_ENOUGH_BITS;
  }

  bs.mpeg_id = FDKreadBits(hBs, Adts_Length_Id);
  bs.layer = FDKreadBits(hBs, Adts_Length_Layer);
  bs.protection_absent = FDKreadBits(hBs, Adts_Length_ProtectionAbsent);
  bs.profile = FDKreadBits(hBs, Adts_Length_Profile);
  bs.sample_freq_index = FDKreadBits(hBs, Adts_Length_SamplingFrequencyIndex);
  bs.private_bit = FDKreadBits(hBs, Adts_Length_PrivateBit);
  bs.channel_config = FDKreadBits(hBs, Adts_Length_ChannelConfiguration);
  bs.original = FDKreadBits(hBs, Adts_Length_OriginalCopy);
  bs.home = FDKreadBits(hBs, Adts_Length_Home);

  bs.copyright_id = FDKreadBits(hBs, Adts_Length_CopyrightIdentificationBit);
  bs.copyright_start =
      FDKreadBits(hBs, Adts_Length_CopyrightIdentificationStart);
  bs.frame_length = FDKreadBits(hBs, Adts_Length_FrameLength);
  bs.adts_fullness = FDKreadBits(hBs, Adts_Length_BufferFullness);
  bs.num_raw_blocks =
      FDKreadBits(hBs, Adts_Length_NumberOfRawDataBlocksInFrame);
  bs.num_pce_bits = 0;

  adtsHeaderLength = ADTS_HEADERLENGTH;

  if (valBits < bs.frame_length * 8) {
    goto bail;
  }

  FDKcrcReset(&pAdts->crcInfo);
  if (!bs.protection_absent) {
    FDKpushBack(hBs, 56);
    crcReg = FDKcrcStartReg(&pAdts->crcInfo, hBs, 0);
    FDKpushFor(hBs, 56);
  }

  if (!bs.protection_absent && bs.num_raw_blocks > 0) {
    if ((INT)FDKgetValidBits(hBs) < bs.num_raw_blocks * 16) {
      goto bail;
    }
    for (i = 0; i < bs.num_raw_blocks; i++) {
      pAdts->rawDataBlockDist[i] = (USHORT)FDKreadBits(hBs, 16);
      adtsHeaderLength += 16;
    }
    pAdts->rawDataBlockDist[bs.num_raw_blocks] =
        bs.frame_length - 7 - bs.num_raw_blocks * 2 - 2;
    for (i = bs.num_raw_blocks; i > 0; i--) {
      pAdts->rawDataBlockDist[i] -= pAdts->rawDataBlockDist[i - 1];
    }
  }

  if (!bs.protection_absent) {
    USHORT crc_check;

    FDKcrcEndReg(&pAdts->crcInfo, hBs, crcReg);

    if ((INT)FDKgetValidBits(hBs) < Adts_Length_CrcCheck) {
      goto bail;
    }

    crc_check = FDKreadBits(hBs, Adts_Length_CrcCheck);
    adtsHeaderLength += Adts_Length_CrcCheck;

    pAdts->crcReadValue = crc_check;
    if (bs.num_raw_blocks > 0) {
      if (pAdts->crcReadValue != FDKcrcGetCRC(&pAdts->crcInfo)) {
        return TRANSPORTDEC_CRC_ERROR;
      }
      FDKcrcReset(&pAdts->crcInfo);
    }
  }

  if ((bs.layer != 0) ||
      (bs.sample_freq_index >= 13)
  ) {
    FDKpushFor(hBs, bs.frame_length * 8);
    return TRANSPORTDEC_UNSUPPORTED_FORMAT;
  }

  if ((bs.mpeg_id == 0) && (pAdts->decoderCanDoMpeg4 == 0)) {

    FDKpushFor(hBs, bs.frame_length * 8);
    return TRANSPORTDEC_UNSUPPORTED_FORMAT;
  }

  if (!ignoreBufferFullness) {
    cmp_buffer_fullness =
        bs.frame_length * 8 +
        bs.adts_fullness * 32 * getNumberOfEffectiveChannels(bs.channel_config);

    if (bs.adts_fullness != 0x7FF) {
      if (pAdts->BufferFullnesStartFlag) {
        if (valBits < cmp_buffer_fullness) {
          FDKpushBack(hBs, adtsHeaderLength);

          if ((cmp_buffer_fullness + adtsHeaderLength) >
              (((8192 * 4) << 3) - 7)) {
            return TRANSPORTDEC_SYNC_ERROR;
          } else {
            return TRANSPORTDEC_NOT_ENOUGH_BITS;
          }
        } else {
          pAdts->BufferFullnesStartFlag = 0;
        }
      }
    }
  }

  AudioSpecificConfig_Init(pAsc);
  pAsc->m_aot = (AUDIO_OBJECT_TYPE)(bs.profile + 1);
  pAsc->m_samplingFrequencyIndex = bs.sample_freq_index;
  pAsc->m_samplingFrequency = SamplingRateTable[bs.sample_freq_index];
  pAsc->m_channelConfiguration = bs.channel_config;
  pAsc->m_samplesPerFrame = 1024;

  if (bs.channel_config == 0) {
    int pceBits = 0;
    UINT alignAnchor = FDKgetValidBits(hBs);
    CProgramConfig tmpPce;

    if (FDKreadBits(hBs, 3) == ID_PCE) {
      crcReg = adtsRead_CrcStartReg(pAdts, hBs, 0);

      CProgramConfig_Init(&tmpPce);
      CProgramConfig_Read(&tmpPce, hBs, alignAnchor);

      if (CProgramConfig_IsValid(&tmpPce)) {
        if (CProgramConfig_IsValid(&oldPce)) {
          switch (CProgramConfig_Compare(&tmpPce, &oldPce)) {
            case 0:
            case 1:
              FDKmemcpy(&pAsc->m_progrConfigElement, &tmpPce,
                        sizeof(CProgramConfig));
              break;
            case 2:
            case -1:
            default:
              FDKmemcpy(&pAsc->m_progrConfigElement, &oldPce,
                        sizeof(CProgramConfig));
              FDKpushBack(hBs, adtsHeaderLength);
              return TRANSPORTDEC_PARSE_ERROR;
          }
        } else {
          FDKmemcpy(&pAsc->m_progrConfigElement, &tmpPce,
                    sizeof(CProgramConfig));
        }
      } else {
        if (CProgramConfig_IsValid(&oldPce)) {
          FDKmemcpy(&pAsc->m_progrConfigElement, &oldPce,
                    sizeof(CProgramConfig));
        } else {
          FDKpushBack(hBs, adtsHeaderLength);
          return TRANSPORTDEC_PARSE_ERROR;
        }
      }

      adtsRead_CrcEndReg(pAdts, hBs, crcReg);
      pceBits = (INT)alignAnchor - (INT)FDKgetValidBits(hBs);
      adtsHeaderLength += pceBits;

      if (pceBits > (INT)alignAnchor) {
        goto bail;
      }

      bs.num_pce_bits = pceBits;
    } else {
      FDKpushBack(hBs, 3);

      if (oldPce.isValid &&
          (bs.sample_freq_index ==
           pAdts->bs.sample_freq_index)
          && (bs.channel_config == pAdts->bs.channel_config)
          &&
          (bs.mpeg_id ==
           pAdts->bs.mpeg_id)) {
        FDKmemcpy(&pAsc->m_progrConfigElement, &oldPce, sizeof(CProgramConfig));
      } else if (bs.mpeg_id == 0) {
        FDKpushFor(hBs, (bs.frame_length << 3) - adtsHeaderLength - 3);
        return TRANSPORTDEC_UNSUPPORTED_FORMAT;
      }
    }
  }

  FDKmemcpy(&pAdts->bs, &bs, sizeof(STRUCT_ADTS_BS));

  return TRANSPORTDEC_OK;

bail:
  FDKpushBack(hBs, adtsHeaderLength);
  return TRANSPORTDEC_NOT_ENOUGH_BITS;
}

int adtsRead_GetRawDataBlockLength(HANDLE_ADTS pAdts, INT blockNum) {
  int length;

  if (pAdts->bs.num_raw_blocks == 0) {
    length =
        (pAdts->bs.frame_length - 7)
        << 3;
    if (pAdts->bs.protection_absent == 0)
      length -= 16;
  } else {
    if (pAdts->bs.protection_absent) {
      length = -1;
    } else {
      if (blockNum < 0 || blockNum > 3) {
        length = -1;
      } else {
        length = (pAdts->rawDataBlockDist[blockNum] << 3) - 16;
      }
    }
  }
  if (blockNum == 0 && length > 0) {
    length -= pAdts->bs.num_pce_bits;
  }
  return length;
}
