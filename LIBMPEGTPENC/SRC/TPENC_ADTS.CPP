/******************* MPEG transport format encoder library *********************

   Author(s):   Alex Groeschel

   Description: ADTS Transport Headers support

*******************************************************************************/

#include "tpenc_adts.h"
#include "tpenc_lib.h"
#include "tpenc_asc.h"

int adtsWrite_CrcStartReg(
    HANDLE_ADTS pAdts,
    HANDLE_FDK_BITSTREAM hBs,
    int mBits
) {
  if (pAdts->protection_absent) {
    return 0;
  }
  return (FDKcrcStartReg(&pAdts->crcInfo, hBs, mBits));
}

void adtsWrite_CrcEndReg(
    HANDLE_ADTS pAdts,
    HANDLE_FDK_BITSTREAM hBs,
    int reg
) {
  if (pAdts->protection_absent == 0) {
    FDKcrcEndReg(&pAdts->crcInfo, hBs, reg);
  }
}

int adtsWrite_GetHeaderBits(HANDLE_ADTS hAdts) {
  int bits = 0;

  if (hAdts->currentBlock == 0) {
    bits = 56;
    if (!hAdts->protection_absent) {
      bits += 16;
      if (hAdts->num_raw_blocks > 0) {
        bits += (hAdts->num_raw_blocks) * 16;
      }
    }
  }
  if (!hAdts->protection_absent && hAdts->num_raw_blocks > 0) {
    bits += 16;
  }

  hAdts->headerBits = bits;

  return bits;
}

INT adtsWrite_Init(HANDLE_ADTS hAdts, CODER_CONFIG *config) {
  if (config->nSubFrames < 1 || config->nSubFrames > 4 ||
      (int)config->aot > 4 || (int)config->aot < 1) {
    return -1;
  }

  if (config->flags & CC_MPEG_ID) {
    hAdts->mpeg_id = 0;
  } else {
    hAdts->mpeg_id = 1;
  }
  hAdts->layer = 0;
  hAdts->protection_absent = !(config->flags & CC_PROTECTION);
  hAdts->profile = ((int)config->aot) - 1;
  hAdts->sample_freq_index = getSamplingRateIndex(config->samplingRate, 4);
  hAdts->sample_freq = config->samplingRate;
  hAdts->private_bit = 0;
  hAdts->channel_mode = config->channelMode;
  hAdts->original = 0;
  hAdts->home = 0;
  hAdts->copyright_id = 0;
  hAdts->copyright_start = 0;

  hAdts->num_raw_blocks = config->nSubFrames - 1;

  hAdts->channel_config_zero = config->channelConfigZero;

  FDKcrcInit(&hAdts->crcInfo, 0x8005, 0xFFFF, 16);

  hAdts->currentBlock = 0;

  return 0;
}

int adtsWrite_EncodeHeader(HANDLE_ADTS hAdts, HANDLE_FDK_BITSTREAM hBitStream,
                           int buffer_fullness, int frame_length) {
  INT crcIndex = 0;

  hAdts->headerBits = adtsWrite_GetHeaderBits(hAdts);

  FDK_ASSERT(((frame_length + hAdts->headerBits) / 8) < 0x2000);
  FDK_ASSERT(buffer_fullness < 0x800);

  if (!hAdts->protection_absent) {
    FDKcrcReset(&hAdts->crcInfo);
  }

  if (hAdts->currentBlock == 0) {
    FDKresetBitbuffer(hBitStream, BS_WRITER);
  }

  hAdts->subFrameStartBit = FDKgetValidBits(hBitStream);

  if (hAdts->currentBlock == 0) {
    FDKresetBitbuffer(hBitStream, BS_WRITER);

    if (hAdts->num_raw_blocks == 0) {
      crcIndex = adtsWrite_CrcStartReg(hAdts, hBitStream, 0);
    }

    FDKwriteBits(hBitStream, 0xFFF, 12);
    FDKwriteBits(hBitStream, hAdts->mpeg_id, 1);
    FDKwriteBits(hBitStream, hAdts->layer, 2);
    FDKwriteBits(hBitStream, hAdts->protection_absent, 1);
    FDKwriteBits(hBitStream, hAdts->profile, 2);
    FDKwriteBits(hBitStream, hAdts->sample_freq_index, 4);
    FDKwriteBits(hBitStream, hAdts->private_bit, 1);
    FDKwriteBits(
        hBitStream,
        getChannelConfig(hAdts->channel_mode, hAdts->channel_config_zero), 3);
    FDKwriteBits(hBitStream, hAdts->original, 1);
    FDKwriteBits(hBitStream, hAdts->home, 1);
    FDKwriteBits(hBitStream, hAdts->copyright_id, 1);
    FDKwriteBits(hBitStream, hAdts->copyright_start, 1);
    FDKwriteBits(hBitStream, (frame_length + hAdts->headerBits) >> 3, 13);
    FDKwriteBits(hBitStream, buffer_fullness, 11);
    FDKwriteBits(hBitStream, hAdts->num_raw_blocks, 2);

    if (!hAdts->protection_absent) {
      int i;

      if (hAdts->num_raw_blocks == 0) {
        adtsWrite_CrcEndReg(hAdts, hBitStream, crcIndex);
      } else {
        for (i = 0; i < hAdts->num_raw_blocks; i++) {
          FDKwriteBits(hBitStream, 0, 16);
        }
      }
      FDKwriteBits(hBitStream, 0, 16);
    }
  }

  return 0;
}

void adtsWrite_EndRawDataBlock(HANDLE_ADTS hAdts, HANDLE_FDK_BITSTREAM hBs,
                               int *pBits) {
  if (!hAdts->protection_absent) {
    FDK_BITSTREAM bsWriter;

    FDKinitBitStream(&bsWriter, hBs->hBitBuf.Buffer, hBs->hBitBuf.bufSize, 0,
                     BS_WRITER);
    FDKpushFor(&bsWriter, 56);

    if (hAdts->num_raw_blocks == 0) {
      FDKwriteBits(&bsWriter, FDKcrcGetCRC(&hAdts->crcInfo), 16);
    } else {
      int distance;

      FDKwriteBits(hBs, FDKcrcGetCRC(&hAdts->crcInfo), 16);

      if (hAdts->currentBlock < hAdts->num_raw_blocks) {
        FDKpushFor(&bsWriter, hAdts->currentBlock * 16);
        distance =
            FDKgetValidBits(hBs) - (56 + (hAdts->num_raw_blocks) * 16 + 16);
        FDKwriteBits(&bsWriter, distance >> 3, 16);
      }
    }
    FDKsyncCache(&bsWriter);
  }

  if (hAdts->num_raw_blocks > 0 &&
      hAdts->currentBlock == hAdts->num_raw_blocks) {
    FDK_BITSTREAM bsWriter;
    int crcIndex = 0;

    FDKinitBitStream(&bsWriter, hBs->hBitBuf.Buffer, hBs->hBitBuf.bufSize, 0,
                     BS_WRITER);

    if (!hAdts->protection_absent) {
      FDKcrcReset(&hAdts->crcInfo);
      crcIndex = FDKcrcStartReg(&hAdts->crcInfo, &bsWriter, 0);
    }
    FDKpushFor(&bsWriter, 56 - 28 + 2);
    FDKwriteBits(&bsWriter, FDKgetValidBits(hBs) >> 3, 13);

    if (!hAdts->protection_absent) {
      FDKpushFor(&bsWriter, 11 + 2 + (hAdts->num_raw_blocks) * 16);
      FDKcrcEndReg(&hAdts->crcInfo, &bsWriter, crcIndex);
      FDKwriteBits(&bsWriter, FDKcrcGetCRC(&hAdts->crcInfo), 16);
    }
    FDKsyncCache(&bsWriter);
  }

  *pBits -= hAdts->subFrameStartBit;
  if (!hAdts->protection_absent && hAdts->num_raw_blocks > 0) {
    *pBits += 16;
  }
  hAdts->currentBlock++;
}
