/******************* MPEG transport format encoder library *********************

   Author(s):   Alex Groeschel

   Description: ADTS Transport writer

*******************************************************************************/

#ifndef TPENC_ADTS_H
#define TPENC_ADTS_H

#include "tp_data.h"

#include "FDK_crc.h"

typedef struct {
  INT sample_freq;
  CHANNEL_MODE channel_mode;
  UCHAR decoderCanDoMpeg4;
  UCHAR mpeg_id;
  UCHAR layer;
  UCHAR protection_absent;
  UCHAR profile;
  UCHAR sample_freq_index;
  UCHAR private_bit;
  UCHAR original;
  UCHAR home;
  UCHAR copyright_id;
  UCHAR copyright_start;
  USHORT frame_length;
  UCHAR num_raw_blocks;
  UCHAR BufferFullnesStartFlag;
  UCHAR channel_config_zero;
  int headerBits;
  int currentBlock;
  int subFrameStartBit;
  FDK_CRCINFO crcInfo;
} STRUCT_ADTS;

typedef STRUCT_ADTS *HANDLE_ADTS;

INT adtsWrite_Init(HANDLE_ADTS hAdts, CODER_CONFIG *config);

int adtsWrite_GetHeaderBits(HANDLE_ADTS hAdts);

INT adtsWrite_EncodeHeader(HANDLE_ADTS hAdts, HANDLE_FDK_BITSTREAM hBitStream,
                           int bufferFullness, int frame_length);

void adtsWrite_EndRawDataBlock(HANDLE_ADTS hAdts, HANDLE_FDK_BITSTREAM hBs,
                               int *bits);

int adtsWrite_CrcStartReg(HANDLE_ADTS pAdts, HANDLE_FDK_BITSTREAM hBs,
                          int mBits);

void adtsWrite_CrcEndReg(HANDLE_ADTS pAdts, HANDLE_FDK_BITSTREAM hBs, int reg);

#endif /* TPENC_ADTS_H */
