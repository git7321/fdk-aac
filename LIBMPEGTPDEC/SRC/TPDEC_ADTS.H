/******************* MPEG transport format decoder library *********************

   Author(s):   Josef Hoepfl

   Description: ADTS interface

*******************************************************************************/

#ifndef TPDEC_ADTS_H
#define TPDEC_ADTS_H

#include "tpdec_lib.h"

#define ADTS_SYNCWORD (0xfff)
#define ADTS_SYNCLENGTH (12)
#define ADTS_HEADERLENGTH (56)
#define ADTS_FIXED_HEADERLENGTH (28)
#define ADTS_VARIABLE_HEADERLENGTH (ADTS_HEADERLENGTH - ADTS_FIXED_HEADERLENGTH)

#ifdef CHECK_TWO_SYNCS
#define ADTS_MIN_TP_BUF_SIZE (8191 + 2)
#else
#define ADTS_MIN_TP_BUF_SIZE (8191)
#endif

#include "FDK_crc.h"

typedef struct {
  UCHAR mpeg_id;
  UCHAR layer;
  UCHAR protection_absent;
  UCHAR profile;
  UCHAR sample_freq_index;
  UCHAR private_bit;
  UCHAR channel_config;
  UCHAR original;
  UCHAR home;
  UCHAR copyright_id;
  UCHAR copyright_start;
  USHORT frame_length;
  USHORT adts_fullness;
  UCHAR num_raw_blocks;
  UCHAR num_pce_bits;
} STRUCT_ADTS_BS;

struct STRUCT_ADTS {
  STRUCT_ADTS_BS bs;

  UCHAR decoderCanDoMpeg4;
  UCHAR BufferFullnesStartFlag;

  FDK_CRCINFO crcInfo;
  USHORT crcReadValue;
  USHORT rawDataBlockDist[4];
};

typedef struct STRUCT_ADTS *HANDLE_ADTS;

void adtsRead_CrcInit(HANDLE_ADTS pAdts);

int adtsRead_CrcStartReg(HANDLE_ADTS pAdts, HANDLE_FDK_BITSTREAM hBs,
                         int mBits);

void adtsRead_CrcEndReg(HANDLE_ADTS pAdts, HANDLE_FDK_BITSTREAM hBs, int reg);

TRANSPORTDEC_ERROR adtsRead_CrcCheck(HANDLE_ADTS pAdts);
TRANSPORTDEC_ERROR adtsRead_DecodeHeader(HANDLE_ADTS pAdts,
                                         CSAudioSpecificConfig *pAsc,
                                         HANDLE_FDK_BITSTREAM bs,
                                         const INT ignoreBufferFullness);

int adtsRead_GetRawDataBlockLength(HANDLE_ADTS pAdts, INT blockNum);

#endif /* TPDEC_ADTS_H */
