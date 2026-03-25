/******************* MPEG transport format encoder library *********************

   Author(s):

   Description:

*******************************************************************************/

#ifndef TPENC_LATM_H
#define TPENC_LATM_H

#include "tpenc_lib.h"
#include "FDK_bitstream.h"

#define DEFAULT_LATM_NR_OF_SUBFRAMES 1
#define DEFAULT_LATM_SMC_REPEAT 8

#define MAX_AAC_LAYERS 9

#define LATM_MAX_PROGRAMS 1
#define LATM_MAX_STREAM_ID 16

#define LATM_MAX_LAYERS 1

#define MAX_NR_OF_SUBFRAMES \
  2

typedef enum { LATMVAR_SIMPLE_SEQUENCE } LATM_VAR_MODE;

typedef struct {
  signed int frameLengthType;
  signed int frameLengthBits;
  signed int varFrameLengthTable[4];
  signed int streamID;
} LATM_LAYER_INFO;

typedef struct {
  LATM_LAYER_INFO m_linfo[LATM_MAX_PROGRAMS][LATM_MAX_LAYERS];
  CODER_CONFIG *config[LATM_MAX_PROGRAMS][LATM_MAX_LAYERS];

  LATM_VAR_MODE varMode;
  TRANSPORT_TYPE tt;

  int audioMuxLengthBytes;

  int audioMuxLengthBytesPos;
  int taraBufferFullness;
  int varStreamCnt;

  UCHAR
  latmFrameCounter;
  UCHAR muxConfigPeriod;

  UCHAR
  audioMuxVersion;
  UCHAR audioMuxVersionA;

  UCHAR noProgram;
  UCHAR noLayer[LATM_MAX_PROGRAMS];
  UCHAR fractDelayPresent;

  UCHAR allStreamsSameTimeFraming;
  UCHAR subFrameCnt;
  UCHAR noSubframes;
  UINT latmSubframeStart;
  UCHAR noSubframes_next;

  UCHAR otherDataLenBits;
  UCHAR fillBits;
  UINT streamMuxConfigBits;

} LATM_STREAM;

typedef LATM_STREAM *HANDLE_LATM_STREAM;

TRANSPORTENC_ERROR transportEnc_Latm_Init(HANDLE_LATM_STREAM hLatmStreamInfo,
                                          HANDLE_FDK_BITSTREAM hBs,
                                          CODER_CONFIG *layerConfig,
                                          UINT audioMuxVersion,
                                          TRANSPORT_TYPE tt, CSTpCallBacks *cb);

TRANSPORTENC_ERROR transportEnc_LatmAddOtherDataBits(HANDLE_LATM_STREAM hAss,
                                                     const int otherDataBits);

unsigned int transportEnc_LatmCountTotalBitDemandHeader(
    HANDLE_LATM_STREAM hAss, unsigned int streamDataLength);

TRANSPORTENC_ERROR
transportEnc_LatmWrite(HANDLE_LATM_STREAM hAss, HANDLE_FDK_BITSTREAM hBitstream,
                       int auBits, int bufferFullness, CSTpCallBacks *cb);

void transportEnc_LatmAdjustSubframeBits(HANDLE_LATM_STREAM hAss, int *pBits);

TRANSPORTENC_ERROR transportEnc_LatmGetFrame(HANDLE_LATM_STREAM hAss,
                                             HANDLE_FDK_BITSTREAM hBs,
                                             int *pBytes);

TRANSPORTENC_ERROR
CreateStreamMuxConfig(HANDLE_LATM_STREAM hAss, HANDLE_FDK_BITSTREAM hBs,
                      int bufferFullness, CSTpCallBacks *cb);

#endif /* TPENC_LATM_H */
