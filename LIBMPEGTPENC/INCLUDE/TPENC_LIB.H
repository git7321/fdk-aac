/******************* MPEG transport format encoder library *********************

   Author(s):   Manuel Jander

   Description: MPEG Transport encode

*******************************************************************************/

#ifndef TPENC_LIB_H
#define TPENC_LIB_H

#include "tp_data.h"
#include "FDK_bitstream.h"

#define TRANSPORTENC_INBUF_SIZE 8192

typedef enum {
  TRANSPORTENC_OK = 0,
  TRANSPORTENC_NO_MEM,
  TRANSPORTENC_UNKOWN_ERROR = 1,
  TRANSPORTENC_INVALID_PARAMETER,
  TRANSPORTENC_UNSUPPORTED_FORMAT,
  TRANSPORTENC_NOT_ENOUGH_BITS,
  TRANSPORTENC_INVALID_CONFIG,
  TRANSPORTENC_LATM_INVALID_NR_OF_SUBFRAMES,
  TRANSPORTENC_LOAS_NOT_AVAILABLE,
  TRANSPORTENC_INVALID_LATM_ALIGNMENT,

  TRANSPORTENC_INVALID_TRANSMISSION_FRAME_LENGTH,
  TRANSPORTENC_INVALID_CELP_FRAME_LENGTH,
  TRANSPORTENC_INVALID_FRAME_BITS,
  TRANSPORTENC_INVALID_AOT,
  TRANSPORTENC_INVALID_AU_LENGTH
} TRANSPORTENC_ERROR;

typedef struct TRANSPORTENC *HANDLE_TRANSPORTENC;

CHANNEL_MODE transportEnc_GetChannelMode(int noChannels);

int transportEnc_RegisterSbrCallback(HANDLE_TRANSPORTENC hTpEnc,
                                     const cbSbr_t cbSbr, void *user_data);
int transportEnc_RegisterUsacCallback(HANDLE_TRANSPORTENC hTpEnc,
                                      const cbUsac_t cbUsac, void *user_data);
int transportEnc_RegisterSscCallback(HANDLE_TRANSPORTENC hTpEnc,
                                     const cbSsc_t cbSsc, void *user_data);
int transportEnc_writeASC(HANDLE_FDK_BITSTREAM asc, CODER_CONFIG *config,
                          CSTpCallBacks *cb);

#define TP_FLAG_MPEG4 1
#define TP_FLAG_LATM_AMV 2
#define TP_FLAG_LATM_AMVA 4

TRANSPORTENC_ERROR transportEnc_Open(HANDLE_TRANSPORTENC *phTpEnc);
TRANSPORTENC_ERROR transportEnc_Init(HANDLE_TRANSPORTENC hTpEnc,
                                     UCHAR *bsBuffer, INT bsBufferSize,
                                     TRANSPORT_TYPE transportFmt,
                                     CODER_CONFIG *config, UINT flags);
TRANSPORTENC_ERROR transportEnc_AddOtherDataBits(HANDLE_TRANSPORTENC hTpEnc,
                                                 const int nBits);

HANDLE_FDK_BITSTREAM transportEnc_GetBitstream(HANDLE_TRANSPORTENC hTp);

INT transportEnc_GetStaticBits(HANDLE_TRANSPORTENC hTp, int auBits);

void transportEnc_Close(HANDLE_TRANSPORTENC *phTp);

TRANSPORTENC_ERROR transportEnc_WriteAccessUnit(HANDLE_TRANSPORTENC hTp,
                                                INT total_bits,
                                                int bufferFullness,
                                                int noConsideredChannels);
TRANSPORTENC_ERROR transportEnc_EndAccessUnit(HANDLE_TRANSPORTENC hTp,
                                              int *pBits);
TRANSPORTENC_ERROR transportEnc_GetFrame(HANDLE_TRANSPORTENC hTpEnc,
                                         int *nbytes);

int transportEnc_CrcStartReg(HANDLE_TRANSPORTENC hTpEnc, int mBits);

void transportEnc_CrcEndReg(HANDLE_TRANSPORTENC hTpEnc, int reg);

TRANSPORTENC_ERROR transportEnc_GetConf(HANDLE_TRANSPORTENC hTpEnc,
                                        CODER_CONFIG *cc,
                                        FDK_BITSTREAM *dataBuffer,
                                        UINT *confType);
TRANSPORTENC_ERROR transportEnc_GetLibInfo(LIB_INFO *info);

#endif /* #ifndef TPENC_LIB_H */
