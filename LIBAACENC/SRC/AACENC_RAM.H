/**************************** AAC encoder library ******************************

   Author(s):   M. Lohwasser, M. Gayer

   Description:

*******************************************************************************/

#ifndef AACENC_RAM_H
#define AACENC_RAM_H

#include "common_fix.h"

#include "aacenc.h"
#include "psy_data.h"
#include "interface.h"
#include "psy_main.h"
#include "bitenc.h"
#include "bit_cnt.h"
#include "psy_const.h"

#define OUTPUTBUFFER_SIZE \
  (8192)

struct AAC_ENC {
  AACENC_CONFIG *config;

  INT ancillaryBitsPerFrame;

  CHANNEL_MAPPING channelMapping;

  QC_STATE *qcKernel;
  QC_OUT *qcOut[(1)];

  PSY_OUT *psyOut[(1)];
  PSY_INTERNAL *psyKernel;

  CHANNEL_MODE encoderMode;
  INT bandwidth90dB;
  AACENC_BITRATE_MODE bitrateMode;

  INT dontWriteAdif;
  FIXP_DBL *dynamic_RAM;

  INT maxChannels;
  INT maxElements;
  INT maxFrames;

  AUDIO_OBJECT_TYPE aot;
};

#define maxSize(a, b) (((a) > (b)) ? (a) : (b))

#define BIT_LOOK_UP_SIZE \
  (sizeof(INT) * (MAX_SFB_LONG * (CODE_BOOK_ESC_NDX + 1)))
#define MERGE_GAIN_LOOK_UP_SIZE (sizeof(INT) * MAX_SFB_LONG)

#define ADJ_THR_AH_FLAG_SIZE (sizeof(UCHAR) * ((8)) * (2) * MAX_GROUPED_SFB)
#define ADJ_THR_THR_EXP_SIZE (sizeof(FIXP_DBL) * ((8)) * (2) * MAX_GROUPED_SFB)
#define ADJ_THR_ACT_LIN_LD_DATA_SIZE \
  (sizeof(FIXP_DBL) * ((8)) * (2) * MAX_GROUPED_SFB)
#define ADJ_THR_SIZE \
  (ADJ_THR_AH_FLAG_SIZE + ADJ_THR_THR_EXP_SIZE + ADJ_THR_ACT_LIN_LD_DATA_SIZE)

#define BUF_SIZE_0 (ALIGN_SIZE(sizeof(QC_OUT_CHANNEL)) * (8))
#define BUF_SIZE_1                                                           \
  (ALIGN_SIZE(maxSize(maxSize(sizeof(PSY_DYNAMIC),                           \
                              (BIT_LOOK_UP_SIZE + MERGE_GAIN_LOOK_UP_SIZE)), \
                      ADJ_THR_SIZE)))

#define P_BUF_0 (0)
#define P_BUF_1 (P_BUF_0 + BUF_SIZE_0)

#define AAC_ENC_DYN_RAM_SIZE (BUF_SIZE_0 + BUF_SIZE_1)

H_ALLOC_MEM(AACdynamic_RAM, FIXP_DBL)

H_ALLOC_MEM(Ram_aacEnc_AacEncoder, AAC_ENC)

H_ALLOC_MEM(Ram_aacEnc_PsyElement, PSY_ELEMENT)

H_ALLOC_MEM(Ram_aacEnc_PsyInternal, PSY_INTERNAL)
H_ALLOC_MEM(Ram_aacEnc_PsyStatic, PSY_STATIC)
H_ALLOC_MEM(Ram_aacEnc_PsyInputBuffer, INT_PCM)

PSY_DYNAMIC *GetRam_aacEnc_PsyDynamic(int n, UCHAR *dynamic_RAM);

H_ALLOC_MEM(Ram_aacEnc_PsyOutChannel, PSY_OUT_CHANNEL)

H_ALLOC_MEM(Ram_aacEnc_PsyOut, PSY_OUT)
H_ALLOC_MEM(Ram_aacEnc_PsyOutElements, PSY_OUT_ELEMENT)

H_ALLOC_MEM(Ram_aacEnc_QCstate, QC_STATE)
H_ALLOC_MEM(Ram_aacEnc_AdjustThreshold, ADJ_THR_STATE)

H_ALLOC_MEM(Ram_aacEnc_AdjThrStateElement, ATS_ELEMENT)
H_ALLOC_MEM(Ram_aacEnc_ElementBits, ELEMENT_BITS)
H_ALLOC_MEM(Ram_aacEnc_BitCntrState, BITCNTR_STATE)

INT *GetRam_aacEnc_BitLookUp(int n, UCHAR *dynamic_RAM);
INT *GetRam_aacEnc_MergeGainLookUp(int n, UCHAR *dynamic_RAM);
QC_OUT_CHANNEL *GetRam_aacEnc_QCchannel(int n, UCHAR *dynamic_RAM);

H_ALLOC_MEM(Ram_aacEnc_QCout, QC_OUT)
H_ALLOC_MEM(Ram_aacEnc_QCelement, QC_OUT_ELEMENT)

#endif /* #ifndef AACENC_RAM_H */
