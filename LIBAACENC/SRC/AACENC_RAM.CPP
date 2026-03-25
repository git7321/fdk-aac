/**************************** AAC encoder library ******************************

   Author(s):   M. Lohwasser, M. Gayer

   Description:

*******************************************************************************/

#include "aacEnc_ram.h"

C_AALLOC_MEM(AACdynamic_RAM, FIXP_DBL, AAC_ENC_DYN_RAM_SIZE / sizeof(FIXP_DBL))

C_ALLOC_MEM(Ram_aacEnc_AacEncoder, struct AAC_ENC, 1)

C_ALLOC_MEM2(Ram_aacEnc_PsyElement, PSY_ELEMENT, 1, ((8)))

C_ALLOC_MEM(Ram_aacEnc_PsyInternal, PSY_INTERNAL, 1)
C_ALLOC_MEM2(Ram_aacEnc_PsyStatic, PSY_STATIC, 1, (8))

C_ALLOC_MEM2(Ram_aacEnc_PsyInputBuffer, INT_PCM, MAX_INPUT_BUFFER_SIZE, (8))

PSY_DYNAMIC *GetRam_aacEnc_PsyDynamic(int n, UCHAR *dynamic_RAM) {
  FDK_ASSERT(dynamic_RAM != 0);
  return reinterpret_cast<PSY_DYNAMIC *>(reinterpret_cast<void *>(
      dynamic_RAM + P_BUF_1 + n * sizeof(PSY_DYNAMIC)));
}

C_ALLOC_MEM2(Ram_aacEnc_PsyOut, PSY_OUT, 1, (1))
C_ALLOC_MEM2(Ram_aacEnc_PsyOutElements, PSY_OUT_ELEMENT, 1, (1) * ((8)))
C_ALLOC_MEM2(Ram_aacEnc_PsyOutChannel, PSY_OUT_CHANNEL, 1, (1) * (8))

C_ALLOC_MEM(Ram_aacEnc_QCstate, QC_STATE, 1)
C_ALLOC_MEM(Ram_aacEnc_AdjustThreshold, ADJ_THR_STATE, 1)

C_ALLOC_MEM2(Ram_aacEnc_AdjThrStateElement, ATS_ELEMENT, 1, ((8)))
C_ALLOC_MEM2(Ram_aacEnc_ElementBits, ELEMENT_BITS, 1, ((8)))
C_ALLOC_MEM(Ram_aacEnc_BitCntrState, struct BITCNTR_STATE, 1)

INT *GetRam_aacEnc_BitLookUp(int n, UCHAR *dynamic_RAM) {
  FDK_ASSERT(dynamic_RAM != 0);
  return reinterpret_cast<INT *>(
      reinterpret_cast<void *>(dynamic_RAM + P_BUF_1));
}
INT *GetRam_aacEnc_MergeGainLookUp(int n, UCHAR *dynamic_RAM) {
  FDK_ASSERT(dynamic_RAM != 0);
  return reinterpret_cast<INT *>(reinterpret_cast<void *>(
      dynamic_RAM + P_BUF_1 +
      sizeof(INT) * (MAX_SFB_LONG * (CODE_BOOK_ESC_NDX + 1))));
}

C_ALLOC_MEM2(Ram_aacEnc_QCout, QC_OUT, 1, (1))
C_ALLOC_MEM2(Ram_aacEnc_QCelement, QC_OUT_ELEMENT, 1, (1) * ((8)))
QC_OUT_CHANNEL *GetRam_aacEnc_QCchannel(int n, UCHAR *dynamic_RAM) {
  FDK_ASSERT(dynamic_RAM != 0);
  return reinterpret_cast<QC_OUT_CHANNEL *>(reinterpret_cast<void *>(
      dynamic_RAM + P_BUF_0 + n * ALIGN_SIZE(sizeof(QC_OUT_CHANNEL))));
}
