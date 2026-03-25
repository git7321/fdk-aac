/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "sbrenc_ram.h"
#include "sbr.h"
#include "genericStds.h"

C_AALLOC_MEM(Ram_SbrDynamic_RAM, FIXP_DBL,
             ((SBR_ENC_DYN_RAM_SIZE) / sizeof(FIXP_DBL)))

C_ALLOC_MEM(Ram_SbrEncoder, SBR_ENCODER, 1)
C_ALLOC_MEM2(Ram_SbrChannel, SBR_CHANNEL, 1, (8))
C_ALLOC_MEM2(Ram_SbrElement, SBR_ELEMENT, 1, (8))

C_AALLOC_MEM2_L(Ram_Sbr_QmfStatesAnalysis, FIXP_QAS, 640, (8), SECT_DATA_L1)
C_ALLOC_MEM2_L(Ram_Sbr_quotaMatrix, FIXP_DBL, (MAX_NO_OF_ESTIMATES * 64), (8),
               SECT_DATA_L1)
C_ALLOC_MEM2(Ram_Sbr_signMatrix, INT, (MAX_NO_OF_ESTIMATES * 64), (8))
C_ALLOC_MEM2(Ram_Sbr_freqBandTableLO, UCHAR, (MAX_FREQ_COEFFS / 2 + 1), (8))
C_ALLOC_MEM2(Ram_Sbr_freqBandTableHI, UCHAR, (MAX_FREQ_COEFFS + 1), (8))
C_ALLOC_MEM2(Ram_Sbr_v_k_master, UCHAR, (MAX_FREQ_COEFFS + 1), (8))
C_ALLOC_MEM2(Ram_Sbr_detectionVectors, UCHAR,
             (MAX_NO_OF_ESTIMATES * MAX_FREQ_COEFFS), (8))
C_ALLOC_MEM2(Ram_Sbr_prevEnvelopeCompensation, UCHAR, MAX_FREQ_COEFFS, (8))
C_ALLOC_MEM2(Ram_Sbr_guideScfb, UCHAR, MAX_FREQ_COEFFS, (8))
C_ALLOC_MEM2(Ram_Sbr_guideVectorDetected, UCHAR,
             (MAX_NO_OF_ESTIMATES * MAX_FREQ_COEFFS), (8))
C_ALLOC_MEM2(Ram_Sbr_guideVectorDiff, FIXP_DBL,
             (MAX_NO_OF_ESTIMATES * MAX_FREQ_COEFFS), (8))
C_ALLOC_MEM2(Ram_Sbr_guideVectorOrig, FIXP_DBL,
             (MAX_NO_OF_ESTIMATES * MAX_FREQ_COEFFS), (8))
C_AALLOC_MEM_L(Ram_PsQmfStatesSynthesis, FIXP_DBL, 640 / 2, SECT_DATA_L1)
C_ALLOC_MEM_L(Ram_PsEncode, PS_ENCODE, 1, SECT_DATA_L1)
C_ALLOC_MEM(Ram_ParamStereo, PARAMETRIC_STEREO, 1)
C_ALLOC_MEM2(Ram_Sbr_envYBuffer, FIXP_DBL, (32 / 2 * 64), (8))

FIXP_DBL* GetRam_Sbr_envYBuffer(int n, UCHAR* dynamic_RAM) {
  FDK_ASSERT(dynamic_RAM != 0);
  return reinterpret_cast<FIXP_DBL*>(
      reinterpret_cast<void*>(dynamic_RAM + OFFSET_NRG + (n * Y_2_BUF_BYTE)));
}

FIXP_DBL* GetRam_Sbr_envRBuffer(int n, UCHAR* dynamic_RAM) {
  FDK_ASSERT(dynamic_RAM != 0);
  return reinterpret_cast<FIXP_DBL*>(reinterpret_cast<void*>(
      dynamic_RAM + OFFSET_QMF + (n * (ENV_R_BUFF_BYTE + ENV_I_BUFF_BYTE))));
}
FIXP_DBL* GetRam_Sbr_envIBuffer(int n, UCHAR* dynamic_RAM) {
  FDK_ASSERT(dynamic_RAM != 0);
  return reinterpret_cast<FIXP_DBL*>(
      reinterpret_cast<void*>(dynamic_RAM + OFFSET_QMF + (ENV_R_BUFF_BYTE) +
                              (n * (ENV_R_BUFF_BYTE + ENV_I_BUFF_BYTE))));
}
