/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef SBRENC_RAM_H
#define SBRENC_RAM_H

#include "sbr_def.h"
#include "env_est.h"
#include "sbr_encoder.h"
#include "sbr.h"

#include "ps_main.h"
#include "ps_encode.h"

#define ENV_TRANSIENTS_BYTE ((sizeof(FIXP_DBL) * (MAX_NUM_CHANNELS * 3 * 32)))

#define ENV_R_BUFF_BYTE ((sizeof(FIXP_DBL) * ((32) * MAX_HYBRID_BANDS)))
#define ENV_I_BUFF_BYTE ((sizeof(FIXP_DBL) * ((32) * MAX_HYBRID_BANDS)))
#define Y_BUF_CH_BYTE \
  ((2 * sizeof(FIXP_DBL) * (((32) - (32 / 2)) * MAX_HYBRID_BANDS)))

#define ENV_R_BUF_PS_BYTE ((sizeof(FIXP_DBL) * 32 * 64 / 2))
#define ENV_I_BUF_PS_BYTE ((sizeof(FIXP_DBL) * 32 * 64 / 2))

#define TON_BUF_CH_BYTE \
  ((sizeof(FIXP_DBL) * (MAX_NO_OF_ESTIMATES * MAX_FREQ_COEFFS)))

#define Y_2_BUF_BYTE (Y_BUF_CH_BYTE)

#define BUF_NRG_SIZE ((MAX_NUM_CHANNELS * Y_2_BUF_BYTE))
#define BUF_QMF_SIZE (ENV_R_BUFF_BYTE + ENV_I_BUFF_BYTE)

#define SBR_ENC_DYN_RAM_SIZE (BUF_QMF_SIZE + BUF_NRG_SIZE)

#define OFFSET_QMF (0)
#define OFFSET_NRG (OFFSET_QMF + BUF_QMF_SIZE)

H_ALLOC_MEM(Ram_SbrDynamic_RAM, FIXP_DBL)

H_ALLOC_MEM(Ram_SbrEncoder, SBR_ENCODER)
H_ALLOC_MEM(Ram_SbrChannel, SBR_CHANNEL)
H_ALLOC_MEM(Ram_SbrElement, SBR_ELEMENT)

H_ALLOC_MEM(Ram_Sbr_quotaMatrix, FIXP_DBL)
H_ALLOC_MEM(Ram_Sbr_signMatrix, INT)

H_ALLOC_MEM(Ram_Sbr_QmfStatesAnalysis, FIXP_QAS)

H_ALLOC_MEM(Ram_Sbr_freqBandTableLO, UCHAR)
H_ALLOC_MEM(Ram_Sbr_freqBandTableHI, UCHAR)
H_ALLOC_MEM(Ram_Sbr_v_k_master, UCHAR)

H_ALLOC_MEM(Ram_Sbr_detectionVectors, UCHAR)
H_ALLOC_MEM(Ram_Sbr_prevEnvelopeCompensation, UCHAR)
H_ALLOC_MEM(Ram_Sbr_guideScfb, UCHAR)
H_ALLOC_MEM(Ram_Sbr_guideVectorDetected, UCHAR)

H_ALLOC_MEM(Ram_Sbr_envYBuffer, FIXP_DBL)
FIXP_DBL* GetRam_Sbr_envYBuffer(int n, UCHAR* dynamic_RAM);
FIXP_DBL* GetRam_Sbr_envRBuffer(int n, UCHAR* dynamic_RAM);
FIXP_DBL* GetRam_Sbr_envIBuffer(int n, UCHAR* dynamic_RAM);

H_ALLOC_MEM(Ram_Sbr_guideVectorDiff, FIXP_DBL)
H_ALLOC_MEM(Ram_Sbr_guideVectorOrig, FIXP_DBL)

H_ALLOC_MEM(Ram_PsQmfStatesSynthesis, FIXP_DBL)

H_ALLOC_MEM(Ram_PsEncode, PS_ENCODE)

FIXP_DBL* FDKsbrEnc_SliceRam_PsRqmf(FIXP_DBL* rQmfData, UCHAR* dynamic_RAM,
                                    int n, int i, int qmfSlots);
FIXP_DBL* FDKsbrEnc_SliceRam_PsIqmf(FIXP_DBL* iQmfData, UCHAR* dynamic_RAM,
                                    int n, int i, int qmfSlots);

H_ALLOC_MEM(Ram_ParamStereo, PARAMETRIC_STEREO)
#endif
