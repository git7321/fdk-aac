/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Dec QMF processing

*******************************************************************************/

#include "sac_qmf.h"
#include "FDK_matrixCalloc.h"
#include "sac_dec_interface.h"
#include "sac_rom.h"

#include "qmf.h"

SACDEC_ERROR CalculateSpaceSynthesisQmf(
    const HANDLE_FDK_QMF_DOMAIN_OUT hQmfDomainOutCh, const FIXP_DBL *Sr,
    const FIXP_DBL *Si, const INT stride, PCM_MPS *timeSig) {
  SACDEC_ERROR err = MPS_OK;

  if (hQmfDomainOutCh == NULL) {
    err = MPS_INVALID_HANDLE;
  } else {
    HANDLE_SPACE_SYNTHESIS_QMF hSpaceSynthesisQmf = &hQmfDomainOutCh->fb;
#if (QMF_MAX_SYNTHESIS_BANDS <= 64)
    C_AALLOC_SCRATCH_START(pWorkBuffer, FIXP_DBL,
                           (QMF_MAX_SYNTHESIS_BANDS << 1));
#else
    C_AALLOC_STACK_START(pWorkBuffer, FIXP_DBL, (QMF_MAX_SYNTHESIS_BANDS << 1));
#endif

    qmfSynthesisFilteringSlot(hSpaceSynthesisQmf, Sr, Si, 0, 0, timeSig, stride,
                              pWorkBuffer);

#if (QMF_MAX_SYNTHESIS_BANDS <= 64)
    C_AALLOC_SCRATCH_END(pWorkBuffer, FIXP_DBL, (QMF_MAX_SYNTHESIS_BANDS << 1));
#else
    C_AALLOC_STACK_END(pWorkBuffer, FIXP_DBL, (QMF_MAX_SYNTHESIS_BANDS << 1));
#endif
  }

  return err;
}

SACDEC_ERROR CalculateSpaceAnalysisQmf(
    HANDLE_SPACE_ANALYSIS_QMF hSpaceAnalysisQmf, const PCM_MPS *timeSig,
    FIXP_DBL *Sr, FIXP_DBL *Si) {
  SACDEC_ERROR err = MPS_OK;

  if (hSpaceAnalysisQmf == NULL) {
    err = MPS_INVALID_HANDLE;
  } else {
    C_AALLOC_SCRATCH_START(pWorkBuffer, FIXP_DBL, (64 << 1));

    qmfAnalysisFilteringSlot(hSpaceAnalysisQmf, Sr, Si, timeSig, 1,
                             pWorkBuffer);
    C_AALLOC_SCRATCH_END(pWorkBuffer, FIXP_DBL, (64 << 1));
  }

  return err;
}
