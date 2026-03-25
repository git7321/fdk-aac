/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Dec QMF processing

*******************************************************************************/

#ifndef SAC_QMF_H
#define SAC_QMF_H

#include "common_fix.h"
#include "sac_dec_interface.h"

#include "FDK_qmf_domain.h"
#define HANDLE_SPACE_ANALYSIS_QMF HANDLE_QMF_FILTER_BANK
#define HANDLE_SPACE_SYNTHESIS_QMF HANDLE_QMF_FILTER_BANK

SACDEC_ERROR CalculateSpaceSynthesisQmf(
    const HANDLE_FDK_QMF_DOMAIN_OUT hQmfDomainOutCh, const FIXP_DBL *Sr,
    const FIXP_DBL *Si, const INT stride, PCM_MPS *timeSig);

SACDEC_ERROR CalculateSpaceAnalysisQmf(
    HANDLE_SPACE_ANALYSIS_QMF hSpaceAnalysisQmf, const PCM_MPS *timeSig,
    FIXP_DBL *Sr, FIXP_DBL *Si);

#endif /* SAC_QMF_H */
