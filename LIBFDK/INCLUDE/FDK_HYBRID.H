/******************* Library for basic calculation routines ********************

   Author(s):   Markus Lohwasser

   Description: FDK Tools Hybrid Filterbank

*******************************************************************************/

#ifndef FDK_HYBRID_H
#define FDK_HYBRID_H

#include "common_fix.h"

typedef enum {
  THREE_TO_TEN,
  THREE_TO_TWELVE,
  THREE_TO_SIXTEEN

} FDK_HYBRID_MODE;

typedef const struct FDK_HYBRID_SETUP *HANDLE_FDK_HYBRID_SETUP;

typedef struct {
  FIXP_DBL *bufferLFReal[3];
  FIXP_DBL *bufferLFImag[3];
  FIXP_DBL *bufferHFReal[13];
  FIXP_DBL *bufferHFImag[13];

  INT bufferLFpos;
  INT bufferHFpos;
  INT nrBands;
  INT cplxBands;
  UCHAR hfMode;

  FIXP_DBL *pLFmemory;
  FIXP_DBL *pHFmemory;

  UINT LFmemorySize;
  UINT HFmemorySize;

  HANDLE_FDK_HYBRID_SETUP pSetup;

} FDK_ANA_HYB_FILTER;

typedef struct {
  INT nrBands;
  INT cplxBands;

  HANDLE_FDK_HYBRID_SETUP pSetup;

} FDK_SYN_HYB_FILTER;

typedef FDK_ANA_HYB_FILTER *HANDLE_FDK_ANA_HYB_FILTER;
typedef FDK_SYN_HYB_FILTER *HANDLE_FDK_SYN_HYB_FILTER;

INT FDKhybridAnalysisOpen(HANDLE_FDK_ANA_HYB_FILTER hAnalysisHybFilter,
                          FIXP_DBL *const pLFmemory, const UINT LFmemorySize,
                          FIXP_DBL *const pHFmemory, const UINT HFmemorySize);

INT FDKhybridAnalysisInit(HANDLE_FDK_ANA_HYB_FILTER hAnalysisHybFilter,
                          const FDK_HYBRID_MODE mode, const INT qmfBands,
                          const INT cplxBands, const INT initStatesFlag);

INT FDKhybridAnalysisScaleStates(HANDLE_FDK_ANA_HYB_FILTER hAnalysisHybFilter,
                                 const INT scalingValue);

INT FDKhybridAnalysisApply(HANDLE_FDK_ANA_HYB_FILTER hAnalysisHybFilter,
                           const FIXP_DBL *const pQmfReal,
                           const FIXP_DBL *const pQmfImag,
                           FIXP_DBL *const pHybridReal,
                           FIXP_DBL *const pHybridImag);

INT FDKhybridAnalysisClose(HANDLE_FDK_ANA_HYB_FILTER hAnalysisHybFilter);

INT FDKhybridSynthesisInit(HANDLE_FDK_SYN_HYB_FILTER hSynthesisHybFilter,
                           const FDK_HYBRID_MODE mode, const INT qmfBands,
                           const INT cplxBands);

void FDKhybridSynthesisApply(HANDLE_FDK_SYN_HYB_FILTER hSynthesisHybFilter,
                             const FIXP_DBL *const pHybridReal,
                             const FIXP_DBL *const pHybridImag,
                             FIXP_DBL *const pQmfReal,
                             FIXP_DBL *const pQmfImag);

#endif /* FDK_HYBRID_H */
