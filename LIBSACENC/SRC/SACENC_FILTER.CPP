/*********************** MPEG surround encoder library *************************

   Author(s):   M. Multrus

   Description: Encoder Library
                Filter functions

*******************************************************************************/

/* Includes ******************************************************************/
#include "sacenc_filter.h"

typedef struct T_DC_FILTER {
  FIXP_DBL c__FDK;
  FIXP_DBL state__FDK;

} DC_FILTER;

FDK_SACENC_ERROR fdk_sacenc_createDCFilter(HANDLE_DC_FILTER *hDCFilter) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == hDCFilter) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDK_ALLOCATE_MEMORY_1D(*hDCFilter, 1, DC_FILTER);
  }
  return error;

bail:
  fdk_sacenc_destroyDCFilter(hDCFilter);
  return ((SACENC_OK == error) ? SACENC_MEMORY_ERROR : error);
}

FDK_SACENC_ERROR fdk_sacenc_initDCFilter(HANDLE_DC_FILTER hDCFilter,
                                         const UINT sampleRate) {
  FDK_SACENC_ERROR error = SACENC_OK;

  FIXP_DBL expC;
  int s;

  expC = fDivNormHighPrec((FIXP_DBL)20, (FIXP_DBL)sampleRate, &s);
  expC = fMultDiv2(FL2FXCONST_DBL(-1.44269504089 * 0.5), expC) >>
         (LD_DATA_SHIFT - 1 - 1);

  if (s < 0)
    expC = expC >> (-s);
  else
    expC = expC << (s);

  expC = CalcInvLdData(expC);

  hDCFilter->c__FDK = expC;
  hDCFilter->state__FDK = FL2FXCONST_DBL(0.0f);

  return error;
}

FDK_SACENC_ERROR fdk_sacenc_destroyDCFilter(HANDLE_DC_FILTER *hDCFilter) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((hDCFilter != NULL) && (*hDCFilter != NULL)) {
    FDKfree(*hDCFilter);

    *hDCFilter = NULL;
  }

  return error;
}

FDK_SACENC_ERROR fdk_sacenc_applyDCFilter(HANDLE_DC_FILTER hDCFilter,
                                          const INT_PCM *const signalIn,
                                          INT_PCM *const signalOut,
                                          const INT signalLength) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((hDCFilter == NULL) || (signalIn == NULL) || (signalOut == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    const INT_PCM *const x = signalIn;
    INT_PCM *const y = signalOut;
    const FIXP_DBL c = hDCFilter->c__FDK;
    FIXP_DBL *const state = &hDCFilter->state__FDK;
    int i;
    FIXP_DBL x0, x1, y1;

    x1 = x0 = FX_PCM2FX_DBL(x[0]) >> DC_FILTER_SF;
    y1 = x0 + (*state);

    for (i = 1; i < signalLength; i++) {
      x0 = FX_PCM2FX_DBL(x[i]) >> DC_FILTER_SF;
      y[i - 1] = FX_DBL2FX_PCM(y1);
      y1 = x0 - x1 + fMult(c, y1);
      x1 = x0;
    }

    *state = fMult(c, y1) - x1;
    y[i - 1] = FX_DBL2FX_PCM(y1);
  }

  return error;
}
