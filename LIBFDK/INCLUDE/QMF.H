/******************* Library for basic calculation routines ********************

   Author(s):

   Description:

*******************************************************************************/

#ifndef QMF_H
#define QMF_H

#include "common_fix.h"
#include "FDK_tools_rom.h"
#include "dct.h"

#define FIXP_QAS FIXP_PCM
#define QAS_BITS SAMPLE_BITS
#define INT_PCM_QMFIN INT_PCM

#define FIXP_QSS FIXP_DBL
#define QSS_BITS DFRACT_BITS

#define QMF_FLAG_LP 1
#define QMF_FLAG_NONSYMMETRIC 2
#define QMF_FLAG_CLDFB 4
#define QMF_FLAG_KEEP_STATES 8
#define QMF_FLAG_MPSLDFB 16
#define QMF_FLAG_MPSLDFB_OPTIMIZE_MODULATION 32
#define QMF_FLAG_DOWNSAMPLED 64

#define QMF_MAX_SYNTHESIS_BANDS (64)

#define ALGORITHMIC_SCALING_IN_ANALYSIS_FILTERBANK 7
#define ALGORITHMIC_SCALING_IN_SYNTHESIS_FILTERBANK 1

typedef struct {
  int lb_scale;
  int ov_lb_scale;
  int hb_scale;
  int ov_hb_scale;
} QMF_SCALE_FACTOR;

struct QMF_FILTER_BANK {
  const FIXP_PFT *p_filter;

  void *FilterStates;
  int FilterSize;
  const FIXP_QTW *t_cos;
  const FIXP_QTW *t_sin;
  int filterScale;

  int no_channels;
  int no_col;
  int lsb;
  int usb;

  int synScalefactor;
  int outScalefactor;
  FIXP_DBL outGain_m;
  int outGain_e;

  UINT flags;
  UCHAR p_stride;
};

typedef struct QMF_FILTER_BANK *HANDLE_QMF_FILTER_BANK;

int qmfInitAnalysisFilterBank(
    HANDLE_QMF_FILTER_BANK h_Qmf,
    FIXP_QAS *pFilterStates,
    int noCols,
    int lsb,
    int usb,
    int no_channels,
    int flags);
#if SAMPLE_BITS == 16

int qmfInitAnalysisFilterBank(
    HANDLE_QMF_FILTER_BANK h_Qmf,
    FIXP_DBL *pFilterStates,
    int noCols,
    int lsb,
    int usb,
    int no_channels,
    int flags);
#endif

void qmfAnalysisFiltering(
    HANDLE_QMF_FILTER_BANK anaQmf,
    FIXP_DBL **qmfReal,
    FIXP_DBL **qmfImag,
    QMF_SCALE_FACTOR *scaleFactor,
    const INT_PCM *timeIn,
    const int timeIn_e,
    const int stride,
    FIXP_DBL *pWorkBuffer
);
#if SAMPLE_BITS == 16

void qmfAnalysisFiltering(
    HANDLE_QMF_FILTER_BANK anaQmf,
    FIXP_DBL **qmfReal,
    FIXP_DBL **qmfImag,
    QMF_SCALE_FACTOR *scaleFactor,
    const LONG *timeIn,
    const int timeIn_e,
    const int stride,
    FIXP_DBL *pWorkBuffer
);
#endif

void qmfAnalysisFilteringSlot(
    HANDLE_QMF_FILTER_BANK anaQmf,
    FIXP_DBL *qmfReal,
    FIXP_DBL *qmfImag,
    const INT_PCM *timeIn,
    const int stride,
    FIXP_DBL *pWorkBuffer
);
#if SAMPLE_BITS == 16

void qmfAnalysisFilteringSlot(
    HANDLE_QMF_FILTER_BANK anaQmf,
    FIXP_DBL *qmfReal,
    FIXP_DBL *qmfImag,
    const LONG *timeIn,
    const int stride,
    FIXP_DBL *pWorkBuffer
);
#endif

int qmfInitSynthesisFilterBank(
    HANDLE_QMF_FILTER_BANK h_Qmf,
    FIXP_QSS *pFilterStates,
    int noCols,
    int lsb,
    int usb,
    int no_channels,
    int flags);

void qmfSynthesisFiltering(
    HANDLE_QMF_FILTER_BANK synQmf,
    FIXP_DBL **QmfBufferReal,
    FIXP_DBL **QmfBufferImag,
    const QMF_SCALE_FACTOR *scaleFactor,
    const int ov_len,
    INT_PCM *timeOut,
    const INT stride,
    FIXP_DBL *pWorkBuffer
);
#if SAMPLE_BITS == 16

void qmfSynthesisFiltering(
    HANDLE_QMF_FILTER_BANK synQmf,
    FIXP_DBL **QmfBufferReal,
    FIXP_DBL **QmfBufferImag,
    const QMF_SCALE_FACTOR *scaleFactor,
    const int ov_len,
    LONG *timeOut,
    const int timeOut_e,
    FIXP_DBL *pWorkBuffer
);
#endif

void qmfSynthesisFilteringSlot(HANDLE_QMF_FILTER_BANK synQmf,
                               const FIXP_DBL *realSlot,
                               const FIXP_DBL *imagSlot,
                               const int scaleFactorLowBand,
                               const int scaleFactorHighBand, INT_PCM *timeOut,
                               const int timeOut_e, FIXP_DBL *pWorkBuffer);
#if SAMPLE_BITS == 16

void qmfSynthesisFilteringSlot(HANDLE_QMF_FILTER_BANK synQmf,
                               const FIXP_DBL *realSlot,
                               const FIXP_DBL *imagSlot,
                               const int scaleFactorLowBand,
                               const int scaleFactorHighBand, LONG *timeOut,
                               const int timeOut_e, FIXP_DBL *pWorkBuffer);
#endif

void qmfChangeOutScalefactor(
    HANDLE_QMF_FILTER_BANK synQmf,
    int outScalefactor
);

int qmfGetOutScalefactor(
    HANDLE_QMF_FILTER_BANK synQmf
);

void qmfChangeOutGain(
    HANDLE_QMF_FILTER_BANK synQmf,
    FIXP_DBL outputGain,
    int outputGainScale
);

#endif /*ifndef QMF_H       */
