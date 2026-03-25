/*********************** MPEG surround encoder library *************************

   Author(s):   Josef Hoepfl

   Description: Encoder Library Interface
                vector functions

*******************************************************************************/

#ifndef SACENC_VECTORFUNCTIONS_H
#define SACENC_VECTORFUNCTIONS_H

/* Includes ******************************************************************/
#include "common_fix.h"

/* Defines *******************************************************************/
#define SUM_UP_STATIC_SCALE 0
#define SUM_UP_DYNAMIC_SCALE 1

FIXP_DBL sumUpCplxPow2(const FIXP_DPK *const x, const INT scaleMode,
                       const INT inScaleFactor, INT *const outScaleFactor,
                       const INT n);

FIXP_DBL sumUpCplxPow2Dim2(const FIXP_DPK *const *const x, const INT scaleMode,
                           const INT inScaleFactor, INT *const outScaleFactor,
                           const INT sDim1, const INT nDim1, const INT sDim2,
                           const INT nDim2);

void copyCplxVec(FIXP_DPK *const Z, const FIXP_DPK *const X, const INT n);

void setCplxVec(FIXP_DPK *const Z, const FIXP_DBL a, const INT n);

void cplx_cplxScalarProduct(FIXP_DPK *const Z, const FIXP_DPK *const *const X,
                            const FIXP_DPK *const *const Y, const INT scaleX,
                            const INT scaleY, INT *const scaleZ,
                            const INT sDim1, const INT nDim1, const INT sDim2,
                            const INT nDim2);

void FDKcalcCorrelationVec(FIXP_DBL *const z, const FIXP_DBL *const pr12,
                           const FIXP_DBL *const p1, const FIXP_DBL *const p2,
                           const INT n);

void calcCoherenceVec(FIXP_DBL *const z, const FIXP_DBL *const p12r,
                      const FIXP_DBL *const p12i, const FIXP_DBL *const p1,
                      const FIXP_DBL *const p2, const INT scaleP12,
                      const INT scaleP, const INT n);

void addWeightedCplxVec(FIXP_DPK *const *const Z, const FIXP_DBL *const a,
                        const FIXP_DPK *const *const X, const FIXP_DBL *const b,
                        const FIXP_DPK *const *const Y, const INT scale,
                        INT *const scaleCh1, const INT scaleCh2,
                        const UCHAR *const pParameterBand2HybridBandOffset,
                        const INT nParameterBands, const INT nTimeSlots,
                        const INT startTimeSlot);

void FDKcalcPbScaleFactor(const FIXP_DPK *const *const x,
                          const UCHAR *const pParameterBand2HybridBandOffset,
                          INT *const outScaleFactor, const INT startTimeSlot,
                          const INT nTimeSlots, const INT nParamBands);

INT FDKcalcScaleFactor(const FIXP_DBL *const x, const FIXP_DBL *const y,
                       const INT n);

INT FDKcalcScaleFactorDPK(const FIXP_DPK *RESTRICT x, const INT startBand,
                          const INT bands);

template <class T>
inline void FDKmemcpy_flex(T *const dst, const INT dstStride,
                           const T *const src, const INT srcStride,
                           const INT nSamples) {
  int i;

  for (i = 0; i < nSamples; i++) {
    dst[i * dstStride] = src[i * srcStride];
  }
}

template <class T>
inline void FDKmemset_flex(T *const x, const T c, const INT nSamples) {
  int i;

  for (i = 0; i < nSamples; i++) {
    x[i] = c;
  }
}

#endif /* SACENC_VECTORFUNCTIONS_H */
