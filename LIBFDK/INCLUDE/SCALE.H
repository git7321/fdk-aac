/******************* Library for basic calculation routines ********************

   Author(s):

   Description: Scaling operations

*******************************************************************************/

#ifndef SCALE_H
#define SCALE_H

#include "common_fix.h"
#include "genericStds.h"
#include "fixminmax.h"

#if defined(__arm__)
#include "arm/scale_arm.h"

#elif defined(__mips__)
#include "mips/scale_mips.h"

#endif

void scaleValues(FIXP_SGL *vector, INT len, INT scalefactor);
void scaleValues(FIXP_DBL *vector, INT len, INT scalefactor);
void scaleValues(FIXP_DBL *dst, const FIXP_DBL *src, INT len, INT scalefactor);
#if (SAMPLE_BITS == 16)
void scaleValues(FIXP_PCM *dst, const FIXP_DBL *src, INT len, INT scalefactor);
#endif
void scaleValues(FIXP_SGL *dst, const FIXP_SGL *src, INT len, INT scalefactor);
void scaleCplxValues(FIXP_DBL *r_dst, FIXP_DBL *i_dst, const FIXP_DBL *r_src,
                     const FIXP_DBL *i_src, INT len, INT scalefactor);
void scaleValuesWithFactor(FIXP_DBL *vector, FIXP_DBL factor, INT len,
                           INT scalefactor);
void scaleValuesSaturate(FIXP_DBL *vector, INT len, INT scalefactor);
void scaleValuesSaturate(FIXP_DBL *dst, const FIXP_DBL *src, INT len,
                         INT scalefactor);
void scaleValuesSaturate(FIXP_SGL *dst, const FIXP_DBL *src, INT len,
                         INT scalefactor);
void scaleValuesSaturate(FIXP_SGL *vector, INT len, INT scalefactor);
void scaleValuesSaturate(FIXP_SGL *dst, const FIXP_SGL *src, INT len,
                         INT scalefactor);
INT getScalefactorShort(const SHORT *vector, INT len);
INT getScalefactorPCM(const INT_PCM *vector, INT len, INT stride);
INT getScalefactor(const FIXP_DBL *vector, INT len);
INT getScalefactor(const FIXP_SGL *vector, INT len);

#ifndef FUNCTION_scaleValue
#define FUNCTION_scaleValue
inline FIXP_DBL scaleValue(const FIXP_DBL value,
                           INT scalefactor
) {
  if (scalefactor > 0)
    return (value << scalefactor);
  else
    return (value >> (-scalefactor));
}
inline FIXP_SGL scaleValue(const FIXP_SGL value,
                           INT scalefactor
) {
  if (scalefactor > 0)
    return (value << scalefactor);
  else
    return (value >> (-scalefactor));
}
#endif

#ifndef FUNCTION_scaleValueSaturate
#define FUNCTION_scaleValueSaturate
inline FIXP_DBL scaleValueSaturate(const FIXP_DBL value,
                                   INT scalefactor
) {
  int headroom = fixnormz_D(
      (INT)value ^ (INT)((value >> 31)));
  if (scalefactor >= 0) {
    if (headroom <= scalefactor) {
      if (value > (FIXP_DBL)0)
        return (FIXP_DBL)MAXVAL_DBL;
      else
        return (FIXP_DBL)MINVAL_DBL + (FIXP_DBL)1;
    } else {
      return fMax((value << scalefactor), (FIXP_DBL)MINVAL_DBL + (FIXP_DBL)1);
    }
  } else {
    scalefactor = -scalefactor;
    if ((DFRACT_BITS - headroom) <= scalefactor) {
      return (FIXP_DBL)0;
    } else {
      return fMax((value >> scalefactor), (FIXP_DBL)MINVAL_DBL + (FIXP_DBL)1);
    }
  }
}
#endif

#ifndef FUNCTION_scaleValueInPlace
#define FUNCTION_scaleValueInPlace
inline void scaleValueInPlace(FIXP_DBL *value,
                              INT scalefactor
) {
  INT newscale;
  if ((newscale = (scalefactor)) >= 0) {
    *(value) <<= newscale;
  } else {
    *(value) >>= -newscale;
  }
}
#endif

#ifndef SATURATE_RIGHT_SHIFT
#define SATURATE_RIGHT_SHIFT(src, scale, dBits)                            \
  ((((LONG)(src) >> (scale)) > (LONG)(((1U) << ((dBits)-1)) - 1))          \
       ? (LONG)(((1U) << ((dBits)-1)) - 1)                                 \
       : (((LONG)(src) >> (scale)) < ~((LONG)(((1U) << ((dBits)-1)) - 1))) \
             ? ~((LONG)(((1U) << ((dBits)-1)) - 1))                        \
             : ((LONG)(src) >> (scale)))
#endif

#ifndef SATURATE_LEFT_SHIFT
#define SATURATE_LEFT_SHIFT(src, scale, dBits)                           \
  (((LONG)(src) > ((LONG)(((1U) << ((dBits)-1)) - 1) >> (scale)))        \
       ? (LONG)(((1U) << ((dBits)-1)) - 1)                               \
       : ((LONG)(src) < ~((LONG)(((1U) << ((dBits)-1)) - 1) >> (scale))) \
             ? ~((LONG)(((1U) << ((dBits)-1)) - 1))                      \
             : ((LONG)(src) << (scale)))
#endif

#ifndef SATURATE_SHIFT
#define SATURATE_SHIFT(src, scale, dBits)                        \
  (((scale) < 0) ? SATURATE_LEFT_SHIFT((src), -(scale), (dBits)) \
                 : SATURATE_RIGHT_SHIFT((src), (scale), (dBits)))
#endif

#ifndef SATURATE_LEFT_SHIFT_ALT
#define SATURATE_LEFT_SHIFT_ALT(src, scale, dBits)                        \
  (((LONG)(src) > ((LONG)(((1U) << ((dBits)-1)) - 1) >> (scale)))         \
       ? (LONG)(((1U) << ((dBits)-1)) - 1)                                \
       : ((LONG)(src) <= ~((LONG)(((1U) << ((dBits)-1)) - 1) >> (scale))) \
             ? ~((LONG)(((1U) << ((dBits)-1)) - 2))                       \
             : ((LONG)(src) << (scale)))
#endif

#ifndef SATURATE_RIGHT_SHIFT_ALT
#define SATURATE_RIGHT_SHIFT_ALT(src, scale, dBits)                        \
  ((((LONG)(src) >> (scale)) > (LONG)(((1U) << ((dBits)-1)) - 1))          \
       ? (LONG)(((1U) << ((dBits)-1)) - 1)                                 \
       : (((LONG)(src) >> (scale)) < ~((LONG)(((1U) << ((dBits)-1)) - 2))) \
             ? ~((LONG)(((1U) << ((dBits)-1)) - 2))                        \
             : ((LONG)(src) >> (scale)))
#endif

#ifndef SATURATE_INT_PCM_RIGHT_SHIFT
#define SATURATE_INT_PCM_RIGHT_SHIFT(src, scale) \
  SATURATE_RIGHT_SHIFT(src, scale, SAMPLE_BITS)
#endif

#ifndef SATURATE_INT_PCM_LEFT_SHIFT
#define SATURATE_INT_PCM_LEFT_SHIFT(src, scale) \
  SATURATE_LEFT_SHIFT(src, scale, SAMPLE_BITS)
#endif

#endif /* #ifndef SCALE_H */
