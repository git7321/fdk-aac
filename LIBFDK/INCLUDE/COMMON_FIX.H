#ifndef COMMON_FIX_H
#define COMMON_FIX_H

#include "FDK_archdef.h"
#include "machine_type.h"

#define FRACT_BITS 16
#define DFRACT_BITS 32
#define ACCU_BITS 40

#if defined(SAMPLE_BITS)
#if (SAMPLE_BITS == DFRACT_BITS)
#define FIXP_PCM FIXP_DBL
#define MAXVAL_FIXP_PCM MAXVAL_DBL
#define MINVAL_FIXP_PCM MINVAL_DBL
#define FX_PCM2FX_DBL(x) ((FIXP_DBL)(x))
#define FX_DBL2FX_PCM(x) ((INT_PCM)(x))
#elif (SAMPLE_BITS == FRACT_BITS)
#define FIXP_PCM FIXP_SGL
#define MAXVAL_FIXP_PCM MAXVAL_SGL
#define MINVAL_FIXP_PCM MINVAL_SGL
#define FX_PCM2FX_DBL(x) FX_SGL2FX_DBL((FIXP_SGL)(x))
#define FX_DBL2FX_PCM(x) FX_DBL2FX_SGL(x)
#else
#error SAMPLE_BITS different from FRACT_BITS or DFRACT_BITS not implemented!
#endif
#endif

#define SGL_MASK ((1UL << FRACT_BITS) - 1)

#define MAX_SHIFT_SGL \
  (FRACT_BITS - 1)
#define MAX_SHIFT_DBL \
  (DFRACT_BITS - 1)

#define FRACT_FIX_SCALE ((INT64(1) << (FRACT_BITS - 1)))
#define DFRACT_FIX_SCALE ((INT64(1) << (DFRACT_BITS - 1)))

#define MAXVAL_SGL \
  ((signed)0x00007FFF)
#define MINVAL_SGL \
  ((signed)0xFFFF8000)
#define MAXVAL_DBL \
  ((signed)0x7FFFFFFF)
#define MINVAL_DBL \
  ((signed)0x80000000)

#define FX_DBL2FXCONST_SGL(val)                                               \
  ((((((val) >> (DFRACT_BITS - FRACT_BITS - 1)) + 1) >                        \
     (((LONG)1 << FRACT_BITS) - 1)) &&                                        \
    ((LONG)(val) > 0))                                                        \
       ? (FIXP_SGL)(SHORT)(((LONG)1 << (FRACT_BITS - 1)) - 1)                 \
       : (FIXP_SGL)(SHORT)((((val) >> (DFRACT_BITS - FRACT_BITS - 1)) + 1) >> \
                           1))

#define shouldBeUnion union

typedef SHORT FIXP_SGL;
typedef LONG FIXP_DBL;

#define FL2FXCONST_SPC FL2FXCONST_DBL

#define MINVAL_DBL_CONST MINVAL_DBL
#define MINVAL_SGL_CONST MINVAL_SGL

#define FL2FXCONST_SGL(val)                                                  \
  (FIXP_SGL)(                                                                \
      ((val) >= 0)                                                           \
          ? ((((double)(val) * (FRACT_FIX_SCALE) + 0.5) >=                   \
              (double)(MAXVAL_SGL))                                          \
                 ? (SHORT)(MAXVAL_SGL)                                       \
                 : (SHORT)((double)(val) * (double)(FRACT_FIX_SCALE) + 0.5)) \
          : ((((double)(val) * (FRACT_FIX_SCALE)-0.5) <=                     \
              (double)(MINVAL_SGL_CONST))                                    \
                 ? (SHORT)(MINVAL_SGL_CONST)                                 \
                 : (SHORT)((double)(val) * (double)(FRACT_FIX_SCALE)-0.5)))

#define FL2FXCONST_DBL(val)                                                  \
  (FIXP_DBL)(                                                                \
      ((val) >= 0)                                                           \
          ? ((((double)(val) * (DFRACT_FIX_SCALE) + 0.5) >=                  \
              (double)(MAXVAL_DBL))                                          \
                 ? (LONG)(MAXVAL_DBL)                                        \
                 : (LONG)((double)(val) * (double)(DFRACT_FIX_SCALE) + 0.5)) \
          : ((((double)(val) * (DFRACT_FIX_SCALE)-0.5) <=                    \
              (double)(MINVAL_DBL_CONST))                                    \
                 ? (LONG)(MINVAL_DBL_CONST)                                  \
                 : (LONG)((double)(val) * (double)(DFRACT_FIX_SCALE)-0.5)))

#define FL2FX_SPC FL2FX_DBL
#define FL2FX_SGL(val)                                             \
  ((val) > 0.0f ? (SHORT)((val) * (float)(FRACT_FIX_SCALE) + 0.5f) \
                : (SHORT)((val) * (float)(FRACT_FIX_SCALE)-0.5f))
#define FL2FX_DBL(val)                                             \
  ((val) > 0.0f ? (LONG)((val) * (float)(DFRACT_FIX_SCALE) + 0.5f) \
                : (LONG)((val) * (float)(DFRACT_FIX_SCALE)-0.5f))

#define FX_ACC2FX_SGL(val) ((FIXP_SGL)((val) >> (ACCU_BITS - FRACT_BITS)))
#define FX_ACC2FX_DBL(val) ((FIXP_DBL)((val) >> (ACCU_BITS - DFRACT_BITS)))
#define FX_SGL2FX_ACC(val) ((FIXP_ACC)((LONG)(val) << (ACCU_BITS - FRACT_BITS)))
#define FX_SGL2FX_DBL(val) \
  ((FIXP_DBL)((LONG)(val) << (DFRACT_BITS - FRACT_BITS)))
#define FX_DBL2FX_SGL(val) ((FIXP_SGL)((val) >> (DFRACT_BITS - FRACT_BITS)))

#define FX_DBL2FL(val) \
  ((float)((double)(val) / (double)DFRACT_FIX_SCALE))
#define FX_DBL2DOUBLE(val) (((double)(val) / (double)DFRACT_FIX_SCALE))

#include "fixmul.h"

FDK_INLINE LONG fMult(SHORT a, SHORT b) { return fixmul_SS(a, b); }
FDK_INLINE LONG fMult(SHORT a, LONG b) { return fixmul_SD(a, b); }
FDK_INLINE LONG fMult(LONG a, SHORT b) { return fixmul_DS(a, b); }
FDK_INLINE LONG fMult(LONG a, LONG b) { return fixmul_DD(a, b); }
FDK_INLINE LONG fPow2(LONG a) { return fixpow2_D(a); }
FDK_INLINE LONG fPow2(SHORT a) { return fixpow2_S(a); }

FDK_INLINE LONG fMultDiv2(SHORT a, SHORT b) { return fixmuldiv2_SS(a, b); }
FDK_INLINE LONG fMultDiv2(SHORT a, LONG b) { return fixmuldiv2_SD(a, b); }
FDK_INLINE LONG fMultDiv2(LONG a, SHORT b) { return fixmuldiv2_DS(a, b); }
FDK_INLINE LONG fMultDiv2(LONG a, LONG b) { return fixmuldiv2_DD(a, b); }
FDK_INLINE LONG fPow2Div2(LONG a) { return fixpow2div2_D(a); }
FDK_INLINE LONG fPow2Div2(SHORT a) { return fixpow2div2_S(a); }

FDK_INLINE LONG fMultDiv2BitExact(LONG a, LONG b) {
  return fixmuldiv2BitExact_DD(a, b);
}
FDK_INLINE LONG fMultDiv2BitExact(SHORT a, LONG b) {
  return fixmuldiv2BitExact_SD(a, b);
}
FDK_INLINE LONG fMultDiv2BitExact(LONG a, SHORT b) {
  return fixmuldiv2BitExact_DS(a, b);
}
FDK_INLINE LONG fMultBitExact(LONG a, LONG b) {
  return fixmulBitExact_DD(a, b);
}
FDK_INLINE LONG fMultBitExact(SHORT a, LONG b) {
  return fixmulBitExact_SD(a, b);
}
FDK_INLINE LONG fMultBitExact(LONG a, SHORT b) {
  return fixmulBitExact_DS(a, b);
}

#include "abs.h"

FDK_INLINE FIXP_DBL fAbs(FIXP_DBL x) { return fixabs_D(x); }
FDK_INLINE FIXP_SGL fAbs(FIXP_SGL x) { return fixabs_S(x); }

#if !defined(__LP64__)
FDK_INLINE INT fAbs(INT x) { return fixabs_I(x); }
#endif

#include "clz.h"

FDK_INLINE INT fNormz(INT64 x) {
  INT clz = fixnormz_D((INT)(x >> 32));
  if (clz == 32) clz += fixnormz_D((INT)x);
  return clz;
}
FDK_INLINE INT fNormz(FIXP_DBL x) { return fixnormz_D(x); }
FDK_INLINE INT fNormz(FIXP_SGL x) { return fixnormz_S(x); }
FDK_INLINE INT fNorm(FIXP_DBL x) { return fixnorm_D(x); }
FDK_INLINE INT fNorm(FIXP_SGL x) { return fixnorm_S(x); }

#include "clz.h"
#define fixp_abs(x) fAbs(x)
#define fixMin(a, b) fMin(a, b)
#define fixMax(a, b) fMax(a, b)
#define CntLeadingZeros(x) fixnormz_D(x)
#define CountLeadingBits(x) fixnorm_D(x)

#include "fixmadd.h"

FDK_INLINE FIXP_DBL fMultAddDiv2(FIXP_DBL x, FIXP_DBL a, FIXP_DBL b) {
  return fixmadddiv2_DD(x, a, b);
}
FDK_INLINE FIXP_DBL fMultAddDiv2(FIXP_DBL x, FIXP_SGL a, FIXP_DBL b) {
  return fixmadddiv2_SD(x, a, b);
}
FDK_INLINE FIXP_DBL fMultAddDiv2(FIXP_DBL x, FIXP_DBL a, FIXP_SGL b) {
  return fixmadddiv2_DS(x, a, b);
}
FDK_INLINE FIXP_DBL fMultAddDiv2(FIXP_DBL x, FIXP_SGL a, FIXP_SGL b) {
  return fixmadddiv2_SS(x, a, b);
}

FDK_INLINE FIXP_DBL fPow2AddDiv2(FIXP_DBL x, FIXP_DBL a) {
  return fixpadddiv2_D(x, a);
}
FDK_INLINE FIXP_DBL fPow2AddDiv2(FIXP_DBL x, FIXP_SGL a) {
  return fixpadddiv2_S(x, a);
}

FDK_INLINE FIXP_DBL fMultAdd(FIXP_DBL x, FIXP_DBL a, FIXP_DBL b) {
  return fixmadd_DD(x, a, b);
}
inline FIXP_DBL fMultAdd(FIXP_DBL x, FIXP_SGL a, FIXP_DBL b) {
  return fixmadd_SD(x, a, b);
}
inline FIXP_DBL fMultAdd(FIXP_DBL x, FIXP_DBL a, FIXP_SGL b) {
  return fixmadd_DS(x, a, b);
}
inline FIXP_DBL fMultAdd(FIXP_DBL x, FIXP_SGL a, FIXP_SGL b) {
  return fixmadd_SS(x, a, b);
}

inline FIXP_DBL fPow2Add(FIXP_DBL x, FIXP_DBL a) { return fixpadd_D(x, a); }
inline FIXP_DBL fPow2Add(FIXP_DBL x, FIXP_SGL a) { return fixpadd_S(x, a); }

inline FIXP_DBL fMultSubDiv2(FIXP_DBL x, FIXP_DBL a, FIXP_DBL b) {
  return fixmsubdiv2_DD(x, a, b);
}
inline FIXP_DBL fMultSubDiv2(FIXP_DBL x, FIXP_SGL a, FIXP_DBL b) {
  return fixmsubdiv2_SD(x, a, b);
}
inline FIXP_DBL fMultSubDiv2(FIXP_DBL x, FIXP_DBL a, FIXP_SGL b) {
  return fixmsubdiv2_DS(x, a, b);
}
inline FIXP_DBL fMultSubDiv2(FIXP_DBL x, FIXP_SGL a, FIXP_SGL b) {
  return fixmsubdiv2_SS(x, a, b);
}

FDK_INLINE FIXP_DBL fMultSub(FIXP_DBL x, FIXP_DBL a, FIXP_DBL b) {
  return fixmsub_DD(x, a, b);
}
inline FIXP_DBL fMultSub(FIXP_DBL x, FIXP_SGL a, FIXP_DBL b) {
  return fixmsub_SD(x, a, b);
}
inline FIXP_DBL fMultSub(FIXP_DBL x, FIXP_DBL a, FIXP_SGL b) {
  return fixmsub_DS(x, a, b);
}
inline FIXP_DBL fMultSub(FIXP_DBL x, FIXP_SGL a, FIXP_SGL b) {
  return fixmsub_SS(x, a, b);
}

FDK_INLINE FIXP_DBL fMultAddDiv2BitExact(FIXP_DBL x, FIXP_DBL a, FIXP_DBL b) {
  return fixmadddiv2BitExact_DD(x, a, b);
}
FDK_INLINE FIXP_DBL fMultAddDiv2BitExact(FIXP_DBL x, FIXP_SGL a, FIXP_DBL b) {
  return fixmadddiv2BitExact_SD(x, a, b);
}
FDK_INLINE FIXP_DBL fMultAddDiv2BitExact(FIXP_DBL x, FIXP_DBL a, FIXP_SGL b) {
  return fixmadddiv2BitExact_DS(x, a, b);
}
FDK_INLINE FIXP_DBL fMultSubDiv2BitExact(FIXP_DBL x, FIXP_DBL a, FIXP_DBL b) {
  return fixmsubdiv2BitExact_DD(x, a, b);
}
FDK_INLINE FIXP_DBL fMultSubDiv2BitExact(FIXP_DBL x, FIXP_SGL a, FIXP_DBL b) {
  return fixmsubdiv2BitExact_SD(x, a, b);
}
FDK_INLINE FIXP_DBL fMultSubDiv2BitExact(FIXP_DBL x, FIXP_DBL a, FIXP_SGL b) {
  return fixmsubdiv2BitExact_DS(x, a, b);
}

#include "fixminmax.h"

FDK_INLINE FIXP_DBL fMin(FIXP_DBL a, FIXP_DBL b) { return fixmin_D(a, b); }
FDK_INLINE FIXP_DBL fMax(FIXP_DBL a, FIXP_DBL b) { return fixmax_D(a, b); }

FDK_INLINE FIXP_SGL fMin(FIXP_SGL a, FIXP_SGL b) { return fixmin_S(a, b); }
FDK_INLINE FIXP_SGL fMax(FIXP_SGL a, FIXP_SGL b) { return fixmax_S(a, b); }

#if !defined(__LP64__)
FDK_INLINE INT fMax(INT a, INT b) { return fixmax_I(a, b); }
FDK_INLINE INT fMin(INT a, INT b) { return fixmin_I(a, b); }
#endif

inline UINT fMax(UINT a, UINT b) { return fixmax_UI(a, b); }
inline UINT fMin(UINT a, UINT b) { return fixmin_UI(a, b); }

inline UCHAR fMax(UCHAR a, UCHAR b) {
  return (UCHAR)fixmax_UI((UINT)a, (UINT)b);
}
inline UCHAR fMin(UCHAR a, UCHAR b) {
  return (UCHAR)fixmin_UI((UINT)a, (UINT)b);
}

typedef shouldBeUnion {
  struct {
    FIXP_SGL re;
    FIXP_SGL im;
  } v;
  LONG w;
}
FIXP_SPK;

typedef shouldBeUnion {
  struct {
    FIXP_DBL re;
    FIXP_DBL im;
  } v;
  INT64 w;
}
FIXP_DPK;

#include "fixmul.h"
#include "fixmadd.h"
#include "cplx_mul.h"
#include "fixpoint_math.h"

#endif
