/******************* Library for basic calculation routines ********************

   Author(s):

   Description: fixed point intrinsics

*******************************************************************************/

#if !defined(CPLX_MUL_H)
#define CPLX_MUL_H

#include "common_fix.h"

#if defined(__arm__) || defined(_M_ARM)
#include "arm/cplx_mul_arm.h"

#elif defined(__GNUC__) && defined(__mips__) && __mips_isa_rev < 6
#include "mips/cplx_mul_mips.h"

#endif

#if !defined(FUNCTION_cplxMultDiv2_32x16X2)
#define FUNCTION_cplxMultDiv2_32x16X2

inline void cplxMultDiv2(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_DBL a_Re,
                         const FIXP_DBL a_Im, const FIXP_SGL b_Re,
                         const FIXP_SGL b_Im) {
  *c_Re = fMultDiv2(a_Re, b_Re) - fMultDiv2(a_Im, b_Im);
  *c_Im = fMultDiv2(a_Re, b_Im) + fMultDiv2(a_Im, b_Re);
}
#endif

#if !defined(FUNCTION_cplxMultDiv2_16x16X2)
#define FUNCTION_cplxMultDiv2_16x16X2

inline void cplxMultDiv2(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_SGL a_Re,
                         const FIXP_SGL a_Im, const FIXP_SGL b_Re,
                         const FIXP_SGL b_Im) {
  *c_Re = fMultDiv2(a_Re, b_Re) - fMultDiv2(a_Im, b_Im);
  *c_Im = fMultDiv2(a_Re, b_Im) + fMultDiv2(a_Im, b_Re);
}

inline void cplxMultDiv2(FIXP_SGL *c_Re, FIXP_SGL *c_Im, const FIXP_SGL a_Re,
                         const FIXP_SGL a_Im, const FIXP_SGL b_Re,
                         const FIXP_SGL b_Im) {
  *c_Re = FX_DBL2FX_SGL(fMultDiv2(a_Re, b_Re) - fMultDiv2(a_Im, b_Im));
  *c_Im = FX_DBL2FX_SGL(fMultDiv2(a_Re, b_Im) + fMultDiv2(a_Im, b_Re));
}
#endif

#if !defined(FUNCTION_cplxMultDiv2_32x16)
#define FUNCTION_cplxMultDiv2_32x16

inline void cplxMultDiv2(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_DBL a_Re,
                         const FIXP_DBL a_Im, const FIXP_SPK w) {
  cplxMultDiv2(c_Re, c_Im, a_Re, a_Im, w.v.re, w.v.im);
}
#endif

#if !defined(FUNCTION_cplxMultDiv2_16x16)
#define FUNCTION_cplxMultDiv2_16x16

inline void cplxMultDiv2(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_SGL a_Re,
                         const FIXP_SGL a_Im, const FIXP_SPK w) {
  cplxMultDiv2(c_Re, c_Im, a_Re, a_Im, w.v.re, w.v.im);
}

inline void cplxMultDiv2(FIXP_SGL *c_Re, FIXP_SGL *c_Im, const FIXP_SGL a_Re,
                         const FIXP_SGL a_Im, const FIXP_SPK w) {
  cplxMultDiv2(c_Re, c_Im, a_Re, a_Im, w.v.re, w.v.im);
}
#endif

#if !defined(FUNCTION_cplxMultSubDiv2_32x16X2)
#define FUNCTION_cplxMultSubDiv2_32x16X2

inline void cplxMultSubDiv2(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_DBL a_Re,
                            const FIXP_DBL a_Im, const FIXP_SGL b_Re,
                            const FIXP_SGL b_Im) {
  *c_Re -= fMultDiv2(a_Re, b_Re) - fMultDiv2(a_Im, b_Im);
  *c_Im -= fMultDiv2(a_Re, b_Im) + fMultDiv2(a_Im, b_Re);
}
#endif

#if !defined(FUNCTION_cplxMultDiv2_32x32X2)
#define FUNCTION_cplxMultDiv2_32x32X2

inline void cplxMultDiv2(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_DBL a_Re,
                         const FIXP_DBL a_Im, const FIXP_DBL b_Re,
                         const FIXP_DBL b_Im) {
  *c_Re = fMultDiv2(a_Re, b_Re) - fMultDiv2(a_Im, b_Im);
  *c_Im = fMultDiv2(a_Re, b_Im) + fMultDiv2(a_Im, b_Re);
}
#endif

#if !defined(FUNCTION_cplxMultDiv2_32x32)
#define FUNCTION_cplxMultDiv2_32x32

inline void cplxMultDiv2(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_DBL a_Re,
                         const FIXP_DBL a_Im, const FIXP_DPK w) {
  cplxMultDiv2(c_Re, c_Im, a_Re, a_Im, w.v.re, w.v.im);
}
#endif

#if !defined(FUNCTION_cplxMultSubDiv2_32x32X2)
#define FUNCTION_cplxMultSubDiv2_32x32X2

inline void cplxMultSubDiv2(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_DBL a_Re,
                            const FIXP_DBL a_Im, const FIXP_DBL b_Re,
                            const FIXP_DBL b_Im) {
  *c_Re -= fMultDiv2(a_Re, b_Re) - fMultDiv2(a_Im, b_Im);
  *c_Im -= fMultDiv2(a_Re, b_Im) + fMultDiv2(a_Im, b_Re);
}
#endif

#if !defined(FUNCTION_cplxMult_32x16X2)
#define FUNCTION_cplxMult_32x16X2

inline void cplxMult(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_DBL a_Re,
                     const FIXP_DBL a_Im, const FIXP_SGL b_Re,
                     const FIXP_SGL b_Im) {
  *c_Re = fMult(a_Re, b_Re) - fMult(a_Im, b_Im);
  *c_Im = fMult(a_Re, b_Im) + fMult(a_Im, b_Re);
}
inline void cplxMult(FIXP_SGL *c_Re, FIXP_SGL *c_Im, const FIXP_SGL a_Re,
                     const FIXP_SGL a_Im, const FIXP_SGL b_Re,
                     const FIXP_SGL b_Im) {
  *c_Re = FX_DBL2FX_SGL(fMult(a_Re, b_Re) - fMult(a_Im, b_Im));
  *c_Im = FX_DBL2FX_SGL(fMult(a_Re, b_Im) + fMult(a_Im, b_Re));
}
#endif

#if !defined(FUNCTION_cplxMult_32x16)
#define FUNCTION_cplxMult_32x16

inline void cplxMult(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_DBL a_Re,
                     const FIXP_DBL a_Im, const FIXP_SPK w) {
  cplxMult(c_Re, c_Im, a_Re, a_Im, w.v.re, w.v.im);
}
#endif

#if !defined(FUNCTION_cplxMult_32x32X2)
#define FUNCTION_cplxMult_32x32X2

inline void cplxMult(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_DBL a_Re,
                     const FIXP_DBL a_Im, const FIXP_DBL b_Re,
                     const FIXP_DBL b_Im) {
  *c_Re = fMult(a_Re, b_Re) - fMult(a_Im, b_Im);
  *c_Im = fMult(a_Re, b_Im) + fMult(a_Im, b_Re);
}
#endif

#if !defined(FUNCTION_cplxMult_32x32)
#define FUNCTION_cplxMult_32x32
inline void cplxMult(FIXP_DBL *c_Re, FIXP_DBL *c_Im, const FIXP_DBL a_Re,
                     const FIXP_DBL a_Im, const FIXP_DPK w) {
  cplxMult(c_Re, c_Im, a_Re, a_Im, w.v.re, w.v.im);
}
#endif

#endif /* CPLX_MUL_H */
