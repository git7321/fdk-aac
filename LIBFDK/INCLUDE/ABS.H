/******************* Library for basic calculation routines ********************

   Author(s):   M. Lohwasser

   Description: fixed point abs definitions

*******************************************************************************/

#if !defined(ABS_H)
#define ABS_H

#if defined(__mips__)
#include "mips/abs_mips.h"

#elif defined(__x86__)
#include "x86/abs_x86.h"

#endif /* all cores */

#if !defined(FUNCTION_fixabs_D)
inline FIXP_DBL fixabs_D(FIXP_DBL x) {
  return ((x) > (FIXP_DBL)(0)) ? (x) : -(x);
}
#endif

#if !defined(FUNCTION_fixabs_I)
inline INT fixabs_I(INT x) { return ((x) > (INT)(0)) ? (x) : -(x); }
#endif

#if !defined(FUNCTION_fixabs_S)
inline FIXP_SGL fixabs_S(FIXP_SGL x) {
  return ((x) > (FIXP_SGL)(0)) ? (x) : -(x);
}
#endif

#endif /* ABS_H */
