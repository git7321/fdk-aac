/******************* Library for basic calculation routines ********************

   Author(s):   Marc Gayer

   Description: fixed point intrinsics

*******************************************************************************/

#if !defined(CLZ_H)
#define CLZ_H

#include "FDK_archdef.h"
#include "machine_type.h"

#if defined(__arm__)
#include "arm/clz_arm.h"

#elif defined(__mips__)
#include "mips/clz_mips.h"

#elif defined(__x86__)
#include "x86/clz_x86.h"

#elif defined(__powerpc__)
#include "ppc/clz_ppc.h"

#endif /* all cores */

#if !defined(FUNCTION_fixnormz_S)
#ifdef FUNCTION_fixnormz_D
inline INT fixnormz_S(SHORT a) {
  if (a < 0) {
    return 0;
  }
  return fixnormz_D((INT)(a)) - 16;
}
#else
inline INT fixnormz_S(SHORT a) {
  int leadingBits = 0;
  a = ~a;
  while (a & 0x8000) {
    leadingBits++;
    a <<= 1;
  }

  return (leadingBits);
}
#endif
#endif

#if !defined(FUNCTION_fixnormz_D)
inline INT fixnormz_D(LONG a) {
  INT leadingBits = 0;
  a = ~a;
  while (a & 0x80000000) {
    leadingBits++;
    a <<= 1;
  }

  return (leadingBits);
}
#endif

#if !defined(FUNCTION_fixnorm_S)
#ifdef FUNCTION_fixnorm_D
inline INT fixnorm_S(FIXP_SGL val) {
  if (val == (FIXP_SGL)0) {
    return 0;
  }
  return fixnorm_D((INT)(val)) - 16;
}
#else
inline INT fixnorm_S(FIXP_SGL val) {
  INT leadingBits = 0;
  if (val != (FIXP_SGL)0) {
    if (val < (FIXP_SGL)0) {
      val = ~val;
    }
    leadingBits = fixnormz_S(val) - 1;
  }
  return (leadingBits);
}
#endif
#endif

#if !defined(FUNCTION_fixnorm_D)
inline INT fixnorm_D(FIXP_DBL val) {
  INT leadingBits = 0;
  if (val != (FIXP_DBL)0) {
    if (val < (FIXP_DBL)0) {
      val = ~val;
    }
    leadingBits = fixnormz_D(val) - 1;
  }
  return (leadingBits);
}
#endif

#endif /* CLZ_H */
