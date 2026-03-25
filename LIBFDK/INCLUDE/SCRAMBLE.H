/******************* Library for basic calculation routines ********************

   Author(s):

   Description:

*******************************************************************************/

#ifndef SCRAMBLE_H
#define SCRAMBLE_H

#include "common_fix.h"

#if defined(__arm__)
#include "arm/scramble_arm.h"

#elif defined(__mips__) && defined(__mips_dsp) && !defined(__mips16)
#include "mips/scramble_mips.h"

#endif

#if !defined(FUNCTION_scramble)

inline void scramble(FIXP_DBL *x, INT length) {
  INT m, k, j;
  FDK_ASSERT(!(((INT)(INT64)x) & (ALIGNMENT_DEFAULT - 1)));
  C_ALLOC_ALIGNED_CHECK(x);

  for (m = 1, j = 0; m < length - 1; m++) {
    {
      for (k = length >> 1; (!((j ^= k) & k)); k >>= 1)
        ;
    }

    if (j > m) {
      FIXP_DBL tmp;
      tmp = x[2 * m];
      x[2 * m] = x[2 * j];
      x[2 * j] = tmp;

      tmp = x[2 * m + 1];
      x[2 * m + 1] = x[2 * j + 1];
      x[2 * j + 1] = tmp;
    }
  }
}
#endif /* !defined(FUNCTION_scramble) */

#endif /* SCRAMBLE_H */
