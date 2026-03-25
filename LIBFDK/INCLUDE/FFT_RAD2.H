/******************* Library for basic calculation routines ********************

   Author(s):   M. Lohwasser, M. Gayer

   Description:

*******************************************************************************/

#ifndef FFT_RAD2_H
#define FFT_RAD2_H

#include "common_fix.h"

void dit_fft(FIXP_DBL *x, const INT ldn, const FIXP_STP *trigdata,
             const INT trigDataSize);

#endif /* FFT_RAD2_H */
