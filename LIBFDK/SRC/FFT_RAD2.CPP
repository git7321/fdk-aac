/******************* Library for basic calculation routines ********************

   Author(s):   M. Lohwasser, M. Gayer

   Description:

*******************************************************************************/

#include "fft_rad2.h"

#include "scramble.h"

#define __FFT_RAD2_CPP__

#if defined(__arm__)
#include "arm/fft_rad2_arm.cpp"

#elif defined(__GNUC__) && defined(__mips__) && defined(__mips_dsp) && !defined(__mips16)
#include "mips/fft_rad2_mips.cpp"

#endif

#ifndef FUNCTION_dit_fft

void dit_fft(FIXP_DBL *x, const INT ldn, const FIXP_STP *trigdata,
             const INT trigDataSize) {
  const INT n = 1 << ldn;
  INT trigstep, i, ldm;

  C_ALLOC_ALIGNED_CHECK(x);

  scramble(x, n);

  for (i = 0; i < n * 2; i += 8) {
    FIXP_DBL a00, a10, a20, a30;
    a00 = (x[i + 0] + x[i + 2]) >> 1; /* Re A + Re B */
    a10 = (x[i + 4] + x[i + 6]) >> 1; /* Re C + Re D */
    a20 = (x[i + 1] + x[i + 3]) >> 1; /* Im A + Im B */
    a30 = (x[i + 5] + x[i + 7]) >> 1; /* Im C + Im D */

    x[i + 0] = a00 + a10; /* Re A' = Re A + Re B + Re C + Re D */
    x[i + 4] = a00 - a10; /* Re C' = Re A + Re B - Re C - Re D */
    x[i + 1] = a20 + a30; /* Im A' = Im A + Im B + Im C + Im D */
    x[i + 5] = a20 - a30; /* Im C' = Im A + Im B - Im C - Im D */

    a00 = a00 - x[i + 2]; /* Re A - Re B */
    a10 = a10 - x[i + 6]; /* Re C - Re D */
    a20 = a20 - x[i + 3]; /* Im A - Im B */
    a30 = a30 - x[i + 7]; /* Im C - Im D */

    x[i + 2] = a00 + a30; /* Re B' = Re A - Re B + Im C - Im D */
    x[i + 6] = a00 - a30; /* Re D' = Re A - Re B - Im C + Im D */
    x[i + 3] = a20 - a10; /* Im B' = Im A - Im B - Re C + Re D */
    x[i + 7] = a20 + a10; /* Im D' = Im A - Im B + Re C - Re D */
  }

  for (ldm = 3; ldm <= ldn; ++ldm) {
    INT m = (1 << ldm);
    INT mh = (m >> 1);
    INT j, r;

    trigstep = ((trigDataSize << 2) >> ldm);

    FDK_ASSERT(trigstep > 0);

    {

      j = 0;

      for (r = 0; r < n; r += m) {
        INT t1 = (r + j) << 1;
        INT t2 = t1 + (mh << 1);
        FIXP_DBL vr, vi, ur, ui;

        vi = x[t2 + 1] >> 1;
        vr = x[t2] >> 1;

        ur = x[t1] >> 1;
        ui = x[t1 + 1] >> 1;

        x[t1] = ur + vr;
        x[t1 + 1] = ui + vi;

        x[t2] = ur - vr;
        x[t2 + 1] = ui - vi;

        t1 += mh;
        t2 = t1 + (mh << 1);

        vr = x[t2 + 1] >> 1;
        vi = x[t2] >> 1;

        ur = x[t1] >> 1;
        ui = x[t1 + 1] >> 1;

        x[t1] = ur + vr;
        x[t1 + 1] = ui - vi;

        x[t2] = ur - vr;
        x[t2 + 1] = ui + vi;
      }
    }

    for (j = 1; j < mh / 4; ++j) {
      FIXP_STP cs;

      cs = trigdata[j * trigstep];

      for (r = 0; r < n; r += m) {
        INT t1 = (r + j) << 1;
        INT t2 = t1 + (mh << 1);
        FIXP_DBL vr, vi, ur, ui;

        cplxMultDiv2(&vi, &vr, x[t2 + 1], x[t2], cs);

        ur = x[t1] >> 1;
        ui = x[t1 + 1] >> 1;

        x[t1] = ur + vr;
        x[t1 + 1] = ui + vi;

        x[t2] = ur - vr;
        x[t2 + 1] = ui - vi;

        t1 += mh;
        t2 = t1 + (mh << 1);

        cplxMultDiv2(&vr, &vi, x[t2 + 1], x[t2], cs);

        ur = x[t1] >> 1;
        ui = x[t1 + 1] >> 1;

        x[t1] = ur + vr;
        x[t1 + 1] = ui - vi;

        x[t2] = ur - vr;
        x[t2 + 1] = ui + vi;

        t1 = (r + mh / 2 - j) << 1;
        t2 = t1 + (mh << 1);

        cplxMultDiv2(&vi, &vr, x[t2], x[t2 + 1], cs);

        ur = x[t1] >> 1;
        ui = x[t1 + 1] >> 1;

        x[t1] = ur + vr;
        x[t1 + 1] = ui - vi;

        x[t2] = ur - vr;
        x[t2 + 1] = ui + vi;

        t1 += mh;
        t2 = t1 + (mh << 1);

        cplxMultDiv2(&vr, &vi, x[t2], x[t2 + 1], cs);

        ur = x[t1] >> 1;
        ui = x[t1 + 1] >> 1;

        x[t1] = ur - vr;
        x[t1 + 1] = ui - vi;

        x[t2] = ur + vr;
        x[t2 + 1] = ui + vi;
      }
    }

    {
      j = mh / 4;

      for (r = 0; r < n; r += m) {
        INT t1 = (r + j) << 1;
        INT t2 = t1 + (mh << 1);
        FIXP_DBL vr, vi, ur, ui;

        cplxMultDiv2(&vi, &vr, x[t2 + 1], x[t2], STC(0x5a82799a),
                     STC(0x5a82799a));

        ur = x[t1] >> 1;
        ui = x[t1 + 1] >> 1;

        x[t1] = ur + vr;
        x[t1 + 1] = ui + vi;

        x[t2] = ur - vr;
        x[t2 + 1] = ui - vi;

        t1 += mh;
        t2 = t1 + (mh << 1);

        cplxMultDiv2(&vr, &vi, x[t2 + 1], x[t2], STC(0x5a82799a),
                     STC(0x5a82799a));

        ur = x[t1] >> 1;
        ui = x[t1 + 1] >> 1;

        x[t1] = ur + vr;
        x[t1 + 1] = ui - vi;

        x[t2] = ur - vr;
        x[t2 + 1] = ui + vi;
      }
    }
  }
}

#endif
