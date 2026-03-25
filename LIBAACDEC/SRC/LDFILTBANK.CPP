/**************************** AAC decoder library ******************************

   Author(s):

   Description: low delay filterbank

*******************************************************************************/

#include "ldfiltbank.h"

#include "aac_rom.h"
#include "dct.h"
#include "FDK_tools_rom.h"
#include "mdct.h"

#define LDFB_HEADROOM 2

#if defined(__arm__)
#endif

static void multE2_DinvF_fdk(PCM_DEC *output, FIXP_DBL *x, const FIXP_WTB *fb,
                             FIXP_DBL *z, const int N) {
  int i;

  const int scale = (DFRACT_BITS - PCM_OUT_BITS) - LDFB_HEADROOM + (3);

#if ((DFRACT_BITS - PCM_OUT_BITS - LDFB_HEADROOM + (3) - 1) > 0)
  FIXP_DBL rnd_val_wts0 = (FIXP_DBL)0;
  FIXP_DBL rnd_val_wts1 = (FIXP_DBL)0;
#if ((DFRACT_BITS - PCM_OUT_BITS - LDFB_HEADROOM + (3) - WTS0 - 1) > 0)
  if (-WTS0 - 1 + scale)
    rnd_val_wts0 = (FIXP_DBL)(1 << (-WTS0 - 1 + scale - 1));
#endif
  if (-WTS1 - 1 + scale)
    rnd_val_wts1 = (FIXP_DBL)(1 << (-WTS1 - 1 + scale - 1));
#endif

  for (i = 0; i < N / 4; i++) {
    FIXP_DBL z0, z2, tmp;

    z2 = x[N / 2 + i];
    z0 = fAddSaturate(z2,
                      (fMultDiv2(z[N / 2 + i], fb[2 * N + i]) >> (-WTS2 - 1)));

    z[N / 2 + i] = fAddSaturate(
        x[N / 2 - 1 - i],
        (fMultDiv2(z[N + i], fb[2 * N + N / 2 + i]) >> (-WTS2 - 1)));

    tmp = (fMultDiv2(z[N / 2 + i], fb[N + N / 2 - 1 - i]) +
           fMultDiv2(z[i], fb[N + N / 2 + i]));

#if ((DFRACT_BITS - PCM_OUT_BITS - LDFB_HEADROOM + (3) - 1) > 0)
    FDK_ASSERT((-WTS1 - 1 + scale) >= 0);
    FDK_ASSERT(tmp <= ((FIXP_DBL)0x7FFFFFFF -
                       rnd_val_wts1));
    output[(N * 3 / 4 - 1 - i)] = (PCM_DEC)SATURATE_RIGHT_SHIFT(
        tmp + rnd_val_wts1, -WTS1 - 1 + scale, PCM_OUT_BITS);
#else
    FDK_ASSERT((WTS1 + 1 - scale) >= 0);
    output[(N * 3 / 4 - 1 - i)] =
        (PCM_DEC)SATURATE_LEFT_SHIFT(tmp, WTS1 + 1 - scale, PCM_OUT_BITS);
#endif

    z[i] = z0;
    z[N + i] = z2;
  }

  for (i = N / 4; i < N / 2; i++) {
    FIXP_DBL z0, z2, tmp0, tmp1;

    z2 = x[N / 2 + i];
    z0 = fAddSaturate(z2,
                      (fMultDiv2(z[N / 2 + i], fb[2 * N + i]) >> (-WTS2 - 1)));

    z[N / 2 + i] = fAddSaturate(
        x[N / 2 - 1 - i],
        (fMultDiv2(z[N + i], fb[2 * N + N / 2 + i]) >> (-WTS2 - 1)));

    tmp0 = (fMultDiv2(z[N / 2 + i], fb[N / 2 - 1 - i]) +
            fMultDiv2(z[i], fb[N / 2 + i]));
    tmp1 = (fMultDiv2(z[N / 2 + i], fb[N + N / 2 - 1 - i]) +
            fMultDiv2(z[i], fb[N + N / 2 + i]));

#if ((DFRACT_BITS - PCM_OUT_BITS - LDFB_HEADROOM + (3) - 1) > 0)
    FDK_ASSERT((-WTS0 - 1 + scale) >= 0);
    FDK_ASSERT(tmp0 <= ((FIXP_DBL)0x7FFFFFFF -
                        rnd_val_wts0));
    FDK_ASSERT(tmp1 <= ((FIXP_DBL)0x7FFFFFFF -
                        rnd_val_wts1));
    output[(i - N / 4)] = (PCM_DEC)SATURATE_RIGHT_SHIFT(
        tmp0 + rnd_val_wts0, -WTS0 - 1 + scale, PCM_OUT_BITS);
    output[(N * 3 / 4 - 1 - i)] = (PCM_DEC)SATURATE_RIGHT_SHIFT(
        tmp1 + rnd_val_wts1, -WTS1 - 1 + scale, PCM_OUT_BITS);
#else
    FDK_ASSERT((WTS0 + 1 - scale) >= 0);
    output[(i - N / 4)] =
        (PCM_DEC)SATURATE_LEFT_SHIFT(tmp0, WTS0 + 1 - scale, PCM_OUT_BITS);
    output[(N * 3 / 4 - 1 - i)] =
        (PCM_DEC)SATURATE_LEFT_SHIFT(tmp1, WTS1 + 1 - scale, PCM_OUT_BITS);
#endif
    z[i] = z0;
    z[N + i] = z2;
  }

  for (i = 0; i < N / 4; i++) {
    FIXP_DBL tmp0 = fMultDiv2(z[i], fb[N / 2 + i]);

#if ((DFRACT_BITS - PCM_OUT_BITS - LDFB_HEADROOM + (3) - 1) > 0)
    FDK_ASSERT((-WTS0 - 1 + scale) >= 0);
    FDK_ASSERT(tmp0 <= ((FIXP_DBL)0x7FFFFFFF -
                        rnd_val_wts0));
    output[(N * 3 / 4 + i)] = (PCM_DEC)SATURATE_RIGHT_SHIFT(
        tmp0 + rnd_val_wts0, -WTS0 - 1 + scale, PCM_OUT_BITS);
#else
    FDK_ASSERT((WTS0 + 1 - scale) >= 0);
    output[(N * 3 / 4 + i)] =
        (PCM_DEC)SATURATE_LEFT_SHIFT(tmp0, WTS0 + 1 - scale, PCM_OUT_BITS);
#endif
  }
}

int InvMdctTransformLowDelay_fdk(FIXP_DBL *mdctData, const int mdctData_e,
                                 PCM_DEC *output, FIXP_DBL *fs_buffer,
                                 const int N) {
  const FIXP_WTB *coef;
  FIXP_DBL gain = (FIXP_DBL)0;
  int scale = mdctData_e + MDCT_OUT_HEADROOM -
              LDFB_HEADROOM;
  int i;

  switch (N) {
    case 256:
      coef = LowDelaySynthesis256;
      break;
    case 240:
      coef = LowDelaySynthesis240;
      break;
    case 160:
      coef = LowDelaySynthesis160;
      break;
    case 128:
      coef = LowDelaySynthesis128;
      break;
    case 120:
      coef = LowDelaySynthesis120;
      break;
    case 512:
      coef = LowDelaySynthesis512;
      break;
    case 480:
    default:
      coef = LowDelaySynthesis480;
      break;
  }

  imdct_gain(&gain, &scale, N);
  dct_IV(mdctData, N, &scale);

  if (N == 256 || N == 240 || N == 160) {
    scale -= 1;
  } else if (N == 128 || N == 120) {
    scale -= 2;
  }

  if (gain != (FIXP_DBL)0) {
    for (i = 0; i < N; i++) {
      mdctData[i] = fMult(mdctData[i], gain);
    }
  }
  scaleValuesSaturate(mdctData, N, scale);

  multE2_DinvF_fdk(output, mdctData, coef, fs_buffer, N);

  return (1);
}
