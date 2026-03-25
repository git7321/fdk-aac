/**************************** AAC decoder library ******************************

   Author(s):   Matthias Hildenbrand

   Description: USAC ACELP LTP filter

*******************************************************************************/

#include "usacdec_ace_ltp.h"

#include "genericStds.h"
#include "common_fix.h"

#define UP_SAMP 4
#define L_INTERPOL2 16
#define L_SUBFR 64

#define A2 FL2FX_SGL(2 * 0.18f)
#define B FL2FX_SGL(0.64f)

static const LONG Pred_lt4_inter4_2[UP_SAMP][L_INTERPOL2] = {
    {(LONG)0x0000FFFC, (LONG)0x0008FFFC, (LONG)0xFFEB004C, (LONG)0xFF50014A,
     (LONG)0xFDD90351, (LONG)0xFB2A06CD, (LONG)0xF6920D46, (LONG)0xEBB42B35,
     (LONG)0x6D9EEF39, (LONG)0x0618FE0F, (LONG)0xFFE00131, (LONG)0xFE5501C5,
     (LONG)0xFE5E015D, (LONG)0xFEF700B6, (LONG)0xFF920037, (LONG)0xFFEC0003},
    {(LONG)0x0002FFF2, (LONG)0x0026FFBD, (LONG)0x005DFF98, (LONG)0x0055FFEF,
     (LONG)0xFF89015F, (LONG)0xFD3A04E5, (LONG)0xF7D90DAA, (LONG)0xE67A50EE,
     (LONG)0x50EEE67A, (LONG)0x0DAAF7D9, (LONG)0x04E5FD3A, (LONG)0x015FFF89,
     (LONG)0xFFEF0055, (LONG)0xFF98005D, (LONG)0xFFBD0026, (LONG)0xFFF20002},
    {(LONG)0x0003FFEC, (LONG)0x0037FF92, (LONG)0x00B6FEF7, (LONG)0x015DFE5E,
     (LONG)0x01C5FE55, (LONG)0x0131FFE0, (LONG)0xFE0F0618, (LONG)0xEF396D9E,
     (LONG)0x2B35EBB4, (LONG)0x0D46F692, (LONG)0x06CDFB2A, (LONG)0x0351FDD9,
     (LONG)0x014AFF50, (LONG)0x004CFFEB, (LONG)0xFFFC0008, (LONG)0xFFFC0000},
    {(LONG)0x0002FFF2, (LONG)0x002BFF9E, (LONG)0x00B9FECE, (LONG)0x01CFFD75,
     (LONG)0x035EFBC1, (LONG)0x0521FA0C, (LONG)0x06AAF8C9, (LONG)0x07907852,
     (LONG)0x0790F8C9, (LONG)0x06AAFA0C, (LONG)0x0521FBC1, (LONG)0x035EFD75,
     (LONG)0x01CFFECE, (LONG)0x00B9FF9E, (LONG)0x002BFFF2, (LONG)0x00020000}};

void Pred_lt4(FIXP_DBL exc[],
              int T0,
              int frac
) {
  int j;
  FIXP_DBL *x;
  const LONG *interpol;
  FIXP_DBL L_sumb, L_sumt;

  x = &exc[-T0 - L_INTERPOL2 + 1];

  if (--frac < 0)
    frac += UP_SAMP;
  else
    x--;

  j = L_SUBFR + 1;
  do {
    LONG filt;
    FIXP_DBL x0, x1;
    FIXP_DBL *xi = x++;
    interpol = Pred_lt4_inter4_2[frac];
    int i = 3;

    filt = *interpol++;
    x0 = *xi++;
    x1 = *xi++;
    L_sumt = fMultDiv2(x0, (FIXP_SGL)((SHORT)(filt >> 16)));
    L_sumb = fMultDiv2(x1, (FIXP_SGL)((SHORT)filt));
    do {
      filt = *interpol++;
      x0 = *xi++;
      x1 = *xi++;
      L_sumt = fMultAddDiv2(L_sumt, x0, (FIXP_SGL)((SHORT)(filt >> 16)));
      L_sumb = fMultAddDiv2(L_sumb, x1, (FIXP_SGL)((SHORT)filt));

      filt = *interpol++;
      x0 = *xi++;
      x1 = *xi++;
      L_sumt = fMultAddDiv2(L_sumt, x0, (FIXP_SGL)((SHORT)(filt >> 16)));
      L_sumb = fMultAddDiv2(L_sumb, x1, (FIXP_SGL)((SHORT)filt));

      filt = *interpol++;
      x0 = *xi++;
      x1 = *xi++;
      L_sumt = fMultAddDiv2(L_sumt, x0, (FIXP_SGL)((SHORT)(filt >> 16)));
      L_sumb = fMultAddDiv2(L_sumb, x1, (FIXP_SGL)((SHORT)filt));

      filt = *interpol++;
      x0 = *xi++;
      x1 = *xi++;
      L_sumt = fMultAddDiv2(L_sumt, x0, (FIXP_SGL)((SHORT)(filt >> 16)));
      L_sumb = fMultAddDiv2(L_sumb, x1, (FIXP_SGL)((SHORT)filt));

      filt = *interpol++;
      x0 = *xi++;
      x1 = *xi++;
      L_sumt = fMultAddDiv2(L_sumt, x0, (FIXP_SGL)((SHORT)(filt >> 16)));
      L_sumb = fMultAddDiv2(L_sumb, x1, (FIXP_SGL)((SHORT)filt));
    } while (--i != 0);

    L_sumb <<= 1;
    L_sumb = fAddSaturate(L_sumt << 1, L_sumb);
    *exc++ = L_sumb;
  } while (--j != 0);
  return;
}

void Pred_lt4_postfilter(FIXP_DBL exc[]
) {
  int i;
  FIXP_DBL sum0, sum1, a_exc0, a_exc1;
  a_exc0 = fMultDiv2(A2, exc[-1]);
  a_exc1 = fMultDiv2(A2, exc[0]);

  for (i = 0; i < L_SUBFR; i += 2) {
    sum0 = a_exc0 + fMult(B, exc[i]);
    sum1 = a_exc1 + fMult(B, exc[i + 1]);
    a_exc0 = fMultDiv2(A2, exc[i + 1]);
    a_exc1 = fMultDiv2(A2, exc[i + 2]);
    exc[i] = sum0 + a_exc0;
    exc[i + 1] = sum1 + a_exc1;
  }
  return;
}
