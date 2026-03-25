/**************************** AAC decoder library ******************************

   Author(s):   Matthias Hildenbrand, Manuel Jander

   Description: USAC LPC/AVQ decode

*******************************************************************************/

#include "usacdec_lpc.h"
#include "usacdec_rom.h"
#include "FDK_trigFcts.h"

#define NQ_MAX 36

static int get_vlclbf(HANDLE_FDK_BITSTREAM hBs) {
  int result = 0;

  while (FDKreadBits(hBs, 1) && result <= NQ_MAX) {
    result++;
  }
  return result;
}

static int get_vlclbf_n(HANDLE_FDK_BITSTREAM hBs, int n) {
  int result = 0;

  while (FDKreadBits(hBs, 1)) {
    result++;
    n--;
    if (n <= 0) {
      break;
    }
  }

  return result;
}

#define ZF_SCALE ((DFRACT_BITS / 2))
#define FIXP_ZF FIXP_DBL
#define INT2ZF(x, s) (FIXP_ZF)((x) << (ZF_SCALE - (s)))
#define ZF2INT(x) (INT)((x) >> ZF_SCALE)
#define ONEZF ((FIXP_ZF)INT2ZF(1, 0))

void nearest_neighbor_2D8(FIXP_ZF x[8], int y[8]) {
  FIXP_ZF s, em, e[8];
  int i, j, sum;

  sum = 0;
  for (i = 0; i < 8; i++) {
    FIXP_ZF tmp;
    if (x[i] < (FIXP_ZF)0) {
      tmp = ONEZF - x[i];
      y[i] = -2 * ((ZF2INT(tmp)) >> 1);
    } else {
      tmp = ONEZF + x[i];
      y[i] = 2 * ((ZF2INT(tmp)) >> 1);
    }
    sum += y[i];
  }
  if (sum % 4) {
    em = (FIXP_SGL)0;
    j = 0;
    for (i = 0; i < 8; i++) {
      e[i] = x[i] - INT2ZF(y[i], 0);
    }
    for (i = 0; i < 8; i++) {
      if (e[i] < (FIXP_ZF)0) {
        s = -e[i];
      } else {
        s = e[i];
      }
      if (em < s) {
        em = s;
        j = i;
      }
    }
    if (e[j] < (FIXP_ZF)0) {
      y[j] -= 2;
    } else {
      y[j] += 2;
    }
  }
}

void RE8_PPV(FIXP_ZF x[], SHORT y[], int r) {
  int i, y0[8], y1[8];
  FIXP_ZF x1[8], tmp;
  INT64 e;

  nearest_neighbor_2D8(x, y0);
  for (i = 0; i < 8; i++) {
    x1[i] = x[i] - ONEZF;
  }
  nearest_neighbor_2D8(x1, y1);
  for (i = 0; i < 8; i++) {
    y1[i] += 1;
  }

  e = 0;
  for (i = 0; i < 8; i++) {
    tmp = x[i] - INT2ZF(y0[i], 0);
    e += (INT64)fPow2Div2(
        tmp << r);
    tmp = x[i] - INT2ZF(y1[i], 0);
    e -= (INT64)fPow2Div2(tmp << r);
  }
  if (e < 0) {
    for (i = 0; i < 8; i++) {
      y[i] = y0[i];
    }
  } else {
    for (i = 0; i < 8; i++) {
      y[i] = y1[i];
    }
  }
}

static int table_lookup(const USHORT *table, unsigned int index, int range) {
  int i;

  for (i = 4; i < range; i += 4) {
    if (index < table[i]) {
      break;
    }
  }
  if (i > range) {
    i = range;
  }

  if (index < table[i - 2]) {
    i -= 2;
  }
  if (index < table[i - 1]) {
    i--;
  }
  i--;

  return (i);
}

static void re8_decode_rank_of_permutation(int rank, int *xs, SHORT x[8]) {
  INT a[8], w[8], B, fac, fac_B, target;
  int i, j;

  j = 0;
  w[j] = 1;
  a[j] = xs[0];
  B = 1;
  for (i = 1; i < 8; i++) {
    if (xs[i] != xs[i - 1]) {
      j++;
      w[j] = 1;
      a[j] = xs[i];
    } else {
      w[j]++;
      B *= w[j];
    }
  }

  if (w[0] == 8) {
    for (i = 0; i < 8; i++) {
      x[i] = a[0];
    }
  } else {
    target = rank * B;
    fac_B = 1;
    for (i = 0; i < 8; i++) {
      fac = fac_B * fdk_dec_tab_factorial[i];
      j = -1;
      do {
        target -= w[++j] * fac;
      } while (target >= 0);
      x[i] = a[j];
      target += w[j] * fac;
      fac_B *= w[j];
      w[j]--;
    }
  }
}

static void re8_decode_base_index(int *n, UINT index, SHORT y[8]) {
  int i, im, t, sign_code, ka, ks, rank, leader[8];

  if (*n < 2) {
    for (i = 0; i < 8; i++) {
      y[i] = 0;
    }
  } else {
    switch (*n) {
      case 2:
      case 3:
        i = table_lookup(fdk_dec_I3, index, NB_LDQ3);
        ka = fdk_dec_A3[i];
        break;
      case 4:
        i = table_lookup(fdk_dec_I4, index, NB_LDQ4);
        ka = fdk_dec_A4[i];
        break;
      default:
        FDK_ASSERT(0);
        return;
    }
    for (i = 0; i < 8; i++) {
      leader[i] = fdk_dec_Da[ka][i];
    }
    t = fdk_dec_Ia[ka];
    im = fdk_dec_Ns[ka];
    ks = table_lookup(fdk_dec_Is + t, index, im);

    sign_code = 2 * fdk_dec_Ds[t + ks];
    for (i = 7; i >= 0; i--) {
      leader[i] *= (1 - (sign_code & 2));
      sign_code >>= 1;
    }

    rank = index - fdk_dec_Is[t + ks];
    re8_decode_rank_of_permutation(rank, leader, y);
  }
  return;
}

static void re8_k2y(int *k, int r, SHORT *y) {
  int i, tmp, sum;
  SHORT v[8];
  FIXP_ZF zf[8];

  FDK_ASSERT(r <= ZF_SCALE);

  for (i = 0; i < 8; i++) {
    y[i] = k[7];
  }
  zf[7] = INT2ZF(y[7], r);
  sum = 0;
  for (i = 6; i >= 1; i--) {
    tmp = 2 * k[i];
    sum += tmp;
    y[i] += tmp;
    zf[i] = INT2ZF(y[i], r);
  }
  y[0] += (4 * k[0] + sum);
  zf[0] = INT2ZF(y[0] - 2, r);
  RE8_PPV(zf, v, r);
  for (i = 0; i < 8; i++) {
    y[i] -= (SHORT)(v[i] << r);
  }
}

static int RE8_dec(int n, int I, int *k, FIXP_DBL *y) {
  SHORT v[8];
  SHORT _y[8];
  UINT r;
  int i;

  if (n > NQ_MAX) {
    return -1;
  }

  if (n <= 4) {
    re8_decode_base_index(&n, I, _y);
    for (i = 0; i < 8; i++) {
      y[i] = (LONG)_y[i];
    }
  } else {
    r = ((n - 3) >> 1);

    while (n > 4) {
      n -= 2;
    }
    re8_decode_base_index(&n, I, _y);
    re8_k2y(k, r, v);
    for (i = 0; i < 8; i++) {
      y[i] = (LONG)((_y[i] << r) + v[i]);
    }
  }
  return 0;
}

#define FREQ_MAX 6400.0f
#define FREQ_DIV 400.0f
#define LSF_GAP 50.0f

static void lsf_weight_2st(FIXP_LPC *lsfq, FIXP_DBL *xq, int nk_mode) {
  FIXP_LPC d[M_LP_FILTER_ORDER + 1];
  FIXP_SGL factor;
  LONG w;
  int i;

  d[0] = lsfq[0];
  d[M_LP_FILTER_ORDER] =
      FL2FXCONST_LPC(FREQ_MAX / (1 << LSF_SCALE)) - lsfq[M_LP_FILTER_ORDER - 1];
  for (i = 1; i < M_LP_FILTER_ORDER; i++) {
    d[i] = lsfq[i] - lsfq[i - 1];
  }

  switch (nk_mode) {
    case 0:
      factor = FL2FXCONST_SGL(2.0f * 60.0f / FREQ_DIV);
      break;
    case 1:
      factor = FL2FXCONST_SGL(2.0f * 65.0f / FREQ_DIV);
      break;
    case 2:
      factor = FL2FXCONST_SGL(2.0f * 64.0f / FREQ_DIV);
      break;
    default:
      factor = FL2FXCONST_SGL(2.0f * 63.0f / FREQ_DIV);
      break;
  }
  for (i = 0; i < M_LP_FILTER_ORDER; i++) {
    w = (LONG)fMultDiv2(factor, sqrtFixp(fMult(d[i], d[i + 1])));
    lsfq[i] = fAddSaturate(lsfq[i],
                           FX_DBL2FX_LPC((FIXP_DBL)((INT64)w * (LONG)xq[i])));
  }

  return;
}

static void decode_qn(HANDLE_FDK_BITSTREAM hBs, int nk_mode, int nqn,
                      int qn[]) {
  int n;

  if (nk_mode == 1) {
    for (n = 0; n < nqn; n++) {
      qn[n] = get_vlclbf(hBs);
      if (qn[n] > 0) {
        qn[n]++;
      }
    }
  } else {
    for (n = 0; n < nqn; n++) {
      qn[n] = 2 + FDKreadBits(hBs, 2);
    }
    if (nk_mode == 2) {
      for (n = 0; n < nqn; n++) {
        if (qn[n] > 4) {
          qn[n] = get_vlclbf(hBs);
          if (qn[n] > 0) qn[n] += 4;
        }
      }
    } else {
      for (n = 0; n < nqn; n++) {
        if (qn[n] > 4) {
          qn[n] = get_vlclbf(hBs);
          switch (qn[n]) {
            case 0:
              qn[n] = 5;
              break;
            case 1:
              qn[n] = 6;
              break;
            case 2:
              qn[n] = 0;
              break;
            default:
              qn[n] += 4;
              break;
          }
        }
      }
    }
  }
}

static void reorder_lsf(FIXP_LPC *lsf, FIXP_LPC min_dist, int n) {
  FIXP_LPC lsf_min;
  int i;

  lsf_min = min_dist;
  for (i = 0; i < n; i++) {
    if (lsf[i] < lsf_min) {
      lsf[i] = lsf_min;
    }
    lsf_min = fAddSaturate(lsf[i], min_dist);
  }

  lsf_min = FL2FXCONST_LPC(FREQ_MAX / (1 << LSF_SCALE)) - min_dist;
  for (i = n - 1; i >= 0; i--) {
    if (lsf[i] > lsf_min) {
      lsf[i] = lsf_min;
    }

    lsf_min = lsf[i] - min_dist;
  }
}

static void vlpc_1st_dec(
    HANDLE_FDK_BITSTREAM hBs,
    FIXP_LPC *lsfq
) {
  const FIXP_LPC *p_dico;
  int i, index;

  index = FDKreadBits(hBs, 8);
  p_dico = &fdk_dec_dico_lsf_abs_8b[index * M_LP_FILTER_ORDER];
  for (i = 0; i < M_LP_FILTER_ORDER; i++) {
    lsfq[i] = p_dico[i];
  }
}

static int vlpc_2st_dec(
    HANDLE_FDK_BITSTREAM hBs,
    FIXP_LPC *lsfq,
    int nk_mode
) {
  int err;
  FIXP_DBL xq[M_LP_FILTER_ORDER];

  { err = CLpc_DecodeAVQ(hBs, xq, nk_mode, 2, 8); }
  if (err != 0) {
    return -1;
  }

  lsf_weight_2st(lsfq, xq, nk_mode);
  reorder_lsf(lsfq, FL2FXCONST_LPC(LSF_GAP / (1 << LSF_SCALE)),
              M_LP_FILTER_ORDER);

  return 0;
}

int CLpc_DecodeAVQ(HANDLE_FDK_BITSTREAM hBs, FIXP_DBL *pOutput, int nk_mode,
                   int no_qn, int length) {
  int i, l;

  for (i = 0; i < length; i += 8 * no_qn) {
    int qn[2], nk, n, I;
    int kv[8] = {0};

    decode_qn(hBs, nk_mode, no_qn, qn);

    for (l = 0; l < no_qn; l++) {
      if (qn[l] == 0) {
        FDKmemclear(&pOutput[i + l * 8], 8 * sizeof(FIXP_DBL));
      }

      nk = 0;
      n = qn[l];
      if (qn[l] > 4) {
        nk = (qn[l] - 3) >> 1;
        n = qn[l] - nk * 2;
      }

      I = FDKreadBits(hBs, 4 * n);

      if (nk > 0) {
        int j;

        for (j = 0; j < 8; j++) {
          kv[j] = FDKreadBits(hBs, nk);
        }
      }

      if (RE8_dec(qn[l], I, kv, &pOutput[i + l * 8]) != 0) {
        return -1;
      }
    }
  }
  return 0;
}

int CLpc_Read(HANDLE_FDK_BITSTREAM hBs, FIXP_LPC lsp[][M_LP_FILTER_ORDER],
              FIXP_LPC lpc4_lsf[M_LP_FILTER_ORDER],
              FIXP_LPC lsf_adaptive_mean_cand[M_LP_FILTER_ORDER],
              FIXP_SGL pStability[], UCHAR *mod, int first_lpd_flag,
              int last_lpc_lost, int last_frame_ok) {
  int i, k, err;
  int mode_lpc_bin = 0;
  int lpc_present[5] = {0, 0, 0, 0, 0};
  int lpc0_available = 1;
  int s = 0;
  int l = 3;
  const int nbDiv = NB_DIV;

  lpc_present[4 >> s] = 1;

  vlpc_1st_dec(hBs, lsp[4 >> s]);
  err = vlpc_2st_dec(hBs, lsp[4 >> s], 0);
  if (err != 0) {
    return err;
  }

  k = 0;
  if (!first_lpd_flag) {
    lpc_present[0] = 1;
    lpc0_available = !last_lpc_lost;
    for (i = 0; i < M_LP_FILTER_ORDER; i++) {
      lsp[0][i] = lpc4_lsf[i];
    }
    k = 2;
  }

  for (; k < l; k += 2) {
    int nk_mode = 0;

    if ((k == 2) && (mod[0] == 3)) {
      break;
    }

    lpc_present[k >> s] = 1;

    mode_lpc_bin = FDKreadBit(hBs);

    if (mode_lpc_bin == 0) {
      vlpc_1st_dec(hBs, lsp[k >> s]);
    } else {
      for (i = 0; i < M_LP_FILTER_ORDER; i++) {
        lsp[k >> s][i] = lsp[4 >> s][i];
      }
      nk_mode = 3;
    }

    err = vlpc_2st_dec(hBs, lsp[k >> s], nk_mode);
    if (err != 0) {
      return err;
    }
  }

  if (mod[0] < 2) {
    lpc_present[1] = 1;
    mode_lpc_bin = get_vlclbf_n(hBs, 2);

    switch (mode_lpc_bin) {
      case 1:
        vlpc_1st_dec(hBs, lsp[1]);
        err = vlpc_2st_dec(hBs, lsp[1], 0);
        if (err != 0) {
          return err;
        }
        break;
      case 2:
        if (lpc0_available) {
          for (i = 0; i < M_LP_FILTER_ORDER; i++) {
            lsp[1][i] = (lsp[0][i] >> 1) + (lsp[2][i] >> 1);
          }
        } else {
          for (i = 0; i < M_LP_FILTER_ORDER; i++) {
            lsp[1][i] = lsp[2][i];
          }
        }
        break;
      case 0:
        for (i = 0; i < M_LP_FILTER_ORDER; i++) {
          lsp[1][i] = lsp[2][i];
        }
        err = vlpc_2st_dec(hBs, lsp[1], 2 << s);
        if (err != 0) {
          return err;
        }
        break;
    }
  }

  if ((mod[2] < 2)) {
    int nk_mode = 0;
    lpc_present[3] = 1;

    mode_lpc_bin = get_vlclbf_n(hBs, 3);

    switch (mode_lpc_bin) {
      case 1:
        vlpc_1st_dec(hBs, lsp[3]);
        break;
      case 0:
        for (i = 0; i < M_LP_FILTER_ORDER; i++) {
          lsp[3][i] = (lsp[2][i] >> 1) + (lsp[4][i] >> 1);
        }
        nk_mode = 1;
        break;
      case 2:
        for (i = 0; i < M_LP_FILTER_ORDER; i++) {
          lsp[3][i] = lsp[2][i];
        }
        nk_mode = 2;
        break;
      case 3:
        for (i = 0; i < M_LP_FILTER_ORDER; i++) {
          lsp[3][i] = lsp[4][i];
        }
        nk_mode = 2;
        break;
    }
    err = vlpc_2st_dec(hBs, lsp[3], nk_mode);
    if (err != 0) {
      return err;
    }
  }

  if (!lpc0_available && !last_frame_ok) {
    for (k = 1; k < (nbDiv + 1); k++) {
      if (lpc_present[k]) {
        for (i = 0; i < M_LP_FILTER_ORDER; i++) {
#define LSF_INIT_TILT (0.25f)
          if (mod[0] > 0) {
            lsp[0][i] = FX_DBL2FX_LPC(
                fMult(lsp[k][i], FL2FXCONST_SGL(1.0f - LSF_INIT_TILT)) +
                fMult(fdk_dec_lsf_init[i], FL2FXCONST_SGL(LSF_INIT_TILT)));
          } else {
            lsp[0][i] = lsp[k][i];
          }
        }
        break;
      }
    }
  }

  for (i = 0; i < M_LP_FILTER_ORDER; i++) {
    lpc4_lsf[i] = lsp[4 >> s][i];
  }

  {
    FIXP_DBL divFac;
    int last, numLpc = 0;

    i = nbDiv;
    do {
      numLpc += lpc_present[i--];
    } while (i >= 0 && numLpc < 3);

    last = i;

    switch (numLpc) {
      case 3:
        divFac = FL2FXCONST_DBL(1.0f / 3.0f);
        break;
      case 2:
        divFac = FL2FXCONST_DBL(1.0f / 2.0f);
        break;
      default:
        divFac = FL2FXCONST_DBL(1.0f);
        break;
    }

    for (k = 0; k < M_LP_FILTER_ORDER; k++) {
      FIXP_DBL tmp = (FIXP_DBL)0;
      for (i = nbDiv; i > last; i--) {
        if (lpc_present[i]) {
          tmp = fMultAdd(tmp >> 1, lsp[i][k], divFac);
        }
      }
      lsf_adaptive_mean_cand[k] = FX_DBL2FX_LPC(tmp);
    }
  }
  {
    FIXP_LPC *lsf_prev, *lsf_curr;
    k = 0;

    FDK_ASSERT(lpc_present[0] == 1 && lpc_present[4 >> s] == 1);
    lsf_prev = lsp[0];
    for (i = 1; i < (nbDiv + 1); i++) {
      if (lpc_present[i]) {
        FIXP_DBL tmp = (FIXP_DBL)0;
        int j;
        lsf_curr = lsp[i];

        for (j = 0; j < M_LP_FILTER_ORDER; j++) {
          tmp += fPow2Div2((FIXP_SGL)(lsf_curr[j] - lsf_prev[j])) >> 3;
        }

        tmp = FL2FXCONST_DBL(1.25f / (1 << LSF_SCALE)) -
              fMult(tmp, FL2FXCONST_DBL((1 << (LSF_SCALE + 4)) / 400000.0f));
        if (tmp >= FL2FXCONST_DBL(1.0f / (1 << LSF_SCALE))) {
          pStability[k] = FL2FXCONST_SGL(1.0f / 2.0f);
        } else if (tmp < FL2FXCONST_DBL(0.0f)) {
          pStability[k] = FL2FXCONST_SGL(0.0f);
        } else {
          pStability[k] = FX_DBL2FX_SGL(tmp << (LSF_SCALE - 1));
        }

        lsf_prev = lsf_curr;
        k = i;
      } else {
        pStability[i] = (FIXP_SGL)-1;
      }
    }
  }

  for (i = 0; i < (nbDiv + 1); i++) {
    if (lpc_present[i]) {
      for (k = 0; k < M_LP_FILTER_ORDER; k++) {
        lsp[i][k] = FX_DBL2FX_LPC(
            fixp_cos(fMult(lsp[i][k],
                           FL2FXCONST_SGL((1 << LSPARG_SCALE) * M_PI / 6400.0)),
                     LSF_SCALE - LSPARG_SCALE));
      }
    }
  }

  return 0;
}

void CLpc_Conceal(FIXP_LPC lsp[][M_LP_FILTER_ORDER],
                  FIXP_LPC lpc4_lsf[M_LP_FILTER_ORDER],
                  FIXP_LPC lsf_adaptive_mean[M_LP_FILTER_ORDER],
                  const int first_lpd_flag) {
  int i, j;

#define BETA (FL2FXCONST_SGL(0.25f))
#define ONE_BETA (FL2FXCONST_SGL(0.75f))
#define BFI_FAC (FL2FXCONST_SGL(0.90f))
#define ONE_BFI_FAC (FL2FXCONST_SGL(0.10f))

  if (first_lpd_flag) {
    for (i = 0; i < M_LP_FILTER_ORDER; i++) {
      lsp[0][i] = lpc4_lsf[i] = fdk_dec_lsf_init[i];
    }
  } else {
    for (i = 0; i < M_LP_FILTER_ORDER; i++) {
      lsp[0][i] = lpc4_lsf[i];
    }
  }

  for (i = 0; i < M_LP_FILTER_ORDER; i++) {
    FIXP_LPC lsf_mean = FX_DBL2FX_LPC(fMult(BETA, fdk_dec_lsf_init[i]) +
                                      fMult(ONE_BETA, lsf_adaptive_mean[i]));

    lsp[1][i] = FX_DBL2FX_LPC(fMult(BFI_FAC, lpc4_lsf[i]) +
                              fMult(ONE_BFI_FAC, lsf_mean));
  }

  for (j = 2; j <= 4; j++) {
    for (i = 0; i < M_LP_FILTER_ORDER; i++) {

      FIXP_LPC lsf_mean = FX_DBL2FX_LPC(
          fMult((FIXP_SGL)(BETA + (FIXP_SGL)(j * (INT)FL2FXCONST_SGL(0.1f))),
                (FIXP_SGL)fdk_dec_lsf_init[i]) +
          fMult(
              (FIXP_SGL)(ONE_BETA - (FIXP_SGL)(j * (INT)FL2FXCONST_SGL(0.1f))),
              lsf_adaptive_mean[i]));

      lsp[j][i] = FX_DBL2FX_LPC(fMult(BFI_FAC, lsp[j - 1][i]) +
                                fMult(ONE_BFI_FAC, lsf_mean));
    }
  }

  for (i = 0; i < M_LP_FILTER_ORDER; i++) {
    lpc4_lsf[i] = lsp[4][i];
  }

  for (j = 0; j < 5; j++) {
    for (i = 0; i < M_LP_FILTER_ORDER; i++) {
      lsp[j][i] = FX_DBL2FX_LPC(fixp_cos(
          fMult(lsp[j][i], FL2FXCONST_SGL((1 << LSPARG_SCALE) * M_PI / 6400.0)),
          LSF_SCALE - LSPARG_SCALE));
    }
  }
}

void E_LPC_a_weight(FIXP_LPC *wA, const FIXP_LPC *A, int m) {
  FIXP_DBL f;
  int i;

  f = FL2FXCONST_DBL(0.92f);
  for (i = 0; i < m; i++) {
    wA[i] = FX_DBL2FX_LPC(fMult(A[i], f));
    f = fMult(f, FL2FXCONST_DBL(0.92f));
  }
}

void CLpd_DecodeGain(FIXP_DBL *gain, INT *gain_e, int gain_code) {
  *gain = fLdPow(
      FL2FXCONST_DBL(3.3219280948873623478703194294894 / 4.0),
      2,
      fMultDiv2((FIXP_DBL)gain_code << (DFRACT_BITS - 1 - 7),
                FL2FXCONST_DBL(2.0f / 28.0f)),
      7, gain_e);
}

#define SF_F 8

static void get_lsppol(FIXP_LPC lsp[], FIXP_DBL f[], int n, int flag) {
  FIXP_DBL b;
  FIXP_LPC *plsp;
  int i, j;

  plsp = lsp + flag - 1;
  f[0] = FL2FXCONST_DBL(1.0f / (1 << SF_F));
  b = -FX_LPC2FX_DBL(*plsp);
  f[1] = b >> (SF_F - 1);
  for (i = 2; i <= n; i++) {
    plsp += 2;
    b = -FX_LPC2FX_DBL(*plsp);
    f[i] = SATURATE_LEFT_SHIFT((fMultDiv2(b, f[i - 1]) + (f[i - 2] >> 1)), 2,
                               DFRACT_BITS);
    for (j = i - 1; j > 1; j--) {
      f[j] = SATURATE_LEFT_SHIFT(
          ((f[j] >> 2) + fMultDiv2(b, f[j - 1]) + (f[j - 2] >> 2)), 2,
          DFRACT_BITS);
    }
    f[1] = f[1] + (b >> (SF_F - 1));
  }
  return;
}

#define NC M_LP_FILTER_ORDER / 2

void E_LPC_f_lsp_a_conversion(FIXP_LPC *lsp, FIXP_LPC *a, INT *a_exp) {
  FIXP_DBL f1[NC + 1], f2[NC + 1];
  int i, k;

  get_lsppol(lsp, f1, NC, 1);
  get_lsppol(lsp, f2, NC, 2);

  scaleValues(f1, NC + 1, -2);
  scaleValues(f2, NC + 1, -2);

  for (i = NC; i > 0; i--) {
    f1[i] += f1[i - 1];
    f2[i] -= f2[i - 1];
  }

  FIXP_DBL aDBL[M_LP_FILTER_ORDER];

  for (i = 1, k = M_LP_FILTER_ORDER - 1; i <= NC; i++, k--) {
    aDBL[i - 1] = f1[i] + f2[i];
    aDBL[k] = f1[i] - f2[i];
  }

  int headroom_a = getScalefactor(aDBL, M_LP_FILTER_ORDER);

  for (i = 0; i < M_LP_FILTER_ORDER; i++) {
    a[i] = FX_DBL2FX_LPC(aDBL[i] << headroom_a);
  }

  *a_exp = SF_F + (2 - 1) - headroom_a;
}
