/**************************** AAC decoder library ******************************

   Author(s):   Matthias Hildenbrand

   Description: USAC ACELP frame decoder

*******************************************************************************/

#include "usacdec_acelp.h"

#include "usacdec_ace_d4t64.h"
#include "usacdec_ace_ltp.h"
#include "usacdec_rom.h"
#include "usacdec_lpc.h"
#include "genericStds.h"

#define PIT_FR2_12k8 128
#define PIT_FR1_12k8 160
#define TILT_CODE2 \
  FL2FXCONST_SGL(0.3f * 2.0f)
#define PIT_SHARP \
  FL2FXCONST_SGL(0.85f)
#define PREEMPH_FAC \
  FL2FXCONST_SGL(0.68f)

#define ACELP_HEADROOM 1
#define ACELP_OUTSCALE (MDCT_OUT_HEADROOM - ACELP_HEADROOM)

void E_UTIL_preemph(const FIXP_DBL *in, FIXP_DBL *out, INT L) {
  int i;

  for (i = 0; i < L; i++) {
    out[i] = fAddSaturate(in[i], -fMult(PREEMPH_FAC, in[i - 1]));
  }

  return;
}

static void Preemph_code(
    FIXP_COD x[]
) {
  int i;
  FIXP_DBL L_tmp;

  for (i = L_SUBFR - 1; i > 0; i--) {
    L_tmp = FX_COD2FX_DBL(x[i]);
    L_tmp -= fMultDiv2(x[i - 1], TILT_CODE2);
    x[i] = FX_DBL2FX_COD(L_tmp);
  }
}

static void Pit_shrp(
    FIXP_COD x[],
    int pit_lag
) {
  int i;
  FIXP_DBL L_tmp;

  for (i = pit_lag; i < L_SUBFR; i++) {
    L_tmp = FX_COD2FX_DBL(x[i]);
    L_tmp += fMult(x[i - pit_lag], PIT_SHARP);
    x[i] = FX_DBL2FX_COD(L_tmp);
  }

  return;
}

#define SF_MEAN_ENER_LG10 9

static const FIXP_DBL pow_10_mean_energy[4] = {0x01fc5ebd, 0x07e7db92,
                                               0x1f791f65, 0x7d4bfba3};

static void D_gain2_plus(int index, FIXP_COD code[], FIXP_SGL *gain_pit,
                         FIXP_DBL *gain_code, int mean_ener_bits, int bfi,
                         FIXP_SGL *past_gpit, FIXP_DBL *past_gcode,
                         FIXP_DBL *pEner_code, int *pEner_code_e) {
  FIXP_DBL Ltmp;
  FIXP_DBL gcode0, gcode_inov;
  INT gcode0_e, gcode_inov_e;
  int i;

  FIXP_DBL ener_code;
  INT ener_code_e;

  ener_code = FIXP_DBL(0);
  for (i = 0; i < L_SUBFR; i++) {
    ener_code += fPow2Div2(code[i]);
  }

  ener_code_e = fMax(fNorm(ener_code) - 1, 0);
  ener_code <<= ener_code_e;
  ener_code_e = 2 * SF_CODE + 1 - ener_code_e;

  *pEner_code = ener_code;
  *pEner_code_e = ener_code_e;

  ener_code += scaleValue(FL2FXCONST_DBL(0.01f), -ener_code_e);

  if (ener_code_e & 1) {
    ener_code_e -= 5;
    ener_code >>= 1;
  } else {
    ener_code_e -= 6;
  }
  gcode_inov = invSqrtNorm2(ener_code, &gcode0_e);
  gcode_inov_e = gcode0_e - (ener_code_e >> 1);

  if (bfi) {
    FIXP_DBL tgcode;
    FIXP_SGL tgpit;

    tgpit = *past_gpit;

    if (tgpit > FL2FXCONST_SGL(0.95f / (1 << SF_GAIN_P))) {
      tgpit = FL2FXCONST_SGL(0.95f / (1 << SF_GAIN_P));
    } else if (tgpit < FL2FXCONST_SGL(0.5f / (1 << SF_GAIN_P))) {
      tgpit = FL2FXCONST_SGL(0.5f / (1 << SF_GAIN_P));
    }
    *gain_pit = tgpit;
    tgpit = FX_DBL2FX_SGL(fMult(tgpit, FL2FXCONST_DBL(0.95f)));
    *past_gpit = tgpit;

    tgpit = FL2FXCONST_SGL(1.4f / (1 << SF_GAIN_P)) - tgpit;
    tgcode = fMult(*past_gcode, tgpit) << SF_GAIN_P;
    *gain_code = scaleValue(fMult(tgcode, gcode_inov), gcode_inov_e);
    *past_gcode = tgcode;

    return;
  }

  gcode0 = pow_10_mean_energy[mean_ener_bits];
  gcode0 = fMultDiv2(gcode0, gcode_inov);
  gcode0_e = gcode0_e + SF_MEAN_ENER_LG10 - (ener_code_e >> 1) + 1;

  i = index << 1;
  *gain_pit = fdk_t_qua_gain7b[i];
  Ltmp = fMult(fdk_t_qua_gain7b[i + 1], gcode0);
  *gain_code = scaleValue(Ltmp, gcode0_e - SF_GAIN_C + SF_QUA_GAIN7B);

  *past_gpit = *gain_pit;

  {
    FIXP_DBL gcode_m;
    INT gcode_e;

    gcode_m = fDivNormHighPrec(Ltmp, gcode_inov, &gcode_e);
    gcode_e += (gcode0_e - SF_GAIN_C + SF_QUA_GAIN7B) - (gcode_inov_e);
    *past_gcode = scaleValue(gcode_m, gcode_e);
  }
}

static FIXP_DBL calc_period_factor(FIXP_DBL exc[], FIXP_SGL gain_pit,
                                   FIXP_DBL gain_code, FIXP_DBL ener_code,
                                   int ener_code_e) {
  int ener_exc_e, L_tmp_e, s = 0;
  FIXP_DBL ener_exc, L_tmp;
  FIXP_DBL period_fac;

  ener_exc = (FIXP_DBL)0;
  for (int i = 0; i < L_SUBFR; i++) {
    ener_exc += fPow2Div2(exc[i]) >> s;
    if (ener_exc >= FL2FXCONST_DBL(0.5f)) {
      ener_exc >>= 1;
      s++;
    }
  }

  ener_exc_e = fNorm(ener_exc);
  ener_exc = fMult(ener_exc << ener_exc_e, fPow2(gain_pit));
  if (ener_exc != (FIXP_DBL)0) {
    ener_exc_e = 2 * SF_EXC + 1 + 2 * SF_GAIN_P - ener_exc_e + s;
  } else {
    ener_exc_e = 0;
  }

  L_tmp_e = fNorm(gain_code);
  L_tmp = fPow2(gain_code << L_tmp_e);
  L_tmp = fMult(ener_code, L_tmp);
  L_tmp_e = 2 * SF_GAIN_C + ener_code_e - 2 * L_tmp_e;

  {
    FIXP_DBL num, den;
    int exp_diff;

    exp_diff = ener_exc_e - L_tmp_e;
    if (exp_diff >= 0) {
      ener_exc >>= 1;
      if (exp_diff <= DFRACT_BITS - 2) {
        L_tmp >>= exp_diff + 1;
      } else {
        L_tmp = (FIXP_DBL)0;
      }
      den = ener_exc + L_tmp;
      if (ener_exc_e < DFRACT_BITS - 1) {
        den += scaleValue(FL2FXCONST_DBL(0.01f), -ener_exc_e - 1);
      }
    } else {
      if (exp_diff >= -(DFRACT_BITS - 2)) {
        ener_exc >>= 1 - exp_diff;
      } else {
        ener_exc = (FIXP_DBL)0;
      }
      L_tmp >>= 1;
      den = ener_exc + L_tmp;
      if (L_tmp_e < DFRACT_BITS - 1) {
        den += scaleValue(FL2FXCONST_DBL(0.01f), -L_tmp_e - 1);
      }
    }
    num = (ener_exc - L_tmp);
    num >>= SF_PFAC;

    if (den > (FIXP_DBL)0) {
      if (ener_exc > L_tmp) {
        period_fac = schur_div(num, den, 16);
      } else {
        period_fac = -schur_div(-num, den, 16);
      }
    } else {
      period_fac = (FIXP_DBL)MAXVAL_DBL;
    }
  }

  return period_fac;
}

static FIXP_DBL
noise_enhancer(
               FIXP_DBL gain_code,
               FIXP_DBL period_fac,
               FIXP_SGL stab_fac,
               FIXP_DBL
                   *p_gc_threshold)
{
  FIXP_DBL fac, L_tmp, gc_thres;

  gc_thres = *p_gc_threshold;

  L_tmp = gain_code;
  if (L_tmp < gc_thres) {
    L_tmp += fMultDiv2(gain_code,
                       FL2FXCONST_SGL(2.0 * 0.19f));
    if (L_tmp > gc_thres) {
      L_tmp = gc_thres;
    }
  } else {
    L_tmp = fMult(gain_code,
                  FL2FXCONST_SGL(1.0f / 1.19f));
    if (L_tmp < gc_thres) {
      L_tmp = gc_thres;
    }
  }
  *p_gc_threshold = L_tmp;

  fac = (FX_SGL2FX_DBL(stab_fac) >> (SF_PFAC + 1)) -
        fMultDiv2(stab_fac, period_fac);
  FDK_ASSERT(fac >= (FIXP_DBL)0);

  gain_code = fMult(fac, L_tmp) -
              fMult(FL2FXCONST_DBL(-1.0f / (1 << (SF_PFAC + SF_STAB))) + fac,
                    gain_code);
  gain_code <<= (SF_PFAC + SF_STAB);

  return gain_code;
}

void BuildAdaptiveExcitation(
    FIXP_COD code[],
    FIXP_DBL exc[],
    FIXP_SGL gain_pit,
    FIXP_DBL gain_code,
    FIXP_DBL gain_code_smoothed,
    FIXP_DBL period_fac,
    FIXP_DBL exc2[]
) {
#define SF_HEADROOM (1)
#define SF (SF_CODE + SF_GAIN_C + 1 - SF_EXC - SF_HEADROOM)
#define SF_GAIN_P2 (SF_GAIN_P - SF_HEADROOM)

  int i;
  FIXP_DBL tmp, cpe, code_smooth_prev, code_smooth;

  FIXP_COD code_i;
  FIXP_DBL cpe_code_smooth, cpe_code_smooth_prev;

  cpe = (period_fac >> (2 - SF_PFAC)) + FL2FXCONST_DBL(0.25f);

  tmp = fMultDiv2(*exc, gain_pit) << (SF_GAIN_P2 + 1);
  *exc++ = (tmp + (fMultDiv2(code[0], gain_code) << SF)) << SF_HEADROOM;

  code_smooth_prev = fMultDiv2(*code++, gain_code_smoothed)
                     << SF;
  code_i = *code++;
  code_smooth = fMultDiv2(code_i, gain_code_smoothed) << SF;
  tmp += code_smooth_prev;
  cpe_code_smooth = fMultDiv2(cpe, code_smooth);
  *exc2++ = (tmp - cpe_code_smooth) << SF_HEADROOM;
  cpe_code_smooth_prev = fMultDiv2(cpe, code_smooth_prev);

  i = L_SUBFR - 2;
  do
  {
    tmp = fMultDiv2(*exc, gain_pit) << (SF_GAIN_P2 + 1);
    *exc++ = (tmp + (fMultDiv2(code_i, gain_code) << SF)) << SF_HEADROOM;
    tmp += code_smooth;
    tmp -= cpe_code_smooth_prev;
    cpe_code_smooth_prev = cpe_code_smooth;
    code_i = *code++;
    code_smooth = fMultDiv2(code_i, gain_code_smoothed) << SF;
    cpe_code_smooth = fMultDiv2(cpe, code_smooth);
    *exc2++ = (tmp - cpe_code_smooth)
              << SF_HEADROOM;
  } while (--i != 0);

  tmp = fMultDiv2(*exc, gain_pit) << (SF_GAIN_P2 + 1);
  *exc = (tmp + (fMultDiv2(code_i, gain_code) << SF)) << SF_HEADROOM;
  tmp += code_smooth;
  tmp -= cpe_code_smooth_prev;
  *exc2++ = tmp << SF_HEADROOM;

  return;
}

void int_lpc_acelp(
    const FIXP_LPC lsp_old[],
    const FIXP_LPC lsp_new[],
    int subfr_nr, int nb_subfr,
    FIXP_LPC
        A[],
    INT *A_exp) {
  int i;
  FIXP_LPC lsp_interpol[M_LP_FILTER_ORDER];
  FIXP_SGL fac_old, fac_new;

  FDK_ASSERT((nb_subfr == 3) || (nb_subfr == 4));

  fac_old = lsp_interpol_factor[nb_subfr & 0x1][(nb_subfr - 1) - subfr_nr];
  fac_new = lsp_interpol_factor[nb_subfr & 0x1][subfr_nr];
  for (i = 0; i < M_LP_FILTER_ORDER; i++) {
    lsp_interpol[i] = FX_DBL2FX_LPC(
        (fMultDiv2(lsp_old[i], fac_old) + fMultDiv2(lsp_new[i], fac_new)) << 1);
  }

  E_LPC_f_lsp_a_conversion(lsp_interpol, A, A_exp);

  return;
}

void Syn_filt(const FIXP_LPC a[],
              const INT a_exp,
              INT length,
              FIXP_DBL x[],
              FIXP_DBL y[]
) {
  int i, j;
  FIXP_DBL L_tmp;

  for (i = 0; i < length; i++) {
    L_tmp = (FIXP_DBL)0;

    for (j = 0; j < M_LP_FILTER_ORDER; j++) {
      L_tmp -= fMultDiv2(a[j], y[i - (j + 1)]) >> (LP_FILTER_SCALE - 1);
    }

    L_tmp = scaleValue(L_tmp, a_exp + LP_FILTER_SCALE);
    y[i] = fAddSaturate(L_tmp, x[i]);
  }

  return;
}

void Deemph(FIXP_DBL *x, FIXP_DBL *y, int L, FIXP_DBL *mem) {
  int i;
  FIXP_DBL yi = *mem;

  for (i = 0; i < L; i++) {
    FIXP_DBL xi = x[i] >> 1;
    xi = fMultAddDiv2(xi, PREEMPH_FAC, yi);
    yi = SATURATE_LEFT_SHIFT(xi, 1, 32);
    y[i] = yi;
  }
  *mem = yi;
  return;
}

void E_UTIL_residu(const FIXP_LPC *a, const INT a_exp, FIXP_DBL *x, FIXP_DBL *y,
                   INT l) {
  FIXP_DBL s;
  INT i, j;

  for (i = 0; i < l; i++) {
    s = (FIXP_DBL)0;

    for (j = 0; j < M_LP_FILTER_ORDER; j++) {
      s += fMultDiv2(a[j], x[i - j - 1]) >> (LP_FILTER_SCALE - 1);
    }

    s = scaleValue(s, a_exp + LP_FILTER_SCALE);
    y[i] = fAddSaturate(s, x[i]);
  }

  return;
}

static const UCHAR num_acb_idx_bits_table[2][NB_SUBFR] = {
    {9, 6, 9, 6},
    {9, 6, 6, 0}
};

static int DecodePitchLag(HANDLE_FDK_BITSTREAM hBs,
                          const UCHAR num_acb_idx_bits,
                          const int PIT_MIN,
                          const int PIT_FR2,
                          const int PIT_FR1,
                          const int PIT_MAX,
                          int *pT0, int *pT0_frac, int *pT0_min, int *pT0_max) {
  int acb_idx;
  int error = 0;
  int T0, T0_frac;

  FDK_ASSERT((num_acb_idx_bits == 9) || (num_acb_idx_bits == 6));

  acb_idx = FDKreadBits(hBs, num_acb_idx_bits);

  if (num_acb_idx_bits == 6) {
    T0 = *pT0_min + acb_idx / 4;
    T0_frac = acb_idx & 0x3;
  } else {
    int T0_min, T0_max;

    if (acb_idx < (PIT_FR2 - PIT_MIN) * 4) {
      T0 = PIT_MIN + (acb_idx / 4);
      T0_frac = acb_idx & 0x3;
    } else if (acb_idx < ((PIT_FR2 - PIT_MIN) * 4 + (PIT_FR1 - PIT_FR2) * 2)) {
      acb_idx -= (PIT_FR2 - PIT_MIN) * 4;
      T0 = PIT_FR2 + (acb_idx / 2);
      T0_frac = (acb_idx & 0x1) * 2;
    } else {
      T0 = acb_idx + PIT_FR1 - ((PIT_FR2 - PIT_MIN) * 4) -
           ((PIT_FR1 - PIT_FR2) * 2);
      T0_frac = 0;
    }
    T0_min = T0 - 8;
    if (T0_min < PIT_MIN) {
      T0_min = PIT_MIN;
    }
    T0_max = T0_min + 15;
    if (T0_max > PIT_MAX) {
      T0_max = PIT_MAX;
      T0_min = T0_max - 15;
    }
    *pT0_min = T0_min;
    *pT0_max = T0_max;
  }
  *pT0 = T0;
  *pT0_frac = T0_frac;

  return error;
}
static void ConcealPitchLag(CAcelpStaticMem *acelp_mem, const int PIT_MAX,
                            int *pT0, int *pT0_frac) {
  USHORT *pold_T0 = &acelp_mem->old_T0;
  UCHAR *pold_T0_frac = &acelp_mem->old_T0_frac;

  if ((int)*pold_T0 >= PIT_MAX) {
    *pold_T0 = (USHORT)(PIT_MAX - 5);
  }
  *pT0 = (int)*pold_T0;
  *pT0_frac = (int)*pold_T0_frac;
}

static UCHAR tab_coremode2nbits[8] = {20, 28, 36, 44, 52, 64, 12, 16};

static int MapCoreMode2NBits(int core_mode) {
  return (int)tab_coremode2nbits[core_mode];
}

void CLpd_AcelpDecode(CAcelpStaticMem *acelp_mem, INT i_offset,
                      const FIXP_LPC lsp_old[M_LP_FILTER_ORDER],
                      const FIXP_LPC lsp_new[M_LP_FILTER_ORDER],
                      FIXP_SGL stab_fac, CAcelpChannelData *pAcelpData,
                      INT numLostSubframes, int lastLpcLost, int frameCnt,
                      FIXP_DBL synth[], int pT[], FIXP_DBL *pit_gain,
                      INT coreCoderFrameLength) {
  int i_subfr, subfr_nr, l_div, T;
  int T0 = -1, T0_frac = -1;

  int pit_gain_index = 0;

  const int PIT_MAX = PIT_MAX_12k8 + (6 * i_offset);

  FIXP_COD *code;
  FIXP_DBL *exc2;
  FIXP_DBL *syn;
  FIXP_DBL *exc;
  FIXP_LPC A[M_LP_FILTER_ORDER];
  INT A_exp;

  FIXP_DBL period_fac;
  FIXP_SGL gain_pit;
  FIXP_DBL gain_code, gain_code_smooth, Ener_code;
  int Ener_code_e;
  int n;
  int bfi = (numLostSubframes > 0) ? 1 : 0;

  C_ALLOC_SCRATCH_START(
      exc_buf, FIXP_DBL,
      PIT_MAX_MAX + L_INTERPOL + L_DIV + 1);
  C_ALLOC_SCRATCH_START(syn_buf, FIXP_DBL,
                        M_LP_FILTER_ORDER + L_DIV);
  C_ALLOC_SCRATCH_START(tmp_buf, FIXP_DBL, L_SUBFR);
#if (COD_BITS == FRACT_BITS)
  code = (FIXP_COD *)(tmp_buf + L_SUBFR / 2);
#elif (COD_BITS == DFRACT_BITS)
  code = (FIXP_COD *)tmp_buf;
#endif
  exc2 = (FIXP_DBL *)tmp_buf;

  syn = syn_buf + M_LP_FILTER_ORDER;
  exc = exc_buf + PIT_MAX_MAX + L_INTERPOL;

  FDKmemcpy(syn_buf, acelp_mem->old_syn_mem,
            M_LP_FILTER_ORDER * sizeof(FIXP_DBL));
  FDKmemcpy(exc_buf, acelp_mem->old_exc_mem,
            (PIT_MAX_MAX + L_INTERPOL) * sizeof(FIXP_DBL));

  FDKmemclear(exc_buf + (PIT_MAX_MAX + L_INTERPOL),
              (L_DIV + 1) * sizeof(FIXP_DBL));

  l_div = coreCoderFrameLength / NB_DIV;

  for (i_subfr = 0, subfr_nr = 0; i_subfr < l_div;
       i_subfr += L_SUBFR, subfr_nr++) {
    if (bfi) {
      ConcealPitchLag(acelp_mem, PIT_MAX, &T0, &T0_frac);
    } else {
      T0 = (int)pAcelpData->T0[subfr_nr];
      T0_frac = (int)pAcelpData->T0_frac[subfr_nr];
    }

    Pred_lt4(&exc[i_subfr], T0, T0_frac);

    if ((!bfi && pAcelpData->ltp_filtering_flag[subfr_nr] == 0) ||
        (bfi && numLostSubframes == 1 && stab_fac < FL2FXCONST_SGL(0.25f))) {
      Pred_lt4_postfilter(&exc[i_subfr]);
    }

    if (bfi) {
      for (n = 0; n < L_SUBFR; n++) {
        code[n] =
            FX_SGL2FX_COD((FIXP_SGL)E_UTIL_random(&acelp_mem->seed_ace)) >> 4;
      }
    } else {
      int nbits = MapCoreMode2NBits((int)pAcelpData->acelp_core_mode);
      D_ACELP_decode_4t64(pAcelpData->icb_index[subfr_nr], nbits, &code[0]);
    }

    T = T0;
    if (T0_frac > 2) {
      T += 1;
    }

    Preemph_code(code);
    Pit_shrp(code, T);

    if (T > PIT_MAX) {
      pT[subfr_nr] = PIT_MAX;
    } else {
      pT[subfr_nr] = T;
    }
    D_gain2_plus(
        pAcelpData->gains[subfr_nr],
        code,
        &gain_pit,
        &gain_code,
        pAcelpData
            ->mean_energy,
        bfi, &acelp_mem->past_gpit, &acelp_mem->past_gcode,
        &Ener_code,
        &Ener_code_e);

    pit_gain[pit_gain_index++] = FX_SGL2FX_DBL(gain_pit);

    period_fac =
        calc_period_factor(
                           &exc[i_subfr],
                           gain_pit,
                           gain_code,
                           Ener_code,
                           Ener_code_e);

    if (lastLpcLost && frameCnt == 0) {
      if (gain_pit > FL2FXCONST_SGL(1.0f / (1 << SF_GAIN_P))) {
        gain_pit = FL2FXCONST_SGL(1.0f / (1 << SF_GAIN_P));
      }
    }

    gain_code_smooth =
        noise_enhancer(
                       gain_code,
                       period_fac,
                       stab_fac,
                       &acelp_mem->gc_threshold);

    BuildAdaptiveExcitation(code, exc + i_subfr, gain_pit, gain_code,
                            gain_code_smooth, period_fac, exc2);

    int_lpc_acelp(lsp_old,
                  lsp_new,
                  subfr_nr,
                  coreCoderFrameLength / L_DIV,
                  A,
                  &A_exp);

    Syn_filt(A,
             A_exp, L_SUBFR,
             exc2,
             &syn[i_subfr]
    );

  }

  acelp_mem->old_T0_frac = T0_frac;
  acelp_mem->old_T0 = T0;

  FDKmemcpy(acelp_mem->old_exc_mem, exc + l_div - (PIT_MAX_MAX + L_INTERPOL),
            sizeof(FIXP_DBL) * (PIT_MAX_MAX + L_INTERPOL));
  FDKmemcpy(acelp_mem->old_syn_mem, syn_buf + l_div,
            sizeof(FIXP_DBL) * M_LP_FILTER_ORDER);

  Deemph(syn, synth, l_div,
         &acelp_mem->de_emph_mem);

  scaleValues(synth, l_div, -ACELP_OUTSCALE);
  acelp_mem->deemph_mem_wsyn = acelp_mem->de_emph_mem;

  C_ALLOC_SCRATCH_END(tmp_buf, FIXP_DBL, L_SUBFR);
  C_ALLOC_SCRATCH_END(syn_buf, FIXP_DBL, M_LP_FILTER_ORDER + L_DIV);
  C_ALLOC_SCRATCH_END(exc_buf, FIXP_DBL, PIT_MAX_MAX + L_INTERPOL + L_DIV + 1);
  return;
}

void CLpd_AcelpReset(CAcelpStaticMem *acelp) {
  acelp->gc_threshold = (FIXP_DBL)0;

  acelp->past_gpit = (FIXP_SGL)0;
  acelp->past_gcode = (FIXP_DBL)0;
  acelp->old_T0 = 64;
  acelp->old_T0_frac = 0;
  acelp->deemph_mem_wsyn = (FIXP_DBL)0;
  acelp->wsyn_rms = (FIXP_DBL)0;
  acelp->seed_ace = 0;
}

void CLpd_TcxTDConceal(CAcelpStaticMem *acelp_mem, SHORT *pitch,
                       const FIXP_LPC lsp_old[M_LP_FILTER_ORDER],
                       const FIXP_LPC lsp_new[M_LP_FILTER_ORDER],
                       const FIXP_SGL stab_fac, INT nLostSf, FIXP_DBL synth[],
                       INT coreCoderFrameLength, UCHAR last_tcx_noise_factor) {
  C_ALLOC_SCRATCH_START(
      exc_buf, FIXP_DBL,
      PIT_MAX_MAX + L_INTERPOL + L_DIV);
  C_ALLOC_SCRATCH_START(syn_buf, FIXP_DBL,
                        M_LP_FILTER_ORDER + L_DIV);
  FIXP_DBL ns_buf[L_DIV + 1];
  FIXP_DBL *syn = syn_buf + M_LP_FILTER_ORDER;
  FIXP_DBL *exc = exc_buf + PIT_MAX_MAX + L_INTERPOL;
  FIXP_DBL *ns = ns_buf + 1;
  FIXP_DBL tmp, fact_exc;
  INT T = fMin(*pitch, (SHORT)PIT_MAX_MAX);
  int i, i_subfr, subfr_nr;
  int lDiv = coreCoderFrameLength / NB_DIV;

  FDKmemcpy(syn_buf, acelp_mem->old_syn_mem,
            M_LP_FILTER_ORDER * sizeof(FIXP_DBL));
  FDKmemcpy(exc_buf, acelp_mem->old_exc_mem,
            (PIT_MAX_MAX + L_INTERPOL) * sizeof(FIXP_DBL));

  if (nLostSf < 2) {
    fact_exc = FL2FXCONST_DBL(0.8f);
  } else {
    fact_exc = FL2FXCONST_DBL(0.4f);
  }

  for (i = 0; i < lDiv; i++) {
    exc[i] = fMult(fact_exc, exc[i - T]);
  }

  tmp = fMult(fact_exc, acelp_mem->wsyn_rms);
  acelp_mem->wsyn_rms = tmp;

  acelp_mem->deemph_mem_wsyn = exc[-1];

  ns[-1] = acelp_mem->deemph_mem_wsyn;

  for (i_subfr = 0, subfr_nr = 0; i_subfr < lDiv;
       i_subfr += L_SUBFR, subfr_nr++) {
    FIXP_DBL tRes[L_SUBFR];
    FIXP_LPC A[M_LP_FILTER_ORDER];
    INT A_exp;

    int_lpc_acelp(lsp_old, lsp_new, subfr_nr, lDiv / L_SUBFR, A, &A_exp);

    Syn_filt(A,
             A_exp, L_SUBFR,
             &exc[i_subfr],
             &syn[i_subfr]
    );

    E_LPC_a_weight(
        A, A,
        M_LP_FILTER_ORDER);

    E_UTIL_residu(A, A_exp, &syn[i_subfr], tRes, L_SUBFR);

    Deemph(tRes, &ns[i_subfr], L_SUBFR, &acelp_mem->deemph_mem_wsyn);

    for (i = i_subfr; i < i_subfr + L_SUBFR; i++) {
      if (ns[i] > tmp) {
        ns[i] = tmp;
      } else {
        if (ns[i] < -tmp) {
          ns[i] = -tmp;
        }
      }
    }

    E_UTIL_preemph(&ns[i_subfr], tRes, L_SUBFR);

    Syn_filt(A,
             A_exp, L_SUBFR,
             tRes,
             &syn[i_subfr]
    );

    FDKmemmove(&synth[i_subfr], &syn[i_subfr], L_SUBFR * sizeof(FIXP_DBL));
  }

  FDKmemcpy(acelp_mem->old_exc_mem, exc + lDiv - (PIT_MAX_MAX + L_INTERPOL),
            sizeof(FIXP_DBL) * (PIT_MAX_MAX + L_INTERPOL));
  FDKmemcpy(acelp_mem->old_syn_mem, syn_buf + lDiv,
            sizeof(FIXP_DBL) * M_LP_FILTER_ORDER);
  acelp_mem->de_emph_mem = acelp_mem->deemph_mem_wsyn;

  C_ALLOC_SCRATCH_END(syn_buf, FIXP_DBL, M_LP_FILTER_ORDER + L_DIV);
  C_ALLOC_SCRATCH_END(exc_buf, FIXP_DBL, PIT_MAX_MAX + L_INTERPOL + L_DIV);
}

void Acelp_PreProcessing(FIXP_DBL *synth_buf, FIXP_DBL *old_synth, INT *pitch,
                         INT *old_T_pf, FIXP_DBL *pit_gain,
                         FIXP_DBL *old_gain_pf, INT samplingRate, INT *i_offset,
                         INT coreCoderFrameLength, INT synSfd,
                         INT nbSubfrSuperfr) {
  int n;

  FDKmemcpy(synth_buf, old_synth, sizeof(FIXP_DBL) * (PIT_MAX_MAX - BPF_DELAY));

  *i_offset =
      (samplingRate * PIT_MIN_12k8 + (FSCALE_DENOM / 2)) / FSCALE_DENOM -
      PIT_MIN_12k8;

  for (n = 0; n < synSfd; n++) {
    pitch[n] = old_T_pf[n];
    pit_gain[n] = old_gain_pf[n];
  }
  for (n = 0; n < nbSubfrSuperfr; n++) {
    pitch[n + synSfd] = L_SUBFR;
    pit_gain[n + synSfd] = (FIXP_DBL)0;
  }
}

void Acelp_PostProcessing(FIXP_DBL *synth_buf, FIXP_DBL *old_synth, INT *pitch,
                          INT *old_T_pf, INT coreCoderFrameLength, INT synSfd,
                          INT nbSubfrSuperfr) {
  int n;

  FDKmemcpy(old_synth, synth_buf + coreCoderFrameLength,
            sizeof(FIXP_DBL) * (PIT_MAX_MAX - BPF_DELAY));

  for (n = 0; n < synSfd; n++) {
    old_T_pf[n] = pitch[nbSubfrSuperfr + n];
  }
}

#define L_FAC_ZIR (LFAC)

void CLpd_Acelp_Zir(const FIXP_LPC A[], const INT A_exp,
                    CAcelpStaticMem *acelp_mem, const INT length,
                    FIXP_DBL zir[], int doDeemph) {
  C_ALLOC_SCRATCH_START(tmp_buf, FIXP_DBL, L_FAC_ZIR + M_LP_FILTER_ORDER);
  FDK_ASSERT(length <= L_FAC_ZIR);

  FDKmemcpy(tmp_buf, acelp_mem->old_syn_mem,
            M_LP_FILTER_ORDER * sizeof(FIXP_DBL));
  FDKmemset(tmp_buf + M_LP_FILTER_ORDER, 0, L_FAC_ZIR * sizeof(FIXP_DBL));

  Syn_filt(A, A_exp, length, &tmp_buf[M_LP_FILTER_ORDER],
           &tmp_buf[M_LP_FILTER_ORDER]);
  if (!doDeemph) {
    FDKmemcpy(zir, tmp_buf, length * sizeof(*zir));
  } else {
    Deemph(&tmp_buf[M_LP_FILTER_ORDER], &zir[0], length,
           &acelp_mem->de_emph_mem);
    scaleValues(zir, length, -ACELP_OUTSCALE);
  }
  C_ALLOC_SCRATCH_END(tmp_buf, FIXP_DBL, L_FAC_ZIR + M_LP_FILTER_ORDER);
}

void CLpd_AcelpPrepareInternalMem(const FIXP_DBL *synth, UCHAR last_lpd_mode,
                                  UCHAR last_last_lpd_mode,
                                  const FIXP_LPC *A_new, const INT A_new_exp,
                                  const FIXP_LPC *A_old, const INT A_old_exp,
                                  CAcelpStaticMem *acelp_mem,
                                  INT coreCoderFrameLength, INT clearOldExc,
                                  UCHAR lpd_mode) {
  int l_div =
      coreCoderFrameLength / NB_DIV;
  int l_div_partial;
  FIXP_DBL *syn, *old_exc_mem;

  C_ALLOC_SCRATCH_START(synth_buf, FIXP_DBL,
                        PIT_MAX_MAX + L_INTERPOL + M_LP_FILTER_ORDER);
  syn = &synth_buf[M_LP_FILTER_ORDER];

  l_div_partial = PIT_MAX_MAX + L_INTERPOL - l_div;
  old_exc_mem = acelp_mem->old_exc_mem;

  if (lpd_mode == 4) {
    FDKmemcpy(
        synth_buf, &synth[-(PIT_MAX_MAX + L_INTERPOL + M_LP_FILTER_ORDER)],
        (PIT_MAX_MAX + L_INTERPOL + M_LP_FILTER_ORDER) * sizeof(FIXP_DBL));
    acelp_mem->deemph_mem_wsyn = scaleValueSaturate(synth[-1], ACELP_OUTSCALE);
  } else {
    E_UTIL_preemph(&synth[-(PIT_MAX_MAX + L_INTERPOL + M_LP_FILTER_ORDER)],
                   synth_buf, PIT_MAX_MAX + L_INTERPOL + M_LP_FILTER_ORDER);
    scaleValuesSaturate(synth_buf, PIT_MAX_MAX + L_INTERPOL + M_LP_FILTER_ORDER,
                        ACELP_OUTSCALE);
  }

  acelp_mem->de_emph_mem = scaleValueSaturate(synth[-1], ACELP_OUTSCALE);

  FDKmemcpy(acelp_mem->old_syn_mem,
            &syn[PIT_MAX_MAX + L_INTERPOL - M_LP_FILTER_ORDER],
            M_LP_FILTER_ORDER * sizeof(FIXP_DBL));

  if (clearOldExc) {
    FDKmemclear(old_exc_mem, (PIT_MAX_MAX + L_INTERPOL) * sizeof(FIXP_DBL));
    C_ALLOC_SCRATCH_END(synth_buf, FIXP_DBL,
                        PIT_MAX_MAX + L_INTERPOL + M_LP_FILTER_ORDER);
    return;
  }

  if (last_lpd_mode == 1) {
    if (last_last_lpd_mode == 0) {
      FDKmemmove(old_exc_mem, old_exc_mem + l_div,
                 sizeof(FIXP_DBL) * l_div_partial);
    } else if (last_last_lpd_mode > 0) {
      E_UTIL_residu(A_old, A_old_exp, syn, old_exc_mem, l_div_partial);
    }
    E_UTIL_residu(A_new, A_new_exp, syn + l_div_partial,
                  old_exc_mem + l_div_partial, l_div);
  } else {
    int exc_A_new_length = (coreCoderFrameLength / 2 > PIT_MAX_MAX + L_INTERPOL)
                               ? PIT_MAX_MAX + L_INTERPOL
                               : coreCoderFrameLength / 2;
    int exc_A_old_length = PIT_MAX_MAX + L_INTERPOL - exc_A_new_length;
    E_UTIL_residu(A_old, A_old_exp, syn, old_exc_mem, exc_A_old_length);
    E_UTIL_residu(A_new, A_new_exp, &syn[exc_A_old_length],
                  &old_exc_mem[exc_A_old_length], exc_A_new_length);
  }
  C_ALLOC_SCRATCH_END(synth_buf, FIXP_DBL,
                      PIT_MAX_MAX + L_INTERPOL + M_LP_FILTER_ORDER);

  return;
}

FIXP_DBL *CLpd_ACELP_GetFreeExcMem(CAcelpStaticMem *acelp_mem, INT length) {
  FDK_ASSERT(length <= PIT_MAX_MAX + L_INTERPOL);
  return acelp_mem->old_exc_mem;
}

INT CLpd_AcelpRead(HANDLE_FDK_BITSTREAM hBs, CAcelpChannelData *acelp,
                   INT acelp_core_mode, INT coreCoderFrameLength,
                   INT i_offset) {
  int nb_subfr = coreCoderFrameLength / L_DIV;
  const UCHAR *num_acb_index_bits =
      (nb_subfr == 4) ? num_acb_idx_bits_table[0] : num_acb_idx_bits_table[1];
  int nbits;
  int error = 0;

  const int PIT_MIN = PIT_MIN_12k8 + i_offset;
  const int PIT_FR2 = PIT_FR2_12k8 - i_offset;
  const int PIT_FR1 = PIT_FR1_12k8;
  const int PIT_MAX = PIT_MAX_12k8 + (6 * i_offset);
  int T0, T0_frac, T0_min = 0, T0_max;

  if (PIT_MAX > PIT_MAX_MAX) {
    error = AAC_DEC_DECODE_FRAME_ERROR;
    goto bail;
  }

  acelp->acelp_core_mode = acelp_core_mode;

  nbits = MapCoreMode2NBits(acelp_core_mode);

  acelp->mean_energy = FDKreadBits(hBs, 2);

  for (int sfr = 0; sfr < nb_subfr; sfr++) {
    error = DecodePitchLag(hBs, num_acb_index_bits[sfr], PIT_MIN, PIT_FR2,
                           PIT_FR1, PIT_MAX, &T0, &T0_frac, &T0_min, &T0_max);
    if (error) {
      goto bail;
    }
    acelp->T0[sfr] = (USHORT)T0;
    acelp->T0_frac[sfr] = (UCHAR)T0_frac;
    acelp->ltp_filtering_flag[sfr] = FDKreadBits(hBs, 1);
    switch (nbits) {
      case 12:
        acelp->icb_index[sfr][0] = FDKreadBits(hBs, 1);
        acelp->icb_index[sfr][1] = FDKreadBits(hBs, 5);
        acelp->icb_index[sfr][2] = FDKreadBits(hBs, 1);
        acelp->icb_index[sfr][3] = FDKreadBits(hBs, 5);
        break;
      case 16:
        acelp->icb_index[sfr][0] = FDKreadBits(hBs, 1);
        acelp->icb_index[sfr][1] = FDKreadBits(hBs, 5);
        acelp->icb_index[sfr][2] = FDKreadBits(hBs, 5);
        acelp->icb_index[sfr][3] = FDKreadBits(hBs, 5);
        break;
      case 20:
        acelp->icb_index[sfr][0] = FDKreadBits(hBs, 5);
        acelp->icb_index[sfr][1] = FDKreadBits(hBs, 5);
        acelp->icb_index[sfr][2] = FDKreadBits(hBs, 5);
        acelp->icb_index[sfr][3] = FDKreadBits(hBs, 5);
        break;
      case 28:
        acelp->icb_index[sfr][0] = FDKreadBits(hBs, 9);
        acelp->icb_index[sfr][1] = FDKreadBits(hBs, 9);
        acelp->icb_index[sfr][2] = FDKreadBits(hBs, 5);
        acelp->icb_index[sfr][3] = FDKreadBits(hBs, 5);
        break;
      case 36:
        acelp->icb_index[sfr][0] = FDKreadBits(hBs, 9);
        acelp->icb_index[sfr][1] = FDKreadBits(hBs, 9);
        acelp->icb_index[sfr][2] = FDKreadBits(hBs, 9);
        acelp->icb_index[sfr][3] = FDKreadBits(hBs, 9);
        break;
      case 44:
        acelp->icb_index[sfr][0] = FDKreadBits(hBs, 13);
        acelp->icb_index[sfr][1] = FDKreadBits(hBs, 13);
        acelp->icb_index[sfr][2] = FDKreadBits(hBs, 9);
        acelp->icb_index[sfr][3] = FDKreadBits(hBs, 9);
        break;
      case 52:
        acelp->icb_index[sfr][0] = FDKreadBits(hBs, 13);
        acelp->icb_index[sfr][1] = FDKreadBits(hBs, 13);
        acelp->icb_index[sfr][2] = FDKreadBits(hBs, 13);
        acelp->icb_index[sfr][3] = FDKreadBits(hBs, 13);
        break;
      case 64:
        acelp->icb_index[sfr][0] = FDKreadBits(hBs, 2);
        acelp->icb_index[sfr][1] = FDKreadBits(hBs, 2);
        acelp->icb_index[sfr][2] = FDKreadBits(hBs, 2);
        acelp->icb_index[sfr][3] = FDKreadBits(hBs, 2);
        acelp->icb_index[sfr][4] = FDKreadBits(hBs, 14);
        acelp->icb_index[sfr][5] = FDKreadBits(hBs, 14);
        acelp->icb_index[sfr][6] = FDKreadBits(hBs, 14);
        acelp->icb_index[sfr][7] = FDKreadBits(hBs, 14);
        break;
      default:
        FDK_ASSERT(0);
        break;
    }
    acelp->gains[sfr] = FDKreadBits(hBs, 7);
  }

bail:
  return error;
}
