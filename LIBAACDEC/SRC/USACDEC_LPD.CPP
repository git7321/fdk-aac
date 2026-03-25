/**************************** AAC decoder library ******************************

   Author(s):   Manuel Jander

   Description: USAC Linear Prediction Domain coding

*******************************************************************************/

#include "usacdec_lpd.h"
#include "usacdec_rom.h"
#include "usacdec_fac.h"
#include "usacdec_lpc.h"
#include "FDK_tools_rom.h"
#include "fft.h"
#include "mdct.h"
#include "usacdec_acelp.h"
#include "overlapadd.h"
#include "conceal.h"
#include "block.h"

#define SF_PITCH_TRACK 6
#define SF_GAIN 3
#define MIN_VAL FL2FXCONST_DBL(0.0f)
#define MAX_VAL (FIXP_DBL) MAXVAL_DBL

#include "ac_arith_coder.h"

void filtLP(const FIXP_DBL *syn, PCM_DEC *syn_out, FIXP_DBL *noise,
            const FIXP_SGL *filt, const INT aacOutDataHeadroom, INT stop,
            int len) {
  INT i, j;
  FIXP_DBL tmp;

  FDK_ASSERT((aacOutDataHeadroom - 1) >= -(MDCT_OUTPUT_SCALE));

  for (i = 0; i < stop; i++) {
    tmp = fMultDiv2(noise[i], filt[0]);
    for (j = 1; j <= len; j++) {
      tmp += fMult((noise[i - j] >> 1) + (noise[i + j] >> 1), filt[j]);
    }
    syn_out[i] = (PCM_DEC)(
        IMDCT_SCALE((syn[i] >> 1) - (tmp >> 1), aacOutDataHeadroom - 1));
  }
}

void bass_pf_1sf_delay(
    FIXP_DBL *syn,
    const INT *T_sf,
    FIXP_DBL *pit_gain,
    const int frame_length,
    const INT l_frame,
    const INT l_next,
    PCM_DEC *synth_out,
    const INT aacOutDataHeadroom,
    FIXP_DBL mem_bpf[])
{
  INT i, sf, i_subfr, T, T2, lg;

  FIXP_DBL tmp, ener, corr, gain;
  FIXP_DBL *noise, *noise_in;
  FIXP_DBL
  noise_buf[L_FILT + (2 * L_SUBFR)];
  const FIXP_DBL *x, *y;

  {
    noise = noise_buf + L_FILT;
    noise_in = noise_buf + L_FILT + L_SUBFR;
    scaleValues(mem_bpf, (L_FILT + L_SUBFR), 1);
  }

  int gain_exp = 17;

  sf = 0;
  for (i_subfr = 0; i_subfr < l_frame; i_subfr += L_SUBFR, sf++) {
    T = T_sf[sf];
    gain = pit_gain[sf];

    if (gain > (FIXP_DBL)(1 << 14)) gain = (FIXP_DBL)(1 << 14);

    if (gain < (FIXP_DBL)0) gain = (FIXP_DBL)0;

    if (gain > (FIXP_DBL)0) {
      T2 = T >> 1;
      x = &syn[i_subfr - L_EXTRA];
      y = &syn[i_subfr - T2 - L_EXTRA];

      ener = (FIXP_DBL)0;
      corr = (FIXP_DBL)0;
      tmp = (FIXP_DBL)0;

      int headroom_x = getScalefactor(x, L_SUBFR + L_EXTRA);
      int headroom_y = getScalefactor(y, L_SUBFR + L_EXTRA);

      int width_shift = 7;

      for (i = 0; i < (L_SUBFR + L_EXTRA); i++) {
        ener += fPow2Div2((x[i] << headroom_x)) >> width_shift;
        corr += fMultDiv2((x[i] << headroom_x), (y[i] << headroom_y)) >>
                width_shift;
        tmp += fPow2Div2((y[i] << headroom_y)) >> width_shift;
      }

      int exp_ener = ((17 - headroom_x) << 1) + width_shift + 1;
      int exp_corr = (17 - headroom_x) + (17 - headroom_y) + width_shift + 1;
      int exp_tmp = ((17 - headroom_y) << 1) + width_shift + 1;

      FIXP_DBL point_zero_one = (FIXP_DBL)0x51eb851f;
      int diff;
      ener = fAddNorm(ener, exp_ener, point_zero_one, -6, &exp_ener);
      corr = fAddNorm(corr, exp_corr, point_zero_one, -6, &exp_corr);
      tmp = fAddNorm(tmp, exp_tmp, point_zero_one, -6, &exp_tmp);

      INT s1, s2;
      s1 = CntLeadingZeros(ener) - 1;
      s2 = CntLeadingZeros(tmp) - 1;

      FIXP_DBL ener_by_tmp = fMultDiv2(ener << s1, tmp << s2);
      int ener_by_tmp_exp = (exp_ener - s1) + (exp_tmp - s2) + 1;

      if (ener_by_tmp_exp & 1) {
        ener_by_tmp <<= 1;
        ener_by_tmp_exp -= 1;
      }

      int temp_exp = 0;

      FIXP_DBL temp1 = invSqrtNorm2(ener_by_tmp, &temp_exp);

      int temp1_exp = temp_exp - (ener_by_tmp_exp >> 1);

      FIXP_DBL tmp_result = fMult(corr, temp1);

      int tmp_result_exp = exp_corr + temp1_exp;

      diff = tmp_result_exp - 0;
      FIXP_DBL point95 = FL2FXCONST_DBL(0.95f);
      if (diff >= 0) {
        diff = fMin(diff, 31);
        point95 = FL2FXCONST_DBL(0.95f) >> diff;
      } else {
        diff = fMax(diff, -31);
        tmp_result >>= (-diff);
      }

      if (tmp_result > point95) T = T2;

      lg = l_frame + l_next - T - i_subfr;

      if (lg > L_SUBFR)
        lg = L_SUBFR;
      else if (lg < 0)
        lg = 0;

      if (lg > 0) {
        FIXP_DBL tmp1;

        s1 = getScalefactor(&syn[i_subfr], lg);
        s2 = getScalefactor(&syn[i_subfr + T], lg);
        INT s = fixMin(s1, s2);

        tmp = (FIXP_DBL)0;
        ener = (FIXP_DBL)0;
        for (i = 0; i < lg; i++) {
          tmp += fPow2Div2(syn[i + i_subfr] << s1) >> (SF_PITCH_TRACK);
          ener += fPow2Div2(syn[i + i_subfr + T] << s2) >> (SF_PITCH_TRACK);
        }
        tmp = tmp >> fMin(DFRACT_BITS - 1, (2 * (s1 - s)));
        ener = ener >> fMin(DFRACT_BITS - 1, (2 * (s2 - s)));

        if (ener <= (FIXP_DBL)0) {
          ener = (FIXP_DBL)1;
        }
        if (tmp <= (FIXP_DBL)0) {
          tmp = (FIXP_DBL)1;
        }
        FDK_ASSERT(ener > (FIXP_DBL)0);

        int result_e = 0;
        tmp1 = fDivNorm(tmp, ener, &result_e);
        if (result_e & 1) {
          tmp1 >>= 1;
          result_e += 1;
        }
        tmp = sqrtFixp(tmp1);
        result_e >>= 1;

        gain_exp = 17;

        diff = result_e - gain_exp;

        FIXP_DBL gain1 = gain;

        if (diff >= 0) {
          diff = fMin(diff, 31);
          gain1 >>= diff;
        } else {
          result_e += (-diff);
          diff = fMax(diff, -31);
          tmp >>= (-diff);
        }

        if (tmp < gain1) {
          gain = tmp;
          gain_exp = result_e;
        }
      }

      diff = gain_exp - 17;
      if (diff >= 0) {
        gain <<= diff;
      } else {
        gain >>= (-diff);
      }

      s1 = CntLeadingZeros(gain) - 1;
      s1 -= 16;

      FIXP_SGL gainSGL = FX_DBL2FX_SGL(gain << 16);

      gainSGL = gainSGL << s1;

      {
        for (i = 0; i < lg; i++) {
          noise_in[i] =
              scaleValue(fMultDiv2(gainSGL, (syn[i + i_subfr] >> 1) -
                                                (syn[i + i_subfr - T] >> 2) -
                                                (syn[i + i_subfr + T] >> 2)),
                         2 - s1);
        }

        for (i = lg; i < L_SUBFR; i++) {
          noise_in[i] =
              scaleValue(fMultDiv2(gainSGL, (syn[i + i_subfr] >> 1) -
                                                (syn[i + i_subfr - T] >> 1)),
                         2 - s1);
        }
      }
    } else {
      FDKmemset(noise_in, (FIXP_DBL)0, L_SUBFR * sizeof(FIXP_DBL));
    }

    {
      FDKmemcpy(noise_buf, mem_bpf, (L_FILT + L_SUBFR) * sizeof(FIXP_DBL));

      FDKmemcpy(mem_bpf, noise_buf + L_SUBFR,
                (L_FILT + L_SUBFR) * sizeof(FIXP_DBL));
    }

    {
      filtLP(&syn[i_subfr - L_SUBFR], &synth_out[i_subfr], noise,
             fdk_dec_filt_lp, aacOutDataHeadroom, L_SUBFR, L_FILT);
    }
  }

  {
    scaleValues(mem_bpf, (L_FILT + L_SUBFR), -1);
    scaleValuesSaturate(
        (PCM_DEC *)&synth_out[l_frame], (FIXP_DBL *)&syn[l_frame - L_SUBFR],
        (frame_length - l_frame), MDCT_OUT_HEADROOM - aacOutDataHeadroom);
  }

  return;
}

#define ALFDPOW2_SCALE 3
void CLpd_AdaptLowFreqDeemph(FIXP_DBL x[], int lg, FIXP_DBL alfd_gains[],
                             INT s) {
  {
    int i, j, k, i_max;
    FIXP_DBL max, fac;
    FIXP_DBL tmp_pow2[32];

    s = s * 2 + ALFDPOW2_SCALE;
    s = fMin(31, s);

    k = 8;
    i_max = lg / 4;

    max = FL2FX_DBL(0.01f) >> s;
    for (i = 0; i < i_max; i += k) {
      FIXP_DBL tmp;

      tmp = FIXP_DBL(0);
      FIXP_DBL *pX = &x[i];

      j = 8;
      do {
        FIXP_DBL x0 = *pX++;
        FIXP_DBL x1 = *pX++;
        x0 = fPow2Div2(x0);
        x1 = fPow2Div2(x1);
        tmp = tmp + (x0 >> (ALFDPOW2_SCALE - 1));
        tmp = tmp + (x1 >> (ALFDPOW2_SCALE - 1));
      } while ((j = j - 2) != 0);
      tmp = fMax(tmp, (FL2FX_DBL(0.01f) >> s));
      tmp_pow2[i >> 3] = tmp;
      if (tmp > max) {
        max = tmp;
      }
    }

    fac = FL2FX_DBL(0.1f) >> 1;
    for (i = 0; i < i_max; i += k) {
      FIXP_DBL tmp;
      INT shifti;

      tmp = tmp_pow2[i >> 3];

      {
        INT sd;

        if (tmp != (FIXP_DBL)0) {
          tmp = fDivNorm(max, tmp, &sd);
          if (sd & 1) {
            sd++;
            tmp >>= 1;
          }
        } else {
          tmp = (FIXP_DBL)MAXVAL_DBL;
          sd = 0;
        }
        tmp = invSqrtNorm2(tmp, &shifti);
        tmp = scaleValue(tmp, shifti - 1 - (sd / 2));
      }
      if (tmp > fac) {
        fac = tmp;
      }
      FIXP_DBL *pX = &x[i];

      j = 8;
      do {
        FIXP_DBL x0 = pX[0];
        FIXP_DBL x1 = pX[1];
        x0 = fMultDiv2(x0, fac);
        x1 = fMultDiv2(x1, fac);
        x0 = x0 << 2;
        x1 = x1 << 2;
        *pX++ = x0;
        *pX++ = x1;

      } while ((j = j - 2) != 0);
      *alfd_gains++ = fac;
    }
  }
}

#define NSHAPE_SCALE (4)

#define LPC2MDCT_CALC (1)
#define LPC2MDCT_GAIN_LOAD (2)
#define LPC2MDCT_GAIN_SAVE (4)
#define LPC2MDCT_APPLY_NSHAPE (8)

void lpc2mdctAndNoiseShaping(FIXP_DBL *r, SHORT *pScale, const INT lg,
                             const INT fdns_npts, const FIXP_LPC *A1,
                             const INT A1_exp, const FIXP_LPC *A2,
                             const INT A2_exp) {
  FIXP_DBL *tmp2 = NULL;
  FIXP_DBL rr_minus_one;
  int i, k, s, step;

  C_AALLOC_SCRATCH_START(tmp1, FIXP_DBL, FDNS_NPTS * 8)

  {
    tmp2 = tmp1 + fdns_npts * 4;

    FIXP_DBL f = FL2FXCONST_DBL(0.92f);

    const FIXP_STP *SinTab;
    int k_step;
    switch (fdns_npts) {
      case 64:
        SinTab = SineTable512;
        k_step = (512 / 64);
        FDK_ASSERT(512 >= 64);
        break;
      case 48:
        SinTab = SineTable384;
        k_step = 384 / 48;
        FDK_ASSERT(384 >= 48);
        break;
      default:
        FDK_ASSERT(0);
        return;
    }

    for (i = 0, k = k_step; i < M_LP_FILTER_ORDER; i++, k += k_step) {
      FIXP_STP cs = SinTab[k];
      FIXP_DBL wA1, wA2;

      wA1 = fMult(A1[i], f);
      wA2 = fMult(A2[i], f);

      tmp1[2 + i * 2] = fMult(wA1, cs.v.re);
      tmp2[2 + i * 2] = fMult(wA2, cs.v.re);
      tmp1[3 + i * 2] = -fMult(wA1, cs.v.im);
      tmp2[3 + i * 2] = -fMult(wA2, cs.v.im);

      f = fMult(f, FL2FXCONST_DBL(0.92f));
    }

    int A1_exp_fix = fMax(3, A1_exp + 2);
    int A2_exp_fix = fMax(3, A2_exp + 2);

    tmp1[0] = (FIXP_DBL)(INT)((ULONG)0x80000000 >> A1_exp_fix);
    tmp2[0] = (FIXP_DBL)(INT)((ULONG)0x80000000 >> A2_exp_fix);

    tmp1[1] = tmp2[1] = (FIXP_DBL)0;

    FDKmemclear(
        tmp1 + 2 * (M_LP_FILTER_ORDER + 1),
        2 * (fdns_npts * 2 - (M_LP_FILTER_ORDER + 1)) * sizeof(FIXP_DBL));
    FDKmemclear(
        tmp2 + 2 * (M_LP_FILTER_ORDER + 1),
        2 * (fdns_npts * 2 - (M_LP_FILTER_ORDER + 1)) * sizeof(FIXP_DBL));

    scaleValues(&tmp1[2], (2 * M_LP_FILTER_ORDER), (A1_exp - A1_exp_fix));
    scaleValues(&tmp2[2], (2 * M_LP_FILTER_ORDER), (A2_exp - A2_exp_fix));

    INT s2;
    s = A1_exp_fix;
    s2 = A2_exp_fix;

    fft(2 * fdns_npts, tmp1, &s);
    fft(2 * fdns_npts, tmp2, &s2);

    if (s > s2) {
      scaleValues(tmp2, 2 * fdns_npts, s2 - s);
      s2 = s;
    } else if (s < s2) {
      scaleValues(tmp1, 2 * fdns_npts, s - s2);
      s = s2;
    }

    FDK_ASSERT(s == s2);
  }

  step = lg / fdns_npts;
  rr_minus_one = (FIXP_DBL)0;

  for (k = 0; k < fdns_npts; k++) {
    FIXP_DBL g1, g2, inv_g1_g2, a, b;
    INT inv_g1_g2_e;
    int g_e, shift;

    {
      FIXP_DBL real, imag;
      int si1, si2, sInput;

      real = tmp1[k * 2];
      imag = tmp1[k * 2 + 1];
      sInput = fMax(fMin(fNorm(real), fNorm(imag)) - 1, 0);
      real <<= sInput;
      imag <<= sInput;
      g1 = invSqrtNorm2(fPow2(real) + fPow2(imag), &si1);
      si1 += sInput;

      real = tmp2[k * 2];
      imag = tmp2[k * 2 + 1];
      sInput = fMax(fMin(fNorm(real), fNorm(imag)) - 1, 0);
      real <<= sInput;
      imag <<= sInput;
      g2 = invSqrtNorm2(fPow2(real) + fPow2(imag), &si2);
      si2 += sInput;

      if (si1 > si2) {
        g2 >>= si1 - si2;
        g_e = si1 - s;
      } else {
        g1 >>= si2 - si1;
        g_e = si2 - s;
      }
    }

    FDK_ASSERT(g1 >= (FIXP_DBL)0);
    FDK_ASSERT(g2 >= (FIXP_DBL)0);

    {
      inv_g1_g2 = (g1 >> 1) + (g2 >> 1);
      if (inv_g1_g2 != (FIXP_DBL)0) {
        inv_g1_g2 = fDivNorm(FL2FXCONST_DBL(0.5f), inv_g1_g2, &inv_g1_g2_e);
        inv_g1_g2_e = inv_g1_g2_e - g_e;
      } else {
        inv_g1_g2 = (FIXP_DBL)MAXVAL_DBL;
        inv_g1_g2_e = 0;
      }

      if (g_e < 0) {
        a = scaleValue(fMult(fMult(g1, g2), inv_g1_g2), g_e);
        b = fMult(g2 - g1, inv_g1_g2);
        shift = g_e + inv_g1_g2_e + 1 - NSHAPE_SCALE;
      } else {
        a = fMult(fMult(g1, g2), inv_g1_g2);
        b = scaleValue(fMult(g2 - g1, inv_g1_g2), -g_e);
        shift = (g_e + g_e) + inv_g1_g2_e + 1 - NSHAPE_SCALE;
      }

      for (i = k * step; i < (k + 1) * step; i++) {
        FIXP_DBL tmp;

        tmp = fMult(a, r[i]);
        tmp += scaleValue(fMultDiv2(b, rr_minus_one), NSHAPE_SCALE);
        tmp = scaleValueSaturate(tmp, shift);
        rr_minus_one = tmp;
        r[i] = tmp;
      }
    }
  }

  { *pScale += NSHAPE_SCALE; }

  C_AALLOC_SCRATCH_END(tmp1, FIXP_DBL, FDNS_NPTS * 8)
}

static FIXP_DBL calcEnergy(const FIXP_DBL *r, const SHORT rs, const INT lg,
                           INT *rms_e) {
  int headroom = getScalefactor(r, lg);

  FIXP_DBL rms_m = 0;

  INT shift = (INT)(fNormz((FIXP_DBL)lg));
  shift = 31 - shift;

  const FIXP_DBL value0_01 = 0x51eb851e;
  const INT value0_01_exp = -6;

  *rms_e = ((rs - headroom) << 1) + shift + 1;

  INT delta = *rms_e - value0_01_exp;
  if (delta > 0) {
    delta = fMin(31, delta);
    rms_m = value0_01 >> delta;
  } else {
    rms_m = value0_01;
    *rms_e = value0_01_exp;
    shift = shift - delta;
    shift = fMin(31, shift);
  }

  for (int i = 0; i < lg; i++) {
    rms_m += fPow2Div2(r[i] << headroom) >> shift;
  }

  return rms_m;
}

static void calcTCXGain(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                        FIXP_DBL *r, FIXP_DBL rms_m, INT rms_e, const INT frame,
                        const INT lg) {
  if ((rms_m != (FIXP_DBL)0)) {
    FIXP_DBL tcx_gain_m;
    INT tcx_gain_e;

    CLpd_DecodeGain(&tcx_gain_m, &tcx_gain_e,
                    pAacDecoderChannelInfo->pDynData->specificTo.usac
                        .tcx_global_gain[frame]);

    if (rms_e & 1) {
      rms_m >>= 1;
      rms_e++;
    }

    {
      FIXP_DBL fx_lg;
      INT fx_lg_e, s;
      INT inv_e;

      s = fNorm((FIXP_DBL)lg);
      fx_lg = (FIXP_DBL)lg << s;
      fx_lg_e = DFRACT_BITS - 1 - s;
      rms_m = invSqrtNorm2(rms_m, &inv_e);
      rms_m = fMult(rms_m, fx_lg);
      rms_e = inv_e - (rms_e >> 1) + fx_lg_e;
    }

    {
      int s = fNorm(tcx_gain_m);
      tcx_gain_m = tcx_gain_m << s;
      tcx_gain_e -= s;
    }

    tcx_gain_m = fMultDiv2(tcx_gain_m, rms_m);
    tcx_gain_e = tcx_gain_e + rms_e;

    {
      { tcx_gain_e += 1; }
    }

    pAacDecoderChannelInfo->data.usac.tcx_gain[frame] = tcx_gain_m;
    pAacDecoderChannelInfo->data.usac.tcx_gain_e[frame] = tcx_gain_e;

    pAacDecoderChannelInfo->specScale[frame] += tcx_gain_e;
  }
}

void CLpd_FdnsDecode(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                     CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
                     FIXP_DBL r[], const INT lg, const INT frame, SHORT *pScale,
                     const FIXP_LPC A1[M_LP_FILTER_ORDER], const INT A1_exp,
                     const FIXP_LPC A2[M_LP_FILTER_ORDER], const INT A2_exp,
                     FIXP_DBL pAlfdGains[LFAC / 4], const INT fdns_npts) {
  CLpd_AdaptLowFreqDeemph(r, lg, pAlfdGains, *pScale);

  FIXP_DBL rms_m = (FIXP_DBL)0;
  INT rms_e = 0;
  {
    rms_m = calcEnergy(r, *pScale, lg, &rms_e);
  }

  calcTCXGain(pAacDecoderChannelInfo, r, rms_m, rms_e, frame, lg);

  lpc2mdctAndNoiseShaping(r, pScale, lg, fdns_npts, A1, A1_exp, A2, A2_exp);
}

static int find_mpitch(FIXP_DBL xri[], int lg) {
  FIXP_DBL max, pitch;
  INT pitch_e;
  int i, n;

  max = (FIXP_DBL)0;
  n = 2;

  for (i = 2; i < (lg >> 4); i += 2) {
    FIXP_DBL tmp = fPow2Div2(xri[i]) + fPow2Div2(xri[i + 1]);
    if (tmp > max) {
      max = tmp;
      n = i;
    }
  }

  pitch = fDivNorm((FIXP_DBL)lg << 1, (FIXP_DBL)n, &pitch_e);
  pitch >>= fixMax(0, DFRACT_BITS - 1 - pitch_e - 16);

  if (pitch >= (FIXP_DBL)((256 << 16) - 1)) {
    n = 256;
  } else {
    FIXP_DBL mpitch = pitch;
    while (mpitch < (FIXP_DBL)(255 << 16)) {
      mpitch += pitch;
    }
    n = (int)(mpitch - pitch) >> 16;
  }

  return (n);
}

static const int lg_table_ccfl[2][4] = {
    {256, 256, 512, 1024},
    {192, 192, 384, 768}
};

static void CLpd_TcxDecode(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo, UINT flags,
    int mod, int last_mod, int frame, int frameOk) {
  FIXP_DBL *pAlfd_gains = pAacDecoderStaticChannelInfo->last_alfd_gains;
  ULONG *pSeed = &pAacDecoderStaticChannelInfo->nfRandomSeed;
  int lg = (pAacDecoderChannelInfo->granuleLength == 128)
               ? lg_table_ccfl[0][mod + 0]
               : lg_table_ccfl[1][mod + 0];
  int next_frame = frame + (1 << (mod - 1));
  int isFullBandLpd = 0;

  {
    FIXP_DBL noise_level;
    FIXP_DBL *coeffs =
        SPEC_TCX(pAacDecoderChannelInfo->pSpectralCoefficient, frame,
                 pAacDecoderChannelInfo->granuleLength, isFullBandLpd);
    int scale = pAacDecoderChannelInfo->specScale[frame];
    int i, nfBgn, nfEnd;
    UCHAR tcx_noise_factor = pAacDecoderChannelInfo->pDynData->specificTo.usac
                                 .tcx_noise_factor[frame];

    pAacDecoderStaticChannelInfo->last_tcx_pitch = find_mpitch(coeffs, lg);

    if (frameOk) {
      pAacDecoderStaticChannelInfo->last_tcx_noise_factor = tcx_noise_factor;
    } else {
      tcx_noise_factor = pAacDecoderStaticChannelInfo->last_tcx_noise_factor;
    }

    noise_level =
        (FIXP_DBL)((LONG)FL2FXCONST_DBL(0.0625f) * (8 - tcx_noise_factor));
    noise_level = scaleValue(noise_level, -scale);

    const FIXP_DBL neg_noise_level = -noise_level;

    {
      nfBgn = lg / 6;
      nfEnd = lg;
    }

    for (i = nfBgn; i < nfEnd - 7; i += 8) {
      LONG tmp;

      tmp = coeffs[i + 0] | coeffs[i + 1] | coeffs[i + 2] | coeffs[i + 3] |
            coeffs[i + 4] | coeffs[i + 5] | coeffs[i + 6] | coeffs[i + 7];

      if (tmp == 0) {
        for (int k = i; k < i + 8; k++) {
          UsacRandomSign(pSeed) ? (coeffs[k] = neg_noise_level)
                                : (coeffs[k] = noise_level);
        }
      }
    }
    if ((nfEnd - i) >
        0) {
      LONG tmp = (LONG)coeffs[i];
      int k;

      FDK_ASSERT((nfEnd - i) < 8);
      for (k = 1; k < (nfEnd - i); k++) {
        tmp |= (LONG)coeffs[i + k];
      }
      if (tmp == 0) {
        for (k = i; k < nfEnd; k++) {
          UsacRandomSign(pSeed) ? (coeffs[k] = neg_noise_level)
                                : (coeffs[k] = noise_level);
        }
      }
    }
  }

  {
    if (last_mod == 0) {
      E_LPC_f_lsp_a_conversion(
          pAacDecoderChannelInfo->data.usac.lsp_coeff[frame],
          pAacDecoderChannelInfo->data.usac.lp_coeff[frame],
          &pAacDecoderChannelInfo->data.usac.lp_coeff_exp[frame]);
    }

    E_LPC_f_lsp_a_conversion(
        pAacDecoderChannelInfo->data.usac.lsp_coeff[next_frame],
        pAacDecoderChannelInfo->data.usac.lp_coeff[next_frame],
        &pAacDecoderChannelInfo->data.usac.lp_coeff_exp[next_frame]);

    CLpd_FdnsDecode(
        pAacDecoderChannelInfo, pAacDecoderStaticChannelInfo,
        SPEC_TCX(pAacDecoderChannelInfo->pSpectralCoefficient, frame,
                 pAacDecoderChannelInfo->granuleLength, isFullBandLpd),
        lg, frame, pAacDecoderChannelInfo->specScale + frame,
        pAacDecoderChannelInfo->data.usac.lp_coeff[frame],
        pAacDecoderChannelInfo->data.usac.lp_coeff_exp[frame],
        pAacDecoderChannelInfo->data.usac.lp_coeff[next_frame],
        pAacDecoderChannelInfo->data.usac.lp_coeff_exp[next_frame], pAlfd_gains,
        pAacDecoderChannelInfo->granuleLength / 2
    );
  }
}

static AAC_DECODER_ERROR CLpd_TCX_Read(
    HANDLE_FDK_BITSTREAM hBs, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo, int lg,
    int first_tcx_flag, int frame, UINT flags) {
  AAC_DECODER_ERROR errorAAC = AAC_DEC_OK;
  ARITH_CODING_ERROR error = ARITH_CODER_OK;
  FIXP_DBL *pSpec;
  int arith_reset_flag = 0;
  int isFullBandLpd = 0;

  pSpec = SPEC_TCX(pAacDecoderChannelInfo->pSpectralCoefficient, frame,
                   pAacDecoderChannelInfo->granuleLength, isFullBandLpd);
  {
    pAacDecoderChannelInfo->pDynData->specificTo.usac.tcx_noise_factor[frame] =
        FDKreadBits(hBs, 3);
  }
  pAacDecoderChannelInfo->pDynData->specificTo.usac.tcx_global_gain[frame] =
      FDKreadBits(hBs, 7);

  if (first_tcx_flag) {
    if (flags & AC_INDEP) {
      arith_reset_flag = 1;
    } else {
      arith_reset_flag = FDKreadBits(hBs, 1);
    }
  }

  error = CArco_DecodeArithData(pAacDecoderStaticChannelInfo->hArCo, hBs, pSpec,
                                lg, lg, arith_reset_flag);

  {
    int scale = getScalefactor(pSpec, lg) - 2;

    scaleValues(pSpec, lg, scale);
    scale = DFRACT_BITS - 1 - scale;

    pAacDecoderChannelInfo->specScale[frame] = scale;
  }

  if (error == ARITH_CODER_ERROR) errorAAC = AAC_DEC_UNKNOWN;

  return errorAAC;
}

static AAC_DECODER_ERROR CLpd_ReadAndMapLpdModeToModArray(
    UCHAR mod[4], HANDLE_FDK_BITSTREAM hBs, UINT elFlags) {
  int lpd_mode;

  {
    lpd_mode = FDKreadBits(hBs, 5);

    if (lpd_mode > 25 || lpd_mode < 0) {
      return AAC_DEC_PARSE_ERROR;
    }

    switch (lpd_mode) {
      case 25:
        mod[0] = mod[1] = mod[2] = mod[3] = 3;
        break;
      case 24:
        mod[0] = mod[1] = mod[2] = mod[3] = 2;
        break;
      default:
        switch (lpd_mode >> 2) {
          case 4:
            mod[0] = mod[1] = 2;
            mod[2] = (lpd_mode & 1) ? 1 : 0;
            mod[3] = (lpd_mode & 2) ? 1 : 0;
            break;
          case 5:
            mod[2] = mod[3] = 2;
            mod[0] = (lpd_mode & 1) ? 1 : 0;
            mod[1] = (lpd_mode & 2) ? 1 : 0;
            break;
          default:
            mod[0] = (lpd_mode & 1) ? 1 : 0;
            mod[1] = (lpd_mode & 2) ? 1 : 0;
            mod[2] = (lpd_mode & 4) ? 1 : 0;
            mod[3] = (lpd_mode & 8) ? 1 : 0;
            break;
        }
        break;
    }
  }
  return AAC_DEC_OK;
}

static void CLpd_Reset(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    int keep_past_signal) {
  int i;

  if (!keep_past_signal) {
    FDKmemclear(pAacDecoderStaticChannelInfo->old_synth,
                sizeof(pAacDecoderStaticChannelInfo->old_synth));
  }

  for (i = 0; i < M_LP_FILTER_ORDER; i++) {
    pAacDecoderStaticChannelInfo->lpc4_lsf[i] = fdk_dec_lsf_init[i];
  }

  FDKmemclear(pAacDecoderStaticChannelInfo->mem_bpf,
              sizeof(pAacDecoderStaticChannelInfo->mem_bpf));

  pAacDecoderStaticChannelInfo->old_bpf_control_info = 0;
  for (i = 0; i < SYN_SFD; i++) {
    pAacDecoderStaticChannelInfo->old_T_pf[i] = 64;
    pAacDecoderStaticChannelInfo->old_gain_pf[i] = (FIXP_DBL)0;
  }

  CLpd_AcelpReset(&pAacDecoderStaticChannelInfo->acelp);

  pAacDecoderStaticChannelInfo->last_lpc_lost = 0;
  pAacDecoderStaticChannelInfo->last_tcx_pitch = L_DIV;
  pAacDecoderStaticChannelInfo->numLostLpdFrames = 0;
}

AAC_DECODER_ERROR CLpdChannelStream_Read(
    HANDLE_FDK_BITSTREAM hBs, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    const SamplingRateInfo *pSamplingRateInfo, UINT flags) {
  AAC_DECODER_ERROR error = AAC_DEC_OK;
  int first_tcx_flag;
  int k, nbDiv, fFacDataPresent, first_lpd_flag, acelp_core_mode,
      facGetMemState = 0;
  UCHAR *mod = pAacDecoderChannelInfo->data.usac.mod;
  int lpd_mode_last, prev_frame_was_lpd;
  USAC_COREMODE core_mode_last;
  const int lg_table_offset = 0;
  const int *lg_table = (pAacDecoderChannelInfo->granuleLength == 128)
                            ? &lg_table_ccfl[0][lg_table_offset]
                            : &lg_table_ccfl[1][lg_table_offset];
  int last_lpc_lost = pAacDecoderStaticChannelInfo->last_lpc_lost;

  int last_frame_ok = CConcealment_GetLastFrameOk(
      &pAacDecoderStaticChannelInfo->concealmentInfo, 1);

  INT i_offset;
  UINT samplingRate;

  samplingRate = pSamplingRateInfo->samplingRate;

  i_offset =
      (INT)(samplingRate * PIT_MIN_12k8 + (FSCALE_DENOM / 2)) / FSCALE_DENOM -
      (INT)PIT_MIN_12k8;

  if ((samplingRate < FAC_FSCALE_MIN) || (samplingRate > FAC_FSCALE_MAX)) {
    error = AAC_DEC_PARSE_ERROR;
    goto bail;
  }

  acelp_core_mode = FDKreadBits(hBs, 3);

  error = CLpd_ReadAndMapLpdModeToModArray(mod, hBs, 0);
  if (error != AAC_DEC_OK) {
    goto bail;
  }

  pAacDecoderChannelInfo->data.usac.bpf_control_info = FDKreadBit(hBs);

  prev_frame_was_lpd = FDKreadBit(hBs);
  fFacDataPresent = FDKreadBit(hBs);

  pAacDecoderChannelInfo->data.usac.core_mode_last =
      pAacDecoderStaticChannelInfo->last_core_mode;
  lpd_mode_last = pAacDecoderChannelInfo->data.usac.lpd_mode_last =
      pAacDecoderStaticChannelInfo->last_lpd_mode;

  if (prev_frame_was_lpd == 0) {
    pAacDecoderChannelInfo->data.usac.core_mode_last = FD_LONG;
    pAacDecoderChannelInfo->data.usac.lpd_mode_last = 255;
  } else {
    pAacDecoderChannelInfo->data.usac.core_mode_last = LPD;
    if (((mod[0] == 0) && fFacDataPresent) ||
        ((mod[0] != 0) && !fFacDataPresent)) {
      if (lpd_mode_last == 0) {
        pAacDecoderChannelInfo->data.usac.lpd_mode_last = 1;
      }
    } else {
      pAacDecoderChannelInfo->data.usac.lpd_mode_last = 0;
    }
  }

  first_lpd_flag = (pAacDecoderChannelInfo->data.usac.core_mode_last !=
                    LPD);
  first_tcx_flag = 1;

  if (pAacDecoderStaticChannelInfo->last_core_mode !=
      LPD) {
    CLpd_Reset(pAacDecoderChannelInfo, pAacDecoderStaticChannelInfo, 0);

    if (!last_frame_ok) {
      last_lpc_lost |= (first_lpd_flag) ? 0 : 1;
    }
  }

  core_mode_last = pAacDecoderChannelInfo->data.usac.core_mode_last;
  lpd_mode_last = pAacDecoderChannelInfo->data.usac.lpd_mode_last;

  nbDiv = NB_DIV;

  k = 0;
  while (k < nbDiv) {
    pAacDecoderChannelInfo->data.usac.fac_data[k] = NULL;

    if ((k == 0 && core_mode_last == LPD && fFacDataPresent) ||
        (lpd_mode_last == 0 && mod[k] > 0) ||
        ((lpd_mode_last != 255) && lpd_mode_last > 0 && mod[k] == 0)) {
      int err;

      pAacDecoderChannelInfo->data.usac.fac_data[k] =
          CLpd_FAC_GetMemory(pAacDecoderChannelInfo, mod, &facGetMemState);

      err = CLpd_FAC_Read(
          hBs, pAacDecoderChannelInfo->data.usac.fac_data[k],
          pAacDecoderChannelInfo->data.usac.fac_data_e,
          pAacDecoderChannelInfo->granuleLength,
          0, k);
      if (err != 0) {
        error = AAC_DEC_PARSE_ERROR;
        goto bail;
      }
    }

    if (mod[k] == 0)
    {
      int err;
      err = CLpd_AcelpRead(
          hBs, &pAacDecoderChannelInfo->data.usac.acelp[k], acelp_core_mode,
          pAacDecoderChannelInfo->granuleLength * 8,
          i_offset);
      if (err != 0) {
        error = AAC_DEC_PARSE_ERROR;
        goto bail;
      }

      lpd_mode_last = 0;
      k++;
    } else
    {
      error = CLpd_TCX_Read(hBs, pAacDecoderChannelInfo,
                            pAacDecoderStaticChannelInfo, lg_table[mod[k]],
                            first_tcx_flag, k, flags);

      lpd_mode_last = mod[k];
      first_tcx_flag = 0;
      k += 1 << (mod[k] - 1);
    }
    if (error != AAC_DEC_OK) {
      error = AAC_DEC_PARSE_ERROR;
      goto bail;
    }
  }

  {
    int err;

    err = CLpc_Read(
        hBs, pAacDecoderChannelInfo->data.usac.lsp_coeff,
        pAacDecoderStaticChannelInfo->lpc4_lsf,
        pAacDecoderChannelInfo->data.usac.lsf_adaptive_mean_cand,
        pAacDecoderChannelInfo->data.usac.aStability, mod, first_lpd_flag,
        (mod[0] & 0x3) ? 0
                       : (last_lpc_lost &&
                          pAacDecoderStaticChannelInfo->last_core_mode != LPD),
        last_frame_ok);
    if (err != 0) {
      error = AAC_DEC_PARSE_ERROR;
      goto bail;
    }
  }

  if (last_lpc_lost && !last_frame_ok) {
    int k_next;
    k = 0;
    while (k < nbDiv) {
      int i;
      k_next = k + (((mod[k] & 0x3) == 0) ? 1 : (1 << (mod[k] - 1)));
      FIXP_LPC *lsp_old = pAacDecoderChannelInfo->data.usac.lsp_coeff[k];
      FIXP_LPC *lsp_new = pAacDecoderChannelInfo->data.usac.lsp_coeff[k_next];

      for (i = 0; i < M_LP_FILTER_ORDER; i++) {
        if (lsp_new[i] < lsp_old[i]) {
          lsp_old[i] = lsp_new[i];
        }
      }
      k = k_next;
    }
  }

  if (!CConcealment_GetLastFrameOk(
          &pAacDecoderStaticChannelInfo->concealmentInfo, 1)) {
    E_LPC_f_lsp_a_conversion(
        pAacDecoderChannelInfo->data.usac.lsp_coeff[0],
        pAacDecoderChannelInfo->data.usac.lp_coeff[0],
        &pAacDecoderChannelInfo->data.usac.lp_coeff_exp[0]);
  } else if (pAacDecoderStaticChannelInfo->last_lpd_mode != 0) {
    if (pAacDecoderStaticChannelInfo->last_lpd_mode == 255) {
      E_LPC_f_lsp_a_conversion(
          pAacDecoderChannelInfo->data.usac.lsp_coeff[0],
          pAacDecoderChannelInfo->data.usac.lp_coeff[0],
          &pAacDecoderChannelInfo->data.usac.lp_coeff_exp[0]);
    } else {
      FDKmemcpy(pAacDecoderChannelInfo->data.usac.lp_coeff[0],
                pAacDecoderStaticChannelInfo->lp_coeff_old[0],
                M_LP_FILTER_ORDER * sizeof(FIXP_LPC));
      pAacDecoderChannelInfo->data.usac.lp_coeff_exp[0] =
          pAacDecoderStaticChannelInfo->lp_coeff_old_exp[0];
    }
  }

  if (fFacDataPresent && (core_mode_last != LPD)) {
    int prev_frame_was_short;

    prev_frame_was_short = FDKreadBit(hBs);

    if (prev_frame_was_short) {
      core_mode_last = pAacDecoderChannelInfo->data.usac.core_mode_last =
          FD_SHORT;
      pAacDecoderChannelInfo->data.usac.lpd_mode_last = 255;

      if ((pAacDecoderStaticChannelInfo->last_core_mode != FD_SHORT) &&
          CConcealment_GetLastFrameOk(
              &pAacDecoderStaticChannelInfo->concealmentInfo, 1)) {
        error = AAC_DEC_PARSE_ERROR;
        goto bail;
      }
    }

    pAacDecoderChannelInfo->data.usac.fac_data[0] =
        CLpd_FAC_GetMemory(pAacDecoderChannelInfo, mod, &facGetMemState);

    {
      int err;

      err = CLpd_FAC_Read(
          hBs, pAacDecoderChannelInfo->data.usac.fac_data[0],
          pAacDecoderChannelInfo->data.usac.fac_data_e,
          CLpd_FAC_getLength(core_mode_last != FD_SHORT,
                             pAacDecoderChannelInfo->granuleLength),
          1, 0);
      if (err != 0) {
        error = AAC_DEC_PARSE_ERROR;
        goto bail;
      }
    }
  }

bail:
  if (error == AAC_DEC_OK) {
    if ((pAacDecoderChannelInfo->data.usac.core_mode_last !=
         pAacDecoderStaticChannelInfo->last_core_mode) &&
        (pAacDecoderStaticChannelInfo->last_lpc_lost == 0)) {
    } else if ((pAacDecoderChannelInfo->data.usac.core_mode_last == LPD) &&
               (pAacDecoderChannelInfo->data.usac.lpd_mode_last !=
                pAacDecoderStaticChannelInfo->last_lpd_mode) &&
               (pAacDecoderStaticChannelInfo->last_lpc_lost == 0)) {
    }
  }

  return error;
}

void CLpdChannelStream_Decode(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo,
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo, UINT flags) {
  UCHAR *mod = pAacDecoderChannelInfo->data.usac.mod;
  int k;
  UCHAR last_lpd_mode;
  int nbDiv = NB_DIV;

  k = 0;
  last_lpd_mode =
      pAacDecoderChannelInfo->data.usac
          .lpd_mode_last;
  while (k < nbDiv) {
    if (mod[k] == 0) {
      if (last_lpd_mode > 0 && last_lpd_mode != 255 &&
          pAacDecoderChannelInfo->data.usac.fac_data[k]) {
        CFac_ApplyGains(pAacDecoderChannelInfo->data.usac.fac_data[k],
                        pAacDecoderChannelInfo->granuleLength,
                        pAacDecoderStaticChannelInfo->last_tcx_gain,
                        pAacDecoderStaticChannelInfo->last_alfd_gains,
                        (last_lpd_mode < 4) ? last_lpd_mode : 3);

        pAacDecoderChannelInfo->data.usac.fac_data_e[k] +=
            pAacDecoderStaticChannelInfo->last_tcx_gain_e;
      }
    } else {
      /* TCX */
      CLpd_TcxDecode(pAacDecoderChannelInfo, pAacDecoderStaticChannelInfo,
                     flags, mod[k], last_lpd_mode, k, 1
      );

      pAacDecoderStaticChannelInfo->last_tcx_gain =
          pAacDecoderChannelInfo->data.usac.tcx_gain[k];
      pAacDecoderStaticChannelInfo->last_tcx_gain_e =
          pAacDecoderChannelInfo->data.usac.tcx_gain_e[k];

      if (last_lpd_mode == 0 && pAacDecoderChannelInfo->data.usac.fac_data[k]) {
        CFac_ApplyGains(
            pAacDecoderChannelInfo->data.usac.fac_data[k],
            pAacDecoderChannelInfo->granuleLength,
            pAacDecoderChannelInfo->data.usac.tcx_gain[k],
            pAacDecoderStaticChannelInfo->last_alfd_gains, mod[k]);

        pAacDecoderChannelInfo->data.usac.fac_data_e[k] +=
            pAacDecoderChannelInfo->data.usac.tcx_gain_e[k];
      }
    }

    last_lpd_mode = mod[k];
    k += (mod[k] == 0) ? 1 : (1 << (mod[k] - 1));
  }
}

AAC_DECODER_ERROR CLpd_RenderTimeSignal(
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo,
    CAacDecoderChannelInfo *pAacDecoderChannelInfo, PCM_DEC *pTimeData,
    INT lFrame, SamplingRateInfo *pSamplingRateInfo, UINT frameOk,
    const INT aacOutDataHeadroom, UINT flags, UINT strmFlags) {
  UCHAR *mod = pAacDecoderChannelInfo->data.usac.mod;
  AAC_DECODER_ERROR error = AAC_DEC_OK;
  int k, i_offset;
  int last_k;
  int nrSamples = 0;
  int facFB = 1;
  int nbDiv = NB_DIV;
  int lDiv = lFrame / nbDiv;
  int lFac = lDiv / 2;
  int nbSubfr =
      lFrame / (nbDiv * L_SUBFR);
  int nbSubfrSuperfr = nbDiv * nbSubfr;
  int synSfd = (nbSubfrSuperfr / 2) - BPF_SFD;
  int SynDelay = synSfd * L_SUBFR;
  int aacDelay = lFrame / 2;

  FIXP_DBL *synth_buf =
      pAacDecoderChannelInfo->pComStaticData->pWorkBufferCore1->synth_buf;
  FIXP_DBL *synth = synth_buf + PIT_MAX_MAX - BPF_DELAY;
  UCHAR last_lpd_mode, last_last_lpd_mode, last_lpc_lost, last_frame_lost;

  INT pitch[NB_SUBFR_SUPERFR + SYN_SFD];
  FIXP_DBL pit_gain[NB_SUBFR_SUPERFR + SYN_SFD];

  const int *lg_table;
  int lg_table_offset = 0;

  UINT samplingRate = pSamplingRateInfo->samplingRate;

  FDKmemclear(pitch, (NB_SUBFR_SUPERFR + SYN_SFD) * sizeof(INT));

  if (flags & AACDEC_FLUSH) {
    CLpd_Reset(pAacDecoderChannelInfo, pAacDecoderStaticChannelInfo,
               flags & AACDEC_FLUSH);
    frameOk = 0;
  }

  switch (lFrame) {
    case 1024:
      lg_table = &lg_table_ccfl[0][lg_table_offset];
      break;
    case 768:
      lg_table = &lg_table_ccfl[1][lg_table_offset];
      break;
    default:
      FDK_ASSERT(0);
      return AAC_DEC_UNKNOWN;
  }

  last_frame_lost = !CConcealment_GetLastFrameOk(
      &pAacDecoderStaticChannelInfo->concealmentInfo, 0);

  if ((pAacDecoderStaticChannelInfo->last_core_mode == FD_LONG) ||
      (pAacDecoderStaticChannelInfo->last_core_mode == FD_SHORT)) {
    pAacDecoderStaticChannelInfo->last_lpd_mode = 255;
  }

  if (!frameOk) {
    FIXP_DBL old_tcx_gain;
    FIXP_SGL old_stab;
    SCHAR old_tcx_gain_e;
    int nLostSf;

    last_lpd_mode = pAacDecoderStaticChannelInfo->last_lpd_mode;
    old_tcx_gain = pAacDecoderStaticChannelInfo->last_tcx_gain;
    old_tcx_gain_e = pAacDecoderStaticChannelInfo->last_tcx_gain_e;
    old_stab = pAacDecoderStaticChannelInfo->oldStability;
    nLostSf = pAacDecoderStaticChannelInfo->numLostLpdFrames;

    pAacDecoderChannelInfo->data.usac.lpd_mode_last = last_lpd_mode;

    switch (last_lpd_mode) {
      case 0:
        mod[0] = mod[1] = mod[2] = mod[3] = 0;
        break;
      case 3:
        mod[0] = mod[1] = mod[2] = mod[3] = 3;
        break;
      case 2:
        mod[0] = mod[1] = mod[2] = mod[3] = 2;
        break;
      case 1:
      default:
        mod[0] = mod[1] = mod[2] = mod[3] = 4;
        break;
    }

    CLpc_Conceal(pAacDecoderChannelInfo->data.usac.lsp_coeff,
                 pAacDecoderStaticChannelInfo->lpc4_lsf,
                 pAacDecoderStaticChannelInfo->lsf_adaptive_mean,
                 (last_lpd_mode == 255));

    if ((last_lpd_mode > 0) && (last_lpd_mode < 255)) {
      FDKmemcpy(pAacDecoderChannelInfo->data.usac.lp_coeff[0],
                pAacDecoderStaticChannelInfo->lp_coeff_old[0],
                M_LP_FILTER_ORDER * sizeof(FIXP_LPC));
      pAacDecoderChannelInfo->data.usac.lp_coeff_exp[0] =
          pAacDecoderStaticChannelInfo->lp_coeff_old_exp[0];
    }
    k = 0;
    while (k < nbDiv) {
      pAacDecoderChannelInfo->data.usac.tcx_gain[k] = old_tcx_gain;
      pAacDecoderChannelInfo->data.usac.tcx_gain_e[k] = old_tcx_gain_e;
      pAacDecoderChannelInfo->data.usac.aStability[k] = old_stab;

      k += ((mod[k] & 0x3) == 0) ? 1 : (1 << ((mod[k] & 0x3) - 1));

      nLostSf++;
    }
  } else {
    if ((pAacDecoderStaticChannelInfo->last_lpd_mode == 4) && (mod[0] > 0)) {
      FDKmemcpy(pAacDecoderChannelInfo->data.usac.lp_coeff[0],
                pAacDecoderStaticChannelInfo->lp_coeff_old[0],
                M_LP_FILTER_ORDER * sizeof(FIXP_LPC));
      pAacDecoderChannelInfo->data.usac.lp_coeff_exp[0] =
          pAacDecoderStaticChannelInfo->lp_coeff_old_exp[0];
    }
  }

  Acelp_PreProcessing(synth_buf, pAacDecoderStaticChannelInfo->old_synth, pitch,
                      pAacDecoderStaticChannelInfo->old_T_pf, pit_gain,
                      pAacDecoderStaticChannelInfo->old_gain_pf, samplingRate,
                      &i_offset, lFrame, synSfd, nbSubfrSuperfr);

  k = 0;
  last_k = -1;
  last_lpd_mode = pAacDecoderStaticChannelInfo->last_lpd_mode;
  last_last_lpd_mode = pAacDecoderStaticChannelInfo->last_last_lpd_mode;
  last_lpc_lost = pAacDecoderStaticChannelInfo->last_lpc_lost | last_frame_lost;

  while (k < nbDiv) {
    if (frameOk == 0) {
      pAacDecoderStaticChannelInfo->numLostLpdFrames++;
    } else {
      last_frame_lost |=
          (pAacDecoderStaticChannelInfo->numLostLpdFrames > 0) ? 1 : 0;
      pAacDecoderStaticChannelInfo->numLostLpdFrames = 0;
    }
    if (mod[k] == 0 || mod[k] == 4) {
      FIXP_DBL *acelp_out;

      if ((last_lpd_mode != 0) && (last_lpd_mode != 4))
      {
        FIXP_DBL *pFacData = NULL;

        if (frameOk && !last_frame_lost) {
          pFacData = pAacDecoderChannelInfo->data.usac.fac_data[k];
        }

        nrSamples += CLpd_FAC_Mdct2Acelp(
            &pAacDecoderStaticChannelInfo->IMdct, synth + nrSamples, pFacData,
            pAacDecoderChannelInfo->data.usac.fac_data_e[k],
            pAacDecoderChannelInfo->data.usac.lp_coeff[k],
            pAacDecoderChannelInfo->data.usac.lp_coeff_exp[k],
            lFrame - nrSamples,
            CLpd_FAC_getLength(
                (pAacDecoderStaticChannelInfo->last_core_mode != FD_SHORT) ||
                    (k > 0),
                lFac),
            (pAacDecoderStaticChannelInfo->last_core_mode != LPD) && (k == 0),
            0);

        FDKmemcpy(
            synth + nrSamples, pAacDecoderStaticChannelInfo->IMdct.overlap.time,
            pAacDecoderStaticChannelInfo->IMdct.ov_offset * sizeof(FIXP_DBL));
        {
          FIXP_LPC *lp_prev =
              pAacDecoderChannelInfo->data.usac
                  .lp_coeff[0];
          INT lp_prev_exp = pAacDecoderChannelInfo->data.usac.lp_coeff_exp[0];

          if (last_lpd_mode != 255) {
            last_k = k - (1 << (last_lpd_mode - 1));
            if (last_k < 0) {
              lp_prev = pAacDecoderStaticChannelInfo->lp_coeff_old[1];
              lp_prev_exp = pAacDecoderStaticChannelInfo->lp_coeff_old_exp[1];
            } else {
              lp_prev = pAacDecoderChannelInfo->data.usac.lp_coeff[last_k];
              lp_prev_exp =
                  pAacDecoderChannelInfo->data.usac.lp_coeff_exp[last_k];
            }
          }

          CLpd_AcelpPrepareInternalMem(
              synth + aacDelay + k * lDiv, last_lpd_mode,
              (last_last_lpd_mode == 4) ? 0 : last_last_lpd_mode,
              pAacDecoderChannelInfo->data.usac.lp_coeff[k],
              pAacDecoderChannelInfo->data.usac.lp_coeff_exp[k], lp_prev,
              lp_prev_exp, &pAacDecoderStaticChannelInfo->acelp, lFrame,
              (last_frame_lost && k < 2), mod[k]);
        }
      } else {
        if (k == 0 && pAacDecoderStaticChannelInfo->IMdct.ov_offset !=
                          lFrame / facFB / 2) {
          pAacDecoderStaticChannelInfo->IMdct.ov_offset = lFrame / facFB / 2;
        }
        nrSamples += imdct_drain(&pAacDecoderStaticChannelInfo->IMdct,
                                 synth + nrSamples, lFrame / facFB - nrSamples);
      }

      if (nrSamples >= lFrame / facFB) {
        acelp_out = pAacDecoderStaticChannelInfo->IMdct.overlap.time +
                    pAacDecoderStaticChannelInfo->IMdct.ov_offset;

        pAacDecoderStaticChannelInfo->IMdct.ov_offset += lDiv;
      } else {
        acelp_out = synth + nrSamples;

        nrSamples += lDiv;
      }

      if (mod[k] == 4) {
        pAacDecoderStaticChannelInfo->acelp.wsyn_rms = scaleValue(
            pAacDecoderChannelInfo->data.usac.tcx_gain[k],
            fixMin(0,
                   pAacDecoderChannelInfo->data.usac.tcx_gain_e[k] - SF_EXC));
        CLpd_TcxTDConceal(&pAacDecoderStaticChannelInfo->acelp,
                          &pAacDecoderStaticChannelInfo->last_tcx_pitch,
                          pAacDecoderChannelInfo->data.usac.lsp_coeff[k],
                          pAacDecoderChannelInfo->data.usac.lsp_coeff[k + 1],
                          pAacDecoderChannelInfo->data.usac.aStability[k],
                          pAacDecoderStaticChannelInfo->numLostLpdFrames,
                          acelp_out, lFrame,
                          pAacDecoderStaticChannelInfo->last_tcx_noise_factor);

      } else {
        FDK_ASSERT(pAacDecoderChannelInfo->data.usac.aStability[k] >=
                   (FIXP_SGL)0);
        CLpd_AcelpDecode(&pAacDecoderStaticChannelInfo->acelp, i_offset,
                         pAacDecoderChannelInfo->data.usac.lsp_coeff[k],
                         pAacDecoderChannelInfo->data.usac.lsp_coeff[k + 1],
                         pAacDecoderChannelInfo->data.usac.aStability[k],
                         &pAacDecoderChannelInfo->data.usac.acelp[k],
                         pAacDecoderStaticChannelInfo->numLostLpdFrames,
                         last_lpc_lost, k, acelp_out,
                         &pitch[(k * nbSubfr) + synSfd],
                         &pit_gain[(k * nbSubfr) + synSfd], lFrame);
      }

      if (mod[k] != 4) {
        if (last_lpd_mode != 0 &&
            pAacDecoderChannelInfo->data.usac
                .bpf_control_info) {
          int currentSf = synSfd + k * nbSubfr;

          if ((k > 0) || (pAacDecoderStaticChannelInfo->last_core_mode !=
                          FD_SHORT)) {
            pitch[currentSf - 2] = pitch[currentSf - 1] = pitch[currentSf];
            pit_gain[currentSf - 2] = pit_gain[currentSf - 1] =
                pit_gain[currentSf];
          } else {
            pitch[currentSf - 1] = pitch[currentSf];
            pit_gain[currentSf - 1] = pit_gain[currentSf];
          }
        }
      }
    } else {
      int lg = lg_table[mod[k]];
      int isFullBandLpd = 0;

      if ((last_lpd_mode == 0) || (last_lpd_mode == 4))
      {
        C_AALLOC_SCRATCH_START(fac_buf, FIXP_DBL, 1024 / 8);

        if (last_frame_lost == 1 ||
            pAacDecoderChannelInfo->data.usac.fac_data[k] == NULL) {
          FDKmemclear(fac_buf, 1024 / 8 * sizeof(FIXP_DBL));
          pAacDecoderChannelInfo->data.usac.fac_data[k] = fac_buf;
          pAacDecoderChannelInfo->data.usac.fac_data_e[k] = 0;
        }

        nrSamples += CLpd_FAC_Acelp2Mdct(
            &pAacDecoderStaticChannelInfo->IMdct, synth + nrSamples,
            SPEC_TCX(pAacDecoderChannelInfo->pSpectralCoefficient, k,
                     pAacDecoderChannelInfo->granuleLength, isFullBandLpd),
            pAacDecoderChannelInfo->specScale + k, 1,
            pAacDecoderChannelInfo->data.usac.fac_data[k],
            pAacDecoderChannelInfo->data.usac.fac_data_e[k],
            pAacDecoderChannelInfo->granuleLength,
            lFrame - nrSamples, lg,
            FDKgetWindowSlope(lDiv,
                              GetWindowShape(&pAacDecoderChannelInfo->icsInfo)),
            lDiv, pAacDecoderChannelInfo->data.usac.lp_coeff[k],
            pAacDecoderChannelInfo->data.usac.lp_coeff_exp[k],
            &pAacDecoderStaticChannelInfo->acelp,
            pAacDecoderChannelInfo->data.usac.tcx_gain[k],
            (last_frame_lost || !frameOk), 0
            ,
            last_lpd_mode, k,
            pAacDecoderChannelInfo
                ->currAliasingSymmetry
        );

        pitch[(k * nbSubfr) + synSfd + 1] = pitch[(k * nbSubfr) + synSfd] =
            pitch[(k * nbSubfr) + synSfd - 1];
        pit_gain[(k * nbSubfr) + synSfd + 1] =
            pit_gain[(k * nbSubfr) + synSfd] =
                pit_gain[(k * nbSubfr) + synSfd - 1];

        C_AALLOC_SCRATCH_END(fac_buf, FIXP_DBL, 1024 / 8);
      } else {
        int tl = lg;
        int fl = lDiv;
        int fr = lDiv;

        nrSamples += imlt_block(
            &pAacDecoderStaticChannelInfo->IMdct, synth + nrSamples,
            SPEC_TCX(pAacDecoderChannelInfo->pSpectralCoefficient, k,
                     pAacDecoderChannelInfo->granuleLength, isFullBandLpd),
            pAacDecoderChannelInfo->specScale + k, 1, lFrame - nrSamples, tl,
            FDKgetWindowSlope(fl,
                              GetWindowShape(&pAacDecoderChannelInfo->icsInfo)),
            fl,
            FDKgetWindowSlope(fr,
                              GetWindowShape(&pAacDecoderChannelInfo->icsInfo)),
            fr, pAacDecoderChannelInfo->data.usac.tcx_gain[k],
            pAacDecoderChannelInfo->currAliasingSymmetry
                ? MLT_FLAG_CURR_ALIAS_SYMMETRY
                : 0);
      }
    }
    last_last_lpd_mode = last_lpd_mode;
    last_lpd_mode = mod[k];
    last_lpc_lost = (frameOk == 0) ? 1 : 0;

    last_k = k;
    k += ((mod[k] & 0x3) == 0) ? 1 : (1 << (mod[k] - 1));
  }

  if (frameOk) {
    FDK_ASSERT(pAacDecoderChannelInfo->data.usac.aStability[last_k] >=
               (FIXP_SGL)0);
    pAacDecoderStaticChannelInfo->oldStability =
        pAacDecoderChannelInfo->data.usac.aStability[last_k];
    FDKmemcpy(pAacDecoderStaticChannelInfo->lsf_adaptive_mean,
              pAacDecoderChannelInfo->data.usac.lsf_adaptive_mean_cand,
              M_LP_FILTER_ORDER * sizeof(FIXP_LPC));
  }

  if (last_lpd_mode > 0) {
    FDKmemcpy(pAacDecoderStaticChannelInfo->lp_coeff_old[0],
              pAacDecoderChannelInfo->data.usac.lp_coeff[nbDiv],
              M_LP_FILTER_ORDER * sizeof(FIXP_LPC));
    pAacDecoderStaticChannelInfo->lp_coeff_old_exp[0] =
        pAacDecoderChannelInfo->data.usac.lp_coeff_exp[nbDiv];
    FDKmemcpy(pAacDecoderStaticChannelInfo->lp_coeff_old[1],
              pAacDecoderChannelInfo->data.usac.lp_coeff[last_k],
              M_LP_FILTER_ORDER * sizeof(FIXP_LPC));
    pAacDecoderStaticChannelInfo->lp_coeff_old_exp[1] =
        pAacDecoderChannelInfo->data.usac.lp_coeff_exp[last_k];
  }

  FDK_ASSERT(nrSamples == lFrame);

  if (!(pAacDecoderChannelInfo->data.usac.bpf_control_info)) {
    if (mod[0] != 0 && (pAacDecoderStaticChannelInfo->old_bpf_control_info)) {
      for (int i = 2; i < nbSubfrSuperfr; i++)
        pit_gain[synSfd + i] = (FIXP_DBL)0;
    } else {
      for (int i = 0; i < nbSubfrSuperfr; i++)
        pit_gain[synSfd + i] = (FIXP_DBL)0;
    }
  }

  for (int n = 0; n < synSfd; n++) {
    pAacDecoderStaticChannelInfo->old_T_pf[n] = pitch[nbSubfrSuperfr + n];
    pAacDecoderStaticChannelInfo->old_gain_pf[n] = pit_gain[nbSubfrSuperfr + n];
  }

  pAacDecoderStaticChannelInfo->old_bpf_control_info =
      pAacDecoderChannelInfo->data.usac.bpf_control_info;

  {
    INT lookahead = -BPF_DELAY;
    int copySamp = (mod[nbDiv - 1] == 0) ? (aacDelay) : (aacDelay - lFac);

    lookahead += imdct_copy_ov_and_nr(&pAacDecoderStaticChannelInfo->IMdct,
                                      synth + nrSamples, copySamp);

    FDK_ASSERT(lookahead == copySamp - BPF_DELAY);

    FIXP_DBL *p2_synth = synth + BPF_DELAY;

    for (int i = 0; i < nbSubfrSuperfr; i++) {
      int T = pitch[i];
      FIXP_DBL gain = pit_gain[i];

      if (gain > (FIXP_DBL)0) {
        gain = get_gain(&p2_synth[i * L_SUBFR], &p2_synth[(i * L_SUBFR) - T],
                        L_SUBFR);
        pit_gain[i] = gain;
      }
    }

    {
      bass_pf_1sf_delay(p2_synth, pitch, pit_gain, lFrame, lFrame / facFB,
                        mod[nbDiv - 1] ? (SynDelay - (lDiv / 2)) : SynDelay,
                        pTimeData, aacOutDataHeadroom,
                        pAacDecoderStaticChannelInfo->mem_bpf);
    }
  }

  Acelp_PostProcessing(synth_buf, pAacDecoderStaticChannelInfo->old_synth,
                       pitch, pAacDecoderStaticChannelInfo->old_T_pf, lFrame,
                       synSfd, nbSubfrSuperfr);

  { pAacDecoderStaticChannelInfo->last_core_mode = LPD; }
  pAacDecoderStaticChannelInfo->last_lpd_mode = last_lpd_mode;
  pAacDecoderStaticChannelInfo->last_last_lpd_mode = last_last_lpd_mode;
  pAacDecoderStaticChannelInfo->last_lpc_lost = last_lpc_lost;

  return error;
}
