/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "sbrenc_freq_sca.h"
#include "sbr_misc.h"
#include "genericStds.h"

static INT getStartFreq(INT fsCore, const INT start_freq);

static INT getStopFreq(INT fsCore, const INT stop_freq);

static INT numberOfBands(INT b_p_o, INT start, INT stop, FIXP_DBL warp_factor);
static void CalcBands(INT *diff, INT start, INT stop, INT num_bands);
static INT modifyBands(INT max_band, INT *diff, INT length);
static void cumSum(INT start_value, INT *diff, INT length, UCHAR *start_adress);

INT FDKsbrEnc_getSbrStartFreqRAW(INT startFreq, INT fsCore) {
  INT result;

  if (startFreq < 0 || startFreq > 15) {
    return -1;
  }
  result = getStartFreq(fsCore, startFreq);

  result =
      (result * (fsCore >> 5) + 1) >> 1;

  return (result);
}

INT FDKsbrEnc_getSbrStopFreqRAW(INT stopFreq, INT fsCore) {
  INT result;

  if (stopFreq < 0 || stopFreq > 13) return -1;

  result = getStopFreq(fsCore, stopFreq);
  result =
      (result * (fsCore >> 5) + 1) >> 1;

  return (result);
}

static INT getStartFreq(INT fsCore, const INT start_freq) {
  INT k0_min;

  switch (fsCore) {
    case 8000:
      k0_min = 24; /* (3000 * nQmfChannels / fsSBR ) + 0.5 */
      break;
    case 11025:
      k0_min = 17; /* (3000 * nQmfChannels / fsSBR ) + 0.5 */
      break;
    case 12000:
      k0_min = 16; /* (3000 * nQmfChannels / fsSBR ) + 0.5 */
      break;
    case 16000:
      k0_min = 16; /* (4000 * nQmfChannels / fsSBR ) + 0.5 */
      break;
    case 22050:
      k0_min = 12; /* (4000 * nQmfChannels / fsSBR ) + 0.5 */
      break;
    case 24000:
      k0_min = 11; /* (4000 * nQmfChannels / fsSBR ) + 0.5 */
      break;
    case 32000:
      k0_min = 10; /* (5000 * nQmfChannels / fsSBR ) + 0.5 */
      break;
    case 44100:
      k0_min = 7; /* (5000 * nQmfChannels / fsSBR ) + 0.5 */
      break;
    case 48000:
      k0_min = 7; /* (5000 * nQmfChannels / fsSBR ) + 0.5 */
      break;
    case 96000:
      k0_min = 3; /* (5000 * nQmfChannels / fsSBR ) + 0.5 */
      break;
    default:
      k0_min = 11; /* illegal fs */
  }

  switch (fsCore) {
    case 8000: {
      INT v_offset[] = {-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7};
      return (k0_min + v_offset[start_freq]);
    }
    case 11025: {
      INT v_offset[] = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13};
      return (k0_min + v_offset[start_freq]);
    }
    case 12000: {
      INT v_offset[] = {-5, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16};
      return (k0_min + v_offset[start_freq]);
    }
    case 16000: {
      INT v_offset[] = {-6, -4, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16};
      return (k0_min + v_offset[start_freq]);
    }
    case 22050:
    case 24000:
    case 32000: {
      INT v_offset[] = {-4, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16, 20};
      return (k0_min + v_offset[start_freq]);
    }
    case 44100:
    case 48000:
    case 96000: {
      INT v_offset[] = {-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16, 20, 24};
      return (k0_min + v_offset[start_freq]);
    }
    default: {
      INT v_offset[] = {0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 16, 20, 24, 28, 33};
      return (k0_min + v_offset[start_freq]);
    }
  }
}

static INT getStopFreq(INT fsCore, const INT stop_freq) {
  INT result, i;
  INT k1_min;
  INT v_dstop[13];

  INT *v_stop_freq = NULL;
  INT v_stop_freq_16[14] = {48, 49, 50, 51, 52, 54, 55,
                            56, 57, 59, 60, 61, 63, 64};
  INT v_stop_freq_22[14] = {35, 37, 38, 40, 42, 44, 46,
                            48, 51, 53, 56, 58, 61, 64};
  INT v_stop_freq_24[14] = {32, 34, 36, 38, 40, 42, 44,
                            46, 49, 52, 55, 58, 61, 64};
  INT v_stop_freq_32[14] = {32, 34, 36, 38, 40, 42, 44,
                            46, 49, 52, 55, 58, 61, 64};
  INT v_stop_freq_44[14] = {23, 25, 27, 29, 32, 34, 37,
                            40, 43, 47, 51, 55, 59, 64};
  INT v_stop_freq_48[14] = {21, 23, 25, 27, 30, 32, 35,
                            38, 42, 45, 49, 54, 59, 64};
  INT v_stop_freq_64[14] = {20, 22, 24, 26, 29, 31, 34,
                            37, 41, 45, 49, 54, 59, 64};
  INT v_stop_freq_88[14] = {15, 17, 19, 21, 23, 26, 29,
                            33, 37, 41, 46, 51, 57, 64};
  INT v_stop_freq_96[14] = {13, 15, 17, 19, 21, 24, 27,
                            31, 35, 39, 44, 50, 57, 64};
  INT v_stop_freq_192[14] = {7,  8,  10, 12, 14, 16, 19,
                             23, 27, 32, 38, 46, 54, 64};

  switch (fsCore) {
    case 8000:
      k1_min = 48;
      v_stop_freq = v_stop_freq_16;
      break;
    case 11025:
      k1_min = 35;
      v_stop_freq = v_stop_freq_22;
      break;
    case 12000:
      k1_min = 32;
      v_stop_freq = v_stop_freq_24;
      break;
    case 16000:
      k1_min = 32;
      v_stop_freq = v_stop_freq_32;
      break;
    case 22050:
      k1_min = 23;
      v_stop_freq = v_stop_freq_44;
      break;
    case 24000:
      k1_min = 21;
      v_stop_freq = v_stop_freq_48;
      break;
    case 32000:
      k1_min = 20;
      v_stop_freq = v_stop_freq_64;
      break;
    case 44100:
      k1_min = 15;
      v_stop_freq = v_stop_freq_88;
      break;
    case 48000:
      k1_min = 13;
      v_stop_freq = v_stop_freq_96;
      break;
    case 96000:
      k1_min = 7;
      v_stop_freq = v_stop_freq_192;
      break;
    default:
      k1_min = 21;
  }

  for (i = 0; i <= 12; i++) {
    v_dstop[i] = v_stop_freq[i + 1] - v_stop_freq[i];
  }

  FDKsbrEnc_Shellsort_int(v_dstop, 13);

  result = k1_min;
  for (i = 0; i < stop_freq; i++) {
    result = result + v_dstop[i];
  }

  return (result);
}

INT FDKsbrEnc_FindStartAndStopBand(const INT srSbr, const INT srCore,
                                   const INT noChannels, const INT startFreq,
                                   const INT stopFreq, INT *k0, INT *k2) {
  *k0 = getStartFreq(srCore, startFreq);

  if (srSbr * noChannels < *k0 * srCore) {
    return (
        1);
  }

  if (stopFreq < 14) {
    *k2 = getStopFreq(srCore, stopFreq);
  } else if (stopFreq == 14) {
    *k2 = 2 * *k0;
  } else {
    *k2 = 3 * *k0;
  }

  if (*k2 > noChannels) {
    *k2 = noChannels;
  }

  if ((srCore == 22050) && ((*k2 - *k0) > MAX_FREQ_COEFFS_FS44100))
    return (1);

  if ((srCore >= 24000) && ((*k2 - *k0) > MAX_FREQ_COEFFS_FS48000))
    return (1);

  if ((*k2 - *k0) > MAX_FREQ_COEFFS)
    return (1);

  if ((*k2 - *k0) < 0) return (1);

  return (0);
}

INT FDKsbrEnc_UpdateFreqScale(UCHAR *v_k_master, INT *h_num_bands, const INT k0,
                              const INT k2, const INT freqScale,
                              const INT alterScale)

{
  INT b_p_o = 0;
  FIXP_DBL warp = FL2FXCONST_DBL(0.0f);
  INT dk = 0;

  INT k1 = 0, i;
  INT num_bands0;
  INT num_bands1;
  INT diff_tot[MAX_OCTAVE + MAX_SECOND_REGION];
  INT *diff0 = diff_tot;
  INT *diff1 = diff_tot + MAX_OCTAVE;
  INT k2_achived;
  INT k2_diff;
  INT incr = 0;

  if (freqScale == 1) b_p_o = 12;
  if (freqScale == 2) b_p_o = 10;
  if (freqScale == 3) b_p_o = 8;

  if (freqScale > 0)
  {
    if (alterScale == 0)
      warp = FL2FXCONST_DBL(0.5f);
    else
      warp = FL2FXCONST_DBL(1.0f / 2.6f);

    if (4 * k2 >= 9 * k0)
    {
      k1 = 2 * k0;

      num_bands0 = numberOfBands(b_p_o, k0, k1, FL2FXCONST_DBL(0.5f));
      num_bands1 = numberOfBands(b_p_o, k1, k2, warp);

      CalcBands(diff0, k0, k1, num_bands0);
      FDKsbrEnc_Shellsort_int(diff0, num_bands0);

      if (diff0[0] == 0)
      {
        return (1);
      }

      cumSum(k0, diff0, num_bands0, v_k_master);

      CalcBands(diff1, k1, k2, num_bands1);
      FDKsbrEnc_Shellsort_int(diff1, num_bands1);
      if (diff0[num_bands0 - 1] > diff1[0])
      {
        if (modifyBands(diff0[num_bands0 - 1], diff1, num_bands1)) return (1);
      }

      cumSum(k1, diff1, num_bands1, &v_k_master[num_bands0]);
      *h_num_bands = num_bands0 + num_bands1;

    } else
    {
      k1 = k2;

      num_bands0 = numberOfBands(b_p_o, k0, k1, FL2FXCONST_DBL(0.5f));
      CalcBands(diff0, k0, k1, num_bands0);
      FDKsbrEnc_Shellsort_int(diff0, num_bands0);

      if (diff0[0] == 0)
      {
        return (1);
      }

      cumSum(k0, diff0, num_bands0, v_k_master);
      *h_num_bands = num_bands0;
    }
  } else
  {
    if (alterScale == 0) {
      dk = 1;
      num_bands0 = 2 * ((k2 - k0) / 2);
    } else {
      dk = 2;
      num_bands0 =
          2 * (((k2 - k0) / dk + 1) / 2);
    }

    k2_achived = k0 + num_bands0 * dk;
    k2_diff = k2 - k2_achived;

    for (i = 0; i < num_bands0; i++) diff_tot[i] = dk;

    if (k2_diff < 0) {
      incr = 1;
      i = 0;
    }

    if (k2_diff > 0) {
      incr = -1;
      i = num_bands0 - 1;
    }

    while (k2_diff != 0) {
      diff_tot[i] = diff_tot[i] - incr;
      i = i + incr;
      k2_diff = k2_diff + incr;
    }

    cumSum(k0, diff_tot, num_bands0, v_k_master);
    *h_num_bands = num_bands0;
  }

  if (*h_num_bands < 1) return (1);

  return (0);
}

static INT numberOfBands(INT b_p_o, INT start, INT stop, FIXP_DBL warp_factor) {
  INT result = 0;
  result = ((b_p_o * fMult((CalcLdInt(stop) - CalcLdInt(start)), warp_factor) +
             (FL2FX_DBL(0.5f) >> LD_DATA_SHIFT)) >>
            ((DFRACT_BITS - 1) - LD_DATA_SHIFT))
           << 1;

  return (result);
}

static void CalcBands(INT *diff, INT start, INT stop, INT num_bands) {
  INT i, qb, qe, qtmp;
  INT previous;
  INT current;
  FIXP_DBL base, exp, tmp;

  previous = start;
  for (i = 1; i <= num_bands; i++) {
    base = fDivNorm((FIXP_DBL)stop, (FIXP_DBL)start, &qb);
    exp = fDivNorm((FIXP_DBL)i, (FIXP_DBL)num_bands, &qe);
    tmp = fPow(base, qb, exp, qe, &qtmp);
    tmp = fMult(tmp, (FIXP_DBL)(start << 24));
    current = (INT)scaleValue(tmp, qtmp - 23);
    current = (current + 1) >> 1;
    diff[i - 1] = current - previous;
    previous = current;
  }
}

static void cumSum(INT start_value, INT *diff, INT length,
                   UCHAR *start_adress) {
  INT i;
  start_adress[0] = start_value;
  for (i = 1; i <= length; i++)
    start_adress[i] = start_adress[i - 1] + diff[i - 1];
}

static INT modifyBands(INT max_band_previous, INT *diff, INT length) {
  INT change = max_band_previous - diff[0];

  if (change > (diff[length - 1] - diff[0]) / 2)
    change = (diff[length - 1] - diff[0]) / 2;

  diff[0] += change;
  diff[length - 1] -= change;
  FDKsbrEnc_Shellsort_int(diff, length);

  return (0);
}

INT FDKsbrEnc_UpdateHiRes(UCHAR *h_hires, INT *num_hires, UCHAR *v_k_master,
                          INT num_master, INT *xover_band) {
  INT i;
  INT max1, max2;

  if ((v_k_master[*xover_band] >
       32) ||
      (*xover_band > num_master)) {

    max1 = 0;
    max2 = num_master;
    while ((v_k_master[max1 + 1] < 32) &&
           ((max1 + 1) < max2)) {
      max1++;
    }

    *xover_band = max1;
  }

  *num_hires = num_master - *xover_band;
  for (i = *xover_band; i <= num_master; i++) {
    h_hires[i - *xover_band] = v_k_master[i];
  }

  return (0);
}

void FDKsbrEnc_UpdateLoRes(UCHAR *h_lores, INT *num_lores, UCHAR *h_hires,
                           INT num_hires) {
  INT i;

  if (num_hires % 2 == 0)
  {
    *num_lores = num_hires / 2;
    for (i = 0; i <= *num_lores; i++) h_lores[i] = h_hires[i * 2];

  } else
  {
    *num_lores = (num_hires + 1) / 2;

    h_lores[0] = h_hires[0];
    for (i = 1; i <= *num_lores; i++) {
      h_lores[i] = h_hires[i * 2 - 1];
    }
  }
}
