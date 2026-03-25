/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "sbrdec_freq_sca.h"
#include "transcendent.h"
#include "sbr_rom.h"
#include "env_extr.h"
#include "genericStds.h"

#define MAX_OCTAVE 29
#define MAX_SECOND_REGION 50

static int numberOfBands(FIXP_SGL bpo_div16, int start, int stop, int warpFlag);
static void CalcBands(UCHAR *diff, UCHAR start, UCHAR stop, UCHAR num_bands);
static SBR_ERROR modifyBands(UCHAR max_band, UCHAR *diff, UCHAR length);
static void cumSum(UCHAR start_value, UCHAR *diff, UCHAR length,
                   UCHAR *start_adress);

static UCHAR getStartBand(
    UINT fs,
    UCHAR startFreq,
    UINT headerDataFlags)
{
  INT band;
  UINT fsMapped = fs;
  SBR_RATE rate = DUAL;

  if (headerDataFlags & (SBRDEC_SYNTAX_USAC | SBRDEC_SYNTAX_RSVD50)) {
    if (headerDataFlags & SBRDEC_QUAD_RATE) {
      rate = QUAD;
    }
    fsMapped = sbrdec_mapToStdSampleRate(fs, 1);
  }

  FDK_ASSERT(2 * (rate + 1) <= (4));

  switch (fsMapped) {
    case 192000:
      band = FDK_sbrDecoder_sbr_start_freq_192[startFreq];
      break;
    case 176400:
      band = FDK_sbrDecoder_sbr_start_freq_176[startFreq];
      break;
    case 128000:
      band = FDK_sbrDecoder_sbr_start_freq_128[startFreq];
      break;
    case 96000:
    case 88200:
      band = FDK_sbrDecoder_sbr_start_freq_88[rate][startFreq];
      break;
    case 64000:
      band = FDK_sbrDecoder_sbr_start_freq_64[rate][startFreq];
      break;
    case 48000:
      band = FDK_sbrDecoder_sbr_start_freq_48[rate][startFreq];
      break;
    case 44100:
      band = FDK_sbrDecoder_sbr_start_freq_44[rate][startFreq];
      break;
    case 40000:
      band = FDK_sbrDecoder_sbr_start_freq_40[rate][startFreq];
      break;
    case 32000:
      band = FDK_sbrDecoder_sbr_start_freq_32[rate][startFreq];
      break;
    case 24000:
      band = FDK_sbrDecoder_sbr_start_freq_24[rate][startFreq];
      break;
    case 22050:
      band = FDK_sbrDecoder_sbr_start_freq_22[rate][startFreq];
      break;
    case 16000:
      band = FDK_sbrDecoder_sbr_start_freq_16[rate][startFreq];
      break;
    default:
      band = 255;
  }

  return band;
}

static UCHAR getStopBand(
    UINT fs,
    UCHAR stopFreq,
    UINT headerDataFlags,
    UCHAR k0)
{
  UCHAR k2;

  if (stopFreq < 14) {
    INT stopMin;
    INT num = 2 * (64);
    UCHAR diff_tot[MAX_OCTAVE + MAX_SECOND_REGION];
    UCHAR *diff0 = diff_tot;
    UCHAR *diff1 = diff_tot + MAX_OCTAVE;

    if (headerDataFlags & SBRDEC_QUAD_RATE) {
      num >>= 1;
    }

    if (fs < 32000) {
      stopMin = (((2 * 6000 * num) / fs) + 1) >> 1;
    } else {
      if (fs < 64000) {
        stopMin = (((2 * 8000 * num) / fs) + 1) >> 1;
      } else {
        stopMin = (((2 * 10000 * num) / fs) + 1) >> 1;
      }
    }

    stopMin = fMin(stopMin, 64);

    CalcBands(diff0, stopMin, 64, 13);
    shellsort(diff0, 13);
    cumSum(stopMin, diff0, 13, diff1);
    k2 = diff1[stopFreq];
  } else if (stopFreq == 14)
    k2 = 2 * k0;
  else
    k2 = 3 * k0;

  if (k2 > (64)) k2 = (64);

  {
    UCHAR max_freq_coeffs = (headerDataFlags & SBRDEC_QUAD_RATE)
                                ? MAX_FREQ_COEFFS_QUAD_RATE
                                : MAX_FREQ_COEFFS;
    if (((k2 - k0) > max_freq_coeffs) || (k2 <= k0)) {
      return 255;
    }
  }

  if (headerDataFlags & SBRDEC_QUAD_RATE) {
    return k2;
  }
  if (headerDataFlags & (SBRDEC_SYNTAX_USAC | SBRDEC_SYNTAX_RSVD50)) {
    if ((fs >= 42000) && ((k2 - k0) > MAX_FREQ_COEFFS_FS44100)) {
      return 255;
    }
    if ((fs >= 46009) && ((k2 - k0) > MAX_FREQ_COEFFS_FS48000)) {
      return 255;
    }
  } else {
    if ((fs == 44100) && ((k2 - k0) > MAX_FREQ_COEFFS_FS44100)) {
      return 255;
    }
    if ((fs >= 48000) && ((k2 - k0) > MAX_FREQ_COEFFS_FS48000)) {
      return 255;
    }
  }

  return k2;
}

SBR_ERROR
sbrdecUpdateFreqScale(
    UCHAR *v_k_master,
    UCHAR *numMaster,
    UINT fs,
    HANDLE_SBR_HEADER_DATA hHeaderData,
    UINT flags) {
  FIXP_SGL bpo_div16;
  INT dk = 0;

  UCHAR k0, k2, i;
  UCHAR num_bands0 = 0;
  UCHAR num_bands1 = 0;
  UCHAR diff_tot[MAX_OCTAVE + MAX_SECOND_REGION];
  UCHAR *diff0 = diff_tot;
  UCHAR *diff1 = diff_tot + MAX_OCTAVE;
  INT k2_achived;
  INT k2_diff;
  INT incr = 0;

  if (flags & SBRDEC_QUAD_RATE) {
    fs >>= 1;
  }

  k0 = getStartBand(fs, hHeaderData->bs_data.startFreq, flags);
  if (k0 == 255) {
    return SBRDEC_UNSUPPORTED_CONFIG;
  }

  k2 = getStopBand(fs, hHeaderData->bs_data.stopFreq, flags, k0);
  if (k2 == 255) {
    return SBRDEC_UNSUPPORTED_CONFIG;
  }

  if (hHeaderData->bs_data.freqScale > 0) {
    INT k1;

    if (hHeaderData->bs_data.freqScale == 1) {
      bpo_div16 = FL2FXCONST_SGL(12.0f / 16.0f);
    } else if (hHeaderData->bs_data.freqScale == 2) {
      bpo_div16 = FL2FXCONST_SGL(10.0f / 16.0f);
    } else {
      bpo_div16 = FL2FXCONST_SGL(8.0f / 16.0f);
    }

    if (flags & SBRDEC_QUAD_RATE) {
      if ((SHORT)k0 < (SHORT)(bpo_div16 >> ((FRACT_BITS - 1) - 4))) {
        bpo_div16 = (FIXP_SGL)(k0 & (UCHAR)0xfe)
                    << ((FRACT_BITS - 1) - 4);
      }
    }

    if (1000 * k2 > 2245 * k0) {
      k1 = 2 * k0;

      num_bands0 = numberOfBands(bpo_div16, k0, k1, 0);
      num_bands1 =
          numberOfBands(bpo_div16, k1, k2, hHeaderData->bs_data.alterScale);
      if (num_bands0 < 1) {
        return SBRDEC_UNSUPPORTED_CONFIG;
      }
      if (num_bands1 < 1) {
        return SBRDEC_UNSUPPORTED_CONFIG;
      }

      CalcBands(diff0, k0, k1, num_bands0);
      shellsort(diff0, num_bands0);
      if (diff0[0] == 0) {
        return SBRDEC_UNSUPPORTED_CONFIG;
      }

      cumSum(k0, diff0, num_bands0, v_k_master);

      CalcBands(diff1, k1, k2, num_bands1);
      shellsort(diff1, num_bands1);
      if (diff0[num_bands0 - 1] > diff1[0]) {
        SBR_ERROR err;

        err = modifyBands(diff0[num_bands0 - 1], diff1, num_bands1);
        if (err) return SBRDEC_UNSUPPORTED_CONFIG;
      }

      cumSum(k1, diff1, num_bands1, &v_k_master[num_bands0]);
      *numMaster = num_bands0 + num_bands1;

    } else {
      k1 = k2;

      num_bands0 = numberOfBands(bpo_div16, k0, k1, 0);
      if (num_bands0 < 1) {
        return SBRDEC_UNSUPPORTED_CONFIG;
      }
      CalcBands(diff0, k0, k1, num_bands0);
      shellsort(diff0, num_bands0);
      if (diff0[0] == 0) {
        return SBRDEC_UNSUPPORTED_CONFIG;
      }

      cumSum(k0, diff0, num_bands0, v_k_master);
      *numMaster = num_bands0;
    }
  } else {
    if (hHeaderData->bs_data.alterScale == 0) {
      dk = 1;
      num_bands0 = (k2 - k0) & 254;
    } else {
      dk = 2;
      num_bands0 = (((k2 - k0) >> 1) + 1) & 254;
    }

    if (num_bands0 < 1) {
      return SBRDEC_UNSUPPORTED_CONFIG;
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
    *numMaster = num_bands0;
  }

  if (*numMaster < 1) {
    return SBRDEC_UNSUPPORTED_CONFIG;
  }

  if (flags & SBRDEC_QUAD_RATE) {
    int k;
    for (k = 1; k < *numMaster; k++) {
      if (!(v_k_master[k] - v_k_master[k - 1] <= k0 - 2)) {
        return SBRDEC_UNSUPPORTED_CONFIG;
      }
    }
  }

  return SBRDEC_OK;
}

static FIXP_SGL calcFactorPerBand(int k_start, int k_stop, int num_bands) {
  FIXP_DBL bandfactor = FL2FXCONST_DBL(0.25f);
  FIXP_DBL step = FL2FXCONST_DBL(0.125f);

  int direction = 1;

  FIXP_DBL start = k_start << (DFRACT_BITS - 8);
  FIXP_DBL stop = k_stop << (DFRACT_BITS - 8);

  FIXP_DBL temp;

  int j, i = 0;

  while (step > FL2FXCONST_DBL(0.0f)) {
    i++;
    temp = stop;

    for (j = 0; j < num_bands; j++)
      temp = fMultDiv2(temp, bandfactor) << 2;

    if (temp < start) {
      if (direction == 0)
        step = (FIXP_DBL)((LONG)step >> 1);
      direction = 1;
      bandfactor = bandfactor + step;
    } else {
      if (direction == 1) step = (FIXP_DBL)((LONG)step >> 1);
      direction = 0;
      bandfactor = bandfactor - step;
    }

    if (i > 100) {
      step = FL2FXCONST_DBL(0.0f);
    }
  }
  return (bandfactor >= FL2FXCONST_DBL(0.5)) ? (FIXP_SGL)MAXVAL_SGL
                                             : FX_DBL2FX_SGL(bandfactor << 1);
}

static int numberOfBands(
    FIXP_SGL bpo_div16,
    int start,
    int stop,
    int warpFlag)
{
  FIXP_SGL num_bands_div128;
  int num_bands;

  num_bands_div128 =
      FX_DBL2FX_SGL(fMult(FDK_getNumOctavesDiv8(start, stop), bpo_div16));

  if (warpFlag) {
    num_bands_div128 = FX_DBL2FX_SGL(
        fMult(num_bands_div128, FL2FXCONST_SGL(25200.0 / 32768.0)));
  }

  num_bands_div128 = num_bands_div128 + FL2FXCONST_SGL(1.0f / 128.0f);
  num_bands = 2 * ((LONG)num_bands_div128 >> (FRACT_BITS - 7));

  return (num_bands);
}

static void CalcBands(UCHAR *diff,
                      UCHAR start,
                      UCHAR stop,
                      UCHAR num_bands)
{
  int i;
  int previous;
  int current;
  FIXP_SGL exact, temp;
  FIXP_SGL bandfactor = calcFactorPerBand(start, stop, num_bands);

  previous = stop;
  exact = (FIXP_SGL)(
      stop << (FRACT_BITS - 8));

  for (i = num_bands - 1; i >= 0; i--) {
    exact = FX_DBL2FX_SGL(fMult(exact, bandfactor));

    temp = exact + FL2FXCONST_SGL(128.0 / 32768.0);

    current = (LONG)temp >> (FRACT_BITS - 8);

    diff[i] = previous - current;
    previous = current;
  }
}

static void cumSum(UCHAR start_value, UCHAR *diff, UCHAR length,
                   UCHAR *start_adress) {
  int i;
  start_adress[0] = start_value;
  for (i = 1; i <= length; i++)
    start_adress[i] = start_adress[i - 1] + diff[i - 1];
}

static SBR_ERROR modifyBands(UCHAR max_band_previous, UCHAR *diff,
                             UCHAR length) {
  int change = max_band_previous - diff[0];

  if (change > (diff[length - 1] - diff[0]) >> 1)
    change = (diff[length - 1] - diff[0]) >> 1;

  diff[0] += change;
  diff[length - 1] -= change;
  shellsort(diff, length);

  return SBRDEC_OK;
}

static void sbrdecUpdateHiRes(UCHAR *h_hires, UCHAR *num_hires,
                              UCHAR *v_k_master, UCHAR num_bands,
                              UCHAR xover_band) {
  UCHAR i;

  *num_hires = num_bands - xover_band;

  for (i = xover_band; i <= num_bands; i++) {
    h_hires[i - xover_band] = v_k_master[i];
  }
}

static void sbrdecUpdateLoRes(UCHAR *h_lores, UCHAR *num_lores, UCHAR *h_hires,
                              UCHAR num_hires) {
  UCHAR i;

  if ((num_hires & 1) == 0) {
    *num_lores = num_hires >> 1;
    for (i = 0; i <= *num_lores; i++) h_lores[i] = h_hires[i * 2];
  } else {
    *num_lores = (num_hires + 1) >> 1;
    h_lores[0] = h_hires[0];
    for (i = 1; i <= *num_lores; i++) {
      h_lores[i] = h_hires[i * 2 - 1];
    }
  }
}

void sbrdecDownSampleLoRes(UCHAR *v_result, UCHAR num_result,
                           UCHAR *freqBandTableRef, UCHAR num_Ref) {
  int step;
  int i, j;
  int org_length, result_length;
  int v_index[MAX_FREQ_COEFFS >> 1];

  org_length = num_Ref;
  result_length = num_result;

  v_index[0] = 0;
  i = 0;
  while (org_length > 0) {
    i++;
    step = org_length / result_length;
    org_length = org_length - step;
    result_length--;
    v_index[i] = v_index[i - 1] + step;
  }

  for (j = 0; j <= i; j++) {
    v_result[j] = freqBandTableRef[v_index[j]];
  }
}

void shellsort(UCHAR *in, UCHAR n) {
  int i, j, v, w;
  int inc = 1;

  do
    inc = 3 * inc + 1;
  while (inc <= n);

  do {
    inc = inc / 3;
    for (i = inc; i < n; i++) {
      v = in[i];
      j = i;
      while ((w = in[j - inc]) > v) {
        in[j] = w;
        j -= inc;
        if (j < inc) break;
      }
      in[j] = v;
    }
  } while (inc > 1);
}

SBR_ERROR
resetFreqBandTables(HANDLE_SBR_HEADER_DATA hHeaderData, const UINT flags) {
  SBR_ERROR err = SBRDEC_OK;
  int k2, kx, lsb, usb;
  int intTemp;
  UCHAR nBandsLo, nBandsHi;
  HANDLE_FREQ_BAND_DATA hFreq = &hHeaderData->freqBandData;

  err = sbrdecUpdateFreqScale(hFreq->v_k_master, &hFreq->numMaster,
                              hHeaderData->sbrProcSmplRate, hHeaderData, flags);

  if (err || (hHeaderData->bs_info.xover_band > hFreq->numMaster)) {
    return SBRDEC_UNSUPPORTED_CONFIG;
  }

  sbrdecUpdateHiRes(hFreq->freqBandTable[1], &nBandsHi, hFreq->v_k_master,
                    hFreq->numMaster, hHeaderData->bs_info.xover_band);
  sbrdecUpdateLoRes(hFreq->freqBandTable[0], &nBandsLo, hFreq->freqBandTable[1],
                    nBandsHi);

  if (!(nBandsLo > 0) ||
      (nBandsLo > (((hHeaderData->numberOfAnalysisBands == 16)
                        ? MAX_FREQ_COEFFS_QUAD_RATE
                        : MAX_FREQ_COEFFS_DUAL_RATE) >>
                   1))) {
    return SBRDEC_UNSUPPORTED_CONFIG;
  }

  hFreq->nSfb[0] = nBandsLo;
  hFreq->nSfb[1] = nBandsHi;

  lsb = hFreq->freqBandTable[0][0];
  usb = hFreq->freqBandTable[0][nBandsLo];

  if ((lsb > ((flags & SBRDEC_QUAD_RATE) ? 16 : (32))) || (lsb >= usb)) {
    return SBRDEC_UNSUPPORTED_CONFIG;
  }

  k2 = hFreq->freqBandTable[1][nBandsHi];
  kx = hFreq->freqBandTable[1][0];

  if (hHeaderData->bs_data.noise_bands == 0) {
    hFreq->nNfb = 1;
  } else
  {
    intTemp = (LONG)FDK_getNumOctavesDiv8(kx, k2) >> 2;
    intTemp = intTemp * hHeaderData->bs_data.noise_bands;
    intTemp = intTemp + (LONG)FL2FXCONST_SGL(0.5f / 32.0f);
    intTemp = intTemp >> (FRACT_BITS - 1 - 5);

    if (intTemp == 0) intTemp = 1;

    if (intTemp > MAX_NOISE_COEFFS) {
      return SBRDEC_UNSUPPORTED_CONFIG;
    }

    hFreq->nNfb = intTemp;
  }

  hFreq->nInvfBands = hFreq->nNfb;

  sbrdecDownSampleLoRes(hFreq->freqBandTableNoise, hFreq->nNfb,
                        hFreq->freqBandTable[0], nBandsLo);

  hFreq->ov_highSubband = hFreq->highSubband;
  hFreq->lowSubband = lsb;
  hFreq->highSubband = usb;

  return SBRDEC_OK;
}
