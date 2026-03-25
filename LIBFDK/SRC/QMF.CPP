/******************* Library for basic calculation routines ********************

   Author(s):   Markus Lohwasser, Josef Hoepfl, Manuel Jander

   Description: QMF filterbank

*******************************************************************************/

#include "qmf.h"

#include "FDK_trigFcts.h"
#include "fixpoint_math.h"
#include "dct.h"

#define QSSCALE (0)
#define FX_DBL2FX_QSS(x) (x)
#define FX_QSS2FX_DBL(x) (x)

static void qmfForwardModulationLP_even(
    HANDLE_QMF_FILTER_BANK anaQmf,
    FIXP_DBL *timeIn,
    FIXP_DBL *rSubband)
{
  int i;
  int L = anaQmf->no_channels;
  int M = L >> 1;
  int scale = 0;
  FIXP_DBL accu;

  const FIXP_DBL *timeInTmp1 = (FIXP_DBL *)&timeIn[3 * M];
  const FIXP_DBL *timeInTmp2 = timeInTmp1;
  FIXP_DBL *rSubbandTmp = rSubband;

  rSubband[0] = timeIn[3 * M] >> 1;

  for (i = M - 1; i != 0; i--) {
    accu = ((*--timeInTmp1) >> 1) + ((*++timeInTmp2) >> 1);
    *++rSubbandTmp = accu;
  }

  timeInTmp1 = &timeIn[2 * M];
  timeInTmp2 = &timeIn[0];
  rSubbandTmp = &rSubband[M];

  for (i = L - M; i != 0; i--) {
    accu = ((*timeInTmp1--) >> 1) - ((*timeInTmp2++) >> 1);
    *rSubbandTmp++ = accu;
  }

  dct_III(rSubband, timeIn, L, &scale);
}

#if !defined(FUNCTION_qmfForwardModulationLP_odd)
static void qmfForwardModulationLP_odd(
    HANDLE_QMF_FILTER_BANK anaQmf,
    const FIXP_DBL *timeIn,
    FIXP_DBL *rSubband)
{
  int i;
  int L = anaQmf->no_channels;
  int M = L >> 1;
  int shift = (anaQmf->no_channels >> 6) + 1;

  for (i = 0; i < M; i++) {
    rSubband[M + i] = (timeIn[L - 1 - i] >> 1) - (timeIn[i] >> shift);
    rSubband[M - 1 - i] =
        (timeIn[L + i] >> 1) + (timeIn[2 * L - 1 - i] >> shift);
  }

  dct_IV(rSubband, L, &shift);
}
#endif

#if !defined(FUNCTION_qmfForwardModulationHQ)
static void qmfForwardModulationHQ(
    HANDLE_QMF_FILTER_BANK anaQmf,
    const FIXP_DBL *RESTRICT timeIn,
    FIXP_DBL *RESTRICT rSubband,
    FIXP_DBL *RESTRICT iSubband
) {
  int i;
  int L = anaQmf->no_channels;
  int L2 = L << 1;
  int shift = 0;

  if ((L == 64) && !(anaQmf->flags & (QMF_FLAG_CLDFB | QMF_FLAG_MPSLDFB))) {
    FIXP_DBL x, y;

    x = timeIn[1] >> 1;
    y = timeIn[0];
    rSubband[0] = x + (y >> 1);
    iSubband[0] = x - (y >> 1);

    for (i = 1; i < L; i++) {
      x = timeIn[i + 1] >> 1;
      y = timeIn[L2 - i];
      rSubband[i] = x - (y >> 1);
      iSubband[i] = x + (y >> 1);
    }
  } else {
    for (i = 0; i < L; i += 2) {
      FIXP_DBL x0, x1, y0, y1;

      x0 = timeIn[i + 0] >> 1;
      x1 = timeIn[i + 1] >> 1;
      y0 = timeIn[L2 - 1 - i];
      y1 = timeIn[L2 - 2 - i];

      rSubband[i + 0] = x0 - (y0 >> 1);
      rSubband[i + 1] = x1 - (y1 >> 1);
      iSubband[i + 0] = x0 + (y0 >> 1);
      iSubband[i + 1] = x1 + (y1 >> 1);
    }
  }

  dct_IV(rSubband, L, &shift);
  dst_IV(iSubband, L, &shift);

  if ((L != 64) || (anaQmf->flags & (QMF_FLAG_CLDFB | QMF_FLAG_MPSLDFB))) {
    if (anaQmf->flags & QMF_FLAG_MPSLDFB_OPTIMIZE_MODULATION) {
      FIXP_DBL iBand;
      for (i = 0; i < fMin(anaQmf->lsb, L); i += 2) {
        iBand = rSubband[i];
        rSubband[i] = -iSubband[i];
        iSubband[i] = iBand;

        iBand = -rSubband[i + 1];
        rSubband[i + 1] = iSubband[i + 1];
        iSubband[i + 1] = iBand;
      }
    } else {
      const FIXP_QTW *sbr_t_cos;
      const FIXP_QTW *sbr_t_sin;
      const int len = L;
      sbr_t_cos = anaQmf->t_cos;
      sbr_t_sin = anaQmf->t_sin;

      for (i = 0; i < len; i++) {
        cplxMult(&iSubband[i], &rSubband[i], iSubband[i], rSubband[i],
                 sbr_t_cos[i], sbr_t_sin[i]);
      }
    }
  }
}
#endif

inline static void qmfInverseModulationLP_even(
    HANDLE_QMF_FILTER_BANK synQmf,
    const FIXP_DBL *qmfReal,
    const int scaleFactorLowBand,
    const int scaleFactorHighBand,
    FIXP_DBL *pTimeOut
) {
  int i;
  int L = synQmf->no_channels;
  int M = L >> 1;
  int scale = 0;
  FIXP_DBL tmp;
  FIXP_DBL *RESTRICT tReal = pTimeOut;
  FIXP_DBL *RESTRICT tImag = pTimeOut + L;

  scaleValuesSaturate(&tReal[0], &qmfReal[0], synQmf->lsb, scaleFactorLowBand);
  scaleValuesSaturate(&tReal[0 + synQmf->lsb], &qmfReal[0 + synQmf->lsb],
                      synQmf->usb - synQmf->lsb, scaleFactorHighBand);
  FDKmemclear(&tReal[0 + synQmf->usb], (L - synQmf->usb) * sizeof(FIXP_DBL));

  dct_II(tReal, tImag, L, &scale);

  tImag[0] = tReal[M];
  tImag[M] = (FIXP_DBL)0;
  tmp = tReal[0];
  tReal[0] = tReal[M];
  tReal[M] = tmp;

  for (i = 1; i < M / 2; i++) {
    tmp = tReal[L - i];
    tImag[M - i] = tmp;
    tImag[i + M] = -tmp;

    tmp = tReal[M + i];
    tImag[i] = tmp;
    tImag[L - i] = -tmp;

    tReal[M + i] = tReal[i];
    tReal[L - i] = tReal[M - i];
    tmp = tReal[i];
    tReal[i] = tReal[M - i];
    tReal[M - i] = tmp;
  }
  tmp = tReal[M + M / 2];
  tImag[M / 2] = tmp;
  tImag[M / 2 + M] = -tmp;

  tReal[M + M / 2] = tReal[M / 2];
}

inline static void qmfInverseModulationLP_odd(
    HANDLE_QMF_FILTER_BANK synQmf,
    const FIXP_DBL *qmfReal,
    const int scaleFactorLowBand,
    const int scaleFactorHighBand,
    FIXP_DBL *pTimeOut
) {
  int i;
  int L = synQmf->no_channels;
  int M = L >> 1;
  int shift = 0;

  scaleValuesSaturate(pTimeOut + M, qmfReal, synQmf->lsb, scaleFactorLowBand);
  scaleValuesSaturate(pTimeOut + M + synQmf->lsb, qmfReal + synQmf->lsb,
                      synQmf->usb - synQmf->lsb, scaleFactorHighBand);
  FDKmemclear(pTimeOut + M + synQmf->usb, (L - synQmf->usb) * sizeof(FIXP_DBL));

  dct_IV(pTimeOut + M, L, &shift);
  for (i = 0; i < M; i++) {
    pTimeOut[i] = pTimeOut[L - 1 - i];
    pTimeOut[2 * L - 1 - i] = -pTimeOut[L + i];
  }
}

#ifndef FUNCTION_qmfInverseModulationHQ
inline static void qmfInverseModulationHQ(
    HANDLE_QMF_FILTER_BANK synQmf,
    const FIXP_DBL *qmfReal,
    const FIXP_DBL *qmfImag,
    const int scaleFactorLowBand,
    const int scaleFactorHighBand,
    FIXP_DBL *pWorkBuffer
) {
  int i;
  int L = synQmf->no_channels;
  int M = L >> 1;
  int shift = 0;
  FIXP_DBL *RESTRICT tReal = pWorkBuffer;
  FIXP_DBL *RESTRICT tImag = pWorkBuffer + L;

  if (synQmf->flags & QMF_FLAG_CLDFB) {
    for (i = 0; i < synQmf->usb; i++) {
      cplxMultDiv2(&tImag[i], &tReal[i], qmfImag[i], qmfReal[i],
                   synQmf->t_cos[i], synQmf->t_sin[i]);
    }
    scaleValuesSaturate(&tReal[0], synQmf->lsb, scaleFactorLowBand + 1);
    scaleValuesSaturate(&tReal[0 + synQmf->lsb], synQmf->usb - synQmf->lsb,
                        scaleFactorHighBand + 1);
    scaleValuesSaturate(&tImag[0], synQmf->lsb, scaleFactorLowBand + 1);
    scaleValuesSaturate(&tImag[0 + synQmf->lsb], synQmf->usb - synQmf->lsb,
                        scaleFactorHighBand + 1);
  }

  if ((synQmf->flags & QMF_FLAG_CLDFB) == 0) {
    scaleValuesSaturate(&tReal[0], &qmfReal[0], synQmf->lsb,
                        scaleFactorLowBand);
    scaleValuesSaturate(&tReal[0 + synQmf->lsb], &qmfReal[0 + synQmf->lsb],
                        synQmf->usb - synQmf->lsb, scaleFactorHighBand);
    scaleValuesSaturate(&tImag[0], &qmfImag[0], synQmf->lsb,
                        scaleFactorLowBand);
    scaleValuesSaturate(&tImag[0 + synQmf->lsb], &qmfImag[0 + synQmf->lsb],
                        synQmf->usb - synQmf->lsb, scaleFactorHighBand);
  }

  FDKmemclear(&tReal[synQmf->usb],
              (synQmf->no_channels - synQmf->usb) * sizeof(FIXP_DBL));
  FDKmemclear(&tImag[synQmf->usb],
              (synQmf->no_channels - synQmf->usb) * sizeof(FIXP_DBL));

  dct_IV(tReal, L, &shift);
  dst_IV(tImag, L, &shift);

  if (synQmf->flags & QMF_FLAG_CLDFB) {
    for (i = 0; i < M; i++) {
      FIXP_DBL r1, i1, r2, i2;
      r1 = tReal[i];
      i2 = tImag[L - 1 - i];
      r2 = tReal[L - i - 1];
      i1 = tImag[i];

      tReal[i] = (r1 - i1) >> 1;
      tImag[L - 1 - i] = -(r1 + i1) >> 1;
      tReal[L - i - 1] = (r2 - i2) >> 1;
      tImag[i] = -(r2 + i2) >> 1;
    }
  } else {
    for (i = 0; i < M; i++) {
      FIXP_DBL r1, i1, r2, i2;
      r1 = -tReal[i];
      i2 = -tImag[L - 1 - i];
      r2 = -tReal[L - i - 1];
      i1 = -tImag[i];

      tReal[i] = (r1 - i1) >> 1;
      tImag[L - 1 - i] = -(r1 + i1) >> 1;
      tReal[L - i - 1] = (r2 - i2) >> 1;
      tImag[i] = -(r2 + i2) >> 1;
    }
  }
}
#endif

static int qmfInitFilterBank(
    HANDLE_QMF_FILTER_BANK h_Qmf,
    void *pFilterStates,
    int noCols,
    int lsb,
    int usb,
    int no_channels,
    UINT flags,
    int synflag)
{
  FDKmemclear(h_Qmf, sizeof(QMF_FILTER_BANK));

  if (flags & QMF_FLAG_MPSLDFB) {
    flags |= QMF_FLAG_NONSYMMETRIC;
    flags |= QMF_FLAG_MPSLDFB_OPTIMIZE_MODULATION;

    h_Qmf->t_cos = NULL;
    h_Qmf->t_sin = NULL;
    h_Qmf->filterScale = QMF_MPSLDFB_PFT_SCALE;
    h_Qmf->p_stride = 1;

    switch (no_channels) {
      case 64:
        h_Qmf->p_filter = qmf_mpsldfb_640;
        h_Qmf->FilterSize = 640;
        break;
      case 32:
        h_Qmf->p_filter = qmf_mpsldfb_320;
        h_Qmf->FilterSize = 320;
        break;
      default:
        return -1;
    }
  }

  if (!(flags & QMF_FLAG_MPSLDFB) && (flags & QMF_FLAG_CLDFB)) {
    flags |= QMF_FLAG_NONSYMMETRIC;
    h_Qmf->filterScale = QMF_CLDFB_PFT_SCALE;

    h_Qmf->p_stride = 1;
    switch (no_channels) {
      case 64:
        h_Qmf->t_cos = qmf_phaseshift_cos64_cldfb;
        h_Qmf->t_sin = qmf_phaseshift_sin64_cldfb;
        h_Qmf->p_filter = qmf_cldfb_640;
        h_Qmf->FilterSize = 640;
        break;
      case 32:
        h_Qmf->t_cos = (synflag) ? qmf_phaseshift_cos32_cldfb_syn
                                 : qmf_phaseshift_cos32_cldfb_ana;
        h_Qmf->t_sin = qmf_phaseshift_sin32_cldfb;
        h_Qmf->p_filter = qmf_cldfb_320;
        h_Qmf->FilterSize = 320;
        break;
      case 16:
        h_Qmf->t_cos = (synflag) ? qmf_phaseshift_cos16_cldfb_syn
                                 : qmf_phaseshift_cos16_cldfb_ana;
        h_Qmf->t_sin = qmf_phaseshift_sin16_cldfb;
        h_Qmf->p_filter = qmf_cldfb_160;
        h_Qmf->FilterSize = 160;
        break;
      case 8:
        h_Qmf->t_cos = (synflag) ? qmf_phaseshift_cos8_cldfb_syn
                                 : qmf_phaseshift_cos8_cldfb_ana;
        h_Qmf->t_sin = qmf_phaseshift_sin8_cldfb;
        h_Qmf->p_filter = qmf_cldfb_80;
        h_Qmf->FilterSize = 80;
        break;
      default:
        return -1;
    }
  }

  if (!(flags & QMF_FLAG_MPSLDFB) && ((flags & QMF_FLAG_CLDFB) == 0)) {
    switch (no_channels) {
      case 64:
        h_Qmf->p_filter = qmf_pfilt640;
        h_Qmf->t_cos = qmf_phaseshift_cos64;
        h_Qmf->t_sin = qmf_phaseshift_sin64;
        h_Qmf->p_stride = 1;
        h_Qmf->FilterSize = 640;
        h_Qmf->filterScale = 0;
        break;
      case 40:
        if (synflag) {
          break;
        } else {
          h_Qmf->p_filter = qmf_pfilt400;
          h_Qmf->t_cos = qmf_phaseshift_cos40;
          h_Qmf->t_sin = qmf_phaseshift_sin40;
          h_Qmf->filterScale = 1;
          h_Qmf->p_stride = 1;
          h_Qmf->FilterSize = no_channels * 10;
        }
        break;
      case 32:
        h_Qmf->p_filter = qmf_pfilt640;
        if (flags & QMF_FLAG_DOWNSAMPLED) {
          h_Qmf->t_cos = qmf_phaseshift_cos_downsamp32;
          h_Qmf->t_sin = qmf_phaseshift_sin_downsamp32;
        } else {
          h_Qmf->t_cos = qmf_phaseshift_cos32;
          h_Qmf->t_sin = qmf_phaseshift_sin32;
        }
        h_Qmf->p_stride = 2;
        h_Qmf->FilterSize = 640;
        h_Qmf->filterScale = 0;
        break;
      case 20:
        h_Qmf->p_filter = qmf_pfilt200;
        h_Qmf->p_stride = 1;
        h_Qmf->FilterSize = 200;
        h_Qmf->filterScale = 0;
        break;
      case 12:
        h_Qmf->p_filter = qmf_pfilt120;
        h_Qmf->p_stride = 1;
        h_Qmf->FilterSize = 120;
        h_Qmf->filterScale = 0;
        break;
      case 8:
        h_Qmf->p_filter = qmf_pfilt640;
        h_Qmf->p_stride = 8;
        h_Qmf->FilterSize = 640;
        h_Qmf->filterScale = 0;
        break;
      case 16:
        h_Qmf->p_filter = qmf_pfilt640;
        h_Qmf->t_cos = qmf_phaseshift_cos16;
        h_Qmf->t_sin = qmf_phaseshift_sin16;
        h_Qmf->p_stride = 4;
        h_Qmf->FilterSize = 640;
        h_Qmf->filterScale = 0;
        break;
      case 24:
        h_Qmf->p_filter = qmf_pfilt240;
        h_Qmf->t_cos = qmf_phaseshift_cos24;
        h_Qmf->t_sin = qmf_phaseshift_sin24;
        h_Qmf->p_stride = 1;
        h_Qmf->FilterSize = 240;
        h_Qmf->filterScale = 1;
        break;
      default:
        return -1;
    }
  }

  h_Qmf->synScalefactor = h_Qmf->filterScale;
  switch (no_channels) {
    case 128:
      h_Qmf->synScalefactor += ALGORITHMIC_SCALING_IN_SYNTHESIS_FILTERBANK + 1;
      break;
    case 40: {
      h_Qmf->synScalefactor += ALGORITHMIC_SCALING_IN_SYNTHESIS_FILTERBANK - 1;
    } break;
    case 64:
      h_Qmf->synScalefactor += ALGORITHMIC_SCALING_IN_SYNTHESIS_FILTERBANK;
      break;
    case 8:
      h_Qmf->synScalefactor += ALGORITHMIC_SCALING_IN_SYNTHESIS_FILTERBANK - 3;
      break;
    case 12:
      h_Qmf->synScalefactor += ALGORITHMIC_SCALING_IN_SYNTHESIS_FILTERBANK;
      break;
    case 20:
      h_Qmf->synScalefactor += ALGORITHMIC_SCALING_IN_SYNTHESIS_FILTERBANK + 1;
      break;
    case 32:
      h_Qmf->synScalefactor += ALGORITHMIC_SCALING_IN_SYNTHESIS_FILTERBANK - 1;
      break;
    case 16:
      h_Qmf->synScalefactor += ALGORITHMIC_SCALING_IN_SYNTHESIS_FILTERBANK - 2;
      break;
    case 24:
      h_Qmf->synScalefactor += ALGORITHMIC_SCALING_IN_SYNTHESIS_FILTERBANK - 1;
      break;
    default:
      return -1;
  }

  h_Qmf->flags = flags;

  h_Qmf->no_channels = no_channels;
  h_Qmf->no_col = noCols;

  h_Qmf->lsb = fMin(lsb, h_Qmf->no_channels);
  h_Qmf->usb = synflag
                   ? fMin(usb, h_Qmf->no_channels)
                   : usb;

  h_Qmf->FilterStates = (void *)pFilterStates;

  h_Qmf->outScalefactor =
      (ALGORITHMIC_SCALING_IN_ANALYSIS_FILTERBANK + h_Qmf->filterScale) +
      h_Qmf->synScalefactor;

  h_Qmf->outGain_m =
      (FIXP_DBL)0x80000000;
  h_Qmf->outGain_e = 0;

  return (0);
}

static inline void qmfAdaptFilterStates(
    HANDLE_QMF_FILTER_BANK synQmf,
    int scaleFactorDiff)
{
  if (synQmf == NULL || synQmf->FilterStates == NULL) {
    return;
  }
  if (scaleFactorDiff > 0) {
    scaleValuesSaturate((FIXP_QSS *)synQmf->FilterStates,
                        synQmf->no_channels * (QMF_NO_POLY * 2 - 1),
                        scaleFactorDiff);
  } else {
    scaleValues((FIXP_QSS *)synQmf->FilterStates,
                synQmf->no_channels * (QMF_NO_POLY * 2 - 1), scaleFactorDiff);
  }
}

int qmfInitSynthesisFilterBank(
    HANDLE_QMF_FILTER_BANK h_Qmf,
    FIXP_QSS *pFilterStates,
    int noCols,
    int lsb,
    int usb,
    int no_channels,
    int flags)
{
  int oldOutScale = h_Qmf->outScalefactor;
  int err = qmfInitFilterBank(h_Qmf, pFilterStates, noCols, lsb, usb,
                              no_channels, flags, 1);
  if (h_Qmf->FilterStates != NULL) {
    if (!(flags & QMF_FLAG_KEEP_STATES)) {
      FDKmemclear(
          h_Qmf->FilterStates,
          (2 * QMF_NO_POLY - 1) * h_Qmf->no_channels * sizeof(FIXP_QSS));
    } else {
      qmfAdaptFilterStates(h_Qmf, oldOutScale - h_Qmf->outScalefactor);
    }
  }

  FDK_ASSERT(h_Qmf->no_channels >= h_Qmf->lsb);
  FDK_ASSERT(h_Qmf->no_channels >= h_Qmf->usb);

  return err;
}

void qmfChangeOutScalefactor(
    HANDLE_QMF_FILTER_BANK synQmf,
    int outScalefactor
) {
  if (synQmf == NULL) {
    return;
  }

  outScalefactor +=
      (ALGORITHMIC_SCALING_IN_ANALYSIS_FILTERBANK + synQmf->filterScale) +
      synQmf->synScalefactor;

  if (synQmf->outScalefactor != outScalefactor) {
    int diff;

    diff = synQmf->outScalefactor - outScalefactor;

    qmfAdaptFilterStates(synQmf, diff);

    synQmf->outScalefactor = outScalefactor;
  }
}

int qmfGetOutScalefactor(
    HANDLE_QMF_FILTER_BANK synQmf)
{
  int scaleFactor = synQmf->outScalefactor
                        ? (synQmf->outScalefactor -
                           (ALGORITHMIC_SCALING_IN_ANALYSIS_FILTERBANK +
                            synQmf->filterScale + synQmf->synScalefactor))
                        : 0;
  return scaleFactor;
}

void qmfChangeOutGain(
    HANDLE_QMF_FILTER_BANK synQmf,
    FIXP_DBL outputGain,
    int outputGainScale
) {
  synQmf->outGain_m = outputGain;
  synQmf->outGain_e = outputGainScale;
}

#define INT_PCM_QMFOUT INT_PCM
#define SAMPLE_BITS_QMFOUT SAMPLE_BITS
#include "qmf_pcm.h"
#if SAMPLE_BITS == 16
#undef INT_PCM_QMFOUT
#undef SAMPLE_BITS_QMFOUT
#undef QMF_PCM_H
#undef FIXP_QAS
#undef QAS_BITS
#undef INT_PCM_QMFIN
#define INT_PCM_QMFOUT LONG
#define SAMPLE_BITS_QMFOUT 32
#define FIXP_QAS FIXP_DBL
#define QAS_BITS 32
#define INT_PCM_QMFIN LONG
#include "qmf_pcm.h"
#endif
