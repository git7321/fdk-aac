/******************* Library for basic calculation routines ********************

   Author(s):   Markus Lohwasser, Josef Hoepfl, Manuel Jander

   Description: QMF filterbank

*******************************************************************************/

#ifndef QMF_PCM_H
#define QMF_PCM_H

#define QSSCALE (0)
#define FX_DBL2FX_QSS(x) (x)
#define FX_QSS2FX_DBL(x) (x)

#ifndef FUNCTION_qmfSynPrototypeFirSlot
void qmfSynPrototypeFirSlot(
#else
void qmfSynPrototypeFirSlot_fallback(
#endif
    HANDLE_QMF_FILTER_BANK qmf,
    FIXP_DBL *RESTRICT realSlot,
    FIXP_DBL *RESTRICT imagSlot,
    INT_PCM_QMFOUT *RESTRICT timeOut,
    int stride) {
  FIXP_QSS *FilterStates = (FIXP_QSS *)qmf->FilterStates;
  int no_channels = qmf->no_channels;
  const FIXP_PFT *p_Filter = qmf->p_filter;
  int p_stride = qmf->p_stride;
  int j;
  FIXP_QSS *RESTRICT sta = FilterStates;
  const FIXP_PFT *RESTRICT p_flt, *RESTRICT p_fltm;
  int scale = (DFRACT_BITS - SAMPLE_BITS_QMFOUT) - 1 - qmf->outScalefactor -
              qmf->outGain_e;

  p_flt = p_Filter + p_stride * QMF_NO_POLY;
  p_fltm = p_Filter + (qmf->FilterSize / 2) -
           p_stride * QMF_NO_POLY;

  FIXP_SGL gain = FX_DBL2FX_SGL(qmf->outGain_m);

  FIXP_DBL rnd_val = 0;

  if (scale > 0) {
    if (scale < (DFRACT_BITS - 1))
      rnd_val = FIXP_DBL(1 << (scale - 1));
    else
      scale = (DFRACT_BITS - 1);
  } else {
    scale = fMax(scale, -(DFRACT_BITS - 1));
  }

  for (j = no_channels - 1; j >= 0; j--) {
    FIXP_DBL imag = imagSlot[j];
    FIXP_DBL real = realSlot[j];
    {
      INT_PCM_QMFOUT tmp;
      FIXP_DBL Are = fMultAddDiv2(FX_QSS2FX_DBL(sta[0]), p_fltm[0], real);

      if (gain != (FIXP_SGL)(-32768))
      {
        Are = fMult(Are, gain);
      }
      if (scale >= 0) {
        FDK_ASSERT(
            Are <=
            (Are + rnd_val));
        tmp = (INT_PCM_QMFOUT)(
            SATURATE_RIGHT_SHIFT(Are + rnd_val, scale, SAMPLE_BITS_QMFOUT));
      } else {
        tmp = (INT_PCM_QMFOUT)(
            SATURATE_LEFT_SHIFT(Are, -scale, SAMPLE_BITS_QMFOUT));
      }

      { timeOut[(j)*stride] = tmp; }
    }

    sta[0] = FX_DBL2FX_QSS(fMultAddDiv2(FX_QSS2FX_DBL(sta[1]), p_flt[4], imag));
    sta[1] =
        FX_DBL2FX_QSS(fMultAddDiv2(FX_QSS2FX_DBL(sta[2]), p_fltm[1], real));
    sta[2] = FX_DBL2FX_QSS(fMultAddDiv2(FX_QSS2FX_DBL(sta[3]), p_flt[3], imag));
    sta[3] =
        FX_DBL2FX_QSS(fMultAddDiv2(FX_QSS2FX_DBL(sta[4]), p_fltm[2], real));
    sta[4] = FX_DBL2FX_QSS(fMultAddDiv2(FX_QSS2FX_DBL(sta[5]), p_flt[2], imag));
    sta[5] =
        FX_DBL2FX_QSS(fMultAddDiv2(FX_QSS2FX_DBL(sta[6]), p_fltm[3], real));
    sta[6] = FX_DBL2FX_QSS(fMultAddDiv2(FX_QSS2FX_DBL(sta[7]), p_flt[1], imag));
    sta[7] =
        FX_DBL2FX_QSS(fMultAddDiv2(FX_QSS2FX_DBL(sta[8]), p_fltm[4], real));
    sta[8] = FX_DBL2FX_QSS(fMultDiv2(p_flt[0], imag));
    p_flt += (p_stride * QMF_NO_POLY);
    p_fltm -= (p_stride * QMF_NO_POLY);
    sta += 9;
  }
}

#ifndef FUNCTION_qmfSynPrototypeFirSlot_NonSymmetric
static void qmfSynPrototypeFirSlot_NonSymmetric(
    HANDLE_QMF_FILTER_BANK qmf,
    FIXP_DBL *RESTRICT realSlot,
    FIXP_DBL *RESTRICT imagSlot,
    INT_PCM_QMFOUT *RESTRICT timeOut,
    int stride) {
  FIXP_QSS *FilterStates = (FIXP_QSS *)qmf->FilterStates;
  int no_channels = qmf->no_channels;
  const FIXP_PFT *p_Filter = qmf->p_filter;
  int p_stride = qmf->p_stride;
  int j;
  FIXP_QSS *RESTRICT sta = FilterStates;
  const FIXP_PFT *RESTRICT p_flt, *RESTRICT p_fltm;
  int scale = (DFRACT_BITS - SAMPLE_BITS_QMFOUT) - 1 - qmf->outScalefactor -
              qmf->outGain_e;

  p_flt = p_Filter;
  p_fltm =
      &p_flt[qmf->FilterSize / 2];

  FIXP_SGL gain = FX_DBL2FX_SGL(qmf->outGain_m);

  FIXP_DBL rnd_val = (FIXP_DBL)0;

  if (scale > 0) {
    if (scale < (DFRACT_BITS - 1))
      rnd_val = FIXP_DBL(1 << (scale - 1));
    else
      scale = (DFRACT_BITS - 1);
  } else {
    scale = fMax(scale, -(DFRACT_BITS - 1));
  }

  for (j = no_channels - 1; j >= 0; j--) {
    FIXP_DBL imag = imagSlot[j];
    FIXP_DBL real = realSlot[j];
    {
      INT_PCM_QMFOUT tmp;
      FIXP_DBL Are = sta[0] + FX_DBL2FX_QSS(fMultDiv2(p_fltm[4], real));

      if (gain != (FIXP_SGL)(-32768))
      {
        Are = fMult(Are, gain);
      }
      if (scale > 0) {
        FDK_ASSERT(Are <
                   (Are + rnd_val));
        tmp = (INT_PCM_QMFOUT)(
            SATURATE_RIGHT_SHIFT(Are + rnd_val, scale, SAMPLE_BITS_QMFOUT));
      } else {
        tmp = (INT_PCM_QMFOUT)(
            SATURATE_LEFT_SHIFT(Are, -scale, SAMPLE_BITS_QMFOUT));
      }
      timeOut[j * stride] = tmp;
    }

    sta[0] = sta[1] + FX_DBL2FX_QSS(fMultDiv2(p_flt[4], imag));
    sta[1] = sta[2] + FX_DBL2FX_QSS(fMultDiv2(p_fltm[3], real));
    sta[2] = sta[3] + FX_DBL2FX_QSS(fMultDiv2(p_flt[3], imag));

    sta[3] = sta[4] + FX_DBL2FX_QSS(fMultDiv2(p_fltm[2], real));
    sta[4] = sta[5] + FX_DBL2FX_QSS(fMultDiv2(p_flt[2], imag));
    sta[5] = sta[6] + FX_DBL2FX_QSS(fMultDiv2(p_fltm[1], real));
    sta[6] = sta[7] + FX_DBL2FX_QSS(fMultDiv2(p_flt[1], imag));

    sta[7] = sta[8] + FX_DBL2FX_QSS(fMultDiv2(p_fltm[0], real));
    sta[8] = FX_DBL2FX_QSS(fMultDiv2(p_flt[0], imag));

    p_flt += (p_stride * QMF_NO_POLY);
    p_fltm += (p_stride * QMF_NO_POLY);
    sta += 9;
  }
}
#endif /* FUNCTION_qmfSynPrototypeFirSlot_NonSymmetric */

void qmfSynthesisFilteringSlot(HANDLE_QMF_FILTER_BANK synQmf,
                               const FIXP_DBL *realSlot,
                               const FIXP_DBL *imagSlot,
                               const int scaleFactorLowBand,
                               const int scaleFactorHighBand,
                               INT_PCM_QMFOUT *timeOut, const int stride,
                               FIXP_DBL *pWorkBuffer) {
  if (!(synQmf->flags & QMF_FLAG_LP))
    qmfInverseModulationHQ(synQmf, realSlot, imagSlot, scaleFactorLowBand,
                           scaleFactorHighBand, pWorkBuffer);
  else {
    if (synQmf->flags & QMF_FLAG_CLDFB) {
      qmfInverseModulationLP_odd(synQmf, realSlot, scaleFactorLowBand,
                                 scaleFactorHighBand, pWorkBuffer);
    } else {
      qmfInverseModulationLP_even(synQmf, realSlot, scaleFactorLowBand,
                                  scaleFactorHighBand, pWorkBuffer);
    }
  }

  if (synQmf->flags & QMF_FLAG_NONSYMMETRIC) {
    qmfSynPrototypeFirSlot_NonSymmetric(synQmf, pWorkBuffer,
                                        pWorkBuffer + synQmf->no_channels,
                                        timeOut, stride);
  } else {
    qmfSynPrototypeFirSlot(synQmf, pWorkBuffer,
                           pWorkBuffer + synQmf->no_channels, timeOut, stride);
  }
}

void qmfSynthesisFiltering(
    HANDLE_QMF_FILTER_BANK synQmf,
    FIXP_DBL **QmfBufferReal,
    FIXP_DBL **QmfBufferImag,
    const QMF_SCALE_FACTOR *scaleFactor,
    const INT ov_len,
    INT_PCM_QMFOUT *timeOut,
    const INT stride,
    FIXP_DBL *pWorkBuffer
) {
  int i;
  int L = synQmf->no_channels;
  int scaleFactorHighBand;
  int scaleFactorLowBand_ov, scaleFactorLowBand_no_ov;

  FDK_ASSERT(synQmf->no_channels >= synQmf->lsb);
  FDK_ASSERT(synQmf->no_channels >= synQmf->usb);

  scaleFactorHighBand = -ALGORITHMIC_SCALING_IN_ANALYSIS_FILTERBANK -
                        scaleFactor->hb_scale - synQmf->filterScale;
  scaleFactorLowBand_ov = -ALGORITHMIC_SCALING_IN_ANALYSIS_FILTERBANK -
                          scaleFactor->ov_lb_scale - synQmf->filterScale;
  scaleFactorLowBand_no_ov = -ALGORITHMIC_SCALING_IN_ANALYSIS_FILTERBANK -
                             scaleFactor->lb_scale - synQmf->filterScale;

  for (i = 0; i < synQmf->no_col; i++)
  {
    const FIXP_DBL *QmfBufferImagSlot = NULL;

    int scaleFactorLowBand =
        (i < ov_len) ? scaleFactorLowBand_ov : scaleFactorLowBand_no_ov;

    if (!(synQmf->flags & QMF_FLAG_LP)) QmfBufferImagSlot = QmfBufferImag[i];

    qmfSynthesisFilteringSlot(synQmf, QmfBufferReal[i], QmfBufferImagSlot,
                              scaleFactorLowBand, scaleFactorHighBand,
                              timeOut + (i * L * stride), stride, pWorkBuffer);
  }
}

int qmfInitAnalysisFilterBank(
    HANDLE_QMF_FILTER_BANK h_Qmf,
    FIXP_QAS *pFilterStates,
    int noCols,
    int lsb,
    int usb,
    int no_channels,
    int flags)
{
  int err = qmfInitFilterBank(h_Qmf, pFilterStates, noCols, lsb, usb,
                              no_channels, flags, 0);
  if (!(flags & QMF_FLAG_KEEP_STATES) && (h_Qmf->FilterStates != NULL)) {
    FDKmemclear(h_Qmf->FilterStates,
                (2 * QMF_NO_POLY - 1) * h_Qmf->no_channels * sizeof(FIXP_QAS));
  }

  FDK_ASSERT(h_Qmf->no_channels >= h_Qmf->lsb);

  return err;
}

#ifndef FUNCTION_qmfAnaPrototypeFirSlot
static void qmfAnaPrototypeFirSlot(
    FIXP_DBL *analysisBuffer,
    INT no_channels,
    const FIXP_PFT *p_filter, INT p_stride,
    FIXP_QAS *RESTRICT pFilterStates) {
  INT k;

  FIXP_DBL accu;
  const FIXP_PFT *RESTRICT p_flt = p_filter;
  FIXP_DBL *RESTRICT pData_0 = analysisBuffer + 2 * no_channels - 1;
  FIXP_DBL *RESTRICT pData_1 = analysisBuffer;

  FIXP_QAS *RESTRICT sta_0 = (FIXP_QAS *)pFilterStates;
  FIXP_QAS *RESTRICT sta_1 =
      (FIXP_QAS *)pFilterStates + (2 * QMF_NO_POLY * no_channels) - 1;
  INT pfltStep = QMF_NO_POLY * (p_stride);
  INT staStep1 = no_channels << 1;
  INT staStep2 = (no_channels << 3) - 1;

  for (k = 0; k < no_channels; k++) {
    accu = fMultDiv2(p_flt[0], *sta_1);
    sta_1 -= staStep1;
    accu += fMultDiv2(p_flt[1], *sta_1);
    sta_1 -= staStep1;
    accu += fMultDiv2(p_flt[2], *sta_1);
    sta_1 -= staStep1;
    accu += fMultDiv2(p_flt[3], *sta_1);
    sta_1 -= staStep1;
    accu += fMultDiv2(p_flt[4], *sta_1);
    *pData_1++ = (accu << 1);
    sta_1 += staStep2;

    p_flt += pfltStep;
    accu = fMultDiv2(p_flt[0], *sta_0);
    sta_0 += staStep1;
    accu += fMultDiv2(p_flt[1], *sta_0);
    sta_0 += staStep1;
    accu += fMultDiv2(p_flt[2], *sta_0);
    sta_0 += staStep1;
    accu += fMultDiv2(p_flt[3], *sta_0);
    sta_0 += staStep1;
    accu += fMultDiv2(p_flt[4], *sta_0);
    *pData_0-- = (accu << 1);
    sta_0 -= staStep2;
  }
}
#endif /* !defined(FUNCTION_qmfAnaPrototypeFirSlot) */

#ifndef FUNCTION_qmfAnaPrototypeFirSlot_NonSymmetric
static void qmfAnaPrototypeFirSlot_NonSymmetric(
    FIXP_DBL *analysisBuffer,
    int no_channels,
    const FIXP_PFT *p_filter, int p_stride,
    FIXP_QAS *RESTRICT pFilterStates) {
  const FIXP_PFT *RESTRICT p_flt = p_filter;
  int p, k;

  for (k = 0; k < 2 * no_channels; k++) {
    FIXP_DBL accu = (FIXP_DBL)0;

    p_flt += QMF_NO_POLY * (p_stride - 1);

    for (p = 0; p < QMF_NO_POLY; p++) {
      accu += fMultDiv2(*p_flt++, pFilterStates[2 * no_channels * p]);
    }
    analysisBuffer[2 * no_channels - 1 - k] = (accu << 1);
    pFilterStates++;
  }
}
#endif /* FUNCTION_qmfAnaPrototypeFirSlot_NonSymmetric */

void qmfAnalysisFilteringSlot(
    HANDLE_QMF_FILTER_BANK anaQmf,
    FIXP_DBL *qmfReal,
    FIXP_DBL *qmfImag,
    const INT_PCM_QMFIN *RESTRICT timeIn,
    const int stride,
    FIXP_DBL *pWorkBuffer
) {
  int offset = anaQmf->no_channels * (QMF_NO_POLY * 2 - 1);
  {
    FIXP_QAS *FilterStatesAnaTmp = ((FIXP_QAS *)anaQmf->FilterStates) + offset;

    for (int i = anaQmf->no_channels >> 1; i != 0; i--) {
      *FilterStatesAnaTmp++ = (FIXP_QAS)*timeIn;
      timeIn += stride;
      *FilterStatesAnaTmp++ = (FIXP_QAS)*timeIn;
      timeIn += stride;
    }
  }

  if (anaQmf->flags & QMF_FLAG_NONSYMMETRIC) {
    qmfAnaPrototypeFirSlot_NonSymmetric(pWorkBuffer, anaQmf->no_channels,
                                        anaQmf->p_filter, anaQmf->p_stride,
                                        (FIXP_QAS *)anaQmf->FilterStates);
  } else {
    qmfAnaPrototypeFirSlot(pWorkBuffer, anaQmf->no_channels, anaQmf->p_filter,
                           anaQmf->p_stride, (FIXP_QAS *)anaQmf->FilterStates);
  }

  if (anaQmf->flags & QMF_FLAG_LP) {
    if (anaQmf->flags & QMF_FLAG_CLDFB)
      qmfForwardModulationLP_odd(anaQmf, pWorkBuffer, qmfReal);
    else
      qmfForwardModulationLP_even(anaQmf, pWorkBuffer, qmfReal);

  } else {
    qmfForwardModulationHQ(anaQmf, pWorkBuffer, qmfReal, qmfImag);
  }
  FDKmemmove(anaQmf->FilterStates,
             (FIXP_QAS *)anaQmf->FilterStates + anaQmf->no_channels,
             offset * sizeof(FIXP_QAS));
}

void qmfAnalysisFiltering(
    HANDLE_QMF_FILTER_BANK anaQmf,
    FIXP_DBL **qmfReal,
    FIXP_DBL **qmfImag,
    QMF_SCALE_FACTOR *scaleFactor,
    const INT_PCM_QMFIN *timeIn,
    const int timeIn_e, const int stride,
    FIXP_DBL *pWorkBuffer
) {
  int i;
  int no_channels = anaQmf->no_channels;

  scaleFactor->lb_scale =
      -ALGORITHMIC_SCALING_IN_ANALYSIS_FILTERBANK - timeIn_e;
  scaleFactor->lb_scale -= anaQmf->filterScale;

  for (i = 0; i < anaQmf->no_col; i++) {
    FIXP_DBL *qmfImagSlot = NULL;

    if (!(anaQmf->flags & QMF_FLAG_LP)) {
      qmfImagSlot = qmfImag[i];
    }

    qmfAnalysisFilteringSlot(anaQmf, qmfReal[i], qmfImagSlot, timeIn, stride,
                             pWorkBuffer);

    timeIn += no_channels * stride;
  }
}
#endif /* QMF_PCM_H */
