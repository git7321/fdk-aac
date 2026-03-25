/**************************** SBR encoder library ******************************

   Author(s):   M. Multrus

   Description: PS Wrapper, Downmix

*******************************************************************************/

#include "ps_main.h"
#include "ps_bitenc.h"
#include "sbrenc_ram.h"

static void psFindBestScaling(
    HANDLE_PARAMETRIC_STEREO hParametricStereo,
    FIXP_DBL *hybridData[HYBRID_FRAMESIZE][MAX_PS_CHANNELS][2],
    UCHAR *dynBandScale, FIXP_DBL *maxBandValue, SCHAR *dmxScale);

FDK_PSENC_ERROR PSEnc_Create(HANDLE_PARAMETRIC_STEREO *phParametricStereo) {
  FDK_PSENC_ERROR error = PSENC_OK;
  HANDLE_PARAMETRIC_STEREO hParametricStereo = NULL;

  if (phParametricStereo == NULL) {
    error = PSENC_INVALID_HANDLE;
  } else {
    int i;

    if (NULL == (hParametricStereo = GetRam_ParamStereo())) {
      error = PSENC_MEMORY_ERROR;
      goto bail;
    }
    FDKmemclear(hParametricStereo, sizeof(PARAMETRIC_STEREO));

    if (PSENC_OK !=
        (error = FDKsbrEnc_CreatePSEncode(&hParametricStereo->hPsEncode))) {
      error = PSENC_MEMORY_ERROR;
      goto bail;
    }

    for (i = 0; i < MAX_PS_CHANNELS; i++) {
      if (FDKhybridAnalysisOpen(
              &hParametricStereo->fdkHybAnaFilter[i],
              hParametricStereo->__staticHybAnaStatesLF[i],
              sizeof(hParametricStereo->__staticHybAnaStatesLF[i]),
              hParametricStereo->__staticHybAnaStatesHF[i],
              sizeof(hParametricStereo->__staticHybAnaStatesHF[i])) != 0) {
        error = PSENC_MEMORY_ERROR;
        goto bail;
      }
    }
  }

bail:
  if (phParametricStereo != NULL) {
    *phParametricStereo = hParametricStereo;
  }

  if (error != PSENC_OK) {
    PSEnc_Destroy(phParametricStereo);
  }
  return error;
}

FDK_PSENC_ERROR PSEnc_Init(HANDLE_PARAMETRIC_STEREO hParametricStereo,
                           const HANDLE_PSENC_CONFIG hPsEncConfig,
                           INT noQmfSlots, INT noQmfBands, UCHAR *dynamic_RAM) {
  FDK_PSENC_ERROR error = PSENC_OK;

  if ((NULL == hParametricStereo) || (NULL == hPsEncConfig)) {
    error = PSENC_INVALID_HANDLE;
  } else {
    int ch, i;

    hParametricStereo->initPS = 1;
    hParametricStereo->noQmfSlots = noQmfSlots;
    hParametricStereo->noQmfBands = noQmfBands;

    FDKmemclear(hParametricStereo->qmfDelayLines,
                sizeof(hParametricStereo->qmfDelayLines));

    hParametricStereo->qmfDelayScale = FRACT_BITS - 1;

    for (ch = 0; ch < MAX_PS_CHANNELS; ch++) {
      FDKhybridAnalysisInit(&hParametricStereo->fdkHybAnaFilter[ch],
                            THREE_TO_TEN, 64, 64, 1);
    }

    FDKhybridSynthesisInit(&hParametricStereo->fdkHybSynFilter, THREE_TO_TEN,
                           64, 64);

    hParametricStereo->psDelay =
        (HYBRID_FILTER_DELAY * hParametricStereo->noQmfBands);

    if ((hPsEncConfig->maxEnvelopes < PSENC_NENV_1) ||
        (hPsEncConfig->maxEnvelopes > PSENC_NENV_MAX)) {
      hPsEncConfig->maxEnvelopes = PSENC_NENV_DEFAULT;
    }
    hParametricStereo->maxEnvelopes = hPsEncConfig->maxEnvelopes;

    if (PSENC_OK !=
        (error = FDKsbrEnc_InitPSEncode(
             hParametricStereo->hPsEncode, (PS_BANDS)hPsEncConfig->nStereoBands,
             hPsEncConfig->iidQuantErrorThreshold))) {
      goto bail;
    }

    for (ch = 0; ch < MAX_PS_CHANNELS; ch++) {
      FIXP_DBL *pDynReal = GetRam_Sbr_envRBuffer(ch, dynamic_RAM);
      FIXP_DBL *pDynImag = GetRam_Sbr_envIBuffer(ch, dynamic_RAM);

      for (i = 0; i < HYBRID_FRAMESIZE; i++) {
        hParametricStereo->pHybridData[i + HYBRID_READ_OFFSET][ch][0] =
            &pDynReal[i * MAX_HYBRID_BANDS];
        hParametricStereo->pHybridData[i + HYBRID_READ_OFFSET][ch][1] =
            &pDynImag[i * MAX_HYBRID_BANDS];
        ;
      }

      for (i = 0; i < HYBRID_READ_OFFSET; i++) {
        hParametricStereo->pHybridData[i][ch][0] =
            hParametricStereo->__staticHybridData[i][ch][0];
        hParametricStereo->pHybridData[i][ch][1] =
            hParametricStereo->__staticHybridData[i][ch][1];
      }
    }

    FDKmemclear(hParametricStereo->__staticHybridData,
                sizeof(hParametricStereo->__staticHybridData));
    FDKmemclear(hParametricStereo->psOut, sizeof(hParametricStereo->psOut));

    hParametricStereo->psOut[0].enablePSHeader =
        1;

    FDKmemclear(hParametricStereo->dynBandScale, sizeof(UCHAR) * PS_MAX_BANDS);
    FDKmemclear(hParametricStereo->maxBandValue,
                sizeof(FIXP_DBL) * PS_MAX_BANDS);
  }
bail:
  return error;
}

FDK_PSENC_ERROR PSEnc_Destroy(HANDLE_PARAMETRIC_STEREO *phParametricStereo) {
  FDK_PSENC_ERROR error = PSENC_OK;

  if (NULL != phParametricStereo) {
    HANDLE_PARAMETRIC_STEREO hParametricStereo = *phParametricStereo;
    if (hParametricStereo != NULL) {
      FDKsbrEnc_DestroyPSEncode(&hParametricStereo->hPsEncode);
      FreeRam_ParamStereo(phParametricStereo);
    }
  }

  return error;
}

static FDK_PSENC_ERROR ExtractPSParameters(
    HANDLE_PARAMETRIC_STEREO hParametricStereo, const int sendHeader,
    FIXP_DBL *hybridData[HYBRID_FRAMESIZE][MAX_PS_CHANNELS][2]) {
  FDK_PSENC_ERROR error = PSENC_OK;

  if (hParametricStereo == NULL) {
    error = PSENC_INVALID_HANDLE;
  } else {
    if (hParametricStereo->initPS) {
      hParametricStereo->psOut[1] = hParametricStereo->psOut[0];
    }
    hParametricStereo->psOut[0] = hParametricStereo->psOut[1];

    if (PSENC_OK !=
        (error = FDKsbrEnc_PSEncode(
             hParametricStereo->hPsEncode, &hParametricStereo->psOut[1],
             hParametricStereo->dynBandScale, hParametricStereo->maxEnvelopes,
             hybridData, hParametricStereo->noQmfSlots, sendHeader))) {
      goto bail;
    }

    if (hParametricStereo->initPS) {
      hParametricStereo->psOut[0] = hParametricStereo->psOut[1];
      hParametricStereo->initPS = 0;
    }
  }
bail:
  return error;
}

static FDK_PSENC_ERROR DownmixPSQmfData(
    HANDLE_PARAMETRIC_STEREO hParametricStereo,
    HANDLE_QMF_FILTER_BANK sbrSynthQmf, FIXP_DBL **RESTRICT mixRealQmfData,
    FIXP_DBL **RESTRICT mixImagQmfData, INT_PCM *downsampledOutSignal,
    const UINT downsampledOutSignalBufSize,
    FIXP_DBL *hybridData[HYBRID_FRAMESIZE][MAX_PS_CHANNELS][2],
    const INT noQmfSlots, const INT psQmfScale[MAX_PS_CHANNELS],
    SCHAR *qmfScale) {
  FDK_PSENC_ERROR error = PSENC_OK;

  if (hParametricStereo == NULL) {
    error = PSENC_INVALID_HANDLE;
  } else {
    int n, k;
    C_AALLOC_SCRATCH_START(pWorkBuffer, FIXP_DBL, 2 * 64)

    int dynQmfScale = fixMax(
        0, hParametricStereo->dmxScale -
               1);
    int downmixScale = psQmfScale[0] - dynQmfScale;
    const FIXP_DBL maxStereoScaleFactor = MAXVAL_DBL;

    for (n = 0; n < noQmfSlots; n++) {
      FIXP_DBL tmpHybrid[2][MAX_HYBRID_BANDS];

      for (k = 0; k < 71; k++) {
        int dynScale, sc;
        FIXP_DBL tmpLeftReal, tmpRightReal, tmpLeftImag, tmpRightImag;
        FIXP_DBL tmpScaleFactor, stereoScaleFactor;

        tmpLeftReal = hybridData[n][0][0][k];
        tmpLeftImag = hybridData[n][0][1][k];
        tmpRightReal = hybridData[n][1][0][k];
        tmpRightImag = hybridData[n][1][1][k];

        sc = fixMax(
            0, CntLeadingZeros(fixMax(
                   fixMax(fixp_abs(tmpLeftReal), fixp_abs(tmpLeftImag)),
                   fixMax(fixp_abs(tmpRightReal), fixp_abs(tmpRightImag)))) -
                   2);

        tmpLeftReal <<= sc;
        tmpLeftImag <<= sc;
        tmpRightReal <<= sc;
        tmpRightImag <<= sc;
        dynScale = fixMin(sc - dynQmfScale, DFRACT_BITS - 1);

        stereoScaleFactor = fPow2Div2(tmpLeftReal) + fPow2Div2(tmpLeftImag) +
                            fPow2Div2(tmpRightReal) + fPow2Div2(tmpRightImag);

        tmpScaleFactor =
            fixp_abs(stereoScaleFactor + fMult(tmpLeftReal, tmpRightReal) +
                     fMult(tmpLeftImag, tmpRightImag));

        if ((stereoScaleFactor >> 1) <
            fMult(maxStereoScaleFactor, tmpScaleFactor)) {
          int sc_num = CountLeadingBits(stereoScaleFactor);
          int sc_denum = CountLeadingBits(tmpScaleFactor);
          sc = -(sc_num - sc_denum);

          tmpScaleFactor = schur_div((stereoScaleFactor << (sc_num)) >> 1,
                                     tmpScaleFactor << sc_denum, 16);

          if (sc & 0x1) {
            sc++;
            tmpScaleFactor >>= 1;
          }
          stereoScaleFactor = sqrtFixp(tmpScaleFactor);
          stereoScaleFactor <<= (sc >> 1);
        } else {
          stereoScaleFactor = maxStereoScaleFactor;
        }

        tmpHybrid[0][k] = fMultDiv2(stereoScaleFactor,
                                    (FIXP_DBL)(tmpLeftReal + tmpRightReal)) >>
                          dynScale;
        tmpHybrid[1][k] = fMultDiv2(stereoScaleFactor,
                                    (FIXP_DBL)(tmpLeftImag + tmpRightImag)) >>
                          dynScale;
      }

      FDKhybridSynthesisApply(&hParametricStereo->fdkHybSynFilter, tmpHybrid[0],
                              tmpHybrid[1], mixRealQmfData[n],
                              mixImagQmfData[n]);

      qmfSynthesisFilteringSlot(
          sbrSynthQmf, mixRealQmfData[n], mixImagQmfData[n], downmixScale - 7,
          downmixScale - 7,
          downsampledOutSignal + (n * sbrSynthQmf->no_channels), 1,
          pWorkBuffer);
    }

    *qmfScale = -downmixScale + 7;

    C_AALLOC_SCRATCH_END(pWorkBuffer, FIXP_DBL, 2 * 64)

    {
      const INT noQmfSlots2 = hParametricStereo->noQmfSlots >> 1;
      const int noQmfBands = hParametricStereo->noQmfBands;

      INT scale, i, j, slotOffset;

      FIXP_DBL tmp[2][64];

      for (i = 0; i < noQmfSlots2; i++) {
        FDKmemcpy(tmp[0], hParametricStereo->qmfDelayLines[0][i],
                  noQmfBands * sizeof(FIXP_DBL));
        FDKmemcpy(tmp[1], hParametricStereo->qmfDelayLines[1][i],
                  noQmfBands * sizeof(FIXP_DBL));

        FDKmemcpy(hParametricStereo->qmfDelayLines[0][i],
                  mixRealQmfData[i + noQmfSlots2],
                  noQmfBands * sizeof(FIXP_DBL));
        FDKmemcpy(hParametricStereo->qmfDelayLines[1][i],
                  mixImagQmfData[i + noQmfSlots2],
                  noQmfBands * sizeof(FIXP_DBL));

        FDKmemcpy(mixRealQmfData[i + noQmfSlots2], mixRealQmfData[i],
                  noQmfBands * sizeof(FIXP_DBL));
        FDKmemcpy(mixImagQmfData[i + noQmfSlots2], mixImagQmfData[i],
                  noQmfBands * sizeof(FIXP_DBL));

        FDKmemcpy(mixRealQmfData[i], tmp[0], noQmfBands * sizeof(FIXP_DBL));
        FDKmemcpy(mixImagQmfData[i], tmp[1], noQmfBands * sizeof(FIXP_DBL));
      }

      if (hParametricStereo->qmfDelayScale > *qmfScale) {
        scale = hParametricStereo->qmfDelayScale - *qmfScale;
        slotOffset = 0;
      } else {
        scale = *qmfScale - hParametricStereo->qmfDelayScale;
        slotOffset = noQmfSlots2;
      }

      for (i = 0; i < noQmfSlots2; i++) {
        for (j = 0; j < noQmfBands; j++) {
          mixRealQmfData[i + slotOffset][j] >>= scale;
          mixImagQmfData[i + slotOffset][j] >>= scale;
        }
      }

      scale = *qmfScale;
      *qmfScale = fMin(*qmfScale, hParametricStereo->qmfDelayScale);
      hParametricStereo->qmfDelayScale = scale;
    }
  }

  return error;
}

INT FDKsbrEnc_PSEnc_WritePSData(HANDLE_PARAMETRIC_STEREO hParametricStereo,
                                HANDLE_FDK_BITSTREAM hBitstream) {
  return (
      (hParametricStereo != NULL)
          ? FDKsbrEnc_WritePSBitstream(&hParametricStereo->psOut[0], hBitstream)
          : 0);
}

FDK_PSENC_ERROR FDKsbrEnc_PSEnc_ParametricStereoProcessing(
    HANDLE_PARAMETRIC_STEREO hParametricStereo, INT_PCM *samples[2],
    UINT samplesBufSize, QMF_FILTER_BANK **hQmfAnalysis,
    FIXP_DBL **RESTRICT downmixedRealQmfData,
    FIXP_DBL **RESTRICT downmixedImagQmfData, INT_PCM *downsampledOutSignal,
    HANDLE_QMF_FILTER_BANK sbrSynthQmf, SCHAR *qmfScale, const int sendHeader) {
  FDK_PSENC_ERROR error = PSENC_OK;
  INT psQmfScale[MAX_PS_CHANNELS] = {0};
  int psCh, i;
  C_AALLOC_SCRATCH_START(pWorkBuffer, FIXP_DBL, 4 * 64)

  for (psCh = 0; psCh < MAX_PS_CHANNELS; psCh++) {
    for (i = 0; i < hQmfAnalysis[psCh]->no_col; i++) {
      qmfAnalysisFilteringSlot(
          hQmfAnalysis[psCh], &pWorkBuffer[2 * 64],
          &pWorkBuffer[3 * 64],
          samples[psCh] + i * hQmfAnalysis[psCh]->no_channels, 1,
          &pWorkBuffer[0 * 64]
      );

      FDKhybridAnalysisApply(
          &hParametricStereo->fdkHybAnaFilter[psCh],
          &pWorkBuffer[2 * 64],
          &pWorkBuffer[3 * 64],
          hParametricStereo->pHybridData[i + HYBRID_READ_OFFSET][psCh][0],
          hParametricStereo->pHybridData[i + HYBRID_READ_OFFSET][psCh][1]);
    }

    psQmfScale[psCh] = hQmfAnalysis[psCh]->outScalefactor;
  }

  C_AALLOC_SCRATCH_END(pWorkBuffer, FIXP_DBL, 4 * 64)

  psFindBestScaling(
      hParametricStereo, &hParametricStereo->pHybridData[HYBRID_READ_OFFSET],
      hParametricStereo->dynBandScale, hParametricStereo->maxBandValue,
      &hParametricStereo->dmxScale);

  if (PSENC_OK !=
      (error = ExtractPSParameters(hParametricStereo, sendHeader,
                                   &hParametricStereo->pHybridData[0]))) {
    goto bail;
  }

  for (i = 0; i < HYBRID_READ_OFFSET; i++) {
    FDKmemcpy(
        hParametricStereo->pHybridData[i][0][0],
        hParametricStereo->pHybridData[hParametricStereo->noQmfSlots + i][0][0],
        MAX_HYBRID_BANDS * sizeof(FIXP_DBL)); /* left, real */
    FDKmemcpy(
        hParametricStereo->pHybridData[i][0][1],
        hParametricStereo->pHybridData[hParametricStereo->noQmfSlots + i][0][1],
        MAX_HYBRID_BANDS * sizeof(FIXP_DBL)); /* left, imag */
    FDKmemcpy(
        hParametricStereo->pHybridData[i][1][0],
        hParametricStereo->pHybridData[hParametricStereo->noQmfSlots + i][1][0],
        MAX_HYBRID_BANDS * sizeof(FIXP_DBL)); /* right, real */
    FDKmemcpy(
        hParametricStereo->pHybridData[i][1][1],
        hParametricStereo->pHybridData[hParametricStereo->noQmfSlots + i][1][1],
        MAX_HYBRID_BANDS * sizeof(FIXP_DBL)); /* right, imag */
  }

  if (PSENC_OK !=
      (error = DownmixPSQmfData(
           hParametricStereo, sbrSynthQmf, downmixedRealQmfData,
           downmixedImagQmfData, downsampledOutSignal, samplesBufSize,
           &hParametricStereo->pHybridData[HYBRID_READ_OFFSET],
           hParametricStereo->noQmfSlots, psQmfScale, qmfScale))) {
    goto bail;
  }

bail:

  return error;
}

static void psFindBestScaling(
    HANDLE_PARAMETRIC_STEREO hParametricStereo,
    FIXP_DBL *hybridData[HYBRID_FRAMESIZE][MAX_PS_CHANNELS][2],
    UCHAR *dynBandScale, FIXP_DBL *maxBandValue, SCHAR *dmxScale) {
  HANDLE_PS_ENCODE hPsEncode = hParametricStereo->hPsEncode;

  INT group, bin, col, band;
  const INT frameSize = hParametricStereo->noQmfSlots;
  const INT psBands = (INT)hPsEncode->psEncMode;
  const INT nIidGroups = hPsEncode->nQmfIidGroups + hPsEncode->nSubQmfIidGroups;

  FIXP_DBL maxVal[2][PS_MAX_BANDS];
  FIXP_DBL maxValue = FL2FXCONST_DBL(0.f);

  FDKmemclear(maxVal, sizeof(maxVal));

  for (group = 0; group < nIidGroups; group++) {
    bin = hPsEncode->subband2parameterIndex[group];

    if (hPsEncode->psEncMode == PS_BANDS_COARSE) {
      bin >>= 1;
    }

    for (col = 0; col < frameSize; col++) {
      int i, section = (col < frameSize - HYBRID_READ_OFFSET) ? 0 : 1;
      FIXP_DBL tmp = maxVal[section][bin];
      for (i = hPsEncode->iidGroupBorders[group];
           i < hPsEncode->iidGroupBorders[group + 1]; i++) {
        tmp = fixMax(tmp, (FIXP_DBL)fixp_abs(hybridData[col][0][0][i]));
        tmp = fixMax(tmp, (FIXP_DBL)fixp_abs(hybridData[col][0][1][i]));
        tmp = fixMax(tmp, (FIXP_DBL)fixp_abs(hybridData[col][1][0][i]));
        tmp = fixMax(tmp, (FIXP_DBL)fixp_abs(hybridData[col][1][1][i]));
      }
      maxVal[section][bin] = tmp;
    }
  }

  for (band = 0; band < psBands; band++) {
#ifndef MULT_16x16
    dynBandScale[band] =
        CountLeadingBits(fixMax(maxVal[0][band], maxBandValue[band]));
#else
    dynBandScale[band] = fixMax(
        0, CountLeadingBits(fixMax(maxVal[0][band], maxBandValue[band])) -
               FRACT_BITS);
#endif
    maxValue = fixMax(maxValue, fixMax(maxVal[0][band], maxVal[1][band]));
    maxBandValue[band] = fixMax(maxVal[0][band], maxVal[1][band]);
  }

#ifndef MULT_16x16
  *dmxScale = fixMin(DFRACT_BITS, CountLeadingBits(maxValue));
#else
  *dmxScale = fixMax(0, fixMin(FRACT_BITS, CountLeadingBits((maxValue))));
#endif
}
