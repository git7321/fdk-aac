/************************* MPEG-D DRC decoder library **************************

   Author(s):

   Description:

*******************************************************************************/

#include "drcDec_types.h"
#include "drcDec_gainDecoder.h"
#include "drcGainDec_process.h"

#define E_TGAINSTEP 12

static DRC_ERROR _prepareLnbIndex(ACTIVE_DRC* pActiveDrc,
                                  const int channelOffset,
                                  const int drcChannelOffset,
                                  const int numChannelsProcessed,
                                  const int lnbPointer) {
  int g, c;
  DRC_INSTRUCTIONS_UNI_DRC* pInst = pActiveDrc->pInst;

  if ((channelOffset + numChannelsProcessed) > 8) return DE_NOT_OK;

  if ((channelOffset + drcChannelOffset + numChannelsProcessed) > 8)
    return DE_NOT_OK;

  if ((channelOffset + drcChannelOffset) < 0) return DE_NOT_OK;

  for (c = channelOffset; c < channelOffset + numChannelsProcessed; c++) {
    if (pInst->drcSetId > 0) {
      int drcChannel = c + drcChannelOffset;
      if (drcChannel >= pInst->drcChannelCount) drcChannel = 0;
      g = pActiveDrc->channelGroupForChannel[drcChannel];
      if ((g >= 0) && !pActiveDrc->channelGroupIsParametricDrc[g]) {
        pActiveDrc->lnbIndexForChannel[c][lnbPointer] =
            pActiveDrc->activeDrcOffset + pActiveDrc->gainElementForGroup[g];
      }
    }
  }

  return DE_OK;
}

static DRC_ERROR _interpolateDrcGain(
    const GAIN_INTERPOLATION_TYPE gainInterpolationType,
    const SHORT timePrev,
    const SHORT tGainStep,
    const SHORT start, const SHORT stop, const SHORT stepsize,
    const FIXP_DBL gainLeft, const FIXP_DBL gainRight, const FIXP_DBL slopeLeft,
    const FIXP_DBL slopeRight, FIXP_DBL* buffer) {
  int n, n_buf;
  int start_modulo, start_offset;

  if (tGainStep < 0) {
    return DE_NOT_OK;
  }
  if (tGainStep == 0) {
    return DE_OK;
  }

  start_modulo = (start + timePrev) & (stepsize - 1);
  start_offset = (start_modulo ? stepsize - start_modulo : 0);
  n_buf = (start + timePrev + start_offset) >> (15 - fixnormz_S(stepsize));

  {
    LONG a;
    INT runs = (INT)(stop - start - start_offset + stepsize - 1) >>
               (30 - CountLeadingBits(stepsize));
    INT n_min = fMin(
        fMin(CntLeadingZeros(gainRight), CntLeadingZeros(gainLeft)) - 1, 8);
    a = (LONG)((gainRight << n_min) - (gainLeft << n_min)) / tGainStep;
    LONG a_step = a * stepsize;
    n = start + start_offset;
    a = a * n + (LONG)(gainLeft << n_min);
    buffer += n_buf;
#if defined(FUNCTION_interpolateDrcGain_func1)
    interpolateDrcGain_func1(buffer, a, a_step, n_min, runs);
#else
    a -= a_step;
    n_min = 8 - n_min;
    for (int i = 0; i < runs; i++) {
      a += a_step;
      buffer[i] = fMultDiv2(buffer[i], (FIXP_DBL)a) << n_min;
    }
#endif
  }
  return DE_OK;
}

static DRC_ERROR _processNodeSegments(
    const int frameSize, const GAIN_INTERPOLATION_TYPE gainInterpolationType,
    const int nNodes, const NODE_LIN* pNodeLin, const int offset,
    const SHORT stepsize,
    const NODE_LIN nodePrevious,
    const FIXP_DBL channelGain, FIXP_DBL* buffer) {
  DRC_ERROR err = DE_OK;
  SHORT timePrev, duration, start, stop, time;
  int n;
  FIXP_DBL gainLin = FL2FXCONST_DBL(1.0f / (float)(1 << 7)), gainLinPrev;
  FIXP_DBL slopeLin = (FIXP_DBL)0, slopeLinPrev = (FIXP_DBL)0;

  timePrev = nodePrevious.time + offset;
  gainLinPrev = nodePrevious.gainLin;
  for (n = 0; n < nNodes; n++) {
    time = pNodeLin[n].time + offset;
    duration = time - timePrev;
    gainLin = pNodeLin[n].gainLin;
    if (channelGain != FL2FXCONST_DBL(1.0f / (float)(1 << 8)))
      gainLin =
          SATURATE_LEFT_SHIFT(fMultDiv2(gainLin, channelGain), 9, DFRACT_BITS);

    if ((timePrev >= (frameSize - 1)) ||
        (time < 0)) {
      timePrev = time;
      gainLinPrev = gainLin;
      slopeLinPrev = slopeLin;
      continue;
    }

    start = fMax(-timePrev, 1);
    stop = fMin(time, (SHORT)(frameSize - 1)) - timePrev + 1;

    err = _interpolateDrcGain(gainInterpolationType, timePrev, duration, start,
                              stop, stepsize, gainLinPrev, gainLin,
                              slopeLinPrev, slopeLin, buffer);
    if (err) return err;

    timePrev = time;
    gainLinPrev = gainLin;
  }
  return err;
}

DRC_ERROR
processDrcTime(HANDLE_DRC_GAIN_DECODER hGainDec, const int activeDrcIndex,
               const int delaySamples, const int channelOffset,
               const int drcChannelOffset, const int numChannelsProcessed,
               const int timeDataChannelOffset, FIXP_DBL* deinterleavedAudio) {
  DRC_ERROR err = DE_OK;
  int c, b, i;
  ACTIVE_DRC* pActiveDrc = &(hGainDec->activeDrc[activeDrcIndex]);
  DRC_GAIN_BUFFERS* pDrcGainBuffers = &(hGainDec->drcGainBuffers);
  int lnbPointer = pDrcGainBuffers->lnbPointer, lnbIx;
  LINEAR_NODE_BUFFER* pLinearNodeBuffer = pDrcGainBuffers->linearNodeBuffer;
  LINEAR_NODE_BUFFER* pDummyLnb = &(pDrcGainBuffers->dummyLnb);
  int offset = 0;

  if (hGainDec->delayMode == DM_REGULAR_DELAY) {
    offset = hGainDec->frameSize;
  }

  if ((delaySamples + offset) >
      (NUM_LNB_FRAMES - 2) *
          hGainDec->frameSize)
    return DE_NOT_OK;

  err = _prepareLnbIndex(pActiveDrc, channelOffset, drcChannelOffset,
                         numChannelsProcessed, lnbPointer);
  if (err) return err;

  deinterleavedAudio +=
      channelOffset * timeDataChannelOffset;

  for (c = channelOffset; c < channelOffset + numChannelsProcessed; c++) {
    if (activeDrcIndex == hGainDec->channelGainActiveDrcIndex)
      pDrcGainBuffers->channelGain[c][lnbPointer] = hGainDec->channelGain[c];

    b = 0;
    {
      LINEAR_NODE_BUFFER *pLnb, *pLnbPrevious;
      NODE_LIN nodePrevious;
      int lnbPointerDiff;
      FIXP_DBL channelGain;
      lnbIx = lnbPointer + 1 - NUM_LNB_FRAMES;
      while (lnbIx < 0) lnbIx += NUM_LNB_FRAMES;

      if (activeDrcIndex == hGainDec->channelGainActiveDrcIndex)
        channelGain = pDrcGainBuffers->channelGain[c][lnbIx];
      else
        channelGain = FL2FXCONST_DBL(1.0f / (float)(1 << 8));

      for (i = 0; i < NUM_LNB_FRAMES - 1; i++) {
        if (pActiveDrc->lnbIndexForChannel[c][lnbIx] >= 0)
          pLnbPrevious = &(
              pLinearNodeBuffer[pActiveDrc->lnbIndexForChannel[c][lnbIx] + b]);
        else
          pLnbPrevious = pDummyLnb;
        nodePrevious =
            pLnbPrevious->linearNode[lnbIx][pLnbPrevious->nNodes[lnbIx] - 1];
        nodePrevious.time -= hGainDec->frameSize;
        if (channelGain != FL2FXCONST_DBL(1.0f / (float)(1 << 8)))
          nodePrevious.gainLin = SATURATE_LEFT_SHIFT(
              fMultDiv2(nodePrevious.gainLin,
                        pDrcGainBuffers->channelGain[c][lnbIx]),
              9, DFRACT_BITS);

        lnbIx++;
        if (lnbIx >= NUM_LNB_FRAMES) lnbIx = 0;

        if (pActiveDrc->lnbIndexForChannel[c][lnbIx] >= 0)
          pLnb = &(
              pLinearNodeBuffer[pActiveDrc->lnbIndexForChannel[c][lnbIx] + b]);
        else
          pLnb = pDummyLnb;

        if (activeDrcIndex == hGainDec->channelGainActiveDrcIndex)
          channelGain = pDrcGainBuffers->channelGain[c][lnbIx];

        lnbPointerDiff = i - (NUM_LNB_FRAMES - 2);

        err = _processNodeSegments(
            hGainDec->frameSize, pLnb->gainInterpolationType,
            pLnb->nNodes[lnbIx], pLnb->linearNode[lnbIx],
            lnbPointerDiff * hGainDec->frameSize + delaySamples + offset, 1,
            nodePrevious, channelGain, deinterleavedAudio);
        if (err) return err;
      }
      deinterleavedAudio += timeDataChannelOffset;
    }
  }
  return DE_OK;
}

DRC_ERROR
processDrcSubband(HANDLE_DRC_GAIN_DECODER hGainDec, const int activeDrcIndex,
                  const int delaySamples, const int channelOffset,
                  const int drcChannelOffset, const int numChannelsProcessed,
                  const int processSingleTimeslot,
                  FIXP_DBL* deinterleavedAudioReal[],
                  FIXP_DBL* deinterleavedAudioImag[]) {
  DRC_ERROR err = DE_OK;
  int b, c, g, m, m_start, m_stop, s, i;
  FIXP_DBL gainSb;
  DRC_INSTRUCTIONS_UNI_DRC* pInst = hGainDec->activeDrc[activeDrcIndex].pInst;
  DRC_GAIN_BUFFERS* pDrcGainBuffers = &(hGainDec->drcGainBuffers);
  ACTIVE_DRC* pActiveDrc = &(hGainDec->activeDrc[activeDrcIndex]);
  int activeDrcOffset = pActiveDrc->activeDrcOffset;
  int lnbPointer = pDrcGainBuffers->lnbPointer, lnbIx;
  LINEAR_NODE_BUFFER* pLinearNodeBuffer = pDrcGainBuffers->linearNodeBuffer;
  FIXP_DBL(*subbandGains)[4 * 1024 / 256] = hGainDec->subbandGains;
  FIXP_DBL* dummySubbandGains = hGainDec->dummySubbandGains;
  SUBBAND_DOMAIN_MODE subbandDomainMode = hGainDec->subbandDomainSupported;
  int signalIndex = 0;
  int frameSizeSb = 0;
  int nDecoderSubbands;
  SHORT L = 0;
  int offset = 0;
  FIXP_DBL *audioReal = NULL, *audioImag = NULL;

  if (hGainDec->delayMode == DM_REGULAR_DELAY) {
    offset = hGainDec->frameSize;
  }

  if ((delaySamples + offset) >
      (NUM_LNB_FRAMES - 2) *
          hGainDec->frameSize)
    return DE_NOT_OK;

  switch (subbandDomainMode) {
#if ((1024 / 256) >= (4096 / SUBBAND_DOWNSAMPLING_FACTOR_QMF64))
    case SDM_QMF64:
      nDecoderSubbands = SUBBAND_NUM_BANDS_QMF64;
      L = SUBBAND_DOWNSAMPLING_FACTOR_QMF64;
      break;
    case SDM_QMF71:
      nDecoderSubbands = SUBBAND_NUM_BANDS_QMF71;
      L = SUBBAND_DOWNSAMPLING_FACTOR_QMF71;
      break;
#else
    case SDM_QMF64:
    case SDM_QMF71:
      return DE_NOT_OK;
#endif
    case SDM_STFT256:
      nDecoderSubbands = SUBBAND_NUM_BANDS_STFT256;
      L = SUBBAND_DOWNSAMPLING_FACTOR_STFT256;
      break;
    default:
      return DE_NOT_OK;
  }

  frameSizeSb =
      hGainDec->frameSize >> (15 - fixnormz_S(L));

  if ((processSingleTimeslot < 0) || (processSingleTimeslot >= frameSizeSb)) {
    m_start = 0;
    m_stop = frameSizeSb;
  } else {
    m_start = processSingleTimeslot;
    m_stop = m_start + 1;
  }

  err = _prepareLnbIndex(pActiveDrc, channelOffset, drcChannelOffset,
                         numChannelsProcessed, lnbPointer);
  if (err) return err;

  if (!pActiveDrc->subbandGainsReady)
  {
    for (g = 0; g < pInst->nDrcChannelGroups; g++) {
      b = 0;
      {
        LINEAR_NODE_BUFFER* pLnb =
            &(pLinearNodeBuffer[activeDrcOffset +
                                pActiveDrc->gainElementForGroup[g] + b]);
        NODE_LIN nodePrevious;
        int lnbPointerDiff;

        for (m = 0; m < frameSizeSb; m++) {
          subbandGains[activeDrcOffset + g][b * frameSizeSb + m] =
              FL2FXCONST_DBL(1.0f / (float)(1 << 7));
        }

        lnbIx = lnbPointer - (NUM_LNB_FRAMES - 1);
        while (lnbIx < 0) lnbIx += NUM_LNB_FRAMES;

        for (i = 0; i < NUM_LNB_FRAMES - 1; i++) {
          nodePrevious = pLnb->linearNode[lnbIx][pLnb->nNodes[lnbIx] - 1];
          nodePrevious.time -= hGainDec->frameSize;

          lnbIx++;
          if (lnbIx >= NUM_LNB_FRAMES) lnbIx = 0;

          lnbPointerDiff = i - (NUM_LNB_FRAMES - 2);

          err = _processNodeSegments(
              hGainDec->frameSize, pLnb->gainInterpolationType,
              pLnb->nNodes[lnbIx], pLnb->linearNode[lnbIx],
              lnbPointerDiff * hGainDec->frameSize + delaySamples + offset -
                  (L - 1) / 2,
              L, nodePrevious, FL2FXCONST_DBL(1.0f / (float)(1 << 8)),
              &(subbandGains[activeDrcOffset + g][b * frameSizeSb]));
          if (err) return err;
        }
      }
    }
    pActiveDrc->subbandGainsReady = 1;
  }

  for (c = channelOffset; c < channelOffset + numChannelsProcessed; c++) {
    FIXP_DBL* thisSubbandGainsBuffer;
    if (pInst->drcSetId > 0)
      g = pActiveDrc->channelGroupForChannel[c + drcChannelOffset];
    else
      g = -1;

    audioReal = deinterleavedAudioReal[signalIndex];
    if (subbandDomainMode != SDM_STFT256) {
      audioImag = deinterleavedAudioImag[signalIndex];
    }

    if ((g >= 0) && !pActiveDrc->channelGroupIsParametricDrc[g]) {
      thisSubbandGainsBuffer = subbandGains[activeDrcOffset + g];
    } else {
      thisSubbandGainsBuffer = dummySubbandGains;
    }

    for (m = m_start; m < m_stop; m++) {
      INT n_min = 8;
      {
        gainSb = thisSubbandGainsBuffer[m];
        if (activeDrcIndex == hGainDec->channelGainActiveDrcIndex)
          gainSb = SATURATE_LEFT_SHIFT(
              fMultDiv2(gainSb, hGainDec->channelGain[c]), 9, DFRACT_BITS);
        n_min = fMin(CntLeadingZeros(gainSb) - 1, n_min);
        gainSb <<= n_min;
        n_min = 8 - n_min;
        if (subbandDomainMode ==
            SDM_STFT256) {
          for (s = 0; s < nDecoderSubbands; s++) {
            *audioReal = fMultDiv2(*audioReal, gainSb) << n_min;
            audioReal++;
            *audioReal = fMultDiv2(*audioReal, gainSb) << n_min;
            audioReal++;
          }
        } else {
          for (s = 0; s < nDecoderSubbands; s++) {
            *audioReal = fMultDiv2(*audioReal, gainSb) << n_min;
            audioReal++;
            *audioImag = fMultDiv2(*audioImag, gainSb) << n_min;
            audioImag++;
          }
        }
      }
    }
    signalIndex++;
  }
  return DE_OK;
}
