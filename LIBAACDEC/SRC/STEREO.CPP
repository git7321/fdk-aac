/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: joint stereo processing

*******************************************************************************/

#include "stereo.h"

#include "aac_rom.h"
#include "FDK_bitstream.h"
#include "channelinfo.h"
#include "FDK_audio.h"

enum { L = 0, R = 1 };

#include "block.h"

int CJointStereo_Read(HANDLE_FDK_BITSTREAM bs,
                      CJointStereoData *pJointStereoData,
                      const int windowGroups,
                      const int scaleFactorBandsTransmitted,
                      const int max_sfb_ste_clear,
                      CJointStereoPersistentData *pJointStereoPersistentData,
                      CCplxPredictionData *cplxPredictionData,
                      int cplxPredictionActiv, int scaleFactorBandsTotal,
                      int windowSequence, const UINT flags) {
  int group, band;

  pJointStereoData->MsMaskPresent = (UCHAR)FDKreadBits(bs, 2);

  FDKmemclear(pJointStereoData->MsUsed,
              scaleFactorBandsTransmitted * sizeof(UCHAR));

  pJointStereoData->cplx_pred_flag = 0;
  if (cplxPredictionActiv) {
    cplxPredictionData->pred_dir = 0;
    cplxPredictionData->complex_coef = 0;
    cplxPredictionData->use_prev_frame = 0;
    cplxPredictionData->igf_pred_dir = 0;
  }

  switch (pJointStereoData->MsMaskPresent) {
    case 0:
      break;

    case 1:
      for (group = 0; group < windowGroups; group++) {
        for (band = 0; band < scaleFactorBandsTransmitted; band++) {
          pJointStereoData->MsUsed[band] |= (FDKreadBits(bs, 1) << group);
        }
      }
      break;

    case 2:
      for (band = 0; band < scaleFactorBandsTransmitted; band++) {
        pJointStereoData->MsUsed[band] = 255;
      }
      break;

    case 3:
      if (flags & (AC_USAC | AC_RSVD50 | AC_RSV603DA)) {
        if (cplxPredictionActiv) {

          pJointStereoData->cplx_pred_flag = 1;

          int cplx_pred_all = 0;
          cplx_pred_all = FDKreadBits(bs, 1);

          if (cplx_pred_all) {
            for (group = 0; group < windowGroups; group++) {
              UCHAR groupmask = ((UCHAR)1 << group);
              for (band = 0; band < scaleFactorBandsTransmitted; band++) {
                pJointStereoData->MsUsed[band] |= groupmask;
              }
            }
          } else {
            for (group = 0; group < windowGroups; group++) {
              for (band = 0; band < scaleFactorBandsTransmitted;
                   band += SFB_PER_PRED_BAND) {
                pJointStereoData->MsUsed[band] |= (FDKreadBits(bs, 1) << group);
                if ((band + 1) < scaleFactorBandsTotal) {
                  pJointStereoData->MsUsed[band + 1] |=
                      (pJointStereoData->MsUsed[band] & ((UCHAR)1 << group));
                }
              }
            }
          }
        } else {
          return -1;
        }
      }
      break;
  }

  if (cplxPredictionActiv) {
    if (pJointStereoData->MsMaskPresent == 3) {
      if (pJointStereoData->cplx_pred_flag) {
        int delta_code_time = 0;

        const CodeBookDescription *hcb = &AACcodeBookDescriptionTable[BOOKSCL];
        if (((windowSequence == BLOCK_SHORT) &&
             (pJointStereoPersistentData->winSeqPrev != BLOCK_SHORT)) ||
            ((pJointStereoPersistentData->winSeqPrev == BLOCK_SHORT) &&
             (windowSequence != BLOCK_SHORT))) {
          FDKmemclear(pJointStereoPersistentData->alpha_q_re_prev,
                      JointStereoMaximumGroups * JointStereoMaximumBands *
                          sizeof(SHORT));
          FDKmemclear(pJointStereoPersistentData->alpha_q_im_prev,
                      JointStereoMaximumGroups * JointStereoMaximumBands *
                          sizeof(SHORT));
        }
        {
          FDKmemclear(cplxPredictionData->alpha_q_re,
                      JointStereoMaximumGroups * JointStereoMaximumBands *
                          sizeof(SHORT));
          FDKmemclear(cplxPredictionData->alpha_q_im,
                      JointStereoMaximumGroups * JointStereoMaximumBands *
                          sizeof(SHORT));
        }

        cplxPredictionData->pred_dir = FDKreadBits(bs, 1);
        cplxPredictionData->complex_coef = FDKreadBits(bs, 1);

        if (cplxPredictionData->complex_coef) {
          if (flags & AC_INDEP) {
            cplxPredictionData->use_prev_frame = 0;
          } else {
            cplxPredictionData->use_prev_frame = FDKreadBits(bs, 1);
          }
        }

        if (flags & AC_INDEP) {
          delta_code_time = 0;
        } else {
          delta_code_time = FDKreadBits(bs, 1);
        }

        {
          int last_alpha_q_re = 0, last_alpha_q_im = 0;

          for (group = 0; group < windowGroups; group++) {
            for (band = 0; band < scaleFactorBandsTransmitted;
                 band += SFB_PER_PRED_BAND) {
              if (delta_code_time == 1) {
                if (group > 0) {
                  last_alpha_q_re =
                      cplxPredictionData->alpha_q_re[group - 1][band];
                  last_alpha_q_im =
                      cplxPredictionData->alpha_q_im[group - 1][band];
                } else if ((windowSequence == BLOCK_SHORT) &&
                           (pJointStereoPersistentData->winSeqPrev ==
                            BLOCK_SHORT)) {
                  if (pJointStereoPersistentData->winGroupsPrev == 0) return -1;

                  last_alpha_q_re =
                      pJointStereoPersistentData->alpha_q_re_prev
                          [pJointStereoPersistentData->winGroupsPrev - 1][band];
                  last_alpha_q_im =
                      pJointStereoPersistentData->alpha_q_im_prev
                          [pJointStereoPersistentData->winGroupsPrev - 1][band];
                } else {
                  last_alpha_q_re =
                      pJointStereoPersistentData->alpha_q_re_prev[group][band];
                  last_alpha_q_im =
                      pJointStereoPersistentData->alpha_q_im_prev[group][band];
                }

              } else {
                if (band > 0) {
                  last_alpha_q_re =
                      cplxPredictionData->alpha_q_re[group][band - 1];
                  last_alpha_q_im =
                      cplxPredictionData->alpha_q_im[group][band - 1];
                } else {
                  last_alpha_q_re = 0;
                  last_alpha_q_im = 0;
                }
              }

              if (pJointStereoData->MsUsed[band] & ((UCHAR)1 << group)) {
                int dpcm_alpha_re, dpcm_alpha_im;

                dpcm_alpha_re = CBlock_DecodeHuffmanWord(bs, hcb);
                dpcm_alpha_re -= 60;
                dpcm_alpha_re *= -1;

                cplxPredictionData->alpha_q_re[group][band] =
                    dpcm_alpha_re + last_alpha_q_re;

                if (cplxPredictionData->complex_coef) {
                  dpcm_alpha_im = CBlock_DecodeHuffmanWord(bs, hcb);
                  dpcm_alpha_im -= 60;
                  dpcm_alpha_im *= -1;

                  cplxPredictionData->alpha_q_im[group][band] =
                      dpcm_alpha_im + last_alpha_q_im;
                } else {
                  cplxPredictionData->alpha_q_im[group][band] = 0;
                }

              } else {
                cplxPredictionData->alpha_q_re[group][band] = 0;
                cplxPredictionData->alpha_q_im[group][band] = 0;
              }

              if ((band + 1) <
                  scaleFactorBandsTransmitted) {
                cplxPredictionData->alpha_q_re[group][band + 1] =
                    cplxPredictionData->alpha_q_re[group][band];
                cplxPredictionData->alpha_q_im[group][band + 1] =
                    cplxPredictionData->alpha_q_im[group][band];
              }

              pJointStereoPersistentData->alpha_q_re_prev[group][band] =
                  cplxPredictionData->alpha_q_re[group][band];
              pJointStereoPersistentData->alpha_q_im_prev[group][band] =
                  cplxPredictionData->alpha_q_im[group][band];
            }

            for (band = scaleFactorBandsTransmitted; band < max_sfb_ste_clear;
                 band++) {
              cplxPredictionData->alpha_q_re[group][band] = 0;
              cplxPredictionData->alpha_q_im[group][band] = 0;
              pJointStereoPersistentData->alpha_q_re_prev[group][band] = 0;
              pJointStereoPersistentData->alpha_q_im_prev[group][band] = 0;
            }
          }
        }
      }
    } else {
      for (group = 0; group < windowGroups; group++) {
        for (band = 0; band < max_sfb_ste_clear; band++) {
          pJointStereoPersistentData->alpha_q_re_prev[group][band] = 0;
          pJointStereoPersistentData->alpha_q_im_prev[group][band] = 0;
        }
      }
    }

    pJointStereoPersistentData->winGroupsPrev = windowGroups;
  }

  return 0;
}

static void CJointStereo_filterAndAdd(
    FIXP_DBL *in, int len, int windowLen, const FIXP_FILT *coeff, FIXP_DBL *out,
    UCHAR isCurrent
) {
  int i, j;

  int indices_1[] = {2, 1, 0, 1, 2, 3};
  int indices_2[] = {1, 0, 0, 2, 3, 4};
  int indices_3[] = {0, 0, 1, 3, 4, 5};

  int subtr_1[] = {6, 5, 4, 2, 1, 1};
  int subtr_2[] = {5, 4, 3, 1, 1, 2};
  int subtr_3[] = {4, 3, 2, 1, 2, 3};

  if (isCurrent == 1) {

    for (i = 0; i < 3; i++) {
      out[0] -= (FIXP_DBL)fMultDiv2(coeff[i], in[indices_1[i]]) >> SR_FNA_OUT;
      out[0] +=
          (FIXP_DBL)fMultDiv2(coeff[i], in[indices_1[5 - i]]) >> SR_FNA_OUT;
    }

    for (i = 0; i < 3; i++) {
      out[1] -= (FIXP_DBL)fMultDiv2(coeff[i], in[indices_2[i]]) >> SR_FNA_OUT;
      out[1] +=
          (FIXP_DBL)fMultDiv2(coeff[i], in[indices_2[5 - i]]) >> SR_FNA_OUT;
    }

    for (i = 0; i < 3; i++) {
      out[2] -= (FIXP_DBL)fMultDiv2(coeff[i], in[indices_3[i]]) >> SR_FNA_OUT;
      out[2] +=
          (FIXP_DBL)fMultDiv2(coeff[i], in[indices_3[5 - i]]) >> SR_FNA_OUT;
    }

    for (j = 3; j < (len - 3); j++) {
      for (i = 0; i < 3; i++) {
        out[j] -= (FIXP_DBL)fMultDiv2(coeff[i], in[j - 3 + i]) >> SR_FNA_OUT;
        out[j] += (FIXP_DBL)fMultDiv2(coeff[i], in[j + 3 - i]) >> SR_FNA_OUT;
      }
    }

    for (i = 0; i < 3; i++) {
      out[len - 3] -=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_1[i]]) >> SR_FNA_OUT;
      out[len - 3] +=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_1[5 - i]]) >> SR_FNA_OUT;
    }

    for (i = 0; i < 3; i++) {
      out[len - 2] -=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_2[i]]) >> SR_FNA_OUT;
      out[len - 2] +=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_2[5 - i]]) >> SR_FNA_OUT;
    }

    for (i = 0; i < 3; i++) {
      out[len - 1] -=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_3[i]]) >> SR_FNA_OUT;
      out[len - 1] +=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_3[5 - i]]) >> SR_FNA_OUT;
    }

  } else {

    for (i = 0; i < 3; i++) {
      out[0] -= (FIXP_DBL)fMultDiv2(coeff[i], in[indices_1[i]] >> SR_FNA_OUT);
      out[0] -=
          (FIXP_DBL)fMultDiv2(coeff[i], in[indices_1[5 - i]] >> SR_FNA_OUT);
    }
    out[0] -= (FIXP_DBL)fMultDiv2(coeff[3], in[0] >> SR_FNA_OUT);

    for (i = 0; i < 3; i++) {
      out[1] += (FIXP_DBL)fMultDiv2(coeff[i], in[indices_2[i]] >> SR_FNA_OUT);
      out[1] +=
          (FIXP_DBL)fMultDiv2(coeff[i], in[indices_2[5 - i]] >> SR_FNA_OUT);
    }
    out[1] += (FIXP_DBL)fMultDiv2(coeff[3], in[1] >> SR_FNA_OUT);

    for (i = 0; i < 3; i++) {
      out[2] -= (FIXP_DBL)fMultDiv2(coeff[i], in[indices_3[i]] >> SR_FNA_OUT);
      out[2] -=
          (FIXP_DBL)fMultDiv2(coeff[i], in[indices_3[5 - i]] >> SR_FNA_OUT);
    }
    out[2] -= (FIXP_DBL)fMultDiv2(coeff[3], in[2] >> SR_FNA_OUT);

    for (j = 3; j < (len - 4); j++) {
      for (i = 0; i < 3; i++) {
        out[j] += (FIXP_DBL)fMultDiv2(coeff[i], in[j - 3 + i] >> SR_FNA_OUT);
        out[j] += (FIXP_DBL)fMultDiv2(coeff[i], in[j + 3 - i] >> SR_FNA_OUT);
      }
      out[j] += (FIXP_DBL)fMultDiv2(coeff[3], in[j] >> SR_FNA_OUT);

      j++;

      for (i = 0; i < 3; i++) {
        out[j] -= (FIXP_DBL)fMultDiv2(coeff[i], in[j - 3 + i] >> SR_FNA_OUT);
        out[j] -= (FIXP_DBL)fMultDiv2(coeff[i], in[j + 3 - i] >> SR_FNA_OUT);
      }
      out[j] -= (FIXP_DBL)fMultDiv2(coeff[3], in[j] >> SR_FNA_OUT);
    }

    for (i = 0; i < 3; i++) {
      out[len - 3] +=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_1[i]] >> SR_FNA_OUT);
      out[len - 3] +=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_1[5 - i]] >> SR_FNA_OUT);
    }
    out[len - 3] += (FIXP_DBL)fMultDiv2(coeff[3], in[len - 3] >> SR_FNA_OUT);

    for (i = 0; i < 3; i++) {
      out[len - 2] -=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_2[i]] >> SR_FNA_OUT);
      out[len - 2] -=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_2[5 - i]] >> SR_FNA_OUT);
    }
    out[len - 2] -= (FIXP_DBL)fMultDiv2(coeff[3], in[len - 2] >> SR_FNA_OUT);

    for (i = 0; i < 3; i++) {
      out[len - 1] +=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_3[i]] >> SR_FNA_OUT);
      out[len - 1] +=
          (FIXP_DBL)fMultDiv2(coeff[i], in[len - subtr_3[5 - i]] >> SR_FNA_OUT);
    }
    out[len - 1] += (FIXP_DBL)fMultDiv2(coeff[3], in[len - 1] >> SR_FNA_OUT);
  }
}

static inline void CJointStereo_GenerateMSOutput(FIXP_DBL *pSpecLCurrBand,
                                                 FIXP_DBL *pSpecRCurrBand,
                                                 UINT leftScale,
                                                 UINT rightScale,
                                                 UINT nSfbBands) {
  unsigned int i;

  FIXP_DBL leftCoefficient0;
  FIXP_DBL leftCoefficient1;
  FIXP_DBL leftCoefficient2;
  FIXP_DBL leftCoefficient3;

  FIXP_DBL rightCoefficient0;
  FIXP_DBL rightCoefficient1;
  FIXP_DBL rightCoefficient2;
  FIXP_DBL rightCoefficient3;

  for (i = nSfbBands; i > 0; i -= 4) {
    leftCoefficient0 = pSpecLCurrBand[i - 4];
    leftCoefficient1 = pSpecLCurrBand[i - 3];
    leftCoefficient2 = pSpecLCurrBand[i - 2];
    leftCoefficient3 = pSpecLCurrBand[i - 1];

    rightCoefficient0 = pSpecRCurrBand[i - 4];
    rightCoefficient1 = pSpecRCurrBand[i - 3];
    rightCoefficient2 = pSpecRCurrBand[i - 2];
    rightCoefficient3 = pSpecRCurrBand[i - 1];

    leftCoefficient0 >>= leftScale;
    leftCoefficient1 >>= leftScale;
    leftCoefficient2 >>= leftScale;
    leftCoefficient3 >>= leftScale;

    rightCoefficient0 >>= rightScale;
    rightCoefficient1 >>= rightScale;
    rightCoefficient2 >>= rightScale;
    rightCoefficient3 >>= rightScale;

    pSpecLCurrBand[i - 4] = leftCoefficient0 + rightCoefficient0;
    pSpecLCurrBand[i - 3] = leftCoefficient1 + rightCoefficient1;
    pSpecLCurrBand[i - 2] = leftCoefficient2 + rightCoefficient2;
    pSpecLCurrBand[i - 1] = leftCoefficient3 + rightCoefficient3;

    pSpecRCurrBand[i - 4] = leftCoefficient0 - rightCoefficient0;
    pSpecRCurrBand[i - 3] = leftCoefficient1 - rightCoefficient1;
    pSpecRCurrBand[i - 2] = leftCoefficient2 - rightCoefficient2;
    pSpecRCurrBand[i - 1] = leftCoefficient3 - rightCoefficient3;
  }
}

void CJointStereo_ApplyMS(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo[2],
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[2],
    FIXP_DBL *spectrumL, FIXP_DBL *spectrumR, SHORT *SFBleftScale,
    SHORT *SFBrightScale, SHORT *specScaleL, SHORT *specScaleR,
    const SHORT *pScaleFactorBandOffsets, const UCHAR *pWindowGroupLength,
    const int windowGroups, const int max_sfb_ste_outside,
    const int scaleFactorBandsTransmittedL,
    const int scaleFactorBandsTransmittedR, FIXP_DBL *store_dmx_re_prev,
    SHORT *store_dmx_re_prev_e, const int mainband_flag) {
  int window, group, band;
  UCHAR groupMask;
  CJointStereoData *pJointStereoData =
      &pAacDecoderChannelInfo[L]->pComData->jointStereoData;
  CCplxPredictionData *cplxPredictionData =
      pAacDecoderChannelInfo[L]->pComStaticData->cplxPredictionData;

  int max_sfb_ste =
      fMax(scaleFactorBandsTransmittedL, scaleFactorBandsTransmittedR);
  int min_sfb_ste =
      fMin(scaleFactorBandsTransmittedL, scaleFactorBandsTransmittedR);
  int scaleFactorBandsTransmitted = min_sfb_ste;

  if (pJointStereoData->cplx_pred_flag) {
    int windowLen, groupwin, frameMaxScale;
    CJointStereoPersistentData *pJointStereoPersistentData =
        &pAacDecoderStaticChannelInfo[L]
             ->pCpeStaticData->jointStereoPersistentData;
    FIXP_DBL *const staticSpectralCoeffsL =
        pAacDecoderStaticChannelInfo[L]
            ->pCpeStaticData->jointStereoPersistentData.spectralCoeffs[L];
    FIXP_DBL *const staticSpectralCoeffsR =
        pAacDecoderStaticChannelInfo[L]
            ->pCpeStaticData->jointStereoPersistentData.spectralCoeffs[R];
    SHORT *const staticSpecScaleL =
        pAacDecoderStaticChannelInfo[L]
            ->pCpeStaticData->jointStereoPersistentData.specScale[L];
    SHORT *const staticSpecScaleR =
        pAacDecoderStaticChannelInfo[L]
            ->pCpeStaticData->jointStereoPersistentData.specScale[R];

    FIXP_DBL *dmx_re =
        pAacDecoderStaticChannelInfo[L]
            ->pCpeStaticData->jointStereoPersistentData.scratchBuffer;
    FIXP_DBL *dmx_re_prev =
        pAacDecoderStaticChannelInfo[L]
            ->pCpeStaticData->jointStereoPersistentData.scratchBuffer +
        1024;

    SHORT dmx_re_prev_e = *store_dmx_re_prev_e;

    const FIXP_FILT *pCoeff;
    const FIXP_FILT *pCoeffPrev;
    int coeffPointerOffset;

    int previousShape = (int)pJointStereoPersistentData->winShapePrev;
    int currentShape = (int)pAacDecoderChannelInfo[L]->icsInfo.WindowShape;

    C_AALLOC_SCRATCH_START(dmx_im, FIXP_DBL, 1024);

    windowLen =
        pAacDecoderChannelInfo[L]->granuleLength;

    if (pAacDecoderChannelInfo[L]->icsInfo.WindowSequence != BLOCK_SHORT) {
      windowLen *= 8;
    }

    pCoeffPrev = mdst_filt_coef_prev[previousShape];

    if (cplxPredictionData->complex_coef == 1) {
      switch (pAacDecoderChannelInfo[L]
                  ->icsInfo.WindowSequence) {
        case BLOCK_SHORT:
        case BLOCK_LONG:
          pCoeffPrev = mdst_filt_coef_prev[previousShape];
          break;

        case BLOCK_START:
          if ((pJointStereoPersistentData->winSeqPrev == BLOCK_SHORT) ||
              (pJointStereoPersistentData->winSeqPrev == BLOCK_START)) {
            pCoeffPrev = mdst_filt_coef_prev[2 + previousShape];
          } else {
            pCoeffPrev = mdst_filt_coef_prev[previousShape];
          }
          break;

        case BLOCK_STOP:
          pCoeffPrev = mdst_filt_coef_prev[2 + previousShape];
          break;

        default:
          pCoeffPrev = mdst_filt_coef_prev[previousShape];
          break;
      }
    }

    if ((previousShape == SHAPE_SINE) && (currentShape == SHAPE_SINE)) {
      coeffPointerOffset = 0;
    } else if ((previousShape == SHAPE_SINE) && (currentShape == SHAPE_KBD)) {
      coeffPointerOffset = 2;
    } else if ((previousShape == SHAPE_KBD) && (currentShape == SHAPE_KBD)) {
      coeffPointerOffset = 1;
    } else
    {
      coeffPointerOffset = 3;
    }

    switch (pAacDecoderChannelInfo[L]
                ->icsInfo.WindowSequence) {
      case BLOCK_SHORT:
      case BLOCK_LONG:
        pCoeff = mdst_filt_coef_curr[coeffPointerOffset];
        break;

      case BLOCK_START:
        if ((pJointStereoPersistentData->winSeqPrev == BLOCK_SHORT) ||
            (pJointStereoPersistentData->winSeqPrev == BLOCK_START)) {
          pCoeff = mdst_filt_coef_curr[12 + coeffPointerOffset];
        } else {
          pCoeff = mdst_filt_coef_curr[4 + coeffPointerOffset];
        }
        break;

      case BLOCK_STOP:
        pCoeff = mdst_filt_coef_curr[8 + coeffPointerOffset];
        break;

      default:
        pCoeff = mdst_filt_coef_curr[coeffPointerOffset];
    }

    frameMaxScale = 0;
    for (window = 0, group = 0; group < windowGroups; group++) {
      for (groupwin = 0; groupwin < pWindowGroupLength[group];
           groupwin++, window++) {
        SHORT *leftScale = &SFBleftScale[window * 16];
        SHORT *rightScale = &SFBrightScale[window * 16];
        int windowMaxScale = 0;

        for (band = 0; band < min_sfb_ste; band++) {
          int lScale = leftScale[band];
          int rScale = rightScale[band];
          int commonScale = ((lScale > rScale) ? lScale : rScale);
          windowMaxScale =
              (windowMaxScale < commonScale) ? commonScale : windowMaxScale;
        }
        if (scaleFactorBandsTransmittedL >
            min_sfb_ste) {
          for (; band < max_sfb_ste; band++) {
            int lScale = leftScale[band];
            windowMaxScale =
                (windowMaxScale < lScale) ? lScale : windowMaxScale;
          }
        } else {
          if (scaleFactorBandsTransmittedR >
              min_sfb_ste) {
            for (; band < max_sfb_ste; band++) {
              int rScale = rightScale[band];
              windowMaxScale =
                  (windowMaxScale < rScale) ? rScale : windowMaxScale;
            }
          }
        }

        frameMaxScale =
            (frameMaxScale < windowMaxScale) ? windowMaxScale : frameMaxScale;
      }
    }

    frameMaxScale += 2;

    for (window = 0, group = 0; group < windowGroups; group++) {
      groupMask = 1 << group;

      for (groupwin = 0; groupwin < pWindowGroupLength[group];
           groupwin++, window++) {
        FDKmemclear(&dmx_im[windowLen * window], windowLen * sizeof(FIXP_DBL));

        if (cplxPredictionData->complex_coef == 1) {
          if ((cplxPredictionData->use_prev_frame == 1) && (mainband_flag)) {
            if ((pAacDecoderChannelInfo[L]->icsInfo.WindowSequence !=
                 BLOCK_SHORT) ||
                (window == 0)) {
              int index_offset = 0;
              int srLeftChan = 0;
              int srRightChan = 0;
              if (pAacDecoderChannelInfo[L]->icsInfo.WindowSequence ==
                  BLOCK_SHORT) {
                index_offset = windowLen * 7;
                if (staticSpecScaleL[7] > staticSpecScaleR[7]) {
                  srRightChan = staticSpecScaleL[7] - staticSpecScaleR[7];
                  dmx_re_prev_e = staticSpecScaleL[7];
                } else {
                  srLeftChan = staticSpecScaleR[7] - staticSpecScaleL[7];
                  dmx_re_prev_e = staticSpecScaleR[7];
                }
              } else {
                if (staticSpecScaleL[0] > staticSpecScaleR[0]) {
                  srRightChan = staticSpecScaleL[0] - staticSpecScaleR[0];
                  dmx_re_prev_e = staticSpecScaleL[0];
                } else {
                  srLeftChan = staticSpecScaleR[0] - staticSpecScaleL[0];
                  dmx_re_prev_e = staticSpecScaleR[0];
                }
              }

              if (pAacDecoderStaticChannelInfo[L]
                      ->pCpeStaticData->jointStereoPersistentData
                      .clearSpectralCoeffs == 1) {
                FDKmemclear(dmx_re_prev, windowLen * sizeof(FIXP_DBL));
                dmx_re_prev_e = 0;
              } else {
                if (cplxPredictionData->pred_dir == 0) {
                  for (int i = 0; i < windowLen; i++) {
                    dmx_re_prev[i] =
                        ((staticSpectralCoeffsL[index_offset + i] >>
                          fMin(DFRACT_BITS - 1, srLeftChan + 1)) +
                         (staticSpectralCoeffsR[index_offset + i] >>
                          fMin(DFRACT_BITS - 1, srRightChan + 1)));
                  }
                } else {
                  for (int i = 0; i < windowLen; i++) {
                    dmx_re_prev[i] =
                        ((staticSpectralCoeffsL[index_offset + i] >>
                          fMin(DFRACT_BITS - 1, srLeftChan + 1)) -
                         (staticSpectralCoeffsR[index_offset + i] >>
                          fMin(DFRACT_BITS - 1, srRightChan + 1)));
                  }
                }
              }

              FDKmemcpy(store_dmx_re_prev, &dmx_re_prev[0],
                        windowLen * sizeof(FIXP_DBL));

              *store_dmx_re_prev_e = dmx_re_prev_e;
            }
          }
        }

        SHORT *leftScale = &SFBleftScale[window * 16];
        SHORT *rightScale = &SFBrightScale[window * 16];

        specScaleL[window] = specScaleR[window] = frameMaxScale;

        if (cplxPredictionData->use_prev_frame == 1) {
          if (window == 0) {
            if (dmx_re_prev_e < frameMaxScale) {
              if (mainband_flag == 0) {
                scaleValues(
                    dmx_re_prev, store_dmx_re_prev, windowLen,
                    -fMin(DFRACT_BITS - 1, (frameMaxScale - dmx_re_prev_e)));
              } else {
                scaleValues(
                    dmx_re_prev, windowLen,
                    -fMin(DFRACT_BITS - 1, (frameMaxScale - dmx_re_prev_e)));
              }
            } else {
              if (mainband_flag == 0) {
                FDKmemcpy(dmx_re_prev, store_dmx_re_prev,
                          windowLen * sizeof(FIXP_DBL));
              }
              specScaleL[0] = dmx_re_prev_e;
              specScaleR[0] = dmx_re_prev_e;
            }
          } else {
            FDK_ASSERT(pAacDecoderChannelInfo[L]->icsInfo.WindowSequence ==
                       BLOCK_SHORT);
            if (specScaleL[window - 1] < frameMaxScale) {
              scaleValues(&dmx_re[windowLen * (window - 1)], windowLen,
                          -fMin(DFRACT_BITS - 1,
                                (frameMaxScale - specScaleL[window - 1])));
            } else {
              specScaleL[window] = specScaleL[window - 1];
              specScaleR[window] = specScaleR[window - 1];
            }
          }
        }

        FDK_ASSERT(specScaleL[window] == specScaleR[window]);

        for (band = 0; band < max_sfb_ste; band++) {
          int lScale = leftScale[band];
          int rScale = rightScale[band];

          lScale = fMin(DFRACT_BITS - 1, specScaleL[window] - lScale);
          rScale = fMin(DFRACT_BITS - 1,
                        specScaleL[window] - rScale);

          leftScale[band] = rightScale[band] = specScaleL[window];

          for (int i = pScaleFactorBandOffsets[band];
               i < pScaleFactorBandOffsets[band + 1]; i++) {
            spectrumL[windowLen * window + i] >>= lScale;
            spectrumR[windowLen * window + i] >>= rScale;
          }

          if (pJointStereoData->MsUsed[band] & groupMask) {
            for (int i = pScaleFactorBandOffsets[band];
                 i < pScaleFactorBandOffsets[band + 1]; i++) {
              dmx_re[windowLen * window + i] =
                  spectrumL[windowLen * window + i];
            }
          } else {
            if (cplxPredictionData->pred_dir == 0) {
              for (int i = pScaleFactorBandOffsets[band];
                   i < pScaleFactorBandOffsets[band + 1]; i++) {
                dmx_re[windowLen * window + i] =
                    (spectrumL[windowLen * window + i] +
                     spectrumR[windowLen * window + i]) >>
                    1;
              }
            } else {
              for (int i = pScaleFactorBandOffsets[band];
                   i < pScaleFactorBandOffsets[band + 1]; i++) {
                dmx_re[windowLen * window + i] =
                    (spectrumL[windowLen * window + i] -
                     spectrumR[windowLen * window + i]) >>
                    1;
              }
            }
          }
        }
        for (int i = pScaleFactorBandOffsets[max_sfb_ste_outside];
             i < windowLen; i++) {
          dmx_re[windowLen * window + i] = (FIXP_DBL)0;
        }

        if (cplxPredictionData->complex_coef == 1) {
          {
            if ((pAacDecoderChannelInfo[L]->icsInfo.WindowSequence ==
                 BLOCK_SHORT) &&
                (window != 0)) {
              pCoeff = mdst_filt_coef_curr[currentShape];
              pCoeffPrev = mdst_filt_coef_prev[currentShape];
            }

            int extended_band = fMin(
                pScaleFactorBandOffsets[max_sfb_ste_outside] + 7, windowLen);

            if ((pAacDecoderChannelInfo[L]->icsInfo.WindowSequence ==
                 BLOCK_SHORT) &&
                (window != 0)) {
              CJointStereo_filterAndAdd(&dmx_re[windowLen * window],
                                        extended_band, windowLen, pCoeff,
                                        &dmx_im[windowLen * window], 1);

              CJointStereo_filterAndAdd(&dmx_re[windowLen * (window - 1)],
                                        extended_band, windowLen, pCoeffPrev,
                                        &dmx_im[windowLen * window], 0);
            } else {
              CJointStereo_filterAndAdd(dmx_re, extended_band, windowLen,
                                        pCoeff, dmx_im, 1);

              if (cplxPredictionData->use_prev_frame == 1) {
                CJointStereo_filterAndAdd(dmx_re_prev, extended_band, windowLen,
                                          pCoeffPrev,
                                          &dmx_im[windowLen * window], 0);
              }
            }
          }
        }

        LONG pred_dir = cplxPredictionData->pred_dir ? -1 : 1;
        const FIXP_DBL pointOne = 0x66666666;
        const INT shift_dmx = SF_FNA_COEFFS + 1;

        for (band = 0; band < max_sfb_ste_outside; band++) {
          if (pJointStereoData->MsUsed[band] & groupMask) {
            FIXP_SGL tempRe =
                (FIXP_SGL)cplxPredictionData->alpha_q_re[group][band];
            FIXP_SGL tempIm =
                (FIXP_SGL)cplxPredictionData->alpha_q_im[group][band];

            int alpha_re_headroom = CountLeadingBits((INT)tempRe) - 16;
            if (tempRe == (FIXP_SGL)0) alpha_re_headroom = 15;
            int alpha_im_headroom = CountLeadingBits((INT)tempIm) - 16;
            if (tempIm == (FIXP_SGL)0) alpha_im_headroom = 15;
            int val = fMin(alpha_re_headroom, alpha_im_headroom);

            FDK_ASSERT(val >= 0);
            FIXP_DBL alpha_re_tmp = fMult((FIXP_SGL)(tempRe << val), pointOne);
            FIXP_DBL alpha_im_tmp = fMult((FIXP_SGL)(tempIm << val), pointOne);

            int alpha_re_exp = -3 + 15 - val;
            int help3_shift = alpha_re_exp + 1;

            FIXP_DBL *p2CoeffL = &(
                spectrumL[windowLen * window + pScaleFactorBandOffsets[band]]);
            FIXP_DBL *p2CoeffR = &(
                spectrumR[windowLen * window + pScaleFactorBandOffsets[band]]);
            FIXP_DBL *p2dmxIm =
                &(dmx_im[windowLen * window + pScaleFactorBandOffsets[band]]);
            FIXP_DBL *p2dmxRe =
                &(dmx_re[windowLen * window + pScaleFactorBandOffsets[band]]);

            for (int i = pScaleFactorBandOffsets[band];
                 i < pScaleFactorBandOffsets[band + 1]; i++) {
              FIXP_DBL side, left, right;

              side = fMultAddDiv2(fMultDiv2(alpha_re_tmp, *p2dmxRe++),
                                  alpha_im_tmp, (*p2dmxIm++) << shift_dmx);
              side = ((*p2CoeffR) >> 2) -
                     (FIXP_DBL)SATURATE_SHIFT(side, -(help3_shift - 2),
                                              DFRACT_BITS - 2);

              left = ((*p2CoeffL) >> 2) + side;
              right = ((*p2CoeffL) >> 2) - side;
              right = (FIXP_DBL)((LONG)right * pred_dir);

              *p2CoeffL++ = SATURATE_LEFT_SHIFT_ALT(left, 2, DFRACT_BITS);
              *p2CoeffR++ = SATURATE_LEFT_SHIFT_ALT(right, 2, DFRACT_BITS);
            }
          }
        }
      }
    }

    C_AALLOC_SCRATCH_END(dmx_im, FIXP_DBL, 1024);

  } else {

    for (window = 0, group = 0; group < windowGroups; group++) {
      groupMask = 1 << group;

      for (int groupwin = 0; groupwin < pWindowGroupLength[group];
           groupwin++, window++) {
        FIXP_DBL *leftSpectrum, *rightSpectrum;
        SHORT *leftScale = &SFBleftScale[window * 16];
        SHORT *rightScale = &SFBrightScale[window * 16];

        leftSpectrum =
            SPEC(spectrumL, window, pAacDecoderChannelInfo[L]->granuleLength);
        rightSpectrum =
            SPEC(spectrumR, window, pAacDecoderChannelInfo[R]->granuleLength);

        for (band = 0; band < max_sfb_ste_outside; band++) {
          if (pJointStereoData->MsUsed[band] & groupMask) {
            int lScale = leftScale[band];
            int rScale = rightScale[band];
            int commonScale = lScale > rScale ? lScale : rScale;
            unsigned int offsetCurrBand, offsetNextBand;

            FDK_ASSERT(GetWindowSequence(&pAacDecoderChannelInfo[L]->icsInfo) ==
                       GetWindowSequence(&pAacDecoderChannelInfo[R]->icsInfo));
            FDK_ASSERT(GetWindowShape(&pAacDecoderChannelInfo[L]->icsInfo) ==
                       GetWindowShape(&pAacDecoderChannelInfo[R]->icsInfo));

            commonScale++;
            leftScale[band] = commonScale;
            rightScale[band] = commonScale;

            lScale = fMin(DFRACT_BITS - 1, commonScale - lScale);
            rScale = fMin(DFRACT_BITS - 1, commonScale - rScale);

            FDK_ASSERT(lScale >= 0 && rScale >= 0);

            offsetCurrBand = pScaleFactorBandOffsets[band];
            offsetNextBand = pScaleFactorBandOffsets[band + 1];

            CJointStereo_GenerateMSOutput(&(leftSpectrum[offsetCurrBand]),
                                          &(rightSpectrum[offsetCurrBand]),
                                          lScale, rScale,
                                          offsetNextBand - offsetCurrBand);
          }
        }
        if (scaleFactorBandsTransmittedL > scaleFactorBandsTransmitted) {
          for (; band < scaleFactorBandsTransmittedL; band++) {
            if (pJointStereoData->MsUsed[band] & groupMask) {
              rightScale[band] = leftScale[band];

              for (int index = pScaleFactorBandOffsets[band];
                   index < pScaleFactorBandOffsets[band + 1]; index++) {
                FIXP_DBL leftCoefficient = leftSpectrum[index];
                rightSpectrum[index] = leftCoefficient;
              }
            }
          }
        } else if (scaleFactorBandsTransmittedR > scaleFactorBandsTransmitted) {
          for (; band < scaleFactorBandsTransmittedR; band++) {
            if (pJointStereoData->MsUsed[band] & groupMask) {
              leftScale[band] = rightScale[band];

              for (int index = pScaleFactorBandOffsets[band];
                   index < pScaleFactorBandOffsets[band + 1]; index++) {
                FIXP_DBL rightCoefficient = rightSpectrum[index];

                leftSpectrum[index] = rightCoefficient;
                rightSpectrum[index] = -rightCoefficient;
              }
            }
          }
        }
      }
    }

    if (pJointStereoData->MsMaskPresent == 2) {
      FDKmemclear(pJointStereoData->MsUsed,
                  JointStereoMaximumBands * sizeof(UCHAR));
    }
  }
}

void CJointStereo_ApplyIS(CAacDecoderChannelInfo *pAacDecoderChannelInfo[2],
                          const SHORT *pScaleFactorBandOffsets,
                          const UCHAR *pWindowGroupLength,
                          const int windowGroups,
                          const int scaleFactorBandsTransmitted) {
  CJointStereoData *pJointStereoData =
      &pAacDecoderChannelInfo[L]->pComData->jointStereoData;

  for (int window = 0, group = 0; group < windowGroups; group++) {
    UCHAR *CodeBook;
    SHORT *ScaleFactor;
    UCHAR groupMask = 1 << group;

    CodeBook = &pAacDecoderChannelInfo[R]->pDynData->aCodeBook[group * 16];
    ScaleFactor =
        &pAacDecoderChannelInfo[R]->pDynData->aScaleFactor[group * 16];

    for (int groupwin = 0; groupwin < pWindowGroupLength[group];
         groupwin++, window++) {
      FIXP_DBL *leftSpectrum, *rightSpectrum;
      SHORT *leftScale =
          &pAacDecoderChannelInfo[L]->pDynData->aSfbScale[window * 16];
      SHORT *rightScale =
          &pAacDecoderChannelInfo[R]->pDynData->aSfbScale[window * 16];
      int band;

      leftSpectrum = SPEC(pAacDecoderChannelInfo[L]->pSpectralCoefficient,
                          window, pAacDecoderChannelInfo[L]->granuleLength);
      rightSpectrum = SPEC(pAacDecoderChannelInfo[R]->pSpectralCoefficient,
                           window, pAacDecoderChannelInfo[R]->granuleLength);

      for (band = 0; band < scaleFactorBandsTransmitted; band++) {
        if ((CodeBook[band] == INTENSITY_HCB) ||
            (CodeBook[band] == INTENSITY_HCB2)) {
          int bandScale = -(ScaleFactor[band] + 100);

          int msb = bandScale >> 2;
          int lsb = bandScale & 0x03;

          FIXP_DBL scale = MantissaTable[lsb][0];

          FDK_ASSERT(GetWindowSequence(&pAacDecoderChannelInfo[L]->icsInfo) ==
                     GetWindowSequence(&pAacDecoderChannelInfo[R]->icsInfo));
          FDK_ASSERT(GetWindowShape(&pAacDecoderChannelInfo[L]->icsInfo) ==
                     GetWindowShape(&pAacDecoderChannelInfo[R]->icsInfo));

          rightScale[band] = leftScale[band] + msb + 1;

          if (pJointStereoData->MsUsed[band] & groupMask) {
            if (CodeBook[band] == INTENSITY_HCB)
            {
              scale = -scale;
            }
          } else {
            if (CodeBook[band] == INTENSITY_HCB2)
            {
              scale = -scale;
            }
          }

          for (int index = pScaleFactorBandOffsets[band];
               index < pScaleFactorBandOffsets[band + 1]; index++) {
            rightSpectrum[index] = fMult(leftSpectrum[index], scale);
          }
        }
      }
    }
  }
}
