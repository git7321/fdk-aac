/**************************** AAC decoder library ******************************

   Author(s):   Manuel Jander

   Description: USAC FAC

*******************************************************************************/

#include "usacdec_fac.h"

#include "usacdec_const.h"
#include "usacdec_lpc.h"
#include "usacdec_acelp.h"
#include "usacdec_rom.h"
#include "dct.h"
#include "FDK_tools_rom.h"
#include "mdct.h"

#define SPEC_FAC(ptr, i, gl) ((ptr) + ((i) * (gl)))

FIXP_DBL *CLpd_FAC_GetMemory(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                             UCHAR mod[NB_DIV], int *pState) {
  FIXP_DBL *ptr;
  int i;
  int k = 0;
  int max_windows = 8;

  FDK_ASSERT(*pState >= 0 && *pState < max_windows);

  for (i = *pState; i < max_windows; i++) {
    if (mod[i >> 1] == 0) {
      break;
    }
  }

  *pState = i + 1;

  if (i == max_windows) {
    ptr = pAacDecoderChannelInfo->data.usac.fac_data0;
  } else {
    FDK_ASSERT(mod[(i >> 1)] == 0);
    ptr = SPEC_FAC(pAacDecoderChannelInfo->pSpectralCoefficient, i,
                   pAacDecoderChannelInfo->granuleLength << k);
  }

  return ptr;
}

int CLpd_FAC_Read(HANDLE_FDK_BITSTREAM hBs, FIXP_DBL *pFac, SCHAR *pFacScale,
                  int length, int use_gain, int frame) {
  FIXP_DBL fac_gain;
  int fac_gain_e = 0;

  if (use_gain) {
    CLpd_DecodeGain(&fac_gain, &fac_gain_e, FDKreadBits(hBs, 7));
  }

  if (CLpc_DecodeAVQ(hBs, pFac, 1, 1, length) != 0) {
    return -1;
  }

  {
    int scale;

    scale = getScalefactor(pFac, length);
    scaleValues(pFac, length, scale);
    pFacScale[frame] = DFRACT_BITS - 1 - scale;
  }

  if (use_gain) {
    int i;

    pFacScale[frame] += fac_gain_e;

    for (i = 0; i < length; i++) {
      pFac[i] = fMult(pFac[i], fac_gain);
    }
  }
  return 0;
}

static void Syn_filt_zero(const FIXP_LPC a[], const INT a_exp, INT length,
                          FIXP_DBL x[]) {
  int i, j;
  FIXP_DBL L_tmp;

  for (i = 0; i < length; i++) {
    L_tmp = (FIXP_DBL)0;

    for (j = 0; j < fMin(i, M_LP_FILTER_ORDER); j++) {
      L_tmp -= fMultDiv2(a[j], x[i - (j + 1)]) >> (LP_FILTER_SCALE - 1);
    }

    L_tmp = scaleValue(L_tmp, a_exp + LP_FILTER_SCALE);
    x[i] = fAddSaturate(x[i], L_tmp);
  }
}

static const FIXP_DBL gainFac[4] = {0x40000000, 0x2d413ccd, 0x20000000,
                                    0x16a09e66};

void CFac_ApplyGains(FIXP_DBL fac_data[LFAC], const INT fac_length,
                     const FIXP_DBL tcx_gain, const FIXP_DBL alfd_gains[],
                     const INT mod) {
  FIXP_DBL facFactor;
  int i;

  FDK_ASSERT((fac_length == 128) || (fac_length == 96));

  facFactor = fMult(gainFac[mod], tcx_gain);
  for (i = 0; i < fac_length; i++) {
    fac_data[i] = fMult(fac_data[i], facFactor);
  }

  for (i = 0; i < fac_length / 4; i++) {
    int k;

    k = i >> (3 - mod);
    fac_data[i] = fMult(fac_data[i], alfd_gains[k])
                  << 1;
  }
}

static void CFac_CalcFacSignal(FIXP_DBL *pOut, FIXP_DBL *pFac,
                               const int fac_scale, const int fac_length,
                               const FIXP_LPC A[M_LP_FILTER_ORDER],
                               const INT A_exp, const int fAddZir,
                               const int isFdFac) {
  FIXP_LPC wA[M_LP_FILTER_ORDER];
  FIXP_DBL tf_gain = (FIXP_DBL)0;
  int wlength;
  int scale = fac_scale;

  imdct_gain(&tf_gain, &scale, isFdFac ? 0 : fac_length);

  dct_IV(pFac, fac_length, &scale);
  if (tf_gain != (FIXP_DBL)0) {
    int i;

    for (i = 0; i < fac_length; i++) {
      pFac[i] = fMult(tf_gain, pFac[i]);
    }
  }
  scaleValuesSaturate(pOut, pFac, fac_length,
                      scale);

  E_LPC_a_weight(wA, A, M_LP_FILTER_ORDER);

  FDKmemclear(pOut + fac_length, fac_length * sizeof(FIXP_DBL));
  wlength = 2 * fac_length;

  Syn_filt_zero(wA, A_exp, wlength, pOut);
}

INT CLpd_FAC_Mdct2Acelp(H_MDCT hMdct, FIXP_DBL *output, FIXP_DBL *pFac,
                        const int fac_scale, FIXP_LPC *A, INT A_exp,
                        INT nrOutSamples, const INT fac_length,
                        const INT isFdFac, UCHAR prevWindowShape) {
  FIXP_DBL *pOvl;
  FIXP_DBL *pOut0;
  const FIXP_WTP *pWindow;
  int i, fl, nrSamples = 0;

  FDK_ASSERT(fac_length <= 1024 / (4 * 2));

  fl = fac_length * 2;

  pWindow = FDKgetWindowSlope(fl, prevWindowShape);

  if (hMdct->prev_fr != fl) {
    int nl = 0;
    imdct_adapt_parameters(hMdct, &fl, &nl, fac_length, pWindow, nrOutSamples);
    FDK_ASSERT(nl == 0);
  }

  if (nrSamples < nrOutSamples) {
    pOut0 = output;
    nrSamples += hMdct->ov_offset;
    FDKmemcpy(pOut0, hMdct->overlap.time, hMdct->ov_offset * sizeof(pOut0[0]));
    hMdct->ov_offset = 0;
  }

  pOvl = hMdct->overlap.freq + hMdct->ov_size - 1;

  if (nrSamples >= nrOutSamples) {
    pOut0 = hMdct->overlap.time + hMdct->ov_offset;
    hMdct->ov_offset += hMdct->prev_nr + fl / 2;
  } else {
    pOut0 = output + nrSamples;
    nrSamples += hMdct->prev_nr + fl / 2;
  }
  if (hMdct->prevPrevAliasSymmetry == 0) {
    for (i = 0; i < hMdct->prev_nr; i++) {
      FIXP_DBL x = -(*pOvl--);
      *pOut0 = IMDCT_SCALE_DBL(x);
      pOut0++;
    }
  } else {
    for (i = 0; i < hMdct->prev_nr; i++) {
      FIXP_DBL x = (*pOvl--);
      *pOut0 = IMDCT_SCALE_DBL(x);
      pOut0++;
    }
  }
  hMdct->prev_nr = 0;

  {
    if (pFac != NULL) {
      CFac_CalcFacSignal(pOut0, pFac, fac_scale, fac_length, A, A_exp, 0,
                         isFdFac);
    } else {
      FDKmemclear(pOut0, fac_length * sizeof(FIXP_DBL));
    }
  }

  i = 0;

  if (hMdct->prevPrevAliasSymmetry == 0) {
    for (; i < fl / 2; i++) {
      FIXP_DBL x0;

      x0 = -fMult(*pOvl--, pWindow[i].v.re);

      *pOut0 = fAddSaturate(*pOut0, IMDCT_SCALE_DBL(x0));
      pOut0++;
    }
  } else {
    for (; i < fl / 2; i++) {
      FIXP_DBL x0;

      x0 = fMult(*pOvl--, pWindow[i].v.re);

      *pOut0 = fAddSaturate(*pOut0, IMDCT_SCALE_DBL(x0));
      pOut0++;
    }
  }
  if (hMdct->pFacZir !=
      0) {
    FIXP_DBL *pOut = pOut0 - fl / 2;
    for (i = 0; i < fl / 2; i++) {
      pOut[i] = fAddSaturate(pOut[i], IMDCT_SCALE_DBL(hMdct->pFacZir[i]));
    }
    hMdct->pFacZir = NULL;
  }

  hMdct->prev_fr = 0;
  hMdct->prev_nr = 0;
  hMdct->prev_tl = 0;
  hMdct->prevPrevAliasSymmetry = hMdct->prevAliasSymmetry;

  return nrSamples;
}

INT CLpd_FAC_Acelp2Mdct(H_MDCT hMdct, FIXP_DBL *output, FIXP_DBL *_pSpec,
                        const SHORT spec_scale[], const int nSpec,
                        FIXP_DBL *pFac, const int fac_scale,
                        const INT fac_length, INT noOutSamples, const INT tl,
                        const FIXP_WTP *wrs, const INT fr, FIXP_LPC A[16],
                        INT A_exp, CAcelpStaticMem *acelp_mem,
                        const FIXP_DBL gain, const int last_frame_lost,
                        const int isFdFac, const UCHAR last_lpd_mode,
                        const int k, int currAliasingSymmetry) {
  FIXP_DBL *pCurr, *pOvl, *pSpec;
  const FIXP_WTP *pWindow;
  const FIXP_WTB *FacWindowZir_conceal;
  UCHAR doFacZirConceal = 0;
  int doDeemph = 1;
  const FIXP_WTB *FacWindowZir, *FacWindowSynth;
  FIXP_DBL *pOut0 = output, *pOut1;
  int w, i, fl, nl, nr, f_len, nrSamples = 0, s = 0, scale, total_gain_e;
  FIXP_DBL *pF, *pFAC_and_FAC_ZIR = NULL;
  FIXP_DBL total_gain = gain;

  FDK_ASSERT(fac_length <= 1024 / (4 * 2));
  switch (fac_length) {
    case 128:
      pWindow = SineWindow256;
      FacWindowZir = FacWindowZir128;
      FacWindowSynth = FacWindowSynth128;
      break;
    case 64:
      pWindow = SineWindow128;
      FacWindowZir = FacWindowZir64;
      FacWindowSynth = FacWindowSynth64;
      break;
    case 32:
      pWindow = SineWindow64;
      FacWindowZir = FacWindowZir32;
      FacWindowSynth = FacWindowSynth32;
      break;
    case 96:
      pWindow = SineWindow192;
      FacWindowZir = FacWindowZir96;
      FacWindowSynth = FacWindowSynth96;
      break;
    case 48:
      pWindow = SineWindow96;
      FacWindowZir = FacWindowZir48;
      FacWindowSynth = FacWindowSynth48;
      break;
    default:
      FDK_ASSERT(0);
      return 0;
  }

  FacWindowZir_conceal = FacWindowSynth;
  fl = fac_length * 2;
  nl = (tl - fl) >> 1;
  nr = (tl - fr) >> 1;

  if (noOutSamples > nrSamples) {
    FDKmemcpy(pOut0, hMdct->overlap.time, hMdct->ov_offset * sizeof(pOut0[0]));
    nrSamples = hMdct->ov_offset;
    hMdct->ov_offset = 0;
  }

  if (nrSamples >= noOutSamples) {
    pOut1 = hMdct->overlap.time + hMdct->ov_offset;
    if (hMdct->ov_offset < fac_length) {
      pOut0 = output + nrSamples;
    } else {
      pOut0 = pOut1;
    }
    hMdct->ov_offset += fac_length + nl;
  } else {
    pOut1 = output + nrSamples;
    pOut0 = output + nrSamples;
  }

  {
    pFAC_and_FAC_ZIR = CLpd_ACELP_GetFreeExcMem(acelp_mem, 2 * fac_length);
    {
      const FIXP_DBL *pTmp1, *pTmp2;

      doFacZirConceal |= ((last_frame_lost != 0) && (k == 0));
      doDeemph &= (last_lpd_mode != 4);
      if (doFacZirConceal) {
        FDKmemclear(pFAC_and_FAC_ZIR, 2 * fac_length * sizeof(FIXP_DBL));
        FacWindowSynth = (FIXP_WTB *)pFAC_and_FAC_ZIR;
        FacWindowZir = FacWindowZir_conceal;
      } else {
        CFac_CalcFacSignal(pFAC_and_FAC_ZIR, pFac, fac_scale + s, fac_length, A,
                           A_exp, 1, isFdFac);
      }

      {
        {
          CLpd_Acelp_Zir(A, A_exp, acelp_mem, fac_length, pFac, doDeemph);

          pTmp1 = pOut0;
          pTmp2 = pFac;
        }

        for (i = 0, w = 0; i < fac_length; i++) {
          FIXP_DBL x;
          x = fMultDiv2(pTmp2[i], FacWindowZir[w]);
          x += fMultDiv2(pTmp1[-i - 1], FacWindowSynth[w]);
          pOut1[i] = fAddSaturate(x, pFAC_and_FAC_ZIR[i]);
          w++;
        }
      }

      if (doFacZirConceal) {
        scaleValues(pOut1, fac_length, 1);
      }
    }
  }

  if (nrSamples < noOutSamples) {
    nrSamples += fac_length + nl;
  }

  total_gain = gain;
  total_gain_e = 0;
  imdct_gain(&total_gain, &total_gain_e, tl);

  scale = total_gain_e;
  pSpec = _pSpec;

  if (currAliasingSymmetry == 0) {
    dct_IV(pSpec, tl, &scale);
  } else {
    FIXP_DBL _tmp[1024 + ALIGNMENT_DEFAULT / sizeof(FIXP_DBL)];
    FIXP_DBL *tmp = (FIXP_DBL *)ALIGN_PTR(_tmp);
    C_ALLOC_ALIGNED_REGISTER(tmp, sizeof(_tmp));
    dst_III(pSpec, tmp, tl, &scale);
    C_ALLOC_ALIGNED_UNREGISTER(tmp);
  }

  if (total_gain != (FIXP_DBL)0) {
    for (i = 0; i < tl; i++) {
      pSpec[i] = fMult(pSpec[i], total_gain);
    }
  }
  int loc_scale = fixmin_I(spec_scale[0] + scale, (INT)DFRACT_BITS - 1);
  scaleValuesSaturate(pSpec, tl, loc_scale);

  pOut1 += fl / 2 - 1;
  pCurr = pSpec + tl - fl / 2;

  for (i = 0; i < fl / 2; i++) {
    FIXP_DBL x1;

    x1 = fMult(*pCurr++, pWindow[i].v.re);
    FDK_ASSERT((pOut1 >= hMdct->overlap.time &&
                pOut1 < hMdct->overlap.time + hMdct->ov_size) ||
               (pOut1 >= output && pOut1 < output + 1024));
    *pOut1 = fAddSaturate(*pOut1, IMDCT_SCALE_DBL(-x1));
    pOut1--;
  }

  pOut1 += (fl / 2) + 1;

  pFAC_and_FAC_ZIR += fac_length;

  if (nl == 0) {
    hMdct->pFacZir = pFAC_and_FAC_ZIR;
  } else {
    FDK_ASSERT(nl >= fac_length);
    hMdct->pFacZir = NULL;
  }

  pF = pFAC_and_FAC_ZIR;
  f_len = fac_length;

  pCurr = pSpec + tl - fl / 2 - 1;
  for (i = 0; i < nl; i++) {
    FIXP_DBL x = -(*pCurr--);
    if (i < f_len) {
      x = fAddSaturate(x, *pF++);
    }

    FDK_ASSERT((pOut1 >= hMdct->overlap.time &&
                pOut1 < hMdct->overlap.time + hMdct->ov_size) ||
               (pOut1 >= output && pOut1 < output + 1024));
    *pOut1 = IMDCT_SCALE_DBL(x);
    pOut1++;
  }

  hMdct->prev_nr = nr;
  hMdct->prev_fr = fr;
  hMdct->prev_wrs = wrs;
  hMdct->prev_tl = tl;
  hMdct->prevPrevAliasSymmetry = hMdct->prevAliasSymmetry;
  hMdct->prevAliasSymmetry = currAliasingSymmetry;
  fl = fr;
  nl = nr;

  pOvl = pSpec + tl / 2 - 1;
  pOut0 = pOut1;

  for (w = 1; w < nSpec; w++)
  {
    const FIXP_WTP *pWindow_prev;

    pWindow_prev = hMdct->prev_wrs;

    pSpec = _pSpec + w * tl;

    scale = total_gain_e;

    if (currAliasingSymmetry == 0) {
      dct_IV(pSpec, tl, &scale);
    } else {
      dst_IV(pSpec, tl, &scale);
    }

    if (total_gain != (FIXP_DBL)0) {
      for (i = 0; i < tl; i++) {
        pSpec[i] = fMult(pSpec[i], total_gain);
      }
    }
    loc_scale = fixmin_I(spec_scale[w] + scale, (INT)DFRACT_BITS - 1);
    scaleValuesSaturate(pSpec, tl, loc_scale);

    if (noOutSamples <= nrSamples) {
      pOut0 = hMdct->overlap.time + hMdct->ov_offset;
      hMdct->ov_offset += hMdct->prev_nr + fl / 2;
    } else {
      nrSamples += hMdct->prev_nr + fl / 2;
    }

    for (i = 0; i < hMdct->prev_nr; i++) {
      FIXP_DBL x = -(*pOvl--);
      *pOut0 = IMDCT_SCALE_DBL(x);
      pOut0++;
    }

    if (noOutSamples <= nrSamples) {
      pOut1 = hMdct->overlap.time + hMdct->ov_offset + fl / 2 - 1;
      hMdct->ov_offset += fl / 2 + nl;
    } else {
      pOut1 = pOut0 + (fl - 1);
      nrSamples += fl / 2 + nl;
    }

    pCurr = pSpec + tl - fl / 2;
    if (currAliasingSymmetry == 0) {
      for (i = 0; i < fl / 2; i++) {
        FIXP_DBL x0, x1;

        cplxMultDiv2(&x1, &x0, *pCurr++, -*pOvl--, pWindow_prev[i]);
        *pOut0 = IMDCT_SCALE_DBL_LSH1(x0);
        *pOut1 = IMDCT_SCALE_DBL_LSH1(-x1);
        pOut0++;
        pOut1--;
      }
    } else {
      if (hMdct->prevPrevAliasSymmetry == 0) {
        for (i = 0; i < fl / 2; i++) {
          FIXP_DBL x0, x1;

          cplxMultDiv2(&x1, &x0, *pCurr++, -*pOvl--, pWindow_prev[i]);
          *pOut0 = IMDCT_SCALE_DBL_LSH1(x0);
          *pOut1 = IMDCT_SCALE_DBL_LSH1(x1);
          pOut0++;
          pOut1--;
        }
      } else {
        for (i = 0; i < fl / 2; i++) {
          FIXP_DBL x0, x1;

          cplxMultDiv2(&x1, &x0, *pCurr++, *pOvl--, pWindow_prev[i]);
          *pOut0 = IMDCT_SCALE_DBL_LSH1(x0);
          *pOut1 = IMDCT_SCALE_DBL_LSH1(x1);
          pOut0++;
          pOut1--;
        }
      }
    }

    if (hMdct->pFacZir != 0) {
      FIXP_DBL *pOut = pOut0 - fl / 2;
      FDK_ASSERT(fl / 2 <= 128);
      for (i = 0; i < fl / 2; i++) {
        pOut[i] = fAddSaturate(pOut[i], IMDCT_SCALE_DBL(hMdct->pFacZir[i]));
      }
      hMdct->pFacZir = NULL;
    }
    pOut0 += (fl / 2);

    pOut1 += (fl / 2) + 1;
    pCurr = pSpec + tl - fl / 2 - 1;
    for (i = 0; i < nl; i++) {
      FIXP_DBL x = -(*pCurr--);
      *pOut1 = IMDCT_SCALE_DBL(x);
      pOut1++;
    }

    pOvl = pSpec + tl / 2 - 1;

    hMdct->prev_nr = nr;
    hMdct->prev_fr = fr;
    hMdct->prev_tl = tl;
    hMdct->prev_wrs = pWindow_prev;
    hMdct->prevPrevAliasSymmetry = hMdct->prevAliasSymmetry;
    hMdct->prevAliasSymmetry = currAliasingSymmetry;
  }

  pOvl = hMdct->overlap.freq + hMdct->ov_size - tl / 2;
  FDK_ASSERT(pOvl >= hMdct->overlap.time + hMdct->ov_offset);
  FDK_ASSERT(tl / 2 <= hMdct->ov_size);
  for (i = 0; i < tl / 2; i++) {
    pOvl[i] = _pSpec[i + (w - 1) * tl];
  }

  return nrSamples;
}
