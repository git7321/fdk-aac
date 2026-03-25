/******************* Library for basic calculation routines ********************

   Author(s):   Josef Hoepfl, Manuel Jander, Youliy Ninov, Daniel Hagel

   Description: MDCT/MDST routines

*******************************************************************************/

#include "mdct.h"

#include "FDK_tools_rom.h"
#include "dct.h"
#include "fixpoint_math.h"

void mdct_init(H_MDCT hMdct, FIXP_DBL *overlap, INT overlapBufferSize) {
  hMdct->overlap.freq = overlap;
  hMdct->prev_fr = 0;
  hMdct->prev_nr = 0;
  hMdct->prev_tl = 0;
  hMdct->ov_size = overlapBufferSize;
  hMdct->prevAliasSymmetry = 0;
  hMdct->prevPrevAliasSymmetry = 0;
  hMdct->pFacZir = NULL;
  hMdct->pAsymOvlp = NULL;
}

INT mdct_block(H_MDCT hMdct, const INT_PCM *RESTRICT timeData,
               const INT noInSamples, FIXP_DBL *RESTRICT mdctData,
               const INT nSpec, const INT tl, const FIXP_WTP *pRightWindowPart,
               const INT fr, SHORT *pMdctData_e) {
  int i, n;
  int fl, nl, nr;
  const FIXP_WTP *wls, *wrs;

  wrs = pRightWindowPart;

  if (hMdct->prev_fr ==
      0) {
    hMdct->prev_fr = fr;
    hMdct->prev_wrs = wrs;
    hMdct->prev_tl = tl;
  }

  nr = (tl - fr) >> 1;
  timeData += (noInSamples - tl) >> 1;

  for (n = 0; n < nSpec; n++) {
    INT mdctData_e = 1 + 1;

    wls = hMdct->prev_wrs;
    fl = hMdct->prev_fr;
    nl = (tl - fl) >> 1;

    for (i = 0; i < nl; i++) {
#if SAMPLE_BITS == DFRACT_BITS
      mdctData[(tl / 2) + i] = -((FIXP_DBL)timeData[tl - i - 1] >> (1));
#else
      mdctData[(tl / 2) + i] = -(FIXP_DBL)timeData[tl - i - 1]
                               << (DFRACT_BITS - SAMPLE_BITS - 1);
#endif
    }

    for (i = 0; i < fl / 2; i++) {
      FIXP_DBL tmp0;
      tmp0 = fMultDiv2((FIXP_PCM)timeData[i + nl], wls[i].v.im);
      mdctData[(tl / 2) + i + nl] =
          fMultSubDiv2(tmp0, (FIXP_PCM)timeData[tl - nl - i - 1],
                       wls[i].v.re);
    }

    for (i = 0; i < nr; i++) {
#if SAMPLE_BITS == \
    DFRACT_BITS
      mdctData[(tl / 2) - 1 - i] = -((FIXP_DBL)timeData[tl + i] >> (1));
#else
      mdctData[(tl / 2) - 1 - i] =
          -(FIXP_DBL)timeData[tl + i]
          << (DFRACT_BITS - SAMPLE_BITS - 1);
#endif
    }

    for (i = 0; i < fr / 2; i++) {
      FIXP_DBL tmp1;
      tmp1 = fMultDiv2((FIXP_PCM)timeData[tl + nr + i],
                       wrs[i].v.re);
      mdctData[(tl / 2) - nr - i - 1] =
          -fMultAddDiv2(tmp1, (FIXP_PCM)timeData[(tl * 2) - nr - i - 1],
                        wrs[i].v.im);
    }

    dct_IV(mdctData, tl, &mdctData_e);

    pMdctData_e[n] = (SHORT)mdctData_e;

    timeData += tl;
    mdctData += tl;

    hMdct->prev_wrs = wrs;
    hMdct->prev_fr = fr;
    hMdct->prev_tl = tl;
  }

  return nSpec * tl;
}

void imdct_gain(FIXP_DBL *pGain_m, int *pGain_e, int tl) {
  FIXP_DBL gain_m = *pGain_m;
  int gain_e = *pGain_e;
  int log2_tl;

  gain_e += -MDCT_OUTPUT_GAIN - MDCT_OUT_HEADROOM + 1;
  if (tl == 0) {
    *pGain_e = gain_e;
    return;
  }

  log2_tl = DFRACT_BITS - 1 - fNormz((FIXP_DBL)tl);
  gain_e += -log2_tl;

  switch ((tl) >> (log2_tl - 2)) {
    case 0x7:
      if (gain_m == (FIXP_DBL)0) {
        gain_m = FL2FXCONST_DBL(0.53333333333333333333f);
      } else {
        gain_m = fMult(gain_m, FL2FXCONST_DBL(0.53333333333333333333f));
      }
      break;
    case 0x6:
      if (gain_m == (FIXP_DBL)0) {
        gain_m = FL2FXCONST_DBL(2.0 / 3.0f);
      } else {
        gain_m = fMult(gain_m, FL2FXCONST_DBL(2.0 / 3.0f));
      }
      break;
    case 0x5:
      if (gain_m == (FIXP_DBL)0) {
        gain_m = FL2FXCONST_DBL(0.53333333333333333333f);
      } else {
        gain_m = fMult(gain_m, FL2FXCONST_DBL(0.53333333333333333333f));
      }
      break;
    case 0x4:
      break;
    default:
      FDK_ASSERT(0);
      break;
  }

  *pGain_m = gain_m;
  *pGain_e = gain_e;
}

INT imdct_drain(H_MDCT hMdct, FIXP_DBL *output, INT nrSamplesRoom) {
  int buffered_samples = 0;

  if (nrSamplesRoom > 0) {
    buffered_samples = hMdct->ov_offset;

    FDK_ASSERT(buffered_samples <= nrSamplesRoom);

    if (buffered_samples > 0) {
      FDKmemcpy(output, hMdct->overlap.time,
                buffered_samples * sizeof(FIXP_DBL));
      hMdct->ov_offset = 0;
    }
  }
  return buffered_samples;
}

INT imdct_copy_ov_and_nr(H_MDCT hMdct, FIXP_DBL *pTimeData, INT nrSamples) {
  FIXP_DBL *pOvl;
  int nt, nf, i;

  nt = fMin(hMdct->ov_offset, nrSamples);
  nrSamples -= nt;
  nf = fMin(hMdct->prev_nr, nrSamples);
  FDKmemcpy(pTimeData, hMdct->overlap.time, nt * sizeof(FIXP_DBL));
  pTimeData += nt;

  pOvl = hMdct->overlap.freq + hMdct->ov_size - 1;
  if (hMdct->prevPrevAliasSymmetry == 0) {
    for (i = 0; i < nf; i++) {
      FIXP_DBL x = -(*pOvl--);
      *pTimeData = IMDCT_SCALE_DBL(x);
      pTimeData++;
    }
  } else {
    for (i = 0; i < nf; i++) {
      FIXP_DBL x = (*pOvl--);
      *pTimeData = IMDCT_SCALE_DBL(x);
      pTimeData++;
    }
  }

  return (nt + nf);
}

void imdct_adapt_parameters(H_MDCT hMdct, int *pfl, int *pnl, int tl,
                            const FIXP_WTP *wls, int noOutSamples) {
  int fl = *pfl, nl = *pnl;
  int window_diff, use_current = 0, use_previous = 0;
  if (hMdct->prev_tl == 0) {
    hMdct->prev_wrs = wls;
    hMdct->prev_fr = fl;
    hMdct->prev_nr = (noOutSamples - fl) >> 1;
    hMdct->prev_tl = noOutSamples;
    hMdct->ov_offset = 0;
    use_current = 1;
  }

  window_diff = (hMdct->prev_fr - fl) >> 1;

  if (hMdct->prev_nr + window_diff > 0) {
    use_current = 1;
  }
  if (nl - window_diff > 0) {
    use_previous = 1;
  }

  if (use_current && use_previous) {
    if (fl < hMdct->prev_fr) {
      use_current = 0;
    }
  }
  if (use_current) {
    hMdct->prev_nr += window_diff;
    hMdct->prev_fr = fl;
    hMdct->prev_wrs = wls;
  } else {
    nl -= window_diff;
    fl = hMdct->prev_fr;
  }

  *pfl = fl;
  *pnl = nl;
}

INT imlt_block(H_MDCT hMdct, FIXP_DBL *output, FIXP_DBL *spectrum,
               const SHORT scalefactor[], const INT nSpec,
               const INT noOutSamples, const INT tl, const FIXP_WTP *wls,
               INT fl, const FIXP_WTP *wrs, const INT fr, FIXP_DBL gain,
               int flags) {
  FIXP_DBL *pOvl;
  FIXP_DBL *pOut0 = output, *pOut1;
  INT nl, nr;
  int w, i, nrSamples = 0, specShiftScale, transform_gain_e = 0;
  int currAliasSymmetry = (flags & MLT_FLAG_CURR_ALIAS_SYMMETRY);

  nr = (tl - fr) >> 1;
  nl = (tl - fl) >> 1;

  imdct_gain(&gain, &transform_gain_e, tl);

  if (hMdct->prev_fr != fl) {
    imdct_adapt_parameters(hMdct, &fl, &nl, tl, wls, noOutSamples);
  }

  pOvl = hMdct->overlap.freq + hMdct->ov_size - 1;

  if (noOutSamples > nrSamples) {
    for (i = 0; i < hMdct->ov_offset; i++) {
      *pOut0 = hMdct->overlap.time[i];
      pOut0++;
    }
    nrSamples = hMdct->ov_offset;
    hMdct->ov_offset = 0;
  }

  for (w = 0; w < nSpec; w++) {
    FIXP_DBL *pSpec, *pCurr;
    const FIXP_WTP *pWindow;

    if (hMdct->prev_fr != fl) {
      imdct_adapt_parameters(hMdct, &fl, &nl, tl, wls, noOutSamples);
    }

    specShiftScale = transform_gain_e;

    pWindow = hMdct->prev_wrs;

    pSpec = spectrum + w * tl;

    if (currAliasSymmetry == 0) {
      if (hMdct->prevAliasSymmetry == 0) {
        dct_IV(pSpec, tl, &specShiftScale);
      } else {
        FIXP_DBL _tmp[1024 + ALIGNMENT_DEFAULT / sizeof(FIXP_DBL)];
        FIXP_DBL *tmp = (FIXP_DBL *)ALIGN_PTR(_tmp);
        C_ALLOC_ALIGNED_REGISTER(tmp, sizeof(_tmp));
        dct_III(pSpec, tmp, tl, &specShiftScale);
        C_ALLOC_ALIGNED_UNREGISTER(tmp);
      }
    } else {
      if (hMdct->prevAliasSymmetry == 0) {
        FIXP_DBL _tmp[1024 + ALIGNMENT_DEFAULT / sizeof(FIXP_DBL)];
        FIXP_DBL *tmp = (FIXP_DBL *)ALIGN_PTR(_tmp);
        C_ALLOC_ALIGNED_REGISTER(tmp, sizeof(_tmp));
        dst_III(pSpec, tmp, tl, &specShiftScale);
        C_ALLOC_ALIGNED_UNREGISTER(tmp);
      } else {
        dst_IV(pSpec, tl, &specShiftScale);
      }
    }

    if (gain != (FIXP_DBL)0) {
      for (i = 0; i < tl; i++) {
        pSpec[i] = fMult(pSpec[i], gain);
      }
    }

    {
      int loc_scale =
          fixmin_I(scalefactor[w] + specShiftScale, (INT)DFRACT_BITS - 1);
      DWORD_ALIGNED(pSpec);
      scaleValuesSaturate(pSpec, tl, loc_scale);
    }

    if (noOutSamples <= nrSamples) {
      pOut0 = hMdct->overlap.time + hMdct->ov_offset;
      hMdct->ov_offset += hMdct->prev_nr + fl / 2;
    } else {
      nrSamples += hMdct->prev_nr + fl / 2;
    }

    if ((hMdct->pFacZir != 0) && (hMdct->prev_nr == fl / 2)) {
      for (i = 0; i < hMdct->prev_nr; i++) {
        FIXP_DBL x = -(*pOvl--);
        *pOut0 = fAddSaturate(x, IMDCT_SCALE_DBL(hMdct->pFacZir[i]));
        pOut0++;
      }
      hMdct->pFacZir = NULL;
    } else {
      if (hMdct->prevPrevAliasSymmetry == 0) {
        for (i = 0; i < hMdct->prev_nr; i++) {
          FIXP_DBL x = -(*pOvl--);
          *pOut0 = IMDCT_SCALE_DBL(x);
          pOut0++;
        }
      } else {
        for (i = 0; i < hMdct->prev_nr; i++) {
          FIXP_DBL x = *pOvl--;
          *pOut0 = IMDCT_SCALE_DBL(x);
          pOut0++;
        }
      }
    }

    if (noOutSamples <= nrSamples) {
      pOut1 = hMdct->overlap.time + hMdct->ov_offset + fl / 2 - 1;
      hMdct->ov_offset += fl / 2 + nl;
    } else {
      pOut1 = pOut0 + (fl - 1);
      nrSamples += fl / 2 + nl;
    }

    pCurr = pSpec + tl - fl / 2;
    DWORD_ALIGNED(pCurr);
    C_ALLOC_ALIGNED_REGISTER(pWindow, fl);
    DWORD_ALIGNED(pWindow);
    C_ALLOC_ALIGNED_UNREGISTER(pWindow);

    if (hMdct->prevPrevAliasSymmetry == 0) {
      if (hMdct->prevAliasSymmetry == 0) {
        if (!hMdct->pAsymOvlp) {
          for (i = 0; i < fl / 2; i++) {
            FIXP_DBL x0, x1;
            cplxMultDiv2(&x1, &x0, *pCurr++, -*pOvl--, pWindow[i]);
            *pOut0 = IMDCT_SCALE_DBL_LSH1(x0);
            *pOut1 = IMDCT_SCALE_DBL_LSH1(-x1);
            pOut0++;
            pOut1--;
          }
        } else {
          FIXP_DBL *pAsymOvl = hMdct->pAsymOvlp + fl / 2 - 1;
          for (i = 0; i < fl / 2; i++) {
            FIXP_DBL x0, x1;
            x1 = -fMultDiv2(*pCurr, pWindow[i].v.re) +
                 fMultDiv2(*pAsymOvl, pWindow[i].v.im);
            x0 = fMultDiv2(*pCurr, pWindow[i].v.im) -
                 fMultDiv2(*pOvl, pWindow[i].v.re);
            pCurr++;
            pOvl--;
            pAsymOvl--;
            *pOut0++ = IMDCT_SCALE_DBL_LSH1(x0);
            *pOut1-- = IMDCT_SCALE_DBL_LSH1(x1);
          }
          hMdct->pAsymOvlp = NULL;
        }
      } else {
        for (i = 0; i < fl / 2; i++) {
          FIXP_DBL x0, x1;
          cplxMultDiv2(&x1, &x0, *pCurr++, -*pOvl--, pWindow[i]);
          *pOut0 = IMDCT_SCALE_DBL_LSH1(x0);
          *pOut1 = IMDCT_SCALE_DBL_LSH1(x1);
          pOut0++;
          pOut1--;
        }
      }
    } else {
      if (hMdct->prevAliasSymmetry == 0) {
        for (i = 0; i < fl / 2; i++) {
          FIXP_DBL x0, x1;
          cplxMultDiv2(&x1, &x0, *pCurr++, *pOvl--, pWindow[i]);
          *pOut0 = IMDCT_SCALE_DBL_LSH1(x0);
          *pOut1 = IMDCT_SCALE_DBL_LSH1(-x1);
          pOut0++;
          pOut1--;
        }
      } else {
        for (i = 0; i < fl / 2; i++) {
          FIXP_DBL x0, x1;
          cplxMultDiv2(&x1, &x0, *pCurr++, *pOvl--, pWindow[i]);
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
    pOut0 += (fl / 2) + nl;

    pOut1 += (fl / 2) + 1;
    pCurr = pSpec + tl - fl / 2 - 1;
    if (hMdct->prevAliasSymmetry == 0) {
      for (i = 0; i < nl; i++) {
        FIXP_DBL x = -(*pCurr--);
        *pOut1++ = IMDCT_SCALE_DBL(x);
      }
    } else {
      for (i = 0; i < nl; i++) {
        FIXP_DBL x = *pCurr--;
        *pOut1++ = IMDCT_SCALE_DBL(x);
      }
    }

    pOvl = pSpec + tl / 2 - 1;

    hMdct->prev_nr = nr;
    hMdct->prev_fr = fr;
    hMdct->prev_tl = tl;
    hMdct->prev_wrs = wrs;

    hMdct->prevPrevAliasSymmetry = hMdct->prevAliasSymmetry;
    hMdct->prevAliasSymmetry = currAliasSymmetry;
  }

  pOvl = hMdct->overlap.freq + hMdct->ov_size - tl / 2;
  FDKmemcpy(pOvl, &spectrum[(nSpec - 1) * tl], (tl / 2) * sizeof(FIXP_DBL));

  return nrSamples;
}
