/*********************** MPEG surround encoder library *************************

   Author(s):   Max Neuendorf

   Description: Encoder Library Interface
                Get windows for framing

*******************************************************************************/

#include "sacenc_framewindowing.h"
#include "sacenc_vectorfunctions.h"

typedef struct T_FRAMEWINDOW {
  INT nTimeSlotsMax;
  INT bFrameKeep;
  INT startSlope;
  INT stopSlope;
  INT startRect;
  INT stopRect;

  INT taperAnaLen;
  INT taperSynLen;
  FIXP_WIN pTaperAna__FDK[MAX_TIME_SLOTS];
  FIXP_WIN pTaperSyn__FDK[MAX_TIME_SLOTS];

} FRAMEWINDOW;

typedef enum {
  FIX_INVALID = -1,
  FIX_RECT_SMOOTH = 0,
  FIX_SMOOTH_RECT = 1,
  FIX_LARGE_SMOOTH = 2,
  FIX_RECT_TRIANG = 3

} FIX_TYPE;

typedef enum {
  VAR_INVALID = -1,
  VAR_HOLD = 0,
  VAR_ISOLATE = 1

} VAR_TYPE;

static void calcTaperWin(FIXP_WIN *pTaperWin, INT timeSlots) {
  FIXP_DBL x;
  int i, scale;

  for (i = 0; i < timeSlots; i++) {
    x = fDivNormHighPrec((FIXP_DBL)i, (FIXP_DBL)timeSlots, &scale);

    if (scale < 0) {
      pTaperWin[i] = FX_DBL2FX_WIN(x >> (-scale));
    } else {
      pTaperWin[i] = FX_DBL2FX_WIN(x << (scale));
    }
  }
  pTaperWin[timeSlots] = FX_DBL2FX_WIN((FIXP_DBL)MAXVAL_DBL);
}

FDK_SACENC_ERROR fdk_sacenc_frameWindow_Create(
    HANDLE_FRAMEWINDOW *phFrameWindow) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == phFrameWindow) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDK_ALLOCATE_MEMORY_1D(*phFrameWindow, 1, FRAMEWINDOW);
  }
  return error;

bail:
  fdk_sacenc_frameWindow_Destroy(phFrameWindow);
  return ((SACENC_OK == error) ? SACENC_MEMORY_ERROR : error);
}

FDK_SACENC_ERROR fdk_sacenc_frameWindow_Init(
    HANDLE_FRAMEWINDOW hFrameWindow,
    const FRAMEWINDOW_CONFIG *const pFrameWindowConfig) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((hFrameWindow == NULL) || (pFrameWindowConfig == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else if (pFrameWindowConfig->nTimeSlotsMax < 0) {
    error = SACENC_INIT_ERROR;
  } else {
    int ts;
    hFrameWindow->bFrameKeep = pFrameWindowConfig->bFrameKeep;
    hFrameWindow->nTimeSlotsMax = pFrameWindowConfig->nTimeSlotsMax;

    FIXP_WIN winMaxVal = FX_DBL2FX_WIN((FIXP_DBL)MAXVAL_DBL);
    int timeSlots = pFrameWindowConfig->nTimeSlotsMax;
    {
      hFrameWindow->startSlope = 0;
      hFrameWindow->stopSlope = ((3 * timeSlots) >> 1) - 1;
      hFrameWindow->startRect = timeSlots >> 1;
      hFrameWindow->stopRect = timeSlots;
      calcTaperWin(hFrameWindow->pTaperSyn__FDK, timeSlots >> 1);
      hFrameWindow->taperSynLen = timeSlots >> 1;
    }

    hFrameWindow->taperAnaLen =
        hFrameWindow->startRect - hFrameWindow->startSlope;
    for (ts = 0; ts < hFrameWindow->taperAnaLen; ts++) {
      { hFrameWindow->pTaperAna__FDK[ts] = winMaxVal; }
    }
  }

  return error;
}

FDK_SACENC_ERROR fdk_sacenc_frameWindow_Destroy(
    HANDLE_FRAMEWINDOW *phFrameWindow) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((NULL != phFrameWindow) && (NULL != *phFrameWindow)) {
    FDKfree(*phFrameWindow);
    *phFrameWindow = NULL;
  }
  return error;
}

static FDK_SACENC_ERROR FrameWinList_Reset(FRAMEWIN_LIST *const pFrameWinList) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == pFrameWinList) {
    error = SACENC_INVALID_HANDLE;
  } else {
    int k = 0;
    for (k = 0; k < MAX_NUM_PARAMS; k++) {
      pFrameWinList->dat[k].slot = -1;
      pFrameWinList->dat[k].hold = FW_INTP;
    }
    pFrameWinList->n = 0;
  }
  return error;
}

static FDK_SACENC_ERROR FrameWindowList_Add(FRAMEWIN_LIST *const pFrameWinList,
                                            const INT slot,
                                            const FW_SLOTTYPE hold) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == pFrameWinList) {
    error = SACENC_INVALID_HANDLE;
  } else {
    if (pFrameWinList->n >= MAX_NUM_PARAMS) {
      error = SACENC_PARAM_ERROR;
    } else if (pFrameWinList->n > 0 &&
               pFrameWinList->dat[pFrameWinList->n - 1].slot - slot > 0) {
      error = SACENC_PARAM_ERROR;
    } else {
      pFrameWinList->dat[pFrameWinList->n].slot = slot;
      pFrameWinList->dat[pFrameWinList->n].hold = hold;
      pFrameWinList->n++;
    }
  }
  return error;
}

static FDK_SACENC_ERROR FrameWindowList_Remove(
    FRAMEWIN_LIST *const pFrameWinList, const INT idx) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == pFrameWinList) {
    error = SACENC_INVALID_HANDLE;
  } else {
    int k = 0;
    if (idx < 0 || idx >= MAX_NUM_PARAMS) {
      error = SACENC_PARAM_ERROR;
    } else if (pFrameWinList->n > 0) {
      if (idx == MAX_NUM_PARAMS - 1) {
        pFrameWinList->dat[idx].slot = -1;
        pFrameWinList->dat[idx].hold = FW_INTP;
      } else {
        for (k = idx; k < MAX_NUM_PARAMS - 1; k++) {
          pFrameWinList->dat[k] = pFrameWinList->dat[k + 1];
        }
      }
      pFrameWinList->n--;
    }
  }
  return error;
}

static FDK_SACENC_ERROR FrameWindowList_Limit(
    FRAMEWIN_LIST *const pFrameWinList, const INT ll,
    const INT ul
) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == pFrameWinList) {
    error = SACENC_INVALID_HANDLE;
  } else {
    int k = 0;
    for (k = 0; k < pFrameWinList->n; k++) {
      if (pFrameWinList->dat[k].slot < ll || pFrameWinList->dat[k].slot > ul) {
        FrameWindowList_Remove(pFrameWinList, k);
        --k;
      }
    }
  }
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_frameWindow_GetWindow(
    HANDLE_FRAMEWINDOW hFrameWindow, INT tr_pos[MAX_NUM_PARAMS],
    const INT timeSlots, FRAMINGINFO *const pFramingInfo,
    FIXP_WIN *pWindowAna__FDK[MAX_NUM_PARAMS],
    FRAMEWIN_LIST *const pFrameWinList, const INT avoid_keep) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((hFrameWindow == NULL) || (tr_pos == NULL) || (pFramingInfo == NULL) ||
      (pFrameWinList == NULL) || (pWindowAna__FDK == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    const VAR_TYPE varType = VAR_HOLD;
    const int tranL = 4;
    int winCnt = 0;
    int w, ps;

    int startSlope = hFrameWindow->startSlope;
    int stopSlope = hFrameWindow->stopSlope;
    int startRect = hFrameWindow->startRect;
    int stopRect = hFrameWindow->stopRect;
    int taperAnaLen = hFrameWindow->taperAnaLen;

    FIXP_WIN winMaxVal = FX_DBL2FX_WIN((FIXP_DBL)MAXVAL_DBL);
    FIXP_WIN applyRightWindowGain__FDK[MAX_NUM_PARAMS];
    FIXP_WIN *pTaperAna__FDK = hFrameWindow->pTaperAna__FDK;

    for (ps = 0; ps < MAX_NUM_PARAMS; ps++) {
      if (pWindowAna__FDK[ps] == NULL) {
        error = SACENC_INVALID_HANDLE;
        goto bail;
      }
    }

    if ((timeSlots > hFrameWindow->nTimeSlotsMax) || (timeSlots < 0)) {
      error = SACENC_INVALID_CONFIG;
      goto bail;
    }

    if (SACENC_OK != (error = FrameWinList_Reset(pFrameWinList))) goto bail;

    FDKmemclear(applyRightWindowGain__FDK, sizeof(applyRightWindowGain__FDK));

    if (tr_pos[0] > -1) {
      int p_l = tr_pos[0];
      winCnt = 0;

      switch (varType) {
        case VAR_HOLD:
          if (SACENC_OK !=
              (error = FrameWindowList_Add(pFrameWinList, p_l - 1, FW_HOLD)))
            goto bail;
          if (SACENC_OK !=
              (error = FrameWindowList_Add(pFrameWinList, p_l, FW_INTP)))
            goto bail;
          break;
        case VAR_ISOLATE:
          if (SACENC_OK !=
              (error = FrameWindowList_Add(pFrameWinList, p_l - 1, FW_HOLD)))
            goto bail;
          if (SACENC_OK !=
              (error = FrameWindowList_Add(pFrameWinList, p_l, FW_INTP)))
            goto bail;
          if (SACENC_OK != (error = FrameWindowList_Add(pFrameWinList,
                                                        p_l + tranL, FW_HOLD)))
            goto bail;
          if (SACENC_OK != (error = FrameWindowList_Add(
                                pFrameWinList, p_l + tranL + 1, FW_INTP)))
            goto bail;
          break;
        default:
          error = SACENC_INVALID_CONFIG;
          break;
      }

      if (SACENC_OK !=
          (error = FrameWindowList_Limit(pFrameWinList, 0, timeSlots - 1)))
        goto bail;

      if (SACENC_OK !=
          (error = FrameWindowList_Add(pFrameWinList, timeSlots - 1, FW_HOLD)))
        goto bail;

      for (ps = 0; ps < pFrameWinList->n - 1; ps++) {
        if (FW_HOLD != pFrameWinList->dat[ps].hold) {
          int const start = pFrameWinList->dat[ps].slot;
          int const stop = pFrameWinList->dat[ps + 1].slot;

          FDKmemset_flex(pWindowAna__FDK[winCnt], FX_DBL2FX_WIN((FIXP_DBL)0),
                         start);
          FDKmemset_flex(&pWindowAna__FDK[winCnt][start], winMaxVal,
                         stop - start + 1);
          FDKmemset_flex(&pWindowAna__FDK[winCnt][stop + 1],
                         FX_DBL2FX_WIN((FIXP_DBL)0), timeSlots - stop - 1);

          applyRightWindowGain__FDK[winCnt] =
              pWindowAna__FDK[winCnt][timeSlots - 1];
          winCnt++;
        }
      }

      if (SACENC_OK !=
          (error = FrameWindowList_Remove(pFrameWinList, pFrameWinList->n - 1)))
        goto bail;
    } else {
      winCnt = 0;

      if (SACENC_OK !=
          (error = FrameWindowList_Add(pFrameWinList, timeSlots - 1, FW_INTP)))
        goto bail;
      FDKmemset_flex(pWindowAna__FDK[winCnt], FX_DBL2FX_WIN((FIXP_DBL)0),
                     startSlope);
      FDKmemcpy_flex(&pWindowAna__FDK[winCnt][startSlope], 1, pTaperAna__FDK, 1,
                     taperAnaLen);
      FDKmemset_flex(&pWindowAna__FDK[winCnt][startRect], winMaxVal,
                     timeSlots - startRect);

      applyRightWindowGain__FDK[winCnt] = winMaxVal;
      winCnt++;
    }

    for (w = 0; w < winCnt; w++) {
      if (applyRightWindowGain__FDK[w] > (FIXP_WIN)0) {
        if (tr_pos[1] > -1) {
          int p_r = tr_pos[1];

          FDKmemset_flex(&pWindowAna__FDK[w][timeSlots], winMaxVal,
                         p_r - timeSlots);
          FDKmemset_flex(&pWindowAna__FDK[w][p_r], FX_DBL2FX_WIN((FIXP_DBL)0),
                         2 * timeSlots - p_r);

        } else {
          FDKmemset_flex(&pWindowAna__FDK[w][timeSlots], winMaxVal,
                         stopRect - timeSlots + 1);
          FDKmemcpy_flex(&pWindowAna__FDK[w][stopRect], 1,
                         &pTaperAna__FDK[taperAnaLen - 1], -1, taperAnaLen);
          FDKmemset_flex(&pWindowAna__FDK[w][stopSlope + 1],
                         FX_DBL2FX_WIN((FIXP_DBL)0),
                         2 * timeSlots - stopSlope - 1);
        }

        if (applyRightWindowGain__FDK[w] < winMaxVal) {
          int ts;
          for (ts = 0; ts < timeSlots; ts++) {
            pWindowAna__FDK[w][timeSlots + ts] =
                FX_DBL2FX_WIN(fMult(pWindowAna__FDK[w][timeSlots + ts],
                                    applyRightWindowGain__FDK[w]));
          }
        }
      }
      else {
        FDKmemset_flex(&pWindowAna__FDK[w][timeSlots],
                       FX_DBL2FX_WIN((FIXP_DBL)0), timeSlots);
      }
    }

    if (hFrameWindow->bFrameKeep == 1) {
      FDKmemcpy_flex(&pWindowAna__FDK[0][2 * timeSlots], 1,
                     &pWindowAna__FDK[0][timeSlots], 1, timeSlots);
      FDKmemcpy_flex(&pWindowAna__FDK[0][timeSlots], 1, pWindowAna__FDK[0], 1,
                     timeSlots);

      if (avoid_keep != 0) {
        FDKmemset_flex(pWindowAna__FDK[0], FX_DBL2FX_WIN((FIXP_DBL)0),
                       timeSlots);
      } else {
        FDKmemset_flex(pWindowAna__FDK[0], winMaxVal, timeSlots);
      }
    }

    pFramingInfo->numParamSets = pFrameWinList->n;
    pFramingInfo->bsFramingType = 1;
    for (ps = 0; ps < pFramingInfo->numParamSets; ps++) {
      pFramingInfo->bsParamSlots[ps] = pFrameWinList->dat[ps].slot;
    }

    if ((pFramingInfo->numParamSets == 1) &&
        (pFramingInfo->bsParamSlots[0] == timeSlots - 1)) {
      pFramingInfo->bsFramingType = 0;
    }
  }

bail:

  return error;
}

FDK_SACENC_ERROR fdk_sacenc_analysisWindowing(
    const INT nTimeSlots, const INT startTimeSlot,
    FIXP_WIN *pFrameWindowAna__FDK, const FIXP_DPK *const *const ppDataIn__FDK,
    FIXP_DPK *const *const ppDataOut__FDK, const INT nHybridBands,
    const INT dim) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if ((pFrameWindowAna__FDK == NULL) || (ppDataIn__FDK == NULL) ||
      (ppDataOut__FDK == NULL)) {
    error = SACENC_INVALID_HANDLE;
  } else {
    int i, ts;
    FIXP_WIN maxVal = FX_DBL2FX_WIN((FIXP_DBL)MAXVAL_DBL);

    if (dim == FW_CHANGE_DIM) {
      for (ts = startTimeSlot; ts < nTimeSlots; ts++) {
        FIXP_WIN win = pFrameWindowAna__FDK[ts];
        if (win == maxVal) {
          for (i = 0; i < nHybridBands; i++) {
            ppDataOut__FDK[i][ts].v.re = ppDataIn__FDK[ts][i].v.re;
            ppDataOut__FDK[i][ts].v.im = ppDataIn__FDK[ts][i].v.im;
          }
        } else {
          for (i = 0; i < nHybridBands; i++) {
            ppDataOut__FDK[i][ts].v.re = fMult(win, ppDataIn__FDK[ts][i].v.re);
            ppDataOut__FDK[i][ts].v.im = fMult(win, ppDataIn__FDK[ts][i].v.im);
          }
        }
      }
    } else {
      for (ts = startTimeSlot; ts < nTimeSlots; ts++) {
        FIXP_WIN win = pFrameWindowAna__FDK[ts];
        if (win == maxVal) {
          for (i = 0; i < nHybridBands; i++) {
            ppDataOut__FDK[ts][i].v.re = ppDataIn__FDK[ts][i].v.re;
            ppDataOut__FDK[ts][i].v.im = ppDataIn__FDK[ts][i].v.im;
          }
        } else {
          for (i = 0; i < nHybridBands; i++) {
            ppDataOut__FDK[ts][i].v.re = fMult(win, ppDataIn__FDK[ts][i].v.re);
            ppDataOut__FDK[ts][i].v.im = fMult(win, ppDataIn__FDK[ts][i].v.im);
          }
        }
      }
    }
  }

  return error;
}
