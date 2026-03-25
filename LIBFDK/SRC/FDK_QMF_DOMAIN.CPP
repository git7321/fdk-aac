/******************* Library for basic calculation routines ********************

   Author(s):   Matthias Hildenbrand

   Description: Module to efficiently handle QMF data for multiple channels and
                to share the data between e.g. SBR and MPS

*******************************************************************************/

#include "FDK_qmf_domain.h"

#include "common_fix.h"

#define WORKBUFFER1_TAG 0
#define WORKBUFFER3_TAG 4
#define WORKBUFFER4_TAG 5
#define WORKBUFFER6_TAG 7
#define WORKBUFFER7_TAG 8

C_ALLOC_MEM_OVERLAY(QmfWorkBufferCore1, FIXP_DBL, QMF_WB_SECTION_SIZE,
                    SECT_DATA_L1, WORKBUFFER1_TAG)
C_ALLOC_MEM_OVERLAY(QmfWorkBufferCore3, FIXP_DBL, QMF_WB_SECTION_SIZE,
                    SECT_DATA_L2, WORKBUFFER3_TAG)
C_ALLOC_MEM_OVERLAY(QmfWorkBufferCore4, FIXP_DBL, QMF_WB_SECTION_SIZE,
                    SECT_DATA_L2, WORKBUFFER4_TAG)
C_ALLOC_MEM_OVERLAY(QmfWorkBufferCore6, FIXP_DBL, QMF_WB_SECTION_SIZE,
                    SECT_DATA_L2, WORKBUFFER6_TAG)
C_ALLOC_MEM_OVERLAY(QmfWorkBufferCore7, FIXP_DBL, QMF_WB_SECTION_SIZE,
                    SECT_DATA_L2, WORKBUFFER7_TAG)

C_AALLOC_MEM2(AnaQmfStates, FIXP_DBL, 10 * QMF_DOMAIN_MAX_ANALYSIS_QMF_BANDS,
              ((8) + (1)))

C_AALLOC_MEM2(SynQmfStates, FIXP_QSS, 9 * QMF_DOMAIN_MAX_SYNTHESIS_QMF_BANDS,
              ((8) + (1)))

C_ALLOC_MEM2(QmfSlotsReal, FIXP_DBL *,
             QMF_DOMAIN_MAX_TIMESLOTS + QMF_DOMAIN_MAX_OV_TIMESLOTS,
             ((8) + (1)))

C_ALLOC_MEM2(QmfSlotsImag, FIXP_DBL *,
             QMF_DOMAIN_MAX_TIMESLOTS + QMF_DOMAIN_MAX_OV_TIMESLOTS,
             ((8) + (1)))

C_AALLOC_MEM2(QmfOverlapBuffer, FIXP_DBL,
              2 * QMF_DOMAIN_MAX_OV_TIMESLOTS * QMF_DOMAIN_MAX_QMF_PROC_BANDS,
              ((8) + (1)))

C_AALLOC_MEM2(AnaQmfStates16, FIXP_DBL, 10 * QMF_DOMAIN_ANALYSIS_QMF_BANDS_16,
              ((8) + (1)))
C_AALLOC_MEM2(AnaQmfStates24, FIXP_DBL, 10 * QMF_DOMAIN_ANALYSIS_QMF_BANDS_24,
              ((8) + (1)))

C_AALLOC_MEM2(AnaQmfStates32, FIXP_DBL, 10 * QMF_DOMAIN_ANALYSIS_QMF_BANDS_32,
              ((8) + (1)))

C_ALLOC_MEM2(QmfSlotsReal16, FIXP_DBL *,
             QMF_DOMAIN_TIMESLOTS_16 + QMF_DOMAIN_OV_TIMESLOTS_16, ((8) + (1)))

C_ALLOC_MEM2(QmfSlotsReal32, FIXP_DBL *,
             QMF_DOMAIN_TIMESLOTS_32 + QMF_DOMAIN_OV_TIMESLOTS_32, ((8) + (1)))

C_ALLOC_MEM2(QmfSlotsImag16, FIXP_DBL *,
             QMF_DOMAIN_TIMESLOTS_16 + QMF_DOMAIN_OV_TIMESLOTS_16, ((8) + (1)))

C_ALLOC_MEM2(QmfSlotsImag32, FIXP_DBL *,
             QMF_DOMAIN_TIMESLOTS_32 + QMF_DOMAIN_OV_TIMESLOTS_32, ((8) + (1)))

C_AALLOC_MEM2(QmfOverlapBuffer16, FIXP_DBL,
              2 * QMF_DOMAIN_OV_TIMESLOTS_16 * QMF_DOMAIN_MAX_QMF_PROC_BANDS,
              ((8) + (1)))

C_AALLOC_MEM2(QmfOverlapBuffer32, FIXP_DBL,
              2 * QMF_DOMAIN_OV_TIMESLOTS_32 * QMF_DOMAIN_MAX_QMF_PROC_BANDS,
              ((8) + (1)))

static int FDK_QmfDomain_FreePersistentMemory(HANDLE_FDK_QMF_DOMAIN qd) {
  int err = 0;
  int ch;

  for (ch = 0; ch < ((8) + (1)); ch++) {
    if (qd->QmfDomainIn[ch].pAnaQmfStates) {
      if (qd->globalConf.nBandsAnalysis == QMF_DOMAIN_ANALYSIS_QMF_BANDS_16) {
        FreeAnaQmfStates16(&qd->QmfDomainIn[ch].pAnaQmfStates);
      } else if (qd->globalConf.nBandsAnalysis ==
                 QMF_DOMAIN_ANALYSIS_QMF_BANDS_24) {
        FreeAnaQmfStates24(&qd->QmfDomainIn[ch].pAnaQmfStates);
      } else if (qd->globalConf.nBandsAnalysis ==
                 QMF_DOMAIN_ANALYSIS_QMF_BANDS_32) {
        FreeAnaQmfStates32(&qd->QmfDomainIn[ch].pAnaQmfStates);
      } else {
        FreeAnaQmfStates(&qd->QmfDomainIn[ch].pAnaQmfStates);
      }
    }

    if (qd->QmfDomainIn[ch].pOverlapBuffer) {
      if (qd->globalConf.nQmfOvTimeSlots == QMF_DOMAIN_OV_TIMESLOTS_16) {
        FreeQmfOverlapBuffer16(&qd->QmfDomainIn[ch].pOverlapBuffer);
      } else if (qd->globalConf.nQmfOvTimeSlots == QMF_DOMAIN_OV_TIMESLOTS_32) {
        FreeQmfOverlapBuffer32(&qd->QmfDomainIn[ch].pOverlapBuffer);
      } else {
        FreeQmfOverlapBuffer(&qd->QmfDomainIn[ch].pOverlapBuffer);
      }
    }

    if (qd->QmfDomainIn[ch].hQmfSlotsReal) {
      if (qd->globalConf.nQmfTimeSlots == QMF_DOMAIN_TIMESLOTS_16) {
        FreeQmfSlotsReal16(&qd->QmfDomainIn[ch].hQmfSlotsReal);
      } else if (qd->globalConf.nQmfTimeSlots == QMF_DOMAIN_TIMESLOTS_32) {
        FreeQmfSlotsReal32(&qd->QmfDomainIn[ch].hQmfSlotsReal);
      } else {
        FreeQmfSlotsReal(&qd->QmfDomainIn[ch].hQmfSlotsReal);
      }
    }

    if (qd->QmfDomainIn[ch].hQmfSlotsImag) {
      if (qd->globalConf.nQmfTimeSlots == QMF_DOMAIN_TIMESLOTS_16) {
        FreeQmfSlotsImag16(&qd->QmfDomainIn[ch].hQmfSlotsImag);
      }
      if (qd->globalConf.nQmfTimeSlots == QMF_DOMAIN_TIMESLOTS_32) {
        FreeQmfSlotsImag32(&qd->QmfDomainIn[ch].hQmfSlotsImag);
      } else {
        FreeQmfSlotsImag(&qd->QmfDomainIn[ch].hQmfSlotsImag);
      }
    }
  }

  for (ch = 0; ch < ((8) + (1)); ch++) {
    if (qd->QmfDomainOut[ch].pSynQmfStates) {
      FreeSynQmfStates(&qd->QmfDomainOut[ch].pSynQmfStates);
    }
  }

  return err;
}

static int FDK_QmfDomain_AllocatePersistentMemory(HANDLE_FDK_QMF_DOMAIN qd) {
  int err = 0;
  int ch;
  HANDLE_FDK_QMF_DOMAIN_GC gc = &qd->globalConf;

  if ((gc->nInputChannels > ((8) + (1))) || (gc->nOutputChannels > ((8) + (1))))
    return err = 1;
  for (ch = 0; ch < gc->nInputChannels; ch++) {
    int size;

    size = gc->nBandsAnalysis * 10;
    if (size > 0) {
      if (gc->nBandsAnalysis == QMF_DOMAIN_ANALYSIS_QMF_BANDS_16) {
        if (qd->QmfDomainIn[ch].pAnaQmfStates == NULL) {
          if (NULL ==
              (qd->QmfDomainIn[ch].pAnaQmfStates = GetAnaQmfStates16(ch)))
            goto bail;
        }
      } else if (gc->nBandsAnalysis == QMF_DOMAIN_ANALYSIS_QMF_BANDS_24) {
        if (qd->QmfDomainIn[ch].pAnaQmfStates == NULL) {
          if (NULL ==
              (qd->QmfDomainIn[ch].pAnaQmfStates = GetAnaQmfStates24(ch)))
            goto bail;
        }
      } else if (gc->nBandsAnalysis == QMF_DOMAIN_ANALYSIS_QMF_BANDS_32) {
        if (qd->QmfDomainIn[ch].pAnaQmfStates == NULL) {
          if (NULL ==
              (qd->QmfDomainIn[ch].pAnaQmfStates = GetAnaQmfStates32(ch)))
            goto bail;
        }
      } else {
        if (qd->QmfDomainIn[ch].pAnaQmfStates == NULL) {
          if (NULL == (qd->QmfDomainIn[ch].pAnaQmfStates = GetAnaQmfStates(ch)))
            goto bail;
        }
      }
    } else {
      qd->QmfDomainIn[ch].pAnaQmfStates = NULL;
    }

    size = gc->nQmfOvTimeSlots + gc->nQmfTimeSlots;
    if (size > 0) {
      if (gc->nQmfTimeSlots == QMF_DOMAIN_TIMESLOTS_16) {
        if (qd->QmfDomainIn[ch].hQmfSlotsReal == NULL) {
          if (NULL ==
              (qd->QmfDomainIn[ch].hQmfSlotsReal = GetQmfSlotsReal16(ch)))
            goto bail;
        }
        if (qd->QmfDomainIn[ch].hQmfSlotsImag == NULL) {
          if (NULL ==
              (qd->QmfDomainIn[ch].hQmfSlotsImag = GetQmfSlotsImag16(ch)))
            goto bail;
        }
      } else if (gc->nQmfTimeSlots == QMF_DOMAIN_TIMESLOTS_32) {
        if (qd->QmfDomainIn[ch].hQmfSlotsReal == NULL) {
          if (NULL ==
              (qd->QmfDomainIn[ch].hQmfSlotsReal = GetQmfSlotsReal32(ch)))
            goto bail;
        }
        if (qd->QmfDomainIn[ch].hQmfSlotsImag == NULL) {
          if (NULL ==
              (qd->QmfDomainIn[ch].hQmfSlotsImag = GetQmfSlotsImag32(ch)))
            goto bail;
        }
      } else {
        if (qd->QmfDomainIn[ch].hQmfSlotsReal == NULL) {
          if (NULL == (qd->QmfDomainIn[ch].hQmfSlotsReal = GetQmfSlotsReal(ch)))
            goto bail;
        }
        if (qd->QmfDomainIn[ch].hQmfSlotsImag == NULL) {
          if (NULL == (qd->QmfDomainIn[ch].hQmfSlotsImag = GetQmfSlotsImag(ch)))
            goto bail;
        }
      }
    } else {
      qd->QmfDomainIn[ch].hQmfSlotsReal = NULL;
      qd->QmfDomainIn[ch].hQmfSlotsImag = NULL;
    }

    size = gc->nQmfOvTimeSlots * gc->nQmfProcBands * CMPLX_MOD;
    if (size > 0) {
      if (gc->nQmfOvTimeSlots == QMF_DOMAIN_OV_TIMESLOTS_16) {
        if (qd->QmfDomainIn[ch].pOverlapBuffer == NULL) {
          if (NULL ==
              (qd->QmfDomainIn[ch].pOverlapBuffer = GetQmfOverlapBuffer16(ch)))
            goto bail;
        }
      } else if (gc->nQmfOvTimeSlots == QMF_DOMAIN_OV_TIMESLOTS_32) {
        if (qd->QmfDomainIn[ch].pOverlapBuffer == NULL) {
          if (NULL ==
              (qd->QmfDomainIn[ch].pOverlapBuffer = GetQmfOverlapBuffer32(ch)))
            goto bail;
        }
      } else {
        if (qd->QmfDomainIn[ch].pOverlapBuffer == NULL) {
          if (NULL ==
              (qd->QmfDomainIn[ch].pOverlapBuffer = GetQmfOverlapBuffer(ch)))
            goto bail;
        }
      }
    } else {
      qd->QmfDomainIn[ch].pOverlapBuffer = NULL;
    }
  }

  for (ch = 0; ch < gc->nOutputChannels; ch++) {
    int size = gc->nBandsSynthesis * 9;
    if (size > 0) {
      if (qd->QmfDomainOut[ch].pSynQmfStates == NULL) {
        if (NULL == (qd->QmfDomainOut[ch].pSynQmfStates = GetSynQmfStates(ch)))
          goto bail;
      }
    } else {
      qd->QmfDomainOut[ch].pSynQmfStates = NULL;
    }
  }

  return err;

bail:
  FDK_QmfDomain_FreePersistentMemory(qd);
  return -1;
}

QMF_DOMAIN_ERROR FDK_QmfDomain_ClearPersistentMemory(
    HANDLE_FDK_QMF_DOMAIN hqd) {
  QMF_DOMAIN_ERROR err = QMF_DOMAIN_OK;
  int ch, size;
  if (hqd) {
    HANDLE_FDK_QMF_DOMAIN_GC gc = &hqd->globalConf;

    size = gc->nQmfOvTimeSlots * gc->nQmfProcBands * CMPLX_MOD;
    for (ch = 0; ch < gc->nInputChannels; ch++) {
      if (hqd->QmfDomainIn[ch].pOverlapBuffer) {
        FDKmemclear(hqd->QmfDomainIn[ch].pOverlapBuffer,
                    size * sizeof(FIXP_DBL));
      }
    }
    if (FDK_QmfDomain_InitFilterBank(hqd, 0)) {
      err = QMF_DOMAIN_INIT_ERROR;
    }
  } else {
    err = QMF_DOMAIN_INIT_ERROR;
  }
  return err;
}

static FIXP_DBL *FDK_getWorkBuffer(FIXP_DBL **pWorkBuffer,
                                   USHORT workBufferOffset,
                                   USHORT workBufferSectSize, USHORT memSize) {
  int idx1;
  int idx2;
  FIXP_DBL *pwb;

  FDK_ASSERT((workBufferSectSize % 64) == 0);
  idx2 = workBufferOffset % workBufferSectSize;
  idx1 = (workBufferOffset - idx2) / workBufferSectSize;
  FDK_ASSERT(idx1 < QMF_MAX_WB_SECTIONS);
  FDK_ASSERT(pWorkBuffer[idx1] != NULL);
  FDK_ASSERT((idx2 + memSize) <= workBufferSectSize);
  pwb = &pWorkBuffer[idx1][idx2];

  return pwb;
}

static int FDK_QmfDomain_FeedWorkBuffer(HANDLE_FDK_QMF_DOMAIN qd, int ch,
                                        FIXP_DBL **pWorkBuffer,
                                        USHORT workBufferOffset,
                                        USHORT workBufferSectSize, int size) {
  int err = 0;
  int mem_needed;

  mem_needed = qd->QmfDomainIn[ch].workBuf_nBands *
               qd->QmfDomainIn[ch].workBuf_nTimeSlots * CMPLX_MOD;
  if (mem_needed > size) {
    return (err = 1);
  }
  qd->QmfDomainIn[ch].pWorkBuffer = pWorkBuffer;
  qd->QmfDomainIn[ch].workBufferOffset = workBufferOffset;
  qd->QmfDomainIn[ch].workBufferSectSize = workBufferSectSize;

  return err;
}

int FDK_QmfDomain_IsInitialized(const HANDLE_FDK_QMF_DOMAIN qd) {
  FDK_ASSERT(qd != NULL);
  return ((qd->QmfDomainIn[0].pAnaQmfStates == NULL) &&
          (qd->QmfDomainOut[0].pSynQmfStates == NULL))
             ? 0
             : 1;
}

int FDK_QmfDomain_InitFilterBank(HANDLE_FDK_QMF_DOMAIN qd, UINT extra_flags) {
  FDK_ASSERT(qd != NULL);
  int err = 0;
  int ch, ts;
  HANDLE_FDK_QMF_DOMAIN_GC gc = &qd->globalConf;
  int noCols = gc->nQmfTimeSlots;
  int lsb = gc->nBandsAnalysis;
  int usb = fMin((INT)gc->nBandsSynthesis, 64);
  int nProcBands = gc->nQmfProcBands;
  FDK_ASSERT(nProcBands % ALIGNMENT_DEFAULT == 0);

  if (extra_flags & QMF_FLAG_MPSLDFB) {
    gc->flags &= ~QMF_FLAG_CLDFB;
    gc->flags |= QMF_FLAG_MPSLDFB;
  }
  for (ch = 0; ch < gc->nInputChannels; ch++) {
    FIXP_DBL *ptrOv =
        qd->QmfDomainIn[ch].pOverlapBuffer;
    if ((ptrOv == NULL) && (gc->nQmfOvTimeSlots != 0)) {
      err = 1;
      return err;
    }
    FIXP_DBL **ptr =
        qd->QmfDomainIn[fMin(ch, fMax((INT)gc->nQmfProcChannels - 1, 0))]
            .pWorkBuffer;
    USHORT workBufferOffset =
        qd->QmfDomainIn[fMin(ch, fMax((INT)gc->nQmfProcChannels - 1, 0))]
            .workBufferOffset;
    USHORT workBufferSectSize =
        qd->QmfDomainIn[fMin(ch, fMax((INT)gc->nQmfProcChannels - 1, 0))]
            .workBufferSectSize;

    if ((ptr == NULL) && (gc->nQmfTimeSlots != 0)) {
      err = 1;
      return err;
    }

    qd->QmfDomainIn[ch].pGlobalConf = gc;
    for (ts = 0; ts < gc->nQmfOvTimeSlots; ts++) {
      qd->QmfDomainIn[ch].hQmfSlotsReal[ts] = ptrOv;
      ptrOv += nProcBands;
      qd->QmfDomainIn[ch].hQmfSlotsImag[ts] = ptrOv;
      ptrOv += nProcBands;
    }
    for (; ts < (gc->nQmfOvTimeSlots + gc->nQmfTimeSlots); ts++) {
      qd->QmfDomainIn[ch].hQmfSlotsReal[ts] = FDK_getWorkBuffer(
          ptr, workBufferOffset, workBufferSectSize, nProcBands);
      workBufferOffset += nProcBands;
      qd->QmfDomainIn[ch].hQmfSlotsImag[ts] = FDK_getWorkBuffer(
          ptr, workBufferOffset, workBufferSectSize, nProcBands);
      workBufferOffset += nProcBands;
    }
    err |= qmfInitAnalysisFilterBank(
        &qd->QmfDomainIn[ch].fb, qd->QmfDomainIn[ch].pAnaQmfStates, noCols,
        (qd->QmfDomainIn[ch].fb.lsb == 0) ? lsb : qd->QmfDomainIn[ch].fb.lsb,
        (qd->QmfDomainIn[ch].fb.usb == 0) ? usb : qd->QmfDomainIn[ch].fb.usb,
        gc->nBandsAnalysis, gc->flags | extra_flags);
  }

  for (ch = 0; ch < gc->nOutputChannels; ch++) {
    FIXP_DBL outGain_m = qd->QmfDomainOut[ch].fb.outGain_m;
    int outGain_e = qd->QmfDomainOut[ch].fb.outGain_e;
    int outScale = qmfGetOutScalefactor(&qd->QmfDomainOut[ch].fb);
    err |= qmfInitSynthesisFilterBank(
        &qd->QmfDomainOut[ch].fb, qd->QmfDomainOut[ch].pSynQmfStates, noCols,
        (qd->QmfDomainOut[ch].fb.lsb == 0) ? lsb : qd->QmfDomainOut[ch].fb.lsb,
        (qd->QmfDomainOut[ch].fb.usb == 0) ? usb : qd->QmfDomainOut[ch].fb.usb,
        gc->nBandsSynthesis, gc->flags | extra_flags);
    if (outGain_m != (FIXP_DBL)0) {
      qmfChangeOutGain(&qd->QmfDomainOut[ch].fb, outGain_m, outGain_e);
    }
    if (outScale) {
      qmfChangeOutScalefactor(&qd->QmfDomainOut[ch].fb, outScale);
    }
  }

  return err;
}

void FDK_QmfDomain_SaveOverlap(HANDLE_FDK_QMF_DOMAIN_IN qd_ch, int offset) {
  FDK_ASSERT(qd_ch != NULL);
  int ts;
  HANDLE_FDK_QMF_DOMAIN_GC gc = qd_ch->pGlobalConf;
  int ovSlots = gc->nQmfOvTimeSlots;
  int nCols = gc->nQmfTimeSlots;
  int nProcBands = gc->nQmfProcBands;
  FIXP_DBL **qmfReal = qd_ch->hQmfSlotsReal;
  FIXP_DBL **qmfImag = qd_ch->hQmfSlotsImag;
  QMF_SCALE_FACTOR *pScaling = &qd_ch->scaling;

  if (qmfImag != NULL) {
    for (ts = offset; ts < ovSlots; ts++) {
      FDKmemcpy(qmfReal[ts], qmfReal[nCols + ts],
                sizeof(FIXP_DBL) * nProcBands);
      FDKmemcpy(qmfImag[ts], qmfImag[nCols + ts],
                sizeof(FIXP_DBL) * nProcBands);
    }
  } else {
    for (ts = 0; ts < ovSlots; ts++) {
      FDKmemcpy(qmfReal[ts], qmfReal[nCols + ts],
                sizeof(FIXP_DBL) * nProcBands);
    }
  }
  pScaling->ov_lb_scale = pScaling->lb_scale;
}

#define SCALE2EXP(s) (15 - (s))
#define EXP2SCALE(e) (15 - (e))

void FDK_QmfDomain_GetSlot(const HANDLE_FDK_QMF_DOMAIN_IN qd_ch, const int ts,
                           const int start_band, const int stop_band,
                           FIXP_DBL *pQmfOutReal, FIXP_DBL *pQmfOutImag,
                           const int exp_out) {
  FDK_ASSERT(qd_ch != NULL);
  FDK_ASSERT(pQmfOutReal != NULL);
  HANDLE_FDK_QMF_DOMAIN_GC gc = qd_ch->pGlobalConf;
  const FIXP_DBL *real = qd_ch->hQmfSlotsReal[ts];
  const FIXP_DBL *imag = qd_ch->hQmfSlotsImag[ts];
  const int ovSlots = gc->nQmfOvTimeSlots;
  const int exp_lb = SCALE2EXP((ts < ovSlots) ? qd_ch->scaling.ov_lb_scale
                                              : qd_ch->scaling.lb_scale);
  const int exp_hb = SCALE2EXP(qd_ch->scaling.hb_scale);
  const int lsb = qd_ch->fb.lsb;
  const int usb = qd_ch->fb.usb;
  int b = start_band;
  int lb_sf, hb_sf;

  int target_exp =
      ALGORITHMIC_SCALING_IN_ANALYSIS_FILTERBANK + qd_ch->fb.filterScale;

  FDK_ASSERT(ts < (gc->nQmfTimeSlots + gc->nQmfOvTimeSlots));
  FDK_ASSERT(start_band >= 0);
  FDK_ASSERT(stop_band <= gc->nQmfProcBands);

  if (qd_ch->fb.no_channels == 24) {
    target_exp -= 1;
  }

  lb_sf = fMax(exp_lb - target_exp - exp_out, -31);
  FDK_ASSERT(lb_sf < 32);
  hb_sf = fMax(exp_hb - target_exp - exp_out, -31);
  FDK_ASSERT(hb_sf < 32);

  if (pQmfOutImag == NULL) {
    for (; b < fMin(lsb, stop_band); b++) {
      pQmfOutReal[b] = scaleValueSaturate(real[b], lb_sf);
    }
    for (; b < fMin(usb, stop_band); b++) {
      pQmfOutReal[b] = scaleValueSaturate(real[b], hb_sf);
    }
    for (; b < stop_band; b++) {
      pQmfOutReal[b] = (FIXP_DBL)0;
    }
  } else {
    FDK_ASSERT(imag != NULL);
    for (; b < fMin(lsb, stop_band); b++) {
      pQmfOutReal[b] = scaleValueSaturate(real[b], lb_sf);
      pQmfOutImag[b] = scaleValueSaturate(imag[b], lb_sf);
    }
    for (; b < fMin(usb, stop_band); b++) {
      pQmfOutReal[b] = scaleValueSaturate(real[b], hb_sf);
      pQmfOutImag[b] = scaleValueSaturate(imag[b], hb_sf);
    }
    for (; b < stop_band; b++) {
      pQmfOutReal[b] = (FIXP_DBL)0;
      pQmfOutImag[b] = (FIXP_DBL)0;
    }
  }
}

void FDK_QmfDomain_GetWorkBuffer(const HANDLE_FDK_QMF_DOMAIN_IN qd_ch,
                                 const int ts, FIXP_DBL **ppQmfReal,
                                 FIXP_DBL **ppQmfImag) {
  FDK_ASSERT(qd_ch != NULL);
  FDK_ASSERT(ppQmfReal != NULL);
  FDK_ASSERT(ppQmfImag != NULL);
  const int bands = qd_ch->workBuf_nBands;
  FIXP_DBL **pWorkBuf = qd_ch->pWorkBuffer;
  USHORT workBufferOffset = qd_ch->workBufferOffset;
  USHORT workBufferSectSize = qd_ch->workBufferSectSize;

  FDK_ASSERT(bands > 0);
  FDK_ASSERT(ts < qd_ch->workBuf_nTimeSlots);

  *ppQmfReal = FDK_getWorkBuffer(
      pWorkBuf, workBufferOffset + (ts * CMPLX_MOD + 0) * bands,
      workBufferSectSize, bands);
  *ppQmfImag = FDK_getWorkBuffer(
      pWorkBuf, workBufferOffset + (ts * CMPLX_MOD + 1) * bands,
      workBufferSectSize, bands);
}

void FDK_QmfDomain_WorkBuffer2ProcChannel(
    const HANDLE_FDK_QMF_DOMAIN_IN qd_ch) {
  FDK_ASSERT(qd_ch != NULL);
  HANDLE_FDK_QMF_DOMAIN_GC gc = qd_ch->pGlobalConf;
  FIXP_DBL **pWorkBuf = qd_ch->pWorkBuffer;
  USHORT workBufferOffset = qd_ch->workBufferOffset;
  USHORT workBufferSectSize = qd_ch->workBufferSectSize;

  if (FDK_getWorkBuffer(pWorkBuf, workBufferOffset, workBufferSectSize,
                        qd_ch->workBuf_nBands) ==
      qd_ch->hQmfSlotsReal[gc->nQmfOvTimeSlots]) {
    return;
  } else {
    const int bands = qd_ch->workBuf_nBands;
    const int slots = qd_ch->workBuf_nTimeSlots;
    int ts;
    for (ts = 0; ts < slots; ts++) {
      FDKmemcpy(qd_ch->hQmfSlotsReal[gc->nQmfOvTimeSlots + ts],
                FDK_getWorkBuffer(pWorkBuf, workBufferOffset,
                                  workBufferSectSize, bands),
                sizeof(FIXP_DBL) * bands);
      workBufferOffset += bands;
      FDKmemcpy(qd_ch->hQmfSlotsImag[gc->nQmfOvTimeSlots + ts],
                FDK_getWorkBuffer(pWorkBuf, workBufferOffset,
                                  workBufferSectSize, bands),
                sizeof(FIXP_DBL) * bands);
      workBufferOffset += bands;
    }
  }
}

void FDK_QmfDomain_QmfData2HBE(HANDLE_FDK_QMF_DOMAIN_IN qd_ch,
                               FIXP_DBL **ppQmfReal, FIXP_DBL **ppQmfImag) {
  FDK_ASSERT(qd_ch != NULL);
  FDK_ASSERT(ppQmfReal != NULL);
  FDK_ASSERT(ppQmfImag != NULL);
  HANDLE_FDK_QMF_DOMAIN_GC gc = qd_ch->pGlobalConf;
  FIXP_DBL **pWorkBuf = qd_ch->pWorkBuffer;
  USHORT workBufferOffset = qd_ch->workBufferOffset;
  USHORT workBufferSectSize = qd_ch->workBufferSectSize;

  if (FDK_getWorkBuffer(pWorkBuf, workBufferOffset, workBufferSectSize,
                        qd_ch->workBuf_nBands) ==
      qd_ch->hQmfSlotsReal[gc->nQmfOvTimeSlots]) {
    int ts;
    const int bands = gc->nBandsAnalysis;
    const int slots = qd_ch->workBuf_nTimeSlots;
    FDK_ASSERT(bands <= 64);
    for (ts = 0; ts < slots; ts++) {
      FIXP_DBL tmp[64];
      FDKmemcpy(tmp, qd_ch->hQmfSlotsReal[gc->nQmfOvTimeSlots + ts],
                sizeof(FIXP_DBL) * bands);
      FDKmemcpy(qd_ch->hQmfSlotsReal[gc->nQmfOvTimeSlots + ts], ppQmfReal[ts],
                sizeof(FIXP_DBL) * bands);
      FDKmemcpy(ppQmfReal[ts], tmp, sizeof(FIXP_DBL) * bands);
      FDKmemcpy(tmp, qd_ch->hQmfSlotsImag[gc->nQmfOvTimeSlots + ts],
                sizeof(FIXP_DBL) * bands);
      FDKmemcpy(qd_ch->hQmfSlotsImag[gc->nQmfOvTimeSlots + ts], ppQmfImag[ts],
                sizeof(FIXP_DBL) * bands);
      FDKmemcpy(ppQmfImag[ts], tmp, sizeof(FIXP_DBL) * bands);
    }
  } else {
    const int bands = qd_ch->workBuf_nBands;
    const int slots = qd_ch->workBuf_nTimeSlots;
    int ts;
    FDK_ASSERT(qd_ch->workBuf_nBands == gc->nBandsAnalysis);
    for (ts = 0; ts < slots; ts++) {
      FDKmemcpy(qd_ch->hQmfSlotsReal[gc->nQmfOvTimeSlots + ts], ppQmfReal[ts],
                sizeof(FIXP_DBL) * bands);
      FDKmemcpy(qd_ch->hQmfSlotsImag[gc->nQmfOvTimeSlots + ts], ppQmfImag[ts],
                sizeof(FIXP_DBL) * bands);
      FDKmemcpy(ppQmfReal[ts],
                FDK_getWorkBuffer(pWorkBuf, workBufferOffset,
                                  workBufferSectSize, bands),
                sizeof(FIXP_DBL) * bands);
      workBufferOffset += bands;
      FDKmemcpy(ppQmfImag[ts],
                FDK_getWorkBuffer(pWorkBuf, workBufferOffset,
                                  workBufferSectSize, bands),
                sizeof(FIXP_DBL) * bands);
      workBufferOffset += bands;
    }
  }
}

void FDK_QmfDomain_ClearRequested(HANDLE_FDK_QMF_DOMAIN_GC hgc) {
  hgc->qmfDomainExplicitConfig = 0;
  hgc->flags_requested = 0;
  hgc->nInputChannels_requested = 0;
  hgc->nOutputChannels_requested = 0;
  hgc->parkChannel_requested = 0;
  hgc->nBandsAnalysis_requested = 0;
  hgc->nBandsSynthesis_requested = 0;
  hgc->nQmfTimeSlots_requested = 0;
  hgc->nQmfOvTimeSlots_requested = 0;
  hgc->nQmfProcBands_requested = 0;
  hgc->nQmfProcChannels_requested = 0;
}

static void FDK_QmfDomain_ClearConfigured(HANDLE_FDK_QMF_DOMAIN_GC hgc) {
  hgc->flags = 0;
  hgc->nInputChannels = 0;
  hgc->nOutputChannels = 0;
  hgc->parkChannel = 0;
  hgc->nBandsAnalysis = 0;
  hgc->nBandsSynthesis = 0;
  hgc->nQmfTimeSlots = 0;
  hgc->nQmfOvTimeSlots = 0;
  hgc->nQmfProcBands = 0;
  hgc->nQmfProcChannels = 0;
}

static void FDK_QmfDomain_ClearFilterBank(HANDLE_FDK_QMF_DOMAIN hqd) {
  int ch;

  for (ch = 0; ch < ((8) + (1)); ch++) {
    FDKmemclear(&hqd->QmfDomainIn[ch].fb, sizeof(hqd->QmfDomainIn[ch].fb));
  }

  for (ch = 0; ch < ((8) + (1)); ch++) {
    FDKmemclear(&hqd->QmfDomainOut[ch].fb, sizeof(hqd->QmfDomainIn[ch].fb));
  }
}

QMF_DOMAIN_ERROR FDK_QmfDomain_Configure(HANDLE_FDK_QMF_DOMAIN hqd) {
  FDK_ASSERT(hqd != NULL);
  QMF_DOMAIN_ERROR err = QMF_DOMAIN_OK;
  int i, size_main, size, size_temp = 0;

  HANDLE_FDK_QMF_DOMAIN_GC hgc = &hqd->globalConf;
  FIXP_DBL **pWorkBuffer = hgc->pWorkBuffer;

  int hasChanged = 0;

  if ((hgc->nQmfProcChannels_requested > 0) &&
      (hgc->nQmfProcBands_requested != 64)) {
    return QMF_DOMAIN_INIT_ERROR;
  }
  if (hgc->nBandsAnalysis_requested > hgc->nQmfProcBands_requested) {
    return QMF_DOMAIN_INIT_ERROR;
  }

  if ((hgc->nQmfProcChannels != hgc->nQmfProcChannels_requested) ||
      (hgc->nQmfProcBands != hgc->nQmfProcBands_requested) ||
      (hgc->nQmfTimeSlots != hgc->nQmfTimeSlots_requested)) {
    for (i = 0; i < hgc->nQmfProcChannels_requested; i++) {
      hqd->QmfDomainIn[i].workBuf_nBands = hgc->nQmfProcBands_requested;
      hgc->nQmfProcBands = hgc->nQmfProcBands_requested;

      hqd->QmfDomainIn[i].workBuf_nTimeSlots = hgc->nQmfTimeSlots_requested;
    }

    hgc->nQmfProcChannels =
        hgc->nQmfProcChannels_requested;

    hasChanged = 1;
  }

  if ((hgc->nInputChannels != hgc->nInputChannels_requested) ||
      (hgc->nBandsAnalysis != hgc->nBandsAnalysis_requested) ||
      (hgc->nQmfTimeSlots != hgc->nQmfTimeSlots_requested) ||
      (hgc->nQmfOvTimeSlots != hgc->nQmfOvTimeSlots_requested) ||
      (hgc->nOutputChannels != hgc->nOutputChannels_requested) ||
      (hgc->nBandsSynthesis != hgc->nBandsSynthesis_requested) ||
      (hgc->parkChannel != hgc->parkChannel_requested)) {
    hgc->nInputChannels = hgc->nInputChannels_requested;
    hgc->nBandsAnalysis = hgc->nBandsAnalysis_requested;
    hgc->nQmfTimeSlots = hgc->nQmfTimeSlots_requested;
    hgc->nQmfOvTimeSlots = hgc->nQmfOvTimeSlots_requested;
    hgc->nOutputChannels = hgc->nOutputChannels_requested;
    hgc->nBandsSynthesis = hgc->nBandsSynthesis_requested;
    hgc->parkChannel = hgc->parkChannel_requested;

    if (FDK_QmfDomain_AllocatePersistentMemory(hqd)) {
      err = QMF_DOMAIN_OUT_OF_MEMORY;
      goto bail;
    }

    if ((hgc->nBandsAnalysis == 32) && (hgc->nBandsSynthesis == 32) &&
        !(hgc->flags & (QMF_FLAG_CLDFB | QMF_FLAG_MPSLDFB))) {
      hgc->flags_requested |= QMF_FLAG_DOWNSAMPLED;
    }

    hasChanged = 1;
  }

  if (hgc->flags != hgc->flags_requested) {
    if ((hgc->flags_requested & QMF_FLAG_MPSLDFB) &&
        (hgc->flags_requested & QMF_FLAG_CLDFB)) {
      hgc->flags_requested &= ~QMF_FLAG_CLDFB;
    }
    hgc->flags = hgc->flags_requested;
    hasChanged = 1;
  }

  if (hasChanged) {

    if (hgc->parkChannel && (hqd->globalConf.nQmfProcChannels == 1)) {
      hqd->QmfDomainIn[1].workBuf_nBands = hqd->globalConf.nBandsAnalysis;
      hqd->QmfDomainIn[1].workBuf_nTimeSlots = hqd->globalConf.nQmfTimeSlots;
      size_temp = hqd->QmfDomainIn[1].workBuf_nBands *
                  hqd->QmfDomainIn[1].workBuf_nTimeSlots * CMPLX_MOD;
    }

    size_main = hqd->QmfDomainIn[0].workBuf_nBands *
                hqd->QmfDomainIn[0].workBuf_nTimeSlots * CMPLX_MOD;

    size = size_main * hgc->nQmfProcChannels + size_temp;

    if (size > (QMF_MAX_WB_SECTIONS * QMF_WB_SECTION_SIZE)) {
      err = QMF_DOMAIN_OUT_OF_MEMORY;
      goto bail;
    }

    if ((size > 0 ) && (pWorkBuffer[0] == NULL)) {
      pWorkBuffer[0] = GetQmfWorkBufferCore6();
    }

    if ((size > 1 * QMF_WB_SECTION_SIZE) && (pWorkBuffer[1] == NULL)) {
      pWorkBuffer[1] = GetQmfWorkBufferCore1();
    }

    if ((size > 2 * QMF_WB_SECTION_SIZE) && (pWorkBuffer[2] == NULL)) {
      pWorkBuffer[2] = GetQmfWorkBufferCore3();
    }

    if ((size > 3 * QMF_WB_SECTION_SIZE) && (pWorkBuffer[3] == NULL)) {
      pWorkBuffer[3] = GetQmfWorkBufferCore4();
    }

    if ((size > 4 * QMF_WB_SECTION_SIZE) && (pWorkBuffer[4] == NULL)) {
      pWorkBuffer[4] = GetQmfWorkBufferCore7();
    }

    for (i = 0; i < hgc->nQmfProcChannels; i++) {
      FDK_QmfDomain_FeedWorkBuffer(hqd, i, pWorkBuffer, size_main * i,
                                   QMF_WB_SECTION_SIZE, size_main);
    }
    if (hgc->parkChannel) {
      for (; i < hgc->nInputChannels; i++) {
        FDK_QmfDomain_FeedWorkBuffer(hqd, 1, pWorkBuffer,
                                     size_main * hgc->nQmfProcChannels,
                                     QMF_WB_SECTION_SIZE, size_temp);
      }
    }

    for (i = 0; i < hgc->nOutputChannels; i++) {
      if ((hqd->QmfDomainOut[i].fb.lsb == 0) &&
          (hqd->QmfDomainOut[i].fb.usb == 0)) {
        hqd->QmfDomainOut[i].fb.lsb = hgc->nBandsAnalysis_requested;
        hqd->QmfDomainOut[i].fb.usb =
            fMin((INT)hgc->nBandsSynthesis_requested, 64);
      }
    }
    if (FDK_QmfDomain_InitFilterBank(hqd, 0)) {
      err = QMF_DOMAIN_INIT_ERROR;
    }
  }

bail:
  if (err) {
    FDK_QmfDomain_FreeMem(hqd);
  }
  return err;
}

static void FDK_QmfDomain_FreeWorkBuffer(HANDLE_FDK_QMF_DOMAIN hqd) {
  FIXP_DBL **pWorkBuffer = hqd->globalConf.pWorkBuffer;

  if (pWorkBuffer[0]) FreeQmfWorkBufferCore6(&pWorkBuffer[0]);
  if (pWorkBuffer[1]) FreeQmfWorkBufferCore1(&pWorkBuffer[1]);
  if (pWorkBuffer[2]) FreeQmfWorkBufferCore3(&pWorkBuffer[2]);
  if (pWorkBuffer[3]) FreeQmfWorkBufferCore4(&pWorkBuffer[3]);
  if (pWorkBuffer[4]) FreeQmfWorkBufferCore7(&pWorkBuffer[4]);
}

void FDK_QmfDomain_FreeMem(HANDLE_FDK_QMF_DOMAIN hqd) {
  FDK_QmfDomain_FreeWorkBuffer(hqd);

  FDK_QmfDomain_FreePersistentMemory(hqd);

  FDK_QmfDomain_ClearFilterBank(hqd);

  FDK_QmfDomain_ClearConfigured(&hqd->globalConf);

  FDK_QmfDomain_ClearRequested(&hqd->globalConf);
}

void FDK_QmfDomain_Close(HANDLE_FDK_QMF_DOMAIN hqd) {
  FDK_QmfDomain_FreeWorkBuffer(hqd);

  FDK_QmfDomain_FreePersistentMemory(hqd);
}
