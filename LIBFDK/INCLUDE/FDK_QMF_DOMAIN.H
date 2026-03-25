/******************* Library for basic calculation routines ********************

   Author(s):   Matthias Hildenbrand

   Description: Module to efficiently handle QMF data for multiple channels and
                to share the data between e.g. SBR and MPS

*******************************************************************************/

#ifndef FDK_QMF_DOMAIN_H
#define FDK_QMF_DOMAIN_H

#include "qmf.h"

typedef enum {
  QMF_DOMAIN_OK = 0x0,
  QMF_DOMAIN_OUT_OF_MEMORY =
      0x1,
  QMF_DOMAIN_INIT_ERROR =
      0x2,
  QMF_DOMAIN_RESAMPLER_INIT_ERROR =
      0x3
} QMF_DOMAIN_ERROR;

#define CMPLX_MOD (2)

#define QMF_MAX_WB_SECTIONS (5)
#define QMF_WB_SECTION_SIZE (1024 * 2)

H_ALLOC_MEM_OVERLAY(QmfWorkBufferCore1, FIXP_DBL)
H_ALLOC_MEM_OVERLAY(QmfWorkBufferCore3, FIXP_DBL)
H_ALLOC_MEM_OVERLAY(QmfWorkBufferCore4, FIXP_DBL)
H_ALLOC_MEM_OVERLAY(QmfWorkBufferCore6, FIXP_DBL)
H_ALLOC_MEM_OVERLAY(QmfWorkBufferCore7, FIXP_DBL)

#define QMF_DOMAIN_MAX_ANALYSIS_QMF_BANDS (64)
#define QMF_DOMAIN_MAX_SYNTHESIS_QMF_BANDS (QMF_MAX_SYNTHESIS_BANDS)
#define QMF_DOMAIN_MAX_QMF_PROC_BANDS (64)
#define QMF_DOMAIN_MAX_TIMESLOTS (64)
#define QMF_DOMAIN_MAX_OV_TIMESLOTS (12)

#define QMF_DOMAIN_ANALYSIS_QMF_BANDS_16 (16)
#define QMF_DOMAIN_ANALYSIS_QMF_BANDS_24 (24)
#define QMF_DOMAIN_ANALYSIS_QMF_BANDS_32 (32)

#define QMF_DOMAIN_TIMESLOTS_16 (16)
#define QMF_DOMAIN_TIMESLOTS_32 (32)

#define QMF_DOMAIN_OV_TIMESLOTS_16 (3)
#define QMF_DOMAIN_OV_TIMESLOTS_32 (6)

H_ALLOC_MEM(AnaQmfStates, FIXP_DBL)
H_ALLOC_MEM(SynQmfStates, FIXP_QSS)
H_ALLOC_MEM(QmfSlotsReal, FIXP_DBL *)
H_ALLOC_MEM(QmfSlotsImag, FIXP_DBL *)
H_ALLOC_MEM(QmfOverlapBuffer, FIXP_DBL)

H_ALLOC_MEM(AnaQmfStates16, FIXP_DBL)
H_ALLOC_MEM(AnaQmfStates24, FIXP_DBL)
H_ALLOC_MEM(AnaQmfStates32, FIXP_DBL)
H_ALLOC_MEM(QmfSlotsReal16, FIXP_DBL *)
H_ALLOC_MEM(QmfSlotsReal32, FIXP_DBL *)
H_ALLOC_MEM(QmfSlotsImag16, FIXP_DBL *)
H_ALLOC_MEM(QmfSlotsImag32, FIXP_DBL *)
H_ALLOC_MEM(QmfOverlapBuffer16, FIXP_DBL)
H_ALLOC_MEM(QmfOverlapBuffer32, FIXP_DBL)

typedef struct {
  UCHAR qmfDomainExplicitConfig;
  UCHAR nInputChannels;
  UCHAR nInputChannels_requested;
  UCHAR nOutputChannels;
  UCHAR nOutputChannels_requested;
  UCHAR
  parkChannel;
  UCHAR parkChannel_requested;
  FIXP_DBL *
      pWorkBuffer[QMF_MAX_WB_SECTIONS];
  UINT flags;
  UINT flags_requested;
  UCHAR nBandsAnalysis;
  UCHAR nBandsAnalysis_requested;
  USHORT nBandsSynthesis;
  USHORT
  nBandsSynthesis_requested;
  UCHAR nQmfTimeSlots;
  UCHAR nQmfTimeSlots_requested;
  UCHAR
  nQmfOvTimeSlots;
  UCHAR nQmfOvTimeSlots_requested;
  UCHAR nQmfProcBands;
  UCHAR nQmfProcBands_requested;
  UCHAR
  nQmfProcChannels;
  UCHAR
  nQmfProcChannels_requested;
} FDK_QMF_DOMAIN_GC;
typedef FDK_QMF_DOMAIN_GC *HANDLE_FDK_QMF_DOMAIN_GC;

typedef struct {
  HANDLE_FDK_QMF_DOMAIN_GC
  pGlobalConf;
  QMF_FILTER_BANK fb;
  QMF_SCALE_FACTOR scaling;
  UCHAR workBuf_nTimeSlots;
  UCHAR workBuf_nBands;
  USHORT workBufferOffset;
  USHORT workBufferSectSize;
  FIXP_DBL *
      pAnaQmfStates;
  FIXP_DBL
  *pOverlapBuffer;
  FIXP_DBL **pWorkBuffer;
  FIXP_DBL *
      *hQmfSlotsReal;
  FIXP_DBL **hQmfSlotsImag;
} FDK_QMF_DOMAIN_IN;
typedef FDK_QMF_DOMAIN_IN *HANDLE_FDK_QMF_DOMAIN_IN;

typedef struct {
  QMF_FILTER_BANK fb;
  FIXP_QSS *pSynQmfStates;
} FDK_QMF_DOMAIN_OUT;
typedef FDK_QMF_DOMAIN_OUT *HANDLE_FDK_QMF_DOMAIN_OUT;

typedef struct {
  FDK_QMF_DOMAIN_GC globalConf;
  FDK_QMF_DOMAIN_IN
  QmfDomainIn[((8) + (1))];
  FDK_QMF_DOMAIN_OUT
  QmfDomainOut[((8) + (1))];
} FDK_QMF_DOMAIN;
typedef FDK_QMF_DOMAIN *HANDLE_FDK_QMF_DOMAIN;

int FDK_QmfDomain_IsInitialized(const HANDLE_FDK_QMF_DOMAIN qd);

int FDK_QmfDomain_InitFilterBank(HANDLE_FDK_QMF_DOMAIN qd, UINT extra_flags);

void FDK_QmfDomain_SaveOverlap(HANDLE_FDK_QMF_DOMAIN_IN qd_ch, int offset);

void FDK_QmfDomain_GetSlot(const HANDLE_FDK_QMF_DOMAIN_IN qd_ch, const int ts,
                           const int start_band, const int stop_band,
                           FIXP_DBL *pQmfOutReal, FIXP_DBL *pQmfOutImag,
                           const int exp_out);

void FDK_QmfDomain_GetWorkBuffer(const HANDLE_FDK_QMF_DOMAIN_IN qd_ch,
                                 const int ts, FIXP_DBL **ppQmfReal,
                                 FIXP_DBL **ppQmfImag);

void FDK_QmfDomain_WorkBuffer2ProcChannel(const HANDLE_FDK_QMF_DOMAIN_IN qd_ch);

void FDK_QmfDomain_QmfData2HBE(HANDLE_FDK_QMF_DOMAIN_IN qd_ch,
                               FIXP_DBL **ppQmfReal, FIXP_DBL **ppQmfImag);

void FDK_QmfDomain_ClearRequested(HANDLE_FDK_QMF_DOMAIN_GC hgc);

QMF_DOMAIN_ERROR FDK_QmfDomain_Configure(HANDLE_FDK_QMF_DOMAIN hqd);

void FDK_QmfDomain_FreeMem(HANDLE_FDK_QMF_DOMAIN hqd);

QMF_DOMAIN_ERROR FDK_QmfDomain_ClearPersistentMemory(HANDLE_FDK_QMF_DOMAIN hqd);

void FDK_QmfDomain_Close(HANDLE_FDK_QMF_DOMAIN hqd);

#endif /* FDK_QMF_DOMAIN_H */
