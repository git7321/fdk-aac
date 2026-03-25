/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef SBR_DEC_H
#define SBR_DEC_H

#include "sbrdecoder.h"

#include "lpp_tran.h"
#include "qmf.h"
#include "env_calc.h"
#include "FDK_audio.h"

#include "sbrdec_drc.h"

#include "pvc_dec.h"

#include "hbe.h"

enum SBRDEC_QMF_SKIP {
  qmfSkipNothing = 0,
  qmfSkipAnalysis = 1 << 0,
  qmfSkipSynthesis = 1 << 1
};

typedef struct {
  SBR_CALCULATE_ENVELOPE SbrCalculateEnvelope;
  SBR_LPP_TRANS LppTrans;
  PVC_STATIC_DATA PvcStaticData;

  SHORT scale_ov;
  SHORT scale_lb;
  SHORT scale_hbe;

  SHORT prev_frame_lSbr;
  SHORT prev_frame_hbeSbr;

  int codecFrameSize;

  HANDLE_HBE_TRANSPOSER hHBE;

  HANDLE_FDK_QMF_DOMAIN_IN qmfDomainInCh;
  HANDLE_FDK_QMF_DOMAIN_OUT qmfDomainOutCh;

  SBRDEC_DRC_CHANNEL sbrDrcChannel;

#if (SBRDEC_MAX_HB_FADE_FRAMES > 0)
  INT highBandFadeCnt;

#endif
  FIXP_DBL **tmp_memory;

  FIXP_DBL **hQmfHBESlotsReal;
  FIXP_DBL **hQmfHBESlotsImag;

  FIXP_DBL **codecQMFBufferReal;
  FIXP_DBL **codecQMFBufferImag;
  UCHAR savedStates;
  int applySbrProc_old;
} SBR_DEC;

typedef SBR_DEC *HANDLE_SBR_DEC;

typedef struct {
  SBR_FRAME_DATA frameData[(1) + 1];
  SBR_PREV_FRAME_DATA prevFrameData;
  SBR_DEC SbrDec;
} SBR_CHANNEL;

typedef SBR_CHANNEL *HANDLE_SBR_CHANNEL;

void sbr_dec(
    HANDLE_SBR_DEC hSbrDec,             /*!< handle to Decoder channel */
    LONG *timeIn,                       /*!< pointer to input time signal */
    LONG *timeOut,                      /*!< pointer to output time signal */
    HANDLE_SBR_DEC hSbrDecRight,        /*!< handle to Decoder channel right */
    LONG *timeOutRight,                 /*!< pointer to output time signal */
    INT strideOut,                      /*!< Time data traversal strideOut */
    HANDLE_SBR_HEADER_DATA hHeaderData, /*!< Static control data */
    HANDLE_SBR_FRAME_DATA hFrameData,   /*!< Control data of current frame */
    HANDLE_SBR_PREV_FRAME_DATA
        hPrevFrameData,        /*!< Some control data of last frame */
    const int applyProcessing, /*!< Flag for SBR operation */
    HANDLE_PS_DEC h_ps_d, const UINT flags, const int codecFrameSize,
    const INT sbrInDataHeadroom);

SBR_ERROR
createSbrDec(SBR_CHANNEL *hSbrChannel, HANDLE_SBR_HEADER_DATA hHeaderData,
             TRANSPOSER_SETTINGS *pSettings, const int downsampleFac,
             const UINT qmfFlags, const UINT flags, const int overlap, int chan,
             int codecFrameSize);

int deleteSbrDec(SBR_CHANNEL *hSbrChannel);

SBR_ERROR
resetSbrDec(HANDLE_SBR_DEC hSbrDec, HANDLE_SBR_HEADER_DATA hHeaderData,
            HANDLE_SBR_PREV_FRAME_DATA hPrevFrameData, const int downsampleFac,
            const UINT flags, HANDLE_SBR_FRAME_DATA hFrameData);

#endif
