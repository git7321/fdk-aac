/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Quantizing & coding

*******************************************************************************/

#ifndef QC_MAIN_H
#define QC_MAIN_H

#include "aacenc.h"
#include "qc_data.h"
#include "interface.h"
#include "psy_main.h"
#include "tpenc_lib.h"

AAC_ENCODER_ERROR FDKaacEnc_QCOutNew(QC_OUT **phQC, const INT nElements,
                                     const INT nChannels, const INT nSubFrames,
                                     UCHAR *dynamic_RAM);

AAC_ENCODER_ERROR FDKaacEnc_QCOutInit(QC_OUT *phQC[(1)], const INT nSubFrames,
                                      const CHANNEL_MAPPING *cm);

AAC_ENCODER_ERROR FDKaacEnc_QCNew(QC_STATE **phQC, INT nElements,
                                  UCHAR *dynamic_RAM);

AAC_ENCODER_ERROR FDKaacEnc_QCInit(QC_STATE *hQC, struct QC_INIT *init,
                                   const ULONG initFlags);

AAC_ENCODER_ERROR FDKaacEnc_QCMainPrepare(
    ELEMENT_INFO *elInfo, ATS_ELEMENT *RESTRICT adjThrStateElement,
    PSY_OUT_ELEMENT *RESTRICT psyOutElement,
    QC_OUT_ELEMENT *RESTRICT qcOutElement,
    AUDIO_OBJECT_TYPE aot, UINT syntaxFlags, SCHAR epConfig);

AAC_ENCODER_ERROR FDKaacEnc_QCMain(QC_STATE *RESTRICT hQC, PSY_OUT **psyOut,
                                   QC_OUT **qcOut, INT avgTotalBits,
                                   CHANNEL_MAPPING *cm, AUDIO_OBJECT_TYPE aot,
                                   UINT syntaxFlags, SCHAR epConfig);

AAC_ENCODER_ERROR FDKaacEnc_updateFillBits(CHANNEL_MAPPING *cm,
                                           QC_STATE *qcKernel,
                                           ELEMENT_BITS *RESTRICT elBits[((8))],
                                           QC_OUT **qcOut);

void FDKaacEnc_updateBitres(CHANNEL_MAPPING *cm, QC_STATE *qcKernel,
                            QC_OUT **qcOut);

AAC_ENCODER_ERROR FDKaacEnc_FinalizeBitConsumption(
    CHANNEL_MAPPING *cm, QC_STATE *hQC, QC_OUT *qcOut,
    QC_OUT_ELEMENT **qcElement, HANDLE_TRANSPORTENC hTpEnc,
    AUDIO_OBJECT_TYPE aot, UINT syntaxFlags, SCHAR epConfig);

AAC_ENCODER_ERROR FDKaacEnc_AdjustBitrate(QC_STATE *RESTRICT hQC,
                                          CHANNEL_MAPPING *RESTRICT cm,
                                          INT *avgTotalBits, INT bitRate,
                                          INT sampleRate, INT granuleLength);

void FDKaacEnc_QCClose(QC_STATE **phQCstate, QC_OUT **phQC);

#endif /* QC_MAIN_H */
