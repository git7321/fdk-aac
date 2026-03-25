/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Threshold compensation

*******************************************************************************/

#ifndef ADJ_THR_H
#define ADJ_THR_H

#include "common_fix.h"
#include "adj_thr_data.h"
#include "qc_data.h"
#include "line_pe.h"
#include "interface.h"

/*****************************************************************************
    functionname: FDKaacEnc_peCalculation
    description:
*****************************************************************************/
void FDKaacEnc_peCalculation(PE_DATA *const peData,
                             const PSY_OUT_CHANNEL *const psyOutChannel[(2)],
                             QC_OUT_CHANNEL *const qcOutChannel[(2)],
                             const struct TOOLSINFO *const toolsInfo,
                             ATS_ELEMENT *const adjThrStateElement,
                             const INT nChannels);

/*****************************************************************************
functionname: FDKaacEnc_AdjThrNew
description:  allocate ADJ_THR_STATE
*****************************************************************************/
INT FDKaacEnc_AdjThrNew(ADJ_THR_STATE **phAdjThr, INT nElements);

/*****************************************************************************
functionname: FDKaacEnc_AdjThrInit
description:  initialize ADJ_THR_STATE
*****************************************************************************/
void FDKaacEnc_AdjThrInit(
    ADJ_THR_STATE *const hAdjThr, const INT meanPe, const INT invQuant,
    const CHANNEL_MAPPING *const channelMapping, const INT sampleRate,
    const INT totalBitrate, const INT isLowDelay,
    const AACENC_BITRES_MODE bitResMode, const INT dZoneQuantEnable,
    const INT bitDistributionMode, const FIXP_DBL vbrQualFactor);

/*****************************************************************************
functionname: FDKaacEnc_DistributeBits
description:
*****************************************************************************/
void FDKaacEnc_DistributeBits(
    ADJ_THR_STATE *adjThrState, ATS_ELEMENT *AdjThrStateElement,
    PSY_OUT_CHANNEL *psyOutChannel[(2)], PE_DATA *peData, INT *grantedPe,
    INT *grantedPeCorr, const INT nChannels, const INT commonWindow,
    const INT avgBits, const INT bitresBits, const INT maxBitresBits,
    const FIXP_DBL maxBitFac, const AACENC_BITRES_MODE bitResMode);

/*****************************************************************************
functionname: FDKaacEnc_AdjustThresholds
description:  adjust thresholds
*****************************************************************************/
void FDKaacEnc_AdjustThresholds(
    ADJ_THR_STATE *const hAdjThr, QC_OUT_ELEMENT *const qcElement[((8))],
    QC_OUT *const qcOut, const PSY_OUT_ELEMENT *const psyOutElement[((8))],
    const INT CBRbitrateMode, const CHANNEL_MAPPING *const cm);

/*****************************************************************************
functionname: FDKaacEnc_AdjThrClose
description:
*****************************************************************************/
void FDKaacEnc_AdjThrClose(ADJ_THR_STATE **hAdjThr);

#endif
