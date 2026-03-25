/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Processing

*******************************************************************************/

#ifndef SAC_PROCESS_H
#define SAC_PROCESS_H

#include "sac_dec.h"

void SpatialDecApplyPhase(spatialDec *self, FIXP_SGL alpha,
                          int lastSlotOfParamSet);

SACDEC_ERROR SpatialDecQMFAnalysis(spatialDec *self, const PCM_MPS *inData,
                                   const INT ts, const INT bypassMode,
                                   FIXP_DBL **qmfReal, FIXP_DBL **qmfImag,
                                   const int numInputChannels);

SACDEC_ERROR SpatialDecFeedQMF(spatialDec *self, FIXP_DBL **qmfInDataReal,
                               FIXP_DBL **qmfInDataImag, const INT ts,
                               const INT bypassMode, FIXP_DBL **qmfReal,
                               FIXP_DBL **qmfImag, const INT numInputChannels);

SACDEC_ERROR SpatialDecHybridAnalysis(spatialDec *self, FIXP_DBL **qmfInputReal,
                                      FIXP_DBL **qmfInputImag,
                                      FIXP_DBL **hybOutputReal,
                                      FIXP_DBL **hybOutputImag, const INT ts,
                                      const INT numInputChannels);

SACDEC_ERROR SpatialDecCreateX(spatialDec *self, FIXP_DBL **hybInputReal,
                               FIXP_DBL **hybInputImag, FIXP_DBL **pxReal,
                               FIXP_DBL **pxImag);

SACDEC_ERROR SpatialDecApplyM1_CreateW_Mode212(
    spatialDec *self, const SPATIAL_BS_FRAME *frame, FIXP_DBL **xReal,
    FIXP_DBL **xImag, FIXP_DBL **vReal, FIXP_DBL **vImag);

SACDEC_ERROR SpatialDecApplyM2(spatialDec *self, INT ps, const FIXP_SGL alpha,
                               FIXP_DBL **wReal, FIXP_DBL **wImag,
                               FIXP_DBL **hybOutputRealDry,
                               FIXP_DBL **hybOutputImagDry,
                               FIXP_DBL **hybOutputRealWet,
                               FIXP_DBL **hybOutputImagWet);

SACDEC_ERROR SpatialDecApplyM2_Mode212_ResidualsPlusPhaseCoding(
    spatialDec *self, INT ps, const FIXP_SGL alpha, FIXP_DBL **wReal,
    FIXP_DBL **wImag, FIXP_DBL **hybOutputRealDry, FIXP_DBL **hybOutputImagDry);

SACDEC_ERROR SpatialDecApplyM2_Mode212(spatialDec *self, INT ps,
                                       const FIXP_SGL alpha, FIXP_DBL **wReal,
                                       FIXP_DBL **wImag,
                                       FIXP_DBL **hybOutputRealDry,
                                       FIXP_DBL **hybOutputImagDry);

SACDEC_ERROR SpatialDecSynthesis(spatialDec *self, const INT ts,
                                 FIXP_DBL **hybOutputReal,
                                 FIXP_DBL **hybOutputImag, PCM_MPS *timeOut,
                                 const INT numInputChannels,
                                 const FDK_channelMapDescr *const mapDescr);

void SpatialDecBufferMatrices(spatialDec *self);

FIXP_DBL getChGain(spatialDec *self, UINT ch, INT *scale);

#endif
