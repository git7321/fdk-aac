/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Dec bitstream decoder

*******************************************************************************/

#ifndef SAC_BITDEC_H
#define SAC_BITDEC_H

#include "sac_dec.h"

typedef struct {
  SCHAR numInputChannels;
  SCHAR numOutputChannels;
  SCHAR numOttBoxes;
  SCHAR numTttBoxes;
  SCHAR ottModeLfe[MAX_NUM_OTT];
} TREEPROPERTIES;

enum { TREE_212 = 7, TREE_DUMMY = 255 };

enum { QUANT_FINE = 0, QUANT_EBQ1 = 1, QUANT_EBQ2 = 2 };

SACDEC_ERROR SpatialDecParseSpecificConfigHeader(
    HANDLE_FDK_BITSTREAM bitstream,
    SPATIAL_SPECIFIC_CONFIG *pSpatialSpecificConfig,
    AUDIO_OBJECT_TYPE coreCodec, SPATIAL_DEC_UPMIX_TYPE upmixType);

SACDEC_ERROR SpatialDecParseMps212Config(
    HANDLE_FDK_BITSTREAM bitstream,
    SPATIAL_SPECIFIC_CONFIG *pSpatialSpecificConfig, int samplingRate,
    AUDIO_OBJECT_TYPE coreCodec, INT stereoConfigIndex,
    INT coreSbrFrameLengthIndex);

SACDEC_ERROR SpatialDecParseSpecificConfig(
    HANDLE_FDK_BITSTREAM bitstream,
    SPATIAL_SPECIFIC_CONFIG *pSpatialSpecificConfig, int sacHeaderLen,
    AUDIO_OBJECT_TYPE coreCodec);

int SpatialDecDefaultSpecificConfig(
    SPATIAL_SPECIFIC_CONFIG *pSpatialSpecificConfig,
    AUDIO_OBJECT_TYPE coreCodec, int samplingFreq, int nTimeSlots,
    int sacDecoderLevel, int isBlind, int coreChannels);

SACDEC_ERROR SpatialDecCreateBsFrame(SPATIAL_BS_FRAME *bsFrame,
                                     BS_LL_STATE *llState);

void SpatialDecCloseBsFrame(SPATIAL_BS_FRAME *bsFrame);

SACDEC_ERROR SpatialDecParseFrameData(
    spatialDec *self, SPATIAL_BS_FRAME *frame, HANDLE_FDK_BITSTREAM bitstream,
    const SPATIAL_SPECIFIC_CONFIG *pSpatialSpecificConfig, UPMIXTYPE upmixType,
    int fGlobalIndependencyFlag);

SACDEC_ERROR SpatialDecDecodeFrame(spatialDec *self, SPATIAL_BS_FRAME *frame);

SACDEC_ERROR SpatialDecDecodeHeader(
    spatialDec *self, SPATIAL_SPECIFIC_CONFIG *pSpatialSpecificConfig);

#endif
