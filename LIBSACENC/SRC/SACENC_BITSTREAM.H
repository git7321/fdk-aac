/*********************** MPEG surround encoder library *************************

   Author(s):

   Description: Encoder Library Interface
                Bitstream Writer

*******************************************************************************/

#ifndef SACENC_BITSTREAM_H
#define SACENC_BITSTREAM_H

#include "FDK_bitstream.h"
#include "FDK_matrixCalloc.h"
#include "sacenc_lib.h"
#include "sacenc_const.h"

#define MAX_NUM_BINS 23
#define MAX_NUM_PARAMS 2
#define MAX_NUM_OUTPUTCHANNELS SACENC_MAX_OUTPUT_CHANNELS
#define MAX_TIME_SLOTS 32

typedef enum {
  TREE_212 = 7,
  TREE_ESCAPE = 15
} TREECONFIG;

typedef enum {
  FREQ_RES_40 = 0,
  FREQ_RES_20 = 1,
  FREQ_RES_10 = 2,
  FREQ_RES_5 = 3
} FREQ;

typedef enum {
  QUANTMODE_INVALID = -1,
  QUANTMODE_FINE = 0,
  QUANTMODE_EBQ1 = 1,
  QUANTMODE_EBQ2 = 2
} QUANTMODE;

typedef enum {
  TEMPSHAPE_OFF = 0
} TEMPSHAPECONFIG;

typedef enum {
  FIXEDGAINDMX_INVALID = -1,
  FIXEDGAINDMX_0 = 0,
  FIXEDGAINDMX_1 = 1,
  FIXEDGAINDMX_2 = 2,
  FIXEDGAINDMX_3 = 3,
  FIXEDGAINDMX_4 = 4,
  FIXEDGAINDMX_5 = 5,
  FIXEDGAINDMX_6 = 6,
  FIXEDGAINDMX_7 = 7
} FIXEDGAINDMXCONFIG;

typedef enum {
  DECORR_INVALID = -1,
  DECORR_QMFSPLIT0 = 0,
  DECORR_QMFSPLIT1 = 1,
  DECORR_QMFSPLIT2 = 2
} DECORRCONFIG;

typedef enum {
  DEFAULT = 0,
  KEEP = 1,
  INTERPOLATE = 2,
  FINECOARSE = 3
} DATA_MODE;

typedef enum {
  READ_SPATIALFRAME = 0,
  WRITE_SPATIALFRAME = 1
} SPATIALFRAME_TYPE;

typedef struct {
  INT numOttBoxes;
  INT numInChan;
  INT numOutChan;

} TREEDESCRIPTION;

typedef struct {
  INT bsOttBands;

} OTTCONFIG;

typedef struct {
  INT bsSamplingFrequency;
  INT bsFrameLength;
  INT numBands;
  TREECONFIG bsTreeConfig;
  QUANTMODE bsQuantMode;
  FIXEDGAINDMXCONFIG bsFixedGainDMX;
  int bsEnvQuantMode;
  DECORRCONFIG bsDecorrConfig;
  TREEDESCRIPTION treeDescription;
  OTTCONFIG ottConfig[SACENC_MAX_NUM_BOXES];
} SPATIALSPECIFICCONFIG;

typedef struct {
  UCHAR bsFramingType;
  UCHAR numParamSets;
  UCHAR bsParamSlots[MAX_NUM_PARAMS];
} FRAMINGINFO;

typedef struct {
  SCHAR cld[SACENC_MAX_NUM_BOXES][MAX_NUM_PARAMS][MAX_NUM_BINS];
  SCHAR icc[SACENC_MAX_NUM_BOXES][MAX_NUM_PARAMS][MAX_NUM_BINS];
} OTTDATA;

typedef struct {
  UCHAR bsSmoothMode[MAX_NUM_PARAMS];
  UCHAR bsSmoothTime[MAX_NUM_PARAMS];
  UCHAR bsFreqResStride[MAX_NUM_PARAMS];
  UCHAR bsSmgData[MAX_NUM_PARAMS][MAX_NUM_BINS];
} SMGDATA;

typedef struct {
  UCHAR bsEnvShapeChannel[MAX_NUM_OUTPUTCHANNELS];
  UCHAR bsEnvShapeData[MAX_NUM_OUTPUTCHANNELS][MAX_TIME_SLOTS];
} TEMPSHAPEDATA;

typedef struct {
  UCHAR bsXXXDataMode[SACENC_MAX_NUM_BOXES][MAX_NUM_PARAMS];
  UCHAR bsDataPair[SACENC_MAX_NUM_BOXES][MAX_NUM_PARAMS];
  UCHAR bsQuantCoarseXXX[SACENC_MAX_NUM_BOXES][MAX_NUM_PARAMS];
  UCHAR bsFreqResStrideXXX[SACENC_MAX_NUM_BOXES][MAX_NUM_PARAMS];
} LOSSLESSDATA;

typedef struct {
  FRAMINGINFO framingInfo;
  UCHAR bsIndependencyFlag;
  OTTDATA ottData;
  SMGDATA smgData;
  TEMPSHAPEDATA tempShapeData;
  LOSSLESSDATA CLDLosslessData;
  LOSSLESSDATA ICCLosslessData;
  UCHAR bUseBBCues;
} SPATIALFRAME;

typedef struct BSF_INSTANCE *HANDLE_BSF_INSTANCE;

FDK_SACENC_ERROR fdk_sacenc_destroySpatialBitstreamEncoder(
    HANDLE_BSF_INSTANCE *selfPtr);
FDK_SACENC_ERROR fdk_sacenc_createSpatialBitstreamEncoder(
    HANDLE_BSF_INSTANCE *selfPtr);
FDK_SACENC_ERROR fdk_sacenc_initSpatialBitstreamEncoder(
    HANDLE_BSF_INSTANCE selfPtr);

SPATIALSPECIFICCONFIG *fdk_sacenc_getSpatialSpecificConfig(
    HANDLE_BSF_INSTANCE selfPtr);

FDK_SACENC_ERROR fdk_sacenc_writeSpatialSpecificConfig(
    SPATIALSPECIFICCONFIG *const spatialSpecificConfig,
    UCHAR *const pOutputBuffer, const INT outputBufferSize,
    INT *const pnOutputBits);

SPATIALFRAME *fdk_sacenc_getSpatialFrame(HANDLE_BSF_INSTANCE selfPtr,
                                         const SPATIALFRAME_TYPE frameType);

FDK_SACENC_ERROR fdk_sacenc_writeSpatialFrame(UCHAR *const pOutputBuffer,
                                              const INT outputBufferSize,
                                              INT *const pnOutputBits,
                                              HANDLE_BSF_INSTANCE selfPtr);
FDK_SACENC_ERROR fdk_sacenc_duplicateParameterSet(
    const SPATIALFRAME *const hFrom, const INT setFrom, SPATIALFRAME *const hTo,
    const INT setTo);

#endif /* SACENC_BITSTREAM_H */
