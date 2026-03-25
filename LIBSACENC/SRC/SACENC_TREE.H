/*********************** MPEG surround encoder library *************************

   Author(s):   Max Neuendorf

   Description: Encoder Library Interface
                Tree Structure for Space Encoder

*******************************************************************************/

#ifndef SACENC_TREE_H
#define SACENC_TREE_H

#include "sacenc_framewindowing.h"
#include "sacenc_lib.h"
#include "sacenc_bitstream.h"

typedef enum {
  SPACETREE_INVALID_MODE = 0,
  SPACETREE_212 = 8

} SPACETREE_MODE;

typedef struct SPACE_TREE *HANDLE_SPACE_TREE;

typedef struct {
  UCHAR nParamBands;
  UCHAR bUseCoarseQuantTtoCld;
  UCHAR bUseCoarseQuantTtoIcc;
  QUANTMODE quantMode;
  SPACETREE_MODE mode;

  UCHAR nChannelsInMax;
  UCHAR nHybridBandsMax;

} SPACE_TREE_SETUP;

typedef struct {
  UCHAR nOttBoxes;
  UCHAR nInChannels;
  UCHAR nOutChannels;

} SPACE_TREE_DESCRIPTION;

FDK_SACENC_ERROR fdk_sacenc_spaceTree_Open(HANDLE_SPACE_TREE *phSpaceTree);

FDK_SACENC_ERROR fdk_sacenc_spaceTree_Init(
    HANDLE_SPACE_TREE hST, const SPACE_TREE_SETUP *const hSetup,
    UCHAR *pParameterBand2HybridBandOffset, const INT bFrameKeep);

FDK_SACENC_ERROR fdk_sacenc_spaceTree_Apply(
    HANDLE_SPACE_TREE hST, const INT paramSet, const INT nChannelsIn,
    const INT nTimeSlots, const INT startTimeSlot, const INT nHybridBands,
    FIXP_WIN *pFrameWindowAna__FDK,
    FIXP_DPK *const *const *const pppHybrid__FDK,
    FIXP_DPK *const *const *const pppHybridIn__FDK, SPATIALFRAME *const hSTOut,
    const INT avoid_keep, INT *pEncoderInputChScale);

FDK_SACENC_ERROR fdk_sacenc_spaceTree_Close(HANDLE_SPACE_TREE *phSpaceTree);

FDK_SACENC_ERROR fdk_sacenc_spaceTree_GetDescription(
    const HANDLE_SPACE_TREE hSpaceTree,
    SPACE_TREE_DESCRIPTION *pSpaceTreeDescription);

INT fdk_sacenc_spaceTree_Hybrid2ParamBand(const INT nParamBands,
                                          const INT nHybridBand);

#endif /* SACENC_TREE_H */
