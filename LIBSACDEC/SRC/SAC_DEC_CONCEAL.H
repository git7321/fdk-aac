/*********************** MPEG surround decoder library *************************

   Author(s):   Christian Ertel, Christian Griebel

   Description: SAC Dec error concealment

*******************************************************************************/

#ifndef SAC_DEC_CONCEAL_H
#define SAC_DEC_CONCEAL_H

#include "sac_dec_interface.h"

typedef enum {
  SAC_DEC_CONCEAL_METHOD = 0,
  SAC_DEC_CONCEAL_NUM_KEEP_FRAMES,
  SAC_DEC_CONCEAL_FADE_OUT_SLOPE_LENGTH,
  SAC_DEC_CONCEAL_FADE_IN_SLOPE_LENGTH,
  SAC_DEC_CONCEAL_NUM_RELEASE_FRAMES

} SAC_DEC_CONCEAL_PARAM;

typedef enum {
  SAC_DEC_CONCEAL_WITH_ZERO_VALUED_OUTPUT = 0,
  SAC_DEC_CONCEAL_BY_FADING_PARAMETERS = 1

} SpatialDecConcealmentMethod;

#define MPEGS_CONCEAL_DEFAULT_METHOD SAC_DEC_CONCEAL_BY_FADING_PARAMETERS
#define MPEGS_CONCEAL_DEFAULT_NUM_KEEP_FRAMES (10)
#define MPEGS_CONCEAL_DEFAULT_FADE_OUT_SLOPE_LENGTH (5)
#define MPEGS_CONCEAL_DEFAULT_FADE_IN_SLOPE_LENGTH (5)
#define MPEGS_CONCEAL_DEFAULT_NUM_RELEASE_FRAMES (3)

typedef enum {
  SpatialDecConcealState_Init = 0,
  SpatialDecConcealState_Ok,
  SpatialDecConcealState_Keep,
  SpatialDecConcealState_FadeToDefault,
  SpatialDecConcealState_Default,
  SpatialDecConcealState_FadeFromDefault

} SpatialDecConcealmentState;

typedef struct {
  SpatialDecConcealmentMethod method;

  UINT numKeepFrames;
  UINT numFadeOutFrames;
  UINT numFadeInFrames;
  UINT numReleaseFrames;

} SpatialDecConcealmentParams;

typedef struct {
  SpatialDecConcealmentParams concealParams;
  SpatialDecConcealmentState
      concealState;

  UINT cntStateFrames;
  UINT cntValidFrames;

} SpatialDecConcealmentInfo;

#define MPEGS_CONCEAL_RESET_STATE (0x01)
#define MPEGS_CONCEAL_RESET_PARAMETER (0x02)
#define MPEGS_CONCEAL_RESET_ALL (0xFF)

void SpatialDecConcealment_Init(SpatialDecConcealmentInfo *info,
                                const UINT resetFlags);

int SpatialDecConcealment_Apply(SpatialDecConcealmentInfo *info,
                                const SCHAR (*cmpIdxData)[MAX_PARAMETER_BANDS],
                                SCHAR **diffIdxData, SCHAR *idxPrev,
                                SCHAR *bsXXXDataMode, const int startBand,
                                const int stopBand, const SCHAR defaultValue,
                                const int paramType, const int numParamSets);

void SpatialDecConcealment_UpdateState(SpatialDecConcealmentInfo *info,
                                       const int frameOk);

SACDEC_ERROR SpatialDecConcealment_SetParam(SpatialDecConcealmentInfo *info,
                                            const SAC_DEC_CONCEAL_PARAM param,
                                            const INT value);

#endif /* SAC_DEC_CONCEAL_H */
