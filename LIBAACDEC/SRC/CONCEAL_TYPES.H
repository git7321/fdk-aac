/**************************** AAC decoder library ******************************

   Author(s):   Christian Griebel

   Description: Error concealment structs and types

*******************************************************************************/

#ifndef CONCEAL_TYPES_H
#define CONCEAL_TYPES_H

#include "machine_type.h"
#include "common_fix.h"

#include "rvlc_info.h"

#include "usacdec_lpc.h"

#define CONCEAL_MAX_NUM_FADE_FACTORS (32)

#define FIXP_CNCL FIXP_DBL
#define FL2FXCONST_CNCL FL2FXCONST_DBL
#define FX_DBL2FX_CNCL
#define FX_CNCL2FX_DBL
#define CNCL_FRACT_BITS DFRACT_BITS

typedef enum {
  ConcealMethodNone = -1,
  ConcealMethodMute = 0,
  ConcealMethodNoise = 1,
  ConcealMethodInter = 2,
  ConcealMethodTonal = 3

} CConcealmentMethod;

typedef enum {
  ConcealState_Ok,
  ConcealState_Single,
  ConcealState_FadeIn,
  ConcealState_Mute,
  ConcealState_FadeOut

} CConcealmentState;

typedef struct {
  FIXP_SGL fadeOutFactor[CONCEAL_MAX_NUM_FADE_FACTORS];
  FIXP_SGL fadeInFactor[CONCEAL_MAX_NUM_FADE_FACTORS];

  CConcealmentMethod method;

  int numFadeOutFrames;
  int numFadeInFrames;
  int numMuteReleaseFrames;
  FIXP_DBL comfortNoiseLevel;

} CConcealParams;

typedef enum {
  FADE_TIMEDOMAIN_TOSPECTRALMUTE = 1,
  FADE_TIMEDOMAIN_FROMSPECTRALMUTE,
  FADE_TIMEDOMAIN
} TDfadingType;

typedef struct {
  CConcealParams *pConcealParams;

  FIXP_CNCL spectralCoefficient[1024];
  SHORT specScale[8];

  INT iRandomPhase;
  INT prevFrameOk[2];
  INT cntValidFrames;
  INT cntFadeFrames;
  int winGrpOffset[2];
  int attGrpOffset[2];

  SCHAR lastRenderMode;

  UCHAR windowShape;
  BLOCK_TYPE windowSequence;
  UCHAR lastWinGrpLen;

  CConcealmentState concealState;
  CConcealmentState concealState_old;
  FIXP_DBL fade_old;
  TDfadingType lastFadingType;

  SHORT aRvlcPreviousScaleFactor[RVLC_MAX_SFB];
  UCHAR aRvlcPreviousCodebook[RVLC_MAX_SFB];
  SCHAR rvlcPreviousScaleFactorOK;
  SCHAR rvlcPreviousBlockType;

  FIXP_LPC lsf4[M_LP_FILTER_ORDER];
  FIXP_DBL last_tcx_gain;
  INT last_tcx_gain_e;
  ULONG TDNoiseSeed;
  PCM_DEC TDNoiseStates[3];
  FIXP_SGL TDNoiseCoef[3];
  FIXP_SGL TDNoiseAtt;

} CConcealmentInfo;

#endif /* #ifndef CONCEAL_TYPES_H */
