/**************************** PCM utility library ******************************

   Author(s):   Matthias Neusinger

   Description: Hard limiter for clipping prevention

*******************************************************************************/

#ifndef LIMITER_H
#define LIMITER_H

#include "common_fix.h"
#include "FDK_audio.h"

#define TDL_ATTACK_DEFAULT_MS (15)
#define TDL_RELEASE_DEFAULT_MS (50)

#ifdef __cplusplus
extern "C" {
#endif

struct TDLimiter {
  unsigned int attack;
  FIXP_DBL attackConst, releaseConst;
  unsigned int attackMs, releaseMs, maxAttackMs;
  FIXP_DBL threshold;
  unsigned int channels, maxChannels;
  UINT sampleRate, maxSampleRate;
  FIXP_DBL cor, max;
  FIXP_DBL* maxBuf;
  FIXP_DBL* delayBuf;
  unsigned int maxBufIdx, delayBufIdx;
  FIXP_DBL smoothState0;
  FIXP_DBL minGain;
  INT scaling;
};

typedef enum {
  TDLIMIT_OK = 0,
  TDLIMIT_UNKNOWN = -1,

  __error_codes_start = -100,

  TDLIMIT_INVALID_HANDLE,
  TDLIMIT_INVALID_PARAMETER,

  __error_codes_end
} TDLIMITER_ERROR;

struct TDLimiter;
typedef struct TDLimiter* TDLimiterPtr;

#define PCM_LIM LONG
#define FIXP_DBL2PCM_LIM(x) (x)
#define PCM_LIM2FIXP_DBL(x) (x)
#define PCM_LIM_BITS 32
#define FIXP_PCM_LIM FIXP_DBL

#define SAMPLE_BITS_LIM DFRACT_BITS

TDLIMITER_ERROR pcmLimiter_Reset(TDLimiterPtr limiter);

TDLIMITER_ERROR pcmLimiter_Destroy(TDLimiterPtr limiter);

unsigned int pcmLimiter_GetDelay(TDLimiterPtr limiter);

INT pcmLimiter_GetMaxGainReduction(TDLimiterPtr limiter);

TDLIMITER_ERROR pcmLimiter_SetNChannels(TDLimiterPtr limiter,
                                        unsigned int nChannels);

TDLIMITER_ERROR pcmLimiter_SetSampleRate(TDLimiterPtr limiter, UINT sampleRate);

TDLIMITER_ERROR pcmLimiter_SetAttack(TDLimiterPtr limiter,
                                     unsigned int attackMs);

TDLIMITER_ERROR pcmLimiter_SetRelease(TDLimiterPtr limiter,
                                      unsigned int releaseMs);

TDLIMITER_ERROR pcmLimiter_GetLibInfo(LIB_INFO* info);

#ifdef __cplusplus
}
#endif

TDLimiterPtr pcmLimiter_Create(unsigned int maxAttackMs, unsigned int releaseMs,
                               FIXP_DBL threshold, unsigned int maxChannels,
                               UINT maxSampleRate);

TDLIMITER_ERROR pcmLimiter_SetThreshold(TDLimiterPtr limiter,
                                        FIXP_DBL threshold);

TDLIMITER_ERROR pcmLimiter_Apply(TDLimiterPtr limiter, PCM_LIM* samplesIn,
                                 INT_PCM* samplesOut, FIXP_DBL* pGainPerSample,
                                 const INT scaling, const UINT nSamples);

#endif /* #ifndef LIMITER_H */
