/**************************** AAC encoder library ******************************

   Author(s):   M. Lohwasser

   Description: PNS parameters depending on bitrate and bandwidth

*******************************************************************************/

#ifndef PNSPARAM_H
#define PNSPARAM_H

#include "aacenc.h"
#include "common_fix.h"
#include "psy_const.h"

#define NUM_PNSINFOTAB 4
#define PNS_TABLE_ERROR -1

#define USE_POWER_DISTRIBUTION (1 << 0)
#define USE_PSYCH_TONALITY (1 << 1)
#define USE_TNS_GAIN_THR (1 << 2)
#define USE_TNS_PNS (1 << 3)
#define JUST_LONG_WINDOW (1 << 4)
#define IS_LOW_COMPLEXITY (1 << 5)

typedef struct {
  short startSfb;

  USHORT detectionAlgorithmFlags;

  FIXP_DBL refPower;
  FIXP_DBL refTonality;
  INT tnsGainThreshold;
  INT tnsPNSGainThreshold;
  INT minSfbWidth;
  FIXP_SGL powDistPSDcurve[MAX_GROUPED_SFB];
  FIXP_SGL gapFillThr;
} NOISEPARAMS;

int FDKaacEnc_lookUpPnsUse(int bitRate, int sampleRate, int numChan,
                           const int isLC);

AAC_ENCODER_ERROR FDKaacEnc_GetPnsParam(NOISEPARAMS *np, INT bitRate,
                                        INT sampleRate, INT sfbCnt,
                                        const INT *sfbOffset, INT *usePns,
                                        INT numChan, const INT isLC);

#endif /* PNSPARAM_H */
