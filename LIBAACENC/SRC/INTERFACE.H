/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Interface psychoaccoustic/quantizer

*******************************************************************************/

#ifndef INTERFACE_H
#define INTERFACE_H

#include "common_fix.h"
#include "FDK_audio.h"

#include "psy_data.h"
#include "aacenc_tns.h"

enum { MS_NONE = 0, MS_SOME = 1, MS_ALL = 2 };

enum { MS_ON = 1 };

struct TOOLSINFO {
  INT msDigest;
  INT msMask[MAX_GROUPED_SFB];
};

typedef struct {
  INT sfbCnt;
  INT sfbPerGroup;
  INT maxSfbPerGroup;
  INT lastWindowSequence;
  INT windowShape;
  INT groupingMask;
  INT sfbOffsets[MAX_GROUPED_SFB + 1];

  INT mdctScale;
  INT groupLen[MAX_NO_OF_GROUPS];

  TNS_INFO tnsInfo;
  INT noiseNrg[MAX_GROUPED_SFB];
  INT isBook[MAX_GROUPED_SFB];
  INT isScale[MAX_GROUPED_SFB];

  FIXP_DBL *mdctSpectrum;
  FIXP_DBL *sfbEnergy;
  FIXP_DBL *sfbSpreadEnergy;
  FIXP_DBL *sfbThresholdLdData;
  FIXP_DBL *sfbMinSnrLdData;
  FIXP_DBL *sfbEnergyLdData;

} PSY_OUT_CHANNEL;

typedef struct {
  PSY_OUT_CHANNEL *psyOutChannel[(2)];
  INT commonWindow;
  struct TOOLSINFO toolsInfo;

} PSY_OUT_ELEMENT;

typedef struct {
  PSY_OUT_ELEMENT *psyOutElement[((8))];
  PSY_OUT_CHANNEL *pPsyOutChannels[(8)];

} PSY_OUT;

inline int isLowDelay(AUDIO_OBJECT_TYPE aot) {
  return (aot == AOT_ER_AAC_LD || aot == AOT_ER_AAC_ELD);
}

#endif /* INTERFACE_H */
