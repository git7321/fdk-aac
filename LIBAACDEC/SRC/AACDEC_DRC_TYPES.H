/**************************** AAC decoder library ******************************

   Author(s):   Christian Griebel

   Description: Dynamic range control (DRC) global data types

*******************************************************************************/

#ifndef AACDEC_DRC_TYPES_H
#define AACDEC_DRC_TYPES_H

#include "common_fix.h"

#define MAX_DRC_THREADS \
  ((8) + 1)
#define MAX_DRC_BANDS (16)

typedef enum {
  UNKNOWN_PAYLOAD = 0,
  MPEG_DRC_EXT_DATA = 1,
  DVB_DRC_ANC_DATA = 2

} AACDEC_DRC_PAYLOAD_TYPE;

typedef enum {
  DISABLED_PARAMETER_HANDLING = -1,
  ENABLED_PARAMETER_HANDLING =
      0,
  DRC_PRESENTATION_MODE_1 = 1,
  DRC_PRESENTATION_MODE_2 = 2

} AACDEC_DRC_PARAMETER_HANDLING;

typedef struct {
  UINT expiryCount;
  UINT numBands;
  USHORT bandTop[MAX_DRC_BANDS];
  SHORT drcInterpolationScheme;
  UCHAR drcValue[MAX_DRC_BANDS];
  SCHAR drcDataType;

} CDrcChannelData;

typedef struct {
  UINT excludedChnsMask;
  SCHAR progRefLevel;
  SCHAR presMode;
  SCHAR pceInstanceTag;

  CDrcChannelData channelData;

} CDrcPayload;

typedef struct {
  FIXP_DBL usrCut;
  FIXP_DBL usrBoost;
  UCHAR usrApplyHeavyCompression;

  FIXP_DBL cut;
  FIXP_DBL boost;
  SCHAR targetRefLevel;
  UCHAR applyHeavyCompression;

  UINT expiryFrame;
  UCHAR bsDelayEnable;

  AACDEC_DRC_PARAMETER_HANDLING defaultPresentationMode;
  UCHAR encoderTargetLevel;

} CDrcParams;

typedef struct {
  CDrcParams
      params;

  UCHAR enable;
  UCHAR digitalNorm;

  UCHAR update;
  INT numOutChannels;
  INT prevAacNumChannels;

  USHORT numPayloads;
  USHORT
  numThreads;
  SCHAR progRefLevel;
  UCHAR progRefLevelPresent;

  UINT prlExpiryCount;

  SCHAR presMode;
  UCHAR dvbAncDataAvailable;
  UINT dvbAncDataPosition;
  UINT drcPayloadPosition[MAX_DRC_THREADS];
  UCHAR applyExtGain;

  FIXP_DBL additionalGainPrev;
  FIXP_DBL additionalGainFilterState;
  FIXP_DBL additionalGainFilterState1;

} CDrcInfo;

typedef CDrcInfo *HANDLE_AAC_DRC;

#endif /* AACDEC_DRC_TYPES_H */
