/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef SBR_H
#define SBR_H

#include "fram_gen.h"
#include "bit_sbr.h"
#include "tran_det.h"
#include "code_env.h"
#include "env_est.h"
#include "cmondata.h"
#include "qmf.h"
#include "resampler.h"
#include "ton_corr.h"

#define MAX_DELAY_FRAMES 2

typedef enum { SBRENC_DS_NONE, SBRENC_DS_TIME, SBRENC_DS_QMF } SBRENC_DS_TYPE;

typedef struct SBR_CHANNEL {
  struct ENV_CHANNEL hEnvChannel;
  DOWNSAMPLER downSampler;

} SBR_CHANNEL;
typedef SBR_CHANNEL* HANDLE_SBR_CHANNEL;

typedef struct SBR_ELEMENT {
  HANDLE_SBR_CHANNEL sbrChannel[2];
  QMF_FILTER_BANK* hQmfAnalysis[2];
  SBR_CONFIG_DATA sbrConfigData;
  SBR_HEADER_DATA sbrHeaderData;
  SBR_BITSTREAM_DATA sbrBitstreamData;
  COMMON_DATA CmonData;
  INT dynXOverFreqDelay[5];
  SBR_ELEMENT_INFO elInfo;

  UCHAR payloadDelayLine[1 + MAX_DELAY_FRAMES][MAX_PAYLOAD_SIZE];
  UINT payloadDelayLineSize[1 + MAX_DELAY_FRAMES];

} SBR_ELEMENT, *HANDLE_SBR_ELEMENT;

typedef struct SBR_ENCODER {
  HANDLE_SBR_ELEMENT sbrElement[(8)];
  HANDLE_SBR_CHANNEL pSbrChannel[(8)];
  QMF_FILTER_BANK QmfAnalysis[(8)];
  DOWNSAMPLER lfeDownSampler;
  int lfeChIdx;
  int noElements;
  int nChannels;
  int frameSize;
  int bufferOffset;
  int downsampledOffset;
  int downmixSize;
  INT downSampleFactor;
  SBRENC_DS_TYPE
  downsamplingMethod;
  int nBitstrDelay;
  int sbrDecDelay;
  INT estimateBitrate;
  INT inputDataDelay;

  UCHAR* dynamicRam;
  UCHAR* pSBRdynamic_RAM;

  HANDLE_PARAMETRIC_STEREO hParametricStereo;
  QMF_FILTER_BANK qmfSynthesisPS;

  INT maxElements;
  INT maxChannels;
  INT supportPS;
} SBR_ENCODER;

#endif /* SBR_H */
