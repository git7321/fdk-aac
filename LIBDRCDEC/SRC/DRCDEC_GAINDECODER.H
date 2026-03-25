/************************* MPEG-D DRC decoder library **************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef DRCDEC_GAINDECODER_H
#define DRCDEC_GAINDECODER_H

#include "drcDecoder.h"

#define NUM_LNB_FRAMES \
  5

#define SUBBAND_NUM_BANDS_QMF64 64
#define SUBBAND_DOWNSAMPLING_FACTOR_QMF64 64
#define SUBBAND_ANALYSIS_DELAY_QMF64 320

#define SUBBAND_NUM_BANDS_QMF71 71
#define SUBBAND_DOWNSAMPLING_FACTOR_QMF71 64
#define SUBBAND_ANALYSIS_DELAY_QMF71 320 + 384

#define SUBBAND_NUM_BANDS_STFT256 256
#define SUBBAND_DOWNSAMPLING_FACTOR_STFT256 256
#define SUBBAND_ANALYSIS_DELAY_STFT256 256

typedef enum {
  GAIN_DEC_DRC1,
  GAIN_DEC_DRC1_DRC2,
  GAIN_DEC_DRC2,
  GAIN_DEC_DRC3,
  GAIN_DEC_DRC2_DRC3
} GAIN_DEC_LOCATION;

typedef enum { GAIN_DEC_FRAME_SIZE, GAIN_DEC_SAMPLE_RATE } GAIN_DEC_PARAM;

typedef struct {
  FIXP_DBL gainLin;
  SHORT time;
} NODE_LIN;

typedef struct {
  GAIN_INTERPOLATION_TYPE gainInterpolationType;
  int nNodes[NUM_LNB_FRAMES];
  NODE_LIN linearNode[NUM_LNB_FRAMES][16];
} LINEAR_NODE_BUFFER;

typedef struct {
  int lnbPointer;
  LINEAR_NODE_BUFFER linearNodeBuffer[12];
  LINEAR_NODE_BUFFER dummyLnb;
  FIXP_DBL channelGain[8][NUM_LNB_FRAMES];
} DRC_GAIN_BUFFERS;

typedef struct {
  int activeDrcOffset;
  DRC_INSTRUCTIONS_UNI_DRC* pInst;
  DRC_COEFFICIENTS_UNI_DRC* pCoef;

  DUCKING_MODIFICATION duckingModificationForChannelGroup[8];
  SCHAR channelGroupForChannel[8];

  UCHAR bandCountForChannelGroup[8];
  UCHAR gainElementForGroup[8];
  UCHAR channelGroupIsParametricDrc[8];
  UCHAR gainElementCount;
  int lnbIndexForChannel[8][NUM_LNB_FRAMES];
  int subbandGainsReady;
} ACTIVE_DRC;

typedef struct {
  int deltaTminDefault;
  INT frameSize;
  FIXP_DBL loudnessNormalisationGainDb;
  DELAY_MODE delayMode;

  int nActiveDrcs;
  ACTIVE_DRC activeDrc[MAX_ACTIVE_DRCS];
  int multiBandActiveDrcIndex;
  int channelGainActiveDrcIndex;
  FIXP_DBL channelGain[8];

  DRC_GAIN_BUFFERS drcGainBuffers;
  FIXP_DBL subbandGains[12][4 * 1024 / 256];
  FIXP_DBL dummySubbandGains[4 * 1024 / 256];

  int status;
  int timeDomainSupported;
  SUBBAND_DOMAIN_MODE subbandDomainSupported;
} DRC_GAIN_DECODER, *HANDLE_DRC_GAIN_DECODER;

DRC_ERROR
drcDec_GainDecoder_Open(HANDLE_DRC_GAIN_DECODER* phGainDec);

DRC_ERROR
drcDec_GainDecoder_Init(HANDLE_DRC_GAIN_DECODER hGainDec);

DRC_ERROR
drcDec_GainDecoder_SetParam(HANDLE_DRC_GAIN_DECODER hGainDec,
                            const GAIN_DEC_PARAM paramType,
                            const int paramValue);

DRC_ERROR
drcDec_GainDecoder_SetCodecDependentParameters(
    HANDLE_DRC_GAIN_DECODER hGainDec, const DELAY_MODE delayMode,
    const int timeDomainSupported,
    const SUBBAND_DOMAIN_MODE subbandDomainSupported);

DRC_ERROR
drcDec_GainDecoder_Config(HANDLE_DRC_GAIN_DECODER hGainDec,
                          HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
                          const UCHAR numSelectedDrcSets,
                          const SCHAR* selectedDrcSetIds,
                          const UCHAR* selectedDownmixIds);

DRC_ERROR
drcDec_GainDecoder_Close(HANDLE_DRC_GAIN_DECODER* phGainDec);

DRC_ERROR
drcDec_GainDecoder_Preprocess(HANDLE_DRC_GAIN_DECODER hGainDec,
                              HANDLE_UNI_DRC_GAIN hUniDrcGain,
                              const FIXP_DBL loudnessNormalizationGainDb,
                              const FIXP_SGL boost, const FIXP_SGL compress);

DRC_ERROR
drcDec_GainDecoder_ProcessTimeDomain(
    HANDLE_DRC_GAIN_DECODER hGainDec, const int delaySamples,
    const GAIN_DEC_LOCATION drcLocation, const int channelOffset,
    const int drcChannelOffset, const int numChannelsProcessed,
    const int timeDataChannelOffset, FIXP_DBL* audioIOBuffer);

DRC_ERROR
drcDec_GainDecoder_ProcessSubbandDomain(
    HANDLE_DRC_GAIN_DECODER hGainDec, const int delaySamples,
    GAIN_DEC_LOCATION drcLocation, const int channelOffset,
    const int drcChannelOffset, const int numChannelsProcessed,
    const int processSingleTimeslot, FIXP_DBL* audioIOBufferReal[],
    FIXP_DBL* audioIOBufferImag[]);

DRC_ERROR
drcDec_GainDecoder_Conceal(HANDLE_DRC_GAIN_DECODER hGainDec,
                           HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
                           HANDLE_UNI_DRC_GAIN hUniDrcGain);

DRC_ERROR
drcDec_GainDecoder_SetLoudnessNormalizationGainDb(
    HANDLE_DRC_GAIN_DECODER hGainDec, FIXP_DBL loudnessNormalizationGainDb);

int drcDec_GainDecoder_GetFrameSize(HANDLE_DRC_GAIN_DECODER hGainDec);

int drcDec_GainDecoder_GetDeltaTminDefault(HANDLE_DRC_GAIN_DECODER hGainDec);

void drcDec_GainDecoder_SetChannelGains(HANDLE_DRC_GAIN_DECODER hGainDec,
                                        const int numChannels,
                                        const int frameSize,
                                        const FIXP_DBL* channelGainDb,
                                        const int audioBufferChannelOffset,
                                        FIXP_DBL* audioBuffer);

#endif
