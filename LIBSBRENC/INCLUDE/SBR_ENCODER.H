/**************************** SBR encoder library ******************************

   Author(s):

   Description: SBR encoder top level processing prototype

*******************************************************************************/

#ifndef SBR_ENCODER_H
#define SBR_ENCODER_H

#include "common_fix.h"
#include "FDK_audio.h"
#include "FDK_bitstream.h"

#define MAX_TRANS_FAC 8
#define MAX_CODEC_FRAME_RATIO 2
#define MAX_PAYLOAD_SIZE 256

typedef enum codecType {
  CODEC_AAC = 0,
  CODEC_AACLD = 1,
  CODEC_UNSPECIFIED = 99
} CODEC_TYPE;

typedef struct {
  INT bitRate;
  INT nChannels;
  INT sampleFreq;
  INT transFac;
  INT standardBitrate;
} CODEC_PARAM;

typedef enum {
  SBR_MONO,
  SBR_LEFT_RIGHT,
  SBR_COUPLING,
  SBR_SWITCH_LRC
} SBR_STEREO_MODE;

enum {
  SBR_SYNTAX_LOW_DELAY = 0x0001,
  SBR_SYNTAX_SCALABLE = 0x0002,
  SBR_SYNTAX_CRC = 0x0004,
  SBR_SYNTAX_DRM_CRC = 0x0008,
  SBR_SYNTAX_ELD_REDUCED_DELAY = 0x0010
};

typedef enum { FREQ_RES_LOW = 0, FREQ_RES_HIGH } FREQ_RES;

typedef struct {
  CODEC_TYPE coreCoder;
  UINT bitrateFrom;
  UINT bitrateTo;

  UINT sampleRate;
  UCHAR numChannels;

  UCHAR startFreq;
  UCHAR startFreqSpeech;
  UCHAR stopFreq;
  UCHAR stopFreqSpeech;

  UCHAR numNoiseBands;
  UCHAR noiseFloorOffset;
  SCHAR noiseMaxLevel;
  SBR_STEREO_MODE stereoMode;
  UCHAR freqScale;
} sbrTuningTable_t;

typedef struct sbrConfiguration {
  CODEC_PARAM
  codecSettings;
  INT SendHeaderDataTime;
  INT useWaveCoding;
  INT crcSbr;
  INT dynBwSupported;
  INT parametricCoding;
  INT downSampleFactor;
  FREQ_RES freq_res_fixfix[2];
  UCHAR fResTransIsLow;

  INT tran_thr;
  INT noiseFloorOffset;
  UINT useSpeechConfig;

  INT sbrFrameSize;
  INT sbr_data_extra;
  INT amp_res;
  INT ana_max_level;
  INT tran_fc;
  INT tran_det_mode;
  INT spread;
  INT stat;
  INT e;
  SBR_STEREO_MODE stereoMode;
  INT deltaTAcrossFrames;
  FIXP_DBL dF_edge_1stEnv;
  FIXP_DBL dF_edge_incr;
  INT sbr_invf_mode;
  INT sbr_xpos_mode;
  INT sbr_xpos_ctrl;
  INT sbr_xpos_level;
  INT startFreq;
  INT stopFreq;
  INT useSaPan;
  INT dynBwEnabled;
  INT bParametricStereo;

  UCHAR freqScale;
  INT alterScale;
  INT sbr_noise_bands;

  INT sbr_limiter_bands;
  INT sbr_limiter_gains;
  INT sbr_interpol_freq;
  INT sbr_smoothing_length;
  UCHAR init_amp_res_FF;
  FIXP_DBL threshold_AmpRes_FF_m;
  SCHAR threshold_AmpRes_FF_e;
} sbrConfiguration, *sbrConfigurationPtr;

typedef struct SBR_CONFIG_DATA {
  UINT sbrSyntaxFlags;
  INT nChannels;

  INT nSfb[2];
  INT num_Master;
  INT sampleFreq;
  INT frameSize;
  INT xOverFreq;
  INT dynXOverFreq;

  INT noQmfBands;
  INT noQmfSlots;

  UCHAR *freqBandTable[2];
  UCHAR *v_k_master;

  SBR_STEREO_MODE stereoMode;
  INT noEnvChannels;

  INT useWaveCoding;
  INT useParametricCoding;
  INT xposCtrlSwitch;
  INT switchTransposers;
  UCHAR initAmpResFF;
  FIXP_DBL thresholdAmpResFF_m;
  SCHAR thresholdAmpResFF_e;
} SBR_CONFIG_DATA, *HANDLE_SBR_CONFIG_DATA;

typedef struct {
  MP4_ELEMENT_ID elType;
  INT bitRate;
  int instanceTag;
  UCHAR fParametricStereo;
  UCHAR fDualMono;
  UCHAR nChannelsInEl;
  UCHAR ChannelIndex[2];
} SBR_ELEMENT_INFO;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SBR_ENCODER *HANDLE_SBR_ENCODER;

INT sbrEncoder_GetInBufferSize(int noChannels);
INT sbrEncoder_Open(HANDLE_SBR_ENCODER *phSbrEncoder, INT nElements,
                    INT nChannels, INT supportPS);

UINT sbrEncoder_LimitBitRate(UINT bitRate, UINT numChannels,
                             UINT coreSampleRate, AUDIO_OBJECT_TYPE aot);
UINT sbrEncoder_IsSingleRatePossible(AUDIO_OBJECT_TYPE aot);

INT sbrEncoder_Init(HANDLE_SBR_ENCODER hSbrEncoder,
                    SBR_ELEMENT_INFO elInfo[(8)], int noElements,
                    INT_PCM *inputBuffer, UINT inputBufferBufSize,
                    INT *coreBandwidth, INT *inputBufferOffset,
                    INT *numChannels, const UINT syntaxFlags, INT *sampleRate,
                    UINT *downSampleFactor, INT *frameLength,
                    AUDIO_OBJECT_TYPE aot, int *delay, int transformFactor,
                    const int headerPeriod, ULONG statesInitFlag);
INT sbrEncoder_UpdateBuffers(HANDLE_SBR_ENCODER hEnvEnc, INT_PCM *timeBuffer,
                             UINT timeBufferBufSIze);

void sbrEncoder_Close(HANDLE_SBR_ENCODER *phEbrEncoder);

INT sbrEncoder_EncodeFrame(HANDLE_SBR_ENCODER hEnvEncoder, INT_PCM *samples,
                           UINT samplesBufSize, UINT sbrDataBits[(8)],
                           UCHAR sbrData[(8)][MAX_PAYLOAD_SIZE]);

void sbrEncoder_GetHeader(HANDLE_SBR_ENCODER sbrEncoder,
                          HANDLE_FDK_BITSTREAM hBs, INT element_index,
                          int fSendHeaders);

INT sbrEncoder_SendHeader(HANDLE_SBR_ENCODER hSbrEncoder);
INT sbrEncoder_ContainsHeader(HANDLE_SBR_ENCODER hSbrEncoder);
INT sbrEncoder_GetHeaderDelay(HANDLE_SBR_ENCODER hSbrEncoder);
INT sbrEncoder_GetBsDelay(HANDLE_SBR_ENCODER hSbrEncoder);
INT sbrEncoder_SAPPrepare(HANDLE_SBR_ENCODER hSbrEncoder);
INT sbrEncoder_GetEstimateBitrate(HANDLE_SBR_ENCODER hSbrEncoder);
INT sbrEncoder_GetInputDataDelay(HANDLE_SBR_ENCODER hSbrEncoder);
INT sbrEncoder_GetSbrDecDelay(HANDLE_SBR_ENCODER hSbrEncoder);
INT sbrEncoder_GetLibInfo(LIB_INFO *info);

void sbrPrintRAM(void);
void sbrPrintROM(void);

#ifdef __cplusplus
}
#endif

#endif /* ifndef __SBR_MAIN_H */
