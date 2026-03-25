/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef ENV_EXTR_H
#define ENV_EXTR_H

#include "sbrdecoder.h"
#include "FDK_bitstream.h"
#include "lpp_tran.h"
#include "psdec.h"
#include "pvc_dec.h"

#define ENV_EXP_FRACT 0
#define EXP_BITS 6

#define MASK_M \
  (((1 << (FRACT_BITS - EXP_BITS)) - 1) \
   << EXP_BITS)
#define MASK_E \
  ((1 << EXP_BITS) - 1)

#define SIGN_EXT \
  (((SCHAR)-1) ^ \
   MASK_E)
#define ROUNDING \
  ((FIXP_SGL)( \
      1 << (EXP_BITS - 1)))
#define NRG_EXP_OFFSET \
  16
#define NOISE_EXP_OFFSET \
  38

#define ADD_HARMONICS_FLAGS_SIZE 2

typedef enum {
  HEADER_NOT_PRESENT,
  HEADER_ERROR,
  HEADER_OK,
  HEADER_RESET
} SBR_HEADER_STATUS;

typedef enum {
  SBR_NOT_INITIALIZED = 0,
  UPSAMPLING = 1,
  SBR_HEADER = 2,
  SBR_ACTIVE = 3
} SBR_SYNC_STATE;

typedef enum { COUPLING_OFF = 0, COUPLING_LEVEL, COUPLING_BAL } COUPLING_MODE;

typedef struct {
  UCHAR nSfb[2];
  UCHAR nNfb;
  UCHAR numMaster;
  UCHAR lowSubband;
  UCHAR highSubband;
  UCHAR ov_highSubband;
  UCHAR limiterBandTable[MAX_NUM_LIMITERS + 1];
  UCHAR noLimiterBands;
  UCHAR nInvfBands;
  UCHAR *freqBandTable[2];
  UCHAR freqBandTableLo[MAX_FREQ_COEFFS / 2 + 1];
  UCHAR freqBandTableHi[MAX_FREQ_COEFFS + 1];
  UCHAR freqBandTableNoise[MAX_NOISE_COEFFS + 1];
  UCHAR v_k_master[MAX_FREQ_COEFFS + 1];
} FREQ_BAND_DATA;

typedef FREQ_BAND_DATA *HANDLE_FREQ_BAND_DATA;

#define SBRDEC_ELD_GRID 1
#define SBRDEC_SYNTAX_SCAL 2
#define SBRDEC_SYNTAX_USAC 4
#define SBRDEC_SYNTAX_RSVD50 8
#define SBRDEC_USAC_INDEP \
  16
#define SBRDEC_LOW_POWER \
  32
#define SBRDEC_PS_DECODED \
  64
#define SBRDEC_QUAD_RATE \
  128
#define SBRDEC_USAC_HARMONICSBR \
  256
#define SBRDEC_LD_MPS_QMF \
  512
#define SBRDEC_USAC_ITES \
  1024
#define SBRDEC_SYNTAX_DRM \
  2048
#define SBRDEC_ELD_DOWNSCALE \
  4096
#define SBRDEC_DOWNSAMPLE \
  8192
#define SBRDEC_FLUSH 16384
#define SBRDEC_FORCE_RESET \
  32768
#define SBRDEC_SKIP_QMF_ANA \
  (1 << 21)
#define SBRDEC_SKIP_QMF_SYN \
  (1 << 22)

#define SBRDEC_HDR_STAT_RESET 1
#define SBRDEC_HDR_STAT_UPDATE 2

typedef struct {
  UCHAR ampResolution;
  UCHAR
  xover_band;
  UCHAR sbr_preprocessing;
  UCHAR pvc_mode;
} SBR_HEADER_DATA_BS_INFO;

typedef struct {
  UCHAR startFreq;
  UCHAR stopFreq;
  UCHAR freqScale;
  UCHAR alterScale;
  UCHAR noise_bands;

  UCHAR limiterBands;
  UCHAR limiterGains;
  UCHAR interpolFreq;
  UCHAR smoothingLength;
} SBR_HEADER_DATA_BS;

typedef struct {
  SBR_SYNC_STATE
  syncState;

  UCHAR status;
  UCHAR
  frameErrorFlag;
  UCHAR numberTimeSlots;
  UCHAR numberOfAnalysisBands;
  UCHAR timeStep;
  UINT sbrProcSmplRate;

  SBR_HEADER_DATA_BS bs_data;
  SBR_HEADER_DATA_BS bs_dflt;
  SBR_HEADER_DATA_BS_INFO bs_info;

  FREQ_BAND_DATA freqBandData;
  UCHAR pvcIDprev;
} SBR_HEADER_DATA;

typedef SBR_HEADER_DATA *HANDLE_SBR_HEADER_DATA;

typedef struct {
  UCHAR frameClass;
  UCHAR nEnvelopes;
  UCHAR borders[MAX_ENVELOPES + 1];
  UCHAR freqRes[MAX_ENVELOPES];
  SCHAR tranEnv;
  UCHAR nNoiseEnvelopes;
  UCHAR
  bordersNoise[MAX_NOISE_ENVELOPES + 1];
  UCHAR pvcBorders[MAX_PVC_ENVELOPES + 1];
  UCHAR noisePosition;
  UCHAR varLength;
} FRAME_INFO;

typedef struct {
  FIXP_SGL sfb_nrg_prev[MAX_FREQ_COEFFS];
  FIXP_SGL
  prevNoiseLevel[MAX_NOISE_COEFFS];
  COUPLING_MODE coupling;
  INVF_MODE sbr_invf_mode[MAX_INVF_BANDS];
  UCHAR ampRes;
  UCHAR stopPos;
  UCHAR frameErrorFlag;
  UCHAR prevSbrPitchInBins;
  FRAME_INFO prevFrameInfo;
} SBR_PREV_FRAME_DATA;

typedef SBR_PREV_FRAME_DATA *HANDLE_SBR_PREV_FRAME_DATA;

typedef struct {
  int nScaleFactors;

  FRAME_INFO frameInfo;
  UCHAR domain_vec[MAX_ENVELOPES];
  UCHAR domain_vec_noise
      [MAX_NOISE_ENVELOPES];

  INVF_MODE
  sbr_invf_mode[MAX_INVF_BANDS];
  COUPLING_MODE coupling;
  int ampResolutionCurrentFrame;

  ULONG addHarmonics[ADD_HARMONICS_FLAGS_SIZE];

  FIXP_SGL iEnvelope[MAX_NUM_ENVELOPE_VALUES];
  FIXP_SGL sbrNoiseFloorLevel[MAX_NUM_NOISE_VALUES];
  UCHAR iTESactive;
  UCHAR
  interTempShapeMode[MAX_ENVELOPES];
  UCHAR pvcID[PVC_NTIMESLOT];
  UCHAR ns;
  UCHAR sinusoidal_position;

  UCHAR sbrPatchingMode;
  UCHAR sbrOversamplingFlag;
  UCHAR sbrPitchInBins;
} SBR_FRAME_DATA;

typedef SBR_FRAME_DATA *HANDLE_SBR_FRAME_DATA;

UINT sbrdec_mapToStdSampleRate(UINT fs,
                               UINT isUsac);

void initSbrPrevFrameData(HANDLE_SBR_PREV_FRAME_DATA h_prev_data,
                          int timeSlots);

int sbrGetChannelElement(HANDLE_SBR_HEADER_DATA hHeaderData,
                         HANDLE_SBR_FRAME_DATA hFrameDataLeft,
                         HANDLE_SBR_FRAME_DATA hFrameDataRight,
                         HANDLE_SBR_PREV_FRAME_DATA hFrameDataLeftPrev,
                         UCHAR pvc_mode_last, HANDLE_FDK_BITSTREAM hBitBuf,
                         HANDLE_PS_DEC hParametricStereoDec, const UINT flags,
                         const int overlap);

SBR_HEADER_STATUS
sbrGetHeaderData(HANDLE_SBR_HEADER_DATA headerData,
                 HANDLE_FDK_BITSTREAM hBitBuf, const UINT flags,
                 const int fIsSbrData, const UCHAR configMode);

SBR_ERROR
initHeaderData(HANDLE_SBR_HEADER_DATA hHeaderData, const int sampleRateIn,
               const int sampleRateOut, const INT downscaleFactor,
               const int samplesPerFrame, const UINT flags,
               const int setDefaultHdr);
#endif

#define SCALE2EXP(s) (15 - (s))
#define EXP2SCALE(e) (15 - (e))
