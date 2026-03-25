/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef TRAN_DET_H
#define TRAN_DET_H

#include "sbr_encoder.h"
#include "sbr_def.h"

typedef struct {
  FIXP_DBL transients[32 + (32 / 2)];
  FIXP_DBL thresholds[64];
  FIXP_DBL tran_thr;    /* Master threshold for transient signals */
  FIXP_DBL split_thr_m; /* Threshold for splitting FIXFIX-frames into 2 env */
  INT split_thr_e;      /* Scale for splitting threshold */
  FIXP_DBL prevLowBandEnergy;  /* Energy of low band */
  FIXP_DBL prevHighBandEnergy; /* Energy of high band */
  INT tran_fc;                 /* Number of lowband subbands to discard  */
  INT no_cols;
  INT no_rows;
  INT mode;

  int frameShift;
  int tran_off; /* Offset for reading energy values. */
} SBR_TRANSIENT_DETECTOR;

typedef SBR_TRANSIENT_DETECTOR *HANDLE_SBR_TRANSIENT_DETECTOR;

#define TRAN_DET_LOOKAHEAD 2
#define TRAN_DET_START_FREQ 4500 /*start frequency for transient detection*/
#define TRAN_DET_STOP_FREQ 13500 /*stop frequency for transient detection*/
#define TRAN_DET_MIN_QMFBANDS \
  4

#define QMF_HP_dBd_SLOPE_FIX \
  FL2FXCONST_DBL(0.00075275f) /* 0.002266f/10 * log2(10) */
#define TRAN_DET_THRSHLD FL2FXCONST_DBL(5.0f / 8.0f)
#define TRAN_DET_THRSHLD_SCALE (3)

typedef struct {
  INT transientCandidates[32 + TRAN_DET_LOOKAHEAD];
  INT nTimeSlots;
  INT lookahead;
  INT startBand;
  INT stopBand;

  FIXP_DBL dBf_m[64];
  INT dBf_e[64];

  FIXP_DBL energy_timeSlots[32 + TRAN_DET_LOOKAHEAD];
  INT energy_timeSlots_scale[32 + TRAN_DET_LOOKAHEAD];

  FIXP_DBL delta_energy[32 + TRAN_DET_LOOKAHEAD];
  INT delta_energy_scale[32 + TRAN_DET_LOOKAHEAD];

  FIXP_DBL lowpass_energy[32 + TRAN_DET_LOOKAHEAD];
  INT lowpass_energy_scale[32 + TRAN_DET_LOOKAHEAD];
} FAST_TRAN_DETECTOR;
typedef FAST_TRAN_DETECTOR *HANDLE_FAST_TRAN_DET;

INT FDKsbrEnc_InitSbrFastTransientDetector(
    HANDLE_FAST_TRAN_DET h_sbrFastTransientDetector,
    const INT time_slots_per_frame, const INT bandwidth_qmf_slot,
    const INT no_qmf_channels, const INT sbr_qmf_1st_band);

void FDKsbrEnc_fastTransientDetect(
    const HANDLE_FAST_TRAN_DET h_sbrFastTransientDetector,
    const FIXP_DBL *const *Energies, const int *const scaleEnergies,
    const INT YBufferWriteOffset, UCHAR *const tran_vector);

void FDKsbrEnc_transientDetect(
    HANDLE_SBR_TRANSIENT_DETECTOR h_sbrTransientDetector, FIXP_DBL **Energies,
    INT *scaleEnergies, UCHAR *tran_vector, int YBufferWriteOffset,
    int YBufferSzShift, int timeStep, int frameMiddleBorder);

int FDKsbrEnc_InitSbrTransientDetector(
    HANDLE_SBR_TRANSIENT_DETECTOR h_sbrTransientDetector,
    UINT sbrSyntaxFlags, /* SBR syntax flags derived from AOT. */
    INT frameSize, INT sampleFreq, sbrConfigurationPtr params, int tran_fc,
    int no_cols, int no_rows, int YBufferWriteOffset, int YBufferSzShift,
    int frameShift, int tran_off);

void FDKsbrEnc_frameSplitter(
    FIXP_DBL **Energies, INT *scaleEnergies,
    HANDLE_SBR_TRANSIENT_DETECTOR h_sbrTransientDetector, UCHAR *freqBandTable,
    UCHAR *tran_vector, int YBufferWriteOffset, int YBufferSzShift, int nSfb,
    int timeStep, int no_cols, FIXP_DBL *tonality);
#endif
