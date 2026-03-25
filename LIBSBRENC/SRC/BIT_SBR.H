/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef BIT_SBR_H
#define BIT_SBR_H

#include "sbr_def.h"
#include "cmondata.h"
#include "fram_gen.h"

struct SBR_ENV_DATA;

struct SBR_BITSTREAM_DATA {
  INT TotalBits;
  INT PayloadBits;
  INT FillBits;
  INT HeaderActive;
  INT HeaderActiveDelay;
  INT NrSendHeaderData;
  INT CountSendHeaderData;
  INT rightBorderFIX;
};

typedef struct SBR_BITSTREAM_DATA *HANDLE_SBR_BITSTREAM_DATA;

struct SBR_HEADER_DATA {
  AMP_RES sbr_amp_res;
  INT sbr_start_frequency;
  INT sbr_stop_frequency;
  INT sbr_xover_band;
  INT sbr_noise_bands;
  INT sbr_data_extra;
  INT header_extra_1;
  INT header_extra_2;
  INT sbr_lc_stereo_mode;
  INT sbr_limiter_bands;
  INT sbr_limiter_gains;
  INT sbr_interpol_freq;
  INT sbr_smoothing_length;
  INT alterScale;
  INT freqScale;

  INT coupling;
  INT prev_coupling;
};
typedef struct SBR_HEADER_DATA *HANDLE_SBR_HEADER_DATA;

struct SBR_ENV_DATA {
  INT sbr_xpos_ctrl;
  FREQ_RES freq_res_fixfix[2];
  UCHAR fResTransIsLow;

  INVF_MODE sbr_invf_mode;
  INVF_MODE sbr_invf_mode_vec[MAX_NUM_NOISE_VALUES];

  XPOS_MODE sbr_xpos_mode;

  INT ienvelope[MAX_ENVELOPES][MAX_FREQ_COEFFS];

  INT codeBookScfLavBalance;
  INT codeBookScfLav;
  const INT *hufftableTimeC;
  const INT *hufftableFreqC;
  const UCHAR *hufftableTimeL;
  const UCHAR *hufftableFreqL;

  const INT *hufftableLevelTimeC;
  const INT *hufftableBalanceTimeC;
  const INT *hufftableLevelFreqC;
  const INT *hufftableBalanceFreqC;
  const UCHAR *hufftableLevelTimeL;
  const UCHAR *hufftableBalanceTimeL;
  const UCHAR *hufftableLevelFreqL;
  const UCHAR *hufftableBalanceFreqL;

  const UCHAR *hufftableNoiseTimeL;
  const INT *hufftableNoiseTimeC;
  const UCHAR *hufftableNoiseFreqL;
  const INT *hufftableNoiseFreqC;

  const UCHAR *hufftableNoiseLevelTimeL;
  const INT *hufftableNoiseLevelTimeC;
  const UCHAR *hufftableNoiseBalanceTimeL;
  const INT *hufftableNoiseBalanceTimeC;
  const UCHAR *hufftableNoiseLevelFreqL;
  const INT *hufftableNoiseLevelFreqC;
  const UCHAR *hufftableNoiseBalanceFreqL;
  const INT *hufftableNoiseBalanceFreqC;

  HANDLE_SBR_GRID hSbrBSGrid;

  INT noHarmonics;
  INT addHarmonicFlag;
  UCHAR addHarmonic[MAX_FREQ_COEFFS];

  INT si_sbr_start_env_bits_balance;
  INT si_sbr_start_env_bits;
  INT si_sbr_start_noise_bits_balance;
  INT si_sbr_start_noise_bits;

  INT noOfEnvelopes;
  INT noScfBands[MAX_ENVELOPES];
  INT domain_vec[MAX_ENVELOPES];
  INT domain_vec_noise[MAX_ENVELOPES];
  SCHAR sbr_noise_levels[MAX_FREQ_COEFFS];
  INT noOfnoisebands;

  INT balance;
  AMP_RES init_sbr_amp_res;
  AMP_RES currentAmpResFF;
  FIXP_DBL
  ton_HF[SBR_GLOBAL_TONALITY_VALUES];
  FIXP_DBL global_tonality;

  INT extended_data;
  INT extension_size;
  INT extension_id;
  UCHAR extended_data_buffer[SBR_EXTENDED_DATA_MAX_CNT];

  UCHAR ldGrid;
};
typedef struct SBR_ENV_DATA *HANDLE_SBR_ENV_DATA;

INT FDKsbrEnc_WriteEnvSingleChannelElement(
    struct SBR_HEADER_DATA *sbrHeaderData,
    struct T_PARAMETRIC_STEREO *hParametricStereo,
    struct SBR_BITSTREAM_DATA *sbrBitstreamData,
    struct SBR_ENV_DATA *sbrEnvData, struct COMMON_DATA *cmonData,
    UINT sbrSyntaxFlags);

INT FDKsbrEnc_WriteEnvChannelPairElement(
    struct SBR_HEADER_DATA *sbrHeaderData,
    struct T_PARAMETRIC_STEREO *hParametricStereo,
    struct SBR_BITSTREAM_DATA *sbrBitstreamData,
    struct SBR_ENV_DATA *sbrEnvDataLeft, struct SBR_ENV_DATA *sbrEnvDataRight,
    struct COMMON_DATA *cmonData, UINT sbrSyntaxFlags);

INT FDKsbrEnc_CountSbrChannelPairElement(
    struct SBR_HEADER_DATA *sbrHeaderData,
    struct T_PARAMETRIC_STEREO *hParametricStereo,
    struct SBR_BITSTREAM_DATA *sbrBitstreamData,
    struct SBR_ENV_DATA *sbrEnvDataLeft, struct SBR_ENV_DATA *sbrEnvDataRight,
    struct COMMON_DATA *cmonData, UINT sbrSyntaxFlags);

#endif
