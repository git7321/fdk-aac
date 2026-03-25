/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: individual channel stream info

*******************************************************************************/

#ifndef CHANNELINFO_H
#define CHANNELINFO_H

#include "common_fix.h"
#include "aac_rom.h"
#include "aacdecoder_lib.h"
#include "FDK_bitstream.h"
#include "overlapadd.h"
#include "mdct.h"
#include "stereo.h"
#include "pulsedata.h"
#include "aacdec_tns.h"
#include "aacdec_pns.h"
#include "aacdec_hcr_types.h"
#include "rvlc_info.h"
#include "usacdec_acelp.h"
#include "usacdec_const.h"
#include "usacdec_rom.h"
#include "ac_arith_coder.h"
#include "conceal_types.h"
#include "aacdec_drc_types.h"

#define WB_SECTION_SIZE (1024 * 2)
#define DRM_BS_BUFFER_SIZE \
  (512)

typedef enum {
  AACDEC_RENDER_INVALID = 0,
  AACDEC_RENDER_IMDCT,
  AACDEC_RENDER_ELDFB,
  AACDEC_RENDER_LPD,
  AACDEC_RENDER_INTIMDCT
} AACDEC_RENDER_MODE;

enum { MAX_QUANTIZED_VALUE = 8191 };

typedef enum { FD_LONG, FD_SHORT, LPD } USAC_COREMODE;

typedef struct {
  const SHORT *ScaleFactorBands_Long;
  const SHORT *ScaleFactorBands_Short;
  UCHAR NumberOfScaleFactorBands_Long;
  UCHAR NumberOfScaleFactorBands_Short;
  UINT samplingRateIndex;
  UINT samplingRate;
} SamplingRateInfo;

typedef struct {
  UCHAR CommonWindow;
  UCHAR GlobalGain;

} CRawDataInfo;

typedef struct {
  UCHAR WindowGroupLength[8];
  UCHAR WindowGroups;
  UCHAR Valid;

  UCHAR WindowShape;
  BLOCK_TYPE WindowSequence;
  UCHAR MaxSfBands;
  UCHAR max_sfb_ste;
  UCHAR ScaleFactorGrouping;

  UCHAR TotalSfBands;

} CIcsInfo;

enum {
  ZERO_HCB = 0,
  ESCBOOK = 11,
  NSPECBOOKS = ESCBOOK + 1,
  BOOKSCL = NSPECBOOKS,
  NOISE_HCB = 13,
  INTENSITY_HCB2 = 14,
  INTENSITY_HCB = 15,
  LAST_HCB
};

typedef struct {
  CJointStereoPersistentData jointStereoPersistentData;
} CpePersistentData;

typedef struct {
  FIXP_DBL *pOverlapBuffer;
  mdct_t IMdct;

  CArcoData *hArCo;
  INT pnsCurrentSeed;

  FIXP_DBL old_synth[PIT_MAX_MAX - L_SUBFR];
  INT old_T_pf[SYN_SFD];
  FIXP_DBL old_gain_pf[SYN_SFD];
  FIXP_DBL mem_bpf[L_FILT + L_SUBFR];
  UCHAR
  old_bpf_control_info;

  USAC_COREMODE last_core_mode;
  UCHAR last_lpd_mode;
  UCHAR last_last_lpd_mode;
  UCHAR last_lpc_lost;

  FIXP_LPC
  lpc4_lsf[M_LP_FILTER_ORDER];
  FIXP_LPC lsf_adaptive_mean[M_LP_FILTER_ORDER];
  FIXP_LPC lp_coeff_old[2][M_LP_FILTER_ORDER];
  INT lp_coeff_old_exp[2];

  FIXP_SGL
  oldStability;
  UINT numLostLpdFrames;

  FIXP_DBL last_tcx_gain;
  INT last_tcx_gain_e;
  FIXP_DBL last_alfd_gains[32];
  SHORT last_tcx_pitch;
  UCHAR last_tcx_noise_factor;

  CAcelpStaticMem acelp;
  ULONG nfRandomSeed;

  CDrcChannelData drcData;
  CConcealmentInfo concealmentInfo;

  CpePersistentData *pCpeStaticData;
} CAacDecoderStaticChannelInfo;

typedef struct {
  SHORT aScaleFactor[(
      8 * 16)];
  SHORT aSfbScale[(8 * 16)];
  UCHAR
  aCodeBook[(8 * 16)];
  UCHAR band_is_noise[(8 * 16)];
  CTnsData TnsData;
  CRawDataInfo RawDataInfo;

  shouldBeUnion {
    struct {
      CPulseData PulseData;
      SHORT aNumLineInSec4Hcr[MAX_SFB_HCR];
      UCHAR
      aCodeBooks4Hcr[MAX_SFB_HCR];
      SHORT lenOfReorderedSpectralData;
      SCHAR lenOfLongestCodeword;
      SCHAR numberSection;
      SCHAR rvlcCurrentScaleFactorOK;
      SCHAR rvlcIntensityUsed;
    } aac;
    struct {
      UCHAR fd_noise_level_and_offset;
      UCHAR tns_active;
      UCHAR tns_on_lr;
      UCHAR tcx_noise_factor[4];
      UCHAR tcx_global_gain[4];
    } usac;
  }
  specificTo;

} CAacDecoderDynamicData;

typedef shouldBeUnion {
  UCHAR DrmBsBuffer[DRM_BS_BUFFER_SIZE];

  FIXP_DBL mdctOutTemp[1024];
  FIXP_DBL synth_buf[(PIT_MAX_MAX + SYN_DELAY + L_FRAME_PLUS)];
  FIXP_DBL workBuffer[WB_SECTION_SIZE];
}
CWorkBufferCore1;

typedef struct {
  CAacDecoderDynamicData pAacDecoderDynamicData[2];

  CPnsInterChannelData pnsInterChannelData;
  INT pnsRandomSeed[(8 * 16)];

  CJointStereoData jointStereoData;

  shouldBeUnion {
    struct {
      CErHcrInfo erHcrInfo;
      CErRvlcInfo erRvlcInfo;
      SHORT aRvlcScfEsc[RVLC_MAX_SFB];
      SHORT aRvlcScfFwd[RVLC_MAX_SFB];
      SHORT aRvlcScfBwd[RVLC_MAX_SFB];
    } aac;
  }
  overlay;

} CAacDecoderCommonData;

typedef struct {
  CWorkBufferCore1 *pWorkBufferCore1;
  CCplxPredictionData *cplxPredictionData;
} CAacDecoderCommonStaticData;

typedef struct {
  shouldBeUnion {
    struct {
      FIXP_DBL fac_data0[LFAC];
      SCHAR fac_data_e[4];
      FIXP_DBL
      *fac_data[4];

      UCHAR core_mode;
      USAC_COREMODE
      core_mode_last;
      UCHAR lpd_mode_last;
      UCHAR mod[4];
      UCHAR bpf_control_info;

      FIXP_LPC lsp_coeff[5][M_LP_FILTER_ORDER];
      FIXP_LPC
      lp_coeff[5][M_LP_FILTER_ORDER];
      INT lp_coeff_exp[5];
      FIXP_LPC lsf_adaptive_mean_cand
          [M_LP_FILTER_ORDER];
      FIXP_SGL aStability[4];

      CAcelpChannelData acelp[4];

      FIXP_DBL tcx_gain[4];
      SCHAR tcx_gain_e[4];
    } usac;

    struct {
      CPnsData PnsData;
    } aac;
  }
  data;

  SPECTRAL_PTR pSpectralCoefficient;
  SHORT specScale[8];
  CIcsInfo icsInfo;
  INT granuleLength;
  UCHAR ElementInstanceTag;

  AACDEC_RENDER_MODE renderMode;

  CAacDecoderDynamicData *
      pDynData;
  CAacDecoderCommonData
      *pComData;
  CAacDecoderCommonStaticData *pComStaticData;

  int currAliasingSymmetry;
} CAacDecoderChannelInfo;

AAC_DECODER_ERROR getSamplingRateInfo(SamplingRateInfo *t, UINT samplesPerFrame,
                                      UINT samplingRateIndex,
                                      UINT samplingRate);

AAC_DECODER_ERROR IcsReadMaxSfb(HANDLE_FDK_BITSTREAM hBs, CIcsInfo *pIcsInfo,
                                const SamplingRateInfo *pSamplingRateInfo);

AAC_DECODER_ERROR IcsRead(HANDLE_FDK_BITSTREAM bs, CIcsInfo *pIcsInfo,
                          const SamplingRateInfo *SamplingRateInfoTable,
                          const UINT flags);

void CJointStereo_ApplyMS(
    CAacDecoderChannelInfo *pAacDecoderChannelInfo[2],
    CAacDecoderStaticChannelInfo *pAacDecoderStaticChannelInfo[2],
    FIXP_DBL *spectrumL, FIXP_DBL *spectrumR, SHORT *SFBleftScale,
    SHORT *SFBrightScale, SHORT *specScaleL, SHORT *specScaleR,
    const SHORT *pScaleFactorBandOffsets, const UCHAR *pWindowGroupLength,
    const int windowGroups, const int max_sfb_ste_outside,
    const int scaleFactorBandsTransmittedL,
    const int scaleFactorBandsTransmittedR, FIXP_DBL *store_dmx_re_prev,
    SHORT *store_dmx_re_prev_e, const int mainband_flag);

void CJointStereo_ApplyIS(CAacDecoderChannelInfo *pAacDecoderChannelInfo[2],
                          const short *pScaleFactorBandOffsets,
                          const UCHAR *pWindowGroupLength,
                          const int windowGroups,
                          const int scaleFactorBandsTransmitted);

int CPns_IsPnsUsed(const CPnsData *pPnsData, const int group, const int band);

void CPns_SetCorrelation(CPnsData *pPnsData, const int group, const int band,
                         const int outofphase);

inline UCHAR IsValid(const CIcsInfo *pIcsInfo) { return pIcsInfo->Valid; }

inline UCHAR IsLongBlock(const CIcsInfo *pIcsInfo) {
  return (pIcsInfo->WindowSequence != BLOCK_SHORT);
}

inline UCHAR GetWindowShape(const CIcsInfo *pIcsInfo) {
  return pIcsInfo->WindowShape;
}

inline BLOCK_TYPE GetWindowSequence(const CIcsInfo *pIcsInfo) {
  return pIcsInfo->WindowSequence;
}

inline const SHORT *GetScaleFactorBandOffsets(
    const CIcsInfo *pIcsInfo, const SamplingRateInfo *samplingRateInfo) {
  if (IsLongBlock(pIcsInfo)) {
    return samplingRateInfo->ScaleFactorBands_Long;
  } else {
    return samplingRateInfo->ScaleFactorBands_Short;
  }
}

inline UCHAR GetNumberOfScaleFactorBands(
    const CIcsInfo *pIcsInfo, const SamplingRateInfo *samplingRateInfo) {
  if (IsLongBlock(pIcsInfo)) {
    return samplingRateInfo->NumberOfScaleFactorBands_Long;
  } else {
    return samplingRateInfo->NumberOfScaleFactorBands_Short;
  }
}

inline int GetWindowsPerFrame(const CIcsInfo *pIcsInfo) {
  return (pIcsInfo->WindowSequence == BLOCK_SHORT) ? 8 : 1;
}

inline UCHAR GetWindowGroups(const CIcsInfo *pIcsInfo) {
  return pIcsInfo->WindowGroups;
}

inline UCHAR GetWindowGroupLength(const CIcsInfo *pIcsInfo, const INT index) {
  return pIcsInfo->WindowGroupLength[index];
}

inline const UCHAR *GetWindowGroupLengthTable(const CIcsInfo *pIcsInfo) {
  return pIcsInfo->WindowGroupLength;
}

inline UCHAR GetScaleFactorBandsTransmitted(const CIcsInfo *pIcsInfo) {
  return pIcsInfo->MaxSfBands;
}

inline UCHAR GetScaleMaxFactorBandsTransmitted(const CIcsInfo *pIcsInfo0,
                                               const CIcsInfo *pIcsInfo1) {
  return fMax(pIcsInfo0->MaxSfBands, pIcsInfo1->MaxSfBands);
}

inline UCHAR GetScaleFactorBandsTotal(const CIcsInfo *pIcsInfo) {
  return pIcsInfo->TotalSfBands;
}

inline UCHAR GetMaximumTnsBands(const CIcsInfo *pIcsInfo,
                                const int samplingRateIndex) {
  return tns_max_bands_tbl[samplingRateIndex][!IsLongBlock(pIcsInfo)];
}

#endif /* #ifndef CHANNELINFO_H */
