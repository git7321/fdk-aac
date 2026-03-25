/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef PSDEC_H
#define PSDEC_H

#include "sbrdecoder.h"
#include "FDK_hybrid.h"

#include "FDK_decorrelate.h"
#include "FDK_bitstream.h"

#define SCAL_HEADROOM (2)

#define PS_EXTENSION_SIZE_BITS (4)
#define PS_EXTENSION_ESC_COUNT_BITS (8)

#define NO_QMF_CHANNELS (64)
#define MAX_NUM_COL (32)

#define NO_QMF_BANDS_HYBRID20 (3)
#define NO_SUB_QMF_CHANNELS (12)
#define HYBRID_FILTER_DELAY (6)

#define MAX_NO_PS_ENV (4 + 1)

#define NO_HI_RES_BINS (34)
#define NO_MID_RES_BINS (20)
#define NO_LOW_RES_BINS (10)

#define NO_HI_RES_IID_BINS (NO_HI_RES_BINS)
#define NO_HI_RES_ICC_BINS (NO_HI_RES_BINS)

#define NO_MID_RES_IID_BINS (NO_MID_RES_BINS)
#define NO_MID_RES_ICC_BINS (NO_MID_RES_BINS)

#define NO_LOW_RES_IID_BINS (NO_LOW_RES_BINS)
#define NO_LOW_RES_ICC_BINS (NO_LOW_RES_BINS)

#define SUBQMF_GROUPS (10)
#define QMF_GROUPS (12)

#define NO_IID_GROUPS (SUBQMF_GROUPS + QMF_GROUPS)

#define NO_IID_STEPS (7)       /* 1 .. + 7 */
#define NO_IID_STEPS_FINE (15) /* 1 .. +15 */
#define NO_ICC_STEPS (8)       /* 0 .. + 7 */

#define NO_IID_LEVELS (2 * NO_IID_STEPS + 1)           /* - 7 ..  + 7 */
#define NO_IID_LEVELS_FINE (2 * NO_IID_STEPS_FINE + 1) /* -15 ..  +15 */
#define NO_ICC_LEVELS (NO_ICC_STEPS)                   /*   0 ..  + 7 */

#define FIXP_SQRT05 ((FIXP_DBL)0x5a827980) /* 1/SQRT2 */

struct PS_DEC_COEFFICIENTS {
  FIXP_DBL H11r[NO_IID_GROUPS]; /*!< coefficients of the sub-subband groups */
  FIXP_DBL H12r[NO_IID_GROUPS]; /*!< coefficients of the sub-subband groups */
  FIXP_DBL H21r[NO_IID_GROUPS]; /*!< coefficients of the sub-subband groups */
  FIXP_DBL H22r[NO_IID_GROUPS]; /*!< coefficients of the sub-subband groups */

  FIXP_DBL
  DeltaH11r[NO_IID_GROUPS]; /*!< coefficients of the sub-subband groups */
  FIXP_DBL
  DeltaH12r[NO_IID_GROUPS]; /*!< coefficients of the sub-subband groups */
  FIXP_DBL
  DeltaH21r[NO_IID_GROUPS]; /*!< coefficients of the sub-subband groups */
  FIXP_DBL
  DeltaH22r[NO_IID_GROUPS]; /*!< coefficients of the sub-subband groups */

  SCHAR
  aaIidIndexMapped[MAX_NO_PS_ENV]
                  [NO_HI_RES_IID_BINS]; /*!< The mapped IID index for all
                                           envelopes and all IID bins */
  SCHAR
  aaIccIndexMapped[MAX_NO_PS_ENV]
                  [NO_HI_RES_ICC_BINS]; /*!< The mapped ICC index for all
                                           envelopes and all ICC bins */
};

typedef enum { ppt_none = 0, ppt_mpeg = 1, ppt_drm = 2 } PS_PAYLOAD_TYPE;

typedef struct {
  UCHAR bPsHeaderValid; /*!< set if new header is available from bitstream */

  UCHAR bEnableIid; /*!< One bit denoting the presence of IID parameters */
  UCHAR bEnableIcc; /*!< One bit denoting the presence of ICC parameters */
  UCHAR bEnableExt; /*!< The PS extension layer is enabled using the enable_ext
                       bit. If it is set to %1 the IPD and OPD parameters are
                       sent. If it is disabled, i.e. %0, the extension layer is
                       skipped.   */

  UCHAR
  modeIid;       /*!< The configuration of IID parameters (number of bands and
                      quantisation grid, iid_quant) is determined by iid_mode.   */
  UCHAR modeIcc; /*!< The configuration of Inter-channel Coherence parameters
                      (number of bands and quantisation grid) is determined by
                      icc_mode. */

  UCHAR freqResIid; /*!< 0=low, 1=mid or 2=high frequency resolution for iid */
  UCHAR freqResIcc; /*!< 0=low, 1=mid or 2=high frequency resolution for icc */

  UCHAR bFineIidQ; /*!< Use fine Iid quantisation. */

  UCHAR bFrameClass; /*!< The frame_class bit determines whether the parameter
                          positions of the current frame are uniformly spaced
                          accross the frame or they are defined using the
                        positions described by border_position.
                      */

  UCHAR noEnv; /*!< The number of envelopes per frame */
  UCHAR aEnvStartStop[MAX_NO_PS_ENV + 1]; /*!< In case of variable parameter
                                             spacing the parameter positions are
                                             determined by border_position */

  SCHAR abIidDtFlag[MAX_NO_PS_ENV]; /*!< Deltacoding time/freq flag for IID, 0
                                       => freq           */
  SCHAR abIccDtFlag[MAX_NO_PS_ENV]; /*!< Deltacoding time/freq flag for ICC, 0
                                       => freq           */

  SCHAR
  aaIidIndex[MAX_NO_PS_ENV]
            [NO_HI_RES_IID_BINS]; /*!< The IID index for all envelopes and
                                     all IID bins        */
  SCHAR
  aaIccIndex[MAX_NO_PS_ENV]
            [NO_HI_RES_ICC_BINS]; /*!< The ICC index for all envelopes and
                                     all ICC bins        */

} MPEG_PS_BS_DATA;

struct PS_DEC {
  SCHAR noSubSamples;
  SCHAR noChannels;

  SCHAR procFrameBased; /*!< Helper to detected switching from frame based to
                           slot based processing
                         */

  PS_PAYLOAD_TYPE
  bPsDataAvail[(1) + 1]; /*!< set if new data available from bitstream */
  UCHAR psDecodedPrv;    /*!< set if PS has been processed in the last frame */

  /* helpers for frame delay line */
  UCHAR bsLastSlot;  /*!< Index of last read slot.  */
  UCHAR bsReadSlot;  /*!< Index of current read slot for additional delay.  */
  UCHAR processSlot; /*!< Index of current slot for processing (need for add.
                        delay).   */

  union { /* Bitstream data */
    MPEG_PS_BS_DATA
    mpeg; /*!< Struct containing all MPEG specific PS data from bitstream.
           */
  } bsData[(1) + 1];

  shouldBeUnion { /* Static data */
    struct {
      SCHAR aIidPrevFrameIndex[NO_HI_RES_IID_BINS]; /*!< The IID index for
                                                       previous frame */
      SCHAR aIccPrevFrameIndex[NO_HI_RES_ICC_BINS]; /*!< The ICC index for
                                                       previous frame */
      UCHAR
      bPrevFrameFineIidQ;   /*!< The IID quantization of the previous frame */
      UCHAR prevFreqResIid; /*!< Frequency resolution for IID of the previous
                               frame            */
      UCHAR prevFreqResIcc; /*!< Frequency resolution for ICC of the previous
                               frame            */
      UCHAR lastUsb; /*!< uppermost WMF delay band of last frame          */

      FIXP_DBL pHybridAnaStatesLFdmx
          [2 * 13 * NO_QMF_BANDS_HYBRID20]; /*!< Memory used in hybrid analysis
                                                 for filter states. */
      FDK_ANA_HYB_FILTER hybridAnalysis;
      FDK_SYN_HYB_FILTER hybridSynthesis[2];

      DECORR_DEC apDecor; /*!< Decorrelator instance. */
      FIXP_DBL decorrBufferCplx[(2 * ((825) + (373)))];

      FIXP_DBL h11rPrev[NO_IID_GROUPS]; /*!< previous calculated h(xy)
                                           coefficients */
      FIXP_DBL h12rPrev[NO_IID_GROUPS]; /*!< previous calculated h(xy)
                                           coefficients */
      FIXP_DBL h21rPrev[NO_IID_GROUPS]; /*!< previous calculated h(xy)
                                           coefficients */
      FIXP_DBL h22rPrev[NO_IID_GROUPS]; /*!< previous calculated h(xy)
                                           coefficients */

      PS_DEC_COEFFICIENTS
      *pCoef; /*!< temporal coefficients are on reusable scratch memory */

    } mpeg;
  }
  specificTo;
};

typedef struct PS_DEC *HANDLE_PS_DEC;

int CreatePsDec(HANDLE_PS_DEC *h_PS_DEC, int aacSamplesPerFrame);

int DeletePsDec(HANDLE_PS_DEC *h_PS_DEC);

void PreparePsProcessing(HANDLE_PS_DEC h_ps_d,
                         const FIXP_DBL *const *const rIntBufferLeft,
                         const FIXP_DBL *const *const iIntBufferLeft,
                         const int scaleFactorLowBand);

void initSlotBasedRotation(HANDLE_PS_DEC h_ps_d, int env, int usb);

void ApplyPsSlot(
    HANDLE_PS_DEC h_ps_d,      /* parametric stereo decoder handle    */
    FIXP_DBL **rIntBufferLeft, /* real values of left qmf timeslot    */
    FIXP_DBL **iIntBufferLeft, /* imag values of left qmf timeslot    */
    FIXP_DBL *rIntBufferRight, /* real values of right qmf timeslot   */
    FIXP_DBL *iIntBufferRight, /* imag values of right qmf timeslot   */
    const int scaleFactorLowBand_no_ov, const int scaleFactorLowBand,
    const int scaleFactorHighBand, const int lsb, const int usb);

#endif /* PSDEC_H */
