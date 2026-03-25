/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef LPP_TRAN_H
#define LPP_TRAN_H

#include "sbrdecoder.h"
#include "hbe.h"
#include "qmf.h"

#define QMF_OUT_SCALE 8

#define MAX_NOISE_ENVELOPES 2
#define MAX_NOISE_COEFFS 5
#define MAX_NUM_NOISE_VALUES (MAX_NOISE_ENVELOPES * MAX_NOISE_COEFFS)
#define MAX_NUM_LIMITERS 12

#define MAX_ENVELOPES_LEGACY 5
#define MAX_ENVELOPES_USAC 8
#define MAX_ENVELOPES MAX_ENVELOPES_USAC

#define MAX_FREQ_COEFFS_DUAL_RATE 48
#define MAX_FREQ_COEFFS_QUAD_RATE 56
#define MAX_FREQ_COEFFS MAX_FREQ_COEFFS_QUAD_RATE

#define MAX_FREQ_COEFFS_FS44100 35
#define MAX_FREQ_COEFFS_FS48000 32

#define MAX_NUM_ENVELOPE_VALUES (MAX_ENVELOPES * MAX_FREQ_COEFFS)

#define MAX_GAIN_EXP 34
#define MAX_GAIN_CONCEAL_EXP 1

#define LPC_ORDER 2

#define MAX_INVF_BANDS MAX_NOISE_COEFFS

#define MAX_NUM_PATCHES 6
#define SHIFT_START_SB 1 /*!< lowest subband of source range */

typedef enum {
  INVF_OFF = 0,
  INVF_LOW_LEVEL,
  INVF_MID_LEVEL,
  INVF_HIGH_LEVEL,
  INVF_SWITCHED /* not a real choice but used here to control behaviour */
} INVF_MODE;

/** parameter set for one single patch */
typedef struct {
  UCHAR sourceStartBand; /*!< first band in lowbands where to take the samples
                            from */
  UCHAR
  sourceStopBand;       /*!< first band in lowbands which is not included in the
                           patch anymore */
  UCHAR guardStartBand; /*!< first band in highbands to be filled with zeros in
                           order to reduce interferences between patches */
  UCHAR
  targetStartBand;       /*!< first band in highbands to be filled with whitened
                            lowband signal */
  UCHAR targetBandOffs;  /*!< difference between 'startTargetBand' and
                            'startSourceBand' */
  UCHAR numBandsInPatch; /*!< number of consecutive bands in this one patch */
} PATCH_PARAM;

typedef struct {
  FIXP_DBL off; /*!< bw factor for signal OFF */
  FIXP_DBL transitionLevel;
  FIXP_DBL lowLevel;  /*!< bw factor for signal LOW_LEVEL */
  FIXP_DBL midLevel;  /*!< bw factor for signal MID_LEVEL */
  FIXP_DBL highLevel; /*!< bw factor for signal HIGH_LEVEL */
} WHITENING_FACTORS;

typedef struct {
  UCHAR nCols;           /*!< number subsamples of a codec frame */
  UCHAR noOfPatches;     /*!< number of patches */
  UCHAR lbStartPatching; /*!< first band of lowbands that will be patched */
  UCHAR lbStopPatching;  /*!< first band that won't be patched anymore*/
  UCHAR bwBorders[MAX_NUM_NOISE_VALUES]; /*!< spectral bands with different
                                            inverse filtering levels */

  PATCH_PARAM
  patchParam[MAX_NUM_PATCHES + 1]; /*!< new parameter set for patching */
  WHITENING_FACTORS
  whFactors;     /*!< the pole moving factors for certain
                    whitening levels as indicated     in the bitstream
                    depending on the crossover frequency */
  UCHAR overlap; /*!< Overlap size */
} TRANSPOSER_SETTINGS;

typedef struct {
  TRANSPOSER_SETTINGS *pSettings; /*!< Common settings for both channels */
  FIXP_DBL
  bwVectorOld[MAX_NUM_PATCHES]; /*!< pole moving factors of past frame */
  FIXP_DBL lpcFilterStatesRealLegSBR[LPC_ORDER + (3 * (4))][(
      32)]; /*!< pointer array to save filter states */

  FIXP_DBL lpcFilterStatesImagLegSBR[LPC_ORDER + (3 * (4))][(
      32)]; /*!< pointer array to save filter states */

  FIXP_DBL lpcFilterStatesRealHBE[LPC_ORDER + (3 * (4))][(
      64)]; /*!< pointer array to save filter states */
  FIXP_DBL lpcFilterStatesImagHBE[LPC_ORDER + (3 * (4))][(
      64)]; /*!< pointer array to save filter states */
} SBR_LPP_TRANS;

typedef SBR_LPP_TRANS *HANDLE_SBR_LPP_TRANS;

void lppTransposer(HANDLE_SBR_LPP_TRANS hLppTrans,
                   QMF_SCALE_FACTOR *sbrScaleFactor, FIXP_DBL **qmfBufferReal,

                   FIXP_DBL *degreeAlias, FIXP_DBL **qmfBufferImag,
                   const int useLP, const int fPreWhitening,
                   const int v_k_master0, const int timeStep,
                   const int firstSlotOffset, const int lastSlotOffset,
                   const int nInvfBands, INVF_MODE *sbr_invf_mode,
                   INVF_MODE *sbr_invf_mode_prev);

void lppTransposerHBE(
    HANDLE_SBR_LPP_TRANS hLppTrans, /*!< Handle of lpp transposer  */
    HANDLE_HBE_TRANSPOSER hQmfTransposer,
    QMF_SCALE_FACTOR *sbrScaleFactor, /*!< Scaling factors */
    FIXP_DBL **qmfBufferReal, /*!< Pointer to pointer to real part of subband
                                 samples (source) */
    FIXP_DBL **qmfBufferImag, /*!< Pointer to pointer to imaginary part of
                                 subband samples (source) */
    const int timeStep,       /*!< Time step of envelope */
    const int firstSlotOffs,  /*!< Start position in time */
    const int lastSlotOffs,   /*!< Number of overlap-slots into next frame */
    const int nInvfBands,     /*!< Number of bands for inverse filtering */
    INVF_MODE *sbr_invf_mode, /*!< Current inverse filtering modes */
    INVF_MODE *sbr_invf_mode_prev /*!< Previous inverse filtering modes */
);

SBR_ERROR
createLppTransposer(HANDLE_SBR_LPP_TRANS hLppTrans,
                    TRANSPOSER_SETTINGS *pSettings, const int highBandStartSb,
                    UCHAR *v_k_master, const int numMaster, const int usb,
                    const int timeSlots, const int nCols, UCHAR *noiseBandTable,
                    const int noNoiseBands, UINT fs, const int chan,
                    const int overlap);

SBR_ERROR
resetLppTransposer(HANDLE_SBR_LPP_TRANS hLppTrans, UCHAR highBandStartSb,
                   UCHAR *v_k_master, UCHAR numMaster, UCHAR *noiseBandTable,
                   UCHAR noNoiseBands, UCHAR usb, UINT fs);

#endif /* LPP_TRAN_H */
