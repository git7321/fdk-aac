/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef MH_DET_H
#define MH_DET_H

#include "sbr_encoder.h"
#include "fram_gen.h"

typedef struct {
  FIXP_DBL thresHoldDiff;      /*!< threshold for tonality difference */
  FIXP_DBL thresHoldDiffGuide; /*!< threshold for tonality difference for the
                                  guide */
  FIXP_DBL thresHoldTone;      /*!< threshold for tonality for a sine */
  FIXP_DBL invThresHoldTone;
  FIXP_DBL thresHoldToneGuide; /*!< threshold for tonality for a sine for the
                                  guide */
  FIXP_DBL sfmThresSbr;    /*!< tonality flatness measure threshold for the SBR
                              signal.*/
  FIXP_DBL sfmThresOrig;   /*!< tonality flatness measure threshold for the
                              original signal.*/
  FIXP_DBL decayGuideOrig; /*!< decay value of the tonality value of the guide
                              for the tone. */
  FIXP_DBL decayGuideDiff; /*!< decay value of the tonality value of the guide
                              for the tonality difference. */
  FIXP_DBL derivThresMaxLD64;   /*!< threshold for detecting LP character in a
                                   signal. */
  FIXP_DBL derivThresBelowLD64; /*!< threshold for detecting LP character in a
                                   signal. */
  FIXP_DBL derivThresAboveLD64; /*!< threshold for detecting LP character in a
                                   signal. */
} THRES_HOLDS;

typedef struct {
  INT deltaTime; /*!< maximum allowed transient distance (from frame border in
                    number of qmf subband sample) for a frame to be considered a
                    transient frame.*/
  THRES_HOLDS thresHolds; /*!< the thresholds used for detection. */
  INT maxComp; /*!< maximum alllowed compensation factor for the envelope data.
                */
} DETECTOR_PARAMETERS_MH;

typedef struct {
  FIXP_DBL *guideVectorDiff;
  FIXP_DBL *guideVectorOrig;
  UCHAR *guideVectorDetected;
} GUIDE_VECTORS;

typedef struct {
  INT qmfNoChannels;
  INT nSfb;
  INT sampleFreq;
  INT previousTransientFlag;
  INT previousTransientFrame;
  INT previousTransientPos;

  INT noVecPerFrame;
  INT transientPosOffset;

  INT move;
  INT totNoEst;
  INT noEstPerFrame;
  INT timeSlots;

  UCHAR *guideScfb;
  UCHAR *prevEnvelopeCompensation;
  UCHAR *detectionVectors[MAX_NO_OF_ESTIMATES];
  FIXP_DBL tonalityDiff[MAX_NO_OF_ESTIMATES / 2][MAX_FREQ_COEFFS];
  FIXP_DBL sfmOrig[MAX_NO_OF_ESTIMATES / 2][MAX_FREQ_COEFFS];
  FIXP_DBL sfmSbr[MAX_NO_OF_ESTIMATES / 2][MAX_FREQ_COEFFS];
  const DETECTOR_PARAMETERS_MH *mhParams;
  GUIDE_VECTORS guideVectors[MAX_NO_OF_ESTIMATES];
} SBR_MISSING_HARMONICS_DETECTOR;

typedef SBR_MISSING_HARMONICS_DETECTOR *HANDLE_SBR_MISSING_HARMONICS_DETECTOR;

void FDKsbrEnc_SbrMissingHarmonicsDetectorQmf(
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR h_sbrMissingHarmonicsDetector,
    FIXP_DBL **pQuotaBuffer, INT **pSignBuffer, SCHAR *indexVector,
    const SBR_FRAME_INFO *pFrameInfo, const UCHAR *pTranInfo,
    INT *pAddHarmonicsFlag, UCHAR *pAddHarmonicsScaleFactorBands,
    const UCHAR *freqBandTable, INT nSfb, UCHAR *envelopeCompensation,
    FIXP_DBL *pNrgVector);

INT FDKsbrEnc_CreateSbrMissingHarmonicsDetector(
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR hSbrMHDet, INT chan);

INT FDKsbrEnc_InitSbrMissingHarmonicsDetector(
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR h_sbrMissingHarmonicsDetector,
    INT sampleFreq, INT frameSize, INT nSfb, INT qmfNoChannels, INT totNoEst,
    INT move, INT noEstPerFrame, UINT sbrSyntaxFlags);

void FDKsbrEnc_DeleteSbrMissingHarmonicsDetector(
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR h_sbrMissingHarmonicsDetector);

INT FDKsbrEnc_ResetSbrMissingHarmonicsDetector(
    HANDLE_SBR_MISSING_HARMONICS_DETECTOR hSbrMissingHarmonicsDetector,
    INT nSfb);

#endif
