/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef RESAMPLER_H
#define RESAMPLER_H

#include "common_fix.h"

#define MAXNR_SECTIONS (15)

typedef FIXP_DBL FIXP_BQS;

typedef struct {
  FIXP_BQS states[MAXNR_SECTIONS + 1][2]; /*! state buffer */
  const FIXP_SGL *coeffa;                 /*! pointer to filter coeffs */
  FIXP_DBL gain;                          /*! overall gain factor */
  int Wc;                                 /*! normalized cutoff freq * 1000 */
  int noCoeffs;                           /*! number of filter coeffs sets */
  int ptr;                                /*! index to rinbuffers */
} LP_FILTER;

typedef struct {
  LP_FILTER downFilter; /*! filter instance */
  int ratio;            /*! downsampling ration */
  int delay;            /*! downsampling delay (source fs)   */
  int pending;          /*! number of pending output samples */
} DOWNSAMPLER;

INT FDKaacEnc_InitDownsampler(
    DOWNSAMPLER *DownSampler, /*!< pointer to downsampler instance */
    INT Wc,                   /*!< normalized cutoff freq * 1000 */
    INT ratio);               /*!< downsampler ratio */

INT FDKaacEnc_Downsample(
    DOWNSAMPLER *DownSampler, /*!< pointer to downsampler instance */
    INT_PCM *inSamples,       /*!< pointer to input samples */
    INT numInSamples,         /*!< number  of input samples  */
    INT_PCM *outSamples,      /*!< pointer to output samples */
    INT *numOutSamples);      /*!< pointer tp number of output samples */

#endif /* RESAMPLER_H */
