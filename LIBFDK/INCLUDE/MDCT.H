/******************* Library for basic calculation routines ********************

   Author(s):   Manuel Jander, Josef Hoepfl, Youliy Ninov, Daniel Hagel

   Description: MDCT/MDST routines

*******************************************************************************/

#ifndef MDCT_H
#define MDCT_H

#include "common_fix.h"

#define MDCT_OUT_HEADROOM 2

#define PCM_OUT_BITS DFRACT_BITS
#define PCM_OUT_HEADROOM 8

#define MDCT_OUTPUT_SCALE (-MDCT_OUT_HEADROOM + (DFRACT_BITS - PCM_OUT_BITS))
#define MDCT_OUTPUT_GAIN 16

#define IMDCT_SCALE(x, s) \
  SATURATE_RIGHT_SHIFT((x), ((s) + MDCT_OUTPUT_SCALE), PCM_OUT_BITS)
#define IMDCT_SCALE_DBL(x) (FIXP_DBL)(x)
#define IMDCT_SCALE_DBL_LSH1(x) SATURATE_LEFT_SHIFT_ALT((x), 1, DFRACT_BITS)

#define MLT_FLAG_CURR_ALIAS_SYMMETRY 1

typedef enum {
  BLOCK_LONG = 0,
  BLOCK_START,
  BLOCK_SHORT,
  BLOCK_STOP
} BLOCK_TYPE;

typedef enum { SHAPE_SINE = 0, SHAPE_KBD, SHAPE_LOL } WINDOW_SHAPE;

typedef struct {
  union {
    FIXP_DBL *freq;
    FIXP_DBL *time;
  } overlap;

  const FIXP_WTP *prev_wrs;
  int prev_tl;
  int prev_nr;
  int prev_fr;
  int ov_offset;
  int ov_size;

  int prevAliasSymmetry;
  int prevPrevAliasSymmetry;

  FIXP_DBL *pFacZir;
  FIXP_DBL *pAsymOvlp;
} mdct_t;

typedef mdct_t *H_MDCT;

void mdct_init(H_MDCT hMdct, FIXP_DBL *overlap, INT overlapBufferSize);

INT mdct_block(H_MDCT hMdct, const INT_PCM *pTimeData, const INT noInSamples,
               FIXP_DBL *RESTRICT mdctData, const INT nSpec, const INT tl,
               const FIXP_WTP *pRightWindowPart, const INT fr,
               SHORT *pMdctData_e);

void imdct_gain(FIXP_DBL *pGain, int *pExponent, int tl);

INT imdct_drain(H_MDCT hMdct, FIXP_DBL *pTimeData, INT nrSamplesRoom);

INT imdct_copy_ov_and_nr(H_MDCT hMdct, FIXP_DBL *pTimeData, INT nrSamples);

void imdct_adapt_parameters(H_MDCT hMdct, int *pfl, int *pnl, int tl,
                            const FIXP_WTP *wls, int noOutSamples);

INT imlt_block(H_MDCT hMdct, FIXP_DBL *output, FIXP_DBL *spectrum,
               const SHORT scalefactor[], const INT nSpec,
               const INT noOutSamples, const INT tl, const FIXP_WTP *wls,
               INT fl, const FIXP_WTP *wrs, const INT fr, FIXP_DBL gain,
               int flags);

#endif /* MDCT_H */
