/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef CODE_ENV_H
#define CODE_ENV_H

#include "sbr_def.h"
#include "bit_sbr.h"
#include "fram_gen.h"

typedef struct {
  INT offset;
  INT upDate;
  INT nSfb[2];
  SCHAR sfb_nrg_prev[MAX_FREQ_COEFFS];
  INT deltaTAcrossFrames;
  FIXP_DBL dF_edge_1stEnv;
  FIXP_DBL dF_edge_incr;
  INT dF_edge_incr_fac;

  INT codeBookScfLavTime;
  INT codeBookScfLavFreq;

  INT codeBookScfLavLevelTime;
  INT codeBookScfLavLevelFreq;
  INT codeBookScfLavBalanceTime;
  INT codeBookScfLavBalanceFreq;

  INT start_bits;
  INT start_bits_balance;

  const UCHAR *hufftableTimeL;
  const UCHAR *hufftableFreqL;

  const UCHAR *hufftableLevelTimeL;
  const UCHAR *hufftableBalanceTimeL;
  const UCHAR *hufftableLevelFreqL;
  const UCHAR *hufftableBalanceFreqL;
} SBR_CODE_ENVELOPE;
typedef SBR_CODE_ENVELOPE *HANDLE_SBR_CODE_ENVELOPE;

void FDKsbrEnc_codeEnvelope(SCHAR *sfb_nrg, const FREQ_RES *freq_res,
                            SBR_CODE_ENVELOPE *h_sbrCodeEnvelope,
                            INT *directionVec, INT coupling, INT nEnvelopes,
                            INT channel, INT headerActive);

INT FDKsbrEnc_InitSbrCodeEnvelope(HANDLE_SBR_CODE_ENVELOPE h_sbrCodeEnvelope,
                                  INT *nSfb, INT deltaTAcrossFrames,
                                  FIXP_DBL dF_edge_1stEnv,
                                  FIXP_DBL dF_edge_incr);

INT FDKsbrEnc_InitSbrHuffmanTables(struct SBR_ENV_DATA *sbrEnvData,
                                   HANDLE_SBR_CODE_ENVELOPE henv,
                                   HANDLE_SBR_CODE_ENVELOPE hnoise,
                                   AMP_RES amp_res);

#endif
