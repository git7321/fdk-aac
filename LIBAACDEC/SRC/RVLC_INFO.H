/**************************** AAC decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef RVLC_INFO_H
#define RVLC_INFO_H

#define FWD 0
#define BWD 1

#define MAX_RVL 7
#define MIN_RVL -7
#define MAX_ALLOWED_DPCM_INDEX \
  14
#define TABLE_OFFSET \
  7
#define MAX_LEN_RVLC_CODE_WORD 9
#define MAX_LEN_RVLC_ESCAPE_WORD \
  20

#define DPCM_NOISE_NRG_BITS 9
#define SF_OFFSET 100

#define CONCEAL_MAX_INIT 1311
#define CONCEAL_MIN_INIT -1311

#define RVLC_MAX_SFB ((8) * (16))

typedef struct {
  INT sf_concealment;
  INT rev_global_gain;
  SHORT length_of_rvlc_sf;
  INT dpcm_noise_nrg;
  INT sf_escapes_present;
  SHORT length_of_rvlc_escapes;
  INT dpcm_noise_last_position;

  INT dpcm_is_last_position;

  SHORT length_of_rvlc_sf_fwd;
  SHORT
  length_of_rvlc_sf_bwd;

  SHORT *pRvlBitCnt_RVL;
  INT *pBitstrIndxRvl_RVL;

  UCHAR numWindowGroups;
  UCHAR maxSfbTransmitted;
  UCHAR first_noise_group;
  UCHAR first_noise_band;
  UCHAR direction;

  INT bsAnchor;
  INT bitstreamIndexRvlFwd;
  INT bitstreamIndexRvlBwd;
  INT bitstreamIndexEsc;

  const UINT *pHuffTreeRvlCodewds;
  const UINT *pHuffTreeRvlcEscape;

  UCHAR numDecodedEscapeWordsFwd;
  UCHAR numDecodedEscapeWordsBwd;
  UCHAR numDecodedEscapeWordsEsc;

  SCHAR noise_used;
  SCHAR intensity_used;
  SCHAR sf_used;

  SHORT firstScf;
  SHORT lastScf;
  SHORT firstNrg;
  SHORT lastNrg;
  SHORT firstIs;
  SHORT lastIs;

  UINT errorLogRvlc;
  SHORT conceal_min;
  SHORT conceal_max;
  SHORT conceal_min_esc;
  SHORT conceal_max_esc;
} CErRvlcInfo;

typedef CErRvlcInfo RVLC_INFO;

#endif /* RVLC_INFO_H */
