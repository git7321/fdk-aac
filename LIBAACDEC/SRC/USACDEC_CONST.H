/**************************** AAC decoder library ******************************

   Author(s):   Manuel Jander

   Description: USAC related constants

*******************************************************************************/

#ifndef USACDEC_CONST_H
#define USACDEC_CONST_H

#define SF_CODE 6
#define SF_GAIN_C 16
#define SF_EXC 16
#define SF_GAIN_P 1
#define SF_PFAC 0
#define SF_SYNTH SF_EXC
#define SF_A_COEFFS 3
#define SF_STAB 1

#define M_LP_FILTER_ORDER 16
#define LP_FILTER_SCALE 4

#define PIT_MIN_12k8 34
#define PIT_MAX_12k8 231
#define FSCALE_DENOM 12800
#define FAC_FSCALE_MIN \
  6000

#if !defined(LPD_MAX_CORE_SR)
#define LPD_MAX_CORE_SR 24000
#endif
#define FAC_FSCALE_MAX \
  LPD_MAX_CORE_SR

#define PIT_MAX_TMP                                                            \
  (PIT_MAX_12k8 +                                                              \
   (6 *                                                                        \
    ((((FAC_FSCALE_MAX * PIT_MIN_12k8) + (FSCALE_DENOM / 2)) / FSCALE_DENOM) - \
     PIT_MIN_12k8)))
#if (PIT_MAX_TMP < \
     256)
#define PIT_MAX_MAX 256
#else
#define PIT_MAX_MAX PIT_MAX_TMP
#endif

#define NB_DIV 4
#define L_SUBFR 64
#define BPF_SFD 1
#define BPF_DELAY (BPF_SFD * L_SUBFR)

#define L_FILT 12
#define L_EXTRA 96
#define L_INTERPOL \
  (16 + 1)

#define L_FRAME_PLUS_1024 1024
#define L_DIV_1024 \
  (L_FRAME_PLUS_1024 / NB_DIV)
#define NB_SUBFR_1024 \
  (L_DIV_1024 / L_SUBFR)
#define NB_SUBFR_SUPERFR_1024 \
  (L_FRAME_PLUS_1024 / L_SUBFR)
#define AAC_SFD_1024 (NB_SUBFR_SUPERFR_1024 / 2)
#define AAC_DELAY_1024 (AAC_SFD_1024 * L_SUBFR)
#define SYN_SFD_1024 (AAC_SFD_1024 - BPF_SFD)
#define SYN_DELAY_1024 \
  (SYN_SFD_1024 * L_SUBFR)
#define LFAC_1024 (L_DIV_1024 / 2)
#define LFAC_SHORT_1024 \
  (L_DIV_1024 / 4)
#define FDNS_NPTS_1024 64

#define L_FRAME_PLUS_768 768
#define L_DIV_768 \
  (L_FRAME_PLUS_768 / NB_DIV)
#define NB_SUBFR_768 \
  (L_DIV_768 / L_SUBFR)
#define NB_SUBFR_SUPERFR_768 \
  (L_FRAME_PLUS_768 / L_SUBFR)
#define AAC_SFD_768 (NB_SUBFR_SUPERFR_768 / 2)
#define AAC_DELAY_768 (AAC_SFD_768 * L_SUBFR)
#define SYN_SFD_768 (AAC_SFD_768 - BPF_SFD)
#define SYN_DELAY_768 (SYN_SFD_768 * L_SUBFR)
#define LFAC_768 (L_DIV_768 / 2)
#define LFAC_SHORT_768 \
  (L_DIV_768 / 4)

#define L_FRAME_PLUS L_FRAME_PLUS_1024
#define L_DIV L_DIV_1024
#define NB_SUBFR NB_SUBFR_1024
#define NB_SUBFR_SUPERFR NB_SUBFR_SUPERFR_1024
#define AAC_SFD AAC_SFD_1024
#define AAC_DELAY AAC_DELAY_1024
#define SYN_SFD SYN_SFD_1024
#define SYN_DELAY SYN_DELAY_1024
#define LFAC LFAC_1024
#define LFAC_SHORT LFAC_SHORT_1024
#define FDNS_NPTS FDNS_NPTS_1024

#endif /* USACDEC_CONST_H */
