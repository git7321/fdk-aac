/**************************** AAC decoder library ******************************

   Author(s):   Matthias Hildenbrand

   Description: USAC ACELP frame decoder

*******************************************************************************/

#ifndef USACDEC_ACELP_H
#define USACDEC_ACELP_H

#include "common_fix.h"
#include "FDK_bitstream.h"
#include "usacdec_const.h"
#include "usacdec_rom.h"

typedef struct {
  FIXP_DBL old_exc_mem[PIT_MAX_MAX + L_INTERPOL];
  FIXP_DBL old_syn_mem[M_LP_FILTER_ORDER];
  FIXP_SGL A[M_LP_FILTER_ORDER];
  INT A_exp;
  FIXP_DBL gc_threshold;
  FIXP_DBL de_emph_mem;
  FIXP_SGL past_gpit;
  FIXP_DBL past_gcode;
  USHORT old_T0;
  UCHAR old_T0_frac;
  FIXP_DBL deemph_mem_wsyn;
  FIXP_DBL wsyn_rms;
  SHORT seed_ace;
} CAcelpStaticMem;

typedef struct {
  UCHAR
  acelp_core_mode;
  UCHAR mean_energy;
  USHORT T0[NB_SUBFR];
  UCHAR T0_frac[NB_SUBFR];
  UCHAR ltp_filtering_flag[NB_SUBFR];
  SHORT icb_index[NB_SUBFR]
                 [8];
  UCHAR gains[NB_SUBFR];
} CAcelpChannelData;

INT CLpd_AcelpRead(HANDLE_FDK_BITSTREAM hBs, CAcelpChannelData *acelpData,
                   INT acelp_core_mode, INT i_offset, INT coreCoderFrameLength);

void Acelp_PreProcessing(FIXP_DBL *synth_buf, FIXP_DBL *old_synth, INT *pitch,
                         INT *old_T_pf, FIXP_DBL *pit_gain,
                         FIXP_DBL *old_gain_pf, INT samplingRate, INT *i_offset,
                         INT coreCoderFrameLength, INT synSfd,
                         INT nbSubfrSuperfr);

void Acelp_PostProcessing(FIXP_DBL *synth_buf, FIXP_DBL *old_synth, INT *pitch,
                          INT *old_T_pf, INT coreCoderFrameLength, INT synSfd,
                          INT nbSubfrSuperfr);

void CLpd_AcelpDecode(CAcelpStaticMem *acelp_mem, INT i_offset,
                      const FIXP_LPC lsp_old[M_LP_FILTER_ORDER],
                      const FIXP_LPC lsp_new[M_LP_FILTER_ORDER],
                      FIXP_SGL stab_fac, CAcelpChannelData *acelpData,
                      INT numLostSubframes, int lastLpcLost, int frameCnt,
                      FIXP_DBL synth[], int pT[], FIXP_DBL *pit_gain,
                      INT coreCoderFrameLength);

void CLpd_AcelpReset(CAcelpStaticMem *acelp_mem);

void CLpd_AcelpPrepareInternalMem(const FIXP_DBL *synth, UCHAR last_lpd_mode,
                                  UCHAR last_last_lpd_mode,
                                  const FIXP_LPC *A_new, const INT A_new_exp,
                                  const FIXP_LPC *A_old, const INT A_old_exp,
                                  CAcelpStaticMem *acelp_mem,
                                  INT coreCoderFrameLength, INT clearOldExc,
                                  UCHAR lpd_mode);

void CLpd_Acelp_Zir(const FIXP_LPC A[], const INT A_exp,
                    CAcelpStaticMem *acelp_mem, const INT length,
                    FIXP_DBL zir[], int doDeemph);

FIXP_DBL *CLpd_ACELP_GetFreeExcMem(CAcelpStaticMem *acelp_mem, INT length);

void CLpd_TcxTDConceal(CAcelpStaticMem *acelp_mem, SHORT *pitch,
                       const FIXP_LPC lsp_old[M_LP_FILTER_ORDER],
                       const FIXP_LPC lsp_new[M_LP_FILTER_ORDER],
                       const FIXP_SGL stab_fac, INT numLostSubframes,
                       FIXP_DBL synth[], INT coreCoderFrameLength,
                       UCHAR last_tcx_noise_factor);

inline SHORT E_UTIL_random(SHORT *seed) {
  *seed = (SHORT)((((LONG)*seed * (LONG)31821) >> 1) + (LONG)13849);
  return (*seed);
}

#endif /* USACDEC_ACELP_H */
