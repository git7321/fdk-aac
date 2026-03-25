/**************************** AAC decoder library ******************************

   Author(s):   Matthias Hildenbrand, Manuel Jander

   Description: USAC LPC/AVQ decode

*******************************************************************************/

#ifndef USACDEC_LPC_H
#define USACDEC_LPC_H

#include "channelinfo.h"
#include "common_fix.h"
#include "FDK_bitstream.h"
#include "usacdec_rom.h"

#define LSPARG_SCALE 10

int CLpc_DecodeAVQ(HANDLE_FDK_BITSTREAM hBs, FIXP_DBL *lsfq, int nk_mode,
                   int nqn, int length);

int CLpc_Read(HANDLE_FDK_BITSTREAM hBs, FIXP_LPC lsp[][M_LP_FILTER_ORDER],
              FIXP_LPC lpc4_lsf[M_LP_FILTER_ORDER],
              FIXP_LPC lsf_adaptive_mean_cand[M_LP_FILTER_ORDER],
              FIXP_SGL pStability[], UCHAR *mod, int first_lpd_flag,
              int last_lpc_lost, int last_frame_ok);

void CLpc_Conceal(FIXP_LPC lsp[][M_LP_FILTER_ORDER],
                  FIXP_LPC lpc4_lsf[M_LP_FILTER_ORDER],
                  FIXP_LPC isf_adaptive_mean[M_LP_FILTER_ORDER],
                  const int first_lpd_flag);

void E_LPC_a_weight(FIXP_LPC *wA, const FIXP_LPC *A, const int m);
void CLpd_DecodeGain(FIXP_DBL *gain, INT *gain_e, int gain_code);
void E_LPC_f_lsp_a_conversion(FIXP_LPC *lsp, FIXP_LPC *a, INT *a_exp);

#endif /* USACDEC_LPC_H */
