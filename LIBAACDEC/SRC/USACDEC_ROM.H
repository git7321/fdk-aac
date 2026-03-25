/**************************** AAC decoder library ******************************

   Author(s):   M. Jander

   Description: re8.h

*******************************************************************************/

#ifndef USACDEC_ROM_H
#define USACDEC_ROM_H

#include "common_fix.h"
#include "FDK_lpc.h"
#include "usacdec_const.h"

#define NB_SPHERE 32
#define NB_LEADER 37
#define NB_LDSIGN 226
#define NB_LDQ3 9
#define NB_LDQ4 28

#define LSF_SCALE 13

extern const UINT fdk_dec_tab_factorial[8];
extern const UCHAR fdk_dec_Ia[NB_LEADER];
extern const UCHAR fdk_dec_Ds[NB_LDSIGN];
extern const USHORT fdk_dec_Is[NB_LDSIGN];
extern const UCHAR fdk_dec_Ns[], fdk_dec_A3[], fdk_dec_A4[];
extern const UCHAR fdk_dec_Da[][8];
extern const USHORT fdk_dec_I3[], fdk_dec_I4[];

extern const FIXP_LPC fdk_dec_lsf_init[16];
extern const FIXP_LPC fdk_dec_dico_lsf_abs_8b[16 * 256];

#define SF_QUA_GAIN7B 4
extern const FIXP_SGL fdk_t_qua_gain7b[128 * 2];
extern const FIXP_SGL lsp_interpol_factor[2][NB_SUBFR];

#define L_FILT 12

extern const FIXP_SGL fdk_dec_filt_lp[1 + L_FILT];

extern const FIXP_WTB FacWindowSynth128[128];
extern const FIXP_WTB FacWindowZir128[128];
extern const FIXP_WTB FacWindowSynth64[64];
extern const FIXP_WTB FacWindowZir64[64];
extern const FIXP_WTB FacWindowSynth32[32];
extern const FIXP_WTB FacWindowZir32[32];
extern const FIXP_WTB FacWindowSynth96[96];
extern const FIXP_WTB FacWindowZir96[96];
extern const FIXP_WTB FacWindowSynth48[48];
extern const FIXP_WTB FacWindowZir48[48];

#endif /* USACDEC_ROM_H */
