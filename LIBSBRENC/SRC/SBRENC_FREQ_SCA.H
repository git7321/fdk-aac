/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef SBRENC_FREQ_SCA_H
#define SBRENC_FREQ_SCA_H

#include "sbr_encoder.h"
#include "sbr_def.h"

#define MAX_OCTAVE 29
#define MAX_SECOND_REGION 50

INT FDKsbrEnc_UpdateFreqScale(UCHAR *v_k_master, INT *h_num_bands, const INT k0,
                              const INT k2, const INT freq_scale,
                              const INT alter_scale);

INT FDKsbrEnc_UpdateHiRes(UCHAR *h_hires, INT *num_hires, UCHAR *v_k_master,
                          INT num_master, INT *xover_band);

void FDKsbrEnc_UpdateLoRes(UCHAR *v_lores, INT *num_lores, UCHAR *v_hires,
                           INT num_hires);

INT FDKsbrEnc_FindStartAndStopBand(const INT srSbr, const INT srCore,
                                   const INT noChannels, const INT startFreq,
                                   const INT stop_freq, INT *k0, INT *k2);

INT FDKsbrEnc_getSbrStartFreqRAW(INT startFreq, INT fsCore);
INT FDKsbrEnc_getSbrStopFreqRAW(INT stopFreq, INT fsCore);
#endif
