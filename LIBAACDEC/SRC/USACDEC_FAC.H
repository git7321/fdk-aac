/**************************** AAC decoder library ******************************

   Author(s):   Manuel Jander

   Description: USAC FAC

*******************************************************************************/

#ifndef USACDEC_FAC_H
#define USACDEC_FAC_H

#include "channelinfo.h"
#include "FDK_bitstream.h"

FIXP_DBL *CLpd_FAC_GetMemory(CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                             UCHAR mod[NB_SUBFR], int *pState);

int CLpd_FAC_Read(HANDLE_FDK_BITSTREAM hBs, FIXP_DBL *pFac, SCHAR *pFacScale,
                  int length, int use_gain, int frame);

void CFac_ApplyGains(FIXP_DBL fac_data[LFAC], const INT fac_length,
                     const FIXP_DBL tcx_gain, const FIXP_DBL alfd_gains[],
                     const INT mod);

INT CLpd_FAC_Mdct2Acelp(H_MDCT hMdct, FIXP_DBL *output, FIXP_DBL *pFac_data,
                        const int fac_data_e, FIXP_LPC *A, INT A_exp,
                        INT nrOutSamples, const INT fac_length,
                        const INT isFdFac, UCHAR prevWindowShape);

INT CLpd_FAC_Acelp2Mdct(H_MDCT hMdct, FIXP_DBL *output, FIXP_DBL *pSpec,
                        const SHORT spec_scale[], const int nSpec,
                        FIXP_DBL *pFac_data, const int fac_data_e,
                        const INT fac_length, INT nrSamples, const INT tl,
                        const FIXP_WTP *wrs, const INT fr, FIXP_LPC A[16],
                        INT A_exp, CAcelpStaticMem *acelp_mem,
                        const FIXP_DBL gain, const int last_frame_lost,
                        const int isFdFac, const UCHAR last_lpd, const int k,
                        int currAliasingSymmetry);

#endif /* USACDEC_FAC_H */
