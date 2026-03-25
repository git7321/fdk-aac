/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef PSBITDEC_H
#define PSBITDEC_H

#include "sbrdecoder.h"
#include "psdec.h"

unsigned int ReadPsData(struct PS_DEC *h_ps_d, HANDLE_FDK_BITSTREAM hBs,
                        int nBitsLeft);

int DecodePs(struct PS_DEC *h_ps_d, const UCHAR frameError,
             PS_DEC_COEFFICIENTS *pCoef);

#endif /* PSBITDEC_H */
