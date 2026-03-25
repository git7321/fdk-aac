/*********************** MPEG surround encoder library *************************

   Author(s):   Karsten Linzmeier

   Description: Noiseless Coding
                Huffman encoder

*******************************************************************************/

#ifndef SACENC_NLC_ENC_H
#define SACENC_NLC_ENC_H

/* Includes ******************************************************************/
#include "sacenc_const.h"
#include "FDK_bitstream.h"
#include "sacenc_bitstream.h"

#define MAXBANDS MAX_NUM_BINS

typedef enum {
  t_CLD,
  t_ICC

} DATA_TYPE;

INT fdk_sacenc_ecDataPairEnc(HANDLE_FDK_BITSTREAM strm,
                             SHORT aaInData[][MAXBANDS],
                             SHORT aHistory[MAXBANDS],
                             const DATA_TYPE data_type, const INT setIdx,
                             const INT startBand, const INT dataBands,
                             const INT coarse_flag,
                             const INT independency_flag);

INT fdk_sacenc_ecDataSingleEnc(HANDLE_FDK_BITSTREAM strm,
                               SHORT aaInData[][MAXBANDS],
                               SHORT aHistory[MAXBANDS],
                               const DATA_TYPE data_type, const INT setIdx,
                               const INT startBand, const INT dataBands,
                               const INT coarse_flag,
                               const INT independency_flag);

#endif /* SACENC_NLC_ENC_H */
