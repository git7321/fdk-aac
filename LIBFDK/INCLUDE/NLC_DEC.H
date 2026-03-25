/******************* Library for basic calculation routines ********************

   Author(s):   Omer Osman

   Description: SAC/SAOC Dec Noiseless Coding

*******************************************************************************/

#ifndef NLC_DEC_H
#define NLC_DEC_H

#include "FDK_bitstream.h"
#include "huff_nodes.h"
#include "common_fix.h"

typedef enum {

  SAC_DECODER,
  SAOC_DECODER,
  USAC_DECODER

} DECODER_TYPE;

typedef enum {
  t_CLD,
  t_ICC,
  t_IPD,
  t_OLD,
  t_IOC,
  t_NRG,
  t_DCLD,
  t_DMG,
  t_PDG

} DATA_TYPE;

typedef enum {

  BACKWARDS = 0x0,
  FORWARDS = 0x1

} DIRECTION;

typedef enum {

  DIFF_FREQ = 0x0,
  DIFF_TIME = 0x1

} DIFF_TYPE;

typedef enum {

  HUFF_1D = 0x0,
  HUFF_2D = 0x1

} CODING_SCHEME;

typedef enum {

  FREQ_PAIR = 0x0,
  TIME_PAIR = 0x1

} PAIRING;

#ifndef HUFFDEC_PARAMS
#define HUFFDEC_PARMS

#define MAX_ENTRIES 168
#define HANDLE_HUFF_NODE const SHORT(*)[MAX_ENTRIES][2]

#endif /* HUFFDECPARAMS */

#define HUFFDEC_OK 0
#define HUFFDEC_NOTOK (-1)

typedef int ERROR_t;

ERROR_t EcDataPairDec(DECODER_TYPE DECODER, HANDLE_FDK_BITSTREAM strm,
                      SCHAR *aaOutData1, SCHAR *aaOutData2, SCHAR *aHistory,
                      DATA_TYPE data_type, int startBand, int dataBands,
                      int pair_flag, int coarse_flag,
                      int allowDiffTimeBack_flag);

ERROR_t huff_dec_reshape(HANDLE_FDK_BITSTREAM strm, int *out_data, int num_val);

extern ERROR_t sym_restoreIPD(HANDLE_FDK_BITSTREAM strm, int lav,
                              SCHAR data[2]);

#endif
