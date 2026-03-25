/**************************** AAC decoder library ******************************

   Author(s):   Robert Weidner (DSP Solutions)

   Description: HCR Decoder: Defines of state-constants, masks and
                state-prototypes

*******************************************************************************/

#ifndef AACDEC_HCRS_H
#define AACDEC_HCRS_H

#include "FDK_bitstream.h"
#include "aacdec_hcr_types.h"

#define STOP_THIS_STATE \
  0 /*                                                                      */
#define BODY_ONLY \
  1 /*   X          X                X                                      */
#define BODY_SIGN__BODY \
  2 /*   X          X         X      X [stop if no sign]                    */
#define BODY_SIGN__SIGN \
  3 /*              X                X [stop if sign bits decoded]          */
#define BODY_SIGN_ESC__BODY \
  4 /*   X          X         X      X [stop if no sign]                    */
#define BODY_SIGN_ESC__SIGN \
  5 /*              X         X      X [stop if no escape sequence]         */
#define BODY_SIGN_ESC__ESC_PREFIX \
  6 /*              X         X                                             */
#define BODY_SIGN_ESC__ESC_WORD \
  7 /*              X         X      X [stop if abs(second qsc) != 16]      */

#define MSB_31_MASK 0x80000000
#define DIMENSION_OF_ESCAPE_CODEBOOK 2
#define ESCAPE_CODEBOOK 11

#define MASK_ESCAPE_PREFIX_UP 0x000F0000
#define LSB_ESCAPE_PREFIX_UP 16

#define MASK_ESCAPE_PREFIX_DOWN 0x0000F000
#define LSB_ESCAPE_PREFIX_DOWN 12

#define MASK_ESCAPE_WORD 0x00000FFF
#define MASK_FLAG_A 0x00200000
#define MASK_FLAG_B 0x00100000

extern void DecodeNonPCWs(HANDLE_FDK_BITSTREAM bs, H_HCR_INFO hHcr);

UINT Hcr_State_BODY_ONLY(HANDLE_FDK_BITSTREAM, void*);
UINT Hcr_State_BODY_SIGN__BODY(HANDLE_FDK_BITSTREAM, void*);
UINT Hcr_State_BODY_SIGN__SIGN(HANDLE_FDK_BITSTREAM, void*);
UINT Hcr_State_BODY_SIGN_ESC__BODY(HANDLE_FDK_BITSTREAM, void*);
UINT Hcr_State_BODY_SIGN_ESC__SIGN(HANDLE_FDK_BITSTREAM, void*);
UINT Hcr_State_BODY_SIGN_ESC__ESC_PREFIX(HANDLE_FDK_BITSTREAM, void*);
UINT Hcr_State_BODY_SIGN_ESC__ESC_WORD(HANDLE_FDK_BITSTREAM, void*);

#endif /* AACDEC_HCRS_H */
