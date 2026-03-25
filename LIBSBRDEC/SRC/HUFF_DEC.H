/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef HUFF_DEC_H
#define HUFF_DEC_H

#include "sbrdecoder.h"
#include "FDK_bitstream.h"

typedef const SCHAR (*Huffman)[2];

int DecodeHuffmanCW(Huffman h, HANDLE_FDK_BITSTREAM hBitBuf);

#endif
