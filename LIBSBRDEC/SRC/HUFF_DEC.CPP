/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "huff_dec.h"

int DecodeHuffmanCW(Huffman h,
                    HANDLE_FDK_BITSTREAM hBs)
{
  SCHAR index = 0;
  int value, bit;

  while (index >= 0) {
    bit = FDKreadBits(hBs, 1);
    index = h[index][bit];
  }

  value = index + 64;

  return value;
}
