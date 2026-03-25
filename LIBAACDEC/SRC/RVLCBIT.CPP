/**************************** AAC decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "rvlcbit.h"

UCHAR rvlcReadBitFromBitstream(HANDLE_FDK_BITSTREAM bs, const INT bsAnchor,
                               INT *pPosition, UCHAR readDirection) {
  UINT bit;
  INT readBitOffset = (INT)FDKgetValidBits(bs) - bsAnchor + *pPosition;

  if (readBitOffset) {
    FDKpushBiDirectional(bs, readBitOffset);
  }

  if (readDirection == FWD) {
    bit = FDKreadBits(bs, 1);

    *pPosition += 1;
  } else {
    bit = FDKreadBits(bs, 1);
    FDKpushBack(bs, 2);

    *pPosition -= 1;
  }

  return (bit);
}
