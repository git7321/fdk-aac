/**************************** AAC decoder library ******************************

   Author(s):   Robert Weidner (DSP Solutions)

   Description: HCR Decoder: Bitstream reading

*******************************************************************************/

#include "aacdec_hcr_bit.h"

UCHAR ToggleReadDirection(UCHAR readDirection) {
  if (readDirection == FROM_LEFT_TO_RIGHT) {
    return FROM_RIGHT_TO_LEFT;
  } else {
    return FROM_LEFT_TO_RIGHT;
  }
}

UINT HcrGetABitFromBitstream(HANDLE_FDK_BITSTREAM bs, const INT bsAnchor,
                             INT *pLeftStartOfSegment,
                             INT *pRightStartOfSegment, UCHAR readDirection) {
  UINT bit;
  INT readBitOffset;

  if (readDirection == FROM_LEFT_TO_RIGHT) {
    readBitOffset = (INT)FDKgetValidBits(bs) - bsAnchor + *pLeftStartOfSegment;
    if (readBitOffset) {
      FDKpushBiDirectional(bs, readBitOffset);
    }

    bit = FDKreadBits(bs, 1);

    *pLeftStartOfSegment += 1;
  } else {
    readBitOffset = (INT)FDKgetValidBits(bs) - bsAnchor + *pRightStartOfSegment;
    if (readBitOffset) {
      FDKpushBiDirectional(bs, readBitOffset);
    }

    bit = FDKreadBits(bs, 1);
    FDKpushBack(bs, 2);

    *pRightStartOfSegment -= 1;
  }

  return (bit);
}
