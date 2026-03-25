/**************************** AAC decoder library ******************************

   Author(s):   Robert Weidner (DSP Solutions)

   Description: HCR Decoder: Bitstream reading prototypes

*******************************************************************************/

#ifndef AACDEC_HCR_BIT_H
#define AACDEC_HCR_BIT_H

#include "aacdec_hcr.h"

UCHAR ToggleReadDirection(UCHAR readDirection);

UINT HcrGetABitFromBitstream(HANDLE_FDK_BITSTREAM bs, const INT bsAnchor,
                             INT *pLeftStartOfSegment,
                             INT *pRightStartOfSegment, UCHAR readDirection);

#endif /* AACDEC_HCR_BIT_H */
