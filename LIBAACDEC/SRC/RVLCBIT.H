/**************************** AAC decoder library ******************************

   Author(s):   Robert Weidner (DSP Solutions)

   Description: RVLC Decoder: Bitstream reading

*******************************************************************************/

#ifndef RVLCBIT_H
#define RVLCBIT_H

#include "rvlc.h"

UCHAR rvlcReadBitFromBitstream(HANDLE_FDK_BITSTREAM bs, const INT bsAnchor,
                               INT *pPosition, UCHAR readDirection);

#endif /* RVLCBIT_H */
