/******************* MPEG transport format decoder library *********************

   Author(s):   Josef Hoepfl

   Description: ADIF reader

*******************************************************************************/

#ifndef TPDEC_ADIF_H
#define TPDEC_ADIF_H

#include "tpdec_lib.h"

#define MIN_ADIF_HEADERLENGTH 63

typedef struct {
  INT NumProgramConfigElements;
  UINT BitRate;
  UCHAR CopyrightIdPresent;
  UCHAR OriginalCopy;
  UCHAR Home;
  UCHAR BitstreamType;
} CAdifHeader;

TRANSPORTDEC_ERROR adifRead_DecodeHeader(CAdifHeader *pAdifHeader,
                                         CProgramConfig *pPce,
                                         HANDLE_FDK_BITSTREAM bs);

#endif /* TPDEC_ADIF_H */
