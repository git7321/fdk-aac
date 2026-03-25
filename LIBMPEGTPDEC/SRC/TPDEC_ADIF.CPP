/******************* MPEG transport format decoder library *********************

   Author(s):   Josef Hoepfl

   Description: ADIF reader

*******************************************************************************/

#include "tpdec_adif.h"

#include "FDK_bitstream.h"
#include "genericStds.h"

TRANSPORTDEC_ERROR adifRead_DecodeHeader(CAdifHeader *pAdifHeader,
                                         CProgramConfig *pPce,
                                         HANDLE_FDK_BITSTREAM bs) {
  int i;
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;
  UINT startAnchor = FDKgetValidBits(bs);

  if ((INT)startAnchor < MIN_ADIF_HEADERLENGTH) {
    return (TRANSPORTDEC_NOT_ENOUGH_BITS);
  }

  if (FDKreadBits(bs, 8) != 'A') {
    return (TRANSPORTDEC_SYNC_ERROR);
  }
  if (FDKreadBits(bs, 8) != 'D') {
    return (TRANSPORTDEC_SYNC_ERROR);
  }
  if (FDKreadBits(bs, 8) != 'I') {
    return (TRANSPORTDEC_SYNC_ERROR);
  }
  if (FDKreadBits(bs, 8) != 'F') {
    return (TRANSPORTDEC_SYNC_ERROR);
  }

  if ((pAdifHeader->CopyrightIdPresent = FDKreadBits(bs, 1)) != 0)
    FDKpushBiDirectional(bs, 72);

  pAdifHeader->OriginalCopy = FDKreadBits(bs, 1);
  pAdifHeader->Home = FDKreadBits(bs, 1);
  pAdifHeader->BitstreamType = FDKreadBits(bs, 1);

  pAdifHeader->BitRate = FDKreadBits(bs, 16);
  pAdifHeader->BitRate <<= 7;
  pAdifHeader->BitRate |= FDKreadBits(bs, 7);

  pAdifHeader->NumProgramConfigElements = FDKreadBits(bs, 4) + 1;

  if (pAdifHeader->BitstreamType == 0) {
    FDKpushBiDirectional(bs, 20);
  }

  for (i = 0; i < pAdifHeader->NumProgramConfigElements; i++) {
    CProgramConfig_Read(pPce, bs, startAnchor);
  }

  FDKbyteAlign(bs, startAnchor);

  return (ErrorStatus);
}
