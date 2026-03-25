/******************* MPEG transport format encoder library *********************

   Author(s):

   Description: ADIF Transport Headers writing

*******************************************************************************/

#include "tpenc_adif.h"

#include "tpenc_lib.h"
#include "tpenc_asc.h"

int adifWrite_EncodeHeader(ADIF_INFO *adif, HANDLE_FDK_BITSTREAM hBs,
                           INT adif_buffer_fullness) {
  const char adifId[5] = "ADIF";
  const int copyRightIdPresent = 0;
  const int originalCopy = 0;
  const int home = 0;
  int err = 0;

  int i;

  INT totalBitRate = adif->bitRate;

  if (adif->headerWritten) return 0;

  UINT alignAnchor = FDKgetValidBits(hBs);

  adif->bVariableRate = (adif_buffer_fullness >= (INT)(0x1 << 20)) ? 1 : 0;

  FDKwriteBits(hBs, adifId[0], 8);
  FDKwriteBits(hBs, adifId[1], 8);
  FDKwriteBits(hBs, adifId[2], 8);
  FDKwriteBits(hBs, adifId[3], 8);

  FDKwriteBits(hBs, copyRightIdPresent ? 1 : 0, 1);

  if (copyRightIdPresent) {
    for (i = 0; i < 72; i++) {
      FDKwriteBits(hBs, 0, 1);
    }
  }
  FDKwriteBits(hBs, originalCopy ? 1 : 0, 1);
  FDKwriteBits(hBs, home ? 1 : 0, 1);
  FDKwriteBits(hBs, adif->bVariableRate ? 1 : 0, 1);
  FDKwriteBits(hBs, totalBitRate, 23);

  FDKwriteBits(hBs, 0, 4);

  if (!adif->bVariableRate) {
    FDKwriteBits(hBs, adif_buffer_fullness, 20);
  }
  transportEnc_writePCE(hBs, adif->cm, adif->samplingRate, adif->instanceTag,
                        adif->profile, adif->matrixMixdownA,
                        (adif->pseudoSurroundEnable) ? 1 : 0, alignAnchor);

  return err;
}

int adifWrite_GetHeaderBits(ADIF_INFO *adif) {
  const int copyRightIdPresent = 0;

  if (adif->headerWritten) return 0;

  int bits = 0;

  bits += 8 * 4;

  bits += 1;

  if (copyRightIdPresent) bits += 72;

  bits += 26;

  bits += 4;

  if (!adif->bVariableRate) {
    bits += 20;
  }

  bits = transportEnc_GetPCEBits(adif->cm, adif->matrixMixdownA, bits);

  return bits;
}
