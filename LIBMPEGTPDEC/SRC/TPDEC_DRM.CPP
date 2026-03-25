/******************* MPEG transport format decoder library *********************

   Author(s):   Christian Griebel

   Description: DRM transport stuff

*******************************************************************************/

#include "tpdec_drm.h"
#include "FDK_bitstream.h"

void drmRead_CrcInit(HANDLE_DRM pDrm)
{
  FDK_ASSERT(pDrm != NULL);

  FDKcrcInit(&pDrm->crcInfo, 0x001d, 0xFFFF, 8);
}

int drmRead_CrcStartReg(
    HANDLE_DRM pDrm,
    HANDLE_FDK_BITSTREAM hBs,
    int mBits
) {
  FDK_ASSERT(pDrm != NULL);

  FDKcrcReset(&pDrm->crcInfo);

  pDrm->crcReadValue = FDKreadBits(hBs, 8);

  return (FDKcrcStartReg(&pDrm->crcInfo, hBs, mBits));
}

void drmRead_CrcEndReg(
    HANDLE_DRM pDrm,
    HANDLE_FDK_BITSTREAM hBs,
    int reg
) {
  FDK_ASSERT(pDrm != NULL);

  FDKcrcEndReg(&pDrm->crcInfo, hBs, reg);
}

TRANSPORTDEC_ERROR drmRead_CrcCheck(HANDLE_DRM pDrm) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;
  USHORT crc;

  crc = FDKcrcGetCRC(&pDrm->crcInfo) ^ 0xFF;
  if (crc != pDrm->crcReadValue) {
    return (TRANSPORTDEC_CRC_ERROR);
  }

  return (ErrorStatus);
}
