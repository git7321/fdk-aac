/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "env_bit.h"
#include "cmondata.h"

#ifndef min
#define min(a, b) (a < b ? a : b)
#endif

#ifndef max
#define max(a, b) (a > b ? a : b)
#endif

static void crcAdvance(USHORT crcPoly, USHORT crcMask, USHORT *crc,
                       ULONG bValue, INT bBits) {
  INT i;
  USHORT flag;

  for (i = bBits - 1; i >= 0; i--) {
    flag = ((*crc) & crcMask) ? (1) : (0);
    flag ^= (bValue & (1 << i)) ? (1) : (0);

    (*crc) <<= 1;
    if (flag) (*crc) ^= crcPoly;
  }
}

INT FDKsbrEnc_InitSbrBitstream(
    HANDLE_COMMON_DATA hCmonData,
    UCHAR *memoryBase,
    INT memorySize,
    HANDLE_FDK_CRCINFO hCrcInfo, UINT sbrSyntaxFlags)
{
  INT crcRegion = 0;

  FDKresetBitbuffer(&hCmonData->sbrBitbuf, BS_WRITER);

  FDKinitBitStream(&hCmonData->tmpWriteBitbuf, memoryBase, memorySize, 0,
                   BS_WRITER);

  if (sbrSyntaxFlags & SBR_SYNTAX_CRC) {
    if (sbrSyntaxFlags & SBR_SYNTAX_DRM_CRC) {
      FDKwriteBits(&hCmonData->sbrBitbuf, 0x0, SI_SBR_DRM_CRC_BITS);
      FDKcrcInit(hCrcInfo, 0x001d, 0xFFFF, SI_SBR_DRM_CRC_BITS);
      crcRegion = FDKcrcStartReg(hCrcInfo, &hCmonData->sbrBitbuf, 0);
    } else {
      FDKwriteBits(&hCmonData->sbrBitbuf, 0x0, SI_SBR_CRC_BITS);
    }
  }

  return (crcRegion);
}

void FDKsbrEnc_AssembleSbrBitstream(HANDLE_COMMON_DATA hCmonData,
                                    HANDLE_FDK_CRCINFO hCrcInfo, INT crcRegion,
                                    UINT sbrSyntaxFlags) {
  USHORT crcReg = SBR_CRCINIT;
  INT numCrcBits, i;

  if (hCmonData == NULL) return;

  hCmonData->sbrFillBits = 0;

  if (sbrSyntaxFlags & SBR_SYNTAX_DRM_CRC) {
    FDKcrcEndReg(hCrcInfo, &hCmonData->sbrBitbuf, crcRegion);
    FDKwriteBits(&hCmonData->tmpWriteBitbuf, FDKcrcGetCRC(hCrcInfo) ^ 0xFF,
                 SI_SBR_DRM_CRC_BITS);
  } else {
    if (!(sbrSyntaxFlags & SBR_SYNTAX_LOW_DELAY)) {
      int sbrLoad = hCmonData->sbrHdrBits + hCmonData->sbrDataBits;

      if (sbrSyntaxFlags & SBR_SYNTAX_CRC) {
        sbrLoad += SI_SBR_CRC_BITS;
      }

      sbrLoad += 4;

      hCmonData->sbrFillBits = (8 - (sbrLoad % 8)) % 8;

      FDKwriteBits(&hCmonData->sbrBitbuf, 0, hCmonData->sbrFillBits);

      FDK_ASSERT(FDKgetValidBits(&hCmonData->sbrBitbuf) % 8 == 4);
    }

    if (sbrSyntaxFlags & SBR_SYNTAX_CRC) {
      FDK_BITSTREAM tmpCRCBuf = hCmonData->sbrBitbuf;
      FDKresetBitbuffer(&tmpCRCBuf, BS_READER);

      numCrcBits = hCmonData->sbrHdrBits + hCmonData->sbrDataBits +
                   hCmonData->sbrFillBits;

      for (i = 0; i < numCrcBits; i++) {
        INT bit;
        bit = FDKreadBits(&tmpCRCBuf, 1);
        crcAdvance(SBR_CRC_POLY, SBR_CRC_MASK, &crcReg, bit, 1);
      }
      crcReg &= (SBR_CRC_RANGE);

      FDKwriteBits(&hCmonData->tmpWriteBitbuf, crcReg, SI_SBR_CRC_BITS);
    }
  }

  FDKsyncCache(&hCmonData->tmpWriteBitbuf);
}
