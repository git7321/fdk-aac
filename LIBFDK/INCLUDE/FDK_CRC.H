/******************* Library for basic calculation routines ********************

   Author(s):

   Description: CRC calculation

*******************************************************************************/

#ifndef FDK_CRC_H
#define FDK_CRC_H

#include "FDK_bitstream.h"

#define MAX_CRC_REGS                                                           \
  3

typedef struct {
  UCHAR isActive;
  INT maxBits;
  INT bitBufCntBits;
  INT validBits;

} CCrcRegData;

typedef struct {
  CCrcRegData crcRegData[MAX_CRC_REGS];
  const USHORT*
      pCrcLookup;

  USHORT crcPoly;
  USHORT crcMask;
  USHORT startValue;
  UCHAR crcLen;

  UINT regStart;
  UINT regStop;

  USHORT crcValue;

} FDK_CRCINFO;

typedef FDK_CRCINFO* HANDLE_FDK_CRCINFO;

void FDKcrcInit(HANDLE_FDK_CRCINFO hCrcInfo, const UINT crcPoly,
                const UINT crcStartValue, const UINT crcLen);

void FDKcrcReset(HANDLE_FDK_CRCINFO hCrcInfo);

INT FDKcrcStartReg(HANDLE_FDK_CRCINFO hCrcInfo, const HANDLE_FDK_BITSTREAM hBs,
                   const INT mBits);

INT FDKcrcEndReg(HANDLE_FDK_CRCINFO hCrcInfo, const HANDLE_FDK_BITSTREAM hBs,
                 const INT reg);

USHORT FDKcrcGetCRC(const HANDLE_FDK_CRCINFO hCrcInfo);

#endif /* FDK_CRC_H */
