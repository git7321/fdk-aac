/******************* Library for basic calculation routines ********************

   Author(s):   M. Lohwasser

   Description: common bitbuffer read/write routines

*******************************************************************************/

#ifndef FDK_BITBUFFER_H
#define FDK_BITBUFFER_H

#include "FDK_archdef.h"
#include "machine_type.h"

#define MAX_BUFSIZE_BYTES (0x10000000)

typedef struct {
  UINT ValidBits;
  UINT ReadOffset;
  UINT WriteOffset;
  UINT BitNdx;

  UCHAR *Buffer;
  UINT bufSize;
  UINT bufBits;
} FDK_BITBUF;

typedef FDK_BITBUF *HANDLE_FDK_BITBUF;

#ifdef __cplusplus
extern "C" {
#endif

extern const UINT BitMask[32 + 1];

void FDK_CreateBitBuffer(HANDLE_FDK_BITBUF *hBitBuffer, UCHAR *pBuffer,
                         UINT bufSize);

void FDK_InitBitBuffer(HANDLE_FDK_BITBUF hBitBuffer, UCHAR *pBuffer,
                       UINT bufSize, UINT validBits);

void FDK_ResetBitBuffer(HANDLE_FDK_BITBUF hBitBuffer);

void FDK_DeleteBitBuffer(HANDLE_FDK_BITBUF hBitBuffer);

INT FDK_get(HANDLE_FDK_BITBUF hBitBuffer, const UINT numberOfBits);

INT FDK_get32(HANDLE_FDK_BITBUF hBitBuf);

void FDK_put(HANDLE_FDK_BITBUF hBitBuffer, UINT value, const UINT numberOfBits);

INT FDK_getBwd(HANDLE_FDK_BITBUF hBitBuffer, const UINT numberOfBits);
void FDK_putBwd(HANDLE_FDK_BITBUF hBitBuffer, UINT value,
                const UINT numberOfBits);

void FDK_pushBack(HANDLE_FDK_BITBUF hBitBuffer, const UINT numberOfBits,
                  UCHAR config);

void FDK_pushForward(HANDLE_FDK_BITBUF hBitBuffer, const UINT numberOfBits,
                     UCHAR config);

UINT FDK_getValidBits(HANDLE_FDK_BITBUF hBitBuffer);

INT FDK_getFreeBits(HANDLE_FDK_BITBUF hBitBuffer);

void FDK_Feed(HANDLE_FDK_BITBUF hBitBuffer, const UCHAR inputBuffer[],
              const UINT bufferSize, UINT *bytesValid);

void FDK_Copy(HANDLE_FDK_BITBUF hBitBufDst, HANDLE_FDK_BITBUF hBitBufSrc,
              UINT *bytesValid);

void FDK_Fetch(HANDLE_FDK_BITBUF hBitBuffer, UCHAR outBuf[], UINT *writeBytes);

#ifdef __cplusplus
}
#endif

#endif
