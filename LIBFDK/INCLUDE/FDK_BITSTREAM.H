/******************* Library for basic calculation routines ********************

   Author(s):   M. Lohwasser

   Description: bitstream interface to bitbuffer routines

*******************************************************************************/

#ifndef FDK_BITSTREAM_H
#define FDK_BITSTREAM_H

#include "FDK_bitbuffer.h"
#include "machine_type.h"
#include "genericStds.h"

#define CACHE_BITS 32

#define BUFSIZE_DUMMY_VALUE MAX_BUFSIZE_BYTES

typedef enum { BS_READER, BS_WRITER } FDK_BS_CFG;

typedef struct {
  UINT CacheWord;
  UINT BitsInCache;
  FDK_BITBUF hBitBuf;
  UINT ConfigCache;
} FDK_BITSTREAM;

typedef FDK_BITSTREAM *HANDLE_FDK_BITSTREAM;

FDK_INLINE
HANDLE_FDK_BITSTREAM FDKcreateBitStream(UCHAR *pBuffer, UINT bufSize,
                                        FDK_BS_CFG config = BS_READER) {
  HANDLE_FDK_BITSTREAM hBitStream =
      (HANDLE_FDK_BITSTREAM)FDKcalloc(1, sizeof(FDK_BITSTREAM));
  if (hBitStream == NULL) return NULL;
  FDK_InitBitBuffer(&hBitStream->hBitBuf, pBuffer, bufSize, 0);

  hBitStream->CacheWord = hBitStream->BitsInCache = 0;
  hBitStream->ConfigCache = config;

  return hBitStream;
}

FDK_INLINE
void FDKinitBitStream(HANDLE_FDK_BITSTREAM hBitStream, UCHAR *pBuffer,
                      UINT bufSize, UINT validBits,
                      FDK_BS_CFG config = BS_READER) {
  FDK_InitBitBuffer(&hBitStream->hBitBuf, pBuffer, bufSize, validBits);

  hBitStream->CacheWord = hBitStream->BitsInCache = 0;
  hBitStream->ConfigCache = config;
}

FDK_INLINE void FDKresetBitbuffer(HANDLE_FDK_BITSTREAM hBitStream,
                                  FDK_BS_CFG config = BS_READER) {
  FDK_ResetBitBuffer(&hBitStream->hBitBuf);

  hBitStream->CacheWord = hBitStream->BitsInCache = 0;
  hBitStream->ConfigCache = config;
}

FDK_INLINE void FDKdeleteBitStream(HANDLE_FDK_BITSTREAM hBitStream) {
  FDK_DeleteBitBuffer(&hBitStream->hBitBuf);
  FDKfree(hBitStream);
}

FDK_INLINE UINT FDKreadBits(HANDLE_FDK_BITSTREAM hBitStream,
                            const UINT numberOfBits) {
  UINT bits = 0;
  INT missingBits = (INT)numberOfBits - (INT)hBitStream->BitsInCache;

  FDK_ASSERT(numberOfBits <= 32);
  if (missingBits > 0) {
    if (missingBits != 32) bits = hBitStream->CacheWord << missingBits;
    hBitStream->CacheWord = FDK_get32(&hBitStream->hBitBuf);
    hBitStream->BitsInCache += CACHE_BITS;
  }

  hBitStream->BitsInCache -= numberOfBits;

  return (bits | (hBitStream->CacheWord >> hBitStream->BitsInCache)) &
         BitMask[numberOfBits];
}

FDK_INLINE UINT FDKreadBit(HANDLE_FDK_BITSTREAM hBitStream) {
  if (!hBitStream->BitsInCache) {
    hBitStream->CacheWord = FDK_get32(&hBitStream->hBitBuf);
    hBitStream->BitsInCache = CACHE_BITS - 1;
    return hBitStream->CacheWord >> 31;
  }
  hBitStream->BitsInCache--;

  return (hBitStream->CacheWord >> hBitStream->BitsInCache) & 1;
}

FDK_INLINE UINT FDKread2Bits(HANDLE_FDK_BITSTREAM hBitStream) {
  UINT bits = 0;
  INT missingBits = 2 - (INT)hBitStream->BitsInCache;
  if (missingBits > 0) {
    bits = hBitStream->CacheWord << missingBits;
    hBitStream->CacheWord = FDK_get32(&hBitStream->hBitBuf);
    hBitStream->BitsInCache += CACHE_BITS;
  }

  hBitStream->BitsInCache -= 2;

  return (bits | (hBitStream->CacheWord >> hBitStream->BitsInCache)) & 0x3;
}

FDK_INLINE UINT FDKreadBitsBwd(HANDLE_FDK_BITSTREAM hBitStream,
                               const UINT numberOfBits) {
  const UINT validMask = BitMask[numberOfBits];

  if (hBitStream->BitsInCache <= numberOfBits) {
    const INT freeBits = (CACHE_BITS - 1) - hBitStream->BitsInCache;

    hBitStream->CacheWord = (hBitStream->CacheWord << freeBits) |
                            FDK_getBwd(&hBitStream->hBitBuf, freeBits);
    hBitStream->BitsInCache += freeBits;
  }

  hBitStream->BitsInCache -= numberOfBits;

  return (hBitStream->CacheWord >> hBitStream->BitsInCache) & validMask;
}

FDK_INLINE UINT escapedValue(HANDLE_FDK_BITSTREAM hBitStream, int nBits1,
                             int nBits2, int nBits3) {
  UINT value = FDKreadBits(hBitStream, nBits1);

  if (value == (UINT)(1 << nBits1) - 1) {
    UINT valueAdd = FDKreadBits(hBitStream, nBits2);
    value += valueAdd;
    if (valueAdd == (UINT)(1 << nBits2) - 1) {
      value += FDKreadBits(hBitStream, nBits3);
    }
  }

  return value;
}

FDK_INLINE UINT FDKgetBits(HANDLE_FDK_BITSTREAM hBitStream, UINT numBits) {
  return FDK_get(&hBitStream->hBitBuf, numBits);
}

FDK_INLINE UCHAR FDKwriteBits(HANDLE_FDK_BITSTREAM hBitStream, UINT value,
                              const UINT numberOfBits) {
  const UINT validMask = BitMask[numberOfBits];

  if (hBitStream == NULL) {
    return numberOfBits;
  }

  if ((hBitStream->BitsInCache + numberOfBits) < CACHE_BITS) {
    hBitStream->BitsInCache += numberOfBits;
    hBitStream->CacheWord =
        (hBitStream->CacheWord << numberOfBits) | (value & validMask);
  } else {
    int missing_bits = CACHE_BITS - hBitStream->BitsInCache;
    int remaining_bits = numberOfBits - missing_bits;
    value = value & validMask;
    UINT CacheWord =
        (missing_bits == 32) ? 0 : (hBitStream->CacheWord << missing_bits);
    CacheWord |= (value >> (remaining_bits));
    FDK_put(&hBitStream->hBitBuf, CacheWord, 32);

    hBitStream->CacheWord = value;
    hBitStream->BitsInCache = remaining_bits;
  }

  return numberOfBits;
}

FDK_INLINE UCHAR FDKwriteBitsBwd(HANDLE_FDK_BITSTREAM hBitStream, UINT value,
                                 const UINT numberOfBits) {
  const UINT validMask = BitMask[numberOfBits];

  if ((hBitStream->BitsInCache + numberOfBits) <= CACHE_BITS) {
    hBitStream->BitsInCache += numberOfBits;
    hBitStream->CacheWord =
        (hBitStream->CacheWord << numberOfBits) | (value & validMask);
  } else {
    FDK_putBwd(&hBitStream->hBitBuf, hBitStream->CacheWord,
               hBitStream->BitsInCache);
    hBitStream->BitsInCache = numberOfBits;
    hBitStream->CacheWord = (value & validMask);
  }

  return numberOfBits;
}

FDK_INLINE UCHAR FDKwriteEscapedValue(HANDLE_FDK_BITSTREAM hBitStream,
                                      UINT value, UINT nBits1, UINT nBits2,
                                      UINT nBits3) {
  UCHAR nbits = 0;
  UINT tmp = (1 << nBits1) - 1;

  if (value < tmp) {
    nbits += FDKwriteBits(hBitStream, value, nBits1);
  } else {
    nbits += FDKwriteBits(hBitStream, tmp, nBits1);
    value -= tmp;
    tmp = (1 << nBits2) - 1;

    if (value < tmp) {
      nbits += FDKwriteBits(hBitStream, value, nBits2);
    } else {
      nbits += FDKwriteBits(hBitStream, tmp, nBits2);
      value -= tmp;

      nbits += FDKwriteBits(hBitStream, value, nBits3);
    }
  }

  return nbits;
}

FDK_INLINE void FDKsyncCache(HANDLE_FDK_BITSTREAM hBitStream) {
  if (hBitStream->ConfigCache == BS_READER)
    FDK_pushBack(&hBitStream->hBitBuf, hBitStream->BitsInCache,
                 hBitStream->ConfigCache);
  else if (hBitStream->BitsInCache)
    FDK_put(&hBitStream->hBitBuf, hBitStream->CacheWord,
            hBitStream->BitsInCache);

  hBitStream->BitsInCache = 0;
  hBitStream->CacheWord = 0;
}

FDK_INLINE void FDKsyncCacheBwd(HANDLE_FDK_BITSTREAM hBitStream) {
  if (hBitStream->ConfigCache == BS_READER) {
    FDK_pushForward(&hBitStream->hBitBuf, hBitStream->BitsInCache,
                    hBitStream->ConfigCache);
  } else {
    FDK_putBwd(&hBitStream->hBitBuf, hBitStream->CacheWord,
               hBitStream->BitsInCache);
  }

  hBitStream->BitsInCache = 0;
  hBitStream->CacheWord = 0;
}

FDK_INLINE void FDKbyteAlign(HANDLE_FDK_BITSTREAM hBitStream,
                             UINT alignmentAnchor) {
  FDKsyncCache(hBitStream);
  if (hBitStream->ConfigCache == BS_READER) {
    FDK_pushForward(
        &hBitStream->hBitBuf,
        (UINT)((INT)8 - (((INT)alignmentAnchor -
                          (INT)FDK_getValidBits(&hBitStream->hBitBuf)) &
                         0x07)) &
            0x07,
        hBitStream->ConfigCache);
  } else {
    FDK_put(&hBitStream->hBitBuf, 0,
            (8 - ((FDK_getValidBits(&hBitStream->hBitBuf) - alignmentAnchor) &
                  0x07)) &
                0x07);
  }
}

FDK_INLINE void FDKpushBackCache(HANDLE_FDK_BITSTREAM hBitStream,
                                 const UINT numberOfBits) {
  FDK_ASSERT((hBitStream->BitsInCache + numberOfBits) <= CACHE_BITS);
  hBitStream->BitsInCache += numberOfBits;
}

FDK_INLINE void FDKpushBack(HANDLE_FDK_BITSTREAM hBitStream,
                            const UINT numberOfBits) {
  if ((hBitStream->BitsInCache + numberOfBits) < CACHE_BITS &&
      (hBitStream->ConfigCache == BS_READER)) {
    hBitStream->BitsInCache += numberOfBits;
    FDKsyncCache(hBitStream);
  } else {
    FDKsyncCache(hBitStream);
    FDK_pushBack(&hBitStream->hBitBuf, numberOfBits, hBitStream->ConfigCache);
  }
}

FDK_INLINE void FDKpushFor(HANDLE_FDK_BITSTREAM hBitStream,
                           const UINT numberOfBits) {
  if ((hBitStream->BitsInCache > numberOfBits) &&
      (hBitStream->ConfigCache == BS_READER)) {
    hBitStream->BitsInCache -= numberOfBits;
  } else {
    FDKsyncCache(hBitStream);
    FDK_pushForward(&hBitStream->hBitBuf, numberOfBits,
                    hBitStream->ConfigCache);
  }
}

FDK_INLINE void FDKpushBiDirectional(HANDLE_FDK_BITSTREAM hBitStream,
                                     const INT numberOfBits) {
  if (numberOfBits >= 0)
    FDKpushFor(hBitStream, numberOfBits);
  else
    FDKpushBack(hBitStream, -numberOfBits);
}

FDK_INLINE UINT FDKgetValidBits(HANDLE_FDK_BITSTREAM hBitStream) {
  FDKsyncCache(hBitStream);
  return FDK_getValidBits(&hBitStream->hBitBuf);
}

FDK_INLINE INT FDKgetFreeBits(HANDLE_FDK_BITSTREAM hBitStream) {
  return FDK_getFreeBits(&hBitStream->hBitBuf);
}

FDK_INLINE void FDKfeedBuffer(HANDLE_FDK_BITSTREAM hBitStream,
                              const UCHAR inputBuffer[], const UINT bufferSize,
                              UINT *bytesValid) {
  FDKsyncCache(hBitStream);
  FDK_Feed(&hBitStream->hBitBuf, inputBuffer, bufferSize, bytesValid);
}

FDK_INLINE void FDKcopyBuffer(HANDLE_FDK_BITSTREAM hBSDst,
                              HANDLE_FDK_BITSTREAM hBSSrc, UINT *bytesValid) {
  FDKsyncCache(hBSSrc);
  FDK_Copy(&hBSDst->hBitBuf, &hBSSrc->hBitBuf, bytesValid);
}

FDK_INLINE void FDKfetchBuffer(HANDLE_FDK_BITSTREAM hBitStream,
                               UCHAR *outputBuffer, UINT *writeBytes) {
  FDKsyncCache(hBitStream);
  FDK_Fetch(&hBitStream->hBitBuf, outputBuffer, writeBytes);
}

#endif
