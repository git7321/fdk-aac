/*********************** MPEG surround decoder library *************************

   Author(s):   Matthias Hildenbrand

   Description: USAC MPS212 Transient Steering Decorrelator (TSD)

*******************************************************************************/

#ifndef SAC_TSD_H
#define SAC_TSD_H

#include "FDK_bitstream.h"
#include "common_fix.h"

#define MAX_TSD_TIME_SLOTS (64)

typedef struct {
  UCHAR bsTsdEnable;
  UCHAR numSlots;
  SCHAR
  bsTsdTrPhaseData[MAX_TSD_TIME_SLOTS];
} TSD_DATA;

FDK_INLINE int isTsdActive(const TSD_DATA *pTsdData) {
  return (int)pTsdData->bsTsdEnable;
}

int TsdRead(HANDLE_FDK_BITSTREAM hBs, const int numSlots, TSD_DATA *pTsdData);

void TsdGenerateNonTr(const int numHybridBands, const TSD_DATA *pTsdData,
                      const int ts, FIXP_DBL *pVdirectReal,
                      FIXP_DBL *pVdirectImag, FIXP_DBL *pVnonTrReal,
                      FIXP_DBL *pVnonTrImag, FIXP_DBL **ppDecorrInReal,
                      FIXP_DBL **ppDecorrInImag);
void TsdApply(const int numHybridBands, const TSD_DATA *pTsdData, int *pTsdTs,
              const FIXP_DBL *pVdirectReal, const FIXP_DBL *pVdirectImag,
              FIXP_DBL *pDnonTrReal, FIXP_DBL *pDnonTrImag);

#endif /* #ifndef SAC_TSD_H */
