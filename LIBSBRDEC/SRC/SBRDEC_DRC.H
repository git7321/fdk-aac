/**************************** SBR decoder library ******************************

   Author(s):   Christian Griebel

   Description: Dynamic range control (DRC) decoder tool for SBR

*******************************************************************************/

#ifndef SBRDEC_DRC_H
#define SBRDEC_DRC_H

#include "sbrdecoder.h"

#define SBRDEC_MAX_DRC_CHANNELS (8)
#define SBRDEC_MAX_DRC_BANDS (16)

typedef struct {
  FIXP_DBL prevFact_mag[(64)];
  INT prevFact_exp;

  FIXP_DBL currFact_mag[SBRDEC_MAX_DRC_BANDS];
  FIXP_DBL nextFact_mag[SBRDEC_MAX_DRC_BANDS];
  INT currFact_exp;
  INT nextFact_exp;

  UINT numBandsCurr;
  UINT numBandsNext;
  USHORT bandTopCurr[SBRDEC_MAX_DRC_BANDS];
  USHORT bandTopNext[SBRDEC_MAX_DRC_BANDS];

  SHORT drcInterpolationSchemeCurr;
  SHORT drcInterpolationSchemeNext;

  SHORT enable;

  UCHAR winSequenceCurr;
  UCHAR winSequenceNext;

} SBRDEC_DRC_CHANNEL;

typedef SBRDEC_DRC_CHANNEL *HANDLE_SBR_DRC_CHANNEL;

void sbrDecoder_drcInitChannel(HANDLE_SBR_DRC_CHANNEL hDrcData);

void sbrDecoder_drcUpdateChannel(HANDLE_SBR_DRC_CHANNEL hDrcData);

void sbrDecoder_drcApplySlot(HANDLE_SBR_DRC_CHANNEL hDrcData,
                             FIXP_DBL *qmfRealSlot, FIXP_DBL *qmfImagSlot,
                             int col, int numQmfSubSamples, int maxShift);

void sbrDecoder_drcApply(HANDLE_SBR_DRC_CHANNEL hDrcData,
                         FIXP_DBL **QmfBufferReal, FIXP_DBL **QmfBufferImag,
                         int numQmfSubSamples, int *scaleFactor);

#endif /* SBRDEC_DRC_H */
