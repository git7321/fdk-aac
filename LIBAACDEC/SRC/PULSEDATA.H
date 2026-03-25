/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: pulse data tool

*******************************************************************************/

#ifndef PULSEDATA_H
#define PULSEDATA_H

#include "common_fix.h"
#include "FDK_bitstream.h"

#define N_MAX_LINES 4

typedef struct {
  UCHAR PulseDataPresent;
  UCHAR NumberPulse;
  UCHAR PulseStartBand;
  UCHAR PulseOffset[N_MAX_LINES];
  UCHAR PulseAmp[N_MAX_LINES];
} CPulseData;

INT CPulseData_Read(const HANDLE_FDK_BITSTREAM bs, CPulseData *const PulseData,
                    const SHORT *sfb_startlines, const void *pIcsInfo,
                    const SHORT frame_length);

void CPulseData_Apply(CPulseData *PulseData,
                      const short *pScaleFactorBandOffsets, FIXP_DBL *coef);

#endif /* #ifndef PULSEDATA_H */
