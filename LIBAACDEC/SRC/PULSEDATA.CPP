/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: pulse data tool

*******************************************************************************/

#include "pulsedata.h"

#include "channelinfo.h"

INT CPulseData_Read(HANDLE_FDK_BITSTREAM bs, CPulseData *const PulseData,
                    const SHORT *sfb_startlines, const void *pIcsInfo,
                    const SHORT frame_length) {
  int i, k = 0;
  const UINT MaxSfBands =
      GetScaleFactorBandsTransmitted((const CIcsInfo *)pIcsInfo);

  PulseData->PulseDataPresent = 0;

  if ((PulseData->PulseDataPresent = (UCHAR)FDKreadBit(bs)) != 0) {
    if (!IsLongBlock((const CIcsInfo *)pIcsInfo)) {
      return AAC_DEC_DECODE_FRAME_ERROR;
    }

    PulseData->NumberPulse = (UCHAR)FDKreadBits(bs, 2);
    PulseData->PulseStartBand = (UCHAR)FDKreadBits(bs, 6);

    if (PulseData->PulseStartBand >= MaxSfBands) {
      return AAC_DEC_DECODE_FRAME_ERROR;
    }

    k = sfb_startlines[PulseData->PulseStartBand];

    for (i = 0; i <= PulseData->NumberPulse; i++) {
      PulseData->PulseOffset[i] = (UCHAR)FDKreadBits(bs, 5);
      PulseData->PulseAmp[i] = (UCHAR)FDKreadBits(bs, 4);
      k += PulseData->PulseOffset[i];
    }

    if (k >= frame_length) {
      return AAC_DEC_DECODE_FRAME_ERROR;
    }
  }

  return 0;
}

void CPulseData_Apply(
    CPulseData *PulseData,
    const short
        *pScaleFactorBandOffsets,
    FIXP_DBL *coef)
{
  int i, k;

  if (PulseData->PulseDataPresent) {
    k = pScaleFactorBandOffsets[PulseData->PulseStartBand];

    for (i = 0; i <= PulseData->NumberPulse; i++) {
      k += PulseData->PulseOffset[i];
      if (coef[k] > (FIXP_DBL)0)
        coef[k] += (FIXP_DBL)(int)PulseData->PulseAmp[i];
      else
        coef[k] -= (FIXP_DBL)(int)PulseData->PulseAmp[i];
    }
  }
}
