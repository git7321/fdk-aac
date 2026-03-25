/**************************** AAC encoder library ******************************

   Author(s):   Tobias Chalupka

   Description: FDKaacLdEnc_MdctTransform480:
                The module FDKaacLdEnc_MdctTransform will perform the MDCT.
                The MDCT supports the sine window and
                the zero padded window. The algorithm of the MDCT
                can be divided in  Windowing, PreModulation, Fft and
                PostModulation.

*******************************************************************************/

#include "transform.h"
#include "dct.h"
#include "psy_const.h"
#include "aacEnc_rom.h"
#include "FDK_tools_rom.h"

#if defined(__arm__)
#endif

INT FDKaacEnc_Transform_Real(const INT_PCM *pTimeData,
                             FIXP_DBL *RESTRICT mdctData, const INT blockType,
                             const INT windowShape, INT *prevWindowShape,
                             H_MDCT mdctPers, const INT frameLength,
                             INT *pMdctData_e, INT filterType) {
  const INT_PCM *RESTRICT timeData;

  UINT numSpec;
  UINT numMdctLines;
  UINT offset;
  int fr;
  SHORT mdctData_e[8];

  timeData = pTimeData;

  if (blockType == SHORT_WINDOW) {
    numSpec = 8;
    numMdctLines = frameLength >> 3;
  } else {
    numSpec = 1;
    numMdctLines = frameLength;
  }

  offset = (windowShape == LOL_WINDOW) ? ((frameLength * 3) >> 2) : 0;
  switch (blockType) {
    case LONG_WINDOW:
    case STOP_WINDOW:
      fr = frameLength - offset;
      break;
    case START_WINDOW:
    case SHORT_WINDOW:
      fr = frameLength >> 3;
      break;
    default:
      FDK_ASSERT(0);
      return -1;
  }

  mdct_block(mdctPers, timeData, frameLength, mdctData, numSpec, numMdctLines,
             FDKgetWindowSlope(fr, windowShape), fr, mdctData_e);

  if (blockType == SHORT_WINDOW) {
    if (!(mdctData_e[0] == mdctData_e[1] && mdctData_e[1] == mdctData_e[2] &&
          mdctData_e[2] == mdctData_e[3] && mdctData_e[3] == mdctData_e[4] &&
          mdctData_e[4] == mdctData_e[5] && mdctData_e[5] == mdctData_e[6] &&
          mdctData_e[6] == mdctData_e[7])) {
      return -1;
    }
  }
  *prevWindowShape = windowShape;
  *pMdctData_e = mdctData_e[0];

  return 0;
}

INT FDKaacEnc_Transform_Real_Eld(const INT_PCM *pTimeData,
                                 FIXP_DBL *RESTRICT mdctData,
                                 const INT blockType, const INT windowShape,
                                 INT *prevWindowShape, const INT frameLength,
                                 INT *mdctData_e, INT filterType,
                                 FIXP_DBL *RESTRICT overlapAddBuffer) {
  const INT_PCM *RESTRICT timeData;

  INT i;

  const FIXP_WTB *pWindowELD = NULL;
  int N = frameLength;
  int L = frameLength;

  timeData = pTimeData;

  if (blockType != LONG_WINDOW) {
    return -1;
  }

  *mdctData_e = 1 + 1;

  switch (frameLength) {
    case 512:
      pWindowELD = ELDAnalysis512;
      break;
    case 480:
      pWindowELD = ELDAnalysis480;
      break;
    case 256:
      pWindowELD = ELDAnalysis256;
      *mdctData_e += 1;
      break;
    case 240:
      pWindowELD = ELDAnalysis240;
      *mdctData_e += 1;
      break;
    case 128:
      pWindowELD = ELDAnalysis128;
      *mdctData_e += 2;
      break;
    case 120:
      pWindowELD = ELDAnalysis120;
      *mdctData_e += 2;
      break;
    default:
      FDK_ASSERT(0);
      return -1;
  }

  for (i = 0; i < N / 4; i++) {
    FIXP_DBL z0, outval;

    z0 = (fMult((FIXP_PCM)timeData[L + N * 3 / 4 - 1 - i],
                pWindowELD[N / 2 - 1 - i])
          << (WTS0 - 1)) +
         (fMult((FIXP_PCM)timeData[L + N * 3 / 4 + i], pWindowELD[N / 2 + i])
          << (WTS0 - 1));

    outval = (fMultDiv2((FIXP_PCM)timeData[L + N * 3 / 4 - 1 - i],
                        pWindowELD[N + N / 2 - 1 - i]) >>
              (-WTS1));
    outval += (fMultDiv2((FIXP_PCM)timeData[L + N * 3 / 4 + i],
                         pWindowELD[N + N / 2 + i]) >>
               (-WTS1));
    outval += (fMultDiv2(overlapAddBuffer[N / 2 + i], pWindowELD[2 * N + i]) >>
               (-WTS2 - 1));

    overlapAddBuffer[N / 2 + i] = overlapAddBuffer[i];

    overlapAddBuffer[i] = z0;
    mdctData[i] = overlapAddBuffer[N / 2 + i] +
                  (fMultDiv2(overlapAddBuffer[N + N / 2 - 1 - i],
                             pWindowELD[2 * N + N / 2 + i]) >>
                   (-WTS2 - 1));

    mdctData[N - 1 - i] = outval;
    overlapAddBuffer[N + N / 2 - 1 - i] = outval;
  }

  for (i = N / 4; i < N / 2; i++) {
    FIXP_DBL z0, outval;

    z0 = fMult((FIXP_PCM)timeData[L + N * 3 / 4 - 1 - i],
               pWindowELD[N / 2 - 1 - i])
         << (WTS0 - 1);

    outval = (fMultDiv2((FIXP_PCM)timeData[L + N * 3 / 4 - 1 - i],
                        pWindowELD[N + N / 2 - 1 - i]) >>
              (-WTS1));
    outval += (fMultDiv2(overlapAddBuffer[N / 2 + i], pWindowELD[2 * N + i]) >>
               (-WTS2 - 1));

    overlapAddBuffer[N / 2 + i] =
        overlapAddBuffer[i] +
        (fMult((FIXP_PCM)timeData[L - N / 4 + i], pWindowELD[N / 2 + i])
         << (WTS0 - 1));

    overlapAddBuffer[i] = z0;
    mdctData[i] = overlapAddBuffer[N / 2 + i] +
                  (fMultDiv2(overlapAddBuffer[N + N / 2 - 1 - i],
                             pWindowELD[2 * N + N / 2 + i]) >>
                   (-WTS2 - 1));

    mdctData[N - 1 - i] = outval;
    overlapAddBuffer[N + N / 2 - 1 - i] = outval;
  }
  dct_IV(mdctData, frameLength, mdctData_e);

  *prevWindowShape = windowShape;

  return 0;
}
