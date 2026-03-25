/*********************** MPEG surround encoder library *************************

   Author(s):   Christian Goettlinger

   Description: Encoder Library Interface
                delay management of the encoder

*******************************************************************************/

#include "sacenc_delay.h"
#include "sacenc_const.h"
#include "FDK_matrixCalloc.h"

struct DELAY {
  struct DELAY_CONFIG {
    INT bDmxAlign;
    INT bTimeDomDmx;
    INT bMinimizeDelay;
    INT bSacTimeAlignmentDynamicOut;

    INT nQmfLen;
    INT nFrameLen;
    INT nSurroundDelay;
    INT nArbDmxDelay;
    INT nLimiterDelay;
    INT nCoreCoderDelay;
    INT nSacStreamMuxDelay;
    INT nSacTimeAlignment;
  } config;

  INT nDmxAlignBuffer;
  INT nSurroundAnalysisBuffer;
  INT nArbDmxAnalysisBuffer;
  INT nOutputAudioBuffer;
  INT nBitstreamFrameBuffer;
  INT nOutputAudioQmfFrameBuffer;
  INT nDiscardOutFrames;

  INT nBitstreamFrameBufferSize;

  INT nInfoDmxDelay;
  INT nInfoCodecDelay;
  INT nInfoDecoderDelay;
};

FDK_SACENC_ERROR fdk_sacenc_delay_Open(HANDLE_DELAY *phDelay) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == phDelay) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDK_ALLOCATE_MEMORY_1D(*phDelay, 1, struct DELAY);
  }
  return error;

bail:
  fdk_sacenc_delay_Close(phDelay);
  return ((SACENC_OK == error) ? SACENC_MEMORY_ERROR : error);
}

FDK_SACENC_ERROR fdk_sacenc_delay_Close(HANDLE_DELAY *phDelay) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == phDelay) {
    error = SACENC_INVALID_HANDLE;
  } else {
    if (NULL != *phDelay) {
      FDK_FREE_MEMORY_1D(*phDelay);
    }
  }
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_delay_Init(HANDLE_DELAY hDelay, const INT nQmfLen,
                                       const INT nFrameLen,
                                       const INT nCoreCoderDelay,
                                       const INT nSacStreamMuxDelay) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == hDelay) {
    error = SACENC_INVALID_HANDLE;
  } else {
    FDKmemclear(&hDelay->config, sizeof(hDelay->config));

    hDelay->config.nQmfLen = nQmfLen;
    hDelay->config.nFrameLen = nFrameLen;
    hDelay->config.nCoreCoderDelay = nCoreCoderDelay;
    hDelay->config.nSacStreamMuxDelay = nSacStreamMuxDelay;
  }
  return error;
}

FDK_SACENC_ERROR fdk_sacenc_delay_SubCalulateBufferDelays(HANDLE_DELAY hDel) {
  FDK_SACENC_ERROR error = SACENC_OK;

  if (NULL == hDel) {
    error = SACENC_INVALID_HANDLE;
  } else {
    int nEncoderAnDelay, nEncoderSynDelay, nEncoderWinDelay, nDecoderAnDelay,
        nDecoderSynDelay, nResidualCoderFrameDelay,
        nArbDmxResidualCoderFrameDelay;

    if (hDel->config.bSacTimeAlignmentDynamicOut > 0) {
      hDel->config.nSacTimeAlignment = 0;
    }

    {
      nEncoderAnDelay =
          2 * hDel->config.nQmfLen +
          hDel->config.nQmfLen / 2;
      nEncoderSynDelay = 1 * hDel->config.nQmfLen + hDel->config.nQmfLen / 2;
      nDecoderAnDelay = 2 * hDel->config.nQmfLen + hDel->config.nQmfLen / 2;
      nDecoderSynDelay = 1 * hDel->config.nQmfLen + hDel->config.nQmfLen / 2;
      nEncoderWinDelay =
          hDel->config.nFrameLen / 2;
    }

    { nResidualCoderFrameDelay = 0; }

    { nArbDmxResidualCoderFrameDelay = 0; }

    if (hDel->config.bTimeDomDmx == 0) {
      int tempDelay;

      hDel->nSurroundAnalysisBuffer = 0;
      hDel->nArbDmxAnalysisBuffer = 0;
      tempDelay = nEncoderSynDelay + hDel->config.nLimiterDelay +
                  hDel->config.nCoreCoderDelay +
                  hDel->config.nSacTimeAlignment + nDecoderAnDelay;
      tempDelay = (nResidualCoderFrameDelay * hDel->config.nFrameLen) +
                  hDel->config.nSacStreamMuxDelay - tempDelay;

      if (tempDelay > 0) {
        hDel->nBitstreamFrameBuffer = 0;
        hDel->nOutputAudioBuffer = tempDelay;
      } else {
        tempDelay = -tempDelay;
        hDel->nBitstreamFrameBuffer =
            (tempDelay + hDel->config.nFrameLen - 1) / hDel->config.nFrameLen;
        hDel->nOutputAudioBuffer =
            (hDel->nBitstreamFrameBuffer * hDel->config.nFrameLen) - tempDelay;
      }

      hDel->nOutputAudioQmfFrameBuffer =
          (hDel->nOutputAudioBuffer + (hDel->config.nQmfLen / 2) - 1) /
          hDel->config.nQmfLen;

      if (hDel->config.bDmxAlign > 0) {
        tempDelay = nEncoderWinDelay + nEncoderAnDelay + nEncoderSynDelay +
                    hDel->nOutputAudioBuffer + hDel->config.nLimiterDelay +
                    hDel->config.nCoreCoderDelay;
        hDel->nDiscardOutFrames =
            (tempDelay + hDel->config.nFrameLen - 1) / hDel->config.nFrameLen;
        hDel->nDmxAlignBuffer =
            hDel->nDiscardOutFrames * hDel->config.nFrameLen - tempDelay;
      } else {
        hDel->nDiscardOutFrames = 0;
        hDel->nDmxAlignBuffer = 0;
      }

      hDel->nInfoDmxDelay = hDel->nSurroundAnalysisBuffer + nEncoderAnDelay +
                            nEncoderWinDelay + nEncoderSynDelay +
                            hDel->nOutputAudioBuffer +
                            hDel->config.nLimiterDelay;
      hDel->nInfoCodecDelay =
          hDel->nInfoDmxDelay + hDel->config.nCoreCoderDelay +
          hDel->config.nSacTimeAlignment + nDecoderAnDelay + nDecoderSynDelay;

    } else {
      int tempDelay1, tempDelay2, tempDelay12, tempDelay3;

      tempDelay1 = hDel->config.nArbDmxDelay - hDel->config.nSurroundDelay;

      if (tempDelay1 >= 0) {
        hDel->nSurroundAnalysisBuffer = tempDelay1;
        hDel->nArbDmxAnalysisBuffer = 0;
      } else {
        hDel->nSurroundAnalysisBuffer = 0;
        hDel->nArbDmxAnalysisBuffer = -tempDelay1;
      }

      tempDelay1 = nEncoderWinDelay + hDel->config.nSurroundDelay +
                   hDel->nSurroundAnalysisBuffer +
                   nEncoderAnDelay;
      tempDelay2 = nEncoderWinDelay + hDel->config.nArbDmxDelay +
                   hDel->nArbDmxAnalysisBuffer +
                   nEncoderAnDelay;
      tempDelay3 = hDel->config.nArbDmxDelay + hDel->config.nLimiterDelay +
                   hDel->config.nCoreCoderDelay +
                   hDel->config.nSacTimeAlignment +
                   nDecoderAnDelay;

      tempDelay12 =
          FDKmax(nResidualCoderFrameDelay, nArbDmxResidualCoderFrameDelay) *
          hDel->config.nFrameLen;
      tempDelay12 += hDel->config.nSacStreamMuxDelay;

      if (tempDelay1 > tempDelay2) {
        tempDelay12 += tempDelay1;
      } else {
        tempDelay12 += tempDelay2;
      }

      if (tempDelay3 > tempDelay12) {
        if (hDel->config.bMinimizeDelay > 0) {
          hDel->nBitstreamFrameBuffer =
              (tempDelay3 - tempDelay12) / hDel->config.nFrameLen;
          hDel->nOutputAudioBuffer = 0;
          hDel->nSurroundAnalysisBuffer +=
              (tempDelay3 - tempDelay12 -
               (hDel->nBitstreamFrameBuffer * hDel->config.nFrameLen));
          hDel->nArbDmxAnalysisBuffer +=
              (tempDelay3 - tempDelay12 -
               (hDel->nBitstreamFrameBuffer * hDel->config.nFrameLen));
        } else {
          hDel->nBitstreamFrameBuffer =
              ((tempDelay3 - tempDelay12) + hDel->config.nFrameLen - 1) /
              hDel->config.nFrameLen;
          hDel->nOutputAudioBuffer =
              hDel->nBitstreamFrameBuffer * hDel->config.nFrameLen +
              tempDelay12 - tempDelay3;
        }
      } else {
        hDel->nBitstreamFrameBuffer = 0;
        hDel->nOutputAudioBuffer = tempDelay12 - tempDelay3;
      }

      if (hDel->config.bDmxAlign > 0) {
        int tempDelay = hDel->config.nArbDmxDelay + hDel->nOutputAudioBuffer +
                        hDel->config.nLimiterDelay +
                        hDel->config.nCoreCoderDelay;
        hDel->nDiscardOutFrames =
            (tempDelay + hDel->config.nFrameLen - 1) / hDel->config.nFrameLen;
        hDel->nDmxAlignBuffer =
            hDel->nDiscardOutFrames * hDel->config.nFrameLen - tempDelay;
      } else {
        hDel->nDiscardOutFrames = 0;
        hDel->nDmxAlignBuffer = 0;
      }

      hDel->nInfoDmxDelay = hDel->config.nArbDmxDelay +
                            hDel->nOutputAudioBuffer +
                            hDel->config.nLimiterDelay;
      hDel->nInfoCodecDelay =
          hDel->nInfoDmxDelay + hDel->config.nCoreCoderDelay +
          hDel->config.nSacTimeAlignment + nDecoderAnDelay + nDecoderSynDelay;
      hDel->nInfoDecoderDelay = nDecoderAnDelay + nDecoderSynDelay;
    }

    hDel->nBitstreamFrameBufferSize = hDel->nBitstreamFrameBuffer + 1;
  }

  return error;
}

static FDK_SACENC_ERROR assignParameterInRange(
    const INT startRange,
    const INT stopRange,
    const INT value,
    INT *const ptr
) {
  FDK_SACENC_ERROR error = SACENC_INVALID_CONFIG;

  if ((startRange <= value) && (value <= stopRange)) {
    *ptr = value;
    error = SACENC_OK;
  }

  return error;
}

FDK_SACENC_ERROR fdk_sacenc_delay_SetDmxAlign(HANDLE_DELAY hDelay,
                                              const INT bDmxAlignIn) {
  return (assignParameterInRange(0, 1, bDmxAlignIn, &hDelay->config.bDmxAlign));
}

FDK_SACENC_ERROR fdk_sacenc_delay_SetTimeDomDmx(HANDLE_DELAY hDelay,
                                                const INT bTimeDomDmxIn) {
  return (
      assignParameterInRange(0, 1, bTimeDomDmxIn, &hDelay->config.bTimeDomDmx));
}

FDK_SACENC_ERROR fdk_sacenc_delay_SetSacTimeAlignmentDynamicOut(
    HANDLE_DELAY hDelay, const INT bSacTimeAlignmentDynamicOutIn) {
  return (assignParameterInRange(0, 1, bSacTimeAlignmentDynamicOutIn,
                                 &hDelay->config.bSacTimeAlignmentDynamicOut));
}

FDK_SACENC_ERROR fdk_sacenc_delay_SetNSacTimeAlignment(
    HANDLE_DELAY hDelay, const INT nSacTimeAlignmentIn) {
  return (assignParameterInRange(-32768, 32767, nSacTimeAlignmentIn,
                                 &hDelay->config.nSacTimeAlignment));
}

FDK_SACENC_ERROR fdk_sacenc_delay_SetMinimizeDelay(HANDLE_DELAY hDelay,
                                                   const INT bMinimizeDelay) {
  return (assignParameterInRange(0, 1, bMinimizeDelay,
                                 &hDelay->config.bMinimizeDelay));
}

INT fdk_sacenc_delay_GetOutputAudioBufferDelay(HANDLE_DELAY hDelay) {
  return (hDelay->nOutputAudioBuffer);
}

INT fdk_sacenc_delay_GetSurroundAnalysisBufferDelay(HANDLE_DELAY hDelay) {
  return (hDelay->nSurroundAnalysisBuffer);
}

INT fdk_sacenc_delay_GetArbDmxAnalysisBufferDelay(HANDLE_DELAY hDelay) {
  return (hDelay->nArbDmxAnalysisBuffer);
}

INT fdk_sacenc_delay_GetBitstreamFrameBufferSize(HANDLE_DELAY hDelay) {
  return (hDelay->nBitstreamFrameBufferSize);
}

INT fdk_sacenc_delay_GetDmxAlignBufferDelay(HANDLE_DELAY hDelay) {
  return (hDelay->nDmxAlignBuffer);
}

INT fdk_sacenc_delay_GetDiscardOutFrames(HANDLE_DELAY hDelay) {
  return (hDelay->nDiscardOutFrames);
}

INT fdk_sacenc_delay_GetInfoDmxDelay(HANDLE_DELAY hDelay) {
  return (hDelay->nInfoDmxDelay);
}

INT fdk_sacenc_delay_GetInfoCodecDelay(HANDLE_DELAY hDelay) {
  return (hDelay->nInfoCodecDelay);
}

INT fdk_sacenc_delay_GetInfoDecoderDelay(HANDLE_DELAY hDelay) {
  return (hDelay->nInfoDecoderDelay);
}
