/******************* MPEG transport format decoder library *********************

   Author(s):   Daniel Homm

   Description:

*******************************************************************************/

#ifndef TPDEC_LATM_H
#define TPDEC_LATM_H

#include "tpdec_lib.h"

#include "FDK_bitstream.h"

#define MIN_LATM_HEADERLENGTH 9
#define MIN_LOAS_HEADERLENGTH MIN_LATM_HEADERLENGTH + 24
#define MIN_TP_BUF_SIZE_LOAS (8194)

enum {
  LATM_MAX_PROG = 1,
  LATM_MAX_LAYER = 1,
  LATM_MAX_VAR_CHUNKS = 16,
  LATM_MAX_ID = 16
};

typedef struct {
  UINT m_frameLengthType;
  UINT m_bufferFullness;
  UINT m_streamID;
  UINT m_frameLengthInBits;
} LATM_LAYER_INFO;

typedef struct {
  LATM_LAYER_INFO m_linfo[LATM_MAX_PROG][LATM_MAX_LAYER];
  UINT m_taraBufferFullness;
  UINT m_otherDataLength;
  UINT m_audioMuxLengthBytes;

  UCHAR m_useSameStreamMux;
  UCHAR m_AudioMuxVersion;
  UCHAR m_AudioMuxVersionA;
  UCHAR m_allStreamsSameTimeFraming;
  UCHAR m_noSubFrames;
  UCHAR m_numProgram;
  UCHAR m_numLayer[LATM_MAX_PROG];

  UCHAR m_otherDataPresent;
  UCHAR m_crcCheckPresent;

  SCHAR BufferFullnessAchieved;
  UCHAR
  usacExplicitCfgChanged;
  UCHAR applyAsc;
  UCHAR newCfgHasAudioPreRoll;
} CLatmDemux;

TRANSPORTDEC_ERROR CLatmDemux_Read(HANDLE_FDK_BITSTREAM bs,
                                   CLatmDemux *pLatmDemux, TRANSPORT_TYPE tt,
                                   CSTpCallBacks *pTpDecCallbacks,
                                   CSAudioSpecificConfig *pAsc,
                                   int *pfConfigFound,
                                   const INT ignoreBufferFullness);

TRANSPORTDEC_ERROR CLatmDemux_ReadStreamMuxConfig(
    HANDLE_FDK_BITSTREAM bs, CLatmDemux *pLatmDemux,
    CSTpCallBacks *pTpDecCallbacks, CSAudioSpecificConfig *pAsc,
    int *pfConfigFound, UCHAR configMode, UCHAR configChanged);

TRANSPORTDEC_ERROR CLatmDemux_ReadPayloadLengthInfo(HANDLE_FDK_BITSTREAM bs,
                                                    CLatmDemux *pLatmDemux);

UINT CLatmDemux_GetFrameLengthInBits(CLatmDemux *pLatmDemux, const UINT prog,
                                     const UINT layer);
UINT CLatmDemux_GetOtherDataPresentFlag(CLatmDemux *pLatmDemux);
UINT CLatmDemux_GetOtherDataLength(CLatmDemux *pLatmDemux);
UINT CLatmDemux_GetNrOfSubFrames(CLatmDemux *pLatmDemux);
UINT CLatmDemux_GetNrOfLayers(CLatmDemux *pLatmDemux, const UINT program);

#endif /* TPDEC_LATM_H */
