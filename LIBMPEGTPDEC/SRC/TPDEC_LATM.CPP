/******************* MPEG transport format decoder library *********************

   Author(s):   Daniel Homm

   Description:

*******************************************************************************/

#include "tpdec_latm.h"
#include "FDK_bitstream.h"
#define TPDEC_TRACKINDEX(p, l) (1 * (p) + (l))

static UINT CLatmDemux_GetValue(HANDLE_FDK_BITSTREAM bs) {
  UCHAR bytesForValue = 0, tmp = 0;
  int value = 0;

  bytesForValue = (UCHAR)FDKreadBits(bs, 2);

  for (UINT i = 0; i <= bytesForValue; i++) {
    value <<= 8;
    tmp = (UCHAR)FDKreadBits(bs, 8);
    value += tmp;
  }

  return value;
}

static TRANSPORTDEC_ERROR CLatmDemux_ReadAudioMuxElement(
    HANDLE_FDK_BITSTREAM bs, CLatmDemux *pLatmDemux, int m_muxConfigPresent,
    CSTpCallBacks *pTpDecCallbacks, CSAudioSpecificConfig *pAsc,
    int *pfConfigFound) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;

  if (m_muxConfigPresent) {
    pLatmDemux->m_useSameStreamMux = FDKreadBits(bs, 1);

    if (!pLatmDemux->m_useSameStreamMux) {
      int i;
      UCHAR configChanged = 0;
      UCHAR configMode = 0;

      FDK_BITSTREAM bsAnchor;

      FDK_BITSTREAM bsAnchorDummyParse;

      if (!pLatmDemux->applyAsc) {
        bsAnchorDummyParse = *bs;
        pLatmDemux->newCfgHasAudioPreRoll = 0;
        configMode |= AC_CM_DET_CFG_CHANGE;
        if (TRANSPORTDEC_OK !=
            (ErrorStatus = CLatmDemux_ReadStreamMuxConfig(
                 bs, pLatmDemux, pTpDecCallbacks, pAsc, pfConfigFound,
                 configMode, configChanged))) {
          goto bail;
        }

        if (pAsc->m_sc.m_usacConfig.element[0]
                .extElement.usacExtElementHasAudioPreRoll &&
            pLatmDemux->newCfgHasAudioPreRoll) {
          pLatmDemux->newCfgHasAudioPreRoll = 0;
          pLatmDemux->applyAsc = 0;
        } else {
          *bs = bsAnchorDummyParse;
          pLatmDemux->applyAsc = 1;
        }
      }

      if (pLatmDemux->applyAsc) {
        for (i = 0; i < 2; i++) {
          configMode = 0;

          if (i == 0) {
            configMode |= AC_CM_DET_CFG_CHANGE;
            bsAnchor = *bs;
          } else {
            configMode |= AC_CM_ALLOC_MEM;
            *bs = bsAnchor;
          }

          if (TRANSPORTDEC_OK !=
              (ErrorStatus = CLatmDemux_ReadStreamMuxConfig(
                   bs, pLatmDemux, pTpDecCallbacks, pAsc, pfConfigFound,
                   configMode, configChanged))) {
            goto bail;
          }

          if (ErrorStatus == TRANSPORTDEC_OK) {
            if ((i == 0) && (pAsc->AacConfigChanged || pAsc->SbrConfigChanged ||
                             pAsc->SacConfigChanged)) {
              int errC;

              configChanged = 1;
              errC = pTpDecCallbacks->cbFreeMem(pTpDecCallbacks->cbFreeMemData,
                                                pAsc);
              if (errC != 0) {
                ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
                goto bail;
              }
            }
          }
        }
      }
    }
  }

  if (!*pfConfigFound) {
    ErrorStatus = TRANSPORTDEC_SYNC_ERROR;
    goto bail;
  }

  if (pLatmDemux->m_AudioMuxVersionA == 0) {
    if (TRANSPORTDEC_OK !=
        (ErrorStatus = CLatmDemux_ReadPayloadLengthInfo(bs, pLatmDemux))) {
      *pfConfigFound = 0;
      goto bail;
    }
  } else {
    ErrorStatus = TRANSPORTDEC_UNSUPPORTED_FORMAT;
    *pfConfigFound = 0;
    goto bail;
  }

bail:
  if (ErrorStatus != TRANSPORTDEC_OK) {
    pLatmDemux->applyAsc = 1;
  }

  return (ErrorStatus);
}

TRANSPORTDEC_ERROR CLatmDemux_Read(HANDLE_FDK_BITSTREAM bs,
                                   CLatmDemux *pLatmDemux, TRANSPORT_TYPE tt,
                                   CSTpCallBacks *pTpDecCallbacks,
                                   CSAudioSpecificConfig *pAsc,
                                   int *pfConfigFound,
                                   const INT ignoreBufferFullness) {
  UINT cntBits;
  UINT cmpBufferFullness;
  UINT audioMuxLengthBytesLast = 0;
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;

  cntBits = FDKgetValidBits(bs);

  if ((INT)cntBits < MIN_LATM_HEADERLENGTH) {
    return TRANSPORTDEC_NOT_ENOUGH_BITS;
  }

  if (TRANSPORTDEC_OK != (ErrorStatus = CLatmDemux_ReadAudioMuxElement(
                              bs, pLatmDemux, (tt != TT_MP4_LATM_MCP0),
                              pTpDecCallbacks, pAsc, pfConfigFound)))
    return (ErrorStatus);

  if (!ignoreBufferFullness) {
    cmpBufferFullness =
        24 + audioMuxLengthBytesLast * 8 +
        pLatmDemux->m_linfo[0][0].m_bufferFullness *
            pAsc[TPDEC_TRACKINDEX(0, 0)].m_channelConfiguration * 32;

    if (pLatmDemux->m_linfo[0][0].m_bufferFullness != 0xFF) {
      if (!pLatmDemux->BufferFullnessAchieved) {
        if (cntBits < cmpBufferFullness) {
          return TRANSPORTDEC_NOT_ENOUGH_BITS;
        } else {
          pLatmDemux->BufferFullnessAchieved = 1;
        }
      }
    }
  }

  return (ErrorStatus);
}

TRANSPORTDEC_ERROR CLatmDemux_ReadStreamMuxConfig(
    HANDLE_FDK_BITSTREAM bs, CLatmDemux *pLatmDemux,
    CSTpCallBacks *pTpDecCallbacks, CSAudioSpecificConfig *pAsc,
    int *pfConfigFound, UCHAR configMode, UCHAR configChanged) {
  CSAudioSpecificConfig ascDummy;
  CSAudioSpecificConfig *pAscDummy;
  pAscDummy = &ascDummy;
  pLatmDemux->usacExplicitCfgChanged = 0;
  LATM_LAYER_INFO *p_linfo = NULL;
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;
  UCHAR updateConfig[1 * 1] = {0};

  pLatmDemux->m_AudioMuxVersion = FDKreadBits(bs, 1);

  if (pLatmDemux->m_AudioMuxVersion == 0) {
    pLatmDemux->m_AudioMuxVersionA = 0;
  } else {
    pLatmDemux->m_AudioMuxVersionA = FDKreadBits(bs, 1);
  }

  if (pLatmDemux->m_AudioMuxVersionA == 0) {
    if (pLatmDemux->m_AudioMuxVersion == 1) {
      pLatmDemux->m_taraBufferFullness = CLatmDemux_GetValue(bs);
    }
    pLatmDemux->m_allStreamsSameTimeFraming = FDKreadBits(bs, 1);
    pLatmDemux->m_noSubFrames = FDKreadBits(bs, 6) + 1;
    pLatmDemux->m_numProgram = FDKreadBits(bs, 4) + 1;

    if (pLatmDemux->m_numProgram > LATM_MAX_PROG) {
      ErrorStatus = TRANSPORTDEC_UNSUPPORTED_FORMAT;
      goto bail;
    }

    int idCnt = 0;
    for (UINT prog = 0; prog < pLatmDemux->m_numProgram; prog++) {
      pLatmDemux->m_numLayer[prog] = FDKreadBits(bs, 3) + 1;
      if (pLatmDemux->m_numLayer[prog] > LATM_MAX_LAYER) {
        ErrorStatus = TRANSPORTDEC_UNSUPPORTED_FORMAT;
        goto bail;
      }

      for (UINT lay = 0; lay < pLatmDemux->m_numLayer[prog]; lay++) {
        int useSameConfig;
        p_linfo = &pLatmDemux->m_linfo[prog][lay];

        p_linfo->m_streamID = idCnt++;
        p_linfo->m_frameLengthInBits = 0;

        if ((prog == 0) && (lay == 0)) {
          useSameConfig = 0;
        } else {
          useSameConfig = FDKreadBits(bs, 1);
        }

        if (useSameConfig) {
          if (lay > 0) {
            FDKmemcpy(&pAsc[TPDEC_TRACKINDEX(prog, lay)],
                      &pAsc[TPDEC_TRACKINDEX(prog, lay - 1)],
                      sizeof(CSAudioSpecificConfig));
          } else {
            ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
            goto bail;
          }
        } else {
          UINT usacConfigLengthPrev = 0;
          UCHAR usacConfigPrev[TP_USAC_MAX_CONFIG_LEN];

          if (!(pLatmDemux->applyAsc) &&
              (pAsc[TPDEC_TRACKINDEX(prog, lay)].m_aot == AOT_USAC)) {
            usacConfigLengthPrev =
                (UINT)(pAsc[TPDEC_TRACKINDEX(prog, lay)]
                           .m_sc.m_usacConfig.UsacConfigBits +
                       7) >>
                3;
            if (usacConfigLengthPrev > TP_USAC_MAX_CONFIG_LEN) {
              ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
              goto bail;
            }
            FDKmemclear(usacConfigPrev, TP_USAC_MAX_CONFIG_LEN);
            FDKmemcpy(
                usacConfigPrev,
                &pAsc[TPDEC_TRACKINDEX(prog, lay)].m_sc.m_usacConfig.UsacConfig,
                usacConfigLengthPrev);
          }
          if (pLatmDemux->m_AudioMuxVersion == 1) {
            FDK_BITSTREAM tmpBs;
            INT ascLen = 0;
            ascLen = CLatmDemux_GetValue(bs);
            if (ascLen < 0 || ascLen > (INT)FDKgetValidBits(bs)) {
              ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
              goto bail;
            }
            FDKsyncCache(bs);
            tmpBs = *bs;
            tmpBs.hBitBuf.ValidBits = ascLen;

            if (pLatmDemux->applyAsc) {
              if (TRANSPORTDEC_OK !=
                  (ErrorStatus = AudioSpecificConfig_Parse(
                       &pAsc[TPDEC_TRACKINDEX(prog, lay)], &tmpBs, 1,
                       pTpDecCallbacks, configMode, configChanged,
                       AOT_NULL_OBJECT)))
                goto bail;
            } else {
              if (TRANSPORTDEC_OK !=
                  (ErrorStatus = AudioSpecificConfig_Parse(
                       pAscDummy, &tmpBs, 1, pTpDecCallbacks, configMode,
                       configChanged, AOT_NULL_OBJECT)))
                goto bail;
            }

            if (0 > (INT)FDKgetValidBits(&tmpBs)) {
              ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
              goto bail;
            }
            FDKpushFor(bs, ascLen);
          } else {
            if (pLatmDemux->applyAsc) {
              if (TRANSPORTDEC_OK != (ErrorStatus = AudioSpecificConfig_Parse(
                                          &pAsc[TPDEC_TRACKINDEX(prog, lay)],
                                          bs, 0, pTpDecCallbacks, configMode,
                                          configChanged, AOT_NULL_OBJECT)))
                goto bail;
            } else {
              if (TRANSPORTDEC_OK !=
                  (ErrorStatus = AudioSpecificConfig_Parse(
                       pAscDummy, bs, 0, pTpDecCallbacks, configMode,
                       configChanged, AOT_NULL_OBJECT)))
                goto bail;
            }
          }
          if (!pLatmDemux->applyAsc) {
            updateConfig[TPDEC_TRACKINDEX(prog, lay)] = 0;
          } else {
            updateConfig[TPDEC_TRACKINDEX(prog, lay)] = 1;
          }

          if (!pLatmDemux->applyAsc) {
            if (pAscDummy[TPDEC_TRACKINDEX(prog, lay)].m_aot ==
                AOT_USAC) {
              const UINT usacConfigLength =
                  (UINT)(pAscDummy->m_sc.m_usacConfig.UsacConfigBits + 7) >> 3;
              if (usacConfigLength > TP_USAC_MAX_CONFIG_LEN) {
                ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
                goto bail;
              }
              if (usacConfigLength != usacConfigLengthPrev) {
                FDKmemclear(&pAsc[TPDEC_TRACKINDEX(prog, lay)]
                                 .m_sc.m_usacConfig.UsacConfig,
                            TP_USAC_MAX_CONFIG_LEN);
                FDKmemcpy(&pAsc[TPDEC_TRACKINDEX(prog, lay)]
                               .m_sc.m_usacConfig.UsacConfig,
                          &pAscDummy->m_sc.m_usacConfig.UsacConfig,
                          usacConfigLength);
                pAsc[TPDEC_TRACKINDEX(prog, lay)]
                    .m_sc.m_usacConfig.UsacConfigBits =
                    pAscDummy->m_sc.m_usacConfig.UsacConfigBits;
                pLatmDemux->usacExplicitCfgChanged = 1;
              } else {
                if (FDKmemcmp(usacConfigPrev,
                              pAscDummy->m_sc.m_usacConfig.UsacConfig,
                              usacConfigLengthPrev)) {
                  FDKmemclear(&pAsc[TPDEC_TRACKINDEX(prog, lay)]
                                   .m_sc.m_usacConfig.UsacConfig,
                              TP_USAC_MAX_CONFIG_LEN);
                  FDKmemcpy(&pAsc[TPDEC_TRACKINDEX(prog, lay)]
                                 .m_sc.m_usacConfig.UsacConfig,
                            &pAscDummy->m_sc.m_usacConfig.UsacConfig,
                            usacConfigLength);
                  pAsc[TPDEC_TRACKINDEX(prog, lay)]
                      .m_sc.m_usacConfig.UsacConfigBits =
                      pAscDummy->m_sc.m_usacConfig.UsacConfigBits;
                  pLatmDemux->usacExplicitCfgChanged = 1;
                }
              }

              if (pAscDummy[TPDEC_TRACKINDEX(prog, lay)]
                      .m_sc.m_usacConfig.m_usacNumElements) {
                if (pAscDummy[TPDEC_TRACKINDEX(prog, lay)]
                        .m_sc.m_usacConfig.element[0]
                        .extElement.usacExtElementHasAudioPreRoll) {
                  pLatmDemux->newCfgHasAudioPreRoll =
                      1;
                }
              }
            }
          }
        }

        p_linfo->m_frameLengthType = FDKreadBits(bs, 3);
        switch (p_linfo->m_frameLengthType) {
          case 0:
            p_linfo->m_bufferFullness = FDKreadBits(bs, 8);

            if (!pLatmDemux->m_allStreamsSameTimeFraming) {
              if ((lay > 0) &&
                  (pAsc[TPDEC_TRACKINDEX(prog, lay)].m_aot == AOT_AAC_SCAL ||
                   pAsc[TPDEC_TRACKINDEX(prog, lay)].m_aot ==
                       AOT_ER_AAC_SCAL) &&
                  (pAsc[TPDEC_TRACKINDEX(prog, lay - 1)].m_aot == AOT_CELP ||
                   pAsc[TPDEC_TRACKINDEX(prog, lay - 1)].m_aot ==
                       AOT_ER_CELP)) {
                FDKreadBits(bs, 6);
              }
            }
            break;
          case 1:
            p_linfo->m_frameLengthInBits = FDKreadBits(bs, 9);
            break;
          case 3:
          case 4:
          case 5:
          case 6:
          case 7:
          default:
            ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
            goto bail;
        }
      }
    }

    pLatmDemux->m_otherDataPresent = FDKreadBits(bs, 1);
    pLatmDemux->m_otherDataLength = 0;

    if (pLatmDemux->m_otherDataPresent) {
      if (pLatmDemux->m_AudioMuxVersion == 1) {
        pLatmDemux->m_otherDataLength = CLatmDemux_GetValue(bs);
      } else {
        int otherDataLenEsc = 0;
        do {
          pLatmDemux->m_otherDataLength <<= 8;
          otherDataLenEsc = FDKreadBits(bs, 1);
          pLatmDemux->m_otherDataLength += FDKreadBits(bs, 8);
        } while (otherDataLenEsc);
      }
      if (pLatmDemux->m_audioMuxLengthBytes <
          (pLatmDemux->m_otherDataLength >> 3)) {
        ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
        goto bail;
      }
    }

    pLatmDemux->m_crcCheckPresent = FDKreadBits(bs, 1);

    if (pLatmDemux->m_crcCheckPresent) {
      FDKreadBits(bs, 8);
    }

  } else {
    ErrorStatus = TRANSPORTDEC_UNSUPPORTED_FORMAT;
  }

  if (ErrorStatus == TRANSPORTDEC_OK) {
    UINT prog;
    for (prog = 0; prog < pLatmDemux->m_numProgram; prog++) {
      UINT lay;
      for (lay = 0; lay < pLatmDemux->m_numLayer[prog]; lay++) {
        if (updateConfig[TPDEC_TRACKINDEX(prog, lay)] != 0) {
          int cbError;
          cbError = pTpDecCallbacks->cbUpdateConfig(
              pTpDecCallbacks->cbUpdateConfigData,
              &pAsc[TPDEC_TRACKINDEX(prog, lay)],
              pAsc[TPDEC_TRACKINDEX(prog, lay)].configMode,
              &pAsc[TPDEC_TRACKINDEX(prog, lay)].AacConfigChanged);
          if (cbError == TRANSPORTDEC_NEED_TO_RESTART) {
            *pfConfigFound = 0;
            ErrorStatus = TRANSPORTDEC_NEED_TO_RESTART;
            goto bail;
          }
          if (cbError != 0) {
            *pfConfigFound = 0;
            if (lay == 0) {
              ErrorStatus = TRANSPORTDEC_SYNC_ERROR;
              goto bail;
            }
          } else {
            *pfConfigFound = 1;
          }
        } else {
          *pfConfigFound = 1;
        }
      }
    }
  }

bail:
  if (ErrorStatus != TRANSPORTDEC_OK) {
    UCHAR applyAsc = pLatmDemux->applyAsc;
    FDKmemclear(pLatmDemux, sizeof(CLatmDemux));
    pLatmDemux->applyAsc = applyAsc;
  } else {
    if (configMode == AC_CM_ALLOC_MEM) pLatmDemux->applyAsc = 0;
  }

  return (ErrorStatus);
}

static int CLatmDemux_ReadAuChunkLengthInfo(HANDLE_FDK_BITSTREAM bs) {
  int len = 0, tmp = 255;
  int validBytes = (int)FDKgetValidBits(bs) >> 3;

  while (tmp == 255 && validBytes-- > 0) {
    tmp = (int)FDKreadBits(bs, 8);
    len += tmp;
  }

  return ((tmp == 255) ? -1 : (len << 3));
}

TRANSPORTDEC_ERROR CLatmDemux_ReadPayloadLengthInfo(HANDLE_FDK_BITSTREAM bs,
                                                    CLatmDemux *pLatmDemux) {
  TRANSPORTDEC_ERROR ErrorStatus = TRANSPORTDEC_OK;
  int totalPayloadBits = 0;

  if (pLatmDemux->m_allStreamsSameTimeFraming == 1) {
    FDK_ASSERT(pLatmDemux->m_numProgram <= LATM_MAX_PROG);
    for (UINT prog = 0; prog < pLatmDemux->m_numProgram; prog++) {
      FDK_ASSERT(pLatmDemux->m_numLayer[prog] <= LATM_MAX_LAYER);
      for (UINT lay = 0; lay < pLatmDemux->m_numLayer[prog]; lay++) {
        LATM_LAYER_INFO *p_linfo = &pLatmDemux->m_linfo[prog][lay];
        int auChunkLengthInfo = 0;

        switch (p_linfo->m_frameLengthType) {
          case 0:
            auChunkLengthInfo = CLatmDemux_ReadAuChunkLengthInfo(bs);
            if (auChunkLengthInfo >= 0) {
              p_linfo->m_frameLengthInBits = (UINT)auChunkLengthInfo;
              totalPayloadBits += p_linfo->m_frameLengthInBits;
            } else {
              return TRANSPORTDEC_PARSE_ERROR;
            }
            break;
          case 3:
          case 5:
          case 7:
          default:
            return TRANSPORTDEC_PARSE_ERROR;
        }
      }
    }
  } else {
    ErrorStatus = TRANSPORTDEC_PARSE_ERROR;
  }
  if (pLatmDemux->m_audioMuxLengthBytes > (UINT)0 &&
      totalPayloadBits > (int)pLatmDemux->m_audioMuxLengthBytes * 8) {
    return TRANSPORTDEC_PARSE_ERROR;
  }

  return (ErrorStatus);
}

UINT CLatmDemux_GetFrameLengthInBits(CLatmDemux *pLatmDemux, const UINT prog,
                                     const UINT layer) {
  UINT nFrameLenBits = 0;
  if (prog < pLatmDemux->m_numProgram) {
    if (layer < pLatmDemux->m_numLayer[prog]) {
      nFrameLenBits = pLatmDemux->m_linfo[prog][layer].m_frameLengthInBits;
    }
  }
  return nFrameLenBits;
}

UINT CLatmDemux_GetOtherDataPresentFlag(CLatmDemux *pLatmDemux) {
  return pLatmDemux->m_otherDataPresent ? 1 : 0;
}

UINT CLatmDemux_GetOtherDataLength(CLatmDemux *pLatmDemux) {
  return pLatmDemux->m_otherDataLength;
}

UINT CLatmDemux_GetNrOfSubFrames(CLatmDemux *pLatmDemux) {
  return pLatmDemux->m_noSubFrames;
}

UINT CLatmDemux_GetNrOfLayers(CLatmDemux *pLatmDemux, const UINT prog) {
  UINT numLayer = 0;
  if (prog < pLatmDemux->m_numProgram) {
    numLayer = pLatmDemux->m_numLayer[prog];
  }
  return numLayer;
}
