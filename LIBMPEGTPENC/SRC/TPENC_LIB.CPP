/******************* MPEG transport format encoder library *********************

   Author(s):   Manuel Jander

   Description: MPEG Transport encode

*******************************************************************************/

#include "tpenc_lib.h"
#include "tp_version.h"

#define MODULE_NAME "transportEnc"

#include "tpenc_asc.h"
#include "tpenc_adts.h"
#include "tpenc_adif.h"
#include "tpenc_latm.h"

typedef struct {
  int curSubFrame;
  int nSubFrames;
  int prevBits;
} RAWPACKETS_INFO;

struct TRANSPORTENC {
  CODER_CONFIG config;
  TRANSPORT_TYPE transportFmt;

  FDK_BITSTREAM bitStream;
  UCHAR *bsBuffer;
  INT bsBufferSize;

  INT pceFrameCounter;
  union {
    STRUCT_ADTS adts;

    ADIF_INFO adif;

    LATM_STREAM latm;

    RAWPACKETS_INFO raw;

  } writer;

  CSTpCallBacks callbacks;
};

typedef struct _TRANSPORTENC_STRUCT TRANSPORTENC_STRUCT;

C_ALLOC_MEM(Ram_TransportEncoder, struct TRANSPORTENC, 1)

TRANSPORTENC_ERROR transportEnc_Open(HANDLE_TRANSPORTENC *phTpEnc) {
  HANDLE_TRANSPORTENC hTpEnc;

  if (phTpEnc == NULL) {
    return TRANSPORTENC_INVALID_PARAMETER;
  }

  hTpEnc = GetRam_TransportEncoder(0);

  if (hTpEnc == NULL) {
    return TRANSPORTENC_NO_MEM;
  }

  *phTpEnc = hTpEnc;
  return TRANSPORTENC_OK;
}

static INT getPceRepetitionRate(const CHANNEL_MODE channelMode,
                                const int channelConfigZero,
                                const TRANSPORT_TYPE transportFmt,
                                const int headerPeriod,
                                const int matrixMixdownA) {
  INT pceFrameCounter = -1;

  if (headerPeriod > 0) {
    switch (getChannelConfig(channelMode, channelConfigZero)) {
      case 0:
        switch (transportFmt) {
          case TT_MP4_ADTS:
          case TT_MP4_LATM_MCP0:
          case TT_MP4_RAW:
            pceFrameCounter = headerPeriod;
            break;
          case TT_MP4_ADIF:
            if ((channelMode == MODE_1_2_2) || (channelMode == MODE_1_2_2_1)) {
              pceFrameCounter = headerPeriod;
              break;
            }
            FDK_FALLTHROUGH;
          case TT_MP4_LOAS:
          case TT_MP4_LATM_MCP1:
          default:
            pceFrameCounter = -1;
        }
        break;
      case 5:
      case 6:
        if (matrixMixdownA != 0) {
          switch (transportFmt) {
            case TT_MP4_ADIF:
            case TT_MP4_ADTS:
            case TT_MP4_LOAS:
            case TT_MP4_LATM_MCP1:
            case TT_MP4_LATM_MCP0:
            case TT_MP4_RAW:
              pceFrameCounter = headerPeriod;
              break;
            default:
              pceFrameCounter = -1;
          }
        }
        break;
      default:
        pceFrameCounter = -1;
    }
  }
  else {
    pceFrameCounter = -1;
  }

  return pceFrameCounter;
}

TRANSPORTENC_ERROR transportEnc_Init(HANDLE_TRANSPORTENC hTpEnc,
                                     UCHAR *bsBuffer, INT bsBufferSize,
                                     TRANSPORT_TYPE transportFmt,
                                     CODER_CONFIG *cconfig, UINT flags) {
  FDKmemcpy(&hTpEnc->config, cconfig, sizeof(CODER_CONFIG));

  hTpEnc->transportFmt = transportFmt;

  hTpEnc->bsBuffer = bsBuffer;
  hTpEnc->bsBufferSize = bsBufferSize;

  FDKinitBitStream(&hTpEnc->bitStream, hTpEnc->bsBuffer, hTpEnc->bsBufferSize,
                   0, BS_WRITER);

  switch (transportFmt) {
    case TT_MP4_ADIF:
      if ((hTpEnc->config.aot != AOT_AAC_LC) ||
          (hTpEnc->config.samplesPerFrame != 1024)) {
        return TRANSPORTENC_INVALID_PARAMETER;
      }
      hTpEnc->writer.adif.headerWritten = 0;
      hTpEnc->writer.adif.samplingRate = hTpEnc->config.samplingRate;
      hTpEnc->writer.adif.bitRate = hTpEnc->config.bitRate;
      hTpEnc->writer.adif.profile = ((int)hTpEnc->config.aot) - 1;
      hTpEnc->writer.adif.cm = hTpEnc->config.channelMode;
      hTpEnc->writer.adif.bVariableRate = 0;
      hTpEnc->writer.adif.instanceTag = 0;
      hTpEnc->writer.adif.matrixMixdownA = hTpEnc->config.matrixMixdownA;
      hTpEnc->writer.adif.pseudoSurroundEnable =
          (hTpEnc->config.flags & CC_PSEUDO_SURROUND) ? 1 : 0;
      break;

    case TT_MP4_ADTS:
      if ((hTpEnc->config.aot != AOT_AAC_LC) ||
          (hTpEnc->config.samplesPerFrame != 1024)) {
        return TRANSPORTENC_INVALID_PARAMETER;
      }
      if (adtsWrite_Init(&hTpEnc->writer.adts, &hTpEnc->config) != 0) {
        return TRANSPORTENC_INVALID_PARAMETER;
      }
      break;

    case TT_MP4_LOAS:
    case TT_MP4_LATM_MCP0:
    case TT_MP4_LATM_MCP1: {
      TRANSPORTENC_ERROR error;

      error = transportEnc_Latm_Init(&hTpEnc->writer.latm, &hTpEnc->bitStream,
                                     &hTpEnc->config, flags & TP_FLAG_LATM_AMV,
                                     transportFmt, &hTpEnc->callbacks);
      if (error != TRANSPORTENC_OK) {
        return error;
      }
    } break;

    case TT_MP4_RAW:
      hTpEnc->writer.raw.curSubFrame = 0;
      hTpEnc->writer.raw.nSubFrames = hTpEnc->config.nSubFrames;
      break;

    default:
      return TRANSPORTENC_INVALID_PARAMETER;
  }

  hTpEnc->pceFrameCounter = getPceRepetitionRate(
      hTpEnc->config.channelMode, hTpEnc->config.channelConfigZero,
      transportFmt, hTpEnc->config.headerPeriod, hTpEnc->config.matrixMixdownA);

  return TRANSPORTENC_OK;
}

TRANSPORTENC_ERROR transportEnc_AddOtherDataBits(HANDLE_TRANSPORTENC hTpEnc,
                                                 const int nBits) {
  TRANSPORTENC_ERROR tpErr = TRANSPORTENC_OK;

  switch (hTpEnc->transportFmt) {
    case TT_MP4_LATM_MCP0:
    case TT_MP4_LATM_MCP1:
    case TT_MP4_LOAS:
      tpErr = transportEnc_LatmAddOtherDataBits(&hTpEnc->writer.latm, nBits);
      break;
    case TT_MP4_ADTS:
    case TT_MP4_ADIF:
    case TT_MP4_RAW:
    default:
      tpErr = TRANSPORTENC_UNKOWN_ERROR;
  }

  return tpErr;
}

HANDLE_FDK_BITSTREAM transportEnc_GetBitstream(HANDLE_TRANSPORTENC hTp) {
  return &hTp->bitStream;
}

int transportEnc_RegisterSbrCallback(HANDLE_TRANSPORTENC hTpEnc,
                                     const cbSbr_t cbSbr, void *user_data) {
  if (hTpEnc == NULL) {
    return -1;
  }
  hTpEnc->callbacks.cbSbr = cbSbr;
  hTpEnc->callbacks.cbSbrData = user_data;
  return 0;
}
int transportEnc_RegisterUsacCallback(HANDLE_TRANSPORTENC hTpEnc,
                                      const cbUsac_t cbUsac, void *user_data) {
  if (hTpEnc == NULL) {
    return -1;
  }
  hTpEnc->callbacks.cbUsac = cbUsac;
  hTpEnc->callbacks.cbUsacData = user_data;
  return 0;
}

int transportEnc_RegisterSscCallback(HANDLE_TRANSPORTENC hTpEnc,
                                     const cbSsc_t cbSsc, void *user_data) {
  if (hTpEnc == NULL) {
    return -1;
  }
  hTpEnc->callbacks.cbSsc = cbSsc;
  hTpEnc->callbacks.cbSscData = user_data;
  return 0;
}

TRANSPORTENC_ERROR transportEnc_WriteAccessUnit(HANDLE_TRANSPORTENC hTp,
                                                INT frameUsedBits,
                                                int bufferFullness, int ncc) {
  TRANSPORTENC_ERROR err = TRANSPORTENC_OK;

  if (!hTp) {
    return TRANSPORTENC_INVALID_PARAMETER;
  }
  HANDLE_FDK_BITSTREAM hBs = &hTp->bitStream;

  if (hTp->pceFrameCounter >= hTp->config.headerPeriod) {
    frameUsedBits += transportEnc_GetPCEBits(
        hTp->config.channelMode, hTp->config.matrixMixdownA,
        3);
  }

  switch (hTp->transportFmt) {
    case TT_MP4_ADIF:
      FDKinitBitStream(&hTp->bitStream, hTp->bsBuffer, hTp->bsBufferSize, 0,
                       BS_WRITER);
      if (0 != adifWrite_EncodeHeader(&hTp->writer.adif, hBs, bufferFullness)) {
        err = TRANSPORTENC_INVALID_CONFIG;
      }
      break;
    case TT_MP4_ADTS:
      bufferFullness /= ncc;
      bufferFullness /= 32;
      bufferFullness = FDKmin(0x7FF, bufferFullness);
      adtsWrite_EncodeHeader(&hTp->writer.adts, &hTp->bitStream, bufferFullness,
                             frameUsedBits);
      break;
    case TT_MP4_LOAS:
    case TT_MP4_LATM_MCP0:
    case TT_MP4_LATM_MCP1:
      bufferFullness /= ncc;
      bufferFullness /= 32;
      bufferFullness = FDKmin(0xFF, bufferFullness);
      transportEnc_LatmWrite(&hTp->writer.latm, hBs, frameUsedBits,
                             bufferFullness, &hTp->callbacks);
      break;
    case TT_MP4_RAW:
      if (hTp->writer.raw.curSubFrame >= hTp->writer.raw.nSubFrames) {
        hTp->writer.raw.curSubFrame = 0;
        FDKinitBitStream(&hTp->bitStream, hTp->bsBuffer, hTp->bsBufferSize, 0,
                         BS_WRITER);
      }
      hTp->writer.raw.prevBits = FDKgetValidBits(hBs);
      break;
    default:
      err = TRANSPORTENC_UNSUPPORTED_FORMAT;
      break;
  }

  if (hTp->pceFrameCounter >= hTp->config.headerPeriod) {
    INT crcIndex = 0;
    UINT alignAnchor = FDKgetValidBits(&hTp->bitStream);

    FDKwriteBits(&hTp->bitStream, ID_PCE, 3);

    if ((hTp->transportFmt == TT_MP4_ADTS) &&
        !hTp->writer.adts.protection_absent) {
      crcIndex = adtsWrite_CrcStartReg(&hTp->writer.adts, &hTp->bitStream, 0);
    }

    transportEnc_writePCE(
        &hTp->bitStream, hTp->config.channelMode, hTp->config.samplingRate, 0,
        1, hTp->config.matrixMixdownA,
        (hTp->config.flags & CC_PSEUDO_SURROUND) ? 1 : 0, alignAnchor);

    if ((hTp->transportFmt == TT_MP4_ADTS) &&
        !hTp->writer.adts.protection_absent) {
      adtsWrite_CrcEndReg(&hTp->writer.adts, &hTp->bitStream, crcIndex);
    }
    hTp->pceFrameCounter = 0;
  }

  if (hTp->pceFrameCounter != -1) {
    hTp->pceFrameCounter++;
  }

  return err;
}

TRANSPORTENC_ERROR transportEnc_EndAccessUnit(HANDLE_TRANSPORTENC hTp,
                                              int *bits) {
  switch (hTp->transportFmt) {
    case TT_MP4_LATM_MCP0:
    case TT_MP4_LATM_MCP1:
    case TT_MP4_LOAS:
      transportEnc_LatmAdjustSubframeBits(&hTp->writer.latm, bits);
      break;
    case TT_MP4_ADTS:
      adtsWrite_EndRawDataBlock(&hTp->writer.adts, &hTp->bitStream, bits);
      break;
    case TT_MP4_ADIF:
      *bits -= adifWrite_GetHeaderBits(&hTp->writer.adif);
      hTp->writer.adif.headerWritten = 1;
      break;
    case TT_MP4_RAW:
      *bits -= hTp->writer.raw.prevBits;
      break;
    default:
      break;
  }

  return TRANSPORTENC_OK;
}

TRANSPORTENC_ERROR transportEnc_GetFrame(HANDLE_TRANSPORTENC hTpEnc,
                                         int *nbytes) {
  TRANSPORTENC_ERROR tpErr = TRANSPORTENC_OK;
  HANDLE_FDK_BITSTREAM hBs = &hTpEnc->bitStream;

  switch (hTpEnc->transportFmt) {
    case TT_MP4_LATM_MCP0:
    case TT_MP4_LATM_MCP1:
    case TT_MP4_LOAS:
      *nbytes = hTpEnc->bsBufferSize;
      tpErr = transportEnc_LatmGetFrame(&hTpEnc->writer.latm, hBs, nbytes);
      break;
    case TT_MP4_ADTS:
      if (hTpEnc->writer.adts.currentBlock >=
          hTpEnc->writer.adts.num_raw_blocks + 1) {
        *nbytes = (FDKgetValidBits(hBs) + 7) >> 3;
        hTpEnc->writer.adts.currentBlock = 0;
      } else {
        *nbytes = 0;
      }
      break;
    case TT_MP4_ADIF:
      FDK_ASSERT((INT)FDKgetValidBits(hBs) >= 0);
      *nbytes = (FDKgetValidBits(hBs) + 7) >> 3;
      break;
    case TT_MP4_RAW:
      FDKsyncCache(hBs);
      hTpEnc->writer.raw.curSubFrame++;
      *nbytes = ((FDKgetValidBits(hBs) - hTpEnc->writer.raw.prevBits) + 7) >> 3;
      break;
    default:
      break;
  }

  return tpErr;
}

INT transportEnc_GetStaticBits(HANDLE_TRANSPORTENC hTp, int auBits) {
  INT nbits = 0, nPceBits = 0;

  if (hTp->pceFrameCounter >= hTp->config.headerPeriod) {
    nPceBits = transportEnc_GetPCEBits(
        hTp->config.channelMode, hTp->config.matrixMixdownA,
        3);
    auBits += nPceBits;
  }

  switch (hTp->transportFmt) {
    case TT_MP4_ADIF:
    case TT_MP4_RAW:
      nbits = 0;
      break;
    case TT_MP4_ADTS:
      nbits = adtsWrite_GetHeaderBits(&hTp->writer.adts);
      break;
    case TT_MP4_LOAS:
    case TT_MP4_LATM_MCP0:
    case TT_MP4_LATM_MCP1:
      nbits =
          transportEnc_LatmCountTotalBitDemandHeader(&hTp->writer.latm, auBits);
      break;
    default:
      nbits = 0;
      break;
  }

  nbits += nPceBits;

  return nbits;
}

void transportEnc_Close(HANDLE_TRANSPORTENC *phTp) {
  if (phTp != NULL) {
    if (*phTp != NULL) {
      FreeRam_TransportEncoder(phTp);
    }
  }
}

int transportEnc_CrcStartReg(HANDLE_TRANSPORTENC hTpEnc, int mBits) {
  int crcReg = 0;

  switch (hTpEnc->transportFmt) {
    case TT_MP4_ADTS:
      crcReg = adtsWrite_CrcStartReg(&hTpEnc->writer.adts, &hTpEnc->bitStream,
                                     mBits);
      break;
    default:
      break;
  }

  return crcReg;
}

void transportEnc_CrcEndReg(HANDLE_TRANSPORTENC hTpEnc, int reg) {
  switch (hTpEnc->transportFmt) {
    case TT_MP4_ADTS:
      adtsWrite_CrcEndReg(&hTpEnc->writer.adts, &hTpEnc->bitStream, reg);
      break;
    default:
      break;
  }
}

TRANSPORTENC_ERROR transportEnc_GetConf(HANDLE_TRANSPORTENC hTpEnc,
                                        CODER_CONFIG *cc,
                                        FDK_BITSTREAM *dataBuffer,
                                        UINT *confType) {
  TRANSPORTENC_ERROR tpErr = TRANSPORTENC_OK;
  HANDLE_LATM_STREAM hLatmConfig = &hTpEnc->writer.latm;

  *confType = 0;

  switch (hTpEnc->transportFmt) {
    case TT_MP4_LATM_MCP0:
    case TT_MP4_LATM_MCP1:
    case TT_MP4_LOAS:
      tpErr =
          CreateStreamMuxConfig(hLatmConfig, dataBuffer, 0, &hTpEnc->callbacks);
      *confType = 1;
      break;
    default:
      if (transportEnc_writeASC(dataBuffer, cc, &hTpEnc->callbacks) != 0) {
        tpErr = TRANSPORTENC_UNKOWN_ERROR;
      }
  }

  return tpErr;
}

TRANSPORTENC_ERROR transportEnc_GetLibInfo(LIB_INFO *info) {
  int i;

  if (info == NULL) {
    return TRANSPORTENC_INVALID_PARAMETER;
  }
  for (i = 0; i < FDK_MODULE_LAST; i++) {
    if (info[i].module_id == FDK_NONE) break;
  }
  if (i == FDK_MODULE_LAST) {
    return TRANSPORTENC_UNKOWN_ERROR;
  }
  info += i;

  info->module_id = FDK_TPENC;
  info->version = LIB_VERSION(TP_LIB_VL0, TP_LIB_VL1, TP_LIB_VL2);
  LIB_VERSION_STRING(info);
  info->build_date = "";
  info->build_time = "";
  info->title = TP_LIB_TITLE;

  info->flags =
      0 | CAPF_ADIF | CAPF_ADTS | CAPF_LATM | CAPF_LOAS | CAPF_RAWPACKETS;

  return TRANSPORTENC_OK;
}
