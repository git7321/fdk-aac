/************************* MPEG-D DRC decoder library **************************

   Author(s):

   Description:

*******************************************************************************/

#include "fixpoint_math.h"
#include "drcDec_reader.h"
#include "drcDec_tools.h"
#include "drcDec_rom.h"
#include "drcDecoder.h"

#define UNIDRCCONFEXT_PARAM_DRC 0x1
#define UNIDRCCONFEXT_V1 0x2
#define UNIDRCLOUDEXT_EQ 0x1

#define UNIDRCGAINEXT_TERM 0x0
#define UNIDRCLOUDEXT_TERM 0x0
#define UNIDRCCONFEXT_TERM 0x0

static int _getZ(const int nNodesMax) {
  int Z = 1;
  while ((1 << Z) < (2 * nNodesMax)) {
    Z++;
  }
  return Z;
}

static int _getTimeDeltaMin(const GAIN_SET* pGset, const int deltaTminDefault) {
  if (pGset->timeDeltaMinPresent) {
    return pGset->timeDeltaMin;
  } else {
    return deltaTminDefault;
  }
}

static inline int _compAssign(UCHAR* dest, const UCHAR src) {
  int diff = 0;
  if (*dest != src) diff = 1;
  *dest = src;
  return diff;
}

static inline int _compAssign(ULONG* dest, const ULONG src) {
  int diff = 0;
  if (*dest != src) diff = 1;
  *dest = src;
  return diff;
}

typedef const SCHAR (*Huffman)[2];

int _decodeHuffmanCW(Huffman h,
                     HANDLE_FDK_BITSTREAM hBs)
{
  SCHAR index = 0;
  int value, bit;

  while (index >= 0) {
    bit = FDKreadBits(hBs, 1);
    index = h[index][bit];
  }

  value = index + 64;

  return value;
}

DRC_ERROR
drcDec_readUniDrc(HANDLE_FDK_BITSTREAM hBs, HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
                  HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet,
                  const int frameSize, const int deltaTminDefault,
                  HANDLE_UNI_DRC_GAIN hUniDrcGain) {
  DRC_ERROR err = DE_OK;
  int loudnessInfoSetPresent, uniDrcConfigPresent;

  loudnessInfoSetPresent = FDKreadBits(hBs, 1);
  if (loudnessInfoSetPresent) {
    uniDrcConfigPresent = FDKreadBits(hBs, 1);
    if (uniDrcConfigPresent) {
      err = drcDec_readUniDrcConfig(hBs, hUniDrcConfig);
      if (err) {
        FDKmemclear(hUniDrcConfig, sizeof(UNI_DRC_CONFIG));
        hUniDrcConfig->diff = 1;
      }
    }
    err = drcDec_readLoudnessInfoSet(hBs, hLoudnessInfoSet);
    if (err) {
      FDKmemclear(hLoudnessInfoSet, sizeof(LOUDNESS_INFO_SET));
      hLoudnessInfoSet->diff = 1;
    }
  }

  err = drcDec_readUniDrcGain(hBs, hUniDrcConfig, frameSize, deltaTminDefault,
                              hUniDrcGain);

  return err;
}

static FIXP_SGL _decodeGainInitial(
    HANDLE_FDK_BITSTREAM hBs, const GAIN_CODING_PROFILE gainCodingProfile) {
  int sign, magn;
  FIXP_SGL gainInitial = (FIXP_SGL)0;
  switch (gainCodingProfile) {
    case GCP_REGULAR:
      sign = FDKreadBits(hBs, 1);
      magn = FDKreadBits(hBs, 8);

      gainInitial =
          (FIXP_SGL)(magn << (FRACT_BITS - 1 - 3 - 7));
      if (sign) gainInitial = -gainInitial;
      break;
    case GCP_FADING:
      sign = FDKreadBits(hBs, 1);
      if (sign == 0)
        gainInitial = (FIXP_SGL)0;
      else {
        magn = FDKreadBits(hBs, 10);
        gainInitial = -(FIXP_SGL)(
            (magn + 1) << (FRACT_BITS - 1 - 3 - 7));
      }
      break;
    case GCP_CLIPPING_DUCKING:
      sign = FDKreadBits(hBs, 1);
      if (sign == 0)
        gainInitial = (FIXP_SGL)0;
      else {
        magn = FDKreadBits(hBs, 8);
        gainInitial = -(FIXP_SGL)(
            (magn + 1) << (FRACT_BITS - 1 - 3 - 7));
      }
      break;
    case GCP_CONSTANT:
      break;
  }
  return gainInitial;
}

static int _decodeNNodes(HANDLE_FDK_BITSTREAM hBs) {
  int nNodes = 0, endMarker = 0;

  while (endMarker != 1) {
    nNodes++;
    if (nNodes >= 128) break;
    endMarker = FDKreadBits(hBs, 1);
  }
  return nNodes;
}

static void _decodeGains(HANDLE_FDK_BITSTREAM hBs,
                         const GAIN_CODING_PROFILE gainCodingProfile,
                         const int nNodes, GAIN_NODE* pNodes) {
  int k, deltaGain;
  Huffman deltaGainCodebook;

  pNodes[0].gainDb = _decodeGainInitial(hBs, gainCodingProfile);

  if (gainCodingProfile == GCP_CLIPPING_DUCKING) {
    deltaGainCodebook = (Huffman)&deltaGain_codingProfile_2_huffman;
  } else {
    deltaGainCodebook = (Huffman)&deltaGain_codingProfile_0_1_huffman;
  }

  for (k = 1; k < nNodes; k++) {
    deltaGain = _decodeHuffmanCW(deltaGainCodebook, hBs);
    if (k >= 16) continue;
    pNodes[k].gainDb =
        pNodes[k - 1].gainDb +
        (FIXP_SGL)(deltaGain << (FRACT_BITS - 1 - 7 -
                                 3));
  }
}

static void _decodeSlopes(HANDLE_FDK_BITSTREAM hBs,
                          const GAIN_INTERPOLATION_TYPE gainInterpolationType,
                          const int nNodes, GAIN_NODE* pNodes) {
  int k = 0;

  if (gainInterpolationType == GIT_SPLINE) {
    for (k = 0; k < nNodes; k++) {
      _decodeHuffmanCW((Huffman)&slopeSteepness_huffman, hBs);
    }
  }
}

static int _decodeTimeDelta(HANDLE_FDK_BITSTREAM hBs, const int Z) {
  int prefix, mu;

  prefix = FDKreadBits(hBs, 2);
  switch (prefix) {
    case 0x0:
      return 1;
    case 0x1:
      mu = FDKreadBits(hBs, 2);
      return mu + 2;
    case 0x2:
      mu = FDKreadBits(hBs, 3);
      return mu + 6;
    case 0x3:
      mu = FDKreadBits(hBs, Z);
      return mu + 14;
    default:
      return 0;
  }
}

static void _decodeTimes(HANDLE_FDK_BITSTREAM hBs, const int deltaTmin,
                         const int frameSize, const int fullFrame,
                         const int timeOffset, const int Z, const int nNodes,
                         GAIN_NODE* pNodes) {
  int timeDelta, k;
  int timeOffs = timeOffset;
  int frameEndFlag, nodeTimeTmp, nodeResFlag;

  if (fullFrame == 0) {
    frameEndFlag = FDKreadBits(hBs, 1);
  } else {
    frameEndFlag = 1;
  }

  if (frameEndFlag ==
      1) {
    nodeResFlag = 0;
    for (k = 0; k < nNodes - 1; k++) {
      timeDelta = _decodeTimeDelta(hBs, Z);
      if (k >= (16 - 1)) continue;
      nodeTimeTmp = timeOffs + timeDelta * deltaTmin;
      if (nodeTimeTmp > frameSize + timeOffset) {
        if (nodeResFlag == 0) {
          pNodes[k].time = frameSize + timeOffset;
          nodeResFlag = 1;
        }
        pNodes[k + 1].time = nodeTimeTmp;
      } else {
        pNodes[k].time = nodeTimeTmp;
      }
      timeOffs = nodeTimeTmp;
    }
    if (nodeResFlag == 0) {
      k = fMin(k, 16 - 1);
      pNodes[k].time = frameSize + timeOffset;
    }
  } else {
    for (k = 0; k < nNodes; k++) {
      timeDelta = _decodeTimeDelta(hBs, Z);
      if (k >= 16) continue;
      pNodes[k].time = timeOffs + timeDelta * deltaTmin;
      timeOffs = pNodes[k].time;
    }
  }
}

static void _readNodes(HANDLE_FDK_BITSTREAM hBs, GAIN_SET* gainSet,
                       const int frameSize, const int timeDeltaMin,
                       UCHAR* pNNodes, GAIN_NODE* pNodes) {
  int timeOffset, drcGainCodingMode, nNodes;
  int Z = _getZ(frameSize / timeDeltaMin);
  if (gainSet->timeAlignment == 0) {
    timeOffset = -1;
  } else {
    timeOffset = -timeDeltaMin +
                 (timeDeltaMin - 1) /
                     2;
  }

  drcGainCodingMode = FDKreadBits(hBs, 1);
  if (drcGainCodingMode == 0) {
    nNodes = 1;
    pNodes[0].gainDb = _decodeGainInitial(
        hBs, (GAIN_CODING_PROFILE)gainSet->gainCodingProfile);
    pNodes[0].time = frameSize + timeOffset;
  } else {
    nNodes = _decodeNNodes(hBs);

    _decodeSlopes(hBs, (GAIN_INTERPOLATION_TYPE)gainSet->gainInterpolationType,
                  nNodes, pNodes);

    _decodeTimes(hBs, timeDeltaMin, frameSize, gainSet->fullFrame, timeOffset,
                 Z, nNodes, pNodes);

    _decodeGains(hBs, (GAIN_CODING_PROFILE)gainSet->gainCodingProfile, nNodes,
                 pNodes);
  }
  *pNNodes = (UCHAR)nNodes;
}

static void _readDrcGainSequence(HANDLE_FDK_BITSTREAM hBs, GAIN_SET* gainSet,
                                 const int frameSize, const int timeDeltaMin,
                                 UCHAR* pNNodes, GAIN_NODE pNodes[16]) {
  SHORT timeBufPrevFrame[16], timeBufCurFrame[16];
  int nNodesNodeRes, nNodesCur, k, m;

  if (gainSet->gainCodingProfile == GCP_CONSTANT) {
    *pNNodes = 1;
    pNodes[0].time = frameSize - 1;
    pNodes[0].gainDb = (FIXP_SGL)0;
  } else {
    _readNodes(hBs, gainSet, frameSize, timeDeltaMin, pNNodes, pNodes);

    nNodesNodeRes = 0;
    nNodesCur = 0;
    for (k = 0; k < *pNNodes; k++) {
      if (k >= 16) continue;
      if (pNodes[k].time >= frameSize) {
        timeBufPrevFrame[nNodesNodeRes] = pNodes[k].time;
        nNodesNodeRes++;
      } else {
        timeBufCurFrame[nNodesCur] = pNodes[k].time;
        nNodesCur++;
      }
    }
    for (k = 0; k < nNodesNodeRes; k++) {
      pNodes[k].time = timeBufPrevFrame[k] - 2 * frameSize;
    }
    for (m = 0; m < nNodesCur; m++, k++) {
      pNodes[k].time = timeBufCurFrame[m];
    }
  }
}

static DRC_ERROR _readUniDrcGainExtension(HANDLE_FDK_BITSTREAM hBs,
                                          UNI_DRC_GAIN_EXTENSION* pExt) {
  DRC_ERROR err = DE_OK;
  int k, bitSizeLen, extSizeBits, bitSize;

  k = 0;
  pExt->uniDrcGainExtType[k] = FDKreadBits(hBs, 4);
  while (pExt->uniDrcGainExtType[k] != UNIDRCGAINEXT_TERM) {
    if (k >= (8 - 1)) return DE_MEMORY_ERROR;
    bitSizeLen = FDKreadBits(hBs, 3);
    extSizeBits = bitSizeLen + 4;

    bitSize = FDKreadBits(hBs, extSizeBits);
    pExt->extBitSize[k] = bitSize + 1;

    switch (pExt->uniDrcGainExtType[k]) {
      default:
        FDKpushFor(hBs, pExt->extBitSize[k]);
        break;
    }
    k++;
    pExt->uniDrcGainExtType[k] = FDKreadBits(hBs, 4);
  }

  return err;
}

DRC_ERROR
drcDec_readUniDrcGain(HANDLE_FDK_BITSTREAM hBs,
                      HANDLE_UNI_DRC_CONFIG hUniDrcConfig, const int frameSize,
                      const int deltaTminDefault,
                      HANDLE_UNI_DRC_GAIN hUniDrcGain) {
  DRC_ERROR err = DE_OK;
  int seq, gainSequenceCount;
  DRC_COEFFICIENTS_UNI_DRC* pCoef =
      selectDrcCoefficients(hUniDrcConfig, LOCATION_SELECTED);
  if (hUniDrcGain == NULL) return DE_NOT_OK;
  hUniDrcGain->status = 0;
  if (pCoef) {
    gainSequenceCount = fMin(pCoef->gainSequenceCount, (UCHAR)12);
  } else {
    gainSequenceCount = 0;
  }

  for (seq = 0; seq < gainSequenceCount; seq++) {
    UCHAR index = pCoef->gainSetIndexForGainSequence[seq];
    GAIN_SET* gainSet;
    int timeDeltaMin;
    UCHAR tmpNNodes = 0;
    GAIN_NODE tmpNodes[16];

    if ((index >= pCoef->gainSetCount) || (index >= 12)) return DE_NOT_OK;
    gainSet = &(pCoef->gainSet[index]);

    timeDeltaMin = _getTimeDeltaMin(gainSet, deltaTminDefault);

    _readDrcGainSequence(hBs, gainSet, frameSize, timeDeltaMin, &tmpNNodes,
                         tmpNodes);

    hUniDrcGain->nNodes[seq] = tmpNNodes;
    FDKmemcpy(hUniDrcGain->gainNode[seq], tmpNodes,
              fMin(tmpNNodes, (UCHAR)16) * sizeof(GAIN_NODE));
  }

  if (pCoef && (gainSequenceCount ==
                pCoef->gainSequenceCount)) {
    hUniDrcGain->uniDrcGainExtPresent = FDKreadBits(hBs, 1);
    if (hUniDrcGain->uniDrcGainExtPresent == 1) {
      err = _readUniDrcGainExtension(hBs, &(hUniDrcGain->uniDrcGainExtension));
      if (err) return err;
    }
  }

  if (err == DE_OK && gainSequenceCount > 0) {
    hUniDrcGain->status = 1;
  }
  return err;
}

static void _decodeDuckingModification(HANDLE_FDK_BITSTREAM hBs,
                                       DUCKING_MODIFICATION* pDMod, int isBox) {
  int bsDuckingScaling, sigma, mu;

  if (isBox) FDKpushFor(hBs, 7);
  pDMod->duckingScalingPresent = FDKreadBits(hBs, 1);

  if (pDMod->duckingScalingPresent) {
    if (isBox) FDKpushFor(hBs, 4);
    bsDuckingScaling = FDKreadBits(hBs, 4);
    sigma = bsDuckingScaling >> 3;
    mu = bsDuckingScaling & 0x7;

    if (sigma) {
      pDMod->duckingScaling = (FIXP_SGL)(
          (7 - mu) << (FRACT_BITS - 1 - 3 - 2));
    } else {
      pDMod->duckingScaling = (FIXP_SGL)(
          (9 + mu) << (FRACT_BITS - 1 - 3 - 2));
    }
  } else {
    pDMod->duckingScaling = (FIXP_SGL)(1 << (FRACT_BITS - 1 - 2));
  }
}

static void _decodeGainModification(HANDLE_FDK_BITSTREAM hBs, const int version,
                                    int bandCount, GAIN_MODIFICATION* pGMod,
                                    int isBox) {
  int sign, bsGainOffset, bsAttenuationScaling, bsAmplificationScaling;

  if (version > 0) {
    int b, shapeFilterPresent;

    if (isBox) {
      FDKpushFor(hBs, 4);
      bandCount = FDKreadBits(hBs, 4);
    }

    for (b = 0; b < bandCount; b++) {
      if (isBox) {
        FDKpushFor(hBs, 4);
        pGMod[b].targetCharacteristicLeftPresent = FDKreadBits(hBs, 1);
        pGMod[b].targetCharacteristicRightPresent = FDKreadBits(hBs, 1);
        pGMod[b].gainScalingPresent = FDKreadBits(hBs, 1);
        pGMod[b].gainOffsetPresent = FDKreadBits(hBs, 1);
      }

      if (!isBox)
        pGMod[b].targetCharacteristicLeftPresent = FDKreadBits(hBs, 1);
      if (pGMod[b].targetCharacteristicLeftPresent) {
        if (isBox) FDKpushFor(hBs, 4);
        pGMod[b].targetCharacteristicLeftIndex = FDKreadBits(hBs, 4);
      }
      if (!isBox)
        pGMod[b].targetCharacteristicRightPresent = FDKreadBits(hBs, 1);
      if (pGMod[b].targetCharacteristicRightPresent) {
        if (isBox) FDKpushFor(hBs, 4);
        pGMod[b].targetCharacteristicRightIndex = FDKreadBits(hBs, 4);
      }
      if (!isBox) pGMod[b].gainScalingPresent = FDKreadBits(hBs, 1);
      if (pGMod[b].gainScalingPresent) {
        bsAttenuationScaling = FDKreadBits(hBs, 4);
        pGMod[b].attenuationScaling = (FIXP_SGL)(
            bsAttenuationScaling
            << (FRACT_BITS - 1 - 3 - 2));
        bsAmplificationScaling = FDKreadBits(hBs, 4);
        pGMod[b].amplificationScaling = (FIXP_SGL)(
            bsAmplificationScaling
            << (FRACT_BITS - 1 - 3 - 2));
      }
      if (!isBox) pGMod[b].gainOffsetPresent = FDKreadBits(hBs, 1);
      if (pGMod[b].gainOffsetPresent) {
        if (isBox) FDKpushFor(hBs, 2);
        sign = FDKreadBits(hBs, 1);
        bsGainOffset = FDKreadBits(hBs, 5);
        pGMod[b].gainOffset = (FIXP_SGL)(
            (1 + bsGainOffset)
            << (FRACT_BITS - 1 - 2 - 4));
        if (sign) {
          pGMod[b].gainOffset = -pGMod[b].gainOffset;
        }
      }
    }
    if (bandCount == 1) {
      shapeFilterPresent = FDKreadBits(hBs, 1);
      if (shapeFilterPresent) {
        if (isBox) FDKpushFor(hBs, 3);
        FDKpushFor(hBs, 4);
      } else {
        if (isBox) FDKpushFor(hBs, 7);
      }
    }
  } else {
    int b, gainScalingPresent, gainOffsetPresent;
    FIXP_SGL attenuationScaling = FL2FXCONST_SGL(1.0f / (float)(1 << 2)),
             amplificationScaling = FL2FXCONST_SGL(1.0f / (float)(1 << 2)),
             gainOffset = (FIXP_SGL)0;
    if (isBox) FDKpushFor(hBs, 7);
    gainScalingPresent = FDKreadBits(hBs, 1);
    if (gainScalingPresent) {
      bsAttenuationScaling = FDKreadBits(hBs, 4);
      attenuationScaling = (FIXP_SGL)(
          bsAttenuationScaling
          << (FRACT_BITS - 1 - 3 - 2));
      bsAmplificationScaling = FDKreadBits(hBs, 4);
      amplificationScaling = (FIXP_SGL)(
          bsAmplificationScaling
          << (FRACT_BITS - 1 - 3 - 2));
    }
    if (isBox) FDKpushFor(hBs, 7);
    gainOffsetPresent = FDKreadBits(hBs, 1);
    if (gainOffsetPresent) {
      if (isBox) FDKpushFor(hBs, 2);
      sign = FDKreadBits(hBs, 1);
      bsGainOffset = FDKreadBits(hBs, 5);
      gainOffset =
          (FIXP_SGL)((1 + bsGainOffset) << (FRACT_BITS - 1 - 2 -
                                            4));
      if (sign) {
        gainOffset = -gainOffset;
      }
    }
    for (b = 0; b < 4; b++) {
      pGMod[b].targetCharacteristicLeftPresent = 0;
      pGMod[b].targetCharacteristicRightPresent = 0;
      pGMod[b].gainScalingPresent = gainScalingPresent;
      pGMod[b].attenuationScaling = attenuationScaling;
      pGMod[b].amplificationScaling = amplificationScaling;
      pGMod[b].gainOffsetPresent = gainOffsetPresent;
      pGMod[b].gainOffset = gainOffset;
    }
  }
}

static void _readDrcCharacteristic(HANDLE_FDK_BITSTREAM hBs, const int version,
                                   DRC_CHARACTERISTIC* pDChar, int isBox) {
  if (version == 0) {
    if (isBox) FDKpushFor(hBs, 1);
    pDChar->cicpIndex = FDKreadBits(hBs, 7);
    if (pDChar->cicpIndex > 0) {
      pDChar->present = 1;
      pDChar->isCICP = 1;
    } else {
      pDChar->present = 0;
    }
  } else {
    pDChar->present = FDKreadBits(hBs, 1);
    if (isBox) pDChar->isCICP = FDKreadBits(hBs, 1);
    if (pDChar->present) {
      if (!isBox) pDChar->isCICP = FDKreadBits(hBs, 1);
      if (pDChar->isCICP) {
        if (isBox) FDKpushFor(hBs, 1);
        pDChar->cicpIndex = FDKreadBits(hBs, 7);
      } else {
        pDChar->custom.left = FDKreadBits(hBs, 4);
        pDChar->custom.right = FDKreadBits(hBs, 4);
      }
    }
  }
}

static void _readBandBorder(HANDLE_FDK_BITSTREAM hBs, BAND_BORDER* pBBord,
                            int drcBandType, int isBox) {
  if (drcBandType) {
    if (isBox) FDKpushFor(hBs, 4);
    pBBord->crossoverFreqIndex = FDKreadBits(hBs, 4);
  } else {
    if (isBox) FDKpushFor(hBs, 6);
    pBBord->startSubBandIndex = FDKreadBits(hBs, 10);
  }
}

static DRC_ERROR _readGainSet(HANDLE_FDK_BITSTREAM hBs, const int version,
                              int* gainSequenceIndex, GAIN_SET* pGSet,
                              int isBox) {
  if (isBox) FDKpushFor(hBs, 2);
  pGSet->gainCodingProfile = FDKreadBits(hBs, 2);
  pGSet->gainInterpolationType = FDKreadBits(hBs, 1);
  pGSet->fullFrame = FDKreadBits(hBs, 1);
  pGSet->timeAlignment = FDKreadBits(hBs, 1);
  pGSet->timeDeltaMinPresent = FDKreadBits(hBs, 1);

  if (pGSet->timeDeltaMinPresent) {
    int bsTimeDeltaMin;
    if (isBox) FDKpushFor(hBs, 5);
    bsTimeDeltaMin = FDKreadBits(hBs, 11);
    pGSet->timeDeltaMin = bsTimeDeltaMin + 1;
  }

  if (pGSet->gainCodingProfile != GCP_CONSTANT) {
    int i;
    if (isBox) FDKpushFor(hBs, 3);
    pGSet->bandCount = FDKreadBits(hBs, 4);
    if (pGSet->bandCount > 4) return DE_MEMORY_ERROR;

    if ((pGSet->bandCount > 1) || isBox) {
      pGSet->drcBandType = FDKreadBits(hBs, 1);
    }

    for (i = 0; i < pGSet->bandCount; i++) {
      if (version == 0) {
        *gainSequenceIndex = (*gainSequenceIndex) + 1;
      } else {
        int indexPresent;
        indexPresent = (isBox) ? 1 : FDKreadBits(hBs, 1);
        if (indexPresent) {
          int bsIndex;
          bsIndex = FDKreadBits(hBs, 6);
          *gainSequenceIndex = bsIndex;
        } else {
          *gainSequenceIndex = (*gainSequenceIndex) + 1;
        }
      }
      pGSet->gainSequenceIndex[i] = *gainSequenceIndex;
      _readDrcCharacteristic(hBs, version, &(pGSet->drcCharacteristic[i]),
                             isBox);
    }
    for (i = 1; i < pGSet->bandCount; i++) {
      _readBandBorder(hBs, &(pGSet->bandBorder[i]), pGSet->drcBandType, isBox);
    }
  } else {
    pGSet->bandCount = 1;
    *gainSequenceIndex = (*gainSequenceIndex) + 1;
    pGSet->gainSequenceIndex[0] = *gainSequenceIndex;
  }

  return DE_OK;
}

static DRC_ERROR _readCustomDrcCharacteristic(HANDLE_FDK_BITSTREAM hBs,
                                              const CHARACTERISTIC_SIDE side,
                                              UCHAR* pCharacteristicFormat,
                                              CUSTOM_DRC_CHAR* pCChar,
                                              int isBox) {
  if (isBox) FDKpushFor(hBs, 7);
  *pCharacteristicFormat = FDKreadBits(hBs, 1);
  if (*pCharacteristicFormat == CF_SIGMOID) {
    int bsGain, bsIoRatio, bsExp;
    if (isBox) FDKpushFor(hBs, 1);
    bsGain = FDKreadBits(hBs, 6);
    if (side == CS_LEFT) {
      pCChar->sigmoid.gain = (FIXP_SGL)(bsGain << (FRACT_BITS - 1 - 6));
    } else {
      pCChar->sigmoid.gain = (FIXP_SGL)(-bsGain << (FRACT_BITS - 1 - 6));
    }
    bsIoRatio = FDKreadBits(hBs, 4);
    pCChar->sigmoid.ioRatio =
        FL2FXCONST_SGL(0.05f / (float)(1 << 2)) +
        (FIXP_SGL)((((3 * bsIoRatio) << (FRACT_BITS - 1)) / 5) >> 4);
    bsExp = FDKreadBits(hBs, 4);
    if (bsExp < 15) {
      pCChar->sigmoid.exp = (FIXP_SGL)((1 + 2 * bsExp) << (FRACT_BITS - 1 - 5));
    } else {
      pCChar->sigmoid.exp = (FIXP_SGL)MAXVAL_SGL;
    }
    pCChar->sigmoid.flipSign = FDKreadBits(hBs, 1);
  } else {
    int i, bsCharacteristicNodeCount, bsNodeLevelDelta, bsNodeGain;
    if (isBox) FDKpushFor(hBs, 6);
    bsCharacteristicNodeCount = FDKreadBits(hBs, 2);
    pCChar->nodes.characteristicNodeCount = bsCharacteristicNodeCount + 1;
    if (pCChar->nodes.characteristicNodeCount > 4) return DE_MEMORY_ERROR;
    pCChar->nodes.nodeLevel[0] = DRC_INPUT_LOUDNESS_TARGET_SGL;
    pCChar->nodes.nodeGain[0] = (FIXP_SGL)0;
    for (i = 0; i < pCChar->nodes.characteristicNodeCount; i++) {
      if (isBox) FDKpushFor(hBs, 3);
      bsNodeLevelDelta = FDKreadBits(hBs, 5);
      if (side == CS_LEFT) {
        pCChar->nodes.nodeLevel[i + 1] =
            pCChar->nodes.nodeLevel[i] -
            (FIXP_SGL)((1 + bsNodeLevelDelta) << (FRACT_BITS - 1 - 7));
      } else {
        pCChar->nodes.nodeLevel[i + 1] =
            pCChar->nodes.nodeLevel[i] +
            (FIXP_SGL)((1 + bsNodeLevelDelta) << (FRACT_BITS - 1 - 7));
      }
      bsNodeGain = FDKreadBits(hBs, 8);
      pCChar->nodes.nodeGain[i + 1] = (FIXP_SGL)(
          (bsNodeGain - 128)
          << (FRACT_BITS - 1 - 1 - 7));
    }
  }
  return DE_OK;
}

static void _skipLoudEqInstructions(HANDLE_FDK_BITSTREAM hBs) {
  int i;
  int downmixIdPresent, additionalDownmixIdPresent,
      additionalDownmixIdCount = 0;
  int drcSetIdPresent, additionalDrcSetIdPresent, additionalDrcSetIdCount = 0;
  int eqSetIdPresent, additionalEqSetIdPresent, additionalEqSetIdCount = 0;
  int loudEqGainSequenceCount, drcCharacteristicFormatIsCICP;

  FDKpushFor(hBs, 4);
  FDKpushFor(hBs, 4);
  downmixIdPresent = FDKreadBits(hBs, 1);
  if (downmixIdPresent) {
    FDKpushFor(hBs, 7);
    additionalDownmixIdPresent = FDKreadBits(hBs, 1);
    if (additionalDownmixIdPresent) {
      additionalDownmixIdCount = FDKreadBits(hBs, 7);
      for (i = 0; i < additionalDownmixIdCount; i++) {
        FDKpushFor(hBs, 7);
      }
    }
  }

  drcSetIdPresent = FDKreadBits(hBs, 1);
  if (drcSetIdPresent) {
    FDKpushFor(hBs, 6);
    additionalDrcSetIdPresent = FDKreadBits(hBs, 1);
    if (additionalDrcSetIdPresent) {
      additionalDrcSetIdCount = FDKreadBits(hBs, 6);
      for (i = 0; i < additionalDrcSetIdCount; i++) {
        FDKpushFor(hBs, 6);
      }
    }
  }

  eqSetIdPresent = FDKreadBits(hBs, 1);
  if (eqSetIdPresent) {
    FDKpushFor(hBs, 6);
    additionalEqSetIdPresent = FDKreadBits(hBs, 1);
    if (additionalEqSetIdPresent) {
      additionalEqSetIdCount = FDKreadBits(hBs, 6);
      for (i = 0; i < additionalEqSetIdCount; i++) {
        FDKpushFor(hBs, 6);
      }
    }
  }

  FDKpushFor(hBs, 1);
  FDKpushFor(hBs, 1);
  loudEqGainSequenceCount = FDKreadBits(hBs, 6);
  for (i = 0; i < loudEqGainSequenceCount; i++) {
    FDKpushFor(hBs, 6);
    drcCharacteristicFormatIsCICP = FDKreadBits(hBs, 1);
    if (drcCharacteristicFormatIsCICP) {
      FDKpushFor(hBs, 7);
    } else {
      FDKpushFor(hBs, 4);
      FDKpushFor(hBs, 4);
    }
    FDKpushFor(hBs, 6);
    FDKpushFor(hBs, 3);
    FDKpushFor(hBs, 5);
  }
}

static void _skipEqSubbandGainSpline(HANDLE_FDK_BITSTREAM hBs) {
  int nEqNodes, k, bits;
  nEqNodes = FDKreadBits(hBs, 5);
  nEqNodes += 2;
  for (k = 0; k < nEqNodes; k++) {
    bits = FDKreadBits(hBs, 1);
    if (!bits) {
      FDKpushFor(hBs, 4);
    }
  }
  FDKpushFor(hBs, 4 * (nEqNodes - 1));
  bits = FDKreadBits(hBs, 2);
  switch (bits) {
    case 0:
      FDKpushFor(hBs, 5);
      break;
    case 1:
    case 2:
      FDKpushFor(hBs, 4);
      break;
    case 3:
      FDKpushFor(hBs, 3);
      break;
  }
  FDKpushFor(hBs, 5 * (nEqNodes - 1));
}

static void _skipEqCoefficients(HANDLE_FDK_BITSTREAM hBs) {
  int j, k;
  int eqDelayMaxPresent;
  int uniqueFilterBlockCount, filterElementCount, filterElementGainPresent;
  int uniqueTdFilterElementCount, eqFilterFormat, bsRealZeroRadiusOneCount,
      realZeroCount, genericZeroCount, realPoleCount, complexPoleCount,
      firFilterOrder;
  int uniqueEqSubbandGainsCount, eqSubbandGainRepresentation,
      eqSubbandGainCount;
  int eqSubbandGainFormat;

  eqDelayMaxPresent = FDKreadBits(hBs, 1);
  if (eqDelayMaxPresent) {
    FDKpushFor(hBs, 8);
  }

  uniqueFilterBlockCount = FDKreadBits(hBs, 6);
  for (j = 0; j < uniqueFilterBlockCount; j++) {
    filterElementCount = FDKreadBits(hBs, 6);
    for (k = 0; k < filterElementCount; k++) {
      FDKpushFor(hBs, 6);
      filterElementGainPresent = FDKreadBits(hBs, 1);
      if (filterElementGainPresent) {
        FDKpushFor(hBs, 10);
      }
    }
  }
  uniqueTdFilterElementCount = FDKreadBits(hBs, 6);
  for (j = 0; j < uniqueTdFilterElementCount; j++) {
    eqFilterFormat = FDKreadBits(hBs, 1);
    if (eqFilterFormat == 0) {
      bsRealZeroRadiusOneCount = FDKreadBits(hBs, 3);
      realZeroCount = FDKreadBits(hBs, 6);
      genericZeroCount = FDKreadBits(hBs, 6);
      realPoleCount = FDKreadBits(hBs, 4);
      complexPoleCount = FDKreadBits(hBs, 4);
      FDKpushFor(hBs, 2 * bsRealZeroRadiusOneCount * 1);
      FDKpushFor(hBs, realZeroCount * 8);
      FDKpushFor(hBs, genericZeroCount * 14);
      FDKpushFor(hBs, realPoleCount * 8);
      FDKpushFor(hBs, complexPoleCount * 14);
    } else {
      firFilterOrder = FDKreadBits(hBs, 7);
      FDKpushFor(hBs, 1);
      FDKpushFor(hBs, (firFilterOrder / 2 + 1) * 11);
    }
  }
  uniqueEqSubbandGainsCount = FDKreadBits(hBs, 6);
  if (uniqueEqSubbandGainsCount > 0) {
    eqSubbandGainRepresentation = FDKreadBits(hBs, 1);
    eqSubbandGainFormat = FDKreadBits(hBs, 4);
    switch (eqSubbandGainFormat) {
      case GF_QMF32:
        eqSubbandGainCount = 32;
        break;
      case GF_QMFHYBRID39:
        eqSubbandGainCount = 39;
        break;
      case GF_QMF64:
        eqSubbandGainCount = 64;
        break;
      case GF_QMFHYBRID71:
        eqSubbandGainCount = 71;
        break;
      case GF_QMF128:
        eqSubbandGainCount = 128;
        break;
      case GF_QMFHYBRID135:
        eqSubbandGainCount = 135;
        break;
      case GF_UNIFORM:
      default:
        eqSubbandGainCount = FDKreadBits(hBs, 8);
        eqSubbandGainCount++;
        break;
    }
    for (k = 0; k < uniqueEqSubbandGainsCount; k++) {
      if (eqSubbandGainRepresentation == 1) {
        _skipEqSubbandGainSpline(hBs);
      } else {
        FDKpushFor(hBs, eqSubbandGainCount * 9);
      }
    }
  }
}

static void _skipTdFilterCascade(HANDLE_FDK_BITSTREAM hBs,
                                 const int eqChannelGroupCount) {
  int i, eqCascadeGainPresent, filterBlockCount, eqPhaseAlignmentPresent;
  for (i = 0; i < eqChannelGroupCount; i++) {
    eqCascadeGainPresent = FDKreadBits(hBs, 1);
    if (eqCascadeGainPresent) {
      FDKpushFor(hBs, 10);
    }
    filterBlockCount = FDKreadBits(hBs, 4);
    FDKpushFor(hBs, filterBlockCount * 7);
  }
  eqPhaseAlignmentPresent = FDKreadBits(hBs, 1);
  {
    if (eqPhaseAlignmentPresent) {
      for (i = 0; i < eqChannelGroupCount; i++) {
        FDKpushFor(hBs, (eqChannelGroupCount - i - 1) * 1);
      }
    }
  }
}

static DRC_ERROR _skipEqInstructions(HANDLE_FDK_BITSTREAM hBs,
                                     HANDLE_UNI_DRC_CONFIG hUniDrcConfig) {
  DRC_ERROR err = DE_OK;
  int c, i, k, channelCount;
  int downmixIdPresent, downmixId, eqApplyToDownmix, additionalDownmixIdPresent,
      additionalDownmixIdCount = 0;
  int additionalDrcSetIdPresent, additionalDrcSetIdCount;
  int dependsOnEqSetPresent, eqChannelGroupCount, tdFilterCascadePresent,
      subbandGainsPresent, eqTransitionDurationPresent;
  UCHAR eqChannelGroupForChannel[8];

  FDKpushFor(hBs, 6);
  FDKpushFor(hBs, 4);
  downmixIdPresent = FDKreadBits(hBs, 1);
  if (downmixIdPresent) {
    downmixId = FDKreadBits(hBs, 7);
    eqApplyToDownmix = FDKreadBits(hBs, 1);
    additionalDownmixIdPresent = FDKreadBits(hBs, 1);
    if (additionalDownmixIdPresent) {
      additionalDownmixIdCount = FDKreadBits(hBs, 7);
      FDKpushFor(hBs, additionalDownmixIdCount * 7);
    }
  } else {
    downmixId = 0;
    eqApplyToDownmix = 0;
  }
  FDKpushFor(hBs, 6);
  additionalDrcSetIdPresent = FDKreadBits(hBs, 1);
  if (additionalDrcSetIdPresent) {
    additionalDrcSetIdCount = FDKreadBits(hBs, 6);
    for (i = 0; i < additionalDrcSetIdCount; i++) {
      FDKpushFor(hBs, 6);
    }
  }
  FDKpushFor(hBs, 16);
  dependsOnEqSetPresent = FDKreadBits(hBs, 1);
  if (dependsOnEqSetPresent) {
    FDKpushFor(hBs, 6);
  } else {
    FDKpushFor(hBs, 1);
  }

  channelCount = hUniDrcConfig->channelLayout.baseChannelCount;
  if ((downmixIdPresent == 1) && (eqApplyToDownmix == 1) && (downmixId != 0) &&
      (downmixId != DOWNMIX_ID_ANY_DOWNMIX) &&
      (additionalDownmixIdCount == 0)) {
    DOWNMIX_INSTRUCTIONS* pDown =
        selectDownmixInstructions(hUniDrcConfig, downmixId);
    if (pDown == NULL) return DE_NOT_OK;

    channelCount =
        pDown->targetChannelCount;
  } else if ((downmixId == DOWNMIX_ID_ANY_DOWNMIX) ||
             (additionalDownmixIdCount > 1)) {
    channelCount = 1;
  }

  eqChannelGroupCount = 0;
  for (c = 0; c < channelCount; c++) {
    int newGroup = 1;
    if (c >= 8) return DE_MEMORY_ERROR;
    eqChannelGroupForChannel[c] = FDKreadBits(hBs, 7);
    for (k = 0; k < c; k++) {
      if (eqChannelGroupForChannel[c] == eqChannelGroupForChannel[k]) {
        newGroup = 0;
      }
    }
    if (newGroup == 1) {
      eqChannelGroupCount += 1;
    }
  }
  tdFilterCascadePresent = FDKreadBits(hBs, 1);
  if (tdFilterCascadePresent) {
    _skipTdFilterCascade(hBs, eqChannelGroupCount);
  }
  subbandGainsPresent = FDKreadBits(hBs, 1);
  if (subbandGainsPresent) {
    FDKpushFor(hBs, eqChannelGroupCount * 6);
  }
  eqTransitionDurationPresent = FDKreadBits(hBs, 1);
  if (eqTransitionDurationPresent) {
    FDKpushFor(hBs, 5);
  }
  return err;
}

static void _skipDrcCoefficientsBasic(HANDLE_FDK_BITSTREAM hBs) {
  FDKpushFor(hBs, 4);
  FDKpushFor(hBs, 7);
}

static DRC_ERROR _readDrcCoefficientsUniDrc(HANDLE_FDK_BITSTREAM hBs,
                                            const int version,
                                            DRC_COEFFICIENTS_UNI_DRC* pCoef) {
  DRC_ERROR err = DE_OK;
  int i, bsDrcFrameSize;
  int gainSequenceIndex = -1;

  pCoef->drcLocation = FDKreadBits(hBs, 4);
  pCoef->drcFrameSizePresent = FDKreadBits(hBs, 1);

  if (pCoef->drcFrameSizePresent == 1) {
    bsDrcFrameSize = FDKreadBits(hBs, 15);
    pCoef->drcFrameSize = bsDrcFrameSize + 1;
  }
  if (version == 0) {
    int gainSequenceCount = 0, gainSetCount;
    pCoef->characteristicLeftCount = 0;
    pCoef->characteristicRightCount = 0;
    gainSetCount = FDKreadBits(hBs, 6);
    pCoef->gainSetCount = fMin(gainSetCount, 12);
    for (i = 0; i < gainSetCount; i++) {
      GAIN_SET tmpGset;
      FDKmemclear(&tmpGset, sizeof(GAIN_SET));
      err = _readGainSet(hBs, version, &gainSequenceIndex, &tmpGset, 0);
      if (err) return err;
      gainSequenceCount += tmpGset.bandCount;

      if (i >= 12) continue;
      pCoef->gainSet[i] = tmpGset;
    }
    pCoef->gainSequenceCount = gainSequenceCount;
  } else {
    UCHAR drcCharacteristicLeftPresent, drcCharacteristicRightPresent;
    UCHAR shapeFiltersPresent, shapeFilterCount, tmpPresent;
    int gainSetCount;
    drcCharacteristicLeftPresent = FDKreadBits(hBs, 1);
    if (drcCharacteristicLeftPresent) {
      pCoef->characteristicLeftCount = FDKreadBits(hBs, 4);
      if ((pCoef->characteristicLeftCount + 1) > 16) return DE_MEMORY_ERROR;
      for (i = 0; i < pCoef->characteristicLeftCount; i++) {
        err = _readCustomDrcCharacteristic(
            hBs, CS_LEFT, &(pCoef->characteristicLeftFormat[i + 1]),
            &(pCoef->customCharacteristicLeft[i + 1]), 0);
        if (err) return err;
      }
    }
    drcCharacteristicRightPresent = FDKreadBits(hBs, 1);
    if (drcCharacteristicRightPresent) {
      pCoef->characteristicRightCount = FDKreadBits(hBs, 4);
      if ((pCoef->characteristicRightCount + 1) > 16) return DE_MEMORY_ERROR;
      for (i = 0; i < pCoef->characteristicRightCount; i++) {
        err = _readCustomDrcCharacteristic(
            hBs, CS_RIGHT, &(pCoef->characteristicRightFormat[i + 1]),
            &(pCoef->customCharacteristicRight[i + 1]), 0);
        if (err) return err;
      }
    }
    shapeFiltersPresent = FDKreadBits(hBs, 1);
    if (shapeFiltersPresent) {
      shapeFilterCount = FDKreadBits(hBs, 4);
      for (i = 0; i < shapeFilterCount; i++) {
        tmpPresent = FDKreadBits(hBs, 1);
        if (tmpPresent)
          FDKpushFor(hBs, 5);

        tmpPresent = FDKreadBits(hBs, 1);
        if (tmpPresent)
          FDKpushFor(hBs, 5);

        tmpPresent = FDKreadBits(hBs, 1);
        if (tmpPresent)
          FDKpushFor(hBs, 5);

        tmpPresent = FDKreadBits(hBs, 1);
        if (tmpPresent)
          FDKpushFor(hBs, 5);
      }
    }
    pCoef->gainSequenceCount = FDKreadBits(hBs, 6);
    gainSetCount = FDKreadBits(hBs, 6);
    pCoef->gainSetCount = fMin(gainSetCount, 12);
    for (i = 0; i < gainSetCount; i++) {
      GAIN_SET tmpGset;
      FDKmemclear(&tmpGset, sizeof(GAIN_SET));
      err = _readGainSet(hBs, version, &gainSequenceIndex, &tmpGset, 0);
      if (err) return err;

      if (i >= 12) continue;
      pCoef->gainSet[i] = tmpGset;
    }
  }
  for (i = 0; i < 12; i++) {
    pCoef->gainSetIndexForGainSequence[i] = 255;
  }
  for (i = 0; i < pCoef->gainSetCount; i++) {
    int b;
    for (b = 0; b < pCoef->gainSet[i].bandCount; b++) {
      if (pCoef->gainSet[i].gainSequenceIndex[b] >= 12) continue;
      pCoef->gainSetIndexForGainSequence[pCoef->gainSet[i]
                                             .gainSequenceIndex[b]] = i;
    }
  }

  return err;
}

static void _skipDrcInstructionsBasic(HANDLE_FDK_BITSTREAM hBs) {
  int drcSetEffect;
  int additionalDownmixIdPresent, additionalDownmixIdCount,
      limiterPeakTargetPresent;
  int drcSetTargetLoudnessPresent, drcSetTargetLoudnessValueLowerPresent;

  FDKpushFor(hBs, 6);
  FDKpushFor(hBs, 4);
  FDKpushFor(hBs, 7);
  additionalDownmixIdPresent = FDKreadBits(hBs, 1);
  if (additionalDownmixIdPresent) {
    additionalDownmixIdCount = FDKreadBits(hBs, 3);
    FDKpushFor(hBs, 7 * additionalDownmixIdCount);
  }

  drcSetEffect = FDKreadBits(hBs, 16);
  if (!(drcSetEffect & (EB_DUCK_OTHER | EB_DUCK_SELF))) {
    limiterPeakTargetPresent = FDKreadBits(hBs, 1);
    if (limiterPeakTargetPresent) {
      FDKpushFor(hBs, 8);
    }
  }

  drcSetTargetLoudnessPresent = FDKreadBits(hBs, 1);
  if (drcSetTargetLoudnessPresent) {
    FDKpushFor(hBs, 6);
    drcSetTargetLoudnessValueLowerPresent = FDKreadBits(hBs, 1);
    if (drcSetTargetLoudnessValueLowerPresent) {
      FDKpushFor(hBs, 6);
    }
  }
}

static DRC_ERROR _readDrcInstructionsUniDrc(HANDLE_FDK_BITSTREAM hBs,
                                            const int version,
                                            HANDLE_UNI_DRC_CONFIG hUniDrcConfig,
                                            DRC_INSTRUCTIONS_UNI_DRC* pInst) {
  DRC_ERROR err = DE_OK;
  int i, g, c;
  int downmixIdPresent, additionalDownmixIdPresent, additionalDownmixIdCount;
  int bsLimiterPeakTarget, channelCount;
  DRC_COEFFICIENTS_UNI_DRC* pCoef = NULL;
  int repeatParameters, bsRepeatParametersCount;
  int repeatSequenceIndex, bsRepeatSequenceCount;
  SCHAR* gainSetIndex = pInst->gainSetIndex;
  SCHAR channelGroupForChannel[8];
  DUCKING_MODIFICATION duckingModificationForChannelGroup[8];

  pInst->drcSetId = FDKreadBits(hBs, 6);
  if (version == 0) {
    pInst->drcSetComplexityLevel = 2;
  } else {
    pInst->drcSetComplexityLevel = FDKreadBits(hBs, 4);
  }
  pInst->drcLocation = FDKreadBits(hBs, 4);
  if (version == 0) {
    downmixIdPresent = 1;
  } else {
    downmixIdPresent = FDKreadBits(hBs, 1);
  }
  if (downmixIdPresent) {
    pInst->downmixId[0] = FDKreadBits(hBs, 7);
    if (version == 0) {
      if (pInst->downmixId[0] == 0)
        pInst->drcApplyToDownmix = 0;
      else
        pInst->drcApplyToDownmix = 1;
    } else {
      pInst->drcApplyToDownmix = FDKreadBits(hBs, 1);
    }

    additionalDownmixIdPresent = FDKreadBits(hBs, 1);
    if (additionalDownmixIdPresent) {
      additionalDownmixIdCount = FDKreadBits(hBs, 3);
      if ((1 + additionalDownmixIdCount) > 8) return DE_MEMORY_ERROR;
      for (i = 0; i < additionalDownmixIdCount; i++) {
        pInst->downmixId[i + 1] = FDKreadBits(hBs, 7);
      }
      pInst->downmixIdCount = 1 + additionalDownmixIdCount;
    } else {
      pInst->downmixIdCount = 1;
    }
  } else {
    pInst->downmixId[0] = 0;
    pInst->downmixIdCount = 1;
  }

  pInst->drcSetEffect = FDKreadBits(hBs, 16);

  if ((pInst->drcSetEffect & (EB_DUCK_OTHER | EB_DUCK_SELF)) == 0) {
    pInst->limiterPeakTargetPresent = FDKreadBits(hBs, 1);
    if (pInst->limiterPeakTargetPresent) {
      bsLimiterPeakTarget = FDKreadBits(hBs, 8);
      pInst->limiterPeakTarget = -(FIXP_SGL)(
          bsLimiterPeakTarget
          << (FRACT_BITS - 1 - 3 - 5));
    }
  }

  pInst->drcSetTargetLoudnessPresent = FDKreadBits(hBs, 1);

  pInst->drcSetTargetLoudnessValueUpper = 0;
  pInst->drcSetTargetLoudnessValueLower = -63;

  if (pInst->drcSetTargetLoudnessPresent == 1) {
    int bsDrcSetTargetLoudnessValueUpper, bsDrcSetTargetLoudnessValueLower;
    int drcSetTargetLoudnessValueLowerPresent;
    bsDrcSetTargetLoudnessValueUpper = FDKreadBits(hBs, 6);
    pInst->drcSetTargetLoudnessValueUpper =
        bsDrcSetTargetLoudnessValueUpper - 63;
    drcSetTargetLoudnessValueLowerPresent = FDKreadBits(hBs, 1);
    if (drcSetTargetLoudnessValueLowerPresent == 1) {
      bsDrcSetTargetLoudnessValueLower = FDKreadBits(hBs, 6);
      pInst->drcSetTargetLoudnessValueLower =
          bsDrcSetTargetLoudnessValueLower - 63;
    }
  }

  pInst->dependsOnDrcSetPresent = FDKreadBits(hBs, 1);

  pInst->noIndependentUse = 0;
  if (pInst->dependsOnDrcSetPresent) {
    pInst->dependsOnDrcSet = FDKreadBits(hBs, 6);
  } else {
    pInst->noIndependentUse = FDKreadBits(hBs, 1);
  }

  if (version == 0) {
    pInst->requiresEq = 0;
  } else {
    pInst->requiresEq = FDKreadBits(hBs, 1);
  }

  pCoef = selectDrcCoefficients(hUniDrcConfig, pInst->drcLocation);

  pInst->drcChannelCount = channelCount =
      hUniDrcConfig->channelLayout.baseChannelCount;

  if (pInst->drcSetEffect & (EB_DUCK_OTHER | EB_DUCK_SELF)) {
    DUCKING_MODIFICATION* pDModForChannel =
        pInst->duckingModificationForChannel;
    c = 0;
    while (c < channelCount) {
      int bsGainSetIndex;
      bsGainSetIndex = FDKreadBits(hBs, 6);
      if (c >= 8) return DE_MEMORY_ERROR;
      gainSetIndex[c] = bsGainSetIndex - 1;
      _decodeDuckingModification(hBs, &(pDModForChannel[c]), 0);

      c++;
      repeatParameters = FDKreadBits(hBs, 1);
      if (repeatParameters == 1) {
        bsRepeatParametersCount = FDKreadBits(hBs, 5);
        bsRepeatParametersCount += 1;
        for (i = 0; i < bsRepeatParametersCount; i++) {
          if (c >= 8) return DE_MEMORY_ERROR;
          gainSetIndex[c] = gainSetIndex[c - 1];
          pDModForChannel[c] = pDModForChannel[c - 1];
          c++;
        }
      }
    }
    if (c > channelCount) {
      return DE_NOT_OK;
    }

    err = deriveDrcChannelGroups(
        pInst->drcSetEffect, pInst->drcChannelCount, gainSetIndex,
        pDModForChannel, &pInst->nDrcChannelGroups,
        pInst->gainSetIndexForChannelGroup, channelGroupForChannel,
        duckingModificationForChannelGroup);
    if (err) return (err);
  } else {
    int deriveChannelCount = 0;
    if (((version == 0) || (pInst->drcApplyToDownmix != 0)) &&
        (pInst->downmixId[0] != DOWNMIX_ID_BASE_LAYOUT) &&
        (pInst->downmixId[0] != DOWNMIX_ID_ANY_DOWNMIX) &&
        (pInst->downmixIdCount == 1)) {
      if (hUniDrcConfig->downmixInstructionsCount != 0) {
        DOWNMIX_INSTRUCTIONS* pDown =
            selectDownmixInstructions(hUniDrcConfig, pInst->downmixId[0]);
        if (pDown == NULL) return DE_NOT_OK;
        pInst->drcChannelCount = channelCount =
            pDown->targetChannelCount;
      } else {
        deriveChannelCount = 1;
        channelCount = 1;
      }
    } else if (((version == 0) || (pInst->drcApplyToDownmix != 0)) &&
               ((pInst->downmixId[0] == DOWNMIX_ID_ANY_DOWNMIX) ||
                (pInst->downmixIdCount > 1))) {
      pInst->drcChannelCount = 8;
      channelCount = 1;
    }

    c = 0;
    while (c < channelCount) {
      int bsGainSetIndex;
      bsGainSetIndex = FDKreadBits(hBs, 6);
      if (c >= 8) return DE_MEMORY_ERROR;
      gainSetIndex[c] = bsGainSetIndex - 1;
      c++;
      repeatSequenceIndex = FDKreadBits(hBs, 1);

      if (repeatSequenceIndex == 1) {
        bsRepeatSequenceCount = FDKreadBits(hBs, 5);
        bsRepeatSequenceCount += 1;
        if (deriveChannelCount) {
          channelCount = 1 + bsRepeatSequenceCount;
        }
        for (i = 0; i < bsRepeatSequenceCount; i++) {
          if (c >= 8) return DE_MEMORY_ERROR;
          gainSetIndex[c] = bsGainSetIndex - 1;
          c++;
        }
      }
    }
    if (c > channelCount) {
      return DE_NOT_OK;
    }
    if (deriveChannelCount) {
      pInst->drcChannelCount = channelCount;
    }

    if ((pInst->downmixId[0] == DOWNMIX_ID_ANY_DOWNMIX) ||
        (pInst->downmixIdCount > 1)) {
      for (c = 1; c < pInst->drcChannelCount; c++) {
        gainSetIndex[c] = gainSetIndex[0];
      }
    }

    err = deriveDrcChannelGroups(pInst->drcSetEffect, pInst->drcChannelCount,
                                 gainSetIndex, NULL, &pInst->nDrcChannelGroups,
                                 pInst->gainSetIndexForChannelGroup,
                                 channelGroupForChannel, NULL);
    if (err) return (err);

    for (g = 0; g < pInst->nDrcChannelGroups; g++) {
      int set, bandCount;
      set = pInst->gainSetIndexForChannelGroup[g];

      if (pCoef != NULL && set < pCoef->gainSetCount) {
        bandCount = pCoef->gainSet[set].bandCount;
      } else {
        bandCount = 1;
      }

      _decodeGainModification(hBs, version, bandCount,
                              pInst->gainModificationForChannelGroup[g], 0);
    }
  }

  return err;
}

static DRC_ERROR _readChannelLayout(HANDLE_FDK_BITSTREAM hBs,
                                    CHANNEL_LAYOUT* pChan) {
  DRC_ERROR err = DE_OK;

  pChan->baseChannelCount = FDKreadBits(hBs, 7);

  if (pChan->baseChannelCount > 8) return DE_NOT_OK;

  pChan->layoutSignalingPresent = FDKreadBits(hBs, 1);

  if (pChan->layoutSignalingPresent) {
    pChan->definedLayout = FDKreadBits(hBs, 8);

    if (pChan->definedLayout == 0) {
      int i;
      for (i = 0; i < pChan->baseChannelCount; i++) {
        if (i < 8) {
          pChan->speakerPosition[i] = FDKreadBits(hBs, 7);
        } else {
          FDKpushFor(hBs, 7);
        }
      }
    }
  }
  return err;
}

static DRC_ERROR _readDownmixInstructions(HANDLE_FDK_BITSTREAM hBs,
                                          const int version,
                                          CHANNEL_LAYOUT* pChan,
                                          DOWNMIX_INSTRUCTIONS* pDown) {
  DRC_ERROR err = DE_OK;

  pDown->downmixId = FDKreadBits(hBs, 7);
  pDown->targetChannelCount = FDKreadBits(hBs, 7);
  pDown->targetLayout = FDKreadBits(hBs, 8);
  pDown->downmixCoefficientsPresent = FDKreadBits(hBs, 1);

  if (pDown->downmixCoefficientsPresent) {
    int nDownmixCoeffs = pDown->targetChannelCount * pChan->baseChannelCount;
    int i;
    if (nDownmixCoeffs > 8 * 8) return DE_NOT_OK;
    if (version == 0) {
      pDown->bsDownmixOffset = 0;
      for (i = 0; i < nDownmixCoeffs; i++) {
        pDown->downmixCoefficient[i] = downmixCoeff[FDKreadBits(hBs, 4)];
      }
    } else {
      pDown->bsDownmixOffset = FDKreadBits(hBs, 4);
      for (i = 0; i < nDownmixCoeffs; i++) {
        pDown->downmixCoefficient[i] = downmixCoeffV1[FDKreadBits(hBs, 5)];
      }
    }
  }
  return err;
}

static DRC_ERROR _readDrcExtensionV1(HANDLE_FDK_BITSTREAM hBs,
                                     HANDLE_UNI_DRC_CONFIG hUniDrcConfig) {
  DRC_ERROR err = DE_OK;
  int downmixInstructionsV1Present;
  int drcCoeffsAndInstructionsUniDrcV1Present;
  int loudEqInstructionsPresent, loudEqInstructionsCount;
  int eqPresent, eqInstructionsCount;
  int i, offset;
  int diff = hUniDrcConfig->diff;

  downmixInstructionsV1Present = FDKreadBits(hBs, 1);
  if (downmixInstructionsV1Present == 1) {
    diff |= _compAssign(&hUniDrcConfig->downmixInstructionsCountV1,
                        FDKreadBits(hBs, 7));
    offset = hUniDrcConfig->downmixInstructionsCountV0;
    hUniDrcConfig->downmixInstructionsCount = fMin(
        (UCHAR)(offset + hUniDrcConfig->downmixInstructionsCountV1), (UCHAR)6);
    for (i = 0; i < hUniDrcConfig->downmixInstructionsCountV1; i++) {
      DOWNMIX_INSTRUCTIONS tmpDown;
      FDKmemclear(&tmpDown, sizeof(DOWNMIX_INSTRUCTIONS));
      err = _readDownmixInstructions(hBs, 1, &hUniDrcConfig->channelLayout,
                                     &tmpDown);
      if (err) return err;
      if ((offset + i) >= 6) continue;
      if (!diff)
        diff |= (FDKmemcmp(&tmpDown,
                           &(hUniDrcConfig->downmixInstructions[offset + i]),
                           sizeof(DOWNMIX_INSTRUCTIONS)) != 0);
      hUniDrcConfig->downmixInstructions[offset + i] = tmpDown;
    }
  } else {
    diff |= _compAssign(&hUniDrcConfig->downmixInstructionsCountV1, 0);
  }

  drcCoeffsAndInstructionsUniDrcV1Present = FDKreadBits(hBs, 1);
  if (drcCoeffsAndInstructionsUniDrcV1Present == 1) {
    diff |= _compAssign(&hUniDrcConfig->drcCoefficientsUniDrcCountV1,
                        FDKreadBits(hBs, 3));
    offset = hUniDrcConfig->drcCoefficientsUniDrcCountV0;
    hUniDrcConfig->drcCoefficientsUniDrcCount =
        fMin((UCHAR)(offset + hUniDrcConfig->drcCoefficientsUniDrcCountV1),
             (UCHAR)2);
    for (i = 0; i < hUniDrcConfig->drcCoefficientsUniDrcCountV1; i++) {
      DRC_COEFFICIENTS_UNI_DRC tmpCoef;
      FDKmemclear(&tmpCoef, sizeof(DRC_COEFFICIENTS_UNI_DRC));
      err = _readDrcCoefficientsUniDrc(hBs, 1, &tmpCoef);
      if (err) return err;
      if ((offset + i) >= 2) continue;
      if (!diff)
        diff |= (FDKmemcmp(&tmpCoef,
                           &(hUniDrcConfig->drcCoefficientsUniDrc[offset + i]),
                           sizeof(DRC_COEFFICIENTS_UNI_DRC)) != 0);
      hUniDrcConfig->drcCoefficientsUniDrc[offset + i] = tmpCoef;
    }

    diff |= _compAssign(&hUniDrcConfig->drcInstructionsUniDrcCountV1,
                        FDKreadBits(hBs, 6));
    offset = hUniDrcConfig->drcInstructionsUniDrcCount;
    hUniDrcConfig->drcInstructionsUniDrcCount =
        fMin((UCHAR)(offset + hUniDrcConfig->drcInstructionsUniDrcCountV1),
             (UCHAR)12);
    for (i = 0; i < hUniDrcConfig->drcInstructionsUniDrcCount; i++) {
      DRC_INSTRUCTIONS_UNI_DRC tmpInst;
      FDKmemclear(&tmpInst, sizeof(DRC_INSTRUCTIONS_UNI_DRC));
      err = _readDrcInstructionsUniDrc(hBs, 1, hUniDrcConfig, &tmpInst);
      if (err) return err;
      if ((offset + i) >= 12) continue;
      if (!diff)
        diff |= (FDKmemcmp(&tmpInst,
                           &(hUniDrcConfig->drcInstructionsUniDrc[offset + i]),
                           sizeof(DRC_INSTRUCTIONS_UNI_DRC)) != 0);
      hUniDrcConfig->drcInstructionsUniDrc[offset + i] = tmpInst;
    }
  } else {
    diff |= _compAssign(&hUniDrcConfig->drcCoefficientsUniDrcCountV1, 0);
    diff |= _compAssign(&hUniDrcConfig->drcInstructionsUniDrcCountV1, 0);
  }

  loudEqInstructionsPresent = FDKreadBits(hBs, 1);
  if (loudEqInstructionsPresent == 1) {
    loudEqInstructionsCount = FDKreadBits(hBs, 4);
    for (i = 0; i < loudEqInstructionsCount; i++) {
      _skipLoudEqInstructions(hBs);
    }
  }

  eqPresent = FDKreadBits(hBs, 1);
  if (eqPresent == 1) {
    _skipEqCoefficients(hBs);
    eqInstructionsCount = FDKreadBits(hBs, 4);
    for (i = 0; i < eqInstructionsCount; i++) {
      _skipEqInstructions(hBs, hUniDrcConfig);
    }
  }

  hUniDrcConfig->diff = diff;

  return err;
}

static DRC_ERROR _readUniDrcConfigExtension(
    HANDLE_FDK_BITSTREAM hBs, HANDLE_UNI_DRC_CONFIG hUniDrcConfig) {
  DRC_ERROR err = DE_OK;
  int k, bitSizeLen, extSizeBits, bitSize;
  INT nBitsRemaining;
  UNI_DRC_CONFIG_EXTENSION* pExt = &(hUniDrcConfig->uniDrcConfigExt);

  k = 0;
  pExt->uniDrcConfigExtType[k] = FDKreadBits(hBs, 4);
  while (pExt->uniDrcConfigExtType[k] != UNIDRCCONFEXT_TERM) {
    if (k >= (8 - 1)) return DE_MEMORY_ERROR;
    bitSizeLen = FDKreadBits(hBs, 4);
    extSizeBits = bitSizeLen + 4;

    bitSize = FDKreadBits(hBs, extSizeBits);
    pExt->extBitSize[k] = bitSize + 1;
    nBitsRemaining = (INT)FDKgetValidBits(hBs);

    switch (pExt->uniDrcConfigExtType[k]) {
      case UNIDRCCONFEXT_V1:
        err = _readDrcExtensionV1(hBs, hUniDrcConfig);
        if (err) return err;
        if (nBitsRemaining !=
            ((INT)pExt->extBitSize[k] + (INT)FDKgetValidBits(hBs)))
          return DE_NOT_OK;
        break;
      case UNIDRCCONFEXT_PARAM_DRC:
      default:
        FDKpushFor(hBs, pExt->extBitSize[k]);
        break;
    }
    k++;
    pExt->uniDrcConfigExtType[k] = FDKreadBits(hBs, 4);
  }

  return err;
}

DRC_ERROR
drcDec_readUniDrcConfig(HANDLE_FDK_BITSTREAM hBs,
                        HANDLE_UNI_DRC_CONFIG hUniDrcConfig) {
  DRC_ERROR err = DE_OK;
  int i, diff = 0;
  int drcDescriptionBasicPresent, drcCoefficientsBasicCount,
      drcInstructionsBasicCount;
  CHANNEL_LAYOUT tmpChan;
  FDKmemclear(&tmpChan, sizeof(CHANNEL_LAYOUT));
  if (hUniDrcConfig == NULL) return DE_NOT_OK;

  diff |= _compAssign(&hUniDrcConfig->sampleRatePresent, FDKreadBits(hBs, 1));

  if (hUniDrcConfig->sampleRatePresent == 1) {
    diff |=
        _compAssign(&hUniDrcConfig->sampleRate, FDKreadBits(hBs, 18) + 1000);
  }

  diff |= _compAssign(&hUniDrcConfig->downmixInstructionsCountV0,
                      FDKreadBits(hBs, 7));

  drcDescriptionBasicPresent = FDKreadBits(hBs, 1);
  if (drcDescriptionBasicPresent == 1) {
    drcCoefficientsBasicCount = FDKreadBits(hBs, 3);
    drcInstructionsBasicCount = FDKreadBits(hBs, 4);
  } else {
    drcCoefficientsBasicCount = 0;
    drcInstructionsBasicCount = 0;
  }

  diff |= _compAssign(&hUniDrcConfig->drcCoefficientsUniDrcCountV0,
                      FDKreadBits(hBs, 3));
  diff |= _compAssign(&hUniDrcConfig->drcInstructionsUniDrcCountV0,
                      FDKreadBits(hBs, 6));

  err = _readChannelLayout(hBs, &tmpChan);
  if (err) return err;

  if (!diff)
    diff |= (FDKmemcmp(&tmpChan, &hUniDrcConfig->channelLayout,
                       sizeof(CHANNEL_LAYOUT)) != 0);
  hUniDrcConfig->channelLayout = tmpChan;

  hUniDrcConfig->downmixInstructionsCount =
      fMin(hUniDrcConfig->downmixInstructionsCountV0, (UCHAR)6);
  for (i = 0; i < hUniDrcConfig->downmixInstructionsCountV0; i++) {
    DOWNMIX_INSTRUCTIONS tmpDown;
    FDKmemclear(&tmpDown, sizeof(DOWNMIX_INSTRUCTIONS));
    err = _readDownmixInstructions(hBs, 0, &hUniDrcConfig->channelLayout,
                                   &tmpDown);
    if (err) return err;
    if (i >= 6) continue;
    if (!diff)
      diff |= (FDKmemcmp(&tmpDown, &(hUniDrcConfig->downmixInstructions[i]),
                         sizeof(DOWNMIX_INSTRUCTIONS)) != 0);
    hUniDrcConfig->downmixInstructions[i] = tmpDown;
  }

  for (i = 0; i < drcCoefficientsBasicCount; i++) {
    _skipDrcCoefficientsBasic(hBs);
  }
  for (i = 0; i < drcInstructionsBasicCount; i++) {
    _skipDrcInstructionsBasic(hBs);
  }

  hUniDrcConfig->drcCoefficientsUniDrcCount =
      fMin(hUniDrcConfig->drcCoefficientsUniDrcCountV0, (UCHAR)2);
  for (i = 0; i < hUniDrcConfig->drcCoefficientsUniDrcCountV0; i++) {
    DRC_COEFFICIENTS_UNI_DRC tmpCoef;
    FDKmemclear(&tmpCoef, sizeof(DRC_COEFFICIENTS_UNI_DRC));
    err = _readDrcCoefficientsUniDrc(hBs, 0, &tmpCoef);
    if (err) return err;
    if (i >= 2) continue;
    if (!diff)
      diff |= (FDKmemcmp(&tmpCoef, &(hUniDrcConfig->drcCoefficientsUniDrc[i]),
                         sizeof(DRC_COEFFICIENTS_UNI_DRC)) != 0);
    hUniDrcConfig->drcCoefficientsUniDrc[i] = tmpCoef;
  }

  hUniDrcConfig->drcInstructionsUniDrcCount =
      fMin(hUniDrcConfig->drcInstructionsUniDrcCountV0, (UCHAR)12);
  for (i = 0; i < hUniDrcConfig->drcInstructionsUniDrcCountV0; i++) {
    DRC_INSTRUCTIONS_UNI_DRC tmpInst;
    FDKmemclear(&tmpInst, sizeof(DRC_INSTRUCTIONS_UNI_DRC));
    err = _readDrcInstructionsUniDrc(hBs, 0, hUniDrcConfig, &tmpInst);
    if (err) return err;
    if (i >= 12) continue;
    if (!diff)
      diff |= (FDKmemcmp(&tmpInst, &(hUniDrcConfig->drcInstructionsUniDrc[i]),
                         sizeof(DRC_INSTRUCTIONS_UNI_DRC)) != 0);
    hUniDrcConfig->drcInstructionsUniDrc[i] = tmpInst;
  }

  diff |=
      _compAssign(&hUniDrcConfig->uniDrcConfigExtPresent, FDKreadBits(hBs, 1));
  hUniDrcConfig->diff = diff;

  if (hUniDrcConfig->uniDrcConfigExtPresent == 1) {
    err = _readUniDrcConfigExtension(hBs, hUniDrcConfig);
    if (err) return err;
  }

  return err;
}

static DRC_ERROR _decodeMethodValue(HANDLE_FDK_BITSTREAM hBs,
                                    const UCHAR methodDefinition,
                                    FIXP_DBL* methodValue, INT isBox) {
  int tmp;
  FIXP_DBL val;
  switch (methodDefinition) {
    case MD_UNKNOWN_OTHER:
    case MD_PROGRAM_LOUDNESS:
    case MD_ANCHOR_LOUDNESS:
    case MD_MAX_OF_LOUDNESS_RANGE:
    case MD_MOMENTARY_LOUDNESS_MAX:
    case MD_SHORT_TERM_LOUDNESS_MAX:
      tmp = FDKreadBits(hBs, 8);
      val = FL2FXCONST_DBL(-57.75f / (float)(1 << 7)) +
            (FIXP_DBL)(
                tmp << (DFRACT_BITS - 1 - 2 - 7));
      break;
    case MD_LOUDNESS_RANGE:
      tmp = FDKreadBits(hBs, 8);
      if (tmp == 0)
        val = (FIXP_DBL)0;
      else if (tmp <= 128)
        val = (FIXP_DBL)(tmp << (DFRACT_BITS - 1 - 2 - 7));
      else if (tmp <= 204) {
        val = (FIXP_DBL)(tmp << (DFRACT_BITS - 1 - 1 - 7)) -
              FL2FXCONST_DBL(32.0f / (float)(1 << 7));
      } else {
        val = (FIXP_DBL)(tmp << (DFRACT_BITS - 1 - 8)) -
              FL2FXCONST_DBL(134.0f / (float)(1 << 8));
        val <<= 1;
      }
      break;
    case MD_MIXING_LEVEL:
      tmp = FDKreadBits(hBs, isBox ? 8 : 5);
      val = (FIXP_DBL)(tmp << (DFRACT_BITS - 1 - 7)) +
            FL2FXCONST_DBL(80.0f / (float)(1 << 7));
      break;
    case MD_ROOM_TYPE:
      tmp = FDKreadBits(hBs, isBox ? 8 : 2);
      val = (FIXP_DBL)(tmp << (DFRACT_BITS - 1 - 7));
      break;
    case MD_SHORT_TERM_LOUDNESS:
      tmp = FDKreadBits(hBs, 8);
      val = FL2FXCONST_DBL(-116.0f / (float)(1 << 7)) +
            (FIXP_DBL)(
                tmp << (DFRACT_BITS - 1 - 1 - 7));
      break;
    default:
      return DE_NOT_OK;
  }
  *methodValue = val;
  return DE_OK;
}

static DRC_ERROR _readLoudnessMeasurement(HANDLE_FDK_BITSTREAM hBs,
                                          LOUDNESS_MEASUREMENT* pMeas) {
  DRC_ERROR err = DE_OK;

  pMeas->methodDefinition = FDKreadBits(hBs, 4);
  err =
      _decodeMethodValue(hBs, pMeas->methodDefinition, &pMeas->methodValue, 0);
  if (err) return err;
  pMeas->measurementSystem = FDKreadBits(hBs, 4);
  pMeas->reliability = FDKreadBits(hBs, 2);

  return err;
}

static DRC_ERROR _readLoudnessInfo(HANDLE_FDK_BITSTREAM hBs, const int version,
                                   LOUDNESS_INFO* loudnessInfo) {
  DRC_ERROR err = DE_OK;
  int bsSamplePeakLevel, bsTruePeakLevel, i;
  int measurementCount;

  loudnessInfo->drcSetId = FDKreadBits(hBs, 6);
  if (version >= 1) {
    loudnessInfo->eqSetId = FDKreadBits(hBs, 6);
  } else {
    loudnessInfo->eqSetId = 0;
  }
  loudnessInfo->downmixId = FDKreadBits(hBs, 7);

  loudnessInfo->samplePeakLevelPresent = FDKreadBits(hBs, 1);
  if (loudnessInfo->samplePeakLevelPresent) {
    bsSamplePeakLevel = FDKreadBits(hBs, 12);
    if (bsSamplePeakLevel == 0) {
      loudnessInfo->samplePeakLevelPresent = 0;
      loudnessInfo->samplePeakLevel = (FIXP_DBL)0;
    } else {
      loudnessInfo->samplePeakLevel =
          FL2FXCONST_DBL(20.0f / (float)(1 << 7)) -
          (FIXP_DBL)(bsSamplePeakLevel << (DFRACT_BITS - 1 - 5 - 7));
    }
  }

  loudnessInfo->truePeakLevelPresent = FDKreadBits(hBs, 1);
  if (loudnessInfo->truePeakLevelPresent) {
    bsTruePeakLevel = FDKreadBits(hBs, 12);
    if (bsTruePeakLevel == 0) {
      loudnessInfo->truePeakLevelPresent = 0;
      loudnessInfo->truePeakLevel = (FIXP_DBL)0;
    } else {
      loudnessInfo->truePeakLevel =
          FL2FXCONST_DBL(20.0f / (float)(1 << 7)) -
          (FIXP_DBL)(bsTruePeakLevel << (DFRACT_BITS - 1 - 5 - 7));
    }
    loudnessInfo->truePeakLevelMeasurementSystem = FDKreadBits(hBs, 4);
    loudnessInfo->truePeakLevelReliability = FDKreadBits(hBs, 2);
  }

  measurementCount = FDKreadBits(hBs, 4);
  loudnessInfo->measurementCount = fMin(measurementCount, 8);
  for (i = 0; i < measurementCount; i++) {
    LOUDNESS_MEASUREMENT tmpMeas;
    FDKmemclear(&tmpMeas, sizeof(LOUDNESS_MEASUREMENT));
    err = _readLoudnessMeasurement(hBs, &tmpMeas);
    if (err) return err;
    if (i >= 8) continue;
    loudnessInfo->loudnessMeasurement[i] = tmpMeas;
  }

  return err;
}

static DRC_ERROR _readLoudnessInfoSetExtEq(
    HANDLE_FDK_BITSTREAM hBs, HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet) {
  DRC_ERROR err = DE_OK;
  int i, offset;
  int diff = hLoudnessInfoSet->diff;

  diff |= _compAssign(&hLoudnessInfoSet->loudnessInfoAlbumCountV1,
                      FDKreadBits(hBs, 6));
  diff |=
      _compAssign(&hLoudnessInfoSet->loudnessInfoCountV1, FDKreadBits(hBs, 6));

  offset = hLoudnessInfoSet->loudnessInfoAlbumCountV0;
  hLoudnessInfoSet->loudnessInfoAlbumCount = fMin(
      (UCHAR)(offset + hLoudnessInfoSet->loudnessInfoAlbumCountV1), (UCHAR)12);
  for (i = 0; i < hLoudnessInfoSet->loudnessInfoAlbumCountV1; i++) {
    LOUDNESS_INFO tmpLoud;
    FDKmemclear(&tmpLoud, sizeof(LOUDNESS_INFO));
    err = _readLoudnessInfo(hBs, 1, &tmpLoud);
    if (err) return err;
    if ((offset + i) >= 12) continue;
    if (!diff)
      diff |= (FDKmemcmp(&tmpLoud,
                         &(hLoudnessInfoSet->loudnessInfoAlbum[offset + i]),
                         sizeof(LOUDNESS_INFO)) != 0);
    hLoudnessInfoSet->loudnessInfoAlbum[offset + i] = tmpLoud;
  }

  offset = hLoudnessInfoSet->loudnessInfoCountV0;
  hLoudnessInfoSet->loudnessInfoCount =
      fMin((UCHAR)(offset + hLoudnessInfoSet->loudnessInfoCountV1), (UCHAR)12);
  for (i = 0; i < hLoudnessInfoSet->loudnessInfoCountV1; i++) {
    LOUDNESS_INFO tmpLoud;
    FDKmemclear(&tmpLoud, sizeof(LOUDNESS_INFO));
    err = _readLoudnessInfo(hBs, 1, &tmpLoud);
    if (err) return err;
    if ((offset + i) >= 12) continue;
    if (!diff)
      diff |=
          (FDKmemcmp(&tmpLoud, &(hLoudnessInfoSet->loudnessInfo[offset + i]),
                     sizeof(LOUDNESS_INFO)) != 0);
    hLoudnessInfoSet->loudnessInfo[offset + i] = tmpLoud;
  }
  hLoudnessInfoSet->diff = diff;
  return err;
}

static DRC_ERROR _readLoudnessInfoSetExtension(
    HANDLE_FDK_BITSTREAM hBs, HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet) {
  DRC_ERROR err = DE_OK;
  int k, bitSizeLen, extSizeBits, bitSize;
  INT nBitsRemaining;
  LOUDNESS_INFO_SET_EXTENSION* pExt = &(hLoudnessInfoSet->loudnessInfoSetExt);

  k = 0;
  pExt->loudnessInfoSetExtType[k] = FDKreadBits(hBs, 4);
  while (pExt->loudnessInfoSetExtType[k] != UNIDRCLOUDEXT_TERM) {
    if (k >= (8 - 1)) return DE_MEMORY_ERROR;
    bitSizeLen = FDKreadBits(hBs, 4);
    extSizeBits = bitSizeLen + 4;

    bitSize = FDKreadBits(hBs, extSizeBits);
    pExt->extBitSize[k] = bitSize + 1;
    nBitsRemaining = (INT)FDKgetValidBits(hBs);

    switch (pExt->loudnessInfoSetExtType[k]) {
      case UNIDRCLOUDEXT_EQ:
        err = _readLoudnessInfoSetExtEq(hBs, hLoudnessInfoSet);
        if (err) return err;
        if (nBitsRemaining !=
            ((INT)pExt->extBitSize[k] + (INT)FDKgetValidBits(hBs)))
          return DE_NOT_OK;
        break;
      default:
        FDKpushFor(hBs, pExt->extBitSize[k]);
        break;
    }
    k++;
    pExt->loudnessInfoSetExtType[k] = FDKreadBits(hBs, 4);
  }

  return err;
}

DRC_ERROR
drcDec_readLoudnessInfoSet(HANDLE_FDK_BITSTREAM hBs,
                           HANDLE_LOUDNESS_INFO_SET hLoudnessInfoSet) {
  DRC_ERROR err = DE_OK;
  int i, diff = 0;
  if (hLoudnessInfoSet == NULL) return DE_NOT_OK;

  diff |= _compAssign(&hLoudnessInfoSet->loudnessInfoAlbumCountV0,
                      FDKreadBits(hBs, 6));
  diff |=
      _compAssign(&hLoudnessInfoSet->loudnessInfoCountV0, FDKreadBits(hBs, 6));

  hLoudnessInfoSet->loudnessInfoAlbumCount =
      fMin(hLoudnessInfoSet->loudnessInfoAlbumCountV0, (UCHAR)12);
  for (i = 0; i < hLoudnessInfoSet->loudnessInfoAlbumCountV0; i++) {
    LOUDNESS_INFO tmpLoud;
    FDKmemclear(&tmpLoud, sizeof(LOUDNESS_INFO));
    err = _readLoudnessInfo(hBs, 0, &tmpLoud);
    if (err) return err;
    if (i >= 12) continue;
    if (!diff)
      diff |= (FDKmemcmp(&tmpLoud, &(hLoudnessInfoSet->loudnessInfoAlbum[i]),
                         sizeof(LOUDNESS_INFO)) != 0);
    hLoudnessInfoSet->loudnessInfoAlbum[i] = tmpLoud;
  }

  hLoudnessInfoSet->loudnessInfoCount =
      fMin(hLoudnessInfoSet->loudnessInfoCountV0, (UCHAR)12);
  for (i = 0; i < hLoudnessInfoSet->loudnessInfoCountV0; i++) {
    LOUDNESS_INFO tmpLoud;
    FDKmemclear(&tmpLoud, sizeof(LOUDNESS_INFO));
    err = _readLoudnessInfo(hBs, 0, &tmpLoud);
    if (err) return err;
    if (i >= 12) continue;
    if (!diff)
      diff |= (FDKmemcmp(&tmpLoud, &(hLoudnessInfoSet->loudnessInfo[i]),
                         sizeof(LOUDNESS_INFO)) != 0);
    hLoudnessInfoSet->loudnessInfo[i] = tmpLoud;
  }

  diff |= _compAssign(&hLoudnessInfoSet->loudnessInfoSetExtPresent,
                      FDKreadBits(hBs, 1));
  hLoudnessInfoSet->diff = diff;

  if (hLoudnessInfoSet->loudnessInfoSetExtPresent) {
    err = _readLoudnessInfoSetExtension(hBs, hLoudnessInfoSet);
    if (err) return err;
  }

  return err;
}
