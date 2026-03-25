/**************************** PCM utility library ******************************

   Author(s):   Christian Griebel

   Description: Defines functions that perform downmixing or a simple channel
                expansion in the PCM time domain.

*******************************************************************************/

#include "pcmdmx_lib.h"
#include "genericStds.h"
#include "fixpoint_math.h"
#include "FDK_core.h"
#include "version.h"

#define PCMDMX_LIB_TITLE "PCM Downmix Lib"

#define FALSE 0
#define TRUE 1

#define FIXP_DMX FIXP_SGL
#define FX_DMX2FX_DBL(x) FX_SGL2FX_DBL((FIXP_SGL)(x))
#define FX_DBL2FX_DMX(x) FX_DBL2FX_SGL(x)
#define FL2FXCONST_DMX(x) FL2FXCONST_SGL(x)
#define MAXVAL_DMX MAXVAL_SGL
#define FX_DMX2SHRT(x) ((SHORT)(x))
#define FX_DMX2FL(x) FX_DBL2FL(FX_DMX2FX_DBL(x))

#define CH_GROUP_FRONT (0)
#define CH_GROUP_SIDE (1)
#define CH_GROUP_REAR (2)
#define CH_GROUP_LFE (3)

#define CH_PLAIN_NORMAL (0)
#define CH_PLAIN_TOP (1)
#define CH_PLAIN_BOTTOM (2)

#define CENTER_FRONT_CHANNEL (0)
#define LEFT_FRONT_CHANNEL (1)
#define RIGHT_FRONT_CHANNEL (2)
#define LEFT_REAR_CHANNEL \
  (3)
#define RIGHT_REAR_CHANNEL (4)
#define LOW_FREQUENCY_CHANNEL (5)
#define LEFT_MULTIPRPS_CHANNEL (6)
#define RIGHT_MULTIPRPS_CHANNEL (7)

#define LEFT_SIDE_CHANNEL (8)
#define RIGHT_SIDE_CHANNEL (9)
#define CENTER_REAR_CHANNEL (10)
#define CENTER_FRONT_CHANNEL_TOP (11)
#define LEFT_FRONT_CHANNEL_TOP (12)
#define RIGHT_FRONT_CHANNEL_TOP (13)
#define LEFT_SIDE_CHANNEL_TOP (14)
#define RIGHT_SIDE_CHANNEL_TOP (15)
#define CENTER_SIDE_CHANNEL_TOP (16)
#define LEFT_REAR_CHANNEL_TOP (17)
#define RIGHT_REAR_CHANNEL_TOP (18)
#define CENTER_REAR_CHANNEL_TOP (19)
#define CENTER_FRONT_CHANNEL_BOTTOM (20)
#define LEFT_FRONT_CHANNEL_BOTTOM (21)
#define RIGHT_FRONT_CHANNEL_BOTTOM (22)
#define LOW_FREQUENCY_CHANNEL_2 (23)

#define ONE_CHANNEL (1)
#define TWO_CHANNEL (2)
#define SIX_CHANNEL (6)
#define EIGHT_CHANNEL (8)
#define TWENTY_FOUR_CHANNEL (24)

#define PCMDMX_THRESHOLD_MAP_HEAT_1 (0)
#define PCMDMX_THRESHOLD_MAP_HEAT_2 (20)
#define PCMDMX_THRESHOLD_MAP_HEAT_3 \
  (256)

#define SP_Z_NRM (0)
#define SP_Z_TOP (2)
#define SP_Z_BOT (-2)
#define SP_Z_LFE (-18)
#define SP_Z_MUL (8)

typedef struct {
  SCHAR x;
  SCHAR y;
  SCHAR z;
} PCM_DMX_SPEAKER_POSITION;

static const PCM_DMX_SPEAKER_POSITION spkrSlotPos[] = {
    {0, 0, SP_Z_NRM},
    {-2, 0, SP_Z_NRM},
    {2, 0, SP_Z_NRM},
    {-3, 4, SP_Z_NRM},
    {3, 4, SP_Z_NRM},
    {0, 0, SP_Z_LFE},
    {-2, 2, SP_Z_MUL},
    {2, 2, SP_Z_MUL}
};

typedef enum {
               CH_MODE_UNDEFINED = 0x0000,
               CH_MODE_1_0_0_0 = 0x0001,
               CH_MODE_2_0_0_0 = 0x0002,
               CH_MODE_3_0_0_0 = 0x0003,
               CH_MODE_2_0_1_0 = 0x0102,
               CH_MODE_2_0_0_1 = 0x1002,
               CH_MODE_3_0_1_0 = 0x0103,
               CH_MODE_2_0_2_0 = 0x0202,
               CH_MODE_2_0_1_1 = 0x1102,
               CH_MODE_4_0_0_0 = 0x0004,
               CH_MODE_3_0_2_0 = 0x0203,
               CH_MODE_2_0_2_1 = 0x1202,
               CH_MODE_3_0_1_1 = 0x1103,
               CH_MODE_3_2_0_0 = 0x0023,
               CH_MODE_5_0_0_0 = 0x0005,
               CH_MODE_3_0_2_1 = 0x1203,
               CH_MODE_3_2_0_1 = 0x1023,
               CH_MODE_3_2_1_0 = 0x0123,
               CH_MODE_5_0_1_0 = 0x0105,
               CH_MODE_6_0_0_0 = 0x0006,
               CH_MODE_2_2_2_1 = 0x1222,
               CH_MODE_3_0_3_1 = 0x1303,
               CH_MODE_3_2_1_1 = 0x1123,
               CH_MODE_3_2_2_0 = 0x0223,
               CH_MODE_3_0_2_2 = 0x2203,
               CH_MODE_5_0_2_0 = 0x0205,
               CH_MODE_5_0_1_1 = 0x1105,
               CH_MODE_7_0_0_0 = 0x0007,
               CH_MODE_3_2_2_1 = 0x1223,
               CH_MODE_3_0_4_1 = 0x1403,
               CH_MODE_5_0_2_1 = 0x1205,
               CH_MODE_5_2_1_0 = 0x0125,
               CH_MODE_3_2_1_2 = 0x2123,
               CH_MODE_2_2_2_2 = 0x2222,
               CH_MODE_3_0_3_2 = 0x2303,
               CH_MODE_8_0_0_0 = 0x0008

} PCM_DMX_CHANNEL_MODE;

static const PCM_DMX_CHANNEL_MODE outChModeTable[(8) + 1] = {
    CH_MODE_UNDEFINED,
    CH_MODE_1_0_0_0,
    CH_MODE_2_0_0_0
    ,
    CH_MODE_3_0_0_0,
    CH_MODE_3_0_1_0,
    CH_MODE_3_0_2_0,
    CH_MODE_3_0_2_1
    ,
    CH_MODE_3_0_3_1,
    CH_MODE_3_0_4_1
};

static const FIXP_DMX abMixLvlValueTab[8] = {
    FL2FXCONST_DMX(0.500f),
    FL2FXCONST_DMX(0.841f), FL2FXCONST_DMX(0.707f), FL2FXCONST_DMX(0.596f),
    FL2FXCONST_DMX(0.500f), FL2FXCONST_DMX(0.422f), FL2FXCONST_DMX(0.355f),
    FL2FXCONST_DMX(0.0f)};

static const FIXP_DMX lfeMixLvlValueTab[16] = {
    FL2FXCONST_DMX(0.7905f), /*     2 */
    FL2FXCONST_DMX(0.5000f), /*     2 */
    FL2FXCONST_DMX(0.8395f), /*     1 */
    FL2FXCONST_DMX(0.7065f), /*     1 */
    FL2FXCONST_DMX(0.5945f), /*     1 */
    FL2FXCONST_DMX(0.500f),  /*     1 */
    FL2FXCONST_DMX(0.841f),  /*     0 */
    FL2FXCONST_DMX(0.707f),  /*     0 */
    FL2FXCONST_DMX(0.596f),  /*     0 */
    FL2FXCONST_DMX(0.500f),  /*     0 */
    FL2FXCONST_DMX(0.316f),  /*     0 */
    FL2FXCONST_DMX(0.178f),  /*     0 */
    FL2FXCONST_DMX(0.100f),  /*     0 */
    FL2FXCONST_DMX(0.032f),  /*     0 */
    FL2FXCONST_DMX(0.010f),  /*     0 */
    FL2FXCONST_DMX(0.000f)   /*     0 */
};

static const FIXP_DMX mpegMixDownIdx2Coef[4] = {
    FL2FXCONST_DMX(0.70710678f), FL2FXCONST_DMX(0.5f),
    FL2FXCONST_DMX(0.35355339f), FL2FXCONST_DMX(0.0f)};

static const FIXP_DMX mpegMixDownIdx2PreFact[3][4] = {
    {/* Set 1: */
     FL2FXCONST_DMX(0.4142135623730950f), FL2FXCONST_DMX(0.4530818393219728f),
     FL2FXCONST_DMX(0.4852813742385703f), FL2FXCONST_DMX(0.5857864376269050f)},
    {/* Set 2: */
     FL2FXCONST_DMX(0.3203772410170407f), FL2FXCONST_DMX(0.3693980625181293f),
     FL2FXCONST_DMX(0.4142135623730950f), FL2FXCONST_DMX(0.5857864376269050f)},
    {/* Mono DMX set: */
     FL2FXCONST_DMX(0.2265409196609864f), FL2FXCONST_DMX(0.25f),
     FL2FXCONST_DMX(0.2697521433898179f), FL2FXCONST_DMX(0.3333333333333333f)}};

#define TYPE_NONE (0x00)
#define TYPE_PCE_DATA (0x01)
#define TYPE_DSE_CLEV_DATA (0x02)
#define TYPE_DSE_SLEV_DATA (0x04)
#define TYPE_DSE_DMIX_AB_DATA (0x08)
#define TYPE_DSE_DMIX_LFE_DATA (0x10)
#define TYPE_DSE_DMX_GAIN_DATA (0x20)
#define TYPE_DSE_DMX_CGL_DATA (0x40)
#define TYPE_DSE_DATA (0x7E)

typedef struct {
  UINT typeFlags;
  UCHAR cLevIdx;
  UCHAR sLevIdx;
  UCHAR dmixIdxA;
  UCHAR dmixIdxB;
  UCHAR dmixIdxLfe;
  UCHAR dmxGainIdx2;
  UCHAR dmxGainIdx5;
  UCHAR matrixMixdownIdx;
  SCHAR pseudoSurround;
  UINT expiryCount;

} DMX_BS_META_DATA;

static const DMX_BS_META_DATA dfltMetaData = {0, 2, 2, 2,  2, 15,
                                              0, 0, 0, -1, 0};

typedef struct {
  DMX_PROFILE_TYPE dmxProfile;
  UINT expiryFrame;
  DUAL_CHANNEL_MODE dualChannelMode;
  PSEUDO_SURROUND_MODE
  pseudoSurrMode;
  SHORT numOutChannelsMin;
  SHORT numOutChannelsMax;
  UCHAR frameDelay;

} PCM_DMX_USER_PARAMS;

struct PCM_DMX_INSTANCE {
  DMX_BS_META_DATA bsMetaData[(1) + 1];
  PCM_DMX_USER_PARAMS userParams;
  UCHAR applyProcessing;
};

C_ALLOC_MEM(PcmDmxInstance, struct PCM_DMX_INSTANCE, 1)

static UINT getSpeakerDistance(PCM_DMX_SPEAKER_POSITION posA,
                               PCM_DMX_SPEAKER_POSITION posB) {
  PCM_DMX_SPEAKER_POSITION diff;

  diff.x = posA.x - posB.x;
  diff.y = posA.y - posB.y;
  diff.z = posA.z - posB.z;

  return ((diff.x * diff.x) + (diff.y * diff.y) + (diff.z * diff.z));
}

static PCM_DMX_SPEAKER_POSITION getSpeakerPos(AUDIO_CHANNEL_TYPE chType,
                                              UCHAR chIndex, UCHAR numChInGrp) {
#define PCMDMX_SPKR_POS_X_MAX_WIDTH (3)
#define PCMDMX_SPKR_POS_Y_SPREAD (2)
#define PCMDMX_SPKR_POS_Z_SPREAD (2)

  PCM_DMX_SPEAKER_POSITION spkrPos = {0, 0, 0};
  AUDIO_CHANNEL_TYPE chGrp = (AUDIO_CHANNEL_TYPE)(chType & 0x0F);
  unsigned fHasCenter = numChInGrp & 0x1;
  unsigned chGrpWidth = numChInGrp >> 1;
  unsigned fIsCenter = 0;
  unsigned fIsLfe = (chType == ACT_LFE) ? 1 : 0;
  int offset = 0;

  FDK_ASSERT(chIndex < numChInGrp);

  if ((chGrp == ACT_FRONT) && fHasCenter) {
    if (chIndex == 0) fIsCenter = 1;
    chIndex = (UCHAR)fMax(0, chIndex - 1);
  } else if (fHasCenter && (chIndex == numChInGrp - 1)) {
    fIsCenter = 1;
  }
  if (!fIsCenter) {
    offset = chIndex >> 1;
    if ((chGrp > ACT_FRONT) && (chType != ACT_SIDE) && !fIsLfe) {
      offset = chGrpWidth - fHasCenter - offset;
    }
    if ((chIndex & 0x1) == 0) {
      offset = -(offset + 1);
    } else {
      offset += 1;
    }
  }
  if (chType == ACT_SIDE) {
    spkrPos.x = (offset < 0) ? -PCMDMX_SPKR_POS_X_MAX_WIDTH
                             : PCMDMX_SPKR_POS_X_MAX_WIDTH;
    spkrPos.y = PCMDMX_SPKR_POS_Y_SPREAD + (SCHAR)fAbs(offset) - 1;
    spkrPos.z = 0;
  } else {
    unsigned spread =
        ((chGrpWidth == 1) && (!fIsLfe)) ? PCMDMX_SPKR_POS_X_MAX_WIDTH - 1 : 1;
    spkrPos.x = (SCHAR)offset * (SCHAR)spread;
    if (fIsLfe) {
      spkrPos.y = 0;
      spkrPos.z = SP_Z_LFE;
    } else {
      spkrPos.y = (SCHAR)fMax((SCHAR)chGrp - 1, 0) * PCMDMX_SPKR_POS_Y_SPREAD;
      spkrPos.z = (SCHAR)chType >> 4;
      if (spkrPos.z == 2) {
        spkrPos.z = -1;
      }
      spkrPos.z *= PCMDMX_SPKR_POS_Z_SPREAD;
    }
  }
  return spkrPos;
}

static PCM_DMX_CHANNEL_MODE getChMode4Plain(
    const int plainIndex, const PCM_DMX_CHANNEL_MODE totChMode,
    const int chCfg) {
  PCM_DMX_CHANNEL_MODE plainChMode = totChMode;

  switch (totChMode) {
    case CH_MODE_5_0_2_1:
      if (chCfg == 14) {
        switch (plainIndex) {
          case CH_PLAIN_BOTTOM:
            plainChMode = (PCM_DMX_CHANNEL_MODE)0x0000;
            break;
          case CH_PLAIN_TOP:
            plainChMode = CH_MODE_2_0_0_0;
            break;
          case CH_PLAIN_NORMAL:
          default:
            plainChMode = CH_MODE_3_0_2_1;
            break;
        }
      }
      break;
    default:
      break;
  }

  return plainChMode;
}

static UINT validateIndices(UINT numChannels, UINT numChannelsPlaneAndGrp,
                            AUDIO_CHANNEL_TYPE aChType,
                            const AUDIO_CHANNEL_TYPE channelType[],
                            const UCHAR channelIndices[]) {
  for (UINT reqValue = 0; reqValue < numChannelsPlaneAndGrp; reqValue++) {
    int found = FALSE;
    for (UINT i = 0; i < numChannels; i++) {
      if (channelType[i] == aChType) {
        if (channelIndices[i] == reqValue) {
          if (found == TRUE) {
            return 0;
          } else {
            found = TRUE;
          }
        }
      }
    }
    if (found == FALSE) {
      return 0;
    }
  }
  return 1;
}

static PCMDMX_ERROR getChannelMode(
    const UINT numChannels,
    const AUDIO_CHANNEL_TYPE channelType[],
    UCHAR channelIndices[],
    UCHAR offsetTable[(8)],
    PCM_DMX_CHANNEL_MODE *chMode
) {
  UCHAR numCh[(3)][(4)];
  UCHAR mapped[(8)];
  PCM_DMX_SPEAKER_POSITION spkrPos[(8)];
  PCMDMX_ERROR err = PCMDMX_OK;
  unsigned ch, numMappedInChs = 0;
  unsigned startSlot;
  unsigned stopSlot = LOW_FREQUENCY_CHANNEL;

  FDK_ASSERT(channelType != NULL);
  FDK_ASSERT(channelIndices != NULL);
  FDK_ASSERT(offsetTable != NULL);
  FDK_ASSERT(chMode != NULL);

  FDKmemclear(numCh, (3) * (4) * sizeof(UCHAR));
  FDKmemclear(mapped, (8) * sizeof(UCHAR));
  FDKmemclear(spkrPos, (8) * sizeof(PCM_DMX_SPEAKER_POSITION));
  FDKmemset(offsetTable, 255, (8) * sizeof(UCHAR));
  *chMode = CH_MODE_UNDEFINED;

  for (ch = 0; ch < numChannels; ch += 1) {
    unsigned chGrp = fMax(
        (channelType[ch] & 0x0F) - 1,
        0);
    numCh[channelType[ch] >> 4][chGrp] += 1;
  }

  {
    int chGrp;
    for (chGrp = 0; chGrp < (4); chGrp += 1) {
      int plane;
      for (plane = 0; plane < (3); plane += 1) {
        if (numCh[plane][chGrp] == 0) continue;
        AUDIO_CHANNEL_TYPE aChType =
            (AUDIO_CHANNEL_TYPE)((plane << 4) | ((chGrp + 1) & 0xF));
        if (!validateIndices(numChannels, numCh[plane][chGrp], aChType,
                             channelType, channelIndices)) {
          unsigned idxCnt = 0;
          for (ch = 0; ch < numChannels; ch += 1) {
            if (channelType[ch] == aChType) {
              channelIndices[ch] = idxCnt++;
            }
          }
          err = PCMDMX_INVALID_CH_CONFIG;
        }
      }
    }
  }
  for (ch = 0; ch < numChannels; ch += 1) {
    UINT mapDist = (unsigned)-1;
    unsigned mapCh, mapPos = (unsigned)-1;
    unsigned chGrp = fMax(
        (channelType[ch] & 0x0F) - 1,
        0);

    spkrPos[ch] = getSpeakerPos(channelType[ch], channelIndices[ch],
                                numCh[channelType[ch] >> 4][chGrp]);

    for (mapCh = 0; mapCh <= stopSlot; mapCh += 1) {
      if (offsetTable[mapCh] == 255) {
        UINT dist = getSpeakerDistance(spkrPos[ch], spkrSlotPos[mapCh]);
        if (dist < mapDist) {
          mapPos = mapCh;
          mapDist = dist;
        }
      }
    }
    if (mapDist <= PCMDMX_THRESHOLD_MAP_HEAT_1) {
      offsetTable[mapPos] = (UCHAR)ch;
      mapped[ch] = 1;
      numMappedInChs += 1;
    }
  }

  startSlot =
      ((numCh[CH_PLAIN_NORMAL][CH_GROUP_FRONT] & 0x1) || (numChannels >= (8)))
          ? 0
          : 1;
  for (ch = 0; ch < (unsigned)numChannels; ch += 1) {
    if (!mapped[ch]) {
      UINT mapDist = (unsigned)-1;
      unsigned mapCh, mapPos = (unsigned)-1;

      for (mapCh = startSlot; mapCh <= stopSlot; mapCh += 1) {
        if (offsetTable[mapCh] == 255) {
          UINT dist = getSpeakerDistance(spkrPos[ch], spkrSlotPos[mapCh]);
          if (dist < mapDist) {
            mapPos = mapCh;
            mapDist = dist;
          }
        }
      }
      if ((mapPos <= stopSlot) && (mapDist < PCMDMX_THRESHOLD_MAP_HEAT_2) &&
          (((spkrPos[ch].x != 0) && (spkrSlotPos[mapPos].x != 0))
           || ((spkrPos[ch].x == 0) &&
               (spkrSlotPos[mapPos].x ==
                0)))) {
        offsetTable[mapPos] = (UCHAR)ch;
        mapped[ch] = 1;
        numMappedInChs += 1;
      }
    }
  }

  for (ch = startSlot; (ch < (8)) && (numMappedInChs < numChannels); ch += 1) {
    if (offsetTable[ch] == 255) {
      UINT mapDist = (unsigned)-1;
      unsigned mapCh, mapPos = (unsigned)-1;

      for (mapCh = 0; mapCh < (unsigned)numChannels; mapCh += 1) {
        if (!mapped[mapCh]) {
          UINT dist = getSpeakerDistance(spkrPos[mapCh], spkrSlotPos[ch]);
          if (dist < mapDist) {
            mapPos = mapCh;
            mapDist = dist;
          }
        }
      }
      if (mapDist < PCMDMX_THRESHOLD_MAP_HEAT_3) {
        offsetTable[ch] = (UCHAR)mapPos;
        mapped[mapPos] = 1;
        numMappedInChs += 1;
        if ((spkrPos[mapPos].x == 0) && (spkrSlotPos[ch].x != 0) &&
            (numChannels <
             (8))) {
          ch += 1;
        }
      }
    }
  }

  for (ch = 0; ch < (4); ch += 1) {
    int plane, numChInGrp = 0;
    for (plane = 0; plane < (3); plane += 1) {
      numChInGrp += numCh[plane][ch];
    }
    *chMode = (PCM_DMX_CHANNEL_MODE)(*chMode | (numChInGrp << (ch * 4)));
  }

  return err;
}

static void getChannelDescription(
    const PCM_DMX_CHANNEL_MODE chMode,
    const FDK_channelMapDescr *const mapDescr,
    AUDIO_CHANNEL_TYPE channelType[],
    UCHAR channelIndices[],
    UCHAR offsetTable[(8)]
) {
  int grpIdx, plainIdx, numPlains = 1, numTotalChannels = 0;
  int chCfg, ch = 0;

  FDK_ASSERT(channelType != NULL);
  FDK_ASSERT(channelIndices != NULL);
  FDK_ASSERT(mapDescr != NULL);
  FDK_ASSERT(offsetTable != NULL);

  FDKmemclear(channelType, (8) * sizeof(AUDIO_CHANNEL_TYPE));
  FDKmemclear(channelIndices, (8) * sizeof(UCHAR));
  FDKmemset(offsetTable, 255, (8) * sizeof(UCHAR));

  for (grpIdx = 0; grpIdx < (4); grpIdx += 1) {
    numTotalChannels += (chMode >> (grpIdx * 4)) & 0xF;
  }

  switch (chMode) {
    case CH_MODE_1_0_0_0:
    case CH_MODE_2_0_0_0:
    case CH_MODE_3_0_0_0:
    case CH_MODE_3_0_1_0:
    case CH_MODE_3_0_2_0:
    case CH_MODE_3_0_2_1:
      chCfg = numTotalChannels;
      break;
    case CH_MODE_3_0_3_1:
      chCfg = 11;
      break;
    case CH_MODE_3_0_4_1:
      chCfg = 12;
      break;
    case CH_MODE_5_0_2_1:
      chCfg = 7;
      break;
    default:
      chCfg = 0;
      break;
  }

  for (plainIdx = 0; plainIdx < numPlains; plainIdx += 1) {
    PCM_DMX_CHANNEL_MODE plainChMode;
    UCHAR numChInGrp[(4)];

    plainChMode = getChMode4Plain(plainIdx, chMode, chCfg);

    numChInGrp[CH_GROUP_FRONT] = plainChMode & 0xF;
    numChInGrp[CH_GROUP_SIDE] = (plainChMode >> 4) & 0xF;
    numChInGrp[CH_GROUP_REAR] = (plainChMode >> 8) & 0xF;
    numChInGrp[CH_GROUP_LFE] = (plainChMode >> 12) & 0xF;

    if ((numChInGrp[CH_GROUP_FRONT] & 0x1) && (plainIdx == CH_PLAIN_NORMAL)) {
      int mappedIdx = FDK_chMapDescr_getMapValue(mapDescr, (UCHAR)ch, chCfg);
      offsetTable[CENTER_FRONT_CHANNEL] = (UCHAR)mappedIdx;
      channelType[mappedIdx] = ACT_FRONT;
      channelIndices[mappedIdx] = 0;
      ch += 1;
    }

    for (grpIdx = 0; grpIdx < (4); grpIdx += 1) {
      AUDIO_CHANNEL_TYPE type = ACT_NONE;
      int chMapPos = 0, maxChannels = 0;
      int chIdx = 0;

      switch (grpIdx) {
        case CH_GROUP_FRONT:
          type = (AUDIO_CHANNEL_TYPE)((plainIdx << 4) | ACT_FRONT);
          switch (plainIdx) {
            default:
              chMapPos = LEFT_FRONT_CHANNEL;
              chIdx = numChInGrp[grpIdx] & 0x1;
              break;
          }
          maxChannels = 3;
          break;
        case CH_GROUP_SIDE:
          type = (AUDIO_CHANNEL_TYPE)((plainIdx << 4) | ACT_SIDE);
          if (plainIdx == CH_PLAIN_TOP) {
            chMapPos = LEFT_SIDE_CHANNEL_TOP;
            maxChannels = 3;
          } else {
            chMapPos = LEFT_MULTIPRPS_CHANNEL;
            maxChannels = 2;
          }
          break;
        case CH_GROUP_REAR:
          type = (AUDIO_CHANNEL_TYPE)((plainIdx << 4) | ACT_BACK);
          if (plainIdx == CH_PLAIN_TOP) {
            chMapPos = LEFT_REAR_CHANNEL_TOP;
            maxChannels = 3;
          } else {
            chMapPos = LEFT_REAR_CHANNEL;
            maxChannels = 2;
          }
          break;
        case CH_GROUP_LFE:
          if (plainIdx == CH_PLAIN_NORMAL) {
            type = ACT_LFE;
            chMapPos = LOW_FREQUENCY_CHANNEL;
            maxChannels = 1;
          }
          break;
        default:
          break;
      }

      for (; chIdx < numChInGrp[grpIdx]; chIdx += 1) {
        int mappedIdx = FDK_chMapDescr_getMapValue(mapDescr, (UCHAR)ch, chCfg);
        if ((chIdx == maxChannels) || (offsetTable[chMapPos] < 255)) {
          if (offsetTable[LEFT_MULTIPRPS_CHANNEL] ==
              255) {
            chMapPos = LEFT_MULTIPRPS_CHANNEL;
          } else {
            FDK_ASSERT(0);
          }
        }
        offsetTable[chMapPos] = (UCHAR)mappedIdx;
        channelType[mappedIdx] = type;
        channelIndices[mappedIdx] = (UCHAR)chIdx;
        chMapPos += 1;
        ch += 1;
      }
    }
  }
}

static void dmxInitChannel(FIXP_DMX mixFactors[(8)][(8)],
                           INT mixScales[(8)][(8)], const unsigned int outCh) {
  unsigned int inCh;
  for (inCh = 0; inCh < (8); inCh += 1) {
    if (inCh == outCh) {
      mixFactors[outCh][inCh] = FL2FXCONST_DMX(0.5f);
      mixScales[outCh][inCh] = 1;
    } else {
      mixFactors[outCh][inCh] = FL2FXCONST_DMX(0.0f);
      mixScales[outCh][inCh] = 0;
    }
  }
}

static void dmxClearChannel(FIXP_DMX mixFactors[(8)][(8)],
                            INT mixScales[(8)][(8)], const unsigned int outCh) {
  FDK_ASSERT((outCh >= 0) && (outCh < (8)));
  FDKmemclear(&mixFactors[outCh], (8) * sizeof(FIXP_DMX));
  FDKmemclear(&mixScales[outCh], (8) * sizeof(INT));
}

static void dmxSetChannel(FIXP_DMX mixFactors[(8)][(8)],
                          INT mixScales[(8)][(8)], const unsigned int dstCh,
                          const unsigned int srcCh, const FIXP_DMX factor,
                          const INT scale) {
  int ch;
  for (ch = 0; ch < (8); ch += 1) {
    if (mixFactors[srcCh][ch] != (FIXP_DMX)0) {
      mixFactors[dstCh][ch] =
          FX_DBL2FX_DMX(fMult(mixFactors[srcCh][ch], factor));
      mixScales[dstCh][ch] = mixScales[srcCh][ch] + scale;
    }
  }
}

static void dmxAddChannel(FIXP_DMX mixFactors[(8)][(8)],
                          INT mixScales[(8)][(8)], const unsigned int dstCh,
                          const unsigned int srcCh, const FIXP_DMX factor,
                          const INT scale) {
  int ch;
  for (ch = 0; ch < (8); ch += 1) {
    FIXP_DBL addFact = fMult(mixFactors[srcCh][ch], factor);
    if (addFact != (FIXP_DMX)0) {
      INT newScale = mixScales[srcCh][ch] + scale;
      if (mixFactors[dstCh][ch] != (FIXP_DMX)0) {
        if (newScale > mixScales[dstCh][ch]) {
          mixFactors[dstCh][ch] >>= newScale - mixScales[dstCh][ch];
        } else {
          addFact >>= mixScales[dstCh][ch] - newScale;
          newScale = mixScales[dstCh][ch];
        }
      }
      mixFactors[dstCh][ch] += FX_DBL2FX_DMX(addFact);
      mixScales[dstCh][ch] = newScale;
    }
  }
}

static PCMDMX_ERROR getMixFactors(const UCHAR inModeIsCfg,
                                  PCM_DMX_CHANNEL_MODE inChMode,
                                  const PCM_DMX_CHANNEL_MODE outChMode,
                                  const PCM_DMX_USER_PARAMS *pParams,
                                  const DMX_BS_META_DATA *pMetaData,
                                  FIXP_DMX mixFactors[(8)][(8)],
                                  INT *pOutScale) {
  PCMDMX_ERROR err = PCMDMX_OK;
  INT mixScales[(8)][(8)];
  INT maxScale = 0;
  int numInChannel;
  int numOutChannel;
  int dmxMethod;
  unsigned int outCh, inChCfg = 0;
  unsigned int valid[(8)] = {0};

  FDK_ASSERT(pMetaData != NULL);
  FDK_ASSERT(mixFactors != NULL);
  if (!((outChMode == CH_MODE_1_0_0_0) || (outChMode == CH_MODE_2_0_0_0) ||
        (outChMode == CH_MODE_3_0_2_1) || (outChMode == CH_MODE_3_0_4_1) ||
        (outChMode == CH_MODE_5_0_2_1))) {
    FDK_ASSERT(0);
  }

  if (inModeIsCfg) {
    inChCfg = (unsigned int)inChMode;
    switch (inChCfg) {
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
        inChMode = outChModeTable[inChCfg];
        break;
      case 11:
        inChMode = CH_MODE_3_0_3_1;
        break;
      case 12:
        inChMode = CH_MODE_3_0_4_1;
        break;
      case 7:
      case 14:
        inChMode = CH_MODE_5_0_2_1;
        break;
      default:
        FDK_ASSERT(0);
    }
  }

  numInChannel = (inChMode & 0xF) + ((inChMode >> 4) & 0xF) +
                 ((inChMode >> 8) & 0xF) + ((inChMode >> 12) & 0xF);
  numOutChannel = (outChMode & 0xF) + ((outChMode >> 4) & 0xF) +
                  ((outChMode >> 8) & 0xF) + ((outChMode >> 12) & 0xF);

  for (outCh = 0; outCh < (8); outCh += 1) {
    dmxInitChannel(mixFactors, mixScales, outCh);
  }
  if (((inChMode >> 12) & 0xF) == 0) {
    dmxClearChannel(mixFactors, mixScales, LOW_FREQUENCY_CHANNEL);
  }

  if (numInChannel > SIX_CHANNEL) {
    FIXP_DMX dMixFactA, dMixFactB;
    INT dMixScaleA, dMixScaleB;
    int isValidCfg = TRUE;

    dMixFactA = abMixLvlValueTab[pMetaData->dmixIdxA];
    dMixScaleA = (pMetaData->dmixIdxA == 0) ? 1 : 0;
    dMixFactB = abMixLvlValueTab[pMetaData->dmixIdxB];
    dMixScaleB = (pMetaData->dmixIdxB == 0) ? 1 : 0;

    switch (inChMode) {
      case CH_MODE_3_2_1_1:
      case CH_MODE_3_2_1_0:
        isValidCfg = FALSE;
        err = PCMDMX_INVALID_MODE;
        FDK_FALLTHROUGH;
      case CH_MODE_3_0_3_1:
        dmxClearChannel(
            mixFactors, mixScales,
            RIGHT_MULTIPRPS_CHANNEL);
        dmxSetChannel(mixFactors, mixScales, LEFT_REAR_CHANNEL,
                      LEFT_REAR_CHANNEL, dMixFactA, dMixScaleA);
        dmxSetChannel(mixFactors, mixScales, LEFT_REAR_CHANNEL,
                      LEFT_MULTIPRPS_CHANNEL, dMixFactB, dMixScaleB);
        dmxSetChannel(mixFactors, mixScales, RIGHT_REAR_CHANNEL,
                      RIGHT_REAR_CHANNEL, dMixFactA, dMixScaleA);
        dmxSetChannel(mixFactors, mixScales, RIGHT_REAR_CHANNEL,
                      LEFT_MULTIPRPS_CHANNEL, dMixFactB, dMixScaleB);
        break;
      case CH_MODE_3_0_4_1:
        dmxSetChannel(mixFactors, mixScales, LEFT_REAR_CHANNEL,
                      LEFT_REAR_CHANNEL, dMixFactA, dMixScaleA);
        dmxSetChannel(mixFactors, mixScales, LEFT_REAR_CHANNEL,
                      LEFT_MULTIPRPS_CHANNEL, dMixFactB, dMixScaleB);
        dmxSetChannel(mixFactors, mixScales, RIGHT_REAR_CHANNEL,
                      RIGHT_REAR_CHANNEL, dMixFactA, dMixScaleA);
        dmxSetChannel(mixFactors, mixScales, RIGHT_REAR_CHANNEL,
                      RIGHT_MULTIPRPS_CHANNEL, dMixFactB, dMixScaleB);
        break;
      case CH_MODE_5_0_1_0:
      case CH_MODE_5_0_1_1:
        dmxClearChannel(mixFactors, mixScales,
                        RIGHT_REAR_CHANNEL);
        dmxSetChannel(mixFactors, mixScales, RIGHT_REAR_CHANNEL,
                      LEFT_REAR_CHANNEL, FL2FXCONST_DMX(0.5f), 1);
        dmxSetChannel(mixFactors, mixScales, LEFT_REAR_CHANNEL,
                      LEFT_REAR_CHANNEL, FL2FXCONST_DMX(0.5f), 1);
        FDK_FALLTHROUGH;
      case CH_MODE_5_2_1_0:
        isValidCfg = FALSE;
        err = PCMDMX_INVALID_MODE;
        FDK_FALLTHROUGH;
      case CH_MODE_5_0_2_1:
        if (inChCfg == 14) {
          dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                        LEFT_FRONT_CHANNEL, dMixFactA, dMixScaleA);
          dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                        LEFT_MULTIPRPS_CHANNEL, dMixFactB, dMixScaleB);
          dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                        RIGHT_FRONT_CHANNEL, dMixFactA, dMixScaleA);
          dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                        RIGHT_MULTIPRPS_CHANNEL, dMixFactB, dMixScaleB);
        } else {
          dmxSetChannel(mixFactors, mixScales, CENTER_FRONT_CHANNEL,
                        LEFT_MULTIPRPS_CHANNEL, dMixFactA, dMixScaleA);
          dmxSetChannel(mixFactors, mixScales, CENTER_FRONT_CHANNEL,
                        RIGHT_MULTIPRPS_CHANNEL, dMixFactA, dMixScaleA);
          dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                        LEFT_MULTIPRPS_CHANNEL, dMixFactB, dMixScaleB);
          dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                        LEFT_FRONT_CHANNEL, FL2FXCONST_DMX(0.5f), 1);
          dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                        RIGHT_MULTIPRPS_CHANNEL, dMixFactB, dMixScaleB);
          dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                        RIGHT_FRONT_CHANNEL, FL2FXCONST_DMX(0.5f), 1);
        }
        break;
      default:
        isValidCfg = FALSE;
        err = PCMDMX_INVALID_MODE;
        break;
    }

    if ((isValidCfg == TRUE) &&
        (pMetaData->dmxGainIdx5 != 0)) {
      FIXP_DMX dmxGain;
      INT dmxScale;
      INT sign = (pMetaData->dmxGainIdx5 & 0x40) ? -1 : 1;
      INT val = pMetaData->dmxGainIdx5 & 0x3F;

      dmxGain = FX_DBL2FX_DMX(
          fLdPow(FL2FXCONST_DBL(0.830482023721841f), 2,
                 (FIXP_DBL)(sign * val * (LONG)FL2FXCONST_DBL(0.0125f)), 0,
                 &dmxScale));
      if (dmxScale < 0) {
        dmxGain >>= -dmxScale;
        dmxScale = 0;
      }

      dmxSetChannel(mixFactors, mixScales, CENTER_FRONT_CHANNEL,
                    CENTER_FRONT_CHANNEL, dmxGain, dmxScale);
      dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                    LEFT_FRONT_CHANNEL, dmxGain, dmxScale);
      dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                    RIGHT_FRONT_CHANNEL, dmxGain, dmxScale);
      dmxSetChannel(mixFactors, mixScales, LEFT_REAR_CHANNEL, LEFT_REAR_CHANNEL,
                    dmxGain, dmxScale);
      dmxSetChannel(mixFactors, mixScales, RIGHT_REAR_CHANNEL,
                    RIGHT_REAR_CHANNEL, dmxGain, dmxScale);
      dmxSetChannel(mixFactors, mixScales, LOW_FREQUENCY_CHANNEL,
                    LOW_FREQUENCY_CHANNEL, dmxGain, dmxScale);
    }

    valid[CENTER_FRONT_CHANNEL] = 1;
    valid[LEFT_FRONT_CHANNEL] = 1;
    valid[RIGHT_FRONT_CHANNEL] = 1;
    valid[LEFT_REAR_CHANNEL] = 1;
    valid[RIGHT_REAR_CHANNEL] = 1;
    valid[LOW_FREQUENCY_CHANNEL] = 1;

    inChMode = CH_MODE_3_0_2_1;
  }

#define DMX_METHOD_MPEG_AMD4 1
#define DMX_METHOD_MPEG_LEGACY 2
#define DMX_METHOD_ARIB_JAPAN 4
#define DMX_METHOD_ITU_RECOM 8
#define DMX_METHOD_CUSTOM 16

  dmxMethod = DMX_METHOD_MPEG_AMD4;

  if ((pParams->dmxProfile == DMX_PRFL_FORCE_MATRIX_MIX) &&
      (pMetaData->typeFlags & TYPE_PCE_DATA)) {
    dmxMethod = DMX_METHOD_MPEG_LEGACY;
  } else if (!(pMetaData->typeFlags &
               (TYPE_DSE_CLEV_DATA | TYPE_DSE_SLEV_DATA))) {
    switch (pParams->dmxProfile) {
      default:
      case DMX_PRFL_STANDARD:
        break;
      case DMX_PRFL_MATRIX_MIX:
      case DMX_PRFL_FORCE_MATRIX_MIX:
        if (pMetaData->typeFlags & TYPE_PCE_DATA) {
          dmxMethod = DMX_METHOD_MPEG_LEGACY;
        }
        break;
      case DMX_PRFL_ARIB_JAPAN:
        dmxMethod = DMX_METHOD_ARIB_JAPAN;
        break;
    }
  }

  if (numOutChannel <= TWO_CHANNEL) {
    switch (inChMode) {
      case CH_MODE_2_0_0_0:
        switch (pParams->dualChannelMode) {
          case CH1_MODE:
            dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                          LEFT_FRONT_CHANNEL, FL2FXCONST_DMX(0.707f), 0);
            dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                          LEFT_FRONT_CHANNEL, FL2FXCONST_DMX(0.707f), 0);
            break;
          case CH2_MODE:
            dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                          RIGHT_FRONT_CHANNEL, FL2FXCONST_DMX(0.707f), 0);
            dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                          RIGHT_FRONT_CHANNEL, FL2FXCONST_DMX(0.707f), 0);
            break;
          case MIXED_MODE:
            dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                          LEFT_FRONT_CHANNEL, FL2FXCONST_DMX(0.5f), 0);
            dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                          RIGHT_FRONT_CHANNEL, FL2FXCONST_DMX(0.5f), 0);
            dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                          LEFT_FRONT_CHANNEL, FL2FXCONST_DMX(0.5f), 0);
            dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                          RIGHT_FRONT_CHANNEL, FL2FXCONST_DMX(0.5f), 0);
            break;
          default:
          case STEREO_MODE:
            break;
        }
        break;
      case CH_MODE_2_0_1_0: {
        FIXP_DMX sMixLvl;
        if (dmxMethod == DMX_METHOD_ARIB_JAPAN) {
          dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                        LEFT_FRONT_CHANNEL, FL2FXCONST_DMX(0.707f), 0);
          dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                        RIGHT_FRONT_CHANNEL, FL2FXCONST_DMX(0.707f), 0);
          sMixLvl = FL2FXCONST_DMX(0.5f);
        } else {
          sMixLvl = FL2FXCONST_DMX(0.707f);
        }
        dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                      LEFT_REAR_CHANNEL, sMixLvl, 0);
        dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                      LEFT_REAR_CHANNEL, sMixLvl, 0);
      } break;
      case CH_MODE_3_0_0_0:
      {
        FIXP_DMX cMixLvl;
        if (dmxMethod == DMX_METHOD_ARIB_JAPAN) {
          dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                        LEFT_FRONT_CHANNEL, FL2FXCONST_DMX(0.707f), 0);
          dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                        RIGHT_FRONT_CHANNEL, FL2FXCONST_DMX(0.707f), 0);
          cMixLvl = FL2FXCONST_DMX(0.5f);
        } else {
          cMixLvl = FL2FXCONST_DMX(0.707f);
        }
        dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                      CENTER_FRONT_CHANNEL, cMixLvl, 0);
        dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                      CENTER_FRONT_CHANNEL, cMixLvl, 0);
      } break;
      case CH_MODE_3_0_1_0:
      {
        FIXP_DMX csMixLvl;
        if (dmxMethod == DMX_METHOD_ARIB_JAPAN) {
          dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                        LEFT_FRONT_CHANNEL, FL2FXCONST_DMX(0.707f), 0);
          dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                        RIGHT_FRONT_CHANNEL, FL2FXCONST_DMX(0.707f), 0);
          csMixLvl = FL2FXCONST_DMX(0.5f);
        } else {
          csMixLvl = FL2FXCONST_DMX(0.707f);
        }
        dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                      CENTER_FRONT_CHANNEL, csMixLvl, 0);
        dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                      LEFT_REAR_CHANNEL, csMixLvl, 0);
        dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                      CENTER_FRONT_CHANNEL, csMixLvl, 0);
        dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                      LEFT_REAR_CHANNEL, csMixLvl, 0);
      } break;
      case CH_MODE_3_0_2_0:
      case CH_MODE_3_0_2_1:
      {
        switch (dmxMethod) {
          default:
          case DMX_METHOD_MPEG_AMD4: {
            FIXP_DMX cMixLvl, sMixLvl, lMixLvl;
            INT cMixScale, sMixScale, lMixScale;

            cMixLvl = abMixLvlValueTab[pMetaData->cLevIdx];
            cMixScale = (pMetaData->cLevIdx == 0) ? 1 : 0;
            sMixLvl = abMixLvlValueTab[pMetaData->sLevIdx];
            sMixScale = (pMetaData->sLevIdx == 0) ? 1 : 0;
            lMixLvl = lfeMixLvlValueTab[pMetaData->dmixIdxLfe];
            if (pMetaData->dmixIdxLfe <= 1) {
              lMixScale = 2;
            } else if (pMetaData->dmixIdxLfe <= 5) {
              lMixScale = 1;
            } else {
              lMixScale = 0;
            }
            if ((pParams->pseudoSurrMode == FORCE_PS_DMX) ||
                ((pParams->pseudoSurrMode == AUTO_PS_DMX) &&
                 (pMetaData->pseudoSurround ==
                  1))) {
              dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                            CENTER_FRONT_CHANNEL, cMixLvl, cMixScale);
              dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                            LEFT_REAR_CHANNEL, -sMixLvl, sMixScale);
              dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                            RIGHT_REAR_CHANNEL, -sMixLvl, sMixScale);
              dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                            LOW_FREQUENCY_CHANNEL, lMixLvl, lMixScale);
              dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                            CENTER_FRONT_CHANNEL, cMixLvl, cMixScale);
              dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                            LEFT_REAR_CHANNEL, sMixLvl, sMixScale);
              dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                            RIGHT_REAR_CHANNEL, sMixLvl, sMixScale);
              dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                            LOW_FREQUENCY_CHANNEL, lMixLvl, lMixScale);
            } else {
              dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                            CENTER_FRONT_CHANNEL, cMixLvl, cMixScale);
              dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                            LEFT_REAR_CHANNEL, sMixLvl, sMixScale);
              dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                            LOW_FREQUENCY_CHANNEL, lMixLvl, lMixScale);
              dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                            CENTER_FRONT_CHANNEL, cMixLvl, cMixScale);
              dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                            RIGHT_REAR_CHANNEL, sMixLvl, sMixScale);
              dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                            LOW_FREQUENCY_CHANNEL, lMixLvl, lMixScale);
            }

            if (pMetaData->dmxGainIdx2 != 0) {
              FIXP_DMX dmxGain;
              INT dmxScale;
              INT sign = (pMetaData->dmxGainIdx2 & 0x40) ? -1 : 1;
              INT val = pMetaData->dmxGainIdx2 & 0x3F;

              dmxGain = FX_DBL2FX_DMX(
                  fLdPow(FL2FXCONST_DBL(0.830482023721841f), 2,
                         (FIXP_DBL)(sign * val * (LONG)FL2FXCONST_DBL(0.0125f)),
                         0, &dmxScale));
              if (dmxScale < 0) {
                dmxGain >>= -dmxScale;
                dmxScale = 0;
              }

              dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                            LEFT_FRONT_CHANNEL, dmxGain, dmxScale);
              dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                            RIGHT_FRONT_CHANNEL, dmxGain, dmxScale);
            }
          } break;
          case DMX_METHOD_ARIB_JAPAN:
          case DMX_METHOD_MPEG_LEGACY: {
            FIXP_DMX flev, clev, slevLL, slevLR, slevRL, slevRR;
            FIXP_DMX mtrxMixDwnCoef =
                mpegMixDownIdx2Coef[pMetaData->matrixMixdownIdx];

            if ((pParams->pseudoSurrMode == FORCE_PS_DMX) ||
                ((pParams->pseudoSurrMode == AUTO_PS_DMX) &&
                 (pMetaData->pseudoSurround == 1))) {
              if (dmxMethod == DMX_METHOD_ARIB_JAPAN) {
                flev = mpegMixDownIdx2Coef[0];
              } else {
                flev = mpegMixDownIdx2PreFact[1][pMetaData->matrixMixdownIdx];
              }
              slevRR = slevRL = FX_DBL2FX_DMX(fMult(flev, mtrxMixDwnCoef));
              slevLL = slevLR = -slevRL;
            } else {
              if (dmxMethod == DMX_METHOD_ARIB_JAPAN) {
                flev = mpegMixDownIdx2Coef[0];
              } else {
                flev = mpegMixDownIdx2PreFact[0][pMetaData->matrixMixdownIdx];
              }
              slevRR = slevLL = FX_DBL2FX_DMX(fMult(flev, mtrxMixDwnCoef));
              slevLR = slevRL = (FIXP_DMX)0;
            }
            clev =
                FX_DBL2FX_DMX(fMult(flev, mpegMixDownIdx2Coef[0]));

            dmxSetChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                          LEFT_FRONT_CHANNEL, flev, 0);
            dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                          CENTER_FRONT_CHANNEL, clev, 0);
            dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                          LEFT_REAR_CHANNEL, slevLL, 0);
            dmxAddChannel(mixFactors, mixScales, LEFT_FRONT_CHANNEL,
                          RIGHT_REAR_CHANNEL, slevLR, 0);

            dmxSetChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                          RIGHT_FRONT_CHANNEL, flev, 0);
            dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                          CENTER_FRONT_CHANNEL, clev, 0);
            dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                          LEFT_REAR_CHANNEL, slevRL, 0);
            dmxAddChannel(mixFactors, mixScales, RIGHT_FRONT_CHANNEL,
                          RIGHT_REAR_CHANNEL, slevRR, 0);
          } break;
        }
      } break;
      default:
        err = PCMDMX_INVALID_MODE;
        break;
    }

    FDKmemclear(valid, (8) * sizeof(unsigned int));
    valid[LEFT_FRONT_CHANNEL] = 1;
    valid[RIGHT_FRONT_CHANNEL] = 1;
  }

  if (numOutChannel == ONE_CHANNEL) {
    FIXP_DMX monoMixLevel;
    INT monoMixScale = 0;

    dmxClearChannel(mixFactors, mixScales,
                    CENTER_FRONT_CHANNEL);

    if (dmxMethod ==
        DMX_METHOD_MPEG_LEGACY) {
      monoMixLevel = mpegMixDownIdx2PreFact[2][pMetaData->matrixMixdownIdx];

      mixFactors[CENTER_FRONT_CHANNEL][CENTER_FRONT_CHANNEL] = monoMixLevel;
      mixFactors[CENTER_FRONT_CHANNEL][LEFT_FRONT_CHANNEL] = monoMixLevel;
      mixFactors[CENTER_FRONT_CHANNEL][RIGHT_FRONT_CHANNEL] = monoMixLevel;
      monoMixLevel = FX_DBL2FX_DMX(fMult(
          monoMixLevel, mpegMixDownIdx2Coef[pMetaData->matrixMixdownIdx]));
      mixFactors[CENTER_FRONT_CHANNEL][LEFT_REAR_CHANNEL] = monoMixLevel;
      mixFactors[CENTER_FRONT_CHANNEL][RIGHT_REAR_CHANNEL] = monoMixLevel;
    } else {
      switch (dmxMethod) {
        case DMX_METHOD_MPEG_AMD4:
          monoMixLevel = FL2FXCONST_DMX(0.5f);
          monoMixScale = 1;
          break;
        default:
          monoMixLevel = FL2FXCONST_DMX(0.5f);
          monoMixScale = 0;
          break;
      }
      dmxSetChannel(mixFactors, mixScales, CENTER_FRONT_CHANNEL,
                    LEFT_FRONT_CHANNEL, monoMixLevel, monoMixScale);
      dmxAddChannel(mixFactors, mixScales, CENTER_FRONT_CHANNEL,
                    RIGHT_FRONT_CHANNEL, monoMixLevel, monoMixScale);
    }

    FDKmemclear(valid, (8) * sizeof(unsigned int));
    valid[CENTER_FRONT_CHANNEL] = 1;
  }

#define MAX_SEARCH_START_VAL (-7)

  {
    LONG chSum[(8)];
    INT chSumMax = MAX_SEARCH_START_VAL;

    for (outCh = 0; outCh < (8); outCh += 1) {
      if (valid[outCh] != 0) {
        unsigned int inCh;
        for (inCh = 0; inCh < (8); inCh += 1) {
          if (mixScales[outCh][inCh] > maxScale) {
            maxScale = mixScales[outCh][inCh];
          }
        }
      }
    }

    for (outCh = 0; outCh < (8); outCh += 1) {
      chSum[outCh] = MAX_SEARCH_START_VAL;
      if (valid[outCh] != 0) {
        int ovrflwProtScale = 0;
        unsigned int inCh;

        chSum[outCh] = 0;
        for (inCh = 0; inCh < (8); inCh += 1) {
          SHORT addFact = FX_DMX2SHRT(mixFactors[outCh][inCh]);
          if (mixScales[outCh][inCh] <= maxScale) {
            addFact >>= maxScale - mixScales[outCh][inCh];
          } else {
            addFact <<= mixScales[outCh][inCh] - maxScale;
          }
          chSum[outCh] += addFact;
        }
        if (chSum[outCh] > (LONG)MAXVAL_SGL) {
          while (chSum[outCh] > (LONG)MAXVAL_SGL) {
            ovrflwProtScale += 1;
            chSum[outCh] >>= 1;
          }
        } else if (chSum[outCh] > 0) {
          while ((chSum[outCh] << 1) <= (LONG)MAXVAL_SGL) {
            ovrflwProtScale -= 1;
            chSum[outCh] <<= 1;
          }
        }
        chSum[outCh] = ovrflwProtScale;
      }
    }

    for (outCh = 0; outCh < (8); outCh += 1) {
      if ((valid[outCh] != 0) &&
          (chSum[outCh] > chSumMax)) {
        chSumMax = chSum[outCh];
      }
    }
    maxScale = fMax(maxScale + chSumMax, 0);

    for (outCh = 0; outCh < (8); outCh += 1) {
      if (valid[outCh] != 0) {
        unsigned int inCh;
        for (inCh = 0; inCh < (8); inCh += 1) {
          if (mixFactors[outCh][inCh] != (FIXP_DMX)0) {
            if (mixScales[outCh][inCh] <= maxScale) {
              mixFactors[outCh][inCh] >>= maxScale - mixScales[outCh][inCh];
            } else {
              mixFactors[outCh][inCh] <<= mixScales[outCh][inCh] - maxScale;
            }
            mixScales[outCh][inCh] = maxScale;
          }
        }
      }
    }
  }

  *pOutScale = maxScale;

  return (err);
}

PCMDMX_ERROR pcmDmx_Open(HANDLE_PCM_DOWNMIX *pSelf) {
  HANDLE_PCM_DOWNMIX self;

  if (pSelf == NULL) {
    return (PCMDMX_INVALID_HANDLE);
  }

  *pSelf = NULL;

  self = (HANDLE_PCM_DOWNMIX)GetPcmDmxInstance(0);
  if (self == NULL) {
    return (PCMDMX_OUT_OF_MEMORY);
  }

  pcmDmx_Reset(self, PCMDMX_RESET_FULL);

  *pSelf = self;

  return (PCMDMX_OK);
}

PCMDMX_ERROR pcmDmx_Reset(HANDLE_PCM_DOWNMIX self, UINT flags) {
  if (self == NULL) {
    return (PCMDMX_INVALID_HANDLE);
  }

  if (flags & PCMDMX_RESET_PARAMS) {
    PCM_DMX_USER_PARAMS *pParams = &self->userParams;

    pParams->dualChannelMode = STEREO_MODE;
    pParams->pseudoSurrMode = NEVER_DO_PS_DMX;
    pParams->numOutChannelsMax = (6);
    pParams->numOutChannelsMin = (0);
    pParams->frameDelay = 0;
    pParams->expiryFrame = (0);

    self->applyProcessing = 0;
  }

  if (flags & PCMDMX_RESET_BS_DATA) {
    int slot;
    for (slot = 0; slot <= (1); slot += 1) {
      FDKmemcpy(&self->bsMetaData[slot], &dfltMetaData,
                sizeof(DMX_BS_META_DATA));
    }
  }

  return (PCMDMX_OK);
}

PCMDMX_ERROR pcmDmx_SetParam(HANDLE_PCM_DOWNMIX self, const PCMDMX_PARAM param,
                             const INT value) {
  switch (param) {
    case DMX_PROFILE_SETTING:
      switch ((DMX_PROFILE_TYPE)value) {
        case DMX_PRFL_STANDARD:
        case DMX_PRFL_MATRIX_MIX:
        case DMX_PRFL_FORCE_MATRIX_MIX:
        case DMX_PRFL_ARIB_JAPAN:
          break;
        default:
          return (PCMDMX_UNABLE_TO_SET_PARAM);
      }
      if (self == NULL) return (PCMDMX_INVALID_HANDLE);
      self->userParams.dmxProfile = (DMX_PROFILE_TYPE)value;
      break;

    case DMX_BS_DATA_EXPIRY_FRAME:
      if (self == NULL) return (PCMDMX_INVALID_HANDLE);
      self->userParams.expiryFrame = (value > 0) ? (UINT)value : 0;
      break;

    case DMX_BS_DATA_DELAY:
      if ((value > (1)) || (value < 0)) {
        return (PCMDMX_UNABLE_TO_SET_PARAM);
      }
      if (self == NULL) {
        return (PCMDMX_INVALID_HANDLE);
      }
      self->userParams.frameDelay = (UCHAR)value;
      break;

    case MIN_NUMBER_OF_OUTPUT_CHANNELS:
      switch (value) {
        case -1:
        case 0:
        case ONE_CHANNEL:
        case TWO_CHANNEL:
        case SIX_CHANNEL:
        case EIGHT_CHANNEL:
          break;
        default:
          return (PCMDMX_UNABLE_TO_SET_PARAM);
      }
      if (self == NULL) return (PCMDMX_INVALID_HANDLE);
      self->userParams.numOutChannelsMin = (value > 0) ? (SHORT)value : -1;
      if ((value > 0) && (self->userParams.numOutChannelsMax > 0) &&
          (value > self->userParams
                       .numOutChannelsMax)) {
        self->userParams.numOutChannelsMax = self->userParams.numOutChannelsMin;
      }
      break;

    case MAX_NUMBER_OF_OUTPUT_CHANNELS:
      switch (value) {
        case -1:
        case 0:
        case ONE_CHANNEL:
        case TWO_CHANNEL:
        case SIX_CHANNEL:
        case EIGHT_CHANNEL:
          break;
        default:
          return (PCMDMX_UNABLE_TO_SET_PARAM);
      }
      if (self == NULL) return (PCMDMX_INVALID_HANDLE);
      self->userParams.numOutChannelsMax = (value > 0) ? (SHORT)value : -1;
      if ((value > 0) &&
          (value < self->userParams
                       .numOutChannelsMin)) {
        self->userParams.numOutChannelsMin = self->userParams.numOutChannelsMax;
      }
      break;

    case DMX_DUAL_CHANNEL_MODE:
      switch ((DUAL_CHANNEL_MODE)value) {
        case STEREO_MODE:
        case CH1_MODE:
        case CH2_MODE:
        case MIXED_MODE:
          break;
        default:
          return (PCMDMX_UNABLE_TO_SET_PARAM);
      }
      if (self == NULL) return (PCMDMX_INVALID_HANDLE);
      self->userParams.dualChannelMode = (DUAL_CHANNEL_MODE)value;
      self->applyProcessing = ((DUAL_CHANNEL_MODE)value != STEREO_MODE)
                                  ? 1
                                  : 0;
      break;

    case DMX_PSEUDO_SURROUND_MODE:
      switch ((PSEUDO_SURROUND_MODE)value) {
        case NEVER_DO_PS_DMX:
        case AUTO_PS_DMX:
        case FORCE_PS_DMX:
          break;
        default:
          return (PCMDMX_UNABLE_TO_SET_PARAM);
      }
      if (self == NULL) return (PCMDMX_INVALID_HANDLE);
      self->userParams.pseudoSurrMode = (PSEUDO_SURROUND_MODE)value;
      break;

    default:
      return (PCMDMX_UNKNOWN_PARAM);
  }

  return (PCMDMX_OK);
}

PCMDMX_ERROR pcmDmx_GetParam(HANDLE_PCM_DOWNMIX self, const PCMDMX_PARAM param,
                             INT *const pValue) {
  PCM_DMX_USER_PARAMS *pUsrParams;

  if ((self == NULL) || (pValue == NULL)) {
    return (PCMDMX_INVALID_HANDLE);
  }
  pUsrParams = &self->userParams;

  switch (param) {
    case DMX_PROFILE_SETTING:
      *pValue = (INT)pUsrParams->dmxProfile;
      break;
    case DMX_BS_DATA_EXPIRY_FRAME:
      *pValue = (INT)pUsrParams->expiryFrame;
      break;
    case DMX_BS_DATA_DELAY:
      *pValue = (INT)pUsrParams->frameDelay;
      break;
    case MIN_NUMBER_OF_OUTPUT_CHANNELS:
      *pValue = (INT)pUsrParams->numOutChannelsMin;
      break;
    case MAX_NUMBER_OF_OUTPUT_CHANNELS:
      *pValue = (INT)pUsrParams->numOutChannelsMax;
      break;
    case DMX_DUAL_CHANNEL_MODE:
      *pValue = (INT)pUsrParams->dualChannelMode;
      break;
    case DMX_PSEUDO_SURROUND_MODE:
      *pValue = (INT)pUsrParams->pseudoSurrMode;
      break;
    default:
      return (PCMDMX_UNKNOWN_PARAM);
  }

  return (PCMDMX_OK);
}

PCMDMX_ERROR pcmDmx_Parse(HANDLE_PCM_DOWNMIX self, HANDLE_FDK_BITSTREAM hBs,
                          UINT ancDataBits, int isMpeg2) {
  PCMDMX_ERROR errorStatus = PCMDMX_OK;

#define MAX_DSE_ANC_BYTES (16)
#define ANC_DATA_SYNC_BYTE (0xBC)

  DMX_BS_META_DATA *pBsMetaData;

  int skip4Dmx = 0, skip4Ext = 0;
  int dmxLvlAvail = 0, extDataAvail = 0;
  UINT foundNewData = 0;
  UINT minAncBits = ((isMpeg2) ? 5 : 3) * 8;

  if ((self == NULL) || (hBs == NULL)) {
    return (PCMDMX_INVALID_HANDLE);
  }

  if ((ancDataBits < minAncBits) || (ancDataBits > FDKgetValidBits(hBs))) {
    return (PCMDMX_CORRUPT_ANC_DATA);
  }

  pBsMetaData = &self->bsMetaData[0];

  if (isMpeg2) {
    FDKpushFor(hBs, 16);
  }

  if (FDKreadBits(hBs, 8) != ANC_DATA_SYNC_BYTE) {
    return (PCMDMX_CORRUPT_ANC_DATA);
  }

  FDKpushFor(hBs, 4);

  if (isMpeg2) {
    FDKreadBits(hBs, 4);
    if (FDKreadBit(hBs)) skip4Dmx += 24;
    if (FDKreadBit(hBs)) skip4Dmx += 8;
    if (FDKreadBit(hBs)) skip4Dmx += 8;
  } else {
    FDKpushFor(hBs, 2);
    pBsMetaData->pseudoSurround = (SCHAR)FDKreadBit(hBs);
    FDKpushFor(hBs, 4);
  }

  dmxLvlAvail = FDKreadBit(hBs);

  if (isMpeg2) {
    if (FDKreadBit(hBs)) skip4Ext += 16;
  } else {
    extDataAvail = FDKreadBit(hBs);
  }

  if (FDKreadBit(hBs)) skip4Ext += 16;
  if (FDKreadBit(hBs)) skip4Ext += 16;
  if (FDKreadBit(hBs)) skip4Ext += 16;

  FDKpushFor(hBs, skip4Dmx);

  if (dmxLvlAvail) {
    if (FDKreadBit(hBs)) {
      pBsMetaData->cLevIdx = (UCHAR)FDKreadBits(hBs, 3);
      foundNewData |= TYPE_DSE_CLEV_DATA;
    } else {
      FDKreadBits(hBs, 3);
    }
    if (FDKreadBit(hBs)) {
      pBsMetaData->sLevIdx = (UCHAR)FDKreadBits(hBs, 3);
      foundNewData |= TYPE_DSE_SLEV_DATA;
    } else {
      FDKreadBits(hBs, 3);
    }
  }

  FDKpushFor(hBs, skip4Ext);

  if (extDataAvail) {
    int extDmxLvlSt, extDmxGainSt, extDmxLfeSt;

    FDKreadBit(hBs);
    extDmxLvlSt = FDKreadBit(hBs);
    extDmxGainSt = FDKreadBit(hBs);
    extDmxLfeSt = FDKreadBit(hBs);
    FDKreadBits(hBs, 4);

    if (extDmxLvlSt) {
      pBsMetaData->dmixIdxA = (UCHAR)FDKreadBits(hBs, 3);
      pBsMetaData->dmixIdxB = (UCHAR)FDKreadBits(hBs, 3);
      FDKreadBits(hBs, 2);
      foundNewData |= TYPE_DSE_DMIX_AB_DATA;
    }
    if (extDmxGainSt) {
      pBsMetaData->dmxGainIdx5 = (UCHAR)FDKreadBits(hBs, 7);
      FDKreadBit(hBs);
      pBsMetaData->dmxGainIdx2 = (UCHAR)FDKreadBits(hBs, 7);
      FDKreadBit(hBs);
      foundNewData |= TYPE_DSE_DMX_GAIN_DATA;
    }
    if (extDmxLfeSt) {
      pBsMetaData->dmixIdxLfe = (UCHAR)FDKreadBits(hBs, 4);
      FDKreadBits(hBs, 4);
      foundNewData |= TYPE_DSE_DMIX_LFE_DATA;
    }
  }

  if ((INT)FDKgetValidBits(hBs) < 0) {
    errorStatus = PCMDMX_CORRUPT_ANC_DATA;
  }

  if ((errorStatus == PCMDMX_OK) && (foundNewData != 0)) {
    pBsMetaData->typeFlags |= foundNewData;
    pBsMetaData->expiryCount = 0;
  }

  return (errorStatus);
}

PCMDMX_ERROR pcmDmx_ReadDvbAncData(HANDLE_PCM_DOWNMIX self, UCHAR *pAncDataBuf,
                                   UINT ancDataBytes, int isMpeg2) {
  PCMDMX_ERROR errorStatus = PCMDMX_OK;
  FDK_BITSTREAM bs;
  HANDLE_FDK_BITSTREAM hBs = &bs;

  if (self == NULL) {
    return (PCMDMX_INVALID_HANDLE);
  }

  if ((pAncDataBuf == NULL) || (ancDataBytes == 0)) {
    return (PCMDMX_CORRUPT_ANC_DATA);
  }

  FDKinitBitStream(hBs, pAncDataBuf, MAX_DSE_ANC_BYTES, ancDataBytes * 8,
                   BS_READER);

  errorStatus = pcmDmx_Parse(self, hBs, ancDataBytes * 8, isMpeg2);

  return (errorStatus);
}

PCMDMX_ERROR pcmDmx_SetMatrixMixdownFromPce(HANDLE_PCM_DOWNMIX self,
                                            int matrixMixdownPresent,
                                            int matrixMixdownIdx,
                                            int pseudoSurroundEnable) {
  if (self == NULL) {
    return (PCMDMX_INVALID_HANDLE);
  }

  {
    DMX_BS_META_DATA *pBsMetaData = &self->bsMetaData[0];

    if (matrixMixdownPresent) {
      pBsMetaData->pseudoSurround = (pseudoSurroundEnable) ? 1 : 0;
      pBsMetaData->matrixMixdownIdx = matrixMixdownIdx & 0x03;
      pBsMetaData->typeFlags |= TYPE_PCE_DATA;
      pBsMetaData->expiryCount = 0;
    }
  }

  return (PCMDMX_OK);
}

PCMDMX_ERROR pcmDmx_ApplyFrame(HANDLE_PCM_DOWNMIX self, DMX_PCM *pPcmBuf,
                               const int pcmBufSize, UINT frameSize,
                               INT *nChannels, INT fInterleaved,
                               AUDIO_CHANNEL_TYPE channelType[],
                               UCHAR channelIndices[],
                               const FDK_channelMapDescr *const mapDescr,
                               INT *pDmxOutScale) {
  PCM_DMX_USER_PARAMS *pParam = NULL;
  PCMDMX_ERROR errorStatus = PCMDMX_OK;
  DUAL_CHANNEL_MODE dualChannelMode;
  PCM_DMX_CHANNEL_MODE inChMode;
  PCM_DMX_CHANNEL_MODE outChMode;
  INT devNull;
  int numOutChannels, numInChannels;
  int inStride, outStride, offset;
  int dmxMaxScale, dmxScale;
  int slot;
  UCHAR inOffsetTable[(8)];

  DMX_BS_META_DATA bsMetaData;

  if ((self == NULL) || (nChannels == NULL) || (channelType == NULL) ||
      (channelIndices == NULL) || (!FDK_chMapDescr_isValid(mapDescr))) {
    return (PCMDMX_INVALID_HANDLE);
  }

  dmxScale = 0;
  if (pDmxOutScale != NULL) {
    *pDmxOutScale = 0;
    dmxMaxScale = (3);
  } else {
    pDmxOutScale = &devNull;
    dmxMaxScale = 0;
  }

  pParam = &self->userParams;
  numInChannels = *nChannels;

  if (pPcmBuf == NULL) {
    return (PCMDMX_INVALID_ARGUMENT);
  }
  if (frameSize == 0) {
    return (PCMDMX_INVALID_ARGUMENT);
  }
  if (numInChannels == 0) {
    return (PCMDMX_INVALID_ARGUMENT);
  }
  if (numInChannels > (8)) {
    return (PCMDMX_INVALID_CH_CONFIG);
  }

  FDK_ASSERT((pParam->numOutChannelsMax <= 0) ||
             (pParam->numOutChannelsMax >= pParam->numOutChannelsMin));

  if ((self->applyProcessing == 0) &&
      ((pParam->numOutChannelsMax <= 0) ||
       (pParam->numOutChannelsMax >= numInChannels)) &&
      (pParam->numOutChannelsMin <= numInChannels)) {
    return (errorStatus);
  }

  if ((pParam->numOutChannelsMax > 0) &&
      (numInChannels > pParam->numOutChannelsMax)) {
    numOutChannels = pParam->numOutChannelsMax;
  } else if (numInChannels < pParam->numOutChannelsMin) {
    numOutChannels = pParam->numOutChannelsMin;
  } else {
    numOutChannels = numInChannels;
  }

  if ((UINT)pcmBufSize < (UINT)numOutChannels * frameSize) {
    return (PCMDMX_OUTPUT_BUFFER_TOO_SMALL);
  }

  dualChannelMode = pParam->dualChannelMode;

  errorStatus = getChannelMode(numInChannels, channelType, channelIndices,
                               inOffsetTable, &inChMode);
  if (PCMDMX_IS_FATAL_ERROR(errorStatus) || (inChMode == CH_MODE_UNDEFINED)) {
    return (PCMDMX_INVALID_CH_CONFIG);
  }

  if (fInterleaved) {
    inStride = numInChannels;
    offset = 1;
  } else {
    inStride = 1;
    offset = frameSize;
  }

  if ((pParam->expiryFrame > 0) &&
      (++self->bsMetaData[0].expiryCount >
       pParam
           ->expiryFrame)) {
#ifdef FDK_ASSERT_ENABLE
    PCMDMX_ERROR err = pcmDmx_Reset(self, PCMDMX_RESET_BS_DATA);
    FDK_ASSERT(err == PCMDMX_OK);
#else
    pcmDmx_Reset(self, PCMDMX_RESET_BS_DATA);
#endif
  }
  FDKmemcpy(&bsMetaData, &self->bsMetaData[pParam->frameDelay],
            sizeof(DMX_BS_META_DATA));
  for (slot = pParam->frameDelay; slot > 0; slot -= 1) {
    FDKmemcpy(&self->bsMetaData[slot], &self->bsMetaData[slot - 1],
              sizeof(DMX_BS_META_DATA));
  }

  if (numInChannels > numOutChannels) {
    DMX_PCM *pInPcm[(8)] = {NULL};
    DMX_PCM *pOutPcm[(8)] = {NULL};
    FIXP_DMX mixFactors[(8)][(8)];
    UCHAR outOffsetTable[(8)];
    UINT sample;
    int chCfg = 0;
    int bypScale = 0;

    if (numInChannels > SIX_CHANNEL) {
      AUDIO_CHANNEL_TYPE multiPurposeChType[2];

      multiPurposeChType[0] =
          channelType[inOffsetTable[LEFT_MULTIPRPS_CHANNEL]];
      multiPurposeChType[1] =
          channelType[inOffsetTable[RIGHT_MULTIPRPS_CHANNEL]];

      switch (inChMode) {
        case CH_MODE_5_0_2_1:
          if ((multiPurposeChType[0] == ACT_FRONT_TOP) &&
              (multiPurposeChType[1] == ACT_FRONT_TOP)) {
            chCfg = 14;
          } else {
            chCfg = 7;
          }
          break;
        case CH_MODE_3_0_3_1:
          chCfg = 11;
          break;
        case CH_MODE_3_0_4_1:
          chCfg = 12;
          break;
        default:
          chCfg = 0;
          break;
      }
    }

    outStride = (fInterleaved) ? numOutChannels : 1;
    outChMode = outChModeTable[numOutChannels];
    FDK_ASSERT(outChMode != CH_MODE_UNDEFINED);

    getChannelDescription(outChMode, mapDescr, channelType, channelIndices,
                          outOffsetTable);
    errorStatus =
        getMixFactors((chCfg > 0) ? 1 : 0,
                      (chCfg > 0) ? (PCM_DMX_CHANNEL_MODE)chCfg : inChMode,
                      outChMode, pParam, &bsMetaData, mixFactors, &dmxScale);

    bypScale = fMin(dmxMaxScale, dmxScale);
    *pDmxOutScale += bypScale;
    dmxScale -= bypScale;

    {
      int inCh, outCh, map[(8)];
      int ch = 0;
      for (inCh = 0; inCh < (8); inCh += 1) {
        if (inOffsetTable[inCh] < (UCHAR)numInChannels) {
          pInPcm[ch] = &pPcmBuf[inOffsetTable[inCh] * offset];
          map[ch++] = inCh;
        }
      }
      for (; ch < (8); ch += 1) {
        map[ch] = ch;
      }

      for (inCh = 0; inCh < numInChannels; inCh += 1) {
        if (inCh != map[inCh]) {
          for (outCh = 0; outCh < (8); outCh += 1) {
            mixFactors[outCh][inCh] = mixFactors[outCh][map[inCh]];
          }
        }
      }

      ch = 0;
      for (outCh = 0; outCh < (8); outCh += 1) {
        if (outOffsetTable[outCh] < (UCHAR)numOutChannels) {
          pOutPcm[ch] = &pPcmBuf[outOffsetTable[outCh] * offset];
          map[ch++] = outCh;
        }
      }
      for (; ch < (8); ch += 1) {
        map[ch] = ch;
      }

      for (outCh = 0; outCh < numOutChannels; outCh += 1) {
        if (outCh != map[outCh]) {
          FDKmemcpy(&mixFactors[outCh], &mixFactors[map[outCh]],
                    (8) * sizeof(FIXP_DMX));
        }
      }
    }

    for (sample = 0; sample < frameSize; sample++) {
      DMX_PCM tIn[(8)] = {0};
      FIXP_DBL tOut[(8)] = {(FIXP_DBL)0};
      int inCh, outCh;

      for (inCh = 0; inCh < numInChannels; inCh += 1) {
        if (pInPcm[inCh] != NULL) {
          tIn[inCh] = *pInPcm[inCh];
          pInPcm[inCh] += inStride;
        } else {
          tIn[inCh] = (DMX_PCM)0;
        }
      }
      for (outCh = 0; outCh < numOutChannels; outCh += 1) {
        for (inCh = 0; inCh < numInChannels; inCh += 1) {
          tOut[outCh] += fMult((DMX_PCMF)tIn[inCh], mixFactors[outCh][inCh]);
        }
        FDK_ASSERT(pOutPcm[outCh] >= pPcmBuf);
        FDK_ASSERT(pOutPcm[outCh] < &pPcmBuf[pcmBufSize]);
        *pOutPcm[outCh] = (DMX_PCM)SATURATE_SHIFT(
            tOut[outCh], DFRACT_BITS - DMX_PCM_BITS - dmxScale, DMX_PCM_BITS);
        pOutPcm[outCh] += outStride;
      }
    }

    *nChannels = numOutChannels;

  }
  else if (numInChannels < numOutChannels) {
    UCHAR outOffsetTable[(8)];

    if (numInChannels == ONE_CHANNEL) {
      DMX_PCM *pInPcm[(8)];
      DMX_PCM *pOutLF, *pOutRF;
      UINT sample;

      outStride = (fInterleaved) ? TWO_CHANNEL : 1;
      outChMode = outChModeTable[TWO_CHANNEL];

      getChannelDescription(outChMode, mapDescr, channelType, channelIndices,
                            outOffsetTable);
      pInPcm[CENTER_FRONT_CHANNEL] =
          &pPcmBuf[(frameSize - 1) *
                   inStride];

      pOutLF = &pPcmBuf[outOffsetTable[LEFT_FRONT_CHANNEL] * offset +
                        (frameSize - 1) * outStride];
      pOutRF = &pPcmBuf[outOffsetTable[RIGHT_FRONT_CHANNEL] * offset +
                        (frameSize - 1) * outStride];

      for (sample = 0; sample < frameSize; sample++) {
        *pOutLF = *pOutRF = *pInPcm[CENTER_FRONT_CHANNEL];

        pInPcm[CENTER_FRONT_CHANNEL] -= inStride;
        pOutLF -= outStride;
        pOutRF -= outStride;
      }

      inStride = outStride;
      inChMode = outChMode;
      FDKmemcpy(inOffsetTable, outOffsetTable, (8) * sizeof(UCHAR));
    }

    if (numOutChannels > TWO_CHANNEL) {
      DMX_PCM *pIn[(8)] = {NULL};
      DMX_PCM *pOut[(8)] = {NULL};
      UINT sample;
      AUDIO_CHANNEL_TYPE inChTypes[(8)];
      UCHAR inChIndices[(8)];
      UCHAR numChPerGrp[2][(4)];
      int nContentCh = 0;
      int nEmptyCh = 0;
      int ch, chGrp, isCompatible = 1;

      FDKmemcpy(inChTypes, channelType,
                numInChannels * sizeof(AUDIO_CHANNEL_TYPE));
      FDKmemclear(inChTypes + numInChannels,
                  ((8) - numInChannels) * sizeof(AUDIO_CHANNEL_TYPE));
      FDKmemcpy(inChIndices, channelIndices, numInChannels * sizeof(UCHAR));
      FDKmemclear(inChIndices + numInChannels,
                  ((8) - numInChannels) * sizeof(UCHAR));

      outStride = (fInterleaved) ? numOutChannels : 1;
      outChMode = outChModeTable[numOutChannels];
      FDK_ASSERT(outChMode != CH_MODE_UNDEFINED);

      for (chGrp = 0; chGrp < (4); chGrp += 1) {
        numChPerGrp[0][chGrp] = (inChMode >> (chGrp * 4)) & 0xF;
        numChPerGrp[1][chGrp] = (outChMode >> (chGrp * 4)) & 0xF;

        if (numChPerGrp[0][chGrp] > numChPerGrp[1][chGrp]) {
          isCompatible = 0;
          break;
        }
      }

      if (isCompatible) {
        getChannelDescription(outChMode, mapDescr, channelType, channelIndices,
                              outOffsetTable);
        if ((numChPerGrp[0][CH_GROUP_REAR] % 2) &&
            !(numChPerGrp[1][CH_GROUP_REAR] % 2)) {
          if (numChPerGrp[0][CH_GROUP_REAR] == 1) {
            inOffsetTable[RIGHT_REAR_CHANNEL] =
                inOffsetTable[LEFT_REAR_CHANNEL];
          } else if (numChPerGrp[0][CH_GROUP_REAR] == 3) {
            inOffsetTable[RIGHT_MULTIPRPS_CHANNEL] =
                inOffsetTable[LEFT_MULTIPRPS_CHANNEL];
          }
        }
      } else {
        FDKmemcpy(outOffsetTable, inOffsetTable, (8) * sizeof(UCHAR));
      }

      for (ch = 0; ch < (8); ch += 1) {
        if ((outOffsetTable[ch] < 255) &&
            (inOffsetTable[ch] < 255)) {
          pIn[nContentCh] =
              &pPcmBuf[inOffsetTable[ch] * offset + (frameSize - 1) * inStride];
          pOut[nContentCh] = &pPcmBuf[outOffsetTable[ch] * offset +
                                      (frameSize - 1) * outStride];
          channelType[outOffsetTable[ch]] = inChTypes[inOffsetTable[ch]];
          channelIndices[outOffsetTable[ch]] = inChIndices[inOffsetTable[ch]];
          inOffsetTable[ch] = 255;
          outOffsetTable[ch] = 255;
          nContentCh += 1;
        }
      }
      if (isCompatible) {
        for (ch = 0; ch < (8); ch += 1) {
          if (inOffsetTable[ch] < 255) {
            int outCh;
            for (outCh = 0; outCh < (8); outCh += 1) {
              if (outOffsetTable[outCh] < 255) {
                break;
              }
            }
            if (outCh >= (8)) {
              FDK_ASSERT(0);
              break;
            }
            pIn[nContentCh] = &pPcmBuf[inOffsetTable[ch] * offset +
                                       (frameSize - 1) * inStride];
            pOut[nContentCh] = &pPcmBuf[outOffsetTable[outCh] * offset +
                                        (frameSize - 1) * outStride];
            FDK_ASSERT(inOffsetTable[outCh] < numInChannels);
            FDK_ASSERT(outOffsetTable[outCh] < numOutChannels);
            channelType[outOffsetTable[outCh]] = inChTypes[inOffsetTable[ch]];
            channelIndices[outOffsetTable[outCh]] =
                inChIndices[inOffsetTable[ch]];
            inOffsetTable[ch] = 255;
            outOffsetTable[outCh] = 255;
            nContentCh += 1;
          }
        }
        for (ch = 0; ch < (8); ch += 1) {
          if (outOffsetTable[ch] < 255) {
            pOut[nContentCh + nEmptyCh] = &pPcmBuf[outOffsetTable[ch] * offset +
                                                   (frameSize - 1) * outStride];
            channelType[outOffsetTable[ch]] = ACT_NONE;
            channelIndices[outOffsetTable[ch]] = (UCHAR)nEmptyCh;
            outOffsetTable[ch] = 255;
            nEmptyCh += 1;
          }
        }
      } else {
        for (ch = nContentCh; ch < numOutChannels; ch += 1) {
          pOut[ch] = &pPcmBuf[ch * offset + (frameSize - 1) * outStride];
          channelType[ch] = ACT_NONE;
          channelIndices[ch] = (UCHAR)nEmptyCh;
          nEmptyCh += 1;
        }
      }

      for (sample = 0; sample < frameSize; sample += 1) {
        DMX_PCM tIn[(8)];
        for (ch = 0; ch < nContentCh; ch += 1) {
          tIn[ch] = *pIn[ch];
          pIn[ch] -= inStride;
        }
        for (ch = 0; ch < nContentCh; ch += 1) {
          *pOut[ch] = tIn[ch];
          pOut[ch] -= outStride;
        }
      }

      for (sample = 0; sample < frameSize; sample++) {
        for (ch = nContentCh; ch < numOutChannels; ch += 1) {
          *pOut[ch] = (DMX_PCM)0;
          pOut[ch] -= outStride;
        }
      }
    }

    *nChannels = numOutChannels;
  }
  else if (numInChannels == numOutChannels) {
    switch (numInChannels) {
      case 2: {
        DMX_PCM *pInPcm[(8)];
        DMX_PCM *pOutL, *pOutR;
        FIXP_DMX flev;

        UINT sample;

        if (fInterleaved) {
          inStride = numInChannels;
          outStride =
              2;
          offset = 1;
        } else {
          inStride = 1;
          outStride = 1;
          offset = frameSize;
        }

        pInPcm[LEFT_FRONT_CHANNEL] =
            &pPcmBuf[inOffsetTable[LEFT_FRONT_CHANNEL] * offset];
        pInPcm[RIGHT_FRONT_CHANNEL] =
            &pPcmBuf[inOffsetTable[RIGHT_FRONT_CHANNEL] * offset];

        pOutL = pInPcm[LEFT_FRONT_CHANNEL];
        pOutR = pInPcm[RIGHT_FRONT_CHANNEL];

        flev = FL2FXCONST_DMX(0.70710678f);
        switch (dualChannelMode) {
          case CH1_MODE:
            for (sample = 0; sample < frameSize; sample++) {
              *pOutL = *pOutR = (DMX_PCM)SATURATE_RIGHT_SHIFT(
                  fMult((DMX_PCMF)*pInPcm[LEFT_FRONT_CHANNEL], flev),
                  DFRACT_BITS - DMX_PCM_BITS, DMX_PCM_BITS);

              pInPcm[LEFT_FRONT_CHANNEL] += inStride;
              pOutL += outStride;
              pOutR += outStride;
            }
            break;
          case CH2_MODE:
            for (sample = 0; sample < frameSize; sample++) {
              *pOutL = *pOutR = (DMX_PCM)SATURATE_RIGHT_SHIFT(
                  fMult((DMX_PCMF)*pInPcm[RIGHT_FRONT_CHANNEL], flev),
                  DFRACT_BITS - DMX_PCM_BITS, DMX_PCM_BITS);

              pInPcm[RIGHT_FRONT_CHANNEL] += inStride;
              pOutL += outStride;
              pOutR += outStride;
            }
            break;
          case MIXED_MODE:
            for (sample = 0; sample < frameSize; sample++) {
              *pOutL = *pOutR = (*pInPcm[LEFT_FRONT_CHANNEL] >> 1) +
                                (*pInPcm[RIGHT_FRONT_CHANNEL] >> 1);

              pInPcm[LEFT_FRONT_CHANNEL] += inStride;
              pInPcm[RIGHT_FRONT_CHANNEL] += inStride;
              pOutL += outStride;
              pOutR += outStride;
            }
            break;
          default:
          case STEREO_MODE:
            break;
        }
      } break;

      default:
        break;
    }
  }

  return (errorStatus);
}

PCMDMX_ERROR pcmDmx_Close(HANDLE_PCM_DOWNMIX *pSelf) {
  if (pSelf == NULL) {
    return (PCMDMX_INVALID_HANDLE);
  }

  FreePcmDmxInstance(pSelf);
  *pSelf = NULL;

  return (PCMDMX_OK);
}

PCMDMX_ERROR pcmDmx_GetLibInfo(LIB_INFO *info) {
  int i;

  if (info == NULL) {
    return PCMDMX_INVALID_ARGUMENT;
  }

  for (i = 0; i < FDK_MODULE_LAST; i++) {
    if (info[i].module_id == FDK_NONE) break;
  }
  if (i == FDK_MODULE_LAST) {
    return PCMDMX_INVALID_ARGUMENT;
  }

  info[i].module_id = FDK_PCMDMX;
  info[i].version =
      LIB_VERSION(PCMUTIL_LIB_VL0, PCMUTIL_LIB_VL1, PCMUTIL_LIB_VL2);
  LIB_VERSION_STRING(info + i);
  info[i].build_date = PCMUTIL_LIB_BUILD_DATE;
  info[i].build_time = PCMUTIL_LIB_BUILD_TIME;
  info[i].title = PCMDMX_LIB_TITLE;

  info[i].flags = 0 | CAPF_DMX_BLIND
                  | CAPF_DMX_PCE
                  | CAPF_DMX_ARIB
                  | CAPF_DMX_DVB
                  | CAPF_DMX_CH_EXP
                  | CAPF_DMX_6_CH | CAPF_DMX_8_CH;

  FDK_toolsGetLibInfo(info);

  return PCMDMX_OK;
}
