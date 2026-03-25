#include "syslib_channelMapDescr.h"

#define DFLT_CH_MAP_TAB_LEN \
  (15)

static const UCHAR mapFallback[] = {0,  1,  2,  3,  4,  5,  6,  7,
                                    8,  9,  10, 11, 12, 13, 14, 15,
                                    16, 17, 18, 19, 20, 21, 22, 23};
static const UCHAR mapCfg1[] = {0, 1};
static const UCHAR mapCfg2[] = {0, 1};
static const UCHAR mapCfg3[] = {2, 0, 1};
static const UCHAR mapCfg4[] = {2, 0, 1, 3};
static const UCHAR mapCfg5[] = {2, 0, 1, 3, 4};
static const UCHAR mapCfg6[] = {2, 0, 1, 4, 5, 3};
static const UCHAR mapCfg7[] = {2, 6, 7, 0, 1, 4, 5, 3};
static const UCHAR mapCfg11[] = {2, 0, 1, 4, 5, 6, 3};
static const UCHAR mapCfg12[] = {2, 0, 1, 6, 7, 4, 5, 3};
static const UCHAR mapCfg13[] = {2,  6,  7,  0,  1,  10, 11, 4,
                                 5,  8,  3,  9,  14, 12, 13, 18,
                                 19, 15, 16, 17, 20, 21, 22, 23};
static const UCHAR mapCfg14[] = {2, 0, 1, 4, 5, 3, 6, 7};

static const CHANNEL_MAP_INFO mapInfoTabDflt[DFLT_CH_MAP_TAB_LEN] =
    {
     /*  0 */ {mapFallback, 24},
     /*  1 */ {mapCfg1, 2},
     /*  2 */ {mapCfg2, 2},
     /*  3 */ {mapCfg3, 3},
     /*  4 */ {mapCfg4, 4},
     /*  5 */ {mapCfg5, 5},
     /*  6 */ {mapCfg6, 6},
     /*  7 */ {mapCfg7, 8},
     /*  8 */ {mapFallback, 24},
     /*  9 */ {mapFallback, 24},
     /* 10 */ {mapFallback, 24},
     /* 11 */ {mapCfg11, 7},
     /* 12 */ {mapCfg12, 8},
     /* 13 */ {mapCfg13, 24},
     /* 14 */ {mapCfg14, 8}};

static const UCHAR mapWg4Cfg1[]  = {0, 1};
static const UCHAR mapWg4Cfg2[]  = {0, 1};
static const UCHAR mapWg4Cfg3[]  = {2, 0, 1};
static const UCHAR mapWg4Cfg4[]  = {3, 0, 1, 2};
static const UCHAR mapWg4Cfg5[]  = {4, 0, 1, 2, 3};
static const UCHAR mapWg4Cfg6[]  = {4, 0, 1, 2, 3, 5};
static const UCHAR mapWg4Cfg7[]  = {6, 0, 1, 2, 3, 4, 5, 7};
static const UCHAR mapWg4Cfg14[] = {6, 0, 1, 2, 3, 4, 5, 7};

const CHANNEL_MAP_INFO FDK_mapInfoTabWg4[] =
    {
     /*  0 */ {mapFallback, 24},
     /*  1 */ {mapWg4Cfg1, 2},
     /*  2 */ {mapWg4Cfg2, 2},
     /*  3 */ {mapWg4Cfg3, 3},
     /*  4 */ {mapWg4Cfg4, 4},
     /*  5 */ {mapWg4Cfg5, 5},
     /*  6 */ {mapWg4Cfg6, 6},
     /*  7 */ {mapWg4Cfg7, 8},
     /*  8 */ {mapFallback, 24},
     /*  9 */ {mapFallback, 24},
     /* 10 */ {mapFallback, 24},
     /* 11 */ {mapFallback, 24},
     /* 12 */ {mapFallback, 24},
     /* 13 */ {mapFallback, 24},
     /* 14 */ {mapFallback, 24}};

const UINT FDK_mapInfoTabLenWg4 = sizeof(FDK_mapInfoTabWg4)/sizeof(FDK_mapInfoTabWg4[0]);

UCHAR FDK_chMapDescr_getMapValue(const FDK_channelMapDescr* const pMapDescr,
                                 const UCHAR chIdx, const UINT mapIdx) {
  UCHAR mapValue = chIdx;

  FDK_ASSERT(pMapDescr != NULL);

  if ((pMapDescr->fPassThrough == 0) && (pMapDescr->pMapInfoTab != NULL) &&
      (pMapDescr->mapInfoTabLen > mapIdx)) {
    if (chIdx < pMapDescr->pMapInfoTab[mapIdx].numChannels) {
      mapValue = pMapDescr->pMapInfoTab[mapIdx].pChannelMap[chIdx];
    }
  }
  return mapValue;
}

static int fdk_chMapDescr_isValidMap(const CHANNEL_MAP_INFO* const pMapInfo) {
  int result = 1;
  UINT i;

  if (pMapInfo == NULL) {
    result = 0;
  } else {
    UINT numChannels = pMapInfo->numChannels;

    if (numChannels < 32) {
      UINT mappedChMask = 0x0;
      for (i = 0; i < numChannels; i += 1) {
        mappedChMask |= 1 << pMapInfo->pChannelMap[i];
      }
      if (mappedChMask != (((UINT)1 << numChannels) - 1)) {
        result = 0;
      }
    } else {
      for (i = 0; (i < numChannels) && result; i += 1) {
        UINT j;
        UCHAR value0 = pMapInfo->pChannelMap[i];

        if (value0 > numChannels - 1) {
          result = 0;
        }
        for (j = numChannels - 1; (j > i) && result; j -= 1) {
          if (value0 == pMapInfo->pChannelMap[j]) {
            result = 0;
          }
        }
      }
    }
  }

  return result;
}

int FDK_chMapDescr_isValid(const FDK_channelMapDescr* const pMapDescr) {
  int result = 0;
  UINT i;

  if (pMapDescr != NULL) {
    result = 1;
    for (i = 0; (i < pMapDescr->mapInfoTabLen) && result; i += 1) {
      if (!fdk_chMapDescr_isValidMap(&pMapDescr->pMapInfoTab[i])) {
        result = 0;
      }
    }
  }
  return result;
}

void FDK_chMapDescr_init(FDK_channelMapDescr* const pMapDescr,
                         const CHANNEL_MAP_INFO* const pMapInfoTab,
                         const UINT mapInfoTabLen, const UINT fPassThrough) {
  if (pMapDescr != NULL) {
    int useDefaultTab = 1;

    pMapDescr->fPassThrough = (fPassThrough == 0) ? 0 : 1;

    if ((pMapInfoTab != NULL) && (mapInfoTabLen > 0)) {
      pMapDescr->pMapInfoTab = pMapInfoTab;
      pMapDescr->mapInfoTabLen = mapInfoTabLen;
      useDefaultTab = (FDK_chMapDescr_isValid(pMapDescr) == 0) ? 1 : 0;
    }
    if (useDefaultTab != 0) {
      pMapDescr->pMapInfoTab = mapInfoTabDflt;
      pMapDescr->mapInfoTabLen = DFLT_CH_MAP_TAB_LEN;
    }
  }
}

int FDK_chMapDescr_setPassThrough(FDK_channelMapDescr* const pMapDescr,
                                  UINT fPassThrough) {
  int err = 1;

  if (pMapDescr != NULL) {
    if ((pMapDescr->pMapInfoTab != NULL) && (pMapDescr->mapInfoTabLen > 0)) {
      pMapDescr->fPassThrough = (fPassThrough == 0) ? 0 : 1;
      err = 0;
    }
  }

  return err;
}
