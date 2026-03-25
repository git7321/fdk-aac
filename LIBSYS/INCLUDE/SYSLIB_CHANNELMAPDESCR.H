/************************* System integration library **************************

   Author(s):   Thomas Dietzen

   Description:

*******************************************************************************/

#ifndef SYSLIB_CHANNELMAPDESCR_H
#define SYSLIB_CHANNELMAPDESCR_H

#include "machine_type.h"

typedef struct {
  const UCHAR*
      pChannelMap;
  UCHAR numChannels;
} CHANNEL_MAP_INFO;

typedef struct {
  const CHANNEL_MAP_INFO* pMapInfoTab;
  UINT mapInfoTabLen;
  UINT fPassThrough;
} FDK_channelMapDescr;

#ifdef __cplusplus
extern "C" {
#endif

void FDK_chMapDescr_init(FDK_channelMapDescr* const pMapDescr,
                         const CHANNEL_MAP_INFO* const pMapInfoTab,
                         const UINT mapInfoTabLen, const UINT fPassThrough);

int FDK_chMapDescr_setPassThrough(FDK_channelMapDescr* const pMapDescr,
                                  UINT fPassThrough);

UCHAR FDK_chMapDescr_getMapValue(const FDK_channelMapDescr* const pMapDescr,
                                 const UCHAR chIdx, const UINT mapIdx);

int FDK_chMapDescr_isValid(const FDK_channelMapDescr* const pMapDescr);

extern const CHANNEL_MAP_INFO FDK_mapInfoTabWg4[];
extern const UINT FDK_mapInfoTabLenWg4;

#ifdef __cplusplus
}
#endif

#endif /* !defined(SYSLIB_CHANNELMAPDESCR_H) */
