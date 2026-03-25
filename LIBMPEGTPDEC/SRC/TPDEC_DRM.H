/******************* MPEG transport format decoder library *********************

   Author(s):   Josef Hoepfl

   Description: DRM interface

*******************************************************************************/

#ifndef TPDEC_DRM_H
#define TPDEC_DRM_H

#include "tpdec_lib.h"

#include "FDK_crc.h"

typedef struct {
  FDK_CRCINFO crcInfo;
  USHORT crcReadValue;

} STRUCT_DRM;

typedef STRUCT_DRM *HANDLE_DRM;

void drmRead_CrcInit(HANDLE_DRM pDrm);

int drmRead_CrcStartReg(HANDLE_DRM pDrm, HANDLE_FDK_BITSTREAM hBs, int mBits);

void drmRead_CrcEndReg(HANDLE_DRM pDrm, HANDLE_FDK_BITSTREAM hBs, int reg);

TRANSPORTDEC_ERROR drmRead_CrcCheck(HANDLE_DRM pDrm);

TRANSPORTDEC_ERROR drmRead_DecodeHeader(HANDLE_DRM pDrm,
                                        HANDLE_FDK_BITSTREAM bs);

TRANSPORTDEC_ERROR DrmRawSdcAudioConfig_Parse(CSAudioSpecificConfig *pAsc,
                                              HANDLE_FDK_BITSTREAM hBs,
                                              CSTpCallBacks *cb,
                                              const UCHAR configMode,
                                              const UCHAR configChanged);

#endif /* TPDEC_DRM_H */
