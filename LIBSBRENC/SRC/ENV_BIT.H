/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef ENV_BIT_H
#define ENV_BIT_H

#include "sbr_encoder.h"
#include "FDK_crc.h"

#define SBR_CRC_POLY (0x0233)
#define SBR_CRC_MASK (0x0200)
#define SBR_CRC_RANGE (0x03FF)
#define SBR_CRC_MAXREGS 1
#define SBR_CRCINIT (0x0)

#define SI_SBR_CRC_ENABLE_BITS 0
#define SI_SBR_CRC_BITS 10
#define SI_SBR_DRM_CRC_BITS 8

struct COMMON_DATA;

INT FDKsbrEnc_InitSbrBitstream(struct COMMON_DATA *hCmonData, UCHAR *memoryBase,
                               INT memorySize, HANDLE_FDK_CRCINFO hCrcInfo,
                               UINT sbrSyntaxFlags);

void FDKsbrEnc_AssembleSbrBitstream(struct COMMON_DATA *hCmonData,
                                    HANDLE_FDK_CRCINFO hCrcInfo, INT crcRegion,
                                    UINT sbrSyntaxFlags);

#endif /* #ifndef ENV_BIT_H */
