/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description:

*******************************************************************************/

#include "aac_ram.h"
#include "aac_rom.h"

#define WORKBUFFER1_TAG 0
#define WORKBUFFER2_TAG 1
#define WORKBUFFER5_TAG 6
#define WORKBUFFER6_TAG 7

C_ALLOC_MEM(AacDecoder, struct AAC_DECODER_INSTANCE, 1)
C_ALLOC_MEM2(AacDecoderStaticChannelInfo, CAacDecoderStaticChannelInfo, 1, (8))
C_AALLOC_MEM2(AacDecoderChannelInfo, CAacDecoderChannelInfo, 1, (8))
C_AALLOC_MEM2(OverlapBuffer, FIXP_DBL, OverlapBufferSize, (8))

C_ALLOC_MEM(DrcInfo, CDrcInfo, 1)
C_ALLOC_MEM(CpePersistentData, CpePersistentData, 1)
C_ALLOC_MEM(CplxPredictionData, CCplxPredictionData, 1)
C_ALLOC_MEM2(TimeDataFlush, PCM_DEC, TIME_DATA_FLUSH_SIZE, (8))

C_ALLOC_MEM_OVERLAY(WorkBufferCore2, FIXP_DBL, ((8) * 1024), SECT_DATA_L2,
                    WORKBUFFER2_TAG)
C_ALLOC_MEM_OVERLAY(WorkBufferCore6, SCHAR,
                    fMax((INT)(sizeof(FIXP_DBL) * WB_SECTION_SIZE),
                         (INT)sizeof(CAacDecoderCommonData)),
                    SECT_DATA_L2, WORKBUFFER6_TAG)
C_ALLOC_MEM_OVERLAY(WorkBufferCore1, CWorkBufferCore1, 1, SECT_DATA_L1,
                    WORKBUFFER1_TAG)
C_ALLOC_MEM_OVERLAY(WorkBufferCore5, PCM_DEC, (8) * (1024 * 4) * 2,
                    SECT_DATA_EXTERN, WORKBUFFER5_TAG)
