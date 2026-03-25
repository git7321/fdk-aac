/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description:

*******************************************************************************/

#ifndef AAC_RAM_H
#define AAC_RAM_H
#include "common_fix.h"
#include "aacdecoder.h"
#include "channel.h"
#include "ac_arith_coder.h"
#include "aacdec_hcr_types.h"
#include "aacdec_hcr.h"

#define MAX_SYNCHS 10
#define SAMPL_FREQS 12

H_ALLOC_MEM(AacDecoder, AAC_DECODER_INSTANCE)

H_ALLOC_MEM(DrcInfo, CDrcInfo)

H_ALLOC_MEM(AacDecoderStaticChannelInfo, CAacDecoderStaticChannelInfo)
H_ALLOC_MEM(AacDecoderChannelInfo, CAacDecoderChannelInfo)
H_ALLOC_MEM(OverlapBuffer, FIXP_DBL)

H_ALLOC_MEM(CpePersistentData, CpePersistentData)
H_ALLOC_MEM(CplxPredictionData, CCplxPredictionData)
H_ALLOC_MEM(SpectralCoeffs, FIXP_DBL)
H_ALLOC_MEM(SpecScale, SHORT)

H_ALLOC_MEM(TimeDataFlush, PCM_DEC)

H_ALLOC_MEM_OVERLAY(WorkBufferCore1, CWorkBufferCore1)
H_ALLOC_MEM_OVERLAY(WorkBufferCore2, FIXP_DBL)
H_ALLOC_MEM_OVERLAY(WorkBufferCore5, PCM_DEC)
H_ALLOC_MEM_OVERLAY(WorkBufferCore6, SCHAR)

#endif /* #ifndef AAC_RAM_H */
