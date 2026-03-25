/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Bitstream encoder

*******************************************************************************/

#ifndef BITENC_H
#define BITENC_H

#include "qc_data.h"
#include "aacenc_tns.h"
#include "channel_map.h"
#include "interface.h"
#include "FDK_audio.h"
#include "aacenc.h"

#include "tpenc_lib.h"

typedef enum {
  MAX_ENCODER_CHANNELS = 9,
  MAX_BLOCK_TYPES = 4,
  MAX_AAC_LAYERS = 9,
  MAX_LAYERS = MAX_AAC_LAYERS,
  FIRST_LAY = 1
} _MAX_CONST;

#define BUFFER_MX_HUFFCB_SIZE \
  (32 * sizeof(INT))

#define EL_ID_BITS (3)

AAC_ENCODER_ERROR FDKaacEnc_ChannelElementWrite(
    HANDLE_TRANSPORTENC hTpEnc, ELEMENT_INFO *pElInfo,
    QC_OUT_CHANNEL *qcOutChannel[(2)], PSY_OUT_ELEMENT *psyOutElement,
    PSY_OUT_CHANNEL *psyOutChannel[(2)], UINT syntaxFlags,
    AUDIO_OBJECT_TYPE aot, SCHAR epConfig, INT *pBitDemand, UCHAR minCnt);

AAC_ENCODER_ERROR FDKaacEnc_WriteBitstream(HANDLE_TRANSPORTENC hTpEnc,
                                           CHANNEL_MAPPING *channelMapping,
                                           QC_OUT *qcOut, PSY_OUT *psyOut,
                                           QC_STATE *qcKernel,
                                           AUDIO_OBJECT_TYPE aot,
                                           UINT syntaxFlags, SCHAR epConfig);

INT FDKaacEnc_writeExtensionData(HANDLE_TRANSPORTENC hTpEnc,
                                 QC_OUT_EXTENSION *pExtension,
                                 INT elInstanceTag, UINT alignAnchor,
                                 UINT syntaxFlags, AUDIO_OBJECT_TYPE aot,
                                 SCHAR epConfig);

#endif /* BITENC_H */
