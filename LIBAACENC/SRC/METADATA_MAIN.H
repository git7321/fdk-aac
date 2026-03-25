/**************************** AAC encoder library ******************************

   Author(s):   V. Bacigalupo

   Description: Metadata Encoder library interface functions

*******************************************************************************/

#ifndef METADATA_MAIN_H
#define METADATA_MAIN_H

#include "aacenc_lib.h"
#include "aacenc.h"

typedef enum {
  METADATA_OK = 0x0000,
  METADATA_INVALID_HANDLE =
      0x0020,
  METADATA_MEMORY_ERROR = 0x0021,
  METADATA_INIT_ERROR = 0x0040,
  METADATA_ENCODE_ERROR =
      0x0060

} FDK_METADATA_ERROR;

typedef struct FDK_METADATA_ENCODER *HANDLE_FDK_METADATA_ENCODER;

FDK_METADATA_ERROR FDK_MetadataEnc_Open(
    HANDLE_FDK_METADATA_ENCODER *phMetadataEnc, const UINT maxChannels);

FDK_METADATA_ERROR FDK_MetadataEnc_Init(
    HANDLE_FDK_METADATA_ENCODER hMetadataEnc, const INT resetStates,
    const INT metadataMode, const INT audioDelay, const UINT frameLength,
    const UINT sampleRate, const UINT nChannels, const CHANNEL_MODE channelMode,
    const CHANNEL_ORDER channelOrder);

FDK_METADATA_ERROR FDK_MetadataEnc_Process(
    HANDLE_FDK_METADATA_ENCODER hMetadataEnc, INT_PCM *const pAudioSamples,
    const UINT audioSamplesBufSize, const INT nAudioSamples,
    const AACENC_MetaData *const pMetadata,
    AACENC_EXT_PAYLOAD **ppMetaDataExtPayload, UINT *nMetaDataExtensions,
    INT *matrix_mixdown_idx);

FDK_METADATA_ERROR FDK_MetadataEnc_Close(
    HANDLE_FDK_METADATA_ENCODER *phMetaData);

INT FDK_MetadataEnc_GetDelay(HANDLE_FDK_METADATA_ENCODER hMetadataEnc);

#endif /* METADATA_MAIN_H */
