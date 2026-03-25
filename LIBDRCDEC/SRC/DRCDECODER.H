/************************* MPEG-D DRC decoder library **************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef DRCDECODER_H
#define DRCDECODER_H

#define MAX_ACTIVE_DRCS 3

typedef enum { DM_REGULAR_DELAY = 0, DM_LOW_DELAY = 1 } DELAY_MODE;

typedef enum {
  DE_OK = 0,
  DE_NOT_OK = -100,
  DE_PARAM_OUT_OF_RANGE,
  DE_PARAM_INVALID,
  DE_MEMORY_ERROR
} DRC_ERROR;

typedef enum { SDM_OFF, SDM_QMF64, SDM_QMF71, SDM_STFT256 } SUBBAND_DOMAIN_MODE;

#define DOWNMIX_ID_BASE_LAYOUT 0x0
#define DOWNMIX_ID_ANY_DOWNMIX 0x7F
#define DRC_SET_ID_NO_DRC 0x0
#define DRC_SET_ID_ANY_DRC 0x3F

#define LOCATION_MP4_INSTREAM_UNIDRC 0x1
#define LOCATION_MP4_DYN_RANGE_INFO 0x2
#define LOCATION_MP4_COMPRESSION_VALUE 0x3
#define LOCATION_SELECTED \
  LOCATION_MP4_INSTREAM_UNIDRC

#define MAX_REQUESTS_DOWNMIX_ID 15
#define MAX_REQUESTS_DRC_FEATURE 7
#define MAX_REQUESTS_DRC_EFFECT_TYPE 15

#define DEFAULT_LOUDNESS_DEVIATION_MAX 63

#define DRC_INPUT_LOUDNESS_TARGET FL2FXCONST_DBL(-31.0f / (float)(1 << 7))
#define DRC_INPUT_LOUDNESS_TARGET_SGL FL2FXCONST_SGL(-31.0f / (float)(1 << 7))

#endif
