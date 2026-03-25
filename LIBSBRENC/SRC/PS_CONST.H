/**************************** SBR encoder library ******************************

   Author(s):   N. Rettelbach

   Description: Parametric Stereo constants

*******************************************************************************/

#ifndef PS_CONST_H
#define PS_CONST_H

#define MAX_PS_CHANNELS (2)
#define HYBRID_MAX_QMF_BANDS (3)
#define HYBRID_FILTER_LENGTH (13)
#define HYBRID_FILTER_DELAY ((HYBRID_FILTER_LENGTH - 1) / 2)

#define HYBRID_FRAMESIZE (32)
#define HYBRID_READ_OFFSET (10)

#define MAX_HYBRID_BANDS ((64 - HYBRID_MAX_QMF_BANDS + 10))

typedef enum {
  PS_RES_COARSE = 0,
  PS_RES_MID = 1,
  PS_RES_FINE = 2
} PS_RESOLUTION;

typedef enum {
  PS_BANDS_COARSE = 10,
  PS_BANDS_MID = 20,
  PS_MAX_BANDS = PS_BANDS_MID
} PS_BANDS;

typedef enum { PS_IID_RES_COARSE = 0, PS_IID_RES_FINE } PS_IID_RESOLUTION;

typedef enum { PS_ICC_ROT_A = 0, PS_ICC_ROT_B } PS_ICC_ROTATION_MODE;

typedef enum { PS_DELTA_FREQ, PS_DELTA_TIME } PS_DELTA;

typedef enum {
  PS_MAX_ENVELOPES = 4

} PS_CONSTS;

typedef enum {
  PSENC_OK = 0x0000, /*!< No error happened. All fine. */
  PSENC_INVALID_HANDLE =
      0x0020, /*!< Handle passed to function call was invalid. */
  PSENC_MEMORY_ERROR = 0x0021, /*!< Memory allocation failed. */
  PSENC_INIT_ERROR = 0x0040,   /*!< General initialization error. */
  PSENC_ENCODE_ERROR = 0x0060  /*!< The encoding process was interrupted by an
                                  unexpected error. */

} FDK_PSENC_ERROR;

#endif
