/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef SBR_RAM_H
#define SBR_RAM_H

#include "sbrdecoder.h"

#include "env_extr.h"
#include "sbr_dec.h"

#define SBRDEC_MAX_CH_PER_ELEMENT (2)

#define FRAME_OK (0)
#define FRAME_ERROR (1)
#define FRAME_ERROR_ALLSLOTS (2)

typedef struct {
  SBR_CHANNEL *pSbrChannel[SBRDEC_MAX_CH_PER_ELEMENT];
  TRANSPOSER_SETTINGS
  transposerSettings; /* Common transport settings for each individual
                         channel of an element */
  HANDLE_FDK_BITSTREAM hBs;

  MP4_ELEMENT_ID
  elementID;     /* Element ID set during initialization. Can be used for
                    concealment */
  int nChannels; /* Number of elements output channels (=2 in case of PS) */

  UCHAR frameErrorFlag[(1) + 1]; /* Frame error status (for every slot in the
                                    delay line). Will be copied into header at
                                    the very beginning of decodeElement()
                                    routine. */

  UCHAR useFrameSlot; /* Index which defines which slot will be decoded/filled
                         next (used with additional delay) */
  UCHAR useHeaderSlot[(1) + 1]; /* Index array that provides the link between
                                   header and frame data (important when
                                   processing with additional delay). */
} SBR_DECODER_ELEMENT;

struct SBR_DECODER_INSTANCE {
  SBR_DECODER_ELEMENT *pSbrElement[(8)];
  SBR_HEADER_DATA sbrHeader[(
      8)][(1) + 1]; /* Sbr header for each individual channel of an element */

  HANDLE_FDK_QMF_DOMAIN pQmfDomain;

  HANDLE_PS_DEC hParametricStereoDec;

  /* Global parameters */
  AUDIO_OBJECT_TYPE coreCodec; /* AOT of core codec */
  int numSbrElements;
  int numSbrChannels;
  INT sampleRateIn;  /* SBR decoder input sampling rate; might be different than
                        the transposer input sampling rate. */
  INT sampleRateOut; /* Sampling rate of the SBR decoder output audio samples.
                      */
  USHORT codecFrameSize;
  UCHAR synDownsampleFac;
  INT downscaleFactor;
  UCHAR numDelayFrames; /* The current number of additional delay frames used
                           for processing. */
  UCHAR harmonicSBR;
  UCHAR
  numFlushedFrames; /* The variable counts the number of frames which are
                       flushed consecutively. */

  UINT flags;

  INT sbrInDataHeadroom; /* Headroom of the SBR input time signal to prevent
                            clipping */
};

H_ALLOC_MEM(Ram_SbrDecElement, SBR_DECODER_ELEMENT)
H_ALLOC_MEM(Ram_SbrDecChannel, SBR_CHANNEL)
H_ALLOC_MEM(Ram_SbrDecoder, struct SBR_DECODER_INSTANCE)

H_ALLOC_MEM(Ram_sbr_QmfStatesSynthesis, FIXP_QSS)
H_ALLOC_MEM(Ram_sbr_OverlapBuffer, FIXP_DBL)

H_ALLOC_MEM(Ram_sbr_HBEOverlapBuffer, FIXP_DBL)

H_ALLOC_MEM(Ram_ps_dec, PS_DEC)

#endif /* SBR_RAM_H */
