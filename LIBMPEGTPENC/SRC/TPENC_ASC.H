/******************* MPEG transport format encoder library *********************

   Author(s):   Manuel Jander

   Description: Audio Specific Config writer

*******************************************************************************/

#ifndef TPENC_ASC_H
#define TPENC_ASC_H

int getChannelConfig(const CHANNEL_MODE channel_mode,
                     const UCHAR channel_config_zero);

int transportEnc_writePCE(HANDLE_FDK_BITSTREAM hBs, CHANNEL_MODE channelMode,
                          INT sampleRate, int instanceTagPCE, int profile,
                          int matrixMixdownA, int pseudoSurroundEnable,
                          UINT alignAnchor);

int transportEnc_GetPCEBits(CHANNEL_MODE channelMode, int matrixMixdownA,
                            int bits);

#endif /* TPENC_ASC_H */
