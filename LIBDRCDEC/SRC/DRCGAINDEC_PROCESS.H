/************************* MPEG-D DRC decoder library **************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef DRCGAINDEC_PROCESS_H
#define DRCGAINDEC_PROCESS_H

DRC_ERROR
processDrcTime(HANDLE_DRC_GAIN_DECODER hGainDec, const int activeDrcIndex,
               const int delaySamples, const int channelOffset,
               const int drcChannelOffset, const int numChannelsProcessed,
               const int timeDataChannelOffset, FIXP_DBL* deinterleavedAudio);

DRC_ERROR
processDrcSubband(HANDLE_DRC_GAIN_DECODER hGainDec, const int activeDrcIndex,
                  const int delaySamples, const int channelOffset,
                  const int drcChannelOffset, const int numChannelsProcessed,
                  const int processSingleTimeslot,
                  FIXP_DBL* deinterleavedAudioReal[],
                  FIXP_DBL* deinterleavedAudioImag[]);
#endif
