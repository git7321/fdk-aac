/**************************** AAC decoder library ******************************

   Author(s):   Matthias Hildenbrand

   Description:

*******************************************************************************/

#include "FDK_delay.h"

#include "genericStds.h"

#define MAX_FRAME_LENGTH (1024)

INT FDK_Delay_Create(FDK_SignalDelay* data, const USHORT delay,
                     const UCHAR num_channels) {
  FDK_ASSERT(data != NULL);
  FDK_ASSERT(num_channels > 0);

  if (delay > 0) {
    data->delay_line =
        (PCM_DEC*)FDKcalloc(num_channels * delay, sizeof(PCM_DEC));
    if (data->delay_line == NULL) {
      return -1;
    }
  } else {
    data->delay_line = NULL;
  }
  data->num_channels = num_channels;
  data->delay = delay;

  return 0;
}

void FDK_Delay_Apply(FDK_SignalDelay* data, PCM_DEC* time_buffer,
                     const UINT frame_length, const UCHAR channel) {
  FDK_ASSERT(data != NULL);

  if (data->delay > 0) {
    C_ALLOC_SCRATCH_START(tmp, PCM_DEC, MAX_FRAME_LENGTH)
    FDK_ASSERT(frame_length <= MAX_FRAME_LENGTH);
    FDK_ASSERT(channel < data->num_channels);
    FDK_ASSERT(time_buffer != NULL);
    if (frame_length >= data->delay) {
      FDKmemcpy(tmp, &time_buffer[frame_length - data->delay],
                data->delay * sizeof(PCM_DEC));
      FDKmemmove(&time_buffer[data->delay], &time_buffer[0],
                 (frame_length - data->delay) * sizeof(PCM_DEC));
      FDKmemcpy(&time_buffer[0], &data->delay_line[channel * data->delay],
                data->delay * sizeof(PCM_DEC));
      FDKmemcpy(&data->delay_line[channel * data->delay], tmp,
                data->delay * sizeof(PCM_DEC));
    } else {
      FDKmemcpy(tmp, &time_buffer[0], frame_length * sizeof(PCM_DEC));
      FDKmemcpy(&time_buffer[0], &data->delay_line[channel * data->delay],
                frame_length * sizeof(PCM_DEC));
      FDKmemcpy(&data->delay_line[channel * data->delay],
                &data->delay_line[channel * data->delay + frame_length],
                (data->delay - frame_length) * sizeof(PCM_DEC));
      FDKmemcpy(&data->delay_line[channel * data->delay +
                                  (data->delay - frame_length)],
                tmp, frame_length * sizeof(PCM_DEC));
    }
    C_ALLOC_SCRATCH_END(tmp, PCM_DEC, MAX_FRAME_LENGTH)
  }

  return;
}

void FDK_Delay_Destroy(FDK_SignalDelay* data) {
  if (data->delay_line != NULL) {
    FDKfree(data->delay_line);
  }
  data->delay_line = NULL;
  data->delay = 0;
  data->num_channels = 0;

  return;
}
