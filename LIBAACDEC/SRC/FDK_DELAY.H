/**************************** AAC decoder library ******************************

   Author(s):   Matthias Hildenbrand

   Description:

*******************************************************************************/

#ifndef FDK_DELAY_H
#define FDK_DELAY_H

#include "aac_rom.h"

typedef struct {
  PCM_DEC* delay_line;
  USHORT delay;
  UCHAR num_channels;
} FDK_SignalDelay;

INT FDK_Delay_Create(FDK_SignalDelay* data, const USHORT delay,
                     const UCHAR num_channels);

void FDK_Delay_Apply(FDK_SignalDelay* data, PCM_DEC* time_buffer,
                     const UINT frame_length, const UCHAR channel);

void FDK_Delay_Destroy(FDK_SignalDelay* data);

#endif /* #ifndef FDK_DELAY_H */
