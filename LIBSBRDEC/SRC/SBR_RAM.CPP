/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "sbr_ram.h"

C_ALLOC_MEM(Ram_SbrDecoder, struct SBR_DECODER_INSTANCE, 1)
C_ALLOC_MEM2(Ram_SbrDecElement, SBR_DECODER_ELEMENT, 1, (8))
C_ALLOC_MEM2(Ram_SbrDecChannel, SBR_CHANNEL, 1, (8) + 1)
C_ALLOC_MEM(Ram_ps_dec, struct PS_DEC, 1)
