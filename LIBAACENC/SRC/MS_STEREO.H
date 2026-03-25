/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: MS stereo processing

*******************************************************************************/

#ifndef MS_STEREO_H
#define MS_STEREO_H

#include "interface.h"

void FDKaacEnc_MsStereoProcessing(PSY_DATA *RESTRICT psyData[(2)],
                                  PSY_OUT_CHANNEL *psyOutChannel[2],
                                  const INT *isBook, INT *msDigest,
                                  INT *msMask,
                                  const INT allowMS, const INT sfbCnt,
                                  const INT sfbPerGroup,
                                  const INT maxSfbPerGroup,
                                  const INT *sfbOffset);

#endif /* MS_STEREO_H */
