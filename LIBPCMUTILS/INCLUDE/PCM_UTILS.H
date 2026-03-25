/**************************** PCM utility library ******************************

   Author(s):   Alfonso Pino Garcia

   Description: Functions that perform (de)interleaving combined with format
change

*******************************************************************************/

#if !defined(PCM_UTILS_H)
#define PCM_UTILS_H

#include "common_fix.h"

void FDK_interleave(const FIXP_DBL *RESTRICT pIn, LONG *RESTRICT pOut,
                    const UINT channels, const UINT frameSize,
                    const UINT length);
void FDK_interleave(const FIXP_DBL *RESTRICT pIn, SHORT *RESTRICT pOut,
                    const UINT channels, const UINT frameSize,
                    const UINT length);
void FDK_interleave(const FIXP_SGL *RESTRICT pIn, SHORT *RESTRICT pOut,
                    const UINT channels, const UINT frameSize,
                    const UINT length);

void FDK_deinterleave(const LONG *RESTRICT pIn, SHORT *RESTRICT pOut,
                      const UINT channels, const UINT frameSize,
                      const UINT length);
void FDK_deinterleave(const LONG *RESTRICT pIn, LONG *RESTRICT pOut,
                      const UINT channels, const UINT frameSize,
                      const UINT length);
void FDK_deinterleave(const SHORT *RESTRICT pIn, SHORT *RESTRICT pOut,
                      const UINT channels, const UINT frameSize,
                      const UINT length);
void FDK_deinterleave(const SHORT *RESTRICT pIn, LONG *RESTRICT pOut,
                      const UINT channels, const UINT frameSize,
                      const UINT length);

#endif /* !defined(PCM_UTILS_H) */
