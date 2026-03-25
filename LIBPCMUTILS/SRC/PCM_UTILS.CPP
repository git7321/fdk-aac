/**************************** PCM utility library ******************************

   Author(s):   Arthur Tritthart, Alfonso Pino Garcia

   Description: Functions that perform (de)interleaving combined with format
change

*******************************************************************************/

#include "pcm_utils.h"
#include "version.h"

void FDK_interleave(const FIXP_DBL *RESTRICT pIn, LONG *RESTRICT pOut,
                    const UINT channels, const UINT frameSize,
                    const UINT length) {
  for (UINT sample = 0; sample < length; sample++) {
    const FIXP_DBL *In = &pIn[sample];
    for (UINT ch = 0; ch < channels; ch++) {
      *pOut++ = (LONG)In[0];
      In += frameSize;
    }
  }
}

void FDK_interleave(const FIXP_DBL *RESTRICT pIn, SHORT *RESTRICT pOut,
                    const UINT channels, const UINT frameSize,
                    const UINT length) {
  for (UINT sample = 0; sample < length; sample++) {
    const FIXP_DBL *In = &pIn[sample];
    for (UINT ch = 0; ch < channels; ch++) {
      *pOut++ = (SHORT)FX_DBL2FX_SGL(In[0]);
      In += frameSize;
    }
  }
}

void FDK_interleave(const FIXP_SGL *RESTRICT pIn, SHORT *RESTRICT pOut,
                    const UINT channels, const UINT frameSize,
                    const UINT length) {
  for (UINT sample = 0; sample < length; sample++) {
    const FIXP_SGL *In = &pIn[sample];
    for (UINT ch = 0; ch < channels; ch++) {
      *pOut++ = (SHORT)In[0];
      In += frameSize;
    }
  }
}

void FDK_deinterleave(const LONG *RESTRICT pIn, SHORT *RESTRICT _pOut,
                      const UINT channels, const UINT frameSize,
                      const UINT length) {
  for (UINT ch = 0; ch < channels; ch++) {
    SHORT *pOut = _pOut + length * ch;
    const LONG *In = &pIn[ch];
    for (UINT sample = 0; sample < frameSize; sample++) {
      *pOut++ = (SHORT)(In[0] >> 16);
      In += channels;
    }
  }
}

void FDK_deinterleave(const LONG *RESTRICT pIn, LONG *RESTRICT _pOut,
                      const UINT channels, const UINT frameSize,
                      const UINT length) {
  for (UINT ch = 0; ch < channels; ch++) {
    LONG *pOut = _pOut + length * ch;
    const LONG *In = &pIn[ch];
    for (UINT sample = 0; sample < frameSize; sample++) {
      *pOut++ = In[0];
      In += channels;
    }
  }
}

void FDK_deinterleave(const SHORT *RESTRICT pIn, SHORT *RESTRICT _pOut,
                      const UINT channels, const UINT frameSize,
                      const UINT length) {
  for (UINT ch = 0; ch < channels; ch++) {
    SHORT *pOut = _pOut + length * ch;
    const SHORT *In = &pIn[ch];
    for (UINT sample = 0; sample < frameSize; sample++) {
      *pOut++ = In[0];
      In += channels;
    }
  }
}

void FDK_deinterleave(const SHORT *RESTRICT pIn, LONG *RESTRICT _pOut,
                      const UINT channels, const UINT frameSize,
                      const UINT length) {
  for (UINT ch = 0; ch < channels; ch++) {
    LONG *pOut = _pOut + length * ch;
    const SHORT *In = &pIn[ch];
    for (UINT sample = 0; sample < frameSize; sample++) {
      *pOut++ = (LONG)In[0] << 16;
      In += channels;
    }
  }
}
