/******************* Library for basic calculation routines ********************

   Author(s):

   Description: Scaling operations

*******************************************************************************/

#include "common_fix.h"
#include "genericStds.h"
#include "scale.h"

#if defined(__mips__)
#include "mips/scale_mips.cpp"

#elif defined(__arm__)
#include "arm/scale_arm.cpp"

#endif

#ifndef FUNCTION_scaleValues_SGL
#define FUNCTION_scaleValues_SGL
void scaleValues(FIXP_SGL *vector,
                 INT len,
                 INT scalefactor
) {
  INT i;

  if (scalefactor == 0) return;

  if (scalefactor > 0) {
    scalefactor = fixmin_I(scalefactor, (INT)(FRACT_BITS - 1));
    for (i = len & 3; i--;) {
      *(vector++) <<= scalefactor;
    }
    for (i = len >> 2; i--;) {
      *(vector++) <<= scalefactor;
      *(vector++) <<= scalefactor;
      *(vector++) <<= scalefactor;
      *(vector++) <<= scalefactor;
    }
  } else {
    INT negScalefactor = fixmin_I(-scalefactor, (INT)FRACT_BITS - 1);
    for (i = len & 3; i--;) {
      *(vector++) >>= negScalefactor;
    }
    for (i = len >> 2; i--;) {
      *(vector++) >>= negScalefactor;
      *(vector++) >>= negScalefactor;
      *(vector++) >>= negScalefactor;
      *(vector++) >>= negScalefactor;
    }
  }
}
#endif

#ifndef FUNCTION_scaleValues_DBL
#define FUNCTION_scaleValues_DBL
void scaleValues(FIXP_DBL *vector,
                 INT len,
                 INT scalefactor
) {
  INT i;

  if (scalefactor == 0) return;

  if (scalefactor > 0) {
    scalefactor = fixmin_I(scalefactor, (INT)DFRACT_BITS - 1);
    for (i = len & 3; i--;) {
      *(vector++) <<= scalefactor;
    }
    for (i = len >> 2; i--;) {
      *(vector++) <<= scalefactor;
      *(vector++) <<= scalefactor;
      *(vector++) <<= scalefactor;
      *(vector++) <<= scalefactor;
    }
  } else {
    INT negScalefactor = fixmin_I(-scalefactor, (INT)DFRACT_BITS - 1);
    for (i = len & 3; i--;) {
      *(vector++) >>= negScalefactor;
    }
    for (i = len >> 2; i--;) {
      *(vector++) >>= negScalefactor;
      *(vector++) >>= negScalefactor;
      *(vector++) >>= negScalefactor;
      *(vector++) >>= negScalefactor;
    }
  }
}
#endif

#ifndef FUNCTION_scaleValuesSaturate_DBL
#define FUNCTION_scaleValuesSaturate_DBL
void scaleValuesSaturate(FIXP_DBL *vector,
                         INT len,
                         INT scalefactor
) {
  INT i;

  if (scalefactor == 0) return;

  scalefactor = fixmax_I(fixmin_I(scalefactor, (INT)DFRACT_BITS - 1),
                         (INT) - (DFRACT_BITS - 1));

  for (i = 0; i < len; i++) {
    vector[i] = scaleValueSaturate(vector[i], scalefactor);
  }
}
#endif /* FUNCTION_scaleValuesSaturate_DBL */

#ifndef FUNCTION_scaleValuesSaturate_DBL_DBL
#define FUNCTION_scaleValuesSaturate_DBL_DBL
void scaleValuesSaturate(FIXP_DBL *dst,
                         const FIXP_DBL *src,
                         INT len,
                         INT scalefactor
) {
  INT i;

  if (scalefactor == 0) {
    FDKmemmove(dst, src, len * sizeof(FIXP_DBL));
    return;
  }

  scalefactor = fixmax_I(fixmin_I(scalefactor, (INT)DFRACT_BITS - 1),
                         (INT) - (DFRACT_BITS - 1));

  for (i = 0; i < len; i++) {
    dst[i] = scaleValueSaturate(src[i], scalefactor);
  }
}
#endif /* FUNCTION_scaleValuesSaturate_DBL_DBL */

#ifndef FUNCTION_scaleValuesSaturate_SGL_DBL
#define FUNCTION_scaleValuesSaturate_SGL_DBL
void scaleValuesSaturate(FIXP_SGL *dst,
                         const FIXP_DBL *src,
                         INT len,
                         INT scalefactor)
{
  INT i;
  scalefactor = fixmax_I(fixmin_I(scalefactor, (INT)DFRACT_BITS - 1),
                         (INT) - (DFRACT_BITS - 1));

  for (i = 0; i < len; i++) {
    dst[i] = FX_DBL2FX_SGL(fAddSaturate(scaleValueSaturate(src[i], scalefactor),
                                        (FIXP_DBL)0x8000));
  }
}
#endif /* FUNCTION_scaleValuesSaturate_SGL_DBL */

#ifndef FUNCTION_scaleValuesSaturate_SGL
#define FUNCTION_scaleValuesSaturate_SGL
void scaleValuesSaturate(FIXP_SGL *vector,
                         INT len,
                         INT scalefactor
) {
  INT i;

  if (scalefactor == 0) return;

  scalefactor = fixmax_I(fixmin_I(scalefactor, (INT)DFRACT_BITS - 1),
                         (INT) - (DFRACT_BITS - 1));

  for (i = 0; i < len; i++) {
    vector[i] = FX_DBL2FX_SGL(
        scaleValueSaturate(FX_SGL2FX_DBL(vector[i]), scalefactor));
  }
}
#endif /* FUNCTION_scaleValuesSaturate_SGL */

#ifndef FUNCTION_scaleValuesSaturate_SGL_SGL
#define FUNCTION_scaleValuesSaturate_SGL_SGL
void scaleValuesSaturate(FIXP_SGL *dst,
                         const FIXP_SGL *src,
                         INT len,
                         INT scalefactor
) {
  INT i;

  if (scalefactor == 0) {
    FDKmemmove(dst, src, len * sizeof(FIXP_SGL));
    return;
  }

  scalefactor = fixmax_I(fixmin_I(scalefactor, (INT)DFRACT_BITS - 1),
                         (INT) - (DFRACT_BITS - 1));

  for (i = 0; i < len; i++) {
    dst[i] =
        FX_DBL2FX_SGL(scaleValueSaturate(FX_SGL2FX_DBL(src[i]), scalefactor));
  }
}
#endif /* FUNCTION_scaleValuesSaturate_SGL_SGL */

#ifndef FUNCTION_scaleValues_DBLDBL
#define FUNCTION_scaleValues_DBLDBL
void scaleValues(FIXP_DBL *dst,
                 const FIXP_DBL *src,
                 INT len,
                 INT scalefactor
) {
  INT i;

  if (scalefactor == 0) {
    if (dst != src) FDKmemmove(dst, src, len * sizeof(FIXP_DBL));
  } else {
    if (scalefactor > 0) {
      scalefactor = fixmin_I(scalefactor, (INT)DFRACT_BITS - 1);
      for (i = len & 3; i--;) {
        *(dst++) = *(src++) << scalefactor;
      }
      for (i = len >> 2; i--;) {
        *(dst++) = *(src++) << scalefactor;
        *(dst++) = *(src++) << scalefactor;
        *(dst++) = *(src++) << scalefactor;
        *(dst++) = *(src++) << scalefactor;
      }
    } else {
      INT negScalefactor = fixmin_I(-scalefactor, (INT)DFRACT_BITS - 1);
      for (i = len & 3; i--;) {
        *(dst++) = *(src++) >> negScalefactor;
      }
      for (i = len >> 2; i--;) {
        *(dst++) = *(src++) >> negScalefactor;
        *(dst++) = *(src++) >> negScalefactor;
        *(dst++) = *(src++) >> negScalefactor;
        *(dst++) = *(src++) >> negScalefactor;
      }
    }
  }
}
#endif

#if (SAMPLE_BITS == 16)
#ifndef FUNCTION_scaleValues_PCMDBL
#define FUNCTION_scaleValues_PCMDBL
void scaleValues(FIXP_PCM *dst,
                 const FIXP_DBL *src,
                 INT len,
                 INT scalefactor
) {
  INT i;

  scalefactor -= DFRACT_BITS - SAMPLE_BITS;

  {
    if (scalefactor > 0) {
      scalefactor = fixmin_I(scalefactor, (INT)DFRACT_BITS - 1);
      for (i = len & 3; i--;) {
        *(dst++) = (FIXP_PCM)(*(src++) << scalefactor);
      }
      for (i = len >> 2; i--;) {
        *(dst++) = (FIXP_PCM)(*(src++) << scalefactor);
        *(dst++) = (FIXP_PCM)(*(src++) << scalefactor);
        *(dst++) = (FIXP_PCM)(*(src++) << scalefactor);
        *(dst++) = (FIXP_PCM)(*(src++) << scalefactor);
      }
    } else {
      INT negScalefactor = fixmin_I(-scalefactor, (INT)DFRACT_BITS - 1);
      for (i = len & 3; i--;) {
        *(dst++) = (FIXP_PCM)(*(src++) >> negScalefactor);
      }
      for (i = len >> 2; i--;) {
        *(dst++) = (FIXP_PCM)(*(src++) >> negScalefactor);
        *(dst++) = (FIXP_PCM)(*(src++) >> negScalefactor);
        *(dst++) = (FIXP_PCM)(*(src++) >> negScalefactor);
        *(dst++) = (FIXP_PCM)(*(src++) >> negScalefactor);
      }
    }
  }
}
#endif
#endif /* (SAMPLE_BITS == 16) */

#ifndef FUNCTION_scaleValues_SGLSGL
#define FUNCTION_scaleValues_SGLSGL
void scaleValues(FIXP_SGL *dst,
                 const FIXP_SGL *src,
                 INT len,
                 INT scalefactor
) {
  INT i;

  if (scalefactor == 0) {
    if (dst != src) FDKmemmove(dst, src, len * sizeof(FIXP_DBL));
  } else {
    if (scalefactor > 0) {
      scalefactor = fixmin_I(scalefactor, (INT)DFRACT_BITS - 1);
      for (i = len & 3; i--;) {
        *(dst++) = *(src++) << scalefactor;
      }
      for (i = len >> 2; i--;) {
        *(dst++) = *(src++) << scalefactor;
        *(dst++) = *(src++) << scalefactor;
        *(dst++) = *(src++) << scalefactor;
        *(dst++) = *(src++) << scalefactor;
      }
    } else {
      INT negScalefactor = fixmin_I(-scalefactor, (INT)DFRACT_BITS - 1);
      for (i = len & 3; i--;) {
        *(dst++) = *(src++) >> negScalefactor;
      }
      for (i = len >> 2; i--;) {
        *(dst++) = *(src++) >> negScalefactor;
        *(dst++) = *(src++) >> negScalefactor;
        *(dst++) = *(src++) >> negScalefactor;
        *(dst++) = *(src++) >> negScalefactor;
      }
    }
  }
}
#endif

#ifndef FUNCTION_scaleValuesWithFactor_DBL
#define FUNCTION_scaleValuesWithFactor_DBL
void scaleValuesWithFactor(FIXP_DBL *vector, FIXP_DBL factor, INT len,
                           INT scalefactor) {
  INT i;

  scalefactor++;

  if (scalefactor > 0) {
    scalefactor = fixmin_I(scalefactor, (INT)DFRACT_BITS - 1);
    for (i = len & 3; i--;) {
      *vector = fMultDiv2(*vector, factor) << scalefactor;
      vector++;
    }
    for (i = len >> 2; i--;) {
      *vector = fMultDiv2(*vector, factor) << scalefactor;
      vector++;
      *vector = fMultDiv2(*vector, factor) << scalefactor;
      vector++;
      *vector = fMultDiv2(*vector, factor) << scalefactor;
      vector++;
      *vector = fMultDiv2(*vector, factor) << scalefactor;
      vector++;
    }
  } else {
    INT negScalefactor = fixmin_I(-scalefactor, (INT)DFRACT_BITS - 1);
    for (i = len & 3; i--;) {
      *vector = fMultDiv2(*vector, factor) >> negScalefactor;
      vector++;
    }
    for (i = len >> 2; i--;) {
      *vector = fMultDiv2(*vector, factor) >> negScalefactor;
      vector++;
      *vector = fMultDiv2(*vector, factor) >> negScalefactor;
      vector++;
      *vector = fMultDiv2(*vector, factor) >> negScalefactor;
      vector++;
      *vector = fMultDiv2(*vector, factor) >> negScalefactor;
      vector++;
    }
  }
}
#endif /* FUNCTION_scaleValuesWithFactor_DBL */

#ifndef FUNCTION_getScalefactorShort
#define FUNCTION_getScalefactorShort
INT getScalefactorShort(const SHORT *vector,
                        INT len
) {
  INT i;
  SHORT temp, maxVal = 0;

  for (i = len; i != 0; i--) {
    temp = (SHORT)(*vector++);
    maxVal |= (temp ^ (temp >> (SHORT_BITS - 1)));
  }

  return fixmax_I((INT)0, (INT)(fixnormz_D((INT)maxVal) - (INT)1 -
                                (INT)(DFRACT_BITS - SHORT_BITS)));
}
#endif

#ifndef FUNCTION_getScalefactorPCM
#define FUNCTION_getScalefactorPCM
INT getScalefactorPCM(const INT_PCM *vector,
                      INT len,
                      INT stride) {
  INT i;
  INT_PCM temp, maxVal = 0;

  for (i = len; i != 0; i--) {
    temp = (INT_PCM)(*vector);
    vector += stride;
    maxVal |= (temp ^ (temp >> ((sizeof(INT_PCM) * 8) - 1)));
  }
  return fixmax_I((INT)0, (INT)(fixnormz_D((INT)maxVal) - (INT)1 -
                                (INT)(DFRACT_BITS - SAMPLE_BITS)));
}
#endif

#ifndef FUNCTION_getScalefactorShort
#define FUNCTION_getScalefactorShort
INT getScalefactorShort(const SHORT *vector,
                        INT len,
                        INT stride) {
  INT i;
  SHORT temp, maxVal = 0;

  for (i = len; i != 0; i--) {
    temp = (SHORT)(*vector);
    vector += stride;
    maxVal |= (temp ^ (temp >> (SHORT_BITS - 1)));
  }

  return fixmax_I((INT)0, (INT)(fixnormz_D((INT)maxVal) - (INT)1 -
                                (INT)(DFRACT_BITS - SHORT_BITS)));
}
#endif

#ifndef FUNCTION_getScalefactor_DBL
#define FUNCTION_getScalefactor_DBL
INT getScalefactor(const FIXP_DBL *vector,
                   INT len)
{
  INT i;
  FIXP_DBL temp, maxVal = (FIXP_DBL)0;

  for (i = len; i != 0; i--) {
    temp = (LONG)(*vector++);
    maxVal |= (FIXP_DBL)((LONG)temp ^ (LONG)(temp >> (DFRACT_BITS - 1)));
  }

  return fixmax_I((INT)0, (INT)(fixnormz_D(maxVal) - 1));
}
#endif

#ifndef FUNCTION_getScalefactor_SGL
#define FUNCTION_getScalefactor_SGL
INT getScalefactor(const FIXP_SGL *vector,
                   INT len)
{
  INT i;
  SHORT temp, maxVal = (FIXP_SGL)0;

  for (i = len; i != 0; i--) {
    temp = (SHORT)(*vector++);
    maxVal |= (temp ^ (temp >> (FRACT_BITS - 1)));
  }

  return fixmax_I((INT)0, (INT)(fixnormz_S((FIXP_SGL)maxVal)) - 1);
}
#endif
