/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description:

*******************************************************************************/

#ifndef OVERLAPADD_H
#define OVERLAPADD_H

#include "common_fix.h"

#define OverlapBufferSize (768)

typedef FIXP_DBL SPECTRUM[1024];
typedef FIXP_DBL* SPECTRAL_PTR;

#define SPEC_LONG(ptr) (ptr)
#define SPEC(ptr, w, gl) ((ptr) + ((w) * (gl)))

#define SPEC_TCX(ptr, f, gl, fb) \
  ((ptr) + ((f) * (gl * 2) * (((fb) == 0) ? 1 : 2)))

#endif /* #ifndef OVERLAPADD_H */
