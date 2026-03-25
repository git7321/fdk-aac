/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef SBR_MISC_H
#define SBR_MISC_H

#include "sbr_encoder.h"

void FDKsbrEnc_Shellsort_fract(FIXP_DBL *in, INT n);
void FDKsbrEnc_Shellsort_int(INT *in, INT n);

void FDKsbrEnc_AddLeft(INT *vector, INT *length_vector, INT value);
void FDKsbrEnc_AddRight(INT *vector, INT *length_vector, INT value);
void FDKsbrEnc_AddVecLeft(INT *dst, INT *length_dst, INT *src, INT length_src);
void FDKsbrEnc_AddVecRight(INT *dst, INT *length_vector_dst, INT *src,
                           INT length_src);

FIXP_DBL FDKsbrEnc_LSI_divide_scale_fract(FIXP_DBL num, FIXP_DBL denom,
                                          FIXP_DBL scale);

#endif
