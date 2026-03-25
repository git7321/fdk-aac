/******************* Library for basic calculation routines ********************

   Author(s):

   Description: Library functions to calculate standard DCTs. This will most
                likely be replaced by hand-optimized functions for the specific
                target processor.

*******************************************************************************/

#ifndef DCT_H
#define DCT_H

#include "common_fix.h"

void dct_getTables(const FIXP_WTP **ptwiddle, const FIXP_STP **sin_twiddle,
                   int *sin_step, int length);

void dct_II(FIXP_DBL *pDat, FIXP_DBL *tmp, int size, int *pDat_e);
void dct_III(FIXP_DBL *pDat, FIXP_DBL *tmp, int size, int *pDat_e);
void dst_III(FIXP_DBL *pDat, FIXP_DBL *tmp, int size, int *pDat_e);
void dct_IV(FIXP_DBL *pDat, int size, int *pDat_e);
void dst_IV(FIXP_DBL *pDat, int size, int *pDat_e);

#endif
