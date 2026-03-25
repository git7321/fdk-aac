/**************************** AAC decoder library ******************************

   Author(s):   Matthias Hildenbrand

   Description: USAC ACELP LTP filter

*******************************************************************************/

#ifndef USACDEC_ACE_LTP_H
#define USACDEC_ACE_LTP_H

#include "common_fix.h"

void Pred_lt4(FIXP_DBL exc[],
              int T0,
              int frac
);

void Pred_lt4_postfilter(FIXP_DBL exc[]
);

#endif /* USACDEC_ACE_LTP_H */
