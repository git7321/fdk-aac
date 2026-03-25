/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: temporal noise shaping tool

*******************************************************************************/

#ifndef AACDEC_TNS_H
#define AACDEC_TNS_H

#include "common_fix.h"

enum {
  TNS_MAX_WINDOWS = 8,
  TNS_MAXIMUM_FILTERS = 3
};

#define TNS_MAXIMUM_ORDER (20)

#if (TNS_MAXIMUM_ORDER < 15)
#error USAC: TNS filter order up 15 can be signaled!
#endif

typedef struct {
  SCHAR Coeff[TNS_MAXIMUM_ORDER];

  UCHAR StartBand;
  UCHAR StopBand;

  SCHAR Direction;
  SCHAR Resolution;

  UCHAR Order;
} CFilter;

typedef struct {
  CFilter Filter[TNS_MAX_WINDOWS][TNS_MAXIMUM_FILTERS];
  UCHAR NumberOfFilters[TNS_MAX_WINDOWS];
  UCHAR DataPresent;
  UCHAR Active;

  UCHAR GainLd;
} CTnsData;

void CTns_Reset(CTnsData *pTnsData);

#endif /* #ifndef AACDEC_TNS_H */
