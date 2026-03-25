/************************** Arithmetic coder library ***************************

   Author(s):   Oliver Weiss

   Description: Interface for Spectral Noiseless Coding Scheme based on an
                Arithmetic Coder in Conjunction with an Adaptive Context

*******************************************************************************/

#ifndef AC_ARITH_CODER_H
#define AC_ARITH_CODER_H

#include "common_fix.h"
#include "FDK_bitstream.h"

#include "FDK_audio.h"

typedef enum { ARITH_CODER_OK = 0, ARITH_CODER_ERROR = 5 } ARITH_CODING_ERROR;

typedef struct {
  SHORT m_numberLinesPrev;
  UCHAR c_prev[(1024 / 2) + 4];
} CArcoData;

CArcoData *CArco_Create(void);

void CArco_Destroy(CArcoData *pArcoData);

ARITH_CODING_ERROR CArco_DecodeArithData(CArcoData *pArcoData,
                                         HANDLE_FDK_BITSTREAM hBs,
                                         FIXP_DBL *RESTRICT spectrum, int lg,
                                         int lg_max, int arith_reset_flag);

#endif /* AC_ARITH_CODER_H */
