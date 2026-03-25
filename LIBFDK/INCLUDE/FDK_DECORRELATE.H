/******************* Library for basic calculation routines ********************

   Author(s):   Markus Lohwasser

   Description: FDK Tools Decorrelator

*******************************************************************************/

#ifndef FDK_DECORRELATE_H
#define FDK_DECORRELATE_H

#include "common_fix.h"

#define FIXP_MPS FIXP_DBL

#ifndef ARCH_PREFER_MULT_32x32
#define FIXP_DECORR FIXP_SGL
#define FX_DECORR2FX_DBL FX_SGL2FX_DBL
#define FX_DECORR2FX_SGL
#define FX_DBL2FX_DECORR FX_DBL2FX_SGL
#define FX_SGL2FX_DECORR
#define DECORR(a) (FX_DBL2FXCONST_SGL(a))
#define FL2FXCONST_DECORR FL2FXCONST_SGL
#else
#define FIXP_DECORR FIXP_DBL
#define FX_DECORR2FX_DBL
#define FX_DECORR2FX_SGL FX_DBL2FX_SGL
#define FX_DBL2FX_DECORR
#define FX_SGL2FX_DECORR FX_SGL2FX_DBL
#define DECORR(a) FIXP_DBL(a)
#define FL2FXCONST_DECORR FL2FXCONST_DBL
#endif

typedef enum {
  DECORR_MPS,
  DECORR_PS,
  DECORR_USAC,
  DECORR_LD
} FDK_DECORR_TYPE;

typedef enum {
  DUCKER_AUTOMATIC,
  DUCKER_MPS,
  DUCKER_PS
} FDK_DUCKER_TYPE;

typedef enum {
  NOT_EXIST,
  DELAY,
  COMMON_REAL,
  COMMON_CPLX,
  INDEP_CPLX,
  INDEP_CPLX_PS
} REVBAND_FILT_TYPE;

typedef struct DECORR_DEC *HANDLE_DECORR_DEC;

typedef struct DUCKER_INSTANCE {
  int hybridBands;
  int parameterBands;
  int partiallyComplex;
  FDK_DUCKER_TYPE duckerType;

  const UCHAR *qs_next;
  const UCHAR *mapProcBands2HybBands;
  const UCHAR *mapHybBands2ProcBands;
  FIXP_MPS SmoothDirRevNrg[2 * (28)];

  FIXP_MPS peakDecay[(28)];
  FIXP_MPS peakDiff[(28)];
  FIXP_DBL maxValDirectData;
  FIXP_DBL maxValReverbData;
  SCHAR scaleDirectNrg;
  SCHAR scaleReverbNrg;
  SCHAR scaleSmoothDirRevNrg;
  SCHAR headroomSmoothDirRevNrg;

} DUCKER_INSTANCE;

typedef struct DECORR_FILTER_INSTANCE {
  FIXP_MPS *stateCplx;
  FIXP_DBL *DelayBufferCplx;

  const FIXP_DECORR *numeratorReal;
  const FIXP_STP *coeffsPacked;
  const FIXP_DECORR *denominatorReal;
} DECORR_FILTER_INSTANCE;

typedef struct DECORR_DEC {
  INT L_stateBufferCplx;
  FIXP_DBL *stateBufferCplx;
  INT L_delayBufferCplx;
  FIXP_DBL *delayBufferCplx;

  const REVBAND_FILT_TYPE *REV_filtType;
  const UCHAR *REV_bandOffset;
  const UCHAR *REV_delay;
  const SCHAR *REV_filterOrder;
  INT reverbBandDelayBufferIndex[(4)];
  UCHAR stateBufferOffset[(3)];

  DECORR_FILTER_INSTANCE Filter[(71)];
  DUCKER_INSTANCE ducker;

  int numbins;
  int partiallyComplex;
} DECORR_DEC;

INT FDKdecorrelateOpen(HANDLE_DECORR_DEC hDecorrDec, FIXP_DBL *bufferCplx,
                       const INT bufLen);

INT FDKdecorrelateInit(HANDLE_DECORR_DEC hDecorrDec, const INT nrHybBands,
                       const FDK_DECORR_TYPE decorrType,
                       const FDK_DUCKER_TYPE duckerType, const INT decorrConfig,
                       const INT seed, const INT partiallyComplex,
                       const INT useFractDelay, const INT isLegacyPS,
                       const INT initStatesFlag);

INT FDKdecorrelateApply(HANDLE_DECORR_DEC hDecorrDec, FIXP_DBL *dataRealIn,
                        FIXP_DBL *dataImagIn, FIXP_DBL *dataRealOut,
                        FIXP_DBL *dataImagOut, const INT startHybBand);

INT FDKdecorrelateClose(HANDLE_DECORR_DEC hDecorrDec);

FIXP_DBL *getAddrDirectSignalMaxVal(HANDLE_DECORR_DEC hDecorrDec);

#endif /* FDK_DECORRELATE_H */
