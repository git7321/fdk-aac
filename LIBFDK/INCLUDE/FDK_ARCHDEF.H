/******************* Library for basic calculation routines ********************

   Author(s):

   Description:

*******************************************************************************/

#ifndef FDK_ARCHDEF_H
#define FDK_ARCHDEF_H

#if defined(__i386__) || defined(__i486__) || defined(__i586__) ||  \
    defined(__i686__) || (defined(_MSC_VER) && defined(_M_IX86)) || \
    (defined(_MSC_VER) && defined(_M_X64)) || defined(__x86_64__)
#define __x86__
#endif

#if defined(_M_ARM) && !defined(__arm__) || defined(__aarch64__) || defined(_M_ARM64)
#define __arm__
#endif

#if defined(_ARCH_PPC) && !defined(__powerpc__)
#define __powerpc__ 1
#endif

#if (__TARGET_ARCH_ARM == 5) || defined(__TARGET_FEATURE_DSPMUL) || \
    (_M_ARM == 5) || defined(__ARM_ARCH_5TEJ__) || defined(__ARM_ARCH_7EM__)
#define __ARM_ARCH_5TE__
#endif

#if (__TARGET_ARCH_ARM == 6) || defined(__ARM_ARCH_6J__) || \
    defined(__ARM_ARCH_6ZK__)
#define __ARM_ARCH_5TE__
#define __ARM_ARCH_6__
#endif

#if defined(__TARGET_ARCH_7_R) || defined(__ARM_ARCH_7R__)
#define __ARM_ARCH_5TE__
#define __ARM_ARCH_6__
#define __ARM_ARCH_7_R__
#endif

#if defined(__TARGET_ARCH_7_A) || defined(__ARM_ARCH_7A__) || \
    ((__ARM_ARCH == 8) && (__ARM_32BIT_STATE == 1))
#define __ARM_ARCH_5TE__
#define __ARM_ARCH_6__
#define __ARM_ARCH_7_A__
#endif

#if defined(__TARGET_ARCH_7_M) || defined(__ARM_ARCH_7_M__)
#define __ARM_ARCH_7M__
#endif

#if defined(__TARGET_ARCH_7E_M) || defined(__ARM_ARCH_7E_M__)
#define __ARM_ARCH_7EM__
#endif

#if defined(__aarch64__) || defined(_M_ARM64)
#define __ARM_ARCH_8__
#endif

#ifdef _M_ARM
#include "armintr.h"
#endif

#if defined(__mips__)
#define ARCH_PREFER_MULT_16x16
#undef SINETABLE_16BIT
#undef POW2COEFF_16BIT
#undef LDCOEFF_16BIT
#undef WINDOWTABLE_16BIT

#elif defined(__arm__) && defined(__ARM_ARCH_8__)
#define ARCH_PREFER_MULT_32x16
#define SINETABLE_16BIT
#define POW2COEFF_16BIT
#define LDCOEFF_16BIT
#define WINDOWTABLE_16BIT

#elif defined(__arm__) && defined(__ARM_ARCH_5TE__)
#define ARCH_PREFER_MULT_32x16
#define SINETABLE_16BIT
#define POW2COEFF_16BIT
#define LDCOEFF_16BIT
#define WINDOWTABLE_16BIT

#elif defined(__arm__) && defined(__ARM_ARCH_7M__)
#define ARCH_PREFER_MULT_32x16
#define SINETABLE_16BIT
#define POW2COEFF_16BIT
#define LDCOEFF_16BIT
#define WINDOWTABLE_16BIT

#elif defined(__arm__) && defined(__ARM_ARCH_7EM__)
#define ARCH_PREFER_MULT_32x32
#define ARCH_PREFER_MULT_32x16
#define SINETABLE_16BIT
#define POW2COEFF_16BIT
#define LDCOEFF_16BIT
#define WINDOWTABLE_16BIT

#elif defined(__arm__) && !defined(__ARM_ARCH_5TE__)
#define ARCH_PREFER_MULT_16x16
#undef SINETABLE_16BIT
#undef WINDOWTABLE_16BIT
#undef POW2COEFF_16BIT
#undef LDCOEFF_16BIT

#elif defined(__x86__)
#define ARCH_PREFER_MULT_32x16
#define SINETABLE_16BIT
#define WINDOWTABLE_16BIT
#define POW2COEFF_16BIT
#define LDCOEFF_16BIT

#elif defined(__riscv)
#define ARCH_PREFER_MULT_32x16
#define SINETABLE_16BIT
#define POW2COEFF_16BIT
#define LDCOEFF_16BIT
#define WINDOWTABLE_16BIT

#elif defined(__powerpc__)
#define ARCH_PREFER_MULT_32x32

#elif defined(__s390x__)
#define ARCH_PREFER_MULT_32x32
#define ARCH_PREFER_MULT_32x16
#define SINETABLE_16BIT
#define POW2COEFF_16BIT
#define LDCOEFF_16BIT
#define WINDOWTABLE_16BIT

#else
#warning >>>> Please set architecture characterization defines for your platform (FDK_HIGH_PERFORMANCE)! <<<<

#endif

#ifdef SINETABLE_16BIT
#define FIXP_STB FIXP_SGL
#define FIXP_STP FIXP_SPK
#define STC(a) (FX_DBL2FXCONST_SGL(a))
#else
#define FIXP_STB FIXP_DBL
#define FIXP_STP FIXP_DPK
#define STC(a) ((FIXP_DBL)(LONG)(a))
#endif

#define STCP(cos, sin)     \
  {                        \
    { STC(cos), STC(sin) } \
  }

#ifdef WINDOWTABLE_16BIT
#define FIXP_WTB FIXP_SGL
#define FX_DBL2FX_WTB(x) FX_DBL2FX_SGL(x)
#define FIXP_WTP FIXP_SPK
#define WTC(a) FX_DBL2FXCONST_SGL(a)
#else
#define FIXP_WTB FIXP_DBL
#define FX_DBL2FX_WTB(x) (x)
#define FIXP_WTP FIXP_DPK
#define WTC(a) (FIXP_DBL)(a)
#endif

#define WTCP(a, b)     \
  {                    \
    { WTC(a), WTC(b) } \
  }

#endif /* FDK_ARCHDEF_H */
