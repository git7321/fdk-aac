#if !defined(FIXPOINT_MATH_X86_H)
#define FIXPOINT_MATH_X86_H

#include <math.h>

#define FUNCTION_sqrtFixp

#ifdef FUNCTION_sqrtFixp
static inline FIXP_DBL sqrtFixp(const FIXP_DBL op) {
  FIXP_DBL result;
  result = (FIXP_DBL)(INT)(sqrt((float)(INT)op) * 46340.9492f);
  FDK_ASSERT(result >= (FIXP_DBL)0);
  return result;
}
#endif

#define FUNCTION_invSqrtNorm2
#ifdef FUNCTION_invSqrtNorm2
inline FIXP_DBL invSqrtNorm2(FIXP_DBL op_m, INT *result_e) {
  float result;
  if (op_m == (FIXP_DBL)0) {
    *result_e = 16;
    return ((LONG)0x7fffffff);
  }
  result = (float)(1.0 / sqrt(0.5f * (float)(INT)op_m));
  result = (float)ldexp(frexpf(result, result_e), DFRACT_BITS - 1);
  *result_e += 15;

  FDK_ASSERT(result >= 0);
  return (FIXP_DBL)(INT)result;
}
#endif

#define FUNCTION_invFixp
#ifdef FUNCTION_invFixp
inline FIXP_DBL invFixp(FIXP_DBL op) {
  float result;
  INT result_e;
  if ((op == (FIXP_DBL)0) || (op == (FIXP_DBL)1)) {
    return ((LONG)0x7fffffff);
  }
  result = (float)(1.0 / (float)(INT)op);
  result = frexpf(result, &result_e);
  result = ldexpf(result, 31 + result_e);

  return (FIXP_DBL)(INT)result;
}

inline FIXP_DBL invFixp(FIXP_DBL op_m, int *op_e) {
  float result;
  INT result_e;
  if ((op_m == (FIXP_DBL)0x00000000) || (op_m == (FIXP_DBL)0x00000001)) {
    *op_e = 31 - *op_e;
    return ((LONG)0x7fffffff);
  }
  result = (float)(1.0 / (float)(INT)op_m);
  result = ldexpf(frexpf(result, &result_e), DFRACT_BITS - 1);
  *op_e = result_e - *op_e + 31;
  return (FIXP_DBL)(INT)result;
}
#endif

#define FUNCTION_schur_div
#ifdef FUNCTION_schur_div
inline FIXP_DBL schur_div(FIXP_DBL num, FIXP_DBL denum, INT count) {
  (void)count;
  FDK_ASSERT(num >= (FIXP_DBL)0);
  FDK_ASSERT(denum > (FIXP_DBL)0);
  FDK_ASSERT(num <= denum);

  return (num == denum) ? (FIXP_DBL)MAXVAL_DBL
                        : (FIXP_DBL)(INT)(((INT64)(INT)num << 31) / (INT)denum);
}
#endif /* FUNCTION_schur_div */
#endif /* !defined(FIXPOINT_MATH_X86_H) */
