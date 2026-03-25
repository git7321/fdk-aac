#if !defined(CLZ_X86_H)
#define CLZ_X86_H

#if defined(__GNUC__) && (defined(__x86__) || defined(__x86_64__))

#define FUNCTION_fixnormz_D
#define FUNCTION_fixnorm_D

inline INT fixnormz_D(LONG value) {
  INT result;

  if (value != 0) {
    result = __builtin_clz(value);
  } else {
    result = 32;
  }
  return result;
}

inline INT fixnorm_D(LONG value) {
  INT result;
  if (value == 0) {
    return 0;
  }
  if (value < 0) {
    value = ~value;
  }
  result = fixnormz_D(value);
  return result - 1;
}

#elif (_MSC_VER > 1300) && (defined(_M_IX86) || defined(_M_X64))

#define FUNCTION_fixnormz_D
#define FUNCTION_fixnorm_D
extern "C" unsigned char _BitScanReverse(unsigned long* Index, unsigned long Mask);

inline INT fixnormz_D(LONG value) {
  unsigned long result = 0;
  unsigned char err;
  err = _BitScanReverse(&result, value);
  if (err) {
    return 31 - result;
  } else {
    return 32;
  }
}

inline INT fixnorm_D(LONG value) {
  INT result;
  if (value == 0) {
    return 0;
  }
  if (value < 0) {
    value = ~value;
  }
  result = fixnormz_D(value);
  return result - 1;
}

#endif /* toolchain */
#endif /* !defined(CLZ_X86_H) */
