#if !defined(FIXMUL_X86_H)
#define FIXMUL_X86_H

#if defined(__x86__)
#if defined(_MSC_VER) && defined(_M_IX86)

#define FUNCTION_fixmul_DD
#define FUNCTION_fixmuldiv2_DD
#define FUNCTION_fixmuldiv2BitExact_DD
#define fixmuldiv2BitExact_DD(a, b) fixmuldiv2_DD(a, b)
#define FUNCTION_fixmulBitExact_DD
#define fixmulBitExact_DD(a, b) fixmul_DD(a, b)

#define FUNCTION_fixmuldiv2BitExact_DS
#define fixmuldiv2BitExact_DS(a, b) fixmuldiv2_DS(a, b)

#define FUNCTION_fixmulBitExact_DS
#define fixmulBitExact_DS(a, b) fixmul_DS(a, b)

inline INT fixmul_DD(INT a, const INT b) {
  __asm
  {
    mov eax, a
    imul b
    shl edx, 1
    mov a, edx
  }
  return a;
}

inline INT fixmuldiv2_DD(INT a, const INT b) {
  __asm
  {
    mov eax, a
    imul b
    mov a, edx
  }
  return a;
}

#elif (defined(__GNUC__) || defined(__gnu_linux__)) && defined(__x86__)

#define FUNCTION_fixmul_DD
#define FUNCTION_fixmuldiv2_DD

#define FUNCTION_fixmuldiv2BitExact_DD
#define fixmuldiv2BitExact_DD(a, b) fixmuldiv2_DD(a, b)

#define FUNCTION_fixmulBitExact_DD
#define fixmulBitExact_DD(a, b) fixmul_DD(a, b)

#define FUNCTION_fixmuldiv2BitExact_DS
#define fixmuldiv2BitExact_DS(a, b) fixmuldiv2_DS(a, b)

#define FUNCTION_fixmulBitExact_DS
#define fixmulBitExact_DS(a, b) fixmul_DS(a, b)

inline INT fixmul_DD(INT a, const INT b) {
  INT result;

  asm("imul %2;\n"
      "shl $1, %0;\n"
      : "=d"(result), "+a"(a)
      : "r"(b));

  return result;
}

inline INT fixmuldiv2_DD(INT a, const INT b) {
  INT result;

  asm("imul %2;" : "=d"(result), "+a"(a) : "r"(b));

  return result;
}

#endif /* (defined(__GNUC__)||defined(__gnu_linux__)) && defined(__x86__) */
#endif /* __x86__ */

#endif /* !defined(FIXMUL_X86_H) */
