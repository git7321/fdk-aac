#if !defined(ABS_X86_H)
#define ABS_X86_H

#if defined(__x86__)
#if defined(__x86_64__)

inline INT fixabs_D(INT x) { return ((x) > (INT)(0)) ? (x) : -(x); }
inline INT fixabs_S(INT x) { return ((x) > (INT)(0)) ? (x) : -(x); }

#define fixabs_I(x) fixabs_D(x)

#define FUNCTION_fixabs_S
#define FUNCTION_fixabs_D
#define FUNCTION_fixabs_I

#endif /* __x86_64__ */
#endif /*__x86__ */

#endif /* !defined(ABS_X86_H) */
