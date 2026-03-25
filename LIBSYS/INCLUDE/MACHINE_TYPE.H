/************************* System integration library **************************

   Author(s):

   Description:

*******************************************************************************/

#if !defined(MACHINE_TYPE_H)
#define MACHINE_TYPE_H

#include <stddef.h>

#if defined(__ANDROID__) && (__GNUC__ == 4) && (__GNUC_MINOR__ == 4) && \
    (__GNUC_GNU_INLINE__ == 1)
typedef unsigned long long uint64_t;
#include <sys/types.h>
#endif

typedef signed int INT;
typedef unsigned int UINT;
#ifdef __LP64__
#define LONG INT
#define ULONG UINT
#else  /* __LP64__ */
typedef signed long LONG;
typedef unsigned long ULONG;
#endif /* __LP64__ */
typedef signed short SHORT;
typedef unsigned short USHORT;
typedef signed char SCHAR;
typedef unsigned char UCHAR;

#define SHORT_BITS 16
#define CHAR_BITS 8

#ifdef _MSC_VER
typedef __int64 INT64;
typedef unsigned __int64 UINT64;
#else
typedef long long INT64;
typedef unsigned long long UINT64;
#endif

#ifndef NULL
#ifdef __cplusplus
#define NULL 0
#else
#define NULL ((void *)0)
#endif
#endif

#if ((defined(__i686__) || defined(__i586__) || defined(__i386__) ||  \
      defined(__x86_64__)) ||                                         \
     (defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64)))) && \
    !defined(FDK_ASSERT_ENABLE)
#define FDK_ASSERT_ENABLE
#endif

#if defined(FDK_ASSERT_ENABLE)
#include <assert.h>
#define FDK_ASSERT(x) assert(x)
#else
#define FDK_ASSERT(ignore)
#endif

typedef SHORT INT_PCM;
#define MAXVAL_PCM MAXVAL_SGL
#define MINVAL_PCM MINVAL_SGL
#define WAV_BITS 16
#define SAMPLE_BITS 16
#define SAMPLE_MAX ((INT_PCM)(((ULONG)1 << (SAMPLE_BITS - 1)) - 1))
#define SAMPLE_MIN (~SAMPLE_MAX)

#define ALIGNMENT_DEFAULT 8

#if defined(_MSC_VER)
#define RAM_ALIGN __declspec(align(ALIGNMENT_DEFAULT))
#elif defined(__GNUC__)
#define RAM_ALIGN __attribute__((aligned(ALIGNMENT_DEFAULT)))
#else
#define RAM_ALIGN
#endif

#define RESTRICT
#define WORD_ALIGNED(x) C_ALLOC_ALIGNED_CHECK2((const void *)(x), 2);
#define DWORD_ALIGNED(x) C_ALLOC_ALIGNED_CHECK2((const void *)(x), 4);

#define ALIGN_SIZE(a)                                                          \
  ((a) + (((INT)ALIGNMENT_DEFAULT - ((size_t)(a) & (ALIGNMENT_DEFAULT - 1))) & \
          (ALIGNMENT_DEFAULT - 1)))

#define ALIGN_PTR(a)                                      \
  ((void *)((unsigned char *)(a) +                        \
            ((((INT)ALIGNMENT_DEFAULT -                   \
               ((size_t)(a) & (ALIGNMENT_DEFAULT - 1))) & \
              (ALIGNMENT_DEFAULT - 1)))))

#define ALIGNMENT_EXTRES (ALIGNMENT_DEFAULT)
#define ALGN_SIZE_EXTRES(a)                                               \
  ((a) + (((INT)ALIGNMENT_EXTRES - ((INT)(a) & (ALIGNMENT_EXTRES - 1))) & \
          (ALIGNMENT_EXTRES - 1)))

#ifndef FDK_FORCEINLINE
#if defined(__GNUC__) && !defined(__SDE_MIPS__)
#define FDK_FORCEINLINE inline __attribute((always_inline))
#else
#define FDK_FORCEINLINE inline
#endif
#endif

#define FDK_INLINE static inline

#define LNK_SECTION_CODE_L1
#define LNK_SECTION_CODE_L2
#define LNK_SECTION_INITCODE

#define LNK_SECTION_DATA_L1
#define LNK_SECTION_DATA_L2
#define LNK_SECTION_CONSTDATA
#define LNK_SECTION_CONSTDATA_L1

#define LNK_SECTION_L1_DATA_A
#define LNK_SECTION_L1_DATA_B

#ifdef __cplusplus
#if !defined(__has_cpp_attribute)
#define __has_cpp_attribute(x) 0
#endif
#if defined(__clang__) && __has_cpp_attribute(clang::fallthrough)
#define FDK_FALLTHROUGH [[clang::fallthrough]]
#endif
#endif

#ifndef FDK_FALLTHROUGH
#if defined(__GNUC__) && (__GNUC__ >= 7)
#define FDK_FALLTHROUGH __attribute__((fallthrough))
#else
#define FDK_FALLTHROUGH
#endif
#endif

#include <windows.h>

#ifndef FDKCALL
#define FDKCALL __cdecl
#endif

#ifndef LINKIMPORT
 #ifdef FDK_DLL
  #ifdef FDK_EXPORT
   #define LINKIMPORT __declspec(dllexport)
  #else
   #define LINKIMPORT __declspec(dllimport)
  #endif
 #else
  #define LINKIMPORT
 #endif
#endif

#ifndef LINKSPEC
#define LINKSPEC(type) LINKIMPORT type FDKCALL
#endif

#ifdef _MSC_VER
#pragma warning(disable:4065)
#endif

#endif /* MACHINE_TYPE_H */
