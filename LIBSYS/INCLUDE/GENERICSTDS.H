#if !defined(GENERICSTDS_H)
#define GENERICSTDS_H

#include "machine_type.h"

typedef enum {
  SECT_DATA_L1 = 0x2000,
  SECT_DATA_L2,
  SECT_DATA_L1_A,
  SECT_DATA_L1_B,
  SECT_CONSTDATA_L1,

  SECT_DATA_EXTERN = 0x4000,
  SECT_CONSTDATA_EXTERN

} MEMORY_SECTION;

#define H_ALLOC_MEM(name, type) \
  type *Get##name(int n = 0);   \
  void Free##name(type **p);    \
  UINT GetRequiredMem##name(void);

#define H_ALLOC_MEM_OVERLAY(name, type) \
  type *Get##name(int n = 0);           \
  void Free##name(type **p);            \
  UINT GetRequiredMem##name(void);

#define C_ALLOC_MEM(name, type, num)               \
  type *Get##name(int n) {                         \
    FDK_ASSERT((n) == 0);                          \
    return ((type *)FDKcalloc(num, sizeof(type))); \
  }                                                \
  void Free##name(type **p) {                      \
    if (p != NULL) {                               \
      FDKfree(*p);                                 \
      *p = NULL;                                   \
    }                                              \
  }                                                \
  UINT GetRequiredMem##name(void) {                \
    return ALGN_SIZE_EXTRES((num) * sizeof(type)); \
  }

#define C_ALLOC_MEM2(name, type, n1, n2)                 \
  type *Get##name(int n) {                               \
    FDK_ASSERT((n) < (n2));                              \
    return ((type *)FDKcalloc(n1, sizeof(type)));        \
  }                                                      \
  void Free##name(type **p) {                            \
    if (p != NULL) {                                     \
      FDKfree(*p);                                       \
      *p = NULL;                                         \
    }                                                    \
  }                                                      \
  UINT GetRequiredMem##name(void) {                      \
    return ALGN_SIZE_EXTRES((n1) * sizeof(type)) * (n2); \
  }

#define C_AALLOC_MEM(name, type, num)                                  \
  type *Get##name(int n) {                                             \
    type *ap;                                                          \
    FDK_ASSERT((n) == 0);                                              \
    ap = ((type *)FDKaalloc((num) * sizeof(type), ALIGNMENT_DEFAULT)); \
    return ap;                                                         \
  }                                                                    \
  void Free##name(type **p) {                                          \
    if (p != NULL) {                                                   \
      FDKafree(*p);                                                    \
      *p = NULL;                                                       \
    }                                                                  \
  }                                                                    \
  UINT GetRequiredMem##name(void) {                                    \
    return ALGN_SIZE_EXTRES((num) * sizeof(type) + ALIGNMENT_DEFAULT + \
                            sizeof(void *));                           \
  }

#define C_AALLOC_MEM2(name, type, n1, n2)                             \
  type *Get##name(int n) {                                            \
    type *ap;                                                         \
    FDK_ASSERT((n) < (n2));                                           \
    ap = ((type *)FDKaalloc((n1) * sizeof(type), ALIGNMENT_DEFAULT)); \
    return ap;                                                        \
  }                                                                   \
  void Free##name(type **p) {                                         \
    if (p != NULL) {                                                  \
      FDKafree(*p);                                                   \
      *p = NULL;                                                      \
    }                                                                 \
  }                                                                   \
  UINT GetRequiredMem##name(void) {                                   \
    return ALGN_SIZE_EXTRES((n1) * sizeof(type) + ALIGNMENT_DEFAULT + \
                            sizeof(void *)) *                         \
           (n2);                                                      \
  }

#define C_ALLOC_MEM_L(name, type, num, s)               \
  type *Get##name(int n) {                              \
    FDK_ASSERT((n) == 0);                               \
    return ((type *)FDKcalloc_L(num, sizeof(type), s)); \
  }                                                     \
  void Free##name(type **p) {                           \
    if (p != NULL) {                                    \
      FDKfree_L(*p);                                    \
      *p = NULL;                                        \
    }                                                   \
  }                                                     \
  UINT GetRequiredMem##name(void) {                     \
    return ALGN_SIZE_EXTRES((num) * sizeof(type));      \
  }

#define C_ALLOC_MEM2_L(name, type, n1, n2, s)            \
  type *Get##name(int n) {                               \
    FDK_ASSERT((n) < (n2));                              \
    return (type *)FDKcalloc_L(n1, sizeof(type), s);     \
  }                                                      \
  void Free##name(type **p) {                            \
    if (p != NULL) {                                     \
      FDKfree_L(*p);                                     \
      *p = NULL;                                         \
    }                                                    \
  }                                                      \
  UINT GetRequiredMem##name(void) {                      \
    return ALGN_SIZE_EXTRES((n1) * sizeof(type)) * (n2); \
  }

#define C_AALLOC_MEM_L(name, type, num, s)                                  \
  type *Get##name(int n) {                                                  \
    type *ap;                                                               \
    FDK_ASSERT((n) == 0);                                                   \
    ap = ((type *)FDKaalloc_L((num) * sizeof(type), ALIGNMENT_DEFAULT, s)); \
    return ap;                                                              \
  }                                                                         \
  void Free##name(type **p) {                                               \
    if (p != NULL) {                                                        \
      FDKafree_L(*p);                                                       \
      *p = NULL;                                                            \
    }                                                                       \
  }                                                                         \
  UINT GetRequiredMem##name(void) {                                         \
    return ALGN_SIZE_EXTRES((num) * sizeof(type) + ALIGNMENT_DEFAULT +      \
                            sizeof(void *));                                \
  }

#define C_AALLOC_MEM2_L(name, type, n1, n2, s)                             \
  type *Get##name(int n) {                                                 \
    type *ap;                                                              \
    FDK_ASSERT((n) < (n2));                                                \
    ap = ((type *)FDKaalloc_L((n1) * sizeof(type), ALIGNMENT_DEFAULT, s)); \
    return ap;                                                             \
  }                                                                        \
  void Free##name(type **p) {                                              \
    if (p != NULL) {                                                       \
      FDKafree_L(*p);                                                      \
      *p = NULL;                                                           \
    }                                                                      \
  }                                                                        \
  UINT GetRequiredMem##name(void) {                                        \
    return ALGN_SIZE_EXTRES((n1) * sizeof(type) + ALIGNMENT_DEFAULT +      \
                            sizeof(void *)) *                              \
           (n2);                                                           \
  }

#define C_ALLOC_MEM_OVERLAY(name, type, num, sect, tag) \
  C_AALLOC_MEM_L(name, type, num, sect)

#define C_AALLOC_SCRATCH_START(name, type, n)                 \
  type _##name[(n) + (ALIGNMENT_DEFAULT + sizeof(type) - 1)]; \
  type *name = (type *)ALIGN_PTR(_##name);                    \
  C_ALLOC_ALIGNED_REGISTER(name, (n) * sizeof(type));

#define C_ALLOC_SCRATCH_START(name, type, n) type name[n];

#define C_AALLOC_SCRATCH_END(name, type, n) C_ALLOC_ALIGNED_UNREGISTER(name);
#define C_ALLOC_SCRATCH_END(name, type, n)

#define C_AALLOC_STACK_START(name, type, n)                   \
  type _##name[(n) + (ALIGNMENT_DEFAULT + sizeof(type) - 1)]; \
  type *name = (type *)ALIGN_PTR(_##name);                    \
  C_ALLOC_ALIGNED_REGISTER(name, (n) * sizeof(type));

#define C_AALLOC_STACK_END(name, type, n) C_ALLOC_ALIGNED_UNREGISTER(name);

#define C_ALLOC_ALIGNED_REGISTER(x, size)
#define C_ALLOC_ALIGNED_UNREGISTER(x)
#define C_ALLOC_ALIGNED_CHECK(x)
#define C_ALLOC_ALIGNED_CHECK2(x, y)
#define FDK_showBacktrace(a, b)

#define FDK_EXITCODE_OK 0
#define FDK_EXITCODE_USAGE \
  64
#define FDK_EXITCODE_DATAERROR \
  65
#define FDK_EXITCODE_NOINPUT \
  66
#define FDK_EXITCODE_UNAVAILABLE \
  69
#define FDK_EXITCODE_SOFTWARE \
  70
#define FDK_EXITCODE_CANTCREATE \
  73
#define FDK_EXITCODE_IOERROR \
  74

#ifdef __cplusplus
extern "C" {
#endif

void FDKprintf(const char *szFmt, ...);
void FDKprintfErr(const char *szFmt, ...);

int FDKgetchar(void);

INT FDKfprintf(void *stream, const char *format, ...);
INT FDKsprintf(char *str, const char *format, ...);

char *FDKstrchr(char *s, INT c);
const char *FDKstrstr(const char *haystack, const char *needle);
char *FDKstrcpy(char *dest, const char *src);
char *FDKstrncpy(char *dest, const char *src, const UINT n);

#define FDK_MAX_OVERLAYS 8

void *FDKcalloc(const UINT n, const UINT size);
void *FDKmalloc(const UINT size);
void FDKfree(void *ptr);
void *FDKaalloc(const UINT size, const UINT alignment);
void FDKafree(void *ptr);
void *FDKcalloc_L(const UINT n, const UINT size, MEMORY_SECTION s);
void *FDKaalloc_L(const UINT size, const UINT alignment, MEMORY_SECTION s);
void FDKfree_L(void *ptr);
void FDKafree_L(void *ptr);
void FDKmemcpy(void *dst, const void *src, const UINT size);
void FDKmemmove(void *dst, const void *src, const UINT size);
void FDKmemclear(void *memPtr, const UINT size);
void FDKmemset(void *memPtr, const INT value, const UINT size);

INT FDKmemcmp(const void *s1, const void *s2, const UINT size);
INT FDKstrcmp(const char *s1, const char *s2);
INT FDKstrncmp(const char *s1, const char *s2, const UINT size);

UINT FDKstrlen(const char *s);

#define FDKmax(a, b) ((a) > (b) ? (a) : (b))
#define FDKmin(a, b) ((a) < (b) ? (a) : (b))

#define FDK_INT_MAX ((INT)0x7FFFFFFF)
#define FDK_INT_MIN ((INT)0x80000000)

int IS_LITTLE_ENDIAN(void);

UINT TO_LITTLE_ENDIAN(UINT val);

typedef void FDKFILE;
extern const INT FDKSEEK_SET, FDKSEEK_CUR, FDKSEEK_END;

FDKFILE *FDKfopen(const char *filename, const char *mode);
INT FDKfclose(FDKFILE *FP);
INT FDKfseek(FDKFILE *FP, LONG OFFSET, int WHENCE);
INT FDKftell(FDKFILE *FP);
INT FDKfflush(FDKFILE *fp);
UINT FDKfwrite(const void *ptrf, INT size, UINT nmemb, FDKFILE *fp);
UINT FDKfread(void *dst, INT size, UINT nmemb, FDKFILE *fp);
char *FDKfgets(void *dst, INT size, FDKFILE *fp);
void FDKrewind(FDKFILE *fp);
INT FDKfeof(FDKFILE *fp);
UINT FDKfwrite_EL(const void *ptrf, INT size, UINT nmemb, FDKFILE *fp);
UINT FDKfread_EL(void *dst, INT size, UINT nmemb, FDKFILE *fp);

void FDKprintDisclaimer(void);

#ifdef __cplusplus
}
#endif

#endif /* GENERICSTDS_H */
