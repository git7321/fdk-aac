#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#define __GENERICSTDS_CPP__

#include "genericStds.h"

#define SYS_LIB_VL0 2
#define SYS_LIB_VL1 0
#define SYS_LIB_VL2 0
#define SYS_LIB_TITLE "System Integration Library"
#define SYS_LIB_BUILD_DATE ""
#define SYS_LIB_BUILD_TIME ""

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void FDKprintf(const char *szFmt, ...) {
  va_list ap;
  va_start(ap, szFmt);
  vprintf(szFmt, ap);
  va_end(ap);
}

void FDKprintfErr(const char *szFmt, ...) {
  va_list ap;
  va_start(ap, szFmt);
  vfprintf(stderr, szFmt, ap);
  va_end(ap);
}

int FDKgetchar(void) { return getchar(); }

INT FDKfprintf(FDKFILE *stream, const char *format, ...) {
  INT chars = 0;
  va_list ap;
  va_start(ap, format);
  chars += vfprintf((FILE *)stream, format, ap);
  va_end(ap);
  return chars;
}

INT FDKsprintf(char *str, const char *format, ...) {
  INT chars = 0;
  va_list ap;
  va_start(ap, format);
  chars += vsprintf(str, format, ap);
  va_end(ap);
  return chars;
}

char *FDKstrchr(char *s, INT c) { return strchr(s, c); }
const char *FDKstrstr(const char *haystack, const char *needle) {
  return strstr(haystack, needle);
}
char *FDKstrcpy(char *dest, const char *src) { return strcpy(dest, src); }
char *FDKstrncpy(char *dest, const char *src, UINT n) {
  return strncpy(dest, src, n);
}

void *FDKcalloc(const UINT n, const UINT size) {
  void *ptr;

  ptr = calloc(n, size);

  return ptr;
}

void *FDKmalloc(const UINT size) {
  void *ptr;

  ptr = malloc(size);

  return ptr;
}

void FDKfree(void *ptr) { free((INT *)ptr); }

void *FDKaalloc(const UINT size, const UINT alignment) {
  void *addr, *result = NULL;
  addr = FDKcalloc(1, size + alignment +
                          (UINT)sizeof(void *));

  if (addr != NULL) {
    result = ALIGN_PTR((unsigned char *)addr +
                       sizeof(void *));
    *(((void **)result) - 1) = addr;
    C_ALLOC_ALIGNED_REGISTER(result, size);
  }

  return result;
}

void FDKafree(void *ptr) {
  void *addr;
  addr = *(((void **)ptr) - 1);

  C_ALLOC_ALIGNED_UNREGISTER(ptr);

  FDKfree(addr);
}

void *FDKcalloc_L(const UINT dim, const UINT size, MEMORY_SECTION s) {
  return FDKcalloc(dim, size);
}

void FDKfree_L(void *p) { FDKfree(p); }

void *FDKaalloc_L(const UINT size, const UINT alignment, MEMORY_SECTION s) {
  void *addr, *result = NULL;
  addr = FDKcalloc_L(1, size + alignment + (UINT)sizeof(void *),
                     s);

  if (addr != NULL) {
    result = ALIGN_PTR((unsigned char *)addr +
                       sizeof(void *));
    *(((void **)result) - 1) = addr;
    C_ALLOC_ALIGNED_REGISTER(result, size);
  }

  return result;
}

void FDKafree_L(void *ptr) {
  void *addr;

  addr = *(((void **)ptr) - 1);

  C_ALLOC_ALIGNED_UNREGISTER(ptr);

  FDKfree_L(addr);
}

void FDKmemcpy(void *dst, const void *src, const UINT size) {
  FDK_ASSERT(((const unsigned char *)dst - (const unsigned char *)src) >=
                 (ptrdiff_t)size ||
             ((const unsigned char *)src - (const unsigned char *)dst) >=
                 (ptrdiff_t)size);

  memcpy(dst, src, size);
}

void FDKmemmove(void *dst, const void *src, const UINT size) {
  memmove(dst, src, size);
}

void FDKmemset(void *memPtr, const INT value, const UINT size) {
  memset(memPtr, value, size);
}

void FDKmemclear(void *memPtr, const UINT size) { FDKmemset(memPtr, 0, size); }

UINT FDKstrlen(const char *s) { return (UINT)strlen(s); }

INT FDKmemcmp(const void *s1, const void *s2, const UINT size) {
  return memcmp(s1, s2, size);
}
INT FDKstrcmp(const char *s1, const char *s2) { return strcmp(s1, s2); }
INT FDKstrncmp(const char *s1, const char *s2, const UINT size) {
  return strncmp(s1, s2, size);
}

int IS_LITTLE_ENDIAN(void) {
  int __dummy = 1;
  return (*((UCHAR *)(&(__dummy))));
}

UINT TO_LITTLE_ENDIAN(UINT val) {
  return IS_LITTLE_ENDIAN()
             ? val
             : (((val & 0xff) << 24) | ((val & 0xff00) << 8) |
                ((val & 0xff0000) >> 8) | ((val & 0xff000000) >> 24));
}

FDKFILE *FDKfopen(const char *filename, const char *mode) {
  return fopen(filename, mode);
}
INT FDKfclose(FDKFILE *fp) { return fclose((FILE *)fp); }
INT FDKfseek(FDKFILE *fp, LONG OFFSET, int WHENCE) {
  return fseek((FILE *)fp, OFFSET, WHENCE);
}
INT FDKftell(FDKFILE *fp) { return ftell((FILE *)fp); }
INT FDKfflush(FDKFILE *fp) { return fflush((FILE *)fp); }
const INT FDKSEEK_SET = SEEK_SET;
const INT FDKSEEK_CUR = SEEK_CUR;
const INT FDKSEEK_END = SEEK_END;

UINT FDKfwrite(const void *ptrf, INT size, UINT nmemb, FDKFILE *fp) {
  return (UINT)fwrite(ptrf, size, nmemb, (FILE *)fp);
}
UINT FDKfread(void *dst, INT size, UINT nmemb, FDKFILE *fp) {
  return (UINT)fread(dst, size, nmemb, (FILE *)fp);
}
char *FDKfgets(void *dst, INT size, FDKFILE *fp) {
  return fgets((char *)dst, size, (FILE *)fp);
}
void FDKrewind(FDKFILE *fp) { FDKfseek((FILE *)fp, 0, FDKSEEK_SET); }

UINT FDKfwrite_EL(const void *ptrf, INT size, UINT nmemb, FDKFILE *fp) {
  if (IS_LITTLE_ENDIAN()) {
    FDKfwrite(ptrf, size, nmemb, fp);
  } else {
    UINT n;
    INT s;

    const UCHAR *ptr = (const UCHAR *)ptrf;

    for (n = 0; n < nmemb; n++) {
      for (s = size - 1; s >= 0; s--) {
        FDKfwrite(ptr + s, 1, 1, fp);
      }
      ptr = ptr + size;
    }
  }
  return nmemb;
}

UINT FDKfread_EL(void *dst, INT size, UINT nmemb, FDKFILE *fp) {
  UINT n, s0, s1, err;
  UCHAR tmp, *ptr;
  UCHAR tmp24[3];

  if (size == 3) {
    ptr = (UCHAR *)dst;
    for (n = 0; n < nmemb; n++) {
      if ((err = FDKfread(tmp24, 1, 3, fp)) != 3) {
        return err;
      }
      *ptr++ = tmp24[0];
      *ptr++ = tmp24[1];
      *ptr++ = tmp24[2];
      if (tmp24[2] & 0x80) {
        *ptr++ = 0xff;
      } else {
        *ptr++ = 0;
      }
    }
    err = nmemb;
    size = sizeof(LONG);
  } else {
    if ((err = FDKfread(dst, size, nmemb, fp)) != nmemb) {
      return err;
    }
  }
  if (!IS_LITTLE_ENDIAN() && size > 1) {
    ptr = (UCHAR *)dst;
    for (n = 0; n < nmemb; n++) {
      for (s0 = 0, s1 = size - 1; s0 < s1; s0++, s1--) {
        tmp = ptr[s0];
        ptr[s0] = ptr[s1];
        ptr[s1] = tmp;
      }
      ptr += size;
    }
  }
  return err;
}

INT FDKfeof(FDKFILE *fp) { return feof((FILE *)fp); }

void FDKprintDisclaimer(void) {
  FDKprintf(
      "This program is protected by copyright law and international treaties.\n"
      "Any reproduction or distribution of this program, or any portion\n"
      "of it, may result in severe civil and criminal penalties, and will be\n"
      "prosecuted to the maximum extent possible under law.\n\n");
}
