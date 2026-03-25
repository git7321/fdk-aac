/******************* Library for basic calculation routines ********************

   Author(s):

   Description: matrix memory allocation

*******************************************************************************/

#ifndef FDK_MATRIXCALLOC_H
#define FDK_MATRIXCALLOC_H

#include "machine_type.h"
#include "genericStds.h"

void* fdkCallocMatrix1D(UINT dim1, UINT size);
void* fdkCallocMatrix1D_aligned(UINT dim1, UINT size);
void* fdkCallocMatrix1D_int(UINT dim1, UINT size, MEMORY_SECTION s);
void* fdkCallocMatrix1D_int_aligned(UINT dim1, UINT size, MEMORY_SECTION s);
void fdkFreeMatrix1D(void* p);
void fdkFreeMatrix1D_aligned(void* p);

void** fdkCallocMatrix2D(UINT dim1, UINT dim2, UINT size);
void** fdkCallocMatrix2D_aligned(UINT dim1, UINT dim2, UINT size);
void** fdkCallocMatrix2D_int(UINT dim1, UINT dim2, UINT size, MEMORY_SECTION s);
void** fdkCallocMatrix2D_int_aligned(UINT dim1, UINT dim2, UINT size,
                                     MEMORY_SECTION s);
void fdkFreeMatrix2D(void** p);
void fdkFreeMatrix2D_aligned(void** p);

void*** fdkCallocMatrix3D(UINT dim1, UINT dim2, UINT dim3, UINT size);
void*** fdkCallocMatrix3D_int(UINT dim1, UINT dim2, UINT dim3, UINT size,
                              MEMORY_SECTION s);
void fdkFreeMatrix3D(void*** p);

#define FDK_ALLOCATE_MEMORY_1D(a, dim1, type)                           \
  if (((a) = (type*)fdkCallocMatrix1D((dim1), sizeof(type))) == NULL) { \
    goto bail;                                                          \
  }

#define FDK_ALLOCATE_MEMORY_1D_ALIGNED(a, dim1, type)                   \
  if (((a) = (type*)fdkCallocMatrix1D_aligned((dim1), sizeof(type))) == \
      NULL) {                                                           \
    goto bail;                                                          \
  }

#define FDK_ALLOCATE_MEMORY_1D_P(a, dim1, type, ptype)                  \
  if (((a) = (ptype)fdkCallocMatrix1D((dim1), sizeof(type))) == NULL) { \
    goto bail;                                                          \
  }

#define FDK_ALLOCATE_MEMORY_1D_INT(a, dim1, type, s)                     \
  if (((a) = (type*)fdkCallocMatrix1D_int((dim1), sizeof(type), (s))) == \
      NULL) {                                                            \
    goto bail;                                                           \
  }

#define FDK_FREE_MEMORY_1D(a)    \
  do {                           \
    fdkFreeMatrix1D((void*)(a)); \
    (a) = NULL;                  \
  } while (0)

#define FDK_FREE_MEMORY_1D_ALIGNED(a)    \
  do {                                   \
    fdkFreeMatrix1D_aligned((void*)(a)); \
    (a) = NULL;                          \
  } while (0)

#define FDK_ALLOCATE_MEMORY_2D(a, dim1, dim2, type)                      \
  if (((a) = (type**)fdkCallocMatrix2D((dim1), (dim2), sizeof(type))) == \
      NULL) {                                                            \
    goto bail;                                                           \
  }

#define FDK_ALLOCATE_MEMORY_2D_INT(a, dim1, dim2, type, s)               \
  if (((a) = (type**)fdkCallocMatrix2D_int((dim1), (dim2), sizeof(type), \
                                           (s))) == NULL) {              \
    goto bail;                                                           \
  }

#define FDK_ALLOCATE_MEMORY_2D_INT_ALIGNED(a, dim1, dim2, type, s) \
  if (((a) = (type**)fdkCallocMatrix2D_int_aligned(                \
           (dim1), (dim2), sizeof(type), (s))) == NULL) {          \
    goto bail;                                                     \
  }

#define FDK_FREE_MEMORY_2D(a)     \
  do {                            \
    fdkFreeMatrix2D((void**)(a)); \
    (a) = NULL;                   \
  } while (0)

#define FDK_FREE_MEMORY_2D_ALIGNED(a)     \
  do {                                    \
    fdkFreeMatrix2D_aligned((void**)(a)); \
    (a) = NULL;                           \
  } while (0)

#define FDK_ALLOCATE_MEMORY_3D(a, dim1, dim2, dim3, type)         \
  if (((a) = (type***)fdkCallocMatrix3D((dim1), (dim2), (dim3),   \
                                        sizeof(type))) == NULL) { \
    goto bail;                                                    \
  }

#define FDK_ALLOCATE_MEMORY_3D_INT(a, dim1, dim2, dim3, type, s)           \
  if (((a) = (type***)fdkCallocMatrix3D_int((dim1), (dim2), (dim3),        \
                                            sizeof(type), (s))) == NULL) { \
    goto bail;                                                             \
  }

#define FDK_FREE_MEMORY_3D(a)      \
  do {                             \
    fdkFreeMatrix3D((void***)(a)); \
    (a) = NULL;                    \
  } while (0)

#endif /* FDK_MATRIXCALLOC_H */
