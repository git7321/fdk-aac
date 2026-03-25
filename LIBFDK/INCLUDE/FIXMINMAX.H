/******************* Library for basic calculation routines ********************

   Author(s):   M. Lohwasser, M. Gayer

   Description: min/max inline functions and defines

*******************************************************************************/

#ifndef FIXMINMAX_H
#define FIXMINMAX_H

#include "FDK_archdef.h"
#include "machine_type.h"

template <class T>
inline T fixmin(T a, T b) {
  return (a < b ? a : b);
}

template <class T>
inline T fixmax(T a, T b) {
  return (a > b ? a : b);
}

#define fixmax_D(a, b) fixmax(a, b)
#define fixmin_D(a, b) fixmin(a, b)
#define fixmax_S(a, b) fixmax(a, b)
#define fixmin_S(a, b) fixmin(a, b)
#define fixmax_I(a, b) fixmax(a, b)
#define fixmin_I(a, b) fixmin(a, b)
#define fixmax_UI(a, b) fixmax(a, b)
#define fixmin_UI(a, b) fixmin(a, b)

#endif
