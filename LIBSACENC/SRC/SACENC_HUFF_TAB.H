/*********************** MPEG surround encoder library *************************

   Author(s):   Markus Lohwasser

   Description: SAC-Encoder constant huffman tables

*******************************************************************************/

#ifndef SACENC_HUFF_TAB_H
#define SACENC_HUFF_TAB_H

#include "machine_type.h"

#define HUFF_PACK(a, b) \
  { \
    ((((ULONG)a) & 0x00FFFFFF) << 8) | (((ULONG)b) & 0xFF) \
  }
#define HUFF_VALUE(a) \
  (((a.packed >> 8) & 0x00FFFFFF))
#define HUFF_LENGTH(a) \
  ((a.packed & 0xFF))

typedef struct {
  ULONG packed;
} HUFF_ENTRY;

typedef struct {
  HUFF_ENTRY entry[2][2];
  HUFF_ENTRY escape;

} LAV1_2D;

typedef struct {
  HUFF_ENTRY entry[4][4];
  HUFF_ENTRY escape;

} LAV3_2D;

typedef struct {
  HUFF_ENTRY entry[6][6];
  HUFF_ENTRY escape;

} LAV5_2D;

typedef struct {
  HUFF_ENTRY entry[7][7];
  HUFF_ENTRY escape;

} LAV6_2D;

typedef struct {
  HUFF_ENTRY entry[8][8];
  HUFF_ENTRY escape;

} LAV7_2D;

typedef struct {
  HUFF_ENTRY entry[10][10];
  HUFF_ENTRY escape;

} LAV9_2D;

typedef struct {
  HUFF_ENTRY entry[13][13];
  HUFF_ENTRY escape;

} LAV12_2D;

typedef struct {
  LAV3_2D lav3;
  LAV5_2D lav5;
  LAV7_2D lav7;
  LAV9_2D lav9;

} HUFF_CLD_TAB_2D;

typedef struct {
  LAV1_2D lav1;
  LAV3_2D lav3;
  LAV5_2D lav5;
  LAV7_2D lav7;

} HUFF_ICC_TAB_2D;

typedef struct {
  HUFF_ENTRY h1D[2][31];
  HUFF_CLD_TAB_2D h2D[2][2];

} HUFF_CLD_TABLE;

typedef struct {
  HUFF_ENTRY h1D[2][8];
  HUFF_ICC_TAB_2D h2D[2][2];

} HUFF_ICC_TABLE;

typedef struct {
  HUFF_ENTRY cld[31];
  HUFF_ENTRY icc[8];

} HUFF_PT0_TABLE;

typedef HUFF_ENTRY HUFF_RES_TABLE[5][8];

extern const HUFF_CLD_TABLE fdk_sacenc_huffCLDTab;
extern const HUFF_ICC_TABLE fdk_sacenc_huffICCTab;
extern const HUFF_PT0_TABLE fdk_sacenc_huffPart0Tab;

#endif /* SACENC_HUFF_TAB_H */
