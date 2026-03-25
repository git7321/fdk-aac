/******************* Library for basic calculation routines ********************

   Author(s):   Omer Osman

   Description: MPEG-D SAC/USAC/SAOC Huffman Part0 Tables

*******************************************************************************/

#ifndef HUFF_NODES_H
#define HUFF_NODES_H

#include "genericStds.h"

typedef struct {
  SHORT nodeTab[39][2];

} HUFF_RES_NODES;

typedef struct {
  SHORT nodeTab[30][2];

} HUFF_CLD_NOD_1D;

typedef struct {
  SHORT nodeTab[7][2];

} HUFF_ICC_NOD_1D;

typedef struct {
  SHORT nodeTab[50][2];

} HUFF_CPC_NOD_1D;

typedef struct {
  SHORT nodeTab[15][2];

} HUFF_OLD_NOD_1D;

typedef struct {
  SHORT nodeTab[63][2];

} HUFF_NRG_NOD_1D;

typedef struct {
  SHORT lav3[15][2];
  SHORT lav5[35][2];
  SHORT lav7[63][2];
  SHORT lav9[99][2];

} HUFF_CLD_NOD_2D;

typedef struct {
  SHORT lav1[3][2];
  SHORT lav3[15][2];
  SHORT lav5[35][2];
  SHORT lav7[63][2];

} HUFF_ICC_NOD_2D;

typedef struct {
  SHORT lav3[15][2];
  SHORT lav6[48][2];
  SHORT lav9[99][2];
  SHORT lav12[168][2];

} HUFF_OLD_NOD_2D;

typedef struct {
  SHORT lav3[15][2];
  SHORT lav5[35][2];
  SHORT lav7[63][2];
  SHORT lav9[99][2];

} HUFF_NRG_NOD_2D_df;

typedef struct {
  SHORT lav3[15][2];
  SHORT lav6[48][2];
  SHORT lav9[99][2];
  SHORT lav12[168][2];

} HUFF_NRG_NOD_2D_dt;

typedef struct {
  HUFF_NRG_NOD_2D_df df[2];
  HUFF_NRG_NOD_2D_dt dt[2];
  HUFF_NRG_NOD_2D_df dp[2];

} HUFF_NRG_NOD_2D;

typedef struct {
  const HUFF_CLD_NOD_1D *h1D[3];
  const HUFF_CLD_NOD_2D *h2D[3][2];

} HUFF_CLD_NODES;

typedef struct {
  const HUFF_ICC_NOD_1D *h1D[3];
  const HUFF_ICC_NOD_2D *h2D[3][2];

} HUFF_ICC_NODES;

typedef struct {
  const HUFF_OLD_NOD_1D *h1D[3];
  const HUFF_OLD_NOD_2D *h2D[3][2];

} HUFF_OLD_NODES;

typedef struct {
  const HUFF_NRG_NOD_1D *h1D[3];
  const HUFF_NRG_NOD_2D *h2D;

} HUFF_NRG_NODES;

typedef struct {
  SHORT cld[30][2];
  SHORT icc[7][2];
  SHORT ipd[7][2];
  SHORT old[15][2];
  SHORT nrg[63][2];
} HUFF_PT0_NODES;

typedef struct {
  SHORT nodeTab[3][2];

} HUFF_LAV_NODES;

typedef struct {
  SHORT nodeTab[7][2];

} HUFF_IPD_NOD_1D;

typedef struct {
  SHORT lav1[3][2];
  SHORT lav3[15][2];
  SHORT lav5[35][2];
  SHORT lav7[63][2];

} HUFF_IPD_NOD_2D;

typedef struct {
  HUFF_IPD_NOD_1D h1D[3];
  HUFF_IPD_NOD_2D h2D[3][2];

} HUFF_IPD_NODES;

extern const HUFF_PT0_NODES FDK_huffPart0Nodes;
extern const HUFF_LAV_NODES FDK_huffLavIdxNodes;

extern const HUFF_ICC_NODES FDK_huffICCNodes;
extern const HUFF_CLD_NODES FDK_huffCLDNodes;
extern const HUFF_RES_NODES FDK_huffReshapeNodes;

extern const HUFF_OLD_NODES huffOLDNodes;
extern const HUFF_IPD_NODES FDK_huffIPDNodes;

#endif /* HUFF_NODES_H */
