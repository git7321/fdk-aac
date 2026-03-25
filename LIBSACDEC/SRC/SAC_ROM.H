/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: SAC Dec tables

*******************************************************************************/

#ifndef SAC_ROM_H
#define SAC_ROM_H

#include "FDK_archdef.h"
#include "sac_dec_interface.h"
#include "huff_nodes.h"
#include "sac_bitdec.h"
#include "machine_type.h"

#define FIXP_CFG FIXP_DBL
#define FX_CFG2FX_DBL
#define FX_CFG2FX_SGL FX_DBL2FX_SGL
#define CFG(a) FIXP_DBL(a)
#define FL2FXCONST_CFG FL2FXCONST_DBL
#define FX_DBL2FX_CFG(x) ((FIXP_DBL)(x))

#define SCALE_INV_ICC (2)
#define G_dd_SCALE (2)

#define QCC_SCALE 1
#define M1M2_DATA FIXP_DBL
#define M1M2_CDATA FIXP_DBL
#define M1M2_CDATA2FX_DBL(a) (a)
#define FX_DBL2M1M2_CDATA(a) (a)

#define CLIP_PROTECT_GAIN_0(x) FL2FXCONST_CFG(((x) / (float)(1 << 0)))
#define CLIP_PROTECT_GAIN_1(x) FL2FXCONST_CFG(((x) / (float)(1 << 1)))
#define CLIP_PROTECT_GAIN_2(x) FL2FXCONST_CFG(((x) / (float)(1 << 2)))

#define SF_CLD_C1C2 (9)

extern const FIXP_CFG dequantCPC__FDK[];
extern const FIXP_CFG dequantICC__FDK[8];
extern const FIXP_CFG dequantCLD__FDK[31];

#define IPD_SCALE (5)
#define PI__IPD (FL2FXCONST_DBL(3.1415926535897932f / (float)(1 << IPD_SCALE)))
#define PIx2__IPD \
  (FL2FXCONST_DBL(3.1415926535897932f / (float)(1 << (IPD_SCALE - 1))))

extern const FIXP_CFG dequantIPD__FDK[16];
extern const FIXP_DBL dequantIPD_CLD_ICC_splitAngle__FDK[15][31][8];

extern const FIXP_CFG H11_nc[31][8];
extern const FIXP_CFG H12_nc[31][8];

extern const FIXP_DBL dequantCLD_c1[31];

extern const FIXP_CFG BP__FDK[];
extern const FIXP_CFG BP_GF__FDK[];
extern const SCHAR row2channelSTP[][MAX_M2_INPUT];

extern const INT samplingFreqTable[16];
extern const UCHAR freqResTable[];
extern const UCHAR freqResTable_LD[];
extern const UCHAR tempShapeChanTable[2][8];
extern const TREEPROPERTIES treePropertyTable[];

extern const SCHAR kernels_4_to_71[MAX_HYBRID_BANDS];
extern const SCHAR kernels_5_to_71[MAX_HYBRID_BANDS];
extern const SCHAR kernels_7_to_71[MAX_HYBRID_BANDS];
extern const SCHAR kernels_10_to_71[MAX_HYBRID_BANDS];
extern const SCHAR kernels_14_to_71[MAX_HYBRID_BANDS];
extern const SCHAR kernels_20_to_71[MAX_HYBRID_BANDS];
extern const SCHAR kernels_28_to_71[MAX_HYBRID_BANDS];

extern const UCHAR mapping_4_to_28[MAX_PARAMETER_BANDS];
extern const UCHAR mapping_5_to_28[MAX_PARAMETER_BANDS];
extern const UCHAR mapping_7_to_28[MAX_PARAMETER_BANDS];
extern const UCHAR mapping_10_to_28[MAX_PARAMETER_BANDS];
extern const UCHAR mapping_14_to_28[MAX_PARAMETER_BANDS];
extern const UCHAR mapping_20_to_28[MAX_PARAMETER_BANDS];
extern const SCHAR kernels_4_to_64[MAX_HYBRID_BANDS];
extern const SCHAR kernels_5_to_64[MAX_HYBRID_BANDS];
extern const SCHAR kernels_7_to_64[MAX_HYBRID_BANDS];
extern const SCHAR kernels_9_to_64[MAX_HYBRID_BANDS];
extern const SCHAR kernels_12_to_64[MAX_HYBRID_BANDS];
extern const SCHAR kernels_15_to_64[MAX_HYBRID_BANDS];
extern const SCHAR kernels_23_to_64[MAX_HYBRID_BANDS];

extern const UCHAR mapping_15_to_23[MAX_PARAMETER_BANDS_LD];
extern const UCHAR mapping_12_to_23[MAX_PARAMETER_BANDS_LD];
extern const UCHAR mapping_9_to_23[MAX_PARAMETER_BANDS_LD];
extern const UCHAR mapping_7_to_23[MAX_PARAMETER_BANDS_LD];
extern const UCHAR mapping_5_to_23[MAX_PARAMETER_BANDS_LD];
extern const UCHAR mapping_4_to_23[MAX_PARAMETER_BANDS_LD];

extern const FIXP_CFG clipGainTable__FDK[];
extern const UCHAR clipGainSFTable__FDK[];

extern const UCHAR pbStrideTable[];
extern const int smgTimeTable[];

extern const FIXP_CFG envShapeDataTable__FDK[5][2];
extern const SCHAR row2channelGES[][MAX_M2_INPUT];

extern const SCHAR row2residual[][MAX_M2_INPUT];

void SpatialDequantGetCLDValues(int index, FIXP_DBL* cu, FIXP_DBL* cl);
void SpatialDequantGetCLD2Values(int index, FIXP_DBL* x);

static inline int SacGetHybridSubbands(int qmfSubbands) {
  return qmfSubbands - MAX_QMF_BANDS_TO_HYBRID + 10;
}

#endif /* SAC_ROM_H */
