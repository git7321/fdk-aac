/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: Definition of constant tables

*******************************************************************************/

#ifndef AAC_ROM_H
#define AAC_ROM_H

#include "common_fix.h"
#include "FDK_audio.h"
#include "aacdec_hcr_types.h"
#include "aacdec_hcrs.h"

#define PCM_AAC LONG
#define PCM_DEC FIXP_DBL
#define MAXVAL_PCM_DEC MAXVAL_DBL
#define MINVAL_PCM_DEC MINVAL_DBL
#define FIXP_DBL2PCM_DEC(x) (x)
#define PCM_DEC2FIXP_DBL(x) (x)
#define PCM_DEC_BITS DFRACT_BITS
#define PCM_DEC2FX_PCM(x) FX_DBL2FX_PCM(x)
#define FX_PCM2PCM_DEC(x) FX_PCM2FX_DBL(x)

#define AACDEC_MAX_CH_CONF 14
#define AACDEC_CH_ELEMENTS_TAB_SIZE 7

#define AAC_NF_NO_RANDOM_VAL \
  512

#define INV_QUANT_TABLESIZE 256

extern const FIXP_DBL InverseQuantTable[INV_QUANT_TABLESIZE + 1];
extern const FIXP_DBL MantissaTable[4][14];
extern const SCHAR ExponentTable[4][14];

#define NUM_LD_COEF_512 1536
#define NUM_LD_COEF_480 1440
#define WTS0 (1)
#define WTS1 (0)
#define WTS2 (-2)
extern const FIXP_WTB LowDelaySynthesis512[1536];
extern const FIXP_WTB LowDelaySynthesis480[1440];
extern const FIXP_WTB LowDelaySynthesis256[768];
extern const FIXP_WTB LowDelaySynthesis240[720];
extern const FIXP_WTB LowDelaySynthesis160[480];
extern const FIXP_WTB LowDelaySynthesis128[384];
extern const FIXP_WTB LowDelaySynthesis120[360];

typedef struct {
  const SHORT *sfbOffsetLong;
  const SHORT *sfbOffsetShort;
  UCHAR numberOfSfbLong;
  UCHAR numberOfSfbShort;
} SFB_INFO;

extern const SFB_INFO sfbOffsetTables[5][16];

enum { HuffmanBits = 2, HuffmanEntries = (1 << HuffmanBits) };

typedef struct {
  const USHORT (*CodeBook)[HuffmanEntries];
  UCHAR Dimension;
  UCHAR numBits;
  UCHAR Offset;
} CodeBookDescription;

extern const CodeBookDescription AACcodeBookDescriptionTable[13];
extern const CodeBookDescription AACcodeBookDescriptionSCL;

extern const STATEFUNC aStateConstant2State[];

extern const SCHAR aCodebook2StartInt[];

extern const UCHAR aMinOfCbPair[];
extern const UCHAR aMaxOfCbPair[];

extern const UCHAR aMaxCwLen[];
extern const UCHAR aDimCb[];
extern const UCHAR aDimCbShift[];
extern const UCHAR aSignCb[];
extern const UCHAR aCbPriority[];

extern const UINT *aHuffTable[];
extern const SCHAR *aQuantTable[];

extern const USHORT aLargestAbsoluteValue[];

extern const UINT aHuffTreeRvlcEscape[];
extern const UINT aHuffTreeRvlCodewds[];

extern const UCHAR tns_max_bands_tbl[13][2];

extern const UCHAR tns_max_bands_tbl_480[13];
extern const UCHAR tns_max_bands_tbl_512[13];

#define FIXP_TCC FIXP_DBL

extern const FIXP_TCC FDKaacDec_tnsCoeff3[8];
extern const FIXP_TCC FDKaacDec_tnsCoeff4[16];

extern const UCHAR FDKaacDec_tnsCoeff3_gain_ld[];
extern const UCHAR FDKaacDec_tnsCoeff4_gain_ld[];

extern const USHORT AacDec_randomSign[AAC_NF_NO_RANDOM_VAL / 16];

extern const FIXP_DBL pow2_div24minus1[47];
extern const int offsetTab[2][16];

extern const UCHAR channelMappingTablePassthrough[15][8];
extern const UCHAR channelMappingTableWAV[15][8];

extern const MP4_ELEMENT_ID elementsTab[AACDEC_MAX_CH_CONF]
                                       [AACDEC_CH_ELEMENTS_TAB_SIZE];

#define SF_FNA_COEFFS \
  1

#ifdef CPLX_PRED_FILTER_16BIT
#define FIXP_FILT FIXP_SGL
#define FILT(a) ((FL2FXCONST_SGL(a)) >> SF_FNA_COEFFS)
#else
#define FIXP_FILT FIXP_DBL
#define FILT(a) ((FL2FXCONST_DBL(a)) >> SF_FNA_COEFFS)
#endif

extern const FIXP_FILT mdst_filt_coef_curr[20][3];
extern const FIXP_FILT mdst_filt_coef_prev[6][4];

#endif /* #ifndef AAC_ROM_H */
