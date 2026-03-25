/**************************** AAC encoder library ******************************

   Author(s):   M. Lohwasser, M. Gayer

   Description:

*******************************************************************************/

#ifndef AACENC_ROM_H
#define AACENC_ROM_H

#include "common_fix.h"

#include "psy_const.h"
#include "psy_configuration.h"
#include "FDK_tools_rom.h"
#include "FDK_lpc.h"

extern const ULONG FDKaacEnc_huff_ltab1_2[3][3][3][3];
extern const ULONG FDKaacEnc_huff_ltab3_4[3][3][3][3];
extern const ULONG FDKaacEnc_huff_ltab5_6[9][9];
extern const ULONG FDKaacEnc_huff_ltab7_8[8][8];
extern const ULONG FDKaacEnc_huff_ltab9_10[13][13];
extern const UCHAR FDKaacEnc_huff_ltab11[17][17];
extern const UCHAR FDKaacEnc_huff_ltabscf[121];
extern const USHORT FDKaacEnc_huff_ctab1[3][3][3][3];
extern const USHORT FDKaacEnc_huff_ctab2[3][3][3][3];
extern const USHORT FDKaacEnc_huff_ctab3[3][3][3][3];
extern const USHORT FDKaacEnc_huff_ctab4[3][3][3][3];
extern const USHORT FDKaacEnc_huff_ctab5[9][9];
extern const USHORT FDKaacEnc_huff_ctab6[9][9];
extern const USHORT FDKaacEnc_huff_ctab7[8][8];
extern const USHORT FDKaacEnc_huff_ctab8[8][8];
extern const USHORT FDKaacEnc_huff_ctab9[13][13];
extern const USHORT FDKaacEnc_huff_ctab10[13][13];
extern const USHORT FDKaacEnc_huff_ctab11[21][17];
extern const ULONG FDKaacEnc_huff_ctabscf[121];

#define MANT_DIGITS 9
#define MANT_SIZE (1 << MANT_DIGITS)

#if defined(ARCH_PREFER_MULT_32x16)
#define FIXP_QTD FIXP_SGL
#define QTC FX_DBL2FXCONST_SGL
#else
#define FIXP_QTD FIXP_DBL
#define QTC
#endif

extern const FIXP_QTD FDKaacEnc_mTab_3_4[MANT_SIZE];
extern const FIXP_QTD FDKaacEnc_quantTableQ[4];
extern const FIXP_QTD FDKaacEnc_quantTableE[4];

extern const FIXP_DBL FDKaacEnc_mTab_4_3Elc[512];
extern const FIXP_DBL FDKaacEnc_specExpMantTableCombElc[4][14];
extern const UCHAR FDKaacEnc_specExpTableComb[4][14];

extern const SHORT FDKaacEnc_sideInfoTabLong[];
extern const SHORT FDKaacEnc_sideInfoTabShort[];

extern const SFB_PARAM_LONG p_FDKaacEnc_8000_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_8000_short_128;
extern const SFB_PARAM_LONG p_FDKaacEnc_11025_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_11025_short_128;
extern const SFB_PARAM_LONG p_FDKaacEnc_12000_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_12000_short_128;
extern const SFB_PARAM_LONG p_FDKaacEnc_16000_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_16000_short_128;
extern const SFB_PARAM_LONG p_FDKaacEnc_22050_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_22050_short_128;
extern const SFB_PARAM_LONG p_FDKaacEnc_24000_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_24000_short_128;
extern const SFB_PARAM_LONG p_FDKaacEnc_32000_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_32000_short_128;
extern const SFB_PARAM_LONG p_FDKaacEnc_44100_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_44100_short_128;
extern const SFB_PARAM_LONG p_FDKaacEnc_48000_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_48000_short_128;
extern const SFB_PARAM_LONG p_FDKaacEnc_64000_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_64000_short_128;
extern const SFB_PARAM_LONG p_FDKaacEnc_88200_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_88200_short_128;
extern const SFB_PARAM_LONG p_FDKaacEnc_96000_long_1024;
extern const SFB_PARAM_SHORT p_FDKaacEnc_96000_short_128;

extern const FIXP_LPC FDKaacEnc_tnsEncCoeff3[8];
extern const FIXP_LPC FDKaacEnc_tnsCoeff3Borders[8];
extern const FIXP_LPC FDKaacEnc_tnsEncCoeff4[16];
extern const FIXP_LPC FDKaacEnc_tnsCoeff4Borders[16];

#define WTC0 WTC
#define WTC1 WTC
#define WTC2 WTC

extern const FIXP_WTB ELDAnalysis512[1536];
extern const FIXP_WTB ELDAnalysis480[1440];
extern const FIXP_WTB ELDAnalysis256[768];
extern const FIXP_WTB ELDAnalysis240[720];
extern const FIXP_WTB ELDAnalysis128[384];
extern const FIXP_WTB ELDAnalysis120[360];

#endif /* #ifndef AACENC_ROM_H */
