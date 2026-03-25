/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef SBR_ROM_H
#define SBR_ROM_H

#include "sbrdecoder.h"
#include "env_extr.h"
#include "qmf.h"

#define INV_INT_TABLE_SIZE 49
#define SBR_NF_NO_RANDOM_VAL \
  512

extern const UCHAR FDK_sbrDecoder_sbr_start_freq_16[(4) / 2][16];
extern const UCHAR FDK_sbrDecoder_sbr_start_freq_22[(4) / 2][16];
extern const UCHAR FDK_sbrDecoder_sbr_start_freq_24[(4) / 2][16];
extern const UCHAR FDK_sbrDecoder_sbr_start_freq_32[(4) / 2][16];
extern const UCHAR FDK_sbrDecoder_sbr_start_freq_40[(4) / 2][16];
extern const UCHAR FDK_sbrDecoder_sbr_start_freq_44[(4) / 2][16];
extern const UCHAR FDK_sbrDecoder_sbr_start_freq_48[(4) / 2][16];
extern const UCHAR FDK_sbrDecoder_sbr_start_freq_64[(4) / 2][16];
extern const UCHAR FDK_sbrDecoder_sbr_start_freq_88[(4) / 2][16];
extern const UCHAR FDK_sbrDecoder_sbr_start_freq_192[16];
extern const UCHAR FDK_sbrDecoder_sbr_start_freq_176[16];
extern const UCHAR FDK_sbrDecoder_sbr_start_freq_128[16];

#define NUM_WHFACTOR_TABLE_ENTRIES 9
extern const USHORT
    FDK_sbrDecoder_sbr_whFactorsIndex[NUM_WHFACTOR_TABLE_ENTRIES];
extern const FIXP_DBL
    FDK_sbrDecoder_sbr_whFactorsTable[NUM_WHFACTOR_TABLE_ENTRIES][6];

extern const FIXP_SGL FDK_sbrDecoder_sbr_limGains_m[4];
extern const UCHAR FDK_sbrDecoder_sbr_limGains_e[4];
extern const FIXP_SGL FDK_sbrDecoder_sbr_limGainsPvc_m[4];
extern const UCHAR FDK_sbrDecoder_sbr_limGainsPvc_e[4];
extern const FIXP_SGL FDK_sbrDecoder_sbr_limiterBandsPerOctaveDiv4[4];
extern const FIXP_DBL FDK_sbrDecoder_sbr_limiterBandsPerOctaveDiv4_DBL[4];
extern const FIXP_SGL FDK_sbrDecoder_sbr_smoothFilter[4];
extern const FIXP_SGL FDK_sbrDecoder_sbr_randomPhase[SBR_NF_NO_RANDOM_VAL][2];

extern const int FDK_sbrDecoder_envelopeTable_8[8][5];
extern const int FDK_sbrDecoder_envelopeTable_15[15][6];
extern const int FDK_sbrDecoder_envelopeTable_16[16][6];
extern const FRAME_INFO FDK_sbrDecoder_sbr_frame_info1_15;
extern const FRAME_INFO FDK_sbrDecoder_sbr_frame_info2_15;
extern const FRAME_INFO FDK_sbrDecoder_sbr_frame_info4_15;
extern const FRAME_INFO FDK_sbrDecoder_sbr_frame_info8_15;
extern const FRAME_INFO FDK_sbrDecoder_sbr_frame_info1_16;
extern const FRAME_INFO FDK_sbrDecoder_sbr_frame_info2_16;
extern const FRAME_INFO FDK_sbrDecoder_sbr_frame_info4_16;
extern const FRAME_INFO FDK_sbrDecoder_sbr_frame_info8_16;

extern const SCHAR FDK_sbrDecoder_sbr_huffBook_EnvLevel10T[120][2];
extern const SCHAR FDK_sbrDecoder_sbr_huffBook_EnvLevel10F[120][2];
extern const SCHAR FDK_sbrDecoder_sbr_huffBook_EnvBalance10T[48][2];
extern const SCHAR FDK_sbrDecoder_sbr_huffBook_EnvBalance10F[48][2];
extern const SCHAR FDK_sbrDecoder_sbr_huffBook_EnvLevel11T[62][2];
extern const SCHAR FDK_sbrDecoder_sbr_huffBook_EnvLevel11F[62][2];
extern const SCHAR FDK_sbrDecoder_sbr_huffBook_EnvBalance11T[24][2];
extern const SCHAR FDK_sbrDecoder_sbr_huffBook_EnvBalance11F[24][2];
extern const SCHAR FDK_sbrDecoder_sbr_huffBook_NoiseLevel11T[62][2];
extern const SCHAR FDK_sbrDecoder_sbr_huffBook_NoiseBalance11T[24][2];

extern const UCHAR FDK_sbrDecoder_aFixNoEnvDecode[4];

extern const SCHAR aBookPsIidTimeDecode[28][2];
extern const SCHAR aBookPsIidFreqDecode[28][2];
extern const SCHAR aBookPsIccTimeDecode[14][2];
extern const SCHAR aBookPsIccFreqDecode[14][2];

extern const SCHAR aBookPsIidFineTimeDecode[60][2];
extern const SCHAR aBookPsIidFineFreqDecode[60][2];

extern const FIXP_DBL ScaleFactors[NO_IID_LEVELS];
extern const FIXP_DBL ScaleFactorsFine[NO_IID_LEVELS_FINE];
extern const FIXP_DBL Alphas[NO_ICC_LEVELS];

extern const UCHAR bins2groupMap20[NO_IID_GROUPS];
extern const UCHAR FDK_sbrDecoder_aNoIidBins[3];
extern const UCHAR FDK_sbrDecoder_aNoIccBins[3];

#define INV_TABLE_BITS 8
#define INV_TABLE_SIZE (1 << INV_TABLE_BITS)
extern const FIXP_SGL FDK_sbrDecoder_invTable[INV_TABLE_SIZE];

#endif  // SBR_ROM_H
