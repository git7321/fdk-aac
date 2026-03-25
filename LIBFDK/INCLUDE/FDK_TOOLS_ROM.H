/******************* Library for basic calculation routines ********************

   Author(s):   Oliver Moser

   Description: ROM tables used by FDK tools

*******************************************************************************/

#ifndef FDK_TOOLS_ROM_H
#define FDK_TOOLS_ROM_H

#include "common_fix.h"
#include "FDK_audio.h"

extern RAM_ALIGN const FIXP_STB RotVectorReal60[60];
extern RAM_ALIGN const FIXP_STB RotVectorImag60[60];
extern RAM_ALIGN const FIXP_STB RotVectorReal192[192];
extern RAM_ALIGN const FIXP_STB RotVectorImag192[192];
extern RAM_ALIGN const FIXP_STB RotVectorReal240[210];
extern RAM_ALIGN const FIXP_STB RotVectorImag240[210];
extern RAM_ALIGN const FIXP_STB RotVectorReal480[480];
extern RAM_ALIGN const FIXP_STB RotVectorImag480[480];
extern RAM_ALIGN const FIXP_STB RotVectorReal6[6];
extern RAM_ALIGN const FIXP_STB RotVectorImag6[6];
extern RAM_ALIGN const FIXP_STB RotVectorReal12[12];
extern RAM_ALIGN const FIXP_STB RotVectorImag12[12];
extern RAM_ALIGN const FIXP_STB RotVectorReal24[24];
extern RAM_ALIGN const FIXP_STB RotVectorImag24[24];
extern RAM_ALIGN const FIXP_STB RotVectorReal48[48];
extern RAM_ALIGN const FIXP_STB RotVectorImag48[48];
extern RAM_ALIGN const FIXP_STB RotVectorReal80[80];
extern RAM_ALIGN const FIXP_STB RotVectorImag80[80];
extern RAM_ALIGN const FIXP_STB RotVectorReal96[96];
extern RAM_ALIGN const FIXP_STB RotVectorImag96[96];
extern RAM_ALIGN const FIXP_STB RotVectorReal384[384];
extern RAM_ALIGN const FIXP_STB RotVectorImag384[384];
extern RAM_ALIGN const FIXP_STB RotVectorReal20[20];
extern RAM_ALIGN const FIXP_STB RotVectorImag20[20];
extern RAM_ALIGN const FIXP_STB RotVectorReal120[120];
extern RAM_ALIGN const FIXP_STB RotVectorImag120[120];

extern RAM_ALIGN const FIXP_STP SineTable1024[];
extern RAM_ALIGN const FIXP_STP SineTable512[];
extern RAM_ALIGN const FIXP_STP SineTable480[];
extern RAM_ALIGN const FIXP_STP SineTable384[];
extern RAM_ALIGN const FIXP_STP SineTable80[];
#ifdef INCLUDE_SineTable10
extern RAM_ALIGN const FIXP_STP SineTable10[];
#endif

extern RAM_ALIGN const FIXP_WTP SineWindow1024[];
extern RAM_ALIGN const FIXP_WTP KBDWindow1024[];
extern RAM_ALIGN const FIXP_WTP SineWindow128[];
extern RAM_ALIGN const FIXP_WTP KBDWindow128[];

extern RAM_ALIGN const FIXP_WTP SineWindow960[];
extern RAM_ALIGN const FIXP_WTP KBDWindow960[];
extern RAM_ALIGN const FIXP_WTP SineWindow120[];
extern RAM_ALIGN const FIXP_WTP KBDWindow120[];

extern RAM_ALIGN const FIXP_WTP SineWindow512[];
#define LowOverlapWindow512 SineWindow128
extern RAM_ALIGN const FIXP_WTP SineWindow480[];
#define LowOverlapWindow480 SineWindow120

extern RAM_ALIGN const FIXP_WTP SineWindow256[256];
extern RAM_ALIGN const FIXP_WTP SineWindow192[];

extern RAM_ALIGN const FIXP_WTP SineWindow768[];
extern RAM_ALIGN const FIXP_WTP KBDWindow768[];
extern RAM_ALIGN const FIXP_WTP SineWindow96[];
extern RAM_ALIGN const FIXP_WTP KBDWindow96[];

extern RAM_ALIGN const FIXP_WTP SineWindow64[];
extern RAM_ALIGN const FIXP_WTP SineWindow48[];
extern RAM_ALIGN const FIXP_WTP SineWindow32[];
extern RAM_ALIGN const FIXP_WTP SineWindow24[];
extern RAM_ALIGN const FIXP_WTP SineWindow16[];
extern RAM_ALIGN const FIXP_WTP SineWindow8[];

extern const FIXP_WTP *const windowSlopes[2][4][9];

const FIXP_WTP *FDKgetWindowSlope(int length, int shape);

extern const FIXP_WTP sin_twiddle_L64[];

#if defined(ARCH_PREFER_MULT_16x16) || defined(ARCH_PREFER_MULT_32x16)
#define QMF_COEFF_16BIT
#endif

#define QMF_FILTER_PROTOTYPE_SIZE 640
#define QMF_NO_POLY 5

#ifdef QMF_COEFF_16BIT
#define FIXP_PFT FIXP_SGL
#define FIXP_QTW FIXP_SGL
#define FX_DBL2FX_QTW(x) FX_DBL2FX_SGL(x)
#else
#define FIXP_PFT FIXP_DBL
#define FIXP_QTW FIXP_DBL

#define FX_DBL2FX_QTW(x) (x)

#endif

#define QMF640_PFT_TABLE_SIZE (640 / 2 + QMF_NO_POLY)

extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_cos32[32];
extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_sin32[32];
extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_cos_downsamp32[32];
extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_sin_downsamp32[32];
extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_cos64[64];
extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_sin64[64];
extern RAM_ALIGN const FIXP_PFT
    qmf_pfilt640[QMF640_PFT_TABLE_SIZE + QMF_NO_POLY];
extern RAM_ALIGN const FIXP_PFT qmf_pfilt640_vector[640];

extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_cos40[40];
extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_sin40[40];
extern RAM_ALIGN const FIXP_PFT qmf_pfilt400[];
extern RAM_ALIGN const FIXP_PFT qmf_pfilt200[];
extern RAM_ALIGN const FIXP_PFT qmf_pfilt120[];

extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_cos24[24];
extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_sin24[24];
extern RAM_ALIGN const FIXP_PFT qmf_pfilt240[];

extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_cos16[16];
extern RAM_ALIGN const FIXP_QTW qmf_phaseshift_sin16[16];

#define QMF640_CLDFB_PFT_TABLE_SIZE (640)
#define QMF320_CLDFB_PFT_TABLE_SIZE (320)
#define QMF_CLDFB_PFT_SCALE 1

extern const FIXP_QTW qmf_phaseshift_cos32_cldfb_ana[32];
extern const FIXP_QTW qmf_phaseshift_cos32_cldfb_syn[32];
extern const FIXP_QTW qmf_phaseshift_sin32_cldfb[32];

extern const FIXP_QTW qmf_phaseshift_cos16_cldfb_ana[16];
extern const FIXP_QTW qmf_phaseshift_cos16_cldfb_syn[16];
extern const FIXP_QTW qmf_phaseshift_sin16_cldfb[16];

extern const FIXP_QTW qmf_phaseshift_cos8_cldfb_ana[8];
extern const FIXP_QTW qmf_phaseshift_cos8_cldfb_syn[8];
extern const FIXP_QTW qmf_phaseshift_sin8_cldfb[8];

extern const FIXP_QTW qmf_phaseshift_cos64_cldfb[64];
extern const FIXP_QTW qmf_phaseshift_sin64_cldfb[64];

extern RAM_ALIGN const FIXP_PFT qmf_cldfb_640[QMF640_CLDFB_PFT_TABLE_SIZE];
extern RAM_ALIGN const FIXP_PFT qmf_cldfb_320[QMF320_CLDFB_PFT_TABLE_SIZE];
#define QMF160_CLDFB_PFT_TABLE_SIZE (160)
extern RAM_ALIGN const FIXP_PFT qmf_cldfb_160[QMF160_CLDFB_PFT_TABLE_SIZE];
#define QMF80_CLDFB_PFT_TABLE_SIZE (80)
extern RAM_ALIGN const FIXP_PFT qmf_cldfb_80[QMF80_CLDFB_PFT_TABLE_SIZE];

#define QMF320_MPSLDFB_PFT_TABLE_SIZE (320)
#define QMF640_MPSLDFB_PFT_TABLE_SIZE (640)
#define QMF_MPSLDFB_PFT_SCALE 1

extern const FIXP_PFT qmf_mpsldfb_320[QMF320_MPSLDFB_PFT_TABLE_SIZE];
extern RAM_ALIGN const FIXP_PFT qmf_mpsldfb_640[QMF640_MPSLDFB_PFT_TABLE_SIZE];

#define AC_EL_GA_CCE 0x00000001

typedef enum {
  element_instance_tag,
  common_window,
  global_gain,
  ics_info,
  max_sfb,
  ms,
  ltp_data_present,
  ltp_data,
  section_data,
  scale_factor_data,
  pulse,
  tns_data_present,
  tns_data,
  gain_control_data_present,
  gain_control_data,
  esc1_hcr,
  esc2_rvlc,
  spectral_data,

  scale_factor_data_usac,
  core_mode,
  common_tw,
  lpd_channel_stream,
  tw_data,
  noise,
  ac_spectral_data,
  fac_data,
  tns_active,
  tns_data_present_usac,
  common_max_sfb,

  coupled_elements,
  gain_element_lists,

  adtscrc_start_reg1,
  adtscrc_start_reg2,
  adtscrc_end_reg1,
  adtscrc_end_reg2,
  drmcrc_start_reg,
  drmcrc_end_reg,
  next_channel,
  next_channel_loop,
  link_sequence,
  end_of_sequence
} rbd_id_t;

struct element_list {
  const rbd_id_t *id;
  const struct element_list *next[2];
};

typedef struct element_list element_list_t;
const element_list_t *getBitstreamElementList(AUDIO_OBJECT_TYPE aot,
                                              SCHAR epConfig, UCHAR nChannels,
                                              UCHAR layer, UINT elFlags);

typedef enum {
  FDK_FORMAT_1_0 = 1,
  FDK_FORMAT_2_0 = 2,
  FDK_FORMAT_3_0_FC = 3,
  FDK_FORMAT_3_1_0 = 4,
  FDK_FORMAT_5_0 = 5,
  FDK_FORMAT_5_1 = 6,
  FDK_FORMAT_7_1_ALT = 7,
  FDK_FORMAT_3_0_RC = 9,
  FDK_FORMAT_2_2_0 = 10,
  FDK_FORMAT_6_1 = 11,
  FDK_FORMAT_7_1 = 12,
  FDK_FORMAT_22_2 = 13,
  FDK_FORMAT_5_2_1 = 14,
  FDK_FORMAT_5_5_2 = 15,
  FDK_FORMAT_9_1 = 16,
  FDK_FORMAT_6_5_1 = 17,
  FDK_FORMAT_6_7_1 = 18,
  FDK_FORMAT_5_6_1 = 19,
  FDK_FORMAT_7_6_1 = 20,
  FDK_FORMAT_IN_LISTOFCHANNELS = 21,
  FDK_FORMAT_OUT_LISTOFCHANNELS = 22,
  FDK_NFORMATS = 23,
  FDK_FORMAT_FAIL = -1
} FDK_converter_formatid_t;

extern const INT format_nchan[FDK_NFORMATS + 9 - 2];

#endif
