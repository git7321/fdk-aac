/**************************** AAC encoder library ******************************

   Author(s):   M. Neusinger

   Description: Compressor for AAC Metadata Generator

*******************************************************************************/

#ifndef METADATA_COMPRESSOR_H
#define METADATA_COMPRESSOR_H

#include "FDK_audio.h"
#include "common_fix.h"

#include "aacenc.h"

#define LFE_LEV_SCALE 2

typedef enum DRC_PROFILE {
  DRC_NONE = 0,
  DRC_FILMSTANDARD = 1,
  DRC_FILMLIGHT = 2,
  DRC_MUSICSTANDARD = 3,
  DRC_MUSICLIGHT = 4,
  DRC_SPEECH = 5,
  DRC_DELAY_TEST = 6,
  DRC_NOT_PRESENT = -2

} DRC_PROFILE;

typedef struct DRC_COMP DRC_COMP, *HDRC_COMP;

INT FDK_DRC_Generator_Open(HDRC_COMP *phDrcComp);

INT FDK_DRC_Generator_Close(HDRC_COMP *phDrcComp);

INT FDK_DRC_Generator_Initialize(HDRC_COMP drcComp,
                                 const DRC_PROFILE profileLine,
                                 const DRC_PROFILE profileRF,
                                 const INT blockLength, const UINT sampleRate,
                                 const CHANNEL_MODE channelMode,
                                 const CHANNEL_ORDER channelOrder,
                                 const UCHAR useWeighting);

INT FDK_DRC_Generator_Calc(HDRC_COMP drcComp, const INT_PCM *const inSamples,
                           const UINT inSamplesBufSize, const INT dialnorm,
                           const INT drc_TargetRefLevel,
                           const INT comp_TargetRefLevel, const FIXP_DBL clev,
                           const FIXP_DBL slev, const FIXP_DBL ext_leva,
                           const FIXP_DBL ext_levb, const FIXP_DBL lfe_lev,
                           const INT dmxGain5, const INT dmxGain2,
                           INT *const dynrng, INT *const compr);

INT FDK_DRC_Generator_setDrcProfile(HDRC_COMP drcComp,
                                    const DRC_PROFILE profileLine,
                                    const DRC_PROFILE profileRF);

DRC_PROFILE FDK_DRC_Generator_getDrcProfile(const HDRC_COMP drcComp);

DRC_PROFILE FDK_DRC_Generator_getCompProfile(const HDRC_COMP drcComp);

#endif /* METADATA_COMPRESSOR_H */
