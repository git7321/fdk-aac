/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Psychoaccoustic major function block

*******************************************************************************/

#ifndef PSY_MAIN_H
#define PSY_MAIN_H

#include "psy_configuration.h"
#include "qc_data.h"
#include "aacenc_pns.h"

typedef struct {
  PSY_STATIC *psyStatic[(2)];

} PSY_ELEMENT;

typedef struct {
  PSY_DATA psyData[(2)];
  TNS_DATA tnsData[(2)];
  PNS_DATA pnsData[(2)];

} PSY_DYNAMIC;

typedef struct {
  PSY_CONFIGURATION psyConf[2];
  PSY_ELEMENT *psyElement[((8))];
  PSY_STATIC *pStaticChannels[(8)];
  PSY_DYNAMIC *psyDynamic;
  INT granuleLength;

} PSY_INTERNAL;

AAC_ENCODER_ERROR FDKaacEnc_PsyNew(PSY_INTERNAL **phpsy, const INT nElements,
                                   const INT nChannels, UCHAR *dynamic_RAM);

AAC_ENCODER_ERROR FDKaacEnc_PsyOutNew(PSY_OUT **phpsyOut, const INT nElements,
                                      const INT nChannels, const INT nSubFrames,
                                      UCHAR *dynamic_RAM);

AAC_ENCODER_ERROR FDKaacEnc_psyInit(PSY_INTERNAL *hPsy, PSY_OUT **phpsyOut,
                                    const INT nSubFrames,
                                    const INT nMaxChannels,
                                    const AUDIO_OBJECT_TYPE audioObjectType,
                                    CHANNEL_MAPPING *cm);

AAC_ENCODER_ERROR FDKaacEnc_psyMainInit(
    PSY_INTERNAL *hPsy, AUDIO_OBJECT_TYPE audioObjectType, CHANNEL_MAPPING *cm,
    INT sampleRate, INT granuleLength, INT bitRate, INT tnsMask, INT bandwidth,
    INT usePns, INT useIS, INT useMS, UINT syntaxFlags, ULONG initFlags);

AAC_ENCODER_ERROR FDKaacEnc_psyMain(INT channels, PSY_ELEMENT *psyElement,
                                    PSY_DYNAMIC *psyDynamic,
                                    PSY_CONFIGURATION *psyConf,
                                    PSY_OUT_ELEMENT *psyOutElement,
                                    INT_PCM *pInput, const UINT inputBufSize,
                                    INT *chIdx, INT totalChannels);

void FDKaacEnc_PsyClose(PSY_INTERNAL **phPsyInternal, PSY_OUT **phPsyOut);

#endif /* PSY_MAIN_H */
