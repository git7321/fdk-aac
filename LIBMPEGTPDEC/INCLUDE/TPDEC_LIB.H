/******************* MPEG transport format decoder library *********************

   Author(s):   Manuel Jander

   Description: MPEG Transport decoder

*******************************************************************************/

#ifndef TPDEC_LIB_H
#define TPDEC_LIB_H

#include "tp_data.h"
#include "FDK_bitstream.h"

typedef enum {
  TRANSPORTDEC_OK = 0,

  tpdec_sync_error_start = 0x100,
  TRANSPORTDEC_NOT_ENOUGH_BITS,
  TRANSPORTDEC_SYNC_ERROR,
  tpdec_sync_error_end,

  tpdec_decode_error_start = 0x400,
  TRANSPORTDEC_PARSE_ERROR,
  TRANSPORTDEC_UNSUPPORTED_FORMAT,
  TRANSPORTDEC_CRC_ERROR,
  tpdec_decode_error_end,

  tpdec_fatal_error_start = 0x200,
  TRANSPORTDEC_UNKOWN_ERROR,
  TRANSPORTDEC_INVALID_PARAMETER,
  TRANSPORTDEC_NEED_TO_RESTART,
  TRANSPORTDEC_TOO_MANY_BITS,
  tpdec_fatal_error_end

} TRANSPORTDEC_ERROR;

#define TPDEC_IS_DECODE_ERROR(err)                                            \
  (((err >= tpdec_decode_error_start) && (err <= tpdec_decode_error_end)) ? 1 \
                                                                          : 0)
#define TPDEC_IS_FATAL_ERROR(err) \
  (((err >= tpdec_fatal_error_start) && (err <= tpdec_fatal_error_end)) ? 1 : 0)

typedef enum {
  TPDEC_PARAM_MINIMIZE_DELAY = 1,
  TPDEC_PARAM_EARLY_CONFIG,
  TPDEC_PARAM_IGNORE_BUFFERFULLNESS,
  TPDEC_PARAM_SET_BITRATE,
  TPDEC_PARAM_RESET,
  TPDEC_PARAM_BURST_PERIOD,
  TPDEC_PARAM_TARGETLAYOUT,
  TPDEC_PARAM_FORCE_CONFIG_CHANGE,
  TPDEC_PARAM_USE_ELEM_SKIPPING
} TPDEC_PARAM;

void CProgramConfig_Reset(CProgramConfig *pPce);
void CProgramConfig_Init(CProgramConfig *pPce);

int CProgramConfig_IsValid(const CProgramConfig *pPce);

void CProgramConfig_Read(CProgramConfig *pPce, HANDLE_FDK_BITSTREAM bs,
                         UINT alignAnchor);

int CProgramConfig_Compare(const CProgramConfig *const pPce1,
                           const CProgramConfig *const pPce2);

void CProgramConfig_GetDefault(CProgramConfig *pPce, const UINT channelConfig);

int CProgramConfig_LookupElement(CProgramConfig *pPce, UINT chConfig,
                                 const UINT tag, const UINT channelIdx,
                                 UCHAR chMapping[], AUDIO_CHANNEL_TYPE chType[],
                                 UCHAR chIndex[], const UINT chDescrLen,
                                 UCHAR *elMapping, MP4_ELEMENT_ID elList[],
                                 MP4_ELEMENT_ID elType);
int CProgramConfig_GetPceChMap(const CProgramConfig *pPce, UCHAR pceChMap[],
                               const UINT pceChMapLen);
int CProgramConfig_GetElementTable(const CProgramConfig *pPce,
                                   MP4_ELEMENT_ID table[], const INT elListSize,
                                   UCHAR *pChMapIdx);

void CProgramConfig_GetChannelDescription(const UINT chConfig,
                                          const CProgramConfig *pPce,
                                          AUDIO_CHANNEL_TYPE chType[],
                                          UCHAR chIndex[]);
void AudioSpecificConfig_Init(CSAudioSpecificConfig *pAsc);

TRANSPORTDEC_ERROR AudioSpecificConfig_Parse(
    CSAudioSpecificConfig *pAsc, HANDLE_FDK_BITSTREAM hBs,
    int fExplicitBackwardCompatible, CSTpCallBacks *cb, UCHAR configMode,
    UCHAR configChanged, AUDIO_OBJECT_TYPE m_aot);

enum { MPE = 0, RPE = 1, fs8KHz = 0, fs16KHz = 1 };

#define TP_FLAG_MPEG4 1

#define CAPF_TPDEC_ADIF \
  0x00001000
#define CAPF_TPDEC_ADTS \
  0x00002000
#define CAPF_TPDEC_LOAS \
  0x00004000
#define CAPF_TPDEC_LATM \
  0x00008000
#define CAPF_TPDEC_RAWPACKETS \
  0x00010000

typedef struct TRANSPORTDEC *HANDLE_TRANSPORTDEC;

TRANSPORTDEC_ERROR transportDec_OutOfBandConfig(const HANDLE_TRANSPORTDEC hTp,
                                                UCHAR *conf, const UINT length,
                                                const UINT layer);
TRANSPORTDEC_ERROR transportDec_InBandConfig(
    const HANDLE_TRANSPORTDEC hTp, UCHAR *newConfig, const UINT newConfigLength,
    const UCHAR buildUpStatus, UCHAR *configChanged, const UINT layer,
    UCHAR *implicitExplicitCfgDiff);

HANDLE_TRANSPORTDEC transportDec_Open(TRANSPORT_TYPE transportDecFmt,
                                      const UINT flags, const UINT nrOfLayer);

int transportDec_RegisterAscCallback(HANDLE_TRANSPORTDEC hTp,
                                     const cbUpdateConfig_t cbUpdateConfig,
                                     void *user_data);
int transportDec_RegisterFreeMemCallback(HANDLE_TRANSPORTDEC hTp,
                                         const cbFreeMem_t cbFreeMem,
                                         void *user_data);
int transportDec_RegisterCtrlCFGChangeCallback(
    HANDLE_TRANSPORTDEC hTp, const cbCtrlCFGChange_t cbCtrlCFGChange,
    void *user_data);
int transportDec_RegisterSscCallback(HANDLE_TRANSPORTDEC hTp,
                                     const cbSsc_t cbSscParse, void *user_data);
int transportDec_RegisterSbrCallback(HANDLE_TRANSPORTDEC hTpDec,
                                     const cbSbr_t cbSbr, void *user_data);
int transportDec_RegisterUsacCallback(HANDLE_TRANSPORTDEC hTpDec,
                                      const cbUsac_t cbUsac, void *user_data);
int transportDec_RegisterUniDrcConfigCallback(HANDLE_TRANSPORTDEC hTpDec,
                                              const cbUniDrc_t cbUniDrc,
                                              void *user_data,
                                              UINT *pLoudnessInfoSetPosition);

TRANSPORTDEC_ERROR transportDec_FillData(const HANDLE_TRANSPORTDEC hTp,
                                         UCHAR *pBuffer, const UINT bufferSize,
                                         UINT *pBytesValid, const INT layer);

HANDLE_FDK_BITSTREAM transportDec_GetBitstream(const HANDLE_TRANSPORTDEC hTp,
                                               const UINT layer);

TRANSPORT_TYPE transportDec_GetFormat(const HANDLE_TRANSPORTDEC hTp);

INT transportDec_GetBufferFullness(const HANDLE_TRANSPORTDEC hTp);

void transportDec_Close(HANDLE_TRANSPORTDEC *phTp);

TRANSPORTDEC_ERROR transportDec_ReadAccessUnit(const HANDLE_TRANSPORTDEC hTp,
                                               const UINT layer);
TRANSPORTDEC_ERROR transportDec_GetAsc(const HANDLE_TRANSPORTDEC hTp,
                                       const UINT layer,
                                       CSAudioSpecificConfig *asc);

INT transportDec_GetAuBitsRemaining(const HANDLE_TRANSPORTDEC hTp,
                                    const UINT layer);
INT transportDec_GetAuBitsTotal(const HANDLE_TRANSPORTDEC hTp,
                                const UINT layer);

TRANSPORTDEC_ERROR transportDec_EndAccessUnit(const HANDLE_TRANSPORTDEC hTp);
TRANSPORTDEC_ERROR transportDec_GetMissingAccessUnitCount(
    INT *pNAccessUnits, HANDLE_TRANSPORTDEC hTp);
TRANSPORTDEC_ERROR transportDec_SetParam(const HANDLE_TRANSPORTDEC hTp,
                                         const TPDEC_PARAM param,
                                         const INT value);

UINT transportDec_GetNrOfSubFrames(HANDLE_TRANSPORTDEC hTp);

TRANSPORTDEC_ERROR transportDec_GetLibInfo(LIB_INFO *info);

int transportDec_CrcStartReg(const HANDLE_TRANSPORTDEC hTp, const INT mBits);

void transportDec_CrcEndReg(const HANDLE_TRANSPORTDEC hTp, const INT reg);

TRANSPORTDEC_ERROR transportDec_CrcCheck(const HANDLE_TRANSPORTDEC hTp);
TRANSPORTDEC_ERROR transportDec_DrmRawSdcAudioConfig_Check(UCHAR *conf,
                                                           const UINT length);

#endif /* #ifndef TPDEC_LIB_H */
