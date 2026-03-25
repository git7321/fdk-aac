/**************************** SBR encoder library ******************************

   Author(s):   M. Neuendorf, N. Rettelbach, M. Multrus

   Description: PS Parameter extraction, encoding

*******************************************************************************/

#ifndef PS_ENCODE_H
#define PS_ENCODE_H

#include "ps_const.h"
#include "ps_bitenc.h"

#define IID_SCALE_FT (64.f)
#define IID_SCALE 6
#define IID_MAXVAL (1 << IID_SCALE)

#define PS_QUANT_SCALE_FT \
  (64.f)
#define PS_QUANT_SCALE \
  6

#define QMF_GROUPS_LO_RES 12
#define SUBQMF_GROUPS_LO_RES 10
#define QMF_GROUPS_HI_RES 18
#define SUBQMF_GROUPS_HI_RES 30

typedef struct T_PS_DATA {
  INT iidEnable;
  INT iidEnableLast;
  INT iidQuantMode;
  INT iidQuantModeLast;
  INT iidDiffMode[PS_MAX_ENVELOPES];
  INT iidIdx[PS_MAX_ENVELOPES][PS_MAX_BANDS];
  INT iidIdxLast[PS_MAX_BANDS];

  INT iccEnable;
  INT iccEnableLast;
  INT iccQuantMode;
  INT iccQuantModeLast;
  INT iccDiffMode[PS_MAX_ENVELOPES];
  INT iccIdx[PS_MAX_ENVELOPES][PS_MAX_BANDS];
  INT iccIdxLast[PS_MAX_BANDS];

  INT nEnvelopesLast;

  INT headerCnt;
  INT iidTimeCnt;
  INT iccTimeCnt;
  INT noEnvCnt;

} PS_DATA, *HANDLE_PS_DATA;

typedef struct T_PS_ENCODE {
  PS_DATA psData;

  PS_BANDS psEncMode;
  INT nQmfIidGroups;
  INT nSubQmfIidGroups;
  INT iidGroupBorders[QMF_GROUPS_HI_RES + SUBQMF_GROUPS_HI_RES + 1];
  INT subband2parameterIndex[QMF_GROUPS_HI_RES + SUBQMF_GROUPS_HI_RES];
  UCHAR iidGroupWidthLd[QMF_GROUPS_HI_RES + SUBQMF_GROUPS_HI_RES];
  FIXP_DBL iidQuantErrorThreshold;

  UCHAR psBandNrgScale[PS_MAX_BANDS];

} PS_ENCODE;

typedef struct T_PS_ENCODE *HANDLE_PS_ENCODE;

FDK_PSENC_ERROR FDKsbrEnc_CreatePSEncode(HANDLE_PS_ENCODE *phPsEncode);

FDK_PSENC_ERROR FDKsbrEnc_InitPSEncode(HANDLE_PS_ENCODE hPsEncode,
                                       const PS_BANDS psEncMode,
                                       const FIXP_DBL iidQuantErrorThreshold);

FDK_PSENC_ERROR FDKsbrEnc_DestroyPSEncode(HANDLE_PS_ENCODE *phPsEncode);

FDK_PSENC_ERROR FDKsbrEnc_PSEncode(
    HANDLE_PS_ENCODE hPsEncode, HANDLE_PS_OUT hPsOut, UCHAR *dynBandScale,
    UINT maxEnvelopes,
    FIXP_DBL *hybridData[HYBRID_FRAMESIZE][MAX_PS_CHANNELS][2],
    const INT frameSize, const INT sendHeader);

#endif
