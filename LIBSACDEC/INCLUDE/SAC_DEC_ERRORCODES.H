/*********************** MPEG surround decoder library *************************

   Author(s):

   Description: error codes for mpeg surround decoder

*******************************************************************************/

#ifndef SAC_DEC_ERRORCODES_H
#define SAC_DEC_ERRORCODES_H

typedef enum {
  __mps_error_start = -1000,

  MPS_OK = 0,
  MPS_NOTOK = __mps_error_start,

  MPS_OUTOFMEMORY,
  MPS_INVALID_HANDLE,
  MPS_INVALID_PARAMETER,
  MPS_UNSUPPORTED_HRTFMODEL,
  MPS_UNSUPPORTED_HRTFFREQ,

  MPS_UNSUPPORTED_UPMIX_TYPE,
  MPS_UNSUPPORTED_FORMAT,
  MPS_OUTPUT_BUFFER_TOO_SMALL,

  MPS_INVALID_PARAMETERBANDS,
  MPS_INVALID_TREECONFIG,
  MPS_INVALID_HRTFSET,
  MPS_INVALID_TTT,
  MPS_INVALID_BOXIDX,
  MPS_INVALID_SETIDX,
  MPS_INVALID_QUANTMODE,
  MPS_UNEQUAL_SSC,
  MPS_UNSUPPORTED_CONFIG,

  MPS_PARSE_ERROR,
  MPS_INVALID_TEMPSHAPE,

  MPS_WRONG_PARAMETERSETS,
  MPS_WRONG_PARAMETERBANDS,
  MPS_WRONG_TREECONFIG,
  MPS_WRONG_BLINDCONFIG,
  MPS_WRONG_OTT,
  MPS_WRONG_QUANTMODE,
  MPS_RESDEC_ERROR,
  MPS_APPLY_M2_ERROR,

  __mps_error_end
} SACDEC_ERROR;

#endif
