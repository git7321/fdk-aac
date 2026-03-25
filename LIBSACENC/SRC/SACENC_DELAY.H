/*********************** MPEG surround encoder library *************************

   Author(s):   Christian Goettlinger

   Description: Encoder Library Interface
                delay management of the encoder

*******************************************************************************/

#ifndef SACENC_DELAY_H
#define SACENC_DELAY_H

/* Includes ******************************************************************/
#include "sacenc_lib.h"
#include "machine_type.h"
#include "FDK_matrixCalloc.h"

/* Defines *******************************************************************/
#define MAX_DELAY_INPUT 1024
#define MAX_DELAY_OUTPUT 4096
#define MAX_DELAY_SURROUND_ANALYSIS 5
#define MAX_BITSTREAM_DELAY 1

/* Data Types ****************************************************************/
typedef struct DELAY *HANDLE_DELAY;

FDK_SACENC_ERROR fdk_sacenc_delay_Open(HANDLE_DELAY *phDelay);

FDK_SACENC_ERROR fdk_sacenc_delay_Close(HANDLE_DELAY *phDelay);

FDK_SACENC_ERROR fdk_sacenc_delay_Init(HANDLE_DELAY hDelay, const INT nQmfLen,
                                       const INT nFrameLen,
                                       const INT nCoreCoderDelay,
                                       const INT nSacStreamMuxDelay);

FDK_SACENC_ERROR fdk_sacenc_delay_SubCalulateBufferDelays(HANDLE_DELAY hDel);

/* Set Expert Config Parameters */
FDK_SACENC_ERROR fdk_sacenc_delay_SetDmxAlign(HANDLE_DELAY hDelay,
                                              const INT bDmxAlignIn);

FDK_SACENC_ERROR fdk_sacenc_delay_SetTimeDomDmx(HANDLE_DELAY hDelay,
                                                const INT bTimeDomDmxIn);

FDK_SACENC_ERROR fdk_sacenc_delay_SetSacTimeAlignmentDynamicOut(
    HANDLE_DELAY hDelay, const INT bSacTimeAlignmentDynamicOutIn);

FDK_SACENC_ERROR fdk_sacenc_delay_SetNSacTimeAlignment(
    HANDLE_DELAY hDelay, const INT nSacTimeAlignmentIn);

FDK_SACENC_ERROR fdk_sacenc_delay_SetMinimizeDelay(HANDLE_DELAY hDelay,
                                                   const INT bMinimizeDelay);

/* Get Internal Variables */
INT fdk_sacenc_delay_GetOutputAudioBufferDelay(HANDLE_DELAY hDelay);

INT fdk_sacenc_delay_GetSurroundAnalysisBufferDelay(HANDLE_DELAY hDelay);

INT fdk_sacenc_delay_GetArbDmxAnalysisBufferDelay(HANDLE_DELAY hDelay);

INT fdk_sacenc_delay_GetBitstreamFrameBufferSize(HANDLE_DELAY hDelay);

INT fdk_sacenc_delay_GetDmxAlignBufferDelay(HANDLE_DELAY hDelay);

INT fdk_sacenc_delay_GetDiscardOutFrames(HANDLE_DELAY hDelay);

INT fdk_sacenc_delay_GetInfoDmxDelay(HANDLE_DELAY hDelay);

INT fdk_sacenc_delay_GetInfoCodecDelay(HANDLE_DELAY hDelay);

INT fdk_sacenc_delay_GetInfoDecoderDelay(HANDLE_DELAY hDelay);

#endif /* SACENC_DELAY_H */
