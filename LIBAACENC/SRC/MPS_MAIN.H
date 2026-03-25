/**************************** AAC encoder library ******************************

   Author(s):   Markus Lohwasser

   Description: Mpeg Surround library interface functions

*******************************************************************************/

#ifndef MPS_MAIN_H
#define MPS_MAIN_H

/* Includes ******************************************************************/
#include "aacenc.h"
#include "FDK_audio.h"
#include "machine_type.h"

/* Defines *******************************************************************/
typedef enum {
  MPS_ENCODER_OK = 0x0000,
  MPS_ENCODER_INVALID_HANDLE =
      0x0020,
  MPS_ENCODER_MEMORY_ERROR = 0x0021,
  MPS_ENCODER_INIT_ERROR = 0x0040,
  MPS_ENCODER_ENCODE_ERROR =
      0x0060

} MPS_ENCODER_ERROR;

/* Data Types ****************************************************************/

/**
 *  MPEG Surround Encoder interface handle.
 */
typedef struct MPS_ENCODER MPS_ENCODER, *HANDLE_MPS_ENCODER;

/* Function / Class Declarations *********************************************/

/**
 * \brief  Open a Mpeg Surround Encoder instance.
 *
 * \phMpsEnc                    A pointer to a MPS handle to be allocated.
 * Initialized on return.
 *
 * \return
 *          - MPS_ENCODER_OK, on succes.
 *          - MPS_ENCODER_INVALID_HANDLE, MPS_ENCODER_MEMORY_ERROR, on failure.
 */
MPS_ENCODER_ERROR FDK_MpegsEnc_Open(HANDLE_MPS_ENCODER *phMpsEnc);

/**
 * \brief  Close the Mpeg Surround Encoder instance.
 *
 * Deallocate instance and free whole memory.
 *
 * \param phMpsEnc              Pointer to the MPS handle to be deallocated.
 *
 * \return
 *          - MPS_ENCODER_OK, on succes.
 *          - MPS_ENCODER_INVALID_HANDLE, on failure.
 */
MPS_ENCODER_ERROR FDK_MpegsEnc_Close(HANDLE_MPS_ENCODER *phMpsEnc);

/**
 * \brief  Initialize a Mpeg Surround Encoder instance.
 *
 * \param hMpsEnc                   MPS Encoder handle.
 * \param audioObjectType           Audio object type.
 * \param samplingrate              Sampling rate in Hz of audio input signal.
 * \param bitrate                   Encder target bitrate.
 * \param sbrRatio                  SBR sampling rate ratio.
 * \param framelength               Number of samples to be processes within one
 * frame.
 * \param inputBufferSizePerChannel Size of input buffer per channel.
 * \param coreCoderDelay            Core coder delay.
 *
 * \return
 *          - MPS_ENCODER_OK, on succes.
 *          - MPS_ENCODER_INVALID_HANDLE, MPS_ENCODER_ENCODE_ERROR, on failure.
 */
MPS_ENCODER_ERROR FDK_MpegsEnc_Init(HANDLE_MPS_ENCODER hMpsEnc,
                                    const AUDIO_OBJECT_TYPE audioObjectType,
                                    const UINT samplingrate, const UINT bitrate,
                                    const UINT sbrRatio, const UINT framelength,
                                    const UINT inputBufferSizePerChannel,
                                    const UINT coreCoderDelay);

/**
 * \brief  Calculate Mpeg Surround processing.
 *
 * This fuction applies the MPS processing. The MPS side info will be written to
 * extension payload. The input audio data will be overwritten by the calculated
 * downmix.
 *
 * \param hMpsEnc               MPS Encoder handle.
 * \param pAudioSamples         Pointer to audio input/output data.
 * \param nAudioSamples         Number of input audio samples to be prcessed.
 * \param pMpsExtPayload        Pointer to extension payload to be filled on
 * return.
 *
 * \return
 *          - MPS_ENCODER_OK, on succes.
 *          - MPS_ENCODER_INVALID_HANDLE, MPS_ENCODER_ENCODE_ERROR, on failure.
 */
MPS_ENCODER_ERROR FDK_MpegsEnc_Process(HANDLE_MPS_ENCODER hMpsEnc,
                                       INT_PCM *const pAudioSamples,
                                       const INT nAudioSamples,
                                       AACENC_EXT_PAYLOAD *pMpsExtPayload);

/**
 * \brief  Write Spatial Specific Config.
 *
 * This function can be called via call back from the transport library to write
 * the Spatial Specific Config to given bitstream buffer.
 *
 * \param hMpsEnc               MPS Encoder handle.
 * \param hBs                   Bitstream buffer handle.
 *
 * \return                      Number of written bits.
 */
INT FDK_MpegsEnc_WriteSpatialSpecificConfig(HANDLE_MPS_ENCODER hMpsEnc,
                                            HANDLE_FDK_BITSTREAM hBs);

/**
 * \brief  Get closest valid bitrate supported by given config.
 *
 * \param audioObjectType       Audio object type.
 * \param channelMode           Encoder channel mode.
 * \param samplingrate          Sampling rate in Hz of audio input signal.
 * \param sbrRatio              SBR sampling rate ratio.
 * \param bitrate               The desired target bitrate.
 *
 * \return                      Closest valid bitrate to given bitrate..
 */
INT FDK_MpegsEnc_GetClosestBitRate(const AUDIO_OBJECT_TYPE audioObjectType,
                                   const CHANNEL_MODE channelMode,
                                   const UINT samplingrate, const UINT sbrRatio,
                                   const UINT bitrate);

/**
 * \brief  Get codec delay.
 *
 * This function returns delay of the whole en-/decoded signal, including
 * corecoder delay.
 *
 * \param hMpsEnc               MPS Encoder handle.
 *
 * \return                      Codec delay in samples.
 */
INT FDK_MpegsEnc_GetDelay(HANDLE_MPS_ENCODER hMpsEnc);

/**
 * \brief  Get Mpeg Surround Decoder delay.
 *
 * This function returns delay of the Mpeg Surround decoder.
 *
 * \param hMpsEnc               MPS Encoder handle.
 *
 * \return                      Mpeg Surround Decoder delay in samples.
 */
INT FDK_MpegsEnc_GetDecDelay(HANDLE_MPS_ENCODER hMpsEnc);

/**
 * \brief  Get information about encoder library build.
 *
 * Fill a given LIB_INFO structure with library version information.
 *
 * \param info  Pointer to an allocated LIB_INFO struct.
 *
 * \return
 *          - MPS_ENCODER_OK, on succes.
 *          - MPS_ENCODER_INVALID_HANDLE, MPS_ENCODER_INIT_ERROR, on failure.
 */
MPS_ENCODER_ERROR FDK_MpegsEnc_GetLibInfo(LIB_INFO *info);

#endif /* MPS_MAIN_H */
