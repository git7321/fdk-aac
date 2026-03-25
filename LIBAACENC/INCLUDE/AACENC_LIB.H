/**************************** AAC encoder library ******************************

   Author(s):   M. Lohwasser

   Description:

*******************************************************************************/

#ifndef AACENC_LIB_H
#define AACENC_LIB_H

#include "..\..\libsys\include\machine_type.h"
#include "..\..\libsys\include\FDK_audio.h"

#define AACENCODER_LIB_VL0 4
#define AACENCODER_LIB_VL1 0
#define AACENCODER_LIB_VL2 1

/**
 *  AAC encoder error codes.
 */
typedef enum {
  AACENC_OK = 0x0000, /*!< No error happened. All fine. */

  AACENC_INVALID_HANDLE =
      0x0020, /*!< Handle passed to function call was invalid. */
  AACENC_MEMORY_ERROR = 0x0021,          /*!< Memory allocation failed. */
  AACENC_UNSUPPORTED_PARAMETER = 0x0022, /*!< Parameter not available. */
  AACENC_INVALID_CONFIG = 0x0023,        /*!< Configuration not provided. */

  AACENC_INIT_ERROR = 0x0040,     /*!< General initialization error. */
  AACENC_INIT_AAC_ERROR = 0x0041, /*!< AAC library initialization error. */
  AACENC_INIT_SBR_ERROR = 0x0042, /*!< SBR library initialization error. */
  AACENC_INIT_TP_ERROR = 0x0043, /*!< Transport library initialization error. */
  AACENC_INIT_META_ERROR =
      0x0044, /*!< Meta data library initialization error. */
  AACENC_INIT_MPS_ERROR = 0x0045, /*!< MPS library initialization error. */

  AACENC_ENCODE_ERROR = 0x0060, /*!< The encoding process was interrupted by an
                                   unexpected error. */

  AACENC_ENCODE_EOF = 0x0080 /*!< End of file reached. */

} AACENC_ERROR;

/**
 *  AAC encoder buffer descriptors identifier.
 *  This identifier are used within buffer descriptors
 * AACENC_BufDesc::bufferIdentifiers.
 */
typedef enum {
  /* Input buffer identifier. */
  IN_AUDIO_DATA = 0,    /*!< Audio input buffer, interleaved INT_PCM samples. */
  IN_ANCILLRY_DATA = 1, /*!< Ancillary data to be embedded into bitstream. */
  IN_METADATA_SETUP = 2, /*!< Setup structure for embedding meta data. */

  /* Output buffer identifier. */
  OUT_BITSTREAM_DATA = 3, /*!< Buffer holds bitstream output data. */
  OUT_AU_SIZES =
      4 /*!< Buffer contains sizes of each access unit. This information
             is necessary for superframing. */

} AACENC_BufferIdentifier;

/**
 *  AAC encoder handle.
 */
typedef struct AACENCODER *HANDLE_AACENCODER;

/**
 *  Provides some info about the encoder configuration.
 */
typedef struct {
  UINT maxOutBufBytes; /*!< Maximum number of encoder bitstream bytes within one
                          frame. Size depends on maximum number of supported
                          channels in encoder instance. */

  UINT maxAncBytes; /*!< Maximum number of ancillary data bytes which can be
                       inserted into bitstream within one frame. */

  UINT inBufFillLevel; /*!< Internal input buffer fill level in samples per
                          channel. This parameter will automatically be cleared
                          if samplingrate or channel(Mode/Order) changes. */

  UINT inputChannels; /*!< Number of input channels expected in encoding
                         process. */

  UINT frameLength; /*!< Amount of input audio samples consumed each frame per
                       channel, depending on audio object type configuration. */

  UINT nDelay; /*!< Codec delay in PCM samples/channel. Depends on framelength
                  and AOT. Does not include framing delay for filling up encoder
                  PCM input buffer. */

  UINT nDelayCore; /*!< Codec delay in PCM samples/channel, w/o delay caused by
                      the decoder SBR module. This delay is needed to correctly
                      write edit lists for gapless playback. The decoder may not
                      know how much delay is introdcued by SBR, since it may not
                      know if SBR is active at all (implicit signaling),
                      therefore the decoder must take into account any delay
                      caused by the SBR module. */

  UCHAR confBuf[64]; /*!< Configuration buffer in binary format as an
                        AudioSpecificConfig or StreamMuxConfig according to the
                        selected transport type. */

  UINT confSize; /*!< Number of valid bytes in confBuf. */

} AACENC_InfoStruct;

/**
 *  Describes the input and output buffers for an aacEncEncode() call.
 */
typedef struct {
  INT numBufs;            /*!< Number of buffers. */
  void **bufs;            /*!< Pointer to vector containing buffer addresses. */
  INT *bufferIdentifiers; /*!< Identifier of each buffer element. See
                             ::AACENC_BufferIdentifier. */
  INT *bufSizes;          /*!< Size of each buffer in 8-bit bytes. */
  INT *bufElSizes;        /*!< Size of each buffer element in bytes. */

} AACENC_BufDesc;

/**
 *  Defines the input arguments for an aacEncEncode() call.
 */
typedef struct {
  INT numInSamples; /*!< Number of valid input audio samples (multiple of input
                       channels). */
  INT numAncBytes;  /*!< Number of ancillary data bytes to be encoded. */

} AACENC_InArgs;

/**
 *  Defines the output arguments for an aacEncEncode() call.
 */
typedef struct {
  INT numOutBytes;  /*!< Number of valid bitstream bytes generated during
                       aacEncEncode(). */
  INT numInSamples; /*!< Number of input audio samples consumed by the encoder.
                     */
  INT numAncBytes;  /*!< Number of ancillary data bytes consumed by the encoder.
                     */
  INT bitResState;  /*!< State of the bit reservoir in bits. */

} AACENC_OutArgs;

/**
 *  Meta Data Compression Profiles.
 */
typedef enum {
  AACENC_METADATA_DRC_NONE = 0,          /*!< None. */
  AACENC_METADATA_DRC_FILMSTANDARD = 1,  /*!< Film standard. */
  AACENC_METADATA_DRC_FILMLIGHT = 2,     /*!< Film light. */
  AACENC_METADATA_DRC_MUSICSTANDARD = 3, /*!< Music standard. */
  AACENC_METADATA_DRC_MUSICLIGHT = 4,    /*!< Music light. */
  AACENC_METADATA_DRC_SPEECH = 5,        /*!< Speech. */
  AACENC_METADATA_DRC_NOT_PRESENT =
      256 /*!< Disable writing gain factor (used for comp_profile only). */

} AACENC_METADATA_DRC_PROFILE;

/**
 *  Meta Data setup structure.
 */
typedef struct {
  AACENC_METADATA_DRC_PROFILE
  drc_profile; /*!< MPEG DRC compression profile. See
                  ::AACENC_METADATA_DRC_PROFILE. */
  AACENC_METADATA_DRC_PROFILE
  comp_profile; /*!< ETSI heavy compression profile. See
                   ::AACENC_METADATA_DRC_PROFILE. */

  INT drc_TargetRefLevel;  /*!< Used to define expected level to:
                                Scaled with 16 bit. x*2^16. */
  INT comp_TargetRefLevel; /*!< Adjust limiter to avoid overload.
                                Scaled with 16 bit. x*2^16. */

  INT prog_ref_level_present; /*!< Flag, if prog_ref_level is present */
  INT prog_ref_level;         /*!< Programme Reference Level = Dialogue Level:
                                   -31.75dB .. 0 dB ; stepsize: 0.25dB
                                   Scaled with 16 bit. x*2^16.*/

  UCHAR PCE_mixdown_idx_present; /*!< Flag, if dmx-idx should be written in
                                    programme config element */
  UCHAR ETSI_DmxLvl_present;     /*!< Flag, if dmx-lvl should be written in
                                    ETSI-ancData */

  SCHAR centerMixLevel; /*!< Center downmix level (0...7, according to table) */
  SCHAR surroundMixLevel; /*!< Surround downmix level (0...7, according to
                             table) */

  UCHAR
  dolbySurroundMode; /*!< Indication for Dolby Surround Encoding Mode.
                          - 0: Dolby Surround mode not indicated
                          - 1: 2-ch audio part is not Dolby surround encoded
                          - 2: 2-ch audio part is Dolby surround encoded */

  UCHAR drcPresentationMode; /*!< Indicatin for DRC Presentation Mode.
                                  - 0: Presentation mode not inticated
                                  - 1: Presentation mode 1
                                  - 2: Presentation mode 2 */

  struct {
    /* extended ancillary data */
    UCHAR extAncDataEnable; /*< Indicates if MPEG4_ext_ancillary_data() exists.
                                - 0: No MPEG4_ext_ancillary_data().
                                - 1: Insert MPEG4_ext_ancillary_data(). */

    UCHAR
    extDownmixLevelEnable;   /*< Indicates if ext_downmixing_levels() exists.
                                 - 0: No ext_downmixing_levels().
                                 - 1: Insert ext_downmixing_levels(). */
    UCHAR extDownmixLevel_A; /*< Downmix level index A (0...7, according to
                                table) */
    UCHAR extDownmixLevel_B; /*< Downmix level index B (0...7, according to
                                table) */

    UCHAR dmxGainEnable; /*< Indicates if ext_downmixing_global_gains() exists.
                             - 0: No ext_downmixing_global_gains().
                             - 1: Insert ext_downmixing_global_gains(). */
    INT dmxGain5;        /*< Gain factor for downmix to 5 channels.
                              -15.75dB .. -15.75dB; stepsize: 0.25dB
                              Scaled with 16 bit. x*2^16.*/
    INT dmxGain2;        /*< Gain factor for downmix to 2 channels.
                              -15.75dB .. -15.75dB; stepsize: 0.25dB
                              Scaled with 16 bit. x*2^16.*/

    UCHAR lfeDmxEnable; /*< Indicates if ext_downmixing_lfe_level() exists.
                            - 0: No ext_downmixing_lfe_level().
                            - 1: Insert ext_downmixing_lfe_level(). */
    UCHAR lfeDmxLevel;  /*< Downmix level index for LFE (0..15, according to
                           table) */

  } ExtMetaData;

} AACENC_MetaData;

/**
 * AAC encoder control flags.
 *
 * In interaction with the ::AACENC_CONTROL_STATE parameter it is possible to
 * get information about the internal initialization process. It is also
 * possible to overwrite the internal state from extern when necessary.
 */
typedef enum {
  AACENC_INIT_NONE = 0x0000, /*!< Do not trigger initialization. */
  AACENC_INIT_CONFIG =
      0x0001, /*!< Initialize all encoder modules configuration. */
  AACENC_INIT_STATES = 0x0002, /*!< Reset all encoder modules history buffer. */
  AACENC_INIT_TRANSPORT =
      0x1000, /*!< Initialize transport lib with new parameters. */
  AACENC_RESET_INBUFFER =
      0x2000,              /*!< Reset fill level of internal input buffer. */
  AACENC_INIT_ALL = 0xFFFF /*!< Initialize all. */
} AACENC_CTRLFLAGS;

/**
 * \brief  AAC encoder setting parameters.
 *
 * Use aacEncoder_SetParam() function to configure, or use aacEncoder_GetParam()
 * function to read the internal status of the following parameters.
 */
typedef enum {
  AACENC_AOT =
      0x0100, /*!< Audio object type. See ::AUDIO_OBJECT_TYPE in FDK_audio.h.
                   - 2: MPEG-4 AAC Low Complexity.
                   - 5: MPEG-4 AAC Low Complexity with Spectral Band Replication
                 (HE-AAC).
                   - 29: MPEG-4 AAC Low Complexity with Spectral Band
                 Replication and Parametric Stereo (HE-AAC v2). This
                 configuration can be used only with stereo input audio data.
                   - 23: MPEG-4 AAC Low-Delay.
                   - 39: MPEG-4 AAC Enhanced Low-Delay. Since there is no
                 ::AUDIO_OBJECT_TYPE for ELD in combination with SBR defined,
                 enable SBR explicitely by ::AACENC_SBR_MODE parameter. The ELD
                 v2 212 configuration can be configured by ::AACENC_CHANNELMODE
                 parameter.
                   - 129: MPEG-2 AAC Low Complexity.
                   - 132: MPEG-2 AAC Low Complexity with Spectral Band
                 Replication (HE-AAC).

                   Please note that the virtual MPEG-2 AOT's basically disables
                 non-existing Perceptual Noise Substitution tool in AAC encoder
                 and controls the MPEG_ID flag in adts header. The virtual
                 MPEG-2 AOT doesn't prohibit specific transport formats. */

  AACENC_BITRATE = 0x0101, /*!< Total encoder bitrate. This parameter is
                              mandatory and interacts with ::AACENC_BITRATEMODE.
                                - CBR: Bitrate in bits/second.
                                - VBR: Variable bitrate. Bitrate argument will
                              be ignored. See \ref suppBitrates for details. */

  AACENC_BITRATEMODE = 0x0102, /*!< Bitrate mode. Configuration can be different
                                  kind of bitrate configurations:
                                    - 0: Constant bitrate, use bitrate according
                                  to ::AACENC_BITRATE. (default) Within none
                                  LD/ELD ::AUDIO_OBJECT_TYPE, the CBR mode makes
                                  use of full allowed bitreservoir. In contrast,
                                  at Low-Delay ::AUDIO_OBJECT_TYPE the
                                  bitreservoir is kept very small.
                                    - 1: Variable bitrate mode, \ref vbrmode
                                  "very low bitrate".
                                    - 2: Variable bitrate mode, \ref vbrmode
                                  "low bitrate".
                                    - 3: Variable bitrate mode, \ref vbrmode
                                  "medium bitrate".
                                    - 4: Variable bitrate mode, \ref vbrmode
                                  "high bitrate".
                                    - 5: Variable bitrate mode, \ref vbrmode
                                  "very high bitrate". */

  AACENC_SAMPLERATE = 0x0103, /*!< Audio input data sampling rate. Encoder
                                 supports following sampling rates: 8000, 11025,
                                 12000, 16000, 22050, 24000, 32000, 44100,
                                 48000, 64000, 88200, 96000 */

  AACENC_SBR_MODE = 0x0104, /*!< Configure SBR independently of the chosen Audio
                               Object Type ::AUDIO_OBJECT_TYPE. This parameter
                               is for ELD audio object type only.
                                 - -1: Use ELD SBR auto configurator (default).
                                 - 0: Disable Spectral Band Replication.
                                 - 1: Enable Spectral Band Replication. */

  AACENC_GRANULE_LENGTH =
      0x0105, /*!< Core encoder (AAC) audio frame length in samples:
                   - 1024: Default configuration.
                   - 512: Default length in LD/ELD configuration.
                   - 480: Length in LD/ELD configuration.
                   - 256: Length for ELD reduced delay mode (x2).
                   - 240: Length for ELD reduced delay mode (x2).
                   - 128: Length for ELD reduced delay mode (x4).
                   - 120: Length for ELD reduced delay mode (x4). */

  AACENC_CHANNELMODE = 0x0106, /*!< Set explicit channel mode. Channel mode must
                                  match with number of input channels.
                                    - 1-7, 11,12,14 and 33,34: MPEG channel
                                  modes supported, see ::CHANNEL_MODE in
                                  FDK_audio.h. */

  AACENC_CHANNELORDER =
      0x0107, /*!< Input audio data channel ordering scheme:
                   - 0: MPEG channel ordering (e. g. 5.1: C, L, R, SL, SR, LFE).
                 (default)
                   - 1: WAVE file format channel ordering (e. g. 5.1: L, R, C,
                 LFE, SL, SR). */

  AACENC_SBR_RATIO =
      0x0108, /*!<  Controls activation of downsampled SBR. With downsampled
                 SBR, the delay will be shorter. On the other hand, for
                 achieving the same quality level, downsampled SBR needs more
                 bits than dual-rate SBR. With downsampled SBR, the AAC encoder
                 will work at the same sampling rate as the SBR encoder (single
                 rate). Downsampled SBR is supported for AAC-ELD and HE-AACv1.
                    - 1: Downsampled SBR (default for ELD).
                    - 2: Dual-rate SBR   (default for HE-AAC). */

  AACENC_AFTERBURNER =
      0x0200, /*!< This parameter controls the use of the afterburner feature.
                   The afterburner is a type of analysis by synthesis algorithm
                 which increases the audio quality but also the required
                 processing power. It is recommended to always activate this if
                 additional memory consumption and processing power consumption
                   is not a problem. If increased MHz and memory consumption are
                 an issue then the MHz and memory cost of this optional module
                 need to be evaluated against the improvement in audio quality
                 on a case by case basis.
                   - 0: Disable afterburner (default).
                   - 1: Enable afterburner. */

  AACENC_BANDWIDTH = 0x0203, /*!< Core encoder audio bandwidth:
                                  - 0: Determine audio bandwidth internally
                                (default, see chapter \ref BEHAVIOUR_BANDWIDTH).
                                  - 1 to fs/2: Audio bandwidth in Hertz. Limited
                                to 20kHz max. Not usable if SBR is active. This
                                setting is for experts only, better do not touch
                                this value to avoid degraded audio quality. */

  AACENC_PEAK_BITRATE =
      0x0207, /*!< Peak bitrate configuration parameter to adjust maximum bits
                 per audio frame. Bitrate is in bits/second. The peak bitrate
                 will internally be limited to the chosen bitrate
                 ::AACENC_BITRATE as lower limit and the
                 number_of_effective_channels*6144 bit as upper limit.

                   Setting the peak bitrate equal to ::AACENC_BITRATE does not
                 necessarily mean that the audio frames will be of constant
                 size. Since the peak bitate is in bits/second, the frame sizes
                 can vary by one byte in one or the other direction over various
                 frames. However, it is not recommended to reduce the peak
                 pitrate to ::AACENC_BITRATE - it would disable the
                 bitreservoir, which would affect the audio quality by a large
                 amount. */

  AACENC_TRANSMUX = 0x0300, /*!< Transport type to be used. See ::TRANSPORT_TYPE
                               in FDK_audio.h. Following types can be configured
                               in encoder library:
                                 - 0: raw access units
                                 - 1: ADIF bitstream format
                                 - 2: ADTS bitstream format
                                 - 6: Audio Mux Elements (LATM) with
                               muxConfigPresent = 1
                                 - 7: Audio Mux Elements (LATM) with
                               muxConfigPresent = 0, out of band StreamMuxConfig
                                 - 10: Audio Sync Stream (LOAS) */

  AACENC_HEADER_PERIOD =
      0x0301, /*!< Frame count period for sending in-band configuration buffers
                 within LATM/LOAS transport layer. Additionally this parameter
                 configures the PCE repetition period in raw_data_block(). See
                 \ref encPCE.
                   - 0xFF: auto-mode default 10 for TT_MP4_ADTS, TT_MP4_LOAS and
                 TT_MP4_LATM_MCP1, otherwise 0.
                   - n: Frame count period. */

  AACENC_SIGNALING_MODE =
      0x0302, /*!< Signaling mode of the extension AOT:
                   - 0: Implicit backward compatible signaling (default for
                 non-MPEG-4 based AOT's and for the transport formats ADIF and
                 ADTS)
                        - A stream that uses implicit signaling can be decoded
                 by every AAC decoder, even AAC-LC-only decoders
                        - An AAC-LC-only decoder will only decode the
                 low-frequency part of the stream, resulting in a band-limited
                 output
                        - This method works with all transport formats
                        - This method does not work with downsampled SBR
                   - 1: Explicit backward compatible signaling
                        - A stream that uses explicit backward compatible
                 signaling can be decoded by every AAC decoder, even AAC-LC-only
                 decoders
                        - An AAC-LC-only decoder will only decode the
                 low-frequency part of the stream, resulting in a band-limited
                 output
                        - A decoder not capable of decoding PS will only decode
                 the AAC-LC+SBR part. If the stream contained PS, the result
                 will be a a decoded mono downmix
                        - This method does not work with ADIF or ADTS. For
                 LOAS/LATM, it only works with AudioMuxVersion==1
                        - This method does work with downsampled SBR
                   - 2: Explicit hierarchical signaling (default for MPEG-4
                 based AOT's and for all transport formats excluding ADIF and
                 ADTS)
                        - A stream that uses explicit hierarchical signaling can
                 be decoded only by HE-AAC decoders
                        - An AAC-LC-only decoder will not decode a stream that
                 uses explicit hierarchical signaling
                        - A decoder not capable of decoding PS will not decode
                 the stream at all if it contained PS
                        - This method does not work with ADIF or ADTS. It works
                 with LOAS/LATM and the MPEG-4 File format
                        - This method does work with downsampled SBR

                    For making sure that the listener always experiences the
                 best audio quality, explicit hierarchical signaling should be
                 used. This makes sure that only a full HE-AAC-capable decoder
                 will decode those streams. The audio is played at full
                 bandwidth. For best backwards compatibility, it is recommended
                 to encode with implicit SBR signaling. A decoder capable of
                 AAC-LC only will then only decode the AAC part, which means the
                 decoded audio will sound band-limited.

                    For MPEG-2 transport types (ADTS,ADIF), only implicit
                 signaling is possible.

                    For LOAS and LATM, explicit backwards compatible signaling
                 only works together with AudioMuxVersion==1. The reason is
                 that, for explicit backwards compatible signaling, additional
                 information will be appended to the ASC. A decoder that is only
                 capable of decoding AAC-LC will skip this part. Nevertheless,
                 for jumping to the end of the ASC, it needs to know the ASC
                 length. Transmitting the length of the ASC is a feature of
                 AudioMuxVersion==1, it is not possible to transmit the length
                 of the ASC with AudioMuxVersion==0, therefore an AAC-LC-only
                 decoder will not be able to parse a LOAS/LATM stream that was
                 being encoded with AudioMuxVersion==0.

                    For downsampled SBR, explicit signaling is mandatory. The
                 reason for this is that the extension sampling frequency (which
                 is in case of SBR the sampling frequqncy of the SBR part) can
                 only be signaled in explicit mode.

                    For AAC-ELD, the SBR information is transmitted in the
                 ELDSpecific Config, which is part of the AudioSpecificConfig.
                 Therefore, the settings here will have no effect on AAC-ELD.*/

  AACENC_TPSUBFRAMES =
      0x0303, /*!< Number of sub frames in a transport frame for LOAS/LATM or
                 ADTS (default 1).
                   - ADTS: Maximum number of sub frames restricted to 4.
                   - LOAS/LATM: Maximum number of sub frames restricted to 2.*/

  AACENC_AUDIOMUXVER =
      0x0304, /*!< AudioMuxVersion to be used for LATM. (AudioMuxVersionA,
                 currently not implemented):
                   - 0: Default, no transmission of tara Buffer fullness, no ASC
                 length and including actual latm Buffer fullnes.
                   - 1: Transmission of tara Buffer fullness, ASC length and
                 actual latm Buffer fullness.
                   - 2: Transmission of tara Buffer fullness, ASC length and
                 maximum level of latm Buffer fullness. */

  AACENC_PROTECTION = 0x0306, /*!< Configure protection in transport layer:
                                   - 0: No protection. (default)
                                   - 1: CRC active for ADTS transport format. */

  AACENC_ANCILLARY_BITRATE =
      0x0500, /*!< Constant ancillary data bitrate in bits/second.
                   - 0: Either no ancillary data or insert exact number of
                 bytes, denoted via input parameter, numAncBytes in
                 AACENC_InArgs.
                   - else: Insert ancillary data with specified bitrate. */

  AACENC_METADATA_MODE = 0x0600, /*!< Configure Meta Data. See ::AACENC_MetaData
                                    for further details:
                                      - 0: Do not embed any metadata.
                                      - 1: Embed dynamic_range_info metadata.
                                      - 2: Embed dynamic_range_info and
                                    ancillary_data metadata.
                                      - 3: Embed ancillary_data metadata. */

  AACENC_CONTROL_STATE =
      0xFF00, /*!< There is an automatic process which internally reconfigures
                 the encoder instance when a configuration parameter changed or
                 an error occured. This paramerter allows overwriting or getting
                 the control status of this process. See ::AACENC_CTRLFLAGS. */

  AACENC_NONE = 0xFFFF /*!< ------ */

} AACENC_PARAM;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief  Open an instance of the encoder.
 *
 * Allocate memory for an encoder instance with a functional range denoted by
 * the function parameters. Preinitialize encoder instance with default
 * configuration.
 *
 * \param phAacEncoder  A pointer to an encoder handle. Initialized on return.
 * \param encModules    Specify encoder modules to be supported in this encoder
 * instance:
 *                      - 0x0: Allocate memory for all available encoder
 * modules.
 *                      - else: Select memory allocation regarding encoder
 * modules. Following flags are possible and can be combined.
 *                              - 0x01: AAC module.
 *                              - 0x02: SBR module.
 *                              - 0x04: PS module.
 *                              - 0x08: MPS module.
 *                              - 0x10: Metadata module.
 *                              - example: (0x01|0x02|0x04|0x08|0x10) allocates
 * all modules and is equivalent to default configuration denotet by 0x0.
 * \param maxChannels   Number of channels to be allocated. This parameter can
 * be used in different ways:
 *                      - 0: Allocate maximum number of AAC and SBR channels as
 * supported by the library.
 *                      - nChannels: Use same maximum number of channels for
 * allocating memory in AAC and SBR module.
 *                      - nChannels | (nSbrCh<<8): Number of SBR channels can be
 * different to AAC channels to save data memory.
 *
 * \return
 *          - AACENC_OK, on succes.
 *          - AACENC_INVALID_HANDLE, AACENC_MEMORY_ERROR, AACENC_INVALID_CONFIG,
 * on failure.
 */
LINKSPEC(AACENC_ERROR) aacEncOpen(HANDLE_AACENCODER *phAacEncoder, const UINT encModules,
                        const UINT maxChannels);

/**
 * \brief  Close the encoder instance.
 *
 * Deallocate encoder instance and free whole memory.
 *
 * \param phAacEncoder  Pointer to the encoder handle to be deallocated.
 *
 * \return
 *          - AACENC_OK, on success.
 *          - AACENC_INVALID_HANDLE, on failure.
 */
LINKSPEC(AACENC_ERROR) aacEncClose(HANDLE_AACENCODER *phAacEncoder);

/**
 * \brief Encode audio data.
 *
 * This function is mainly for encoding audio data. In addition the function can
 * be used for an encoder (re)configuration process.
 * - PCM input data will be retrieved from external input buffer until the fill
 * level allows encoding a single frame. This functionality allows an external
 * buffer with reduced size in comparison to the AAC or HE-AAC audio frame
 * length.
 * - If the value of the input samples argument is zero, just internal
 * reinitialization will be applied if it is requested.
 * - At the end of a file the flushing process can be triggerd via setting the
 * value of the input samples argument to -1. The encoder delay lines are fully
 * flushed when the encoder returns no valid bitstream data
 * AACENC_OutArgs::numOutBytes. Furthermore the end of file is signaled by the
 * return value AACENC_ENCODE_EOF.
 * - If an error occured in the previous frame or any of the encoder parameters
 * changed, an internal reinitialization process will be applied before encoding
 * the incoming audio samples.
 * - The function can also be used for an independent reconfiguration process
 * without encoding. The first parameter has to be a valid encoder handle and
 * all other parameters can be set to NULL.
 * - If the size of the external bitbuffer in outBufDesc is not sufficient for
 * writing the whole bitstream, an internal error will be the return value and a
 * reconfiguration will be triggered.
 *
 * \param hAacEncoder           A valid AAC encoder handle.
 * \param inBufDesc             Input buffer descriptor, see AACENC_BufDesc:
 *                              - At least one input buffer with audio data is
 * expected.
 *                              - Optionally a second input buffer with
 * ancillary data can be fed.
 * \param outBufDesc            Output buffer descriptor, see AACENC_BufDesc:
 *                              - Provide one output buffer for the encoded
 * bitstream.
 * \param inargs                Input arguments, see AACENC_InArgs.
 * \param outargs               Output arguments, AACENC_OutArgs.
 *
 * \return
 *          - AACENC_OK, on success.
 *          - AACENC_INVALID_HANDLE, AACENC_ENCODE_ERROR, on failure in encoding
 * process.
 *          - AACENC_INVALID_CONFIG, AACENC_INIT_ERROR, AACENC_INIT_AAC_ERROR,
 * AACENC_INIT_SBR_ERROR, AACENC_INIT_TP_ERROR, AACENC_INIT_META_ERROR,
 * AACENC_INIT_MPS_ERROR, on failure in encoder initialization.
 *          - AACENC_UNSUPPORTED_PARAMETER, on incorrect input or output buffer
 * descriptor initialization.
 *          - AACENC_ENCODE_EOF, when flushing fully concluded.
 */
LINKSPEC(AACENC_ERROR) aacEncEncode(const HANDLE_AACENCODER hAacEncoder,
                          const AACENC_BufDesc *inBufDesc,
                          const AACENC_BufDesc *outBufDesc,
                          const AACENC_InArgs *inargs, AACENC_OutArgs *outargs);

/**
 * \brief  Acquire info about present encoder instance.
 *
 * This function retrieves information of the encoder configuration. In addition
 * to informative internal states, a configuration data block of the current
 * encoder settings will be returned. The format is either Audio Specific Config
 * in case of Raw Packets transport format or StreamMuxConfig in case of
 * LOAS/LATM transport format. The configuration data block is binary coded as
 * specified in ISO/IEC 14496-3 (MPEG-4 audio), to be used directly for MPEG-4
 * File Format or RFC3016 or RFC3640 applications.
 *
 * \param hAacEncoder           A valid AAC encoder handle.
 * \param pInfo                 Pointer to AACENC_InfoStruct. Filled on return.
 *
 * \return
 *          - AACENC_OK, on succes.
 *          - AACENC_INVALID_HANDLE, AACENC_INIT_ERROR, on failure.
 */
LINKSPEC(AACENC_ERROR) aacEncInfo(const HANDLE_AACENCODER hAacEncoder,
                        AACENC_InfoStruct *pInfo);

/**
 * \brief  Set one single AAC encoder parameter.
 *
 * This function allows configuration of all encoder parameters specified in
 * ::AACENC_PARAM. Each parameter must be set with a separate function call. An
 * internal validation of the configuration value range will be done and an
 * internal reconfiguration will be signaled. The actual configuration adoption
 * is part of the subsequent aacEncEncode() call.
 *
 * \param hAacEncoder           A valid AAC encoder handle.
 * \param param                 Parameter to be set. See ::AACENC_PARAM.
 * \param value                 Parameter value. See parameter description in
 * ::AACENC_PARAM.
 *
 * \return
 *          - AACENC_OK, on success.
 *          - AACENC_INVALID_HANDLE, AACENC_UNSUPPORTED_PARAMETER,
 * AACENC_INVALID_CONFIG, on failure.
 */
LINKSPEC(AACENC_ERROR) aacEncoder_SetParam(const HANDLE_AACENCODER hAacEncoder,
                                 const AACENC_PARAM param, const UINT value);

/**
 * \brief  Get one single AAC encoder parameter.
 *
 * This function is the complement to aacEncoder_SetParam(). After encoder
 * reinitialization with user defined settings, the internal status can be
 * obtained of each parameter, specified with ::AACENC_PARAM.
 *
 * \param hAacEncoder           A valid AAC encoder handle.
 * \param param                 Parameter to be returned. See ::AACENC_PARAM.
 *
 * \return  Internal configuration value of specifed parameter ::AACENC_PARAM.
 */
LINKSPEC(UINT) aacEncoder_GetParam(const HANDLE_AACENCODER hAacEncoder,
                         const AACENC_PARAM param);

/**
 * \brief  Get information about encoder library build.
 *
 * Fill a given LIB_INFO structure with library version information.
 *
 * \param info  Pointer to an allocated LIB_INFO struct.
 *
 * \return
 *          - AACENC_OK, on success.
 *          - AACENC_INVALID_HANDLE, AACENC_INIT_ERROR, on failure.
 */
LINKSPEC(AACENC_ERROR) aacEncGetLibInfo(LIB_INFO *info);

#ifdef __cplusplus
}
#endif

#endif /* AACENC_LIB_H */
