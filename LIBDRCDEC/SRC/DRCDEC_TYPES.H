/************************* MPEG-D DRC decoder library **************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef DRCDEC_TYPES_H
#define DRCDEC_TYPES_H

#include "common_fix.h"

typedef struct {
  FIXP_SGL gainDb;
  SHORT time;
} GAIN_NODE;

typedef struct {
  UCHAR uniDrcGainExtType[8];
  ULONG extBitSize[8 - 1];
} UNI_DRC_GAIN_EXTENSION;

typedef struct {
  UCHAR nNodes[12];
  GAIN_NODE gainNode[12][16];

  UCHAR uniDrcGainExtPresent;
  UNI_DRC_GAIN_EXTENSION uniDrcGainExtension;

  UCHAR status;
} UNI_DRC_GAIN, *HANDLE_UNI_DRC_GAIN;

typedef enum {
  EB_NIGHT = 0x0001,
  EB_NOISY = 0x0002,
  EB_LIMITED = 0x0004,
  EB_LOWLEVEL = 0x0008,
  EB_DIALOG = 0x0010,
  EB_GENERAL_COMPR = 0x0020,
  EB_EXPAND = 0x0040,
  EB_ARTISTIC = 0x0080,
  EB_CLIPPING = 0x0100,
  EB_FADE = 0x0200,
  EB_DUCK_OTHER = 0x0400,
  EB_DUCK_SELF = 0x0800
} EFFECT_BIT;

typedef enum {
  GCP_REGULAR = 0,
  GCP_FADING = 1,
  GCP_CLIPPING_DUCKING = 2,
  GCP_CONSTANT = 3
} GAIN_CODING_PROFILE;

typedef enum { GIT_SPLINE = 0, GIT_LINEAR = 1 } GAIN_INTERPOLATION_TYPE;

typedef enum { CS_LEFT = 0, CS_RIGHT = 1 } CHARACTERISTIC_SIDE;

typedef enum { CF_SIGMOID = 0, CF_NODES = 1 } CHARACTERISTIC_FORMAT;

typedef enum {
  GF_QMF32 = 0x1,
  GF_QMFHYBRID39 = 0x2,
  GF_QMF64 = 0x3,
  GF_QMFHYBRID71 = 0x4,
  GF_QMF128 = 0x5,
  GF_QMFHYBRID135 = 0x6,
  GF_UNIFORM = 0x7
} EQ_SUBBAND_GAIN_FORMAT;

typedef struct {
  UCHAR duckingScalingPresent;
  FIXP_SGL duckingScaling;
} DUCKING_MODIFICATION;

typedef struct {
  UCHAR targetCharacteristicLeftPresent;
  UCHAR targetCharacteristicLeftIndex;
  UCHAR targetCharacteristicRightPresent;
  UCHAR targetCharacteristicRightIndex;
  UCHAR gainScalingPresent;
  FIXP_SGL attenuationScaling;
  FIXP_SGL amplificationScaling;
  UCHAR gainOffsetPresent;
  FIXP_SGL gainOffset;
} GAIN_MODIFICATION;

typedef union {
  UCHAR crossoverFreqIndex;
  USHORT startSubBandIndex;
} BAND_BORDER;

typedef struct {
  UCHAR left;
  UCHAR right;
} CUSTOM_INDEX;

typedef struct {
  UCHAR present;
  UCHAR isCICP;
  union {
    UCHAR cicpIndex;
    CUSTOM_INDEX custom;
  };
} DRC_CHARACTERISTIC;

typedef struct {
  UCHAR gainCodingProfile;
  UCHAR gainInterpolationType;
  UCHAR fullFrame;
  UCHAR timeAlignment;
  UCHAR timeDeltaMinPresent;
  USHORT timeDeltaMin;
  UCHAR bandCount;
  UCHAR drcBandType;
  UCHAR gainSequenceIndex[4];
  DRC_CHARACTERISTIC drcCharacteristic[4];
  BAND_BORDER bandBorder[4];
} GAIN_SET;

typedef struct {
  FIXP_SGL gain;
  FIXP_SGL ioRatio;
  FIXP_SGL exp;
  UCHAR flipSign;
} CUSTOM_DRC_CHAR_SIGMOID;

typedef struct {
  UCHAR characteristicNodeCount;
  FIXP_SGL nodeLevel[4 + 1];
  FIXP_SGL nodeGain[4 + 1];
} CUSTOM_DRC_CHAR_NODES;

typedef shouldBeUnion {
  CUSTOM_DRC_CHAR_SIGMOID sigmoid;
  CUSTOM_DRC_CHAR_NODES nodes;
}
CUSTOM_DRC_CHAR;

typedef struct {
  UCHAR drcLocation;
  UCHAR drcFrameSizePresent;
  USHORT drcFrameSize;
  UCHAR characteristicLeftCount;
  UCHAR characteristicLeftFormat[16];
  CUSTOM_DRC_CHAR customCharacteristicLeft[16];
  UCHAR characteristicRightCount;
  UCHAR characteristicRightFormat[16];
  CUSTOM_DRC_CHAR customCharacteristicRight[16];
  UCHAR
  gainSequenceCount;
  UCHAR gainSetCount;
  GAIN_SET gainSet[12];
  UCHAR gainSetIndexForGainSequence[12];
} DRC_COEFFICIENTS_UNI_DRC;

typedef struct {
  SCHAR drcSetId;
  UCHAR drcSetComplexityLevel;
  UCHAR drcLocation;
  UCHAR drcApplyToDownmix;
  UCHAR downmixIdCount;
  UCHAR downmixId[8];
  USHORT drcSetEffect;
  UCHAR limiterPeakTargetPresent;
  FIXP_SGL limiterPeakTarget;
  UCHAR drcSetTargetLoudnessPresent;
  SCHAR drcSetTargetLoudnessValueUpper;
  SCHAR drcSetTargetLoudnessValueLower;
  UCHAR dependsOnDrcSetPresent;
  union {
    SCHAR dependsOnDrcSet;
    UCHAR noIndependentUse;
  };
  UCHAR requiresEq;
  shouldBeUnion {
    GAIN_MODIFICATION gainModificationForChannelGroup[8][4];
    DUCKING_MODIFICATION duckingModificationForChannel[8];
  };
  SCHAR gainSetIndex[8];

  UCHAR drcChannelCount;
  UCHAR nDrcChannelGroups;
  SCHAR gainSetIndexForChannelGroup[8];
} DRC_INSTRUCTIONS_UNI_DRC;

typedef struct {
  UCHAR baseChannelCount;
  UCHAR layoutSignalingPresent;
  UCHAR definedLayout;
  UCHAR speakerPosition[8];
} CHANNEL_LAYOUT;

typedef struct {
  UCHAR downmixId;
  UCHAR targetChannelCount;
  UCHAR targetLayout;
  UCHAR downmixCoefficientsPresent;
  UCHAR bsDownmixOffset;
  FIXP_DBL downmixCoefficient[8 * 8];
} DOWNMIX_INSTRUCTIONS;

typedef struct {
  UCHAR uniDrcConfigExtType[8];
  ULONG extBitSize[8 - 1];
} UNI_DRC_CONFIG_EXTENSION;

typedef struct {
  UCHAR sampleRatePresent;
  ULONG sampleRate;
  UCHAR downmixInstructionsCountV0;
  UCHAR downmixInstructionsCountV1;
  UCHAR downmixInstructionsCount;
  UCHAR drcCoefficientsUniDrcCountV0;
  UCHAR drcCoefficientsUniDrcCountV1;
  UCHAR drcCoefficientsUniDrcCount;
  UCHAR drcInstructionsUniDrcCountV0;
  UCHAR drcInstructionsUniDrcCountV1;
  UCHAR drcInstructionsUniDrcCount;
  CHANNEL_LAYOUT channelLayout;
  DOWNMIX_INSTRUCTIONS downmixInstructions[6];
  DRC_COEFFICIENTS_UNI_DRC drcCoefficientsUniDrc[2];
  DRC_INSTRUCTIONS_UNI_DRC drcInstructionsUniDrc[(12 + 1 + 6)];
  UCHAR uniDrcConfigExtPresent;
  UNI_DRC_CONFIG_EXTENSION uniDrcConfigExt;

  UCHAR drcInstructionsCountInclVirtual;
  UCHAR diff;
} UNI_DRC_CONFIG, *HANDLE_UNI_DRC_CONFIG;

typedef enum {
  MD_UNKNOWN_OTHER = 0,
  MD_PROGRAM_LOUDNESS = 1,
  MD_ANCHOR_LOUDNESS = 2,
  MD_MAX_OF_LOUDNESS_RANGE = 3,
  MD_MOMENTARY_LOUDNESS_MAX = 4,
  MD_SHORT_TERM_LOUDNESS_MAX = 5,
  MD_LOUDNESS_RANGE = 6,
  MD_MIXING_LEVEL = 7,
  MD_ROOM_TYPE = 8,
  MD_SHORT_TERM_LOUDNESS = 9
} METHOD_DEFINITION;

typedef enum {
  MS_UNKNOWN_OTHER = 0,
  MS_EBU_R_128 = 1,
  MS_BS_1770_4 = 2,
  MS_BS_1770_4_PRE_PROCESSING = 3,
  MS_USER = 4,
  MS_EXPERT_PANEL = 5,
  MS_BS_1771_1 = 6,
  MS_RESERVED_A = 7,
  MS_RESERVED_B = 8,
  MS_RESERVED_C = 9,
  MS_RESERVED_D = 10,
  MS_RESERVED_E = 11
} MEASUREMENT_SYSTEM;

typedef enum {
  R_UKNOWN = 0,
  R_UNVERIFIED = 1,
  R_CEILING = 2,
  R_ACCURATE = 3
} RELIABILITY;

typedef struct {
  UCHAR methodDefinition;
  FIXP_DBL methodValue;
  UCHAR measurementSystem;
  UCHAR reliability;
} LOUDNESS_MEASUREMENT;

typedef struct {
  SCHAR drcSetId;
  UCHAR eqSetId;
  UCHAR downmixId;
  UCHAR samplePeakLevelPresent;
  FIXP_DBL samplePeakLevel;
  UCHAR truePeakLevelPresent;
  FIXP_DBL truePeakLevel;
  UCHAR truePeakLevelMeasurementSystem;
  UCHAR truePeakLevelReliability;
  UCHAR measurementCount;
  LOUDNESS_MEASUREMENT loudnessMeasurement[8];
} LOUDNESS_INFO;

typedef struct {
  UCHAR loudnessInfoSetExtType[8];
  ULONG extBitSize[8 - 1];
} LOUDNESS_INFO_SET_EXTENSION;

typedef struct {
  UCHAR loudnessInfoAlbumCountV0;
  UCHAR loudnessInfoAlbumCountV1;
  UCHAR loudnessInfoAlbumCount;
  UCHAR loudnessInfoCountV0;
  UCHAR loudnessInfoCountV1;
  UCHAR loudnessInfoCount;
  LOUDNESS_INFO loudnessInfoAlbum[12];
  LOUDNESS_INFO loudnessInfo[12];
  UCHAR loudnessInfoSetExtPresent;
  LOUDNESS_INFO_SET_EXTENSION loudnessInfoSetExt;
  UCHAR diff;
} LOUDNESS_INFO_SET, *HANDLE_LOUDNESS_INFO_SET;

#endif
