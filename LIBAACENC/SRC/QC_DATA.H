/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Quantizing & coding data

*******************************************************************************/

#ifndef QC_DATA_H
#define QC_DATA_H

#include "aacenc.h"
#include "psy_const.h"
#include "dyn_bits.h"
#include "adj_thr_data.h"
#include "line_pe.h"
#include "FDK_audio.h"
#include "interface.h"

typedef enum {
  QCDATA_BR_MODE_INVALID = -1,
  QCDATA_BR_MODE_CBR = 0,
  QCDATA_BR_MODE_VBR_1 = 1,
  QCDATA_BR_MODE_VBR_2 = 2,
  QCDATA_BR_MODE_VBR_3 = 3,
  QCDATA_BR_MODE_VBR_4 = 4,
  QCDATA_BR_MODE_VBR_5 = 5,
  QCDATA_BR_MODE_FF = 6,
  QCDATA_BR_MODE_SFR = 7

} QCDATA_BR_MODE;

typedef struct {
  MP4_ELEMENT_ID elType;
  INT instanceTag;
  INT nChannelsInEl;
  INT ChannelIndex[2];
  FIXP_DBL relativeBits;
} ELEMENT_INFO;

typedef struct {
  CHANNEL_MODE encMode;
  INT nChannels;
  INT nChannelsEff;
  INT nElements;
  ELEMENT_INFO elInfo[((8))];
} CHANNEL_MAPPING;

typedef struct {
  INT paddingRest;
} PADDING;

struct QC_INIT {
  CHANNEL_MAPPING *channelMapping;
  INT sceCpe;
  INT maxBits;
  INT averageBits;
  INT bitRes;
  INT sampleRate;
  INT isLowDelay;
  INT staticBits;
  QCDATA_BR_MODE bitrateMode;
  INT meanPe;
  INT chBitrate;
  INT invQuant;
  INT maxIterations;
  FIXP_DBL maxBitFac;
  INT bitrate;
  INT nSubFrames;
  INT minBits;
  AACENC_BITRES_MODE bitResMode;
  INT bitDistributionMode;

  PADDING padding;
};

typedef struct {
  FIXP_DBL mdctSpectrum[(1024)];

  SHORT quantSpec[(1024)];

  UINT maxValueInSfb[MAX_GROUPED_SFB];
  INT scf[MAX_GROUPED_SFB];
  INT globalGain;
  SECTION_DATA sectionData;

  FIXP_DBL sfbFormFactorLdData[MAX_GROUPED_SFB];

  FIXP_DBL sfbThresholdLdData[MAX_GROUPED_SFB];
  FIXP_DBL sfbMinSnrLdData[MAX_GROUPED_SFB];
  FIXP_DBL sfbEnergyLdData[MAX_GROUPED_SFB];
  FIXP_DBL sfbEnergy[MAX_GROUPED_SFB];
  FIXP_DBL sfbWeightedEnergyLdData[MAX_GROUPED_SFB];

  FIXP_DBL sfbEnFacLd[MAX_GROUPED_SFB];

  FIXP_DBL sfbSpreadEnergy[MAX_GROUPED_SFB];

} QC_OUT_CHANNEL;

typedef struct {
  EXT_PAYLOAD_TYPE type;
  INT nPayloadBits;
  UCHAR *pPayload;

} QC_OUT_EXTENSION;

typedef struct {
  INT staticBitsUsed;
  INT dynBitsUsed;

  INT extBitsUsed;
  INT nExtensions;
  QC_OUT_EXTENSION extension[(1)];

  INT grantedDynBits;

  INT grantedPe;
  INT grantedPeCorr;

  PE_DATA peData;

  QC_OUT_CHANNEL *qcOutChannel[(2)];

  UCHAR
  *dynMem_Ah_Flag;
  UCHAR
  *dynMem_Thr_Exp;
  UCHAR *dynMem_SfbNActiveLinesLdData;

} QC_OUT_ELEMENT;

typedef struct {
  QC_OUT_ELEMENT *qcElement[((8))];
  QC_OUT_CHANNEL *pQcOutChannels[(8)];
  QC_OUT_EXTENSION extension[(2 + 2)];
  INT nExtensions;
  INT maxDynBits;
  INT grantedDynBits;
  INT totFillBits;
  INT elementExtBits;
  INT globalExtBits;
  INT staticBits;

  INT totalNoRedPe;
  INT totalGrantedPeCorr;

  INT usedDynBits;
  INT alignBits;
  INT totalBits;

} QC_OUT;

typedef struct {
  INT chBitrateEl;
  INT maxBitsEl;
  INT bitResLevelEl;
  INT maxBitResBitsEl;
  FIXP_DBL relativeBitsEl;
} ELEMENT_BITS;

typedef struct {

  INT globHdrBits;
  INT maxBitsPerFrame;
  INT minBitsPerFrame;
  INT nElements;
  QCDATA_BR_MODE bitrateMode;
  AACENC_BITRES_MODE bitResMode;
  INT bitResTot;
  INT bitResTotMax;
  INT maxIterations;
  INT invQuant;

  FIXP_DBL vbrQualFactor;
  FIXP_DBL maxBitFac;

  PADDING padding;

  ELEMENT_BITS *elementBits[((8))];
  BITCNTR_STATE *hBitCounter;
  ADJ_THR_STATE *hAdjThr;

  INT dZoneQuantEnable;

} QC_STATE;

#endif /* QC_DATA_H */
