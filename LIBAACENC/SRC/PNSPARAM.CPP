/**************************** AAC encoder library ******************************

   Author(s):   M.Lohwasser

   Description: PNS parameters depending on bitrate and bandwidth

*******************************************************************************/

#include "pnsparam.h"

#include "psy_configuration.h"

typedef struct {
  SHORT startFreq;
  FIXP_SGL refPower;
  FIXP_SGL refTonality;
  SHORT tnsGainThreshold;
  SHORT tnsPNSGainThreshold;
  FIXP_SGL gapFillThr;
  SHORT minSfbWidth;
  USHORT detectionAlgorithmFlags;
} PNS_INFO_TAB;

typedef struct {
  ULONG brFrom;
  ULONG brTo;
  UCHAR S16000;
  UCHAR S22050;
  UCHAR S24000;
  UCHAR S32000;
  UCHAR S44100;
  UCHAR S48000;
} AUTO_PNS_TAB;

static const AUTO_PNS_TAB levelTable_mono[] = {
    {
        0,
        11999,
        0,
        1,
        1,
        1,
        1,
        1,
    },
    {
        12000,
        19999,
        0,
        1,
        1,
        1,
        1,
        1,
    },
    {
        20000,
        28999,
        0,
        2,
        1,
        1,
        1,
        1,
    },
    {
        29000,
        40999,
        0,
        4,
        4,
        4,
        2,
        2,
    },
    {
        41000,
        55999,
        0,
        9,
        9,
        7,
        7,
        7,
    },
    {
        56000,
        61999,
        0,
        0,
        0,
        0,
        9,
        9,
    },
    {
        62000,
        75999,
        0,
        0,
        0,
        0,
        0,
        0,
    },
    {
        76000,
        92999,
        0,
        0,
        0,
        0,
        0,
        0,
    },
    {
        93000,
        999999,
        0,
        0,
        0,
        0,
        0,
        0,
    },
};

static const AUTO_PNS_TAB levelTable_stereo[] = {
    {
        0,
        11999,
        0,
        1,
        1,
        1,
        1,
        1,
    },
    {
        12000,
        19999,
        0,
        3,
        1,
        1,
        1,
        1,
    },
    {
        20000,
        28999,
        0,
        3,
        3,
        3,
        2,
        2,
    },
    {
        29000,
        40999,
        0,
        7,
        6,
        6,
        5,
        5,
    },
    {
        41000,
        55999,
        0,
        9,
        9,
        7,
        7,
        7,
    },
    {
        56000,
        79999,
        0,
        0,
        0,
        0,
        0,
        0,
    },
    {
        80000,
        99999,
        0,
        0,
        0,
        0,
        0,
        0,
    },
    {
        100000,
        999999,
        0,
        0,
        0,
        0,
        0,
        0,
    },
};

static const PNS_INFO_TAB pnsInfoTab[] = {
    {4000, FL2FXCONST_SGL(0.04), FL2FXCONST_SGL(0.06), 1150, 1200,
           FL2FXCONST_SGL(0.02), 8,
           USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
               USE_TNS_PNS},
    {4000, FL2FXCONST_SGL(0.04), FL2FXCONST_SGL(0.07), 1130, 1300,
     FL2FXCONST_SGL(0.05), 8,
     USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
         USE_TNS_PNS},
    {4100, FL2FXCONST_SGL(0.04), FL2FXCONST_SGL(0.07), 1100, 1400,
     FL2FXCONST_SGL(0.10), 8,
     USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
         USE_TNS_PNS},
    {4100, FL2FXCONST_SGL(0.03), FL2FXCONST_SGL(0.10), 1100, 1400,
     FL2FXCONST_SGL(0.15), 8,
     USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
         USE_TNS_PNS},
    {4300, FL2FXCONST_SGL(0.03), FL2FXCONST_SGL(0.10), 1100, 1400,
     FL2FXCONST_SGL(0.15), 8,
     USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
         USE_TNS_PNS | JUST_LONG_WINDOW},
    {5000, FL2FXCONST_SGL(0.03), FL2FXCONST_SGL(0.10), 1100, 1400,
     FL2FXCONST_SGL(0.25), 8,
     USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
         USE_TNS_PNS | JUST_LONG_WINDOW},
    {5500, FL2FXCONST_SGL(0.03), FL2FXCONST_SGL(0.12), 1100, 1400,
     FL2FXCONST_SGL(0.35), 8,
     USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
         USE_TNS_PNS | JUST_LONG_WINDOW},
    {6000, FL2FXCONST_SGL(0.03), FL2FXCONST_SGL(0.12), 1080, 1400,
     FL2FXCONST_SGL(0.40), 8,
     USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
         USE_TNS_PNS | JUST_LONG_WINDOW},
    {6000, FL2FXCONST_SGL(0.03), FL2FXCONST_SGL(0.14), 1070, 1400,
     FL2FXCONST_SGL(0.45), 8,
     USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
         USE_TNS_PNS | JUST_LONG_WINDOW},
};

static const AUTO_PNS_TAB levelTable_lowComplexity[] = {
    {
        0,
        27999,
        0,
        0,
        0,
        0,
        0,
        0,
    },
    {
        28000,
        31999,
        0,
        2,
        2,
        2,
        2,
        2,
    },
    {
        32000,
        47999,
        0,
        3,
        3,
        3,
        3,
        3,
    },
    {
        48000,
        48000,
        0,
        4,
        4,
        4,
        4,
        4,
    },
    {
        48001,
        999999,
        0,
        0,
        0,
        0,
        0,
        0,
    },
};
static const PNS_INFO_TAB pnsInfoTab_lowComplexity[] = {
    {4100, FL2FXCONST_SGL(0.03), FL2FXCONST_SGL(0.16), 1100, 1400,
           FL2FXCONST_SGL(0.5), 16,
           USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
               USE_TNS_PNS | JUST_LONG_WINDOW},
    {4100, FL2FXCONST_SGL(0.05), FL2FXCONST_SGL(0.10), 1410, 1400,
     FL2FXCONST_SGL(0.5), 16,
     USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
         USE_TNS_PNS | JUST_LONG_WINDOW},
    {4100, FL2FXCONST_SGL(0.05), FL2FXCONST_SGL(0.10), 1100, 1400,
     FL2FXCONST_SGL(0.5), 16,
     USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
         USE_TNS_PNS | JUST_LONG_WINDOW},
    {4100, FL2FXCONST_SGL(0.20), FL2FXCONST_SGL(0.10), 1410, 1400,
     FL2FXCONST_SGL(0.5), 16,
     USE_POWER_DISTRIBUTION | USE_PSYCH_TONALITY | USE_TNS_GAIN_THR |
         USE_TNS_PNS | JUST_LONG_WINDOW},
};

int FDKaacEnc_lookUpPnsUse(int bitRate, int sampleRate, int numChan,
                           const int isLC) {
  int hUsePns = 0, size, i;
  const AUTO_PNS_TAB *levelTable;

  if (isLC) {
    levelTable = &levelTable_lowComplexity[0];
    size = sizeof(levelTable_lowComplexity);
  } else {
    levelTable = (numChan > 1) ? &levelTable_stereo[0] : &levelTable_mono[0];
    size = (numChan > 1) ? sizeof(levelTable_stereo) : sizeof(levelTable_mono);
  }

  for (i = 0; i < (int)(size / sizeof(AUTO_PNS_TAB)); i++) {
    if (((ULONG)bitRate >= levelTable[i].brFrom) &&
        ((ULONG)bitRate <= levelTable[i].brTo))
      break;
  }

  if ((int)(sizeof(pnsInfoTab) / sizeof(PNS_INFO_TAB)) < i) {
    return (PNS_TABLE_ERROR);
  }

  switch (sampleRate) {
    case 16000:
      hUsePns = levelTable[i].S16000;
      break;
    case 22050:
      hUsePns = levelTable[i].S22050;
      break;
    case 24000:
      hUsePns = levelTable[i].S24000;
      break;
    case 32000:
      hUsePns = levelTable[i].S32000;
      break;
    case 44100:
      hUsePns = levelTable[i].S44100;
      break;
    case 48000:
      hUsePns = levelTable[i].S48000;
      break;
    default:
      if (isLC) {
        hUsePns = levelTable[i].S48000;
      }
      break;
  }

  return (hUsePns);
}

AAC_ENCODER_ERROR FDKaacEnc_GetPnsParam(NOISEPARAMS *np, INT bitRate,
                                        INT sampleRate, INT sfbCnt,
                                        const INT *sfbOffset, INT *usePns,
                                        INT numChan, const INT isLC) {
  int i, hUsePns;
  const PNS_INFO_TAB *pnsInfo;

  if (*usePns <= 0) return AAC_ENC_OK;

  if (isLC) {
    np->detectionAlgorithmFlags = IS_LOW_COMPLEXITY;

    pnsInfo = pnsInfoTab_lowComplexity;

    hUsePns = FDKaacEnc_lookUpPnsUse(bitRate, sampleRate, numChan, isLC);
    if (hUsePns == 0) {
      *usePns = 0;
      return AAC_ENC_OK;
    }

    if (hUsePns == PNS_TABLE_ERROR) {
      return AAC_ENC_PNS_TABLE_ERROR;
    }

    pnsInfo += hUsePns - 1;

  } else {
    np->detectionAlgorithmFlags = 0;
    pnsInfo = pnsInfoTab;

    hUsePns = FDKaacEnc_lookUpPnsUse(bitRate, sampleRate, numChan, isLC);
    if (hUsePns == 0) {
      *usePns = 0;
      return AAC_ENC_OK;
    }
    if (hUsePns == PNS_TABLE_ERROR) return AAC_ENC_PNS_TABLE_ERROR;

    pnsInfo += hUsePns - 1;
  }

  np->startSfb = FDKaacEnc_FreqToBandWidthRounding(
      pnsInfo->startFreq, sampleRate, sfbCnt, sfbOffset);

  np->detectionAlgorithmFlags |= pnsInfo->detectionAlgorithmFlags;

  np->refPower = FX_SGL2FX_DBL(pnsInfo->refPower);
  np->refTonality = FX_SGL2FX_DBL(pnsInfo->refTonality);
  np->tnsGainThreshold = pnsInfo->tnsGainThreshold;
  np->tnsPNSGainThreshold = pnsInfo->tnsPNSGainThreshold;
  np->minSfbWidth = pnsInfo->minSfbWidth;

  np->gapFillThr =
      pnsInfo->gapFillThr;

  for (i = 0; i < (sfbCnt - 1); i++) {
    INT qtmp, sfbWidth;
    FIXP_DBL tmp;

    sfbWidth = sfbOffset[i + 1] - sfbOffset[i];

    tmp = fPow(np->refPower, 0, sfbWidth, DFRACT_BITS - 1 - 5, &qtmp);
    np->powDistPSDcurve[i] = (FIXP_SGL)((LONG)(scaleValue(tmp, qtmp) >> 16));
  }
  np->powDistPSDcurve[sfbCnt] = np->powDistPSDcurve[sfbCnt - 1];

  return AAC_ENC_OK;
}
