/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Threshold compensation

*******************************************************************************/

#include "adj_thr.h"
#include "sf_estim.h"
#include "aacEnc_ram.h"

#define NUM_NRG_LEVS (8)
#define INV_INT_TAB_SIZE (8)
static const FIXP_DBL invInt[INV_INT_TAB_SIZE] = {
    0x7fffffff, 0x7fffffff, 0x40000000, 0x2aaaaaaa,
    0x20000000, 0x19999999, 0x15555555, 0x12492492};

#define INV_SQRT4_TAB_SIZE (8)
static const FIXP_DBL invSqrt4[INV_SQRT4_TAB_SIZE] = {
    0x7fffffff, 0x7fffffff, 0x6ba27e65, 0x61424bb5,
    0x5a827999, 0x55994845, 0x51c8e33c, 0x4eb160d1};

static const FIXP_DBL SnrLdMin1 =
    (FIXP_DBL)0xfcad0ddf;
static const FIXP_DBL SnrLdMin2 =
    (FIXP_DBL)0x0351e1a2;
static const FIXP_DBL SnrLdFac =
    (FIXP_DBL)0xff5b2c3e;

static const FIXP_DBL SnrLdMin3 =
    (FIXP_DBL)0xfe000000;
static const FIXP_DBL SnrLdMin4 =
    (FIXP_DBL)0x02000000;
static const FIXP_DBL SnrLdMin5 =
    (FIXP_DBL)0xfc000000;

#define AFTERBURNER_STATI 2
#define MAX_ALLOWED_EL_CHANNELS 2

typedef struct {
  INT bitrate;
  FIXP_DBL bits2PeFactor[AFTERBURNER_STATI][MAX_ALLOWED_EL_CHANNELS];
} BIT_PE_SFAC;

typedef struct {
  INT sampleRate;
  const BIT_PE_SFAC *pPeTab;
  INT nEntries;

} BITS2PE_CFG_TAB;

#define FL2B2PE(value) FL2FXCONST_DBL((value) / (1 << 2))

static const BIT_PE_SFAC S_Bits2PeTab16000[] = {
    {10000,
     {{FL2B2PE(1.60f), FL2B2PE(0.00f)}, {FL2B2PE(1.40f), FL2B2PE(0.00f)}}},
    {24000,
     {{FL2B2PE(1.80f), FL2B2PE(1.40f)}, {FL2B2PE(1.60f), FL2B2PE(1.20f)}}},
    {32000,
     {{FL2B2PE(1.80f), FL2B2PE(1.60f)}, {FL2B2PE(1.60f), FL2B2PE(1.40f)}}},
    {48000,
     {{FL2B2PE(1.60f), FL2B2PE(1.80f)}, {FL2B2PE(1.60f), FL2B2PE(1.60f)}}},
    {64000,
     {{FL2B2PE(1.20f), FL2B2PE(1.60f)}, {FL2B2PE(1.20f), FL2B2PE(1.60f)}}},
    {96000,
     {{FL2B2PE(1.40f), FL2B2PE(1.80f)}, {FL2B2PE(1.40f), FL2B2PE(1.60f)}}},
    {128000,
     {{FL2B2PE(1.40f), FL2B2PE(1.80f)}, {FL2B2PE(1.40f), FL2B2PE(1.80f)}}},
    {148000,
     {{FL2B2PE(1.40f), FL2B2PE(1.80f)}, {FL2B2PE(1.40f), FL2B2PE(1.40f)}}}};

static const BIT_PE_SFAC S_Bits2PeTab22050[] = {
    {16000,
     {{FL2B2PE(1.60f), FL2B2PE(1.40f)}, {FL2B2PE(1.20f), FL2B2PE(0.80f)}}},
    {24000,
     {{FL2B2PE(1.60f), FL2B2PE(1.40f)}, {FL2B2PE(1.40f), FL2B2PE(1.00f)}}},
    {32000,
     {{FL2B2PE(1.40f), FL2B2PE(1.40f)}, {FL2B2PE(1.40f), FL2B2PE(1.20f)}}},
    {48000,
     {{FL2B2PE(1.20f), FL2B2PE(1.60f)}, {FL2B2PE(1.20f), FL2B2PE(1.40f)}}},
    {64000,
     {{FL2B2PE(1.60f), FL2B2PE(1.60f)}, {FL2B2PE(1.60f), FL2B2PE(1.40f)}}},
    {96000,
     {{FL2B2PE(1.80f), FL2B2PE(1.60f)}, {FL2B2PE(1.80f), FL2B2PE(1.60f)}}},
    {128000,
     {{FL2B2PE(1.80f), FL2B2PE(1.80f)}, {FL2B2PE(1.60f), FL2B2PE(1.60f)}}},
    {148000,
     {{FL2B2PE(1.40f), FL2B2PE(1.80f)}, {FL2B2PE(1.40f), FL2B2PE(1.60f)}}}};

static const BIT_PE_SFAC S_Bits2PeTab24000[] = {
    {16000,
     {{FL2B2PE(1.40f), FL2B2PE(1.40f)}, {FL2B2PE(1.20f), FL2B2PE(0.80f)}}},
    {24000,
     {{FL2B2PE(1.60f), FL2B2PE(1.20f)}, {FL2B2PE(1.40f), FL2B2PE(1.00f)}}},
    {32000,
     {{FL2B2PE(1.40f), FL2B2PE(1.20f)}, {FL2B2PE(1.40f), FL2B2PE(0.80f)}}},
    {48000,
     {{FL2B2PE(1.40f), FL2B2PE(1.60f)}, {FL2B2PE(1.40f), FL2B2PE(1.40f)}}},
    {64000,
     {{FL2B2PE(1.60f), FL2B2PE(1.60f)}, {FL2B2PE(1.60f), FL2B2PE(1.40f)}}},
    {96000,
     {{FL2B2PE(1.80f), FL2B2PE(1.60f)}, {FL2B2PE(1.80f), FL2B2PE(1.60f)}}},
    {128000,
     {{FL2B2PE(1.40f), FL2B2PE(1.60f)}, {FL2B2PE(1.80f), FL2B2PE(1.80f)}}},
    {148000,
     {{FL2B2PE(1.40f), FL2B2PE(1.60f)}, {FL2B2PE(1.40f), FL2B2PE(1.80f)}}}};

static const BIT_PE_SFAC S_Bits2PeTab32000[] = {
    {16000,
     {{FL2B2PE(1.20f), FL2B2PE(1.40f)}, {FL2B2PE(0.80f), FL2B2PE(0.80f)}}},
    {24000,
     {{FL2B2PE(1.40f), FL2B2PE(1.20f)}, {FL2B2PE(1.00f), FL2B2PE(0.60f)}}},
    {32000,
     {{FL2B2PE(1.20f), FL2B2PE(1.20f)}, {FL2B2PE(1.00f), FL2B2PE(0.80f)}}},
    {48000,
     {{FL2B2PE(1.40f), FL2B2PE(1.40f)}, {FL2B2PE(1.20f), FL2B2PE(1.20f)}}},
    {64000,
     {{FL2B2PE(1.60f), FL2B2PE(1.40f)}, {FL2B2PE(1.60f), FL2B2PE(1.20f)}}},
    {96000,
     {{FL2B2PE(1.60f), FL2B2PE(1.40f)}, {FL2B2PE(1.60f), FL2B2PE(1.40f)}}},
    {128000,
     {{FL2B2PE(1.80f), FL2B2PE(1.60f)}, {FL2B2PE(1.80f), FL2B2PE(1.60f)}}},
    {148000,
     {{FL2B2PE(1.80f), FL2B2PE(1.60f)}, {FL2B2PE(1.80f), FL2B2PE(1.60f)}}},
    {160000,
     {{FL2B2PE(1.80f), FL2B2PE(1.60f)}, {FL2B2PE(1.80f), FL2B2PE(1.60f)}}},
    {200000,
     {{FL2B2PE(1.40f), FL2B2PE(1.60f)}, {FL2B2PE(1.40f), FL2B2PE(1.60f)}}},
    {320000,
     {{FL2B2PE(3.20f), FL2B2PE(1.80f)}, {FL2B2PE(3.20f), FL2B2PE(1.80f)}}}};

static const BIT_PE_SFAC S_Bits2PeTab44100[] = {
    {16000,
     {{FL2B2PE(1.20f), FL2B2PE(1.60f)}, {FL2B2PE(0.80f), FL2B2PE(1.00f)}}},
    {24000,
     {{FL2B2PE(1.00f), FL2B2PE(1.20f)}, {FL2B2PE(1.00f), FL2B2PE(0.80f)}}},
    {32000,
     {{FL2B2PE(1.20f), FL2B2PE(1.20f)}, {FL2B2PE(0.80f), FL2B2PE(0.60f)}}},
    {48000,
     {{FL2B2PE(1.20f), FL2B2PE(1.20f)}, {FL2B2PE(1.20f), FL2B2PE(0.80f)}}},
    {64000,
     {{FL2B2PE(1.40f), FL2B2PE(1.20f)}, {FL2B2PE(1.20f), FL2B2PE(1.00f)}}},
    {96000,
     {{FL2B2PE(1.60f), FL2B2PE(1.20f)}, {FL2B2PE(1.60f), FL2B2PE(1.20f)}}},
    {128000,
     {{FL2B2PE(1.60f), FL2B2PE(1.60f)}, {FL2B2PE(1.60f), FL2B2PE(1.40f)}}},
    {148000,
     {{FL2B2PE(1.60f), FL2B2PE(1.60f)}, {FL2B2PE(1.60f), FL2B2PE(1.60f)}}},
    {160000,
     {{FL2B2PE(1.60f), FL2B2PE(1.60f)}, {FL2B2PE(1.60f), FL2B2PE(1.60f)}}},
    {200000,
     {{FL2B2PE(1.80f), FL2B2PE(1.60f)}, {FL2B2PE(1.60f), FL2B2PE(1.60f)}}},
    {320000,
     {{FL2B2PE(3.20f), FL2B2PE(1.60f)}, {FL2B2PE(3.20f), FL2B2PE(1.60f)}}}};

static const BIT_PE_SFAC S_Bits2PeTab48000[] = {
    {16000,
     {{FL2B2PE(1.40f), FL2B2PE(0.00f)}, {FL2B2PE(0.80f), FL2B2PE(0.00f)}}},
    {24000,
     {{FL2B2PE(1.40f), FL2B2PE(1.20f)}, {FL2B2PE(1.00f), FL2B2PE(0.80f)}}},
    {32000,
     {{FL2B2PE(1.00f), FL2B2PE(1.20f)}, {FL2B2PE(0.60f), FL2B2PE(0.80f)}}},
    {48000,
     {{FL2B2PE(1.20f), FL2B2PE(1.00f)}, {FL2B2PE(0.80f), FL2B2PE(0.80f)}}},
    {64000,
     {{FL2B2PE(1.20f), FL2B2PE(1.20f)}, {FL2B2PE(1.20f), FL2B2PE(1.00f)}}},
    {96000,
     {{FL2B2PE(1.60f), FL2B2PE(1.40f)}, {FL2B2PE(1.60f), FL2B2PE(1.20f)}}},
    {128000,
     {{FL2B2PE(1.60f), FL2B2PE(1.60f)}, {FL2B2PE(1.60f), FL2B2PE(1.40f)}}},
    {148000,
     {{FL2B2PE(1.60f), FL2B2PE(1.60f)}, {FL2B2PE(1.60f), FL2B2PE(1.40f)}}},
    {160000,
     {{FL2B2PE(1.60f), FL2B2PE(1.60f)}, {FL2B2PE(1.60f), FL2B2PE(1.40f)}}},
    {200000,
     {{FL2B2PE(1.20f), FL2B2PE(1.60f)}, {FL2B2PE(1.60f), FL2B2PE(1.40f)}}},
    {320000,
     {{FL2B2PE(3.20f), FL2B2PE(1.60f)}, {FL2B2PE(3.20f), FL2B2PE(1.60f)}}}};

static const BITS2PE_CFG_TAB bits2PeConfigTab[] = {
    {16000, S_Bits2PeTab16000, sizeof(S_Bits2PeTab16000) / sizeof(BIT_PE_SFAC)},
    {22050, S_Bits2PeTab22050, sizeof(S_Bits2PeTab22050) / sizeof(BIT_PE_SFAC)},
    {24000, S_Bits2PeTab24000, sizeof(S_Bits2PeTab24000) / sizeof(BIT_PE_SFAC)},
    {32000, S_Bits2PeTab32000, sizeof(S_Bits2PeTab32000) / sizeof(BIT_PE_SFAC)},
    {44100, S_Bits2PeTab44100, sizeof(S_Bits2PeTab44100) / sizeof(BIT_PE_SFAC)},
    {48000, S_Bits2PeTab48000,
     sizeof(S_Bits2PeTab48000) / sizeof(BIT_PE_SFAC)}};

enum _avoid_hole_state { NO_AH = 0, AH_INACTIVE = 1, AH_ACTIVE = 2 };

#define Q_BITFAC \
  (24)
#define Q_AVGBITS (17)

static void FDKaacEnc_InitBits2PeFactor(
    FIXP_DBL *bits2PeFactor_m, INT *bits2PeFactor_e, const INT bitRate,
    const INT nChannels, const INT sampleRate, const INT advancedBitsToPe,
    const INT dZoneQuantEnable, const INT invQuant) {
  FIXP_DBL bit2PE_m = FL2FXCONST_DBL(1.18f / (1 << (1)));
  INT bit2PE_e = 1;

  if (advancedBitsToPe && nChannels <= (2)) {
    int i;
    const BIT_PE_SFAC *peTab = NULL;
    INT size = 0;

    for (i = 0; i < (INT)(sizeof(bits2PeConfigTab) / sizeof(BITS2PE_CFG_TAB));
         i++) {
      if (sampleRate >= bits2PeConfigTab[i].sampleRate) {
        peTab = bits2PeConfigTab[i].pPeTab;
        size = bits2PeConfigTab[i].nEntries;
      }
    }

    if ((peTab != NULL) && (size != 0)) {
      INT startB = -1;
      INT stopB = -1;
      FIXP_DBL startPF =
          FL2FXCONST_DBL(0.0f);
      FIXP_DBL stopPF = FL2FXCONST_DBL(0.0f);
      FIXP_DBL slope = FL2FXCONST_DBL(
          0.0f);
      const int qualityIdx = (invQuant == 0) ? 0 : 1;

      if (bitRate >= peTab[size - 1].bitrate) {
        startB = peTab[size - 1].bitrate;
        stopB =
            bitRate +
            1;
        startPF = peTab[size - 1].bits2PeFactor[qualityIdx][nChannels - 1];
        stopPF = peTab[size - 1].bits2PeFactor[qualityIdx][nChannels - 1];
      } else {
        for (i = 0; i < size - 1; i++) {
          if ((peTab[i].bitrate <= bitRate) &&
              (peTab[i + 1].bitrate > bitRate)) {
            startB = peTab[i].bitrate;
            stopB = peTab[i + 1].bitrate;
            startPF = peTab[i].bits2PeFactor[qualityIdx][nChannels - 1];
            stopPF = peTab[i + 1].bits2PeFactor[qualityIdx][nChannels - 1];
            break;
          }
        }
      }

      if (startB != -1) {
        FIXP_DBL bit2PE = 0;

        const FIXP_DBL maxBit2PE = FL2FXCONST_DBL(3.f / 4.f);

        slope = fDivNorm(bitRate - startB, stopB - startB);
        bit2PE = fMult(slope, stopPF - startPF) + startPF;

        bit2PE = fMin(maxBit2PE, bit2PE);

        if (bit2PE >= (FL2FXCONST_DBL(0.35f) >> 2)) {
          bit2PE_m = bit2PE;
          bit2PE_e = 2;
        }
      }
    }
  }

  if (dZoneQuantEnable) {
    if (bit2PE_m >= (FL2FXCONST_DBL(0.6f)) >> bit2PE_e) {
      bit2PE_m >>= 1;
      bit2PE_e += 1;
    }

    if ((bitRate / nChannels > 32000) && (bitRate / nChannels <= 40000)) {
      bit2PE_m += (FL2FXCONST_DBL(0.4f)) >> bit2PE_e;
    } else if (bitRate / nChannels > 20000) {
      bit2PE_m += (FL2FXCONST_DBL(0.3f)) >> bit2PE_e;
    } else if (bitRate / nChannels >= 16000) {
      bit2PE_m += (FL2FXCONST_DBL(0.3f)) >> bit2PE_e;
    } else {
      bit2PE_m += (FL2FXCONST_DBL(0.0f)) >> bit2PE_e;
    }
  }

  *bits2PeFactor_m = bit2PE_m;
  *bits2PeFactor_e = bit2PE_e;
}

FDK_INLINE INT FDKaacEnc_bits2pe2(const INT bits, const FIXP_DBL factor_m,
                                  const INT factor_e) {
  return (INT)(fMult(factor_m, (FIXP_DBL)(bits << Q_AVGBITS)) >>
               (Q_AVGBITS - factor_e));
}

static void FDKaacEnc_calcThreshExp(
    FIXP_DBL thrExp[(2)][MAX_GROUPED_SFB],
    const QC_OUT_CHANNEL *const qcOutChannel[(2)],
    const PSY_OUT_CHANNEL *const psyOutChannel[(2)], const INT nChannels) {
  INT ch, sfb, sfbGrp;
  FIXP_DBL thrExpLdData;

  for (ch = 0; ch < nChannels; ch++) {
    for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
         sfbGrp += psyOutChannel[ch]->sfbPerGroup) {
      for (sfb = 0; sfb < psyOutChannel[ch]->maxSfbPerGroup; sfb++) {
        thrExpLdData = psyOutChannel[ch]->sfbThresholdLdData[sfbGrp + sfb] >> 2;
        thrExp[ch][sfbGrp + sfb] = CalcInvLdData(thrExpLdData);
      }
    }
  }
}

static void FDKaacEnc_adaptMinSnr(
    QC_OUT_CHANNEL *const qcOutChannel[(2)],
    const PSY_OUT_CHANNEL *const psyOutChannel[(2)],
    const MINSNR_ADAPT_PARAM *const msaParam, const INT nChannels) {
  INT ch, sfb, sfbGrp, nSfb;
  FIXP_DBL avgEnLD64, dbRatio, minSnrRed;
  FIXP_DBL minSnrLimitLD64 =
      FL2FXCONST_DBL(-0.00503012648262f);
  FIXP_DBL nSfbLD64;
  FIXP_DBL accu;

  FIXP_DBL msaParam_maxRed = msaParam->maxRed;
  FIXP_DBL msaParam_startRatio = msaParam->startRatio;
  FIXP_DBL msaParam_redRatioFac =
      fMult(msaParam->redRatioFac, FL2FXCONST_DBL(0.3010299956f));
  FIXP_DBL msaParam_redOffs = msaParam->redOffs;

  for (ch = 0; ch < nChannels; ch++) {
    nSfb = 0;
    accu = FL2FXCONST_DBL(0.0f);

    DWORD_ALIGNED(psyOutChannel[ch]->sfbEnergy);

    for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
         sfbGrp += psyOutChannel[ch]->sfbPerGroup) {
      int maxSfbPerGroup = psyOutChannel[ch]->maxSfbPerGroup;
      nSfb += maxSfbPerGroup;
      for (sfb = 0; sfb < maxSfbPerGroup; sfb++) {
        accu += psyOutChannel[ch]->sfbEnergy[sfbGrp + sfb] >> 6;
      }
    }

    if ((accu == FL2FXCONST_DBL(0.0f)) || (nSfb == 0)) {
      avgEnLD64 = FL2FXCONST_DBL(-1.0f);
    } else {
      nSfbLD64 = CalcLdInt(nSfb);
      avgEnLD64 = CalcLdData(accu);
      avgEnLD64 = avgEnLD64 + FL2FXCONST_DBL(0.09375f) -
                  nSfbLD64;
    }

    int maxSfbPerGroup = psyOutChannel[ch]->maxSfbPerGroup;
    int sfbCnt = psyOutChannel[ch]->sfbCnt;
    int sfbPerGroup = psyOutChannel[ch]->sfbPerGroup;

    for (sfbGrp = 0; sfbGrp < sfbCnt; sfbGrp += sfbPerGroup) {
      FIXP_DBL *RESTRICT psfbEnergyLdData =
          &qcOutChannel[ch]->sfbEnergyLdData[sfbGrp];
      FIXP_DBL *RESTRICT psfbMinSnrLdData =
          &qcOutChannel[ch]->sfbMinSnrLdData[sfbGrp];
      for (sfb = 0; sfb < maxSfbPerGroup; sfb++) {
        FIXP_DBL sfbEnergyLdData = *psfbEnergyLdData++;
        FIXP_DBL sfbMinSnrLdData = *psfbMinSnrLdData;
        dbRatio = avgEnLD64 - sfbEnergyLdData;
        int update = (msaParam_startRatio < dbRatio) ? 1 : 0;
        minSnrRed = msaParam_redOffs + fMult(msaParam_redRatioFac,
                                             dbRatio);
        minSnrRed =
            fixMax(minSnrRed, msaParam_maxRed);
        minSnrRed = (fMult(sfbMinSnrLdData, minSnrRed)) << 6;
        minSnrRed = fixMin(minSnrLimitLD64, minSnrRed);
        *psfbMinSnrLdData++ = update ? minSnrRed : sfbMinSnrLdData;
      }
    }
  }
}

static void FDKaacEnc_initAvoidHoleFlag(
    QC_OUT_CHANNEL *const qcOutChannel[(2)],
    const PSY_OUT_CHANNEL *const psyOutChannel[(2)],
    UCHAR ahFlag[(2)][MAX_GROUPED_SFB], const struct TOOLSINFO *const toolsInfo,
    const INT nChannels, const AH_PARAM *const ahParam) {
  INT ch, sfb, sfbGrp;
  FIXP_DBL sfbEn, sfbEnm1;
  FIXP_DBL sfbEnLdData;
  FIXP_DBL avgEnLdData;

  for (ch = 0; ch < nChannels; ch++) {
    QC_OUT_CHANNEL *const qcOutChan = qcOutChannel[ch];

    if (psyOutChannel[ch]->lastWindowSequence != SHORT_WINDOW) {
      for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
           sfbGrp += psyOutChannel[ch]->sfbPerGroup)
        for (sfb = 0; sfb < psyOutChannel[ch]->maxSfbPerGroup; sfb++)
          qcOutChan->sfbSpreadEnergy[sfbGrp + sfb] >>= 1;
    } else {
      for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
           sfbGrp += psyOutChannel[ch]->sfbPerGroup)
        for (sfb = 0; sfb < psyOutChannel[ch]->maxSfbPerGroup; sfb++)
          qcOutChan->sfbSpreadEnergy[sfbGrp + sfb] = fMult(
              FL2FXCONST_DBL(0.63f), qcOutChan->sfbSpreadEnergy[sfbGrp + sfb]);
    }
  }

  if (ahParam->modifyMinSnr) {
    for (ch = 0; ch < nChannels; ch++) {
      QC_OUT_CHANNEL *const qcOutChan = qcOutChannel[ch];
      for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
           sfbGrp += psyOutChannel[ch]->sfbPerGroup) {
        for (sfb = 0; sfb < psyOutChannel[ch]->maxSfbPerGroup; sfb++) {
          FIXP_DBL sfbEnp1, avgEn;
          if (sfb > 0)
            sfbEnm1 = qcOutChan->sfbEnergy[sfbGrp + sfb - 1];
          else
            sfbEnm1 = qcOutChan->sfbEnergy[sfbGrp + sfb];

          if (sfb < psyOutChannel[ch]->maxSfbPerGroup - 1)
            sfbEnp1 = qcOutChan->sfbEnergy[sfbGrp + sfb + 1];
          else
            sfbEnp1 = qcOutChan->sfbEnergy[sfbGrp + sfb];

          avgEn = (sfbEnm1 >> 1) + (sfbEnp1 >> 1);
          avgEnLdData = CalcLdData(avgEn);
          sfbEn = qcOutChan->sfbEnergy[sfbGrp + sfb];
          sfbEnLdData = qcOutChan->sfbEnergyLdData[sfbGrp + sfb];
          if (sfbEn > avgEn) {
            FIXP_DBL tmpMinSnrLdData;
            if (psyOutChannel[ch]->lastWindowSequence == LONG_WINDOW)
              tmpMinSnrLdData = SnrLdFac + fixMax(avgEnLdData - sfbEnLdData,
                                                  SnrLdMin1 - SnrLdFac);
            else
              tmpMinSnrLdData = SnrLdFac + fixMax(avgEnLdData - sfbEnLdData,
                                                  SnrLdMin3 - SnrLdFac);

            qcOutChan->sfbMinSnrLdData[sfbGrp + sfb] = fixMin(
                qcOutChan->sfbMinSnrLdData[sfbGrp + sfb], tmpMinSnrLdData);
          }
          if (((sfbEnLdData + (FIXP_DBL)SnrLdMin4) < (FIXP_DBL)avgEnLdData) &&
              (sfbEn > FL2FXCONST_DBL(0.0))) {
            FIXP_DBL tmpMinSnrLdData = avgEnLdData - sfbEnLdData -
                                       (FIXP_DBL)SnrLdMin4 +
                                       qcOutChan->sfbMinSnrLdData[sfbGrp + sfb];
            tmpMinSnrLdData = fixMin((FIXP_DBL)SnrLdFac, tmpMinSnrLdData);
            qcOutChan->sfbMinSnrLdData[sfbGrp + sfb] =
                fixMin(tmpMinSnrLdData,
                       (FIXP_DBL)(qcOutChan->sfbMinSnrLdData[sfbGrp + sfb] +
                                  SnrLdMin2));
          }
        }
      }
    }
  }

  if (nChannels == 2) {
    QC_OUT_CHANNEL *qcOutChanM = qcOutChannel[0];
    QC_OUT_CHANNEL *qcOutChanS = qcOutChannel[1];
    const PSY_OUT_CHANNEL *const psyOutChanM = psyOutChannel[0];
    for (sfbGrp = 0; sfbGrp < psyOutChanM->sfbCnt;
         sfbGrp += psyOutChanM->sfbPerGroup) {
      for (sfb = 0; sfb < psyOutChanM->maxSfbPerGroup; sfb++) {
        if (toolsInfo->msMask[sfbGrp + sfb]) {
          FIXP_DBL maxSfbEnLd =
              fixMax(qcOutChanM->sfbEnergyLdData[sfbGrp + sfb],
                     qcOutChanS->sfbEnergyLdData[sfbGrp + sfb]);
          FIXP_DBL maxThrLd, sfbMinSnrTmpLd;

          if (((SnrLdMin5 >> 1) + (maxSfbEnLd >> 1) +
               (qcOutChanM->sfbMinSnrLdData[sfbGrp + sfb] >> 1)) <=
              FL2FXCONST_DBL(-0.5f))
            maxThrLd = FL2FXCONST_DBL(-1.0f);
          else
            maxThrLd = SnrLdMin5 + maxSfbEnLd +
                       qcOutChanM->sfbMinSnrLdData[sfbGrp + sfb];

          if (qcOutChanM->sfbEnergy[sfbGrp + sfb] > FL2FXCONST_DBL(0.0f))
            sfbMinSnrTmpLd =
                maxThrLd - qcOutChanM->sfbEnergyLdData[sfbGrp + sfb];
          else
            sfbMinSnrTmpLd = FL2FXCONST_DBL(0.0f);

          qcOutChanM->sfbMinSnrLdData[sfbGrp + sfb] =
              fixMax(qcOutChanM->sfbMinSnrLdData[sfbGrp + sfb], sfbMinSnrTmpLd);

          if (qcOutChanM->sfbMinSnrLdData[sfbGrp + sfb] <= FL2FXCONST_DBL(0.0f))
            qcOutChanM->sfbMinSnrLdData[sfbGrp + sfb] = fixMin(
                qcOutChanM->sfbMinSnrLdData[sfbGrp + sfb], (FIXP_DBL)SnrLdFac);

          if (qcOutChanS->sfbEnergy[sfbGrp + sfb] > FL2FXCONST_DBL(0.0f))
            sfbMinSnrTmpLd =
                maxThrLd - qcOutChanS->sfbEnergyLdData[sfbGrp + sfb];
          else
            sfbMinSnrTmpLd = FL2FXCONST_DBL(0.0f);

          qcOutChanS->sfbMinSnrLdData[sfbGrp + sfb] =
              fixMax(qcOutChanS->sfbMinSnrLdData[sfbGrp + sfb], sfbMinSnrTmpLd);

          if (qcOutChanS->sfbMinSnrLdData[sfbGrp + sfb] <= FL2FXCONST_DBL(0.0f))
            qcOutChanS->sfbMinSnrLdData[sfbGrp + sfb] = fixMin(
                qcOutChanS->sfbMinSnrLdData[sfbGrp + sfb], (FIXP_DBL)SnrLdFac);

          if (qcOutChanM->sfbEnergy[sfbGrp + sfb] >
              qcOutChanM->sfbSpreadEnergy[sfbGrp + sfb])
            qcOutChanS->sfbSpreadEnergy[sfbGrp + sfb] = fMult(
                qcOutChanS->sfbEnergy[sfbGrp + sfb], FL2FXCONST_DBL(0.9f));

          if (qcOutChanS->sfbEnergy[sfbGrp + sfb] >
              qcOutChanS->sfbSpreadEnergy[sfbGrp + sfb])
            qcOutChanM->sfbSpreadEnergy[sfbGrp + sfb] = fMult(
                qcOutChanM->sfbEnergy[sfbGrp + sfb], FL2FXCONST_DBL(0.9f));
        }
      }
    }
  }

  for (ch = 0; ch < nChannels; ch++) {
    QC_OUT_CHANNEL *qcOutChan = qcOutChannel[ch];
    const PSY_OUT_CHANNEL *const psyOutChan = psyOutChannel[ch];
    for (sfbGrp = 0; sfbGrp < psyOutChan->sfbCnt;
         sfbGrp += psyOutChan->sfbPerGroup) {
      for (sfb = 0; sfb < psyOutChan->maxSfbPerGroup; sfb++) {
        if ((qcOutChan->sfbSpreadEnergy[sfbGrp + sfb] >
             qcOutChan->sfbEnergy[sfbGrp + sfb]) ||
            (qcOutChan->sfbMinSnrLdData[sfbGrp + sfb] > FL2FXCONST_DBL(0.0f))) {
          ahFlag[ch][sfbGrp + sfb] = NO_AH;
        } else {
          ahFlag[ch][sfbGrp + sfb] = AH_INACTIVE;
        }
      }
    }
  }
}

static void FDKaacEnc_preparePe(PE_DATA *const peData,
                                const PSY_OUT_CHANNEL *const psyOutChannel[(2)],
                                const QC_OUT_CHANNEL *const qcOutChannel[(2)],
                                const INT nChannels, const INT peOffset) {
  INT ch;

  for (ch = 0; ch < nChannels; ch++) {
    const PSY_OUT_CHANNEL *const psyOutChan = psyOutChannel[ch];
    FDKaacEnc_prepareSfbPe(
        &peData->peChannelData[ch], psyOutChan->sfbEnergyLdData,
        psyOutChan->sfbThresholdLdData, qcOutChannel[ch]->sfbFormFactorLdData,
        psyOutChan->sfbOffsets, psyOutChan->sfbCnt, psyOutChan->sfbPerGroup,
        psyOutChan->maxSfbPerGroup);
  }
  peData->offset = peOffset;
}

static void FDKaacEnc_calcWeighting(
    const PE_DATA *const peData,
    const PSY_OUT_CHANNEL *const psyOutChannel[(2)],
    QC_OUT_CHANNEL *const qcOutChannel[(2)],
    const struct TOOLSINFO *const toolsInfo,
    ATS_ELEMENT *const adjThrStateElement, const INT nChannels,
    const INT usePatchTool) {
  int ch, noShortWindowInFrame = TRUE;
  INT exePatchM = 0;

  for (ch = 0; ch < nChannels; ch++) {
    if (psyOutChannel[ch]->lastWindowSequence == SHORT_WINDOW) {
      noShortWindowInFrame = FALSE;
    }
    FDKmemclear(qcOutChannel[ch]->sfbEnFacLd,
                MAX_GROUPED_SFB * sizeof(FIXP_DBL));
  }

  if (usePatchTool == 0) {
    return;
  }

  for (ch = 0; ch < nChannels; ch++) {
    const PSY_OUT_CHANNEL *const psyOutChan = psyOutChannel[ch];

    if (noShortWindowInFrame) {

      FIXP_DBL nrgSum14, nrgSum12, nrgSum34, nrgTotal;
      FIXP_DBL nrgFacLd_14, nrgFacLd_12, nrgFacLd_34;
      INT usePatch, exePatch;
      int sfb, sfbGrp, nLinesSum = 0;

      nrgSum14 = nrgSum12 = nrgSum34 = nrgTotal = FL2FXCONST_DBL(0.f);

      for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
           sfbGrp += psyOutChannel[ch]->sfbPerGroup) {
        for (sfb = 0; sfb < psyOutChannel[ch]->maxSfbPerGroup; sfb++) {
          FIXP_DBL nrgFac12 = CalcInvLdData(
              psyOutChan->sfbEnergyLdData[sfbGrp + sfb] >> 1);
          FIXP_DBL nrgFac14 = CalcInvLdData(
              psyOutChan->sfbEnergyLdData[sfbGrp + sfb] >> 2);

          nLinesSum += peData->peChannelData[ch]
                           .sfbNLines[sfbGrp + sfb];
          nrgTotal +=
              (psyOutChan->sfbEnergy[sfbGrp + sfb] >> 6);
          nrgSum12 += (nrgFac12 >> 6);
          nrgSum14 += (nrgFac14 >> 6);
          nrgSum34 += (fMult(nrgFac14, nrgFac12) >> 6);
        }
      }

      nrgTotal = CalcLdData(nrgTotal);

      nrgFacLd_14 =
          CalcLdData(nrgSum14) - nrgTotal;
      nrgFacLd_12 =
          CalcLdData(nrgSum12) - nrgTotal;
      nrgFacLd_34 =
          CalcLdData(nrgSum34) - nrgTotal;

      adjThrStateElement->chaosMeasureEnFac[ch] =
          fMax(FL2FXCONST_DBL(0.1875f),
               fDivNorm(nLinesSum, psyOutChan->sfbOffsets[psyOutChan->sfbCnt]));

      usePatch = (adjThrStateElement->chaosMeasureEnFac[ch] >
                  FL2FXCONST_DBL(0.78125f));
      exePatch = ((usePatch) && (adjThrStateElement->lastEnFacPatch[ch]));

      for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
           sfbGrp += psyOutChannel[ch]->sfbPerGroup) {
        for (sfb = 0; sfb < psyOutChannel[ch]->maxSfbPerGroup; sfb++) {
          INT sfbExePatch;
          if ((ch == 1) && (toolsInfo->msMask[sfbGrp + sfb])) {
            sfbExePatch = exePatchM;
          } else {
            sfbExePatch = exePatch;
          }

          if ((sfbExePatch) &&
              (psyOutChan->sfbEnergy[sfbGrp + sfb] > FL2FXCONST_DBL(0.f))) {
            if (adjThrStateElement->chaosMeasureEnFac[ch] >
                FL2FXCONST_DBL(0.8125f)) {
              qcOutChannel[ch]->sfbEnFacLd[sfbGrp + sfb] =
                  ((nrgFacLd_14 +
                    (psyOutChan->sfbEnergyLdData[sfbGrp + sfb] +
                     (psyOutChan->sfbEnergyLdData[sfbGrp + sfb] >> 1))) >>
                   1);
            } else if (adjThrStateElement->chaosMeasureEnFac[ch] >
                       FL2FXCONST_DBL(0.796875f)) {
              qcOutChannel[ch]->sfbEnFacLd[sfbGrp + sfb] =
                  ((nrgFacLd_12 + psyOutChan->sfbEnergyLdData[sfbGrp + sfb]) >>
                   1);
            } else {
              qcOutChannel[ch]->sfbEnFacLd[sfbGrp + sfb] =
                  ((nrgFacLd_34 +
                    (psyOutChan->sfbEnergyLdData[sfbGrp + sfb] >> 1)) >>
                   1);
            }
            qcOutChannel[ch]->sfbEnFacLd[sfbGrp + sfb] =
                fixMin(qcOutChannel[ch]->sfbEnFacLd[sfbGrp + sfb], (FIXP_DBL)0);
          }
        }
      }

      adjThrStateElement->lastEnFacPatch[ch] = usePatch;
      exePatchM = exePatch;
    } else {
      adjThrStateElement->chaosMeasureEnFac[ch] = FL2FXCONST_DBL(0.75f);
      adjThrStateElement->lastEnFacPatch[ch] =
          TRUE;
    }
  }
}

static void FDKaacEnc_calcPe(const PSY_OUT_CHANNEL *const psyOutChannel[(2)],
                             const QC_OUT_CHANNEL *const qcOutChannel[(2)],
                             PE_DATA *const peData, const INT nChannels) {
  INT ch;

  peData->pe = peData->offset;
  peData->constPart = 0;
  peData->nActiveLines = 0;
  for (ch = 0; ch < nChannels; ch++) {
    PE_CHANNEL_DATA *peChanData = &peData->peChannelData[ch];

    FDKaacEnc_calcSfbPe(
        peChanData, qcOutChannel[ch]->sfbWeightedEnergyLdData,
        qcOutChannel[ch]->sfbThresholdLdData, psyOutChannel[ch]->sfbCnt,
        psyOutChannel[ch]->sfbPerGroup, psyOutChannel[ch]->maxSfbPerGroup,
        psyOutChannel[ch]->isBook, psyOutChannel[ch]->isScale);

    peData->pe += peChanData->pe;
    peData->constPart += peChanData->constPart;
    peData->nActiveLines += peChanData->nActiveLines;
  }
}

void FDKaacEnc_peCalculation(PE_DATA *const peData,
                             const PSY_OUT_CHANNEL *const psyOutChannel[(2)],
                             QC_OUT_CHANNEL *const qcOutChannel[(2)],
                             const struct TOOLSINFO *const toolsInfo,
                             ATS_ELEMENT *const adjThrStateElement,
                             const INT nChannels) {
  FDKaacEnc_preparePe(peData, psyOutChannel, qcOutChannel, nChannels,
                      adjThrStateElement->peOffset);

  FDKaacEnc_calcWeighting(peData, psyOutChannel, qcOutChannel, toolsInfo,
                          adjThrStateElement, nChannels, 1);
  {
    int ch;
    for (ch = 0; ch < nChannels; ch++) {
      int sfb, sfbGrp;
      QC_OUT_CHANNEL *pQcOutCh = qcOutChannel[ch];

      for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
           sfbGrp += psyOutChannel[ch]->sfbPerGroup) {
        for (sfb = 0; sfb < psyOutChannel[ch]->maxSfbPerGroup; sfb++) {
          pQcOutCh->sfbWeightedEnergyLdData[sfb + sfbGrp] =
              pQcOutCh->sfbEnergyLdData[sfb + sfbGrp] -
              pQcOutCh->sfbEnFacLd[sfb + sfbGrp];
          pQcOutCh->sfbThresholdLdData[sfb + sfbGrp] -=
              pQcOutCh->sfbEnFacLd[sfb + sfbGrp];
        }
      }
    }
  }

  FDKaacEnc_calcPe(psyOutChannel, qcOutChannel, peData, nChannels);
}

#define CONSTPART_HEADROOM 4
static void FDKaacEnc_FDKaacEnc_calcPeNoAH(
    INT *const pe, INT *const constPart, INT *const nActiveLines,
    const PE_DATA *const peData, const UCHAR ahFlag[(2)][MAX_GROUPED_SFB],
    const PSY_OUT_CHANNEL *const psyOutChannel[(2)], const INT nChannels) {
  INT ch, sfb, sfbGrp;

  INT pe_tmp = peData->offset;
  INT constPart_tmp = 0;
  INT nActiveLines_tmp = 0;
  for (ch = 0; ch < nChannels; ch++) {
    const PE_CHANNEL_DATA *const peChanData = &peData->peChannelData[ch];
    for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
         sfbGrp += psyOutChannel[ch]->sfbPerGroup) {
      for (sfb = 0; sfb < psyOutChannel[ch]->maxSfbPerGroup; sfb++) {
        if (ahFlag[ch][sfbGrp + sfb] < AH_ACTIVE) {
          pe_tmp += peChanData->sfbPe[sfbGrp + sfb];
          constPart_tmp +=
              peChanData->sfbConstPart[sfbGrp + sfb] >> CONSTPART_HEADROOM;
          nActiveLines_tmp += peChanData->sfbNActiveLines[sfbGrp + sfb];
        }
      }
    }
  }
  *pe = pe_tmp >> PE_CONSTPART_SHIFT;
  *constPart = constPart_tmp >> (PE_CONSTPART_SHIFT - CONSTPART_HEADROOM);

  *nActiveLines = nActiveLines_tmp;
}

static const FIXP_DBL limitThrReducedLdData =
    (FIXP_DBL)0x00008000;

static void FDKaacEnc_reduceThresholdsCBR(
    QC_OUT_CHANNEL *const qcOutChannel[(2)],
    const PSY_OUT_CHANNEL *const psyOutChannel[(2)],
    UCHAR ahFlag[(2)][MAX_GROUPED_SFB],
    const FIXP_DBL thrExp[(2)][MAX_GROUPED_SFB], const INT nChannels,
    const FIXP_DBL redVal_m, const SCHAR redVal_e) {
  INT ch, sfb, sfbGrp;
  FIXP_DBL sfbEnLdData, sfbThrLdData, sfbThrReducedLdData;
  FIXP_DBL sfbThrExp;

  for (ch = 0; ch < nChannels; ch++) {
    QC_OUT_CHANNEL *qcOutChan = qcOutChannel[ch];
    for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
         sfbGrp += psyOutChannel[ch]->sfbPerGroup) {
      for (sfb = 0; sfb < psyOutChannel[ch]->maxSfbPerGroup; sfb++) {
        sfbEnLdData = qcOutChan->sfbWeightedEnergyLdData[sfbGrp + sfb];
        sfbThrLdData = qcOutChan->sfbThresholdLdData[sfbGrp + sfb];
        sfbThrExp = thrExp[ch][sfbGrp + sfb];
        if ((sfbEnLdData > sfbThrLdData) &&
            (ahFlag[ch][sfbGrp + sfb] != AH_ACTIVE)) {
          int minScale = fixMin(CountLeadingBits(sfbThrExp),
                                CountLeadingBits(redVal_m) - redVal_e) -
                         1;

          sfbThrReducedLdData =
              CalcLdData(fAbs(scaleValue(sfbThrExp, minScale) +
                              scaleValue(redVal_m, redVal_e + minScale))) -
              (FIXP_DBL)(minScale << (DFRACT_BITS - 1 - LD_DATA_SHIFT));
          sfbThrReducedLdData <<= 2;

          if ((sfbThrReducedLdData >
               (qcOutChan->sfbMinSnrLdData[sfbGrp + sfb] + sfbEnLdData)) &&
              (ahFlag[ch][sfbGrp + sfb] != NO_AH)) {
            if (qcOutChan->sfbMinSnrLdData[sfbGrp + sfb] >
                (FL2FXCONST_DBL(-1.0f) - sfbEnLdData)) {
              sfbThrReducedLdData = fixMax(
                  (qcOutChan->sfbMinSnrLdData[sfbGrp + sfb] + sfbEnLdData),
                  sfbThrLdData);
            } else
              sfbThrReducedLdData = sfbThrLdData;
            ahFlag[ch][sfbGrp + sfb] = AH_ACTIVE;
          }

          if ((sfbEnLdData + (FIXP_DBL)MAXVAL_DBL) >
              FL2FXCONST_DBL(9.6336206 / LD_DATA_SCALING)) {
            sfbThrReducedLdData = fixMax(
                sfbThrReducedLdData,
                (sfbEnLdData - FL2FXCONST_DBL(9.6336206 / LD_DATA_SCALING)));
          }

          qcOutChan->sfbThresholdLdData[sfbGrp + sfb] = sfbThrReducedLdData;
        }
      }
    }
  }
}

static FIXP_DBL FDKaacEnc_calcChaosMeasure(
    const PSY_OUT_CHANNEL *const psyOutChannel,
    const FIXP_DBL *const sfbFormFactorLdData) {
#define SCALE_FORM_FAC \
  (4)
#define SCALE_NRGS (8)
#define SCALE_NLINES (16)
#define SCALE_NRGS_SQRT4 (2)
#define SCALE_NLINES_P34 (12)

  INT sfbGrp, sfb;
  FIXP_DBL chaosMeasure;
  INT frameNLines = 0;
  FIXP_DBL frameFormFactor = FL2FXCONST_DBL(0.f);
  FIXP_DBL frameEnergy = FL2FXCONST_DBL(0.f);

  for (sfbGrp = 0; sfbGrp < psyOutChannel->sfbCnt;
       sfbGrp += psyOutChannel->sfbPerGroup) {
    for (sfb = 0; sfb < psyOutChannel->maxSfbPerGroup; sfb++) {
      if (psyOutChannel->sfbEnergyLdData[sfbGrp + sfb] >
          psyOutChannel->sfbThresholdLdData[sfbGrp + sfb]) {
        frameFormFactor += (CalcInvLdData(sfbFormFactorLdData[sfbGrp + sfb]) >>
                            SCALE_FORM_FAC);
        frameNLines += (psyOutChannel->sfbOffsets[sfbGrp + sfb + 1] -
                        psyOutChannel->sfbOffsets[sfbGrp + sfb]);
        frameEnergy += (psyOutChannel->sfbEnergy[sfbGrp + sfb] >> SCALE_NRGS);
      }
    }
  }

  if (frameNLines > 0) {
    chaosMeasure = CalcInvLdData(
        (((CalcLdData(frameFormFactor) >> 1) -
          (CalcLdData(frameEnergy) >> (2 + 1))) -
         (fMultDiv2(FL2FXCONST_DBL(0.75f),
                    CalcLdData((FIXP_DBL)frameNLines
                               << (DFRACT_BITS - 1 - SCALE_NLINES))) -
          (((FIXP_DBL)(-((-SCALE_FORM_FAC + SCALE_NRGS_SQRT4 - FORM_FAC_SHIFT +
                          SCALE_NLINES_P34)
                         << (DFRACT_BITS - 1 - LD_DATA_SHIFT)))) >>
           1)))
        << 1);
  } else {
    chaosMeasure = FL2FXCONST_DBL(1.f);
  }

  return chaosMeasure;
}

static void FDKaacEnc_reduceThresholdsVBR(
    QC_OUT_CHANNEL *const qcOutChannel[(2)],
    const PSY_OUT_CHANNEL *const psyOutChannel[(2)],
    UCHAR ahFlag[(2)][MAX_GROUPED_SFB],
    const FIXP_DBL thrExp[(2)][MAX_GROUPED_SFB], const INT nChannels,
    const FIXP_DBL vbrQualFactor, FIXP_DBL *const chaosMeasureOld) {
  INT ch, sfbGrp, sfb;
  FIXP_DBL chGroupEnergy[TRANS_FAC][2];
  FIXP_DBL chChaosMeasure[2];
  FIXP_DBL frameEnergy = FL2FXCONST_DBL(1e-10f);
  FIXP_DBL chaosMeasure = FL2FXCONST_DBL(0.f);
  FIXP_DBL sfbEnLdData, sfbThrLdData, sfbThrExp;
  FIXP_DBL sfbThrReducedLdData;
  FIXP_DBL chaosMeasureAvg;
  INT groupCnt;
  FIXP_DBL redVal[TRANS_FAC];
  QC_OUT_CHANNEL *qcOutChan = NULL;
  const PSY_OUT_CHANNEL *psyOutChan = NULL;

#define SCALE_GROUP_ENERGY (8)

#define CONST_CHAOS_MEAS_AVG_FAC_0 (FL2FXCONST_DBL(0.25f))
#define CONST_CHAOS_MEAS_AVG_FAC_1 (FL2FXCONST_DBL(1.f - 0.25f))

#define MIN_LDTHRESH (FL2FXCONST_DBL(-0.515625f))

  for (ch = 0; ch < nChannels; ch++) {
    psyOutChan = psyOutChannel[ch];

    FIXP_DBL chEnergy = FL2FXCONST_DBL(0.f);
    groupCnt = 0;

    for (sfbGrp = 0; sfbGrp < psyOutChan->sfbCnt;
         sfbGrp += psyOutChan->sfbPerGroup, groupCnt++) {
      chGroupEnergy[groupCnt][ch] = FL2FXCONST_DBL(0.f);
      for (sfb = 0; sfb < psyOutChan->maxSfbPerGroup; sfb++) {
        chGroupEnergy[groupCnt][ch] +=
            (psyOutChan->sfbEnergy[sfbGrp + sfb] >> SCALE_GROUP_ENERGY);
      }
      chEnergy += chGroupEnergy[groupCnt][ch];
    }
    frameEnergy += chEnergy;

    if (psyOutChannel[0]->lastWindowSequence == SHORT_WINDOW) {
      chChaosMeasure[ch] = FL2FXCONST_DBL(
          0.5f);
    } else {
      chChaosMeasure[ch] = FDKaacEnc_calcChaosMeasure(
          psyOutChannel[ch], qcOutChannel[ch]->sfbFormFactorLdData);
    }
    chaosMeasure += fMult(chChaosMeasure[ch], chEnergy);
  }

  if (frameEnergy > chaosMeasure) {
    INT scale = CntLeadingZeros(frameEnergy) - 1;
    FIXP_DBL num = chaosMeasure << scale;
    FIXP_DBL denum = frameEnergy << scale;
    chaosMeasure = schur_div(num, denum, 16);
  } else {
    chaosMeasure = FL2FXCONST_DBL(1.f);
  }

  chaosMeasureAvg = fMult(CONST_CHAOS_MEAS_AVG_FAC_0, chaosMeasure) +
                    fMult(CONST_CHAOS_MEAS_AVG_FAC_1,
                          *chaosMeasureOld);
  *chaosMeasureOld = chaosMeasure = (fixMin(
      chaosMeasure, chaosMeasureAvg));

  chaosMeasure = ((FL2FXCONST_DBL(0.2f) >> 2) +
                  fMult(FL2FXCONST_DBL(0.7f / (4.f * 0.3f)),
                        (chaosMeasure - FL2FXCONST_DBL(0.2f))));
  chaosMeasure =
      (fixMin((FIXP_DBL)(FL2FXCONST_DBL(1.0f) >> 2),
              fixMax((FIXP_DBL)(FL2FXCONST_DBL(0.1f) >> 2), chaosMeasure)))
      << 2;

  if (psyOutChannel[0]->lastWindowSequence == SHORT_WINDOW) {
    FDK_ASSERT(TRANS_FAC == 8);
#define WIN_TYPE_SCALE (3)

    groupCnt = 0;
    for (sfbGrp = 0; sfbGrp < psyOutChannel[0]->sfbCnt;
         sfbGrp += psyOutChannel[0]->sfbPerGroup, groupCnt++) {
      FIXP_DBL groupEnergy = FL2FXCONST_DBL(0.f);

      for (ch = 0; ch < nChannels; ch++) {
        groupEnergy +=
            chGroupEnergy[groupCnt]
                         [ch];
      }

      FDK_ASSERT(psyOutChannel[0]->groupLen[groupCnt] <= INV_INT_TAB_SIZE);
      groupEnergy = fMult(
          groupEnergy,
          invInt[psyOutChannel[0]->groupLen[groupCnt]]);
      groupEnergy = fixMin(groupEnergy,
                           frameEnergy >> WIN_TYPE_SCALE);

      groupEnergy >>=
          2;

      redVal[groupCnt] =
          fMult(fMult(vbrQualFactor, chaosMeasure),
                CalcInvLdData(CalcLdData(groupEnergy) >> 2))
          << (int)((2 + (2 * WIN_TYPE_SCALE) + SCALE_GROUP_ENERGY) >> 2);
    }
  } else {

    redVal[0] = fMult(fMult(vbrQualFactor, chaosMeasure),
                      CalcInvLdData(CalcLdData(frameEnergy) >> 2))
                << (int)(SCALE_GROUP_ENERGY >> 2);
  }

  for (ch = 0; ch < nChannels; ch++) {
    qcOutChan = qcOutChannel[ch];
    psyOutChan = psyOutChannel[ch];

    for (sfbGrp = 0; sfbGrp < psyOutChan->sfbCnt;
         sfbGrp += psyOutChan->sfbPerGroup) {
      for (sfb = 0; sfb < psyOutChan->maxSfbPerGroup; sfb++) {
        sfbEnLdData = (qcOutChan->sfbWeightedEnergyLdData[sfbGrp + sfb]);
        sfbThrLdData = (qcOutChan->sfbThresholdLdData[sfbGrp + sfb]);
        sfbThrExp = thrExp[ch][sfbGrp + sfb];

        if ((sfbThrLdData >= MIN_LDTHRESH) && (sfbEnLdData > sfbThrLdData) &&
            (ahFlag[ch][sfbGrp + sfb] != AH_ACTIVE)) {
          if (psyOutChannel[ch]->lastWindowSequence == SHORT_WINDOW) {
            const int groupNumber = (int)sfb / psyOutChan->sfbPerGroup;

            FDK_ASSERT(INV_SQRT4_TAB_SIZE > psyOutChan->groupLen[groupNumber]);

            sfbThrExp =
                fMult(sfbThrExp,
                      fMult(FL2FXCONST_DBL(2.82f / 4.f),
                            invSqrt4[psyOutChan->groupLen[groupNumber]]))
                << 2;

            if (sfbThrExp <= (limitThrReducedLdData - redVal[groupNumber])) {
              sfbThrReducedLdData = FL2FXCONST_DBL(-1.0f);
            } else {
              if ((FIXP_DBL)redVal[groupNumber] >=
                  FL2FXCONST_DBL(1.0f) - sfbThrExp)
                sfbThrReducedLdData = FL2FXCONST_DBL(0.0f);
              else {
                sfbThrReducedLdData =
                    CalcLdData(sfbThrExp + redVal[groupNumber]);
                sfbThrReducedLdData <<= 2;
              }
            }
            sfbThrReducedLdData +=
                (CalcLdInt(psyOutChan->groupLen[groupNumber]) -
                 ((FIXP_DBL)6 << (DFRACT_BITS - 1 - LD_DATA_SHIFT)));
          }

          else {
            if ((FIXP_DBL)redVal[0] >= FL2FXCONST_DBL(1.0f) - sfbThrExp) {
              sfbThrReducedLdData = FL2FXCONST_DBL(0.0f);
            } else {
              sfbThrReducedLdData = CalcLdData(sfbThrExp + redVal[0]);
              sfbThrReducedLdData <<= 2;
            }
          }

          if (((sfbThrReducedLdData - sfbEnLdData) >
               qcOutChan->sfbMinSnrLdData[sfbGrp + sfb]) &&
              (ahFlag[ch][sfbGrp + sfb] != NO_AH)) {
            if (qcOutChan->sfbMinSnrLdData[sfbGrp + sfb] >
                (FL2FXCONST_DBL(-1.0f) - sfbEnLdData)) {
              sfbThrReducedLdData = fixMax(
                  (qcOutChan->sfbMinSnrLdData[sfbGrp + sfb] + sfbEnLdData),
                  sfbThrLdData);
            } else
              sfbThrReducedLdData = sfbThrLdData;
            ahFlag[ch][sfbGrp + sfb] = AH_ACTIVE;
          }

          if (sfbThrReducedLdData < FL2FXCONST_DBL(-0.5f))
            sfbThrReducedLdData = FL2FXCONST_DBL(-1.f);

          sfbThrReducedLdData = fixMax(MIN_LDTHRESH, sfbThrReducedLdData);

          qcOutChan->sfbThresholdLdData[sfbGrp + sfb] = sfbThrReducedLdData;
        }
      }
    }
  }
}

static void FDKaacEnc_correctThresh(
    const CHANNEL_MAPPING *const cm, QC_OUT_ELEMENT *const qcElement[((8))],
    const PSY_OUT_ELEMENT *const psyOutElement[((8))],
    UCHAR ahFlag[((8))][(2)][MAX_GROUPED_SFB],
    const FIXP_DBL thrExp[((8))][(2)][MAX_GROUPED_SFB], const FIXP_DBL redVal_m,
    const SCHAR redVal_e, const INT deltaPe, const INT processElements,
    const INT elementOffset) {
  INT ch, sfb, sfbGrp;
  QC_OUT_CHANNEL *qcOutChan;
  PSY_OUT_CHANNEL *psyOutChan;
  PE_CHANNEL_DATA *peChanData;
  FIXP_DBL thrFactorLdData;
  FIXP_DBL sfbEnLdData, sfbThrLdData, sfbThrReducedLdData;
  FIXP_DBL *sfbPeFactorsLdData[((8))][(2)];
  FIXP_DBL(*sfbNActiveLinesLdData)[(2)][MAX_GROUPED_SFB];

  INT normFactorInt;
  FIXP_DBL normFactorLdData;

  INT nElements = elementOffset + processElements;
  INT elementId;

  for (elementId = elementOffset; elementId < nElements; elementId++) {
    for (ch = 0; ch < cm->elInfo[elementId].nChannelsInEl; ch++) {
      sfbPeFactorsLdData[elementId][ch] =
          (FIXP_DBL *)((void *)(
              qcElement[elementId]->qcOutChannel[ch]->quantSpec));
    }
  }
  sfbNActiveLinesLdData = (FIXP_DBL(*)[(2)][MAX_GROUPED_SFB])(
      (void *)(qcElement[0]->dynMem_SfbNActiveLinesLdData));

  normFactorInt = 0;

  for (elementId = elementOffset; elementId < nElements; elementId++) {
    if (cm->elInfo[elementId].elType != ID_DSE) {
      for (ch = 0; ch < cm->elInfo[elementId].nChannelsInEl; ch++) {
        psyOutChan = psyOutElement[elementId]->psyOutChannel[ch];
        peChanData = &qcElement[elementId]->peData.peChannelData[ch];

        for (sfbGrp = 0; sfbGrp < psyOutChan->sfbCnt;
             sfbGrp += psyOutChan->sfbPerGroup) {
          for (sfb = 0; sfb < psyOutChan->maxSfbPerGroup; sfb++) {
            if (peChanData->sfbNActiveLines[sfbGrp + sfb] == 0) {
              sfbNActiveLinesLdData[elementId][ch][sfbGrp + sfb] =
                  FL2FXCONST_DBL(-1.0f);
            } else {
              sfbNActiveLinesLdData[elementId][ch][sfbGrp + sfb] =
                  CalcLdInt(peChanData->sfbNActiveLines[sfbGrp + sfb]);
            }
            if (((ahFlag[elementId][ch][sfbGrp + sfb] < AH_ACTIVE) ||
                 (deltaPe > 0)) &&
                peChanData->sfbNActiveLines[sfbGrp + sfb] != 0) {
              if (thrExp[elementId][ch][sfbGrp + sfb] > -redVal_m) {

                int minScale =
                    fixMin(
                        CountLeadingBits(thrExp[elementId][ch][sfbGrp + sfb]),
                        CountLeadingBits(redVal_m) - redVal_e) -
                    1;

                FIXP_DBL sumLd =
                    CalcLdData(scaleValue(thrExp[elementId][ch][sfbGrp + sfb],
                                          minScale) +
                               scaleValue(redVal_m, redVal_e + minScale)) -
                    (FIXP_DBL)(minScale << (DFRACT_BITS - 1 - LD_DATA_SHIFT));

                if (sumLd < FL2FXCONST_DBL(0.f)) {
                  sfbPeFactorsLdData[elementId][ch][sfbGrp + sfb] =
                      sfbNActiveLinesLdData[elementId][ch][sfbGrp + sfb] -
                      sumLd;
                } else {
                  if (sfbNActiveLinesLdData[elementId][ch][sfbGrp + sfb] >
                      (FL2FXCONST_DBL(-1.f) + sumLd)) {
                    sfbPeFactorsLdData[elementId][ch][sfbGrp + sfb] =
                        sfbNActiveLinesLdData[elementId][ch][sfbGrp + sfb] -
                        sumLd;
                  } else {
                    sfbPeFactorsLdData[elementId][ch][sfbGrp + sfb] =
                        sfbNActiveLinesLdData[elementId][ch][sfbGrp + sfb];
                  }
                }

                normFactorInt += (INT)CalcInvLdData(
                    sfbPeFactorsLdData[elementId][ch][sfbGrp + sfb]);
              } else
                sfbPeFactorsLdData[elementId][ch][sfbGrp + sfb] =
                    FL2FXCONST_DBL(1.0f);
            } else
              sfbPeFactorsLdData[elementId][ch][sfbGrp + sfb] =
                  FL2FXCONST_DBL(-1.0f);
          }
        }
      }
    }
  }

  normFactorLdData =
      CalcLdData((FIXP_DBL)((deltaPe < 0) ? (-deltaPe) : (deltaPe))) -
      CalcLdData((FIXP_DBL)normFactorInt);

  for (elementId = elementOffset; elementId < nElements; elementId++) {
    if (cm->elInfo[elementId].elType != ID_DSE) {
      for (ch = 0; ch < cm->elInfo[elementId].nChannelsInEl; ch++) {
        qcOutChan = qcElement[elementId]->qcOutChannel[ch];
        psyOutChan = psyOutElement[elementId]->psyOutChannel[ch];
        peChanData = &qcElement[elementId]->peData.peChannelData[ch];

        for (sfbGrp = 0; sfbGrp < psyOutChan->sfbCnt;
             sfbGrp += psyOutChan->sfbPerGroup) {
          for (sfb = 0; sfb < psyOutChan->maxSfbPerGroup; sfb++) {
            if (peChanData->sfbNActiveLines[sfbGrp + sfb] > 0) {
              if ((sfbPeFactorsLdData[elementId][ch][sfbGrp + sfb] ==
                   FL2FXCONST_DBL(-1.0f)) ||
                  (deltaPe == 0)) {
                thrFactorLdData = FL2FXCONST_DBL(0.f);
              } else {
                FIXP_DBL tmp = CalcInvLdData(
                    sfbPeFactorsLdData[elementId][ch][sfbGrp + sfb] +
                    normFactorLdData -
                    sfbNActiveLinesLdData[elementId][ch][sfbGrp + sfb] -
                    FL2FXCONST_DBL((float)LD_DATA_SHIFT / LD_DATA_SCALING));

                tmp = (deltaPe < 0) ? tmp : (-tmp);
                thrFactorLdData =
                    fMin(tmp, FL2FXCONST_DBL(20.f / LD_DATA_SCALING));
              }

              sfbThrLdData = qcOutChan->sfbThresholdLdData[sfbGrp + sfb];
              sfbEnLdData = qcOutChan->sfbWeightedEnergyLdData[sfbGrp + sfb];

              if (thrFactorLdData < FL2FXCONST_DBL(0.f)) {
                if (sfbThrLdData > (FL2FXCONST_DBL(-1.f) - thrFactorLdData)) {
                  sfbThrReducedLdData = sfbThrLdData + thrFactorLdData;
                } else {
                  sfbThrReducedLdData = FL2FXCONST_DBL(-1.f);
                }
              } else {
                sfbThrReducedLdData = sfbThrLdData + thrFactorLdData;
              }

              if ((sfbThrReducedLdData - sfbEnLdData >
                   qcOutChan->sfbMinSnrLdData[sfbGrp + sfb]) &&
                  (ahFlag[elementId][ch][sfbGrp + sfb] == AH_INACTIVE)) {
                if (sfbEnLdData >
                    (sfbThrLdData - qcOutChan->sfbMinSnrLdData[sfbGrp + sfb])) {
                  sfbThrReducedLdData =
                      qcOutChan->sfbMinSnrLdData[sfbGrp + sfb] + sfbEnLdData;
                } else {
                  sfbThrReducedLdData = sfbThrLdData;
                }
                ahFlag[elementId][ch][sfbGrp + sfb] = AH_ACTIVE;
              }

              qcOutChan->sfbThresholdLdData[sfbGrp + sfb] = sfbThrReducedLdData;
            }
          }
        }
      }
    }
  }
}

static void FDKaacEnc_reduceMinSnr(
    const CHANNEL_MAPPING *const cm, QC_OUT_ELEMENT *const qcElement[((8))],
    const PSY_OUT_ELEMENT *const psyOutElement[((8))],
    const UCHAR ahFlag[((8))][(2)][MAX_GROUPED_SFB], const INT desiredPe,
    INT *const redPeGlobal, const INT processElements, const INT elementOffset)

{
  INT ch, elementId, globalMaxSfb = 0;
  const INT nElements = elementOffset + processElements;
  INT newGlobalPe = *redPeGlobal;

  if (newGlobalPe <= desiredPe) {
    goto bail;
  }

  for (elementId = elementOffset; elementId < nElements; elementId++) {
    if (cm->elInfo[elementId].elType != ID_DSE) {
      for (ch = 0; ch < cm->elInfo[elementId].nChannelsInEl; ch++) {
        globalMaxSfb =
            fMax(globalMaxSfb,
                 psyOutElement[elementId]->psyOutChannel[ch]->maxSfbPerGroup);
      }
    }
  }

  while ((newGlobalPe > desiredPe) && (--globalMaxSfb >= 0)) {
    for (elementId = elementOffset; elementId < nElements; elementId++) {
      if (cm->elInfo[elementId].elType != ID_DSE) {
        PE_DATA *peData = &qcElement[elementId]->peData;

        for (ch = 0; ch < cm->elInfo[elementId].nChannelsInEl; ch++) {
          QC_OUT_CHANNEL *qcOutChan = qcElement[elementId]->qcOutChannel[ch];
          PSY_OUT_CHANNEL *psyOutChan =
              psyOutElement[elementId]->psyOutChannel[ch];

          if (globalMaxSfb < psyOutChan->maxSfbPerGroup) {
            INT sfb, deltaPe = 0;

            for (sfb = globalMaxSfb; sfb < psyOutChan->sfbCnt;
                 sfb += psyOutChan->sfbPerGroup) {
              if (ahFlag[elementId][ch][sfb] != NO_AH &&
                  qcOutChan->sfbMinSnrLdData[sfb] < SnrLdFac &&
                  (qcOutChan->sfbWeightedEnergyLdData[sfb] >
                   qcOutChan->sfbThresholdLdData[sfb] - SnrLdFac)) {
                qcOutChan->sfbMinSnrLdData[sfb] = SnrLdFac;
                qcOutChan->sfbThresholdLdData[sfb] =
                    qcOutChan->sfbWeightedEnergyLdData[sfb] + SnrLdFac;

                deltaPe -= peData->peChannelData[ch].sfbPe[sfb];

                peData->peChannelData[ch].sfbPe[sfb] =
                    (3 * peData->peChannelData[ch].sfbNLines[sfb])
                    << (PE_CONSTPART_SHIFT - 1);
                deltaPe += peData->peChannelData[ch].sfbPe[sfb];
              }
            }

            deltaPe >>= PE_CONSTPART_SHIFT;
            peData->pe += deltaPe;
            peData->peChannelData[ch].pe += deltaPe;
            newGlobalPe += deltaPe;
          }
          if (newGlobalPe <= desiredPe) {
            goto bail;
          }
        }
      }
    }
  }

bail:
  *redPeGlobal = newGlobalPe;
}

static void FDKaacEnc_allowMoreHoles(
    const CHANNEL_MAPPING *const cm, QC_OUT_ELEMENT *const qcElement[((8))],
    const PSY_OUT_ELEMENT *const psyOutElement[((8))],
    const ATS_ELEMENT *const AdjThrStateElement[((8))],
    UCHAR ahFlag[((8))][(2)][MAX_GROUPED_SFB], const INT desiredPe,
    const INT currentPe, const int processElements, const int elementOffset) {
  INT elementId;
  INT nElements = elementOffset + processElements;
  INT actPe = currentPe;

  if (actPe <= desiredPe) {
    return;
  }

  for (elementId = elementOffset; elementId < nElements; elementId++) {
    if (cm->elInfo[elementId].elType != ID_DSE) {
      INT ch, sfb, sfbGrp;

      PE_DATA *peData = &qcElement[elementId]->peData;
      const INT nChannels = cm->elInfo[elementId].nChannelsInEl;

      QC_OUT_CHANNEL *qcOutChannel[(2)] = {NULL};
      PSY_OUT_CHANNEL *psyOutChannel[(2)] = {NULL};

      for (ch = 0; ch < nChannels; ch++) {
        qcOutChannel[ch] = qcElement[elementId]->qcOutChannel[ch];
        psyOutChannel[ch] = psyOutElement[elementId]->psyOutChannel[ch];

        for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
             sfbGrp += psyOutChannel[ch]->sfbPerGroup) {
          for (sfb = psyOutChannel[ch]->maxSfbPerGroup;
               sfb < psyOutChannel[ch]->sfbPerGroup; sfb++) {
            peData->peChannelData[ch].sfbPe[sfbGrp + sfb] = 0;
          }
        }
      }

      if (nChannels == 2 && psyOutChannel[0]->lastWindowSequence ==
                                psyOutChannel[1]->lastWindowSequence) {
        for (sfb = psyOutChannel[0]->maxSfbPerGroup - 1; sfb >= 0; sfb--) {
          for (sfbGrp = 0; sfbGrp < psyOutChannel[0]->sfbCnt;
               sfbGrp += psyOutChannel[0]->sfbPerGroup) {
            if (psyOutElement[elementId]->toolsInfo.msMask[sfbGrp + sfb]) {
              FIXP_DBL EnergyLd_L =
                  qcOutChannel[0]->sfbWeightedEnergyLdData[sfbGrp + sfb];
              FIXP_DBL EnergyLd_R =
                  qcOutChannel[1]->sfbWeightedEnergyLdData[sfbGrp + sfb];

              if ((ahFlag[elementId][1][sfbGrp + sfb] != NO_AH) &&
                  (((FL2FXCONST_DBL(-0.02065512648f) >> 1) +
                    (qcOutChannel[0]->sfbMinSnrLdData[sfbGrp + sfb] >> 1)) >
                   ((EnergyLd_R >> 1) - (EnergyLd_L >> 1)))) {
                ahFlag[elementId][1][sfbGrp + sfb] = NO_AH;
                qcOutChannel[1]->sfbThresholdLdData[sfbGrp + sfb] =
                    FL2FXCONST_DBL(0.015625f) + EnergyLd_R;
                actPe -= peData->peChannelData[1].sfbPe[sfbGrp + sfb] >>
                         PE_CONSTPART_SHIFT;
              }
              else if ((ahFlag[elementId][0][sfbGrp + sfb] != NO_AH) &&
                       (((FL2FXCONST_DBL(-0.02065512648f) >> 1) +
                         (qcOutChannel[1]->sfbMinSnrLdData[sfbGrp + sfb] >>
                          1)) > ((EnergyLd_L >> 1) - (EnergyLd_R >> 1)))) {
                ahFlag[elementId][0][sfbGrp + sfb] = NO_AH;
                qcOutChannel[0]->sfbThresholdLdData[sfbGrp + sfb] =
                    FL2FXCONST_DBL(0.015625f) + EnergyLd_L;
                actPe -= peData->peChannelData[0].sfbPe[sfbGrp + sfb] >>
                         PE_CONSTPART_SHIFT;
              }
            }
          }
          if (actPe <= desiredPe) {
            return;
          }
        }
      }
    }
  }

  if (actPe > desiredPe) {
    INT ch, sfb, sfbGrp;
    INT minSfb, maxSfb;
    INT enIdx, ahCnt, done;
    INT startSfb[(8)];
    INT sfbCnt[(8)];
    INT sfbPerGroup[(8)];
    INT maxSfbPerGroup[(8)];
    FIXP_DBL avgEn;
    FIXP_DBL minEnLD64;
    FIXP_DBL avgEnLD64;
    FIXP_DBL enLD64[NUM_NRG_LEVS];
    INT avgEn_e;

    maxSfb = 0;
    for (elementId = elementOffset; elementId < nElements; elementId++) {
      if (cm->elInfo[elementId].elType != ID_DSE) {
        for (ch = 0; ch < cm->elInfo[elementId].nChannelsInEl; ch++) {
          for (sfbGrp = 0;
               sfbGrp < psyOutElement[elementId]->psyOutChannel[ch]->sfbCnt;
               sfbGrp +=
               psyOutElement[elementId]->psyOutChannel[ch]->sfbPerGroup) {
            maxSfb +=
                psyOutElement[elementId]->psyOutChannel[ch]->maxSfbPerGroup;
          }
        }
      }
    }
    avgEn_e =
        (DFRACT_BITS - fixnormz_D((LONG)fMax(0, maxSfb - 1)));

    ahCnt = 0;
    maxSfb = 0;
    minSfb = MAX_SFB;
    avgEn = FL2FXCONST_DBL(0.0f);
    minEnLD64 = FL2FXCONST_DBL(0.0f);

    for (elementId = elementOffset; elementId < nElements; elementId++) {
      if (cm->elInfo[elementId].elType != ID_DSE) {
        for (ch = 0; ch < cm->elInfo[elementId].nChannelsInEl; ch++) {
          const INT chIdx = cm->elInfo[elementId].ChannelIndex[ch];
          QC_OUT_CHANNEL *qcOutChannel = qcElement[elementId]->qcOutChannel[ch];
          PSY_OUT_CHANNEL *psyOutChannel =
              psyOutElement[elementId]->psyOutChannel[ch];

          maxSfbPerGroup[chIdx] = psyOutChannel->maxSfbPerGroup;
          sfbCnt[chIdx] = psyOutChannel->sfbCnt;
          sfbPerGroup[chIdx] = psyOutChannel->sfbPerGroup;

          maxSfb = fMax(maxSfb, psyOutChannel->maxSfbPerGroup);

          if (psyOutChannel->lastWindowSequence != SHORT_WINDOW) {
            startSfb[chIdx] = AdjThrStateElement[elementId]->ahParam.startSfbL;
          } else {
            startSfb[chIdx] = AdjThrStateElement[elementId]->ahParam.startSfbS;
          }

          minSfb = fMin(minSfb, startSfb[chIdx]);

          sfbGrp = 0;
          sfb = startSfb[chIdx];

          do {
            for (; sfb < psyOutChannel->maxSfbPerGroup; sfb++) {
              if ((ahFlag[elementId][ch][sfbGrp + sfb] != NO_AH) &&
                  (qcOutChannel->sfbWeightedEnergyLdData[sfbGrp + sfb] >
                   qcOutChannel->sfbThresholdLdData[sfbGrp + sfb])) {
                minEnLD64 = fixMin(minEnLD64,
                                   qcOutChannel->sfbEnergyLdData[sfbGrp + sfb]);
                avgEn += qcOutChannel->sfbEnergy[sfbGrp + sfb] >> avgEn_e;
                ahCnt++;
              }
            }

            sfbGrp += psyOutChannel->sfbPerGroup;
            sfb = startSfb[chIdx];

          } while (sfbGrp < psyOutChannel->sfbCnt);
        }
      }
    }

    if ((avgEn == FL2FXCONST_DBL(0.0f)) || (ahCnt == 0)) {
      avgEnLD64 = FL2FXCONST_DBL(0.0f);
    } else {
      avgEnLD64 = CalcLdData(avgEn) +
                  (FIXP_DBL)(avgEn_e << (DFRACT_BITS - 1 - LD_DATA_SHIFT)) -
                  CalcLdInt(ahCnt);
    }

    enLD64[0] =
        minEnLD64 + fMult((avgEnLD64 - minEnLD64), FL2FXCONST_DBL(0.06666667f));
    enLD64[1] =
        minEnLD64 + fMult((avgEnLD64 - minEnLD64), FL2FXCONST_DBL(0.20000000f));
    enLD64[2] =
        minEnLD64 + fMult((avgEnLD64 - minEnLD64), FL2FXCONST_DBL(0.33333334f));
    enLD64[3] =
        minEnLD64 + fMult((avgEnLD64 - minEnLD64), FL2FXCONST_DBL(0.46666667f));
    enLD64[4] =
        minEnLD64 + fMult((avgEnLD64 - minEnLD64), FL2FXCONST_DBL(0.60000002f));
    enLD64[5] =
        minEnLD64 + fMult((avgEnLD64 - minEnLD64), FL2FXCONST_DBL(0.73333335f));
    enLD64[6] =
        minEnLD64 + fMult((avgEnLD64 - minEnLD64), FL2FXCONST_DBL(0.86666667f));
    enLD64[7] = minEnLD64 + (avgEnLD64 - minEnLD64);

    done = 0;
    enIdx = 0;
    sfb = maxSfb - 1;

    while (!done) {
      for (elementId = elementOffset; elementId < nElements; elementId++) {
        if (cm->elInfo[elementId].elType != ID_DSE) {
          PE_DATA *peData = &qcElement[elementId]->peData;
          for (ch = 0; ch < cm->elInfo[elementId].nChannelsInEl; ch++) {
            const INT chIdx = cm->elInfo[elementId].ChannelIndex[ch];
            QC_OUT_CHANNEL *qcOutChannel =
                qcElement[elementId]->qcOutChannel[ch];
            if (sfb >= startSfb[chIdx] && sfb < maxSfbPerGroup[chIdx]) {
              for (sfbGrp = 0; sfbGrp < sfbCnt[chIdx];
                   sfbGrp += sfbPerGroup[chIdx]) {
                if (ahFlag[elementId][ch][sfbGrp + sfb] != NO_AH &&
                    qcOutChannel->sfbEnergyLdData[sfbGrp + sfb] <
                        enLD64[enIdx]) {
                  ahFlag[elementId][ch][sfbGrp + sfb] = NO_AH;
                  qcOutChannel->sfbThresholdLdData[sfbGrp + sfb] =
                      FL2FXCONST_DBL(0.015625f) +
                      qcOutChannel->sfbWeightedEnergyLdData[sfbGrp + sfb];
                  actPe -= peData->peChannelData[ch].sfbPe[sfbGrp + sfb] >>
                           PE_CONSTPART_SHIFT;
                }
                if (actPe <= desiredPe) {
                  return;
                }
              }
            }
          }
        }
      }

      sfb--;
      if (sfb < minSfb) {
        sfb = maxSfb;
        enIdx++;
        if (enIdx >= NUM_NRG_LEVS) {
          done = 1;
        }
      }
    }
  }
}

static void FDKaacEnc_resetAHFlags(
    UCHAR ahFlag[(2)][MAX_GROUPED_SFB], const INT nChannels,
    const PSY_OUT_CHANNEL *const psyOutChannel[(2)]) {
  int ch, sfb, sfbGrp;

  for (ch = 0; ch < nChannels; ch++) {
    for (sfbGrp = 0; sfbGrp < psyOutChannel[ch]->sfbCnt;
         sfbGrp += psyOutChannel[ch]->sfbPerGroup) {
      for (sfb = 0; sfb < psyOutChannel[ch]->maxSfbPerGroup; sfb++) {
        if (ahFlag[ch][sfbGrp + sfb] == AH_ACTIVE) {
          ahFlag[ch][sfbGrp + sfb] = AH_INACTIVE;
        }
      }
    }
  }
}

static FIXP_DBL CalcRedValPower(FIXP_DBL num, FIXP_DBL denum, INT *scaling) {
  FIXP_DBL value = FL2FXCONST_DBL(0.f);

  if (num >= FL2FXCONST_DBL(0.f)) {
    value = fDivNorm(num, denum, scaling);
  } else {
    value = -fDivNorm(-num, denum, scaling);
  }
  value = f2Pow(value, *scaling, scaling);

  return value;
}

static void FDKaacEnc_adaptThresholdsToPe(
    const CHANNEL_MAPPING *const cm,
    ATS_ELEMENT *const AdjThrStateElement[((8))],
    QC_OUT_ELEMENT *const qcElement[((8))],
    const PSY_OUT_ELEMENT *const psyOutElement[((8))], const INT desiredPe,
    const INT maxIter2ndGuess, const INT processElements,
    const INT elementOffset) {
  FIXP_DBL reductionValue_m;
  SCHAR reductionValue_e;
  UCHAR(*pAhFlag)[(2)][MAX_GROUPED_SFB];
  FIXP_DBL(*pThrExp)[(2)][MAX_GROUPED_SFB];
  int iter;

  INT constPartGlobal, noRedPeGlobal, nActiveLinesGlobal, redPeGlobal;
  constPartGlobal = noRedPeGlobal = nActiveLinesGlobal = redPeGlobal = 0;

  int elementId;

  int nElements = elementOffset + processElements;
  if (nElements > cm->nElements) {
    nElements = cm->nElements;
  }

  pAhFlag = (UCHAR(*)[(2)][MAX_GROUPED_SFB])(
      (void *)(qcElement[0]->dynMem_Ah_Flag));
  pThrExp = (FIXP_DBL(*)[(2)][MAX_GROUPED_SFB])(
      (void *)(qcElement[0]->dynMem_Thr_Exp));

  for (elementId = elementOffset; elementId < nElements; elementId++) {
    if (cm->elInfo[elementId].elType != ID_DSE) {
      INT nChannels = cm->elInfo[elementId].nChannelsInEl;
      PE_DATA *peData = &qcElement[elementId]->peData;

      FDKaacEnc_calcThreshExp(
          pThrExp[elementId], qcElement[elementId]->qcOutChannel,
          psyOutElement[elementId]->psyOutChannel, nChannels);

      FDKaacEnc_adaptMinSnr(qcElement[elementId]->qcOutChannel,
                            psyOutElement[elementId]->psyOutChannel,
                            &AdjThrStateElement[elementId]->minSnrAdaptParam,
                            nChannels);

      FDKaacEnc_initAvoidHoleFlag(
          qcElement[elementId]->qcOutChannel,
          psyOutElement[elementId]->psyOutChannel, pAhFlag[elementId],
          &psyOutElement[elementId]->toolsInfo, nChannels,
          &AdjThrStateElement[elementId]->ahParam);

      constPartGlobal += peData->constPart;
      noRedPeGlobal += peData->pe;
      nActiveLinesGlobal += fixMax((INT)peData->nActiveLines, 1);
    }
  }

  int redVal_e, avgThrExp_e, result_e;
  FIXP_DBL redVal_m, avgThrExp_m;

  redVal_m = CalcRedValPower(constPartGlobal - desiredPe,
                             4 * nActiveLinesGlobal, &redVal_e);
  avgThrExp_m = CalcRedValPower(constPartGlobal - noRedPeGlobal,
                                4 * nActiveLinesGlobal, &avgThrExp_e);
  result_e = fMax(redVal_e, avgThrExp_e) + 1;

  reductionValue_m = fMax(FL2FXCONST_DBL(0.f),
                          scaleValue(redVal_m, redVal_e - result_e) -
                              scaleValue(avgThrExp_m, avgThrExp_e - result_e));
  reductionValue_e = result_e;

  for (elementId = elementOffset; elementId < nElements; elementId++) {
    if (cm->elInfo[elementId].elType != ID_DSE) {
      INT nChannels = cm->elInfo[elementId].nChannelsInEl;
      PE_DATA *peData = &qcElement[elementId]->peData;

      FDKaacEnc_reduceThresholdsCBR(
          qcElement[elementId]->qcOutChannel,
          psyOutElement[elementId]->psyOutChannel, pAhFlag[elementId],
          pThrExp[elementId], nChannels, reductionValue_m, reductionValue_e);

      FDKaacEnc_calcPe(psyOutElement[elementId]->psyOutChannel,
                       qcElement[elementId]->qcOutChannel, peData, nChannels);

      redPeGlobal += peData->pe;
    }
  }

  iter = 0;
  while ((fixp_abs(redPeGlobal - desiredPe) >
          fMultI(FL2FXCONST_DBL(0.05f), desiredPe)) &&
         (iter < maxIter2ndGuess)) {
    INT desiredPeNoAHGlobal;
    INT redPeNoAHGlobal = 0;
    INT constPartNoAHGlobal = 0;
    INT nActiveLinesNoAHGlobal = 0;

    for (elementId = elementOffset; elementId < nElements; elementId++) {
      if (cm->elInfo[elementId].elType != ID_DSE) {
        INT redPeNoAH, constPartNoAH, nActiveLinesNoAH;
        INT nChannels = cm->elInfo[elementId].nChannelsInEl;
        PE_DATA *peData = &qcElement[elementId]->peData;

        FDKaacEnc_FDKaacEnc_calcPeNoAH(
            &redPeNoAH, &constPartNoAH, &nActiveLinesNoAH, peData,
            pAhFlag[elementId], psyOutElement[elementId]->psyOutChannel,
            nChannels);

        redPeNoAHGlobal += redPeNoAH;
        constPartNoAHGlobal += constPartNoAH;
        nActiveLinesNoAHGlobal += nActiveLinesNoAH;
      }
    }

    if (desiredPe < redPeGlobal) {
      desiredPeNoAHGlobal = desiredPe - (redPeGlobal - redPeNoAHGlobal);
      desiredPeNoAHGlobal = fMax(0, desiredPeNoAHGlobal);

      if (nActiveLinesNoAHGlobal > 0) {
        redVal_m = CalcRedValPower(constPartNoAHGlobal - desiredPeNoAHGlobal,
                                   4 * nActiveLinesNoAHGlobal, &redVal_e);
        avgThrExp_m = CalcRedValPower(constPartNoAHGlobal - redPeNoAHGlobal,
                                      4 * nActiveLinesNoAHGlobal, &avgThrExp_e);
        result_e = fMax(reductionValue_e, fMax(redVal_e, avgThrExp_e) + 1) + 1;

        reductionValue_m =
            fMax(FL2FXCONST_DBL(0.f),
                 scaleValue(reductionValue_m, reductionValue_e - result_e) +
                     scaleValue(redVal_m, redVal_e - result_e) -
                     scaleValue(avgThrExp_m, avgThrExp_e - result_e));
        reductionValue_e = result_e;
      }
    } else {
      int sc0, sc1;
      reductionValue_m = fMultNorm(
          reductionValue_m,
          fDivNorm((FIXP_DBL)redPeGlobal, (FIXP_DBL)desiredPe, &sc0), &sc1);
      reductionValue_e += sc0 + sc1;

      for (elementId = elementOffset; elementId < nElements; elementId++) {
        if (cm->elInfo[elementId].elType != ID_DSE) {
          FDKaacEnc_resetAHFlags(pAhFlag[elementId],
                                 cm->elInfo[elementId].nChannelsInEl,
                                 psyOutElement[elementId]->psyOutChannel);
        }
      }
    }

    redPeGlobal = 0;
    for (elementId = elementOffset; elementId < nElements; elementId++) {
      if (cm->elInfo[elementId].elType != ID_DSE) {
        INT nChannels = cm->elInfo[elementId].nChannelsInEl;
        PE_DATA *peData = &qcElement[elementId]->peData;

        FDKaacEnc_reduceThresholdsCBR(
            qcElement[elementId]->qcOutChannel,
            psyOutElement[elementId]->psyOutChannel, pAhFlag[elementId],
            pThrExp[elementId], nChannels, reductionValue_m, reductionValue_e);

        FDKaacEnc_calcPe(psyOutElement[elementId]->psyOutChannel,
                         qcElement[elementId]->qcOutChannel, peData, nChannels);
        redPeGlobal += peData->pe;
      }
    }

    iter++;
  }

  if (redPeGlobal > desiredPe) {
    FDKaacEnc_correctThresh(cm, qcElement, psyOutElement, pAhFlag, pThrExp,
                            reductionValue_m, reductionValue_e,
                            desiredPe - redPeGlobal, processElements,
                            elementOffset);

    redPeGlobal = 0;
    for (elementId = elementOffset; elementId < nElements; elementId++) {
      if (cm->elInfo[elementId].elType != ID_DSE) {
        INT nChannels = cm->elInfo[elementId].nChannelsInEl;
        PE_DATA *peData = &qcElement[elementId]->peData;

        FDKaacEnc_calcPe(psyOutElement[elementId]->psyOutChannel,
                         qcElement[elementId]->qcOutChannel, peData, nChannels);
        redPeGlobal += peData->pe;
      }
    }
  }

  if (redPeGlobal > desiredPe) {
    FDKaacEnc_reduceMinSnr(
        cm, qcElement, psyOutElement, pAhFlag,
        (fMultI(FL2FXCONST_DBL(0.15f), desiredPe) + desiredPe), &redPeGlobal,
        processElements, elementOffset);

    FDKaacEnc_allowMoreHoles(cm, qcElement, psyOutElement, AdjThrStateElement,
                             pAhFlag, desiredPe, redPeGlobal, processElements,
                             elementOffset);
  }
}

static void FDKaacEnc_AdaptThresholdsVBR(
    QC_OUT_CHANNEL *const qcOutChannel[(2)],
    const PSY_OUT_CHANNEL *const psyOutChannel[(2)],
    ATS_ELEMENT *const AdjThrStateElement,
    const struct TOOLSINFO *const toolsInfo, const INT nChannels) {
  UCHAR(*pAhFlag)[MAX_GROUPED_SFB];
  FIXP_DBL(*pThrExp)[MAX_GROUPED_SFB];

  C_ALLOC_SCRATCH_START(_pAhFlag, UCHAR, (2) * MAX_GROUPED_SFB)
  C_ALLOC_SCRATCH_START(_pThrExp, FIXP_DBL, (2) * MAX_GROUPED_SFB)
  pAhFlag = (UCHAR(*)[MAX_GROUPED_SFB])_pAhFlag;
  pThrExp = (FIXP_DBL(*)[MAX_GROUPED_SFB])_pThrExp;

  FDKaacEnc_calcThreshExp(pThrExp, qcOutChannel, psyOutChannel, nChannels);

  FDKaacEnc_adaptMinSnr(qcOutChannel, psyOutChannel,
                        &AdjThrStateElement->minSnrAdaptParam, nChannels);

  FDKaacEnc_initAvoidHoleFlag(qcOutChannel, psyOutChannel, pAhFlag, toolsInfo,
                              nChannels, &AdjThrStateElement->ahParam);

  FDKaacEnc_reduceThresholdsVBR(qcOutChannel, psyOutChannel, pAhFlag, pThrExp,
                                nChannels, AdjThrStateElement->vbrQualFactor,
                                &AdjThrStateElement->chaosMeasureOld);

  C_ALLOC_SCRATCH_END(_pThrExp, FIXP_DBL, (2) * MAX_GROUPED_SFB)
  C_ALLOC_SCRATCH_END(_pAhFlag, UCHAR, (2) * MAX_GROUPED_SFB)
}

static FIXP_DBL FDKaacEnc_calcBitSave(FIXP_DBL fillLevel,
                                      const FIXP_DBL clipLow,
                                      const FIXP_DBL clipHigh,
                                      const FIXP_DBL minBitSave,
                                      const FIXP_DBL maxBitSave,
                                      const FIXP_DBL bitsave_slope) {
  FIXP_DBL bitsave;

  fillLevel = fixMax(fillLevel, clipLow);
  fillLevel = fixMin(fillLevel, clipHigh);

  bitsave = maxBitSave - fMult((fillLevel - clipLow), bitsave_slope);

  return (bitsave);
}

static FIXP_DBL FDKaacEnc_calcBitSpend(FIXP_DBL fillLevel,
                                       const FIXP_DBL clipLow,
                                       const FIXP_DBL clipHigh,
                                       const FIXP_DBL minBitSpend,
                                       const FIXP_DBL maxBitSpend,
                                       const FIXP_DBL bitspend_slope) {
  FIXP_DBL bitspend;

  fillLevel = fixMax(fillLevel, clipLow);
  fillLevel = fixMin(fillLevel, clipHigh);

  bitspend = minBitSpend + fMult(fillLevel - clipLow, bitspend_slope);

  return (bitspend);
}

static void FDKaacEnc_adjustPeMinMax(const INT currPe, INT *peMin, INT *peMax) {
  FIXP_DBL minFacHi = FL2FXCONST_DBL(0.3f), maxFacHi = (FIXP_DBL)MAXVAL_DBL,
           minFacLo = FL2FXCONST_DBL(0.14f), maxFacLo = FL2FXCONST_DBL(0.07f);
  INT diff;

  INT minDiff_fix = fMultI(FL2FXCONST_DBL(0.1666666667f), currPe);

  if (currPe > *peMax) {
    diff = (currPe - *peMax);
    *peMin += fMultI(minFacHi, diff);
    *peMax += fMultI(maxFacHi, diff);
  } else if (currPe < *peMin) {
    diff = (*peMin - currPe);
    *peMin -= fMultI(minFacLo, diff);
    *peMax -= fMultI(maxFacLo, diff);
  } else {
    *peMin += fMultI(minFacHi, (currPe - *peMin));
    *peMax -= fMultI(maxFacLo, (*peMax - currPe));
  }

  if ((*peMax - *peMin) < minDiff_fix) {
    INT peMax_fix = *peMax, peMin_fix = *peMin;
    FIXP_DBL partLo_fix, partHi_fix;

    partLo_fix = (FIXP_DBL)fixMax(0, currPe - peMin_fix);
    partHi_fix = (FIXP_DBL)fixMax(0, peMax_fix - currPe);

    peMax_fix =
        (INT)(currPe + fMultI(fDivNorm(partHi_fix, (partLo_fix + partHi_fix)),
                              minDiff_fix));
    peMin_fix =
        (INT)(currPe - fMultI(fDivNorm(partLo_fix, (partLo_fix + partHi_fix)),
                              minDiff_fix));
    peMin_fix = fixMax(0, peMin_fix);

    *peMax = peMax_fix;
    *peMin = peMin_fix;
  }
}

void FDKaacEnc_bitresCalcBitFac(const INT bitresBits, const INT maxBitresBits,
                                const INT pe, const INT lastWindowSequence,
                                const INT avgBits, const FIXP_DBL maxBitFac,
                                const ADJ_THR_STATE *const AdjThr,
                                ATS_ELEMENT *const adjThrChan,
                                FIXP_DBL *const pBitresFac,
                                INT *const pBitresFac_e) {
  const BRES_PARAM *bresParam;
  INT pex;
  FIXP_DBL fillLevel;
  INT fillLevel_e = 0;

  FIXP_DBL bitresFac;
  INT bitresFac_e;

  FIXP_DBL bitSave, bitSpend;
  FIXP_DBL bitsave_slope, bitspend_slope;
  FIXP_DBL fillLevel_fix = MAXVAL_DBL;

  FIXP_DBL slope = MAXVAL_DBL;

  if (lastWindowSequence != SHORT_WINDOW) {
    bresParam = &(AdjThr->bresParamLong);
    bitsave_slope = FL2FXCONST_DBL(0.466666666);
    bitspend_slope = FL2FXCONST_DBL(0.666666666);
  } else {
    bresParam = &(AdjThr->bresParamShort);
    bitsave_slope = (FIXP_DBL)0x2E8BA2E9;
    bitspend_slope = (FIXP_DBL)0x7fffffff;
  }

  if (bitresBits < maxBitresBits) {
    fillLevel_fix = fDivNorm(bitresBits, maxBitresBits);
  }

  pex = fMax(pe, adjThrChan->peMin);
  pex = fMin(pex, adjThrChan->peMax);

  bitSave = FDKaacEnc_calcBitSave(
      fillLevel_fix, bresParam->clipSaveLow, bresParam->clipSaveHigh,
      bresParam->minBitSave, bresParam->maxBitSave, bitsave_slope);

  bitSpend = FDKaacEnc_calcBitSpend(
      fillLevel_fix, bresParam->clipSpendLow, bresParam->clipSpendHigh,
      bresParam->minBitSpend, bresParam->maxBitSpend, bitspend_slope);

  slope = schur_div((pex - adjThrChan->peMin),
                    (adjThrChan->peMax - adjThrChan->peMin), 31);

  bitresFac = ((FIXP_DBL)(MAXVAL_DBL >> 1) - (bitSave >> 1));
  bitresFac_e = 1;
  bitresFac = fMultAddDiv2(bitresFac, slope, bitSpend + bitSave);

  fillLevel = fDivNorm(bitresBits, avgBits, &fillLevel_e);
  if (fillLevel_e < 0) {
    fillLevel = scaleValue(fillLevel, fillLevel_e);
    fillLevel_e = 0;
  }
  fillLevel >>= 1;
  fillLevel_e += 1;
  fillLevel += scaleValue(FL2FXCONST_DBL(0.7f), -fillLevel_e);
  if (scaleValue(bitresFac, -fillLevel_e + 1) > fillLevel) {
    bitresFac = fillLevel;
    bitresFac_e = fillLevel_e;
  }

  if (scaleValue(bitresFac, bitresFac_e - (DFRACT_BITS - 1 - 24)) > maxBitFac) {
    bitresFac = maxBitFac;
    bitresFac_e = (DFRACT_BITS - 1 - 24);
  }

  FDKaacEnc_adjustPeMinMax(pe, &adjThrChan->peMin, &adjThrChan->peMax);

  *pBitresFac = bitresFac;
  *pBitresFac_e = bitresFac_e;
}

INT FDKaacEnc_AdjThrNew(ADJ_THR_STATE **phAdjThr, INT nElements) {
  INT err = 0;
  INT i;
  ADJ_THR_STATE *hAdjThr = GetRam_aacEnc_AdjustThreshold();
  if (hAdjThr == NULL) {
    err = 1;
    goto bail;
  }

  for (i = 0; i < nElements; i++) {
    hAdjThr->adjThrStateElem[i] = GetRam_aacEnc_AdjThrStateElement(i);
    if (hAdjThr->adjThrStateElem[i] == NULL) {
      err = 1;
      goto bail;
    }
  }

bail:
  *phAdjThr = hAdjThr;
  return err;
}

void FDKaacEnc_AdjThrInit(
    ADJ_THR_STATE *const hAdjThr, const INT meanPe, const INT invQuant,
    const CHANNEL_MAPPING *const channelMapping, const INT sampleRate,
    const INT totalBitrate, const INT isLowDelay,
    const AACENC_BITRES_MODE bitResMode, const INT dZoneQuantEnable,
    const INT bitDistributionMode, const FIXP_DBL vbrQualFactor) {
  INT i;

  FIXP_DBL POINT8 = FL2FXCONST_DBL(0.8f);
  FIXP_DBL POINT6 = FL2FXCONST_DBL(0.6f);

  if (bitDistributionMode == 1) {
    hAdjThr->bitDistributionMode = AACENC_BD_MODE_INTRA_ELEMENT;
  } else {
    hAdjThr->bitDistributionMode = AACENC_BD_MODE_INTER_ELEMENT;
  }

  hAdjThr->maxIter2ndGuess =
      (isLowDelay != 0 || channelMapping->nElements > 1) ? 3 : 1;

  hAdjThr->bresParamLong.clipSaveLow =
      (FIXP_DBL)0x1999999a;
  hAdjThr->bresParamLong.clipSaveHigh =
      (FIXP_DBL)0x7999999a;
  hAdjThr->bresParamLong.minBitSave =
      (FIXP_DBL)0xf999999a;
  hAdjThr->bresParamLong.maxBitSave =
      (FIXP_DBL)0x26666666;
  hAdjThr->bresParamLong.clipSpendLow =
      (FIXP_DBL)0x1999999a;
  hAdjThr->bresParamLong.clipSpendHigh =
      (FIXP_DBL)0x7999999a;
  hAdjThr->bresParamLong.minBitSpend =
      (FIXP_DBL)0xf3333333;
  hAdjThr->bresParamLong.maxBitSpend =
      (FIXP_DBL)0x33333333;

  hAdjThr->bresParamShort.clipSaveLow =
      (FIXP_DBL)0x199999a0;
  hAdjThr->bresParamShort.clipSaveHigh =
      (FIXP_DBL)0x5fffffff;
  hAdjThr->bresParamShort.minBitSave =
      (FIXP_DBL)0x00000000;
  hAdjThr->bresParamShort.maxBitSave =
      (FIXP_DBL)0x199999a0;
  hAdjThr->bresParamShort.clipSpendLow =
      (FIXP_DBL)0x199999a0;
  hAdjThr->bresParamShort.clipSpendHigh =
      (FIXP_DBL)0x5fffffff;
  hAdjThr->bresParamShort.minBitSpend =
      (FIXP_DBL)0xf9999998;
  hAdjThr->bresParamShort.maxBitSpend =
      (FIXP_DBL)0x40000000;

  for (i = 0; i < channelMapping->nElements; i++) {
    const FIXP_DBL relativeBits = channelMapping->elInfo[i].relativeBits;
    const INT nChannelsInElement = channelMapping->elInfo[i].nChannelsInEl;
    const INT bitrateInElement =
        (relativeBits != (FIXP_DBL)MAXVAL_DBL)
            ? (INT)fMultNorm(relativeBits, (FIXP_DBL)totalBitrate)
            : totalBitrate;
    const INT chBitrate = bitrateInElement >> (nChannelsInElement == 1 ? 0 : 1);

    ATS_ELEMENT *atsElem = hAdjThr->adjThrStateElem[i];
    MINSNR_ADAPT_PARAM *msaParam = &atsElem->minSnrAdaptParam;

    if (isLowDelay) {
      atsElem->peMin = fMultI(POINT8, meanPe);
      atsElem->peMax = fMultI(POINT6, meanPe) << 1;
    } else {
      atsElem->peMin = fMultI(POINT8, meanPe) >> 1;
      atsElem->peMax = fMultI(POINT6, meanPe);
    }

    atsElem->chaosMeasureOld = FL2FXCONST_DBL(0.3f);

    atsElem->peOffset = 0;

    atsElem->vbrQualFactor = vbrQualFactor;
    if (chBitrate < 32000) {
      atsElem->peOffset =
          fixMax(50, 100 - fMultI((FIXP_DBL)0x666667, chBitrate));
    }

    if (chBitrate >= 20000) {
      atsElem->ahParam.modifyMinSnr = TRUE;
      atsElem->ahParam.startSfbL = 15;
      atsElem->ahParam.startSfbS = 3;
    } else {
      atsElem->ahParam.modifyMinSnr = FALSE;
      atsElem->ahParam.startSfbL = 0;
      atsElem->ahParam.startSfbS = 0;
    }

    msaParam->maxRed = FL2FXCONST_DBL(0.00390625f);
    msaParam->startRatio = FL2FXCONST_DBL(0.05190512648f);
    msaParam->redRatioFac = FL2FXCONST_DBL(-0.375f);
    msaParam->redOffs = FL2FXCONST_DBL(0.021484375);

    atsElem->peCorrectionFactor_m = FL2FXCONST_DBL(0.5f);
    atsElem->peCorrectionFactor_e = 1;

    atsElem->dynBitsLast = -1;
    atsElem->peLast = 0;

    FDKaacEnc_InitBits2PeFactor(
        &atsElem->bits2PeFactor_m, &atsElem->bits2PeFactor_e, bitrateInElement,
        nChannelsInElement, sampleRate, isLowDelay, dZoneQuantEnable, invQuant);
  }
}

static void FDKaacEnc_FDKaacEnc_calcPeCorrection(
    FIXP_DBL *const correctionFac_m, INT *const correctionFac_e,
    const INT peAct, const INT peLast, const INT bitsLast,
    const FIXP_DBL bits2PeFactor_m, const INT bits2PeFactor_e) {
  if ((bitsLast > 0) && (peAct < 1.5f * peLast) && (peAct > 0.7f * peLast) &&
      (FDKaacEnc_bits2pe2(bitsLast,
                          fMult(FL2FXCONST_DBL(1.2f / 2.f), bits2PeFactor_m),
                          bits2PeFactor_e + 1) > peLast) &&
      (FDKaacEnc_bits2pe2(bitsLast,
                          fMult(FL2FXCONST_DBL(0.65f), bits2PeFactor_m),
                          bits2PeFactor_e) < peLast)) {
    FIXP_DBL corrFac = *correctionFac_m;

    int scaling = 0;
    FIXP_DBL denum = (FIXP_DBL)FDKaacEnc_bits2pe2(bitsLast, bits2PeFactor_m,
                                                  bits2PeFactor_e);
    FIXP_DBL newFac = fDivNorm((FIXP_DBL)peLast, denum, &scaling);

    if ((FIXP_DBL)peLast <= denum) {
      newFac = fixMax(
          scaleValue(fixMin(fMult(FL2FXCONST_DBL(1.1f / 2.f), newFac),
                            scaleValue(FL2FXCONST_DBL(1.f / 2.f), -scaling)),
                     scaling),
          FL2FXCONST_DBL(0.85f / 2.f));
    } else {
      newFac = fixMax(
          fixMin(scaleValue(fMult(FL2FXCONST_DBL(0.9f / 2.f), newFac), scaling),
                 FL2FXCONST_DBL(1.15f / 2.f)),
          FL2FXCONST_DBL(1.f / 2.f));
    }

    if (((newFac > FL2FXCONST_DBL(1.f / 2.f)) &&
         (corrFac < FL2FXCONST_DBL(1.f / 2.f))) ||
        ((newFac < FL2FXCONST_DBL(1.f / 2.f)) &&
         (corrFac > FL2FXCONST_DBL(1.f / 2.f)))) {
      corrFac = FL2FXCONST_DBL(1.f / 2.f);
    }

    if ((corrFac < FL2FXCONST_DBL(1.f / 2.f) && newFac < corrFac) ||
        (corrFac > FL2FXCONST_DBL(1.f / 2.f) && newFac > corrFac)) {
      corrFac = fMult(FL2FXCONST_DBL(0.85f), corrFac) +
                fMult(FL2FXCONST_DBL(0.15f), newFac);
    } else {
      corrFac = fMult(FL2FXCONST_DBL(0.7f), corrFac) +
                fMult(FL2FXCONST_DBL(0.3f), newFac);
    }

    corrFac = fixMax(fixMin(corrFac, FL2FXCONST_DBL(1.15f / 2.f)),
                     FL2FXCONST_DBL(0.85 / 2.f));

    *correctionFac_m = corrFac;
    *correctionFac_e = 1;
  } else {
    *correctionFac_m = FL2FXCONST_DBL(1.f / 2.f);
    *correctionFac_e = 1;
  }
}

static void FDKaacEnc_calcPeCorrectionLowBitRes(
    FIXP_DBL *const correctionFac_m, INT *const correctionFac_e,
    const INT peLast, const INT bitsLast, const INT bitresLevel,
    const INT nChannels, const FIXP_DBL bits2PeFactor_m,
    const INT bits2PeFactor_e) {
  const FIXP_DBL amp = FL2FXCONST_DBL(0.005);
  const FIXP_DBL maxDiff = FL2FXCONST_DBL(0.25f);

  if (bitsLast > 0) {
    const int bitsBalLast =
        peLast - FDKaacEnc_bits2pe2(bitsLast, bits2PeFactor_m, bits2PeFactor_e);

    int headroom = (bitresLevel >= 50 * nChannels) ? 0 : (100 * nChannels);

    headroom = FDKaacEnc_bits2pe2(headroom, bits2PeFactor_m, bits2PeFactor_e);

    FIXP_DBL denominator = (FIXP_DBL)FDKaacEnc_bits2pe2(
                               bitresLevel, bits2PeFactor_m, bits2PeFactor_e) +
                           (FIXP_DBL)headroom;

    int scaling = 0;
    FIXP_DBL diff =
        (bitsBalLast >= headroom)
            ? fMult(amp, fDivNorm((FIXP_DBL)(bitsBalLast - headroom),
                                  denominator, &scaling))
            : -fMult(amp, fDivNorm(-(FIXP_DBL)(bitsBalLast - headroom),
                                   denominator, &scaling));
    scaling -= 1;

    diff = (scaling <= 0)
               ? fMax(fMin(diff >> (-scaling), maxDiff >> 1), -maxDiff >> 1)
               : fMax(fMin(diff, maxDiff >> (1 + scaling)),
                      -maxDiff >> (1 + scaling))
                     << scaling;

    *correctionFac_m =
        fMax(fMin((*correctionFac_m) + diff, FL2FXCONST_DBL(1.0f / 2.f)),
             FL2FXCONST_DBL(0.75f / 2.f));
    *correctionFac_e = 1;
  } else {
    *correctionFac_m = FL2FXCONST_DBL(0.75 / 2.f);
    *correctionFac_e = 1;
  }
}

void FDKaacEnc_DistributeBits(
    ADJ_THR_STATE *adjThrState, ATS_ELEMENT *AdjThrStateElement,
    PSY_OUT_CHANNEL *psyOutChannel[(2)], PE_DATA *peData, INT *grantedPe,
    INT *grantedPeCorr, const INT nChannels, const INT commonWindow,
    const INT grantedDynBits, const INT bitresBits, const INT maxBitresBits,
    const FIXP_DBL maxBitFac, const AACENC_BITRES_MODE bitResMode) {
  FIXP_DBL bitFactor;
  INT bitFactor_e;
  INT noRedPe = peData->pe;

  INT curWindowSequence = LONG_WINDOW;
  if (nChannels == 2) {
    if ((psyOutChannel[0]->lastWindowSequence == SHORT_WINDOW) ||
        (psyOutChannel[1]->lastWindowSequence == SHORT_WINDOW)) {
      curWindowSequence = SHORT_WINDOW;
    }
  } else {
    curWindowSequence = psyOutChannel[0]->lastWindowSequence;
  }

  if (grantedDynBits >= 1) {
    if (bitResMode != AACENC_BR_MODE_FULL) {
      *grantedPe = FDKaacEnc_bits2pe2(grantedDynBits,
                                      AdjThrStateElement->bits2PeFactor_m,
                                      AdjThrStateElement->bits2PeFactor_e);
    } else {
      FDKaacEnc_bitresCalcBitFac(
          bitresBits, maxBitresBits, noRedPe, curWindowSequence, grantedDynBits,
          maxBitFac, adjThrState, AdjThrStateElement, &bitFactor, &bitFactor_e);

      *grantedPe = FDKaacEnc_bits2pe2(
          grantedDynBits, fMult(bitFactor, AdjThrStateElement->bits2PeFactor_m),
          AdjThrStateElement->bits2PeFactor_e + bitFactor_e);
    }
  } else {
    *grantedPe = 0;
  }

  switch (bitResMode) {
    case AACENC_BR_MODE_DISABLED:
    case AACENC_BR_MODE_REDUCED:
      FDKaacEnc_calcPeCorrectionLowBitRes(
          &AdjThrStateElement->peCorrectionFactor_m,
          &AdjThrStateElement->peCorrectionFactor_e, AdjThrStateElement->peLast,
          AdjThrStateElement->dynBitsLast, bitresBits, nChannels,
          AdjThrStateElement->bits2PeFactor_m,
          AdjThrStateElement->bits2PeFactor_e);
      break;
    case AACENC_BR_MODE_FULL:
    default:
      FDKaacEnc_FDKaacEnc_calcPeCorrection(
          &AdjThrStateElement->peCorrectionFactor_m,
          &AdjThrStateElement->peCorrectionFactor_e,
          fixMin(*grantedPe, noRedPe), AdjThrStateElement->peLast,
          AdjThrStateElement->dynBitsLast, AdjThrStateElement->bits2PeFactor_m,
          AdjThrStateElement->bits2PeFactor_e);
      break;
  }

  *grantedPeCorr =
      (INT)(fMult((FIXP_DBL)(*grantedPe << Q_AVGBITS),
                  AdjThrStateElement->peCorrectionFactor_m) >>
            (Q_AVGBITS - AdjThrStateElement->peCorrectionFactor_e));

  AdjThrStateElement->peLast = *grantedPe;
  AdjThrStateElement->dynBitsLast = -1;
}

void FDKaacEnc_AdjustThresholds(
    ADJ_THR_STATE *const hAdjThr, QC_OUT_ELEMENT *const qcElement[((8))],
    QC_OUT *const qcOut, const PSY_OUT_ELEMENT *const psyOutElement[((8))],
    const INT CBRbitrateMode, const CHANNEL_MAPPING *const cm) {
  int i;

  if (CBRbitrateMode) {
    if (hAdjThr->bitDistributionMode == AACENC_BD_MODE_INTRA_ELEMENT) {
      for (i = 0; i < cm->nElements; i++) {
        ELEMENT_INFO elInfo = cm->elInfo[i];

        if ((elInfo.elType == ID_SCE) || (elInfo.elType == ID_CPE) ||
            (elInfo.elType == ID_LFE)) {
          if (qcElement[i]->grantedPeCorr < qcElement[i]->peData.pe) {
            FDKaacEnc_adaptThresholdsToPe(
                cm, hAdjThr->adjThrStateElem, qcElement, psyOutElement,
                qcElement[i]->grantedPeCorr, hAdjThr->maxIter2ndGuess,
                1,
                i
            );
          }
        }
      }
    }
    else if (hAdjThr->bitDistributionMode == AACENC_BD_MODE_INTER_ELEMENT) {
      if (qcOut->totalGrantedPeCorr < qcOut->totalNoRedPe) {
        FDKaacEnc_adaptThresholdsToPe(cm, hAdjThr->adjThrStateElem, qcElement,
                                      psyOutElement, qcOut->totalGrantedPeCorr,
                                      hAdjThr->maxIter2ndGuess,
                                      cm->nElements,
                                      0);
      } else {
        for (i = 0; i < cm->nElements; i++) {
          if ((cm->elInfo[i].elType == ID_SCE) ||
              (cm->elInfo[i].elType == ID_CPE) ||
              (cm->elInfo[i].elType == ID_LFE)) {
            const int maxElementPe = FDKaacEnc_bits2pe2(
                (cm->elInfo[i].nChannelsInEl * MIN_BUFSIZE_PER_EFF_CHAN) -
                    qcElement[i]->staticBitsUsed - qcElement[i]->extBitsUsed,
                hAdjThr->adjThrStateElem[i]->bits2PeFactor_m,
                hAdjThr->adjThrStateElem[i]->bits2PeFactor_e);

            if (maxElementPe < qcElement[i]->peData.pe) {
              FDKaacEnc_adaptThresholdsToPe(
                  cm, hAdjThr->adjThrStateElem, qcElement, psyOutElement,
                  maxElementPe, hAdjThr->maxIter2ndGuess, 1, i);
            }
          }
        }
      }
    }
  } else {
    for (i = 0; i < cm->nElements; i++) {
      ELEMENT_INFO elInfo = cm->elInfo[i];

      if ((elInfo.elType == ID_SCE) || (elInfo.elType == ID_CPE) ||
          (elInfo.elType == ID_LFE)) {
        FDKaacEnc_AdaptThresholdsVBR(
            qcElement[i]->qcOutChannel, psyOutElement[i]->psyOutChannel,
            hAdjThr->adjThrStateElem[i], &psyOutElement[i]->toolsInfo,
            cm->elInfo[i].nChannelsInEl);
      }
    }
  }
  for (i = 0; i < cm->nElements; i++) {
    int ch, sfb, sfbGrp;
    for (ch = 0; ch < cm->elInfo[i].nChannelsInEl; ch++) {
      QC_OUT_CHANNEL *pQcOutCh = qcElement[i]->qcOutChannel[ch];
      for (sfbGrp = 0; sfbGrp < psyOutElement[i]->psyOutChannel[ch]->sfbCnt;
           sfbGrp += psyOutElement[i]->psyOutChannel[ch]->sfbPerGroup) {
        for (sfb = 0; sfb < psyOutElement[i]->psyOutChannel[ch]->maxSfbPerGroup;
             sfb++) {
          pQcOutCh->sfbThresholdLdData[sfb + sfbGrp] +=
              pQcOutCh->sfbEnFacLd[sfb + sfbGrp];
        }
      }
    }
  }
}

void FDKaacEnc_AdjThrClose(ADJ_THR_STATE **phAdjThr) {
  INT i;
  ADJ_THR_STATE *hAdjThr = *phAdjThr;

  if (hAdjThr != NULL) {
    for (i = 0; i < ((8)); i++) {
      if (hAdjThr->adjThrStateElem[i] != NULL) {
        FreeRam_aacEnc_AdjThrStateElement(&hAdjThr->adjThrStateElem[i]);
      }
    }
    FreeRam_aacEnc_AdjustThreshold(phAdjThr);
  }
}
