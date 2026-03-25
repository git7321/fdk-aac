/**************************** AAC encoder library ******************************

   Author(s):   M. Neusinger

   Description: Compressor for AAC Metadata Generator

*******************************************************************************/

#include "metadata_compressor.h"
#include "channel_map.h"

#define LOG2 0.69314718056f
#define ILOG2 1.442695041f
#define FIXP_ILOG2_DIV2 (FL2FXCONST_DBL(ILOG2 / 2))

#define MAX_DRC_CHANNELS (8)
#define DOWNMIX_SHIFT (3)
#define WEIGHTING_FILTER_SHIFT (2)

#define METADATA_INT_BITS 10
#define METADATA_LINT_BITS 20
#define METADATA_INT_SCALE (INT64(1) << (METADATA_INT_BITS))
#define METADATA_FRACT_BITS (DFRACT_BITS - 1 - METADATA_INT_BITS)
#define METADATA_FRACT_SCALE (INT64(1) << (METADATA_FRACT_BITS))

enum { L = 0, R = 1, C = 2, LFE = 3, LS = 4, RS = 5, S = 6, LS2 = 7, RS2 = 8 };

struct WEIGHTING_STATES {
  FIXP_DBL x1;
  FIXP_DBL x2;
  FIXP_DBL y1;
  FIXP_DBL y2;
};

struct DRC_COMP {
  FIXP_DBL maxBoostThr[2];
  FIXP_DBL boostThr[2];
  FIXP_DBL earlyCutThr[2];
  FIXP_DBL cutThr[2];
  FIXP_DBL maxCutThr[2];

  FIXP_DBL boostFac[2];
  FIXP_DBL
  earlyCutFac[2];
  FIXP_DBL cutFac[2];

  FIXP_DBL maxBoost[2];
  FIXP_DBL maxCut[2];
  FIXP_DBL maxEarlyCut[2];

  FIXP_DBL fastAttack[2];
  FIXP_DBL fastDecay[2];
  FIXP_DBL slowAttack[2];
  FIXP_DBL slowDecay[2];
  UINT holdOff[2];

  FIXP_DBL attackThr[2];
  FIXP_DBL decayThr[2];

  DRC_PROFILE profile[2];
  INT blockLength;
  UINT sampleRate;
  CHANNEL_MODE chanConfig;

  UCHAR useWeighting;

  UINT channels;
  UINT fullChannels;
  INT channelIdx[9];

  FIXP_DBL smoothLevel[2];
  FIXP_DBL smoothGain[2];
  UINT holdCnt[2];

  FIXP_DBL limGain[2];
  FIXP_DBL limDecay;
  FIXP_DBL prevPeak[2];

  WEIGHTING_STATES
  filter[MAX_DRC_CHANNELS];
};

static const FIXP_DBL tabMaxBoostThr[] = {
    (FIXP_DBL)(-(43 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(53 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(55 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(65 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(50 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(40 << METADATA_FRACT_BITS))};
static const FIXP_DBL tabBoostThr[] = {
    (FIXP_DBL)(-(31 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(41 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(31 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(41 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(31 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(31 << METADATA_FRACT_BITS))};
static const FIXP_DBL tabEarlyCutThr[] = {
    (FIXP_DBL)(-(26 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(21 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(26 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(21 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(26 << METADATA_FRACT_BITS)),
    (FIXP_DBL)(-(20 << METADATA_FRACT_BITS))};
static const FIXP_DBL tabCutThr[] = {(FIXP_DBL)(-(16 << METADATA_FRACT_BITS)),
                                     (FIXP_DBL)(-(11 << METADATA_FRACT_BITS)),
                                     (FIXP_DBL)(-(16 << METADATA_FRACT_BITS)),
                                     (FIXP_DBL)(-(21 << METADATA_FRACT_BITS)),
                                     (FIXP_DBL)(-(16 << METADATA_FRACT_BITS)),
                                     (FIXP_DBL)(-(10 << METADATA_FRACT_BITS))};
static const FIXP_DBL tabMaxCutThr[] = {
    (FIXP_DBL)(4 << METADATA_FRACT_BITS), (FIXP_DBL)(9 << METADATA_FRACT_BITS),
    (FIXP_DBL)(4 << METADATA_FRACT_BITS), (FIXP_DBL)(9 << METADATA_FRACT_BITS),
    (FIXP_DBL)(4 << METADATA_FRACT_BITS), (FIXP_DBL)(4 << METADATA_FRACT_BITS)};
static const FIXP_DBL tabBoostRatio[] = {
    FL2FXCONST_DBL(((1.f / 2.f) - 1.f)), FL2FXCONST_DBL(((1.f / 2.f) - 1.f)),
    FL2FXCONST_DBL(((1.f / 2.f) - 1.f)), FL2FXCONST_DBL(((1.f / 2.f) - 1.f)),
    FL2FXCONST_DBL(((1.f / 5.f) - 1.f)), FL2FXCONST_DBL(((1.f / 5.f) - 1.f))};
static const FIXP_DBL tabEarlyCutRatio[] = {
    FL2FXCONST_DBL(((1.f / 2.f) - 1.f)), FL2FXCONST_DBL(((1.f / 2.f) - 1.f)),
    FL2FXCONST_DBL(((1.f / 2.f) - 1.f)), FL2FXCONST_DBL(((1.f / 1.f) - 1.f)),
    FL2FXCONST_DBL(((1.f / 2.f) - 1.f)), FL2FXCONST_DBL(((1.f / 2.f) - 1.f))};
static const FIXP_DBL tabCutRatio[] = {
    FL2FXCONST_DBL(((1.f / 20.f) - 1.f)), FL2FXCONST_DBL(((1.f / 20.f) - 1.f)),
    FL2FXCONST_DBL(((1.f / 20.f) - 1.f)), FL2FXCONST_DBL(((1.f / 2.f) - 1.f)),
    FL2FXCONST_DBL(((1.f / 20.f) - 1.f)), FL2FXCONST_DBL(((1.f / 20.f) - 1.f))};
static const FIXP_DBL tabMaxBoost[] = {(FIXP_DBL)(6 << METADATA_FRACT_BITS),
                                       (FIXP_DBL)(6 << METADATA_FRACT_BITS),
                                       (FIXP_DBL)(12 << METADATA_FRACT_BITS),
                                       (FIXP_DBL)(12 << METADATA_FRACT_BITS),
                                       (FIXP_DBL)(15 << METADATA_FRACT_BITS),
                                       (FIXP_DBL)(15 << METADATA_FRACT_BITS)};
static const FIXP_DBL tabMaxCut[] = {(FIXP_DBL)(24 << METADATA_FRACT_BITS),
                                     (FIXP_DBL)(24 << METADATA_FRACT_BITS),
                                     (FIXP_DBL)(24 << METADATA_FRACT_BITS),
                                     (FIXP_DBL)(15 << METADATA_FRACT_BITS),
                                     (FIXP_DBL)(24 << METADATA_FRACT_BITS),
                                     (FIXP_DBL)(24 << METADATA_FRACT_BITS)};
static const FIXP_DBL tabFastAttack[] = {
    FL2FXCONST_DBL((10.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((10.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((10.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((10.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((10.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((0.f / 1000.f) / METADATA_INT_SCALE)};
static const FIXP_DBL tabFastDecay[] = {
    FL2FXCONST_DBL((1000.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((1000.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((1000.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((1000.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((200.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((0.f / 1000.f) / METADATA_INT_SCALE)};
static const FIXP_DBL tabSlowAttack[] = {
    FL2FXCONST_DBL((100.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((100.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((100.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((100.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((100.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((0.f / 1000.f) / METADATA_INT_SCALE)};
static const FIXP_DBL tabSlowDecay[] = {
    FL2FXCONST_DBL((3000.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((3000.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((10000.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((3000.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((1000.f / 1000.f) / METADATA_INT_SCALE),
    FL2FXCONST_DBL((0.f / 1000.f) / METADATA_INT_SCALE)};

static const INT tabHoldOff[] = {10, 10, 10, 10, 10, 0};

static const FIXP_DBL tabAttackThr[] = {(FIXP_DBL)(15 << METADATA_FRACT_BITS),
                                        (FIXP_DBL)(15 << METADATA_FRACT_BITS),
                                        (FIXP_DBL)(15 << METADATA_FRACT_BITS),
                                        (FIXP_DBL)(15 << METADATA_FRACT_BITS),
                                        (FIXP_DBL)(10 << METADATA_FRACT_BITS),
                                        (FIXP_DBL)(0 << METADATA_FRACT_BITS)};
static const FIXP_DBL tabDecayThr[] = {(FIXP_DBL)(20 << METADATA_FRACT_BITS),
                                       (FIXP_DBL)(20 << METADATA_FRACT_BITS),
                                       (FIXP_DBL)(20 << METADATA_FRACT_BITS),
                                       (FIXP_DBL)(20 << METADATA_FRACT_BITS),
                                       (FIXP_DBL)(10 << METADATA_FRACT_BITS),
                                       (FIXP_DBL)(0 << METADATA_FRACT_BITS)};

static const FIXP_DBL b0 = FL2FXCONST_DBL(0.53050662f);
static const FIXP_DBL a1 = FL2FXCONST_DBL(-0.95237983f),
                      a2 = FL2FXCONST_DBL(-0.02248836f);

static UINT getShiftFactor(const UINT length) {
  UINT ldN;
  for (ldN = 1; (((UINT)1) << ldN) < length; ldN++)
    ;

  return ldN;
}

static void fixpAdd(const FIXP_DBL value1, const int q1,
                    FIXP_DBL* const pValue2, int* const pQ2) {
  const int headroom1 = fNormz(fixp_abs(value1)) - 1;
  const int headroom2 = fNormz(fixp_abs(*pValue2)) - 1;
  int resultScale = fixMax(q1 - headroom1, (*pQ2) - headroom2);

  if ((value1 != FL2FXCONST_DBL(0.f)) && (*pValue2 != FL2FXCONST_DBL(0.f))) {
    resultScale++;
  }

  *pValue2 = scaleValue(value1, q1 - resultScale) +
             scaleValue(*pValue2, (*pQ2) - resultScale);
  *pQ2 = (*pValue2 != (FIXP_DBL)0) ? resultScale : DFRACT_BITS - 1;
}

static FIXP_DBL tc2Coeff(const FIXP_DBL t, const INT sampleRate,
                         const INT blockLength) {
  FIXP_DBL sampleRateFract;
  FIXP_DBL blockLengthFract;
  FIXP_DBL f, product;
  FIXP_DBL exponent, result;
  INT e_res;

  sampleRateFract =
      (FIXP_DBL)(sampleRate << (DFRACT_BITS - 1 - METADATA_LINT_BITS));
  blockLengthFract =
      (FIXP_DBL)(blockLength << (DFRACT_BITS - 1 - METADATA_LINT_BITS));
  f = fDivNorm(sampleRateFract, blockLengthFract, &e_res);
  f = scaleValue(f, e_res - METADATA_INT_BITS);

  product = fMultNorm(t, f, &e_res);
  product = scaleValue(
      product, e_res + METADATA_INT_BITS);

  exponent = fDivNorm(METADATA_FRACT_SCALE, product, &e_res);
  exponent = scaleValue(
      exponent, e_res - METADATA_INT_BITS);

  exponent = fMult(exponent, FIXP_ILOG2_DIV2) << 1;

  result = f2Pow(-exponent, DFRACT_BITS - 1 - METADATA_FRACT_BITS, &e_res);

  result = (FIXP_DBL)MAXVAL_DBL - scaleValue(result, e_res);

  return result;
}

static void findPeakLevels(HDRC_COMP drcComp, const INT_PCM* const inSamples,
                           const FIXP_DBL clev, const FIXP_DBL slev,
                           const FIXP_DBL ext_leva, const FIXP_DBL ext_levb,
                           const FIXP_DBL lfe_lev, const FIXP_DBL dmxGain5,
                           const FIXP_DBL dmxGain2, FIXP_DBL peak[2]) {
  int i, c;
  FIXP_DBL tmp = FL2FXCONST_DBL(0.f);
  INT_PCM maxSample = 0;

  peak[0] = peak[1] = FL2FXCONST_DBL(0.f);
  for (i = 0; i < drcComp->blockLength; i++) {
    const INT_PCM* pSamples = &inSamples[i * drcComp->channels];

    for (c = 0; c < (int)drcComp->channels; c++) {
      maxSample = fMax(maxSample, (INT_PCM)fAbs(pSamples[c]));
    }
  }
  peak[0] = fixMax(peak[0], FX_PCM2FX_DBL(maxSample) >> DOWNMIX_SHIFT);

  if (drcComp->fullChannels > 5) {
    for (i = 0; i < drcComp->blockLength; i++) {
      const INT_PCM* pSamples = &inSamples[i * drcComp->channels];

      tmp = FL2FXCONST_DBL(0.f);
      switch (drcComp->chanConfig) {
        case MODE_6_1:
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[4]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[6]]) >>
              (DOWNMIX_SHIFT - 1);
          break;
        case MODE_7_1_BACK:
        case MODE_7_1_REAR_SURROUND:
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[4]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
              (DOWNMIX_SHIFT - 1);
          break;
        case MODE_1_2_2_2_1:
        case MODE_7_1_FRONT_CENTER:
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[0]]) >>
                  DOWNMIX_SHIFT);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
              (DOWNMIX_SHIFT - 1);
          break;
        case MODE_7_1_TOP_FRONT:
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[0]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
              (DOWNMIX_SHIFT - 1);
          break;
        default:
          break;
      }
      peak[0] = fixMax(peak[0], fixp_abs(tmp));

      tmp = FL2FXCONST_DBL(0.f);
      switch (drcComp->chanConfig) {
        case MODE_6_1:
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[5]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[6]]) >>
              (DOWNMIX_SHIFT - 1);
          break;
        case MODE_7_1_BACK:
        case MODE_7_1_REAR_SURROUND:
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[5]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
              (DOWNMIX_SHIFT - 1);
          break;
        case MODE_1_2_2_2_1:
        case MODE_7_1_FRONT_CENTER:
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[1]]) >>
                  DOWNMIX_SHIFT);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
              (DOWNMIX_SHIFT - 1);
          break;
        case MODE_7_1_TOP_FRONT:
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[1]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
              (DOWNMIX_SHIFT - 1);
          break;
        default:
          break;
      }
      peak[0] = fixMax(peak[0], fixp_abs(tmp));

      tmp = FL2FXCONST_DBL(0.f);
      switch (drcComp->chanConfig) {
        case MODE_1_2_2_2_1:
        case MODE_7_1_FRONT_CENTER:
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                  DOWNMIX_SHIFT);
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
              (DOWNMIX_SHIFT - 1);
          break;
        default:
          break;
      }
      peak[0] = fixMax(peak[0], fixp_abs(tmp));
    }

    peak[0] = fMult(dmxGain5, peak[0])
              << (DFRACT_BITS - 1 - METADATA_FRACT_BITS);
  }

  if (drcComp->fullChannels > 2) {
    for (i = 0; i < drcComp->blockLength; i++) {
      const INT_PCM* pSamples = &inSamples[i * drcComp->channels];

      tmp = FL2FXCONST_DBL(0.f);
      if (drcComp->channelIdx[LS] >= 0)
        tmp -= fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[LS]]) >>
               (DOWNMIX_SHIFT - 1);
      if (drcComp->channelIdx[LS2] >= 0)
        tmp -= fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[LS2]]) >>
               (DOWNMIX_SHIFT - 1);
      if (drcComp->channelIdx[RS] >= 0)
        tmp -= fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[RS]]) >>
               (DOWNMIX_SHIFT - 1);
      if (drcComp->channelIdx[RS2] >= 0)
        tmp -= fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[RS2]]) >>
               (DOWNMIX_SHIFT - 1);
      if ((drcComp->channelIdx[LS] >= 0) && (drcComp->channelIdx[LS2] >= 0))
        tmp = fMult(FL2FXCONST_DBL(0.707f), tmp);
      if (drcComp->channelIdx[S] >= 0)
        tmp -= fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[S]]) >>
               (DOWNMIX_SHIFT - 1);
      if (drcComp->channelIdx[C] >= 0)
        tmp += fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[C]]) >>
               (DOWNMIX_SHIFT - 1);
      tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[L]]) >>
              DOWNMIX_SHIFT);

      if (dmxGain2 > FL2FXCONST_DBL(0.f)) {
        if (drcComp->fullChannels > 5) {
          tmp = fMult(dmxGain5, tmp) << (DFRACT_BITS - 1 - METADATA_FRACT_BITS);
        }
        tmp = fMult(dmxGain2, tmp) << (DFRACT_BITS - 1 - METADATA_FRACT_BITS);
      }
      peak[0] = fixMax(peak[0], fixp_abs(tmp));

      tmp = FL2FXCONST_DBL(0.f);
      if (drcComp->channelIdx[LS] >= 0)
        tmp += fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[LS]]) >>
               (DOWNMIX_SHIFT - 1);
      if (drcComp->channelIdx[LS2] >= 0)
        tmp += fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[LS2]]) >>
               (DOWNMIX_SHIFT - 1);
      if (drcComp->channelIdx[RS] >= 0)
        tmp += fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[RS]]) >>
               (DOWNMIX_SHIFT - 1);
      if (drcComp->channelIdx[RS2] >= 0)
        tmp += fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[RS2]]) >>
               (DOWNMIX_SHIFT - 1);
      if ((drcComp->channelIdx[RS] >= 0) && (drcComp->channelIdx[RS2] >= 0))
        tmp = fMult(FL2FXCONST_DBL(0.707f), tmp);
      if (drcComp->channelIdx[S] >= 0)
        tmp += fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[S]]) >>
               (DOWNMIX_SHIFT - 1);
      if (drcComp->channelIdx[C] >= 0)
        tmp += fMultDiv2(FL2FXCONST_DBL(0.707f),
                         (FIXP_PCM)pSamples[drcComp->channelIdx[C]]) >>
               (DOWNMIX_SHIFT - 1);
      tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[R]]) >>
              DOWNMIX_SHIFT);

      if (dmxGain2 > FL2FXCONST_DBL(0.f)) {
        if (drcComp->fullChannels > 5) {
          tmp = fMult(dmxGain5, tmp) << (DFRACT_BITS - 1 - METADATA_FRACT_BITS);
        }
        tmp = fMult(dmxGain2, tmp) << (DFRACT_BITS - 1 - METADATA_FRACT_BITS);
      }
      peak[0] = fixMax(peak[0], fixp_abs(tmp));
    }

    for (i = 0; i < drcComp->blockLength; i++) {
      const INT_PCM* pSamples = &inSamples[i * drcComp->channels];

      tmp = FL2FXCONST_DBL(0.f);
      switch (drcComp->chanConfig) {
        case MODE_6_1:
          tmp += fMultDiv2(fMult(slev, ext_leva),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[4]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(slev, ext_levb),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[6]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[0]]) >>
                  DOWNMIX_SHIFT);
          tmp +=
              fMultDiv2(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
              (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        case MODE_7_1_BACK:
        case MODE_7_1_REAR_SURROUND:
          tmp += fMultDiv2(fMult(slev, ext_leva),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[4]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(slev, ext_levb),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[0]]) >>
                  DOWNMIX_SHIFT);
          tmp +=
              fMultDiv2(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
              (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        case MODE_1_2_2_2_1:
        case MODE_7_1_FRONT_CENTER:
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[0]]) >>
                  DOWNMIX_SHIFT);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(ext_leva, clev),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[4]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
              (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        case MODE_7_1_TOP_FRONT:
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[0]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[4]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
              (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        default:
          if (drcComp->channelIdx[LS] >= 0)
            tmp +=
                fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[LS]]) >>
                (DOWNMIX_SHIFT - 1);
          if (drcComp->channelIdx[LS2] >= 0)
            tmp +=
                fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[LS2]]) >>
                (DOWNMIX_SHIFT - 1);
          if ((drcComp->channelIdx[LS] >= 0) && (drcComp->channelIdx[LS2] >= 0))
            tmp = fMult(FL2FXCONST_DBL(0.707f), tmp);
          if (drcComp->channelIdx[S] >= 0)
            tmp +=
                fMultDiv2(slev,
                          fMult(FL2FXCONST_DBL(0.7f),
                                (FIXP_PCM)pSamples[drcComp->channelIdx[S]])) >>
                (DOWNMIX_SHIFT - 1);
          if (drcComp->channelIdx[C] >= 0)
            tmp +=
                fMultDiv2(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[C]]) >>
                (DOWNMIX_SHIFT - 1);
          if (drcComp->channelIdx[3] >= 0)
            tmp += fMultDiv2(lfe_lev,
                             (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
                   (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[L]]) >>
                  DOWNMIX_SHIFT);
          break;
      }

      if (dmxGain2 > FL2FXCONST_DBL(0.f)) {
        if (drcComp->fullChannels > 5) {
          tmp = fMult(dmxGain5, tmp) << (DFRACT_BITS - 1 - METADATA_FRACT_BITS);
        }
        tmp = fMult(dmxGain2, tmp) << (DFRACT_BITS - 1 - METADATA_FRACT_BITS);
      }
      peak[0] = fixMax(peak[0], fixp_abs(tmp));

      tmp = FL2FXCONST_DBL(0.f);
      switch (drcComp->chanConfig) {
        case MODE_6_1:
          tmp += fMultDiv2(fMult(slev, ext_leva),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[5]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(slev, ext_levb),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[6]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[1]]) >>
                  DOWNMIX_SHIFT);
          tmp +=
              fMultDiv2(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
              (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        case MODE_7_1_BACK:
        case MODE_7_1_REAR_SURROUND:
          tmp += fMultDiv2(fMult(slev, ext_leva),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[5]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(slev, ext_levb),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[1]]) >>
                  DOWNMIX_SHIFT);
          tmp +=
              fMultDiv2(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
              (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        case MODE_1_2_2_2_1:
        case MODE_7_1_FRONT_CENTER:
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[1]]) >>
                  DOWNMIX_SHIFT);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(ext_leva, clev),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[5]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
              (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        case MODE_7_1_TOP_FRONT:
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[1]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[5]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
              (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        default:
          if (drcComp->channelIdx[RS] >= 0)
            tmp +=
                fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[RS]]) >>
                (DOWNMIX_SHIFT - 1);
          if (drcComp->channelIdx[RS2] >= 0)
            tmp +=
                fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[RS2]]) >>
                (DOWNMIX_SHIFT - 1);
          if ((drcComp->channelIdx[RS] >= 0) && (drcComp->channelIdx[RS2] >= 0))
            tmp = fMult(FL2FXCONST_DBL(0.707f), tmp);
          if (drcComp->channelIdx[S] >= 0)
            tmp +=
                fMultDiv2(slev,
                          fMult(FL2FXCONST_DBL(0.7f),
                                (FIXP_PCM)pSamples[drcComp->channelIdx[S]])) >>
                (DOWNMIX_SHIFT - 1);
          if (drcComp->channelIdx[C] >= 0)
            tmp +=
                fMultDiv2(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[C]]) >>
                (DOWNMIX_SHIFT - 1);
          if (drcComp->channelIdx[3] >= 0)
            tmp += fMultDiv2(lfe_lev,
                             (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
                   (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[R]]) >>
                  DOWNMIX_SHIFT);
      }

      if (dmxGain2 > FL2FXCONST_DBL(0.f)) {
        if (drcComp->fullChannels > 5) {
          tmp = fMult(dmxGain5, tmp) << (DFRACT_BITS - 1 - METADATA_FRACT_BITS);
        }
        tmp = fMult(dmxGain2, tmp) << (DFRACT_BITS - 1 - METADATA_FRACT_BITS);
      }
      peak[0] = fixMax(peak[0], fixp_abs(tmp));
    }
  }

  peak[1] = fixMax(peak[0], peak[1]);

  if (drcComp->fullChannels > 1) {
    for (i = 0; i < drcComp->blockLength; i++) {
      const INT_PCM* pSamples = &inSamples[i * drcComp->channels];

      tmp = FL2FXCONST_DBL(0.f);
      switch (drcComp->chanConfig) {
        case MODE_6_1:
          tmp += fMultDiv2(fMult(slev, ext_leva),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[4]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(slev, ext_leva),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[5]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMult(fMult(slev, ext_levb),
                       (FIXP_PCM)pSamples[drcComp->channelIdx[6]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMult(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[0]]) >>
                  DOWNMIX_SHIFT);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[1]]) >>
                  DOWNMIX_SHIFT);
          tmp += fMult(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
                 (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        case MODE_7_1_BACK:
        case MODE_7_1_REAR_SURROUND:
          tmp += fMultDiv2(fMult(slev, ext_leva),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[4]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(slev, ext_leva),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[5]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(slev, ext_levb),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(slev, ext_levb),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMult(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[0]]) >>
                  DOWNMIX_SHIFT);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[1]]) >>
                  DOWNMIX_SHIFT);
          tmp += fMult(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
                 (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        case MODE_1_2_2_2_1:
        case MODE_7_1_FRONT_CENTER:
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[0]]) >>
                  DOWNMIX_SHIFT);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[1]]) >>
                  DOWNMIX_SHIFT);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(ext_leva, clev),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(fMult(ext_leva, clev),
                           (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMult(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[4]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[5]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMult(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
                 (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        case MODE_7_1_TOP_FRONT:
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[0]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_leva, (FIXP_PCM)pSamples[drcComp->channelIdx[1]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[7]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp +=
              fMultDiv2(ext_levb, (FIXP_PCM)pSamples[drcComp->channelIdx[8]]) >>
              (DOWNMIX_SHIFT - 1);
          tmp += fMult(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[2]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[4]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[5]]) >>
                 (DOWNMIX_SHIFT - 1);
          tmp += fMult(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
                 (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          break;
        default:
          if (drcComp->channelIdx[LS] >= 0)
            tmp +=
                fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[LS]]) >>
                (DOWNMIX_SHIFT - 1);
          if (drcComp->channelIdx[LS2] >= 0)
            tmp +=
                fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[LS2]]) >>
                (DOWNMIX_SHIFT - 1);
          if (drcComp->channelIdx[RS] >= 0)
            tmp +=
                fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[RS]]) >>
                (DOWNMIX_SHIFT - 1);
          if (drcComp->channelIdx[RS2] >= 0)
            tmp +=
                fMultDiv2(slev, (FIXP_PCM)pSamples[drcComp->channelIdx[RS2]]) >>
                (DOWNMIX_SHIFT - 1);
          if ((drcComp->channelIdx[LS] >= 0) && (drcComp->channelIdx[LS2] >= 0))
            tmp = fMult(FL2FXCONST_DBL(0.707f), tmp);
          if (drcComp->channelIdx[S] >= 0)
            tmp +=
                fMultDiv2(slev,
                          fMult(FL2FXCONST_DBL(0.7f),
                                (FIXP_PCM)pSamples[drcComp->channelIdx[S]])) >>
                (DOWNMIX_SHIFT - 1);
          if (drcComp->channelIdx[C] >= 0)
            tmp += fMult(clev, (FIXP_PCM)pSamples[drcComp->channelIdx[C]]) >>
                   (DOWNMIX_SHIFT - 1);
          if (drcComp->channelIdx[3] >= 0)
            tmp += fMult(lfe_lev, (FIXP_PCM)pSamples[drcComp->channelIdx[3]]) >>
                   (DOWNMIX_SHIFT - 1 - LFE_LEV_SCALE);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[L]]) >>
                  DOWNMIX_SHIFT);
          tmp += (FX_PCM2FX_DBL((FIXP_PCM)pSamples[drcComp->channelIdx[R]]) >>
                  DOWNMIX_SHIFT);
      }

      if (dmxGain2 > FL2FXCONST_DBL(0.f)) {
        if (drcComp->fullChannels > 5) {
          tmp = fMult(dmxGain5, tmp) << (DFRACT_BITS - 1 - METADATA_FRACT_BITS);
        }
        tmp = fMult(dmxGain2, tmp) << (DFRACT_BITS - 1 - METADATA_FRACT_BITS);
      }
      peak[1] = fixMax(peak[1], fixp_abs(tmp));
    }
  }
}

INT FDK_DRC_Generator_Open(HDRC_COMP* phDrcComp) {
  INT err = 0;
  HDRC_COMP hDcComp = NULL;

  if (phDrcComp == NULL) {
    err = -1;
    goto bail;
  }

  hDcComp = (HDRC_COMP)FDKcalloc(1, sizeof(DRC_COMP));

  if (hDcComp == NULL) {
    err = -1;
    goto bail;
  }

  FDKmemclear(hDcComp, sizeof(DRC_COMP));

  *phDrcComp = hDcComp;
  return err;
bail:
  FDK_DRC_Generator_Close(&hDcComp);
  return err;
}

INT FDK_DRC_Generator_Close(HDRC_COMP* phDrcComp) {
  if (phDrcComp == NULL) {
    return -1;
  }
  if (*phDrcComp != NULL) {
    FDKfree(*phDrcComp);
    *phDrcComp = NULL;
  }
  return 0;
}

INT FDK_DRC_Generator_Initialize(HDRC_COMP drcComp,
                                 const DRC_PROFILE profileLine,
                                 const DRC_PROFILE profileRF,
                                 const INT blockLength, const UINT sampleRate,
                                 const CHANNEL_MODE channelMode,
                                 const CHANNEL_ORDER channelOrder,
                                 const UCHAR useWeighting) {
  int i;
  CHANNEL_MAPPING channelMapping;

  drcComp->limDecay =
      FL2FXCONST_DBL(((0.006f / 256) * blockLength) / METADATA_INT_SCALE);

  drcComp->blockLength = blockLength;
  drcComp->sampleRate = sampleRate;
  drcComp->chanConfig = channelMode;
  drcComp->useWeighting = useWeighting;

  if (FDK_DRC_Generator_setDrcProfile(drcComp, profileLine, profileRF) !=
      0) {
    return (-1);
  }

  if (FDKaacEnc_InitChannelMapping(channelMode, channelOrder,
                                   &channelMapping) != AAC_ENC_OK) {
    return (-2);
  }

  for (i = 0; i < 9; i++) drcComp->channelIdx[i] = -1;

  switch (channelMode) {
    case MODE_1:
      drcComp->channelIdx[C] = channelMapping.elInfo[0].ChannelIndex[0];
      break;
    case MODE_2:
      drcComp->channelIdx[L] = channelMapping.elInfo[0].ChannelIndex[0];
      drcComp->channelIdx[R] = channelMapping.elInfo[0].ChannelIndex[1];
      break;
    case MODE_1_2:
      drcComp->channelIdx[L] = channelMapping.elInfo[1].ChannelIndex[0];
      drcComp->channelIdx[R] = channelMapping.elInfo[1].ChannelIndex[1];
      drcComp->channelIdx[C] = channelMapping.elInfo[0].ChannelIndex[0];
      break;
    case MODE_1_2_1:
      drcComp->channelIdx[L] = channelMapping.elInfo[1].ChannelIndex[0];
      drcComp->channelIdx[R] = channelMapping.elInfo[1].ChannelIndex[1];
      drcComp->channelIdx[C] = channelMapping.elInfo[0].ChannelIndex[0];
      drcComp->channelIdx[S] = channelMapping.elInfo[2].ChannelIndex[0];
      break;
    case MODE_1_2_2:
      drcComp->channelIdx[L] = channelMapping.elInfo[1].ChannelIndex[0];
      drcComp->channelIdx[R] = channelMapping.elInfo[1].ChannelIndex[1];
      drcComp->channelIdx[C] = channelMapping.elInfo[0].ChannelIndex[0];
      drcComp->channelIdx[LS] = channelMapping.elInfo[2].ChannelIndex[0];
      drcComp->channelIdx[RS] = channelMapping.elInfo[2].ChannelIndex[1];
      break;
    case MODE_1_2_2_1:
      drcComp->channelIdx[L] = channelMapping.elInfo[1].ChannelIndex[0];
      drcComp->channelIdx[R] = channelMapping.elInfo[1].ChannelIndex[1];
      drcComp->channelIdx[C] = channelMapping.elInfo[0].ChannelIndex[0];
      drcComp->channelIdx[LFE] = channelMapping.elInfo[3].ChannelIndex[0];
      drcComp->channelIdx[LS] = channelMapping.elInfo[2].ChannelIndex[0];
      drcComp->channelIdx[RS] = channelMapping.elInfo[2].ChannelIndex[1];
      break;
    case MODE_1_2_2_2_1:
    case MODE_7_1_FRONT_CENTER:
      drcComp->channelIdx[L] = channelMapping.elInfo[2].ChannelIndex[0];
      drcComp->channelIdx[R] = channelMapping.elInfo[2].ChannelIndex[1];
      drcComp->channelIdx[C] = channelMapping.elInfo[0].ChannelIndex[0];
      drcComp->channelIdx[LFE] =
          channelMapping.elInfo[4].ChannelIndex[0];
      drcComp->channelIdx[LS] =
          channelMapping.elInfo[3].ChannelIndex[0];
      drcComp->channelIdx[RS] =
          channelMapping.elInfo[3].ChannelIndex[1];
      drcComp->channelIdx[LS2] =
          channelMapping.elInfo[1].ChannelIndex[0];
      drcComp->channelIdx[RS2] =
          channelMapping.elInfo[1].ChannelIndex[1];
      break;
    case MODE_7_1_BACK:
    case MODE_7_1_REAR_SURROUND:
      drcComp->channelIdx[L] = channelMapping.elInfo[1].ChannelIndex[0];
      drcComp->channelIdx[R] = channelMapping.elInfo[1].ChannelIndex[1];
      drcComp->channelIdx[C] = channelMapping.elInfo[0].ChannelIndex[0];
      drcComp->channelIdx[LFE] =
          channelMapping.elInfo[4].ChannelIndex[0];
      drcComp->channelIdx[LS] =
          channelMapping.elInfo[3].ChannelIndex[0];
      drcComp->channelIdx[RS] =
          channelMapping.elInfo[3].ChannelIndex[1];
      drcComp->channelIdx[LS2] =
          channelMapping.elInfo[2].ChannelIndex[0];
      drcComp->channelIdx[RS2] =
          channelMapping.elInfo[2].ChannelIndex[1];
      break;
    case MODE_6_1:
      drcComp->channelIdx[L] = channelMapping.elInfo[1].ChannelIndex[0];
      drcComp->channelIdx[R] = channelMapping.elInfo[1].ChannelIndex[1];
      drcComp->channelIdx[C] = channelMapping.elInfo[0].ChannelIndex[0];
      drcComp->channelIdx[LFE] =
          channelMapping.elInfo[4].ChannelIndex[0];
      drcComp->channelIdx[LS] =
          channelMapping.elInfo[2].ChannelIndex[0];
      drcComp->channelIdx[RS] =
          channelMapping.elInfo[2].ChannelIndex[1];
      drcComp->channelIdx[S] = channelMapping.elInfo[3].ChannelIndex[0];
      break;
    case MODE_7_1_TOP_FRONT:
      drcComp->channelIdx[L] = channelMapping.elInfo[1].ChannelIndex[0];
      drcComp->channelIdx[R] = channelMapping.elInfo[1].ChannelIndex[1];
      drcComp->channelIdx[C] = channelMapping.elInfo[0].ChannelIndex[0];
      drcComp->channelIdx[LFE] =
          channelMapping.elInfo[3].ChannelIndex[0];
      drcComp->channelIdx[LS] =
          channelMapping.elInfo[2].ChannelIndex[0];
      drcComp->channelIdx[RS] =
          channelMapping.elInfo[2].ChannelIndex[1];
      drcComp->channelIdx[LS2] =
          channelMapping.elInfo[4].ChannelIndex[0];
      drcComp->channelIdx[RS2] =
          channelMapping.elInfo[4].ChannelIndex[1];
      break;
    default:
      return (-1);
  }

  drcComp->fullChannels = channelMapping.nChannelsEff;
  drcComp->channels = channelMapping.nChannels;

  drcComp->smoothLevel[0] = drcComp->smoothLevel[1] =
      (FIXP_DBL)(-(135 << METADATA_FRACT_BITS));

  FDKmemclear(drcComp->smoothGain, sizeof(drcComp->smoothGain));
  FDKmemclear(drcComp->holdCnt, sizeof(drcComp->holdCnt));
  FDKmemclear(drcComp->limGain, sizeof(drcComp->limGain));
  FDKmemclear(drcComp->prevPeak, sizeof(drcComp->prevPeak));
  FDKmemclear(drcComp->filter, sizeof(drcComp->filter));

  return (0);
}

INT FDK_DRC_Generator_setDrcProfile(HDRC_COMP drcComp,
                                    const DRC_PROFILE profileLine,
                                    const DRC_PROFILE profileRF) {
  int profileIdx, i;

  drcComp->profile[0] = profileLine;
  drcComp->profile[1] = profileRF;

  for (i = 0; i < 2; i++) {
    switch (drcComp->profile[i]) {
      case DRC_NONE:
      case DRC_NOT_PRESENT:
      case DRC_FILMSTANDARD:
        profileIdx = 0;
        break;
      case DRC_FILMLIGHT:
        profileIdx = 1;
        break;
      case DRC_MUSICSTANDARD:
        profileIdx = 2;
        break;
      case DRC_MUSICLIGHT:
        profileIdx = 3;
        break;
      case DRC_SPEECH:
        profileIdx = 4;
        break;
      case DRC_DELAY_TEST:
        profileIdx = 5;
        break;
      default:
        return (-1);
    }

    if (profileIdx >= 0) {
      drcComp->maxBoostThr[i] = tabMaxBoostThr[profileIdx];
      drcComp->boostThr[i] = tabBoostThr[profileIdx];
      drcComp->earlyCutThr[i] = tabEarlyCutThr[profileIdx];
      drcComp->cutThr[i] = tabCutThr[profileIdx];
      drcComp->maxCutThr[i] = tabMaxCutThr[profileIdx];

      drcComp->boostFac[i] = tabBoostRatio[profileIdx];
      drcComp->earlyCutFac[i] = tabEarlyCutRatio[profileIdx];
      drcComp->cutFac[i] = tabCutRatio[profileIdx];

      drcComp->maxBoost[i] = tabMaxBoost[profileIdx];
      drcComp->maxCut[i] = tabMaxCut[profileIdx];
      drcComp->maxEarlyCut[i] =
          -fMult((drcComp->cutThr[i] - drcComp->earlyCutThr[i]),
                 drcComp->earlyCutFac[i]);

      drcComp->fastAttack[i] = tc2Coeff(
          tabFastAttack[profileIdx], drcComp->sampleRate, drcComp->blockLength);
      drcComp->fastDecay[i] = tc2Coeff(
          tabFastDecay[profileIdx], drcComp->sampleRate, drcComp->blockLength);
      drcComp->slowAttack[i] = tc2Coeff(
          tabSlowAttack[profileIdx], drcComp->sampleRate, drcComp->blockLength);
      drcComp->slowDecay[i] = tc2Coeff(
          tabSlowDecay[profileIdx], drcComp->sampleRate, drcComp->blockLength);
      drcComp->holdOff[i] = tabHoldOff[profileIdx] * 256 / drcComp->blockLength;

      drcComp->attackThr[i] = tabAttackThr[profileIdx];
      drcComp->decayThr[i] = tabDecayThr[profileIdx];
    }

    drcComp->smoothGain[i] = FL2FXCONST_DBL(0.f);
  }
  return (0);
}

INT FDK_DRC_Generator_Calc(HDRC_COMP drcComp, const INT_PCM* const inSamples,
                           const UINT inSamplesBufSize, const INT dialnorm,
                           const INT drc_TargetRefLevel,
                           const INT comp_TargetRefLevel, const FIXP_DBL clev,
                           const FIXP_DBL slev, const FIXP_DBL ext_leva,
                           const FIXP_DBL ext_levb, const FIXP_DBL lfe_lev,
                           const INT dmxGain5, const INT dmxGain2,
                           INT* const pDynrng, INT* const pCompr) {
  int i, c;
  FIXP_DBL peak[2];

  if ((drcComp->profile[0] != DRC_NONE) || (drcComp->profile[1] != DRC_NONE)) {
    FIXP_DBL level_b = FL2FXCONST_DBL(0.f);
    int level_e = DFRACT_BITS - 1;

    const int granuleLength = fixMin(16, drcComp->blockLength);

    if (drcComp->useWeighting) {
      FIXP_DBL x1, x2, y, y1, y2;
      const int granuleShift = getShiftFactor(granuleLength) - 1;

      for (c = 0; c < (int)drcComp->channels; c++) {
        const INT_PCM* pSamples = inSamples + c * inSamplesBufSize;

        if (c == drcComp->channelIdx[LFE]) {
          continue;
        }

        x1 = drcComp->filter[c].x1;
        x2 = drcComp->filter[c].x2;
        y1 = drcComp->filter[c].y1;
        y2 = drcComp->filter[c].y2;

        i = 0;

        do {
          int offset = i;
          FIXP_DBL accu = FL2FXCONST_DBL(0.f);

          for (i = offset;
               i < fixMin(offset + granuleLength, drcComp->blockLength); i++) {
            FIXP_DBL x =
                FX_PCM2FX_DBL((FIXP_PCM)pSamples[i]) >> WEIGHTING_FILTER_SHIFT;

            y = fMult(b0, x - x2) - fMult(a1, y1) - fMult(a2, y2);

            x2 = x1;
            x1 = x;
            y2 = y1;
            y1 = y;

            accu += fPow2Div2(y) >> (granuleShift - 1);
          }

          fixpAdd(accu, granuleShift + 2 * WEIGHTING_FILTER_SHIFT, &level_b,
                  &level_e);

        } while (i < drcComp->blockLength);

        drcComp->filter[c].x1 = x1;
        drcComp->filter[c].x2 = x2;
        drcComp->filter[c].y1 = y1;
        drcComp->filter[c].y2 = y2;
      }
    }
    else {
      const int granuleShift = getShiftFactor(granuleLength);

      for (c = 0; c < (int)drcComp->channels; c++) {
        const INT_PCM* pSamples = inSamples + c * inSamplesBufSize;

        if ((int)c == drcComp->channelIdx[LFE]) {
          continue;
        }

        i = 0;

        do {
          int offset = i;
          FIXP_DBL accu = FL2FXCONST_DBL(0.f);

          for (i = offset;
               i < fixMin(offset + granuleLength, drcComp->blockLength); i++) {
            accu += fPow2Div2((FIXP_PCM)pSamples[i]) >> (granuleShift - 1);
          }

          fixpAdd(accu, granuleShift, &level_b,
                  &level_e);

        } while (i < drcComp->blockLength);
      }
    }

    FIXP_DBL ldLevel =
        CalcLdData(level_b) +
        (FIXP_DBL)((level_e - 12) << (DFRACT_BITS - 1 - LD_DATA_SHIFT)) -
        CalcLdData((FIXP_DBL)(drcComp->blockLength << (DFRACT_BITS - 1 - 12)));

    ldLevel =
        fMax(ldLevel, FL2FXCONST_DBL(-0.51905126482615036685473741085772f));

    FIXP_DBL level = fMult(
        (FIXP_DBL)(10 << (METADATA_FRACT_BITS + LD_DATA_SHIFT)),
        fMult(FL2FXCONST_DBL(0.30102999566398119521373889472449f), ldLevel) +
            (FIXP_DBL)(FL2FXCONST_DBL(0.3f) >> LD_DATA_SHIFT));

    level -= ((FIXP_DBL)(dialnorm << (METADATA_FRACT_BITS - 16)) +
              (FIXP_DBL)(31 << METADATA_FRACT_BITS));

    for (i = 0; i < 2; i++) {
      if (drcComp->profile[i] == DRC_NONE) {
        drcComp->smoothGain[i] = FL2FXCONST_DBL(0.f);
      } else {
        FIXP_DBL gain, alpha, lvl2smthlvl;

        if (level <= drcComp->maxBoostThr[i]) {
          gain = drcComp->maxBoost[i];
        } else if (level < drcComp->boostThr[i]) {
          gain = fMult((level - drcComp->boostThr[i]), drcComp->boostFac[i]);
        } else if (level <= drcComp->earlyCutThr[i]) {
          gain = FL2FXCONST_DBL(0.f);
        } else if (level <= drcComp->cutThr[i]) {
          gain =
              fMult((level - drcComp->earlyCutThr[i]), drcComp->earlyCutFac[i]);
        } else if (level < drcComp->maxCutThr[i]) {
          gain = fMult((level - drcComp->cutThr[i]), drcComp->cutFac[i]) -
                 drcComp->maxEarlyCut[i];
        } else {
          gain = -drcComp->maxCut[i];
        }

        lvl2smthlvl = level - drcComp->smoothLevel[i];
        if (gain < drcComp->smoothGain[i]) {
          if (lvl2smthlvl > drcComp->attackThr[i]) {
            alpha = drcComp->fastAttack[i];
          } else {
            alpha = drcComp->slowAttack[i];
          }
        } else {
          if (lvl2smthlvl < -drcComp->decayThr[i]) {
            alpha = drcComp->fastDecay[i];
          } else {
            alpha = drcComp->slowDecay[i];
          }
        }

        if ((gain < drcComp->smoothGain[i]) ||
            (drcComp->holdCnt[i] ==
             0)) {
          FIXP_DBL accu;

          accu = fMult(((FIXP_DBL)MAXVAL_DBL - alpha), drcComp->smoothLevel[i]);
          accu += fMult(alpha, level);
          drcComp->smoothLevel[i] = accu;

          accu = fMult(((FIXP_DBL)MAXVAL_DBL - alpha), drcComp->smoothGain[i]);
          accu += fMult(alpha, gain);
          drcComp->smoothGain[i] = accu;
        }

        if (drcComp->holdCnt[i]) {
          drcComp->holdCnt[i]--;
        }
        if (gain < drcComp->smoothGain[i]) {
          drcComp->holdCnt[i] = drcComp->holdOff[i];
        }
      }
    }
  } else {
    drcComp->smoothGain[0] = FL2FXCONST_DBL(0.f);
    drcComp->smoothGain[1] = FL2FXCONST_DBL(0.f);
  }

  findPeakLevels(drcComp, inSamples, clev, slev, ext_leva, ext_levb, lfe_lev,
                 (FIXP_DBL)((LONG)(dmxGain5) << (METADATA_FRACT_BITS - 16)),
                 (FIXP_DBL)((LONG)(dmxGain2) << (METADATA_FRACT_BITS - 16)),
                 peak);

  for (i = 0; i < 2; i++) {
    FIXP_DBL tmp = drcComp->prevPeak[i];
    drcComp->prevPeak[i] = peak[i];
    peak[i] = fixMax(peak[i], tmp);

    FIXP_DBL ld_peak =
        CalcLdData(peak[i]) +
        (FIXP_DBL)((LONG)DOWNMIX_SHIFT << (DFRACT_BITS - 1 - LD_DATA_SHIFT));

    ld_peak =
        fMax(ld_peak, FL2FXCONST_DBL(-0.31143075889569022011284244651463f));

    peak[i] = fMult(
        (FIXP_DBL)(10 << (METADATA_FRACT_BITS + LD_DATA_SHIFT)),
        fMult(FL2FX_DBL(2 * 0.30102999566398119521373889472449f), ld_peak));
    peak[i] +=
        (FL2FX_DBL(0.5f) >> METADATA_INT_BITS);
    peak[i] += drcComp->smoothGain[i];
  }

  peak[0] -= (FIXP_DBL)((dialnorm - drc_TargetRefLevel)
                        << (METADATA_FRACT_BITS -
                            16));

  peak[1] -=
      (FIXP_DBL)((dialnorm - comp_TargetRefLevel)
                 << (METADATA_FRACT_BITS -
                     16));

  drcComp->limGain[0] += drcComp->limDecay;
  drcComp->limGain[0] = fixMin(drcComp->limGain[0], -peak[0]);

  drcComp->limGain[1] += 2 * drcComp->limDecay;
  drcComp->limGain[1] = fixMin(drcComp->limGain[1], -peak[1]);

  {
    FIXP_DBL tmp;

    tmp = drcComp->smoothGain[0];
    if (drcComp->limGain[0] < FL2FXCONST_DBL(0.f)) {
      tmp += drcComp->limGain[0];
    }
    *pDynrng = (LONG)scaleValue(tmp, -(METADATA_FRACT_BITS - 16));

    tmp = drcComp->smoothGain[1];
    if (drcComp->limGain[1] < FL2FXCONST_DBL(0.f)) {
      tmp += drcComp->limGain[1];
    }
    *pCompr = (LONG)scaleValue(tmp, -(METADATA_FRACT_BITS - 16));
  }

  return 0;
}

DRC_PROFILE FDK_DRC_Generator_getDrcProfile(const HDRC_COMP drcComp) {
  return drcComp->profile[0];
}

DRC_PROFILE FDK_DRC_Generator_getCompProfile(const HDRC_COMP drcComp) {
  return drcComp->profile[1];
}
