/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Psychoaccoustic configuration

*******************************************************************************/

#include "psy_configuration.h"
#include "adj_thr.h"
#include "aacEnc_rom.h"
#include "genericStds.h"
#include "FDK_trigFcts.h"

typedef struct {
  LONG sampleRate;
  const SFB_PARAM_LONG *paramLong;
  const SFB_PARAM_SHORT *paramShort;
} SFB_INFO_TAB;

static const SFB_INFO_TAB sfbInfoTab[] = {
    {8000, &p_FDKaacEnc_8000_long_1024, &p_FDKaacEnc_8000_short_128},
    {11025, &p_FDKaacEnc_11025_long_1024, &p_FDKaacEnc_11025_short_128},
    {12000, &p_FDKaacEnc_12000_long_1024, &p_FDKaacEnc_12000_short_128},
    {16000, &p_FDKaacEnc_16000_long_1024, &p_FDKaacEnc_16000_short_128},
    {22050, &p_FDKaacEnc_22050_long_1024, &p_FDKaacEnc_22050_short_128},
    {24000, &p_FDKaacEnc_24000_long_1024, &p_FDKaacEnc_24000_short_128},
    {32000, &p_FDKaacEnc_32000_long_1024, &p_FDKaacEnc_32000_short_128},
    {44100, &p_FDKaacEnc_44100_long_1024, &p_FDKaacEnc_44100_short_128},
    {48000, &p_FDKaacEnc_48000_long_1024, &p_FDKaacEnc_48000_short_128},
    {64000, &p_FDKaacEnc_64000_long_1024, &p_FDKaacEnc_64000_short_128},
    {88200, &p_FDKaacEnc_88200_long_1024, &p_FDKaacEnc_88200_short_128},
    {96000, &p_FDKaacEnc_96000_long_1024, &p_FDKaacEnc_96000_short_128}

};

static const SFB_PARAM_LONG p_22050_long_512 = {
    31, {4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  8,  8,  8,  12, 12,
         12, 16, 20, 24, 28, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32}};

static const SFB_PARAM_LONG p_32000_long_512 = {
    37,
    {4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  8,  8,  8,  8, 8,
     12, 12, 12, 12, 16, 16, 16, 20, 24, 24, 28, 32, 32, 32, 32, 32, 32, 32}};

static const SFB_PARAM_LONG p_44100_long_512 = {
    36, {4, 4, 4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  8,  8,  8,
         8, 8, 12, 12, 12, 12, 16, 20, 24, 28, 32, 32, 32, 32, 32, 32, 32, 52}};

static const SFB_INFO_TAB sfbInfoTabLD512[] = {
    {8000, &p_22050_long_512, NULL},   {11025, &p_22050_long_512, NULL},
    {12000, &p_22050_long_512, NULL},  {16000, &p_22050_long_512, NULL},
    {22050, &p_22050_long_512, NULL},  {24000, &p_22050_long_512, NULL},
    {32000, &p_32000_long_512, NULL},  {44100, &p_44100_long_512, NULL},
    {48000, &p_44100_long_512, NULL},  {64000, &p_44100_long_512, NULL},
    {88200, &p_44100_long_512, NULL},  {96000, &p_44100_long_512, NULL},
    {128000, &p_44100_long_512, NULL}, {176400, &p_44100_long_512, NULL},
    {192000, &p_44100_long_512, NULL}, {256000, &p_44100_long_512, NULL},
    {352800, &p_44100_long_512, NULL}, {384000, &p_44100_long_512, NULL},
};

static const SFB_PARAM_LONG p_22050_long_480 = {
    30, {4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  8,  8,  8,  12,
         12, 12, 16, 20, 24, 28, 32, 32, 32, 32, 32, 32, 32, 32, 32}};

static const SFB_PARAM_LONG p_32000_long_480 = {
    37, {4, 4, 4, 4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  8,  8, 8,
         8, 8, 8, 12, 12, 12, 16, 16, 20, 24, 32, 32, 32, 32, 32, 32, 32, 32}};

static const SFB_PARAM_LONG p_44100_long_480 = {
    35, {4, 4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  8,  8,  8, 8,
         8, 12, 12, 12, 12, 12, 16, 16, 24, 28, 32, 32, 32, 32, 32, 32, 48}};

static const SFB_INFO_TAB sfbInfoTabLD480[] = {
    {8000, &p_22050_long_480, NULL},   {11025, &p_22050_long_480, NULL},
    {12000, &p_22050_long_480, NULL},  {16000, &p_22050_long_480, NULL},
    {22050, &p_22050_long_480, NULL},  {24000, &p_22050_long_480, NULL},
    {32000, &p_32000_long_480, NULL},  {44100, &p_44100_long_480, NULL},
    {48000, &p_44100_long_480, NULL},  {64000, &p_44100_long_480, NULL},
    {88200, &p_44100_long_480, NULL},  {96000, &p_44100_long_480, NULL},
    {128000, &p_44100_long_480, NULL}, {176400, &p_44100_long_480, NULL},
    {192000, &p_44100_long_480, NULL}, {256000, &p_44100_long_480, NULL},
    {352800, &p_44100_long_480, NULL}, {384000, &p_44100_long_480, NULL},
};

#define Q_BARCVAL (25)

AAC_ENCODER_ERROR FDKaacEnc_initSfbTable(const LONG sampleRate,
                                         const INT blockType,
                                         const INT granuleLength,
                                         INT *const sfbOffset,
                                         INT *const sfbCnt) {
  INT i, specStartOffset = 0;
  INT granuleLengthWindow = granuleLength;
  const UCHAR *sfbWidth = NULL;
  const SFB_INFO_TAB *sfbInfo = NULL;
  int size;

  switch (granuleLength) {
    case 1024:
    case 960:
      sfbInfo = sfbInfoTab;
      size = (INT)(sizeof(sfbInfoTab) / sizeof(SFB_INFO_TAB));
      break;
    case 512:
      sfbInfo = sfbInfoTabLD512;
      size = sizeof(sfbInfoTabLD512);
      break;
    case 480:
      sfbInfo = sfbInfoTabLD480;
      size = sizeof(sfbInfoTabLD480);
      break;
    default:
      return AAC_ENC_INVALID_FRAME_LENGTH;
  }

  for (i = 0; i < size; i++) {
    if (sfbInfo[i].sampleRate == sampleRate) {
      switch (blockType) {
        case LONG_WINDOW:
        case START_WINDOW:
        case STOP_WINDOW:
          sfbWidth = sfbInfo[i].paramLong->sfbWidth;
          *sfbCnt = sfbInfo[i].paramLong->sfbCnt;
          break;
        case SHORT_WINDOW:
          sfbWidth = sfbInfo[i].paramShort->sfbWidth;
          *sfbCnt = sfbInfo[i].paramShort->sfbCnt;
          granuleLengthWindow /= TRANS_FAC;
          break;
      }
      break;
    }
  }
  if (i == size) {
    return AAC_ENC_UNSUPPORTED_SAMPLINGRATE;
  }

  for (i = 0; i < *sfbCnt; i++) {
    sfbOffset[i] = specStartOffset;
    specStartOffset += sfbWidth[i];
    if (specStartOffset >= granuleLengthWindow) {
      i++;
      break;
    }
  }
  *sfbCnt = fixMin(i, *sfbCnt);
  sfbOffset[*sfbCnt] = fixMin(specStartOffset, granuleLengthWindow);
  return AAC_ENC_OK;
}

static FIXP_DBL FDKaacEnc_BarcLineValue(INT noOfLines, INT fftLine,
                                        LONG samplingFreq) {
  FIXP_DBL FOURBY3EM4 = (FIXP_DBL)0x45e7b273;
  FIXP_DBL PZZZ76 = (FIXP_DBL)0x639d5e4a;
  FIXP_DBL ONE3P3 = (FIXP_DBL)0x35333333;
  FIXP_DBL THREEP5 = (FIXP_DBL)0x1c000000;
  FIXP_DBL INV480 = (FIXP_DBL)0x44444444;

  FIXP_DBL center_freq, x1, x2;
  FIXP_DBL bvalFFTLine, atan1, atan2;

  center_freq = fftLine * samplingFreq;

  switch (noOfLines) {
    case 1024:
      center_freq = center_freq << 2;
      break;
    case 128:
      center_freq = center_freq << 5;
      break;
    case 512:
      center_freq = (fftLine * samplingFreq) << 3;
      break;
    case 480:
      center_freq = fMult(center_freq, INV480) << 4;
      break;
    default:
      center_freq = (FIXP_DBL)0;
  }

  x1 = fMult(center_freq, FOURBY3EM4);
  x2 = fMult(center_freq, PZZZ76)
       << 2;

  atan1 = fixp_atan(x1);
  atan2 = fixp_atan(x2);

  bvalFFTLine = fMult(ONE3P3, atan2) + fMult(THREEP5, fMult(atan1, atan1));
  return (bvalFFTLine);
}

static void FDKaacEnc_InitMinPCMResolution(int numPb, int *pbOffset,
                                           FIXP_DBL *sfbPCMquantThreshold) {
#define PCM_QUANT_NOISE ((FIXP_DBL)0x00547062)

  for (int i = 0; i < numPb; i++) {
    sfbPCMquantThreshold[i] = (pbOffset[i + 1] - pbOffset[i]) * PCM_QUANT_NOISE;
  }
}

static FIXP_DBL getMaskFactor(const FIXP_DBL dbVal_fix, const INT dbVal_e,
                              const FIXP_DBL ten_fix, const INT ten_e) {
  INT q_msk;
  FIXP_DBL mask_factor;

  mask_factor = fPow(ten_fix, DFRACT_BITS - 1 - ten_e, -dbVal_fix,
                     DFRACT_BITS - 1 - dbVal_e, &q_msk);
  q_msk = fixMin(DFRACT_BITS - 1, fixMax(-(DFRACT_BITS - 1), q_msk));

  if ((q_msk > 0) && (mask_factor > (FIXP_DBL)MAXVAL_DBL >> q_msk)) {
    mask_factor = (FIXP_DBL)MAXVAL_DBL;
  } else {
    mask_factor = scaleValue(mask_factor, q_msk);
  }

  return (mask_factor);
}

static void FDKaacEnc_initSpreading(INT numPb, FIXP_DBL *pbBarcValue,
                                    FIXP_DBL *pbMaskLoFactor,
                                    FIXP_DBL *pbMaskHiFactor,
                                    FIXP_DBL *pbMaskLoFactorSprEn,
                                    FIXP_DBL *pbMaskHiFactorSprEn,
                                    const LONG bitrate, const INT blockType)

{
  INT i;
  FIXP_DBL MASKLOWSPREN, MASKHIGHSPREN;

  FIXP_DBL MASKHIGH = (FIXP_DBL)0x30000000;
  FIXP_DBL MASKLOW = (FIXP_DBL)0x60000000;
  FIXP_DBL MASKLOWSPRENLONG = (FIXP_DBL)0x60000000;
  FIXP_DBL MASKHIGHSPRENLONG = (FIXP_DBL)0x40000000;
  FIXP_DBL MASKHIGHSPRENLONGLOWBR = (FIXP_DBL)0x30000000;
  FIXP_DBL MASKLOWSPRENSHORT = (FIXP_DBL)0x40000000;
  FIXP_DBL MASKHIGHSPRENSHORT = (FIXP_DBL)0x30000000;
  FIXP_DBL TEN = (FIXP_DBL)0x50000000;

  if (blockType != SHORT_WINDOW) {
    MASKLOWSPREN = MASKLOWSPRENLONG;
    MASKHIGHSPREN =
        (bitrate > 20000) ? MASKHIGHSPRENLONG : MASKHIGHSPRENLONGLOWBR;
  } else {
    MASKLOWSPREN = MASKLOWSPRENSHORT;
    MASKHIGHSPREN = MASKHIGHSPRENSHORT;
  }

  for (i = 0; i < numPb; i++) {
    if (i > 0) {
      pbMaskHiFactor[i] = getMaskFactor(
          fMult(MASKHIGH, (pbBarcValue[i] - pbBarcValue[i - 1])), 23, TEN, 27);

      pbMaskLoFactor[i - 1] = getMaskFactor(
          fMult(MASKLOW, (pbBarcValue[i] - pbBarcValue[i - 1])), 23, TEN, 27);

      pbMaskHiFactorSprEn[i] = getMaskFactor(
          fMult(MASKHIGHSPREN, (pbBarcValue[i] - pbBarcValue[i - 1])), 23, TEN,
          27);

      pbMaskLoFactorSprEn[i - 1] = getMaskFactor(
          fMult(MASKLOWSPREN, (pbBarcValue[i] - pbBarcValue[i - 1])), 23, TEN,
          27);
    } else {
      pbMaskHiFactor[i] = (FIXP_DBL)0;
      pbMaskLoFactor[numPb - 1] = (FIXP_DBL)0;
      pbMaskHiFactorSprEn[i] = (FIXP_DBL)0;
      pbMaskLoFactorSprEn[numPb - 1] = (FIXP_DBL)0;
    }
  }
}

static void FDKaacEnc_initBarcValues(INT numPb, INT *pbOffset, INT numLines,
                                     INT samplingFrequency, FIXP_DBL *pbBval) {
  INT i;
  FIXP_DBL MAX_BARC = (FIXP_DBL)0x30000000;

  for (i = 0; i < numPb; i++) {
    FIXP_DBL v1, v2, cur_bark;
    v1 = FDKaacEnc_BarcLineValue(numLines, pbOffset[i], samplingFrequency);
    v2 = FDKaacEnc_BarcLineValue(numLines, pbOffset[i + 1], samplingFrequency);
    cur_bark = (v1 >> 1) + (v2 >> 1);
    pbBval[i] = fixMin(cur_bark, MAX_BARC);
  }
}

static void FDKaacEnc_initMinSnr(const LONG bitrate, const LONG samplerate,
                                 const INT numLines, const INT *sfbOffset,
                                 const INT sfbActive, const INT blockType,
                                 FIXP_DBL *sfbMinSnrLdData) {
  INT sfb;

  INT qbfac, qperwin, qdiv, qpeprt_const, qpeprt;
  INT qtmp, qsnr, sfbWidth;

  FIXP_DBL MAX_BARC = (FIXP_DBL)0x30000000;
  FIXP_DBL MAX_BARCP1 = (FIXP_DBL)0x32000000;
  FIXP_DBL BITS2PEFAC = (FIXP_DBL)0x4b851eb8;
  FIXP_DBL PERS2P4 = (FIXP_DBL)0x624dd2f2;
  FIXP_DBL ONEP5 = (FIXP_DBL)0x60000000;
  FIXP_DBL MAX_SNR = (FIXP_DBL)0x33333333;
  FIXP_DBL MIN_SNR = (FIXP_DBL)0x003126e9;

  FIXP_DBL barcFactor, pePerWindow, pePart, barcWidth;
  FIXP_DBL pePart_const, tmp, snr, one_qsnr, one_point5;

  barcFactor = fDivNorm(fixMin(FDKaacEnc_BarcLineValue(
                                   numLines, sfbOffset[sfbActive], samplerate),
                               MAX_BARC),
                        MAX_BARCP1, &qbfac);

  qbfac = DFRACT_BITS - 1 - qbfac;

  pePerWindow = fDivNorm(bitrate, samplerate, &qperwin);
  qperwin = DFRACT_BITS - 1 - qperwin;
  pePerWindow = fMult(pePerWindow, BITS2PEFAC);
  qperwin = qperwin + 30 - (DFRACT_BITS - 1);
  pePerWindow = fMult(pePerWindow, PERS2P4);
  qperwin = qperwin + 36 - (DFRACT_BITS - 1);

  switch (numLines) {
    case 1024:
      qperwin = qperwin - 10;
      break;
    case 128:
      qperwin = qperwin - 7;
      break;
    case 512:
      qperwin = qperwin - 9;
      break;
    case 480:
      qperwin = qperwin - 9;
      pePerWindow = fMult(pePerWindow, FL2FXCONST_DBL(480.f / 512.f));
      break;
  }

  if (blockType == SHORT_WINDOW) {
    pePerWindow = fMult(pePerWindow, ONEP5);
    qperwin = qperwin + 30 - (DFRACT_BITS - 1);
  }
  pePart_const = fDivNorm(pePerWindow, barcFactor, &qdiv);
  qpeprt_const = qperwin - qbfac + DFRACT_BITS - 1 - qdiv;

  for (sfb = 0; sfb < sfbActive; sfb++) {
    barcWidth =
        FDKaacEnc_BarcLineValue(numLines, sfbOffset[sfb + 1], samplerate) -
        FDKaacEnc_BarcLineValue(numLines, sfbOffset[sfb], samplerate);

    pePart = fMult(pePart_const, barcWidth);
    qpeprt = qpeprt_const + 25 - (DFRACT_BITS - 1);

    sfbWidth = (sfbOffset[sfb + 1] - sfbOffset[sfb]);
    pePart = fDivNorm(pePart, sfbWidth, &qdiv);
    qpeprt += DFRACT_BITS - 1 - qdiv;

    tmp = f2Pow(pePart, DFRACT_BITS - 1 - qpeprt, &qtmp);
    qtmp = DFRACT_BITS - 1 - qtmp;

    qsnr = fixMin(qtmp, 30);
    tmp = tmp >> (qtmp - qsnr);

    if ((30 + 1 - qsnr) > (DFRACT_BITS - 1))
      one_point5 = (FIXP_DBL)0;
    else
      one_point5 = (FIXP_DBL)(ONEP5 >> (30 + 1 - qsnr));

    snr = (tmp >> 1) - (one_point5);
    qsnr -= 1;

    if (qsnr > 0)
      one_qsnr = (FIXP_DBL)(1 << qsnr);
    else
      one_qsnr = (FIXP_DBL)0;

    snr = fixMax(one_qsnr, snr);

    snr = fDivNorm(one_qsnr, snr, &qsnr);
    qsnr = DFRACT_BITS - 1 - qsnr;
    snr = (qsnr > 30) ? (snr >> (qsnr - 30)) : snr;

    snr = (snr > MAX_SNR) ? MAX_SNR : snr;

    snr = (snr < MIN_SNR) ? MIN_SNR : snr;
    snr = snr << 1;

    sfbMinSnrLdData[sfb] = CalcLdData(snr);
  }
}

AAC_ENCODER_ERROR FDKaacEnc_InitPsyConfiguration(INT bitrate, INT samplerate,
                                                 INT bandwidth, INT blocktype,
                                                 INT granuleLength, INT useIS,
                                                 INT useMS,
                                                 PSY_CONFIGURATION *psyConf,
                                                 FB_TYPE filterbank) {
  AAC_ENCODER_ERROR ErrorStatus;
  INT sfb;
  FIXP_DBL sfbBarcVal[MAX_SFB];
  const INT frameLengthLong = granuleLength;
  const INT frameLengthShort = granuleLength / TRANS_FAC;
  INT downscaleFactor = 1;

  switch (granuleLength) {
    case 256:
    case 240:
      downscaleFactor = 2;
      break;
    case 128:
    case 120:
      downscaleFactor = 4;
      break;
    default:
      downscaleFactor = 1;
      break;
  }

  FDKmemclear(psyConf, sizeof(PSY_CONFIGURATION));
  psyConf->granuleLength = granuleLength;
  psyConf->filterbank = filterbank;

  psyConf->allowIS = (useIS) && ((bitrate / bandwidth) < 5);
  psyConf->allowMS = useMS;

  ErrorStatus = FDKaacEnc_initSfbTable(samplerate * downscaleFactor, blocktype,
                                       granuleLength * downscaleFactor,
                                       psyConf->sfbOffset, &psyConf->sfbCnt);

  if (ErrorStatus != AAC_ENC_OK) return ErrorStatus;

  FDKaacEnc_initBarcValues(psyConf->sfbCnt, psyConf->sfbOffset,
                           psyConf->sfbOffset[psyConf->sfbCnt], samplerate,
                           sfbBarcVal);

  FDKaacEnc_InitMinPCMResolution(psyConf->sfbCnt, psyConf->sfbOffset,
                                 psyConf->sfbPcmQuantThreshold);

  FDKaacEnc_initSpreading(psyConf->sfbCnt, sfbBarcVal,
                          psyConf->sfbMaskLowFactor, psyConf->sfbMaskHighFactor,
                          psyConf->sfbMaskLowFactorSprEn,
                          psyConf->sfbMaskHighFactorSprEn, bitrate, blocktype);

  psyConf->maxAllowedIncreaseFactor = 2;
  psyConf->minRemainingThresholdFactor = (FIXP_SGL)0x0148;

  psyConf->clipEnergy =
      (FIXP_DBL)0x773593ff;

  if (blocktype != SHORT_WINDOW) {
    psyConf->lowpassLine =
        (INT)((2 * bandwidth * frameLengthLong) / samplerate);
    psyConf->lowpassLineLFE = LFE_LOWPASS_LINE;
  } else {
    psyConf->lowpassLine =
        (INT)((2 * bandwidth * frameLengthShort) / samplerate);
    psyConf->lowpassLineLFE = 0;
    psyConf->clipEnergy >>= 6;
  }

  for (sfb = 0; sfb < psyConf->sfbCnt; sfb++) {
    if (psyConf->sfbOffset[sfb] >= psyConf->lowpassLine) break;
  }
  psyConf->sfbActive = fMax(sfb, 1);

  for (sfb = 0; sfb < psyConf->sfbCnt; sfb++) {
    if (psyConf->sfbOffset[sfb] >= psyConf->lowpassLineLFE) break;
  }
  psyConf->sfbActiveLFE = sfb;
  psyConf->sfbActive = fMax(psyConf->sfbActive, psyConf->sfbActiveLFE);

  FDKaacEnc_initMinSnr(bitrate, samplerate * downscaleFactor,
                       psyConf->sfbOffset[psyConf->sfbCnt], psyConf->sfbOffset,
                       psyConf->sfbActive, blocktype, psyConf->sfbMinSnrLdData);

  return AAC_ENC_OK;
}
