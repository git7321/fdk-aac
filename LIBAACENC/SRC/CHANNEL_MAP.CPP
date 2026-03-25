/**************************** AAC encoder library ******************************

   Author(s):   A. Groeschel

   Description: channel mapping functionality

*******************************************************************************/

#include "channel_map.h"
#include "bitenc.h"
#include "psy_const.h"
#include "qc_data.h"
#include "aacEnc_ram.h"
#include "FDK_tools_rom.h"

static const CHANNEL_MODE_CONFIG_TAB channelModeConfig[] = {
    {MODE_1, 1, 1, 1},
    {MODE_2, 2, 2, 1},
    {MODE_1_2, 3, 3, 2},
    {MODE_1_2_1, 4, 4, 3},
    {MODE_1_2_2, 5, 5, 3},
    {MODE_1_2_2_1, 6, 5, 4},
    {MODE_1_2_2_2_1, 8, 7, 5},
    {MODE_6_1, 7, 6, 5},
    {MODE_7_1_BACK, 8, 7, 5},
    {MODE_7_1_TOP_FRONT, 8, 7, 5},
    {MODE_7_1_REAR_SURROUND, 8, 7,
     5},
    {MODE_7_1_FRONT_CENTER, 8, 7,
     5},
};

AAC_ENCODER_ERROR FDKaacEnc_DetermineEncoderMode(CHANNEL_MODE* mode,
                                                 INT nChannels) {
  INT i;
  CHANNEL_MODE encMode = MODE_INVALID;

  if (*mode == MODE_UNKNOWN) {
    for (i = 0; i < (INT)sizeof(channelModeConfig) /
                        (INT)sizeof(CHANNEL_MODE_CONFIG_TAB);
         i++) {
      if (channelModeConfig[i].nChannels == nChannels) {
        encMode = channelModeConfig[i].encMode;
        break;
      }
    }
    *mode = encMode;
  } else {
    if (FDKaacEnc_GetChannelModeConfiguration(*mode)->nChannels == nChannels) {
      encMode = *mode;
    }
  }

  if (encMode == MODE_INVALID) {
    return AAC_ENC_UNSUPPORTED_CHANNELCONFIG;
  }

  return AAC_ENC_OK;
}

static INT FDKaacEnc_initElement(ELEMENT_INFO* elInfo, MP4_ELEMENT_ID elType,
                                 INT* cnt, FDK_channelMapDescr* mapDescr,
                                 UINT mapIdx, INT* it_cnt,
                                 const FIXP_DBL relBits) {
  INT error = 0;
  INT counter = *cnt;

  elInfo->elType = elType;
  elInfo->relativeBits = relBits;

  switch (elInfo->elType) {
    case ID_SCE:
    case ID_LFE:
    case ID_CCE:
      elInfo->nChannelsInEl = 1;
      elInfo->ChannelIndex[0] =
          FDK_chMapDescr_getMapValue(mapDescr, counter++, mapIdx);
      elInfo->instanceTag = it_cnt[elType]++;
      break;
    case ID_CPE:
      elInfo->nChannelsInEl = 2;
      elInfo->ChannelIndex[0] =
          FDK_chMapDescr_getMapValue(mapDescr, counter++, mapIdx);
      elInfo->ChannelIndex[1] =
          FDK_chMapDescr_getMapValue(mapDescr, counter++, mapIdx);
      elInfo->instanceTag = it_cnt[elType]++;
      break;
    case ID_DSE:
      elInfo->nChannelsInEl = 0;
      elInfo->ChannelIndex[0] = 0;
      elInfo->ChannelIndex[1] = 0;
      elInfo->instanceTag = it_cnt[elType]++;
      break;
    default:
      error = 1;
  };
  *cnt = counter;
  return error;
}

AAC_ENCODER_ERROR FDKaacEnc_InitChannelMapping(CHANNEL_MODE mode,
                                               CHANNEL_ORDER co,
                                               CHANNEL_MAPPING* cm) {
  INT count = 0;
  INT it_cnt[ID_END + 1];
  INT i;
  UINT mapIdx;
  FDK_channelMapDescr mapDescr;

  for (i = 0; i < ID_END; i++) it_cnt[i] = 0;

  FDKmemclear(cm, sizeof(CHANNEL_MAPPING));

  for (i = 0; i < (INT)sizeof(channelModeConfig) /
                      (INT)sizeof(CHANNEL_MODE_CONFIG_TAB);
       i++) {
    if (channelModeConfig[i].encMode == mode) {
      cm->encMode = channelModeConfig[i].encMode;
      cm->nChannels = channelModeConfig[i].nChannels;
      cm->nChannelsEff = channelModeConfig[i].nChannelsEff;
      cm->nElements = channelModeConfig[i].nElements;

      break;
    }
  }

  FDK_chMapDescr_init(&mapDescr, NULL, 0, (co == CH_ORDER_MPEG) ? 1 : 0);
  switch (mode) {
    case MODE_7_1_REAR_SURROUND:
      mapIdx = (INT)MODE_7_1_BACK;
      break;
    case MODE_7_1_FRONT_CENTER:
      mapIdx = (INT)MODE_1_2_2_2_1;
      break;
    default:
      mapIdx =
          (INT)mode > 14
              ? 0
              : (INT)
                    mode;
  }

  switch (mode) {
    case MODE_1:
      FDKaacEnc_initElement(&cm->elInfo[0], ID_SCE, &count, &mapDescr, mapIdx,
                            it_cnt, (FIXP_DBL)MAXVAL_DBL);
      break;
    case MODE_2:
      FDKaacEnc_initElement(&cm->elInfo[0], ID_CPE, &count, &mapDescr, mapIdx,
                            it_cnt, (FIXP_DBL)MAXVAL_DBL);
      break;

    case MODE_1_2:
      FDKaacEnc_initElement(&cm->elInfo[0], ID_SCE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.4f));
      FDKaacEnc_initElement(&cm->elInfo[1], ID_CPE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.6f));
      break;

    case MODE_1_2_1:
      FDKaacEnc_initElement(&cm->elInfo[0], ID_SCE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.3f));
      FDKaacEnc_initElement(&cm->elInfo[1], ID_CPE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.4f));
      FDKaacEnc_initElement(&cm->elInfo[2], ID_SCE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.3f));
      break;

    case MODE_1_2_2:
      FDKaacEnc_initElement(&cm->elInfo[0], ID_SCE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.26f));
      FDKaacEnc_initElement(&cm->elInfo[1], ID_CPE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.37f));
      FDKaacEnc_initElement(&cm->elInfo[2], ID_CPE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.37f));
      break;

    case MODE_1_2_2_1:
      FDKaacEnc_initElement(&cm->elInfo[0], ID_SCE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.24f));
      FDKaacEnc_initElement(&cm->elInfo[1], ID_CPE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.35f));
      FDKaacEnc_initElement(&cm->elInfo[2], ID_CPE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.35f));
      FDKaacEnc_initElement(&cm->elInfo[3], ID_LFE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.06f));
      break;

    case MODE_6_1:
      FDKaacEnc_initElement(&cm->elInfo[0], ID_SCE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.2f));
      FDKaacEnc_initElement(&cm->elInfo[1], ID_CPE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.275f));
      FDKaacEnc_initElement(&cm->elInfo[2], ID_CPE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.275f));
      FDKaacEnc_initElement(&cm->elInfo[3], ID_SCE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.2f));
      FDKaacEnc_initElement(&cm->elInfo[4], ID_LFE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.05f));
      break;

    case MODE_1_2_2_2_1:
    case MODE_7_1_BACK:
    case MODE_7_1_TOP_FRONT:
    case MODE_7_1_REAR_SURROUND:
    case MODE_7_1_FRONT_CENTER: {

      FDKaacEnc_initElement(&cm->elInfo[0], ID_SCE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.18f));
      FDKaacEnc_initElement(&cm->elInfo[1], ID_CPE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.26f));
      FDKaacEnc_initElement(&cm->elInfo[2], ID_CPE, &count, &mapDescr, mapIdx,
                            it_cnt, FL2FXCONST_DBL(0.26f));
      if (mode != MODE_7_1_TOP_FRONT) {
        FDKaacEnc_initElement(&cm->elInfo[3], ID_CPE, &count, &mapDescr, mapIdx,
                              it_cnt, FL2FXCONST_DBL(0.26f));
        FDKaacEnc_initElement(&cm->elInfo[4], ID_LFE, &count, &mapDescr, mapIdx,
                              it_cnt, FL2FXCONST_DBL(0.04f));
      } else {
        FDKaacEnc_initElement(&cm->elInfo[3], ID_LFE, &count, &mapDescr, mapIdx,
                              it_cnt, FL2FXCONST_DBL(0.04f));
        FDKaacEnc_initElement(&cm->elInfo[4], ID_CPE, &count, &mapDescr, mapIdx,
                              it_cnt, FL2FXCONST_DBL(0.26f));
      }
      break;
    }

    default:
      return AAC_ENC_UNSUPPORTED_CHANNELCONFIG;
  };

  FDK_ASSERT(cm->nElements <= ((8)));

  return AAC_ENC_OK;
}

AAC_ENCODER_ERROR FDKaacEnc_InitElementBits(QC_STATE* hQC, CHANNEL_MAPPING* cm,
                                            INT bitrateTot, INT averageBitsTot,
                                            INT maxChannelBits) {
  int sc_brTot = CountLeadingBits(bitrateTot);

  switch (cm->encMode) {
    case MODE_1:
      hQC->elementBits[0]->chBitrateEl = bitrateTot;

      hQC->elementBits[0]->maxBitsEl = maxChannelBits;

      hQC->elementBits[0]->relativeBitsEl = cm->elInfo[0].relativeBits;
      break;

    case MODE_2:
      hQC->elementBits[0]->chBitrateEl = bitrateTot >> 1;

      hQC->elementBits[0]->maxBitsEl = 2 * maxChannelBits;

      hQC->elementBits[0]->relativeBitsEl = cm->elInfo[0].relativeBits;
      break;
    case MODE_1_2: {
      hQC->elementBits[0]->relativeBitsEl = cm->elInfo[0].relativeBits;
      hQC->elementBits[1]->relativeBitsEl = cm->elInfo[1].relativeBits;
      FIXP_DBL sceRate = cm->elInfo[0].relativeBits;
      FIXP_DBL cpeRate = cm->elInfo[1].relativeBits;

      hQC->elementBits[0]->chBitrateEl =
          fMult(sceRate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> sc_brTot;
      hQC->elementBits[1]->chBitrateEl =
          fMult(cpeRate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);

      hQC->elementBits[0]->maxBitsEl = maxChannelBits;
      hQC->elementBits[1]->maxBitsEl = 2 * maxChannelBits;
      break;
    }
    case MODE_1_2_1: {
      hQC->elementBits[0]->relativeBitsEl = cm->elInfo[0].relativeBits;
      hQC->elementBits[1]->relativeBitsEl = cm->elInfo[1].relativeBits;
      hQC->elementBits[2]->relativeBitsEl = cm->elInfo[2].relativeBits;
      FIXP_DBL sce1Rate = cm->elInfo[0].relativeBits;
      FIXP_DBL cpeRate = cm->elInfo[1].relativeBits;
      FIXP_DBL sce2Rate = cm->elInfo[2].relativeBits;

      hQC->elementBits[0]->chBitrateEl =
          fMult(sce1Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> sc_brTot;
      hQC->elementBits[1]->chBitrateEl =
          fMult(cpeRate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);
      hQC->elementBits[2]->chBitrateEl =
          fMult(sce2Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> sc_brTot;

      hQC->elementBits[0]->maxBitsEl = maxChannelBits;
      hQC->elementBits[1]->maxBitsEl = 2 * maxChannelBits;
      hQC->elementBits[2]->maxBitsEl = maxChannelBits;
      break;
    }
    case MODE_1_2_2: {
      hQC->elementBits[0]->relativeBitsEl = cm->elInfo[0].relativeBits;
      hQC->elementBits[1]->relativeBitsEl = cm->elInfo[1].relativeBits;
      hQC->elementBits[2]->relativeBitsEl = cm->elInfo[2].relativeBits;
      FIXP_DBL sceRate = cm->elInfo[0].relativeBits;
      FIXP_DBL cpe1Rate = cm->elInfo[1].relativeBits;
      FIXP_DBL cpe2Rate = cm->elInfo[2].relativeBits;

      hQC->elementBits[0]->chBitrateEl =
          fMult(sceRate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> sc_brTot;
      hQC->elementBits[1]->chBitrateEl =
          fMult(cpe1Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);
      hQC->elementBits[2]->chBitrateEl =
          fMult(cpe2Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);

      hQC->elementBits[0]->maxBitsEl = maxChannelBits;
      hQC->elementBits[1]->maxBitsEl = 2 * maxChannelBits;
      hQC->elementBits[2]->maxBitsEl = 2 * maxChannelBits;
      break;
    }
    case MODE_1_2_2_1: {
      hQC->elementBits[0]->relativeBitsEl = cm->elInfo[0].relativeBits;
      hQC->elementBits[1]->relativeBitsEl = cm->elInfo[1].relativeBits;
      hQC->elementBits[2]->relativeBitsEl = cm->elInfo[2].relativeBits;
      hQC->elementBits[3]->relativeBitsEl = cm->elInfo[3].relativeBits;
      FIXP_DBL sceRate = cm->elInfo[0].relativeBits;
      FIXP_DBL cpe1Rate = cm->elInfo[1].relativeBits;
      FIXP_DBL cpe2Rate = cm->elInfo[2].relativeBits;
      FIXP_DBL lfeRate = cm->elInfo[3].relativeBits;

      int maxBitsTot =
          maxChannelBits * 5;
      int sc = CountLeadingBits(fixMax(maxChannelBits, averageBitsTot));
      int maxLfeBits = (int)fMax(
          (INT)((fMult(lfeRate, (FIXP_DBL)(maxChannelBits << sc)) >> sc) << 1),
          (INT)((fMult(FL2FXCONST_DBL(1.1f / 2.f),
                       fMult(lfeRate, (FIXP_DBL)(averageBitsTot << sc)))
                 << 1) >>
                sc));

      maxChannelBits = (maxBitsTot - maxLfeBits);
      sc = CountLeadingBits(maxChannelBits);

      maxChannelBits =
          fMult((FIXP_DBL)maxChannelBits << sc, GetInvInt(5)) >> sc;

      hQC->elementBits[0]->chBitrateEl =
          fMult(sceRate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> sc_brTot;
      hQC->elementBits[1]->chBitrateEl =
          fMult(cpe1Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);
      hQC->elementBits[2]->chBitrateEl =
          fMult(cpe2Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);
      hQC->elementBits[3]->chBitrateEl =
          fMult(lfeRate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> sc_brTot;

      hQC->elementBits[0]->maxBitsEl = maxChannelBits;
      hQC->elementBits[1]->maxBitsEl = 2 * maxChannelBits;
      hQC->elementBits[2]->maxBitsEl = 2 * maxChannelBits;
      hQC->elementBits[3]->maxBitsEl = maxLfeBits;

      break;
    }
    case MODE_6_1: {
      FIXP_DBL sceRate = hQC->elementBits[0]->relativeBitsEl =
          cm->elInfo[0].relativeBits;
      FIXP_DBL cpe1Rate = hQC->elementBits[1]->relativeBitsEl =
          cm->elInfo[1].relativeBits;
      FIXP_DBL cpe2Rate = hQC->elementBits[2]->relativeBitsEl =
          cm->elInfo[2].relativeBits;
      FIXP_DBL sce2Rate = hQC->elementBits[3]->relativeBitsEl =
          cm->elInfo[3].relativeBits;
      FIXP_DBL lfeRate = hQC->elementBits[4]->relativeBitsEl =
          cm->elInfo[4].relativeBits;

      int maxBitsTot =
          maxChannelBits * 6;
      int sc = CountLeadingBits(fixMax(maxChannelBits, averageBitsTot));
      int maxLfeBits = (int)fMax(
          (INT)((fMult(lfeRate, (FIXP_DBL)(maxChannelBits << sc)) >> sc) << 1),
          (INT)((fMult(FL2FXCONST_DBL(1.1f / 2.f),
                       fMult(lfeRate, (FIXP_DBL)(averageBitsTot << sc)))
                 << 1) >>
                sc));

      maxChannelBits = (maxBitsTot - maxLfeBits) / 6;

      hQC->elementBits[0]->chBitrateEl =
          fMult(sceRate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> sc_brTot;
      hQC->elementBits[1]->chBitrateEl =
          fMult(cpe1Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);
      hQC->elementBits[2]->chBitrateEl =
          fMult(cpe2Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);
      hQC->elementBits[3]->chBitrateEl =
          fMult(sce2Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);
      hQC->elementBits[4]->chBitrateEl =
          fMult(lfeRate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> sc_brTot;

      hQC->elementBits[0]->maxBitsEl = maxChannelBits;
      hQC->elementBits[1]->maxBitsEl = 2 * maxChannelBits;
      hQC->elementBits[2]->maxBitsEl = 2 * maxChannelBits;
      hQC->elementBits[3]->maxBitsEl = maxChannelBits;
      hQC->elementBits[4]->maxBitsEl = maxLfeBits;
      break;
    }
    case MODE_7_1_TOP_FRONT:
    case MODE_7_1_BACK:
    case MODE_7_1_REAR_SURROUND:
    case MODE_7_1_FRONT_CENTER:
    case MODE_1_2_2_2_1: {
      int cpe3Idx = (cm->encMode != MODE_7_1_TOP_FRONT) ? 3 : 4;
      int lfeIdx = (cm->encMode != MODE_7_1_TOP_FRONT) ? 4 : 3;

      FIXP_DBL sceRate = hQC->elementBits[0]->relativeBitsEl =
          cm->elInfo[0].relativeBits;
      FIXP_DBL cpe1Rate = hQC->elementBits[1]->relativeBitsEl =
          cm->elInfo[1].relativeBits;
      FIXP_DBL cpe2Rate = hQC->elementBits[2]->relativeBitsEl =
          cm->elInfo[2].relativeBits;
      FIXP_DBL cpe3Rate = hQC->elementBits[cpe3Idx]->relativeBitsEl =
          cm->elInfo[cpe3Idx].relativeBits;
      FIXP_DBL lfeRate = hQC->elementBits[lfeIdx]->relativeBitsEl =
          cm->elInfo[lfeIdx].relativeBits;

      int maxBitsTot =
          maxChannelBits * 7;
      int sc = CountLeadingBits(fixMax(maxChannelBits, averageBitsTot));
      int maxLfeBits = (int)fMax(
          (INT)((fMult(lfeRate, (FIXP_DBL)(maxChannelBits << sc)) >> sc) << 1),
          (INT)((fMult(FL2FXCONST_DBL(1.1f / 2.f),
                       fMult(lfeRate, (FIXP_DBL)(averageBitsTot << sc)))
                 << 1) >>
                sc));

      maxChannelBits = (maxBitsTot - maxLfeBits) / 7;

      hQC->elementBits[0]->chBitrateEl =
          fMult(sceRate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> sc_brTot;
      hQC->elementBits[1]->chBitrateEl =
          fMult(cpe1Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);
      hQC->elementBits[2]->chBitrateEl =
          fMult(cpe2Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);
      hQC->elementBits[cpe3Idx]->chBitrateEl =
          fMult(cpe3Rate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> (sc_brTot + 1);
      hQC->elementBits[lfeIdx]->chBitrateEl =
          fMult(lfeRate, (FIXP_DBL)(bitrateTot << sc_brTot)) >> sc_brTot;

      hQC->elementBits[0]->maxBitsEl = maxChannelBits;
      hQC->elementBits[1]->maxBitsEl = 2 * maxChannelBits;
      hQC->elementBits[2]->maxBitsEl = 2 * maxChannelBits;
      hQC->elementBits[cpe3Idx]->maxBitsEl = 2 * maxChannelBits;
      hQC->elementBits[lfeIdx]->maxBitsEl = maxLfeBits;
      break;
    }

    default:
      return AAC_ENC_UNSUPPORTED_CHANNELCONFIG;
  }

  return AAC_ENC_OK;
}

ELEMENT_MODE FDKaacEnc_GetMonoStereoMode(const CHANNEL_MODE mode) {
  ELEMENT_MODE monoStereoSetting = EL_MODE_INVALID;

  switch (mode) {
    case MODE_1:
      monoStereoSetting = EL_MODE_MONO;
      break;

    case MODE_2:
    case MODE_1_2:
    case MODE_1_2_1:
    case MODE_1_2_2:
    case MODE_1_2_2_1:
    case MODE_6_1:
    case MODE_1_2_2_2_1:
    case MODE_7_1_REAR_SURROUND:
    case MODE_7_1_FRONT_CENTER:
    case MODE_7_1_BACK:
    case MODE_7_1_TOP_FRONT:
      monoStereoSetting = EL_MODE_STEREO;
      break;

    default:
      monoStereoSetting = EL_MODE_INVALID;
      break;
  }

  return monoStereoSetting;
}

const CHANNEL_MODE_CONFIG_TAB* FDKaacEnc_GetChannelModeConfiguration(
    const CHANNEL_MODE mode) {
  INT i;
  const CHANNEL_MODE_CONFIG_TAB* cm_config = NULL;

  for (i = 0; i < (INT)sizeof(channelModeConfig) /
                      (INT)sizeof(CHANNEL_MODE_CONFIG_TAB);
       i++) {
    if (channelModeConfig[i].encMode == mode) {
      cm_config = &channelModeConfig[i];
      break;
    }
  }
  return cm_config;
}
