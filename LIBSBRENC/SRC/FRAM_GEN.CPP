/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#include "fram_gen.h"
#include "sbr_misc.h"
#include "genericStds.h"

static const SBR_FRAME_INFO frameInfo1_2048 = {1, {0, 16}, {FREQ_RES_HIGH},
                                               0, 1,       {0, 16}};

static const SBR_FRAME_INFO frameInfo2_2048 = {
    2, {0, 8, 16}, {FREQ_RES_HIGH, FREQ_RES_HIGH}, 0, 2, {0, 8, 16}};

static const SBR_FRAME_INFO frameInfo4_2048 = {
    4,
    {0, 4, 8, 12, 16},
    {FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH},
    0,
    2,
    {0, 8, 16}};

static const SBR_FRAME_INFO frameInfo1_2304 = {1, {0, 18}, {FREQ_RES_HIGH},
                                               0, 1,       {0, 18}};

static const SBR_FRAME_INFO frameInfo2_2304 = {
    2, {0, 9, 18}, {FREQ_RES_HIGH, FREQ_RES_HIGH}, 0, 2, {0, 9, 18}};

static const SBR_FRAME_INFO frameInfo4_2304 = {
    4,
    {0, 5, 9, 14, 18},
    {FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH},
    0,
    2,
    {0, 9, 18}};

static const SBR_FRAME_INFO frameInfo1_1920 = {1, {0, 15}, {FREQ_RES_HIGH},
                                               0, 1,       {0, 15}};

static const SBR_FRAME_INFO frameInfo2_1920 = {
    2, {0, 8, 15}, {FREQ_RES_HIGH, FREQ_RES_HIGH}, 0, 2, {0, 8, 15}};

static const SBR_FRAME_INFO frameInfo4_1920 = {
    4,
    {0, 4, 8, 12, 15},
    {FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH},
    0,
    2,
    {0, 8, 15}};

static const SBR_FRAME_INFO frameInfo1_1152 = {1, {0, 9}, {FREQ_RES_HIGH},
                                               0, 1,      {0, 9}};

static const SBR_FRAME_INFO frameInfo2_1152 = {
    2, {0, 5, 9}, {FREQ_RES_HIGH, FREQ_RES_HIGH}, 0, 2, {0, 5, 9}};

static const SBR_FRAME_INFO frameInfo4_1152 = {
    4,
    {0, 2, 5, 7, 9},
    {FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH},
    0,
    2,
    {0, 5, 9}};

static const SBR_FRAME_INFO frameInfo1_512LD = {1, {0, 8}, {FREQ_RES_HIGH},
                                                0, 1,      {0, 8}};

static const SBR_FRAME_INFO frameInfo2_512LD = {
    2, {0, 4, 8}, {FREQ_RES_HIGH, FREQ_RES_HIGH}, 0, 2, {0, 4, 8}};

static const SBR_FRAME_INFO frameInfo4_512LD = {
    4,
    {0, 2, 4, 6, 8},
    {FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH},
    0,
    2,
    {0, 4, 8}};

static int calcFillLengthMax(
    int tranPos,
    int numberTimeSlots
);

static void fillFrameTran(
    const int *v_tuningSegm,
    const int *v_tuningFreq,
    int tran,
    int *v_bord,
    int *length_v_bord,
    int *v_freq,
    int *length_v_freq,
    int *bmin,
    int *bmax
);

static void fillFramePre(INT dmax, INT *v_bord, INT *length_v_bord, INT *v_freq,
                         INT *length_v_freq, INT bmin, INT rest);

static void fillFramePost(INT *parts, INT *d, INT dmax, INT *v_bord,
                          INT *length_v_bord, INT *v_freq, INT *length_v_freq,
                          INT bmax, INT bufferFrameStart, INT numberTimeSlots,
                          INT fmax);

static void fillFrameInter(INT *nL, const int *v_tuningSegm, INT *v_bord,
                           INT *length_v_bord, INT bmin, INT *v_freq,
                           INT *length_v_freq, INT *v_bordFollow,
                           INT *length_v_bordFollow, INT *v_freqFollow,
                           INT *length_v_freqFollow, INT i_fillFollow, INT dmin,
                           INT dmax, INT numberTimeSlots);

static void calcFrameClass(FRAME_CLASS *frameClass, FRAME_CLASS *frameClassOld,
                           INT tranFlag, INT *spreadFlag);

static void specialCase(INT *spreadFlag, INT allowSpread, INT *v_bord,
                        INT *length_v_bord, INT *v_freq, INT *length_v_freq,
                        INT *parts, INT d);

static void calcCmonBorder(INT *i_cmon, INT *i_tran, INT *v_bord,
                           INT *length_v_bord, INT tran, INT bufferFrameStart,
                           INT numberTimeSlots);

static void keepForFollowUp(INT *v_bordFollow, INT *length_v_bordFollow,
                            INT *v_freqFollow, INT *length_v_freqFollow,
                            INT *i_tranFollow, INT *i_fillFollow, INT *v_bord,
                            INT *length_v_bord, INT *v_freq, INT i_cmon,
                            INT i_tran, INT parts, INT numberTimeSlots);

static void calcCtrlSignal(HANDLE_SBR_GRID hSbrGrid, FRAME_CLASS frameClass,
                           INT *v_bord, INT length_v_bord, INT *v_freq,
                           INT length_v_freq, INT i_cmon, INT i_tran,
                           INT spreadFlag, INT nL);

static void ctrlSignal2FrameInfo(HANDLE_SBR_GRID hSbrGrid,
                                 HANDLE_SBR_FRAME_INFO hFrameInfo,
                                 FREQ_RES *freq_res_fixfix);

static const int envelopeTable_8[8][5] = {
    /*[|T-|------]*/ {2, 0, 0, 1, -1},
    /*[|-T-|-----]*/ {2, 0, 0, 2, -1},
    /*[--|T-|----]*/ {3, 1, 1, 2, 4},
    /*[---|T-|---]*/ {3, 1, 1, 3, 5},
    /*[----|T-|--]*/ {3, 1, 1, 4, 6},
    /*[-----|T--|]*/ {2, 1, 1, 5, -1},
    /*[------|T-|]*/ {2, 1, 1, 6, -1},
    /*[-------|T|]*/ {2, 1, 1, 7, -1},
};

static const int envelopeTable_16[16][6] = {
    /*[|T---|------------|]*/ {2, 0, 0, 4, -1, -1},
    /*[|-T---|-----------|]*/ {2, 0, 0, 5, -1, -1},
    /*[|--|T---|----------]*/ {3, 1, 1, 2, 6, -1},
    /*[|---|T---|---------]*/ {3, 1, 1, 3, 7, -1},
    /*[|----|T---|--------]*/ {3, 1, 1, 4, 8, -1},
    /*[|-----|T---|-------]*/ {3, 1, 1, 5, 9, -1},
    /*[|------|T---|------]*/ {3, 1, 1, 6, 10, -1},
    /*[|-------|T---|-----]*/ {3, 1, 1, 7, 11, -1},
    /*[|--------|T---|----]*/ {3, 1, 1, 8, 12, -1},
    /*[|---------|T---|---]*/ {3, 1, 1, 9, 13, -1},
    /*[|----------|T---|--]*/ {3, 1, 1, 10, 14, -1},
    /*[|-----------|T----|]*/ {2, 1, 1, 11, -1, -1},
    /*[|------------|T---|]*/ {2, 1, 1, 12, -1, -1},
    /*[|-------------|T--|]*/ {2, 1, 1, 13, -1, -1},
    /*[|--------------|T-|]*/ {2, 1, 1, 14, -1, -1},
    /*[|---------------|T|]*/ {2, 1, 1, 15, -1, -1},
};

static const int envelopeTable_15[15][6] = {
    /*[|T---|------------]*/ {2, 0, 0, 4, -1, -1},
    /*[|-T---|-----------]*/ {2, 0, 0, 5, -1, -1},
    /*[|--|T---|---------]*/ {3, 1, 1, 2, 6, -1},
    /*[|---|T---|--------]*/ {3, 1, 1, 3, 7, -1},
    /*[|----|T---|-------]*/ {3, 1, 1, 4, 8, -1},
    /*[|-----|T---|------]*/ {3, 1, 1, 5, 9, -1},
    /*[|------|T---|-----]*/ {3, 1, 1, 6, 10, -1},
    /*[|-------|T---|----]*/ {3, 1, 1, 7, 11, -1},
    /*[|--------|T---|---]*/ {3, 1, 1, 8, 12, -1},
    /*[|---------|T---|--]*/ {3, 1, 1, 9, 13, -1},
    /*[|----------|T----|]*/ {2, 1, 1, 10, -1, -1},
    /*[|-----------|T---|]*/ {2, 1, 1, 11, -1, -1},
    /*[|------------|T--|]*/ {2, 1, 1, 12, -1, -1},
    /*[|-------------|T-|]*/ {2, 1, 1, 13, -1, -1},
    /*[|--------------|T|]*/ {2, 1, 1, 14, -1, -1},
};

static const int minFrameTranDistance = 4;

static const FREQ_RES freqRes_table_8[] = {
    FREQ_RES_LOW,  FREQ_RES_LOW,  FREQ_RES_LOW,  FREQ_RES_LOW, FREQ_RES_LOW,
    FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH, FREQ_RES_HIGH};

static const FREQ_RES freqRes_table_16[16] = {
    /* 0-4 */ FREQ_RES_LOW,
    FREQ_RES_LOW,
    FREQ_RES_LOW,
    FREQ_RES_LOW,
    FREQ_RES_LOW,
    /* 5-9 */ FREQ_RES_LOW,
    FREQ_RES_HIGH,
    FREQ_RES_HIGH,
    FREQ_RES_HIGH,
    FREQ_RES_HIGH,
    /* 10-16 */ FREQ_RES_HIGH,
    FREQ_RES_HIGH,
    FREQ_RES_HIGH,
    FREQ_RES_HIGH,
    FREQ_RES_HIGH,
    FREQ_RES_HIGH};

static void generateFixFixOnly(HANDLE_SBR_FRAME_INFO hSbrFrameInfo,
                               HANDLE_SBR_GRID hSbrGrid, int tranPosInternal,
                               int numberTimeSlots, UCHAR fResTransIsLow);

HANDLE_SBR_FRAME_INFO
FDKsbrEnc_frameInfoGenerator(HANDLE_SBR_ENVELOPE_FRAME hSbrEnvFrame,
                             UCHAR *v_transient_info, const INT rightBorderFIX,
                             UCHAR *v_transient_info_pre, int ldGrid,
                             const int *v_tuning) {
  INT numEnv, tranPosInternal = 0, bmin = 0, bmax = 0, parts, d, i_cmon = 0,
              i_tran = 0, nL;
  INT fmax = 0;

  INT *v_bord = hSbrEnvFrame->v_bord;
  INT *v_freq = hSbrEnvFrame->v_freq;
  INT *v_bordFollow = hSbrEnvFrame->v_bordFollow;
  INT *v_freqFollow = hSbrEnvFrame->v_freqFollow;

  INT *length_v_bordFollow = &hSbrEnvFrame->length_v_bordFollow;
  INT *length_v_freqFollow = &hSbrEnvFrame->length_v_freqFollow;
  INT *length_v_bord = &hSbrEnvFrame->length_v_bord;
  INT *length_v_freq = &hSbrEnvFrame->length_v_freq;
  INT *spreadFlag = &hSbrEnvFrame->spreadFlag;
  INT *i_tranFollow = &hSbrEnvFrame->i_tranFollow;
  INT *i_fillFollow = &hSbrEnvFrame->i_fillFollow;
  FRAME_CLASS *frameClassOld = &hSbrEnvFrame->frameClassOld;
  FRAME_CLASS frameClass = FIXFIX;

  INT allowSpread = hSbrEnvFrame->allowSpread;
  INT numEnvStatic = hSbrEnvFrame->numEnvStatic;
  INT staticFraming = hSbrEnvFrame->staticFraming;
  INT dmin = hSbrEnvFrame->dmin;
  INT dmax = hSbrEnvFrame->dmax;

  INT bufferFrameStart = hSbrEnvFrame->SbrGrid.bufferFrameStart;
  INT numberTimeSlots = hSbrEnvFrame->SbrGrid.numberTimeSlots;
  INT frameMiddleSlot = hSbrEnvFrame->frameMiddleSlot;

  INT tranPos = v_transient_info[0];
  INT tranFlag = v_transient_info[1];

  const int *v_tuningSegm = v_tuning;
  const int *v_tuningFreq = v_tuning + 3;

  hSbrEnvFrame->v_tuningSegm = v_tuningSegm;

  if (ldGrid) {
    if (!tranFlag && v_transient_info_pre[1] &&
        (numberTimeSlots - v_transient_info_pre[0] < minFrameTranDistance)) {
      tranFlag = 1;
      tranPos = 0;
    }
  }

  if (staticFraming) {
    frameClass = FIXFIX;
    numEnv = numEnvStatic;
    *frameClassOld = FIXFIX;
    hSbrEnvFrame->SbrGrid.bs_num_env = numEnv;
    hSbrEnvFrame->SbrGrid.frameClass = frameClass;
  } else {
    if (rightBorderFIX) {
      tranFlag = 0;
      *spreadFlag = 0;
    }
    calcFrameClass(&frameClass, frameClassOld, tranFlag, spreadFlag);

    if (tranFlag && ldGrid) {
      frameClass = FIXFIXonly;
      *frameClassOld = FIXFIX;
    }

    if (tranFlag) {
      tranPosInternal =
          frameMiddleSlot + tranPos + bufferFrameStart;

      fillFrameTran(v_tuningSegm, v_tuningFreq, tranPosInternal, v_bord,
                    length_v_bord, v_freq, length_v_freq, &bmin, &bmax);

      fmax = calcFillLengthMax(tranPos, numberTimeSlots);
    }

    switch (frameClass) {
      case FIXFIXonly:
        FDK_ASSERT(ldGrid);
        tranPosInternal = tranPos;
        generateFixFixOnly(&(hSbrEnvFrame->SbrFrameInfo),
                           &(hSbrEnvFrame->SbrGrid), tranPosInternal,
                           numberTimeSlots, hSbrEnvFrame->fResTransIsLow);

        return &(hSbrEnvFrame->SbrFrameInfo);

      case FIXVAR:

        fillFramePre(dmax, v_bord, length_v_bord, v_freq, length_v_freq, bmin,
                     bmin - bufferFrameStart);

        fillFramePost(&parts, &d, dmax, v_bord, length_v_bord, v_freq,
                      length_v_freq, bmax, bufferFrameStart, numberTimeSlots,
                      fmax);

        if (parts == 1 && d < dmin)
          specialCase(spreadFlag, allowSpread, v_bord, length_v_bord, v_freq,
                      length_v_freq, &parts, d);

        calcCmonBorder(&i_cmon, &i_tran, v_bord, length_v_bord, tranPosInternal,
                       bufferFrameStart, numberTimeSlots);

        keepForFollowUp(v_bordFollow, length_v_bordFollow, v_freqFollow,
                        length_v_freqFollow, i_tranFollow, i_fillFollow, v_bord,
                        length_v_bord, v_freq, i_cmon, i_tran, parts,
                        numberTimeSlots);

        calcCtrlSignal(&hSbrEnvFrame->SbrGrid, frameClass, v_bord,
                       *length_v_bord, v_freq, *length_v_freq, i_cmon, i_tran,
                       *spreadFlag, DC);
        break;
      case VARFIX:
        calcCtrlSignal(&hSbrEnvFrame->SbrGrid, frameClass, v_bordFollow,
                       *length_v_bordFollow, v_freqFollow, *length_v_freqFollow,
                       DC, *i_tranFollow, *spreadFlag, DC);
        break;
      case VARVAR:
        if (*spreadFlag) {
          calcCtrlSignal(&hSbrEnvFrame->SbrGrid, frameClass, v_bordFollow,
                         *length_v_bordFollow, v_freqFollow,
                         *length_v_freqFollow, DC, *i_tranFollow, *spreadFlag,
                         DC);

          *spreadFlag = 0;

          v_bordFollow[0] = hSbrEnvFrame->SbrGrid.bs_abs_bord_1 -
                            numberTimeSlots;
          v_freqFollow[0] = 1;
          *length_v_bordFollow = 1;
          *length_v_freqFollow = 1;

          *i_tranFollow = -DC;
          *i_fillFollow = -DC;
        } else {
          fillFrameInter(&nL, v_tuningSegm, v_bord, length_v_bord, bmin, v_freq,
                         length_v_freq, v_bordFollow, length_v_bordFollow,
                         v_freqFollow, length_v_freqFollow, *i_fillFollow, dmin,
                         dmax, numberTimeSlots);

          fillFramePost(&parts, &d, dmax, v_bord, length_v_bord, v_freq,
                        length_v_freq, bmax, bufferFrameStart, numberTimeSlots,
                        fmax);

          if (parts == 1 && d < dmin)
            specialCase(spreadFlag, allowSpread, v_bord, length_v_bord, v_freq,
                        length_v_freq, &parts, d);

          calcCmonBorder(&i_cmon, &i_tran, v_bord, length_v_bord,
                         tranPosInternal, bufferFrameStart, numberTimeSlots);

          keepForFollowUp(v_bordFollow, length_v_bordFollow, v_freqFollow,
                          length_v_freqFollow, i_tranFollow, i_fillFollow,
                          v_bord, length_v_bord, v_freq, i_cmon, i_tran, parts,
                          numberTimeSlots);

          calcCtrlSignal(&hSbrEnvFrame->SbrGrid, frameClass, v_bord,
                         *length_v_bord, v_freq, *length_v_freq, i_cmon, i_tran,
                         0, nL);
        }
        break;
      case FIXFIX:
        if (tranPos == 0)
          numEnv = 1;
        else
          numEnv = 2;

        hSbrEnvFrame->SbrGrid.bs_num_env = numEnv;
        hSbrEnvFrame->SbrGrid.frameClass = frameClass;

        break;
      default:
        FDK_ASSERT(0);
    }
  }

  ctrlSignal2FrameInfo(&hSbrEnvFrame->SbrGrid, &hSbrEnvFrame->SbrFrameInfo,
                       hSbrEnvFrame->freq_res_fixfix);

  return &hSbrEnvFrame->SbrFrameInfo;
}

static void generateFixFixOnly(HANDLE_SBR_FRAME_INFO hSbrFrameInfo,
                               HANDLE_SBR_GRID hSbrGrid, int tranPosInternal,
                               int numberTimeSlots, UCHAR fResTransIsLow) {
  int nEnv, i, k = 0, tranIdx;
  const int *pTable = NULL;
  const FREQ_RES *freqResTable = NULL;

  switch (numberTimeSlots) {
    case 8: {
      pTable = envelopeTable_8[tranPosInternal];
    }
      freqResTable = freqRes_table_8;
      break;
    case 15:
      pTable = envelopeTable_15[tranPosInternal];
      freqResTable = freqRes_table_16;
      break;
    case 16:
      pTable = envelopeTable_16[tranPosInternal];
      freqResTable = freqRes_table_16;
      break;
  }

  nEnv = pTable[0];
  for (i = 1; i < nEnv; i++) hSbrFrameInfo->borders[i] = pTable[i + 2];

  hSbrFrameInfo->borders[0] = 0;
  hSbrFrameInfo->borders[nEnv] = numberTimeSlots;

  for (i = 0; i < nEnv; i++) {
    k = hSbrFrameInfo->borders[i + 1] - hSbrFrameInfo->borders[i];
    if (!fResTransIsLow)
      hSbrFrameInfo->freqRes[i] = freqResTable[k];
    else
      hSbrFrameInfo->freqRes[i] = FREQ_RES_LOW;

    hSbrGrid->v_f[i] = hSbrFrameInfo->freqRes[i];
  }

  hSbrFrameInfo->nEnvelopes = nEnv;
  hSbrFrameInfo->shortEnv = pTable[2];
  tranIdx = pTable[1];

  hSbrFrameInfo->bordersNoise[0] = 0;
  hSbrFrameInfo->bordersNoise[1] =
      hSbrFrameInfo->borders[tranIdx ? tranIdx : 1];
  hSbrFrameInfo->bordersNoise[2] = numberTimeSlots;
  hSbrFrameInfo->nNoiseEnvelopes = 2;

  hSbrGrid->frameClass = FIXFIXonly;
  hSbrGrid->bs_abs_bord = tranPosInternal;
  hSbrGrid->bs_num_env = nEnv;
}

void FDKsbrEnc_initFrameInfoGenerator(HANDLE_SBR_ENVELOPE_FRAME hSbrEnvFrame,
                                      INT allowSpread, INT numEnvStatic,
                                      INT staticFraming, INT timeSlots,
                                      const FREQ_RES *freq_res_fixfix,
                                      UCHAR fResTransIsLow,
                                      INT ldGrid) {

  FDKmemclear(hSbrEnvFrame, sizeof(SBR_ENVELOPE_FRAME));

  hSbrEnvFrame->frameClassOld = FIXFIX;
  hSbrEnvFrame->spreadFlag = 0;

  hSbrEnvFrame->allowSpread = allowSpread;
  hSbrEnvFrame->numEnvStatic = numEnvStatic;
  hSbrEnvFrame->staticFraming = staticFraming;
  hSbrEnvFrame->freq_res_fixfix[0] = freq_res_fixfix[0];
  hSbrEnvFrame->freq_res_fixfix[1] = freq_res_fixfix[1];
  hSbrEnvFrame->fResTransIsLow = fResTransIsLow;

  hSbrEnvFrame->length_v_bord = 0;
  hSbrEnvFrame->length_v_bordFollow = 0;

  hSbrEnvFrame->length_v_freq = 0;
  hSbrEnvFrame->length_v_freqFollow = 0;

  hSbrEnvFrame->i_tranFollow = 0;
  hSbrEnvFrame->i_fillFollow = 0;

  hSbrEnvFrame->SbrGrid.numberTimeSlots = timeSlots;

  if (ldGrid) {
    hSbrEnvFrame->dmin = 2;
    hSbrEnvFrame->dmax = 16;
    hSbrEnvFrame->frameMiddleSlot = FRAME_MIDDLE_SLOT_512LD;
    hSbrEnvFrame->SbrGrid.bufferFrameStart = 0;
  } else
    switch (timeSlots) {
      case NUMBER_TIME_SLOTS_1920:
        hSbrEnvFrame->dmin = 4;
        hSbrEnvFrame->dmax = 12;
        hSbrEnvFrame->SbrGrid.bufferFrameStart = 0;
        hSbrEnvFrame->frameMiddleSlot = FRAME_MIDDLE_SLOT_1920;
        break;
      case NUMBER_TIME_SLOTS_2048:
        hSbrEnvFrame->dmin = 4;
        hSbrEnvFrame->dmax = 12;
        hSbrEnvFrame->SbrGrid.bufferFrameStart = 0;
        hSbrEnvFrame->frameMiddleSlot = FRAME_MIDDLE_SLOT_2048;
        break;
      case NUMBER_TIME_SLOTS_1152:
        hSbrEnvFrame->dmin = 2;
        hSbrEnvFrame->dmax = 8;
        hSbrEnvFrame->SbrGrid.bufferFrameStart = 0;
        hSbrEnvFrame->frameMiddleSlot = FRAME_MIDDLE_SLOT_1152;
        break;
      case NUMBER_TIME_SLOTS_2304:
        hSbrEnvFrame->dmin = 4;
        hSbrEnvFrame->dmax = 15;
        hSbrEnvFrame->SbrGrid.bufferFrameStart = 0;
        hSbrEnvFrame->frameMiddleSlot = FRAME_MIDDLE_SLOT_2304;
        break;
      default:
        FDK_ASSERT(0);
    }
}

static void fillFrameTran(
    const int *v_tuningSegm,
    const int *v_tuningFreq,
    int tran,
    int *v_bord,
    int *length_v_bord,
    int *v_freq,
    int *length_v_freq,
    int *bmin,
    int *bmax
) {
  int bord, i;

  *length_v_bord = 0;
  *length_v_freq = 0;

  if (v_tuningSegm[0]) {
    FDKsbrEnc_AddRight(v_bord, length_v_bord, (tran - v_tuningSegm[0]));
    FDKsbrEnc_AddRight(v_freq, length_v_freq, v_tuningFreq[0]);
  }

  bord = tran;
  FDKsbrEnc_AddRight(v_bord, length_v_bord, tran);

  if (v_tuningSegm[1]) {
    bord += v_tuningSegm[1];

    FDKsbrEnc_AddRight(v_bord, length_v_bord, bord);
    FDKsbrEnc_AddRight(v_freq, length_v_freq, v_tuningFreq[1]);
  }

  if (v_tuningSegm[2] != 0) {
    bord += v_tuningSegm[2];

    FDKsbrEnc_AddRight(v_bord, length_v_bord, bord);
    FDKsbrEnc_AddRight(v_freq, length_v_freq, v_tuningFreq[2]);
  }

  FDKsbrEnc_AddRight(v_freq, length_v_freq, 1);

  *bmin = v_bord[0];
  for (i = 0; i < *length_v_bord; i++)
    if (v_bord[i] < *bmin) *bmin = v_bord[i];

  *bmax = v_bord[0];
  for (i = 0; i < *length_v_bord; i++)
    if (v_bord[i] > *bmax) *bmax = v_bord[i];
}

static void fillFramePre(INT dmax, INT *v_bord, INT *length_v_bord, INT *v_freq,
                         INT *length_v_freq, INT bmin, INT rest) {
  INT parts, d, j, S, s = 0, segm, bord;

  parts = 1;
  d = rest;

  while (d > dmax) {
    parts++;

    segm = rest / parts;
    S = (segm - 2) >> 1;
    s = fixMin(8, 2 * S + 2);
    d = rest - (parts - 1) * s;
  }

  bord = bmin;

  for (j = 0; j <= parts - 2; j++) {
    bord = bord - s;

    FDKsbrEnc_AddLeft(v_bord, length_v_bord, bord);
    FDKsbrEnc_AddLeft(v_freq, length_v_freq, 1);
  }
}

static int calcFillLengthMax(
    int tranPos,
    int numberTimeSlots
) {
  int fmax;

  switch (numberTimeSlots) {
    case NUMBER_TIME_SLOTS_2048:
      if (tranPos < 4)
        fmax = 6;
      else if (tranPos == 4 || tranPos == 5)
        fmax = 4;
      else
        fmax = 8;
      break;

    case NUMBER_TIME_SLOTS_1920:
      if (tranPos < 4)
        fmax = 5;
      else if (tranPos == 4 || tranPos == 5)
        fmax = 3;
      else
        fmax = 7;
      break;

    default:
      fmax = 8;
      break;
  }

  return fmax;
}

static void fillFramePost(INT *parts, INT *d, INT dmax, INT *v_bord,
                          INT *length_v_bord, INT *v_freq, INT *length_v_freq,
                          INT bmax, INT bufferFrameStart, INT numberTimeSlots,
                          INT fmax) {
  INT j, rest, segm, S, s = 0, bord;

  rest = bufferFrameStart + 2 * numberTimeSlots - bmax;
  *d = rest;

  if (*d > 0) {
    *parts = 1;

    while (*d > dmax) {
      *parts = *parts + 1;

      segm = rest / (*parts);
      S = (segm - 2) >> 1;
      s = fixMin(fmax, 2 * S + 2);
      *d = rest - (*parts - 1) * s;
    }

    bord = bmax;
    for (j = 0; j <= *parts - 2; j++) {
      bord += s;

      FDKsbrEnc_AddRight(v_bord, length_v_bord, bord);
      FDKsbrEnc_AddRight(v_freq, length_v_freq, 1);
    }
  } else {
    *parts = 1;

    *length_v_bord = *length_v_bord - 1;
    *length_v_freq = *length_v_freq - 1;
  }
}

static void fillFrameInter(INT *nL, const int *v_tuningSegm, INT *v_bord,
                           INT *length_v_bord, INT bmin, INT *v_freq,
                           INT *length_v_freq, INT *v_bordFollow,
                           INT *length_v_bordFollow, INT *v_freqFollow,
                           INT *length_v_freqFollow, INT i_fillFollow, INT dmin,
                           INT dmax, INT numberTimeSlots) {
  INT middle, b_new, numBordFollow, bordMaxFollow, i;

  if (numberTimeSlots != NUMBER_TIME_SLOTS_1152) {
    if (i_fillFollow >= 1) {
      *length_v_bordFollow = i_fillFollow;
      *length_v_freqFollow = i_fillFollow;
    }

    numBordFollow = *length_v_bordFollow;
    bordMaxFollow = v_bordFollow[numBordFollow - 1];

    middle = bmin - bordMaxFollow;
    while (middle < 0) {
      numBordFollow--;
      bordMaxFollow = v_bordFollow[numBordFollow - 1];
      middle = bmin - bordMaxFollow;
    }

    *length_v_bordFollow = numBordFollow;
    *length_v_freqFollow = numBordFollow;
    *nL = numBordFollow - 1;

    b_new = *length_v_bord;

    if (middle <= dmax) {
      if (middle >= dmin) {
        FDKsbrEnc_AddVecLeft(v_bord, length_v_bord, v_bordFollow,
                             *length_v_bordFollow);
        FDKsbrEnc_AddVecLeft(v_freq, length_v_freq, v_freqFollow,
                             *length_v_freqFollow);
      }

      else {
        if (v_tuningSegm[0] != 0) {
          *length_v_bord = b_new - 1;
          FDKsbrEnc_AddVecLeft(v_bord, length_v_bord, v_bordFollow,
                               *length_v_bordFollow);

          *length_v_freq = b_new - 1;
          FDKsbrEnc_AddVecLeft(v_freq + 1, length_v_freq, v_freqFollow,
                               *length_v_freqFollow);
        } else {
          if (*length_v_bordFollow >
              1) {
            FDKsbrEnc_AddVecLeft(v_bord, length_v_bord, v_bordFollow,
                                 *length_v_bordFollow - 1);
            FDKsbrEnc_AddVecLeft(v_freq, length_v_freq, v_freqFollow,
                                 *length_v_bordFollow - 1);

            *nL = *nL - 1;
          } else {

            for (i = 0; i < *length_v_bord - 1; i++) v_bord[i] = v_bord[i + 1];
            for (i = 0; i < *length_v_freq - 1; i++) v_freq[i] = v_freq[i + 1];

            *length_v_bord = b_new - 1;
            *length_v_freq = b_new - 1;

            FDKsbrEnc_AddVecLeft(v_bord, length_v_bord, v_bordFollow,
                                 *length_v_bordFollow);
            FDKsbrEnc_AddVecLeft(v_freq, length_v_freq, v_freqFollow,
                                 *length_v_freqFollow);
          }
        }
      }
    } else {

      fillFramePre(dmax, v_bord, length_v_bord, v_freq, length_v_freq, bmin,
                   middle);
      FDKsbrEnc_AddVecLeft(v_bord, length_v_bord, v_bordFollow,
                           *length_v_bordFollow);
      FDKsbrEnc_AddVecLeft(v_freq, length_v_freq, v_freqFollow,
                           *length_v_freqFollow);
    }

  } else {

    INT l, m;

    if (i_fillFollow >= 1) {
      *length_v_bordFollow = i_fillFollow;
      *length_v_freqFollow = i_fillFollow;
    }

    numBordFollow = *length_v_bordFollow;
    bordMaxFollow = v_bordFollow[numBordFollow - 1];

    middle = bmin - bordMaxFollow;

    while (middle < 0) {
      if (numBordFollow == 1) break;

      numBordFollow--;
      bordMaxFollow = v_bordFollow[numBordFollow - 1];
      middle = bmin - bordMaxFollow;
    }

    if (middle < 0) {
      for (l = 0, m = 0; l < *length_v_bord; l++) {
        if (v_bord[l] > bordMaxFollow) {
          v_bord[m] = v_bord[l];
          v_freq[m] = v_freq[l];
          m++;
        }
      }

      *length_v_bord = l;
      *length_v_freq = l;

      bmin = v_bord[0];
    }

    *length_v_bordFollow = numBordFollow;
    *length_v_freqFollow = numBordFollow;

    *nL = numBordFollow - 1;

    middle = bmin - bordMaxFollow;

    if (middle <= dmin)
    {
      b_new = *length_v_bord;

      if (v_tuningSegm[0] != 0) {
        *length_v_bord = b_new - 1;
        FDKsbrEnc_AddVecLeft(v_bord, length_v_bord, v_bordFollow,
                             *length_v_bordFollow);

        *length_v_freq = b_new - 1;
        FDKsbrEnc_AddVecLeft(v_freq + 1, length_v_freq, v_freqFollow,
                             *length_v_freqFollow);

      } else if (*length_v_bordFollow > 1) {
        FDKsbrEnc_AddVecLeft(v_bord, length_v_bord, v_bordFollow,
                             *length_v_bordFollow - 1);
        FDKsbrEnc_AddVecLeft(v_freq, length_v_freq, v_freqFollow,
                             *length_v_bordFollow - 1);

        *nL = *nL - 1;
      } else {
        for (i = 0; i < *length_v_bord - 1; i++) v_bord[i] = v_bord[i + 1];

        for (i = 0; i < *length_v_freq - 1; i++) v_freq[i] = v_freq[i + 1];

        *length_v_bord = b_new - 1;
        *length_v_freq = b_new - 1;

        FDKsbrEnc_AddVecLeft(v_bord, length_v_bord, v_bordFollow,
                             *length_v_bordFollow);
        FDKsbrEnc_AddVecLeft(v_freq, length_v_freq, v_freqFollow,
                             *length_v_freqFollow);
      }
    } else if ((middle >= dmin) && (middle <= dmax))
    {
      FDKsbrEnc_AddVecLeft(v_bord, length_v_bord, v_bordFollow,
                           *length_v_bordFollow);
      FDKsbrEnc_AddVecLeft(v_freq, length_v_freq, v_freqFollow,
                           *length_v_freqFollow);
    } else
    {
      fillFramePre(dmax, v_bord, length_v_bord, v_freq, length_v_freq, bmin,
                   middle);
      FDKsbrEnc_AddVecLeft(v_bord, length_v_bord, v_bordFollow,
                           *length_v_bordFollow);
      FDKsbrEnc_AddVecLeft(v_freq, length_v_freq, v_freqFollow,
                           *length_v_freqFollow);
    }
  }
}

static void calcFrameClass(FRAME_CLASS *frameClass, FRAME_CLASS *frameClassOld,
                           INT tranFlag, INT *spreadFlag) {
  switch (*frameClassOld) {
    case FIXFIXonly:
    case FIXFIX:
      if (tranFlag)
        *frameClass = FIXVAR;
      else
        *frameClass = FIXFIX;
      break;
    case FIXVAR:
      if (tranFlag) {
        *frameClass = VARVAR;
        *spreadFlag = 0;
      } else {
        if (*spreadFlag)
          *frameClass = VARVAR;
        else
          *frameClass = VARFIX;
      }
      break;
    case VARFIX:
      if (tranFlag)
        *frameClass = FIXVAR;
      else
        *frameClass = FIXFIX;
      break;
    case VARVAR:
      if (tranFlag) {
        *frameClass = VARVAR;
        *spreadFlag = 0;
      } else {
        if (*spreadFlag)
          *frameClass = VARVAR;
        else
          *frameClass = VARFIX;
      }
      break;
  };

  *frameClassOld = *frameClass;
}

static void specialCase(INT *spreadFlag, INT allowSpread, INT *v_bord,
                        INT *length_v_bord, INT *v_freq, INT *length_v_freq,
                        INT *parts, INT d) {
  INT L;

  L = *length_v_bord;

  if (allowSpread) {
    *spreadFlag = 1;
    FDKsbrEnc_AddRight(v_bord, length_v_bord, v_bord[L - 1] + 8);
    FDKsbrEnc_AddRight(v_freq, length_v_freq, 1);
    (*parts)++;
  } else {
    if (d == 1) {
      *length_v_bord = L - 1;
      *length_v_freq = L - 1;
    } else {
      if ((v_bord[L - 1] - v_bord[L - 2]) > 2) {
        v_bord[L - 1] = v_bord[L - 1] - 2;
        v_freq[*length_v_freq - 1] = 0;
      }
    }
  }
}

static void calcCmonBorder(INT *i_cmon, INT *i_tran, INT *v_bord,
                           INT *length_v_bord, INT tran, INT bufferFrameStart,
                           INT numberTimeSlots) {
  INT i;

  for (i = 0; i < *length_v_bord; i++)
    if (v_bord[i] >= bufferFrameStart + numberTimeSlots) {
      *i_cmon = i;
      break;
    }

  for (i = 0; i < *length_v_bord; i++)
    if (v_bord[i] >= tran) {
      *i_tran = i;
      break;
    } else
      *i_tran = EMPTY;
}

static void keepForFollowUp(INT *v_bordFollow, INT *length_v_bordFollow,
                            INT *v_freqFollow, INT *length_v_freqFollow,
                            INT *i_tranFollow, INT *i_fillFollow, INT *v_bord,
                            INT *length_v_bord, INT *v_freq, INT i_cmon,
                            INT i_tran, INT parts,
                            INT numberTimeSlots) {
  INT L, i, j;

  L = *length_v_bord;

  (*length_v_bordFollow) = 0;
  (*length_v_freqFollow) = 0;

  for (j = 0, i = i_cmon; i < L; i++, j++) {
    v_bordFollow[j] = v_bord[i] - numberTimeSlots;
    v_freqFollow[j] = v_freq[i];
    (*length_v_bordFollow)++;
    (*length_v_freqFollow)++;
  }
  if (i_tran != EMPTY)
    *i_tranFollow = i_tran - i_cmon;
  else
    *i_tranFollow = EMPTY;
  *i_fillFollow = L - (parts - 1) - i_cmon;
}

static void calcCtrlSignal(HANDLE_SBR_GRID hSbrGrid, FRAME_CLASS frameClass,
                           INT *v_bord, INT length_v_bord, INT *v_freq,
                           INT length_v_freq, INT i_cmon, INT i_tran,
                           INT spreadFlag, INT nL) {
  INT i, r, a, n, p, b, aL, aR, ntot, nmax, nR;

  INT *v_f = hSbrGrid->v_f;
  INT *v_fLR = hSbrGrid->v_fLR;
  INT *v_r = hSbrGrid->bs_rel_bord;
  INT *v_rL = hSbrGrid->bs_rel_bord_0;
  INT *v_rR = hSbrGrid->bs_rel_bord_1;

  INT length_v_r = 0;
  INT length_v_rR = 0;
  INT length_v_rL = 0;

  switch (frameClass) {
    case FIXVAR:
      a = v_bord[i_cmon];

      length_v_r = 0;
      i = i_cmon;

      while (i >= 1) {
        r = v_bord[i] - v_bord[i - 1];
        FDKsbrEnc_AddRight(v_r, &length_v_r, r);
        i--;
      }

      n = length_v_r;

      for (i = 0; i < i_cmon; i++) v_f[i] = v_freq[i_cmon - 1 - i];
      v_f[i_cmon] = 1;

      p = (i_cmon >= i_tran && i_tran != EMPTY) ? (i_cmon - i_tran + 1) : (0);

      hSbrGrid->frameClass = frameClass;
      hSbrGrid->bs_abs_bord = a;
      hSbrGrid->n = n;
      hSbrGrid->p = p;

      break;
    case VARFIX:
      a = v_bord[0];

      length_v_r = 0;

      for (i = 1; i < length_v_bord; i++) {
        r = v_bord[i] - v_bord[i - 1];
        FDKsbrEnc_AddRight(v_r, &length_v_r, r);
      }

      n = length_v_r;

      FDKmemcpy(v_f, v_freq, length_v_freq * sizeof(INT));

      p = (i_tran >= 0 && i_tran != EMPTY) ? (i_tran + 1) : (0);

      hSbrGrid->frameClass = frameClass;
      hSbrGrid->bs_abs_bord = a;
      hSbrGrid->n = n;
      hSbrGrid->p = p;

      break;
    case VARVAR:
      if (spreadFlag) {
        b = length_v_bord;

        aL = v_bord[0];
        aR = v_bord[b - 1];

        ntot = b - 2;

        nmax = 2;
        if (ntot > nmax) {
          nL = nmax;
          nR = ntot - nmax;
        } else {
          nL = ntot;
          nR = 0;
        }

        length_v_rL = 0;
        for (i = 1; i <= nL; i++) {
          r = v_bord[i] - v_bord[i - 1];
          FDKsbrEnc_AddRight(v_rL, &length_v_rL, r);
        }

        length_v_rR = 0;
        i = b - 1;
        while (i >= b - nR) {
          r = v_bord[i] - v_bord[i - 1];
          FDKsbrEnc_AddRight(v_rR, &length_v_rR, r);
          i--;
        }

        p = (i_tran > 0 && i_tran != EMPTY) ? (b - i_tran) : (0);

        for (i = 0; i < b - 1; i++) v_fLR[i] = v_freq[i];
      } else {
        length_v_bord = i_cmon + 1;

        b = length_v_bord;

        aL = v_bord[0];
        aR = v_bord[b - 1];

        ntot = b - 2;
        nR = ntot - nL;

        length_v_rL = 0;
        for (i = 1; i <= nL; i++) {
          r = v_bord[i] - v_bord[i - 1];
          FDKsbrEnc_AddRight(v_rL, &length_v_rL, r);
        }

        length_v_rR = 0;
        i = b - 1;
        while (i >= b - nR) {
          r = v_bord[i] - v_bord[i - 1];
          FDKsbrEnc_AddRight(v_rR, &length_v_rR, r);
          i--;
        }

        p = (i_cmon >= i_tran && i_tran != EMPTY) ? (i_cmon - i_tran + 1) : (0);

        for (i = 0; i < b - 1; i++) v_fLR[i] = v_freq[i];
      }

      hSbrGrid->frameClass = frameClass;
      hSbrGrid->bs_abs_bord_0 = aL;
      hSbrGrid->bs_abs_bord_1 = aR;
      hSbrGrid->bs_num_rel_0 = nL;
      hSbrGrid->bs_num_rel_1 = nR;
      hSbrGrid->p = p;

      break;

    default:
      break;
  }
}

static void createDefFrameInfo(HANDLE_SBR_FRAME_INFO hSbrFrameInfo, INT nEnv,
                               INT nTimeSlots) {
  switch (nEnv) {
    case 1:
      switch (nTimeSlots) {
        case NUMBER_TIME_SLOTS_1920:
          FDKmemcpy(hSbrFrameInfo, &frameInfo1_1920, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_2048:
          FDKmemcpy(hSbrFrameInfo, &frameInfo1_2048, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_1152:
          FDKmemcpy(hSbrFrameInfo, &frameInfo1_1152, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_2304:
          FDKmemcpy(hSbrFrameInfo, &frameInfo1_2304, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_512LD:
          FDKmemcpy(hSbrFrameInfo, &frameInfo1_512LD, sizeof(SBR_FRAME_INFO));
          break;
        default:
          FDK_ASSERT(0);
      }
      break;
    case 2:
      switch (nTimeSlots) {
        case NUMBER_TIME_SLOTS_1920:
          FDKmemcpy(hSbrFrameInfo, &frameInfo2_1920, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_2048:
          FDKmemcpy(hSbrFrameInfo, &frameInfo2_2048, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_1152:
          FDKmemcpy(hSbrFrameInfo, &frameInfo2_1152, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_2304:
          FDKmemcpy(hSbrFrameInfo, &frameInfo2_2304, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_512LD:
          FDKmemcpy(hSbrFrameInfo, &frameInfo2_512LD, sizeof(SBR_FRAME_INFO));
          break;
        default:
          FDK_ASSERT(0);
      }
      break;
    case 4:
      switch (nTimeSlots) {
        case NUMBER_TIME_SLOTS_1920:
          FDKmemcpy(hSbrFrameInfo, &frameInfo4_1920, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_2048:
          FDKmemcpy(hSbrFrameInfo, &frameInfo4_2048, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_1152:
          FDKmemcpy(hSbrFrameInfo, &frameInfo4_1152, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_2304:
          FDKmemcpy(hSbrFrameInfo, &frameInfo4_2304, sizeof(SBR_FRAME_INFO));
          break;
        case NUMBER_TIME_SLOTS_512LD:
          FDKmemcpy(hSbrFrameInfo, &frameInfo4_512LD, sizeof(SBR_FRAME_INFO));
          break;
        default:
          FDK_ASSERT(0);
      }
      break;
    default:
      FDK_ASSERT(0);
  }
}

static void ctrlSignal2FrameInfo(
    HANDLE_SBR_GRID hSbrGrid,
    HANDLE_SBR_FRAME_INFO hSbrFrameInfo,
    FREQ_RES
        *freq_res_fixfix
) {
  INT frameSplit = 0;
  INT nEnv = 0, border = 0, i, k, p;
  INT *v_r = hSbrGrid->bs_rel_bord;
  INT *v_f = hSbrGrid->v_f;

  FRAME_CLASS frameClass = hSbrGrid->frameClass;
  INT bufferFrameStart = hSbrGrid->bufferFrameStart;
  INT numberTimeSlots = hSbrGrid->numberTimeSlots;

  switch (frameClass) {
    case FIXFIX:
      createDefFrameInfo(hSbrFrameInfo, hSbrGrid->bs_num_env, numberTimeSlots);

      frameSplit = (hSbrFrameInfo->nEnvelopes > 1);
      for (i = 0; i < hSbrFrameInfo->nEnvelopes; i++) {
        hSbrGrid->v_f[i] = hSbrFrameInfo->freqRes[i] =
            freq_res_fixfix[frameSplit];
      }
      break;

    case FIXVAR:
    case VARFIX:
      nEnv = hSbrGrid->n + 1;
      FDK_ASSERT(nEnv <= MAX_ENVELOPES_FIXVAR_VARFIX);

      hSbrFrameInfo->nEnvelopes = nEnv;

      border = hSbrGrid->bs_abs_bord;

      if (nEnv == 1)
        hSbrFrameInfo->nNoiseEnvelopes = 1;
      else
        hSbrFrameInfo->nNoiseEnvelopes = 2;

      break;

    default:
      break;
  }

  switch (frameClass) {
    case FIXVAR:
      hSbrFrameInfo->borders[0] =
          bufferFrameStart;

      hSbrFrameInfo->borders[nEnv] = border;

      for (k = 0, i = nEnv - 1; k < nEnv - 1; k++, i--) {
        border -= v_r[k];
        hSbrFrameInfo->borders[i] = border;
      }

      p = hSbrGrid->p;
      if (p == 0) {
        hSbrFrameInfo->shortEnv = 0;
      } else {
        hSbrFrameInfo->shortEnv = nEnv + 1 - p;
      }

      for (k = 0, i = nEnv - 1; k < nEnv; k++, i--) {
        hSbrFrameInfo->freqRes[i] = (FREQ_RES)v_f[k];
      }

      if (p == 0 || p == 1) {
        hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[nEnv - 1];
      } else {
        hSbrFrameInfo->bordersNoise[1] =
            hSbrFrameInfo->borders[hSbrFrameInfo->shortEnv];
      }

      break;

    case VARFIX:
      hSbrFrameInfo->borders[0] = border;

      for (k = 0; k < nEnv - 1; k++) {
        border += v_r[k];
        hSbrFrameInfo->borders[k + 1] = border;
      }

      hSbrFrameInfo->borders[nEnv] = bufferFrameStart + numberTimeSlots;

      p = hSbrGrid->p;
      if (p == 0 || p == 1) {
        hSbrFrameInfo->shortEnv = 0;
      } else {
        hSbrFrameInfo->shortEnv = p - 1;
      }

      for (k = 0; k < nEnv; k++) {
        hSbrFrameInfo->freqRes[k] = (FREQ_RES)v_f[k];
      }

      switch (p) {
        case 0:
          hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[1];
          break;
        case 1:
          hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[nEnv - 1];
          break;
        default:
          hSbrFrameInfo->bordersNoise[1] =
              hSbrFrameInfo->borders[hSbrFrameInfo->shortEnv];
          break;
      }
      break;

    case VARVAR:
      nEnv = hSbrGrid->bs_num_rel_0 + hSbrGrid->bs_num_rel_1 + 1;
      FDK_ASSERT(nEnv <= MAX_ENVELOPES_VARVAR);
      hSbrFrameInfo->nEnvelopes = nEnv;

      hSbrFrameInfo->borders[0] = border = hSbrGrid->bs_abs_bord_0;

      for (k = 0, i = 1; k < hSbrGrid->bs_num_rel_0; k++, i++) {
        border += hSbrGrid->bs_rel_bord_0[k];
        hSbrFrameInfo->borders[i] = border;
      }

      border = hSbrGrid->bs_abs_bord_1;
      hSbrFrameInfo->borders[nEnv] = border;

      for (k = 0, i = nEnv - 1; k < hSbrGrid->bs_num_rel_1; k++, i--) {
        border -= hSbrGrid->bs_rel_bord_1[k];
        hSbrFrameInfo->borders[i] = border;
      }

      p = hSbrGrid->p;
      if (p == 0) {
        hSbrFrameInfo->shortEnv = 0;
      } else {
        hSbrFrameInfo->shortEnv = nEnv + 1 - p;
      }

      for (k = 0; k < nEnv; k++) {
        hSbrFrameInfo->freqRes[k] = (FREQ_RES)hSbrGrid->v_fLR[k];
      }

      if (nEnv == 1) {
        hSbrFrameInfo->nNoiseEnvelopes = 1;
        hSbrFrameInfo->bordersNoise[0] = hSbrGrid->bs_abs_bord_0;
        hSbrFrameInfo->bordersNoise[1] = hSbrGrid->bs_abs_bord_1;
      } else {
        hSbrFrameInfo->nNoiseEnvelopes = 2;
        hSbrFrameInfo->bordersNoise[0] = hSbrGrid->bs_abs_bord_0;

        if (p == 0 || p == 1) {
          hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[nEnv - 1];
        } else {
          hSbrFrameInfo->bordersNoise[1] =
              hSbrFrameInfo->borders[hSbrFrameInfo->shortEnv];
        }
        hSbrFrameInfo->bordersNoise[2] = hSbrGrid->bs_abs_bord_1;
      }
      break;

    default:
      break;
  }

  if (frameClass == VARFIX || frameClass == FIXVAR) {
    hSbrFrameInfo->bordersNoise[0] = hSbrFrameInfo->borders[0];
    if (nEnv == 1) {
      hSbrFrameInfo->bordersNoise[1] = hSbrFrameInfo->borders[nEnv];
    } else {
      hSbrFrameInfo->bordersNoise[2] = hSbrFrameInfo->borders[nEnv];
    }
  }
}
