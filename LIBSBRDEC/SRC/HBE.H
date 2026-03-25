/**************************** SBR decoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef HBE_H
#define HBE_H

#include "sbrdecoder.h"

#ifndef QMF_SYNTH_CHANNELS
#define QMF_SYNTH_CHANNELS (64)
#endif

#define HBE_QMF_FILTER_STATE_ANA_SIZE (400)
#define HBE_QMF_FILTER_STATE_SYN_SIZE (200)

#ifndef MAX_NUM_PATCHES_HBE
#define MAX_NUM_PATCHES_HBE (6)
#endif
#define MAX_STRETCH_HBE (4)

typedef enum {
  KEEP_STATES_SYNCED_OFF = 0,
  KEEP_STATES_SYNCED_NORMAL = 1,
  KEEP_STATES_SYNCED_OUTDIFF =
      2,
  KEEP_STATES_SYNCED_NOOUT =
      3
} KEEP_STATES_SYNCED_MODE;

struct hbeTransposer {
  FIXP_DBL anaQmfStates[HBE_QMF_FILTER_STATE_ANA_SIZE];
  FIXP_QSS synQmfStates[HBE_QMF_FILTER_STATE_SYN_SIZE];

  int xOverQmf[MAX_NUM_PATCHES_HBE];

  int maxStretch;
  int timeDomainWinLen;
  int qmfInBufSize;
  int qmfOutBufSize;
  int noCols;
  int noChannels;
  int startBand;
  int stopBand;
  int bSbr41;

  LONG *inBuf_F;
  FIXP_DBL **qmfInBufReal_F;
  FIXP_DBL **qmfInBufImag_F;

  FIXP_DBL *qmfBufferCodecTempSlot_F;

  QMF_FILTER_BANK HBEAnalysiscQMF;
  QMF_FILTER_BANK HBESynthesisQMF;

  FIXP_DBL const *synthesisQmfPreModCos_F;
  FIXP_DBL const *synthesisQmfPreModSin_F;

  FIXP_DBL **qmfHBEBufReal_F;
  FIXP_DBL **qmfHBEBufImag_F;

  int bXProducts[MAX_STRETCH_HBE];

  int kstart;
  int synthSize;

  int highband_exp[2];
  int target_exp[2];
};

typedef struct hbeTransposer *HANDLE_HBE_TRANSPOSER;

SBR_ERROR QmfTransposerCreate(HANDLE_HBE_TRANSPOSER *hQmfTransposer,
                              const int frameSize, int bDisableCrossProducts,
                              int bSbr41);
SBR_ERROR QmfTransposerReInit(HANDLE_HBE_TRANSPOSER hQmfTransposer,
                              UCHAR *FreqBandTable[2], UCHAR NSfb[2]);

void QmfTransposerClose(HANDLE_HBE_TRANSPOSER hQmfTransposer);
void QmfTransposerApply(HANDLE_HBE_TRANSPOSER hQmfTransposer,
                        FIXP_DBL **qmfBufferCodecReal,
                        FIXP_DBL **qmfBufferCodecImag, int nColsIn,
                        FIXP_DBL **ppQmfBufferOutReal_F,
                        FIXP_DBL **ppQmfBufferOutImag_F,
                        FIXP_DBL lpcFilterStatesReal[2 + (3 * (4))][(64)],
                        FIXP_DBL lpcFilterStatesImag[2 + (3 * (4))][(64)],
                        int pitchInBins, int scale_lb, int scale_hbe,
                        int *scale_hb, int timeStep, int firstSlotOffsset,
                        int ov_len,
                        KEEP_STATES_SYNCED_MODE keepStatesSyncedMode);

int *GetxOverBandQmfTransposer(HANDLE_HBE_TRANSPOSER hQmfTransposer);
int Get41SbrQmfTransposer(HANDLE_HBE_TRANSPOSER hQmfTransposer);

#endif /* HBE_H */
