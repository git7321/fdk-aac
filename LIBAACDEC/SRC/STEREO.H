/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: joint stereo processing

*******************************************************************************/

#ifndef STEREO_H
#define STEREO_H

#include "machine_type.h"
#include "FDK_bitstream.h"
#include "common_fix.h"

#define SFB_PER_PRED_BAND 2

#define SR_FNA_OUT \
  0

enum { JointStereoMaximumGroups = 8, JointStereoMaximumBands = 64 };

typedef struct {
  UCHAR pred_dir;
  UCHAR
  igf_pred_dir;
  UCHAR complex_coef;
  UCHAR use_prev_frame;

  SHORT alpha_q_re[JointStereoMaximumGroups][JointStereoMaximumBands];
  SHORT alpha_q_im[JointStereoMaximumGroups][JointStereoMaximumBands];
} CCplxPredictionData;

typedef struct {
  UCHAR MsMaskPresent;
  UCHAR MsUsed[JointStereoMaximumBands];
  UCHAR IGF_MsMaskPresent;

  UCHAR cplx_pred_flag;
  UCHAR igf_cplx_pred_flag;

  FIXP_DBL store_dmx_re_prev[1024];
  SHORT store_dmx_re_prev_e;

} CJointStereoData;

typedef struct {
  UCHAR clearSpectralCoeffs;

  FIXP_DBL *scratchBuffer;

  FIXP_DBL
  *spectralCoeffs[2];
  SHORT *specScale[2];

  SHORT alpha_q_re_prev[JointStereoMaximumGroups][JointStereoMaximumBands];
  SHORT alpha_q_im_prev[JointStereoMaximumGroups][JointStereoMaximumBands];

  UCHAR winSeqPrev;
  UCHAR winShapePrev;
  UCHAR winGroupsPrev;

} CJointStereoPersistentData;

int CJointStereo_Read(HANDLE_FDK_BITSTREAM bs,
                      CJointStereoData *pJointStereoData,
                      const int windowGroups,
                      const int scaleFactorBandsTransmitted,
                      const int max_sfb_ste_clear,
                      CJointStereoPersistentData *pJointStereoPersistentData,
                      CCplxPredictionData *cplxPredictionData,
                      int cplxPredictionActiv, int scaleFactorBandsTotal,
                      int windowSequence, const UINT flags);

#endif /* #ifndef STEREO_H */
