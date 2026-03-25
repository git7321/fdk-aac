/**************************** SBR encoder library ******************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef FRAM_GEN_H
#define FRAM_GEN_H

#include "sbr_def.h"
#include "sbr_encoder.h"

#define MAX_ENVELOPES_VARVAR \
  MAX_ENVELOPES /*!< worst case number of envelopes in a VARVAR frame */
#define MAX_ENVELOPES_FIXVAR_VARFIX \
  4 /*!< worst case number of envelopes in VARFIX and FIXVAR frames */
#define MAX_NUM_REL \
  3 /*!< maximum number of relative borders in any VAR frame */

/* SBR frame class definitions */
typedef enum {
  FIXFIX =
      0,  /*!< bs_frame_class: leading and trailing frame borders are fixed */
  FIXVAR, /*!< bs_frame_class: leading frame border is fixed, trailing frame
             border is variable */
  VARFIX, /*!< bs_frame_class: leading frame border is variable, trailing frame
             border is fixed */
  VARVAR /*!< bs_frame_class: leading and trailing frame borders are variable */
  ,
  FIXFIXonly /*!< bs_frame_class: leading border fixed (0), trailing border
                fixed (nrTimeSlots) and encased borders are dynamically derived
                from the tranPos */
} FRAME_CLASS;

/* helper constants */
#define DC 4711     /*!< helper constant: don't care */
#define EMPTY (-99) /*!< helper constant: empty */

/* system constants: AAC+SBR, DRM Frame-Length */
#define FRAME_MIDDLE_SLOT_1920 4
#define NUMBER_TIME_SLOTS_1920 15

#define LD_PRETRAN_OFF 3
#define FRAME_MIDDLE_SLOT_512LD 4
#define NUMBER_TIME_SLOTS_512LD 8
#define TRANSIENT_OFFSET_LD 0

#define FRAME_MIDDLE_SLOT_2048 4
#define NUMBER_TIME_SLOTS_2048 16

#define FRAME_MIDDLE_SLOT_1152 4
#define NUMBER_TIME_SLOTS_1152 9

#define FRAME_MIDDLE_SLOT_2304 8
#define NUMBER_TIME_SLOTS_2304 18

typedef struct {
  /* system constants */
  INT bufferFrameStart; /*!< frame generator vs analysis buffer time alignment
                           (currently set to 0, offset added elsewhere) */
  INT numberTimeSlots;  /*!< number of SBR timeslots per frame */

  /* will be adjusted for every frame */
  FRAME_CLASS frameClass; /*!< SBR frame class  */
  INT bs_num_env;         /*!< bs_num_env, number of envelopes for FIXFIX */
  INT bs_abs_bord; /*!< bs_abs_bord, absolute border for VARFIX and FIXVAR */
  INT n;           /*!< number of relative borders for VARFIX and FIXVAR   */
  INT p;           /*!< pointer-to-transient-border  */
  INT bs_rel_bord[MAX_NUM_REL]; /*!< bs_rel_bord, relative borders for all VAR
                                 */
  INT v_f[MAX_ENVELOPES_FIXVAR_VARFIX]; /*!< envelope frequency resolutions for
                                           FIXVAR and VARFIX  */

  INT bs_abs_bord_0; /*!< bs_abs_bord_0, leading absolute border for VARVAR */
  INT bs_abs_bord_1; /*!< bs_abs_bord_1, trailing absolute border for VARVAR */
  INT bs_num_rel_0;  /*!< bs_num_rel_0, number of relative borders associated
                        with leading absolute border for VARVAR */
  INT bs_num_rel_1;  /*!< bs_num_rel_1, number of relative borders associated
                        with trailing absolute border for VARVAR */
  INT bs_rel_bord_0[MAX_NUM_REL];  /*!< bs_rel_bord_0, relative borders
                                      associated with  leading absolute border
                                      for  VARVAR */
  INT bs_rel_bord_1[MAX_NUM_REL];  /*!< bs_rel_bord_1, relative borders
                                      associated with trailing absolute border
                                      for VARVAR */
  INT v_fLR[MAX_ENVELOPES_VARVAR]; /*!< envelope frequency resolutions for
                                      VARVAR  */

} SBR_GRID;
typedef SBR_GRID *HANDLE_SBR_GRID;

typedef struct {
  INT nEnvelopes;                  /*!< number of envelopes */
  INT borders[MAX_ENVELOPES + 1];  /*!< envelope borders in SBR timeslots */
  FREQ_RES freqRes[MAX_ENVELOPES]; /*!< frequency resolution of each envelope */
  INT shortEnv; /*!< number of an envelope to be shortened (starting at 1) or 0
                   for no shortened envelope */
  INT nNoiseEnvelopes; /*!< number of noise floors */
  INT bordersNoise[MAX_NOISE_ENVELOPES +
                   1]; /*!< noise floor borders in SBR timeslots */
} SBR_FRAME_INFO;
typedef SBR_FRAME_INFO *HANDLE_SBR_FRAME_INFO;

typedef struct {
  /* system constants */
  INT frameMiddleSlot; /*!< transient detector offset in SBR timeslots */

  /* basic tuning parameters */
  INT staticFraming; /*!< 1: run static framing in time, i.e. exclusive use of
                        bs_frame_class = FIXFIX */
  INT numEnvStatic;  /*!< number of envelopes per frame for static framing */
  FREQ_RES
  freq_res_fixfix[2]; /*!< envelope frequency resolution to use for
                         bs_frame_class = FIXFIX; single env and split */
  UCHAR
  fResTransIsLow; /*!< frequency resolution for transient frames - always
                     low (0) or according to table (1) */

  /* expert tuning parameters */
  const int *v_tuningSegm; /*!< segment lengths to use around transient */
  const int *v_tuningFreq; /*!< frequency resolutions to use around transient */
  INT dmin;                /*!< minimum length of dependent segments */
  INT dmax;                /*!< maximum length of dependent segments */
  INT allowSpread; /*!< 1: allow isolated transient to influence grid of 3
                      consecutive frames */

  /* internally used signals */
  FRAME_CLASS frameClassOld; /*!< frame class used for previous frame */
  INT spreadFlag; /*!< 1: use VARVAR instead of VARFIX to follow up old
                     transient */

  INT v_bord[2 * MAX_ENVELOPES_VARVAR + 1]; /*!< borders for current frame and
                                               preliminary borders for next
                                               frame (fixed borders excluded) */
  INT length_v_bord; /*!< helper variable: length of v_bord */
  INT v_freq[2 * MAX_ENVELOPES_VARVAR + 1]; /*!< frequency resolutions for
                                               current frame and preliminary
                                               resolutions for next frame */
  INT length_v_freq; /*!< helper variable: length of v_freq */

  INT v_bordFollow[MAX_ENVELOPES_VARVAR]; /*!< preliminary borders for current
                                             frame (calculated during previous
                                             frame) */
  INT length_v_bordFollow; /*!< helper variable: length of v_bordFollow */
  INT i_tranFollow; /*!< points to transient border in v_bordFollow (may be
                       negative, see keepForFollowUp()) */
  INT i_fillFollow; /*!< points to first fill border in v_bordFollow */
  INT v_freqFollow[MAX_ENVELOPES_VARVAR]; /*!< preliminary frequency resolutions
                                             for current frame (calculated
                                             during previous frame) */
  INT length_v_freqFollow; /*!< helper variable: length of v_freqFollow */

  /* externally needed signals */
  SBR_GRID
  SbrGrid; /*!< sbr_grid() signals to be converted to bitstream elements */
  SBR_FRAME_INFO
  SbrFrameInfo; /*!< time/frequency grid description for one frame */
} SBR_ENVELOPE_FRAME;
typedef SBR_ENVELOPE_FRAME *HANDLE_SBR_ENVELOPE_FRAME;

void FDKsbrEnc_initFrameInfoGenerator(HANDLE_SBR_ENVELOPE_FRAME hSbrEnvFrame,
                                      INT allowSpread, INT numEnvStatic,
                                      INT staticFraming, INT timeSlots,
                                      const FREQ_RES *freq_res_fixfix,
                                      UCHAR fResTransIsLow, INT ldGrid);

HANDLE_SBR_FRAME_INFO
FDKsbrEnc_frameInfoGenerator(HANDLE_SBR_ENVELOPE_FRAME hSbrEnvFrame,
                             UCHAR *v_transient_info, const INT rightBorderFIX,
                             UCHAR *v_transient_info_pre, int ldGrid,
                             const int *v_tuning);

#endif
