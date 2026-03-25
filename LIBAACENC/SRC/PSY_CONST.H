/**************************** AAC encoder library ******************************

   Author(s):   M.Werner

   Description: Global psychoaccoustic constants

*******************************************************************************/

#ifndef PSY_CONST_H
#define PSY_CONST_H

#define TRUE 1
#define FALSE 0

#define TRANS_FAC 8

#define FRAME_MAXLEN_SHORT ((1024) / TRANS_FAC)
#define FRAME_LEN_SHORT_128 ((1024) / TRANS_FAC)
#define FRAME_LEN_SHORT_120 (FRAME_LEN_LONG_960 / TRANS_FAC)

enum FB_TYPE { FB_LC = 0, FB_LD = 1, FB_ELD = 2 };

#define N_BLOCKTYPES 6
enum {
  LONG_WINDOW = 0,
  START_WINDOW,
  SHORT_WINDOW,
  STOP_WINDOW,
  _LOWOV_WINDOW,
  WRONG_WINDOW
};

enum {
  SINE_WINDOW = 0,
  KBD_WINDOW = 1,
  LOL_WINDOW = 2
};

enum { SI_MS_MASK_NONE = 0, SI_MS_MASK_SOME = 1, SI_MS_MASK_ALL = 2 };

#define MAX_NO_OF_GROUPS 4
#define MAX_SFB_LONG \
  51
#define MAX_SFB_SHORT \
  15

#define MAX_SFB \
  (MAX_SFB_SHORT > MAX_SFB_LONG ? MAX_SFB_SHORT : MAX_SFB_LONG)
#define MAX_GROUPED_SFB                            \
  (MAX_NO_OF_GROUPS * MAX_SFB_SHORT > MAX_SFB_LONG \
       ? MAX_NO_OF_GROUPS * MAX_SFB_SHORT          \
       : MAX_SFB_LONG)

#define MAX_INPUT_BUFFER_SIZE (2 * (1024))

#define PCM_LEVEL 1.0f
#define NORM_PCM (PCM_LEVEL / 32768.0f)
#define NORM_PCM_ENERGY (NORM_PCM * NORM_PCM)
#define LOG_NORM_PCM -15

#define TNS_PREDGAIN_SCALE (1000)

#define LFE_LOWPASS_LINE 12
#define LFE_LOWPASS_LINE_MIN 4

#endif /* PSY_CONST_H */
