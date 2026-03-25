/*********************** MPEG surround encoder library *************************

   Author(s):   Max Neuendorf

   Description: Encoder Library Interface
                Get windows for framing

*******************************************************************************/

#ifndef SACENC_FRAMEWINDOWING_H
#define SACENC_FRAMEWINDOWING_H

#include "genericStds.h"
#include "common_fix.h"
#include "sacenc_lib.h"
#include "sacenc_bitstream.h"

#define FIXP_WIN FIXP_DBL
#define FX_DBL2FX_WIN(x) (x)
#define DALDATATYPE_WIN DALDATATYPE_DFRACT

typedef enum {
  FW_INTP = 0,
  FW_HOLD = 1

} FW_SLOTTYPE;

typedef enum {
  FW_LEAVE_DIM = 0,
  FW_CHANGE_DIM = 1

} FW_DIMENSION;

typedef struct T_FRAMEWINDOW *HANDLE_FRAMEWINDOW;

typedef struct T_FRAMEWINDOW_CONFIG {
  INT nTimeSlotsMax;
  INT bFrameKeep;

} FRAMEWINDOW_CONFIG;

typedef struct {
  INT slot;
  FW_SLOTTYPE hold;

} FRAMEWIN_DATA;

typedef struct {
  FRAMEWIN_DATA dat[MAX_NUM_PARAMS];
  INT n;

} FRAMEWIN_LIST;

FDK_SACENC_ERROR fdk_sacenc_frameWindow_Create(
    HANDLE_FRAMEWINDOW *phFrameWindow);

FDK_SACENC_ERROR fdk_sacenc_frameWindow_Init(
    HANDLE_FRAMEWINDOW hFrameWindow,
    const FRAMEWINDOW_CONFIG *const pFrameWindowConfig);

FDK_SACENC_ERROR fdk_sacenc_frameWindow_Destroy(
    HANDLE_FRAMEWINDOW *phFrameWindow);

FDK_SACENC_ERROR fdk_sacenc_frameWindow_GetWindow(
    HANDLE_FRAMEWINDOW hFrameWindow, INT tr_pos[MAX_NUM_PARAMS],
    const INT timeSlots, FRAMINGINFO *const pFramingInfo,
    FIXP_WIN *pWindowAna__FDK[MAX_NUM_PARAMS],
    FRAMEWIN_LIST *const pFrameWinList, const INT avoid_keep);

FDK_SACENC_ERROR fdk_sacenc_analysisWindowing(
    const INT nTimeSlots, const INT startTimeSlot,
    FIXP_WIN *pFrameWindowAna__FDK, const FIXP_DPK *const *const ppDataIn__FDK,
    FIXP_DPK *const *const ppDataOut__FDK, const INT nHybridBands,
    const INT dim);

#endif /* SACENC_FRAMEWINDOWING_H */
