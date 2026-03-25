/**************************** SBR decoder library ******************************

   Author(s):   Matthias Hildenbrand

   Description: Decode Predictive Vector Coding Data

*******************************************************************************/

#ifndef PVC_DEC_H
#define PVC_DEC_H

#include "common_fix.h"

#define PVC_DIVMODE_BITS 3
#define PVC_REUSEPVCID_BITS 1
#define PVC_PVCID_BITS 7
#define PVC_GRIDINFO_BITS 1

#define MAX_PVC_ENVELOPES 2
#define PVC_NTIMESLOT 16
#define PVC_NBLOW 3

#define PVC_NBHIGH_MODE1 8
#define PVC_NBHIGH_MODE2 6
#define PVC_NBHIGH_MAX (PVC_NBHIGH_MODE1)
#define PVC_NS_MAX 16

typedef struct {
  UCHAR kx_last;        /**< Xover frequency of last frame */
  UCHAR pvc_mode_last;  /**< PVC mode of last frame */
  UCHAR Esg_slot_index; /**< Ring buffer index to current Esg time slot */
  UCHAR pvcBorder0;     /**< Start SBR time slot of PVC frame */
  FIXP_DBL Esg[PVC_NS_MAX][PVC_NBLOW]; /**< Esg(ksg,t) of current and 15
                                          previous time slots (ring buffer) in
                                          logarithmical domain */
} PVC_STATIC_DATA;

typedef struct {
  UCHAR pvc_mode;   /**< PVC mode 1 or 2, 0 means legacy SBR */
  UCHAR pvcBorder0; /**< Start SBR time slot of PVC frame */
  UCHAR kx;         /**< Index of the first QMF subband in the SBR range */
  UCHAR RATE;       /**< Number of QMF subband samples per time slot (2 or 4) */
  UCHAR ns; /**< Number of time slots for time-domain smoothing of Esg(ksg,t) */
  const UCHAR
      *pPvcID; /**< Pointer to prediction coefficient matrix index table */
  UCHAR pastEsgSlotsAvail;   /**< Number of past Esg(ksg,t) which are available
                                for smoothing filter */
  const FIXP_SGL *pSCcoeffs; /**< Pointer to smoothing window table */
  SCHAR
  sg_offset_low[PVC_NBLOW + 1]; /**< Offset table for PVC grouping of SBR
                                   subbands below SBR range */
  SCHAR sg_offset_high_kx[PVC_NBHIGH_MAX + 1]; /**< Offset table for PVC
                                                  grouping of SBR subbands in
                                                  SBR range (relativ to kx) */
  UCHAR nbHigh; /**< Number of grouped QMF subbands in the SBR range */
  const SCHAR *pScalingCoef; /**< Pointer to scaling coeff table */
  const UCHAR *pPVCTab1;     /**< PVC mode 1 table */
  const UCHAR *pPVCTab2;     /**< PVC mode 2 table */
  const UCHAR *pPVCTab1_dp;  /**< Mapping of pvcID to PVC mode 1 table */
  FIXP_DBL predEsg[PVC_NTIMESLOT]
                  [PVC_NBHIGH_MAX]; /**< Predicted Energy in linear domain */
  int predEsg_exp[PVC_NTIMESLOT];   /**< Exponent of predicted Energy in linear
                                       domain */
  int predEsg_expMax;               /**< Maximum of predEsg_exp[] */
} PVC_DYNAMIC_DATA;

int pvcInitFrame(PVC_STATIC_DATA *pPvcStaticData,
                 PVC_DYNAMIC_DATA *pPvcDynamicData, const UCHAR pvcMode,
                 const UCHAR ns, const int RATE, const int kx,
                 const int pvcBorder0, const UCHAR *pPvcID);

void pvcDecodeFrame(PVC_STATIC_DATA *pPvcStaticData,
                    PVC_DYNAMIC_DATA *pPvcDynamicData, FIXP_DBL **qmfBufferReal,
                    FIXP_DBL **qmfBufferImag, const int overlap,
                    const int qmfExponentOverlap, const int qmfExponentCurrent);

void pvcDecodeTimeSlot(PVC_STATIC_DATA *pPvcStaticData,
                       PVC_DYNAMIC_DATA *pPvcDynamicData,
                       FIXP_DBL **qmfSlotReal, FIXP_DBL **qmfSlotImag,
                       const int qmfExponent, const int pvcBorder0,
                       const int timeSlotNumber, FIXP_DBL predictedEsgSlot[],
                       int *predictedEsg_exp);

void pvcEndFrame(PVC_STATIC_DATA *pPvcStaticData,
                 PVC_DYNAMIC_DATA *pPvcDynamicData);

void expandPredEsg(const PVC_DYNAMIC_DATA *pPvcDynamicData, const int timeSlot,
                   const int lengthOutputVector, FIXP_DBL *pOutput,
                   SCHAR *pOutput_exp);

#endif /* PVC_DEC_H*/
