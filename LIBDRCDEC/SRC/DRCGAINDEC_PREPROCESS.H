/************************* MPEG-D DRC decoder library **************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef DRCGAINDEC_PREPROCESS_H
#define DRCGAINDEC_PREPROCESS_H

DRC_ERROR
prepareDrcGain(HANDLE_DRC_GAIN_DECODER hGainDec,
               HANDLE_UNI_DRC_GAIN hUniDrcGain, const FIXP_SGL compress,
               const FIXP_SGL boost, const FIXP_DBL loudnessNormalizationGainDb,
               const int activeDrcIndex);
#endif
