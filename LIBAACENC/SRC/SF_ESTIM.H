/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Scale factor estimation

*******************************************************************************/

#ifndef SF_ESTIM_H
#define SF_ESTIM_H

#include "common_fix.h"
#include "psy_const.h"
#include "qc_data.h"
#include "interface.h"

#define FORM_FAC_SHIFT 6

void FDKaacEnc_CalcFormFactor(QC_OUT_CHANNEL *qcOutChannel[(2)],
                              PSY_OUT_CHANNEL *psyOutChannel[(2)],
                              const INT nChannels);

void FDKaacEnc_EstimateScaleFactors(PSY_OUT_CHANNEL *psyOutChannel[],
                                    QC_OUT_CHANNEL *qcOutChannel[],
                                    const INT invQuant,
                                    const INT dZoneQuantEnable,
                                    const INT nChannels);

#endif
