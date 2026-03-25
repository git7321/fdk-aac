/**************************** AAC encoder library ******************************

   Author(s):   M. Werner

   Description: Band/Line energy calculation

*******************************************************************************/

#ifndef BAND_NRG_H
#define BAND_NRG_H

#include "common_fix.h"

void FDKaacEnc_CalcSfbMaxScaleSpec(const FIXP_DBL *mdctSpectrum,
                                   const INT *bandOffset, INT *sfbMaxScaleSpec,
                                   const INT numBands);

FIXP_DBL
FDKaacEnc_CheckBandEnergyOptim(const FIXP_DBL *const RESTRICT mdctSpectrum,
                               const INT *const RESTRICT sfbMaxScaleSpec,
                               const INT *const RESTRICT bandOffset,
                               const INT numBands,
                               FIXP_DBL *RESTRICT bandEnergy,
                               FIXP_DBL *RESTRICT bandEnergyLdData,
                               const INT minSpecShift);

INT FDKaacEnc_CalcBandEnergyOptimLong(const FIXP_DBL *mdctSpectrum,
                                      INT *sfbMaxScaleSpec,
                                      const INT *bandOffset, const INT numBands,
                                      FIXP_DBL *bandEnergy,
                                      FIXP_DBL *bandEnergyLdData);

void FDKaacEnc_CalcBandEnergyOptimShort(const FIXP_DBL *mdctSpectrum,
                                        INT *sfbMaxScaleSpec,
                                        const INT *bandOffset,
                                        const INT numBands,
                                        FIXP_DBL *bandEnergy);

void FDKaacEnc_CalcBandNrgMSOpt(
    const FIXP_DBL *RESTRICT mdctSpectrumLeft,
    const FIXP_DBL *RESTRICT mdctSpectrumRight,
    INT *RESTRICT sfbMaxScaleSpecLeft, INT *RESTRICT sfbMaxScaleSpecRight,
    const INT *RESTRICT bandOffset, const INT numBands,
    FIXP_DBL *RESTRICT bandEnergyMid, FIXP_DBL *RESTRICT bandEnergySide,
    INT calcLdData, FIXP_DBL *RESTRICT bandEnergyMidLdData,
    FIXP_DBL *RESTRICT bandEnergySideLdData);

#endif
