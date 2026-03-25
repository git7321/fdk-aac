/**************************** AAC decoder library ******************************

   Author(s):   Josef Hoepfl

   Description: perceptual noise substitution tool

*******************************************************************************/

#include "aacdec_pns.h"

#include "aac_ram.h"
#include "aac_rom.h"
#include "channelinfo.h"
#include "block.h"
#include "FDK_bitstream.h"

#include "genericStds.h"

#define NOISE_OFFSET 90

void CPns_ResetData(CPnsData *pPnsData,
                    CPnsInterChannelData *pPnsInterChannelData) {
  FDK_ASSERT(pPnsData != NULL);
  FDK_ASSERT(pPnsInterChannelData != NULL);
  pPnsData->pPnsInterChannelData = pPnsInterChannelData;
  pPnsData->PnsActive = 0;
  pPnsData->CurrentEnergy = 0;

  FDKmemclear(pPnsData->pnsUsed, (8 * 16) * sizeof(UCHAR));
  FDKmemclear(pPnsInterChannelData->correlated, (8 * 16) * sizeof(UCHAR));
}

void CPns_UpdateNoiseState(CPnsData *pPnsData, INT *currentSeed,
                           INT *randomSeed) {
  pPnsData->currentSeed = currentSeed;
  pPnsData->randomSeed = randomSeed;
}

int CPns_IsPnsUsed(const CPnsData *pPnsData, const int group, const int band) {
  unsigned pns_band = group * 16 + band;

  return pPnsData->pnsUsed[pns_band] & (UCHAR)1;
}

void CPns_SetCorrelation(CPnsData *pPnsData, const int group, const int band,
                         const int outofphase) {
  CPnsInterChannelData *pInterChannelData = pPnsData->pPnsInterChannelData;
  unsigned pns_band = group * 16 + band;

  pInterChannelData->correlated[pns_band] = (outofphase) ? 3 : 1;
}

static int CPns_IsCorrelated(const CPnsData *pPnsData, const int group,
                             const int band) {
  CPnsInterChannelData *pInterChannelData = pPnsData->pPnsInterChannelData;
  unsigned pns_band = group * 16 + band;

  return (pInterChannelData->correlated[pns_band] & 0x01) ? 1 : 0;
}

static int CPns_IsOutOfPhase(const CPnsData *pPnsData, const int group,
                             const int band) {
  CPnsInterChannelData *pInterChannelData = pPnsData->pPnsInterChannelData;
  unsigned pns_band = group * 16 + band;

  return (pInterChannelData->correlated[pns_band] & 0x02) ? 1 : 0;
}

void CPns_Read(CPnsData *pPnsData, HANDLE_FDK_BITSTREAM bs,
               const CodeBookDescription *hcb, SHORT *pScaleFactor,
               UCHAR global_gain, int band, int group) {
  int delta;
  UINT pns_band = group * 16 + band;

  if (pPnsData->PnsActive) {
    delta = CBlock_DecodeHuffmanWord(bs, hcb) - 60;
  } else {
    int noiseStartValue = FDKreadBits(bs, 9);

    delta = noiseStartValue - 256;
    pPnsData->PnsActive = 1;
    pPnsData->CurrentEnergy = global_gain - NOISE_OFFSET;
  }

  pPnsData->CurrentEnergy += delta;
  pScaleFactor[pns_band] = pPnsData->CurrentEnergy;

  pPnsData->pnsUsed[pns_band] = 1;
}

static int GenerateRandomVector(FIXP_DBL *RESTRICT spec, int size,
                                int *pRandomState) {
  int i, invNrg_e = 0, nrg_e = 0;
  FIXP_DBL invNrg_m, nrg_m = FL2FXCONST_DBL(0.0f);
  FIXP_DBL *RESTRICT ptr = spec;
  int randomState = *pRandomState;

#define GEN_NOISE_NRG_SCALE 7

  for (i = 0; i < size; i++) {
    randomState =
        (((INT64)1664525 * randomState) + (INT64)1013904223) & 0xFFFFFFFF;
    nrg_m = fPow2AddDiv2(nrg_m, (FIXP_DBL)randomState >> GEN_NOISE_NRG_SCALE);
    *ptr++ = (FIXP_DBL)randomState;
  }
  nrg_e = GEN_NOISE_NRG_SCALE * 2 + 1;

  invNrg_m = invSqrtNorm2(nrg_m << 1, &invNrg_e);
  invNrg_e += -((nrg_e - 1) >> 1);

  for (i = size; i--;) {
    spec[i] = fMult(spec[i], invNrg_m);
  }

  *pRandomState = randomState;

  return invNrg_e;
}

static void ScaleBand(FIXP_DBL *RESTRICT spec, int size, int scaleFactor,
                      int specScale, int noise_e, int out_of_phase) {
  int i, shift, sfExponent;
  FIXP_DBL sfMatissa;

  sfMatissa = MantissaTable[scaleFactor & 0x03][0];
  sfExponent = (scaleFactor >> 2) + 1;

  if (out_of_phase != 0) {
    sfMatissa = -sfMatissa;
  }

  shift = sfExponent - specScale + 1 + noise_e;

  if (shift >= 0) {
    shift = fixMin(shift, DFRACT_BITS - 1);
    for (i = size; i-- != 0;) {
      spec[i] = fMultDiv2(spec[i], sfMatissa) << shift;
    }
  } else {
    shift = fixMin(-shift, DFRACT_BITS - 1);
    for (i = size; i-- != 0;) {
      spec[i] = fMultDiv2(spec[i], sfMatissa) >> shift;
    }
  }
}

void CPns_Apply(const CPnsData *pPnsData, const CIcsInfo *pIcsInfo,
                SPECTRAL_PTR pSpectrum, const SHORT *pSpecScale,
                const SHORT *pScaleFactor,
                const SamplingRateInfo *pSamplingRateInfo,
                const INT granuleLength, const int channel) {
  if (pPnsData->PnsActive) {
    const short *BandOffsets =
        GetScaleFactorBandOffsets(pIcsInfo, pSamplingRateInfo);

    int ScaleFactorBandsTransmitted = GetScaleFactorBandsTransmitted(pIcsInfo);

    for (int window = 0, group = 0; group < GetWindowGroups(pIcsInfo);
         group++) {
      for (int groupwin = 0; groupwin < GetWindowGroupLength(pIcsInfo, group);
           groupwin++, window++) {
        FIXP_DBL *spectrum = SPEC(pSpectrum, window, granuleLength);

        for (int band = 0; band < ScaleFactorBandsTransmitted; band++) {
          if (CPns_IsPnsUsed(pPnsData, group, band)) {
            UINT pns_band = window * 16 + band;

            int bandWidth = BandOffsets[band + 1] - BandOffsets[band];
            int noise_e;

            FDK_ASSERT(bandWidth >= 0);

            if (channel > 0 && CPns_IsCorrelated(pPnsData, group, band)) {
              noise_e =
                  GenerateRandomVector(spectrum + BandOffsets[band], bandWidth,
                                       &pPnsData->randomSeed[pns_band]);
            } else {
              pPnsData->randomSeed[pns_band] = *pPnsData->currentSeed;

              noise_e = GenerateRandomVector(spectrum + BandOffsets[band],
                                             bandWidth, pPnsData->currentSeed);
            }

            int outOfPhase = CPns_IsOutOfPhase(pPnsData, group, band);

            ScaleBand(spectrum + BandOffsets[band], bandWidth,
                      pScaleFactor[group * 16 + band], pSpecScale[window],
                      noise_e, outOfPhase);
          }
        }
      }
    }
  }
}
