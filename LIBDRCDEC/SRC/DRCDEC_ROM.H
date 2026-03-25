/************************* MPEG-D DRC decoder library **************************

   Author(s):

   Description:

*******************************************************************************/

#ifndef DRCDEC_ROM_H
#define DRCDEC_ROM_H

extern const SCHAR deltaGain_codingProfile_0_1_huffman[24][2];
extern const SCHAR deltaGain_codingProfile_2_huffman[48][2];

extern const FIXP_SGL slopeSteepness[];
extern const SCHAR slopeSteepness_huffman[14][2];

extern const FIXP_DBL downmixCoeff[];
extern const FIXP_DBL downmixCoeffV1[];

extern const CUSTOM_DRC_CHAR_SIGMOID cicpDrcCharSigmoidLeft[];
extern const CUSTOM_DRC_CHAR_SIGMOID cicpDrcCharSigmoidRight[];
extern const CUSTOM_DRC_CHAR_NODES cicpDrcCharNodesLeft[];
extern const CUSTOM_DRC_CHAR_NODES cicpDrcCharNodesRight[];

#endif
