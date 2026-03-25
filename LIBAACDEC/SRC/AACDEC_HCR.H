/**************************** AAC decoder library ******************************

   Author(s):   Robert Weidner (DSP Solutions)

   Description: HCR Decoder: Interface function declaration; common defines
                and structures; defines for switching error-generator,
                -detector, and -concealment

*******************************************************************************/

#ifndef AACDEC_HCR_H
#define AACDEC_HCR_H

#include "channelinfo.h"
#include "FDK_bitstream.h"

UINT HcrInit(H_HCR_INFO pHcr, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
             const SamplingRateInfo *pSamplingRateInfo,
             HANDLE_FDK_BITSTREAM bs);
UINT HcrDecoder(H_HCR_INFO hHcr, CAacDecoderChannelInfo *pAacDecoderChannelInfo,
                const SamplingRateInfo *pSamplingRateInfo,
                HANDLE_FDK_BITSTREAM bs);
void CarryBitToBranchValue(UCHAR carryBit, UINT treeNode, UINT *branchValue,
                           UINT *branchNode);

void CHcr_Read(HANDLE_FDK_BITSTREAM bs,
               CAacDecoderChannelInfo *pAacDecoderChannelInfo,
               const MP4_ELEMENT_ID globalHcrType);
void HcrMuteErroneousLines(H_HCR_INFO hHcr);

void setHcrType(H_HCR_INFO hHcr, MP4_ELEMENT_ID type);
INT getHcrType(H_HCR_INFO hHcr);

#endif /* AACDEC_HCR_H */
