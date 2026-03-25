/*********************** MPEG surround encoder library *************************

   Author(s):   M. Luis Valero

   Description: Enhanced Time Domain Downmix

*******************************************************************************/

#ifndef SACENC_DMX_TDOM_ENH_H
#define SACENC_DMX_TDOM_ENH_H

/* Includes ******************************************************************/
#include "sacenc_lib.h"
#include "common_fix.h"

typedef struct T_ENHANCED_TIME_DOMAIN_DMX *HANDLE_ENHANCED_TIME_DOMAIN_DMX;

FDK_SACENC_ERROR fdk_sacenc_open_enhancedTimeDomainDmx(
    HANDLE_ENHANCED_TIME_DOMAIN_DMX *hEnhancedTimeDmx, const INT framelength);

FDK_SACENC_ERROR fdk_sacenc_init_enhancedTimeDomainDmx(
    HANDLE_ENHANCED_TIME_DOMAIN_DMX hEnhancedTimeDmx,
    const FIXP_DBL *const pInputGain_m, const INT inputGain_e,
    const FIXP_DBL outputGain_m, const INT outputGain_e, const INT framelength);

FDK_SACENC_ERROR fdk_sacenc_apply_enhancedTimeDomainDmx(
    HANDLE_ENHANCED_TIME_DOMAIN_DMX hEnhancedTimeDmx,
    const INT_PCM *const *const inputTime, INT_PCM *const outputTimeDmx,
    const INT InputDelay);

FDK_SACENC_ERROR fdk_sacenc_close_enhancedTimeDomainDmx(
    HANDLE_ENHANCED_TIME_DOMAIN_DMX *hEnhancedTimeDmx);

#endif /* SACENC_DMX_TDOM_ENH_H */
