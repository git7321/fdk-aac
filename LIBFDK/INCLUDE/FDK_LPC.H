/******************* Library for basic calculation routines ********************

   Author(s):   Manuel Jander

   Description: LPC related functions

*******************************************************************************/

#ifndef FDK_LPC_H
#define FDK_LPC_H

#include "common_fix.h"

#define LPC_MAX_ORDER 24

typedef FIXP_DBL FIXP_LPC_TNS;
#define FX_DBL2FX_LPC_TNS(x) (x)
#define FX_DBL2FXCONST_LPC_TNS(x) (x)
#define FX_LPC_TNS2FX_DBL(x) (x)
#define FL2FXCONST_LPC_TNS(val) FL2FXCONST_DBL(val)
#define MAXVAL_LPC_TNS MAXVAL_DBL

typedef FIXP_SGL FIXP_LPC;
#define FX_DBL2FX_LPC(x) FX_DBL2FX_SGL((FIXP_DBL)(x))
#define FX_DBL2FXCONST_LPC(x) FX_DBL2FXCONST_SGL(x)
#define FX_LPC2FX_DBL(x) FX_SGL2FX_DBL(x)
#define FL2FXCONST_LPC(val) FL2FXCONST_SGL(val)
#define MAXVAL_LPC MAXVAL_SGL

void CLpc_Analysis(FIXP_DBL signal[], const int signal_size,
                   const FIXP_LPC lpcCoeff_m[], const int lpcCoeff_e,
                   const int order, FIXP_DBL *filtState, int *filtStateIndex);

void CLpc_Synthesis(FIXP_DBL *signal, const int signal_size, const int signal_e,
                    const int inc, const FIXP_LPC_TNS *lpcCoeff_m,
                    const int lpcCoeff_e, const int order, FIXP_DBL *state,
                    int *pStateIndex);
void CLpc_Synthesis(FIXP_DBL *signal, const int signal_size, const int signal_e,
                    const int inc, const FIXP_LPC coeff[], const int coeff_e,
                    const int order, FIXP_DBL *filtState, int *pStateIndex);

void CLpc_SynthesisLattice(FIXP_DBL *signal, const int signal_size,
                           const int signal_e, const int signal_e_out,
                           const int inc, const FIXP_SGL *coeff,
                           const int order, FIXP_DBL *state);

void CLpc_SynthesisLattice(FIXP_DBL *RESTRICT signal, const int signal_size,
                           const int signal_e, const int signal_e_out,
                           const int inc, const FIXP_DBL *RESTRICT coeff,
                           const int order, FIXP_DBL *RESTRICT state);

INT CLpc_ParcorToLpc(const FIXP_LPC_TNS reflCoeff[], FIXP_LPC_TNS LpcCoeff[],
                     INT numOfCoeff, FIXP_DBL workBuffer[]);
INT CLpc_ParcorToLpc(const FIXP_LPC reflCoeff[], FIXP_LPC LpcCoeff[],
                     const int numOfCoeff, FIXP_DBL workBuffer[]);

void CLpc_AutoToParcor(FIXP_DBL acorr[], const int acorr_e,
                       FIXP_LPC reflCoeff[], const int order,
                       FIXP_DBL *pPredictionGain_m, INT *pPredictionGain_e);

#endif /* FDK_LPC_H */
