#ifndef _V_LPC_H_
#define _V_LPC_H_

extern float vorbis_lpc_from_data(short *data,float *lpc,int n,int m,int stride);

extern void vorbis_lpc_predict(float *coeff,short *prime,int m,
                               short *data,long n,int stride);

#endif
