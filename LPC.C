#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lpc.h"
#include "lpcm.h"

float vorbis_lpc_from_data(short *data,float *lpci,int n,int m,int stride){
  double *aut=malloc(sizeof(*aut)*(m+1));
  double *lpc=malloc(sizeof(*lpc)*(m));
  double error;
  double epsilon;
  int i,j;

  j=m+1;
  while(j--){
    double d=0;
    for(i=j;i<n;i++)d+=(double)data[i*stride]*data[(i-j)*stride]/1073741824.0;
    aut[j]=d;
  }

  error=aut[0] * (1. + 1e-10);
  epsilon=1e-9*aut[0]+1e-10;

  for(i=0;i<m;i++){
    double r= -aut[i+1];

    if(error<epsilon){
      memset(lpc+i,0,(m-i)*sizeof(*lpc));
      goto done;
    }

    for(j=0;j<i;j++)r-=lpc[j]*aut[i-j];
    r/=error;

    lpc[i]=r;
    for(j=0;j<i/2;j++){
      double tmp=lpc[j];

      lpc[j]+=r*lpc[i-1-j];
      lpc[i-1-j]+=r*tmp;
    }
    if(i&1)lpc[j]+=lpc[j]*r;

    error*=1.-r*r;

  }

 done:

  {
    double g = .99;
    double damp = g;
    for(j=0;j<m;j++){
      lpc[j]*=damp;
      damp*=g;
    }
  }

  for(j=0;j<m;j++)lpci[j]=(float)lpc[j];

  free(aut);
  free(lpc);
  return (float)error;
}

void vorbis_lpc_predict(float *coeff,short *prime,int m,
                        short *data,long n,int stride){

  long i,j,o,p;
  float y;
  float *work=malloc(sizeof(*work)*(m+n));

  if(!prime)
    for(i=0;i<m;i++)
      work[i]=0.f;
  else
    for(i=0;i<m;i++)
      work[i]=prime[i*stride]/32768.0f;

  for(i=0;i<n;i++){
    y=0;
    o=i;
    p=m;
    for(j=0;j<m;j++)
      y-=work[o++]*coeff[--p];

    work[o]=y;
    data[i*stride]=lrint(pcm_clip(y*32768.0,-32768.0,32767.0));
  }
  free(work);
}
