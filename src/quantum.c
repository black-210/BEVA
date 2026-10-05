#include "beva/quantum.h"
#include <math.h>
beva_complex beva_cadd(beva_complex a,beva_complex b){return(beva_complex){a.real+b.real,a.imag+b.imag};}
beva_complex beva_cmul(beva_complex a,beva_complex b){return(beva_complex){a.real*b.real-a.imag*b.imag,a.real*b.imag+a.imag*b.real};}
double beva_cabs2(beva_complex a){return a.real*a.real+a.imag*a.imag;}
beva_status beva_quantum_normalize(beva_complex*p,size_t n){if(!p)return BEVA_ERR_NULL;double s=0;for(size_t i=0;i<n;i++)s+=beva_cabs2(p[i]);if(s<=0)return BEVA_ERR_DOMAIN;s=sqrt(s);for(size_t i=0;i<n;i++){p[i].real/=s;p[i].imag/=s;}return BEVA_OK;}
beva_status beva_quantum_expectation(size_t n,const beva_complex*o,const beva_complex*p,double*v){if(!o||!p||!v)return BEVA_ERR_NULL;double s=0;for(size_t i=0;i<n;i++)for(size_t j=0;j<n;j++){beva_complex z=beva_cmul(o[i*n+j],p[j]);s+=p[i].real*z.real+p[i].imag*z.imag;}*v=s;return BEVA_OK;}
