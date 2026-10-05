#include "beva/decomposition.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
beva_status beva_cholesky(const beva_matrix*a,beva_matrix*l){
    if(!a||!l)return BEVA_ERR_NULL;if(a->rows!=a->cols)return BEVA_ERR_DIMENSION;size_t n=a->rows;beva_status e=beva_matrix_init(l,n,n);if(e)return e;
    for(size_t i=0;i<n;i++)for(size_t j=0;j<=i;j++){double s=a->data[i*n+j];for(size_t k=0;k<j;k++)s-=l->data[i*n+k]*l->data[j*n+k];if(i==j){if(s<=0){beva_matrix_free(l);return BEVA_ERR_DOMAIN;}l->data[i*n+j]=sqrt(s);}else l->data[i*n+j]=s/l->data[j*n+j];}return BEVA_OK;
}
beva_status beva_power_iteration(const beva_matrix*a,beva_vector*v,double*lambda,size_t max_iter,double tol){
    if(!a||!v||!lambda)return BEVA_ERR_NULL;if(a->rows!=a->cols||!a->rows)return BEVA_ERR_DIMENSION;size_t n=a->rows;if(max_iter==0)max_iter=1000;if(tol<=0)tol=1e-10;
    beva_status e=beva_vector_init(v,n);if(e)return e;for(size_t i=0;i<n;i++)v->data[i]=1.0/sqrt((double)n);double old=0;
    for(size_t it=0;it<max_iter;it++){beva_vector w; e=beva_matrix_vec_mul(a,v,&w);if(e){beva_vector_free(v);return e;}double nn=beva_vector_norm(&w);if(nn==0){beva_vector_free(&w);beva_vector_free(v);return BEVA_ERR_DOMAIN;}for(size_t i=0;i<n;i++)w.data[i]/=nn;double lam=beva_vector_dot(&w,&(beva_vector){n,w.data});if(fabs(lam-old)<tol){*lambda=lam;beva_vector_free(v);*v=w;return BEVA_OK;}old=lam;beva_vector_free(v);*v=w;}*lambda=old;return BEVA_ERR_CONVERGENCE;
}
