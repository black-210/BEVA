#include "beva/vector.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

beva_status beva_vector_init(beva_vector *v, size_t n) {
    if(!v) return BEVA_ERR_NULL;
    v->n=n; v->data=NULL;
    if(n) { v->data=(double*)calloc(n,sizeof(double)); if(!v->data){v->n=0;return BEVA_ERR_ALLOC;} }
    return BEVA_OK;
}
void beva_vector_free(beva_vector *v){ if(v){free(v->data);v->data=NULL;v->n=0;} }
beva_status beva_vector_copy(beva_vector *d,const beva_vector*s){
    if(!d||!s)return BEVA_ERR_NULL; beva_status e=beva_vector_init(d,s->n); if(e)return e;
    if(s->n)memcpy(d->data,s->data,s->n*sizeof(double)); return BEVA_OK;
}
beva_status beva_vector_set(beva_vector*v,size_t i,double x){if(!v||i>=v->n)return BEVA_ERR_DIMENSION;v->data[i]=x;return BEVA_OK;}
beva_status beva_vector_get(const beva_vector*v,size_t i,double*x){if(!v||!x)return BEVA_ERR_NULL;if(i>=v->n)return BEVA_ERR_DIMENSION;*x=v->data[i];return BEVA_OK;}
double beva_vector_dot(const beva_vector*a,const beva_vector*b){
    if(!a||!b||a->n!=b->n)return NAN; double s=0,c=0;
    for(size_t i=0;i<a->n;i++){double y=a->data[i]*b->data[i]-c,t=s+y;c=(t-s)-y;s=t;} return s;
}
double beva_vector_norm2(const beva_vector*v){return v?beva_vector_dot(v,v):NAN;}
double beva_vector_norm(const beva_vector*v){double x=beva_vector_norm2(v);return sqrt(x);}
beva_status beva_vector_scale(beva_vector*v,double a){if(!v)return BEVA_ERR_NULL;for(size_t i=0;i<v->n;i++)v->data[i]*=a;return BEVA_OK;}
beva_status beva_vector_sub_scaled(beva_vector*v,const beva_vector*u,double a){
    if(!v||!u)return BEVA_ERR_NULL;if(v->n!=u->n)return BEVA_ERR_DIMENSION;
    for(size_t i=0;i<v->n;i++)v->data[i]-=a*u->data[i];return BEVA_OK;
}
beva_status beva_vector_gram_schmidt(const beva_vector*in,beva_vector*out,double*mu){
    if(!in||!out)return BEVA_ERR_NULL;if(in->n!=out->n)return BEVA_ERR_DIMENSION;
    size_t n=in->n; double *q=out->data; memset(q,0,n*sizeof(double)); if(mu)memset(mu,0,n*sizeof(double));
    for(size_t i=0;i<n;i++){
        q[i]=in->data[i];
        for(size_t j=0;j<i;j++){
            double den=0,num=0;
            for(size_t k=0;k<n;k++){den+=q[j*n+k]*q[j*n+k];num+=in->data[i]*q[j*n+k];}
            double m=(den>0?num/den:0); if(mu)mu[i*n+j]=m;
            for(size_t k=0;k<n;k++)q[i*n+k]-=m*q[j*n+k];
        }
    }
    return BEVA_OK;
}
