#include "beva/qr.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

beva_status beva_qr_mgs(const beva_matrix*a,beva_matrix*q,beva_matrix*r){
    if(!a||!q||!r)return BEVA_ERR_NULL;if(a->rows<a->cols)return BEVA_ERR_DIMENSION;
    size_t m=a->rows,n=a->cols;beva_status e=beva_matrix_init(q,m,n);if(e)return e;e=beva_matrix_init(r,n,n);if(e){beva_matrix_free(q);return e;}
    double *v=malloc(m*sizeof(double));if(!v){beva_matrix_free(q);beva_matrix_free(r);return BEVA_ERR_ALLOC;}
    for(size_t j=0;j<n;j++){for(size_t i=0;i<m;i++)v[i]=a->data[i*n+j];for(size_t k=0;k<j;k++){double d=0;for(size_t i=0;i<m;i++)d+=q->data[i*n+k]*v[i];r->data[k*n+j]=d;for(size_t i=0;i<m;i++)v[i]-=d*q->data[i*n+k];}double norm=0;for(size_t i=0;i<m;i++)norm+=v[i]*v[i];norm=sqrt(norm);if(norm<1e-14){free(v);beva_matrix_free(q);beva_matrix_free(r);return BEVA_ERR_SINGULAR;}r->data[j*n+j]=norm;for(size_t i=0;i<m;i++)q->data[i*n+j]=v[i]/norm;}free(v);return BEVA_OK;
}
beva_status beva_qr_solve(const beva_matrix*a,const beva_vector*b,beva_vector*x){
    if(!a||!b||!x)return BEVA_ERR_NULL;if(a->rows<a->cols||b->n!=a->rows)return BEVA_ERR_DIMENSION;
    beva_matrix q,r;beva_status e=beva_qr_mgs(a,&q,&r);if(e)return e;size_t n=a->cols;beva_vector y;if(beva_vector_init(&y,n)){beva_matrix_free(&q);beva_matrix_free(&r);return BEVA_ERR_ALLOC;}
    for(size_t j=0;j<n;j++)for(size_t i=0;i<a->rows;i++)y.data[j]+=q.data[i*n+j]*b->data[i];
    beva_vector sol; e=beva_vector_init(&sol,n);if(e){beva_vector_free(&y);beva_matrix_free(&q);beva_matrix_free(&r);return e;}
    for(size_t ii=0;ii<n;ii++){size_t i=n-1-ii;double s=y.data[i];for(size_t j=i+1;j<n;j++)s-=r.data[i*n+j]*sol.data[j];if(fabs(r.data[i*n+i])<1e-14){e=BEVA_ERR_SINGULAR;break;}sol.data[i]=s/r.data[i*n+i];}
    if(e==BEVA_OK){*x=sol;memset(&sol,0,sizeof(sol));}beva_vector_free(&sol);beva_vector_free(&y);beva_matrix_free(&q);beva_matrix_free(&r);return e;
}
