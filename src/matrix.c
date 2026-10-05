#include "beva/matrix.h"
#include "beva/qr.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <float.h>

static int valid(const beva_matrix*m){return m&&m->rows<=SIZE_MAX/m->cols;}
beva_status beva_matrix_init(beva_matrix*m,size_t r,size_t c){
    if(!m)return BEVA_ERR_NULL;m->rows=r;m->cols=c;m->data=NULL;
    if(r&&c){if(r>SIZE_MAX/c)return BEVA_ERR_OVERFLOW;m->data=calloc(r*c,sizeof(double));if(!m->data){m->rows=m->cols=0;return BEVA_ERR_ALLOC;}}
    return BEVA_OK;
}
void beva_matrix_free(beva_matrix*m){if(m){free(m->data);m->data=NULL;m->rows=m->cols=0;}}
beva_status beva_matrix_copy(beva_matrix*d,const beva_matrix*s){
    if(!d||!s)return BEVA_ERR_NULL;beva_status e=beva_matrix_init(d,s->rows,s->cols);if(e)return e;
    if(s->data)memcpy(d->data,s->data,s->rows*s->cols*sizeof(double));return BEVA_OK;
}
beva_status beva_matrix_identity(beva_matrix*m,size_t n){if(!m)return BEVA_ERR_NULL;beva_status e=beva_matrix_init(m,n,n);if(e)return e;for(size_t i=0;i<n;i++)m->data[i*n+i]=1;return BEVA_OK;}
double*beva_matrix_at(beva_matrix*m,size_t r,size_t c){return(m&&r<m->rows&&c<m->cols)?&m->data[r*m->cols+c]:NULL;}
const double*beva_matrix_cat(const beva_matrix*m,size_t r,size_t c){return(m&&r<m->rows&&c<m->cols)?&m->data[r*m->cols+c]:NULL;}
beva_status beva_matrix_set(beva_matrix*m,size_t r,size_t c,double x){double*p=beva_matrix_at(m,r,c);if(!p)return BEVA_ERR_DIMENSION;*p=x;return BEVA_OK;}
beva_status beva_matrix_get(const beva_matrix*m,size_t r,size_t c,double*x){if(!x)return BEVA_ERR_NULL;const double*p=beva_matrix_cat(m,r,c);if(!p)return BEVA_ERR_DIMENSION;*x=*p;return BEVA_OK;}
beva_status beva_matrix_mul(const beva_matrix*a,const beva_matrix*b,beva_matrix*out){
    if(!a||!b||!out)return BEVA_ERR_NULL;if(a->cols!=b->rows)return BEVA_ERR_DIMENSION;beva_status e=beva_matrix_init(out,a->rows,b->cols);if(e)return e;
    for(size_t i=0;i<a->rows;i++)for(size_t k=0;k<a->cols;k++){double x=a->data[i*a->cols+k];for(size_t j=0;j<b->cols;j++)out->data[i*out->cols+j]+=x*b->data[k*b->cols+j];}return BEVA_OK;
}
beva_status beva_matrix_transpose(const beva_matrix*a,beva_matrix*out){if(!a||!out)return BEVA_ERR_NULL;beva_status e=beva_matrix_init(out,a->cols,a->rows);if(e)return e;for(size_t i=0;i<a->rows;i++)for(size_t j=0;j<a->cols;j++)out->data[j*out->cols+i]=a->data[i*a->cols+j];return BEVA_OK;}
beva_status beva_matrix_vec_mul(const beva_matrix*a,const beva_vector*x,beva_vector*y){if(!a||!x||!y)return BEVA_ERR_NULL;if(a->cols!=x->n)return BEVA_ERR_DIMENSION;beva_status e=beva_vector_init(y,a->rows);if(e)return e;for(size_t i=0;i<a->rows;i++)for(size_t j=0;j<a->cols;j++)y->data[i]+=a->data[i*a->cols+j]*x->data[j];return BEVA_OK;}
static double tol_for(const beva_matrix*a){double mx=0;for(size_t i=0;i<a->rows*a->cols;i++)mx=fmax(mx,fabs(a->data[i]));return fmax(1e-14,mx*1e-13);}
beva_status beva_matrix_solve(const beva_matrix*a,const beva_vector*b,beva_vector*x){
    if(!a||!b||!x)return BEVA_ERR_NULL;if(a->rows!=a->cols||b->n!=a->rows)return BEVA_ERR_DIMENSION;
    size_t n=a->rows;double *m=malloc(n*(n+1)*sizeof(double));if(!m)return BEVA_ERR_ALLOC;
    for(size_t i=0;i<n;i++){for(size_t j=0;j<n;j++)m[i*(n+1)+j]=a->data[i*n+j];m[i*(n+1)+n]=b->data[i];}
    double tol=tol_for(a);for(size_t k=0;k<n;k++){size_t p=k;for(size_t i=k+1;i<n;i++)if(fabs(m[i*(n+1)+k])>fabs(m[p*(n+1)+k]))p=i;if(fabs(m[p*(n+1)+k])<=tol){free(m);return BEVA_ERR_SINGULAR;}if(p!=k)for(size_t j=k;j<=n;j++){double t=m[k*(n+1)+j];m[k*(n+1)+j]=m[p*(n+1)+j];m[p*(n+1)+j]=t;}for(size_t i=k+1;i<n;i++){double f=m[i*(n+1)+k]/m[k*(n+1)+k];for(size_t j=k;j<=n;j++)m[i*(n+1)+j]-=f*m[k*(n+1)+j];}}
    beva_status e=beva_vector_init(x,n);if(e){free(m);return e;}for(size_t ii=0;ii<n;ii++){size_t i=n-1-ii;double s=m[i*(n+1)+n];for(size_t j=i+1;j<n;j++)s-=m[i*(n+1)+j]*x->data[j];x->data[i]=s/m[i*(n+1)+i];}free(m);return BEVA_OK;
}
beva_status beva_matrix_inverse(const beva_matrix*a,beva_matrix*out){
    if(!a||!out)return BEVA_ERR_NULL;if(a->rows!=a->cols)return BEVA_ERR_DIMENSION;size_t n=a->rows;beva_status e=beva_matrix_identity(out,n);if(e)return e;double *m=malloc(n*n*sizeof(double));if(!m){beva_matrix_free(out);return BEVA_ERR_ALLOC;}memcpy(m,a->data,n*n*sizeof(double));double tol=tol_for(a);
    for(size_t k=0;k<n;k++){size_t p=k;for(size_t i=k+1;i<n;i++)if(fabs(m[i*n+k])>fabs(m[p*n+k]))p=i;if(fabs(m[p*n+k])<=tol){free(m);beva_matrix_free(out);return BEVA_ERR_SINGULAR;}if(p!=k){for(size_t j=0;j<n;j++){double t=m[k*n+j];m[k*n+j]=m[p*n+j];m[p*n+j]=t;t=out->data[k*n+j];out->data[k*n+j]=out->data[p*n+j];out->data[p*n+j]=t;}}double d=m[k*n+k];for(size_t j=0;j<n;j++){m[k*n+j]/=d;out->data[k*n+j]/=d;}for(size_t i=0;i<n;i++)if(i!=k){double f=m[i*n+k];for(size_t j=0;j<n;j++){m[i*n+j]-=f*m[k*n+j];out->data[i*n+j]-=f*out->data[k*n+j];}}}free(m);return BEVA_OK;
}
double beva_matrix_det(const beva_matrix*a){if(!a||a->rows!=a->cols)return NAN;size_t n=a->rows;if(!n)return 1;double*m=malloc(n*n*sizeof(double));if(!m)return NAN;memcpy(m,a->data,n*n*sizeof(double));double det=1;for(size_t k=0;k<n;k++){size_t p=k;for(size_t i=k+1;i<n;i++)if(fabs(m[i*n+k])>fabs(m[p*n+k]))p=i;if(fabs(m[p*n+k])<tol_for(a)){free(m);return 0;}if(p!=k){for(size_t j=0;j<n;j++){double t=m[k*n+j];m[k*n+j]=m[p*n+j];m[p*n+j]=t;}det=-det;}det*=m[k*n+k];for(size_t i=k+1;i<n;i++){double f=m[i*n+k]/m[k*n+k];for(size_t j=k+1;j<n;j++)m[i*n+j]-=f*m[k*n+j];}}free(m);return det;}
size_t beva_matrix_rank(const beva_matrix*a,double tol){if(!a)return 0;beva_matrix m;if(beva_matrix_copy(&m,a))return 0;size_t r=0;double t=tol>0?tol:tol_for(a);for(size_t c=0;c<m.cols&&r<m.rows;c++){size_t p=r;for(size_t i=r+1;i<m.rows;i++)if(fabs(m.data[i*m.cols+c])>fabs(m.data[p*m.cols+c]))p=i;if(fabs(m.data[p*m.cols+c])<=t)continue;for(size_t j=c;j<m.cols;j++){double z=m.data[r*m.cols+j];m.data[r*m.cols+j]=m.data[p*m.cols+j];m.data[p*m.cols+j]=z;}for(size_t i=r+1;i<m.rows;i++){double f=m.data[i*m.cols+c]/m.data[r*m.cols+c];for(size_t j=c;j<m.cols;j++)m.data[i*m.cols+j]-=f*m.data[r*m.cols+j];}r++;}beva_matrix_free(&m);return r;}
beva_status beva_matrix_qr(const beva_matrix*a,beva_matrix*q,beva_matrix*r){return beva_qr_mgs(a,q,r);}
