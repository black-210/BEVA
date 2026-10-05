#include "beva/lattice.h"
#include "beva/matrix.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>

static int64_t round_i64(double x,int *ok){if(!isfinite(x)||x>(double)INT64_MAX||x<(double)INT64_MIN){*ok=0;return 0;}*ok=1;return(int64_t)llround(x);}
static beva_status safe_row_add(beva_lattice*b,size_t dst,size_t src,int64_t k,beva_lattice*T){
    size_t n=b->cols;for(size_t j=0;j<n;j++){__int128 z=(__int128)b->data[dst*n+j]+(__int128)k*b->data[src*n+j];if(z>INT64_MAX||z<INT64_MIN)return BEVA_ERR_OVERFLOW;T->data[dst*n+j]=(int64_t)z;}return BEVA_OK;
}
beva_lll_options beva_lll_default_options(void){beva_lll_options o={0,0.99,0.5,100000,0,0};return o;}
beva_status beva_lattice_init(beva_lattice*b,size_t r,size_t c){if(!b)return BEVA_ERR_NULL;b->rows=r;b->cols=c;b->data=NULL;if(r&&c){if(r>SIZE_MAX/c)return BEVA_ERR_OVERFLOW;b->data=calloc(r*c,sizeof(int64_t));if(!b->data){b->rows=b->cols=0;return BEVA_ERR_ALLOC;}}return BEVA_OK;}
void beva_lattice_free(beva_lattice*b){if(b){free(b->data);b->data=NULL;b->rows=b->cols=0;}}
beva_status beva_lattice_copy(beva_lattice*d,const beva_lattice*s){if(!d||!s)return BEVA_ERR_NULL;beva_status e=beva_lattice_init(d,s->rows,s->cols);if(e)return e;memcpy(d->data,s->data,s->rows*s->cols*sizeof(int64_t));return BEVA_OK;}
int64_t*beva_lattice_at(beva_lattice*b,size_t r,size_t c){return(b&&r<b->rows&&c<b->cols)?&b->data[r*b->cols+c]:NULL;}
const int64_t*beva_lattice_cat(const beva_lattice*b,size_t r,size_t c){return(b&&r<b->rows&&c<b->cols)?&b->data[r*b->cols+c]:NULL;}
beva_status beva_lattice_set(beva_lattice*b,size_t r,size_t c,int64_t x){int64_t*p=beva_lattice_at(b,r,c);if(!p)return BEVA_ERR_DIMENSION;*p=x;return BEVA_OK;}
beva_status beva_lattice_get(const beva_lattice*b,size_t r,size_t c,int64_t*x){if(!x)return BEVA_ERR_NULL;const int64_t*p=beva_lattice_cat(b,r,c);if(!p)return BEVA_ERR_DIMENSION;*x=*p;return BEVA_OK;}
beva_status beva_lattice_to_double(const beva_lattice*b,beva_matrix*out){if(!b||!out)return BEVA_ERR_NULL;beva_status e=beva_matrix_init(out,b->rows,b->cols);if(e)return e;for(size_t i=0;i<b->rows*b->cols;i++)out->data[i]=(double)b->data[i];return BEVA_OK;}

beva_status beva_lattice_gso(const beva_lattice*b,beva_matrix*bs,beva_matrix*mu,double*n2){
    if(!b||!bs||!mu||!n2)return BEVA_ERR_NULL;if(!b->rows||!b->cols||b->rows>b->cols)return BEVA_ERR_DIMENSION;
    size_t m=b->rows,n=b->cols;beva_status e=beva_matrix_init(bs,m,n);if(e)return e;e=beva_matrix_init(mu,m,m);if(e){beva_matrix_free(bs);return e;}
    for(size_t i=0;i<m;i++){for(size_t k=0;k<n;k++)bs->data[i*n+k]=(double)b->data[i*n+k];for(size_t j=0;j<i;j++){double num=0,den=0;for(size_t k=0;k<n;k++){num+=(double)b->data[i*n+k]*bs->data[j*n+k];den+=bs->data[j*n+k]*bs->data[j*n+k];}if(den<=0){beva_matrix_free(bs);beva_matrix_free(mu);return BEVA_ERR_DOMAIN;}mu->data[i*m+j]=num/den;for(size_t k=0;k<n;k++)bs->data[i*n+k]-=mu->data[i*m+j]*bs->data[j*n+k];}double z=0;for(size_t k=0;k<n;k++)z+=bs->data[i*n+k]*bs->data[i*n+k];n2[i]=z;if(z<=0){beva_matrix_free(bs);beva_matrix_free(mu);return BEVA_ERR_DOMAIN;}}return BEVA_OK;
}

beva_status beva_lattice_lll(beva_lattice*b,const beva_lll_options*opt,beva_lll_stats*st,beva_lattice*T){
    if(!b)return BEVA_ERR_NULL;if(b->rows==0||b->rows>b->cols)return BEVA_ERR_DIMENSION;beva_lll_options o=opt?*opt:beva_lll_default_options();if(o.delta<=0.25||o.delta>=1||o.eta<=0||o.eta>0.5||!o.max_iterations)return BEVA_ERR_INVALID;
    size_t m=b->rows,n=b->cols;if(st)memset(st,0,sizeof(*st));if(st)st->dimension=m;
    if(o.track_transform){if(!T)return BEVA_ERR_NULL;beva_status e=beva_lattice_init(T,m,m);if(e)return e;for(size_t i=0;i<m;i++)T->data[i*m+i]=1;}
    beva_matrix bs,mu;double*n2=calloc(m,sizeof(double));if(!n2)return BEVA_ERR_ALLOC;beva_status e=beva_lattice_gso(b,&bs,&mu,n2);if(e){free(n2);return e;}
    size_t k=1,iter=0;
    while(k<m){
        if(++iter>o.max_iterations){if(st){st->iterations=iter;st->converged=0;}beva_matrix_free(&bs);beva_matrix_free(&mu);free(n2);return BEVA_ERR_CONVERGENCE;}
        for(size_t jj=k;jj>0;jj--){size_t j=jj-1;double x=mu.data[k*m+j];int ok;int64_t q=round_i64(x,&ok);if(!ok){e=BEVA_ERR_OVERFLOW;goto done;}if(llabs(q)>0){
            beva_lattice tmp={0};e=beva_lattice_init(&tmp,m,n);if(e)goto done;e=safe_row_add(b,k,j,-q,&tmp);if(e){beva_lattice_free(&tmp);goto done;}memcpy(b->data+k*n,tmp.data+n*k,n*sizeof(int64_t));beva_lattice_free(&tmp);
            if(T){beva_lattice tt={0};e=beva_lattice_init(&tt,m,m);if(e)goto done;e=safe_row_add(T,k,j,-q,&tt);if(e){beva_lattice_free(&tt);goto done;}memcpy(T->data+k*m,tt.data+m*k,m*sizeof(int64_t));beva_lattice_free(&tt);}
            if(st)st->size_reductions++;e=beva_lattice_gso(b,&bs,&mu,n2);if(e)goto done;
        }}
        double lhs=n2[k], rhs=(o.delta-mu.data[k*m+k-1]*mu.data[k*m+k-1])*n2[k-1];
        int violate=o.use_siegel ? (lhs < (o.delta)*n2[k-1]) : (lhs < rhs);
        if(violate){
            for(size_t j=0;j<n;j++){int64_t t=b->data[k*n+j];b->data[k*n+j]=b->data[(k-1)*n+j];b->data[(k-1)*n+j]=t;}
            if(T)for(size_t j=0;j<m;j++){int64_t t=T->data[k*m+j];T->data[k*m+j]=T->data[(k-1)*m+j];T->data[(k-1)*m+j]=t;}
            if(st)st->swaps++;if(k>1)k--;e=beva_lattice_gso(b,&bs,&mu,n2);if(e)goto done;
        }else k++;
    }
    e=BEVA_OK;
done:
    if(st){st->iterations=iter;st->converged=(e==BEVA_OK);}
    beva_matrix_free(&bs);beva_matrix_free(&mu);free(n2);return e;
}

beva_status beva_lattice_volume(const beva_lattice*b,double*out){
    if(!b||!out)return BEVA_ERR_NULL;if(b->rows>b->cols)return BEVA_ERR_DIMENSION;beva_matrix bs,mu;double*n=calloc(b->rows,sizeof(double));if(!n)return BEVA_ERR_ALLOC;beva_status e=beva_lattice_gso(b,&bs,&mu,n);if(e){free(n);return e;}double v=1;for(size_t i=0;i<b->rows;i++)v*=sqrt(n[i]);*out=v;beva_matrix_free(&bs);beva_matrix_free(&mu);free(n);return BEVA_OK;
}
beva_status beva_lattice_orthogonality_defect(const beva_lattice*b,double*out){
    if(!b||!out)return BEVA_ERR_NULL;double vol;beva_status e=beva_lattice_volume(b,&vol);if(e)return e;if(vol==0)return BEVA_ERR_DOMAIN;double p=1;for(size_t i=0;i<b->rows;i++){double s=0;for(size_t j=0;j<b->cols;j++){double x=(double)b->data[i*b->cols+j];s+=x*x;}p*=sqrt(s);}*out=p/vol;return BEVA_OK;
}
beva_status beva_lattice_is_reduced(const beva_lattice*b,double delta,double eta,int*out){
    if(!b||!out)return BEVA_ERR_NULL;beva_lll_options o=beva_lll_default_options();o.delta=delta;o.eta=eta;o.max_iterations=1;
    beva_matrix bs,mu;double*n=calloc(b->rows,sizeof(double));if(!n)return BEVA_ERR_ALLOC;beva_status e=beva_lattice_gso(b,&bs,&mu,n);if(e){free(n);return e;}*out=1;for(size_t i=1;i<b->rows;i++){for(size_t j=0;j<i;j++)if(fabs(mu.data[i*b->rows+j])>eta){*out=0;break;}if(n[i]<(delta-mu.data[i*b->rows+i-1]*mu.data[i*b->rows+i-1])*n[i-1]){*out=0;break;}}beva_matrix_free(&bs);beva_matrix_free(&mu);free(n);(void)o;return BEVA_OK;
}
beva_status beva_lattice_report_make(const beva_lattice*b,beva_lattice_report*r){if(!b||!r)return BEVA_ERR_NULL;beva_status e=beva_lattice_volume(b,&r->volume);if(e)return e;e=beva_lattice_orthogonality_defect(b,&r->orthogonality_defect);if(e)return e;r->stats.dimension=b->rows;r->stats.swaps=r->stats.size_reductions=r->stats.iterations=0;r->stats.converged=1;return BEVA_OK;}

beva_status beva_lattice_nearest_plane(const beva_lattice*b,const double*target,double*out,double*distance){
    if(!b||!target||!out)return BEVA_ERR_NULL;if(b->rows>b->cols)return BEVA_ERR_DIMENSION;size_t m=b->rows,n=b->cols;beva_matrix bs,mu;double*n2=calloc(m,sizeof(double));if(!n2)return BEVA_ERR_ALLOC;beva_status e=beva_lattice_gso(b,&bs,&mu,n2);if(e){free(n2);return e;}
    double* y=malloc(m*sizeof(double));if(!y){beva_matrix_free(&bs);beva_matrix_free(&mu);free(n2);return BEVA_ERR_ALLOC;}double*res=malloc(n*sizeof(double));if(!res){free(y);beva_matrix_free(&bs);beva_matrix_free(&mu);free(n2);return BEVA_ERR_ALLOC;}memcpy(res,target,n*sizeof(double));
    for(size_t ii=0;ii<m;ii++){size_t i=m-1-ii;double num=0;for(size_t j=0;j<n;j++)num+=res[j]*bs.data[i*n+j];double c=num/n2[i];y[i]=round(c);for(size_t j=0;j<n;j++)res[j]-=y[i]*b->data[i*n+j];}
    for(size_t j=0;j<n;j++)out[j]=target[j]-res[j];if(distance){double d=0;for(size_t j=0;j<n;j++){double z=target[j]-out[j];d+=z*z;}*distance=sqrt(d);}free(y);free(res);beva_matrix_free(&bs);beva_matrix_free(&mu);free(n2);return BEVA_OK;
}
beva_status beva_lattice_babai(const beva_lattice*b,const double*target,int64_t*coeff,double*distance){
    if(!b||!target||!coeff)return BEVA_ERR_NULL;if(b->rows>b->cols)return BEVA_ERR_DIMENSION;size_t m=b->rows,n=b->cols;beva_matrix bs,mu;double*n2=calloc(m,sizeof(double));if(!n2)return BEVA_ERR_ALLOC;beva_status e=beva_lattice_gso(b,&bs,&mu,n2);if(e){free(n2);return e;}double*res=malloc(n*sizeof(double));if(!res){beva_matrix_free(&bs);beva_matrix_free(&mu);free(n2);return BEVA_ERR_ALLOC;}memcpy(res,target,n*sizeof(double));for(size_t i=0;i<m;i++)coeff[i]=0;
    for(size_t ii=0;ii<m;ii++){size_t i=m-1-ii;double num=0;for(size_t j=0;j<n;j++)num+=res[j]*bs.data[i*n+j];double c=round(num/n2[i]);if(c>(double)INT64_MAX||c<(double)INT64_MIN){free(res);beva_matrix_free(&bs);beva_matrix_free(&mu);free(n2);return BEVA_ERR_OVERFLOW;}coeff[i]=(int64_t)c;for(size_t j=0;j<n;j++)res[j]-=c*(double)b->data[i*n+j];}
    if(distance){double d=0;for(size_t j=0;j<n;j++)d+=res[j]*res[j];*distance=sqrt(d);}free(res);beva_matrix_free(&bs);beva_matrix_free(&mu);free(n2);return BEVA_OK;
}
