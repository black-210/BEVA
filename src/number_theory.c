#include "beva/number_theory.h"
#include <stdlib.h>
#include <limits.h>
uint64_t beva_gcd_u64(uint64_t a,uint64_t b){while(b){uint64_t t=a%b;a=b;b=t;}return a;}
uint64_t beva_lcm_u64(uint64_t a,uint64_t b){if(!a||!b)return 0;return(a/beva_gcd_u64(a,b))*b;}
int64_t beva_gcd_i64(int64_t a,int64_t b){uint64_t x=a<0?(uint64_t)(-(a+1))+1:(uint64_t)a,y=b<0?(uint64_t)(-(b+1))+1:(uint64_t)b;uint64_t g=beva_gcd_u64(x,y);return g>INT64_MAX?INT64_MIN:(int64_t)g;}
beva_status beva_extended_gcd_i64(int64_t a,int64_t b,int64_t*g,int64_t*x,int64_t*y){
    if(!g||!x||!y)return BEVA_ERR_NULL;__int128 old_r=a,r=b,old_s=1,s=0,old_t=0,t=1;
    while(r){__int128 q=old_r/r;__int128 z=old_r-q*r;old_r=r;r=z;z=old_s-q*s;old_s=s;s=z;z=old_t-q*t;old_t=t;t=z;}
    if(old_r<0){old_r=-old_r;old_s=-old_s;old_t=-old_t;}if(old_r>INT64_MAX||old_s>INT64_MAX||old_s<INT64_MIN||old_t>INT64_MAX||old_t<INT64_MIN)return BEVA_ERR_OVERFLOW;
    *g=(int64_t)old_r;*x=(int64_t)old_s;*y=(int64_t)old_t;return BEVA_OK;
}
beva_status beva_mod_inverse_i64(int64_t a,int64_t m,int64_t*out){if(!out||m<=0)return BEVA_ERR_INVALID;int64_t g,x,y;beva_status e=beva_extended_gcd_i64(a,m,&g,&x,&y);if(e)return e;if(g!=1)return BEVA_ERR_DOMAIN;int64_t z=x%m;if(z<0)z+=m;*out=z;return BEVA_OK;}
int beva_is_prime_u64(uint64_t n){if(n<2)return 0;if(n%2==0)return n==2;if(n%3==0)return n==3;for(uint64_t d=5;d<=n/d;d+=6)if(n%d==0||n%(d+2)==0)return 0;return 1;}
