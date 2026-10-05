#include <stdio.h>
#include "beva/beva.h"
int main(void){
    beva_lattice b; if(beva_lattice_init(&b,2,2)!=BEVA_OK)return 1;
    b.data[0]=105;b.data[1]=821;b.data[2]=112;b.data[3]=941;
    beva_lll_options o=beva_lll_default_options();o.delta=.99;
    beva_lll_stats s;
    beva_status e=beva_lattice_lll(&b,&o,&s,NULL);
    if(e){fprintf(stderr,"LLL: %s\n",beva_status_string(e));return 2;}
    printf("reduced basis:\n");
    for(size_t i=0;i<b.rows;i++)printf("[%lld %lld]\n",(long long)b.data[i*b.cols],(long long)b.data[i*b.cols+1]);
    printf("iterations=%zu swaps=%zu reductions=%zu\n",s.iterations,s.swaps,s.size_reductions);
    beva_lattice_free(&b);return 0;
}
