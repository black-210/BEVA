#include "beva/beva.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
int main(void){
    beva_matrix a; assert(beva_matrix_init(&a,2,2)==BEVA_OK);
    a.data[0]=4;a.data[1]=7;a.data[2]=2;a.data[3]=6;
    assert(fabs(beva_matrix_det(&a)-10)<1e-9);
    beva_matrix inv; assert(beva_matrix_inverse(&a,&inv)==BEVA_OK);
    assert(fabs(inv.data[0]-.6)<1e-9);
    beva_vector b,x; assert(beva_vector_init(&b,2)==BEVA_OK);b.data[0]=1;b.data[1]=0;
    assert(beva_matrix_solve(&a,&b,&x)==BEVA_OK);assert(fabs(x.data[0]-.6)<1e-9);
    beva_lattice l;assert(beva_lattice_init(&l,2,2)==BEVA_OK);
    l.data[0]=105;l.data[1]=821;l.data[2]=112;l.data[3]=941;
    beva_lll_options o=beva_lll_default_options();beva_lll_stats s;
    assert(beva_lattice_lll(&l,&o,&s,NULL)==BEVA_OK);
    int reduced=0;assert(beva_lattice_is_reduced(&l,.75,.5,&reduced)==BEVA_OK);assert(reduced);
    double v;assert(beva_lattice_volume(&l,&v)==BEVA_OK);assert(fabs(v-6853)<1e-6);
    beva_lattice_free(&l);beva_vector_free(&b);beva_vector_free(&x);beva_matrix_free(&inv);beva_matrix_free(&a);
    puts("BEVA tests: OK");return 0;
}
