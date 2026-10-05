#ifndef BEVA_QUANTUM_H
#define BEVA_QUANTUM_H
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {double real,imag;} beva_complex;
beva_complex beva_cadd(beva_complex,beva_complex);
beva_complex beva_cmul(beva_complex,beva_complex);
double beva_cabs2(beva_complex);
beva_status beva_quantum_normalize(beva_complex*,size_t);
beva_status beva_quantum_expectation(size_t,const beva_complex*,const beva_complex*,double*);
#ifdef __cplusplus
}
#endif
#endif
