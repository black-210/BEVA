#ifndef BEVA_DECOMPOSITION_H
#define BEVA_DECOMPOSITION_H
#include "matrix.h"
beva_status beva_cholesky(const beva_matrix*, beva_matrix*);
beva_status beva_power_iteration(const beva_matrix*, beva_vector*, double*, size_t, double);
#endif
