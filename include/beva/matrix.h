#ifndef BEVA_MATRIX_H
#define BEVA_MATRIX_H
#include <stddef.h>
#include "status.h"
#include "vector.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { size_t rows, cols; double *data; } beva_matrix;
beva_status beva_matrix_init(beva_matrix*, size_t, size_t);
void beva_matrix_free(beva_matrix*);
beva_status beva_matrix_copy(beva_matrix*, const beva_matrix*);
beva_status beva_matrix_identity(beva_matrix*, size_t);
double *beva_matrix_at(beva_matrix*, size_t, size_t);
const double *beva_matrix_cat(const beva_matrix*, size_t, size_t);
beva_status beva_matrix_set(beva_matrix*, size_t, size_t, double);
beva_status beva_matrix_get(const beva_matrix*, size_t, size_t, double*);
beva_status beva_matrix_mul(const beva_matrix*, const beva_matrix*, beva_matrix*);
beva_status beva_matrix_transpose(const beva_matrix*, beva_matrix*);
beva_status beva_matrix_vec_mul(const beva_matrix*, const beva_vector*, beva_vector*);
beva_status beva_matrix_solve(const beva_matrix*, const beva_vector*, beva_vector*);
beva_status beva_matrix_inverse(const beva_matrix*, beva_matrix*);
beva_status beva_matrix_qr(const beva_matrix*, beva_matrix*, beva_matrix*);
double beva_matrix_det(const beva_matrix*);
size_t beva_matrix_rank(const beva_matrix*, double tol);
#ifdef __cplusplus
}
#endif
#endif
