#ifndef BEVA_VECTOR_H
#define BEVA_VECTOR_H
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { size_t n; double *data; } beva_vector;
beva_status beva_vector_init(beva_vector*, size_t);
void beva_vector_free(beva_vector*);
beva_status beva_vector_copy(beva_vector*, const beva_vector*);
beva_status beva_vector_set(beva_vector*, size_t, double);
beva_status beva_vector_get(const beva_vector*, size_t, double*);
double beva_vector_dot(const beva_vector*, const beva_vector*);
double beva_vector_norm2(const beva_vector*);
double beva_vector_norm(const beva_vector*);
beva_status beva_vector_scale(beva_vector*, double);
beva_status beva_vector_sub_scaled(beva_vector*, const beva_vector*, double);
beva_status beva_vector_gram_schmidt(const beva_vector*, beva_vector*, double*);
#ifdef __cplusplus
}
#endif
#endif
