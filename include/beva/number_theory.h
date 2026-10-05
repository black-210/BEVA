#ifndef BEVA_NUMBER_THEORY_H
#define BEVA_NUMBER_THEORY_H
#include <stdint.h>
#include "status.h"
uint64_t beva_gcd_u64(uint64_t a, uint64_t b);
uint64_t beva_lcm_u64(uint64_t a, uint64_t b);
int64_t beva_gcd_i64(int64_t a, int64_t b);
beva_status beva_extended_gcd_i64(int64_t a, int64_t b, int64_t *g, int64_t *x, int64_t *y);
beva_status beva_mod_inverse_i64(int64_t a, int64_t m, int64_t *out);
int beva_is_prime_u64(uint64_t n);
#endif
