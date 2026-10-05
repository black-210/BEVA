#ifndef BEVA_LATTICE_H
#define BEVA_LATTICE_H
#include <stddef.h>
#include <stdint.h>
#include "status.h"
#include "matrix.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { size_t rows, cols; int64_t *data; } beva_lattice;

typedef struct {
    size_t dimension;
    size_t swaps;
    size_t size_reductions;
    size_t iterations;
    int converged;
} beva_lll_stats;

typedef struct {
    size_t dimension;
    double delta;
    double eta;
    size_t max_iterations;
    int use_siegel;
    int track_transform;
} beva_lll_options;

typedef struct {
    beva_lll_stats stats;
    double volume;
    double orthogonality_defect;
} beva_lattice_report;

beva_lll_options beva_lll_default_options(void);
beva_status beva_lattice_init(beva_lattice*, size_t, size_t);
void beva_lattice_free(beva_lattice*);
beva_status beva_lattice_copy(beva_lattice*, const beva_lattice*);
beva_status beva_lattice_set(beva_lattice*, size_t, size_t, int64_t);
beva_status beva_lattice_get(const beva_lattice*, size_t, size_t, int64_t*);
int64_t *beva_lattice_at(beva_lattice*, size_t, size_t);
const int64_t *beva_lattice_cat(const beva_lattice*, size_t, size_t);

beva_status beva_lattice_to_double(const beva_lattice*, beva_matrix*);
beva_status beva_lattice_gso(const beva_lattice*, beva_matrix *bstar, beva_matrix *mu, double *norm2);
beva_status beva_lattice_lll(beva_lattice*, const beva_lll_options*, beva_lll_stats*, beva_lattice *transform);
beva_status beva_lattice_report_make(const beva_lattice*, beva_lattice_report*);
beva_status beva_lattice_is_reduced(const beva_lattice*, double delta, double eta, int *out);
beva_status beva_lattice_volume(const beva_lattice*, double*);
beva_status beva_lattice_orthogonality_defect(const beva_lattice*, double*);
beva_status beva_lattice_nearest_plane(const beva_lattice*, const double*, double*, double*);
beva_status beva_lattice_babai(const beva_lattice*, const double*, int64_t*, double*);
#ifdef __cplusplus
}
#endif
#endif
