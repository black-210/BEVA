#ifndef BEVA_PHYSICS_H
#define BEVA_PHYSICS_H
#include <stddef.h>
#include "status.h"
#include "matrix.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { size_t n; double mass; } beva_particle;
typedef struct { size_t n; double dt, softening, G; } beva_nbody_options;
typedef struct { size_t iterations; double residual; int converged; } beva_solver_stats;
beva_status beva_nbody_step(const beva_particle*,double*,double*,const beva_nbody_options*);
beva_status beva_poisson_jacobi_2d(const double*,size_t,size_t,double,double,size_t,double,double*,beva_solver_stats*);
beva_status beva_wave_1d_step(const double*,const double*,size_t,double,double,double,double*);
beva_status beva_lorentz_gamma(double,double,double*);
beva_status beva_schrodinger_dispersion(double,double,double,double*);
beva_status beva_kg_dispersion(double,double,double*);
beva_status beva_dimensional_check(const int*,const int*,size_t,int*);
#ifdef __cplusplus
}
#endif
#endif
