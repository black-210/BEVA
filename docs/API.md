# BEVA API

## Core

`include/beva/beva.h` is the umbrella header.

### Vector
`beva_vector_dot`, `beva_vector_norm`, `beva_vector_scale`, `beva_vector_sub_scaled`.

### Matrix
`beva_matrix_mul`, `beva_matrix_solve`, `beva_matrix_inverse`, `beva_matrix_det`, `beva_matrix_rank`.

### QR
`beva_qr_mgs`, `beva_qr_solve`.

### Decomposition
`beva_cholesky`, `beva_power_iteration`.

### Number theory
`beva_gcd_u64`, `beva_lcm_u64`, `beva_extended_gcd_i64`, `beva_mod_inverse_i64`, `beva_is_prime_u64`.

### Lattice
`beva_lattice_gso`, `beva_lattice_lll`, `beva_lattice_volume`,
`beva_lattice_orthogonality_defect`, `beva_lattice_nearest_plane`,
`beva_lattice_babai`.

## LLL options

- `delta`: Lovász parameter, normally `0.75..0.999`
- `eta`: size reduction bound, normally just above `0.5`
- `max_iterations`: hard safety bound
- `use_siegel`: use Siegel-style swap condition
- `track_transform`: maintain a unimodular transform matrix
