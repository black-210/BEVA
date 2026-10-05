# BEVA design

BEVA separates exact lattice storage from floating-point numerical analysis.

`int64_t` is used for basis and transformation matrices. GSO uses `double`.
Every integer row update is checked with `__int128` before being committed.

The architecture is intentionally modular:

- status/error layer
- vector/matrix layer
- decomposition layer
- lattice layer
- optional Python ABI layer

The LLL engine is deterministic for a fixed compiler/platform floating-point
environment and has a hard iteration limit to prevent pathological inputs from
running forever.
