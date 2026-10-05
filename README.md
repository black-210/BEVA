# BEVA

**BEVA** is a C-first numerical library for **lattice analysis, linear algebra, Gram-Schmidt orthogonalization, and LLL lattice reduction**, with a thin Python interface.

The computational core is implemented in C11. Python is only a binding layer.

## Highlights

- C11 core with a small, explicit API.
- Integer lattice basis using signed 64-bit coordinates.
- Gram-Schmidt orthogonalization and GSO coefficients.
- Classical floating-point LLL reduction with exact `int64_t` basis transformations.
- Lovász condition and size reduction.
- Lattice volume / fundamental parallelepiped volume.
- Orthogonality defect.
- Reduced-basis validation.
- Dense matrix multiplication, transpose, vector multiplication.
- Gaussian-elimination linear solve with pivoting.
- Matrix inverse, determinant, numerical rank.
- Python `ctypes` bindings: no third-party runtime dependency.
- CMake build, tests, examples, and documentation.
- Overflow checks for integer LLL row operations.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The shared library is produced as `libbeva.so` on Linux.

## Python

After building:

```bash
export BEVA_LIBRARY="$PWD/build/libbeva.so"
PYTHONPATH="$PWD/python" python3 -c "import beva; print(beva.version())"
```

Example:

```python
from beva import Lattice, lll_reduce, gram_schmidt

basis = [
    [105, 821, 12],
    [318, 211, 17],
    [72,  19,  301],
]

l = Lattice(basis)
print("volume:", l.volume())
print("defect:", l.orthogonality_defect())
print("iterations:", l.lll())
print("reduced:", l.to_list())
print("is reduced:", l.is_reduced())

print(gram_schmidt([[1, 1, 0], [1, 0, 1], [0, 1, 1]]))
```

## LLL notes

BEVA uses the classical LLL framework:

- `1/4 < delta < 1`
- size reduction: `|mu(k,j)| <= 1/2`
- Lovász condition:
  `||b*_k||² >= (delta - mu(k,k-1)²) ||b*_(k-1)||²`

Basis vectors are rows.

The C implementation performs integer row operations in `int64_t` and rejects operations that would overflow. The Gram-Schmidt state is recomputed after reductions/swaps for robustness and clarity.

## API

See `include/beva/*.h` and `docs/API.md`.

## Design

BEVA intentionally keeps the numerical core independent of Python. This makes the same C library usable from:

- C/C++
- Python
- Rust FFI
- command-line tools
- embedded/native applications

## License

MIT. See `LICENSE`.


# Expanded Research Edition

## Everything retained

C11 core; Linear Algebra; Matrix / Vector; Gaussian elimination; Inverse / determinant / rank; QR + least-squares; Cholesky; Power iteration; Gram-Schmidt + GSO; LLL reduction; delta and eta; Siegel option; int64_t exact basis updates; overflow protection; LLL statistics; unimodular transform tracking; lattice volume; orthogonality defect; nearest-plane; Babai; GCD / LCM / extended GCD / modular inverse / primality; Python ctypes API; CMake + CTest; GitHub Actions CI; Examples; API + design + security docs.

## Added without removing the above

- N-body numerical dynamics
- Poisson 2-D iterative PDE solver
- 1-D wave equation stepper
- Lorentz gamma calculation
- Klein-Gordon dispersion
- Schrödinger dispersion
- Complex quantum-state normalization and expectation values
- Dimensional-analysis compatibility checks
- Physics research documentation
- AGPL-3.0 licensing and black-210 copyright/notice files

## Scientific scope

BEVA is a research toolkit, not a claim that an open problem has been solved. It provides numerical infrastructure for investigating hard problems in quantum mechanics, many-body dynamics, PDEs, statistical mechanics and future numerical-relativity work.
