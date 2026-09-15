# Mutable matrices in the M2 engine — overview + refactoring direction

## Context
Mike is presenting the mutable matrix code to the Cornell C++ refactoring group, and wants a
plan to make it simpler, at least as fast, and decoupled: today one fat `MutableMatrix`
virtual interface (~60 methods) forwards every operation through `MutableMat<T>` to
DMat/SMat free functions. Part A is the talk overview; Part B is the proposed direction.
(Engine paths are relative to `M2/Macaulay2/e/`.)

---

## Part A — Overview (for the presentation)

### The types involved, from the user down
- **M2 front end:** `Matrix` and `MutableMatrix` are M2 types (hash tables). Each one essentially
  holds a pointer to a raw object: `RawMatrix` or `RawMutableMatrix` (`m2/matrix.m2`, `m2/mutablemat.m2`).
- **Interpreter:** `d/engine.dd` declares `RawMatrix` and `RawMutableMatrix` as pointers to the
  engine classes `const Matrix*` and `MutableMatrix*`.
- **Engine, immutable side:** `Matrix` (`matrices/matrix.hpp`). Not part of this refactor, but
  conversions go both ways.
- **Engine, mutable side:** `MutableMatrix` (`basic-mutable-matrices/mat.hpp`) is an abstract base
  class with about 60 virtual methods: entries, row/column operations, arithmetic, linear algebra, numerics.
- **`MutableMat<Mat>`** (`mutable-matrices/mutablemat-defs.hpp`, `-imp.hpp`) is the only implementation
  of `MutableMatrix`.
  - It is a template over the concrete matrix type: `DMat<RT>` for dense or `SMat<RT>` for sparse,
    where `RT` is an ARing such as `ARingZZpFlint`. So there is one instance per ring type and density.
  - It holds the concrete matrix, checks bounds, converts `ring_elem` to the ARing's element type,
    and forwards each virtual method to a free function.
- **`DMat<RT>` (dense)** comes in these kinds:
  - **Default template** (`dmat.hpp`): a row-major array of ring elements. Used by most rings:
    ZZ/p (M2's own and FFPACK), GF (M2's own), ZZ and QQ via GMP, RR, CC, RRR, CCC, RRi, CCi,
    tower rings, and the generic `CoefficientRingR` fallback.
  - **Flint specializations** (`dmat-*-flint.hpp`): wrap flint's own matrix types for ZZ (`fmpz_mat`),
    QQ (`fmpq_mat`), ZZ/p (`nmod_mat`), GF (`fq_zech_mat`) and big GF (`fq_nmod_mat`).
  - FFPACK, LAPACK and Eigen have no matrix class here. They use the default `DMat` and either use
    its array directly (FFPACK) or copy it in and out (LAPACK, Eigen).
- **`SMat<RT>` (sparse)** (`smat.hpp`): one template for all rings. Each column is a linked list
  of (row, coefficient) nodes.
- **The algorithms** are free functions on `DMat`/`SMat`, not methods:
  - `MatElementaryOps<Mat>` (`mat-elem-ops.hpp`): row and column operations.
  - `MatrixOps::` in `mat-arith.hpp` (add, subtract, transpose, ...) and `mat-linalg.hpp`
    (rank, det, inverse, solve, LU, eigen, SVD, ...).
  - `DMatLinAlg<RT>` (`dmat-lu*.hpp`): M2's own LU-based algorithms for fields, with specializations
    for FFPACK, flint ZZ/p and GMP QQ.

### Layers, top to bottom
Paths starting with `m2/` or `d/` are under `M2/Macaulay2/`; the rest are under `M2/Macaulay2/e/`.

| Layer | Where | What it does |
|---|---|---|
| M2 top level | `m2/mutablemat.m2`, packages (EngineTests, LLLBases, NumericalAlgebraicGeometry, SLPexpressions, …) | User-facing functions |
| Interpreter | `d/interface.dd`, `d/engine.dd` | Exposes the `raw*` functions; `+`, `-`, `*` are inline C code |
| C API | `interface/mutable-matrix.cpp` | About 72 `raw*` functions, all calling virtual methods |
| Abstract class | `basic-mutable-matrices/mat.hpp` | `MutableMatrix`: about 60 pure virtual methods |
| Wrapper | `mutable-matrices/mutablemat-{defs,imp}.hpp` | `MutableMat<Mat>`: bounds checks, element conversion, `dynamic_cast` of matrix arguments |
| Concrete matrices | `dmat.hpp`, `dmat-*-flint.hpp`, `smat.hpp` | `DMat<RT>` (default + 5 flint), `SMat<RT>` |
| Algorithms | `mat-elem-ops.hpp`, `mat-arith.hpp`, `mat-linalg.hpp`, `dmat-lu*.hpp` | Free functions on concrete matrices |
| Backends | `dmat.cpp`, `lapack.cpp`, `eigen.cpp`, flint calls in `mat-linalg.hpp` | FFPACK, LAPACK, Eigen with mpreal, flint |

### Key mechanics
- **Creation:** `MutableMatrix::zero_matrix(R, r, c, dense)` → virtual `Ring::makeMutableMatrix`;
  `ConcreteRing<RT>` (basic-rings/aring-glue.hpp:43) makes `MutableMat<DMat<RT>>` or `<SMat<RT>>`.
  `RingZZ`/`Z_mod` override in mat.cpp; fallback `CoefficientRingR` (generic ring_elem arithmetic).
- **Dispatch of a linear-algebra call** e.g. `rawLinAlgRank(M)` → `M->rank()` (virtual) →
  `MutableMat<T>::rank` → `MatrixOps::rank(mat)` → overload resolution picks, in order:
  concrete inline overload (e.g. `DMatZZpFlint` → flint) > `template<RT> rank(DMat<RT>)` →
  `DMatLinAlg<RT>` (possibly specialized: FFPACK, ZZpFlint, QQ) > `template<Mat>` stub that throws
  "not implemented".
- **Binary ops** (`add`, `mult`, `solve`, `LU(L,U)`): the argument is `MutableMatrix*`;
  `coerce<T>()` = `dynamic_cast<MutableMat<T>*>` — mismatch → error.
- **DMat generic:** row-major `ElementType*` array + row-pointer vector.
  **Flint specializations** (ZZ, QQFlint, ZZpFlint, GFFlint, GFFlintBig) wrap `fmpz_mat_t`, etc.
  — five ~90-line near-copies.
- **SMat:** array of per-column singly-linked lists sorted by decreasing row. Supports
  elementary ops + add/sub/negate only; *every* mat-linalg op throws for SMat; transpose throws.
- **Instantiation:** 18 ring types × {DMat, SMat} = **36 MutableMat classes**, each with a full
  ~60-entry vtable, instantiated implicitly via aring-glue.hpp (which includes all of
  mutablemat.hpp → smat, mat-linalg, flint headers, lapack, eigen decls) in ~19 TUs.

### Backend coverage (who actually does the work)
| op | ZZ | QQ | ZZ/p | GF | RR/CC | RRR/CCC | SMat |
|---|---|---|---|---|---|---|---|
| rank/det/nullspace/solve/inverse | flint (gmp ZZ: copy→flint) | flint | FFPACK / flint / M2 LU | flint LU + M2 | M2 LU w/ LAPACK getrf | M2 LU, abs pivot | ✗ |
| rref | flint(ZZ ✗) | flint | flint only | flint only | ✗ | ✗ | ✗ |
| mult / addMultipleTo | flint | flint | FFPACK/flint/M2 loop | flint/M2 | M2 loop (no BLAS!) | M2 loop | ✗ |
| eigen/SVD/leastSquares | ✗ | ✗ | ✗ | ✗ | LAPACK (Eigen if NO_LAPACK) | Eigen+mpreal | ✗ |
| QR | ✗ | ✗ | ✗ | ✗ | LAPACK | ✗ | ✗ |

### Consumers and coupling
- **Loose (virtual interface only):** LLL, fractionfreeLU, mutablecomplex (pruning), cone.cpp
  (lattice points), f4/res-f4 output, fplll/NTL interfaces, finalize/debug.
- **Tight (concrete types):** SLP/SLP-imp.hpp (`dynamic_cast<MutableMat<DMat<RT>>>`, uses
  `MatElementaryOps`, `SubMatrix` windows), interface/matrix.cpp point arrays (`DMat<ARingCC>`),
  schreyer-resolutions/res-f4-m2-interface.cpp (builds DMat + DMatLinAlg directly), eigen.cpp.
- C++ tests: only `unit-tests/DMatZZpTest.cpp` (2 trivial tests) + LatticePoints; real coverage
  is M2-level `packages/EngineTests`.

### Talking points: pain
1. Adding one operation touches mat.hpp, mutablemat-defs/imp, mat-linalg stub + overloads, raw API, d/.
2. Many combinations of operation × ring × dense/sparse are not implemented, and you only find out
   at runtime. The compiler happily builds all of them, because a catch-all template in
   `mat-linalg.hpp` fills every gap with a function that throws "not implemented". Examples:
   - every linear algebra operation on a sparse matrix, even `rank` or `mult` (only elementary
     row/column operations and add/subtract work for `SMat`; even `transpose` throws);
   - `rowReducedEchelonForm` for RR, CC, M2's own ZZ/p and GF (only flint rings and QQ have it);
   - `inverse`, `nullSpace`, `solve` for ZZ via GMP;
   - `QR` for RRR/CCC.
   So the `MutableMatrix` interface promises about 60 operations for every matrix, but many
   matrix types can only do some of them, and nothing records which.
3. Algorithm duplication: 3× setUpperLower, 5× flint DMat wrappers, 4× clean/norm, 2× rref-solve, 14× Eigen precision save/restore, 4× LAPACK RHS padding.
4. Dead code: dmat-ffpack.cpp (all `#if 0`, still compiled), undefined `Lapack` RRR/CCC decls, rawCharPoly/MinPoly/Smith/Hermite/sort_columns stubs, `promote/lift/eval` virtuals not in base.
5. Header weight: every ring TU parses all of linear algebra.
6. Field algorithms (division-based LU) silently run over ZZ, CoefficientRingR, RRi/CCi.
7. Fragile contracts: `mult` accumulates into its result (correct only because callers pass zero);
   `DMat` has shallow implicit `operator=`; `rowMajorArray()` returns `T*&`.

### Bugs found (verified in source = ✔)
- ✔ `solveLinear(DMatZZ)` passes `(B, A)` to `fmpz_mat_solve` and never sizes X (mat-linalg.hpp:792); `nullSpace(DMatZZ)` never sizes result (:786).
- ✔ `DMatLinAlg::inverse` ignores `solve()` result → always "succeeds" on singular matrices (dmat-lu.hpp:540).
- ✔ `SMat::vec_negate` derefs null on empty column, skips last node (smat.hpp:459) → sparse subtract wrong/crash.
- ✔ `SMat::vec_equals` ignores row indices and length mismatch (smat.hpp:296).
- ✔ FFPACK nullspace swaps `mArray` but not row pointers, `delete[]`s GC memory (dmat.cpp:177).
- `solveLinear(DMatQQFlint)` stub returns false; QQ solve indexes empty profile; DMat copy ctor double-init leaks mpz/mpfr; Lapack::QR with rows<cols; eigenvalues(RR) leaks workspace; Eigen SVD thin vs LAPACK full shapes differ.

---

## Part B — Refactoring direction (recommended, staged; each stage shippable)

**Stage 0 – Safety net & baseline.** C++ gtests per (op × ring × dense/sparse) against a
generic reference; timing benchmarks (rank/det/inverse/solve at 500–2000 over ZZ/p (flint, FFPACK),
QQ, ZZ, RR, RRR; sparse elementary ops; LLL; SLP evaluate). Record numbers before any change.

**Stage 1 – Fix bugs + delete dead code** (list above). No structural change.

**Stage 2 – Shrink the virtual interface.** Keep in `MutableMatrix` only the genuinely
type-erased core: shape, ring, dense?, get/set entry, row/col elementary ops, submatrix,
copy, to_matrix. Move everything else (arith, linear algebra, numerics, SLP factory) to free
functions in a `mutable_linalg` layer that dispatch *once*:
- Closed set of concrete types → `using AnyMat = std::variant<DMat<ARingZZ>*, …, SMat<…>*>`;
  `MutableMatrix` exposes one virtual `AnyMat variant()`.
  (`std::variant` and `std::visit` are C++17, which the build already uses: `CMakeLists.txt:23`.)
- **Matrix types vs. algorithms are separate layers.** The class definitions, including the flint
  specializations `template <> class DMat<M2::ARingZZ>` etc., stay in headers
  (`matrices-mutable/dmat.hpp`, `dmat-flint.hpp`, `smat.hpp`, `mutablemat.hpp`, `matrix-variant.hpp`).
  They must: `MutableMat<Mat>` holds `Mat` by value, dispatch.cpp constructs results, and an explicit
  specialization must be declared before any implicit instantiation. `matrix-variant.hpp` needs only
  forward declarations. Flint C headers in `dmat-flint.hpp` are cheap and already pulled in by the
  ARings; the heavy header is FFPACK, included by `aring-ZZp-ffpack.hpp` itself (ARing-refactor item).
- **Layout: implementations grouped by backend/ring, dispatch grouped by operation.**
  ```
  e/linalg/
    dispatch.hpp/.cpp      std::visit entry points (rank, det, inverse, ...); includes backend *declarations*
    not-implemented.hpp    one catch-all template per operation
    backends/
      flint-zz, flint-qq, flint-zzp, flint-gf, ffpack-zzp,     (.hpp = declarations, .cpp = definitions;
      lapack-rr-cc, eigen-rrr-ccc                               only the .cpp includes flint/ffpack/lapack/eigen)
      generic-field.hpp/.cpp   template <RT> DMat<RT> algorithms, explicitly instantiated only for fields
      sparse.hpp               SMat<RT>
  ```
  "What does flint ZZ/p support?" is one header. "Who implements rank?" is a grep of `backends/`, or a
  unit test that tries every op × type and prints the support table (also catches mistyped overloads).
  Adding an op: declare it in supporting backends, one fallback, one entry point in dispatch.cpp.
- Sketch of the overload set that `dispatch.cpp` sees (declarations come from the backend headers):
  ```cpp
  // rank.hpp
  size_t rank(const DMat<M2::ARingZZpFlint>& A);        // flint
  template <typename RT> size_t rank(const DMat<RT>& A); // M2's generic field LU
  [[noreturn]] inline size_t rankNotImplemented() { throw exc::engine_error("rank: not implemented"); }
  template <typename M> size_t rank(const M&) { return rankNotImplemented(); } // e.g. SMat

  // caller
  size_t rank(const MutableMatrix& M)
  {
    return std::visit([](auto* A) { return rank(*A); }, M.variant());
  }
  ```
  No change to the base class or `MutableMat` is needed to add an operation.
- Binary ops visit two variants; mismatched types are rejected in one place.
- Speed: still a single dispatch per call, with fully typed inner loops, so no change.
- Option to discuss with the group: skip `std::variant` and keep one virtual per operation family
  (e.g. `ExactLinAlg`, `NumericalLinAlg`) implemented only by the matrix types that support it.
  Simpler to read, but still more code to touch when adding an op.

**Stage 3 – Consolidate backends.** One traits-based `FlintDMat<Traits>` replacing 5 specializations;
one LU utility; RAII workspaces for LAPACK; one precision guard for Eigen; route RR/CC `mult` to BLAS
(speed win). Gate division-based algorithms on a "field" ring trait (ZZ → flint only).

**Stage 4 – Compile-time decoupling.** Remove `mutablemat.hpp` from aring-glue.hpp; define
`makeMutableMatrix` in one .cpp with explicit instantiations (`extern template`), so only that TU
and the linalg TUs parse flint/lapack/eigen.

**Stage 5 – SMat decision.** Either replace linked lists with per-column sorted vectors
(better locality, supports transpose/mult) or restrict SMat to its actual users (LLL-type
elementary ops). Decide after Stage 0 benchmarks and a grep of package usage of `Dense => false`.

SLP (tight consumer) moves to use concrete `DMat<RT>` from the variant; d/ and m2/ raw API unchanged
through Stage 4.

## Critical files
mat.hpp, mat.cpp, mutablemat-defs.hpp, mutablemat-imp.hpp, mat-linalg.hpp, mat-elem-ops.hpp,
mat-arith.hpp, dmat*.hpp, smat.hpp, dmat-lu*.hpp, dmat.cpp, lapack.cpp, ../eigen.cpp,
basic-rings/aring-glue.hpp, interface/mutable-matrix.cpp, SLP/SLP-imp.hpp.

## Verification
- `packages/EngineTests` (MutableMatrix.Test.*, LinearAlgebra.Test.RR-CC) via in-tree M2: `check "EngineTests"`.
- e/unit-tests gtest target (new Stage 0 tests).
- Benchmark script run before/after each stage; no regression > noise.
- Also: `check` for LLLBases, NumericalAlgebraicGeometry, SLPexpressions, Complexes, LocalRings, Triangulations.

---

## Appendix — dispatch with `std::variant` (C++17), worked example

### 1. One header listing the concrete types (forward declarations only)
```cpp
// mutable-matrices/matrix-variant.hpp
#include <variant>
namespace M2 { class ARingZZ; class ARingQQFlint; class ARingZZpFlint; class ARingZZpFFPACK;
               class ARingZZp; class ARingGFM2; class ARingRR; class ARingCC; /* ... */ }
template <typename RT> class DMat;
template <typename RT> class SMat;

using MatrixVariant = std::variant<
    DMat<M2::ARingZZ>*,       SMat<M2::ARingZZ>*,
    DMat<M2::ARingQQFlint>*,  SMat<M2::ARingQQFlint>*,
    DMat<M2::ARingZZpFlint>*, SMat<M2::ARingZZpFlint>*,
    DMat<M2::ARingRR>*,       SMat<M2::ARingRR>*,
    /* ... all 36 ... */ >;
using ConstMatrixVariant = std::variant<const DMat<M2::ARingZZ>*, /* ... same, const ... */ >;
```

### 2. Base class gains two virtuals (and loses ~40)
```cpp
class MutableMatrix : public MutableEngineObject {
 public:
  virtual MatrixVariant      variant() = 0;
  virtual ConstMatrixVariant variant() const = 0;
  // core only: n_rows, n_cols, get_ring, get/set_entry, row/col ops, copy, to_matrix ...
};

template <typename Mat>
class MutableMat : public MutableMatrix {
  const Ring* mRing;
  Mat mat;
 public:
  MutableMat(const Ring* R, Mat&& m) : mRing(R), mat(std::move(m)) {}
  MatrixVariant      variant() override       { return &mat; }
  ConstMatrixVariant variant() const override { return &mat; }
};
```
Forgetting to add a new ring type to the variant makes `return &mat;` fail to compile.

### 3. The overload set for `rank`, as seen in dispatch.cpp
In the real layout the specific overloads are *declared* in `backends/*.hpp` (defined in their .cpp),
the generic template comes from `generic-field.hpp`, and the catch-all from `not-implemented.hpp`.
Shown inline here for readability:
```cpp
namespace linalg {
inline size_t rank(const DMat<M2::ARingZZpFlint>& A) { return nmod_mat_rank(A.nmod_mat()); }
inline size_t rank(const DMat<M2::ARingZZ>& A)       { return fmpz_mat_rank(A.fmpz_mat()); }
size_t        rank(const DMat<M2::ARingZZpFFPACK>& A);          // FFPACK, in a .cpp

template <typename RT>
size_t rank(const DMat<RT>& A) { return DMatLinAlg<RT>(A).rank(); }  // generic dense field

template <typename Mat>
size_t rank(const Mat&) { throw exc::engine_error("rank: not implemented for this ring/matrix type"); }

inline size_t rank(const MutableMatrix& M)
{
  return std::visit([](auto* A) -> size_t { return rank(*A); }, M.variant());
}
}
```
Overload resolution per alternative: exact non-template > `DMat<RT>` template > catch-all.

### 4. Other operation shapes
Returning a new matrix:
```cpp
inline MutableMatrix* inverse(const MutableMatrix& M)
{
  return std::visit([&](auto* A) -> MutableMatrix* {
      using Mat = std::remove_const_t<std::remove_pointer_t<decltype(A)>>;
      Mat X(A->ring(), A->numRows(), A->numColumns());
      if (not inverse(*A, X)) return nullptr;
      return new MutableMat<Mat>(M.get_ring(), std::move(X));
    }, M.variant());
}
```
Two matrix arguments: visit one, `get_if` the other (avoid 36×36 instantiations):
```cpp
inline MutableMatrix* mult(const MutableMatrix& A, const MutableMatrix& B)
{
  return std::visit([&](auto* a) -> MutableMatrix* {
      auto bv = B.variant();
      auto* b = std::get_if<decltype(a)>(&bv);
      if (b == nullptr) throw exc::engine_error("mult: matrices must have the same ring and density");
      using Mat = std::remove_const_t<std::remove_pointer_t<decltype(a)>>;
      Mat C(a->ring(), a->numRows(), (*b)->numColumns());
      mult(*a, **b, C);
      return new MutableMat<Mat>(A.get_ring(), std::move(C));
    }, A.variant());
}
```
Result type differs (eigenvalues RR → CC): keep `EigenTypes<RT>`, and build a
`DMat<typename EigenTypes<RT>::EigenvalueType>` inside the lambda.

### 5. What the class definitions look like, and where ARing headers are needed
- `matrices-mutable/mutable-matrix.hpp`: `MutableMatrix` with ~25 core virtuals (shape, ring,
  `variant()`, `copy()`, `makeZero(r,c)`, `to_matrix()`, get/set entry, elementary row/col ops,
  submatrix). No ARing, DMat or flint includes; only `matrix-variant.hpp` (forward declarations).
- `matrices-mutable/mutablemat.hpp`: `template <typename Mat> class MutableMat` with member
  *declarations*; `mutablemat-imp.hpp` has the bodies and is included only by the instantiation files.
- Complete ARing types are needed only in:
  1. `matrices-mutable/instances-{flint,ffpack,numeric,other}.cpp`: `template class MutableMat<DMat<RT>>;`
     and `<SMat<RT>>` for that group's rings, plus the definitions of `ConcreteRing<RT>::makeMutableMatrix`
     (moved out of aring-glue.hpp, which then includes nothing from the matrix code).
  2. each `linalg/backends/*.cpp`, for its own rings;
  3. SLP / point-array code, for the rings they use.
- `dispatch.cpp` needs **no** ARing headers:
  - calls resolve against declarations (incomplete types are fine);
  - generic field algorithms are exposed as plain per-ring overloads declared in `generic-field.hpp`
    and defined in `generic-field.cpp` via the template (explicit list; ZZ/CoefficientRingR fall to
    "not implemented" instead of a link error);
  - result matrices are created with `M.makeZero(r,c)` and filled through the same-type `get_if`.
- Rejected alternative: `MutableMatrix` as a non-virtual class holding a variant *by value*. That forces
  every ARing/flint/FFPACK header into every file that uses `MutableMatrix`.

### Caveats
- Compile cost: each `visit` instantiates 36 lambdas, but only in that operation's .cpp
  (today: 36 full vtables in every TU including aring-glue.hpp).
- Runtime: index + jump, comparable to a virtual call; inner loops unchanged.
- All overloads must be declared before the visiting lambda: dispatch.cpp includes every backend header first.
- Generic `DMat<RT>` templates: explicitly instantiate in generic-field.cpp for exactly the field rings,
  so field algorithms can no longer be applied silently to ZZ or CoefficientRingR.
- The catch-all hides a mistyped overload (e.g. missing `const`): tests must hit each backend.
- Variant alternatives must be unique; `ARingZZp` from `Z_mod` and `ConcreteRing` is one type anyway.
