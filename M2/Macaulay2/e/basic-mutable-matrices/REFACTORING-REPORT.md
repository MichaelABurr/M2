# Basic Mutable Matrices Refactoring Report

Date: 2026-09-16

## Executive summary

The dense-matrix implementation and linear-algebra code have been reorganized
so that generic template interfaces remain in generic headers, while code for
concrete coefficient rings is declared and defined in ring-specific files.
Concrete implementations now live primarily in `.cpp` files, and concrete
ARing headers are included by those implementation files rather than by the
generic matrix headers.

The principal results are:

- `mat-linalg.hpp` was reduced from 1,486 lines to 394 lines.
- Twenty-three ring-specific implementation `.cpp` files now contain DMat,
  linear-algebra, or in-place LU implementation code.
- Special linear-algebra operations are expressed as explicit template
  specializations, not unrelated overloads.
- Generic headers use forward declarations instead of including headers for
  their specializations.
- `DMatLUinPlace::mLU` is now owned through `std::unique_ptr<Mat>`.
- `lapack.hpp` no longer includes concrete ARing headers.
- `smat.hpp` no longer includes `rings/ZZp.hpp` as an accidental umbrella.
- CMake and Autotools both build the reorganized source set.
- The complete CMake test run passed 550 of 550 tests, and the Autotools unit
  test run passed 259 of 259 tests.

## Original problem

`mat-linalg.hpp` previously mixed several distinct responsibilities:

1. generic matrix-operation templates;
2. declarations for concrete ring implementations;
3. implementations specialized for concrete DMat types;
4. inclusion of concrete ARing definitions needed only by those
   implementations.

That arrangement made generic clients compile ring-specific implementation
code and made include order significant. It also encouraged base headers to
include specialization headers, reversing the desired dependency direction.

The intended dependency direction is now:

```text
generic declaration header
        ^
        |
specialization declaration header
        ^
        |
specialization implementation.cpp
        |
        +-- concrete ARing and external-library headers
```

The generic header forward-declares specialized types where necessary. A
specialization header includes its generic base, never the other way around.

## Work completed

### 1. Generic linear-algebra interface separated from implementations

`mat-linalg.hpp` now contains the public generic operations, fallback behavior,
and declarations needed for DMat dispatch. Reusable generic DMat algorithms are
in `dmat-linalg.hpp`.

Concrete implementations were separated into corresponding header/source
pairs:

- `dmat-linalg-cc.hpp/.cpp`
- `dmat-linalg-ccc.hpp/.cpp`
- `dmat-linalg-gf-flint.hpp/.cpp`
- `dmat-linalg-gf-flint-big.hpp/.cpp`
- `dmat-linalg-qq-flint.hpp/.cpp`
- `dmat-linalg-qq-gmp.hpp/.cpp`
- `dmat-linalg-rr.hpp/.cpp`
- `dmat-linalg-rrr.hpp/.cpp`
- `dmat-linalg-zz-flint.hpp/.cpp`
- `dmat-linalg-zz-gmp.hpp/.cpp`
- `dmat-linalg-zzp-ffpack.hpp/.cpp`
- `dmat-linalg-zzp-flint.hpp/.cpp`

Operations such as `rank`, `determinant`, `inverse`, `mult`, `nullSpace`,
`solveLinear`, `rankProfile`, and related routines are declared as explicit
template specializations in the ring-specific headers and defined in the
matching `.cpp` files.

The operation names were kept aligned with the public matrix API, including
`addMultipleTo` and `subtractMultipleTo`.

### 2. DMat specializations separated into header/source pairs

Concrete DMat classes that require FLINT details were split so their headers
declare the specialization and its storage/API, while methods requiring the
complete ARing type are defined in `.cpp` files:

- `dmat-gf-flint.hpp/.cpp`
- `dmat-gf-flint-big.hpp/.cpp`
- `dmat-qq-flint.hpp/.cpp`
- `dmat-zz-flint.hpp/.cpp`
- `dmat-zzp-flint.hpp/.cpp`

`dmat.hpp` no longer includes all of those child specialization headers. It
forward-declares the ARing classes and DMat specializations instead. This
restores the correct base-to-specialization dependency direction.

The former central `dmat.cpp` source was removed from the build because its
implementation responsibilities have been distributed to the corresponding
ring-specific sources.

### 3. In-place LU specializations separated

`dmat-lu-inplace.hpp` retains the generic template implementation. Its size was
reduced from 427 lines to 238 lines, and ring-specific `computeLU`
specializations were moved to:

- `dmat-lu-inplace-cc.hpp/.cpp`
- `dmat-lu-inplace-ccc.hpp/.cpp`
- `dmat-lu-inplace-gf-flint.hpp/.cpp`
- `dmat-lu-inplace-gf-flint-big.hpp/.cpp`
- `dmat-lu-inplace-rr.hpp/.cpp`
- `dmat-lu-inplace-rrr.hpp/.cpp`

Additional LU implementations were separated into their own implementation
files for QQ and the ZZp backends.

`DMatLUinPlace::mLU` was changed from an embedded `Mat` to:

```cpp
std::unique_ptr<Mat> mLU;
```

This permits the class declaration to avoid requiring the complete concrete
matrix/ring representation merely to determine its own object layout.

### 4. Type aliases and declarations localized

The shared `dmat-types.hpp` approach is no longer present, and there are no
remaining inclusions of `dmat-types.hpp` under `Macaulay2/e`. Type aliases and
forward declarations are located with the concrete interfaces that use them.

This reduces the number of unrelated concrete types introduced by including a
single matrix header.

### 5. LAPACK header dependency reduced

Concrete ARing includes for RR, CC, RRR, and CCC were removed from
`lapack.hpp`. The header now uses ARing forward declarations and DMat aliases.
The complete ARing headers are included in `lapack.cpp`, where their definitions
are actually required.

Necessary LAPACK conversion helper declarations were added to the header so
ring-specific implementation files can call them without importing concrete
ARing definitions through `lapack.hpp`.

### 6. Sparse-matrix umbrella dependency removed

`smat.hpp` formerly included `rings/ZZp.hpp` despite containing no ZZp-specific
operation. That include was replaced by the direct dependencies actually used:

- allocation support;
- M2 interface types;
- error reporting;
- `cassert` and `cstddef`;
- forward declarations for `Ring` and `ring_elem`.

This removes the accidental path:

```text
smat.hpp -> rings/ZZp.hpp -> rings/ring.hpp -> basic-rings/aring.hpp
```

Removing that umbrella exposed two pre-existing hidden dependencies. Within
the allowed directory, they were handled by:

- adding `exceptions.hpp` directly to `mat-arith.hpp`, which throws
  `exc::engine_error`; and
- placing the legacy complete-type compatibility includes for `Z_mod` and
  `M2::ARingZZGMP` at the concrete aggregation boundary in
  `dmat-linalg-all.hpp`, rather than restoring them to generic `smat.hpp` or
  changing SLP.

### 7. Build-system integration

Both build descriptions were updated for the new source files:

- `Macaulay2/e/CMakeLists.txt`
- `Macaulay2/e/Makefile.files.in`

The old `basic-mutable-matrices/dmat.cpp` entry was removed, and the new DMat,
linear-algebra, LU, and in-place LU implementation sources were registered.

Call sites and integration headers were adjusted to include the appropriate
new declaration headers. SLP retains its original `DMatLinAlg<RT>` call; its
required specialization declarations and complete `ARingZZGMP` type are
provided by the matrix aggregation header.

## Current header-to-ARing dependency audit

All 39 `.hpp` files in `basic-mutable-matrices` were individually passed
through the compiler dependency generator. Thirty-six have no direct or
indirect dependency on `basic-rings`.

Three headers still reach ARing infrastructure:

| Header | Kind | Dependency path |
| --- | --- | --- |
| `dmat-qq-interface-flint.hpp` | Direct | `basic-rings/aring-QQ.hpp`, `basic-rings/aring-ZZ-gmp.hpp` |
| `dmat-linalg-all.hpp` | Direct and indirect compatibility dependency | `basic-rings/aring-ZZ-gmp.hpp`; `rings/ZZp.hpp -> rings/ring.hpp -> basic-rings/aring.hpp` |
| `mat.hpp` | Indirect | `ring-elements/ring-element.hpp -> rings/ring.hpp -> basic-rings/aring.hpp` |

No other header in this directory currently reaches `basic-rings` when
preprocessed on its own.

## Verification performed

### CMake build (`build-gpt`)

- `M2-engine`: built and linked successfully.
- `M2-unit-tests`: built and linked successfully.
- Auxiliary `memtailor`, `mathic`, and `mathicgb` unit-test executables: built
  successfully.
- CTest result: **550/550 tests passed**.

### Autotools build (`build-gpt-autotool`)

- Engine archive: built successfully.
- Unit-test executable: built and linked successfully.
- Test result: **259/259 tests passed**.

### Source checks

- `git diff --check`: passed.
- All 39 `basic-mutable-matrices` headers completed dependency preprocessing.

Some pre-existing/non-fatal warnings appeared during verification, including
unused variables in `lapack.cpp`, stale coverage checksum messages in the CMake
build directory, a compiler bounds warning in `NewF4Test`, and GC finalization
warnings during a few tests. None caused a build or test failure.

## Benefits

### Reduced compile-time coupling

A translation unit including a generic matrix header no longer automatically
parses most concrete ARing and backend implementation code. This should reduce
front-end work, memory consumption, and unnecessary recompilation after a
ring-specific change.

This benefit is expected from the include graph but has not yet been measured
with a clean-build timing benchmark.

### Better incremental builds

Changing a ZZp/FLINT implementation now primarily invalidates its associated
`.cpp` file and direct dependents instead of every translation unit that
includes the former monolithic header. The addition of 23 focused
implementation sources gives the build system finer-grained recompilation
boundaries.

### Clear ownership of dependencies

Concrete ARing and external-library dependencies are generally located in the
implementation files that use them. Missing includes are detected closer to
their owner instead of being silently supplied by an unrelated umbrella
header.

### More stable public headers

The generic API is smaller and less sensitive to changes in FLINT, FFPACK,
Givaro, GMP, or LAPACK implementation details. This improves maintainability
and makes it easier to reason about which declarations are available from each
header.

### Correct specialization semantics

Representing concrete algorithms as explicit template specializations makes
their relationship to the primary `MatrixOps` templates unambiguous. It also
reduces the risk that overload resolution selects a generic fallback because a
ring-specific overload was not visible at a call site.

### Improved extensibility

Adding another coefficient-ring backend now has a repeatable structure:

1. forward-declare its ARing type;
2. declare its DMat or algorithm specializations in a focused header;
3. define them in a focused `.cpp` containing the concrete dependencies;
4. register the source with CMake and Autotools;
5. add it to the concrete aggregation boundary only when interface templates
   require visibility of its declarations.

### Lower regression risk

The refactor was validated through both supported build systems and their
engine unit tests. That is particularly important because moving template
specializations can otherwise produce failures only at instantiation or link
time.

## Remaining opportunities

### 1. Split `dmat-qq-interface-flint.hpp`

This is the only header with direct `basic-rings` includes. Its methods are
currently inline. Moving the method definitions to a `.cpp` would allow the
header to forward-declare `M2::ARingQQ` and `M2::ARingZZGMP`.

FLINT headers would remain necessary while `fmpq_mat_t` and `fmpz_mat_t` are
stored by value. A PImpl or pointer-owned representation would be needed to
hide those types as well.

### 2. Reduce `mat.hpp`'s indirect dependency

`mat.hpp` includes `ring-elements/ring-element.hpp` but uses `RingElement` only
through pointers. It can likely replace that include with the needed interface
type header, a `ring_elem` forward declaration, and a direct `<vector>` include.
That would remove its path to `basic-rings/aring.hpp`.

### 3. Remove the `dmat-linalg-all.hpp` compatibility includes

The `rings/ZZp.hpp` and `basic-rings/aring-ZZ-gmp.hpp` includes in
`dmat-linalg-all.hpp` are compatibility workarounds for clients outside this
directory that instantiate matrix code requiring the complete `Z_mod` or
`M2::ARingZZGMP` type without including its defining header. In particular,
the latter keeps the existing SLP `DMatLinAlg<RT>` call compiling without an
SLP-layer include or source change.

The architecturally clean fix is for those clients—notably
`basic-rings/vector-arithmetic.hpp`, tests that dynamically cast to `Z_mod`,
and SLP code instantiating matrices over `ARingZZGMP`—to include the concrete
definitions they require. Once those owners are corrected, the compatibility
includes can be removed. This requires changing files outside
`basic-mutable-matrices`.

### 4. Reconsider the aggregation header

`dmat-linalg-all.hpp` is still required because `mutablemat` interface template
code needs all explicit specialization declarations to be visible. A larger
follow-up could explicitly instantiate the relevant `MutableMat` operations in
implementation files, allowing interface headers to include narrower
declaration sets and potentially eliminating the aggregation header.

### 5. Clean working-tree artifacts before committing

The untracked editor backup file
`basic-mutable-matrices/#dmat-lu-gf-flint.hpp#` should not be committed. The
many new `.hpp/.cpp` pairs are also currently untracked and must be deliberately
added to version control. No commit was created during this work.

### 6. Optional warning cleanup and performance measurement

After the structural changes settle, a separate cleanup can address the
non-fatal warnings observed during the builds. Clean-build and incremental-build
timings before and after the refactor would quantify the expected compilation
benefit.

## Recommended next steps

1. Review and add the new specialization files, excluding the editor backup.
2. Split `dmat-qq-interface-flint.hpp` into declarations and implementation.
3. Decouple `mat.hpp` from `ring-elements/ring-element.hpp`.
4. If changes outside this directory are permitted, fix the direct `Z_mod`
   consumers and remove the compatibility include from `dmat-linalg-all.hpp`.
5. Run one clean CMake build and one clean Autotools build before committing,
   followed by the same unit-test suites recorded above.

## Conclusion

The refactor has substantially reduced the monolithic generic headers and
established consistent ownership boundaries for DMat classes, linear-algebra
specializations, and LU implementations. The current tree builds and passes
the available engine tests under both CMake and Autotools. Most remaining
dependency work is confined to three headers and can be handled incrementally
without reversing the new specialization structure.
