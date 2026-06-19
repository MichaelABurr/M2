# ARing quick-start guide

This guide is for contributors adding a new arithmetic ring, or ARing, under
`M2/Macaulay2/e/basic-rings`.  An ARing is the low-level C++ object that knows
how to store and compute with one coefficient type.  `ConcreteRing<ARingType>`
wraps that object and exposes it through the virtual Macaulay2 `Ring`
interface.

The big picture is:

```text
Macaulay2 code
  -> Ring virtual methods
  -> ConcreteRing<ARingType>
  -> ARingType methods
  -> raw ElementType storage, often backed by GMP, FLINT, MPFR, MPFI, or Givaro
```

`ConcreteRing` is the adapter.  It receives a `ring_elem` from the rest of the
engine, asks the ARing to unpack it into the ARing's `ElementType`, calls the
ARing operation, and asks the ARing to pack the result back into a `ring_elem`.
That means a new ARing must define both the mathematical operations and the
conversion boundary between `ring_elem` and its own element representation.

## The shortest path

Start from an existing ARing whose storage model is close to yours.

Use `SimpleARing<YourRing>` when elements can be initialized, copied, and
cleared without carrying extra context.  Examples include many integer,
rational, prime-field, real, and interval rings.

Use `RingInterface` directly when each element needs ring-owned context for
destruction or array management.  Examples include FLINT extension fields and
recursive tower elements; these rings define their own nested `Element` and
`ElementArray` wrappers.

Create `aring-your-ring.hpp` and usually `aring-your-ring.cpp`.  Add a `RingID`
in `aring.hpp`, add the source/header to the build lists, add factory support in
`interface/aring.cpp` if Macaulay2 code should construct the ring directly, and
add tests before opening the PR.

## Minimal class shape

Most simple ARings start with this shape.  Exact signatures may pass small POD
types by value instead of by const reference, but the roles are the same.

```c++
class ARingExample : public SimpleARing<ARingExample>
{
 public:
  static const RingID ringID = ring_example;

  typedef ExampleElement ElementType;
  typedef ElementType elem;

  ARingExample(/* modulus, precision, context, etc. */);

  size_t characteristic() const;
  unsigned int computeHashValue(const ElementType& a) const;
  void text_out(buffer& o) const;

  void init(ElementType& result) const;
  void init_set(ElementType& result, const ElementType& a) const;
  void set(ElementType& result, const ElementType& a) const;
  static void clear(ElementType& result);
  void set_zero(ElementType& result) const;

  void to_ring_elem(ring_elem& result, const ElementType& a) const;
  void from_ring_elem(ElementType& result, const ring_elem& a) const;
  const ElementType& from_ring_elem_const(const ring_elem& a) const;

  void set_from_long(ElementType& result, long a) const;
  void set_from_mpz(ElementType& result, mpz_srcptr a) const;
  bool set_from_mpq(ElementType& result, mpq_srcptr a) const;
  void set_var(ElementType& result, int v) const;

  bool is_unit(const ElementType& a) const;
  bool is_zero(const ElementType& a) const;
  bool is_equal(const ElementType& a, const ElementType& b) const;
  int compare_elems(const ElementType& a, const ElementType& b) const;

  void negate(ElementType& result, const ElementType& a) const;
  void add(ElementType& result, const ElementType& a, const ElementType& b) const;
  void subtract(ElementType& result, const ElementType& a, const ElementType& b) const;
  void mult(ElementType& result, const ElementType& a, const ElementType& b) const;
  void divide(ElementType& result, const ElementType& a, const ElementType& b) const;
  void invert(ElementType& result, const ElementType& a) const;
  void power(ElementType& result, const ElementType& a, int n) const;
  void power_mpz(ElementType& result, const ElementType& a, mpz_srcptr n) const;
  void subtract_multiple(ElementType& result,
                         const ElementType& a,
                         const ElementType& b) const;

  void elem_text_out(buffer& o,
                     const ElementType& a,
                     bool p_one,
                     bool p_plus,
                     bool p_parens) const;
  void random(ElementType& result) const;
  void eval(const RingMap* map,
            const ElementType& f,
            int first_var,
            ring_elem& result) const;
};
```

## Mandatory pieces

These are the pieces a normal `ConcreteRing<YourRing>` needs in order to be a
usable ring.

Ring identity and element storage are mandatory.  Provide a unique `ringID`,
the raw `ElementType`, and any constructor data needed to interpret an element:
modulus, characteristic, precision, coefficient field, FLINT context, original
polynomial ring, primitive element, lookup tables, or tower variables.

Element lifetime is mandatory.  `ConcreteRing` creates temporary ARing elements
before almost every operation.  `init`, `init_set`, `set`, `set_zero`, and
`clear` must leave elements valid and must not leak memory.  If `clear` needs
ring context, do not use `SimpleARing`; define custom wrappers as in the FLINT
field and tower ARings.

`ring_elem` conversion is mandatory.  The rest of the engine does not know your
raw `ElementType`.  Implement `to_ring_elem`, `from_ring_elem`, and
`from_ring_elem_const` carefully, following a nearby ring with similar storage
ownership.

Construction from common values is mandatory for the usual coefficient-ring
paths.  Implement `set_from_long`, `set_from_mpz`, `set_from_mpq`, and
`set_var`.  `set_from_mpq` should return `false` when the rational cannot be
represented, such as a denominator becoming zero in a finite field.

Predicates and comparison are mandatory.  Implement `is_zero`, `is_unit`,
`is_equal`, and `compare_elems`.  The comparison only needs to be consistent for
engine use unless the ring has a mathematically meaningful canonical order.

Core arithmetic is mandatory.  At minimum, implement `negate`, `add`,
`subtract`, `mult`, `divide`, `invert`, `power`, `power_mpz`, and
`subtract_multiple`.  `add` and `subtract` are not optional; many algorithms
assume they are cheap and correct.  Operations that are mathematically invalid
should fail explicitly, usually by throwing an engine exception such as
`exc::division_by_zero_error`.

Printing, random elements, and evaluation are mandatory for a fully integrated
ring.  Implement `text_out` for the ring, `elem_text_out` for elements,
`random`, and `eval`.  For a constant coefficient ring, `eval` is often a
conversion into the target ring.

Hashing is mandatory.  `computeHashValue` is called through `ConcreteRing` and
should be stable for equal elements.

## Important integration work

Build-system updates are important.  Add new `.hpp` and `.cpp` files to
`M2/Macaulay2/e/CMakeLists.txt` and any active source lists used by the build.

Factory support is important when the ring is user-visible.  Add creation code
in `interface/aring.cpp` or the relevant factory path so Macaulay2 can request
the new `ConcreteRing<YourRing>`.

Promotion and lifting are important when there are natural maps to or from
existing rings.  `aring-glue.hpp` dispatches by `RingID`, then
`aring-translate.hpp` performs element-level `mypromote` or `mylift` work.  If
generic `set_from_*` methods already cover your case, keep the wiring small;
otherwise add explicit overloads and tests for success and failure cases.

Matrix and vector support is important for performance-sensitive coefficient
rings.  Generic dense and sparse mutable matrices work through `ConcreteRing`,
but specialized dense arithmetic may require an entry in
`vector-arithmetic.hpp`.

Tests are important, not optional.  Add focused unit tests for the ARing itself
and normal Macaulay2 tests when the ring is user-visible.  Cover construction,
coercion, equality, zero and unit predicates, add/subtract/multiply, division
or division failure, powers, printing, promotion/lifting, and representative
error cases.

## Optional or specialized hooks

Only add these when the ring actually supports the behavior or when a caller
needs the hook.

`cardinality()` is useful for finite rings where the size is meaningful and
cheap to report.

`ElementContainerType` is useful when vector or matrix code needs a standard
container of raw elements.

Approximate and interval setters are useful for real and complex rings:

```c++
bool set_from_double(ElementType& result, double a) const;
bool set_from_BigReal(ElementType& result, gmp_RR a) const;
bool set_from_Interval(ElementType& result, gmp_RRi a) const;
bool set_from_complex_double(ElementType& result, double re, double im) const;
bool set_from_BigComplex(ElementType& result, gmp_CC a) const;
bool set_from_ComplexInterval(ElementType& result, gmp_CCi a) const;
```

Finite-field and extension-field hooks are specialized:

```c++
long coerceToLongInteger(const ElementType& a) const;
void getGenerator(ElementType& result) const;
long discreteLog(const ElementType& a) const;
const PolynomialRing& originalRing() const;
void lift_to_original_ring(ring_elem& result, const ElementType& a) const;
M2_arrayint getModPolynomialCoeffs() const;
M2_arrayint getGeneratorCoeffs() const;
M2_arrayint fieldElementToM2Array(ElementType a) const;
```

`syzygy` is optional algorithmic support, not part of the basic arithmetic
definition in the same way that `add` and `subtract` are.  Implement it when
algorithms targeting the ring need coefficient syzygies.  If the ring cannot
support it, keep the behavior explicit and covered by tests rather than letting
callers silently get nonsense.

Precision hooks such as `get_precision()` are specialized for real, complex,
and interval rings.

## Quick checklist

Choose the closest existing ARing and copy its structure.

Decide between `SimpleARing` and custom `RingInterface` wrappers.

Define `RingID`, `ElementType`, constructor data, and lifetime management.

Implement `ring_elem` conversion, construction from ZZ/QQ, predicates,
comparison, hashing, core arithmetic, printing, random elements, and evaluation.

Wire the ring into build files and factory code.

Add promotion/lift dispatch only for natural maps that should exist.

Add unit tests and user-visible Macaulay2 tests.

Run the ARing tests and at least the relevant engine/Macaulay2 smoke tests
before sending the PR.
