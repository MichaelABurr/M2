# ARing implementation guide

An ARing is the engine-level arithmetic object used by `ConcreteRing` to
implement a Macaulay2 `Ring`.  The ARing owns the representation of one element
type and provides the operations needed to create, convert, compare, print, and
compute with those elements.

Most ARings live in an `aring-*.hpp` and matching `.cpp` file.  Simple rings
inherit from `SimpleARing<ARingType>`, which supplies RAII element wrappers when
the ARing provides `init`, `init_set`, and `clear`.  Rings whose elements need a
ring-specific context for destruction, such as FLINT finite fields, inherit from
`RingInterface` and define their own nested `Element` and `ElementArray`
wrappers.

## Required ring data

Every concrete ARing should provide these basic declarations and ring metadata.

```c++
class ARingExample : public SimpleARing<ARingExample>
{
 public:
  static const RingID ringID = ring_example;

  typedef ExampleElement ElementType;
  typedef ElementType elem;

  size_t characteristic() const;
  unsigned int computeHashValue(const ElementType& a) const;
  void text_out(buffer& o) const;
};
```

`ringID` identifies the ring implementation to dispatch code in
`aring-glue.hpp`; add a new value to `RingID` in `aring.hpp` before using one.
`ElementType` is the raw element representation.  `elem` is the older local
alias used by many ARing files.

Several rings also provide `cardinality()` and
`typedef std::vector<elem> ElementContainerType`.  Add `cardinality()` when the
size of the ring is meaningful and cheap to report.  Add `ElementContainerType`
when matrix or vector code needs a standard container of raw elements.

Constructors should store all data needed to interpret an element: modulus,
precision, FLINT/GMP/Givaro context, original polynomial ring, primitive
element, tower variables, or any lookup tables.  If elements depend on that
context for memory management, define custom `Element` and `ElementArray`
wrappers instead of using `SimpleARing`.

## Element lifetime

For a `SimpleARing`, implement the lifecycle functions for `ElementType`.

```c++
void init(ElementType& result) const;
void init_set(ElementType& result, const ElementType& a) const;
void set(ElementType& result, const ElementType& a) const;
static void clear(ElementType& result);
void set_zero(ElementType& result) const;
```

`init` must create a valid zero-like element.  `init_set` creates a new valid
element initialized from another element.  `set` assigns into an already
initialized element.  `clear` releases resources held by an initialized element.
For POD element types, these functions are often simple assignments and `clear`
can be a no-op.

## Conversion to ring_elem

`ConcreteRing` stores public ring values as `ring_elem`, so an ARing must define
how to pack and unpack its `ElementType`.

```c++
void to_ring_elem(ring_elem& result, const ElementType& a) const;
void from_ring_elem(ElementType& result, const ring_elem& a) const;
const ElementType& from_ring_elem_const(const ring_elem& a) const;
```

`to_ring_elem` creates the engine-facing value.  `from_ring_elem` copies or
reconstructs into an initialized `ElementType`.  `from_ring_elem_const` should
return a read-only view when the stored `ring_elem` representation allows it.
If no stable reference is possible, follow the patterns of existing rings with
similar storage.

## Constructors from common values

The high-level `Ring` interface calls these methods when creating elements from
integers, rationals, variables, and approximate values.

```c++
void set_from_long(ElementType& result, long a) const;
void set_from_mpz(ElementType& result, mpz_srcptr a) const;
bool set_from_mpq(ElementType& result, mpq_srcptr a) const;
void set_var(ElementType& result, int v) const;
```

`set_from_mpq` returns `false` when the rational cannot be represented, for
example when a denominator maps to zero in a finite field.  Approximate and
interval rings should also provide whichever of these conversions are natural:

```c++
bool set_from_double(ElementType& result, double a) const;
bool set_from_BigReal(ElementType& result, gmp_RR a) const;
bool set_from_Interval(ElementType& result, gmp_RRi a) const;
bool set_from_complex_double(ElementType& result, double re, double im) const;
bool set_from_BigComplex(ElementType& result, gmp_CC a) const;
bool set_from_ComplexInterval(ElementType& result, gmp_CCi a) const;
```

`aring-translate.hpp` detects these methods and uses them for promotion and
lifting between ARing implementations.

## Predicates and ordering

Implement the basic predicates and comparison used by `Ring`.

```c++
bool is_unit(const ElementType& a) const;
bool is_zero(const ElementType& a) const;
bool is_equal(const ElementType& a, const ElementType& b) const;
int compare_elems(const ElementType& a, const ElementType& b) const;
```

`compare_elems` should return `-1`, `0`, or `1`.  The ordering only needs to be
consistent for engine use; it does not have to be mathematically canonical
unless callers rely on that for the ring.

## Arithmetic

`ConcreteRing` forwards public arithmetic to methods on the ARing.  The result
argument is already initialized before each call.

```c++
void negate(ElementType& result, const ElementType& a) const;
void add(ElementType& result, const ElementType& a, const ElementType& b) const;
void subtract(ElementType& result, const ElementType& a, const ElementType& b) const;
void mult(ElementType& result, const ElementType& a, const ElementType& b) const;
void divide(ElementType& result, const ElementType& a, const ElementType& b) const;
void invert(ElementType& result, const ElementType& a) const;
void power(ElementType& result, const ElementType& a, int n) const;
void power_mpz(ElementType& result, const ElementType& a, mpz_srcptr n) const;
void subtract_multiple(ElementType& result, const ElementType& a, const ElementType& b) const;
void syzygy(const ElementType& a, const ElementType& b, ElementType& x, ElementType& y) const;
```

`subtract_multiple` updates `result` by subtracting `a * b`.  Division and
inversion should throw the appropriate engine exception, such as
`exc::division_by_zero_error`, when the operation is invalid.

## Output, random elements, and evaluation

These methods connect an ARing to printing, random element generation, and ring
maps.

```c++
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
```

`text_out` prints the ring itself, while `elem_text_out` prints an element.
`eval` is used when applying a `RingMap`; for constant rings this is often just
conversion into the target ring.

## Optional hooks

Some ARings support additional behavior used by finite field, approximate, or
factory code.  Add these only when they are meaningful for the new ring.

```c++
long coerceToLongInteger(const ElementType& a) const;
void getGenerator(ElementType& result) const;
long discreteLog(const ElementType& a) const;
const PolynomialRing& originalRing() const;
void lift_to_original_ring(ring_elem& result, const ElementType& a) const;
M2_arrayint getModPolynomialCoeffs() const;
M2_arrayint getGeneratorCoeffs() const;
M2_arrayint fieldElementToM2Array(ElementType a) const;
unsigned long get_precision() const;
```

Prime finite fields commonly provide integer coercion, generators, and
discrete logarithms.  Extension fields usually keep the original quotient
polynomial ring and provide access to modulus and generator data.  Real and
complex rings provide precision and approximate-value conversions.

## Promotion and lifting

Promotion and lifting happen in two layers.  `ConcreteRing::promote` and
`ConcreteRing::lift` in `aring-glue.hpp` dispatch by `RingID`; the actual
element-level conversions are implemented in `aring-translate.hpp`.

When adding a ring, update `aring-glue.hpp` if the new ring participates in a
new source or target pair.  Update `aring-translate.hpp` with `mypromote` or
`mylift` overloads when conversion cannot be handled by the generic
`set_from_*` helpers.

## Integration checklist

Create the header and implementation files under `basic-rings/`.

Add a `RingID` value in `aring.hpp`.

Add source files and installed headers to `M2/Macaulay2/e/CMakeLists.txt` and
the matching autotools source lists if they are still in use.

Include the new header from `aring-glue.hpp` or `aring-translate.hpp` only when
dispatch or conversion code needs the complete type.

Add a factory entry in `interface/aring.cpp` when Macaulay2 code needs to create
the ring directly.

Add promotion and lifting dispatch in `aring-glue.hpp` and element conversion
helpers in `aring-translate.hpp` when the ring has natural maps to or from
existing ARings.

Check matrix and vector arithmetic support.  Generic mutable matrices work
through `ConcreteRing`, but specialized dense vector arithmetic may require an
entry in `vector-arithmetic.hpp`.

Add focused unit tests under `unit-tests/` and, when the ring is user-visible,
normal Macaulay2 tests for construction, coercion, arithmetic, printing,
promotion, and failure cases.
