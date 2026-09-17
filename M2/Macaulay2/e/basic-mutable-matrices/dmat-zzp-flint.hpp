// Copyright 2013  Michael E. Stillman

#ifndef M2_BASIC_MUTMAT_DMAT_ZZP_FLINT_HPP_
#define M2_BASIC_MUTMAT_DMAT_ZZP_FLINT_HPP_

#include "basic-mutable-matrices/dmat.hpp"

// The following needs to be included before any flint files are included.
#include <M2/gc-include.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <flint/nmod_mat.h>  // for nmod_mat_t, nmod_mat_clear, nmod_mat_init, nmod_m...
#pragma GCC diagnostic pop

//////////////////////////////////////////////////////////////
// Flint: use nmod_mat for implementation of dense matrices //
//////////////////////////////////////////////////////////////

template <>
class DMat<M2::ARingZZpFlint>
{
 public:
  typedef M2::ARingZZpFlint ACoeffRing;
  typedef ACoeffRing CoeffRing;
  typedef mp_limb_t ElementType;
  typedef ElementType elem;

  DMat();
  DMat(const ACoeffRing& R, size_t nrows, size_t ncols);
  DMat(const DMat<ACoeffRing>& M);
  ~DMat();

  // swap the actual matrices of 'this' and 'M'.
  // The rings must be the same.
  void swap(DMat<ACoeffRing>& M);

  const ACoeffRing& ring() const;
  size_t numRows() const;
  size_t numColumns() const;

  ElementType& entry(size_t row, size_t column);
  const ElementType& entry(size_t row, size_t column) const;

  void resize(size_t new_nrows, size_t new_ncols);

  // These are labelled 'unsafe', as it s possible the rows
  // are out of order (which happens in particular if
  // certain flint functions created this.
  const ElementType* unsafeArray() const;
  ElementType*& unsafeArray();
 public:
  // Access routines so that the flint nmod_mat interface may be used
  const nmod_mat_t& nmod_mat() const;
  nmod_mat_t& nmod_mat();
 private:
  const ACoeffRing* mRing;
  nmod_mat_t mArray;
};

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
