// Copyright 2013  Michael E. Stillman

#ifndef M2_BASIC_MUTMAT_DMAT_QQ_FLINT_HPP_
#define M2_BASIC_MUTMAT_DMAT_QQ_FLINT_HPP_

#include "basic-mutable-matrices/dmat.hpp"

// The following needs to be included before any flint files are included.
#include <M2/gc-include.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <flint/fmpq_mat.h>  // for fmpq_mat_t, fmpq_mat_entry, fmpq_mat_init, fmpq_ma...
#pragma GCC diagnostic pop

////////////////////////////////////////////////////
// Dense matrices for Flint type M2::ARingQQFlint //
////////////////////////////////////////////////////

// Warning: objects of this class should *not* go to the front end.
// fmpz_t's might be garbage collected out from under you...
template <>
class DMat<M2::ARingQQFlint>
{
 public:
  typedef M2::ARingQQFlint ACoeffRing;
  typedef ACoeffRing CoeffRing;
  typedef fmpq ElementType;
  // typedef ElementType elem;
  // typedef ACoeffRing::Element Element;

  DMat();
  DMat(const ACoeffRing& R, size_t nrows, size_t ncols);
  DMat(const DMat<ACoeffRing>& M);
  ~DMat();

    // swap the actual matrices of 'this' and 'M'.
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
  // Other routines from flint nmod_mat interface
  const fmpq_mat_t& fmpq_mat() const;
  fmpq_mat_t& fmpq_mat();
 private:
  const ACoeffRing* mRing;
  fmpq_mat_t mArray;
};

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
