// Copyright 2014  Michael E. Stillman

#ifndef M2_BASIC_MUTMAT_DMAT_GF_FLINT__HPP_
#define M2_BASIC_MUTMAT_DMAT_GF_FLINT__HPP_

#include "basic-mutable-matrices/dmat.hpp"

// The following needs to be included before any flint files are included.
#include <M2/gc-include.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <flint/fq_nmod_mat.h>  // for fq_zech_mat_entry, fq_zech_mat_clear
#include <flint/fq_zech_mat.h>  // for fq_zech_mat_t
#pragma GCC diagnostic pop

/////////////////////////////////////////////////////////////////
// Flint: use fq_nmod_mat for implementation of dense matrices //
/////////////////////////////////////////////////////////////////
template <>
class DMat<M2::ARingGFFlint>
{
 public:
  typedef M2::ARingGFFlint ACoeffRing;
  typedef ACoeffRing CoeffRing;
  typedef fq_zech_struct ElementType;
  // typedef ElementType elem;
  // typedef ACoeffRing::Element Element;

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
  // Access routines so that the flint fq_zech_mat interface may be used
  const fq_zech_mat_t& fq_zech_mat() const;
  fq_zech_mat_t& fq_zech_mat();
 private:
  const ACoeffRing* mRing;
  fq_zech_mat_t mArray;
};

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
