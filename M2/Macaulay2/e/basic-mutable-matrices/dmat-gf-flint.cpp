// Copyright 2014  Michael E. Stillman

#include "basic-mutable-matrices/dmat-gf-flint.hpp"

#include "basic-rings/aring-GF-flint.hpp"

#include <utility>

DMat<M2::ARingGFFlint>::DMat() : mRing(nullptr) {}

DMat<M2::ARingGFFlint>::DMat(const ACoeffRing& R,
                             size_t nrows,
                             size_t ncols)
    : mRing(&R)
{
  fq_zech_mat_init(mArray, nrows, ncols, ring().flintContext());
}

DMat<M2::ARingGFFlint>::DMat(const DMat<ACoeffRing>& M)
    : mRing(&M.ring())
{
  fq_zech_mat_init_set(mArray, M.mArray, ring().flintContext());
}

DMat<M2::ARingGFFlint>::~DMat()
{
  fq_zech_mat_clear(mArray, ring().flintContext());
}

void DMat<M2::ARingGFFlint>::swap(DMat<ACoeffRing>& M)
{
  std::swap(mRing, M.mRing);
  std::swap(mArray, M.mArray);
}

const DMat<M2::ARingGFFlint>::ACoeffRing&
DMat<M2::ARingGFFlint>::ring() const
{
  return *mRing;
}

size_t DMat<M2::ARingGFFlint>::numRows() const
{
  return fq_zech_mat_nrows(mArray, ring().flintContext());
}

size_t DMat<M2::ARingGFFlint>::numColumns() const
{
  return fq_zech_mat_ncols(mArray, ring().flintContext());
}

DMat<M2::ARingGFFlint>::ElementType& DMat<M2::ARingGFFlint>::entry(
    size_t row,
    size_t column)
{
  return *fq_zech_mat_entry(mArray, row, column);
}

const DMat<M2::ARingGFFlint>::ElementType&
DMat<M2::ARingGFFlint>::entry(size_t row, size_t column) const
{
  return *fq_zech_mat_entry(mArray, row, column);
}

void DMat<M2::ARingGFFlint>::resize(size_t new_nrows, size_t new_ncols)
{
  DMat newMatrix(ring(), new_nrows, new_ncols);
  swap(newMatrix);
}

const DMat<M2::ARingGFFlint>::ElementType*
DMat<M2::ARingGFFlint>::unsafeArray() const
{
  return mArray->entries;
}

DMat<M2::ARingGFFlint>::ElementType*&
DMat<M2::ARingGFFlint>::unsafeArray()
{
  return mArray->entries;
}

const fq_zech_mat_t& DMat<M2::ARingGFFlint>::fq_zech_mat() const
{
  return mArray;
}

fq_zech_mat_t& DMat<M2::ARingGFFlint>::fq_zech_mat() { return mArray; }

