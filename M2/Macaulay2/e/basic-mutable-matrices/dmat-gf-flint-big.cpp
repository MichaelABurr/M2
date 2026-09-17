// Copyright 2014  Michael E. Stillman

#include "basic-mutable-matrices/dmat-gf-flint-big.hpp"

#include "basic-rings/aring-GF-flint-big.hpp"

#include <utility>

DMat<M2::ARingGFFlintBig>::DMat() : mRing(nullptr) {}

DMat<M2::ARingGFFlintBig>::DMat(const ACoeffRing& R,
                                size_t nrows,
                                size_t ncols)
    : mRing(&R)
{
  fq_nmod_mat_init(mArray, nrows, ncols, ring().flintContext());
}

DMat<M2::ARingGFFlintBig>::DMat(const DMat<ACoeffRing>& M)
    : mRing(&M.ring())
{
  fq_nmod_mat_init_set(mArray, M.mArray, ring().flintContext());
}

DMat<M2::ARingGFFlintBig>::~DMat()
{
  fq_nmod_mat_clear(mArray, ring().flintContext());
}

void DMat<M2::ARingGFFlintBig>::swap(DMat<ACoeffRing>& M)
{
  std::swap(mRing, M.mRing);
  std::swap(mArray, M.mArray);
}

const DMat<M2::ARingGFFlintBig>::ACoeffRing&
DMat<M2::ARingGFFlintBig>::ring() const
{
  return *mRing;
}

size_t DMat<M2::ARingGFFlintBig>::numRows() const
{
  return fq_nmod_mat_nrows(mArray, ring().flintContext());
}

size_t DMat<M2::ARingGFFlintBig>::numColumns() const
{
  return fq_nmod_mat_ncols(mArray, ring().flintContext());
}

DMat<M2::ARingGFFlintBig>::ElementType&
DMat<M2::ARingGFFlintBig>::entry(size_t row, size_t column)
{
  return *fq_nmod_mat_entry(mArray, row, column);
}

const DMat<M2::ARingGFFlintBig>::ElementType&
DMat<M2::ARingGFFlintBig>::entry(size_t row, size_t column) const
{
  return *fq_nmod_mat_entry(mArray, row, column);
}

void DMat<M2::ARingGFFlintBig>::resize(size_t new_nrows,
                                       size_t new_ncols)
{
  DMat newMatrix(ring(), new_nrows, new_ncols);
  swap(newMatrix);
}

const DMat<M2::ARingGFFlintBig>::ElementType*
DMat<M2::ARingGFFlintBig>::unsafeArray() const
{
  return mArray->entries;
}

DMat<M2::ARingGFFlintBig>::ElementType*&
DMat<M2::ARingGFFlintBig>::unsafeArray()
{
  return mArray->entries;
}

const fq_nmod_mat_t& DMat<M2::ARingGFFlintBig>::fq_nmod_mat() const
{
  return mArray;
}

fq_nmod_mat_t& DMat<M2::ARingGFFlintBig>::fq_nmod_mat() { return mArray; }

