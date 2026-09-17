// Copyright 2013  Michael E. Stillman

#include "basic-mutable-matrices/dmat-zzp-flint.hpp"

#include "basic-rings/aring-ZZp-flint.hpp"

#include <utility>

DMat<M2::ARingZZpFlint>::DMat() : mRing(nullptr) {}

DMat<M2::ARingZZpFlint>::DMat(const ACoeffRing& R,
                              size_t nrows,
                              size_t ncols)
    : mRing(&R)
{
  nmod_mat_init(mArray, nrows, ncols, R.characteristic());
}

DMat<M2::ARingZZpFlint>::DMat(const DMat<ACoeffRing>& M)
    : mRing(&M.ring())
{
  nmod_mat_init_set(mArray, M.mArray);
}

DMat<M2::ARingZZpFlint>::~DMat() { nmod_mat_clear(mArray); }

void DMat<M2::ARingZZpFlint>::swap(DMat<ACoeffRing>& M)
{
  std::swap(mRing, M.mRing);
  std::swap(mArray, M.mArray);
}

const DMat<M2::ARingZZpFlint>::ACoeffRing&
DMat<M2::ARingZZpFlint>::ring() const
{
  return *mRing;
}

size_t DMat<M2::ARingZZpFlint>::numRows() const
{
  return nmod_mat_nrows(mArray);
}

size_t DMat<M2::ARingZZpFlint>::numColumns() const
{
  return nmod_mat_ncols(mArray);
}

DMat<M2::ARingZZpFlint>::ElementType& DMat<M2::ARingZZpFlint>::entry(
    size_t row,
    size_t column)
{
  return nmod_mat_entry(mArray, row, column);
}

const DMat<M2::ARingZZpFlint>::ElementType&
DMat<M2::ARingZZpFlint>::entry(size_t row, size_t column) const
{
  return nmod_mat_entry(mArray, row, column);
}

void DMat<M2::ARingZZpFlint>::resize(size_t new_nrows, size_t new_ncols)
{
  DMat newMatrix(ring(), new_nrows, new_ncols);
  swap(newMatrix);
}

const DMat<M2::ARingZZpFlint>::ElementType*
DMat<M2::ARingZZpFlint>::unsafeArray() const
{
  return mArray->entries;
}

DMat<M2::ARingZZpFlint>::ElementType*&
DMat<M2::ARingZZpFlint>::unsafeArray()
{
  return mArray->entries;
}

const nmod_mat_t& DMat<M2::ARingZZpFlint>::nmod_mat() const
{
  return mArray;
}

nmod_mat_t& DMat<M2::ARingZZpFlint>::nmod_mat() { return mArray; }

