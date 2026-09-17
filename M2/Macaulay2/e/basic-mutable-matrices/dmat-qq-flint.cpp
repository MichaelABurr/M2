// Copyright 2013  Michael E. Stillman

#include "basic-mutable-matrices/dmat-qq-flint.hpp"

#include "basic-rings/aring-QQ-flint.hpp"

#include <cassert>
#include <utility>

DMat<M2::ARingQQFlint>::DMat() : mRing(nullptr) {}

DMat<M2::ARingQQFlint>::DMat(const ACoeffRing& R,
                             size_t nrows,
                             size_t ncols)
    : mRing(&R)
{
  fmpq_mat_init(mArray, nrows, ncols);
}

DMat<M2::ARingQQFlint>::DMat(const DMat<ACoeffRing>& M)
    : mRing(&M.ring())
{
  fmpq_mat_init(mArray, M.numRows(), M.numColumns());
  fmpq_mat_set(mArray, M.mArray);
}

DMat<M2::ARingQQFlint>::~DMat() { fmpq_mat_clear(mArray); }

void DMat<M2::ARingQQFlint>::swap(DMat<ACoeffRing>& M)
{
  std::swap(mRing, M.mRing);
  std::swap(*mArray, *M.mArray);
}

const DMat<M2::ARingQQFlint>::ACoeffRing&
DMat<M2::ARingQQFlint>::ring() const
{
  return *mRing;
}

size_t DMat<M2::ARingQQFlint>::numRows() const
{
  return fmpq_mat_nrows(mArray);
}

size_t DMat<M2::ARingQQFlint>::numColumns() const
{
  return fmpq_mat_ncols(mArray);
}

DMat<M2::ARingQQFlint>::ElementType& DMat<M2::ARingQQFlint>::entry(
    size_t row,
    size_t column)
{
  assert(row < numRows());
  assert(column < numColumns());
  return *fmpq_mat_entry(mArray, row, column);
}

const DMat<M2::ARingQQFlint>::ElementType&
DMat<M2::ARingQQFlint>::entry(size_t row, size_t column) const
{
  assert(row < numRows());
  assert(column < numColumns());
  return *fmpq_mat_entry(mArray, row, column);
}

void DMat<M2::ARingQQFlint>::resize(size_t new_nrows, size_t new_ncols)
{
  DMat newMatrix(ring(), new_nrows, new_ncols);
  swap(newMatrix);
}

const DMat<M2::ARingQQFlint>::ElementType*
DMat<M2::ARingQQFlint>::unsafeArray() const
{
  return mArray->entries;
}

DMat<M2::ARingQQFlint>::ElementType*&
DMat<M2::ARingQQFlint>::unsafeArray()
{
  return mArray->entries;
}

const fmpq_mat_t& DMat<M2::ARingQQFlint>::fmpq_mat() const
{
  return mArray;
}

fmpq_mat_t& DMat<M2::ARingQQFlint>::fmpq_mat() { return mArray; }

