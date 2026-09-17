// Copyright 2013  Michael E. Stillman

#include "basic-mutable-matrices/dmat-zz-flint.hpp"

#include "basic-rings/aring-ZZ-flint.hpp"

#include <cassert>
#include <utility>

DMat<M2::ARingZZ>::DMat() : mRing(nullptr) {}

DMat<M2::ARingZZ>::DMat(const ACoeffRing& R,
                        size_t nrows,
                        size_t ncols)
    : mRing(&R)
{
  fmpz_mat_init(mArray, nrows, ncols);
}

DMat<M2::ARingZZ>::DMat(const DMat<ACoeffRing>& M) : mRing(&M.ring())
{
  fmpz_mat_init_set(mArray, M.mArray);
}

DMat<M2::ARingZZ>::~DMat() { fmpz_mat_clear(mArray); }

void DMat<M2::ARingZZ>::swap(DMat<ACoeffRing>& M)
{
  std::swap(mRing, M.mRing);
  fmpz_mat_swap(mArray, M.mArray);
}

const DMat<M2::ARingZZ>::ACoeffRing& DMat<M2::ARingZZ>::ring() const
{
  return *mRing;
}

size_t DMat<M2::ARingZZ>::numRows() const { return fmpz_mat_nrows(mArray); }

size_t DMat<M2::ARingZZ>::numColumns() const
{
  return fmpz_mat_ncols(mArray);
}

DMat<M2::ARingZZ>::ElementType& DMat<M2::ARingZZ>::entry(size_t row,
                                                         size_t column)
{
  assert(row < numRows());
  assert(column < numColumns());
  return *fmpz_mat_entry(mArray, row, column);
}

const DMat<M2::ARingZZ>::ElementType& DMat<M2::ARingZZ>::entry(
    size_t row,
    size_t column) const
{
  assert(row < numRows());
  assert(column < numColumns());
  return *fmpz_mat_entry(mArray, row, column);
}

void DMat<M2::ARingZZ>::resize(size_t new_nrows, size_t new_ncols)
{
  DMat newMatrix(ring(), new_nrows, new_ncols);
  swap(newMatrix);
}

const DMat<M2::ARingZZ>::ElementType* DMat<M2::ARingZZ>::unsafeArray() const
{
  return mArray->entries;
}

DMat<M2::ARingZZ>::ElementType*& DMat<M2::ARingZZ>::unsafeArray()
{
  return mArray->entries;
}

const fmpz_mat_t& DMat<M2::ARingZZ>::fmpz_mat() const { return mArray; }

fmpz_mat_t& DMat<M2::ARingZZ>::fmpz_mat() { return mArray; }

