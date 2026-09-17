// Copyright 2026  The Macaulay2 Authors

#include "basic-mutable-matrices/dmat-lu.hpp"
#include "basic-mutable-matrices/dmat-lu-zzp-flint.hpp"

#include "basic-mutable-matrices/dmat-zzp-flint.hpp"
#include "basic-rings/aring-ZZp-flint.hpp"

#include <cassert>

namespace {
using RingType = M2::ARingZZpFlint;
using Mat = DMat<RingType>;
}  // namespace

DMatLinAlg<RingType>::DMatLinAlg(const Mat& A) : mMatrix(A) {}

size_t DMatLinAlg<RingType>::rank()
{
  return nmod_mat_rank(mMatrix.nmod_mat());
}

void DMatLinAlg<RingType>::determinant(mp_limb_t& result_det)
{
  result_det = nmod_mat_det(mMatrix.nmod_mat());
}

void DMatLinAlg<RingType>::columnRankProfile(std::vector<size_t>& profile)
{
  Mat LU(mMatrix);
  mp_limb_signed_t* perm = newarray_atomic(mp_limb_signed_t, LU.numRows());
  nmod_mat_lu(perm, LU.nmod_mat(), false);
  freemem(perm);
  LUUtil<RingType>::computePivotColumns(LU, profile);
}

void DMatLinAlg<RingType>::matrixPLU(std::vector<size_t>& P,
                                     Mat& L,
                                     Mat& U)
{
  Mat LU(mMatrix);
  mp_limb_signed_t* perm = newarray_atomic(mp_limb_signed_t, LU.numRows());
  nmod_mat_lu(perm, LU.nmod_mat(), false);
  P.clear();
  for (long i = 0; i < LU.numRows(); i++) P.push_back(perm[i]);
  freemem(perm);
  LUUtil<RingType>::setUpperLower(LU, L, U);
}

bool DMatLinAlg<RingType>::solveInvertible(const Mat& B, Mat& X)
{
  assert(mMatrix.numRows() == mMatrix.numColumns());
  assert(mMatrix.numRows() == B.numRows());
  X.resize(mMatrix.numColumns(), B.numColumns());
  return nmod_mat_solve(X.nmod_mat(), mMatrix.nmod_mat(), B.nmod_mat()) != 0;
}

bool DMatLinAlg<RingType>::solve(const Mat& B, Mat& X)
{
  long nrows = mMatrix.numRows();
  long ncols = mMatrix.numColumns();
  std::vector<size_t> profile;
  Mat AB(mMatrix.ring(), nrows, ncols + B.numColumns());
  concatenateMatrices<Mat>(mMatrix, B, AB);
  nmod_mat_rref(AB.nmod_mat());
  LUUtil<RingType>::computePivotColumns(AB, profile);
  if (profile.size() >= 1 and profile[profile.size() - 1] >= ncols)
    return false;

  X.resize(ncols, B.numColumns());
  for (long c = 0; c < B.numColumns(); c++)
    for (long r = 0; r < profile.size(); r++)
      mMatrix.ring().copy(X.entry(profile[r], c), AB.entry(r, ncols + c));
  return true;
}

bool DMatLinAlg<RingType>::inverse(Mat& result_inv)
{
  assert(mMatrix.numRows() == mMatrix.numColumns());
  Mat& A = const_cast<Mat&>(mMatrix);
  result_inv.resize(mMatrix.numRows(), mMatrix.numColumns());
  return nmod_mat_inv(result_inv.nmod_mat(), A.nmod_mat()) != 0;
}

size_t DMatLinAlg<RingType>::kernel(Mat& result_nullspace)
{
  long rank = nmod_mat_rank(mMatrix.nmod_mat());
  result_nullspace.resize(mMatrix.numColumns(), mMatrix.numColumns() - rank);
  long nullity =
      nmod_mat_nullspace(result_nullspace.nmod_mat(), mMatrix.nmod_mat());
  assert(rank == mMatrix.numColumns() - nullity);
  return nullity;
}
