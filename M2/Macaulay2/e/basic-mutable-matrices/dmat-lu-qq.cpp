// Copyright 2026  The Macaulay2 Authors

#include "basic-mutable-matrices/dmat-lu.hpp"
#include "basic-mutable-matrices/dmat-lu-qq.hpp"

#include "basic-mutable-matrices/dmat-qq-interface-flint.hpp"
#include "basic-rings/aring-QQ.hpp"

namespace {
using RingType = M2::ARingQQGMP;
using Mat = DMat<RingType>;

void findColumnRankProfileFromLU(const fmpq_mat_t B,
                                 std::vector<size_t>& profile)
{
  profile.clear();
  long nrows = fmpq_mat_nrows(B);
  long ncols = fmpq_mat_ncols(B);
  long thiscol = 0;
  long thisrow = 0;
  while (thisrow < nrows and thiscol < ncols)
    {
      if (not fmpq_is_zero(fmpq_mat_entry(B, thisrow, thiscol)))
        {
          profile.push_back(thiscol);
          thisrow++;
        }
      thiscol++;
    }
}
}  // namespace

DMatLinAlg<RingType>::DMatLinAlg(const Mat& A) : mInputMatrix(A) {}

size_t DMatLinAlg<RingType>::rank()
{
  FlintQQMat A(mInputMatrix);
  fmpz_mat_t m1;
  fmpz_mat_init(m1, A.numRows(), A.numColumns());
  fmpq_mat_get_fmpz_mat_rowwise(m1, nullptr, A.value());
  size_t result = fmpz_mat_rank(m1);
  fmpz_mat_clear(m1);
  return result;
}

void DMatLinAlg<RingType>::determinant(__mpq_struct& result_det)
{
  FlintQQMat A(mInputMatrix);
  fmpq_t det;
  fmpq_init(det);
  fmpq_mat_det(det, A.value());
  fmpq_get_mpq(&result_det, det);
  fmpq_clear(det);
}

void DMatLinAlg<RingType>::columnRankProfile(std::vector<size_t>& profile)
{
  FlintQQMat A(mInputMatrix);
  FlintQQMat B(A.numRows(), A.numColumns());
  fmpq_mat_rref(B.value(), A.value());
  findColumnRankProfileFromLU(B.value(), profile);
}

void DMatLinAlg<RingType>::matrixPLU(std::vector<size_t>& P,
                                     Mat& L,
                                     Mat& U)
{
  DMatLUinPlace<RingType> luObject(mInputMatrix);
  const Mat& LU = luObject.LUinPlace();
  LUUtil<RingType>::setUpperLower(LU, L, U);
  P = luObject.permutation();
}

bool DMatLinAlg<RingType>::solve(const Mat& B1, Mat& X1)
{
  long nrows = mInputMatrix.numRows();
  long ncols = mInputMatrix.numColumns();
  std::vector<size_t> profile;
  Mat AB1(mInputMatrix.ring(), nrows, ncols + B1.numColumns());
  concatenateMatrices<Mat>(mInputMatrix, B1, AB1);
  FlintQQMat AB(AB1);
  fmpq_mat_rref(AB.value(), AB.value());
  findColumnRankProfileFromLU(AB.value(), profile);
  if (profile[profile.size() - 1] >= ncols)
    return false;

  FlintQQMat X(ncols, B1.numColumns());
  for (long c = 0; c < B1.numColumns(); c++)
    for (long r = 0; r < profile.size(); r++)
      fmpq_set(fmpq_mat_entry(X.value(), profile[r], c),
               fmpq_mat_entry(AB.value(), r, ncols + c));
  X.toDMat(X1);
  return true;
}

bool DMatLinAlg<RingType>::solveInvertible(const Mat& B1, Mat& X1)
{
  FlintQQMat A(mInputMatrix);
  FlintQQMat B(B1);
  FlintQQMat X(mInputMatrix.numColumns(), B1.numColumns());
  int isfullrank = fmpq_mat_solve_dixon(X.value(), A.value(), B.value());
  if (isfullrank == 0) return false;
  X.toDMat(X1);
  return true;
}

bool DMatLinAlg<RingType>::inverse(Mat& result_inv)
{
  FlintQQMat A(mInputMatrix);
  long nrows = A.numRows();
  FlintQQMat inv(nrows, nrows);
  fmpz_t D, den;
  fmpz_mat_t matZZ, invZZ;
  fmpz_init(D);
  fmpz_init(den);
  fmpz_mat_init(matZZ, nrows, nrows);
  fmpz_mat_init(invZZ, nrows, nrows);

  fmpq_mat_get_fmpz_mat_matwise(matZZ, D, A.value());
  int result = fmpz_mat_inv(invZZ, den, matZZ);
  if (result != 0)
    {
      fmpq_mat_set_fmpz_mat_div_fmpz(inv.value(), invZZ, den);
      fmpq_mat_scalar_mul_fmpz(inv.value(), inv.value(), D);
      inv.toDMat(result_inv);
    }
  fmpz_clear(D);
  fmpz_clear(den);
  fmpz_mat_clear(matZZ);
  fmpz_mat_clear(invZZ);
  return result != 0;
}

size_t DMatLinAlg<RingType>::kernel(Mat& result_nullspace)
{
  FlintQQMat A(mInputMatrix);
  fmpz_mat_t matZZ;
  fmpz_mat_t kerZZ;
  fmpz_mat_init(matZZ, A.numRows(), A.numColumns());
  fmpz_mat_init(kerZZ, A.numColumns(), A.numColumns());

  fmpq_mat_get_fmpz_mat_rowwise(matZZ, nullptr, A.value());
  long nullity = fmpz_mat_nullspace(kerZZ, matZZ);
  FlintQQMat kerQQ(A.numColumns(), nullity);
  for (long r = 0; r < A.numColumns(); r++)
    for (long c = 0; c < nullity; c++)
      kerQQ.set_from_fmpz(r, c, fmpz_mat_entry(kerZZ, r, c));
  kerQQ.toDMat(result_nullspace);
  fmpz_mat_clear(matZZ);
  fmpz_mat_clear(kerZZ);
  return nullity;
}
