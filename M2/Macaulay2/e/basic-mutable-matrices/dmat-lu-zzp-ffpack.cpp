// Copyright 2026  The Macaulay2 Authors

#include "basic-mutable-matrices/dmat-lu.hpp"
#include "basic-mutable-matrices/dmat-lu-zzp-ffpack.hpp"

#include "basic-rings/aring-ZZp-ffpack.hpp"

namespace {
using RingType = M2::ARingZZpFFPACK;
using Mat = DMat<RingType>;
}  // namespace

namespace ffpackInterface {
size_t rank(const Mat& A);
void determinant(const Mat& A, double& result_det);
void rankProfile(const Mat& A,
                 bool row_profile,
                 std::vector<size_t>& result_profile);
bool solveLinear(const Mat& A, const Mat& B, Mat& X);
bool inverse(const Mat& A, Mat& result_inv);
size_t nullSpace(const Mat& A, Mat& result_nullspace);
}  // namespace ffpackInterface

DMatLinAlg<RingType>::DMatLinAlg(const Mat& A) : mLU(A) {}

size_t DMatLinAlg<RingType>::rank()
{
  return ffpackInterface::rank(mLU);
}

void DMatLinAlg<RingType>::determinant(double& result)
{
  ffpackInterface::determinant(mLU, result);
}

void DMatLinAlg<RingType>::columnRankProfile(std::vector<size_t>& profile)
{
  ffpackInterface::rankProfile(mLU, false, profile);
}

void DMatLinAlg<RingType>::matrixPLU(std::vector<size_t>& P,
                                     Mat& L,
                                     Mat& U)
{
  DMatLUinPlace<RingType> C(mLU);
  const Mat& LU = C.LUinPlace();
  LUUtil<RingType>::setUpperLower(LU, L, U);
  P = C.permutation();
}

bool DMatLinAlg<RingType>::solve(const Mat& B, Mat& X)
{
  return ffpackInterface::solveLinear(mLU, B, X);
}

bool DMatLinAlg<RingType>::solveInvertible(const Mat& B, Mat& X)
{
  return ffpackInterface::solveLinear(mLU, B, X);
}

bool DMatLinAlg<RingType>::inverse(Mat& X)
{
  return ffpackInterface::inverse(mLU, X);
}

size_t DMatLinAlg<RingType>::kernel(Mat& X)
{
  return ffpackInterface::nullSpace(mLU, X);
}
