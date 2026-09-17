// Copyright 2026  The Macaulay2 Authors

#include "exceptions.hpp"
#include "basic-mutable-matrices/dmat-linalg-zz-flint.hpp"
#include "basic-mutable-matrices/dmat-zz-flint.hpp"
#include "basic-rings/aring-ZZ-flint.hpp"
#include "basic-mutable-matrices/dmat-linalg.hpp"

using DMatZZ = DMat<M2::ARingZZ>;

namespace MatrixOps {
//////////////////////
// ZZFlint ///////////
//////////////////////
// Warning: nullSpace is WRONG, and needs to be rewritten,
// using an algorithm that will compute kernel over ZZ.

// Functions for DMatZZ

template <>
size_t rank(const DMatZZ& A) { return fmpz_mat_rank(A.fmpz_mat()); }
template <>
void determinant<M2::ARingZZ, M2::ARingZZ::ElementType>(
    const DMatZZ& A,
    M2::ARingZZ::ElementType& result_det)
{
  fmpz_mat_det(&result_det, A.fmpz_mat());
}

template <>
bool inverse(const DMatZZ& A, DMatZZ& result_inv)
{
  M2::ARingZZ::Element den(A.ring());
  bool result = fmpz_mat_inv(result_inv.fmpz_mat(), &den.value(), A.fmpz_mat());
  if (!fmpz_is_pm1(&den.value())) result = false;
  return result;
}

template <>
void mult(const DMatZZ& A, const DMatZZ& B, DMatZZ& result_product)
{
  // The A1 and B1 on the next line are switched because the memory layout
  // expected
  // is the transpose of what we have for DMat.
  fmpz_mat_mul(result_product.fmpz_mat(), A.fmpz_mat(), B.fmpz_mat());
}

template <>
size_t nullSpace(const DMatZZ& A, DMatZZ& result_nullspace)
{
  long nullity = fmpz_mat_nullspace(result_nullspace.fmpz_mat(), A.fmpz_mat());
  return nullity;
}

template <>
bool solveLinear(const DMatZZ& A, const DMatZZ& B, DMatZZ& X)
{
  M2::ARingZZ::Element den(A.ring());
  bool result = fmpz_mat_solve(X.fmpz_mat(), &den.value(), B.fmpz_mat(), A.fmpz_mat());
  if (!fmpz_is_pm1(&den.value())) result = false;
  return result;
}

template <>
M2_arrayintOrNull rankProfile(const DMatZZ& A, bool row_profile)
{
  (void) A;
  (void) row_profile;
  throw exc::engine_error(
      "'rankProfile' not implemented for this kind of matrix over this ring");
}

template <>
void addMultipleTo(DMatZZ& C, const DMatZZ& A, const DMatZZ& B)
{
  DMatZZ D(C.ring(), A.numRows(), B.numColumns());
  fmpz_mat_mul(D.fmpz_mat(), A.fmpz_mat(), B.fmpz_mat());
  fmpz_mat_add(C.fmpz_mat(), C.fmpz_mat(), D.fmpz_mat());
}

template <>
void subtractMultipleTo(DMatZZ& C, const DMatZZ& A, const DMatZZ& B)
{
  DMatZZ D(C.ring(), A.numRows(), B.numColumns());
  fmpz_mat_mul(D.fmpz_mat(), A.fmpz_mat(), B.fmpz_mat());
  fmpz_mat_sub(C.fmpz_mat(), C.fmpz_mat(), D.fmpz_mat());
}


}  // namespace MatrixOps
