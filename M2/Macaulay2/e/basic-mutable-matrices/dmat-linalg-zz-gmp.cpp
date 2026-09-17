// Copyright 2026  The Macaulay2 Authors

#include "exceptions.hpp"
#include "basic-mutable-matrices/dmat-linalg-zz-gmp.hpp"
#include "basic-mutable-matrices/dmat-qq-interface-flint.hpp"
#include "basic-mutable-matrices/dmat-linalg.hpp"
namespace MatrixOps {
//////////////////////
// ZZ (ARingZZGMP) ///
//////////////////////

template <>
M2_arrayintOrNull LU(const DMatZZGMP& A, DMatZZGMP& L, DMatZZGMP& U)
{
  (void) A;
  (void) L;
  (void) U;
  throw exc::engine_error(
      "'LU' not implemented for this kind of matrix over this ring");
}

template <>
M2_arrayintOrNull rankProfile(const DMatZZGMP& A, bool row_profile)
{
  (void) A;
  (void) row_profile;
  throw exc::engine_error(
      "'rankProfile' not implemented for this kind of matrix over this ring");
}

template <>
bool inverse(const DMatZZGMP& A, DMatZZGMP& result_inv)
{
  (void) A;
  (void) result_inv;
  throw exc::engine_error(
      "'invert' not implemented for this kind of matrix over this ring");
}

template <>
size_t nullSpace(const DMatZZGMP& A, DMatZZGMP& result_nullspace)
{
  (void) A;
  (void) result_nullspace;
  throw exc::engine_error(
      "'nullSpace' not implemented for this kind of matrix over this ring");
}

template <>
bool solveLinear(const DMatZZGMP& A, const DMatZZGMP& B, DMatZZGMP& X)
{
  (void) A;
  (void) B;
  (void) X;
  throw exc::engine_error(
      "'solveLinear' not implemented for this kind of matrix over this ring");
}

template <>
bool solveInvertible(const DMatZZGMP& A,
                            const DMatZZGMP& B,
                            DMatZZGMP& X)
{
  (void) A;
  (void) B;
  (void) X;
  throw exc::engine_error(
      "'solveInvertible' not implemented for this kind of matrix over this "
      "ring");
}

template <>
void mult(const DMatZZGMP& A,
                 const DMatZZGMP& B,
                 DMatZZGMP& result_product)
{
  FlintZZMat A1(A);
  FlintZZMat B1(B);
  FlintZZMat result1(A.numRows(), B.numColumns());

  fmpz_mat_mul(result1.value(), A1.value(), B1.value());

  result1.toDMat(result_product);
}

template <>
void addMultipleTo(DMatZZGMP& C, const DMatZZGMP& A, const DMatZZGMP& B)
{
  FlintZZMat A1(A);
  FlintZZMat B1(B);
  FlintZZMat C1(C);
  FlintZZMat result1(A.numRows(), B.numColumns());

  FlintZZMat D1(A.numRows(), B.numColumns());
  fmpz_mat_mul(D1.value(), A1.value(), B1.value());
  fmpz_mat_add(C1.value(), C1.value(), D1.value());

  C1.toDMat(C);
}

template <>
void subtractMultipleTo(DMatZZGMP& C,
                               const DMatZZGMP& A,
                               const DMatZZGMP& B)
{
  FlintZZMat A1(A);
  FlintZZMat B1(B);
  FlintZZMat C1(C);
  FlintZZMat result1(A.numRows(), B.numColumns());

  FlintZZMat D1(A.numRows(), B.numColumns());
  fmpz_mat_mul(D1.value(), A1.value(), B1.value());
  fmpz_mat_sub(C1.value(), C1.value(), D1.value());

  C1.toDMat(C);
}

template <>
size_t rank(const DMatZZGMP& A)
{
  FlintZZMat A1(A);
  return fmpz_mat_rank(A1.value());
}

template <>
void determinant<M2::ARingZZGMP, M2::ARingZZGMP::ElementType>(
    const DMatZZGMP& A,
    M2::ARingZZGMP::ElementType& result_det)
{
  FlintZZMat A1(A);
  fmpz_t det;
  fmpz_init(det);
  fmpz_mat_det(det, A1.value());
  fmpz_get_mpz(&result_det, det);
  fmpz_clear(det);
}


}  // namespace MatrixOps
