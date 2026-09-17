// Copyright 2026  The Macaulay2 Authors

#include "exceptions.hpp"
#include "basic-mutable-matrices/dmat-linalg-qq-gmp.hpp"
#include "basic-mutable-matrices/dmat-qq-interface-flint.hpp"
#include "basic-rings/aring-QQ.hpp"
#include "basic-mutable-matrices/dmat-linalg.hpp"

using DMatQQ = DMat<M2::ARingQQ>;

namespace MatrixOps {
//////////////////////
// QQ ////////////////
//////////////////////

template <>
void mult(const DMatQQ& A, const DMatQQ& B, DMatQQ& result_product)
{
  FlintQQMat A1(A);
  FlintQQMat B1(B);
  FlintQQMat result1(A.numRows(), B.numColumns());

  fmpq_mat_mul(result1.value(), A1.value(), B1.value());

  result1.toDMat(result_product);
}

template <>
void addMultipleTo(DMatQQ& C, const DMatQQ& A, const DMatQQ& B)
{
  FlintQQMat A1(A);
  FlintQQMat B1(B);
  FlintQQMat C1(C);
  FlintQQMat result1(A.numRows(), B.numColumns());

  FlintQQMat D1(A.numRows(), B.numColumns());
  fmpq_mat_mul(D1.value(), A1.value(), B1.value());
  fmpq_mat_add(C1.value(), C1.value(), D1.value());

  C1.toDMat(C);
}

template <>
void subtractMultipleTo(DMatQQ& C, const DMatQQ& A, const DMatQQ& B)
{
  FlintQQMat A1(A);
  FlintQQMat B1(B);
  FlintQQMat C1(C);
  FlintQQMat result1(A.numRows(), B.numColumns());

  FlintQQMat D1(A.numRows(), B.numColumns());
  fmpq_mat_mul(D1.value(), A1.value(), B1.value());
  fmpq_mat_sub(C1.value(), C1.value(), D1.value());

  C1.toDMat(C);
}

template <>
size_t rowReducedEchelonForm(const DMatQQ& A, DMatQQ& result_rref)
{
  FlintQQMat A1(A);
  long rank = fmpq_mat_rref(A1.value(), A1.value());
  A1.toDMat(result_rref);
  return rank;
}


}  // namespace MatrixOps
