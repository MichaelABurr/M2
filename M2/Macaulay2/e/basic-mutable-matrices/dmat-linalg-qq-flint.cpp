// Copyright 2026  The Macaulay2 Authors

#include "exceptions.hpp"
#include "basic-mutable-matrices/dmat-linalg-qq-flint.hpp"
#include "basic-mutable-matrices/dmat-qq-flint.hpp"
#include "basic-rings/aring-QQ-flint.hpp"
#include "basic-mutable-matrices/dmat-linalg.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <flint/fmpz_mat.h>
#pragma GCC diagnostic pop

using DMatQQFlint = DMat<M2::ARingQQFlint>;

namespace MatrixOps {
//////////////////////
// QQFlint ///////////
//////////////////////

// Functions for DMatQQFlint

template <>
size_t rank(const DMatQQFlint& A)
{
  // fmpq_mat has no rank function.
  // So we clear denominators row-wise (or column-wise), and compute the rank of
  // that matrix.
  fmpz_mat_t m1;
  fmpz_mat_init(m1, A.numRows(), A.numColumns());
  fmpq_mat_get_fmpz_mat_rowwise(m1, nullptr, A.fmpq_mat());
  // fmpz_mat_print_pretty(m1);
  size_t rk = fmpz_mat_rank(m1);
  fmpz_mat_clear(m1);
  return rk;
}

template <>
void determinant<M2::ARingQQFlint, M2::ARingQQFlint::ElementType>(
    const DMatQQFlint& A,
    M2::ARingQQFlint::ElementType& result_det)
{
  fmpq_mat_det(&result_det, A.fmpq_mat());
}

template <>
bool inverse(const DMatQQFlint& A, DMatQQFlint& result_inv)
{
  return fmpq_mat_inv(result_inv.fmpq_mat(), A.fmpq_mat());
}

template <>
size_t rowReducedEchelonForm(const DMatQQFlint& A,
                                    DMatQQFlint& result_rref)
{
  return fmpq_mat_rref(result_rref.fmpq_mat(), A.fmpq_mat());
}

template <>
size_t nullSpace(const DMatQQFlint& A, DMatQQFlint& result_nullspace)
{
  fmpz_mat_t m1;
  fmpz_mat_t m2;
  fmpz_mat_init(m1, A.numRows(), A.numColumns());
  fmpz_mat_init(m2, A.numColumns(), A.numColumns());
  fmpq_mat_get_fmpz_mat_rowwise(m1, nullptr, A.fmpq_mat());
  // fmpz_mat_print_pretty(m1);
  size_t nullity = fmpz_mat_nullspace(m2, m1);
  // now copy the first 'nullity' columns into result_nullspace
  result_nullspace.resize(A.numColumns(), nullity);
  for (size_t c = 0; c < nullity; c++)
    for (size_t r = 0; r < A.numColumns(); r++)
      fmpz_set(fmpq_numref(&result_nullspace.entry(r, c)),
               fmpz_mat_entry(m2, r, c));
  fmpz_mat_clear(m1);
  fmpz_mat_clear(m2);
  return nullity;
}

template <>
bool solveLinear(const DMatQQFlint& A,
                        const DMatQQFlint& B,
                        DMatQQFlint& X)
{
  // TODO: WRITE ME
  // DMatQQFlint& A1 = const_cast<DMatQQFlint&>(A); // needed because
  // fmpq_mat_solve doesn't declare params const
  // DMatQQFlint& B1 = const_cast<DMatQQFlint&>(B);
  //    return fmpq_mat_solve(X.fmpq_mat(), B1.fmpq_mat(), A1.fmpq_mat());
  (void) A;
  (void) B;
  (void) X;
  return false;
}

template <>
M2_arrayintOrNull rankProfile(const DMatQQFlint& A, bool row_profile)
{
  // TODO: WRITE ME
  (void) A;
  (void) row_profile;
  throw exc::engine_error(
      "'rankProfile' not implemented for this kind of matrix over this ring");
}

template <>
void addMultipleTo(DMatQQFlint& C,
                          const DMatQQFlint& A,
                          const DMatQQFlint& B)
{
  DMatQQFlint D(C.ring(), A.numRows(), B.numColumns());
  fmpq_mat_mul(D.fmpq_mat(), A.fmpq_mat(), B.fmpq_mat());
  fmpq_mat_add(C.fmpq_mat(), C.fmpq_mat(), D.fmpq_mat());
}

template <>
void subtractMultipleTo(DMatQQFlint& C,
                               const DMatQQFlint& A,
                               const DMatQQFlint& B)
{
  DMatQQFlint D(C.ring(), A.numRows(), B.numColumns());
  fmpq_mat_mul(D.fmpq_mat(), A.fmpq_mat(), B.fmpq_mat());
  fmpq_mat_sub(C.fmpq_mat(), C.fmpq_mat(), D.fmpq_mat());
}

template <>
void mult(const DMatQQFlint& A,
                 const DMatQQFlint& B,
                 DMatQQFlint& result_product)
{
  // The A and B on the next line are switched because the memory layout
  // expected
  // is the transpose of what we have for DMat.
  fmpq_mat_mul(result_product.fmpq_mat(), A.fmpq_mat(), B.fmpq_mat());
}


}  // namespace MatrixOps
