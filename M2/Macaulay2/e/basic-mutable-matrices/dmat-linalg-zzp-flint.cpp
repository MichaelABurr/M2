// Copyright 2026  The Macaulay2 Authors

#include "exceptions.hpp"
#include "basic-mutable-matrices/dmat-linalg-zzp-flint.hpp"
#include "basic-mutable-matrices/dmat-lu-zzp-flint.hpp"
#include "basic-mutable-matrices/dmat-zzp-flint.hpp"
#include "basic-mutable-matrices/dmat-linalg.hpp"

using DMatZZpFlint = DMat<M2::ARingZZpFlint>;

namespace MatrixOps {
//////////////////////
// ZZpFlint //////////
//////////////////////

// Functions for DMatZZpFlint

template <>
void addMultipleTo(DMatZZpFlint& C,
                          const DMatZZpFlint& A,
                          const DMatZZpFlint& B)
{
  DMatZZpFlint D(C.ring(), A.numRows(), B.numColumns());
  nmod_mat_mul(D.nmod_mat(), A.nmod_mat(), B.nmod_mat());
  nmod_mat_add(C.nmod_mat(), C.nmod_mat(), D.nmod_mat());
}

template <>
void subtractMultipleTo(DMatZZpFlint& C,
                               const DMatZZpFlint& A,
                               const DMatZZpFlint& B)
{
  DMatZZpFlint D(C.ring(), A.numRows(), B.numColumns());
  nmod_mat_mul(D.nmod_mat(), A.nmod_mat(), B.nmod_mat());
  nmod_mat_sub(C.nmod_mat(), C.nmod_mat(), D.nmod_mat());
}

template <>
void mult(const DMatZZpFlint& A,
                 const DMatZZpFlint& B,
                 DMatZZpFlint& result_product)
{
  //    DMatZZpFlint& A1 = const_cast<DMatZZpFlint&>(A); // needed because
  //    nmod_mat_mul doesn't declare params const
  //    DMatZZpFlint& B1 = const_cast<DMatZZpFlint&>(B);
  // The A1 and B1 on the next line are switched because the memory layout
  // expected
  // is the transpose of what we have for DMat.
  nmod_mat_mul(result_product.nmod_mat(), A.nmod_mat(), B.nmod_mat());
}

template <>
size_t rowReducedEchelonForm(const DMatZZpFlint& A,
                                    DMatZZpFlint& result_rref)
{
  DMatZZpFlint A1(A);
  long rank = nmod_mat_rref(A1.nmod_mat());
  result_rref.swap(A1);
  return rank;
}


}  // namespace MatrixOps
