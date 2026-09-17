// Copyright 2026  The Macaulay2 Authors

#include "exceptions.hpp"
#include "basic-mutable-matrices/dmat-linalg-gf-flint-big.hpp"
#include "basic-mutable-matrices/dmat-lu-inplace-gf-flint-big.hpp"
#include "basic-mutable-matrices/dmat-gf-flint-big.hpp"
#include "basic-rings/aring-GF-flint-big.hpp"
#include "basic-mutable-matrices/dmat-linalg.hpp"

using DMatGFFlintBig = DMat<M2::ARingGFFlintBig>;

namespace MatrixOps {
////////////////////////
// GFFlintBig //////////
////////////////////////

// Functions for DMatGFFlintBig, linear algebra is sent out to LU

template <>
void addMultipleTo(DMatGFFlintBig& C,
                          const DMatGFFlintBig& A,
                          const DMatGFFlintBig& B)
{
  DMatGFFlintBig D(C.ring(), A.numRows(), B.numColumns());
  fq_nmod_mat_mul(D.fq_nmod_mat(),
                  A.fq_nmod_mat(),
                  B.fq_nmod_mat(),
                  A.ring().flintContext());
  fq_nmod_mat_add(C.fq_nmod_mat(),
                  C.fq_nmod_mat(),
                  D.fq_nmod_mat(),
                  A.ring().flintContext());
}

template <>
void subtractMultipleTo(DMatGFFlintBig& C,
                               const DMatGFFlintBig& A,
                               const DMatGFFlintBig& B)
{
  DMatGFFlintBig D(C.ring(), A.numRows(), B.numColumns());
  fq_nmod_mat_mul(D.fq_nmod_mat(),
                  A.fq_nmod_mat(),
                  B.fq_nmod_mat(),
                  A.ring().flintContext());
  fq_nmod_mat_sub(C.fq_nmod_mat(),
                  C.fq_nmod_mat(),
                  D.fq_nmod_mat(),
                  A.ring().flintContext());
}

template <>
void mult(const DMatGFFlintBig& A,
                 const DMatGFFlintBig& B,
                 DMatGFFlintBig& result_product)
{
  //    DMatGFFlintBig& A1 = const_cast<DMatGFFlintBig&>(A); // needed because
  //    nmod_mat_mul doesn't declare params const
  //    DMatGFFlintBig& B1 = const_cast<DMatGFFlintBig&>(B);
  // The A1 and B1 on the next line are switched because the memory layout
  // expected
  // is the transpose of what we have for DMat.
  fq_nmod_mat_mul(result_product.fq_nmod_mat(),
                  A.fq_nmod_mat(),
                  B.fq_nmod_mat(),
                  A.ring().flintContext());
}

template <>
size_t rowReducedEchelonForm(const DMatGFFlintBig& A,
                                    DMatGFFlintBig& result_rref)
{
  DMatGFFlintBig A1(A);
#if __FLINT_RELEASE >= 30100
  long rank = fq_nmod_mat_rref(A1.fq_nmod_mat(), A1.fq_nmod_mat(), A.ring().flintContext());
#else
  long rank = fq_nmod_mat_rref(A1.fq_nmod_mat(), A.ring().flintContext());
#endif
  result_rref.swap(A1);
  return rank;
}


}  // namespace MatrixOps
