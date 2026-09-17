// Copyright 2026  The Macaulay2 Authors

#include "exceptions.hpp"
#include "basic-mutable-matrices/dmat-linalg-gf-flint.hpp"
#include "basic-mutable-matrices/dmat-lu-inplace-gf-flint.hpp"
#include "basic-mutable-matrices/dmat-gf-flint.hpp"
#include "basic-rings/aring-GF-flint.hpp"
#include "basic-mutable-matrices/dmat-linalg.hpp"

using DMatGFFlint = DMat<M2::ARingGFFlint>;

namespace MatrixOps {
////////////////////////
// GFFlint /////////////
////////////////////////

// Functions for DMatGFFlint, linear algebra is sent out to LU

template <>
void addMultipleTo(DMatGFFlint& C,
                          const DMatGFFlint& A,
                          const DMatGFFlint& B)
{
  DMatGFFlint D(C.ring(), A.numRows(), B.numColumns());
  fq_zech_mat_mul(D.fq_zech_mat(),
                  A.fq_zech_mat(),
                  B.fq_zech_mat(),
                  A.ring().flintContext());
  fq_zech_mat_add(C.fq_zech_mat(),
                  C.fq_zech_mat(),
                  D.fq_zech_mat(),
                  A.ring().flintContext());
}

template <>
void subtractMultipleTo(DMatGFFlint& C,
                               const DMatGFFlint& A,
                               const DMatGFFlint& B)
{
  DMatGFFlint D(C.ring(), A.numRows(), B.numColumns());
  fq_zech_mat_mul(D.fq_zech_mat(),
                  A.fq_zech_mat(),
                  B.fq_zech_mat(),
                  A.ring().flintContext());
  fq_zech_mat_sub(C.fq_zech_mat(),
                  C.fq_zech_mat(),
                  D.fq_zech_mat(),
                  A.ring().flintContext());
}

template <>
void mult(const DMatGFFlint& A,
                 const DMatGFFlint& B,
                 DMatGFFlint& result_product)
{
  //    DMatGFFlint& A1 = const_cast<DMatGFFlint&>(A); // needed because
  //    nmod_mat_mul doesn't declare params const
  //    DMatGFFlint& B1 = const_cast<DMatGFFlint&>(B);
  // The A1 and B1 on the next line are switched because the memory layout
  // expected
  // is the transpose of what we have for DMat.
  fq_zech_mat_mul(result_product.fq_zech_mat(),
                  A.fq_zech_mat(),
                  B.fq_zech_mat(),
                  A.ring().flintContext());
}

template <>
size_t rowReducedEchelonForm(const DMatGFFlint& A,
                                    DMatGFFlint& result_rref)
{
  DMatGFFlint A1(A);
#if __FLINT_RELEASE >= 30100
  long rank = fq_zech_mat_rref(A1.fq_zech_mat(), A1.fq_zech_mat(), A.ring().flintContext());
#else
  long rank = fq_zech_mat_rref(A1.fq_zech_mat(), A.ring().flintContext());
#endif
  result_rref.swap(A1);
  return rank;
}


}  // namespace MatrixOps
