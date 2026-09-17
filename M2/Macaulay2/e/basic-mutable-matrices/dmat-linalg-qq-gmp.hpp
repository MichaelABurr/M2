// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_QQ_GMP_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_QQ_GMP_HPP_

#include "basic-mutable-matrices/mat-linalg.hpp"

namespace M2 {
class ARingQQGMP;
}

namespace MatrixOps {
template <>
void mult(const DMat<M2::ARingQQGMP>& A,
          const DMat<M2::ARingQQGMP>& B,
          DMat<M2::ARingQQGMP>& result_product);
template <>
void addMultipleTo(DMat<M2::ARingQQGMP>& C,
                   const DMat<M2::ARingQQGMP>& A,
                   const DMat<M2::ARingQQGMP>& B);
template <>
void subtractMultipleTo(DMat<M2::ARingQQGMP>& C,
                        const DMat<M2::ARingQQGMP>& A,
                        const DMat<M2::ARingQQGMP>& B);
template <>
size_t rowReducedEchelonForm(const DMat<M2::ARingQQGMP>& A,
                             DMat<M2::ARingQQGMP>& result_rref);
}  // namespace MatrixOps

#endif
