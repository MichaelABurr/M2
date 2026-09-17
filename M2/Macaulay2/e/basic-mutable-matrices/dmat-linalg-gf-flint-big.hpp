// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_GF_FLINT_BIG_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_GF_FLINT_BIG_HPP_

#include "basic-mutable-matrices/mat-linalg.hpp"

namespace M2 {
class ARingGFFlintBig;
}

namespace MatrixOps {
template <>
void addMultipleTo(DMat<M2::ARingGFFlintBig>& C,
                   const DMat<M2::ARingGFFlintBig>& A,
                   const DMat<M2::ARingGFFlintBig>& B);
template <>
void subtractMultipleTo(DMat<M2::ARingGFFlintBig>& C,
                        const DMat<M2::ARingGFFlintBig>& A,
                        const DMat<M2::ARingGFFlintBig>& B);
template <>
void mult(const DMat<M2::ARingGFFlintBig>& A,
          const DMat<M2::ARingGFFlintBig>& B,
          DMat<M2::ARingGFFlintBig>& result_product);
template <>
size_t rowReducedEchelonForm(const DMat<M2::ARingGFFlintBig>& A,
                             DMat<M2::ARingGFFlintBig>& result_rref);
}  // namespace MatrixOps

#endif
