// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_GF_FLINT_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_GF_FLINT_HPP_

#include "basic-mutable-matrices/mat-linalg.hpp"

namespace M2 {
class ARingGFFlint;
}

namespace MatrixOps {
template <>
void addMultipleTo(DMat<M2::ARingGFFlint>& C,
                   const DMat<M2::ARingGFFlint>& A,
                   const DMat<M2::ARingGFFlint>& B);
template <>
void subtractMultipleTo(DMat<M2::ARingGFFlint>& C,
                        const DMat<M2::ARingGFFlint>& A,
                        const DMat<M2::ARingGFFlint>& B);
template <>
void mult(const DMat<M2::ARingGFFlint>& A,
          const DMat<M2::ARingGFFlint>& B,
          DMat<M2::ARingGFFlint>& result_product);
template <>
size_t rowReducedEchelonForm(const DMat<M2::ARingGFFlint>& A,
                             DMat<M2::ARingGFFlint>& result_rref);
}  // namespace MatrixOps

#endif
