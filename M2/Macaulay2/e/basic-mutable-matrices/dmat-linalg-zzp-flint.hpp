// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_ZZP_FLINT_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_ZZP_FLINT_HPP_

#include "basic-mutable-matrices/mat-linalg.hpp"

namespace M2 {
class ARingZZpFlint;
}

namespace MatrixOps {
template <>
void addMultipleTo(DMat<M2::ARingZZpFlint>& C,
                   const DMat<M2::ARingZZpFlint>& A,
                   const DMat<M2::ARingZZpFlint>& B);
template <>
void subtractMultipleTo(DMat<M2::ARingZZpFlint>& C,
                        const DMat<M2::ARingZZpFlint>& A,
                        const DMat<M2::ARingZZpFlint>& B);
template <>
void mult(const DMat<M2::ARingZZpFlint>& A,
          const DMat<M2::ARingZZpFlint>& B,
          DMat<M2::ARingZZpFlint>& result_product);
template <>
size_t rowReducedEchelonForm(const DMat<M2::ARingZZpFlint>& A,
                             DMat<M2::ARingZZpFlint>& result_rref);
}  // namespace MatrixOps

#endif
