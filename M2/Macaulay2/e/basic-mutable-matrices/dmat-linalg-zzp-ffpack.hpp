// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_ZZP_FFPACK_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_ZZP_FFPACK_HPP_

#include "basic-mutable-matrices/mat-linalg.hpp"

namespace M2 {
class ARingZZpFFPACK;
}

namespace MatrixOps {
template <>
void mult(const DMat<M2::ARingZZpFFPACK>& A,
          const DMat<M2::ARingZZpFFPACK>& B,
          DMat<M2::ARingZZpFFPACK>& result_product);
template <>
void addMultipleTo(DMat<M2::ARingZZpFFPACK>& C,
                   const DMat<M2::ARingZZpFFPACK>& A,
                   const DMat<M2::ARingZZpFFPACK>& B);
template <>
void subtractMultipleTo(DMat<M2::ARingZZpFFPACK>& C,
                        const DMat<M2::ARingZZpFFPACK>& A,
                        const DMat<M2::ARingZZpFFPACK>& B);
}  // namespace MatrixOps

#endif
