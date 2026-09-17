// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_QQ_FLINT_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_QQ_FLINT_HPP_

#include "basic-mutable-matrices/mat-linalg.hpp"

#include <flint/fmpq.h>

namespace M2 {
class ARingQQFlint;
}

namespace MatrixOps {
template <>
size_t rank(const DMat<M2::ARingQQFlint>& A);
template <>
void determinant<M2::ARingQQFlint, fmpq>(
    const DMat<M2::ARingQQFlint>& A,
    fmpq& result_det);
template <>
bool inverse(const DMat<M2::ARingQQFlint>& A,
             DMat<M2::ARingQQFlint>& result_inv);
template <>
size_t rowReducedEchelonForm(const DMat<M2::ARingQQFlint>& A,
                             DMat<M2::ARingQQFlint>& result_rref);
template <>
size_t nullSpace(const DMat<M2::ARingQQFlint>& A,
                 DMat<M2::ARingQQFlint>& result_nullspace);
template <>
bool solveLinear(const DMat<M2::ARingQQFlint>& A,
                 const DMat<M2::ARingQQFlint>& B,
                 DMat<M2::ARingQQFlint>& X);
template <>
M2_arrayintOrNull rankProfile(const DMat<M2::ARingQQFlint>& A,
                              bool row_profile);
template <>
void addMultipleTo(DMat<M2::ARingQQFlint>& C,
                   const DMat<M2::ARingQQFlint>& A,
                   const DMat<M2::ARingQQFlint>& B);
template <>
void subtractMultipleTo(DMat<M2::ARingQQFlint>& C,
                        const DMat<M2::ARingQQFlint>& A,
                        const DMat<M2::ARingQQFlint>& B);
template <>
void mult(const DMat<M2::ARingQQFlint>& A,
          const DMat<M2::ARingQQFlint>& B,
          DMat<M2::ARingQQFlint>& result_product);
}  // namespace MatrixOps

#endif
