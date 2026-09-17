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
template <>
size_t rank(const DMat<M2::ARingZZpFFPACK>& A);
template <>
void determinant<M2::ARingZZpFFPACK, double>(
    const DMat<M2::ARingZZpFFPACK>& A,
    double& result_det);
template <>
bool inverse(const DMat<M2::ARingZZpFFPACK>& A,
             DMat<M2::ARingZZpFFPACK>& result_inv);
template <>
size_t nullSpace(const DMat<M2::ARingZZpFFPACK>& A,
                 DMat<M2::ARingZZpFFPACK>& result_nullspace);
template <>
bool solveLinear(const DMat<M2::ARingZZpFFPACK>& A,
                 const DMat<M2::ARingZZpFFPACK>& B,
                 DMat<M2::ARingZZpFFPACK>& X);
template <>
bool solveInvertible(const DMat<M2::ARingZZpFFPACK>& A,
                     const DMat<M2::ARingZZpFFPACK>& B,
                     DMat<M2::ARingZZpFFPACK>& X);
template <>
M2_arrayintOrNull rankProfile(const DMat<M2::ARingZZpFFPACK>& A,
                              bool row_profile);
template <>
M2_arrayintOrNull LU(const DMat<M2::ARingZZpFFPACK>& A,
                     DMat<M2::ARingZZpFFPACK>& L,
                     DMat<M2::ARingZZpFFPACK>& U);
}  // namespace MatrixOps

#endif
