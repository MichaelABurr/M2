// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_ZZ_GMP_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_ZZ_GMP_HPP_

#include "basic-mutable-matrices/mat-linalg.hpp"

namespace M2 {
class ARingZZGMP;
}

namespace MatrixOps {
template <>
M2_arrayintOrNull LU(const DMat<M2::ARingZZGMP>& A,
                     DMat<M2::ARingZZGMP>& L,
                     DMat<M2::ARingZZGMP>& U);
template <>
M2_arrayintOrNull rankProfile(const DMat<M2::ARingZZGMP>& A,
                              bool row_profile);
template <>
bool inverse(const DMat<M2::ARingZZGMP>& A,
             DMat<M2::ARingZZGMP>& result_inv);
template <>
size_t nullSpace(const DMat<M2::ARingZZGMP>& A,
                 DMat<M2::ARingZZGMP>& result_nullspace);
template <>
bool solveLinear(const DMat<M2::ARingZZGMP>& A,
                 const DMat<M2::ARingZZGMP>& B,
                 DMat<M2::ARingZZGMP>& X);
template <>
bool solveInvertible(const DMat<M2::ARingZZGMP>& A,
                     const DMat<M2::ARingZZGMP>& B,
                     DMat<M2::ARingZZGMP>& X);
template <>
void mult(const DMat<M2::ARingZZGMP>& A,
          const DMat<M2::ARingZZGMP>& B,
          DMat<M2::ARingZZGMP>& result_product);
template <>
void addMultipleTo(DMat<M2::ARingZZGMP>& C,
                   const DMat<M2::ARingZZGMP>& A,
                   const DMat<M2::ARingZZGMP>& B);
template <>
void subtractMultipleTo(DMat<M2::ARingZZGMP>& C,
                        const DMat<M2::ARingZZGMP>& A,
                        const DMat<M2::ARingZZGMP>& B);
template <>
size_t rank(const DMat<M2::ARingZZGMP>& A);
template <>
void determinant<M2::ARingZZGMP, __mpz_struct>(
    const DMat<M2::ARingZZGMP>& A,
    __mpz_struct& result_det);
}  // namespace MatrixOps

#endif
