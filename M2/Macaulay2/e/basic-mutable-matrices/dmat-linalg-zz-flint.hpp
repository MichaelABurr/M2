// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_ZZ_FLINT_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_ZZ_FLINT_HPP_

#include "basic-mutable-matrices/mat-linalg.hpp"

#include <flint/fmpz.h>

namespace M2 {
class ARingZZ;
}

namespace MatrixOps {
template <>
size_t rank(const DMat<M2::ARingZZ>& A);
template <>
void determinant<M2::ARingZZ, fmpz>(const DMat<M2::ARingZZ>& A,
                                    fmpz& result_det);
template <>
bool inverse(const DMat<M2::ARingZZ>& A, DMat<M2::ARingZZ>& result_inv);
template <>
void mult(const DMat<M2::ARingZZ>& A,
          const DMat<M2::ARingZZ>& B,
          DMat<M2::ARingZZ>& result_product);
template <>
size_t nullSpace(const DMat<M2::ARingZZ>& A,
                 DMat<M2::ARingZZ>& result_nullspace);
template <>
bool solveLinear(const DMat<M2::ARingZZ>& A,
                 const DMat<M2::ARingZZ>& B,
                 DMat<M2::ARingZZ>& X);
template <>
M2_arrayintOrNull rankProfile(const DMat<M2::ARingZZ>& A, bool row_profile);
template <>
void addMultipleTo(DMat<M2::ARingZZ>& C,
                   const DMat<M2::ARingZZ>& A,
                   const DMat<M2::ARingZZ>& B);
template <>
void subtractMultipleTo(DMat<M2::ARingZZ>& C,
                        const DMat<M2::ARingZZ>& A,
                        const DMat<M2::ARingZZ>& B);
}  // namespace MatrixOps

#endif
