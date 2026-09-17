// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_CC_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_CC_HPP_

#include "basic-mutable-matrices/mat-linalg.hpp"

namespace M2 {
class ARingCC;
class ARingRR;
}

namespace MatrixOps {
template <>
bool QR(const DMat<M2::ARingCC>& A,
        DMat<M2::ARingCC>& Q,
        DMat<M2::ARingCC>& R,
        bool return_QR);
template <>
bool eigenvaluesHermitian(const DMat<M2::ARingCC>& A,
                          DMat<M2::ARingRR>& eigenvals);
template <>
bool eigenvalues(const DMat<M2::ARingCC>& A,
                 DMat<M2::ARingCC>& eigenvals);
template <>
bool eigenvectorsHermitian(const DMat<M2::ARingCC>& A,
                           DMat<M2::ARingRR>& eigenvals,
                           DMat<M2::ARingCC>& eigenvecs);
template <>
bool eigenvectors(const DMat<M2::ARingCC>& A,
                  DMat<M2::ARingCC>& eigenvals,
                  DMat<M2::ARingCC>& eigenvecs);
template <>
bool leastSquares(const DMat<M2::ARingCC>& A,
                  const DMat<M2::ARingCC>& B,
                  DMat<M2::ARingCC>& X,
                  bool assume_full_rank);
template <>
bool SVD(const DMat<M2::ARingCC>& A,
         DMat<M2::ARingRR>& Sigma,
         DMat<M2::ARingCC>& U,
         DMat<M2::ARingCC>& Vt,
         int strategy);
template <>
void clean(gmp_RR epsilon, DMat<M2::ARingCC>& mat);
template <>
void increase_norm(gmp_RRmutable norm, const DMat<M2::ARingCC>& mat);
}  // namespace MatrixOps

#endif
