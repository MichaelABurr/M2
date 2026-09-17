// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_RR_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_RR_HPP_

#include "basic-mutable-matrices/mat-linalg.hpp"

namespace M2 {
class ARingCC;
class ARingRR;
}

namespace MatrixOps {
template <>
bool eigenvaluesHermitian(const DMat<M2::ARingRR>& A,
                          DMat<M2::ARingRR>& eigenvals);
template <>
bool eigenvalues(const DMat<M2::ARingRR>& A,
                 DMat<M2::ARingCC>& eigenvals);
template <>
bool eigenvectorsHermitian(const DMat<M2::ARingRR>& A,
                           DMat<M2::ARingRR>& eigenvals,
                           DMat<M2::ARingRR>& eigenvecs);
template <>
bool eigenvectors(const DMat<M2::ARingRR>& A,
                  DMat<M2::ARingCC>& eigenvals,
                  DMat<M2::ARingCC>& eigenvecs);
template <>
bool leastSquares(const DMat<M2::ARingRR>& A,
                  const DMat<M2::ARingRR>& B,
                  DMat<M2::ARingRR>& X,
                  bool assume_full_rank);
template <>
bool SVD(const DMat<M2::ARingRR>& A,
         DMat<M2::ARingRR>& Sigma,
         DMat<M2::ARingRR>& U,
         DMat<M2::ARingRR>& Vt,
         int strategy);
template <>
bool QR(const DMat<M2::ARingRR>& A,
        DMat<M2::ARingRR>& Q,
        DMat<M2::ARingRR>& R,
        bool return_QR);
template <>
void clean(gmp_RR epsilon, DMat<M2::ARingRR>& mat);
template <>
void increase_norm(gmp_RRmutable norm, const DMat<M2::ARingRR>& mat);
}  // namespace MatrixOps

#endif
