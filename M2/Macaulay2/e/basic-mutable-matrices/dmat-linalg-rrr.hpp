// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_RRR_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_RRR_HPP_

#include "basic-mutable-matrices/mat-linalg.hpp"

namespace M2 {
class ARingCCC;
class ARingRRR;
}

namespace MatrixOps {
template <>
bool eigenvaluesHermitian(const DMat<M2::ARingRRR>& A,
                          DMat<M2::ARingRRR>& eigenvals);
template <>
bool eigenvalues(const DMat<M2::ARingRRR>& A,
                 DMat<M2::ARingCCC>& eigenvals);
template <>
bool eigenvectorsHermitian(const DMat<M2::ARingRRR>& A,
                           DMat<M2::ARingRRR>& eigenvals,
                           DMat<M2::ARingRRR>& eigenvecs);
template <>
bool eigenvectors(const DMat<M2::ARingRRR>& A,
                  DMat<M2::ARingCCC>& eigenvals,
                  DMat<M2::ARingCCC>& eigenvecs);
template <>
bool leastSquares(const DMat<M2::ARingRRR>& A,
                  const DMat<M2::ARingRRR>& B,
                  DMat<M2::ARingRRR>& X,
                  bool assume_full_rank);
template <>
bool SVD(const DMat<M2::ARingRRR>& A,
         DMat<M2::ARingRRR>& Sigma,
         DMat<M2::ARingRRR>& U,
         DMat<M2::ARingRRR>& Vt,
         int strategy);
template <>
void clean(gmp_RR epsilon, DMat<M2::ARingRRR>& mat);
template <>
void increase_norm(gmp_RRmutable norm, const DMat<M2::ARingRRR>& mat);
}  // namespace MatrixOps

#endif
