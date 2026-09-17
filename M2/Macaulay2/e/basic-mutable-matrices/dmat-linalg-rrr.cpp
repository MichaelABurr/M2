// Copyright 2026  The Macaulay2 Authors

#include "exceptions.hpp"
#include "basic-mutable-matrices/dmat-linalg-rrr.hpp"
#include "basic-mutable-matrices/dmat-lu-inplace-rrr.hpp"
#include "basic-rings/aring-RRR.hpp"
#include "basic-rings/aring-CCC.hpp"
#include "basic-mutable-matrices/dmat-linalg.hpp"
#include "eigen.hpp"

using DMatRRR = DMat<M2::ARingRRR>;
using DMatCCC = DMat<M2::ARingCCC>;

namespace MatrixOps {
/////////
// RRR //
/////////

template <>
bool eigenvaluesHermitian(const DMatRRR& A, DMatRRR& eigenvals)
{
  return EigenM2::eigenvalues_hermitian(&A, &eigenvals);
}

template <>
bool eigenvalues(const DMatRRR& A, DMatCCC& eigenvals)
{
  return EigenM2::eigenvalues(&A, &eigenvals);
}

template <>
bool eigenvectorsHermitian(const DMatRRR& A,
  DMatRRR& eigenvals,
  DMatRRR& eigenvecs)
{
  return EigenM2::eigenvectors_hermitian(&A, &eigenvals, &eigenvecs);
}

template <>
bool eigenvectors(const DMatRRR& A,
  DMatCCC& eigenvals,
  DMatCCC& eigenvecs)
{
  return EigenM2::eigenvectors(&A, &eigenvals, &eigenvecs);
}

template <>
bool leastSquares(const DMatRRR& A,
  const DMatRRR& B,
  DMatRRR& X,
  bool assume_full_rank)
{
  (void) assume_full_rank;
  return EigenM2::least_squares(&A, &B, &X);
}

template <>
bool SVD(const DMatRRR& A,
  DMatRRR& Sigma,
  DMatRRR& U,
  DMatRRR& Vt,
  int strategy)
{
  if (strategy == 1) return EigenM2::SVD_divide_conquer(&A, &Sigma, &U, &Vt);
  return EigenM2::SVD(&A, &Sigma, &U, &Vt);
}

template <>
void clean(gmp_RR epsilon, DMatRRR& mat)
{
  for (size_t r = 0; r < mat.numRows(); ++r)
    for (size_t c = 0; c < mat.numColumns(); ++c)
      {
        mat.ring().zeroize_tiny(epsilon, mat.entry(r,c));
      }
}

template <>
void increase_norm(gmp_RRmutable norm, const DMatRRR& mat)
{
  for (size_t r = 0; r < mat.numRows(); ++r)
    for (size_t c = 0; c < mat.numColumns(); ++c)
      {
        mat.ring().increase_norm(norm, mat.entry(r,c));
      }
}


}  // namespace MatrixOps
