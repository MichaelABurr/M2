// Copyright 2026  The Macaulay2 Authors

#include "exceptions.hpp"
#include "basic-mutable-matrices/dmat-linalg-rr.hpp"
#include "basic-mutable-matrices/dmat-lu-inplace-rr.hpp"
#include "basic-rings/aring-CC.hpp"
#include "basic-rings/aring-RR.hpp"
#include "basic-mutable-matrices/lapack.hpp"
#include "basic-mutable-matrices/dmat-linalg.hpp"
#include "eigen.hpp"
namespace MatrixOps {
////////
// RR //
////////
template <>
bool eigenvaluesHermitian(const DMatRR& A, DMatRR& eigenvals)
{
#ifndef NO_LAPACK
  return Lapack::eigenvalues_symmetric(&A, &eigenvals);
#else
  return EigenM2::eigenvalues_hermitian(&A, &eigenvals);
#endif
}

template <>
bool eigenvalues(const DMatRR& A, DMatCC& eigenvals)
{
#ifndef NO_LAPACK
  return Lapack::eigenvalues(&A, &eigenvals);
#else
  return EigenM2::eigenvalues(&A, &eigenvals);
#endif
}

template <>
bool eigenvectorsHermitian(const DMatRR& A,
  DMatRR& eigenvals,
  DMatRR& eigenvecs)
{
#ifndef NO_LAPACK
  return Lapack::eigenvectors_symmetric(&A, &eigenvals, &eigenvecs);
#else
  return EigenM2::eigenvectors_hermitian(&A, &eigenvals, &eigenvecs);
#endif
}

template <>
bool eigenvectors(const DMatRR& A, DMatCC& eigenvals, DMatCC& eigenvecs)
{
#ifndef NO_LAPACK
  return Lapack::eigenvectors(&A, &eigenvals, &eigenvecs);
#else
  return EigenM2::eigenvectors(&A, &eigenvals, &eigenvecs);
#endif
}

template <>
bool leastSquares(const DMatRR& A,
  const DMatRR& B,
  DMatRR& X,
  bool assume_full_rank)
{
#ifndef NO_LAPACK
  if (assume_full_rank)
    return Lapack::least_squares(&A, &B, &X);
  else
    return Lapack::least_squares_deficient(&A, &B, &X);
#else
  // place eigen code here
  return EigenM2::least_squares(&A, &B, &X);
  // throw exc::engine_error( "not implemented!!!");
  // return false; // indicates error
#endif
}

template <>
bool SVD(const DMatRR& A,
  DMatRR& Sigma,
  DMatRR& U,
  DMatRR& Vt,
  int strategy)
{
#ifndef NO_LAPACK
  if (strategy == 1) return Lapack::SVD_divide_conquer(&A, &Sigma, &U, &Vt);
  return Lapack::SVD(&A, &Sigma, &U, &Vt);
#else
  if (strategy == 1) return EigenM2::SVD_divide_conquer(&A, &Sigma, &U, &Vt);
  return EigenM2::SVD(&A, &Sigma, &U, &Vt);
#endif
}

template <>
bool QR(const DMatRR& A, DMatRR& Q, DMatRR& R, bool return_QR)
{
  return Lapack::QR(&A, &Q, &R, return_QR);
}

template <>
void clean(gmp_RR epsilon, DMatRR& mat)
{
  for (size_t r = 0; r < mat.numRows(); ++r)
    for (size_t c = 0; c < mat.numColumns(); ++c)
      {
        mat.ring().zeroize_tiny(epsilon, mat.entry(r,c));
      }
}

template <>
void increase_norm(gmp_RRmutable norm, const DMatRR& mat)
{
  for (size_t r = 0; r < mat.numRows(); ++r)
    for (size_t c = 0; c < mat.numColumns(); ++c)
      {
        mat.ring().increase_norm(norm, mat.entry(r,c));
      }
}


}  // namespace MatrixOps
