// Copyright 2013  Michael E. Stillman

#ifndef M2_BASIC_MUTMAT_MAT_LINALG_HPP_
#define M2_BASIC_MUTMAT_MAT_LINALG_HPP_

/**
 * \ingroup matrices
 */

#include "exceptions.hpp"
#include "interface/m2-types.h"

#include <vector>

template <typename RT>
class DMat;

namespace MatrixOps {
// DMat primary overloads. Their definitions are in dmat-linalg.hpp.
template <typename RT>
size_t rank(const DMat<RT>& A);

template <typename RT, typename Element>
void determinant(const DMat<RT>& A, Element& result_det);

template <typename RT>
bool inverse(const DMat<RT>& A, DMat<RT>& result_inv);

template <typename RT>
void mult(const DMat<RT>& A,
          const DMat<RT>& B,
          DMat<RT>& result_product);

template <typename RT>
size_t nullSpace(const DMat<RT>& A, DMat<RT>& result_nullspace);

template <typename RT>
bool solveLinear(const DMat<RT>& A, const DMat<RT>& B, DMat<RT>& X);

template <typename RT>
bool solveInvertible(const DMat<RT>& A,
                     const DMat<RT>& B,
                     DMat<RT>& X);

template <typename RT>
M2_arrayintOrNull rankProfile(const DMat<RT>& A, bool row_profile);

template <typename RT>
void addMultipleTo(DMat<RT>& C, const DMat<RT>& A, const DMat<RT>& B);

template <typename RT>
void subtractMultipleTo(DMat<RT>& C,
                        const DMat<RT>& A,
                        const DMat<RT>& B);

template <typename RT>
M2_arrayintOrNull LU(const DMat<RT>& A, DMat<RT>& L, DMat<RT>& U);

template <typename RT>
M2_arrayintOrNull LUincremental(std::vector<size_t>& P,
                                DMat<RT>& LU,
                                const DMat<RT>& v,
                                int i);

template <typename RT>
void triangularSolve(DMat<RT>& Lv, DMat<RT>& x, int m, int strategy);

/// @brief the rank of a matrix
///
/// throws an engine_error for ring/matrix types where the function is not
/// implemented.
/// This version is deterministic.
template <typename Mat>
size_t rank(const Mat& A)
{
  (void) A;
  throw exc::engine_error(
      "'rank' not implemented for this kind of matrix over this ring");
  return 0;
}

/// @brief the determinant of a square matrix
///
/// result_det should be a previously initialized ElementType.
/// throws an engine_error for ring/matrix types where the function is not
/// implemented.
template <typename Mat>
void determinant(const Mat& A, typename Mat::ElementType& result_det)
{
  (void) A;
  (void) result_det;
  throw exc::engine_error(
      "'determinant' not implemented for this kind of matrix over this ring");
}

/// @brief the inverse of a square matrix
///
/// result_inv is set to the inverse of the square matrix A, if A is invertible.
/// result_inv should be a Mat, with the same ring/type as the input matrix A.
///   result_inv does not need to be the same size as A, it will be resized if
///   needed.
/// returns true exactly when the matrix is invertible, and result_inv has been
/// set.
///
/// throws an engine_error for ring/matrix types where the function is not
/// implemented.
/// throws an error if the matrix is not square.
///
/// Note: the inverse of a 0 x 0 matrix is another 0 x 0 matrix.
template <typename Mat>
bool inverse(const Mat& A, Mat& result_inv)
{
  (void) A;
  (void) result_inv;
  throw exc::engine_error(
      "'invert' not implemented for this kind of matrix over this ring");
}

/// @brief the row reduced echelon form of a matrix over a field, or ZZ.
///
/// result_rref should be a Mat, with the same ring/type as the input matrix A.
///   result_rref does not need to be the same size as A, it will be resized if
///   needed.
/// returns the rank of A.
///
/// throws an engine_error for ring/matrix types where the function is not
/// implemented.
template <typename Mat>
size_t rowReducedEchelonForm(const Mat& A, Mat& result_rref)
{
  (void) A;
  (void) result_rref;
  throw exc::engine_error(
      "'rowReducedEchelonForm' not implemented for this kind of matrix over "
      "this ring");
}

/// @brief the product of two matrices
///
/// result_product is set to the product A*B
/// result_product should be a Mat, with the same ring/type as the input
/// matrices A,B.
///   result_product does not need to be the same size as A*B, it will be
///   resized if needed.
///
/// throws an engine_error for ring/matrix types where the function is not
/// implemented.
/// throws an error if the number of columns of A is not the number of rows of
/// B.
/// result_prod should not be the same as A or B (assertion error).
template <typename Mat>
void mult(const Mat& A, const Mat& B, Mat& result_product)
{
  (void) A;
  (void) B;
  (void) result_product;
  throw exc::engine_error(
      "'mult matrices' not implemented for this kind of matrix over this ring");
}

/// @brief the null space of a matrix
///
///   result_nullspace is set to the matrix whose columns form a basis for {x |
///   Ax = 0}.
/// Returns the dimension of the nullspace.
///
/// result_nullspace should be a Mat, with the same ring/type as the input
/// matrix A.
///   result_nullspace does not need to be the correct size, it will be resized
///   if needed.
///
/// throws an engine_error for ring/matrix types where the function is not
/// implemented.
template <typename Mat>
size_t nullSpace(const Mat& A, Mat& result_nullspace)
{
  (void) A;
  (void) result_nullspace;
  throw exc::engine_error(
      "'nullSpace' not implemented for this kind of matrix over this ring");
}

/// @brief solve AX=B, return true if the system has a solution.
template <typename Mat>
bool solveLinear(const Mat& A, const Mat& B, Mat& X)
{
  (void) A;
  (void) B;
  (void) X;
  throw exc::engine_error(
      "'solveLinear' not implemented for this kind of matrix over this ring");
}

/// @brief solve AX=B, where A is a square (invertible) matrix.
///
/// return true if A is invertible, and in this case, sets X.  If false is
/// returned, X's value is
/// undefined.
template <typename Mat>
bool solveInvertible(const Mat& A, const Mat& B, Mat& X)
{
  (void) A;
  (void) B;
  (void) X;
  throw exc::engine_error(
      "'solveInvertible' not implemented for this kind of matrix over this "
      "ring");
}

/// @brief Returns either the row or column rank profile of A
///
/// if row_profile is true, then row profile is computed, otherwise
/// the column profile is computed.
///
/// The return value is an ascending sequence of non-negative integers
/// with an entry a occurring iff the submatrix of A of the first
/// (a-1) rows (resp columns) has lower rank than the submatrix of the
/// first a rows (resp columns).
///
/// Notice that if the matrix is non-zero and the first row is
/// non-zero, then the first entry will be 0.
template <typename Mat>
M2_arrayintOrNull rankProfile(const Mat& A, bool row_profile)
{
  (void) A;
  (void) row_profile;
  throw exc::engine_error(
      "'rankProfile' not implemented for this kind of matrix over this ring");
}

/// @brief Set C += A*B
///
/// Throws an exception if not yet implemented for this ring/matrix type.
/// The sizes of C,A,B must be compatible.  These are checked only via
/// assertions.
template <typename Mat>
void addMultipleTo(Mat& C, const Mat& A, const Mat& B)
// C = C + A*B
{
  (void) C;
  (void) A;
  (void) B;
  throw exc::engine_error(
      "'addMultipleTo' not implemented for this kind of matrix over this ring");
}

/// @brief Set C -= A*B
///
/// Throws an exception if not yet implemented for this ring/matrix type.
/// The sizes of C,A,B must be compatible.  These are checked only via
/// assertions.
template <typename Mat>
void subtractMultipleTo(Mat& C, const Mat& A, const Mat& B)
// C = C - A*B
{
  (void) C;
  (void) A;
  (void) B;
  throw exc::engine_error(
      "'subtractMultipleTo' not implemented for this kind of matrix over this "
      "ring");
}

template <typename Mat>
M2_arrayintOrNull LU(const Mat& A, Mat& L, Mat& U)
{
  (void) A;
  (void) L;
  (void) U;
  throw exc::engine_error(
      "'LU' not implemented for this kind of matrix over this ring");
}

template <typename Mat>
M2_arrayintOrNull LUincremental(std::vector<size_t>& P, Mat& LU, const Mat& v, int i)
{
  (void) P;
  (void) LU;
  (void) v;
  (void) i;
  throw exc::engine_error(
      "'LUincremental' not implemented for this kind of matrix over this ring");
}

template <typename Mat>
void triangularSolve(Mat& Lv, Mat& x, int m, int strategy)
{
  (void) Lv;
  (void) x;
  (void) m;
  (void) strategy;
  throw exc::engine_error(
      "'triangularSolve' not implemented for this kind of matrix over this "
      "ring");
}

template <typename Mat, typename Mat2>
bool eigenvalues(const Mat& A, Mat2& eigenvals)
{
  (void) A;
  (void) eigenvals;
  throw exc::engine_error(
      "'eigenvalues' not implemented for this kind of matrix over this ring");
}

template <typename Mat, typename Mat2>
bool eigenvaluesHermitian(const Mat& A, Mat2& eigenvals)
{
  (void) A;
  (void) eigenvals;
  throw exc::engine_error(
      "'eigenvalues' not implemented for this kind of matrix over this ring");
}

template <typename Mat, typename Mat2, typename Mat3>
bool eigenvectors(const Mat& A, Mat2& eigenvals, Mat3& eigenvecs)
{
  (void) A;
  (void) eigenvals;
  (void) eigenvecs;
  throw exc::engine_error(
      "'eigenvectors' not implemented for this kind of matrix over this ring");
}

template <typename Mat, typename Mat2, typename Mat3>
bool eigenvectorsHermitian(const Mat& A, Mat2& eigenvals, Mat3& eigenvecs)
{
  (void) A;
  (void) eigenvals;
  (void) eigenvecs;
  throw exc::engine_error(
      "'eigenvectors' not implemented for this kind of matrix over this ring");
}

template <typename Mat>
bool leastSquares(const Mat& A, const Mat& B, Mat& X, bool assume_full_rank)
{
  (void) A;
  (void) B;
  (void) X;
  (void) assume_full_rank;
  throw exc::engine_error(
      "'leastSquares' not implemented for this kind of matrix over this ring");
}

template <typename Mat, typename Mat2>
bool SVD(const Mat& A, Mat2& Sigma, Mat& U, Mat& Vt, int strategy)
{
  (void) A;
  (void) Sigma;
  (void) U;
  (void) Vt;
  (void) strategy;
  throw exc::engine_error(
      "'SVD' not implemented for this kind of matrix over this ring");
}

template <typename Mat, typename Mat2, typename Mat3>
bool QR(const Mat& A, Mat2& Q, Mat3& R, bool return_QR)
{
  (void) A;
  (void) Q;
  (void) R;
  (void) return_QR;
  throw exc::engine_error(
      "'QR' not implemented for this kind of matrix over this ring");
}

template <typename T>
void clean(gmp_RR epsilon, T& mat)
{
  (void) epsilon;
  (void) mat;
  throw exc::engine_error(
      "'clean' not implemented for this kind of matrix over this ring");
}

template <typename T>
void increase_norm(gmp_RRmutable nm, const T& mat)
{
  (void) nm;
  (void) mat;
  throw exc::engine_error(
      "'norm' not implemented for this kind of matrix over this ring");
}

}  // namespace MatrixOps

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
