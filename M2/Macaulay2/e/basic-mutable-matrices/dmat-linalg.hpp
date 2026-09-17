// Copyright 2013  Michael E. Stillman

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_HPP_

#include "util.hpp"
#include "basic-mutable-matrices/dmat.hpp"
#include "basic-mutable-matrices/dmat-lu.hpp"
#include "basic-mutable-matrices/mat-arith.hpp"
#include "basic-mutable-matrices/mat-elem-ops.hpp"
#include "basic-mutable-matrices/mat-linalg.hpp"

#include <algorithm>
#include <cassert>
#include <vector>

namespace MatrixOps {

template <typename RT>
void mult(const DMat<RT>& A, const DMat<RT>& B, DMat<RT>& result_product)
{
  assert(A.numColumns() == B.numRows());
  assert(A.numRows() == result_product.numRows());
  assert(B.numColumns() == result_product.numColumns());

  typename RT::Element tmp(A.ring());
  for (size_t i = 0; i < A.numRows(); i++)
    for (size_t j = 0; j < B.numColumns(); j++)
      {
        auto& val = result_product.entry(i,j);
        for (size_t k = 0; k < A.numColumns(); ++k)
          {
            A.ring().mult(tmp, A.entry(i,k), B.entry(k,j));
            A.ring().add(val, val, tmp);
          }
      }
}

template <typename RT>
void addMultipleTo(DMat<RT>& C, const DMat<RT>& A, const DMat<RT>& B)
{
  mult(A, B, C);
}

template <typename RT>
void subtractMultipleTo(DMat<RT>& C, const DMat<RT>& A, const DMat<RT>& B)
{
  assert(A.numColumns() == B.numRows());
  assert(A.numRows() == C.numRows());
  assert(B.numColumns() == C.numColumns());

  typename RT::Element tmp(A.ring());
  for (size_t i = 0; i < A.numRows(); i++)
    for (size_t j = 0; j < B.numColumns(); j++)
      {
        auto& val = C.entry(i,j);
        for (size_t k = 0; k < A.numColumns(); ++k)
          {
            A.ring().mult(tmp, A.entry(i,k), B.entry(k,j));
            A.ring().subtract(val, val, tmp);
          }
      }
}

// This default version only works for fields. Any other rings must provide
// specializations for these functions.
template <typename RT, typename Element>
inline void determinant(const DMat<RT>& A, Element& result)
{
  DMatLinAlg<RT> LUdecomp(A);
  LUdecomp.determinant(result);
}

template <typename RT>
inline M2_arrayintOrNull LU(const DMat<RT>& A, DMat<RT>& L, DMat<RT>& U)
{
  std::vector<size_t> perm;
  DMatLinAlg<RT> LUdecomp(A);
  LUdecomp.matrixPLU(perm, L, U);
  return stdvector_to_M2_arrayint(perm);
}

/*
  Cases for strategy:
  00 lower triangular (forward substitution)
  01 lower triangular, assume 1 on diagonal
  10 upper triangular (backward substitution)
  11 upper triangular, assume 1 on diagonal
  Note: the rest of the matrix need not be 0 filled.
*/
template <typename RT>
void triangularSolve(DMat<RT>& Lv, DMat<RT>& x, int m, int strategy)
{
  switch (strategy)
    {
      case 0:
        for (size_t i = 0; i < m; i++)
          {
            auto& a = x.entry(i, 0);
            x.ring().divide(a, Lv.entry(i, m), Lv.entry(i, i));
            x.ring().negate(a, a);
            MatElementaryOps<DMat<RT>>::column_op(Lv, m, a, i);
            x.ring().negate(a, a);
          }
        break;
      case 1:
        for (size_t i = 0; i < m; i++)
          {
            auto& a = x.entry(i, 0);
            x.ring().negate(a, Lv.entry(i, m));
            MatElementaryOps<DMat<RT>>::column_op(Lv, m, a, i);
            x.ring().negate(a, a);
          }
        break;
      case 2:
        for (size_t i = 1; i < m + 1; i++)
          {
            auto& a = x.entry(m - i, 0);
            x.ring().divide(a, Lv.entry(m - i, m), Lv.entry(m - i, m - i));
            x.ring().negate(a, a);
            MatElementaryOps<DMat<RT>>::column_op(Lv, m, a, m - i);
            x.ring().negate(a, a);
          }
        break;
      case 3:
        for (size_t i = 1; i < m + 1; i++)
          {
            auto& a = x.entry(m - i, 0);
            x.ring().negate(a, Lv.entry(m - i, m));
            MatElementaryOps<DMat<RT>>::column_op(Lv, m, a, m - i);
            x.ring().negate(a, a);
          }
        break;
    }
}

template <typename RT>
M2_arrayintOrNull LUincremental(std::vector<size_t>& P,
                               DMat<RT>& LU,
                               const DMat<RT>& v,
                               int m)
{
  size_t n = LU.numRows();

  for (size_t j = 0; j < n; j++)
    LU.ring().set(LU.entry(j, m), v.entry(P[j], 0));

  DMat<RT> x{LU.ring(), n, 1};
  triangularSolve(LU, x, m, 1);
  for (size_t i = 0; i < m; i++)
    LU.ring().set(LU.entry(i, m), x.entry(i, 0));

  int pivotPosition = -1;
  for (size_t j = m; j < n; j++)
    if (!LU.ring().is_zero(LU.entry(j, m)))
      {
        pivotPosition = j;
        break;
      }
  if (pivotPosition == -1) return stdvector_to_M2_arrayint(P);

  MatElementaryOps<DMat<RT>>::interchange_rows(LU, pivotPosition, m);
  std::swap(P[pivotPosition], P[m]);

  for (int j = m + 1; j < n; j++)
    LU.ring().divide(LU.entry(j, m), LU.entry(j, m), LU.entry(m, m));

  return stdvector_to_M2_arrayint(P);
}

template <typename RT>
inline size_t rank(const DMat<RT>& A)
{
  DMatLinAlg<RT> LUdecomp(A);
  return LUdecomp.rank();
}

template <typename RT>
inline M2_arrayintOrNull rankProfile(const DMat<RT>& A, bool row_profile)
{
  std::vector<size_t> profile;
  if (row_profile)
    {
      DMat<RT> B(A.ring(), A.numColumns(), A.numRows());
      MatrixOps::transpose(A, B);
      DMatLinAlg<RT> LUdecomp(B);
      LUdecomp.columnRankProfile(profile);
      return stdvector_to_M2_arrayint(profile);
    }

  DMatLinAlg<RT> LUdecomp(A);
  LUdecomp.columnRankProfile(profile);
  return stdvector_to_M2_arrayint(profile);
}

template <typename RT>
inline bool inverse(const DMat<RT>& A, DMat<RT>& result_inv)
{
  DMatLinAlg<RT> LUdecomp(A);
  return LUdecomp.inverse(result_inv);
}

template <typename RT>
inline size_t nullSpace(const DMat<RT>& A, DMat<RT>& result_nullspace)
{
  DMatLinAlg<RT> LUdecomp(A);
  return LUdecomp.kernel(result_nullspace);
}

template <typename RT>
inline bool solveLinear(const DMat<RT>& A, const DMat<RT>& B, DMat<RT>& X)
{
  DMatLinAlg<RT> LUdecomp(A);
  return LUdecomp.solve(B, X);
}

template <typename RT>
inline bool solveInvertible(const DMat<RT>& A,
                            const DMat<RT>& B,
                            DMat<RT>& X)
{
  DMatLinAlg<RT> LUdecomp(A);
  return LUdecomp.solveInvertible(B, X);
}

}  // namespace MatrixOps

#endif
