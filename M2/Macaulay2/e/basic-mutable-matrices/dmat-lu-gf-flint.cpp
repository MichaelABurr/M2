#include "basic-mutable-matrices/dmat-gf-flint.hpp"
#include "basic-mutable-matrices/dmat-lu.hpp"
#include "basic-rings/aring-GF-flint.hpp"

template <>
void solveLowerTriangular<DMat<M2::ARingGFFlint>>(
    const DMat<M2::ARingGFFlint>& LU,
    const DMat<M2::ARingGFFlint>& B,
    DMat<M2::ARingGFFlint>& X)
{
  fq_zech_mat_solve_tril(X.fq_zech_mat(),
                         LU.fq_zech_mat(),
                         B.fq_zech_mat(),
                         1,
                         LU.ring().flintContext());
}

template <>
void solveUpperTriangular<DMat<M2::ARingGFFlint>>(
    const DMat<M2::ARingGFFlint>& LU,
    const DMat<M2::ARingGFFlint>& B,
    DMat<M2::ARingGFFlint>& X)
{
  fq_zech_mat_solve_triu(X.fq_zech_mat(),
                         LU.fq_zech_mat(),
                         B.fq_zech_mat(),
                         0,
                         LU.ring().flintContext());
}
