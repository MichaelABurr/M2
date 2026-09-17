#include "basic-mutable-matrices/dmat-gf-flint-big.hpp"
#include "basic-mutable-matrices/dmat-lu.hpp"
#include "basic-rings/aring-GF-flint-big.hpp"


template <>
void solveLowerTriangular<DMat<M2::ARingGFFlintBig>>(
    const DMat<M2::ARingGFFlintBig>& LU,
    const DMat<M2::ARingGFFlintBig>& B,
    DMat<M2::ARingGFFlintBig>& X)
{
  fq_nmod_mat_solve_tril(X.fq_nmod_mat(),
                         LU.fq_nmod_mat(),
                         B.fq_nmod_mat(),
                         1,
                         LU.ring().flintContext());
}


template <>
void solveUpperTriangular<DMat<M2::ARingGFFlintBig>>(
    const DMat<M2::ARingGFFlintBig>& LU,
    const DMat<M2::ARingGFFlintBig>& B,
    DMat<M2::ARingGFFlintBig>& X)
{
  fq_nmod_mat_solve_triu(X.fq_nmod_mat(),
                         LU.fq_nmod_mat(),
                         B.fq_nmod_mat(),
                         0,
                         LU.ring().flintContext());
}
