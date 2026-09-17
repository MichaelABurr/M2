// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LU_QQ_HPP_
#define M2_BASIC_MUTMAT_DMAT_LU_QQ_HPP_

#include <gmp.h>

#include <cstddef>
#include <vector>

namespace M2 {
class ARingQQGMP;
}

template <typename RingType>
class DMat;

template <class RingType>
class DMatLinAlg;

template <>
class DMatLinAlg<M2::ARingQQGMP>
{
 public:
  using RingType = M2::ARingQQGMP;
  using Mat = DMat<RingType>;

  explicit DMatLinAlg(const Mat& A);

  size_t rank();
  void determinant(__mpq_struct& result_det);
  void columnRankProfile(std::vector<size_t>& profile);
  void matrixPLU(std::vector<size_t>& P, Mat& L, Mat& U);
  bool solve(const Mat& B, Mat& X);
  bool solveInvertible(const Mat& B, Mat& X);
  bool inverse(Mat& result_inv);
  size_t kernel(Mat& result_nullspace);

 private:
  const Mat& mInputMatrix;
};

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
