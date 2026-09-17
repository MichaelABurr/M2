// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LU_ZZP_FFPACK_HPP_
#define M2_BASIC_MUTMAT_DMAT_LU_ZZP_FFPACK_HPP_

#include <cstddef>
#include <vector>

namespace M2 {
class ARingZZpFFPACK;
}

template <typename RingType>
class DMat;

template <class RingType>
class DMatLinAlg;

template <>
class DMatLinAlg<M2::ARingZZpFFPACK>
{
 public:
  using RingType = M2::ARingZZpFFPACK;
  using Mat = DMat<RingType>;

  explicit DMatLinAlg(const Mat& A);

  size_t rank();
  void determinant(double& result);
  void columnRankProfile(std::vector<size_t>& profile);
  void matrixPLU(std::vector<size_t>& P, Mat& L, Mat& U);
  bool solve(const Mat& B, Mat& X);
  bool solveInvertible(const Mat& B, Mat& X);
  bool inverse(Mat& X);
  size_t kernel(Mat& X);

 private:
  const Mat& mLU;
};

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
