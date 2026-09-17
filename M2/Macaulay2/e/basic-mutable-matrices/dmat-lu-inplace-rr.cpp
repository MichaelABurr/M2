// Copyright 2026  The Macaulay2 Authors

#include "basic-mutable-matrices/dmat-lu-inplace-rr.hpp"

#include "basic-rings/aring-RR.hpp"
#include "basic-rings/aring-CC.hpp"
#include "basic-mutable-matrices/lapack.hpp"
#include "exceptions.hpp"

template <>
void DMatLUinPlace<M2::ARingRR>::computeLU()
{
  if (mIsDone) return;

  int rows = static_cast<int>(mLU->numRows());
  int cols = static_cast<int>(mLU->numColumns());
  int info;
  int min = (rows <= cols) ? rows : cols;

  if (min == 0) return;

  int* perm = new int[min];
  std::vector<double> copyA = make_lapack_array(*mLU);

  dgetrf_(&rows, &cols, copyA.data(), &rows, perm, &info);

  if (info < 0)
    {
      delete[] perm;
      throw exc::engine_error("argument passed to dgetrf had an illegal value");
    }

  fill_from_lapack_array(copyA, *mLU);

  for (int i = 0; i < min; i++)
    {
      int thisloc = perm[i] - 1;
      if (i != thisloc)
        {
          mSign = not mSign;
          size_t tmp = mPerm[thisloc];
          mPerm[thisloc] = mPerm[i];
          mPerm[i] = tmp;
        }
    }

  LUUtil<RingType>::computePivotColumns(*mLU, mPivotColumns);
  mIsDone = true;

  delete[] perm;
}
