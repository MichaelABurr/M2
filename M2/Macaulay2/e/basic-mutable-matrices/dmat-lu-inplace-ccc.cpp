// Copyright 2026  The Macaulay2 Authors

#include "basic-mutable-matrices/dmat-lu-inplace-ccc.hpp"

#include "basic-rings/aring-CCC.hpp"
#include "basic-rings/aring-RRR.hpp"

template <>
size_t DMatLUinPlace<M2::ARingCCC>::findPivot(size_t row, size_t col)
{
  const M2::ARingRRR& RR = ring().real_ring();
  M2::ARingRRR::Element largest(RR), abs(RR);
  size_t best_row_so_far = static_cast<size_t>(-1);

  RR.set_zero(largest);
  for (size_t i = row; i < mLU->numRows(); i++)
    {
      ring().abs(abs, mLU->entry(i, col));
      if (RR.compare_elems(abs, largest) > 0)
        {
          best_row_so_far = i;
          RR.set(largest, abs);
        }
    }
  return best_row_so_far;
}
