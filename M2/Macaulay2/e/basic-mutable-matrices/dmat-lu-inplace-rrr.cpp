// Copyright 2026  The Macaulay2 Authors

#include "basic-mutable-matrices/dmat-lu-inplace-rrr.hpp"

#include "basic-rings/aring-RRR.hpp"

template <>
size_t DMatLUinPlace<M2::ARingRRR>::findPivot(size_t row, size_t col)
{
  M2::ARingRRR::Element largest(ring()), abs(ring());
  size_t best_row_so_far = static_cast<size_t>(-1);

  ring().set_zero(largest);
  for (size_t i = row; i < mLU->numRows(); i++)
    {
      ring().abs(abs, mLU->entry(i, col));
      if (ring().compare_elems(abs, largest) > 0)
        {
          best_row_so_far = i;
          ring().copy(largest, abs);
        }
    }
  return best_row_so_far;
}
