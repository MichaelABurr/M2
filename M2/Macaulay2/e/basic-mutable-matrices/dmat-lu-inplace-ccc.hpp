// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LU_INPLACE_CCC_HPP_
#define M2_BASIC_MUTMAT_DMAT_LU_INPLACE_CCC_HPP_

#include "basic-mutable-matrices/dmat-lu-inplace.hpp"

namespace M2 {
class ARingCCC;
}

template <>
size_t DMatLUinPlace<M2::ARingCCC>::findPivot(size_t row, size_t col);

#endif
