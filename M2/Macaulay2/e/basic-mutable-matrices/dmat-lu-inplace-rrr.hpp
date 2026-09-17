// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LU_INPLACE_RRR_HPP_
#define M2_BASIC_MUTMAT_DMAT_LU_INPLACE_RRR_HPP_

#include "basic-mutable-matrices/dmat-lu-inplace.hpp"

namespace M2 {
class ARingRRR;
}

template <>
size_t DMatLUinPlace<M2::ARingRRR>::findPivot(size_t row, size_t col);

#endif
