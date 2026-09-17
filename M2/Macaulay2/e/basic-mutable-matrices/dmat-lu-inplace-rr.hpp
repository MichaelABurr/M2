// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LU_INPLACE_RR_HPP_
#define M2_BASIC_MUTMAT_DMAT_LU_INPLACE_RR_HPP_

#include "basic-mutable-matrices/dmat-lu-inplace.hpp"

namespace M2 {
class ARingRR;
}

template <>
void DMatLUinPlace<M2::ARingRR>::computeLU();

#endif
