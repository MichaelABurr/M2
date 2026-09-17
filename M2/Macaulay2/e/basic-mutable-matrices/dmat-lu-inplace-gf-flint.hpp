// Copyright 2026  The Macaulay2 Authors

#ifndef M2_BASIC_MUTMAT_DMAT_LU_INPLACE_GF_FLINT_HPP_
#define M2_BASIC_MUTMAT_DMAT_LU_INPLACE_GF_FLINT_HPP_

#include "basic-mutable-matrices/dmat-lu-inplace.hpp"

namespace M2 {
class ARingGFFlint;
}

template <>
void DMatLUinPlace<M2::ARingGFFlint>::computeLU();

#endif
