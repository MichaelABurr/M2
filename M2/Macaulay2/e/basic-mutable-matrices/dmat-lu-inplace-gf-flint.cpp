// Copyright 2026  The Macaulay2 Authors

#include "basic-mutable-matrices/dmat-lu-inplace-gf-flint.hpp"

#include <M2/gc-include.h>

#include "basic-mutable-matrices/dmat-gf-flint.hpp"
#include "basic-rings/aring-GF-flint.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <flint/perm.h>
#pragma GCC diagnostic pop

template <>
void DMatLUinPlace<M2::ARingGFFlint>::computeLU()
{
  if (mIsDone) return;

  mp_limb_signed_t* perm = newarray_atomic(mp_limb_signed_t, mLU->numRows());
  fq_zech_mat_lu(perm, mLU->fq_zech_mat(), false, ring().flintContext());

  mPerm.clear();
  for (long i = 0; i < mLU->numRows(); i++) mPerm.push_back(perm[i]);
  mSign = (_perm_parity(perm, mLU->numRows()) == 0);
  freemem(perm);

  LUUtil<RingType>::computePivotColumns(*mLU, mPivotColumns);
  mIsDone = true;
}
