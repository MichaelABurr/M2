// Copyright 2014  Michael E. Stillman

#ifndef M2_BASIC_MUTMAT_DMAT_LINALG_ALL__HPP_
#define M2_BASIC_MUTMAT_DMAT_LINALG_ALL__HPP_

#include "basic-rings/aring-ZZ-gmp.hpp"

#include "basic-mutable-matrices/dmat-gf-flint-big.hpp"
#include "basic-mutable-matrices/dmat-gf-flint.hpp"
#include "basic-mutable-matrices/dmat-qq-flint.hpp"
#include "basic-mutable-matrices/dmat-zz-flint.hpp"
#include "basic-mutable-matrices/dmat-zzp-flint.hpp"

#include "basic-mutable-matrices/dmat-lu-inplace-cc.hpp"
#include "basic-mutable-matrices/dmat-lu-inplace-ccc.hpp"
#include "basic-mutable-matrices/dmat-lu-inplace-gf-flint-big.hpp"
#include "basic-mutable-matrices/dmat-lu-inplace-gf-flint.hpp"
#include "basic-mutable-matrices/dmat-lu-inplace-rr.hpp"
#include "basic-mutable-matrices/dmat-lu-inplace-rrr.hpp"

#include "basic-mutable-matrices/dmat-linalg.hpp"
#include "basic-mutable-matrices/dmat-linalg-cc.hpp"
#include "basic-mutable-matrices/dmat-linalg-ccc.hpp"
#include "basic-mutable-matrices/dmat-linalg-gf-flint-big.hpp"
#include "basic-mutable-matrices/dmat-linalg-gf-flint.hpp"
#include "basic-mutable-matrices/dmat-linalg-qq-flint.hpp"
#include "basic-mutable-matrices/dmat-linalg-qq-gmp.hpp"
#include "basic-mutable-matrices/dmat-linalg-rr.hpp"
#include "basic-mutable-matrices/dmat-linalg-rrr.hpp"
#include "basic-mutable-matrices/dmat-linalg-zz-flint.hpp"
#include "basic-mutable-matrices/dmat-linalg-zz-gmp.hpp"
#include "basic-mutable-matrices/dmat-linalg-zzp-ffpack.hpp"
#include "basic-mutable-matrices/dmat-linalg-zzp-flint.hpp"
#include "basic-mutable-matrices/dmat-lu-qq.hpp"
#include "basic-mutable-matrices/dmat-lu-zzp-ffpack.hpp"
#include "basic-mutable-matrices/dmat-lu-zzp-flint.hpp"

// ConcreteRing's legacy ring_ZZp integration requires the complete Z_mod type.
#include "rings/ZZp.hpp"

#endif
