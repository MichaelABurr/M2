// Copyright 2013 Michael E. Stillman.

/// AI:
/// \file aring-QQ.hpp
/// \brief Selects the default ARing implementation for rational numbers.
///
/// This small coordination header includes both rational ARing backends and
/// defines the `ARingQQ` typedef used throughout the rest of `basic-rings`.
/// Code that dispatches on the abstract rational ARing name, including
/// `aring-glue.hpp` and `aring-translate.hpp`, should include this header
/// rather than committing to a concrete backend directly.
///
/// Changing the typedef changes which rational implementation is used by
/// `RingQQ::create()` and by promotion/lift code that names `ARingQQ`.  Backend
/// headers still remain available for tests or factory paths that need to
/// instantiate a specific GMP or FLINT implementation explicitly.

#ifndef M2_BASIC_RINGS_ARING_QQ_HPP_
#define M2_BASIC_RINGS_ARING_QQ_HPP_

#include "basic-rings/aring-QQ-flint.hpp"
#include "basic-rings/aring-QQ-gmp.hpp"

namespace M2 {
class ARingQQFlint;
class ARingQQGMP;

typedef ARingQQGMP ARingQQ;
};

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
