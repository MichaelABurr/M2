// Copyright 2013 Michael E. Stillman.

/// AI:
/// \file aring-QQ.hpp
/// \brief Selects the default ARing implementation for rational numbers.
///
/// This header is the coordination point for rational ARings.  It includes the
/// available QQ backends and defines the `ARingQQ` typedef used by code that
/// wants the default rational implementation rather than a specific GMP or
/// FLINT class.
///
/// `aring-glue.hpp`, `aring-translate.hpp`, factory code, and tests include
/// this header when they need the selected QQ backend.  Changing the typedef
/// changes the implementation used by `RingQQ::create()` and by promotion/lift
/// code that names `ARingQQ`, while the backend headers remain available for
/// places that need an explicit implementation.

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
