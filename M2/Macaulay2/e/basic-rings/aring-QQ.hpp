// Copyright 2013 Michael E. Stillman.

/// \file aring-QQ.hpp
/// \brief Selects the default ARing implementation for rational numbers.
///
/// This header includes the available FLINT and GMP rational implementations
/// and aliases ARingQQ to the implementation currently used by the engine.

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
