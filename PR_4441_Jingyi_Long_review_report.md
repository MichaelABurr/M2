### Issue: Same-type promotion is not always defined

**Location:** `M2/Macaulay2/e/basic-rings/aring-translate.hpp:209`

The new template allows promotion between any two rings of the same C++ type, but such a promotion is not always defined. For `ZZp`, the prime `p` is a runtime parameter, so two `ZZp` instances may have different primes even though they have the same type. Their elements therefore cannot, in general, be converted by simply copying the underlying value.

#### Reproduction

Build the project in the `build` directory of the repository root, with cmake.

From the repository root, save the following as `.pr4441-promotion-probe.cpp`:

```cpp
#include <iostream>

#include "basic-rings/aring-translate.hpp"

int main()
{
  M2::ARingZZpFFPACK source(5);
  M2::ARingZZpFFPACK target(7);
  M2::ARingZZpFFPACK::ElementType input;
  M2::ARingZZpFFPACK::ElementType output;
  source.set(input, 4);
  const bool promoted = M2::mypromote(source, target, input, output);
  std::cout << promoted << " " << output << "\n";
}
```

Compile it against the configured `build` tree:

```sh
/usr/bin/c++ \
  -std=gnu++17 -O2 \
  -ffunction-sections -fdata-sections \
  -I M2/Macaulay2/e \
  -I M2/include \
  -I build/include \
  -I build/usr-host/include \
  .pr4441-promotion-probe.cpp \
  M2/Macaulay2/e/basic-rings/aring-ZZp-ffpack.cpp \
  -Wl,--gc-sections \
  -L/usr/local/lib \
  -lgivaro -lgmpxx -lgmp -lgc \
  -o /tmp/pr4441-promotion-probe
```

Run it:

```sh
/tmp/pr4441-promotion-probe
```

The output is:

```text
1 4
```

Thus promotion from characteristic 5 to characteristic 7 reports success and copies the underlying value `4`.
