#pragma once

// Everything that was here - popcnt, clz, ctz, bsf, bsr, the zero nibble and byte searches, next_power_of_two and the
// 64-bit atomics - resolves without an include on every supported compiler, and so is in md_common.h now (see THE RULE
// there). This header is kept for its includers, which also got <stdint.h> through it.

#include <core/md_common.h>

#include <stdint.h>
