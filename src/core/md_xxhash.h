#pragma once

// xxHash, inline, as the few files which hash in a loop of their own take it (md_hash.c, md_atom.c). Not for headers:
// it includes the intrinsics.

#define XXH_INLINE_ALL
// Read the input through memcpy. xxHash's default on GCC reads it through a uint64_t (or uint32_t)
// lvalue, which may not alias what the caller stored: hashing a uint32_t[2] that way is undefined,
// and GCC 13 at -O2 removes the stores to it, so the hash is of whatever was on the stack.
#ifndef XXH_FORCE_MEMORY_ACCESS
#define XXH_FORCE_MEMORY_ACCESS 0
#endif
#include <xxhash.h>
