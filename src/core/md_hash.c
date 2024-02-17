#include <core/md_hash.h>

#define XXH_INLINE_ALL
#include <xxhash.h>

uint32_t md_hash32(const void* input, size_t len, uint32_t seed) {
	return XXH32(input, len, seed);
}

uint64_t md_hash64(const void* input, size_t len, uint64_t seed) {
	return XXH64(input, len, seed);
}