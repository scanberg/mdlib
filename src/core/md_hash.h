#pragma once

#include <stddef.h>
#include <stdint.h>

// https://stackoverflow.com/questions/2826559/compile-time-preprocessor-hashing-of-string
// http://lolengine.net/blog/2011/12/20/cpp-constant-string-hash
#define _H1(s,i,x)   (x*65599u+(uint8_t)s[(i)<sizeof(s)?sizeof(s)-1-(i):sizeof(s)])
#define _H4(s,i,x)   H1(s,i,H1(s,i+1,H1(s,i+2,H1(s,i+3,x))))
#define _H16(s,i,x)  H4(s,i,H4(s,i+4,H4(s,i+8,H4(s,i+12,x))))
#define _H64(s,i,x)  H16(s,i,H16(s,i+16,H16(s,i+32,H16(s,i+48,x))))
#define _H256(s,i,x) H64(s,i,H64(s,i+64,H64(s,i+128,H64(s,i+192,x))))
// Hash a string literal (at compile time)
#define HASH_STR_LIT(s)  ((uint32_t)(H256(s"",0,0)^(H256(s"",0,0)>>16)))

#ifdef __cplusplus
extern "C" {
#endif

uint32_t md_hash32(const void* input, size_t len, uint32_t seed);
uint64_t md_hash64(const void* input, size_t len, uint64_t seed);

#ifdef __cplusplus
}
#endif
