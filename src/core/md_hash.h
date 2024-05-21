#pragma once

#include <stddef.h>
#include <stdint.h>
#include <core/md_compiler.h>
#include <core/md_intrinsics.h>
#include <core/md_allocator.h>

#if 0
// This does not seem to work on MSVC
// https://stackoverflow.com/questions/2826559/compile-time-preprocessor-hashing-of-string
// http://lolengine.net/blog/2011/12/20/cpp-constant-string-hash
#define HASH_H1(s,i,x)   (x*65599u+(uint8_t)s[(i)<sizeof(s)?sizeof(s)-1-(i):sizeof(s)])
#define HASH_H4(s,i,x)   HASH_H1 (s,i,HASH_H1 (s,i+1, HASH_H1 (s,i+2,  HASH_H1 (s,i+3,x))))
#define HASH_H16(s,i,x)  HASH_H4 (s,i,HASH_H4 (s,i+4, HASH_H4 (s,i+8,  HASH_H4 (s,i+12,x))))
#define HASH_H64(s,i,x)  HASH_H16(s,i,HASH_H16(s,i+16,HASH_H16(s,i+32, HASH_H16(s,i+48,x))))
#define HASH_H256(s,i,x) HASH_H64(s,i,HASH_H64(s,i+64,HASH_H64(s,i+128,HASH_H64(s,i+192,x))))
// Hash a string literal (at compile time)
#define HASH_STR_LIT(s)  ((uint32_t)(HASH_H256(s"",0,0)^(HASH_H256(s"",0,0)>>16)))
#endif

#ifdef __cplusplus
/*
https://github.com/LordJZ/consthash
The MIT License(MIT)

Copyright(c) 2015

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

static uint32_t constexpr md_crc32_tab[] = {
    0x00000000, 0x77073096, 0xee0e612c, 0x990951ba, 0x076dc419, 0x706af48f, 0xe963a535, 0x9e6495a3, 0x0edb8832, 0x79dcb8a4, 0xe0d5e91e, 0x97d2d988, 0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91,
    0x1db71064, 0x6ab020f2, 0xf3b97148, 0x84be41de, 0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7, 0x136c9856, 0x646ba8c0, 0xfd62f97a, 0x8a65c9ec, 0x14015c4f, 0x63066cd9, 0xfa0f3d63, 0x8d080df5,
    0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172, 0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b, 0x35b5a8fa, 0x42b2986c, 0xdbbbc9d6, 0xacbcf940, 0x32d86ce3, 0x45df5c75, 0xdcd60dcf, 0xabd13d59,
    0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116, 0x21b4f4b5, 0x56b3c423, 0xcfba9599, 0xb8bda50f, 0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924, 0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d,
    0x76dc4190, 0x01db7106, 0x98d220bc, 0xefd5102a, 0x71b18589, 0x06b6b51f, 0x9fbfe4a5, 0xe8b8d433, 0x7807c9a2, 0x0f00f934, 0x9609a88e, 0xe10e9818, 0x7f6a0dbb, 0x086d3d2d, 0x91646c97, 0xe6635c01,
    0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e, 0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457, 0x65b0d9c6, 0x12b7e950, 0x8bbeb8ea, 0xfcb9887c, 0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65,
    0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2, 0x4adfa541, 0x3dd895d7, 0xa4d1c46d, 0xd3d6f4fb, 0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0, 0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9,
    0x5005713c, 0x270241aa, 0xbe0b1010, 0xc90c2086, 0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f, 0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4, 0x59b33d17, 0x2eb40d81, 0xb7bd5c3b, 0xc0ba6cad,
    0xedb88320, 0x9abfb3b6, 0x03b6e20c, 0x74b1d29a, 0xead54739, 0x9dd277af, 0x04db2615, 0x73dc1683, 0xe3630b12, 0x94643b84, 0x0d6d6a3e, 0x7a6a5aa8, 0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1,
    0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe, 0xf762575d, 0x806567cb, 0x196c3671, 0x6e6b06e7, 0xfed41b76, 0x89d32be0, 0x10da7a5a, 0x67dd4acc, 0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5,
    0xd6d6a3e8, 0xa1d1937e, 0x38d8c2c4, 0x4fdff252, 0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b, 0xd80d2bda, 0xaf0a1b4c, 0x36034af6, 0x41047a60, 0xdf60efc3, 0xa867df55, 0x316e8eef, 0x4669be79,
    0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236, 0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f, 0xc5ba3bbe, 0xb2bd0b28, 0x2bb45a92, 0x5cb36a04, 0xc2d7ffa7, 0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d,
    0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x026d930a, 0x9c0906a9, 0xeb0e363f, 0x72076785, 0x05005713, 0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38, 0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7, 0x0bdbdf21,
    0x86d3d2d4, 0xf1d4e242, 0x68ddb3f8, 0x1fda836e, 0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777, 0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c, 0x8f659eff, 0xf862ae69, 0x616bffd3, 0x166ccf45,
    0xa00ae278, 0xd70dd2ee, 0x4e048354, 0x3903b3c2, 0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db, 0xaed16a4a, 0xd9d65adc, 0x40df0b66, 0x37d83bf0, 0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
    0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6, 0xbad03605, 0xcdd70693, 0x54de5729, 0x23d967bf, 0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94, 0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d};

constexpr uint32_t md_crc32impl(uint32_t prevCrc, const char* str, size_t size) { return !size ? prevCrc : md_crc32impl((prevCrc >> 8) ^ md_crc32_tab[(prevCrc ^ *str) & 0xff], str + 1, size - 1); }
constexpr uint32_t md_crc32(const char* ptr, size_t size) { return md_crc32impl(0xffffffff, ptr, size) ^ 0xffffffff; }

#define HASH_STR_LIT(str) md_crc32(str"",sizeof(str))

extern "C" {
#endif

// https://en.wikipedia.org/wiki/Halton_sequence
static inline float md_halton(int index, int base) {
    float f = 1;
    float r = 0;
    const float ifb = 1.f / base;
    while (index > 0) {
        f = f * ifb;
        r = r + f * (float)(index % base);
        index = (int)(index * ifb);
    }
    return r;
}


uint32_t md_hash32(const void* input, size_t len, uint32_t seed);
uint64_t md_hash64(const void* input, size_t len, uint64_t seed);

// This is very inspired by the Ourmachinery hash API

// Sentinel key value that represents a removed key. Note that this value can never be used as an
// actual key, or there will be trouble.
#define MD_HASH_TOMBSTONE 0xfffffffffffffffeULL

// Sentinel key value that represents an unused key. Note that this value can never be used as an
// actual key, or there will be trouble.
#define MD_HASH_UNUSED    0xffffffffffffffffULL

#define MD_HASHMAP_T(V)                   \
    {                                     \
        uint32_t num_buckets;             \
        uint32_t num_used;                \
        uint64_t *keys;                   \
        V *values;                        \
        struct md_allocator_i *allocator; \
    }

typedef struct {
    uint32_t num_buckets;
    uint32_t num_used;
    uint64_t *keys;
    struct md_allocator_i *allocator;
} md_hashset_t;

// Returns *true* if index should be skipped when iterating over the hash table (i.e. if it
// contains a `TOMBSTONE` or `UNUSED` value).
//
// ~~~c
// for (uint32_t i = 0; i < ht.num_buckets; ++i) {
//     if (md_hashmap_skip_index(&ht, i))
//         continue;
//     ....
// }
// ~~~
#define md_hashmap_skip_index(h, index) ((h)->keys[index] >= MD_HASH_TOMBSTONE)

// The opposite of [[md_hashmap_skip_index()]].
#define md_hashmap_use_index(h, index) (!md_hashmap_skip_index(h, index))

// Returns the index in the hash table where the `key` is stored, or `UINT32_MAX` if `key` is
// not in the hash table.
#define md_hashmap_index(h, key) md_hashmap__index((h)->keys, (h)->num_buckets, key)

// Returns *true* if the hash table contains `key`.
#define md_hashmap_has(h, key) (md_hashmap_index(h, key) != UINT32_MAX)

#define md_hashmap_get(h, key) md_hashmap__get((h)->values, sizeof(*(h)->values), (h)->keys, (h)->num_buckets, key)
#define md_hashset_get(h, key) md_hashmap__get((h)->keys,   sizeof(*(h)->keys),   (h)->keys, (h)->num_buckets, key)

// Adds the pair `(key, value)` to the hash table.
#define md_hashmap_add(h, key, value)                                                                                                                                           \
do {                                                                                                                                                                            \
    uint32_t _slot_idx = md_hashmap__add(&(h)->keys, &(h)->num_used, &(h)->num_buckets, key, (void **)&(h)->values, sizeof(*(h)->values), (h)->allocator, __FILE__, __LINE__);  \
    ASSERT(_slot_idx < (h)->num_buckets);                                                                                                                                       \
    (h)->values[_slot_idx] = value;                                                                                                                                             \
} while(0);

// Adds the pair `(key, value)` to the hash set.
#define md_hashset_add(h, key) md_hashmap__add(&(h)->keys, &(h)->num_used, &(h)->num_buckets, key, NULL, 0, (h)->allocator, __FILE__, __LINE__)

// Removes the value with the specified `key`. Removing a value never reallocates the hash table.
#define md_hashmap_remove(h, key) md_hashmap__remove((h)->keys, (h)->num_buckets, key)
#define md_hashset_remove(h, key) md_hashmap__remove((h)->keys, (h)->num_buckets, key)

// Clears the hash table. This does not free any memory. The array will keep the same number of
// buckets, but all the keys and values will be erased.
#define md_hashmap_clear(h) (MEMSET((h)->keys, 0xFF, (h)->num_buckets * sizeof(*(h)->keys)), (h)->num_used = 0)

#define md_hashmap_free(h) (md_hashmap__free((uint64_t **)&(h)->keys, &(h)->num_used, &(h)->num_buckets, (void **)&(h)->values, sizeof(*(h)->values), (h)->allocator, __FILE__, __LINE__))
#define md_hashset_free(h) (md_hashmap__free((uint64_t **)&(h)->keys, &(h)->num_used, &(h)->num_buckets, NULL, 0, (h)->allocator, __FILE__, __LINE__))

// Copies the hash table `from` to `h`.
#define md_hashmap_copy(h, from) \
    (md_hashmap_free(h),         \
        *h = *from,              \
        md_hashmap__copy((uint64_t **)&(h)->keys, (h)->num_buckets, (void **)&(h)->values, sizeof(*(h)->values), (h)->allocator, __FILE__, __LINE__))

#define md_hashset_copy(h, from) \
    (md_hashmap_free(h),         \
        *h = *from,              \
        md_hashmap__copy((uint64_t **)&(h)->keys, (h)->num_buckets, NULL, 0, (h)->allocator, __FILE__, __LINE__))

// Compute the number of buckets required to accommodate `n` elements.
static inline uint32_t md_hashmap__buckets_for_elements(uint32_t n) {
    return next_power_of_two32(MAX(16, (uint32_t)(n / 0.7f + 1)));
}

// Reserves space in the hash table for `num_elements` new elements. The hash table will be resized
// so that its size is a power-of-two and adding `num_elements` elements will keep it below the
// target fillrate (currently 70 %). If the hash table already has room for `num_elements` new
// elements, this function is a NOP.
#define md_hashmap_reserve(h, num_elements) \
    md_hashmap__grow_to((uint64_t **)&(h)->keys, &(h)->num_used, &(h)->num_buckets, (void **)&(h)->values, sizeof(*(h)->values), (h)->allocator, __FILE__, __LINE__, md_hashmap__buckets_for_elements((h)->num_used + (uint32_t)num_elements))

#define md_hashset_reserve(h, num_elements) \
    md_hashmap__grow_to((uint64_t **)&(h)->keys, &(h)->num_used, &(h)->num_buckets, NULL, 0, (h)->allocator, __FILE__, __LINE__, md_hashmap__buckets_for_elements((h)->num_used + (uint32_t)num_elements))


// Standard hash and set types
//
// Defines some commonly used set and hash types.

// Maps from an `uint64_t` key to an `uint64_t` value.
typedef struct MD_HASHMAP_T(uint64_t) md_hashmap64_t;

// Maps from an `uint64_t` key to a `uint32_t` value.
typedef struct MD_HASHMAP_T(uint32_t) md_hashmap32_t;


// Maps a key to a hash table index.
//
// ???: Should we do more to mix up the bits of `key` here? We want to make sure that we don't have
// any bad patterns in the bits, but this function gets called a lot, so we don't want to add too
// many steps here or performance will suffer.
// This is currently not properly tested and evaluated on real data
// The idéa is to increase the entropy by mising in some of the seed
static inline uint32_t md_hashmap__first_index(uint64_t key, uint32_t num_buckets, uint64_t seed) {
    (void)seed;
    uint32_t *v = (uint32_t *)&key;
    return ((v[0] ^ v[1])) & (num_buckets - 1);
}

// Frees the arrays allocated by the hash.
static inline void md_hashmap__free(uint64_t **keys_ptr, uint32_t *num_used_ptr, uint32_t *num_buckets_ptr, void **values_ptr, uint64_t value_bytes, md_allocator_i *allocator, const char *file, uint32_t line) {
    uint64_t *keys = *keys_ptr;
    const uint32_t num_buckets = *num_buckets_ptr;
    const size_t bytes_to_free = num_buckets ? num_buckets * (sizeof(*keys) + value_bytes) : 0;
    if (allocator)
        allocator->realloc(allocator->inst, keys, bytes_to_free, 0, file, line);
    *keys_ptr = 0;
    if (values_ptr)
        *values_ptr = 0;
    *num_used_ptr = 0;
    *num_buckets_ptr = 0;
}

// Returns the index of the `key` in the array `keys` or `UINT32_MAX` if `key` is not in the array.
static inline uint32_t md_hashmap__index(const uint64_t *keys, uint32_t num_buckets, uint64_t key) {
    if (!num_buckets || key >= MD_HASH_TOMBSTONE)
        return UINT32_MAX;

    const uint32_t max_distance = num_buckets;

    uint32_t i = md_hashmap__first_index(key, num_buckets, (uint64_t)keys);
    uint32_t distance = 0;
    while (keys[i] != key) {
        if (distance > max_distance)
            return UINT32_MAX;
        if (keys[i] == MD_HASH_UNUSED)
            return UINT32_MAX;
        i = (i + 1) & (num_buckets - 1);
        ++distance;
    }
    return i;
}

// Removes entry in hash table
static inline void md_hashmap__remove(uint64_t* keys, uint32_t num_buckets, uint64_t key) {
    uint32_t idx = md_hashmap__index(keys, num_buckets, key);
    if (idx != UINT32_MAX) {
        keys[idx] = MD_HASH_TOMBSTONE;
    }
}

// Returns the value address of the `key` if it exists, otherwise NULL
static inline void* md_hashmap__get(void* value_ptr, size_t value_bytes, const uint64_t *keys, uint32_t num_buckets, uint64_t key) {
    uint32_t idx = md_hashmap__index(keys, num_buckets, key);
    return idx == UINT32_MAX ? NULL : (char*)value_ptr + idx * value_bytes;
}

// Tries to add `key` to the `keys` array (without reallocating it). If it succeeds, it will return
// the index of the added key and otherwise `UINT32_MAX`. You should only call this if
// [[md_hashmap_index()]] has failed to find the key.
static inline uint32_t md_hashmap__add_no_grow(const uint64_t *keys, uint32_t num_buckets, uint64_t key) {
    const uint32_t max_distance = num_buckets;

    if (!num_buckets)
        return UINT32_MAX;

    uint32_t i = md_hashmap__first_index(key, num_buckets, (uint64_t)keys);
    uint32_t distance = 0;
    while (keys[i] < MD_HASH_TOMBSTONE) {
        if (distance > max_distance)
            return UINT32_MAX;
        i = (i + 1) & (num_buckets - 1);
        ++distance;
    }
    return i;
}

static inline void md_hashmap__grow(uint64_t **keys_ptr, uint32_t *num_used_ptr, uint32_t *num_buckets_ptr, void **values_ptr, size_t value_bytes, md_allocator_i *allocator, const char *file, uint32_t line);

// Adds `key` to the hash table, growing it as necessary, and returns its index. If the hash table
// needs to grow and `allocator` is NULL, an error will be generated.
static inline uint32_t md_hashmap__add(uint64_t **keys_ptr, uint32_t *num_used_ptr, uint32_t *num_buckets_ptr, uint64_t key, void **values_ptr, size_t value_bytes, md_allocator_i *allocator, const char *file, uint32_t line) {
    uint32_t i = md_hashmap__index(*keys_ptr, *num_buckets_ptr, key);
    if (i != UINT32_MAX)
        return i;

    if ((!*num_buckets_ptr || (float)*num_used_ptr / (float)*num_buckets_ptr > 0.7f) && allocator)
        md_hashmap__grow(keys_ptr, num_used_ptr, num_buckets_ptr, values_ptr, value_bytes, allocator, file, line);

    ASSERT(key < MD_HASH_TOMBSTONE && "Invalid Key");

    i = md_hashmap__add_no_grow(*keys_ptr, *num_buckets_ptr, key);
    ASSERT(i != UINT32_MAX && "Hash table is full!");

    uint64_t *keys = *keys_ptr;
    keys[i] = key;
    ++*num_used_ptr;
    return i;
}

// Grows the hash table to the specified size
static inline void md_hashmap__grow_to(uint64_t **keys_ptr, uint32_t *num_used_ptr, uint32_t *num_buckets_ptr, void **values_ptr, size_t value_bytes, md_allocator_i *allocator, const char *file, uint32_t line, uint32_t new_buckets) {
    const uint64_t *keys = *keys_ptr;
    const uint32_t num_buckets = *num_buckets_ptr;
    const void *values = values_ptr ? *values_ptr : 0;

    if (num_buckets >= new_buckets)
        return;

    // (new_buckets + 1) for the default value and key which is stored at [-1]
    const size_t to_alloc = new_buckets * (sizeof(uint64_t) + value_bytes);
    uint64_t *new_keys = (uint64_t*)allocator->realloc(allocator->inst, 0, 0, to_alloc, file, line);
    MEMSET(new_keys, 0xFF, new_buckets * sizeof(uint64_t));
    void *new_values = new_keys + new_buckets;
    uint32_t new_elements = 0;

    for (uint32_t i = 0; i < num_buckets; ++i) {
        if (keys[i] < MD_HASH_TOMBSTONE) {
            const uint32_t new_i = md_hashmap__add_no_grow(new_keys, new_buckets, keys[i]);
            new_keys[new_i] = keys[i];
            MEMCPY((char *)new_values + new_i * value_bytes, (char *)values + i * value_bytes, value_bytes);
            ++new_elements;
        }
    }
    if (num_buckets)
        md_hashmap__free(keys_ptr, num_used_ptr, num_buckets_ptr, values_ptr, value_bytes, allocator, file, line);

    *num_used_ptr = new_elements;
    *num_buckets_ptr = new_buckets;
    *keys_ptr = new_keys;
    if (values_ptr) {
        *values_ptr = new_values;
    }
}

// Grows the hash table.
static inline void md_hashmap__grow(uint64_t **keys_ptr, uint32_t *num_used_ptr, uint32_t *num_buckets_ptr, void **values_ptr, size_t value_bytes, md_allocator_i *alloc, const char *file, uint32_t line) {
    const uint32_t num_used = *num_used_ptr;
    const uint32_t new_buckets = MAX(16, next_power_of_two32(num_used * 2));
    md_hashmap__grow_to(keys_ptr, num_used_ptr, num_buckets_ptr, values_ptr, value_bytes, alloc, file, line, new_buckets);
}

// Copies the keys and values arrays to new pointers.
static inline void md_hashmap__copy(uint64_t **keys_ptr, uint32_t num_buckets, void **values_ptr, size_t value_bytes, md_allocator_i *alloc, const char *file, uint32_t line) {
    uint64_t *old_keys = *keys_ptr;
    const uint64_t size = num_buckets * (sizeof(old_keys) + value_bytes);
    *keys_ptr = (uint64_t*)alloc->realloc(alloc->inst, 0, 0, size, file, line);

    if (values_ptr) {
        *values_ptr = *keys_ptr + num_buckets;
    }

    MEMCPY(*keys_ptr, old_keys, size);
}

#ifdef __cplusplus
}
#endif
