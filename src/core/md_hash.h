#pragma once

#include <core/md_compiler.h>
#include <core/md_intrinsics.h>
#include <core/md_allocator.h>
#include <core/md_str.h>

#define XXH_INLINE_ALL
#include <xxhash.h>

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus

constexpr uint64_t md_hash_str64(const char* str, size_t len) {
    uint64_t h = 0x9e3779b97f4a7c15ull; // seed

    for (size_t i = 0; i < len; ++i) {
		uint64_t c = static_cast<unsigned char>(str[i]);
        h ^= c;
        h *= 0xbf58476d1ce4e5b9ull;
        h ^= (h >> 27);
    }

    // final avalanche (splitmix64 style)
    h ^= h >> 30;
    h *= 0xbf58476d1ce4e5b9ull;
    h ^= h >> 27;
    h *= 0x94d049bb133111ebull;
    h ^= h >> 31;

    return h;
}

constexpr uint32_t md_hash_str32(const char* str, size_t len) {
    uint64_t h = md_hash_str64(str, len);
    return static_cast<uint32_t>(h ^ (h >> 32));
}

#define HASH_STR_LIT(str)   md_hash_str32(str"",sizeof(str))
#define HASH_STR_LIT32(str) md_hash_str32(str"",sizeof(str))
#define HASH_STR_LIT64(str) md_hash_str64(str"",sizeof(str))

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

// Standard default hash functions.
static inline uint32_t md_hash32(const void* input, size_t len, uint32_t seed) {
    return XXH32(input, len, seed);
}

static inline uint64_t md_hash64(const void* input, size_t len, uint64_t seed) {
    return XXH64(input, len, seed);
}

static inline uint64_t md_hash64_combine(uint64_t a, uint64_t b) {
    return a ^ (b + 0x9e3779b97f4a7c15ull + (a << 6) + (a >> 2));
}

static inline uint64_t md_hash64_str(str_t str, uint64_t seed) { return XXH64(str.ptr, str.len, seed); }


// This is very inspired by the Ourmachinery hash API

// Sentinel key value that represents a removed key. Note that this value can never be used as an
// actual key, or there will be trouble.
#define MD_HASH_TOMBSTONE 0xfffffffffffffffeULL

// Sentinel key value that represents an unused key. Note that this value can never be used as an
// actual key, or there will be trouble.
#define MD_HASH_UNUSED    0xffffffffffffffffULL

#define MD_HASHMAP_T(V)               \
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
#define md_hashset_skip_index(h, index) ((h)->keys[index] >= MD_HASH_TOMBSTONE)
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

// Adds the key to the hash set.
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
// The idea is to increase the entropy by mising in some of the seed
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
