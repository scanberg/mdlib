#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_allocator_i;
struct md_bitblock_t;

typedef struct md_bitfield_t {
    struct md_bitblock_t* bits;
    struct md_allocator_i* alloc;
    uint32_t beg_bit;
    uint32_t end_bit;
    uint32_t magic;
    uint32_t flags;
} md_bitfield_t;

typedef struct md_bitfield_iter_t {
    const md_bitfield_t* bf;
    uint64_t idx;
    uint32_t beg_bit;
    uint32_t end_bit;
} md_bitfield_iter_t;

#ifdef __cplusplus
extern "C" {
#endif

// The bitfield type is a semi-dense bitfield which only needs to hold the range of bits which are set.
// Because there is no explicit number of bits for the bitfield, some operations require explicit range parameters beg, end.

md_bitfield_t md_bitfield_create(struct md_allocator_i* alloc);
void md_bitfield_init           (md_bitfield_t* bf, struct md_allocator_i* alloc);
bool md_bitfield_free           (md_bitfield_t* bf);
bool md_bitfield_validate       (const md_bitfield_t* bf);

// Reset data to empty
void md_bitfield_reset          (md_bitfield_t* bf);

void md_bitfield_reserve_range  (md_bitfield_t* bf, uint64_t beg, uint64_t end);

bool md_bitfield_empty          (const md_bitfield_t* bf);
uint64_t md_bitfield_beg_bit    (const md_bitfield_t* bf);
uint64_t md_bitfield_end_bit    (const md_bitfield_t* bf);

void md_bitfield_set_range      (md_bitfield_t* bf, uint64_t beg, uint64_t end);
void md_bitfield_set_bit        (md_bitfield_t* bf, uint64_t bit_idx);
void md_bitfield_set_indices_u32(md_bitfield_t* bf, const uint32_t* indices, size_t num_indices);

// Clear bits (does not deallocate data or change the range)
void md_bitfield_clear          (md_bitfield_t* bf);
void md_bitfield_clear_range    (md_bitfield_t* bf, uint64_t beg, uint64_t end);
void md_bitfield_clear_bit      (md_bitfield_t* bf, uint64_t bit_idx);

void md_bitfield_or             (md_bitfield_t* dst, const md_bitfield_t* src_a, const md_bitfield_t* src_b);
void md_bitfield_or_inplace     (md_bitfield_t* a,   const md_bitfield_t* b);

void md_bitfield_and            (md_bitfield_t* dst, const md_bitfield_t* src_a, const md_bitfield_t* src_b);
void md_bitfield_and_inplace    (md_bitfield_t* a, const md_bitfield_t* b);

void md_bitfield_andnot         (md_bitfield_t* dst, const md_bitfield_t* src_a, const md_bitfield_t* src_b);
void md_bitfield_andnot_inplace (md_bitfield_t* a, const md_bitfield_t* b);

void md_bitfield_xor            (md_bitfield_t* dst, const md_bitfield_t* src_a, const md_bitfield_t* src_b);
void md_bitfield_xor_inplace    (md_bitfield_t* a, const md_bitfield_t* b);

// We need an explicit range or something here to 
void md_bitfield_not            (md_bitfield_t* dst, const md_bitfield_t* src, uint64_t beg, uint64_t end);
void md_bitfield_not_inplace    (md_bitfield_t* bf, uint64_t beg, uint64_t end);

void md_bitfield_copy           (md_bitfield_t* dst, const md_bitfield_t* src);

// Count number of bits set
size_t md_bitfield_popcount    (const md_bitfield_t* bf);
size_t md_bitfield_popcount_range(const md_bitfield_t* bf, uint64_t beg, uint64_t end);

// Test if bits are set
bool md_bitfield_test_bit   (const md_bitfield_t* bf, uint64_t idx);

// Test if any bits are set
bool md_bitfield_test_all   (const md_bitfield_t* bf);
bool md_bitfield_test_all_range (const md_bitfield_t* bf, uint64_t beg, uint64_t end);

bool md_bitfield_test_any   (const md_bitfield_t* bf);
bool md_bitfield_test_any_range (const md_bitfield_t* bf, uint64_t beg, uint64_t end);

/*
 Bit scan forward, finds the first bit set within a given range (beg, end)
 Returns 0 if no bit is found, otherwise it returns the absolute offset to the bit (indexing starts at 1, posix convention)
 To use this in a loop to visit each bit set, the following operation can be used:

 md_bitfield_t* bitfield;
 int64_t beg_idx = 0;       // Marks the beginning of the range to scan
 int64_t end_idx = 200;     // Marks the end of the range to scan
 while ((beg_idx = md_bitfield_scan(bitfield, beg_idx, end_idx)) != 0) {
     int64_t bit_idx = beg_idx - 1;   // beg_idx has a start index of 1 so we have to subtract one to get the 'real' index
     do_something(bit_idx);
 }
*/
uint64_t md_bitfield_scan(const md_bitfield_t* bf, uint64_t beg, uint64_t end);

// Create default iterator full range
md_bitfield_iter_t md_bitfield_iter_create(const md_bitfield_t* bf);

// Create iterator with range
md_bitfield_iter_t md_bitfield_iter_range_create(const md_bitfield_t* bf, uint64_t beg, uint64_t end);

// Go to the next bit set in the bitfield, returns true if found, false if not found
bool md_bitfield_iter_next(md_bitfield_iter_t* iter);

void md_bitfield_iter_skip_to_idx(md_bitfield_iter_t* iter, uint64_t idx);

// Extract the indices which are represented by the iterator
// Writes the indices of bits set to the given buffer.
// cap is the capacity (num elements) of the buffer.
// returns the number of indices written to the buffer
// Use popcount to determine required capacity of buffer.
size_t md_bitfield_iter_extract_indices(int32_t* buf, size_t cap, md_bitfield_iter_t iter);

// Returns the actual index to the bit pointed to by the iterator
static inline uint64_t md_bitfield_iter_idx(const md_bitfield_iter_t* it) {
    return it->idx - 1;
}

// Get the minimal set range
bool md_bitfield_get_range(uint64_t* first_idx, uint64_t* last_idx, const md_bitfield_t* bf);

// Copy the contents of the bitfield into an external buffer
//bool md_bitfield_extract_bits_u64(uint64_t* dst_ptr, int64_t num_bits, const md_bitfield_t* src);

// Returns the maximum serialization size in bytes of a bitfield
size_t md_bitfield_serialize_size_in_bytes(const md_bitfield_t* bf);

// Serializes a bitfield into a destination buffer
// It is expected that the supplied buffer has the size_in_bytes supplied by bitfield_serialize_size_in_bytes()
size_t md_bitfield_serialize(void* dst, const md_bitfield_t* bf);

// Deserializes a compressed buffer into a bitfield.
// User must ensure that the bitfield is properly initialized with an allocator
bool md_bitfield_deserialize(md_bitfield_t* bf, const void* src, size_t num_bytes);

// Compute hash for the bitfield
uint64_t md_bitfield_hash64(const md_bitfield_t* bf, uint64_t seed);

#ifdef __cplusplus
}
#endif
