#pragma once
#include <stdint.h>
#include <stdbool.h>

struct md_allocator_i;

typedef struct md_exp_bitfield_t {
    uint32_t magic;
    uint32_t flags;
    uint32_t beg_bit;
    uint32_t end_bit;
    void* bits;
    struct md_allocator_i* alloc;
} md_exp_bitfield_t;

#ifdef __cplusplus
extern "C" {
#endif

// The bitfield type is a semi-dense bitfield which only needs to hold the range of bits which are set.
// Because there is no explicit number of bits for the bitfield, some operations require explicit range parameters beg, end.

void md_bitfield_init           (md_exp_bitfield_t* bf, struct md_allocator_i* alloc);
bool md_bitfield_free           (md_exp_bitfield_t* bf);

bool md_bitfield_empty          (const md_exp_bitfield_t* bf);

void md_bitfield_set_range      (md_exp_bitfield_t* bf, int64_t beg, int64_t end);
void md_bitfield_set_bit        (md_exp_bitfield_t* bf, int64_t bit_idx);
void md_bitfield_set_indices_u32(md_exp_bitfield_t* bf, uint32_t* indices, int64_t num_indices);

void md_bitfield_clear          (md_exp_bitfield_t* bf);
void md_bitfield_clear_range    (md_exp_bitfield_t* bf, int64_t beg, int64_t end);
void md_bitfield_clear_bit      (md_exp_bitfield_t* bf, int64_t bit_idx);

void md_bitfield_or             (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);
void md_bitfield_or_inplace     (md_exp_bitfield_t* a,   const md_exp_bitfield_t* b);

void md_bitfield_and            (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);
void md_bitfield_and_inplace    (md_exp_bitfield_t* a, const md_exp_bitfield_t* b);

void md_bitfield_andnot         (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);
void md_bitfield_andnot_inplace (md_exp_bitfield_t* a, const md_exp_bitfield_t* b);

void md_bitfield_xor            (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);
void md_bitfield_xor_inplace    (md_exp_bitfield_t* a, const md_exp_bitfield_t* b);

// We need an explicit range or something here to 
void md_bitfield_not            (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src, int64_t beg, int64_t end);
void md_bitfield_not_inplace    (md_exp_bitfield_t* bf, int64_t beg, int64_t end);

void md_bitfield_copy           (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src);

// Count number of bits set
int64_t md_bitfield_popcount    (const md_exp_bitfield_t* bf);
int64_t md_bitfield_popcount_range(const md_exp_bitfield_t* bf, int64_t beg, int64_t end);

// Test if single bit in field is set
bool md_bitfield_test_bit   (const md_exp_bitfield_t* bf, int64_t idx);

/*
 Bit scan forward, finds the first bit set within a given range (beg, end)
 Returns 0 if no bit is found, otherwise it returns the absolute offset to the bit (indexing starts at 1, posix convention)
 To use this in a loop to visit each bit set, the following operation can be used:

 md_exp_bitfield_t* bitfield;
 int64_t beg_idx = 0;       // Marks the beginning of the range to scan
 int64_t end_idx = 200;     // Marks the end of the range to scan
 while ((beg_idx = md_bitfield_scan(bitfield, beg_idx, end_idx) != 0) {
     int64_t bit_idx = beg_idx - 1;   // beg_idx has a start index of 1 so we have to subtract one to get the 'real' index
     do_something(bit_idx);
 }
*/
int64_t md_bitfield_scan(const md_exp_bitfield_t* bf, int64_t beg, int64_t end);

int64_t md_bitfield_scan_reverse(const md_exp_bitfield_t* bf, int64_t beg, int64_t end);

// Copy the contents of the bitfield into an external buffer
//bool md_bitfield_extract_bits_u64(uint64_t* dst_ptr, int64_t num_bits, const md_exp_bitfield_t* src);

// Will create an md_array of uint32_t containing the indices of all set bits.
// The md_array is allocated using the supplied allocator.
// User is responsible for freeing the memory of the returned array (md_array_free).
uint32_t* md_bitfield_extract_indices_u32(const md_exp_bitfield_t* bf, struct md_allocator_i* alloc);
uint32_t* md_bitfield_extract_bits_u32(const md_exp_bitfield_t* bf, struct md_allocator_i* alloc);

// Returns the maximum serialization size in bytes of a bitfield
int64_t md_bitfield_serialize_size_in_bytes(const md_exp_bitfield_t* bf);

// Serializes a bitfield into a destination buffer
// It is expected that the supplied buffer has the size_in_bytes supplied by bitfield_serialize_size_in_bytes()
int64_t md_bitfield_serialize(void* dst, const md_exp_bitfield_t* bf);

// Deserializes a compressed buffer into a bitfield.
// User must ensure that the bitfield is properly initialized with an allocator
bool md_bitfield_deserialize(md_exp_bitfield_t* bf, const void* src, int64_t num_bytes);

#ifdef __cplusplus
}
#endif
