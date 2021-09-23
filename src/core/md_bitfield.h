#pragma once
#include <stdint.h>
#include <stdbool.h>

typedef struct md_bitfield_t {
    uint64_t *bits;
    int64_t  num_bits;
} md_bitfield_t;

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

void md_bitfield_clear          (md_exp_bitfield_t* bf);
void md_bitfield_clear_range    (md_exp_bitfield_t* bf, int64_t beg, int64_t end);
void md_bitfield_clear_bit      (md_exp_bitfield_t* bf, int64_t bit_idx);

void md_bitfield_or             (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);
void md_bitfield_and            (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);

// We need an explicit range or something here to 
void md_bitfield_not            (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src, int64_t beg, int64_t end);

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

// Copy the contents of the bitfield into an external buffer
bool md_bitfield_extract_u64(uint64_t* dst_ptr, int64_t num_bits, const md_exp_bitfield_t* src);

//void md_bitfield_or            (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);
//void md_bitfield_or_range      (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b, int64_t beg, int64_t end);

/*
void md_bitfield_or_not        (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);
void md_bitfield_or_not_range  (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b, int64_t beg, int64_t end);
*/

//void md_bitfield_and           (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);
//void md_bitfield_and_range     (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b, int64_t beg, int64_t end);

/*
void md_bitfield_and_not       (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);
void md_bitfield_and_not_range (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b, int64_t beg, int64_t end);
*/

//void md_bitfield_xor           (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);
//void md_bitfield_xor_range     (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b, int64_t beg, int64_t end);

//void md_bitfield_not           (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src);
//void md_bitfield_not_range     (md_exp_bitfield_t* dst, const md_exp_bitfield_t* src, int64_t beg, int64_t end);


//void md_bitfield_copy_range(md_exp_bitfield_t* dst, const md_exp_bitfield_t* src, int64_t beg, int64_t end);


// Counts the number of bits set
//int64_t md_bitfield_count_range (const md_exp_bitfield_t* bf, int64_t beg, int64_t end);



// Test if bitfields are equivalent
//bool md_bitfield_cmp        (const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b);
//bool md_bitfield_cmp_range  (const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b, int64_t beg, int64_t end);



#ifdef __cplusplus
}
#endif
