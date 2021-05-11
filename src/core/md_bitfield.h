#pragma once
#include <stdint.h>

typedef struct md_bitfield {
    uint64_t *bits;
    int64_t  num_bits;
} md_bitfield_t;