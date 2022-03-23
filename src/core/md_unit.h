#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <core/md_str.h>

struct md_allocator_i;

/*

This is a very simplistic unit system to cover our needs.
There is a single unit structure which represents a fixed set of units with the corresponding
dimensions. Some of which are dimensionless, but we need to propagate them through our system
anyways.

Since we use fixed units within the structure, we need to deal with unit conversions, in the case
where one would like for example to use nanometer instead of Ångström.

This is currently exposed with a flag which is supplied to the print functions.
Not the most flexible solution, but gets the job done.

*/

// TODO (Robin): Fill in more units over time and cover the equivalent cases in the convert function.
enum {
    UNIT_LENGTH_ANGSTROM   = 0,
    UNIT_LENGTH_NANOMETER  = 1,
    UNIT_LENGTH_COUNT      = 2,
};

enum {
    UNIT_TIME_PIKOSECONDS  = 0,
    UNIT_TIME_NANOSECONDS  = 1,
    UNIT_TIME_COUNT        = 2,
};

enum {
    UNIT_ANGLE_RADIANS     = 0,
    UNIT_ANGLE_DEGREES     = 1,
    UNIT_ANGLE_COUNT       = 2,
};

typedef union md_unit_base_t {
    struct {
        uint32_t length  : 4;
        uint32_t mass    : 4;
        uint32_t time    : 4;
        uint32_t current : 4;
        uint32_t mole    : 4;                    
        uint32_t temp    : 4;
        uint32_t count   : 4;
        uint32_t angle   : 4;
    };
    uint32_t raw_bits;
} md_unit_base_t;

typedef union md_unit_dim_t {
    struct {
        int32_t length  : 4;
        int32_t mass    : 4;
        int32_t time    : 4;
        int32_t current : 4;
        int32_t mole    : 4;                    
        int32_t temp    : 4;
        int32_t count   : 4;
        int32_t angle   : 4;
    };
    uint32_t raw_bits;
} md_unit_dim_t;

typedef struct md_unit_t {
    md_unit_base_t base;
    md_unit_dim_t  dim;
} md_unit_t;

static inline bool unit_empty(md_unit_t unit) {
    return unit.dim.raw_bits == 0 && unit.base.raw_bits == 0;
}

#ifdef __cplusplus
extern "C" {
#endif

md_unit_t unit_mul(md_unit_t a, md_unit_t b);
md_unit_t unit_div(md_unit_t a, md_unit_t b);
md_unit_t unit_add(md_unit_t a, md_unit_t b);
md_unit_t unit_sub(md_unit_t a, md_unit_t b);

void unit_convert_d(double* values, int64_t num_values, md_unit_t* unit, md_unit_base_t new_base);
void unit_convert_f(float*  values, int64_t num_values, md_unit_t* unit, md_unit_base_t new_base);

int unit_print_long (char* buf, int cap, md_unit_t unit);
int unit_print      (char* buf, int cap, md_unit_t unit);

#ifdef __cplusplus
}
#endif